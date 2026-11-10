# Small-strain isotropic elasticity on the straight-sided reference mesh.
# All interior degrees of freedom participate, including primary vertices.
# The original curved boundary is the Dirichlet displacement; a candidate
# is published only after the existing whole-element positivity certificates.
using LinearAlgebra: qr, dot
using ..FEMLinearSolvers: _fem_solve_spd
using ..MeshFunctionSpaces: mesh_basis_functions
using ..Elements: msh_dimension

function _elastic_reference(coords, blocks)
    reference = copy(coords)
    primary = falses(size(coords,2))
    for (msh,cells) in blocks
        fam = _p2opt_family(msh)
        for t in axes(cells,2), k in 1:fam.primary
            primary[cells[k,t]] = true
        end
    end
    for (msh,cells) in blocks
        fam = _p2opt_family(msh)
        for t in axes(cells,2)
            for (k,a,b) in fam.edge_slots
                v = cells[k,t]
                primary[v] && continue
                for d in 1:3
                    reference[d,v] = _midpoint_coordinate(
                        coords[d,cells[a,t]],coords[d,cells[b,t]])
                end
            end
            for (k,slots) in fam.face_slots
                v = cells[k,t]
                primary[v] && continue
                for d in 1:3
                    reference[d,v] = sum(coords[d,cells[a,t]]/length(slots) for a in slots)
                end
            end
            for k in fam.interior_slots
                v = cells[k,t]
                for d in 1:3
                    reference[d,v] = sum(coords[d,cells[a,t]]/fam.primary for a in 1:fam.primary)
                end
            end
        end
    end
    return reference
end

function _elastic_boundary_nodes(blocks)
    facets = Dict{Tuple,Tuple{Int,Vector{Int}}}()
    for (msh,cells) in blocks
        fam = _p2opt_family(msh)
        for t in axes(cells,2), facet in fam.facets
            key = Tuple(sort!([Int(cells[a,t]) for a in facet]))
            nodes = [Int(cells[a,t]) for a in facet]
            for (k,a,b) in fam.edge_slots
                a in facet && b in facet && push!(nodes,Int(cells[k,t]))
            end
            for (k,slots) in fam.face_slots
                all(a -> a in facet,slots) && push!(nodes,Int(cells[k,t]))
            end
            previous = get(facets,key,(0,Int[]))
            facets[key] = (previous[1]+1,nodes)
        end
    end
    boundary = Set{Int}()
    for (count,nodes) in values(facets)
        count == 1 && union!(boundary,nodes)
    end
    return boundary
end

function _elastic_component_root!(parent::Dict{Int,Int}, node::Int)
    root = node
    while parent[root] != root
        root = parent[root]
    end
    while node != root
        following = parent[node]
        parent[node] = root
        node = following
    end
    return root
end

# Disconnected domains have independent linear systems. A local coordinate and
# displacement scale for each prevents a distant large component from erasing
# small finite geometry or underflowing its Dirichlet load.
function _elastic_components(domain)
    parent = Dict{Int,Int}()
    for (_,cells) in domain, t in axes(cells,2)
        firstnode = Int(cells[1,t])
        get!(parent,firstnode,firstnode)
        for node in Base.view(cells,:,t)
            v = Int(node)
            get!(parent,v,v)
            a = _elastic_component_root!(parent,firstnode)
            b = _elastic_component_root!(parent,v)
            parent[b] = a
        end
    end
    groups = Dict{Int,Dict{Int,Vector{Int}}}()
    for (b,(_,cells)) in enumerate(domain), t in axes(cells,2)
        root = _elastic_component_root!(parent,Int(cells[1,t]))
        byblock = get!(groups,root,Dict{Int,Vector{Int}}())
        push!(get!(byblock,b,Int[]),t)
    end
    length(groups) == 1 && return [domain]
    return [[(domain[b][1],domain[b][2][:,byblock[b]])
             for b in sort!(collect(keys(byblock)))]
            for (_,byblock) in sort!(collect(groups);by=first)]
end

function elastic_p2!(coords, allblocks; poisson::Real=0.33,
        ownership=nothing, project=nothing, itermax::Int=100,
        passes::Int=25)
    isfinite(poisson) && -1 < poisson < 0.5 ||
        throw(ArgumentError("High-order elastic Poisson ratio must lie in (-1, 0.5)"))
    itermax >= 0 && passes >= 0 || throw(ArgumentError("negative elastic iteration limit"))
    allblocks = collect(allblocks)
    original = copy(coords)
    blocks = [(Int32(msh),cells) for (msh,cells) in allblocks
              if _p2opt_family(msh) !== nothing && size(cells,2)>0]
    isempty(blocks) && return ((0,NaN,NaN,0,0),Dict{Int,Vector{Float64}}())
    before = minimum(_p2opt_elem_jac(coords,_p2opt_family(msh),cells,t)
                     for (msh,cells) in blocks for t in axes(cells,2))
    (itermax == 0 || passes == 0) &&
        return ((0,before,before,0,0),Dict{Int,Vector{Float64}}())
    reference = _elastic_reference(coords,blocks)
    used = sort!(unique!([Int(v) for (_,cells) in blocks for v in cells]))
    target = copy(original)
    normalized = zeros(size(reference))
    displacement = zeros(size(reference))
    mu = 1/(2*(1+Float64(poisson)))
    lambda = Float64(poisson)/((1+Float64(poisson))*(1-2*Float64(poisson)))
    higher_nodes = Set{Int}()
    accepted_parameters = Dict{Int,Vector{Float64}}()
    for dimension in (3,2)
        domain = [(msh,cells) for (msh,cells) in blocks
                  if (_p2opt_family(msh).volume ? 3 : 2) == dimension]
        isempty(domain) && continue
        active = Set{Int}(Int(v) for (_,cells) in domain for v in cells)
        if dimension == 2
            # Surface blocks attached to a volume are its Dirichlet boundary.
            domain = [(msh,cells[:,[all(v -> !(Int(v) in higher_nodes),cells[:,t])
                       for t in axes(cells,2)]]) for (msh,cells) in domain]
            filter!(block -> size(block[2],2)>0,domain)
            isempty(domain) && continue
            active = Set{Int}(Int(v) for (_,cells) in domain for v in cells)
        else
            union!(higher_nodes,active)
        end
        for domain in _elastic_components(domain)
            component_nodes = sort!(unique!([Int(v) for (_,cells) in domain for v in cells]))
            active = Set(component_nodes)
            origin = copy(reference[:,first(component_nodes)])
            scale = maximum(abs(reference[d,v]-origin[d]) for v in component_nodes for d in 1:3)
            isfinite(scale) && scale > 0 ||
                throw(ArgumentError("invalid elastic reference extent"))
            for v in component_nodes, d in 1:3
                normalized[d,v] = (reference[d,v]-origin[d])/scale
                displacement[d,v] = (original[d,v]-reference[d,v])/scale
            end
            boundary = _elastic_boundary_nodes(domain)
            # A node also used by an unoptimized family cannot move without a
            # positivity certificate for that other element.
            foreign = Set{Int}()
            for (msh,cells) in allblocks
                if _p2opt_family(msh) === nothing && msh_dimension(msh) >= dimension
                    union!(foreign,Int(v) for v in cells)
                end
            end
            unknown = sort!([v for v in active if !(v in boundary) && !(v in foreign) &&
                (ownership === nothing || ownership[v] == dimension)])
            isempty(unknown) && continue
            indices = Dict(v=>i for (i,v) in enumerate(unknown))
            ndofs = 3length(unknown)
            rows=Int[]; columns=Int[]; values=Float64[]; rhs=zeros(ndofs)
            for (msh,cells) in domain
                nn = size(cells,1)
                points,weights = mesh_integration_points(msh,"Gauss6")
                _,flat,_ = mesh_basis_functions(msh,points,"GradLagrange")
                gradients = reshape(flat,3,nn,:)
                localmatrix=zeros(3nn,3nn)
                for t in axes(cells,2)
                    fill!(localmatrix,0.0)
                    conn = cells[:,t]
                    for q in eachindex(weights)
                        g = gradients[1:dimension,:,q]
                        jac = normalized[:,conn]*transpose(g)
                        factorization = qr(jac)
                        triangular = Matrix(factorization.R)
                        diagonal = [abs(triangular[d,d]) for d in 1:dimension]
                        all(v -> isfinite(v) && v > 0,diagonal) ||
                            throw(ArgumentError("degenerate elastic reference element"))
                        # Weight the gradients before division. Forming J'J squares
                        # the aspect ratio, while multiplying detJ by squared
                        # unweighted gradients can underflow/overflow separately.
                        root_measure = prod(sqrt,diagonal)
                        isfinite(root_measure) && root_measure > 0 ||
                            throw(ArgumentError("unrepresentable elastic reference measure"))
                        orthogonal = Matrix(factorization.Q)[:,1:dimension]
                        world = orthogonal*(transpose(triangular)\(g.*root_measure))
                        weight = weights[q]
                        for i in 1:nn, j in 1:nn
                            gg=dot(Base.view(world,:,i),Base.view(world,:,j))
                            for a in 1:3, b in 1:3
                                localmatrix[3(i-1)+a,3(j-1)+b] += weight*(
                                    (a == b ? mu*gg : 0.0) +
                                    mu*world[b,i]*world[a,j] +
                                    lambda*world[a,i]*world[b,j])
                            end
                        end
                    end
                    # Use one rounded coefficient for each symmetric pair; independent
                    # evaluation of cancelling Lamé terms can differ near zero.
                    for a in axes(localmatrix,1), b in a+1:size(localmatrix,2)
                        value = _midpoint_coordinate(localmatrix[a,b],localmatrix[b,a])
                        localmatrix[a,b] = value
                        localmatrix[b,a] = value
                    end
                    for i in 1:nn
                        ii=get(indices,Int(conn[i]),0)
                        ii == 0 && continue
                        for j in 1:nn
                            jj=get(indices,Int(conn[j]),0)
                            for a in 1:3, b in 1:3
                                entry=localmatrix[3(i-1)+a,3(j-1)+b]
                                iszero(entry) && continue
                                row=3(ii-1)+a
                                if jj == 0
                                    rhs[row] -= entry*displacement[b,conn[j]]
                                else
                                    push!(rows,row);push!(columns,3(jj-1)+b);push!(values,entry)
                                end
                            end
                        end
                    end
                end
            end
            solution = _fem_solve_spd(ndofs,rows,columns,values,rhs)
            for (i,v) in enumerate(unknown), d in 1:3
                target[d,v] = reference[d,v]+scale*solution[3(i-1)+d]
            end
        end
    end
    all(isfinite,target) || throw(ArgumentError("nonfinite elastic target"))
    trials = 0
    candidate = copy(original)
    for step in 0:24
        trials += 1
        factor = 2.0^(-step)
        empty!(accepted_parameters)
        copyto!(candidate,original)
        for v in used
            all(d -> target[d,v] == original[d,v],1:3) && continue
            point = [original[d,v]+factor*(target[d,v]-original[d,v]) for d in 1:3]
            if project !== nothing
                projected = project(v,point)
                projected === nothing && continue
                point,parameter = projected
                parameter === nothing || (accepted_parameters[v]=parameter)
            end
            candidate[:,v] .= point
        end
        if all(_p2opt_certify(_p2opt_family(msh).kind,candidate,cells,t)
               for (msh,cells) in blocks for t in axes(cells,2))
            copyto!(coords,candidate)
            moved = count(v -> any(d -> coords[d,v] != original[d,v],1:3),used)
            after = minimum(_p2opt_elem_jac(coords,_p2opt_family(msh),cells,t)
                            for (msh,cells) in blocks for t in axes(cells,2))
            return ((moved,before,after,1,trials),accepted_parameters)
        end
    end
    copyto!(coords,original)
    return ((0,before,before,1,trials),Dict{Int,Vector{Float64}}())
end
