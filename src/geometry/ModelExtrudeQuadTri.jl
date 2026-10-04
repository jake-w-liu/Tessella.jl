# Gmsh 4.15.2 QuadTriExtruded2D / QuadTriExtruded3D: the AddVerts
# transition preserves all unconstrained prism/hex cells and fans each
# constrained cell to its unique-corner centroid. Surface adjacency decides
# whether a lateral quadrangle must stay whole or be divided.

function _extrude_quadtri_regions(m::GeoModel,t::Int,caller::AbstractString)
    regions=Int[]
    for tag in keys(m.volumes)
        for shell in m.volumes[tag]
            if any(s->abs(s)==t,m.surface_loops[abs(shell)])
                push!(regions,tag)
                break
            end
        end
    end
    sort!(regions)
    return regions
end

@inline _extrude_is_quadtri(ep)=ep!==nothing && !isempty(ep.layers) &&
    ep.recombine && ep.quad_to_tri!==:none

# IsInToroidalQuadToTri: copied section faces keep their source quadrangles
# when a chain of structured sweeps closes onto its original source face.
# Follow entity links, rather than the rotation angle: closure is established
# by geometry coherence and also applies to a cycle containing ordinary sweeps.
function _extrude_quadtri_source_root(m::GeoModel,t::Int)
    seen=Set{Int}()
    while !(t in seen)
        push!(seen,t)
        ep=_extrude_gate(m,2,t)
        link=get(m.meshing.extrude_sources,(2,t),nothing)
        (ep===nothing || link===nothing || link[1]!=2) && return t
        t=abs(link[2])
    end
    throw(ArgumentError("QuadToTri toroidal cap has a cyclic source dependency"))
end

function _extrude_quadtri_toroidal_root(m::GeoModel,t::Int,
                                        caller::AbstractString)
    root=_extrude_quadtri_source_root(m,t)
    adjacent=Int[]
    for r in _extrude_quadtri_regions(m,root,caller)
        _extrude_gate(m,3,r)===nothing || push!(adjacent,r)
    end
    length(adjacent)==2 || return 0
    first=false
    last=0
    quadtri=false
    for r in adjacent
        ep=_extrude_gate(m,3,r)
        link=get(m.meshing.extrude_sources,(3,r),nothing)
        link!==nothing && link[1]==2 || return 0
        quadtri|=ep.quad_to_tri!==:none
        if abs(link[2])==root
            first=true
        else
            last=abs(link[2])
        end
    end
    first && last!=0 || return 0
    seen=Set{Int}()
    while !(last in seen)
        push!(seen,last)
        ep=_extrude_gate(m,2,last)
        link=get(m.meshing.extrude_sources,(2,last),nothing)
        ep!==nothing && link!==nothing && link[1]==2 || return 0
        quadtri|=ep.quad_to_tri!==:none
        last=abs(link[2])
        last==root && return quadtri ? root : 0
    end
    throw(ArgumentError("QuadToTri toroidal cap has a cyclic source dependency"))
end

function _extrude_quadtri_closing_root(m::GeoModel,r::Int,
                                       caller::AbstractString)
    link=get(m.meshing.extrude_sources,(3,r),nothing)
    link!==nothing && link[1]==2 || return 0
    src=abs(link[2])
    root=_extrude_quadtri_toroidal_root(m,src,caller)
    root!=0 && root!=src && r in _extrude_quadtri_regions(m,root,caller) || return 0
    return root
end

# Upstream's vertex rtree resolves the closing row to the existing source
# mesh. A spatial hash provides the same coincidence lookup without an O(N²)
# nearest-point pass, and keeps the seam identical in volume and lateral parts.
function _extrude_quadtri_seam_index(m::GeoModel,r::Int,
                                     caller::AbstractString)
    root=_extrude_quadtri_closing_root(m,r,caller)
    root==0 && return nothing
    mesh=mesh_model_surface(m,root)
    tolerance=_coherence_eps(m)
    origin=nnodes(mesh)==0 ? (0.0,0.0,0.0) :
        (mesh.coords[1,1],mesh.coords[2,1],mesh.coords[3,1])
    grid=Dict{NTuple{3,Int},Vector{NTuple{3,Float64}}}()
    for i in 1:nnodes(mesh)
        p=(mesh.coords[1,i],mesh.coords[2,i],mesh.coords[3,i])
        cell=_coherence_spatial_cell(p,origin,tolerance)
        push!(get!(grid,cell,NTuple{3,Float64}[]),p)
    end
    return (grid,tolerance,origin)
end

function _extrude_quadtri_snap_seam(index,p::NTuple{3,Float64},
                                    caller::AbstractString)
    index===nothing && return p
    grid,tolerance,origin=index
    cell=_coherence_spatial_cell(p,origin,tolerance)
    for dx in -1:1,dy in -1:1,dz in -1:1
        bucket=get(grid,(cell[1]+dx,cell[2]+dy,cell[3]+dz),nothing)
        bucket===nothing && continue
        for other in bucket
            _points_close(p,other,tolerance) && return other
        end
    end
    throw(ArgumentError("$caller: closed QuadToTri extrusion does not match its source mesh"))
end

function _extrude_quadtri_lateral_recombine(m::GeoModel,t::Int,params,
                                            caller::AbstractString)
    adjacent=_extrude_quadtri_regions(m,t,caller)
    qt=0
    for r in adjacent
        ep=_extrude_gate(m,3,r)
        _extrude_is_quadtri(ep) || continue
        source=get(m.meshing.extrude_sources,(3,r),nothing)
        source===nothing && continue
        abs(source[2])==t && continue
        link=get(m.meshing.extrude_sources,(2,t),nothing)
        link!==nothing && link[1]==1 && (qt=r;break)
    end
    qt==0 && return params.recombine
    length(adjacent)<=2 || throw(ArgumentError(
        "$caller: too many regions adjacent to QuadToTri Surface[$t]"))
    ep=_extrude_gate(m,3,qt)
    length(adjacent)==1 && return ep.recomb_laterals
    neighbor=adjacent[1]==qt ? adjacent[2] : adjacent[1]
    haskey(m.meshing.transfinite_volumes,neighbor) && return params.recombine
    other=_extrude_gate(m,3,neighbor)
    other===nothing && return false
    source=get(m.meshing.extrude_sources,(3,neighbor),nothing)
    source!==nothing && abs(source[2])==t && return ep.recomb_laterals
    return other.quad_to_tri===:none ? other.recombine : true
end

function _extrude_quadtri_top(m::GeoModel,t::Int,params,
                              caller::AbstractString)
    _extrude_is_quadtri(params) || return false
    _extrude_quadtri_toroidal_root(m,t,caller)!=0 && return false
    for r in _extrude_quadtri_regions(m,t,caller)
        ep=_extrude_gate(m,3,r)
        _extrude_is_quadtri(ep) || continue
        source=get(m.meshing.extrude_sources,(3,r),nothing)
        link=get(m.meshing.extrude_sources,(2,t),nothing)
        source!==nothing && link!==nothing && link[1]==2 &&
            abs(source[2])==abs(link[2]) && return true
    end
    return false
end

# A copied/transfinite root uses its fixed (1,3) diagonal. The unstructured
# root uses MeshQuadToTriTopUnstructured's squared-distance criterion;
# upstream deliberately truncates side lengths to C++ int before squaring.
function _extrude_quadtri_structured_root(m::GeoModel,src::Int)
    seen=Set{Int}()
    while !(src in seen)
        push!(seen,src)
        haskey(m.meshing.transfinite_surfaces,src) && return true
        ep=_extrude_gate(m,2,src)
        ep===nothing && return false
        link=get(m.meshing.extrude_sources,(2,src),nothing)
        link===nothing && return false
        link[1]==1 && return true
        src=abs(link[2])
    end
    throw(ArgumentError("QuadToTri top has a cyclic source dependency"))
end

@inline function _extrude_quadtri_top_diagonal(coords,cell,structured::Bool)
    structured && return true
    average=0.0
    for k in 1:4
        a=coords[cell[k]]; b=coords[cell[mod1(k+1,4)]]
        d=trunc(sqrt((a[1]-b[1])^2+(a[2]-b[2])^2+(a[3]-b[3])^2))
        average+=d*d
    end
    average/=4
    a,b,c,d=(coords[cell[k]] for k in 1:4)
    d1=(a[1]-c[1])^2+(a[2]-c[2])^2+(a[3]-c[3])^2
    d2=(b[1]-d[1])^2+(b[2]-d[2])^2+(b[3]-d[3])^2
    return abs(d1-average)<=abs(d2-average)
end

function _extrude_quadtri_top_part(m::GeoModel,src::Int,srcmesh,
                                   coords,remap,caller::AbstractString;
                                   quad_to_tri::Symbol=:add_verts)
    tris=NTuple{3,Int32}[]
    structured=_extrude_quadtri_structured_root(m,src)
    for block in srcmesh.blocks
        block isa ElementBlock || continue
        if block.msh==2
            for i in axes(block.nodes,2)
                push!(tris,ntuple(k->remap[block.nodes[k,i]],3))
            end
        elseif block.msh==3
            for i in axes(block.nodes,2)
                cell=ntuple(k->remap[block.nodes[k,i]],4)
                a,b,c,d=cell
                if quad_to_tri===:no_new_verts
                    # NoNew uses the original source vertex precedence. The
                    # output remap and AddVerts distance policy do not select
                    # this cap's diagonal.
                    original=ntuple(k->block.nodes[k,i],4)
                    low=1
                    for k in 2:4
                        original[k]<original[low] && (low=k)
                    end
                    p=cell[low];q=cell[mod1(low+1,4)]
                    r=cell[mod1(low+2,4)];s=cell[mod1(low+3,4)]
                    push!(tris,(p,q,r));push!(tris,(p,r,s))
                elseif _extrude_quadtri_top_diagonal(coords,cell,structured)
                    # Structured roots use upstream's reverse pair ordering.
                    if structured
                        push!(tris,(a,c,d)); push!(tris,(a,b,c))
                    else
                        push!(tris,(a,b,c)); push!(tris,(a,c,d))
                    end
                else
                    push!(tris,(b,c,d)); push!(tris,(b,d,a))
                end
            end
        end
    end
    nodes=Matrix{Int32}(undef,3,length(tris))
    for (i,tri) in enumerate(tris),k in 1:3
        nodes[k,i]=tri[k]
    end
    out=Matrix{Float64}(undef,3,length(coords))
    for (i,p) in enumerate(coords),k in 1:3
        out[k,i]=p[k]
    end
    return Mesh(out;tris=nodes)
end

function _extrude_quadtri_boundary_edges(m::GeoModel,t::Int,
                                         caller::AbstractString;edges=nothing)
    diagonals=Set{NTuple{2,NTuple{3,Float64}}}()
    for signed in _model_volume_boundary_surfaces(m,t,caller)
        surface=mesh_model_surface(m,abs(signed);_extrude_edges=edges)
        coords=surface.coords
        if surface isa Mesh
            cells=surface.tris
            for i in axes(cells,2),k in 1:3
                a=cells[k,i];b=cells[mod1(k+1,3),i]
                push!(diagonals,_extrude_ekey(
                    (coords[1,a],coords[2,a],coords[3,a]),
                    (coords[1,b],coords[2,b],coords[3,b])))
            end
        else
            for block in surface.blocks
                block isa ElementBlock && block.msh==2 || continue
                cells=block.nodes
                for i in axes(cells,2),k in 1:3
                    a=cells[k,i];b=cells[mod1(k+1,3),i]
                    push!(diagonals,_extrude_ekey(
                        (coords[1,a],coords[2,a],coords[3,a]),
                        (coords[1,b],coords[2,b],coords[3,b])))
                end
            end
        end
    end
    return diagonals
end

# QtFindVertsCentroid sums in the original corner order, omitting the top
# copy of each fixed column. This order matters for bitwise coordinate parity.
@inline function _extrude_quadtri_centroid(v::NTuple{N,NTuple{3,Float64}}) where N
    n=N÷2
    x=0.0;y=0.0;z=0.0;count=0
    for k in 1:N
        k>n && v[k]==v[k-n] && continue
        p=v[k];x+=p[1];y+=p[2];z+=p[3];count+=1
    end
    return (x/count,y/count,z/count)
end

const _EXTRUDE_QT_HEX_FACES=((1,4,3,2),(5,6,7,8),(1,2,6,5),
    (2,3,7,6),(3,4,8,7),(4,1,5,8))
const _EXTRUDE_QT_PRISM_FACES=((1,3,2),(4,5,6),(1,2,5,4),
    (2,3,6,5),(3,1,4,6))

@inline function _extrude_quadtri_certify_triangle(a,b,c,center,orientation)
    orient3(a,b,c,center)*orientation>0 || throw(ArgumentError(
        "QuadTriAddVerts: centroid fan has a folded or degenerate face"))
    volume=tet_signed_volume(a,c,b,center)*orientation
    isfinite(volume) && volume>0 || throw(ArgumentError(
        "QuadTriAddVerts: centroid tetrahedron volume is not finite and positive"))
    return nothing
end

function _extrude_quadtri_certify(v::NTuple{N,NTuple{3,Float64}},center,
                                  edges) where N
    msh=N==6 ? 6 : 5
    total=0.0
    for (a,b,c,d) in _EXTRUDE_CELL_TETS[msh]
        total+=tet_signed_volume(v[a],v[b],v[c],v[d])
    end
    isfinite(total) && total!=0 || throw(ArgumentError(
        "QuadTriAddVerts: degenerate or nonfinite logical cell volume"))
    orientation=total>0 ? 1 : -1
    faces=N==6 ? _EXTRUDE_QT_PRISM_FACES : _EXTRUDE_QT_HEX_FACES
    for face in faces
        if length(face)==3
            a,b,c=v[face[1]],v[face[2]],v[face[3]]
            _extrude_quadtri_certify_triangle(a,b,c,center,orientation)
            continue
        end
        a,b,c,d=(v[face[k]] for k in 1:4)
        if a==b || b==c || c==d || d==a
            if a==b
                b,c=c,d
            elseif b==c
                c=d
            end
            (a==b || a==c || b==c) && continue
            _extrude_quadtri_certify_triangle(a,b,c,center,orientation)
            continue
        end
        first=_extrude_ein(edges,a,c)
        second=_extrude_ein(edges,b,d)
        first && second && throw(ArgumentError(
            "QuadTriAddVerts: boundary quadrangle has conflicting diagonals"))
        if first || !second
            _extrude_quadtri_certify_triangle(a,b,c,center,orientation)
            _extrude_quadtri_certify_triangle(a,c,d,center,orientation)
        end
        if second || !first
            _extrude_quadtri_certify_triangle(a,b,d,center,orientation)
            _extrude_quadtri_certify_triangle(b,c,d,center,orientation)
        end
    end
    return nothing
end

function _extrude_quadtri_certify_retained(sweep::_ExtrudeVolumeSweep,
        ntets::Int,nhexes::Int,nprisms::Int,npyramids::Int)
    # The factory may replace collapsed logical cells with a different
    # element family. Certify the actual retained maps, as the fan branch
    # already does for its emitted tetrahedra and pyramids.
    caller="QuadTriAddVerts retained cell"
    try
        for index in ntets+1:length(sweep.tets)
            v=sweep.tets[index]
            volume=tet_signed_volume(v...)
            all(p->all(isfinite,p),v) && -orient3(v...)!=0 &&
                isfinite(volume) && volume!=0 ||
                throw(ArgumentError("$caller: degenerate or nonfinite Tet4 map"))
        end
        for index in nhexes+1:length(sweep.hexes)
            v=sweep.hexes[index]
            corner=-orient3(v[1],v[2],v[4],v[5])
            orientation=corner>0 ? 1 : corner<0 ? -1 : 0
            _extrude_nonew_hex_jacobian_certify(v,orientation,caller)
        end
        for index in nprisms+1:length(sweep.prisms)
            _extrude_nonew_prism_jacobian_certify(sweep.prisms[index],caller)
        end
        for index in npyramids+1:length(sweep.pyramids)
            _extrude_nonew_pyramid_jacobian_certify(sweep.pyramids[index],caller)
        end
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError(replace(err.msg,"QuadTriNoNewVerts"=>"QuadTriAddVerts")))
    end
    return nothing
end

function _extrude_quadtri_fan!(sweep::_ExtrudeVolumeSweep,
                               v::NTuple{N,NTuple{3,Float64}},edges) where N
    n=N÷2
    # A one-column collapsed hex has seven distinct corners and always
    # requires a body fan (QuadToTriHexPriPyrTet's m==1 branch).
    divided=n==4 && count(k->v[k]==v[k+n],1:n)==1
    for k in 1:n
        j=mod1(k+1,n)
        (v[k]==v[k+n] || v[j]==v[j+n]) && continue
        if _extrude_ein(edges,v[k],v[j+n]) ||
           _extrude_ein(edges,v[k+n],v[j])
            divided=true;break
        end
    end
    if n==4 && !divided
        divided=_extrude_ein(edges,v[1],v[3]) ||
                _extrude_ein(edges,v[2],v[4]) ||
                _extrude_ein(edges,v[5],v[7]) ||
                _extrude_ein(edges,v[6],v[8])
    end
    if !divided
        ntets=length(sweep.tets);nhexes=length(sweep.hexes)
        nprisms=length(sweep.prisms);npyramids=length(sweep.pyramids)
        if n==3
            _extrude_pripyrtet!(sweep.tets,sweep.pyramids,sweep.prisms,v)
        else
            _extrude_hexpri!(sweep.hexes,sweep.prisms,sweep.degenerate,v)
        end
        _extrude_quadtri_certify_retained(sweep,ntets,nhexes,nprisms,npyramids)
        return nothing
    end
    center=_extrude_quadtri_centroid(v)
    _extrude_quadtri_certify(v,center,edges)
    for k in 1:n
        j=mod1(k+1,n)
        if v[k]==v[k+n] && v[j]==v[j+n]
            continue
        elseif v[k]==v[k+n] || v[j]==v[j+n]
            fixed=v[k]==v[k+n] ? k : j
            other=fixed==k ? j : k
            _extrude_tet!(sweep.tets,v[fixed],v[other],v[other+n],center)
        elseif _extrude_ein(edges,v[k],v[j+n])
            _extrude_tet!(sweep.tets,v[k],v[j],v[j+n],center)
            _extrude_tet!(sweep.tets,v[k],v[j+n],v[k+n],center)
        elseif _extrude_ein(edges,v[k+n],v[j])
            _extrude_tet!(sweep.tets,v[k],v[j],v[k+n],center)
            _extrude_tet!(sweep.tets,v[j],v[j+n],v[k+n],center)
        else
            push!(sweep.pyramids,(v[k],v[j],v[j+n],v[k+n],center))
        end
    end
    if n==3
        _extrude_tet!(sweep.tets,v[1],v[2],v[3],center)
        _extrude_tet!(sweep.tets,v[4],v[6],v[5],center)
    else
        for add in (0,4)
            a,b,c,d=v[1+add],v[2+add],v[3+add],v[4+add]
            if _extrude_ein(edges,a,c)
                _extrude_tet!(sweep.tets,a,b,c,center)
                _extrude_tet!(sweep.tets,a,c,d,center)
            elseif _extrude_ein(edges,b,d)
                _extrude_tet!(sweep.tets,b,c,d,center)
                _extrude_tet!(sweep.tets,b,d,a,center)
            else
                push!(sweep.pyramids,(a,b,c,d,center))
            end
        end
    end
    return nothing
end

function _extrude_quadtri_addverts!(m::GeoModel,sweep::_ExtrudeVolumeSweep,
                                    caller::AbstractString;edges=nothing)
    boundary=_extrude_quadtri_boundary_edges(m,sweep.tag,caller;edges=edges)
    hexes=sweep.hexes;prisms=sweep.prisms
    sweep.hexes=eltype(hexes)[];sweep.prisms=eltype(prisms)[]
    for cell in prisms
        _extrude_quadtri_fan!(sweep,cell,boundary)
    end
    levels=size(sweep.cols,1)-1
    for (i,cell) in enumerate(hexes)
        # QtAddInteriorVerts preallocates a centroid for the final layer of
        # every quadrilateral column, including unchanged toroidal hexes.
        # Such nodes retain volume ownership even if no emitted cell uses them.
        i%levels==0 && push!(sweep.interior,_extrude_quadtri_centroid(cell))
        _extrude_quadtri_fan!(sweep,cell,boundary)
    end
    return nothing
end
