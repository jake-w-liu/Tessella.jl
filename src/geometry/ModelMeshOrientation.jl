# Detached GRegion/orientMeshGFace orientation preparation. The caller owns
# Surface STL lifetime, actual cell/row bindings, CAD evaluation and publication.
# A cell key is a private operation identity, never a public element label.
import Random

const _ModelMeshOrientationPoint = NTuple{3,Float64}
const _ModelMeshOrientationTriangle = NTuple{3,_ModelMeshOrientationPoint}

function _model_orientation_point(value,caller)
    (value isa Tuple || value isa AbstractVector) && length(value)==3 ||
        throw(ArgumentError("$caller: a point or normal must have three entries"))
    for component in value
        component isa Real && !(component isa Bool) && isfinite(component) ||
            throw(ArgumentError("$caller: point and normal entries must be finite real numbers"))
    end
    point=(Float64(value[1]),Float64(value[2]),Float64(value[3]))
    all(isfinite,point) || throw(ArgumentError(
        "$caller: point and normal entries must fit Float64"))
    return point
end

function _model_orientation_permutation(msh::Integer,caller)
    msh isa Bool && throw(ArgumentError("$caller: MSH type must be an integer"))
    msh==2 && return (1,3,2)
    msh==3 && return (1,4,3,2)
    msh==9 && return (1,3,2,6,5,4)
    msh==10 && return (1,4,3,2,8,7,6,5,9)
    msh==16 && return (1,4,3,2,8,7,6,5)
    throw(ArgumentError("$caller: surface orientation does not support MSH type $msh"))
end

struct _ModelMeshOrientationCell
    key::Int
    msh::Int
    nodes::Vector{Int32}
    vertex_normal::Union{Nothing,_ModelMeshOrientationPoint}
    barycenter_normal::Union{Nothing,_ModelMeshOrientationPoint}
    boundary_layer::Bool
end

function _ModelMeshOrientationCell(key::Integer,msh::Integer,nodes;
        vertex_normal=nothing,barycenter_normal=nothing,boundary_layer::Bool=false)
    caller="Model mesh orientation cell"
    !(key isa Bool) && 1<=key<=typemax(Int) || throw(ArgumentError(
        "$caller: private cell key must be a positive Int"))
    permutation=_model_orientation_permutation(msh,caller)
    (nodes isa AbstractVector || nodes isa Tuple) &&
        length(nodes)==length(permutation) || throw(ArgumentError(
            "$caller: MSH type $msh expects $(length(permutation)) node rows"))
    rows=Int32[]
    sizehint!(rows,length(permutation))
    for node in nodes
        node isa Integer && !(node isa Bool) && 1<=node<=typemax(Int32) ||
            throw(ArgumentError("$caller: node rows must be positive Int32 indices"))
        push!(rows,Int32(node))
    end
    normal1=vertex_normal===nothing ? nothing :
        _model_orientation_point(vertex_normal,caller)
    normal2=barycenter_normal===nothing ? nothing :
        _model_orientation_point(barycenter_normal,caller)
    return _ModelMeshOrientationCell(Int(key),Int(msh),rows,normal1,normal2,boundary_layer)
end

struct _ModelMeshOrientationSurface
    tag::Int
    coords::Matrix{Float64}
    cells::Vector{_ModelMeshOrientationCell}
    stl::Vector{_ModelMeshOrientationTriangle}
    mode::Symbol
end

function _ModelMeshOrientationSurface(tag::Integer,coords,cells,stl;
        mode::Symbol=:cad)
    caller="Model mesh orientation Surface"
    !(tag isa Bool) && 1<=tag<=typemax(Int) || throw(ArgumentError(
        "$caller: Surface tag must be a positive Int"))
    coords isa AbstractMatrix && size(coords,1)==3 || throw(ArgumentError(
        "$caller: coordinates must have three rows"))
    mode in (:cad,:periodic,:fully_discrete,:boundary_layer_surface) ||
        throw(ArgumentError("$caller: unknown orientation mode $mode"))
    cells isa AbstractVector{_ModelMeshOrientationCell} || throw(ArgumentError(
        "$caller: cells must carry explicit private identities and CAD-normal bindings"))
    for cell in cells, slot in 1:3
        node=cell.nodes[slot]
        node<=size(coords,2) || throw(ArgumentError(
            "$caller: private cell $(cell.key) references unknown coordinate row $node"))
        _model_orientation_point((coords[1,node],coords[2,node],coords[3,node]),caller)
    end
    coordinates=coords isa Matrix{Float64} ? coords : Matrix{Float64}(coords)
    for cell in cells, node in cell.nodes
        node<=size(coordinates,2) || throw(ArgumentError(
            "$caller: private cell $(cell.key) references unknown coordinate row $node"))
    end
    stl isa AbstractVector || throw(ArgumentError(
        "$caller: STL must be an explicit vector of triangle points"))
    triangles=_ModelMeshOrientationTriangle[]
    sizehint!(triangles,length(stl))
    for triangle in stl
        (triangle isa Tuple || triangle isa AbstractVector) && length(triangle)==3 ||
            throw(ArgumentError("$caller: each STL triangle must have three points"))
        push!(triangles,(_model_orientation_point(triangle[1],caller),
                         _model_orientation_point(triangle[2],caller),
                         _model_orientation_point(triangle[3],caller)))
    end
    return _ModelMeshOrientationSurface(Int(tag),coordinates,copy(cells),triangles,mode)
end

@inline _model_orientation_sub(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _model_orientation_dot(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
@inline _model_orientation_cross(a,b)=
    (a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])

function _model_orientation_unit(vector)
    length=sqrt(_model_orientation_dot(vector,vector))
    length==0 && return (0.0,0.0,0.0)
    return (vector[1]/length,vector[2]/length,vector[3]/length)
end

function _model_orientation_scaled(triangle,scale)
    p1,p2,p3=triangle
    return ((p1[1]/scale,p1[2]/scale,p1[3]/scale),
            (p2[1]/scale,p2[2]/scale,p2[3]/scale),
            (p3[1]/scale,p3[2]/scale,p3[3]/scale))
end

# Numeric.cpp sys3x3_with_tol and GRegion.cpp intersectLineTriangle use the
# normalized edge columns and +N column, including the sign of the solved t.
function _model_orientation_intersection(triangle,point,normal)
    a,b,c=triangle
    e1=_model_orientation_sub(b,a);e2=_model_orientation_sub(c,a)
    rhs=_model_orientation_sub(point,a)
    m00,m10,m20=e1;m01,m11,m21=e2;m02,m12,m22=normal
    determinant=m00*(m11*m22-m12*m21)-m01*(m10*m22-m12*m20)+
                m02*(m10*m21-m11*m20)
    norm=m00*m00+m01*m01+m02*m02+m10*m10+m11*m11+m12*m12+
         m20*m20+m21*m21+m22*m22
    (determinant==0 || norm==0 || abs(determinant)/norm<1e-12) && return 0
    inverse=1/determinant
    u=(rhs[1]*(m11*m22-m12*m21)-m01*(rhs[2]*m22-m12*rhs[3])+
       m02*(rhs[2]*m21-m11*rhs[3]))*inverse
    v=(m00*(rhs[2]*m22-m12*rhs[3])-rhs[1]*(m10*m22-m12*m20)+
       m02*(m10*rhs[3]-rhs[2]*m20))*inverse
    t=(m00*(m11*rhs[3]-rhs[2]*m21)-m01*(m10*rhs[3]-rhs[2]*m20)+
       rhs[1]*(m10*m21-m11*m20))*inverse
    w=1-u-v;epsilon=1e-9
    if epsilon<=u<=1-epsilon && epsilon<=v<=1-epsilon && epsilon<=w<=1-epsilon
        return t>0 ? 1 : 0
    elseif u < -epsilon || u > 1+epsilon || v < -epsilon || v > 1+epsilon ||
           w < -epsilon || w > 1+epsilon
        return 0
    end
    return -10000
end

function _model_orientation_same_triangle(first,second)
    for vertex in 1:3, dimension in 1:3
        abs(first[vertex][dimension]-second[vertex][dimension])<1e-12 || return false
    end
    return true
end

function _model_orientation_perturbation(rng)
    return (1e-4*Random.rand(rng),1e-4*Random.rand(rng),1e-4*Random.rand(rng),
            1e-4*Random.rand(rng),1e-4*Random.rand(rng),1e-4*Random.rand(rng))
end

function _model_orientation_ray_count(triangle,surfaces,scale,perturbation)
    a,b,c=triangle
    point=((a[1]+b[1]+c[1])/3,(a[2]+b[2]+c[2])/3,(a[3]+b[3]+c[3])/3)
    edge1=_model_orientation_sub(a,b);edge2=_model_orientation_sub(c,b)
    normal=_model_orientation_unit(_model_orientation_cross(edge1,edge2))
    unit1=_model_orientation_unit(edge1);unit2=_model_orientation_unit(edge2)
    normal=_model_orientation_unit((
        normal[1]+perturbation[1]*unit1[1]+perturbation[2]*unit2[1],
        normal[2]+perturbation[3]*unit1[2]+perturbation[4]*unit2[2],
        normal[3]+perturbation[5]*unit1[3]+perturbation[6]*unit2[3]))
    count=0
    for surface in surfaces, other in surface.stl
        scaled=_model_orientation_scaled(other,scale)
        _model_orientation_same_triangle(triangle,scaled) && continue
        count+=_model_orientation_intersection(scaled,point,normal)
    end
    return count
end

function _model_orientation_cell_normal(surface,cell)
    first=cell.nodes[1];second=cell.nodes[2];third=cell.nodes[3]
    coordinates=surface.coords
    a=(coordinates[1,first],coordinates[2,first],coordinates[3,first])
    b=(coordinates[1,second],coordinates[2,second],coordinates[3,second])
    c=(coordinates[1,third],coordinates[2,third],coordinates[3,third])
    return _model_orientation_unit(_model_orientation_cross(
        _model_orientation_sub(b,a),_model_orientation_sub(c,a)))
end

function _model_orientation_cad_groups(surface)
    # getGFaceOrientation first tries stored vertex parameters across all cells;
    # barycenter reparametrization is a second pass only for unresolved groups.
    nonbl=0;bl=0;hasbl=false
    for cell in surface.cells
        cell.boundary_layer && (hasbl=true)
    end
    hasbl || (bl=1)
    for pass in 1:2
        for cell in surface.cells
            (cell.boundary_layer ? bl : nonbl)==0 || continue
            normal=pass==1 ? cell.vertex_normal : cell.barycenter_normal
            normal===nothing && continue
            sign=_model_orientation_dot(_model_orientation_cell_normal(surface,cell),normal)>0 ? 1 : -1
            if cell.boundary_layer
                bl=sign
            else
                nonbl=sign
            end
            nonbl!=0 && bl!=0 && break
        end
        nonbl!=0 && bl!=0 && break
    end
    return nonbl,bl,hasbl
end

"""
Prepare the actual GRegion Surface flag and orientMeshGFace cell updates.

`surfaces` follows the region's actual bounding-face order. STL inputs retain
their native cached identity; empty STL stays empty, including quad-only native
faces. `scaling` is the norm of the actual region bounding-box diagonal. An
infinite scale reproduces an empty fully discrete Volume bounding box; zero
uses the upstream unit-scale fallback. Each cell key remains bound to its
original storage column. Changed connectivity only permutes slots, including
P2 support rows. Node/element labels, coordinates and metadata are caller-owned.

Normals must be evaluated from the actual stored vertex parameters and native
barycenter fallback before this call. Periodic slaves, fully discrete faces and
BoundaryLayerSurface skip CAD normalization as in orientMeshGFace. Boundary
layer cell groups on ordinary faces are normalized separately. All inputs are
read-only; a diagnostic or thrown error never publishes any live mutation.
"""
function _model_mesh_orientation_plan(region_tag::Integer,
        surfaces::AbstractVector{_ModelMeshOrientationSurface};
        scaling,reverse_flags=Dict{Int,Bool}(),
        rng::Random.AbstractRNG=Random.Xoshiro(0x4d657368),max_retries::Integer=256)
    caller="Model mesh outward orientation"
    !(region_tag isa Bool) && 1<=region_tag<=typemax(Int) || throw(ArgumentError(
        "$caller: region tag must be a positive Int"))
    scaling isa Real && !(scaling isa Bool) && !isnan(scaling) && scaling>=0 ||
        throw(ArgumentError("$caller: bounding-box scale must be nonnegative"))
    scale=Float64(scaling)
    !(max_retries isa Bool) && 1<=max_retries<=typemax(Int) || throw(ArgumentError(
        "$caller: retry count must be a positive Int"))
    flags=Dict{Int,Bool}()
    for (tag,flag) in reverse_flags
        tag isa Integer && !(tag isa Bool) && 1<=tag<=typemax(Int) && flag isa Bool ||
            throw(ArgumentError("$caller: reverse flags must bind positive Surface tags to Bool values"))
        flags[Int(tag)]=flag
    end
    cells=Dict{Int,Vector{Int32}}();diagnostics=String[]
    tags=Set{Int}();keys=Set{Int}()
    cell_count=0
    for surface in surfaces
        cell_count=Base.checked_add(cell_count,length(surface.cells))
    end
    sizehint!(tags,length(surfaces))
    sizehint!(keys,cell_count)
    sizehint!(flags,Base.checked_add(length(flags),length(surfaces)))
    perturbation=_model_orientation_perturbation(rng)
    for surface in surfaces
        surface.tag in tags && throw(ArgumentError(
            "$caller: Surface $(surface.tag) occurs more than once in the region"))
        push!(tags,surface.tag)
        for cell in surface.cells
            cell.key in keys && throw(ArgumentError(
                "$caller: private cell key $(cell.key) is bound more than once"))
            push!(keys,cell.key)
        end
    end
    if scale==0
        scale=1.0
        push!(diagnostics,"Bad scaling in outward orientation for Volume $region_tag; using unit scale")
    end
    for surface in surfaces
        isempty(surface.stl) || continue
        push!(diagnostics,"No valid STL triangulation found for Surface $(surface.tag); skipping outward orientation for Volume $region_tag")
        return (reverse_flags=flags,cells=cells,diagnostics=diagnostics,completed=false)
    end
    for surface in surfaces
        count=-1;attempt=0
        while count<0
            attempt+=1
            attempt<=max_retries || throw(ArgumentError(
                "$caller: no robust ray intersection for Surface $(surface.tag) " *
                "of Volume $region_tag after $max_retries attempts"))
            count=0
            for triangle in surface.stl
                count+=_model_orientation_ray_count(
                    _model_orientation_scaled(triangle,scale),surfaces,scale,perturbation)
                count>=0 && break
            end
            count<0 && (perturbation=_model_orientation_perturbation(rng))
        end
        reverse=isodd(count)
        flags[surface.tag]=reverse
        nonbl=0;bl=0
        if surface.mode==:cad && !isempty(surface.cells)
            nonbl,bl,_=_model_orientation_cad_groups(surface)
            if nonbl==0 && bl==0
                push!(diagnostics,"Could not orient mesh in Surface $(surface.tag)")
                continue
            end
        end
        changed_count=0
        for cell in surface.cells
            normalize=surface.mode==:cad && (cell.boundary_layer ? bl : nonbl)==-1
            xor(normalize,reverse) && (changed_count+=1)
        end
        changed_count==0 && continue
        sizehint!(cells,Base.checked_add(length(cells),changed_count))
        for cell in surface.cells
            normalize=surface.mode==:cad && (cell.boundary_layer ? bl : nonbl)==-1
            xor(normalize,reverse) || continue
            permutation=_model_orientation_permutation(cell.msh,caller)
            replacement=Vector{Int32}(undef,length(permutation))
            for slot in eachindex(permutation)
                replacement[slot]=cell.nodes[permutation[slot]]
            end
            cells[cell.key]=replacement
        end
    end
    return (reverse_flags=flags,cells=cells,diagnostics=diagnostics,completed=true)
end
