# Point location of the public mesh queries, following Gmsh 4.15.2:
# `MElement::xyz2uvw` with the closed forms of `MTriangle` and `MTetrahedron`,
# the per-family `isInside`, and the `MElementOctree` search (`Octree.cpp`,
# `OctreeInternals.cpp`) with its bucket order, last-hit shortcut and tolerance
# sweep. The standalone locators of `MeshPointLocation` keep their own
# residual-based contract and do not use this file.

# `Mesh.ToleranceReferenceElement` default and the `MElementOctree` bucket size.
const _LOCATION_TOLERANCE=1.0e-6
const _LOCATION_BUCKET_SIZE=100

# --- element table ----------------------------------------------------------

# Every element of the queried cache in the order `GModel::getEntities` and
# `GEntity::getMeshElement` visit them: Points, Curves, Surfaces, Volumes, each
# by entity tag; inside a Surface triangles precede quadrangles and inside a
# Volume tetrahedra precede hexahedra, prisms and pyramids.
struct _LocationElements
    coords::Matrix{Float64}
    cells::Vector{Matrix{Int32}}
    msh::Vector{Int32}
    family::Vector{Symbol}
    dim::Vector{Int8}
    order::Vector{Int8}
    block::Vector{Int32}
    column::Vector{Int32}
    dense::Vector{Int}
    searchable::BitVector
end

function _location_elements(model::GeoModel,mesh::Union{Mesh,MixedMesh},class)
    # An unclassified cache keeps its stored block order inside each dimension.
    catalog,keys=_entity_element_order(mesh,class)
    cells=Matrix{Int32}[entry[4] for entry in catalog]
    # `MElementOctree` leaves out the node of a Point that bounds no Curve.
    bounding_points=Set{Int32}()
    if class!==nothing
        for ((dim,_),children) in class.boundaries
            dim==1 && union!(bounding_points,children)
        end
    end
    for ((dim,_),record) in model.discrete
        dim==1 || continue
        for bound in record.boundary
            bound[1]==0 && push!(bounding_points,Int32(abs(bound[2])))
        end
    end
    for points in values(model.curves), point in points
        push!(bounding_points,Int32(point))
    end
    count=length(keys)
    elements=_LocationElements(mesh.coords,cells,Vector{Int32}(undef,count),
        Vector{Symbol}(undef,count),Vector{Int8}(undef,count),Vector{Int8}(undef,count),
        Vector{Int32}(undef,count),Vector{Int32}(undef,count),Vector{Int}(undef,count),
        trues(count))
    for (element,(dim,owner,_,block,column)) in enumerate(keys)
        msh,_,offset,_,_=catalog[block]
        spec=msh_spec(Int(msh))
        elements.msh[element]=msh;elements.family[element]=spec.family
        elements.dim[element]=dim;elements.order[element]=spec.order
        elements.block[element]=block;elements.column[element]=column
        elements.dense[element]=offset+column
        dim==0 && !(owner in bounding_points) && (elements.searchable[element]=false)
    end
    return elements
end

@inline function _location_node(elements::_LocationElements,element::Int,local_node::Int)
    @inbounds node=elements.cells[elements.block[element]][local_node,elements.column[element]]
    @inbounds return (elements.coords[1,node],elements.coords[2,node],elements.coords[3,node])
end

@inline _location_node_count(elements::_LocationElements,element::Int)=
    size(elements.cells[elements.block[element]],1)

# --- reference coordinates ---------------------------------------------------

# Nodal values and gradients of the element's own Lagrange space at `uvw`.
function _location_shape(msh::Int,family::Symbol,order::Int,uvw::NTuple{3,Float64},
                         caller::AbstractString)
    u,v,w=uvw
    if family===:pnt
        return (1.0,),((0.0,0.0,0.0),)
    elseif msh==1
        return (0.5*(1.0-u),0.5*(1.0+u)),((-0.5,0.0,0.0),(0.5,0.0,0.0))
    elseif msh==3
        return (0.25*(1.0-u)*(1.0-v),0.25*(1.0+u)*(1.0-v),
                0.25*(1.0+u)*(1.0+v),0.25*(1.0-u)*(1.0+v)),
               ((-0.25*(1.0-v),-0.25*(1.0-u),0.0),(0.25*(1.0-v),-0.25*(1.0+u),0.0),
                (0.25*(1.0+v),0.25*(1.0+u),0.0),(-0.25*(1.0+v),0.25*(1.0-u),0.0))
    elseif family in (:tri,:qua,:tet,:hex,:pri) && 1<=order<=10
        return _native_cell_shape(_native_cell_polynomial_cached(msh,caller),uvw)
    end
    coordinates=Float64[u,v,w]
    _,values,_=mesh_basis_functions(msh,coordinates,"Lagrange";caller)
    _,gradients,_=mesh_basis_functions(msh,coordinates,"GradLagrange";caller)
    return values,[(gradients[3i-2],gradients[3i-1],gradients[3i]) for i in eachindex(values)]
end

@inline function _location_unit(a::NTuple{3,Float64})
    modulus=sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3])
    modulus==0.0 && return a
    reciprocal=1.0/modulus
    return (a[1]*reciprocal,a[2]*reciprocal,a[3]*reciprocal)
end

@inline _location_cross(a,b)=(a[2]*b[3]-a[3]*b[2],-a[1]*b[3]+a[3]*b[1],a[1]*b[2]-a[2]*b[1])

# `MElement::getJacobian`: the rows of the Jacobian summed over the nodes in
# order, completed to a 3x3 frame by `_computeDeterminantAndRegularize`.
function _location_jacobian(elements::_LocationElements,element::Int,gradients)
    dim=Int(elements.dim[element])
    a=(0.0,0.0,0.0);b=a;c=a
    @inbounds for node in 1:_location_node_count(elements,element)
        point=_location_node(elements,element,node);gradient=gradients[node]
        dim>=1 && (a=(a[1]+point[1]*gradient[1],a[2]+point[2]*gradient[1],a[3]+point[3]*gradient[1]))
        dim>=2 && (b=(b[1]+point[1]*gradient[2],b[2]+point[2]*gradient[2],b[3]+point[3]*gradient[2]))
        dim>=3 && (c=(c[1]+point[1]*gradient[3],c[2]+point[2]*gradient[3],c[3]+point[3]*gradient[3]))
    end
    if dim==0
        return (1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0),1.0
    elseif dim==1
        determinant=sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3])
        b=(abs(a[1])>=abs(a[2]) && abs(a[1])>=abs(a[3])) ||
          (abs(a[2])>=abs(a[1]) && abs(a[2])>=abs(a[3])) ? (a[2],-a[1],0.0) : (0.0,a[3],-a[2])
        b=_location_unit(b)
        c=_location_unit(_location_cross(a,b))
        return a,b,c,determinant
    elseif dim==2
        first=a[1]*b[2]-a[2]*b[1];second=a[3]*b[1]-a[1]*b[3];third=a[2]*b[3]-a[3]*b[2]
        determinant=sqrt(first*first+second*second+third*third)
        c=_location_unit(_location_cross(a,b))
        return a,b,c,determinant
    end
    determinant=a[1]*b[2]*c[3]+a[3]*b[1]*c[2]+a[2]*b[3]*c[1]-
        a[3]*b[2]*c[1]-a[1]*b[3]*c[2]-a[2]*b[1]*c[3]
    return a,b,c,determinant
end

# `MTriangle::xyz2uvw`: the point is dropped along the coordinate axis the
# triangle faces most, and `w` is zero wherever the point is.
function _location_triangle_uvw(elements::_LocationElements,element::Int,p::NTuple{3,Float64})
    o=_location_node(elements,element,1)
    first=_location_node(elements,element,2);second=_location_node(elements,element,3)
    d=(p[1]-o[1],p[2]-o[2],p[3]-o[3])
    d1=(first[1]-o[1],first[2]-o[2],first[3]-o[3])
    d2=(second[1]-o[1],second[2]-o[2],second[3]-o[3])
    jxy=d1[1]*d2[2]-d1[2]*d2[1]
    jxz=d1[1]*d2[3]-d1[3]*d2[1]
    jyz=d1[2]*d2[3]-d1[3]*d2[2]
    if abs(jxy)>abs(jxz) && abs(jxy)>abs(jyz)
        return ((d[1]*d2[2]-d[2]*d2[1])/jxy,(d[2]*d1[1]-d[1]*d1[2])/jxy,0.0)
    elseif abs(jxz)>abs(jyz)
        return ((d[1]*d2[3]-d[3]*d2[1])/jxz,(d[3]*d1[1]-d[1]*d1[3])/jxz,0.0)
    end
    return ((d[2]*d2[3]-d[3]*d2[2])/jyz,(d[3]*d1[2]-d[2]*d1[3])/jyz,0.0)
end

# `MTetrahedron::xyz2uvw`: `sys3x3` on the three edge vectors.
function _location_tetrahedron_uvw(elements::_LocationElements,element::Int,p::NTuple{3,Float64})
    o=_location_node(elements,element,1)
    e1=_location_node(elements,element,2);e2=_location_node(elements,element,3)
    e3=_location_node(elements,element,4)
    m11=e1[1]-o[1];m12=e2[1]-o[1];m13=e3[1]-o[1]
    m21=e1[2]-o[2];m22=e2[2]-o[2];m23=e3[2]-o[2]
    m31=e1[3]-o[3];m32=e2[3]-o[3];m33=e3[3]-o[3]
    b1=p[1]-o[1];b2=p[2]-o[2];b3=p[3]-o[3]
    determinant=m11*(m22*m33-m23*m32)-m12*(m21*m33-m23*m31)+m13*(m21*m32-m22*m31)
    determinant==0.0 && return (0.0,0.0,0.0)
    reciprocal=1.0/determinant
    u=b1*(m22*m33-m23*m32)-m12*(b2*m33-m23*b3)+m13*(b2*m32-m22*b3)
    v=m11*(b2*m33-m23*b3)-b1*(m21*m33-m23*m31)+m13*(m21*b3-b2*m31)
    w=m11*(m22*b3-b2*m32)-m12*(m21*b3-b2*m31)+b1*(m21*m32-m22*m31)
    return (u*reciprocal,v*reciprocal,w*reciprocal)
end

# `MElement::xyz2uvw`. The Newton loop stops on the size of its step, so a
# point off a Curve or Surface element keeps adding its out-of-element offset
# to the unused coordinates until the nineteenth step; that is the value Gmsh
# reports and tests in `isInside`.
function _location_uvw(elements::_LocationElements,element::Int,p::NTuple{3,Float64},
                       caller::AbstractString)
    msh=Int(elements.msh[element])
    msh==2 && return _location_triangle_uvw(elements,element,p)
    msh==4 && return _location_tetrahedron_uvw(elements,element,p)
    family=elements.family[element];order=Int(elements.order[element])
    uvw=(0.0,0.0,0.0)
    if order>2
        nearest=1;distance=Inf
        for node in 1:_location_node_count(elements,element)
            point=_location_node(elements,element,node)
            candidate=(point[1]-p[1])*(point[1]-p[1])+(point[2]-p[2])*(point[2]-p[2])+
                (point[3]-p[3])*(point[3]-p[3])
            if node==1 || candidate<distance
                nearest=node;distance=candidate
            end
        end
        reference=Elements.lagrange_nodes(msh)
        dim=Int(elements.dim[element])
        uvw=ntuple(axis->axis<=dim ? reference[axis,nearest] : 0.0,3)
    end
    iteration=1;step=1.0
    while step>1.0e-6 && iteration<20
        values,gradients=_location_shape(msh,family,order,uvw,caller)
        a,b,c,determinant=_location_jacobian(elements,element,gradients)
        determinant==0.0 && break
        mapped=(0.0,0.0,0.0)
        @inbounds for node in eachindex(values)
            point=_location_node(elements,element,node);value=values[node]
            mapped=(mapped[1]+point[1]*value,mapped[2]+point[2]*value,mapped[3]+point[3]*value)
        end
        # `inv3x3` of the frame with rows a, b, c.
        frame=a[1]*(b[2]*c[3]-b[3]*c[2])-a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1])
        frame==0.0 && throw(ArgumentError("$caller: Singular matrix 3x3"))
        r=1.0/frame
        i11=(b[2]*c[3]-b[3]*c[2])*r;i21=-(b[1]*c[3]-b[3]*c[1])*r;i31=(b[1]*c[2]-b[2]*c[1])*r
        i12=-(a[2]*c[3]-a[3]*c[2])*r;i22=(a[1]*c[3]-a[3]*c[1])*r;i32=-(a[1]*c[2]-a[2]*c[1])*r
        i13=(a[2]*b[3]-a[3]*b[2])*r;i23=-(a[1]*b[3]-a[3]*b[1])*r;i33=(a[1]*b[2]-a[2]*b[1])*r
        dx=p[1]-mapped[1];dy=p[2]-mapped[2];dz=p[3]-mapped[3]
        un=uvw[1]+i11*dx+i21*dy+i31*dz
        vn=uvw[2]+i12*dx+i22*dy+i32*dz
        wn=uvw[3]+i13*dx+i23*dy+i33*dz
        step=sqrt((un-uvw[1])*(un-uvw[1])+(vn-uvw[2])*(vn-uvw[2])+(wn-uvw[3])*(wn-uvw[3]))
        uvw=(un,vn,wn);iteration+=1
    end
    return uvw
end

# The `isInside` of each element family. The tests are written as upstream
# writes them: a coordinate that is not a number fails none of them.
function _location_inside(family::Symbol,uvw::NTuple{3,Float64},tolerance::Float64)
    u,v,w=uvw
    if family===:pnt
        return !(abs(u)>tolerance || abs(v)>tolerance || abs(w)>tolerance)
    elseif family===:lin
        return !(u< -(1.0+tolerance) || u>(1.0+tolerance) || abs(v)>tolerance || abs(w)>tolerance)
    elseif family===:tri
        return !(u< -tolerance || v< -tolerance || u>((1.0+tolerance)-v) || abs(w)>tolerance)
    elseif family===:qua
        return !(u< -(1.0+tolerance) || v< -(1.0+tolerance) || u>(1.0+tolerance) ||
                 v>(1.0+tolerance) || abs(w)>tolerance)
    elseif family===:tet
        return !(u< -tolerance || v< -tolerance || w< -tolerance || u>((1.0+tolerance)-v-w))
    elseif family===:hex
        return !(u< -(1.0+tolerance) || v< -(1.0+tolerance) || w< -(1.0+tolerance) ||
                 u>(1.0+tolerance) || v>(1.0+tolerance) || w>(1.0+tolerance))
    elseif family===:pri
        return !(w>(1.0+tolerance) || w< -(1.0+tolerance) || u< -tolerance || v< -tolerance ||
                 u>((1.0+tolerance)-v))
    end
    return !(u<(w-(1.0+tolerance)) || u>((1.0+tolerance)-w) || v<(w-(1.0+tolerance)) ||
             v>((1.0+tolerance)-w) || w< -tolerance || w>(1.0+tolerance))
end

@inline _location_contains(elements::_LocationElements,element::Int,p::NTuple{3,Float64},
                           tolerance::Float64,caller::AbstractString)=
    _location_inside(elements.family[element],_location_uvw(elements,element,p,caller),tolerance)

# --- bounding boxes ------------------------------------------------------------

@inline function _location_thicken(lower::NTuple{3,Float64},upper::NTuple{3,Float64})
    dx=upper[1]-lower[1];dy=upper[2]-lower[2];dz=upper[3]-lower[3]
    d=0.01*sqrt(dx*dx+dy*dy+dz*dz)
    return (lower[1]-d,lower[2]-d,lower[3]-d),(upper[1]+d,upper[2]+d,upper[3]+d)
end

# Lagrange values of an element type at the lattice of its complete Bernstein
# space: the sampling of `MElement::getBezierVerticesCoord`.
const _LOCATION_BEZIER_SAMPLES=Dict{Int,Any}()
function _location_bezier_samples(msh::Int,caller::AbstractString)
    return lock(_NATIVE_CELL_POLYNOMIAL_LOCK) do
        get!(_LOCATION_BEZIER_SAMPLES,msh) do
            spec=msh_spec(msh)
            if spec.family===:lin
                points=NTuple{3,Float64}[(2.0*index/spec.order-1.0,0.0,0.0) for index in 0:spec.order]
                space=nothing
            elseif spec.family in (:tri,:qua,:tet,:hex,:pri)
                space=_native_bezier_space(spec.family,spec.order)
                points=space.points
            else
                return nothing
            end
            matrix=Matrix{Float64}(undef,length(points),spec.nnodes)
            for (sample,point) in enumerate(points)
                values,_=_location_shape(msh,spec.family,spec.order,point,caller)
                for node in 1:spec.nnodes
                    matrix[sample,node]=values[node]
                end
            end
            (space,matrix)
        end
    end
end

# `MElementBB`: the box of the vertices of a first-order element, or of the
# Bezier control points of a curved one, made one percent of its diagonal
# thicker. A family without a Bernstein space here (the pyramid) falls back to
# its node box.
function _location_element_box(elements::_LocationElements,element::Int,caller::AbstractString)
    count=_location_node_count(elements,element)
    lower=upper=_location_node(elements,element,1)
    samples=elements.order[element]>1 ?
        _location_bezier_samples(Int(elements.msh[element]),caller) : nothing
    if samples===nothing
        for node in 2:count
            point=_location_node(elements,element,node)
            lower=(min(lower[1],point[1]),min(lower[2],point[2]),min(lower[3],point[3]))
            upper=(max(upper[1],point[1]),max(upper[2],point[2]),max(upper[3],point[3]))
        end
        return _location_thicken(lower,upper)
    end
    space,matrix=samples
    extrema=ntuple(3) do axis
        sampled=zeros(Float64,size(matrix,1))
        for node in 1:count
            coordinate=_location_node(elements,element,node)[axis]
            for sample in axes(matrix,1)
                sampled[sample]+=matrix[sample,node]*coordinate
            end
        end
        control=if space===nothing
            _native_bezier_convert_line!(sampled,length(sampled)-1,1,1,similar(sampled))
        else
            _native_bezier_coefficients(space,sampled)
        end
        low=high=control[1]
        for value in control
            low=min(low,value);high=max(high,value)
        end
        (low,high)
    end
    return _location_thicken((extrema[1][1],extrema[2][1],extrema[3][1]),
        (extrema[1][2],extrema[2][2],extrema[3][2]))
end

function _location_centroid(elements::_LocationElements,element::Int)
    count=_location_node_count(elements,element)
    total=_location_node(elements,element,1)
    for node in 2:count
        point=_location_node(elements,element,node)
        total=(total[1]+point[1],total[2]+point[2],total[3]+point[3])
    end
    reciprocal=1.0/count
    return (total[1]*reciprocal,total[2]*reciprocal,total[3]*reciprocal)
end

# --- octree ----------------------------------------------------------------------

# Bucket 1 is the head; buckets are allocated eight at a time. `head` and
# `next` thread the centroid list of a leaf (last inserted first), `boxes` is
# its `listBB`, and `previous` is `ptrToPrevElement`.
mutable struct _LocationOctree
    source::Union{Mesh,MixedMesh}
    elements::_LocationElements
    lower::Vector{NTuple{3,Float64}}
    upper::Vector{NTuple{3,Float64}}
    centroid::Vector{NTuple{3,Float64}}
    bucket_lower::Vector{NTuple{3,Float64}}
    bucket_upper::Vector{NTuple{3,Float64}}
    bucket_child::Vector{Int}
    bucket_head::Vector{Int}
    bucket_count::Vector{Int}
    bucket_boxes::Vector{Vector{Int}}
    next::Vector{Int}
    inserted::Vector{Int}
    previous::Int
end

function _location_add_children!(tree::_LocationOctree,bucket::Int,half::NTuple{3,Float64})
    first=length(tree.bucket_child)+1
    base=tree.bucket_lower[bucket]
    for k in 0:1,j in 0:1,i in 0:1
        push!(tree.bucket_lower,(base[1]+half[1]*i,base[2]+half[2]*j,base[3]+half[3]*k))
        push!(tree.bucket_upper,(base[1]+half[1]*(i+1),base[2]+half[2]*(j+1),base[3]+half[3]*(k+1)))
        push!(tree.bucket_child,0);push!(tree.bucket_head,0);push!(tree.bucket_count,0)
        push!(tree.bucket_boxes,Int[])
    end
    tree.bucket_child[bucket]=first
    return first
end

# `findElementBucket`: the first child, in index order, whose closed box holds
# the point, followed down to a leaf; 0 when a level has none.
function _location_bucket(tree::_LocationOctree,start::Int,p::NTuple{3,Float64})
    found=0
    children=tree.bucket_child[start]
    @inbounds while children!=0
        hit=0
        for bucket in children:children+7
            lower=tree.bucket_lower[bucket];upper=tree.bucket_upper[bucket]
            if !(lower[1]>p[1] || upper[1]<p[1] || lower[2]>p[2] || upper[2]<p[2] ||
                 lower[3]>p[3] || upper[3]<p[3])
                hit=bucket;break
            end
        end
        hit==0 && return 0
        found=hit;children=tree.bucket_child[hit]
    end
    return found
end

# `addElement2Bucket`.
function _location_insert!(tree::_LocationOctree,element::Int)
    bucket=_location_bucket(tree,1,tree.centroid[element])
    bucket==0 && return nothing
    link=tree.bucket_head[bucket]
    while link!=0
        link==element && return nothing
        link=tree.next[link]
    end
    push!(tree.inserted,element)
    tree.next[element]=tree.bucket_head[bucket]
    tree.bucket_head[bucket]=element
    tree.bucket_count[bucket]+=1
    again=true
    while again
        again=false
        tree.bucket_count[bucket]>_LOCATION_BUCKET_SIZE || break
        lower=tree.bucket_lower[bucket];upper=tree.bucket_upper[bucket]
        _location_add_children!(tree,bucket,((upper[1]-lower[1])/2,(upper[2]-lower[2])/2,
            (upper[3]-lower[3])/2))
        link=tree.bucket_head[bucket]
        while link!=0
            target=_location_bucket(tree,bucket,tree.centroid[link])
            moved=link;link=tree.next[link]
            target==0 && error("API mesh point location: null bucket in octree")
            tree.next[moved]=tree.bucket_head[target]
            tree.bucket_head[target]=moved
            tree.bucket_count[target]+=1
            if tree.bucket_count[target]>_LOCATION_BUCKET_SIZE
                again=true
                tree.bucket_head[bucket]=0
                bucket=target
            end
        end
        again || (tree.bucket_head[bucket]=0)
    end
    return nothing
end

# `insertOneBB`.
function _location_insert_box!(tree::_LocationOctree,element::Int,bucket::Int)
    lower=tree.bucket_lower[bucket];upper=tree.bucket_upper[bucket]
    low=tree.lower[element];high=tree.upper[element]
    (lower[1]>high[1] || upper[1]<low[1] || lower[2]>high[2] || upper[2]<low[2] ||
     lower[3]>high[3] || upper[3]<low[3]) && return nothing
    children=tree.bucket_child[bucket]
    if children==0
        link=tree.bucket_head[bucket]
        while link!=0
            link==element && return nothing
            link=tree.next[link]
        end
        push!(tree.bucket_boxes[bucket],element)
        return nothing
    end
    for child in children:children+7
        _location_insert_box!(tree,element,child)
    end
    return nothing
end

# `MElementOctree::MElementOctree(GModel *)`. The root box is `GModel::bounds`
# — the mesh nodes and the Point coordinates — made one percent thicker.
function _LocationOctree(model::GeoModel,mesh::Union{Mesh,MixedMesh},class,caller::AbstractString)
    elements=_location_elements(model,mesh,class)
    count=length(elements.msh)
    lower=(Inf,Inf,Inf);upper=(-Inf,-Inf,-Inf)
    for node in axes(mesh.coords,2)
        point=(mesh.coords[1,node],mesh.coords[2,node],mesh.coords[3,node])
        lower=(min(lower[1],point[1]),min(lower[2],point[2]),min(lower[3],point[3]))
        upper=(max(upper[1],point[1]),max(upper[2],point[2]),max(upper[3],point[3]))
    end
    if mesh isa MixedMesh && mesh.entity_data!==nothing
        for ((dim,_),entity) in mesh.entity_data.entities
            dim==0 || continue
            point=(entity.bbox[1],entity.bbox[2],entity.bbox[3])
            all(isfinite,point) || continue
            lower=(min(lower[1],point[1]),min(lower[2],point[2]),min(lower[3],point[3]))
            upper=(max(upper[1],point[1]),max(upper[2],point[2]),max(upper[3],point[3]))
        end
    end
    isfinite(lower[1]) || (lower=upper=(0.0,0.0,0.0))
    lower,upper=_location_thicken(lower,upper)
    size=(upper[1]-lower[1],upper[2]-lower[2],upper[3]-lower[3])
    margin=(size[1]*0.01,size[2]*0.01,size[3]*0.01)
    tree=_LocationOctree(mesh,elements,
        Vector{NTuple{3,Float64}}(undef,count),Vector{NTuple{3,Float64}}(undef,count),
        Vector{NTuple{3,Float64}}(undef,count),
        NTuple{3,Float64}[(lower[1]-margin[1],lower[2]-margin[2],lower[3]-margin[3])],
        NTuple{3,Float64}[(size[1]+lower[1]+margin[1],size[2]+lower[2]+margin[2],
            size[3]+lower[3]+margin[3])],
        Int[0],Int[0],Int[0],Vector{Int}[Int[]],zeros(Int,count),Int[],0)
    _location_add_children!(tree,1,((size[1]+2*margin[1])/2,(size[2]+2*margin[2])/2,
        (size[3]+2*margin[3])/2))
    for element in 1:count
        tree.lower[element],tree.upper[element]=_location_element_box(elements,element,caller)
        tree.centroid[element]=_location_centroid(elements,element)
    end
    for element in 1:count
        elements.searchable[element] && _location_insert!(tree,element)
    end
    for element in copy(tree.inserted)
        _location_insert_box!(tree,element,1)
    end
    return tree
end

@inline function _location_in_box(tree::_LocationOctree,element::Int,p::NTuple{3,Float64})
    lower=tree.lower[element];upper=tree.upper[element]
    return !(p[1]>upper[1] || p[1]<lower[1] || p[2]>upper[2] || p[2]<lower[2] ||
             p[3]>upper[3] || p[3]<lower[3])
end

@inline _location_hit(tree::_LocationOctree,element::Int,p::NTuple{3,Float64},
                      caller::AbstractString)=
    _location_in_box(tree,element,p) &&
    _location_contains(tree.elements,element,p,_LOCATION_TOLERANCE,caller)

# `searchElement`: the last element found is tried first, then the leaf's
# centroid list and its box list; the first hit is returned and remembered.
function _location_search(tree::_LocationOctree,p::NTuple{3,Float64},caller::AbstractString)
    previous=tree.previous
    previous!=0 && _location_hit(tree,previous,p,caller) && return previous
    bucket=_location_bucket(tree,1,p)
    bucket==0 && return 0
    link=tree.bucket_head[bucket]
    while link!=0
        if _location_hit(tree,link,p,caller)
            tree.previous=link
            return link
        end
        link=tree.next[link]
    end
    for element in tree.bucket_boxes[bucket]
        if _location_hit(tree,element,p,caller)
            tree.previous=element
            return element
        end
    end
    return 0
end

# `searchAllElements`.
function _location_search_all(tree::_LocationOctree,p::NTuple{3,Float64},caller::AbstractString)
    found=Int[]
    bucket=_location_bucket(tree,1,p)
    bucket==0 && return found
    link=tree.bucket_head[bucket]
    while link!=0
        _location_hit(tree,link,p,caller) && push!(found,link)
        link=tree.next[link]
    end
    for element in tree.bucket_boxes[bucket]
        _location_hit(tree,element,p,caller) && push!(found,element)
    end
    return found
end

# The tolerances of the relaxed sweep: the reference tolerance times ten until
# it reaches one, in the floating-point values that loop produces.
function _location_relaxed_tolerances()
    tolerances=Float64[]
    tolerance=_LOCATION_TOLERANCE
    while tolerance<1.0
        tolerance*=10
        push!(tolerances,tolerance)
    end
    return tolerances
end

# `MElementOctree::find`. A relaxed search that finds nothing in the octree
# sweeps every element of the model, searchable or not, in entity order.
function _location_find(tree::_LocationOctree,p::NTuple{3,Float64},dim::Int,strict::Bool,
                        caller::AbstractString)
    elements=tree.elements
    element=_location_search(tree,p,caller)
    element!=0 && (dim==-1 || elements.dim[element]==dim) && return element
    if element!=0
        for candidate in _location_search_all(tree,p,caller)
            elements.dim[candidate]==dim && return candidate
        end
    end
    strict && return 0
    for tolerance in _location_relaxed_tolerances(), candidate in eachindex(elements.msh)
        (dim==-1 || elements.dim[candidate]==dim) || continue
        _location_contains(elements,candidate,p,tolerance,caller) && return candidate
    end
    return 0
end

# `MElementOctree::findAll`.
function _location_find_all(tree::_LocationOctree,p::NTuple{3,Float64},dim::Int,strict::Bool,
                            caller::AbstractString)
    elements=tree.elements
    found=Int[candidate for candidate in _location_search_all(tree,p,caller)
              if dim==-1 || elements.dim[candidate]==dim]
    (isempty(found) && !strict) || return found
    for tolerance in _location_relaxed_tolerances()
        for candidate in eachindex(elements.msh)
            (dim==-1 || elements.dim[candidate]==dim) || continue
            _location_contains(elements,candidate,p,tolerance,caller) && push!(found,candidate)
        end
        isempty(found) || return found
    end
    return found
end

# --- session cache ------------------------------------------------------------------

# Gmsh keeps one octree per model until `destroyMeshCaches`; the tree and its
# last-hit shortcut live as long as the published cache they were built from.
const LAST_LOCATION_OCTREE=Ref{Union{Nothing,_LocationOctree}}(nothing)

function _location_octree_locked(mesh::Union{Mesh,MixedMesh},class,caller::AbstractString)
    tree=LAST_LOCATION_OCTREE[]
    if tree===nothing || tree.source!==mesh
        tree=_LocationOctree(_model_locked(),mesh,class,caller)
        mesh===LAST_MESH[] && (LAST_LOCATION_OCTREE[]=tree)
    end
    return tree
end

# `MElement::xyz2uvw` of one dense cache element, for
# `get_local_coordinates_in_element`.
function _location_element_uvw(mesh::Union{Mesh,MixedMesh},dense::Int,p::NTuple{3,Float64},
                               caller::AbstractString)
    offset=0
    for (block,(msh,dim,_,cells,_)) in enumerate(_cache_catalog(mesh))
        if offset<dense<=offset+size(cells,2)
            spec=msh_spec(Int(msh))
            elements=_LocationElements(mesh.coords,Matrix{Int32}[cells],Int32[msh],
                Symbol[spec.family],Int8[dim],Int8[spec.order],Int32[1],
                Int32[dense-offset],Int[dense],trues(1))
            return _location_uvw(elements,1,p,caller)
        end
        offset+=size(cells,2)
    end
    throw(ArgumentError("$caller: unknown element tag $dense"))
end
