"""
    Transform

Validated affine transformations for finalized simplex meshes. Connectivity and
physical tags are preserved. Orientation-reversing maps rewind triangles and
tetrahedra so outward surface orientation and positive tetrahedron orientation are
preserved.
"""
module Transform

using ..MeshTypes: Mesh, validate

export affine_transform, translate_mesh, rotate_mesh, dilate_mesh, mirror_mesh

@inline function _transform_float(value,caller::AbstractString,name::AbstractString)
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    value isa Real || throw(ArgumentError("$caller: $name must be real"))
    result=try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be Float64-representable"))
    end
    isfinite(result) || throw(ArgumentError("$caller: $name must be finite"))
    return result
end

function _transform_point3(value,caller::AbstractString,name::AbstractString)
    (value isa Tuple || value isa AbstractVector) || throw(ArgumentError(
        "$caller: $name must be a three-coordinate tuple or vector"))
    length(value)==3 || throw(ArgumentError(
        "$caller: $name must contain three coordinates"))
    return (_transform_float(value[1],caller,"$name[1]"),
            _transform_float(value[2],caller,"$name[2]"),
            _transform_float(value[3],caller,"$name[3]"))
end

function _transform_matrix(value,caller::AbstractString)
    value isa AbstractMatrix || throw(ArgumentError(
        "$caller: matrix must be a 3×3 real matrix"))
    size(value)==(3,3) || throw(ArgumentError(
        "$caller: matrix must be 3×3"))
    return ntuple(9) do flat
        row=mod(flat-1,3)+1;column=(flat-1)÷3+1
        _transform_float(value[row,column],caller,"matrix[$row,$column]")
    end
end

@inline _matrix_entry(matrix,index)=matrix[index]

function _determinant_sign(matrix::NTuple{9,Float64},caller::AbstractString)
    # Conversion to Rational{BigInt} is exact for Float64. This setup-time
    # certificate prevents a rounded determinant from accepting a singular map
    # or choosing the wrong orientation near cancellation.
    a11=Rational{BigInt}(_matrix_entry(matrix,1))
    a21=Rational{BigInt}(_matrix_entry(matrix,2))
    a31=Rational{BigInt}(_matrix_entry(matrix,3))
    a12=Rational{BigInt}(_matrix_entry(matrix,4))
    a22=Rational{BigInt}(_matrix_entry(matrix,5))
    a32=Rational{BigInt}(_matrix_entry(matrix,6))
    a13=Rational{BigInt}(_matrix_entry(matrix,7))
    a23=Rational{BigInt}(_matrix_entry(matrix,8))
    a33=Rational{BigInt}(_matrix_entry(matrix,9))
    determinant=a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+
                a13*(a21*a32-a22*a31)
    determinant==0 && throw(ArgumentError("$caller: matrix must be nonsingular"))
    return determinant>0 ? 1 : -1
end

function _copy_oriented_cells(mesh::Mesh,orientation::Int)
    segments=copy(mesh.segs)
    triangles=copy(mesh.tris)
    tetrahedra=copy(mesh.tets)
    if orientation<0
        @inbounds for cell in axes(triangles,2)
            triangles[2,cell],triangles[3,cell]=triangles[3,cell],triangles[2,cell]
        end
        @inbounds for cell in axes(tetrahedra,2)
            tetrahedra[1,cell],tetrahedra[2,cell]=tetrahedra[2,cell],tetrahedra[1,cell]
        end
    end
    return segments,triangles,tetrahedra
end

@inline function _affine_coordinate_fast(base::Float64,shift::Float64,
                                         a::Float64,b::Float64,c::Float64,
                                         point1::Float64,point2::Float64,
                                         point3::Float64,origin1::Float64,
                                         origin2::Float64,origin3::Float64)
    x=point1-origin1;y=point2-origin2;z=point3-origin3
    (isfinite(x) && isfinite(y) && isfinite(z)) || return NaN
    linear=muladd(c,z,muladd(b,y,a*x))
    return (base+shift)+linear
end

function _affine_coordinate_exact(base::Float64,shift::Float64,
                                  a::Float64,b::Float64,c::Float64,
                                  point1::Float64,point2::Float64,point3::Float64,
                                  origin1::Float64,origin2::Float64,origin3::Float64,
                                  caller::AbstractString,node::Int)
    # Float64 values are dyadic rationals, so this rare path evaluates the affine
    # coordinate exactly even when an intermediate subtraction or product would
    # overflow or when widely separated exponents cancel. A fixed BigFloat
    # precision cannot provide that guarantee for every finite Float64 input.
    R=Rational{BigInt}
    result=R(base)+R(shift)+R(a)*(R(point1)-R(origin1))+
           R(b)*(R(point2)-R(origin2))+R(c)*(R(point3)-R(origin3))
    converted=Float64(result)
    isfinite(converted) || throw(ArgumentError(
        "$caller: transformed node $node is not Float64-representable"))
    return converted
end

@inline function _affine_coordinate(base::Float64,shift::Float64,
                                    a::Float64,b::Float64,c::Float64,
                                    point1::Float64,point2::Float64,point3::Float64,
                                    origin1::Float64,origin2::Float64,origin3::Float64,
                                    caller::AbstractString,node::Int)
    result=_affine_coordinate_fast(base,shift,a,b,c,point1,point2,point3,
                                   origin1,origin2,origin3)
    isfinite(result) && return result
    return _affine_coordinate_exact(base,shift,a,b,c,point1,point2,point3,
                                    origin1,origin2,origin3,caller,node)
end

"""
    affine_transform(mesh, matrix; origin=(0,0,0), translation=(0,0,0),
                     check=true) -> Mesh

Apply `q = origin + matrix * (p - origin) + translation` to every node. `matrix`
must be a finite nonsingular 3×3 matrix. The input is validated before use and, by
default, the result is independently validated. An orientation-reversing matrix
rewinds triangle and tetrahedron connectivity while retaining cell order and tags.
"""
function affine_transform(mesh::Mesh,matrix;origin=(0.,0.,0.),
                          translation=(0.,0.,0.),check::Bool=true)
    caller="affine_transform"
    diagnostic=validate(mesh)
    diagnostic.ok || throw(ArgumentError(
        "$caller: input mesh is invalid — "*join(diagnostic.messages,"; ")))
    coefficients=_transform_matrix(matrix,caller)
    orientation=_determinant_sign(coefficients,caller)
    base=_transform_point3(origin,caller,"origin")
    shift=_transform_point3(translation,caller,"translation")
    a11,a21,a31,a12,a22,a32,a13,a23,a33=coefficients
    coordinates=Matrix{Float64}(undef,3,size(mesh.coords,2))
    @inbounds for node in axes(mesh.coords,2)
        point1=mesh.coords[1,node];point2=mesh.coords[2,node]
        point3=mesh.coords[3,node]
        qx=_affine_coordinate(base[1],shift[1],a11,a12,a13,
                              point1,point2,point3,base[1],base[2],base[3],
                              caller,node)
        qy=_affine_coordinate(base[2],shift[2],a21,a22,a23,
                              point1,point2,point3,base[1],base[2],base[3],
                              caller,node)
        qz=_affine_coordinate(base[3],shift[3],a31,a32,a33,
                              point1,point2,point3,base[1],base[2],base[3],
                              caller,node)
        coordinates[1,node]=qx;coordinates[2,node]=qy;coordinates[3,node]=qz
    end
    segments,triangles,tetrahedra=_copy_oriented_cells(mesh,orientation)
    result=Mesh(coordinates;segs=segments,tris=triangles,tets=tetrahedra,
                seg_tag=mesh.seg_tag,tri_tag=mesh.tri_tag,tet_tag=mesh.tet_tag)
    if check
        output_diagnostic=validate(result)
        output_diagnostic.ok || throw(ArgumentError(
            "$caller: transformation produced an invalid mesh — "*
            join(output_diagnostic.messages,"; ")))
    end
    return result
end

"""Translate `mesh` by a finite three-coordinate displacement."""
function translate_mesh(mesh::Mesh,displacement;check::Bool=true)
    shift=_transform_point3(displacement,"translate_mesh","displacement")
    return affine_transform(mesh,[1.0 0.0 0.0;0.0 1.0 0.0;0.0 0.0 1.0];
                            translation=shift,check=check)
end

@inline function _unit_axis(axis,caller::AbstractString)
    x,y,z=_transform_point3(axis,caller,"axis")
    magnitude=hypot(x,y,z)
    (isfinite(magnitude) && magnitude>0) || throw(ArgumentError(
        "$caller: axis must have finite positive length"))
    return (x/magnitude,y/magnitude,z/magnitude)
end

"""Rotate `mesh` about the oriented line through `center` along `axis`."""
function rotate_mesh(mesh::Mesh,center,axis,angle::Real;check::Bool=true)
    pivot=_transform_point3(center,"rotate_mesh","center")
    ux,uy,uz=_unit_axis(axis,"rotate_mesh")
    theta=_transform_float(angle,"rotate_mesh","angle")
    sine,cosine=sincos(theta);one_minus=1-cosine
    matrix=(
        cosine+ux*ux*one_minus, uy*ux*one_minus+uz*sine,
        uz*ux*one_minus-uy*sine,
        ux*uy*one_minus-uz*sine, cosine+uy*uy*one_minus,
        uz*uy*one_minus+ux*sine,
        ux*uz*one_minus+uy*sine, uy*uz*one_minus-ux*sine,
        cosine+uz*uz*one_minus,
    )
    return affine_transform(mesh,reshape(collect(matrix),3,3);origin=pivot,
                            check=check)
end

"""
    dilate_mesh(mesh, center, factors; check=true) -> Mesh

Scale about `center` by three finite nonzero factors. Negative factors are allowed;
an odd number of them triggers orientation rewinding.
"""
function dilate_mesh(mesh::Mesh,center,factors;check::Bool=true)
    pivot=_transform_point3(center,"dilate_mesh","center")
    sx,sy,sz=_transform_point3(factors,"dilate_mesh","factors")
    (sx!=0 && sy!=0 && sz!=0) || throw(ArgumentError(
        "dilate_mesh: scale factors must be nonzero"))
    return affine_transform(mesh,[sx 0.0 0.0;0.0 sy 0.0;0.0 0.0 sz];
                            origin=pivot,check=check)
end

"""Reflect `mesh` in the plane through `point` with the given normal."""
function mirror_mesh(mesh::Mesh,point,normal;check::Bool=true)
    pivot=_transform_point3(point,"mirror_mesh","point")
    nx,ny,nz=_unit_axis(normal,"mirror_mesh")
    matrix=[1-2nx*nx -2nx*ny -2nx*nz;
            -2ny*nx 1-2ny*ny -2ny*nz;
            -2nz*nx -2nz*ny 1-2nz*nz]
    return affine_transform(mesh,matrix;origin=pivot,check=check)
end

end # module Transform
