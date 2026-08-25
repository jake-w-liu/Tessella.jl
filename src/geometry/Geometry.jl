"""
    Geometry

Stage-5 native constructive primitives (PLAN.md §3 "Geometry", §5 native CSG
path). Each builder returns a **closed, manifold, outward-oriented** triangle
[`Mesh`](@ref) — a boundary surface ready for [`Mesh3D.tetrahedralize`](@ref) /
[`mesh_volume`](@ref). All are verified `Heal.is_meshable` and, when filled, to
reproduce the exact analytic volume.

These are the primitives the ASCENT `solid_model` emits (boxes, cylinders, cones,
geodesic spheres, holes/bores, and axis-aligned cavities via
[`box_shell_surface`](@ref)); general Boolean composition is provided by
`Mesh3D.mesh_boolean`.
"""
module Geometry

using ..MeshTypes: Mesh, validate
using ..Heal: is_meshable

export box_surface, cylinder_surface, sphere_surface, cone_surface
export box_tunnel_surface, box_shell_surface

@inline function _finite_real(x,caller::AbstractString,name::AbstractString)
    x isa Real || throw(ArgumentError("$caller: $name must be real"))
    x isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    y = try
        Float64(x)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be Float64-representable: $(sprint(showerror, err))"))
    end
    isfinite(y) || throw(ArgumentError("$caller: $name must be finite (got $x)"))
    return y
end

@inline function _point3(p, caller::AbstractString, name::AbstractString)
    count=try
        length(p)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name needs three coordinates"))
    end
    count>=3 || throw(ArgumentError("$caller: $name needs three coordinates"))
    values=try
        (p[1],p[2],p[3])
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "$caller: $name coordinates must support integer indexing"))
    end
    return (_finite_real(values[1],caller,"$name[1]"),
            _finite_real(values[2],caller,"$name[2]"),
            _finite_real(values[3],caller,"$name[3]"))
end

# ── box ─────────────────────────────────────────────────────────────────────────
"""
    box_surface(x0,x1, y0,y1, z0,z1) -> Mesh

Closed surface of the axis-aligned box `[x0,x1]×[y0,y1]×[z0,z1]` (12 outward
triangles).
"""
function box_surface(x0,x1,y0,y1,z0,z1)
    x0=_finite_real(x0,"box_surface","x0"); x1=_finite_real(x1,"box_surface","x1")
    y0=_finite_real(y0,"box_surface","y0"); y1=_finite_real(y1,"box_surface","y1")
    z0=_finite_real(z0,"box_surface","z0"); z1=_finite_real(z1,"box_surface","z1")
    (x0<x1 && y0<y1 && z0<z1) || throw(ArgumentError("box_surface: need x0<x1, y0<y1, z0<z1"))
    C = Float64[x0 x1 x1 x0 x0 x1 x1 x0; y0 y0 y1 y1 y0 y0 y1 y1; z0 z0 z0 z0 z1 z1 z1 z1]
    F = [(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),
         (2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]
    t = Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); t[:,k]=Int32[f...]; end
    return _checked_surface(Mesh(C; tris=t),"box_surface")
end

# ── cylinder (solid, watertight lateral + caps) ─────────────────────────────────
"""
    cylinder_surface(center, axis, radius, height; nθ=24, nz=2,
                     max_nodes=10_000_000,
                     max_triangles=20_000_000) -> Mesh

Closed surface of a solid cylinder: base `center`, unit-ish `axis`, `radius`,
`height`. `nθ` circumferential sectors, `nz` axial levels. Lateral wall + two cap
fans share the rim rings (watertight by construction). Resource limits are checked
before coordinate or topology allocation.
"""
function cylinder_surface(center,axis,radius,height;nθ=24,nz=2,
                          max_nodes=10_000_000,
                          max_triangles=20_000_000)
    R=_finite_real(radius,"cylinder_surface","radius")
    H=_finite_real(height,"cylinder_surface","height")
    (R>0 && H>0) || throw(ArgumentError(
        "cylinder_surface: radius and height must be positive"))
    ntheta=_geometry_count(nθ,"cylinder_surface","nθ",3)
    nlevels=_geometry_count(nz,"cylinder_surface","nz",2)
    nodes = try
        Base.checked_add(Base.checked_mul(ntheta,nlevels),2)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("cylinder_surface: requested node count overflows the platform Int limit"))
    end
    triangles=try
        Base.checked_mul(Base.checked_mul(2,ntheta),nlevels)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("cylinder_surface: requested triangle count overflows the platform Int limit"))
    end
    _checked_primitive_counts(nodes,triangles,max_nodes,max_triangles,
                              "cylinder_surface")
    c=_point3(center,"cylinder_surface","center")
    ex,ey,ez=_axis_frame(axis,"cylinder_surface")
    function on(θ,z)
        point=_geometry_point(
            c[1]+R*cos(θ)*ex[1]+R*sin(θ)*ey[1]+z*ez[1],
            c[2]+R*cos(θ)*ex[2]+R*sin(θ)*ey[2]+z*ez[2],
            c[3]+R*cos(θ)*ex[3]+R*sin(θ)*ey[3]+z*ez[3],
            "cylinder_surface")
        return _certify_cylindrical_point(
            point,c,ez,R,z,"cylinder_surface")
    end
    V=Tuple{Float64,Float64,Float64}[]
    sizehint!(V,nodes)
    for j in 0:nlevels-1, i in 0:ntheta-1
        push!(V,on(2π*i/ntheta,H*(j/(nlevels-1))))
    end
    ci=length(V)+1
    push!(V,_certify_cylindrical_point(c,c,ez,0.0,0.0,
                                       "cylinder_surface"))
    cti=length(V)+1
    top=_geometry_point(c[1]+H*ez[1],c[2]+H*ez[2],c[3]+H*ez[3],
                        "cylinder_surface")
    push!(V,_certify_cylindrical_point(top,c,ez,0.0,H,
                                       "cylinder_surface"))
    idx(j,i)=(j-1)*ntheta + mod(i,ntheta) + 1
    Tr=NTuple{3,Int32}[]
    for j in 1:nlevels-1, i in 0:ntheta-1                            # wall (outward)
        a=idx(j,i);b=idx(j,i+1);cc=idx(j+1,i+1);d=idx(j+1,i)
        push!(Tr,(Int32(a),Int32(b),Int32(cc))); push!(Tr,(Int32(a),Int32(cc),Int32(d)))
    end
    for i in 0:ntheta-1; push!(Tr,(Int32(ci),Int32(idx(1,i+1)),Int32(idx(1,i)))); end   # bottom cap
    for i in 0:ntheta-1; push!(Tr,(Int32(cti),Int32(idx(nlevels,i)),Int32(idx(nlevels,i+1)))); end # top cap
    C=Matrix{Float64}(undef,3,length(V)); for (k,p) in enumerate(V); C[:,k]=[p...]; end
    tm=Matrix{Int32}(undef,3,length(Tr)); for (k,f) in enumerate(Tr); tm[:,k]=Int32[f...]; end
    return _checked_surface(Mesh(C; tris=tm),"cylinder_surface")
end

# ── box with a rectangular through-tunnel (genus-1 bore) ────────────────────────
"""
    box_tunnel_surface(ox0,ox1, oy0,oy1, z0,z1, ix0,ix1, iy0,iy1) -> Mesh

Closed surface of the box `[ox0,ox1]×[oy0,oy1]×[z0,z1]` with a rectangular
through-tunnel `[ix0,ix1]×[iy0,iy1]` along `z` (a genus-1 solid — the coax-bore
analogue). The inner `[ix,iy]` rectangle must lie strictly inside the outer one.
"""
function box_tunnel_surface(ox0,ox1, oy0,oy1, z0,z1, ix0,ix1, iy0,iy1)
    raw = (ox0,ox1,oy0,oy1,z0,z1,ix0,ix1,iy0,iy1)
    names = ("ox0","ox1","oy0","oy1","z0","z1","ix0","ix1","iy0","iy1")
    vals = ntuple(i -> _finite_real(raw[i],"box_tunnel_surface",names[i]), 10)
    ox0,ox1,oy0,oy1,z0,z1,ix0,ix1,iy0,iy1 = vals
    (ox0<ix0<ix1<ox1 && oy0<iy0<iy1<oy1 && z0<z1) ||
        throw(ArgumentError("box_tunnel_surface: inner rectangle must be strictly inside the outer"))
    V=Tuple{Float64,Float64,Float64}[]
    for (x,y) in [(ox0,oy0),(ox1,oy0),(ox1,oy1),(ox0,oy1)]; push!(V,(Float64(x),Float64(y),Float64(z0))); end
    for (x,y) in [(ox0,oy0),(ox1,oy0),(ox1,oy1),(ox0,oy1)]; push!(V,(Float64(x),Float64(y),Float64(z1))); end
    for (x,y) in [(ix0,iy0),(ix1,iy0),(ix1,iy1),(ix0,iy1)]; push!(V,(Float64(x),Float64(y),Float64(z0))); end
    for (x,y) in [(ix0,iy0),(ix1,iy0),(ix1,iy1),(ix0,iy1)]; push!(V,(Float64(x),Float64(y),Float64(z1))); end
    Tr=NTuple{3,Int32}[]; q(a,b,c,d)=(push!(Tr,(Int32(a),Int32(b),Int32(c)));push!(Tr,(Int32(a),Int32(c),Int32(d))))
    q(1,9,10,2);q(2,10,11,3);q(3,11,12,4);q(4,12,9,1)      # bottom frame (-z)
    q(5,6,14,13);q(6,7,15,14);q(7,8,16,15);q(8,5,13,16)    # top frame (+z)
    q(1,2,6,5);q(2,3,7,6);q(3,4,8,7);q(4,1,5,8)            # outer walls
    q(9,13,14,10);q(10,14,15,11);q(11,15,16,12);q(12,16,13,9)  # inner tunnel walls
    C=Matrix{Float64}(undef,3,length(V)); for (k,p) in enumerate(V); C[:,k]=[p...]; end
    tm=Matrix{Int32}(undef,3,length(Tr)); for (k,f) in enumerate(Tr); tm[:,k]=Int32[f...]; end
    return _checked_surface(Mesh(C; tris=tm),"box_tunnel_surface")
end

# ── hollow box (axis-aligned Boolean difference: outer box MINUS inner cavity) ──
"""
    box_shell_surface(ox0,ox1, oy0,oy1, oz0,oz1, ix0,ix1, iy0,iy1, iz0,iz1) -> Mesh

Closed boundary surface of the axis-aligned outer box `[ox0,ox1]×[oy0,oy1]×[oz0,oz1]`
MINUS a strictly-interior inner box `[ix0,ix1]×[iy0,iy1]×[iz0,iz1]` — a hollow box
(cavity / "case shell"), the native axis-aligned CSG difference. The inner box must
lie strictly inside the outer one (validated).

The boundary has **two** components: the outer box surface (12 outward triangles,
the [`box_surface`](@ref) pattern) and the inner cavity surface (the same pattern
with **reversed** triangle winding, so its normals face the cavity interior — i.e.
point out of the solid shell). Together they form a closed, manifold, 2-component
surface with zero open edges, ready for [`Mesh3D.tetrahedralize`](@ref) — which
fills the shell and reports the exact volume `outer − inner`.
"""
function box_shell_surface(ox0,ox1,oy0,oy1,oz0,oz1,
                           ix0,ix1,iy0,iy1,iz0,iz1)
    raw = (ox0,ox1,oy0,oy1,oz0,oz1,ix0,ix1,iy0,iy1,iz0,iz1)
    names = ("ox0","ox1","oy0","oy1","oz0","oz1","ix0","ix1","iy0","iy1","iz0","iz1")
    vals = ntuple(i -> _finite_real(raw[i],"box_shell_surface",names[i]), 12)
    ox0,ox1,oy0,oy1,oz0,oz1,ix0,ix1,iy0,iy1,iz0,iz1 = vals
    (ox0<ix0<ix1<ox1 && oy0<iy0<iy1<oy1 && oz0<iz0<iz1<oz1) ||
        throw(ArgumentError("box_shell_surface: inner box must be strictly inside the outer"))
    Co = Float64[ox0 ox1 ox1 ox0 ox0 ox1 ox1 ox0; oy0 oy0 oy1 oy1 oy0 oy0 oy1 oy1; oz0 oz0 oz0 oz0 oz1 oz1 oz1 oz1]
    Ci = Float64[ix0 ix1 ix1 ix0 ix0 ix1 ix1 ix0; iy0 iy0 iy1 iy1 iy0 iy0 iy1 iy1; iz0 iz0 iz0 iz0 iz1 iz1 iz1 iz1]
    C = hcat(Co, Ci)                                     # inner vertices are 9..16
    Fo = [(1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),
          (2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)]   # outer: outward normals
    Fi = [(a+8, c+8, b+8) for (a,b,c) in Fo]             # inner: +8 offset, reversed winding
    F = vcat(Fo, Fi)
    t = Matrix{Int32}(undef,3,length(F)); for (k,f) in enumerate(F); t[:,k]=Int32[f...]; end
    return _checked_surface(Mesh(C; tris=t),"box_shell_surface")
end

@inline _cross(a,b) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline function _unit(a)
    scale=max(abs(a[1]),abs(a[2]),abs(a[3]))
    (isfinite(scale) && scale>0) || throw(ArgumentError(
        "Geometry: axis must have finite positive length"))
    scaled=(a[1]/scale,a[2]/scale,a[3]/scale)
    magnitude=hypot(scaled[1],scaled[2],scaled[3])
    (isfinite(magnitude) && magnitude>0) || throw(ArgumentError(
        "Geometry: axis must have finite positive length"))
    return (scaled[1]/magnitude,scaled[2]/magnitude,scaled[3]/magnitude)
end

@inline function _geometry_count(value,caller::AbstractString,
                                 name::AbstractString,minimum::Int)
    value isa Integer || throw(ArgumentError(
        "$caller: $name must be an integer"))
    value isa Bool && throw(ArgumentError("$caller: $name must not be Bool"))
    converted=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name is outside the platform Int range"))
    end
    converted>=minimum || throw(ArgumentError(
        "$caller: $name must be at least $minimum (got $value)"))
    return converted
end

@inline function _geometry_limit(value,caller::AbstractString,
                                 name::AbstractString)
    return _geometry_count(value,caller,name,1)
end

@inline function _geometry_nonnegative(value,caller::AbstractString,
                                       name::AbstractString)
    result=_finite_real(value,caller,name)
    result>=0 || throw(ArgumentError(
        "$caller: $name must be non-negative (got $value)"))
    return result
end

@inline function _geometry_point(x::Float64,y::Float64,z::Float64,
                                 caller::AbstractString)
    (isfinite(x) && isfinite(y) && isfinite(z)) || throw(ArgumentError(
        "$caller: generated coordinate is not finite"))
    return (x,y,z)
end

function _certify_cylindrical_point(point::NTuple{3,Float64},
                                    origin::NTuple{3,Float64},
                                    axis::NTuple{3,Float64},
                                    expected_radius::Float64,
                                    expected_axial::Float64,
                                    caller::AbstractString)
    displacement=(point[1]-origin[1],point[2]-origin[2],
                  point[3]-origin[3])
    all(isfinite,displacement) || throw(ArgumentError(
        "$caller: generated displacement is not Float64-representable"))
    axis_norm_squared=sum(component*component for component in axis)
    axial=sum(displacement[index]*axis[index] for index in 1:3)/
          axis_norm_squared
    radial=ntuple(index->displacement[index]-axial*axis[index],3)
    radial_distance=hypot(radial...)
    axial_scale=max(abs(expected_axial),expected_radius,nextfloat(0.0))
    radial_scale=expected_radius>0 ? expected_radius : axial_scale
    (isfinite(axial) &&
     abs(axial-expected_axial)<=256eps(Float64)*axial_scale) ||
        throw(ArgumentError(
            "$caller: requested axial position is not representable to round-off"))
    (isfinite(radial_distance) &&
     (expected_radius==0 || radial_distance>0) &&
     abs(radial_distance-expected_radius)<=256eps(Float64)*radial_scale) ||
        throw(ArgumentError(
            "$caller: requested radial position is not representable to round-off"))
    return point
end

function _certify_sphere_point(point::NTuple{3,Float64},
                               center::NTuple{3,Float64},radius::Float64)
    displacement=(point[1]-center[1],point[2]-center[2],
                  point[3]-center[3])
    all(isfinite,displacement) || throw(ArgumentError(
        "sphere_surface: generated displacement is not Float64-representable"))
    represented_radius=hypot(displacement...)
    (isfinite(represented_radius) && represented_radius>0 &&
     abs(represented_radius-radius)<=256eps(Float64)*radius) ||
        throw(ArgumentError(
            "sphere_surface: requested radius is not representable to round-off " *
            "at the supplied center"))
    return point
end

@inline function _axis_frame(axis,caller::AbstractString)
    ez=_unit(_point3(axis,caller,"axis"))
    reference = abs(ez[1])<=abs(ez[2]) ?
        (abs(ez[1])<=abs(ez[3]) ? (1.0,0.0,0.0) : (0.0,0.0,1.0)) :
        (abs(ez[2])<=abs(ez[3]) ? (0.0,1.0,0.0) : (0.0,0.0,1.0))
    ex=_unit(_cross(reference,ez))
    ey=_cross(ez,ex)
    return ex,ey,ez
end

function _checked_primitive_counts(nodes::Int,triangles::Int,max_nodes,
                                   max_triangles,caller::AbstractString)
    node_limit=_geometry_limit(max_nodes,caller,"max_nodes")
    triangle_limit=_geometry_limit(max_triangles,caller,"max_triangles")
    nodes<=typemax(Int32) || throw(ArgumentError(
        "$caller: $nodes nodes exceed the Int32 indexing limit"))
    triangles<=typemax(Int32) || throw(ArgumentError(
        "$caller: $triangles triangles exceed the Int32 topology limit"))
    nodes<=node_limit || throw(ArgumentError(
        "$caller: $nodes nodes exceed max_nodes=$node_limit"))
    triangles<=triangle_limit || throw(ArgumentError(
        "$caller: $triangles triangles exceed max_triangles=$triangle_limit"))
    return nothing
end

function _sphere_midpoint!(vertices::Vector{NTuple{3,Float64}},
                           midpoints::Dict{Tuple{Int32,Int32},Int32},
                           a::Int32,b::Int32,center::NTuple{3,Float64},
                           radius::Float64)
    key=a<b ? (a,b) : (b,a)
    return get!(midpoints,key) do
        pa=vertices[Int(a)];pb=vertices[Int(b)]
        dx=(pa[1]-center[1])/2+(pb[1]-center[1])/2
        dy=(pa[2]-center[2])/2+(pb[2]-center[2])/2
        dz=(pa[3]-center[3])/2+(pb[3]-center[3])/2
        midpoint_norm=hypot(dx,dy,dz)
        (isfinite(midpoint_norm) && midpoint_norm>0) || throw(ArgumentError(
            "sphere_surface: an edge midpoint cannot be projected to the sphere"))
        scale=radius/midpoint_norm
        point=_geometry_point(center[1]+scale*dx,center[2]+scale*dy,
                              center[3]+scale*dz,"sphere_surface")
        _certify_sphere_point(point,center,radius)
        length(vertices)<typemax(Int32) || throw(ArgumentError(
            "sphere_surface: node count exceeds the Int32 indexing limit"))
        push!(vertices,point)
        Int32(length(vertices))
    end
end

"""
    sphere_surface(center, radius; subdivisions=2, max_nodes=10_000_000,
                   max_triangles=20_000_000) -> Mesh

Closed outward-oriented geodesic sphere surface. The six vertices of an octahedron are
projected to the analytical sphere and every subdivision splits one triangle into four,
projecting shared edge midpoints back to the sphere. `subdivisions=0` returns the
octahedron. Resource limits are checked before allocation.
"""
function sphere_surface(center,radius;subdivisions=2,
                        max_nodes=10_000_000,
                        max_triangles=20_000_000)
    r=_finite_real(radius,"sphere_surface","radius")
    r>0 || throw(ArgumentError("sphere_surface: radius must be positive"))
    levels=_geometry_count(subdivisions,"sphere_surface","subdivisions",0)
    power=4
    for _ in 1:levels
        power=try
            Base.checked_mul(power,4)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("sphere_surface: requested subdivision count overflows Int"))
        end
        # Stop before continuing a request that already exceeds either public
        # resource bound; this also bounds work for very large subdivision counts.
        nodes=try Base.checked_add(power,2) catch; typemax(Int) end
        triangles=try Base.checked_mul(2,power) catch; typemax(Int) end
        _checked_primitive_counts(nodes,triangles,max_nodes,max_triangles,
                                  "sphere_surface")
    end
    nodes=Base.checked_add(power,2);triangles=Base.checked_mul(2,power)
    _checked_primitive_counts(nodes,triangles,max_nodes,max_triangles,
                              "sphere_surface")
    c=_point3(center,"sphere_surface","center")

    vertices=NTuple{3,Float64}[
        _geometry_point(c[1]+r,c[2],c[3],"sphere_surface"),
        _geometry_point(c[1]-r,c[2],c[3],"sphere_surface"),
        _geometry_point(c[1],c[2]+r,c[3],"sphere_surface"),
        _geometry_point(c[1],c[2]-r,c[3],"sphere_surface"),
        _geometry_point(c[1],c[2],c[3]+r,"sphere_surface"),
        _geometry_point(c[1],c[2],c[3]-r,"sphere_surface"),
    ]
    foreach(point->_certify_sphere_point(point,c,r),vertices)
    faces=NTuple{3,Int32}[(1,3,5),(1,6,3),(1,5,4),(1,4,6),
                           (2,5,3),(2,3,6),(2,4,5),(2,6,4)]
    sizehint!(vertices,nodes)
    for _ in 1:levels
        midpoints=Dict{Tuple{Int32,Int32},Int32}()
        sizehint!(midpoints,3length(faces)÷2)
        next_faces=NTuple{3,Int32}[]
        sizehint!(next_faces,4length(faces))
        for (a,b,d) in faces
            ab=_sphere_midpoint!(vertices,midpoints,a,b,c,r)
            bd=_sphere_midpoint!(vertices,midpoints,b,d,c,r)
            da=_sphere_midpoint!(vertices,midpoints,d,a,c,r)
            push!(next_faces,(a,ab,da),(ab,b,bd),(da,bd,d),(ab,bd,da))
        end
        faces=next_faces
    end
    length(vertices)==nodes && length(faces)==triangles || throw(ErrorException(
        "sphere_surface: internal subdivision count invariant failed"))
    coordinates=Matrix{Float64}(undef,3,nodes)
    @inbounds for (index,point) in pairs(vertices)
        coordinates[1,index]=point[1];coordinates[2,index]=point[2]
        coordinates[3,index]=point[3]
    end
    topology=Matrix{Int32}(undef,3,triangles)
    @inbounds for (index,face) in pairs(faces)
        topology[1,index]=face[1];topology[2,index]=face[2];topology[3,index]=face[3]
    end
    return _checked_surface(Mesh(coordinates;tris=topology),"sphere_surface")
end

"""
    cone_surface(base, axis, radius1, radius2, height; nθ=24, nz=2,
                 max_nodes=10_000_000, max_triangles=20_000_000) -> Mesh

Closed outward-oriented polygonal cone or conical frustum. `radius1` and `radius2`
are the endpoint radii; either (but not both) can be zero. `nz` is the number of
axial levels, including both endpoints.
"""
function cone_surface(base,axis,radius1,radius2,height;
                      nθ=24,nz=2,
                      max_nodes=10_000_000,
                      max_triangles=20_000_000)
    r1=_geometry_nonnegative(radius1,"cone_surface","radius1")
    r2=_geometry_nonnegative(radius2,"cone_surface","radius2")
    (r1>0 || r2>0) || throw(ArgumentError(
        "cone_surface: at least one endpoint radius must be positive"))
    h=_finite_real(height,"cone_surface","height")
    h>0 || throw(ArgumentError("cone_surface: height must be positive"))
    sectors=_geometry_count(nθ,"cone_surface","nθ",3)
    levels=_geometry_count(nz,"cone_surface","nz",2)
    zero_endpoints=(r1==0 ? 1 : 0)+(r2==0 ? 1 : 0)
    ring_levels=levels-zero_endpoints
    nodes=try
        Base.checked_add(Base.checked_mul(ring_levels,sectors),2)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("cone_surface: requested node count overflows Int"))
    end
    triangles=try
        Base.checked_mul(Base.checked_mul(2,sectors),levels-zero_endpoints)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("cone_surface: requested triangle count overflows Int"))
    end
    _checked_primitive_counts(nodes,triangles,max_nodes,max_triangles,"cone_surface")
    c=_point3(base,"cone_surface","base")
    ex,ey,ez=_axis_frame(axis,"cone_surface")

    coordinates=Matrix{Float64}(undef,3,nodes)
    level_starts=Vector{Int}(undef,levels)
    level_is_apex=falses(levels)
    cursor=1
    @inbounds for level in 1:levels
        fraction=(level-1)/(levels-1)
        radius=level==1 ? r1 : level==levels ? r2 :
            (1-fraction)*r1+fraction*r2
        radius>0 || (r1==0 && level==1) || (r2==0 && level==levels) ||
            throw(ArgumentError(
                "cone_surface: a positive interpolated radius is not " *
                "Float64-representable"))
        axial=fraction*h
        level_starts[level]=cursor
        if radius==0
            point=_geometry_point(c[1]+axial*ez[1],c[2]+axial*ez[2],
                                  c[3]+axial*ez[3],"cone_surface")
            _certify_cylindrical_point(
                point,c,ez,0.0,axial,"cone_surface")
            coordinates[1,cursor]=point[1];coordinates[2,cursor]=point[2]
            coordinates[3,cursor]=point[3]
            level_is_apex[level]=true;cursor+=1
        else
            for sector in 0:sectors-1
                angle=2pi*sector/sectors
                radial_x=cos(angle)*ex[1]+sin(angle)*ey[1]
                radial_y=cos(angle)*ex[2]+sin(angle)*ey[2]
                radial_z=cos(angle)*ex[3]+sin(angle)*ey[3]
                point=_geometry_point(c[1]+axial*ez[1]+radius*radial_x,
                                      c[2]+axial*ez[2]+radius*radial_y,
                                      c[3]+axial*ez[3]+radius*radial_z,
                                      "cone_surface")
                _certify_cylindrical_point(
                    point,c,ez,radius,axial,"cone_surface")
                coordinates[1,cursor]=point[1];coordinates[2,cursor]=point[2]
                coordinates[3,cursor]=point[3];cursor+=1
            end
        end
    end
    bottom_center=0;top_center=0
    if r1>0
        bottom_center=cursor
        _certify_cylindrical_point(c,c,ez,0.0,0.0,"cone_surface")
        coordinates[1,cursor]=c[1];coordinates[2,cursor]=c[2]
        coordinates[3,cursor]=c[3];cursor+=1
    end
    if r2>0
        top_center=cursor
        point=_geometry_point(c[1]+h*ez[1],c[2]+h*ez[2],c[3]+h*ez[3],
                              "cone_surface")
        _certify_cylindrical_point(point,c,ez,0.0,h,"cone_surface")
        coordinates[1,cursor]=point[1];coordinates[2,cursor]=point[2]
        coordinates[3,cursor]=point[3];cursor+=1
    end
    cursor==nodes+1 || throw(ErrorException(
        "cone_surface: internal node count invariant failed"))

    topology=Matrix{Int32}(undef,3,triangles);face=1
    @inbounds for level in 1:levels-1
        lower=level_starts[level];upper=level_starts[level+1]
        if level_is_apex[level]
            for sector in 0:sectors-1
                current=upper+sector;next=upper+mod(sector+1,sectors)
                topology[1,face]=Int32(lower);topology[2,face]=Int32(next)
                topology[3,face]=Int32(current);face+=1
            end
        elseif level_is_apex[level+1]
            for sector in 0:sectors-1
                current=lower+sector;next=lower+mod(sector+1,sectors)
                topology[1,face]=Int32(current);topology[2,face]=Int32(next)
                topology[3,face]=Int32(upper);face+=1
            end
        else
            for sector in 0:sectors-1
                a=lower+sector;b=lower+mod(sector+1,sectors)
                d=upper+sector;e=upper+mod(sector+1,sectors)
                topology[1,face]=Int32(a);topology[2,face]=Int32(b)
                topology[3,face]=Int32(e);face+=1
                topology[1,face]=Int32(a);topology[2,face]=Int32(e)
                topology[3,face]=Int32(d);face+=1
            end
        end
    end
    if bottom_center!=0
        lower=level_starts[1]
        @inbounds for sector in 0:sectors-1
            current=lower+sector;next=lower+mod(sector+1,sectors)
            topology[1,face]=Int32(bottom_center);topology[2,face]=Int32(next)
            topology[3,face]=Int32(current);face+=1
        end
    end
    if top_center!=0
        upper=level_starts[end]
        @inbounds for sector in 0:sectors-1
            current=upper+sector;next=upper+mod(sector+1,sectors)
            topology[1,face]=Int32(top_center);topology[2,face]=Int32(current)
            topology[3,face]=Int32(next);face+=1
        end
    end
    face==triangles+1 || throw(ErrorException(
        "cone_surface: internal triangle count invariant failed"))
    return _checked_surface(Mesh(coordinates;tris=topology),"cone_surface")
end

function _checked_surface(m::Mesh,caller::AbstractString)
    d=validate(m)
    if !d.ok
        numerical=all(message->
            occursin("non-finite computed area",message) ||
            occursin("degenerate (zero-area) triangles",message),d.messages)
        numerical && throw(ArgumentError(
            "$caller: input geometry is not representable as a valid " *
            "Float64 surface — "*join(d.messages,"; ")))
        throw(ErrorException(
            "$caller: constructed surface is invalid — "*join(d.messages,"; ")))
    end
    ok,report=is_meshable(m)
    if !ok
        numerical=report.closed && report.manifold && report.oriented &&
            report.n_duplicate_tris==0 &&
            (report.n_degenerate_tris>0 || report.n_coincident_pairs>0)
        numerical && throw(ArgumentError(
            "$caller: input geometry is not representable as a meshable " *
            "Float64 surface — "*join(report.messages,"; ")))
        throw(ErrorException(
            "$caller: constructed surface is not meshable — "*
            join(report.messages,"; ")))
    end
    return m
end

end # module Geometry
