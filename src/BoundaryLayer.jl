"""
    BoundaryLayer

First-order prismatic boundary-layer extrusion of a triangle surface along
area-weighted vertex normals. The result is a mixed mesh of type-6 prisms.
This is the element-topology counterpart of [`BoundaryLayerField`](@ref).
"""
module BoundaryLayer

using ..MeshTypes: Mesh, nnodes, ntris, triangle_area, tet_volume
using ..Elements: ElementBlock, MixedMesh, validate
using ..Predicates: orient3

export mesh_boundary_layer

function _finite(value, caller, name)
    v=try Float64(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name must be Float64-representable"))
    end
    isfinite(v) || throw(ArgumentError("$caller: $name must be finite"))
    return v
end

function mesh_boundary_layer(surface::Mesh; hwall::Real, ratio::Real, nlayers::Integer,
                             max_prisms::Integer=10_000_000)
    caller="mesh_boundary_layer"
    ntris(surface)>0 || throw(ArgumentError("$caller: surface has no triangles"))
    hw=_finite(hwall,caller,"hwall"); hw>0 || throw(ArgumentError("$caller: hwall must be positive"))
    ra=_finite(ratio,caller,"ratio"); ra>1 || throw(ArgumentError("$caller: ratio must be > 1"))
    nl=Int(nlayers); nl>=1 || throw(ArgumentError("$caller: nlayers must be ≥ 1"))
    npr=nl*ntris(surface)
    npr<=max_prisms || throw(ArgumentError("$caller: $npr prisms exceed max_prisms=$max_prisms"))
    npr<=typemax(Int32) || throw(ArgumentError("$caller: prism count exceeds Int32"))

    nv=nnodes(surface)
    normals=zeros(Float64,3,nv)
    @inbounds for t in 1:ntris(surface)
        i,j,k=Int(surface.tris[1,t]),Int(surface.tris[2,t]),Int(surface.tris[3,t])
        a=(surface.coords[1,i],surface.coords[2,i],surface.coords[3,i])
        b=(surface.coords[1,j],surface.coords[2,j],surface.coords[3,j])
        c=(surface.coords[1,k],surface.coords[2,k],surface.coords[3,k])
        ab=(b[1]-a[1],b[2]-a[2],b[3]-a[3])
        ac=(c[1]-a[1],c[2]-a[2],c[3]-a[3])
        n=(ab[2]*ac[3]-ab[3]*ac[2], ab[3]*ac[1]-ab[1]*ac[3], ab[1]*ac[2]-ab[2]*ac[1])
        area=triangle_area(a,b,c)
        for id in (i,j,k)
            normals[1,id]+=n[1]*area; normals[2,id]+=n[2]*area; normals[3,id]+=n[3]*area
        end
    end
    @inbounds for i in 1:nv
        L=hypot(normals[1,i],normals[2,i],normals[3,i])
        L>0 || throw(ArgumentError("$caller: vertex $i has a zero normal"))
        normals[1,i]/=L; normals[2,i]/=L; normals[3,i]/=L
    end

    # Geometric offsets: h_k = hwall * (ratio^k - 1)/(ratio-1), k=1..nl
    offsets=Vector{Float64}(undef,nl)
    @inbounds for k in 1:nl
        offsets[k]=hw*(ra^k-1)/(ra-1)
        isfinite(offsets[k]) || throw(ArgumentError("$caller: layer offset overflowed"))
    end

    nout=nv*(nl+1)
    nout<=typemax(Int32) || throw(ArgumentError("$caller: node count exceeds Int32"))
    coords=Matrix{Float64}(undef,3,nout)
    @inbounds for i in 1:nv
        coords[1,i]=surface.coords[1,i]; coords[2,i]=surface.coords[2,i]; coords[3,i]=surface.coords[3,i]
    end
    @inbounds for k in 1:nl, i in 1:nv
        id=k*nv+i
        coords[1,id]=surface.coords[1,i]+offsets[k]*normals[1,i]
        coords[2,id]=surface.coords[2,i]+offsets[k]*normals[2,i]
        coords[3,id]=surface.coords[3,i]+offsets[k]*normals[3,i]
        all(isfinite, (coords[1,id],coords[2,id],coords[3,id])) ||
            throw(ArgumentError("$caller: extruded node is non-finite"))
    end

    prisms=Matrix{Int32}(undef,6,npr)
    cursor=0
    @inbounds for k in 0:nl-1, t in 1:ntris(surface)
        cursor+=1
        b1=Int32(k*nv+Int(surface.tris[1,t]))
        b2=Int32(k*nv+Int(surface.tris[2,t]))
        b3=Int32(k*nv+Int(surface.tris[3,t]))
        t1=Int32((k+1)*nv+Int(surface.tris[1,t]))
        t2=Int32((k+1)*nv+Int(surface.tris[2,t]))
        t3=Int32((k+1)*nv+Int(surface.tris[3,t]))
        prisms[:,cursor].=(b1,b2,b3,t1,t2,t3)
    end
    mesh=MixedMesh(coords,[ElementBlock(6,prisms)])
    diag=validate(mesh)
    diag.ok || throw(ErrorException("$caller: invalid mesh — "*join(diag.messages,"; ")))
    return mesh
end

end # module
