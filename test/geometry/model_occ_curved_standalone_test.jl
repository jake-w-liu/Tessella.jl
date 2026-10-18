# Standalone `mesh_model_surface` on parametric OCC faces.
#
# `_model_occ_curved_surface_mesh` tessellates a curved face's (u,v) domain —
# but previously ran only inside a volume PLC, so `mesh_model_surface` (and
# `generate(2)`) rejected every curved OCC face at the planar gate, while
# upstream `meshGFace` takes the parametric path on OCC faces regardless of
# whether the face bounds a volume. `mesh_model_surface` now routes
# `_model_occ_curved_surface` faces — analytic Cylinder/Sphere/Cone/Torus
# records and tensor patches carrying edge pcurves (trimmed NURBS, lofts,
# materialized shells) — through the same kernel, honoring the shared
# Mesh.Algorithm admission. Boundary nodes keep their canonical
# `_periodic_curve_point` coordinates so neighbouring faces still weld
# bitwise; no boundary writeback is needed (or possible — a seam-welded
# patch has no conventional edge chain). uv-space embedding, recombination,
# smoothing and background size fields have no analogue in the lifted
# construction and reject loudly; so do faces bounded by degenerate
# pole/apex edges, whose collapsed uv rows need a pole fan the kernel does
# not build yet.

using Test
using Tessella
const M=Tessella.Model
const MT=Tessella.MeshTypes

# Every mesh node satisfies the analytic cylinder band |ρ|=r, z0≤z≤z1.
function _on_cylinder(mesh,r,z0,z1;atol=1e-9)
    for j in axes(mesh.coords,2)
        abs(hypot(mesh.coords[1,j],mesh.coords[2,j])-r)<=atol || return false
        z0-atol<=mesh.coords[3,j]<=z1+atol || return false
    end
    return true
end

_coord_bits(mesh,j)=ntuple(i->reinterpret(UInt64,mesh.coords[i,j]),3)
_coord_set(mesh)=Set(_coord_bits(mesh,j) for j in axes(mesh.coords,2))

_error(f)=try;f();nothing catch err;err end

@testset "standalone curved OCC faces mesh in their (u,v) domain" begin
    # Materialized cylinder: the lateral is a `:cylinder` face with rim/seam
    # pcurves; the caps stay planar. Standalone meshing used to reject the
    # lateral with "Surface[1] is a Cylinder".
    m=GeoModel()
    M._model_occ_add_cylinder!(m,0.0,0.0,0.0,0.0,0.0,2.0,1.0;tag=1)
    M._model_occ_synchronize!(m)
    barrel=findfirst(==(:cylinder),m.surface_types)
    mesh=M.mesh_model_surface(m,barrel)
    @test MT.validate(mesh).ok
    @test MT.ntris(mesh)>0
    @test _on_cylinder(mesh,1.0,0.0,2.0)
    # The caps still take the planar path and weld the shared rim nodes
    # bitwise — `_periodic_curve_point` canonical coordinates on both sides.
    cap=findfirst(==(:plane),m.surface_types)
    cmesh=M.mesh_model_surface(m,cap)
    @test MT.validate(cmesh).ok
    shared=intersect(_coord_set(mesh),_coord_set(cmesh))
    @test length(shared)>=4          # rim circle births on both faces

    # A half-sector lateral (`angle=π`) carries a trimmed uv window with
    # certified seam edges — the same kernel the volume PLC uses.
    ms=GeoModel()
    M._model_occ_add_cylinder!(ms,0.0,0.0,0.0,0.0,0.0,2.0,1.0;
                               tag=1,angle=pi)
    M._model_occ_synchronize!(ms)
    sector=findfirst(==(:cylinder),ms.surface_types)
    smesh=M.mesh_model_surface(ms,sector)
    @test MT.validate(smesh).ok
    @test MT.ntris(smesh)>0
    @test _on_cylinder(smesh,1.0,0.0,2.0)
    # The θ∈[0,π] sector spans x∈[-r,r] but stays on the y≥0 half.
    for j in axes(smesh.coords,2)
        @test smesh.coords[2,j]>=-1e-9
    end

    # NURBS-converted lateral: `convertToNURBS` rebinds the face to a
    # quadratic-rational patch while the public (angle,height) parameters —
    # and every stored pcurve — stay bitwise. The kernel's `_occ_surface_point`
    # lift evaluates the patch through `_model_nurbs_surface_d2`.
    mn=GeoModel()
    M._model_occ_add_cylinder!(mn,0.0,0.0,0.0,0.0,0.0,2.0,1.0;tag=1)
    M._model_occ_convert_to_nurbs!(mn,[(3,1)],"test")
    M._model_occ_synchronize!(mn)
    nlat=findfirst(==(:cylinder),mn.surface_types)
    ng=mn.surface_geometry[nlat]
    @test hasproperty(ng,:patch) && hasproperty(ng,:reparam)
    nmesh=M.mesh_model_surface(mn,nlat)
    @test MT.validate(nmesh).ok
    @test MT.ntris(nmesh)>0
    @test _on_cylinder(nmesh,1.0,0.0,2.0)

    # A materialized sphere face is bounded by degenerate pole edges whose
    # pcurves span an entire collapsed uv row — an explicit pole fan is not
    # built yet, so the kernel rejects loudly rather than feeding a
    # self-duplicated PSLG to `constrained_delaunay` (the same rejection now
    # covers a sewn sphere shell inside a volume PLC).
    mp=GeoModel()
    M._model_occ_add_sphere!(mp,0.0,0.0,0.0,1.5;tag=1)
    M._model_occ_synchronize!(mp)
    sface=findfirst(==(:sphere),mp.surface_types)
    err=_error(()->M.mesh_model_surface(mp,sface))
    @test err isa ArgumentError
    @test occursin("degenerate",sprint(showerror,err))
    @test occursin("pole",sprint(showerror,err))
end

@testset "trimmed cylinder band meshes standalone" begin
    # `addTrimmedSurface` on a curved base: a quad wire on the barrel's
    # front patch (θ∈[0,½], z∈[2,2.5]) is projected edge-wise — the two
    # chords become spl2d pcurves and the two axial lines lin2d — and the
    # record is marked `trimmed` with `vertex_uv`. Constructed at the
    # kernel/Model level so the check also runs under minimal compilation
    # (the public arc entry points route through an llvmcall intrinsic).
    m=GeoModel()
    M._model_occ_add_cylinder!(m,0.0,0.0,0.0,0.0,0.0,4.0,1.0;tag=1)
    M._model_occ_synchronize!(m)
    barrel=findfirst(==(:cylinder),m.surface_types)
    h=0.5
    a=M._model_occ_add_point!(m,1.0,0.0,2.0)
    b=M._model_occ_add_point!(m,cos(h),sin(h),2.0)
    c=M._model_occ_add_point!(m,1.0,0.0,2.5)
    d=M._model_occ_add_point!(m,cos(h),sin(h),2.5)
    l1=M._model_occ_add_line!(m,a,b)
    l2=M._model_occ_add_line!(m,b,d)
    l3=M._model_occ_add_line!(m,d,c)
    l4=M._model_occ_add_line!(m,c,a)
    w=M._model_occ_add_wire!(m,[l1,l2,l3,l4];check_closed=true)
    ft=M._model_occ_trimmed_surface(m,-1,w,barrel)
    M._model_occ_synchronize!(m)
    sg=m.surface_geometry[ft]
    @test get(sg,:trimmed,false) && hasproperty(sg,:pcurves)
    mesh=M.mesh_model_surface(m,ft)
    @test MT.validate(mesh).ok
    @test MT.ntris(mesh)>0
    # Every emitted node lies on the barrel inside the trimmed band.
    @test _on_cylinder(mesh,1.0,2.0,2.5;atol=1e-8)
end

@testset "unsupported surface attributes on curved faces reject loudly" begin
    function _sector()
        m=GeoModel()
        M._model_occ_add_cylinder!(m,0.0,0.0,0.0,0.0,0.0,2.0,1.0;
                                   tag=1,angle=pi)
        M._model_occ_synchronize!(m)
        return m,findfirst(==(:cylinder),m.surface_types)
    end
    # Mesh.Algorithm admission is shared with the volume path: 7 (BAMG) and
    # the other producer-less selectors reject on curved faces.
    m,lat=_sector()
    M.set_algorithm!(m,2,lat,7)
    err=_error(()->M.mesh_model_surface(m,lat))
    @test err isa ArgumentError
    @test occursin("not implemented",sprint(showerror,err))
    m,lat=_sector()
    m.meshing.algorithm_2d=2   # AUTO resolves to MeshAdapt off-plane
    err=_error(()->M.mesh_model_surface(m,lat))
    @test err isa ArgumentError
    @test occursin("not implemented",sprint(showerror,err))
    # Recombination has no uv analogue in the lifted construction.
    m,lat=_sector()
    m.meshing.recombine_all=true
    err=_error(()->M.mesh_model_surface(m,lat))
    @test err isa ArgumentError
    @test occursin("recombination",sprint(showerror,err))
    # Laplacian smoothing operates on pre-lift plane coordinates only.
    m,lat=_sector()
    M.set_smoothing!(m,2,lat,3)
    err=_error(()->M.mesh_model_surface(m,lat))
    @test err isa ArgumentError
    @test occursin("smoothing",sprint(showerror,err))
    # uv-space embedding is not implemented for the parametric kernel.
    m,lat=_sector()
    m.embeds[(2,lat)]=NTuple{2,Int}[(0,1)]
    err=_error(()->M.mesh_model_surface(m,lat))
    @test err isa ArgumentError
    @test occursin("embedded",sprint(showerror,err))
end
