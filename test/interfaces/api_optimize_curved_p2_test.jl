# Curved quadratic-simplex `API.mesh.optimize` dispatch against Gmsh 4.15.2's
# method gates (src/mesh/Generator.cpp + meshSmoothLaplacian.cpp +
# meshRelocateVertex.cpp), oracle-verified on an OCC cylinder fixture:
#
#   * default / "Gmsh" / "Optimize" visit only region tetrahedra — a
#     surface-only P2 mesh is accepted and unchanged;
#   * "Laplace2D" refuses every order>1 face with upstream's first-order
#     diagnostic — the check precedes its niter early-return, so the error
#     also fires at niter=0;
#   * "RelocateVertices" has no order gate — "Relocate2D" relocates
#     face-owned mid nodes in the owning surface's (u,v) parameters;
#   * "HighOrder"/"HighOrderElastic" run the staged interior relocation —
#     "HighOrderElastic" the unbounded elastic stage, "HighOrder" the
#     Mesh.HighOrderThresholdMin/Max-banded stage;
#   * "Relocate3D" ingests the actual Tet10/Tet4 volume cells unchanged.
#
# Every accepted arm publishes atomically; a rejected arm leaves the cache
# untouched.

module APIOptimizeCurvedP2Tests
using Test,Tessella
const api=Tessella.API

function _oc_nodes()
    tags,flat,_=api.mesh.get_nodes()
    xyz=reshape(flat,3,:)
    return Dict(tag=>Tuple(xyz[:,i]) for (i,tag) in enumerate(tags))
end

function _oc_cylinder(;order=true)
    api.initialize()
    try
        occ=api.model.occ
        occ.add_cylinder(0.,0.,0., 0.,0.,2., 1.)
        occ.synchronize()
        api.mesh.generate(2)
        order && api.mesh.set_order(2)
        return _oc_nodes()
    catch
        api.finalize();rethrow()
    end
end

@testset "Curved P2 surface optimizer dispatch" begin
    before=_oc_cylinder()
    try
        @test api.LAST_MESH[] isa Tessella.Elements.MixedMesh

        # Default and Relocate3D visit regions only — no region cells on a
        # surface mesh, so the mesh is accepted and stays bitwise unchanged.
        for method in ("","Gmsh","Optimize","Relocate3D")
            @test api.mesh.optimize(method)===nothing
            @test _oc_nodes()==before
        end

        # niter=0 Relocate2D returns without work.
        @test api.mesh.optimize("Relocate2D",false,0)===nothing
        @test _oc_nodes()==before

        # Relocate2D relocates face-owned mid nodes parametrically (upstream
        # moved 249 of the same class on the denser oracle fixture).
        @test api.mesh.optimize("Relocate2D")===nothing
        relocated=_oc_nodes()
        moved=count(tag->relocated[tag]!=before[tag],keys(before))
        @test moved>0
        # Every node stays on its owning surface: corner tags keep positions;
        # moved mids keep their (u,v) parametric support.
        for tag in keys(before)
            p=relocated[tag]
            @test isapprox(hypot(p[1],p[2]),1.0,atol=2e-2) ||
                  isapprox(p[3],0.0,atol=2e-2) ||
                  isapprox(p[3],2.0,atol=2e-2)
        end

        # Laplace2D's upstream order gate fires even at niter=0 and publishes
        # nothing.
        for niter in (1,0)
            @test_throws r"only valid for first order mesh" (
                api.mesh.optimize("Laplace2D",false,niter))
            @test _oc_nodes()==relocated
        end

        # The dedicated upstream method names run the staged interior
        # relocation through the same certificate-gated engine.
        @test api.mesh.optimize("HighOrderElastic")===nothing
        elastic=_oc_nodes()
        @test count(tag->elastic[tag]!=relocated[tag],keys(before))>0
        @test api.mesh.optimize("HighOrder")===nothing
        @test api.mesh.optimize("HighOrderElastic",false,0)===nothing

        # Entity selection stays unsupported-but-tolerated, like upstream's
        # "not interfaced yet" warning.
        @test api.mesh.optimize("HighOrder",false,-1,[(2,1)])===nothing

        # "Netgen" is a supported region optimizer: a surface-only cache has
        # no regions, so the call is accepted without work (upstream iterates
        # model regions; there are none to optimize).
        @test api.mesh.optimize("Netgen")===nothing
        @test _oc_nodes()==relocated
        # Unsupported kernels keep their atomic rejection.
        @test_throws r"unknown or unsupported optimizer" (
            api.mesh.optimize("Pachyderm"))
    finally
        api.finalize()
    end
end

@testset "Linear MixedMesh surface optimizer dispatch" begin
    before=_oc_cylinder(;order=false)
    try
        @test api.LAST_MESH[] isa Tessella.Elements.MixedMesh
        # Linear MixedMesh Laplace2D/Relocate2D run upstream's per-face
        # parametric surface kernels — relocated face-owned vertices keep
        # their owning-surface support (lateral r=1 or cap z∈{0,2}).
        for method in ("Laplace2D","Relocate2D")
            @test api.mesh.optimize(method)===nothing
        end
        relocated=_oc_nodes()
        @test keys(relocated)==keys(before)
        for tag in keys(before)
            p=relocated[tag]
            @test isapprox(hypot(p[1],p[2]),1.0,atol=2e-2) ||
                  isapprox(p[3],0.0,atol=2e-2) ||
                  isapprox(p[3],2.0,atol=2e-2)
        end
        # Upstream warns and no-ops when the dedicated optimizer meets a
        # purely linear mesh.
        @test api.mesh.optimize("HighOrder")===nothing
        @test api.mesh.optimize("HighOrderElastic")===nothing
        @test _oc_nodes()==relocated
    finally
        api.finalize()
    end
end
end
