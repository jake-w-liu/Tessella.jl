using Test, Tessella, SHA, TOML

const _AC_B4_MODEL=Tessella.Model
const _AC_B4_FIXTURE=joinpath(@__DIR__,"..","data","native_geo_auto_coherence_v1",
                            "native_geometry_controls_v3.toml")

@testset "Unmerged B4 operation provenance preserves literal CAD" begin
    @test bytes2hex(sha256(read(_AC_B4_FIXTURE)))==
        "991451ac2885f4f647b93abb65f4ff3b79948d35ad48b121c79ceac3a61e5b94"
    for case in TOML.parsefile(_AC_B4_FIXTURE)["cases"]
        case["name"] in ("original_b4_free_mode0","original_b4_recombined_mode0") || continue
        @testset "$(case["name"])" begin
            model=mktemp() do path,io
                write(io,case["source"]);close(io)
                execute_geo(path;mesh_dim=0).model
            end
            volume=only(keys(model.volumes))
            literal=_AC_B4_MODEL._model_volume_boundary_surfaces(model,volume,
                "literal B4 incidence";certify_shells=false)
            before=deepcopy(model)
            @test_throws ArgumentError _AC_B4_MODEL._model_volume_boundary_surfaces(
                model,volume,"strict B4 CAD shell")
            @test _AC_B4_MODEL._extrude_nonew_boundary_surfaces(model,volume,
                "recorded B4 operation")==literal
            @test length(model.points)==44 && length(model.curves)==20
            @test isequal((model.points,model.curves,model.loops,model.surfaces,
                model.surface_loops,model.volumes,model.meshing.extrude_sources,
                model.meshing.extrude_specs,model.next_tag),
                (before.points,before.curves,before.loops,before.surfaces,
                before.surface_loops,before.volumes,before.meshing.extrude_sources,
                before.meshing.extrude_specs,before.next_tag))
            source=abs(model.meshing.extrude_sources[(3,volume)][2])
            top,laterals=_AC_B4_MODEL._extrude_nonew_region_boundary(
                model,volume,source,"recorded B4 operation")
            cap_curve=abs(first(model.loops[only(model.surfaces[top])]))
            lateral=first(laterals)
            lateral_curve=only(curve for curve in
                abs.(model.loops[only(model.surfaces[lateral])]) if
                get(model.meshing.extrude_sources,(1,curve),nothing)!==nothing &&
                model.meshing.extrude_sources[(1,curve)][1]==1)
            mutations=(
                m->delete!(m.meshing.extrude_sources,(2,top)),
                m->(m.meshing.extrude_sources[(2,top)]=(2,source+1000)),
                m->(m.meshing.extrude_sources[(2,lateral)]=m.meshing.extrude_sources[(2,last(laterals))]),
                m->(m.meshing.extrude_specs[(2,top)]=merge(m.meshing.extrude_specs[(2,top)],(T=(0.0,0.0,2.0),))),
                m->(m.meshing.extrude[(2,lateral)]=merge(m.meshing.extrude[(2,lateral)],(layers=[2],))),
                m->delete!(m.meshing.extrude_sources,(1,cap_curve)),
                m->(m.meshing.extrude_sources[(1,cap_curve)]=(1,1000)),
                m->(m.meshing.extrude_sources[(1,lateral_curve)]=(1,1000)),
                m->delete!(m.curves,cap_curve),
                m->push!(m.surface_loops[only(m.volumes[volume])],top),
                m->(m.surface_loops[only(m.volumes[volume])][2]=-top),
            )
            for mutate in mutations
                rejected=deepcopy(model);mutate(rejected)
                snapshot=deepcopy((rejected.points,rejected.curves,rejected.loops,
                    rejected.surfaces,rejected.surface_loops,rejected.volumes,
                    rejected.meshing.extrude_sources,rejected.meshing.extrude_specs,
                    rejected.meshing.extrude,rejected.next_tag))
                @test_throws ArgumentError _AC_B4_MODEL._extrude_nonew_boundary_surfaces(
                    rejected,volume,"corrupt recorded B4 operation")
                @test isequal(snapshot,(rejected.points,rejected.curves,rejected.loops,
                    rejected.surfaces,rejected.surface_loops,rejected.volumes,
                    rejected.meshing.extrude_sources,rejected.meshing.extrude_specs,
                    rejected.meshing.extrude,rejected.next_tag))
            end
            coherent=mktemp() do path,io
                write(io,replace(case["source"],"Geometry.AutoCoherence=0;"=>
                    "Geometry.AutoCoherence=1;"));close(io)
                execute_geo(path;mesh_dim=0).model
            end
            strict=_AC_B4_MODEL._model_volume_boundary_surfaces(
                coherent,volume,"coherent B4 CAD shell")
            @test _AC_B4_MODEL._extrude_nonew_boundary_surfaces(
                coherent,volume,"coherent B4 strict admission")==strict
            coherent_top=abs(coherent.meshing.extrude_sources[(3,volume)][2])
            cap=only(surface for surface in abs.(strict) if
                get(coherent.meshing.extrude_sources,(2,surface),nothing)==
                coherent.meshing.extrude_sources[(3,volume)])
            for curve in abs.(coherent.loops[only(coherent.surfaces[cap])])
                _AC_B4_MODEL.model_set_tag!(coherent,1,curve,1000+curve)
            end
            @test _AC_B4_MODEL._model_volume_boundary_surfaces(
                coherent,volume,"retagged coherent B4")==strict
            @test _AC_B4_MODEL._extrude_nonew_boundary_surfaces(
                coherent,volume,"retagged coherent B4")==strict
        end
    end
end
