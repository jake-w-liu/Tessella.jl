module RetainedSplitMacrofaceCertificateTests
using Test,Tessella,TOML
import Tessella.Elements: nnodes
import Tessella.Model: _extrude_nonew_certify_boundary,_extrude_nonew_split_macroface_certificate
const FIXTURE=joinpath(@__DIR__,"..","data","native_geo_auto_coherence_v1",
    "native_geometry_controls_v3.toml")

@testset "Retained B4 split faces preserve surface coverage" begin
    api=Tessella.API;model_api=Tessella.Model
    source=only(filter(row->row["name"]=="original_b4_recombined_mode0",
        TOML.parsefile(FIXTURE)["cases"]))["source"]
    source=replace(source,"Layers{500}"=>"Layers{2}")
    api.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"retained_scope.geo");write(path,source)
            api.open_geo!(path;mesh_dim=0);api.mesh.generate(2)
            model=api.CURRENT[];volume=only(keys(model.volumes))
            scope=model_api._extrude_nonew_scope(model,"current retained receipt";
                _mesh_input=model.meshing.mesh_authority)
            completed=scope.volumes[volume]
            @test model_api._extrude_nonew_retained_b4_catalog(completed.catalog)
            @test nnodes(completed.volume)==36
            @test isempty(completed.catalog.problem_positions)
            proof=model_api._extrude_nonew_retained_b4_certificate(model,volume,
                completed.volume,completed.surfaces,completed.catalog,"current retained receipt")
            @test Set(values(proof.facet_surfaces))==Set(keys(completed.surfaces))
            surface_signs=Dict(abs(tag)=>sign(tag) for tag in
                model_api._model_volume_boundary_surfaces(model,volume,"retained split signs"))
            @testset "Split surface faces preserve the volume boundary" begin
                let volume=completed.volume,surfaces=completed.surfaces,signs=surface_signs
                    before=deepcopy((volume,surfaces,signs))
                    @test_throws ArgumentError _extrude_nonew_certify_boundary(volume,surfaces,"original typed default";oriented_internal=true)
                    proof=_extrude_nonew_split_macroface_certificate(volume,surfaces,signs,"native retained proof")
                    @test length(proof.boundary)==length(proof.surface_polygons)+length(proof.relations)
                    @test length(proof.surface_polygons)==39
                    @test !isempty(proof.relations)
                    @test Set(keys(proof.facet_surfaces))==proof.boundary
                    @test all(length(unique(vcat(collect(left[1:3]),collect(right[1:3]))))==4 for (left,right) in values(proof.relations))
                    order=collect(nnodes(volume):-1:1)
                    inverse=invperm(order)
                    permuted=MixedMesh(volume.coords[:,order],[ElementBlock(block.msh,
                        Int32.(inverse[block.nodes])) for block in volume.blocks])
                    reordered=_extrude_nonew_split_macroface_certificate(permuted,surfaces,signs,
                        "retained proof binds actual reordered Region rows")
                    @test length(reordered.boundary)==length(proof.boundary) && length(reordered.relations)==length(proof.relations)
                    @test Set(keys(reordered.facet_surfaces))==reordered.boundary
                    @test volume.coords==before[1].coords && all(a.nodes==b.nodes for (a,b) in zip(volume.blocks,before[1].blocks))
                    @test all(surfaces[tag].coords==part.coords && all(a.nodes==b.nodes for (a,b) in zip(surfaces[tag].blocks,part.blocks)) for (tag,part) in before[2])
                    @test signs==before[3]
                    # Rejected inputs leave their meshes and orientations unchanged.
                    for mutation in (
                        (v,s,g)->(s[26]=MixedMesh(s[26].coords,[ElementBlock(2,s[26].blocks[1].nodes[:,2:end])])),
                        (v,s,g)->(g[13]=-g[13]))
                        v,s,g=deepcopy((volume,surfaces,signs));mutation(v,s,g)
                        snapshot=deepcopy((v,s,g))
                        @test_throws ArgumentError _extrude_nonew_split_macroface_certificate(v,s,g,"malformed native retained proof")
                        @test v.coords==snapshot[1].coords && all(a.nodes==b.nodes for (a,b) in zip(v.blocks,snapshot[1].blocks))
                        @test all(s[tag].coords==part.coords && all(a.nodes==b.nodes for (a,b) in zip(s[tag].blocks,part.blocks)) for (tag,part) in snapshot[2])
                        @test g==snapshot[3]
                    end
                end
            end
            before=deepcopy((api.mesh.get_nodes(),api.mesh.get_elements(),
                api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag(),
                model.points,model.curves,model.loops,model.surfaces,model.volumes,model.next_tag))
            for change in (
                m->(m.meshing.extrude_sources[(3,volume)]=(2,999)),
                m->delete!(m.meshing.extrude_sources,(2,13)))
                invalid=deepcopy(model);change(invalid)
                @test_throws ArgumentError model_api._extrude_nonew_retained_b4_certificate(
                    invalid,volume,completed.volume,completed.surfaces,completed.catalog,
                    "stale retained operation")
            end
            @test isequal(before,(api.mesh.get_nodes(),api.mesh.get_elements(),
                api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag(),
                model.points,model.curves,model.loops,model.surfaces,model.volumes,model.next_tag))
        end
    finally
        api.finalize()
    end
end
end
