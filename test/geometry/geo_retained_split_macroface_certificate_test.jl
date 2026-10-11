module RetainedSplitMacrofaceCertificateTests
using Test,Tessella,SHA,TOML
import Tessella.Elements: nnodes
import Tessella.Model: _extrude_nonew_certify_boundary,_extrude_nonew_split_macroface_certificate
const FIXTURE=joinpath(@__DIR__,"..","data","native_geo_auto_coherence_v1","native_retained_recombined.toml")
const WIDTH=Dict(2=>3,3=>4,4=>4,5=>8,6=>6,7=>5)
function native_mesh(row,points)
    used=sort!(unique(vcat(row["connectivity"]...)))
    lookup=Dict(uid=>Int32(index) for (index,uid) in pairs(used))
    coordinates=hcat((points[Int(uid)] for uid in used)...)
    blocks=ElementBlock[]
    for (kind,connectivity) in zip(row["element_types"],row["connectivity"])
        push!(blocks,ElementBlock(kind,reshape(Int32[lookup[uid] for uid in connectivity],WIDTH[kind],:)))
    end
    return MixedMesh(coordinates,blocks)
end
@testset "Actual retained split-macroface proof preserves default typed certificate" begin
    native=TOML.parsefile(FIXTURE)
    coordinates_vector=[reinterpret(Float64,hex2bytes(value))[1] for value in native["global_xyz_bits"]]
    points=[coordinates_vector[index:index+2] for index in 1:3:length(coordinates_vector)]
    region=only(row for row in native["entities"] if row["dimension"]==3)
    volume=native_mesh(region,points)
    surfaces=Dict(Int(row["tag"])=>native_mesh(row,points) for row in native["entities"] if row["dimension"]==2)
    signs=Dict(abs(Int(tag))=>sign(Int(tag)) for (dimension,tag) in region["boundary"])
    before=deepcopy((volume,surfaces,signs))
    @test_throws ArgumentError _extrude_nonew_certify_boundary(volume,surfaces,"original typed default";oriented_internal=true)
    proof=_extrude_nonew_split_macroface_certificate(volume,surfaces,signs,"native retained proof")
    @test length(proof.boundary)==6025
    @test length(proof.surface_polygons)==6015
    @test length(proof.relations)==10
    @test Set(keys(proof.facet_surfaces))==proof.boundary
    @test all(length(unique(vcat(collect(left[1:3]),collect(right[1:3]))))==4 for (left,right) in values(proof.relations))
    order=collect(nnodes(volume):-1:1)
    inverse=invperm(order)
    permuted=MixedMesh(volume.coords[:,order],[ElementBlock(block.msh,
        Int32.(inverse[block.nodes])) for block in volume.blocks])
    reordered=_extrude_nonew_split_macroface_certificate(permuted,surfaces,signs,
        "retained proof binds actual reordered Region rows")
    @test length(reordered.boundary)==6025 && length(reordered.relations)==10
    @test Set(keys(reordered.facet_surfaces))==reordered.boundary
    @test volume.coords==before[1].coords && all(a.nodes==b.nodes for (a,b) in zip(volume.blocks,before[1].blocks))
    @test all(surfaces[tag].coords==part.coords && all(a.nodes==b.nodes for (a,b) in zip(surfaces[tag].blocks,part.blocks)) for (tag,part) in before[2])
    @test signs==before[3]
    # Every real failure keeps all source constructors untouched.
    for mutation in (
        (v,s,g)->delete!(g,13),
        (v,s,g)->(s[13]=MixedMesh(s[13].coords,[ElementBlock(3,hcat(s[13].blocks[1].nodes,s[13].blocks[1].nodes[:,1]))])),
        (v,s,g)->(s[1]=MixedMesh(s[1].coords,[ElementBlock(3,s[1].blocks[1].nodes[:,2:end])])),
        (v,s,g)->(s[26]=MixedMesh(s[26].coords,[ElementBlock(2,s[26].blocks[1].nodes[:,2:end])])),
        (v,s,g)->reverse!(@view s[1].blocks[1].nodes[:,1]),
        (v,s,g)->push!(v.blocks,ElementBlock(2,reshape(Int32[1,2,3],3,1))),
        (v,s,g)->push!(s[13].blocks,ElementBlock(1,reshape(Int32[1,2],2,1))),
        (v,s,g)->(g[13]=-g[13]),
        (v,s,g)->(v.coords[:,2].=v.coords[:,1]),
        (v,s,g)->empty!(v.blocks))
        v,s,g=deepcopy((volume,surfaces,signs));mutation(v,s,g)
        snapshot=deepcopy((v,s,g))
        @test_throws ArgumentError _extrude_nonew_split_macroface_certificate(v,s,g,"malformed native retained proof")
        @test v.coords==snapshot[1].coords && all(a.nodes==b.nodes for (a,b) in zip(v.blocks,snapshot[1].blocks))
        @test all(s[tag].coords==part.coords && all(a.nodes==b.nodes for (a,b) in zip(s[tag].blocks,part.blocks)) for (tag,part) in snapshot[2])
        @test g==snapshot[3]
    end
    @test bytes2hex(sha256(read(FIXTURE)))=="56e9823e8905d1f48edff0b1a0f646a9a463babbcf7e35a2397919c4ec25164b"
end

@testset "Retained B4 macroface receipt binds current operation and lower Source" begin
    api=Tessella.API;model_api=Tessella.Model
    native=TOML.parsefile(FIXTURE)
    api.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"retained_scope.geo");write(path,native["source"])
            api.open_geo!(path;mesh_dim=0);api.mesh.generate(2)
            model=api.CURRENT[];volume=only(keys(model.volumes))
            scope=model_api._extrude_nonew_scope(model,"current retained receipt";
                _mesh_input=model.meshing.mesh_authority)
            completed=scope.volumes[volume]
            @test model_api._extrude_nonew_retained_b4_catalog(completed.catalog)
            @test nnodes(completed.volume)==6012
            @test isempty(completed.catalog.problem_positions)
            proof=model_api._extrude_nonew_retained_b4_certificate(model,volume,
                completed.volume,completed.surfaces,completed.catalog,"current retained receipt")
            @test !isempty(proof.relations)
            @test Set(values(proof.facet_surfaces))==Set(keys(completed.surfaces))
            @test_throws ArgumentError model_api._extrude_nonew_certify_boundary(
                completed.volume,completed.surfaces,"strict default remains";oriented_internal=true)
            before=deepcopy((api.mesh.get_nodes(),api.mesh.get_elements(),
                api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag(),
                model.points,model.curves,model.loops,model.surfaces,model.volumes,model.next_tag))
            for change in (
                m->delete!(m.meshing.extrude_sources,(3,volume)),
                m->(m.meshing.extrude_sources[(3,volume)]=(2,999)),
                m->delete!(m.meshing.extrude_specs,(3,volume)),
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
