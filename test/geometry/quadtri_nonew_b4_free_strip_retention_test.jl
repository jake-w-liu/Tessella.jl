if !isdefined(@__MODULE__,:QuadTriNoNewB4FreeStripCertificates)
    include("quadtri_nonew_b4_free_strip_certificates.jl")
end
module NoNewB4FreeRetentionTests
using Tessella,Test
using Tessella.Elements:MixedMesh,ElementBlock
using Tessella.MeshTypes:nnodes
using ..QuadTriNoNewB4FreeStripCertificates
const B4=QuadTriNoNewB4FreeStripCertificates
const Model=Tessella.Model
const Free=Model._ExtrudeNoNewB4Free

function source_copy(source;kwargs...)
    fields=fieldnames(typeof(source))
    values=NamedTuple{fields}(Tuple(getfield(source,name) for name in fields))
    changed=merge(values,(;kwargs...))
    return typeof(source)(Tuple(getfield(changed,name) for name in fields)...)
end

function source_guards(model,mesh,source,spec,f,caller)
    @test Free.preflight((5,1),caller)==(12,60,0,5)
    for shape in ((4,1),(2,2),(typemax(Int),1))
        @test_throws ArgumentError Free.preflight(shape,caller)
    end
    for levels in ([0.,.5],[0.,NaN,1.],[0.,0.,1.],[1.,2.])
        refs=fill((Int32(1),Int32(1)),length(levels)-1)
        @test_throws ArgumentError Free.plan(source,levels,refs,caller)
    end
    @test_throws ArgumentError Free.plan(source,[0.,1.],NTuple{2,Int32}[],caller)
    for intervals in (0,-1,typemax(Int))
        @test_throws ArgumentError Free.SourcePhase.Phase(source,intervals)
    end
    edges=copy(source.edges);edges[1]=reverse(edges[1])
    @test_throws ArgumentError Free.SourcePhase.Phase(source_copy(source;edges),1)
    category=copy(source.category);category[1]=0x03
    @test_throws ArgumentError Free.SourcePhase.Phase(source_copy(source;category),1)
    @test_throws ArgumentError Free.SourcePhase.Phase(source_copy(source;source_cells=source.source_cells[1:end-1]),1)
    indices=copy(source.edge_indices);indices[1]=(Int32(0),indices[1][2:4]...)
    @test_throws ArgumentError Free.SourcePhase.Phase(source_copy(source;edge_indices=indices),1)
    coordinates=copy(mesh.coords);coordinates[:,2]=coordinates[:,1]
    malformed=MixedMesh(coordinates,[ElementBlock(3,copy(only(mesh.blocks).nodes))])
    @test_throws ArgumentError B4.source_complex(malformed,f)
    @test_throws ArgumentError Model._extrude_nonew_rect_grid_source(model,f.surface,malformed,spec,caller;mode=:b4_strip)
    repeated=copy(only(mesh.blocks).nodes);repeated[2,1]=repeated[1,1]
    malformed=MixedMesh(copy(mesh.coords),[ElementBlock(3,repeated)])
    @test_throws ArgumentError B4.source_complex(malformed,f)
    @test_throws ArgumentError Model._extrude_nonew_rect_grid_source(model,f.surface,malformed,spec,caller;mode=:b4_strip)
    @test_throws ArgumentError Free.FinalFactories.factory_index((0x03,0x00,0x00,0x00,0x00,0x00))
    @test_throws ArgumentError Free.CenterFan.fan((0x00,0x00,0x00,0x00,0x00,0x03))
    # Both rejections precede level/column allocation; no near-limit mesh is built.
    too_many=(Model._EXTRUDE_NONEW_MAX_NODES-12)÷17+1
    @test_throws ArgumentError Model._extrude_nonew_levels((layers=[too_many],heights=[1.]),caller;
        source_nodes=12,extra_nodes=0,extra_nodes_per_interval=5,cells_per_interval=60)
    @test_throws ArgumentError Model._extrude_nonew_levels((layers=[typemax(Int),1],heights=[.5,1.]),caller;
        source_nodes=12,extra_nodes=0,extra_nodes_per_interval=5,cells_per_interval=60)
end

@testset "FREE B4 structural retained problem with supported final mask" begin
    caller="FREE B4 structural retention control"
    f=B4.fixture("supported_prior_problem";strip_length=5,layers=:one,pins=(1,2,3,4))
    f=merge(f,(source=replace(f.source,"Layers{1}"=>"Layers{2}"),levels=[0.,.5,1.],intervals=2))
    execution=B4.execute(f.source;dim=0);model=execution.model
    out=Int.(execution.lists["sweep"]);volume_tag=out[2]
    original=Model._extrude_nonew_plan(model,volume_tag,caller)
    @test original.catalog isa Free.Catalog
    catalog=original.catalog;cols=original.sweep.cols
    phase=Free.SourcePhase.Phase(catalog.source,2);Free.SourcePhase.run!(phase)
    @test !any(phase.problems)
    cell,interval=1,1
    mask=Free.SourcePhase.mask(phase,cell,interval)
    @test mask==(0x01,0x01,0x02,0x02,0x00,0x01)
    @test Free.FinalFactories.factory_index(mask)==169
    # This explicitly injected prior problem tests downstream retention after
    # physical propagation. It is a structural control, not a naturally observed
    # primary C++ case. Separate actual Source witnesses exercise natural C1/C2.
    phase.problems[cell,interval]=true
    retained=Free.catalog(phase,catalog.levels,catalog.layer_refs)
    @test retained.template_indices[cell,interval]==0
    @test retained.problem_positions==[(Int32(cell),Int32(interval))]
    @test retained.problem_masks==[mask]
    spec=model.meshing.extrude_specs[(3,volume_tag)]
    actual=Free.emit(retained,cols,spec,caller)
    cert=B4.certify(actual.volume,f,original.source_mesh)
    # The independent certificate proves positive whole maps, exact per-macro
    # partition, the typed oriented shell, all-node usage and a complete fan.
    @test length(cert.centers)==1
    center=only(cert.centers)
    @test cert.center_parent[center]==(cell,interval)
    @test cert.center_parent[center][2]<f.intervals
    @test Tuple(get(cert.families,msh,0) for msh in 4:7)==retained.cell_counts
    @test length(cert.boundary)+length(cert.internal)==retained.face_capacity
    @test cert.total==cert.source.area*abs(B4.Q(last(cert.heights))-B4.Q(first(cert.heights)))
    surfaces=Dict{Int,Union{Tessella.Mesh,MixedMesh}}(original.source_tag=>original.source_mesh)
    top_params,_,_=Model._extrude_entity_params(model,2,out[1],caller)
    surfaces[out[1]]=Model._extrude_nonew_rect_grid_top(model,out[1],top_params,cols,retained,caller)
    for lateral in out[3:end]
        params,_,link=Model._extrude_entity_params(model,2,lateral,caller)
        surfaces[lateral]=Model._extrude_nonew_rect_grid_lateral(model,lateral,params,link[2],cols,retained,actual.edges,caller)
    end
    @test B4.Quad.surface_faces([(2,tag,mesh) for (tag,mesh) in surfaces],actual.volume)==Set(keys(cert.boundary))
    Model._extrude_nonew_certify_boundary(actual.volume,surfaces,caller;oriented_internal=true,face_capacity=retained.face_capacity)
    completed=Model._ExtrudeNoNewCompletePlan(actual.volume,surfaces,cols,retained)
    scope=Model._ExtrudeNoNewScope(model,Dict(volume_tag=>completed),surfaces)
    projected=Model.model_to_mixed(model,actual.volume,3,volume_tag;_extrude_scope=scope)
    @test nnodes(projected)==nnodes(actual.volume)
    @test projected.coords==actual.volume.coords
    @test projected.entity_data.node_entities[center]==(3,Int32(volume_tag))
    source_guards(model,original.source_mesh,catalog.source,spec,f,caller)
end
end
