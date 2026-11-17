using Test, Tessella, TOML
using Tessella.Elements: MixedMesh, ElementBlock, msh_dimension
using Tessella.MeshTypes: Mesh, nnodes
using Tessella.Model: mesh_model_surface, mesh_model_volume, model_to_mixed
if !isdefined(@__MODULE__,:QuadTriNoNewRectGridCertificates)
    include("quadtri_nonew_rect_grid_certificates.jl")
end
const _RGSC=QuadTriNoNewRectGridCertificates





function _rgsc_quadratic_reference()
    for msh in (11,12,13,14)
        vertices=[Tuple(_RGSC.Q(v)/2 for v in p) for p in _RGSC._P2_REFERENCE_TWICE[msh]]
        map=_RGSC.quadratic_map_certificate(msh,vertices)
        @test map.minimum==1 && map.maximum==1
        @test map.volume==(msh==11 ? 1//6 : msh==12 ? 8 : msh==13 ? 1 : 4//3)
    end
    # Positive primary corners alone do not certify an actual P2 map. These
    # independent nonprimary-node changes leave every primary corner intact.
    for (msh,node,axis,value) in ((11,5,1,2),(12,27,1,4),(13,7,1,2),(14,14,3,2))
        vertices=[Tuple(_RGSC.Q(v)/2 for v in p) for p in _RGSC._P2_REFERENCE_TWICE[msh]]
        vertices[node]=ntuple(k->k==axis ? _RGSC.Q(value) : vertices[node][k],3)
        @test_throws r"whole quadratic reference determinant" _RGSC.quadratic_map_certificate(msh,vertices)
    end
end

function _rgsc_product(f;entry=false)
    geometry=_RGSC.execute(f.source;dim=0)
    out=Int.(geometry.lists["sweep"]);m=geometry.model
    source=mesh_model_surface(m,f.surface)
    volume=mesh_model_volume(m,out[2])
    cert=_RGSC.certify(volume,f,source)
    @test validate(volume).ok
    @test nnodes(volume)==prod(f.grid.+1)*(f.intervals+1)
    parts=[(2,tag,mesh_model_surface(m,tag)) for tag in (f.surface,out[1],out[3:end]...)]
    @test _RGSC.Quad.surface_faces(parts,volume)==Set(keys(cert.boundary))
    # Certify the GEO entity products against their bound Source and surfaces.
    actual=_RGSC.execute(f.source;dim=3)
    @test Int.(actual.lists["sweep"])==out
    actual_source=geo_entity_mesh(actual,2,f.surface)
    actual_volume=geo_entity_mesh(actual,3,out[2])
    actual_cert=_RGSC.certify(actual_volume,f,actual_source;native_top=false)
    @test validate(actual_volume).ok
    @test _RGSC.Quad.surface_faces(actual.actual_mesh_parts,actual_volume)==Set(keys(actual_cert.boundary))
    @test Set(actual_source.entity_data.external_node_tags)<=Set(actual_volume.entity_data.external_node_tags)
    @test nnodes(actual_volume)==nnodes(volume)
    projected=model_to_mixed(actual.model,actual_volume,3,out[2])
    @test validate(projected;reject_duplicate_point_cells=false).ok
    @test _RGSC.Quad.typed_signature(projected)==_RGSC.Quad.typed_signature(actual_volume)
    rows=1:nnodes(actual_volume)
    @test projected.entity_data.external_node_tags[rows]==actual_volume.entity_data.external_node_tags
    @test projected.entity_data.node_entities[rows]==actual_volume.entity_data.node_entities
    entry && @test actual.mesh!==nothing
    return_geometry=entry ? actual : geometry
    return (;geometry=return_geometry,source,volume,cert,projected)
end


function _rgsc_original_frames(f)
    M=Tessella.Model;g=_RGSC.execute(f.source;dim=0);m=g.model;out=Int.(g.lists["sweep"])
    t=out[2];caller="independent rectangular source frames"
    completed=M._extrude_nonew_plan(deepcopy(m),t,caller)
    original=completed.source_mesh;nodes=only(_RGSC.blocks(original)).nodes
    params=M._extrude_gate(m,3,t);spec=m.meshing.extrude_specs[(3,t)]
    baseline=_RGSC.certify(completed.sweep.volume,f,original)
    ncells=size(nodes,2);orders=[collect(1:ncells),collect(ncells:-1:1),circshift(collect(1:ncells),1)]
    for order in orders,rotation in 0:3
        rotated=hcat((circshift(nodes[:,order[column]],mod(rotation+column-1,4)) for column in 1:ncells)...)
        source=MixedMesh(original.coords,[ElementBlock(3,rotated)])
        independent=_RGSC.source_complex(source,f)
        data=M._extrude_nonew_rect_grid_source(m,f.surface,source,spec,caller)
        @test data.source_cells==[Tuple(Int32.(cell)) for cell in eachcol(rotated)]
        @test data.source_coordinates==[Tuple(original.coords[:,node]) for node in 1:nnodes(original)]
        @test Set(findall(data.boundary_vertices))==independent.border
        catalog=M._extrude_nonew_rect_grid_plan(data,completed.catalog.levels,
            completed.catalog.layer_refs,f.laterals,caller)
        M._extrude_nonew_rect_grid_product_certify(completed.sweep.cols,catalog,spec,caller)
        plan=M._extrude_nonew_rect_grid_finish(m,t,params,f.surface,source,catalog,
            completed.sweep.cols,out[1],Tuple(out[3:end]),caller)
        cert=_RGSC.certify(plan.sweep.volume,f,source)
        @test length(cert.domains)==ncells*f.intervals
        @test Set(keys(cert.boundary))==Set(keys(baseline.boundary))
    end
end

function _rgsc_malformed_source(f)
    M=Tessella.Model;g=_RGSC.execute(f.source;dim=0);m=g.model
    source=mesh_model_surface(m,f.surface);spec=m.meshing.extrude_specs[(3,Int(g.lists["sweep"][2]))]
    nodes=only(_RGSC.blocks(source)).nodes
    stable=repr(m);coordinates=copy(source.coords);connectivity=copy(nodes)
    bad=deepcopy(source);bad.coords[:,2]=bad.coords[:,1]
    @test_throws r"physically coincident" M._extrude_nonew_rect_grid_source(m,f.surface,bad,spec,"malformed source")
    @test repr(m)==stable && source.coords==coordinates && nodes==connectivity
    bad=deepcopy(source);only(_RGSC.blocks(bad)).nodes[:,2]=nodes[:,1]
    @test_throws ArgumentError M._extrude_nonew_rect_grid_source(m,f.surface,bad,spec,"malformed source")
    @test repr(m)==stable && source.coords==coordinates && nodes==connectivity
    bad=deepcopy(source);only(_RGSC.blocks(bad)).nodes[:,1]=reverse(nodes[:,1])
    @test_throws ArgumentError M._extrude_nonew_rect_grid_source(m,f.surface,bad,spec,"malformed source")
    @test repr(m)==stable && source.coords==coordinates && nodes==connectivity
    # A simple count/Euler test would accept the untouched incidence here.
    # The actual strict convexity/embedding certificate must reject this
    # displaced interior vertex without modifying native curve parameters.
    cert=_RGSC.source_complex(source,f);bad=deepcopy(source)
    interior=minimum(cert.interior);bad.coords[f.axes[1],interior]=maximum(source.coords[f.axes[1],:])+1.
    @test_throws ArgumentError M._extrude_nonew_rect_grid_source(m,f.surface,bad,spec,"malformed source")
    @test repr(m)==stable && source.coords==coordinates && nodes==connectivity
end

function _rgsc_atomic(f,diagnostic)
    g=_RGSC.execute(f.source;dim=0);m=g.model;out=Int.(g.lists["sweep"])
    before=repr(m);params=deepcopy(m.curve_params)
    @test_throws diagnostic mesh_model_volume(m,out[2])
    @test repr(m)==before && m.curve_params==params
    context=Tessella.IO._GeoNumericContext();sentinel=Mesh(zeros(3,1))
    push!(context.mesh_parts,(0,9001,sentinel));parts=context.mesh_parts;options=copy(context.option_numbers)
    @test_throws diagnostic Tessella.GeoExec._geo_mesh_model(m,3,context)
    @test repr(m)==before && m.curve_params==params
    @test context.mesh_parts===parts && only(parts)[3]===sentinel
    @test context.option_numbers==options
end

@testset "Dynamic rectangular Quad-grid independent certificates" begin
    @testset "Exact quadratic reference spaces and interior folds" begin
        _rgsc_quadratic_reference()
    end
    @testset "Native B3/B2/B0 maps, complete shells and all public surface routes" begin
        first=_rgsc_product(_RGSC.fixture("rectangular";grid=(2,3),plane=:XY,
            height=1.,laterals=false,layers=:one))
        @test_throws r"Points\[[0-9]+\] and \[[0-9]+\] map to mesh node" model_to_mixed(
            first.geometry.model,first.volume,3,Int(first.geometry.lists["sweep"][2]))
        for (grid,plane,height,laterals,layers) in (
                ((2,3),:YZ,-1.,true,:three),
                ((3,3),:ZX,1.,false,:graded),((3,3),:XY,-1.,true,:graded))
            _rgsc_product(_RGSC.fixture("rectangular";grid,plane,height,laterals,layers))
        end
        for grid in ((2,4),(4,3),(4,4)),laterals in (false,true)
            _rgsc_product(_RGSC.fixture("dynamic_sizes";grid,laterals,layers=:graded))
        end
        for shape in (:skew,:rounded),laterals in (false,true)
            _rgsc_product(_RGSC.fixture("sampled_frames";grid=(2,3),shape,laterals,
                plane=:XZ,height=-.75,layers=:nonbinary,winding=-1,pins=:opposite,
                curve_reverse=5,tags=:sparse,offset=(.2,-.3,.4),progression=4.);entry=true)
        end
    end
    @testset "Original source IDs survive cell orders and local cyclic frames" begin
        for grid in ((2,3),(3,3)),laterals in (false,true)
            _rgsc_original_frames(_RGSC.fixture("original_identity";grid,laterals,layers=:graded,height=-1.))
        end
    end
    @testset "Malformed actual complexes and unsupported transforms fail atomically" begin
        _rgsc_malformed_source(_RGSC.fixture("malformed";grid=(3,3),layers=:one))
        for laterals in (false,true)
            f=_RGSC.fixture("collapsed_layers";grid=(2,3),laterals,layers=:four,offset=(0.,0.,2.0^50),height=.25)
            _rgsc_atomic(f,r"strictly ordered|collapse")
            f=_RGSC.fixture("oblique";grid=(2,3),laterals,layers=:one)
            _rgsc_atomic(merge(f,(;source=replace(f.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{0.25,0.0,1.0}"))),r"axis-normal|one finite nonzero")
        end
    end
end
