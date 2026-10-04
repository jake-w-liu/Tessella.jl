using Test,Tessella
using Tessella.Elements:MixedMesh,ElementBlock,msh_spec
using Tessella.MeshTypes:Mesh,nnodes
using Tessella.Model:mesh_model_surface,mesh_model_volume,model_to_mixed
if !isdefined(@__MODULE__,:QuadTriNoNewThreeQuadStripCertificates)
    include("quadtri_nonew_three_quad_strip_certificates.jl")
end
const _TQSC=QuadTriNoNewThreeQuadStripCertificates

function _tqsc_existing_factories(cert)
    templates=Tessella.Model._extrude_nonew_templates()
    for ((parent,layer),cells) in cert.domains
        corners=Tuple(cert.grid[row,node] for row in (layer,layer+1) for node in cert.source.cells[parent])
        all(cell->all(node->node in corners,cell.nodes),cells) || continue
        local_id=Dict(node=>slot for (slot,node) in enumerate(corners))
        actual=Set((cell.msh,Tuple(sort!([local_id[node] for node in cell.nodes]))) for cell in cells)
        # Compare actual typed supports against the complete existing-corner
        # relation, independently of the selected production template index.
        @test any(templates) do template
            expected=Set((Int(cell.msh),Tuple(sort!(Int.(collect(cell.nodes[1:msh_spec(cell.msh).nnodes])))))
                for cell in template.cells[1:template.ncells])
            actual==expected
        end
    end
end

function _tqsc_product(f;entry=false)
    geometry=_TQSC.execute(f.source;dim=entry ? 3 : 0)
    out=Int.(geometry.lists["sweep"]);m=geometry.model
    source=mesh_model_surface(m,f.surface)
    volume=entry ? geo_entity_mesh(geometry,3,out[2]) : mesh_model_volume(m,out[2])
    cert=_TQSC.certify(volume,f,source)
    @test validate(volume).ok
    @test nnodes(source)==8 && length(cert.source.cells)==3
    @test sort(length.(collect(cert.source.adjacency)))==[1,1,2]
    @test sort(length.(collect(cert.source.chains)))==[2,2,4,4]
    @test length(cert.source.shared)==2
    @test length(cert.centers)==(f.laterals ? 3 : 0)
    @test length(cert.boundary)==(f.laterals ? 8f.intervals+9 : 16f.intervals+9)
    _tqsc_existing_factories(cert)
    projected=model_to_mixed(m,volume,3,out[2])
    @test validate(projected).ok
    @test _TQSC.Quad.typed_signature(projected)==_TQSC.Quad.typed_signature(volume)
    owners=projected.entity_data.node_entities
    @test [count(o->o[1]==d,owners) for d in 0:3]==[8,4f.intervals+4,4f.intervals-4,f.laterals ? 3 : 0]
    parts=[(2,s,mesh_model_surface(m,s)) for s in (f.surface,out[1],out[3:end]...)]
    @test _TQSC.Quad.surface_faces(parts,volume)==Set(keys(cert.boundary))
    @test _TQSC.Quad.surface_faces(_TQSC.execute(f.source;dim=2).mesh_parts,volume)==Set(keys(cert.boundary))
    @test _TQSC.Quad.typed_signature(mesh_model_volume(m,out[2]))==_TQSC.Quad.typed_signature(volume)
    return (;geometry,source,volume,cert,projected)
end

function _tqsc_atomic(f,diagnostic)
    g=_TQSC.execute(f.source;dim=0);m=g.model;out=Int.(g.lists["sweep"])
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

@testset "Three actual quadrangles form a conforming NoNew strip" begin
    @testset "Actual source, whole maps, typed shells and projected carriers" begin
        for plane in (:XY,:YZ,:ZX),height in (-1.,1.),laterals in (false,true),layers in (:one,:three,:graded)
            _tqsc_product(_TQSC.fixture("three_strip";plane,height,laterals,layers))
        end
        for plane in (:XY,:YZ,:XZ),laterals in (false,true),shape in (:skew,:rounded)
            f=_TQSC.fixture("graded_frame";plane,laterals,shape,height=-.75,layers=:nonbinary,
                winding=-1,pins=:opposite,curve_reverse=5,tags=:sparse,
                offset=(.2,-.3,.4),progression=4.)
            # Retain the represented interpolation of the declared group;
            # decimal .6 is a different Float64 from .2+(1-.2)/2.
            @test f.levels==[0.,.2,.2+(1.0-.2)/2,1.]
            _tqsc_product(f;entry=true)
        end
        for long_pair in ((1,3),(2,4)),laterals in (false,true)
            _tqsc_product(_TQSC.fixture("long_chain_identity";long_pair,laterals,layers=:one,
                curve_reverse=15,pins=(4,3,2,1),winding=-1))
        end
    end

    @testset "Original cell order and cyclic local frames survive emission" begin
        M=Tessella.Model
        orders=((1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1))
        for laterals in (false,true)
            f=_TQSC.fixture("source_order";laterals,layers=:three,height=-1.)
            g=_TQSC.execute(f.source;dim=0);m=g.model;out=Int.(g.lists["sweep"]);t=out[2]
            original=mesh_model_surface(m,f.surface);nodes=only(_TQSC.blocks(original)).nodes
            params=M._extrude_gate(m,3,t);spec=m.meshing.extrude_specs[(3,t)]
            levels,refs=M._extrude_nonew_levels(params,"independent real source order";
                source_nodes=8,extra_nodes=laterals ? 3 : 0,cells_per_interval=laterals ? 21 : 18)
            for order in orders,rotation in 0:3
                permuted=hcat((circshift(nodes[:,order[column]],rotation+column-1) for column in 1:3)...)
                source=MixedMesh(original.coords,[ElementBlock(3,permuted)])
                independent=_TQSC.source_complex(source,f)
                catalog=M._extrude_nonew_three_quad_strip_catalog(m,f.surface,source,spec,levels,refs,laterals,"source order")
                @test catalog.source_cells==Tuple(Tuple(cell) for cell in eachcol(permuted))
                @test independent.cells==Tuple(Tuple(Int.(cell)) for cell in eachcol(permuted))
                cols=M._extrude_volume_columns(m,t,f.surface,source,params,spec,levels,"source order")
                M._extrude_nonew_three_quad_strip_product_certify(cols,catalog,spec,"source order")
                plan=M._extrude_nonew_three_quad_strip_finish(m,t,params,f.surface,source,catalog,cols,out[1],out[3:end],"source order")
                cert=_TQSC.certify(plan.sweep.volume,f,source)
                @test validate(plan.sweep.volume).ok
                @test cert.total==1 && length(cert.centers)==(laterals ? 3 : 0)
                _tqsc_existing_factories(cert)
            end
        end
    end

    @testset "Independent certificates reject corrupted source and volume" begin
        f=_TQSC.fixture("corruption";layers=:one)
        p=_tqsc_product(f)
        nodes=copy(only(_TQSC.blocks(p.source)).nodes);nodes[:,2]=reverse(nodes[:,2])
        @test_throws r"CAD winding|equal cycles" _TQSC.source_complex(MixedMesh(p.source.coords,[ElementBlock(3,nodes)]),f)
        nodes[:,2]=nodes[:,1]
        @test_throws r"incidence|source|boundary|cycles" _TQSC.source_complex(MixedMesh(p.source.coords,[ElementBlock(3,nodes)]),f)
        coords=copy(p.volume.coords);coords[3,1]=.125
        @test_throws Exception _TQSC.certify(MixedMesh(coords,_TQSC.blocks(p.volume)),f,p.source)
    end

    @testset "Unsupported source and boundary changes reject atomically" begin
        f=_TQSC.fixture("guard";layers=:one)
        larger=replace(f.source,"}=4;"=>"}=6;")
        @test larger!=f.source
        _tqsc_atomic(merge(f,(;source=larger)),r"QuadTriNoNewVerts")
        tilted=replace(f.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{0.1,0.0,1.0}")
        @test tilted!=f.source
        _tqsc_atomic(merge(f,(;source=tilted)),r"normal|translation")
        _tqsc_atomic(merge(f,(;source=f.source*"Transfinite Surface{sweep[0]};")),r"constrained-boundary|override")
        incomplete=replace(f.source,"Layers{1}"=>"Layers{{1},{0.5}}")
        @test incomplete!=f.source
        _tqsc_atomic(merge(f,(;source=incomplete)),r"height|normalized|1.0")
    end
end
