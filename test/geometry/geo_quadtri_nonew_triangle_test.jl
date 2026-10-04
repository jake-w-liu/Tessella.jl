using Test,Tessella
using Tessella.Elements: ElementBlock,MixedMesh,msh_dimension
using Tessella.Model: mesh_model_surface,mesh_model_volume,model_to_mixed
using Tessella.MeshTypes: nnodes

include("quadtri_nonew_triangle_certificates.jl")
const _QTNT=QuadTriNoNewTriangleCertificates
const _QTNTQ=_QTNT.Quad
const _QTNT_PINS=_QTNT.artifact_pins()

function _qtnt_api(source,action)
    api=Tessella.API
    api.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"triangle_nonew.geo");write(path,source)
            geometry=api.open_geo!(path;mesh_dim=0)
            action(api,geometry)
        end
    finally
        api.finalize()
    end
end

function _qtnt_public_volume(api,entity=1)
    tags,xyz,_=api.mesh.get_nodes();coordinates=reshape(xyz,3,:)
    positions=Dict(tag=>index for (index,tag) in enumerate(tags))
    types,_,families=api.mesh.get_elements(3,entity)
    used=sort!(unique(vcat(families...)))
    remap=Dict(tag=>Int32(index) for (index,tag) in enumerate(used))
    blocks=ElementBlock[]
    for (type,nodes) in zip(types,families)
        width=Tessella.Elements.msh_spec(type).nnodes
        push!(blocks,ElementBlock(type,reshape(Int32[remap[tag] for tag in nodes],width,:)))
    end
    return MixedMesh(coordinates[:,[positions[tag] for tag in used]],blocks)
end

function _qtnt_snapshot(api)
    model=api.CURRENT[]
    return (;model,cache=api.LAST_MESH[],class=api.LAST_MESH_CLASS[],
        payload=(api.mesh.get_nodes(),api.mesh.get_elements()),
        model_value=repr(model),params=deepcopy(model.curve_params),
        maximum_tags=(api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag()),
        counters=(api.NODE_TAG_MAX[],api.ELEMENT_TAG_MAX[]))
end

function _qtnt_unchanged(api,before)
    @test api.CURRENT[]===before.model
    @test api.LAST_MESH[]===before.cache && api.LAST_MESH_CLASS[]===before.class
    @test (api.mesh.get_nodes(),api.mesh.get_elements())==before.payload
    @test repr(before.model)==before.model_value && before.model.curve_params==before.params
    @test (api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())==before.maximum_tags
    @test (api.NODE_TAG_MAX[],api.ELEMENT_TAG_MAX[])==before.counters
end

@testset "QuadTriNoNewVerts isolated triangular source" begin
    @test length(_QTNT_PINS)==40
    @testset "actual complex, full P1 Jacobians and classified boundary" begin
        @test length(_QTNT.fixtures())==36
        for fixture in _QTNT.fixtures()
            @testset "$(fixture.name)" begin
                execution=_QTNTQ.execute(fixture.source)
                volume=geo_entity_mesh(execution,3,1)
                certificate=_QTNT.certify(volume,fixture)
                @test validate(execution.mesh).ok && validate(volume).ok
                @test nnodes(volume)==3length(fixture.levels)
                @test certificate.ncell==(fixture.laterals ? 1 : 3)*(length(fixture.levels)-1)
                @test length(certificate.boundary)==2+
                    (fixture.laterals ? 3 : 6)*(length(fixture.levels)-1)
                @test _QTNTQ.surface_faces(execution.mesh_parts,volume)==Set(keys(certificate.boundary))
                @test _QTNT.certify_cap_winding(execution.mesh_parts,volume,fixture,certificate).sigma in (-1,1)
                projected=model_to_mixed(execution.model,volume,3,1)
                @test validate(projected).ok
                @test _QTNTQ.surface_faces([(2,0,projected)],volume)==Set(keys(certificate.boundary))
                @test _QTNTQ.typed_signature(projected)==_QTNTQ.typed_signature(volume)
                owners=projected.entity_data.node_entities
                @test all(owner->owner[1] in (0,1),owners)
                @test count(owner->owner[1]==0,owners)==6
                @test count(owner->owner[1]==1,owners)==3(length(fixture.levels)-2)
                again=_QTNTQ.execute(fixture.source)
                @test _QTNT.crc(execution.mesh)==_QTNT.crc(again.mesh)
                @test _QTNT.global_record(fixture,execution)==_QTNT_PINS[fixture.name]
            end
        end
    end

    @testset "standalone, surface-only GEO and volume-only API agree" begin
        for motion in (:translation,:rotation,:twist),laterals in (false,true),winding in (-1,1)
            fixture=_QTNT.fixture("triangle_entrypoints","Layers{{2,1},{0.2,1.0}}",
                [0.,.1,.2,1.],laterals;motion,winding)
            geometry=_QTNTQ.execute(fixture.source;dim=0)
            surfaces=[1,Int(geometry.lists["sweep"][1]),Int.(geometry.lists["sweep"][3:end])...]
            @test length(surfaces)==5 && allunique(surfaces)
            before=[(2,tag,mesh_model_surface(geometry.model,tag)) for tag in surfaces]
            volume=mesh_model_volume(geometry.model,1)
            certificate=_QTNT.certify(volume,fixture)
            boundary=Set(keys(certificate.boundary))
            @test _QTNTQ.surface_faces(before,volume)==boundary
            after=[(2,tag,mesh_model_surface(geometry.model,tag)) for tag in surfaces]
            @test _QTNTQ.surface_faces(after,volume)==boundary
            surface_only=_QTNTQ.execute(fixture.source;dim=2)
            @test _QTNTQ.surface_faces(surface_only.mesh_parts,volume)==boundary
            _qtnt_api(fixture.source,(api,_)->begin
                surface_product=api.mesh.generate(2)
                @test _QTNTQ.surface_faces([(2,0,surface_product)],volume)==boundary
                generated=api.mesh.generate(3)
                @test _QTNTQ.typed_signature(generated)==_QTNTQ.typed_signature(volume)
                @test _QTNT.certify(generated,fixture).family_counts==certificate.family_counts
                @test _QTNTQ.typed_signature(_qtnt_public_volume(api))==_QTNTQ.typed_signature(volume)
                @test isempty(api.mesh.get_elements(2)[1])
                @test length(api.mesh.get_nodes(3,1,true)[1])==nnodes(volume)
                @test isempty(api.mesh.get_nodes(3,1,false)[1])
                @test sum(length,api.mesh.get_elements(3,1)[2])==certificate.ncell
            end)
        end
    end

    @testset "equivalent graded groups and affine full P2 supports" begin
        for laterals in (false,true)
            a=_QTNT.fixture("triangle_groups_a","Layers{{2,1},{0.2,1.0}}",[0.,.1,.2,1.],laterals)
            b=_QTNT.fixture("triangle_groups_b","Layers{{1,1,1},{0.1,0.2,1.0}}",[0.,.1,.2,1.],laterals)
            @test _QTNT.crc(_QTNTQ.execute(a.source).mesh)==_QTNT.crc(_QTNTQ.execute(b.source).mesh)
        end
        for laterals in (false,true),winding in (-1,1)
            fixture=_QTNT.fixture("triangle_quadratic","Layers{3}",[0.,1/3,2/3,1.],laterals;winding)
            _qtnt_api(fixture.source,(api,geometry)->begin
                api.mesh.generate(3);linear=_qtnt_public_volume(api)
                api.mesh.set_order(2);quadratic=_qtnt_public_volume(api)
                @test _QTNTQ.certify_quadratic(quadratic,linear).ncell==(laterals ? 3 : 9)
                @test Set(b.msh for b in quadratic.blocks)==Set([laterals ? 13 : 11])
                artifact_name="triangle_p2_W$(winding)_R$(laterals)"
                @test _QTNT.api_record(artifact_name,api)==_QTNT_PINS[artifact_name]
                for surface in Int.(geometry.lists["sweep"][3:end])
                    nodes,_,parameters=api.mesh.get_nodes(2,surface,false,true)
                    @test length(nodes)==5 && length(parameters)==2length(nodes)
                    for node in nodes
                        p,uv,dim,owner=api.mesh.get_node(node)
                        @test (dim,owner)==(2,surface) && length(uv)==2
                        @test maximum(abs.(api.model.get_value(2,surface,uv).-p))<=2e-11
                    end
                    nodes,coordinates,parameters=api.mesh.get_nodes(2,surface,true,true)
                    @test length(nodes)==21 && length(parameters)==2length(nodes)
                    @test maximum(abs.(api.model.get_value(2,surface,parameters).-coordinates))<=2e-11
                end
                api.mesh.set_order(1)
                @test _QTNTQ.typed_signature(_qtnt_public_volume(api))==_QTNTQ.typed_signature(linear)
            end)
        end
    end

    @testset "projection accepts actual permutations and rejects changed cells" begin
        for laterals in (false,true)
            fixture=_QTNT.fixture("triangle_project","Layers{3}",[0.,1/3,2/3,1.],laterals;motion=:twist)
            execution=_QTNTQ.execute(fixture.source);volume=geo_entity_mesh(execution,3,1)
            blocks=_QTNTQ.blocks(volume)
            permutation=collect(nnodes(volume):-1:1);remap=invperm(permutation)
            reordered=MixedMesh(volume.coords[:,permutation],
                [ElementBlock(b.msh,Int32.(remap[b.nodes[:,end:-1:1]])) for b in blocks])
            projected=model_to_mixed(execution.model,reordered,3,1)
            @test _QTNTQ.typed_signature(projected)==_QTNTQ.typed_signature(volume)
            @test _QTNTQ.surface_faces([(2,0,projected)],volume)==
                Set(keys(_QTNT.certify(volume,fixture).boundary))
            changed=copy(reordered.coords);changed[1,1]+=1e-5
            parameters=deepcopy(execution.model.curve_params);attributes=repr(execution.model.meshing)
            @test_throws r"QuadTriNoNewVerts" model_to_mixed(execution.model,
                MixedMesh(changed,reordered.blocks),3,1)
            @test execution.model.curve_params==parameters && repr(execution.model.meshing)==attributes
        end
    end

    @testset "physical names and boundary classification survive projection" begin
        for laterals in (false,true)
            fixture=_QTNT.fixture("triangle_physical","Layers{3}",[0.,1/3,2/3,1.],laterals;motion=:twist)
            source=fixture.source*"""
            Physical Surface("source",101)={1};
            Physical Surface("top",102)={sweep[0]};
            Physical Volume("body",103)={sweep[1]};
            """
            execution=_QTNTQ.execute(source);volume=geo_entity_mesh(execution,3,1)
            projected=model_to_mixed(execution.model,volume,3,1)
            @test projected.physical_names==Dict((2,101)=>"source",(2,102)=>"top",(3,103)=>"body")
            @test all(all(==(103),block.tags) for block in projected.blocks if msh_dimension(block.msh)==3)
            @test _QTNTQ.surface_faces([(2,0,projected)],volume)==
                Set(keys(_QTNT.certify(volume,fixture).boundary))
        end
    end

    @testset "unsupported sources and boundaries fail before publication" begin
        fixture=_QTNT.fixture("triangle_guard","Layers{3}",[0.,1/3,2/3,1.],false)
        _qtnt_api(fixture.source,(api,_)->begin
            api.mesh.generate(3)
            for curve in 1:3;api.mesh.set_transfinite_curve(curve,3);end
            before=_qtnt_snapshot(api)
            @test_throws r"multiple source cells" api.mesh.generate(3)
            _qtnt_unchanged(api,before)
        end)
        invalid=(
            (replace(fixture.source,"Transfinite Curve{:}=2"=>"Transfinite Curve{:}=3"),r"multiple source cells"),
            (replace(fixture.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{{0,0,1},{0,0,0},Pi/6}"),r"collapsed"),
            (fixture.source*"Transfinite Surface{sweep[0]};\n",r"transfinite cap/lateral overrides"),
            (fixture.source*"Extrude{0,0,-1}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;}\n",r"shared surfaces"),
            (fixture.source*"Extrude{0,0,1}{Surface{sweep[0]};Layers{3};Recombine;QuadTriNoNewVerts;}\n",r"copied or chained sources"))
        for (source,diagnostic) in invalid
            geometry=_QTNTQ.execute(source;dim=0)
            parameters=deepcopy(geometry.model.curve_params);attributes=repr(geometry.model.meshing)
            volumes=sort!(collect(keys(geometry.model.volumes)))
            target=occursin("chained",diagnostic.pattern) ? last(volumes) : first(volumes)
            # The all-region request encounters the earlier shared cap before
            # reaching the later region's copied-source dependency guard.
            global_diagnostic=target==last(volumes) && length(volumes)>1 &&
                occursin("chained",diagnostic.pattern) ? r"shared surfaces" : diagnostic
            @test_throws diagnostic mesh_model_volume(geometry.model,target)
            @test geometry.model.curve_params==parameters && repr(geometry.model.meshing)==attributes
            @test_throws global_diagnostic _QTNTQ.execute(source)
            context=Tessella.IO._GeoNumericContext()
            old_part=Tessella.MeshTypes.Mesh(zeros(3,1));push!(context.mesh_parts,(0,777,old_part))
            parts=context.mesh_parts;options=deepcopy(context.option_numbers);model_value=repr(geometry.model)
            @test_throws global_diagnostic Tessella.GeoExec._geo_mesh_model(geometry.model,3,context)
            @test repr(geometry.model)==model_value && context.option_numbers==options
            @test context.mesh_parts===parts && only(parts)[3]===old_part
        end
        mixed=MixedMesh(Float64[0 1 0 1;0 0 1 1;0 0 0 0],
            [ElementBlock(2,reshape(Int32[1,2,3],3,1)),ElementBlock(3,reshape(Int32[1,2,4,3],4,1))])
        @test_throws r"multiple source cells" Tessella.Model._extrude_nonew_source_cell(mixed,"triangle-test")
        _qtnt_api(fixture.source,(api,_)->begin
            api.mesh.generate(3)
            model=api.CURRENT[];old=model.meshing.extrude[(3,1)]
            model.meshing.extrude[(3,1)]=merge(old,(;layers=[1],heights=[.5]))
            before=_qtnt_snapshot(api)
            @test_throws r"normalized final layer height of 1.0" api.mesh.generate(3)
            _qtnt_unchanged(api,before)
        end)
    end

    @testset "independent triangular identities remain distinct when coincident" begin
        source=_QTNT.paired_source()
        execution=_QTNTQ.execute(source)
        @test nnodes(execution.mesh)==24
        _qtnt_api(source,(api,geometry)->begin
            for (tag,p) in collect(api.CURRENT[].points)
                p[1]>=3 || continue
                api.model.set_coordinates(tag,p[1]-3,p[2],p[3])
            end
            generated=api.mesh.generate(3)
            left=Int(geometry.lists["left"][2]);right=Int(geometry.lists["right"][2])
            a=api.mesh.get_nodes(3,left,true)[1];b=api.mesh.get_nodes(3,right,true)[1]
            @test length(a)==12 && length(b)==12 && isempty(intersect(a,b))
            @test nnodes(generated)==24
            @test length(unique(Tuple(generated.coords[:,i]) for i in 1:nnodes(generated)))==12
            for (entity,laterals) in ((left,false),(right,true))
                fixture=_QTNT.fixture("triangle_paired","Layers{3}",[0.,1/3,2/3,1.],laterals)
                @test _QTNT.certify(_qtnt_public_volume(api,entity),fixture).ncell==(laterals ? 3 : 9)
            end
            @test validate(generated).ok
        end)
    end
end
