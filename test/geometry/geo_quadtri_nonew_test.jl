using Test, Tessella
using Tessella.Elements: ElementBlock, MixedMesh, msh_spec
using Tessella.Model: mesh_model_surface, mesh_model_volume, model_to_mixed
using Tessella.MeshTypes: nnodes

include("quadtri_nonew_certificates.jl")
const _QTNN=QuadTriNoNewCertificates
const _QTNN_ALLOCATION_RESULT=Ref((0,0))

function _qtnn_api(source,action)
    api=Tessella.API
    api.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"nonew.geo")
            write(path,source)
            execution=api.open_geo!(path;mesh_dim=0)
            action(api,execution)
        end
    finally
        api.finalize()
    end
end

@testset "NoNew independent certificates reject invalid complexes" begin
    coordinates=Float64[0 1 1 0 0 1 1 0;0 0 1 1 0 0 1 1;0 0 0 0 1 1 1 1]
    nodes=reshape(Int32.(1:8),8,1)
    cube=MixedMesh(coordinates,[ElementBlock(5,nodes)])
    certificate=_QTNN.certify_complex(cube)
    @test certificate.total≈1.0 atol=1e-14
    @test length(certificate.boundary)==6
    repeated=copy(nodes);repeated[2]=repeated[1]
    @test_throws r"repeated cell corner" _QTNN.certify_complex(
        MixedMesh(coordinates,[ElementBlock(5,repeated)]))
    reversed=nodes[[2,1,4,3,6,5,8,7],:]
    @test_throws r"nonpositive exact face fan" _QTNN.certify_complex(
        MixedMesh(coordinates,[ElementBlock(5,reversed)]))
    @test_throws r"interior face orientations" _QTNN.certify_complex(
        MixedMesh(coordinates,[ElementBlock(5,hcat(nodes,nodes))]))
end

@testset "QuadTriNoNewVerts single source quadrangle" begin
    @testset "sweep cells, source/cap policies, and classified boundary" begin
        for fixture in _QTNN.fixtures()
            @testset "$(fixture.name)" begin
                execution=_QTNN.execute(fixture.source)
                volume=geo_entity_mesh(execution,3,1;geometry=true)
                certificate=_QTNN.certify_fixture(volume,fixture)
                @test validate(volume).ok
                @test validate(execution.mesh).ok
                @test nnodes(volume)==4length(fixture.levels)+(fixture.laterals ? 1 : 0)
                @test Set(keys(certificate.boundary))==
                    _QTNN.surface_faces(execution.mesh_parts,volume)
                @test length(certificate.boundary)==3+
                    (fixture.laterals ? 4 : 8)*(length(fixture.levels)-1)
                projected=model_to_mixed(execution.model,volume,3,1)
                @test validate(projected).ok
                @test _QTNN.counts(projected)==certificate.family_counts
                @test all(all(==(1),owners) for (block,owners) in
                    zip(projected.blocks,projected.elementary_entities) if
                    block isa ElementBlock && Tessella.Elements.msh_dimension(block.msh)==3)
                again=_QTNN.execute(fixture.source)
                @test mixed_crc(execution.mesh)==mixed_crc(again.mesh)
            end
        end
    end

    @testset "standalone surfaces before and after volume planning" begin
        for laterals in (false,true)
            fixture=_QTNN.fixture("standalone","Layers{{1,2},{0.5,1.0}}",
                [0.,.5,.75,1.],laterals)
            geometry=_QTNN.execute(fixture.source;dim=0)
            top=Int(geometry.lists["sweep"][1])
            sides=Int.(geometry.lists["sweep"][3:end])
            tags=[1,top,sides...]
            before=[(2,tag,mesh_model_surface(geometry.model,tag)) for tag in tags]
            direct=mesh_model_volume(geometry.model,1)
            certificate=_QTNN.certify_fixture(direct,fixture)
            @test _QTNN.surface_faces(before,direct)==Set(keys(certificate.boundary))
            after=[(2,tag,mesh_model_surface(geometry.model,tag)) for tag in tags]
            @test _QTNN.surface_faces(before,direct)==_QTNN.surface_faces(after,direct)
            execution=_QTNN.execute(fixture.source)
            @test _QTNN.typed_signature(direct)==
                _QTNN.typed_signature(geo_entity_mesh(execution,3,1;geometry=true))
            surface_only=_QTNN.execute(fixture.source;dim=2)
            @test _QTNN.surface_faces(surface_only.mesh_parts,direct)==
                Set(keys(certificate.boundary))
            _qtnn_api(fixture.source,(api,_)->begin
                generated=api.mesh.generate(3)
                @test _QTNN.typed_signature(generated)==_QTNN.typed_signature(direct)
                @test _QTNN.certify_fixture(generated,fixture).family_counts==certificate.family_counts
                @test validate(generated).ok
                # The bounding Surfaces keep their triangles and quadrangles
                # after generate(3), as in Gmsh.
                @test api.mesh.get_elements(2)[1]==Int32[2,3]
                @test length(api.mesh.get_nodes(3,1,true)[1])==nnodes(generated)
                @test length(api.mesh.get_nodes(3,1,false)[1])==(laterals ? 1 : 0)
                @test sum(length,api.mesh.get_elements(3,1)[2])==certificate.ncell
            end)
        end
    end

    @testset "curved surface generation and projection preserve the planned boundary" begin
        for motion in (:rotation,:twist),laterals in (false,true)
            fixture=_QTNN.fixture("surface_projection","Layers{3}",
                [0.,1/3,2/3,1.],laterals;motion)
            source=fixture.source*"""
            Physical Surface("source",101)={1};
            Physical Surface("top",102)={sweep[0]};
            Physical Volume("body",103)={sweep[1]};
            """
            execution=_QTNN.execute(source)
            volume=geo_entity_mesh(execution,3,1;geometry=true)
            boundary=Set(keys(_QTNN.certify_fixture(volume,fixture).boundary))
            projected=model_to_mixed(execution.model,volume,3,1)
            @test projected.physical_names==Dict((2,101)=>"source",(2,102)=>"top",(3,103)=>"body")
            @test _QTNN.surface_faces([(2,0,projected)],volume)==boundary
            permutation=collect(nnodes(volume):-1:1);remap=invperm(permutation)
            reordered_blocks=[ElementBlock(b.msh,reshape(Int32[remap[n] for n in b.nodes],
                size(b.nodes))) for b in volume.blocks]
            renumbered=MixedMesh(volume.coords[:,permutation],reordered_blocks)
            reordered=model_to_mixed(execution.model,renumbered,3,1)
            @test reordered.physical_names==projected.physical_names
            @test _QTNN.surface_faces([(2,0,reordered)],renumbered)==
                Set(keys(_QTNN.certify_fixture(renumbered,fixture).boundary))
            _qtnn_api(source,(api,_)->begin
                surface=api.mesh.generate(2)
                @test validate(surface).ok
                @test _QTNN.surface_faces([(2,0,surface)],volume)==boundary
                @test sum(length,api.mesh.get_elements(2)[2])==length(boundary)
                @test length(api.mesh.get_nodes(2,1,true)[1])==4
                top=Int(execution.lists["sweep"][1])
                @test length(api.mesh.get_nodes(2,top,true)[1])==4
                @test length(api.mesh.get_elements(2,top)[2][1])==2
                @test all(Tessella.Elements.msh_dimension(b.msh)==2 for b in surface.blocks)
            end)
            geometry=_QTNN.execute(source;dim=0)
            params=deepcopy(geometry.model.curve_params)
            coordinates=copy(volume.coords);coordinates[1,1]+=.01
            displaced=MixedMesh(coordinates,volume.blocks)
            @test_throws r"coordinates do not match" model_to_mixed(geometry.model,displaced,3,1)
            @test geometry.model.curve_params==params
            reversed_blocks=copy(volume.blocks)
            first=firstindex(reversed_blocks)
            reversed_nodes=copy(reversed_blocks[first].nodes)
            reversed_nodes[[1,2],1]=reversed_nodes[[2,1],1]
            reversed_blocks[first]=ElementBlock(reversed_blocks[first].msh,reversed_nodes)
            reversed=MixedMesh(volume.coords,reversed_blocks)
            @test_throws r"cell topology or winding" model_to_mixed(geometry.model,reversed,3,1)
            @test geometry.model.curve_params==params
        end
    end

    @testset "equivalent positive layer groups preserve the mesh" begin
        for laterals in (false,true)
            first=_QTNN.fixture("first","Layers{{1,2},{0.5,1.0}}",[0.,.5,.75,1.],laterals)
            second=_QTNN.fixture("second","Layers{{1,1,1},{0.5,0.75,1.0}}",[0.,.5,.75,1.],laterals)
            a=_QTNN.execute(first.source);b=_QTNN.execute(second.source)
            @test mixed_crc(a.mesh)==mixed_crc(b.mesh)
        end
    end

    @testset "retained CAD points have stable tag order" begin
        for motion in (:rotation,:twist),laterals in (false,true)
            fixture=_QTNN.fixture("point_order","Layers{{1,2},{0.5,1.0}}",
                [0.,.5,.75,1.],laterals;motion)
            geometry=_QTNN.execute(fixture.source;dim=0)
            expected=mixed_crc(_QTNN.execute(fixture.source).mesh)
            for reverse_order in (false,true),extra_capacity in (0,257)
                model=deepcopy(geometry.model)
                entries=sort!(collect(model.points);by=first,rev=reverse_order)
                model.points=Dict{Int,NTuple{3,Float64}}()
                sizehint!(model.points,length(entries)+extra_capacity)
                for (tag,coordinate) in entries;model.points[tag]=coordinate;end
                context=Tessella.IO._GeoNumericContext()
                actual,_=Tessella.GeoExec._geo_mesh_model(model,3,context)
                @test mixed_crc(actual)==expected
                point_tags=[tag for (dim,tag,_) in context.mesh_parts if dim==0]
                @test point_tags==sort(point_tags)
            end
        end
    end

    @testset "cap precedence follows the original source-node indices" begin
        for corners in ("1,2,3,4","2,3,4,1","3,4,1,2","4,1,2,3")
            fixture=_QTNN.fixture("source_precedence","Layers{1}",[0.,1.],false)
            source=replace(fixture.source,"Transfinite Surface{1};"=>
                           "Transfinite Surface{1}={$corners};")
            execution=_QTNN.execute(source)
            original=mesh_model_surface(execution.model,1)
            cell=only(eachcol(only(original.blocks).nodes))
            first=findmin(cell)[2];opposite=mod1(first+2,4)
            expected=Set(_QTNN.transformed(fixture,_QTNN.point(original,cell[i]),1.)
                         for i in (first,opposite))
            top=geo_entity_mesh(execution,2,Int(execution.lists["sweep"][1]);geometry=true)
            shared=intersect(collect(top.tris[:,1]),collect(top.tris[:,2]))
            @test Set(_QTNN.point(top,node) for node in shared)==expected
            @test _QTNN.certify_fixture(geo_entity_mesh(execution,3,1;geometry=true),fixture).ncell==5
        end
    end

    @testset "full quadratic affine cache and order roundtrip" begin
        for layers in (1,3),laterals in (false,true)
            levels=collect(range(0.,1.;length=layers+1))
            fixture=_QTNN.fixture("quadratic","Layers{$layers}",levels,laterals)
            _qtnn_api(fixture.source,(api,_)->begin
                first=api.mesh.generate(3)
                original=_QTNN.typed_signature(first)
                api.mesh.set_order(2)
                quadratic=api.mesh.get()
                @test validate(quadratic).ok
                @test _QTNN.certify_quadratic(quadratic,first).ncell==
                    _QTNN.certify_fixture(first,fixture).ncell
                @test all(block->msh_spec(block.msh).order==2,quadratic.blocks)
                expected=Dict(4=>11,5=>12,6=>13,7=>14)
                @test Dict(b.msh=>size(b.nodes,2) for b in quadratic.blocks)==
                    Dict(expected[type]=>count for (type,count) in _QTNN.counts(first))
                tags,coordinates,_=api.mesh.get_nodes()
                @test length(tags)==length(unique(tags))
                @test all(isfinite,coordinates)
                api.mesh.set_order(1)
                @test _QTNN.typed_signature(api.mesh.get())==original
                @test _QTNN.certify_fixture(api.mesh.get(),fixture).family_counts==_QTNN.counts(first)
            end)
        end
    end

    @testset "seven source shapes with and without AutoCoherence, as in Gmsh" begin
        # Gmsh 4.15.2 on these sources (n x m cells, recombined or not; node
        # counts without / with RecombLaterals):
        # - Geometry.AutoCoherence=1 meshes all of them, with the `merged`
        #   node counts below.
        # - AutoCoherence=0 keeps 44 Points on the swept square. Without
        #   RecombLaterals Mesh 3 fails on every source with
        #   "FindDiagonalEdgeIndices() encountered unexpected surface
        #   configuration"; with RecombLaterals it returns an unmerged mesh
        #   with the `unmerged` node count.
        # Tessella agrees on every coherent case. Without AutoCoherence and
        # without RecombLaterals it stops on every source (with Gmsh's
        # message for the two it can bind); with RecombLaterals it returns
        # Gmsh's unmerged mesh for the 5x1 and 2x2 sources and still refuses
        # the other five, marked broken: Gmsh's mesh for those comes from
        # its general QuadToTri planner classifying unmerged copies, which
        # Tessella's catalogs do not reproduce.
        shapes=((1,1,true,(8,9),44,true),(2,1,true,(12,14),50,true),
            (3,1,true,(16,19),56,true),(4,1,true,(20,24),62,true),
            (5,1,true,(24,29),68,false),(2,2,true,(18,18),58,false),
            (1,1,false,(8,8),44,true))
        for (n,m,recombined,merged,unmerged,refused) in shapes,laterals in (false,true)
            body="""
            Geometry.OldRuledSurface=0;Mesh.ElementOrder=1;
            Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={1,1,0,1};Point(4)={0,1,0,1};
            Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
            Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
            Transfinite Curve{1,3}=$(n+1);Transfinite Curve{2,4}=$(m+1);
            Transfinite Surface{1} Left;$(recombined ? "Recombine Surface{1};" : "")
            sweep[]=Extrude{0.0,0.0,-1.0}{Surface{1};Layers{1};Recombine;
              QuadTriNoNewVerts$(laterals ? " RecombLaterals" : "");};
            """
            @test nnodes(_QTNN.execute("Geometry.AutoCoherence=1;"*body).mesh)==
                merged[laterals ? 2 : 1]
            separate="Geometry.AutoCoherence=0;"*body
            if !laterals && refused
                @test_throws ArgumentError _QTNN.execute(separate)
            elseif !laterals
                @test_throws r"FindDiagonalEdgeIndices\(\) encountered unexpected surface configuration" _QTNN.execute(separate)
            elseif refused
                @test_broken try nnodes(_QTNN.execute(separate).mesh)==unmerged catch;false end
            else
                @test nnodes(_QTNN.execute(separate).mesh)==unmerged
            end
        end
    end

    @testset "guard failures leave the API cache and model state intact" begin
        fixture=_QTNN.fixture("guard","Layers{3}",[0.,1/3,2/3,1.],false)
        _qtnn_api(fixture.source,(api,_)->begin
            api.mesh.generate(3)
            for curve in 1:4;api.mesh.set_transfinite_curve(curve,4);end
            model=api.CURRENT[];cache=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
            payload=(api.mesh.get_nodes(),api.mesh.get_elements())
            params=deepcopy(model.curve_params)
            attrs=deepcopy(model.meshing.extrude)
            counters=(api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())
            @test_throws r"QuadTriNoNewVerts" api.mesh.generate(3)
            @test api.CURRENT[]===model && api.LAST_MESH[]===cache && api.LAST_MESH_CLASS[]===class
            @test (api.mesh.get_nodes(),api.mesh.get_elements())==payload
            @test model.curve_params==params && model.meshing.extrude==attrs
            @test (api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())==counters
        end)
        # Valid all-boundary strips are supported; unmatched opposite
        # source chains remain malformed.
        sources=(
            replace(fixture.source,"Transfinite Curve{:}=2"=>
                "Transfinite Curve{1,3}=2; Transfinite Curve{2}=6; Transfinite Curve{4}=7"),
            replace(fixture.source,"Recombine Surface{1};"=>"",
                    "Transfinite Curve{:}=2"=>"Transfinite Curve{:}=3"),
            replace(fixture.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{{0,1,0},{0,0,0},Pi/3}"))
        # Gmsh 4.15.2 meshes all three (124, 36 and 10 nodes): a transfinite
        # Surface with unequal opposite sides falls back to the unstructured
        # algorithm, an unrecombined source sweeps to tetrahedra, and a
        # rotation about an edge of the source collapses that edge. Tessella
        # meshes the second and still refuses the first and the third.
        mismatched,triangles,hinged=sources
        for source in (mismatched,hinged)
            geometry=_QTNN.execute(source;dim=0)
            @test_throws r"QuadTriNoNewVerts|non-matching node counts" mesh_model_volume(geometry.model,1)
            @test_throws r"QuadTriNoNewVerts|non-matching node counts" _QTNN.execute(source)
            @test_broken try mesh_model_volume(geometry.model,1);true catch;false end
        end
        geometry=_QTNN.execute(triangles;dim=0)
        @test nnodes(mesh_model_volume(geometry.model,1))==36
        @test _QTNN.execute(triangles).mesh!==nothing
        for (layers,heights,diagnostic) in (([1,1],[.5,.5],r"strictly increasing"),
                ([1],[Inf],r"strictly increasing"),
                ([1,2],[prevfloat(1.),1.],r"levels collapse in Float64"),
                ([2_500_000],[1.],r"node limit"),
                ([typemax(Int),1],[.5,1.],r"overflows Int"))
            geometry=_QTNN.execute(fixture.source;dim=0)
            old=geometry.model.meshing.extrude[(3,1)]
            geometry.model.meshing.extrude[(3,1)]=merge(old,(;layers,heights))
            params=deepcopy(geometry.model.curve_params)
            @test_throws diagnostic mesh_model_volume(geometry.model,1)
            @test geometry.model.curve_params==params
        end
        for height in (.5,2.),laterals in (false,true)
            invalid=_QTNN.fixture("terminal_height","Layers{{1},{$height}}",
                [0.,height],laterals)
            geometry=_QTNN.execute(invalid.source;dim=0)
            params=deepcopy(geometry.model.curve_params)
            @test_throws r"normalized final layer height of 1.0" mesh_model_volume(
                geometry.model,1)
            @test geometry.model.curve_params==params
            @test_throws r"normalized final layer height of 1.0" _QTNN.execute(invalid.source)
            _qtnn_api(fixture.source,(api,_)->begin
                api.mesh.generate(3)
                model=api.CURRENT[];cache=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
                old=model.meshing.extrude[(3,1)]
                model.meshing.extrude[(3,1)]=merge(old,(;layers=[1],heights=[height]))
                attributes=repr(model.meshing)
                params=deepcopy(model.curve_params)
                payload=(api.mesh.get_nodes(),api.mesh.get_elements())
                counters=(api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())
                @test_throws r"normalized final layer height of 1.0" api.mesh.generate(3)
                @test api.CURRENT[]===model && api.LAST_MESH[]===cache &&
                    api.LAST_MESH_CLASS[]===class
                @test (api.mesh.get_nodes(),api.mesh.get_elements())==payload
                @test repr(model.meshing)==attributes && model.curve_params==params
                @test (api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())==counters
            end)
        end
    end

    @testset "curved sweeps are raised to second order with Gmsh's node counts" begin
        # Gmsh 4.15.2, generate(3) then setOrder(2) on the same four sources:
        # 18 -> 65, 19 -> 73, 32 -> 79 and 33 -> 87 nodes.
        gmsh_nodes=Dict((:rotation,false)=>(18,65),(:rotation,true)=>(19,73),
            (:twist,false)=>(32,79),(:twist,true)=>(33,87))
        for motion in (:rotation,:twist),laterals in (false,true)
            fixture=_QTNN.fixture("curved_order_guard","Layers{3}",
                [0.,1/3,2/3,1.],laterals;motion)
            _qtnn_api(fixture.source,(api,_)->begin
                api.mesh.generate(3)
                model=api.CURRENT[];cache=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
                payload=(api.mesh.get_nodes(),api.mesh.get_elements())
                params=deepcopy(model.curve_params)
                # Affine transform records carry mutable step arrays, so the
                # attributes' generic equality does not compare deep copies.
                attributes=repr(model.meshing)
                @test length(payload[1][1])==gmsh_nodes[(motion,laterals)][1]
                api.mesh.set_order(2)
                @test length(api.mesh.get_nodes()[1])==gmsh_nodes[(motion,laterals)][2]
                @test model.curve_params==params && repr(model.meshing)==attributes
            end)
        end
    end

    @testset "full revolutions require global intersection planning" begin
        fixture=_QTNN.fixture("revolution_guard","Layers{100}",
            collect(range(0.,1.;length=101)),true;motion=:rotation)
        source="Geometry.AutoCoherence=1;\n"*replace(fixture.source,"Pi/6"=>"25*Pi/6")
        geometry=_QTNN.execute(source;dim=0)
        params=deepcopy(geometry.model.curve_params)
        @test_throws r"full or multiple revolutions" mesh_model_volume(geometry.model,1)
        @test geometry.model.curve_params==params
        @test_throws r"full or multiple revolutions" _QTNN.execute(source)
        valid=_QTNN.fixture("valid_rotation","Layers{3}",[0.,1/3,2/3,1.],true;
            motion=:rotation)
        _qtnn_api(valid.source,(api,_)->begin
            api.mesh.generate(3)
            model=api.CURRENT[];cache=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
            spec=model.meshing.extrude_specs[(3,1)]
            model.meshing.extrude_specs[(3,1)]=merge(spec,(;angle=25pi/6))
            old=model.meshing.extrude[(3,1)]
            model.meshing.extrude[(3,1)]=merge(old,(;layers=[100],heights=[1.]))
            attributes=repr(model.meshing);params=deepcopy(model.curve_params)
            payload=(api.mesh.get_nodes(),api.mesh.get_elements())
            counters=(api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())
            @test_throws r"full or multiple revolutions" api.mesh.generate(3)
            @test api.CURRENT[]===model && api.LAST_MESH[]===cache &&
                api.LAST_MESH_CLASS[]===class
            @test (api.mesh.get_nodes(),api.mesh.get_elements())==payload
            @test model.curve_params==params && repr(model.meshing)==attributes
            @test (api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())==counters
        end)
    end

    @testset "actual helical columns require global disjointness" begin
        for pitch in (-.1,.1),laterals in (false,true)
            fixture=_QTNN.helical_fixture("thin_helix",pitch,laterals)
            witness=(.4935836460303578,.5+pitch/104,-.09789464416503871)
            for interval in (1,49)
                vertices=Tuple(_QTNN.transformed(fixture,p,fixture.levels[level])
                    for level in (interval,interval+1) for p in _QTNN.CORNERS)
                parameters=_QTNN.trilinear_parameters(vertices,witness)
                @test maximum(abs,parameters)<1-1e-9
            end
            geometry=_QTNN.execute(fixture.source;dim=0)
            params=deepcopy(geometry.model.curve_params)
            @test_throws r"global P1 cell-hull disjointness is not certified" mesh_model_volume(
                geometry.model,1)
            @test geometry.model.curve_params==params
            @test_throws r"global P1 cell-hull disjointness is not certified" _QTNN.execute(fixture.source)
        end
        for pitch in (-3.,3.),laterals in (false,true)
            fixture=_QTNN.helical_fixture("separated_helix",pitch,laterals)
            execution=_QTNN.execute(fixture.source)
            volume=geo_entity_mesh(execution,3,1;geometry=true)
            @test _QTNN.certify_fixture(volume,fixture).ncell>0
            @test nnodes(volume)==212+(laterals ? 1 : 0)
            @test validate(execution.mesh).ok
            _qtnn_api(fixture.source,(api,_)->begin
                generated=api.mesh.generate(3)
                @test _QTNN.typed_signature(generated)==_QTNN.typed_signature(volume)
                model=api.CURRENT[];cache=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
                old=model.meshing.extrude_specs[(3,1)]
                model.meshing.extrude_specs[(3,1)]=merge(old,(;T=(0.,sign(pitch)*.1,0.)))
                params=deepcopy(model.curve_params);attributes=repr(model.meshing)
                payload=(api.mesh.get_nodes(),api.mesh.get_elements())
                counters=(api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())
                @test_throws r"global P1 cell-hull disjointness is not certified" api.mesh.generate(3)
                @test api.CURRENT[]===model && api.LAST_MESH[]===cache &&
                    api.LAST_MESH_CLASS[]===class
                @test (api.mesh.get_nodes(),api.mesh.get_elements())==payload
                @test model.curve_params==params && repr(model.meshing)==attributes
                @test (api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())==counters
            end)
        end
    end

    @testset "a folded Hex8 map is swept as in Gmsh" begin
        vertices=_QTNN.FOLDED_HEX_VERTICES
        coordinates=reduce(hcat,collect.(vertices))
        hex=MixedMesh(coordinates,[ElementBlock(5,reshape(Int32.(1:8),8,1))])
        @test _QTNN.certify_complex(hex).total>0
        determinant=_QTNN.exact_hex_jacobian(vertices,big(63)//64,big(1)//2,big(1)//4096)
        @test determinant<0
        @test Float64(determinant)≈-.015157102010903806 atol=1e-15
        for mode in ("QuadTriNoNewVerts","QuadTriAddVerts")
            source=replace(_QTNN.folded_hex_source(),"QuadTriNoNewVerts"=>mode)
            geometry=_QTNN.execute(source;dim=0)
            model=geometry.model
            params=deepcopy(model.curve_params);attributes=repr(model.meshing)
            # Gmsh 4.15.2 meshes this sweep in both modes: two tetrahedra,
            # five pyramids and one hexahedron whose minimum scaled Jacobian
            # is -0.10, on 13 swept nodes. The Region producer returns the
            # same cells. Gmsh's 33 nodes add one isolated node for each of
            # the 20 Points its screw Extrude puts on the generatrices: those
            # curves end at Points 29, 34, 48 and 62, not at the cap corners,
            # in Gmsh as here, and the entity projection binds curve ends to
            # swept nodes, so the GEO request stops there (broken below;
            # repro in the structured4 scratch directory, p4.jl).
            volume=mesh_model_volume(model,1)
            @test nnodes(volume)==13
            @test _QTNN.counts(volume)==Dict(4=>2,5=>1,7=>5)
            geometry=_QTNN.execute(source;dim=0)
            model=geometry.model
            params=deepcopy(model.curve_params);attributes=repr(model.meshing)
            pipeline=ArgumentError
            @test_throws pipeline _QTNN.execute(source)
            @test_broken try nnodes(_QTNN.execute(source).mesh)==33 catch;false end

            # A failed GEO request also rolls back its preliminary grading,
            # option synchronization, and replacement of existing parts.
            context=Tessella.IO._GeoNumericContext()
            context.option_numbers[("Mesh",0,"CharacteristicLengthFactor")]=2.
            options=deepcopy(context.option_numbers)
            old_part=Tessella.MeshTypes.Mesh(zeros(3,1))
            push!(context.mesh_parts,(0,123,old_part))
            parts=context.mesh_parts;model_value=repr(model)
            @test_throws pipeline Tessella.GeoExec._geo_mesh_model(model,3,context)
            @test repr(model)==model_value
            @test context.option_numbers==options
            @test context.mesh_parts===parts && length(parts)==1 && parts[1][3]===old_part

            _qtnn_api(source,(api,_)->begin
                api.mesh.generate(0)
                model=api.CURRENT[];cache=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
                payload=(api.mesh.get_nodes(),api.mesh.get_elements())
                params=deepcopy(model.curve_params);attributes=repr(model.meshing)
                counters=(api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())
                @test_throws pipeline api.mesh.generate(3)
                @test api.CURRENT[]===model && api.LAST_MESH[]===cache &&
                    api.LAST_MESH_CLASS[]===class
                @test (api.mesh.get_nodes(),api.mesh.get_elements())==payload
                @test model.curve_params==params && repr(model.meshing)==attributes
                @test (api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())==counters
            end)
        end
    end

    @testset "multiple isolated regions retain independent identity" begin
        source=_QTNN.paired_source()
        execution=_QTNN.execute(source)
        left=Int(execution.lists["left"][2]);right=Int(execution.lists["right"][2])
        first=geo_entity_mesh(execution,3,left;geometry=true);second=geo_entity_mesh(execution,3,right;geometry=true)
        free=_QTNN.fixture("uniform3_false","Layers{3}",[0.,1/3,2/3,1.],false)
        recombined=_QTNN.fixture("uniform3_true","Layers{3}",[0.,1/3,2/3,1.],true)
        @test _QTNN.certify_fixture(first,free).ncell>0
        translated=MixedMesh(second.coords.-[3.,0,0],second.blocks)
        @test _QTNN.certify_fixture(translated,recombined).ncell==9
        @test nnodes(execution.mesh)==33
        _qtnn_api(source,(api,geometry)->begin
            # Move the second independent topology with the public point
            # setter, which preserves identity and does not run coherence.
            for (tag,p) in collect(api.CURRENT[].points)
                p[1]>=3 || continue
                api.model.set_coordinates(tag,p[1]-3,p[2],p[3])
            end
            generated=api.mesh.generate(3)
            @test nnodes(generated)==33
            a=api.mesh.get_nodes(3,left,true)[1]
            b=api.mesh.get_nodes(3,right,true)[1]
            @test length(a)==16 && length(b)==17
            @test isempty(intersect(a,b))
            @test length(unique(Tuple(generated.coords[:,i]) for i in 1:nnodes(generated)))==17
            @test sum(length,api.mesh.get_elements(3,left)[2])>0
            @test sum(length,api.mesh.get_elements(3,right)[2])==9
            @test validate(generated).ok
        end)
    end

    @testset "inert lower-dimensional and unrecombined modes" begin
        curve="""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        Transfinite Curve{1}=2;
        Extrude{0,0,1}{Curve{1};Layers{3};Recombine;QuadTriNoNewVerts;}
        """
        execution=_QTNN.execute(curve;dim=2)
        @test sum(size(b.nodes,2) for b in execution.mesh.blocks if b.msh==3)==3
        fixture=_QTNN.fixture("inert","Layers{3}",[0.,1/3,2/3,1.],false)
        source=replace(fixture.source,"Recombine Surface{1};"=>"",
                       "Layers{3};Recombine;"=>"Layers{3};")
        execution=_QTNN.execute(source)
        @test size(geo_entity_mesh(execution,3,1;geometry=true).tets,2)==18
        @test validate(execution.mesh).ok
    end

    @testset "allocation grows with the layer column" begin
        models=[_QTNN.execute(_QTNN.fixture("allocation","Layers{$layers}",
                 collect(range(0.,1.;length=layers+1)),false).source;dim=0).model
                for layers in (8,16)]
        for model in models;mesh_model_volume(model,1);end
        small=minimum(@allocated(mesh_model_volume(models[1],1)) for _ in 1:3)
        large=minimum(@allocated(mesh_model_volume(models[2],1)) for _ in 1:3)
        _QTNN_ALLOCATION_RESULT[]=(small,large)
        @test large<=3small
    end
end
