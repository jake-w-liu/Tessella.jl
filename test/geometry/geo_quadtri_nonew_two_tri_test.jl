using Test,Tessella
using Tessella.Elements: MixedMesh, ElementBlock, msh_dimension
using Tessella.Model: mesh_model_surface,mesh_model_volume,model_to_mixed
using Tessella.MeshTypes: Mesh,nnodes

include("quadtri_nonew_two_tri_certificates.jl")
const _QTT=QuadTriNoNewTwoTriCertificates
const _QTTQ=_QTT.Quad

function _qtt_api(source,action)
    api=Tessella.API;api.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"two_tri_nonew.geo");write(path,source)
            geometry=api.open_geo!(path;mesh_dim=0)
            return action(api,geometry)
        end
    finally
        api.finalize()
    end
end

function _qtt_snapshot(api)
    model=api.CURRENT[]
    return (;model,cache=api.LAST_MESH[],class=api.LAST_MESH_CLASS[],overlay=api.LAST_MESH_HIGH_ORDER[],
        edges=api.LAST_MESH_EDGES[],faces=api.LAST_MESH_FACES[],model_value=repr(model),
        params=deepcopy(model.curve_params),payload=(api.mesh.get_nodes(),api.mesh.get_elements()),
        maxima=(api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag()),
        counters=(api.NODE_TAG_MAX[],api.ELEMENT_TAG_MAX[]))
end

function _qtt_unchanged(api,before)
    @test api.CURRENT[]===before.model
    @test api.LAST_MESH[]===before.cache && api.LAST_MESH_CLASS[]===before.class
    @test api.LAST_MESH_HIGH_ORDER[]===before.overlay
    @test api.LAST_MESH_EDGES[]===before.edges && api.LAST_MESH_FACES[]===before.faces
    @test repr(before.model)==before.model_value && before.model.curve_params==before.params
    @test (api.mesh.get_nodes(),api.mesh.get_elements())==before.payload
    @test (api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())==before.maxima
    @test (api.NODE_TAG_MAX[],api.ELEMENT_TAG_MAX[])==before.counters
end

function _qtt_model_failure(source,diagnostic;mutate=identity)
    geometry=_QTT.execute(source;dim=0)
    mutate(geometry.model)
    parameters=deepcopy(geometry.model.curve_params);value=repr(geometry.model)
    @test_throws diagnostic mesh_model_volume(geometry.model,1)
    @test repr(geometry.model)==value && geometry.model.curve_params==parameters
    context=Tessella.IO._GeoNumericContext()
    sentinel=Mesh(zeros(3,1));push!(context.mesh_parts,(0,777,sentinel))
    parts=context.mesh_parts;options=deepcopy(context.option_numbers)
    @test_throws diagnostic Tessella.GeoExec._geo_mesh_model(geometry.model,3,context)
    @test repr(geometry.model)==value && geometry.model.curve_params==parameters
    @test context.mesh_parts===parts && only(parts)[3]===sentinel
    @test context.option_numbers==options
end

function _qtt_curve_parameters(model,count)
    expected=collect(range(0.,1.;length=count))
    for curve in 1:4
        parameters=model.curve_params[curve]
        @test length(parameters)==count && issorted(parameters) && all(diff(parameters).>0)
        @test first(parameters)==0. && last(parameters)==1.
        @test parameters≈expected atol=8eps(Float64) rtol=0
    end
end

function _qtt_arc_source(plane)
    axes=_QTT.FRAMES[plane]
    planar=((1.,0.),(2.,0.),(0.,2.),(0.,1.),(0.,0.))
    points=join(["Point($tag)={$(join(ntuple(k->k==axes[1] ? p[1] : k==axes[2] ? p[2] : 0.,3),',')),1};" for (tag,p) in enumerate(planar)],"\n")
    return """
    Geometry.AutoCoherence=0;
    $points
    Line(1)={1,2};Circle(2)={2,5,3};Line(3)={3,4};Circle(4)={4,5,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Transfinite Curve{1,2,3,4}=3;Transfinite Surface{1};
    """
end

function _qtt_old_source(category,direction,laterals)
    corners=category===:triangle ? ((0.,0.),(1.,0.),(0.,1.)) :
                                   ((0.,0.),(1.,0.),(1.,1.),(0.,1.))
    count=length(corners);ids=join(1:count,',')
    points=join("Point($tag)={$(p[1]),$(p[2]),0,1};" for (tag,p) in enumerate(corners))
    curves=join("Line($tag)={$tag,$(mod1(tag+1,count))};" for tag in 1:count)
    return "Geometry.AutoCoherence=0;"*points*curves*
        "Curve Loop(1)={$ids};Plane Surface(1)={1};"*
        "Transfinite Curve{$ids}=2;Transfinite Surface{1};"*
        (category===:quadrangle ? "Recombine Surface{1};" : "")*
        "Extrude{0,0,$direction}{Surface{1};Layers{4};Recombine;"*
        "QuadTriNoNewVerts $(laterals ? "RecombLaterals" : "");}"
end

# These affine reference maps are proved with exact actual coordinates.
# Translation to a local frame is used only by the independent face complex
# auditor, whose ordinary Float64 face centers must not round onto a face.
function _qtt_affine_cell_volume(volume)
    total=_QTT.Q(0)
    for block in _QTT.blocks(volume),cell in eachcol(block.nodes)
        msh_dimension(block.msh)==3 || continue
        p=Tuple(_QTT.exact(_QTT.point(volume,node)) for node in cell)
        if block.msh==4
            determinant=_QTT.det(p[2].-p[1],p[3].-p[1],p[4].-p[1])
            @test determinant>0
            total+=determinant/6
        elseif block.msh==6
            shift=p[4].-p[1]
            @test shift==p[5].-p[2]==p[6].-p[3]
            determinant=_QTT.det(p[2].-p[1],p[3].-p[1],shift)
            @test determinant>0
            total+=determinant/2
        elseif block.msh==7
            @test p[3].-p[2]==p[4].-p[1]
            determinant=_QTT.det(p[2].-p[1],p[4].-p[1],p[5].-p[1])
            @test determinant>0
            total+=determinant/3
        else
            error("unexpected old adjacent-plane family $(block.msh)")
        end
    end
    return total
end

@testset "QuadTriNoNewVerts convex two-Tri source" begin
    @testset "packed orientation records distinguish every cyclic face" begin
        labels=Int32[11,37,59,101]
        for width in (3,4)
            permutations=width==3 ?
                [Tuple(labels[[a,b,c]]) for a in 1:3,b in 1:3,c in 1:3 if length(Set((a,b,c)))==3] :
                [Tuple(labels[[a,b,c,d]]) for a in 1:4,b in 1:4,c in 1:4,d in 1:4 if length(Set((a,b,c,d)))==4]
            codes=Dict{Tuple,UInt8}()
            @test length(permutations)==(width==3 ? 6 : 24)
            for face in permutations
                canonical=_QTT.cycle(face)
                recorded=Tessella.Model._extrude_nonew_face_cycle(face)
                code=Tessella.Model._extrude_nonew_face_orientation(recorded,width)
                @test recorded[1:width]==canonical
                @test get(codes,canonical,code)==code
                codes[canonical]=code
                @test all(Tessella.Model._extrude_nonew_face_orientation(
                    Tessella.Model._extrude_nonew_face_cycle(Tuple(circshift(collect(face),shift))),width)==code for shift in 0:width-1)
                reverse_code=Tessella.Model._extrude_nonew_face_orientation(
                    Tessella.Model._extrude_nonew_face_cycle(reverse(face)),width)
                @test reverse_code==xor(code,0x01)
            end
            @test length(codes)==(width==3 ? 2 : 6)
            @test length(unique(values(codes)))==length(codes)
        end
    end

    @testset "opposite oriented incidence across a true shared face" begin
        # Two tetrahedra on opposite sides of ABC have opposite outward
        # cycles. Moving the second apex to the same side keeps both local
        # determinants positive, but cannot create a conforming interior.
        coords=[0. 1. 0. 0. .25;0. 0. 1. 0. .25;0. 0. 0. 1. -1.]
        valid=Int32[1 1;2 3;3 2;4 5]
        for cell in eachcol(valid)
            a,b,c,d=(_QTT.exact(Tuple(coords[:,node])) for node in cell)
            @test _QTT.det(b.-a,c.-a,d.-a)>0
        end
        counts=Dict{NTuple{4,Int32},Tessella.Model._ExtrudeNoNewFaceIncidence}()
        Tessella.Model._extrude_nonew_count_faces!(counts,valid,4,"independent face regression")
        @test count(incidence->incidence.count==2,values(counts))==1
        @test count(incidence->incidence.count==1,values(counts))==6
        invalid=Int32[1 1;2 2;3 3;4 5]
        coords[3,5]=2.
        a,b,c,d=(_QTT.exact(Tuple(coords[:,node])) for node in invalid[:,2])
        @test _QTT.det(b.-a,c.-a,d.-a)>0
        empty!(counts)
        @test_throws r"internal volume face .* equal or inconsistent orientation" Tessella.Model._extrude_nonew_count_faces!(counts,invalid,4,"independent face regression")
    end

    @testset "published immutable packed relations remain unchanged" begin
        @test _QTT.packed_digest(Tessella.Model._extrude_nonew_templates())==
            "3d1c348a08a8645945a6864875ebec7d72c58ca419067279339d8aa0ee079ccc"
        @test _QTT.packed_digest(Tessella.Model._extrude_nonew_prism_templates())==
            "79082312c91d6de90190e7888769ff6f9a1ddf0335d5782e76b2650545d97556"
    end

    @testset "signed actual products on XY/YZ/ZX and grouped layers" begin
        fixtures=_QTT.fixtures()
        @test length(fixtures)==72
        for fixture in fixtures
            @testset "$(fixture.name)" begin
                execution=_QTT.execute(fixture.source)
                volume=geo_entity_mesh(execution,3,1)
                source=only(part for (dim,tag,part) in execution.mesh_parts if dim==2 && tag==1)
                certificate=_QTT.certify(volume,fixture,source)
                n=length(fixture.levels)-1
                @test validate(execution.mesh).ok && validate(volume).ok
                @test nnodes(volume)==4(n+1)
                @test certificate.ncell==(fixture.laterals ? 2n : 6n)
                @test length(certificate.boundary)==(fixture.laterals ? 4n+4 : 8n+4)
                @test length(certificate.internal)==(fixture.laterals ? 3n-2 : 8n-2)
                @test _QTTQ.surface_faces(execution.mesh_parts,volume)==Set(keys(certificate.boundary))
                @test count(part->part[1]==2,execution.mesh_parts)==6
                projected=model_to_mixed(execution.model,volume,3,1)
                @test validate(projected).ok
                @test _QTTQ.typed_signature(projected)==_QTTQ.typed_signature(volume)
                @test _QTTQ.surface_faces([(2,0,projected)],volume)==Set(keys(certificate.boundary))
                owners=projected.entity_data.node_entities
                @test count(owner->owner[1]==0,owners)==8
                @test count(owner->owner[1]==1,owners)==4(n-1)
                @test all(owner->owner[1] in (0,1),owners)
                # Every emitted CAD segment uses a true edge of this volume;
                # the source's interior diagonal has no phantom CAD entity.
                volume_edges=Set{Tuple}()
                for block in _QTT.blocks(volume),cell in eachcol(block.nodes),pattern in _QTT.FACES[Int(block.msh)],i in eachindex(pattern)
                    a=cell[pattern[i]];b=cell[pattern[mod1(i+1,length(pattern))]]
                    push!(volume_edges,_QTTQ.key((_QTT.point(volume,a),_QTT.point(volume,b))))
                end
                line_keys=Tuple[]
                for (dim,_,part) in execution.mesh_parts
                    dim==1 || continue
                    for cell in eachcol(part.segs)
                        push!(line_keys,_QTTQ.key(_QTT.point(part,node) for node in cell))
                    end
                end
                @test all(line->line in volume_edges,line_keys)
            end
        end
    end

    @testset "actual diagonal/frame choices and plane offsets" begin
        for arrangement in ("Left","Right","AlternateLeft","AlternateRight"),pins in (:auto,:opposite),winding in (-1,1),laterals in (false,true)
            fixture=_QTT.fixture("frame";shape=:trapezoid,arrangement,pins,winding,laterals,
                height=-2.,offset=(.125,-.5,2.25),layers=:nonbinary)
            geometry=_QTT.execute(fixture.source;dim=0)
            source=mesh_model_surface(geometry.model,1)
            volume=mesh_model_volume(geometry.model,1)
            certificate=_QTT.certify(volume,fixture,source)
            @test certificate.source.area==2
            @test certificate.total==4
            @test validate(volume).ok
            again=mesh_model_volume(geometry.model,1)
            @test _QTTQ.typed_signature(again)==_QTTQ.typed_signature(volume)
        end
        for plane in (:YZ,:ZX),laterals in (false,true)
            fixture=_QTT.fixture("offset_$(plane)_$(laterals)";plane,shape=:skew,height=-1.25,
                winding=-1,laterals,offset=(1.5,-2.25,.125),pins=:opposite)
            execution=_QTT.execute(fixture.source)
            source=only(part for (dim,tag,part) in execution.mesh_parts if dim==2 && tag==1)
            @test _QTT.certify(geo_entity_mesh(execution,3,1),fixture,source).total==15//4
        end
    end

    @testset "flexible scaling retains the actual four source endpoints" begin
        for (factor,stored_count) in ((2,4),(3,6)),laterals in (false,true)
            fixture=_QTT.fixture("flexible_$(factor)_$(laterals)";laterals,layers=:one)
            source="Mesh.FlexibleTransfinite=1;Mesh.MeshSizeFactor=$factor;\n"*
                replace(fixture.source,"Transfinite Curve{1,2,3,4}=2"=>"Transfinite Curve{1,2,3,4}=$stored_count")
            execution=_QTT.execute(source)
            surface=only(part for (dim,tag,part) in execution.mesh_parts if dim==2 && tag==1)
            @test nnodes(surface)==4 && size(surface.tris,2)==2
            @test _QTT.certify(geo_entity_mesh(execution,3,1),fixture,surface).ncell==(laterals ? 2 : 6)
        end
    end

    @testset "standalone caps/laterals, GEO Mesh2 and API top product agree" begin
        for plane in (:XY,:YZ,:ZX),laterals in (false,true)
            fixture=_QTT.fixture("entry_$(plane)_$(laterals)";plane,shape=:skew,winding=-1,
                height=-1.25,laterals,layers=:graded,pins=:opposite,offset=(.5,-.25,1.))
            geometry=_QTT.execute(fixture.source;dim=0)
            tags=[1,Int(geometry.lists["sweep"][1]),Int.(geometry.lists["sweep"][3:end])...]
            @test length(tags)==6 && allunique(tags)
            before=[(2,tag,mesh_model_surface(geometry.model,tag)) for tag in tags]
            source=only(part for (_,tag,part) in before if tag==1)
            volume=mesh_model_volume(geometry.model,1)
            certificate=_QTT.certify(volume,fixture,source)
            @test _QTTQ.surface_faces(before,volume)==Set(keys(certificate.boundary))
            surface_geo=_QTT.execute(fixture.source;dim=2)
            @test _QTTQ.surface_faces(surface_geo.mesh_parts,volume)==Set(keys(certificate.boundary))
            complete=_QTT.execute(fixture.source)
            @test _QTTQ.typed_signature(geo_entity_mesh(complete,3,1))==_QTTQ.typed_signature(volume)
            _qtt_api(fixture.source,(api,_)->begin
                generated=api.mesh.generate(3)
                actual=_QTT.public_volume(api)
                @test _QTT.certify(actual,fixture,source).total==certificate.total
                @test _QTTQ.typed_signature(actual)==_QTTQ.typed_signature(volume)
                @test validate(generated).ok
                @test api.mesh.get_elements(2)[1]==Int32[]
                @test length(api.mesh.get_nodes(3,1,true)[1])==16
                @test api.mesh.get_elements(3,1)[1]==[laterals ? 6 : 4]
                old=_qtt_snapshot(api)
                api.mesh.generate(3)
                @test _QTTQ.typed_signature(_QTT.public_volume(api))==_QTTQ.typed_signature(actual)
                @test api.mesh.get_nodes()==old.payload[1]
            end)
        end
    end

    @testset "physical metadata and source/top identity" begin
        for laterals in (false,true)
            fixture=_QTT.fixture("physical";laterals,shape=:trapezoid)
            source=fixture.source*"""
            Physical Surface("source",101)={1};
            Physical Surface("top",102)={sweep[0]};
            Physical Volume("body",103)={sweep[1]};
            """
            execution=_QTT.execute(source);volume=geo_entity_mesh(execution,3,1)
            projected=model_to_mixed(execution.model,volume,3,1)
            @test projected.physical_names==Dict((2,101)=>"source",(2,102)=>"top",(3,103)=>"body")
            @test all(all(==(103),block.tags) for block in projected.blocks if msh_dimension(block.msh)==3)
            @test Set(owner[2] for owner in projected.entity_data.node_entities if owner[1]==0)==
                Set(tag for (dim,tag,_) in execution.mesh_parts if dim==0)
        end
    end

    @testset "independent coincident source identities remain separate" begin
        fixture=_QTT.fixture("paired";layers=:three)
        prefix=first(split(fixture.source,"sweep[]="))
        paired=prefix*"""
        Point(101)={3,0,0,1};Point(102)={4,0,0,1};
        Point(103)={4,1,0,1};Point(104)={3,1,0,1};
        Line(101)={101,102};Line(102)={102,103};Line(103)={103,104};Line(104)={104,101};
        Curve Loop(101)={101,102,103,104};Plane Surface(101)={101};
        Transfinite Curve{101,102,103,104}=2;Transfinite Surface{101} Left;
        left[]=Extrude{0,0,1}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{101};Layers{3};Recombine;QuadTriNoNewVerts RecombLaterals;};
        """
        execution=_QTT.execute(paired)
        @test nnodes(execution.mesh)==32
        _qtt_api(paired,(api,geometry)->begin
            for (tag,p) in collect(api.CURRENT[].points)
                p[1]>=3 || continue
                api.model.set_coordinates(tag,p[1]-3,p[2],p[3])
            end
            generated=api.mesh.generate(3)
            left=Int(geometry.lists["left"][2]);right=Int(geometry.lists["right"][2])
            a=api.mesh.get_nodes(3,left,true)[1];b=api.mesh.get_nodes(3,right,true)[1]
            @test length(a)==16 && length(b)==16 && isempty(intersect(a,b))
            @test nnodes(generated)==32
            @test length(unique(_QTT.point(generated,i) for i in 1:nnodes(generated)))==16
            for (entity,source_tag,laterals) in ((left,1,false),(right,101,true))
                source=mesh_model_surface(api.CURRENT[],source_tag)
                expected=_QTT.fixture("coincident";laterals)
                @test _QTT.certify(_QTT.public_volume(api,entity),expected,source).ncell==(laterals ? 6 : 18)
            end
            @test validate(generated).ok
        end)
    end

    @testset "exact translated and small source coordinates retain native parameters" begin
        for plane in (:XY,:YZ,:ZX),direction in (-1,1),laterals in (false,true)
            fixture=_QTT.fixture("exact_offset";plane,laterals,layers=:one)
            geometry=_QTT.execute(fixture.source;dim=0)
            axis=fixture.axes[3];offset=direction*2.0^50
            for (tag,p) in collect(geometry.model.points)
                geometry.model.points[tag]=ntuple(k->p[k]+(k==axis ? offset : 0.),3)
            end
            corners=Tuple(ntuple(k->p[k]+(k==axis ? offset : 0.),3) for p in fixture.corners)
            shifted=merge(fixture,(;corners))
            surface=mesh_model_surface(geometry.model,1)
            @test _QTT.source_complex(surface,shifted).area==1
            @test Set(_QTT.point(surface,i) for i in 1:4)==Set(corners)
            _qtt_curve_parameters(geometry.model,2)
            volume=mesh_model_volume(geometry.model,1)
            @test _QTT.certify(volume,shifted,surface).total==1
            @test validate(volume).ok
        end
        for plane in (:XY,:YZ,:ZX),power in (-40,-48),count in (2,3,4)
            fixture=_QTT.fixture("small_source";plane)
            geometry=_QTT.execute(first(split(fixture.source,"sweep[]="));dim=0)
            for (tag,p) in collect(geometry.model.points)
                geometry.model.points[tag]=Tuple(x*2.0^power for x in p)
            end
            for curve in 1:4;Tessella.Model.set_transfinite_curve!(geometry.model,curve,count);end
            surface=mesh_model_surface(geometry.model,1)
            @test nnodes(surface)==count^2 && size(surface.tris,2)==2(count-1)^2
            @test validate(surface).ok
            _qtt_curve_parameters(geometry.model,count)
            projected=model_to_mixed(geometry.model,surface,2,1)
            @test Set(_QTT.point(projected,i) for i in 1:nnodes(projected))==Set(_QTT.point(surface,i) for i in 1:nnodes(surface))
            again=mesh_model_surface(geometry.model,1)
            @test again.coords==surface.coords && again.tris==surface.tris
            _qtt_curve_parameters(geometry.model,count)
        end
        for plane in (:XY,:YZ,:ZX),direction in (-1,1)
            source=_qtt_arc_source(plane);axes=_QTT.FRAMES[plane]
            baseline=_QTT.execute(source;dim=0)
            ordinary=mesh_model_surface(baseline.model,1)
            geometry=_QTT.execute(source;dim=0)
            for (tag,p) in collect(geometry.model.points)
                geometry.model.points[tag]=ntuple(k->p[k]+(k==axes[3] ? direction*2.0^50 : 0.),3)
            end
            surface=mesh_model_surface(geometry.model,1)
            @test nnodes(surface)==9 && size(surface.tris,2)==8 && validate(surface).ok
            @test surface.coords[collect(axes[1:2]),:]==ordinary.coords[collect(axes[1:2]),:]
            @test all(==(direction*2.0^50),surface.coords[axes[3],:])
            @test surface.tris==ordinary.tris
            _qtt_curve_parameters(geometry.model,3)
            projected=model_to_mixed(geometry.model,surface,2,1)
            @test validate(projected).ok && nnodes(projected)==9
            again=mesh_model_surface(geometry.model,1)
            @test again.coords==surface.coords && again.tris==surface.tris
            _qtt_curve_parameters(geometry.model,3)
        end
        fixture=_QTT.fixture("public_exact_source")
        _qtt_api(first(split(fixture.source,"sweep[]=")),(api,_)->begin
            for tag in 1:4
                p=api.CURRENT[].points[tag]
                api.model.set_coordinates(tag,p[1],p[2],2.0^50)
            end
            generated=api.mesh.generate(2)
            @test nnodes(generated)==4 && size(generated.tris,2)==2
            @test all(==(2.0^50),generated.coords[3,:])
            @test Set(_QTT.point(generated,i) for i in 1:4)==Set((p[1],p[2],2.0^50) for p in fixture.corners)
            _qtt_curve_parameters(api.CURRENT[],2)
        end)
    end

    @testset "graded native samples survive endpoint snapping at every local scale" begin
        for plane in (:XY,:YZ,:ZX),radius in (1.,1e8,2.0^-48)
            source=replace(_qtt_arc_source(plane),
                "Transfinite Curve{1,2,3,4}=3;"=>
                "Transfinite Curve{1,3}=2;Transfinite Curve{2}=3 Using Progression 1000000;Transfinite Curve{4}=3 Using Progression 0.000001;")
            geometry=_QTT.execute(source;dim=0)
            for (tag,p) in collect(geometry.model.points)
                geometry.model.points[tag]=Tuple(radius*x for x in p)
            end
            surface=mesh_model_surface(geometry.model,1)
            @test nnodes(surface)==6 && size(surface.tris,2)==4 && validate(surface).ok
            parameters=deepcopy(geometry.model.curve_params)
            fraction=1/(1+1e6)
            @test parameters[2][2]≈fraction atol=8eps(Float64)*fraction rtol=0
            @test parameters[4][2]≈1-fraction atol=8eps(Float64) rtol=0
            for curve in 1:4
                @test first(parameters[curve])==0 && last(parameters[curve])==1
                @test issorted(parameters[curve]) && all(diff(parameters[curve]).>0)
            end
            projected=model_to_mixed(geometry.model,surface,2,1)
            @test validate(projected).ok && nnodes(projected)==6
            @test count(owner->owner[1]==0,projected.entity_data.node_entities)==4
            @test count(owner->owner[1]==1,projected.entity_data.node_entities)==2
            again=mesh_model_surface(geometry.model,1)
            @test again.coords==surface.coords && again.tris==surface.tris
            @test geometry.model.curve_params==parameters
        end
        for plane in (:XY,:YZ,:ZX),ratio in (1e6,1e15),scale in (1.,2.0^-48)
            fixture=_QTT.fixture("graded_straight";plane)
            source=replace(first(split(fixture.source,"sweep[]=")),
                "Transfinite Curve{1,2,3,4}=2;"=>
                "Transfinite Curve{2,4}=2;Transfinite Curve{1}=3 Using Progression $ratio;Transfinite Curve{3}=3 Using Progression $(1/ratio);")
            geometry=_QTT.execute(source;dim=0)
            for (tag,p) in collect(geometry.model.points)
                geometry.model.points[tag]=Tuple(scale*x for x in p)
            end
            surface=mesh_model_surface(geometry.model,1)
            @test nnodes(surface)==6 && size(surface.tris,2)==4 && validate(surface).ok
            parameters=deepcopy(geometry.model.curve_params)
            # Native F_Transfinite integrates and inverts a sampled density.
            # These saved Gmsh 4.15.2 controls include its extreme endpoint
            # partition behavior; analytic grading remains a separate API.
            forward=ratio==1e6 ? 1.0074495582442956e-6 : 7.450581596923805e-9
            reverse=ratio==1e6 ? 0.9999999925474284 : 0.499999996273405
            @test parameters[1][2]≈forward atol=2e-11 rtol=0
            @test parameters[3][2]≈reverse atol=2e-11 rtol=0
            @test 0<parameters[1][2]<1 && 0<parameters[3][2]<1
            for curve in 1:4
                @test first(parameters[curve])==0 && last(parameters[curve])==1
                @test issorted(parameters[curve]) && all(diff(parameters[curve]).>0)
            end
            projected=model_to_mixed(geometry.model,surface,2,1)
            @test validate(projected).ok && nnodes(projected)==6
            for curve in (1,3)
                owned=findall(==((Int8(1),Int32(curve))),projected.entity_data.node_entities)
                @test length(owned)==1
                @test _QTT.point(projected,only(owned))!=geometry.model.points[geometry.model.curves[curve][1]]
                @test _QTT.point(projected,only(owned))!=geometry.model.points[geometry.model.curves[curve][2]]
            end
            @test count(owner->owner[1]==0,projected.entity_data.node_entities)==4
            again=mesh_model_surface(geometry.model,1)
            @test again.coords==surface.coords && again.tris==surface.tris
            @test geometry.model.curve_params==parameters
        end
    end

    @testset "actual one-ULP layer intervals need no representable centroid" begin
        for plane in (:XY,:YZ,:ZX),direction in (-1,1),laterals in (false,true)
            fixture=_QTT.fixture("one_ulp_layers";plane,laterals,layers=:one,height=direction)
            source=replace(fixture.source,"Layers{1}"=>"Layers{4}")
            geometry=_QTT.execute(source;dim=0)
            axis=fixture.axes[3];offset=2.0^50 + 2.0
            for (tag,p) in collect(geometry.model.points)
                geometry.model.points[tag]=ntuple(k->p[k]+(k==axis ? offset : 0.),3)
            end
            corners=Tuple(ntuple(k->p[k]+(k==axis ? offset : 0.),3) for p in fixture.corners)
            actual=merge(fixture,(;corners,levels=collect(0.:.25:1.),source))
            surface=mesh_model_surface(geometry.model,1)
            volume=mesh_model_volume(geometry.model,1)
            certificate=_QTT.certify(volume,actual,surface)
            @test certificate.total==1 && certificate.ncell==(laterals ? 8 : 24)
            @test nnodes(volume)==20 && validate(volume).ok
            @test all(diff(certificate.heights).==direction*.25)
            # Every exact midpoint lies strictly inside its mathematical
            # interval, but its Float64 rounding equals an endpoint.
            for layer in 1:4
                a,b=certificate.heights[layer:layer+1]
                midpoint=(_QTT.Q(a)+_QTT.Q(b))/2
                @test min(_QTT.Q(a),_QTT.Q(b))<midpoint<max(_QTT.Q(a),_QTT.Q(b))
                @test Float64(midpoint) in (a,b)
            end
            projected=model_to_mixed(geometry.model,volume,3,1)
            @test _QTTQ.typed_signature(projected)==_QTTQ.typed_signature(volume)
            @test count(owner->owner[1]==0,projected.entity_data.node_entities)==8
            @test count(owner->owner[1]==1,projected.entity_data.node_entities)==12
        end
    end

    @testset "closed curved boundaries retain a local extent without samples" begin
        model=GeoModel();add_cylinder!(model,0,0,0,0,0,2,1.0e8)
        curve=1;seam=model.points[model.curves[curve][1]]
        coordinates=hcat(collect(seam),[0.,1e8,2.],[-1e8,0.,2.],
                         [0.,-1e8,2.],[0.,0.,2.])
        surface=Mesh(coordinates;tris=Int32[5 5 5 5;1 2 3 4;2 3 4 1])
        eligible,edges=Tessella.Model._surface_boundary_topology(surface,"closed boundary regression")
        @test validate(surface).ok && isempty(model.curve_params)
        @test model.curves[curve][1]==model.curves[curve][2]
        entries,tolerance=Tessella.Model._curve_parameter_nodes_curved(
            model,surface,curve,eligible,edges,1e-12,"closed boundary regression")
        @test last.(entries)==collect(1:4)
        @test first.(entries)≈[0.,pi/2,pi,3pi/2] atol=8eps(Float64) rtol=0
        @test issorted(first.(entries)) && all(diff(first.(entries)).>0)
        @test 0<tolerance<pi/4
        @test isempty(model.curve_params)
        for (parameter,node) in entries
            point=Tessella.Model._model_curve_point(model,curve,parameter,"closed boundary regression")
            @test maximum(abs.(point.-Tuple(coordinates[:,node])))<1e-6
        end
        again,_=Tessella.Model._curve_parameter_nodes_curved(
            model,surface,curve,eligible,edges,1e-12,"closed boundary repeat")
        @test again==entries && isempty(model.curve_params)
        broken=copy(edges);delete!(broken,(Int32(3),Int32(4)))
        @test_throws r"nodes do not form a mesh-edge chain" Tessella.Model._curve_parameter_nodes_curved(
            model,surface,curve,eligible,broken,1e-12,"closed boundary broken chain")
        # Exact independently supplied native samples close to the seam must
        # remain distinct even when the admission tolerance spans the circle.
        delta=1e-6
        close_coords=hcat(collect(seam),[1e8*cos(delta),1e8*sin(delta),2.],
            [0.,1e8,2.],[-1e8,0.,2.],[0.,-1e8,2.],[0.,0.,2.])
        close_surface=Mesh(close_coords;tris=Int32[6 6 6 6 6;1 2 3 4 5;2 3 4 5 1])
        eligible,edges=Tessella.Model._surface_boundary_topology(close_surface,"closed near-seam regression")
        @test validate(close_surface).ok
        for admission in (1e-12,1e9)
            close_entries,snap_tolerance=Tessella.Model._curve_parameter_nodes_curved(
                model,close_surface,curve,eligible,edges,admission,"closed near-seam regression")
            @test last.(close_entries)==collect(1:5)
            @test close_entries[2][1]≈delta atol=8eps(Float64)*delta rtol=0
            @test all(diff(first.(close_entries)).>0)
            @test 0<=snap_tolerance<delta/4
            @test isempty(model.curve_params)
        end
        drifted_coords=copy(close_coords);drifted_coords[2,1]=1.0
        drifted_surface=Mesh(drifted_coords;tris=close_surface.tris)
        @test validate(drifted_surface).ok
        drifted_entries,disabled_snap=Tessella.Model._curve_parameter_nodes_curved(
            model,drifted_surface,curve,eligible,edges,1e9,"closed admitted endpoint drift")
        @test 0<first(drifted_entries)[1]<delta
        @test disabled_snap==0.0 && all(diff(first.(drifted_entries)).>0)
        cap=only(tag for (tag,loops) in model.surfaces if length(loops)==1 &&
            length(model.loops[only(loops)])==1 && abs(only(model.loops[only(loops)]))==curve)
        projected=model_to_mixed(model,close_surface,2,cap)
        @test validate(projected).ok && nnodes(projected)==6
        @test count(owner->owner[1]==0,projected.entity_data.node_entities)==1
        @test count(owner->owner[1]==1,projected.entity_data.node_entities)==4
    end

    @testset "retained one-cell templates use an exact logical witness" begin
        offset=2.0^50+2.0
        for category in (:triangle,:quadrangle),direction in (-1,1),laterals in (false,true)
            source_text=_qtt_old_source(category,direction,laterals)
            shift_model=model->begin
                for (tag,p) in collect(model.points)
                    model.points[tag]=(p[1],p[2],p[3]+offset)
                end
                model
            end
            geometry=_QTT.execute(source_text;dim=0);shift_model(geometry.model)
            surface=mesh_model_surface(geometry.model,1)
            width=category===:triangle ? 3 : 4
            @test nnodes(surface)==width && validate(surface).ok
            @test all(==(offset),surface.coords[3,:])
            if category===:quadrangle && laterals
                # This branch emits the final centroid. Its adjacent planes
                # have no representable interior vertex, so it must reject.
                lower,upper=sort([offset+.75direction,offset+direction])
                @test nextfloat(lower)==upper
                @test lower<(_QTT.Q(lower)+_QTT.Q(upper))/2<upper
                _qtt_model_failure(source_text,r"centroid fan has a folded or degenerate face";mutate=shift_model)
                _qtt_api(source_text,(api,_)->begin
                    for (tag,p) in collect(api.CURRENT[].points)
                        api.model.set_coordinates(tag,p[1],p[2],p[3]+offset)
                    end
                    before=_qtt_snapshot(api)
                    @test_throws r"centroid fan has a folded or degenerate face" api.mesh.generate(3)
                    _qtt_unchanged(api,before)
                end)
                continue
            end
            volume=mesh_model_volume(geometry.model,1)
            expected_volume=category===:triangle ? _QTT.Q(1)/2 : _QTT.Q(1)
            @test nnodes(volume)==5width && validate(volume).ok
            @test _qtt_affine_cell_volume(volume)==expected_volume
            @test Set(volume.coords[3,:])==Set(offset.+direction.*collect(0.:.25:1.))
            local_coords=copy(volume.coords);local_coords[3,:].-=offset
            local_volume=MixedMesh(local_coords,_QTT.blocks(volume))
            certificate=_QTTQ.certify_complex(local_volume)
            @test certificate.total≈Float64(expected_volume) atol=2e-14 rtol=0
            if category===:triangle
                @test _QTTQ.counts(volume)==Dict((laterals ? 6 : 4)=>(laterals ? 4 : 12))
            end
            projected=model_to_mixed(geometry.model,volume,3,1)
            @test _QTTQ.typed_signature(projected)==_QTTQ.typed_signature(volume)
            @test count(owner->owner[1]==0,projected.entity_data.node_entities)==2width
            @test count(owner->owner[1]==1,projected.entity_data.node_entities)==3width
            context=Tessella.IO._GeoNumericContext()
            merged,_=Tessella.GeoExec._geo_mesh_model(geometry.model,3,context)
            emitted=only(part for (dim,tag,part) in context.mesh_parts if dim==3 && tag==1)
            @test validate(merged).ok
            @test _QTTQ.typed_signature(emitted)==_QTTQ.typed_signature(volume)
            _qtt_api(source_text,(api,_)->begin
                for (tag,p) in collect(api.CURRENT[].points)
                    api.model.set_coordinates(tag,p[1],p[2],p[3]+offset)
                end
                generated=api.mesh.generate(3)
                actual=_QTT.public_volume(api,1)
                @test validate(generated).ok && nnodes(actual)==5width
                @test _qtt_affine_cell_volume(actual)==expected_volume
                @test _QTTQ.typed_signature(actual)==_QTTQ.typed_signature(volume)
                @test length(api.mesh.get_nodes(3,1,true)[1])==5width
            end)
        end
    end

    @testset "malformed products and unsupported constraints reject atomically" begin
        fixture=_QTT.fixture("guard";laterals=true)
        # Unsupported higher grids, cap overrides and non-normal translation
        # must be diagnosed before grading or publishing any output.
        invalid=(
            (replace(fixture.source,"Transfinite Curve{1,2,3,4}=2"=>"Transfinite Curve{1,2,3,4}=3"),r"multiple source cells require the boundary-category planner"),
            (replace(fixture.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{0.25,0,1}"),r"two-triangle source requires exactly one nonzero translation component"),
            (replace(fixture.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{{0,0,1},{-1,0,0},Pi/6}"),r"QuadTriNoNewVerts.*translation"),
            (fixture.source*"Transfinite Surface{sweep[0]};\n",r"transfinite cap/lateral overrides"))
        for (source,diagnostic) in invalid
            _qtt_model_failure(source,diagnostic)
        end
        # The literal CAD construction remains valid. A finite caller edit
        # makes this source concave before the planner is invoked.
        concave=_QTT.fixture("concave_guard";laterals=true,arrangement="Right")
        _qtt_model_failure(concave.source,r"QuadTriNoNewVerts.*convex";mutate=model->(model.points[3]=(.25,.25,0.)))
        # The source at z=1 has four exact endpoints. A positive requested
        # sub-ULP translation rounds all subsequent planes onto that source.
        # Edit the existing shell to preserve its distinct CAD identities.
        rounded=_QTT.fixture("rounded_guard";laterals=true,offset=(0.,0.,1.))
        _qtt_model_failure(rounded.source,r"actual normal planes collapse or reverse after rounding";mutate=model->begin
            for (tag,p) in collect(model.points)
                p[3]>1. && (model.points[tag]=(p[1],p[2],1.0 + 2.0^-54))
            end
            for (key,spec) in collect(model.meshing.extrude_specs)
                spec.type===:translate && (model.meshing.extrude_specs[key]=(type=:translate,T=(0.,0.,2.0^-54)))
            end
        end)
        # A real stacked extrusion has a copied source and a shared cap.
        # These are different unsupported contracts, and each entrypoint
        # reports the first applicable one before it publishes grading.
        chained=fixture.source*"""
        next[]=Extrude{0,0,1}{Surface{sweep[0]};Layers{1};Recombine;QuadTriNoNewVerts;};
        """
        geometry=_QTT.execute(chained;dim=0)
        second_volume=Int(geometry.lists["next"][2])
        before=repr(geometry.model);params=deepcopy(geometry.model.curve_params)
        @test_throws r"copied or chained sources require the dependency planner" mesh_model_volume(geometry.model,second_volume)
        @test repr(geometry.model)==before && geometry.model.curve_params==params
        context=Tessella.IO._GeoNumericContext()
        sentinel=Mesh(zeros(3,1));push!(context.mesh_parts,(0,777,sentinel))
        parts=context.mesh_parts;options=deepcopy(context.option_numbers)
        @test_throws r"shared surfaces require the neighboring-region planner" Tessella.GeoExec._geo_mesh_model(geometry.model,3,context)
        @test repr(geometry.model)==before && geometry.model.curve_params==params
        @test context.mesh_parts===parts && only(parts)[3]===sentinel
        @test context.option_numbers==options

        # The catalog consumes the actual source topology. A malformed or
        # mixed source cannot be silently padded into the bounded product.
        geometry=_QTT.execute(fixture.source;dim=0)
        model=geometry.model;source=mesh_model_surface(model,1)
        levels,refs=Tessella.Model._extrude_nonew_levels(model.meshing.extrude[(3,1)],"independent catalog guard")
        spec=model.meshing.extrude_specs[(3,1)]
        before=repr(model);params=deepcopy(model.curve_params)
        repeated=copy(source.tris);repeated[1,1]=repeated[2,1]
        reversed=copy(source.tris);reversed[1,1],reversed[2,1]=reversed[2,1],reversed[1,1]
        mixed=MixedMesh(source.coords,ElementBlock[ElementBlock(2,source.tris[:,1:1]),ElementBlock(3,reshape(Int32[1,2,3,4],4,1))])
        for (actual,diagnostic) in ((Mesh(source.coords;tris=repeated),r"repeated or foreign triangle nodes"),
                                   (Mesh(source.coords;tris=reversed),r"triangle winding must agree"),
                                   (mixed,r"exactly four nodes and two linear triangles"))
            @test_throws diagnostic Tessella.Model._extrude_nonew_two_tri_catalog(model,1,actual,spec,levels,refs,true,"independent catalog guard")
            @test repr(model)==before && model.curve_params==params
        end
        _qtt_api(fixture.source,(api,_)->begin
            api.mesh.generate(3)
            for curve in 1:4;api.mesh.set_transfinite_curve(curve,3);end
            before=_qtt_snapshot(api)
            @test_throws r"multiple source cells require the boundary-category planner" api.mesh.generate(3)
            _qtt_unchanged(api,before)
        end)
    end
end
