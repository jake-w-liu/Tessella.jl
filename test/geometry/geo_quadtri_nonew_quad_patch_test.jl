using Test,Tessella
using Tessella.Elements: MixedMesh,ElementBlock,msh_dimension
using Tessella.MeshTypes: Mesh,nnodes
using Tessella.Model: mesh_model_surface,mesh_model_volume,model_to_mixed

if !isdefined(@__MODULE__,:QuadTriNoNewQuadPatchCertificates)
    include("quadtri_nonew_quad_patch_certificates.jl")
end
const _QPC=QuadTriNoNewQuadPatchCertificates

function _qpc_api(source,action)
    api=Tessella.API;api.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"quad_patch.geo");write(path,source)
            result=api.open_geo!(path;mesh_dim=0)
            return action(api,result)
        end
    finally
        api.finalize()
    end
end

function _qpc_api_snapshot(api)
    model=api.CURRENT[]
    return (;model,model_value=repr(model),params=deepcopy(model.curve_params),
        cache=api.LAST_MESH[],class=api.LAST_MESH_CLASS[],overlay=api.LAST_MESH_HIGH_ORDER[],
        edges=api.LAST_MESH_EDGES[],faces=api.LAST_MESH_FACES[],
        payload=(api.mesh.get_nodes(),api.mesh.get_elements()),
        history=(api.NODE_TAG_MAX[],api.ELEMENT_TAG_MAX[]))
end

function _qpc_api_unchanged(api,before)
    @test api.CURRENT[]===before.model && repr(before.model)==before.model_value
    @test before.model.curve_params==before.params
    @test api.LAST_MESH[]===before.cache && api.LAST_MESH_CLASS[]===before.class
    @test api.LAST_MESH_HIGH_ORDER[]===before.overlay
    @test api.LAST_MESH_EDGES[]===before.edges && api.LAST_MESH_FACES[]===before.faces
    @test (api.mesh.get_nodes(),api.mesh.get_elements())==before.payload
    @test (api.NODE_TAG_MAX[],api.ELEMENT_TAG_MAX[])==before.history
end

function _qpc_atomic(source,diagnostic)
    geometry=_QPC.execute(source;dim=0)
    model=geometry.model;before=repr(model);params=deepcopy(model.curve_params)
    volume=Int(geometry.lists["sweep"][2])
    @test_throws diagnostic mesh_model_volume(model,volume)
    @test repr(model)==before && model.curve_params==params
    context=Tessella.IO._GeoNumericContext()
    sentinel=Mesh(zeros(3,1));push!(context.mesh_parts,(0,90001,sentinel))
    parts=context.mesh_parts;options=copy(context.option_numbers)
    @test_throws diagnostic Tessella.GeoExec._geo_mesh_model(model,3,context)
    @test repr(model)==before && model.curve_params==params
    @test context.mesh_parts===parts && only(parts)[3]===sentinel
    @test context.option_numbers==options
end

function _qpc_projection(geometry,f,source,volume,certificate)
    out=Int.(geometry.lists["sweep"]);dim=out[2]
    projected=model_to_mixed(geometry.model,volume,3,dim)
    @test Tessella.validate(projected).ok
    @test _QPC.Quad.typed_signature(projected)==_QPC.Quad.typed_signature(volume)
    @test haskey(projected.entity_data.entities,(3,dim))
    @test Set(abs.(projected.entity_data.entities[(3,dim)].boundaries))==Set((f.surface,out[1],out[3:end]...))
    @test all(projected.entity_data.node_entities[v][1] in 0:3 for v in axes(projected.coords,2))
    parts=Tuple{Int,Int,Any}[]
    for surface in (f.surface,out[1],out[3:end]...)
        actual=mesh_model_surface(geometry.model,surface)
        push!(parts,(2,surface,actual))
        expected=surface==f.surface ? Dict(3=>4) : surface==out[1] ? Dict(2=>8) :
            f.laterals ? Dict(3=>2(length(f.levels)-1)) : Dict(2=>4(length(f.levels)-1))
        families=Dict(Int(b.msh)=>size(b.nodes,2) for b in _QPC.blocks(actual))
        @test families==expected
        @test nnodes(actual)==(surface in (f.surface,out[1]) ? 9 : 3length(f.levels))
        projected_surface=model_to_mixed(geometry.model,actual,2,surface)
        @test Tessella.validate(projected_surface).ok
        @test all(all(==(surface),projected_surface.elementary_entities[i]) for
                  (i,b) in enumerate(projected_surface.blocks) if msh_dimension(b.msh)==2)
    end
    @test _QPC.Quad.surface_faces(parts,volume)==Set(keys(certificate.boundary))
    @test _QPC.Quad.surface_faces(geometry.mesh_parts,volume)==Set(keys(certificate.boundary))
    # Every three-node source chain and copied cap chain is present as actual
    # Line1 connectivity. This checks TF3 middles rather than endpoint closure.
    lines=[b for b in projected.blocks if b.msh==1]
    @test length(lines)==1
    line_owners=projected.elementary_entities[findfirst(b->b.msh==1,projected.blocks)]
    for curve in f.curve_tags
        @test count(==(curve),line_owners)==2
        selected=findall(==(curve),line_owners)
        ids=Set(vec(only(lines).nodes[:,selected]))
        @test length(ids)==3
        @test count(v->projected.entity_data.node_entities[v]==(1,curve),ids)==1
    end
    return projected
end

function _qpc_source_progression_case(laterals)
    base=_QPC.fixture("source_progression";laterals,layers=:one)
    fixture=merge(base,(;source=replace(base.source,
        "Transfinite Curve{1,2,3,4}=3;"=>"Transfinite Curve{1,2,3,4}=3 Using Progression 4;")))
    geometry=_QPC.execute(fixture.source;dim=3)
    source=mesh_model_surface(geometry.model,fixture.surface)
    volume=geo_entity_mesh(geometry,3,Int(geometry.lists["sweep"][2]))
    certificate=_QPC.certify(volume,fixture,source)
    @test certificate.total==1
    # Saved Gmsh 4.15.2 native density samples; the stepwise primitive is
    # numerically inverted and does not place this node at the analytic .2.
    samples=(0.2000000039104483,0.2000000039104483,
             0.20000000391155554,0.20000000391155554)
    for curve in fixture.curve_tags
        @test geometry.model.curve_params[curve]≈[0.,samples[curve],1.] atol=2e-11 rtol=0
    end
    _qpc_projection(geometry,fixture,source,volume,certificate)
    complex=_QPC.source_complex(source,fixture)
    _qpc_api(fixture.source,(api,opened)->begin
        api.mesh.generate(3);api.mesh.set_order(2)
        top,body=Int.(opened.lists["sweep"])[1:2]
        tags,xyz,_=api.mesh.get_nodes(3,body,true)
        positions=Dict(Tuple(p)=>tag for (tag,p) in zip(tags,eachcol(reshape(xyz,3,:))))
        faces=reshape(api.mesh.get_element_face_nodes(14,4,body,false),9,:)
        edge_mids=Dict{Tuple,UInt64}()
        for msh in api.mesh.get_element_types(3,body)
            for edge in eachcol(reshape(api.mesh.get_element_edge_nodes(msh,body,false),3,:))
                key=_QPC.key(edge[1:2])
                @test !haskey(edge_mids,key) || edge_mids[key]==edge[3]
                edge_mids[key]=edge[3]
            end
        end
        for cell in complex.cells
            points=[_QPC.point(source,v) for v in cell]
            support=_QPC.key(positions[p] for p in points)
            matches=[face[9] for face in eachcol(faces) if _QPC.key(face[1:4])==support]
            @test length(matches)==1
            face_center=[sum(p[k] for p in points)/4 for k in 1:3]
            xyz,uv,dim,owner=api.mesh.get_node(only(matches))
            @test (dim,owner)==(2,fixture.surface)
            @test maximum(abs.(xyz.-face_center))<=2e-11
            @test maximum(abs.(api.model.get_value(dim,owner,uv).-xyz))<=2e-11
            pivot=findfirst(==(complex.pivot),cell)
            opposite=cell[mod1(pivot+2,4)]
            endpoints=[ntuple(k->_QPC.point(source,v)[k]+(k==3 ? 1. : 0.),3)
                       for v in (complex.pivot,opposite)]
            midpoint=edge_mids[_QPC.key(positions[p] for p in endpoints)]
            xyz,uv,dim,owner=api.mesh.get_node(midpoint)
            center=[sum(p[k] for p in endpoints)/2 for k in 1:3]
            @test (dim,owner)==(2,top)
            @test maximum(abs.(xyz.-center))<=2e-11
            @test maximum(abs.(api.model.get_value(dim,owner,uv).-xyz))<=2e-11
            # Nonuniform native Curve middles distinguish the actual Quad9
            # source-face center from the pivot/opposite cap-edge midpoint.
            @test maximum(abs.(face_center[1:2].-center[1:2]))>.05
        end
    end)
end

@testset "QuadTriNoNewVerts actual 2-by-2 Quad patch" begin
    @testset "Existing immutable relations are unchanged" begin
        @test _QPC.packed_digest(Tessella.Model._extrude_nonew_templates())==
            "3d1c348a08a8645945a6864875ebec7d72c58ca419067279339d8aa0ee079ccc"
        @test _QPC.packed_digest(Tessella.Model._extrude_nonew_prism_templates())==
            "79082312c91d6de90190e7888769ff6f9a1ddf0335d5782e76b2650545d97556"
    end

    @testset "Normal products on all signed coordinate planes" begin
        for plane in (:XY,:YZ,:ZX),winding in (-1,1),direction in (-1,1),laterals in (false,true),
            (profile,shape) in ((:one,:unit),(:three,:skew),(:graded,:rounded))
            fixture=_QPC.fixture("product";plane,winding,height=direction*1.5,laterals,layers=profile,shape)
            geometry=_QPC.execute(fixture.source;dim=3)
            source=mesh_model_surface(geometry.model,fixture.surface)
            volume=geo_entity_mesh(geometry,3,Int(geometry.lists["sweep"][2]))
            cert=_QPC.certify(volume,fixture,source)
            @test cert.ncell==(laterals ? 4(length(fixture.levels)-2)+12 : 24(length(fixture.levels)-1)-4)
            @test Tessella.validate(volume).ok
            @test cert.total==cert.source.area*abs(_QPC.Q(last(cert.heights))-_QPC.Q(first(cert.heights)))
            @test nnodes(volume)==9length(fixture.levels)
        end
    end

    @testset "Source frame, explicit identities and sampled boundary chains" begin
        for winding in (-1,1),laterals in (false,true),reverse in (0,1,2,4,8,15)
            fixture=_QPC.fixture("frame";winding,laterals,curve_reverse=reverse,
                pins=reverse==15 ? :opposite : :auto,tags=reverse==8 ? :sparse : :dense,
                layers=:nonbinary,shape=:rounded)
            geometry=_QPC.execute(fixture.source;dim=3)
            source=mesh_model_surface(geometry.model,fixture.surface)
            cert=_QPC.certify(geo_entity_mesh(geometry,3,Int(geometry.lists["sweep"][2])),fixture,source)
            @test length(cert.source.chains)==4 && length(cert.source.shared)==4
            for curve in fixture.curve_tags
                @test length(geometry.model.curve_params[curve])==3
                @test issorted(geometry.model.curve_params[curve]) && all(diff(geometry.model.curve_params[curve]).>0)
            end
        end
    end

    @testset "Native source Curve progression keeps actual cap supports" begin
        for laterals in (false,true)
            _qpc_source_progression_case(laterals)
        end
    end

    @testset "Complete operation surfaces and lower-dimensional projection" begin
        for plane in (:XY,:YZ,:ZX),laterals in (false,true)
            fixture=_QPC.fixture("entry";plane,laterals,height=-1.,layers=:graded,shape=:rounded,
                                 winding=-1,pins=:opposite,tags=:sparse)
            geometry=_QPC.execute(fixture.source;dim=3)
            source=mesh_model_surface(geometry.model,fixture.surface)
            volume=geo_entity_mesh(geometry,3,Int(geometry.lists["sweep"][2]))
            cert=_QPC.certify(volume,fixture,source)
            _qpc_projection(geometry,fixture,source,volume,cert)
            surface_only=_QPC.execute(fixture.source;dim=2)
            @test _QPC.Quad.surface_faces(surface_only.mesh_parts,volume)==Set(keys(cert.boundary))
            standalone=mesh_model_volume(geometry.model,Int(geometry.lists["sweep"][2]))
            @test _QPC.Quad.typed_signature(standalone)==_QPC.Quad.typed_signature(volume)
            _qpc_api(fixture.source,(api,opened)->begin
                surface_product=api.mesh.generate(2)
                @test _QPC.Quad.surface_faces([(2,0,surface_product)],volume)==Set(keys(cert.boundary))
                actual=api.mesh.generate(3)
                @test _QPC.Quad.typed_signature(actual)==_QPC.Quad.typed_signature(volume)
                @test isempty(api.mesh.get_elements(2)[1]) # Existing API3 actual-cell cache contract.
                @test _QPC.certify(actual,fixture,source).total==cert.total
                payload=api.mesh.get_nodes(3,Int(opened.lists["sweep"][2]),true)
                @test length(payload[1])==nnodes(actual)
                before=copy(actual.coords);fill!(payload[2],NaN)
                @test actual.coords==before
            end)
        end
    end

    @testset "Actual precision rather than rounded surrogate centers" begin
        for plane in (:XY,:YZ,:ZX),direction in (-1,1),laterals in (false,true)
            axes=_QPC.FRAMES[plane]
            offset=ntuple(k->k==axes[3] ? 2.0^50+2. : 0.,3)
            fixture=_QPC.fixture("adjacent";plane,height=direction,laterals,layers=:four,offset)
            geometry=_QPC.execute(fixture.source;dim=3)
            source=mesh_model_surface(geometry.model,fixture.surface)
            volume=geo_entity_mesh(geometry,3,Int(geometry.lists["sweep"][2]))
            cert=_QPC.certify(volume,fixture,source)
            @test abs.(diff(cert.heights))==fill(.25,4)
            @test cert.total==1
            @test Tessella.validate(model_to_mixed(geometry.model,volume,3,Int(geometry.lists["sweep"][2]))).ok
        end
        for exponent in (-40,-48),laterals in (false,true)
            fixture=_QPC.fixture("tiny";scale=2.0^exponent,height=2.0^exponent,laterals,layers=:graded)
            geometry=_QPC.execute(fixture.source;dim=3)
            source=mesh_model_surface(geometry.model,fixture.surface)
            cert=_QPC.certify(geo_entity_mesh(geometry,3,Int(geometry.lists["sweep"][2])),fixture,source)
            @test cert.total==_QPC.Q(2.0^exponent)^3
        end
    end

    @testset "Physical metadata and independent coincident patch identity" begin
        for laterals in (false,true)
            f=_QPC.fixture("physical";laterals,layers=:one)
            text=f.source*"""
            Physical Surface("source",101)={1};Physical Surface("top",102)={sweep[0]};
            Physical Volume("body",103)={sweep[1]};
            """
            geometry=_QPC.execute(text;dim=3)
            volume=geo_entity_mesh(geometry,3,Int(geometry.lists["sweep"][2]))
            classified=model_to_mixed(geometry.model,volume,3,Int(geometry.lists["sweep"][2]))
            @test classified.physical_names==Dict((2,101)=>"source",(2,102)=>"top",(3,103)=>"body")
            @test all(all(==(103),b.tags) for b in classified.blocks if msh_dimension(b.msh)==3)
        end
        first=_QPC.fixture("left";laterals=false,layers=:three)
        second=_QPC.fixture("right";laterals=true,layers=:three,tags=:sparse,offset=(3.,0.,0.))
        pair=replace(first.source,"sweep[]="=>"left[]=")*replace(second.source,"sweep[]="=>"right[]=")
        _qpc_api(pair,(api,geometry)->begin
            # Public point edits preserve the two original CAD identities while
            # making their complete geometric products exactly coincident.
            for (tag,p) in collect(api.CURRENT[].points)
                p[1]>=3 || continue
                api.model.set_coordinates(tag,p[1]-3,p[2],p[3])
            end
            mesh=api.mesh.generate(3)
            left=Int(geometry.lists["left"][2]);right=Int(geometry.lists["right"][2])
            a=api.mesh.get_nodes(3,left,true)[1];b=api.mesh.get_nodes(3,right,true)[1]
            @test length(a)==36 && length(b)==36 && isempty(intersect(a,b))
            @test nnodes(mesh)==72
            @test length(Set(_QPC.point(mesh,i) for i in 1:nnodes(mesh)))==36
            for (entity,source,laterals) in ((left,first.surface,false),(right,second.surface,true))
                expected=_QPC.fixture("coincident";laterals,layers=:three)
                actual=_QPC.public_volume(api,entity)
                src=mesh_model_surface(api.CURRENT[],source)
                @test _QPC.certify(actual,expected,src).total==1
            end
            api.mesh.set_order(2)
            @test length(api.mesh.get_nodes(3,left,true)[1])==175
            @test length(api.mesh.get_nodes(3,right,true)[1])==175
            @test isempty(intersect(api.mesh.get_nodes(3,left,true)[1],api.mesh.get_nodes(3,right,true)[1]))
        end)
    end

    @testset "Independent certificates reject corrupted actual data" begin
        f=_QPC.fixture("corrupt";layers=:one)
        geometry=_QPC.execute(f.source;dim=0)
        src=mesh_model_surface(geometry.model,f.surface)
        pivot=_QPC.source_complex(src,f).pivot
        damaged=deepcopy(src);damaged.coords[:,pivot].=damaged.coords[:,1]
        @test_throws r"source coordinates coincide" _QPC.source_complex(damaged,f)
        damaged=deepcopy(src);damaged.blocks[1].nodes[1,1]=damaged.blocks[1].nodes[2,1]
        @test_throws r"foreign or repeated source node" _QPC.source_complex(damaged,f)
        volume=mesh_model_volume(geometry.model,Int(geometry.lists["sweep"][2]))
        damaged=deepcopy(volume);damaged.blocks[1].nodes[[1,2],1]=damaged.blocks[1].nodes[[2,1],1]
        @test_throws r"nonpositive whole P1 reference map" _QPC.certify(damaged,f,src)
    end

    @testset "Unsupported inputs reject before publishing operation state" begin
        for laterals in (false,true)
            fixture=_QPC.fixture("guard";laterals,layers=:three)
            invalid=(replace(fixture.source,"Transfinite Curve{1,2,3,4}=3"=>
                    "Transfinite Curve{1,3}=2; Transfinite Curve{2,4}=6"),
                replace(fixture.source,"Recombine Surface{1};"=>""),
                replace(fixture.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{0.25,0.0,1.0}"),
                replace(fixture.source,"Layers{3}"=>"Layers{{1},{0.5}}"),
                replace(fixture.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{{0,1,0},{-1,0,0},Pi/6}"),
                fixture.source*"Transfinite Surface{sweep[0]};\n")
            for text in invalid
                _qpc_atomic(text,r"QuadTriNoNewVerts")
            end
            _qpc_api(fixture.source,(api,geometry)->begin
                api.mesh.generate(3);api.mesh.set_order(2)
                for (position,curve) in pairs(fixture.curve_tags)
                    api.mesh.set_transfinite_curve(curve,position in (1,3) ? 2 : 6)
                end
                before=_qpc_api_snapshot(api)
                @test_throws r"QuadTriNoNewVerts" api.mesh.generate(3)
                _qpc_api_unchanged(api,before)
            end)
            collapsed=_QPC.fixture("collapsed";laterals,layers=:four,offset=(0.,0.,2.0^50))
            _qpc_atomic(replace(collapsed.source,"Layers{4}"=>"Layers{8}"),r"QuadTriNoNewVerts")
        end
    end
end
