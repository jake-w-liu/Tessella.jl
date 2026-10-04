using Test,Tessella
using Tessella.Elements:MixedMesh,ElementBlock,msh_dimension
using Tessella.MeshTypes:Mesh,nnodes
using Tessella.Model:mesh_model_surface,mesh_model_volume,model_to_mixed
if !isdefined(@__MODULE__,:QuadTriNoNewQuadStripCertificates)
    include("quadtri_nonew_quad_strip_certificates.jl")
end
const _QSC=QuadTriNoNewQuadStripCertificates

function _qsc_product(f;entry=false)
    geometry=_QSC.execute(f.source;dim=entry ? 3 : 0)
    out=Int.(geometry.lists["sweep"]);m=geometry.model
    source=mesh_model_surface(m,f.surface)
    volume=entry ? geo_entity_mesh(geometry,3,out[2]) : mesh_model_volume(m,out[2])
    cert=_QSC.certify(volume,f,source)
    @test validate(volume).ok
    base=source.coords[f.axes[3],1]
    @test cert.total==sum(cert.source.areas)*abs(_QSC.Q(base+f.height)-_QSC.Q(base))
    @test length(cert.source.chains)==4 && sort!(length.(collect(cert.source.chains)))==[2,2,3,3]
    @test length(cert.centers)==(f.laterals ? 2 : 0)
    for curve in f.curve_tags
        params=m.curve_params[curve]
        @test params[1]==0 && params[end]==1 && all(diff(params).>0)
    end
    projected=model_to_mixed(m,volume,3,out[2])
    @test validate(projected).ok
    @test _QSC.Quad.typed_signature(projected)==_QSC.Quad.typed_signature(volume)
    @test count(o->o[1]==0,projected.entity_data.node_entities)==8
    @test count(o->o[1]==1,projected.entity_data.node_entities)==4(length(f.levels)-1)
    @test count(o->o[1]==2,projected.entity_data.node_entities)==2(length(f.levels)-2)
    @test count(o->o[1]==3,projected.entity_data.node_entities)==(f.laterals ? 2 : 0)
    parts=[(2,s,mesh_model_surface(m,s)) for s in (f.surface,out[1],out[3:end]...)]
    @test _QSC.Quad.surface_faces(parts,volume)==Set(keys(cert.boundary))
    @test _QSC.Quad.surface_faces(_QSC.execute(f.source;dim=2).mesh_parts,volume)==Set(keys(cert.boundary))
    @test _QSC.Quad.typed_signature(mesh_model_volume(m,out[2]))==_QSC.Quad.typed_signature(volume)
    return (;geometry,source,volume,cert,projected)
end

function _qsc_atomic(f,diagnostic)
    g=_QSC.execute(f.source;dim=0);m=g.model;out=Int.(g.lists["sweep"])
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

@testset "NoNew two actual quadrangles form a conforming strip" begin
    @testset "Actual columns, reference maps, signed shells and CAD carriers" begin
        for plane in (:XY,:YZ,:ZX),height in (-1.,1.),laterals in (false,true),layers in (:one,:three,:graded)
            _qsc_product(_QSC.fixture("strip";plane,height,laterals,layers))
        end
        for plane in (:XY,:YZ,:ZX),laterals in (false,true),shape in (:skew,:rounded)
            f=_QSC.fixture("frame";plane,laterals,shape,height=-.75,layers=:graded,
                winding=-1,pins=:opposite,curve_reverse=5,tags=:sparse,
                offset=(.2,-.3,.4),progression=4.)
            _qsc_product(f;entry=true)
        end
    end

    @testset "Independent certificate rejects corrupted actual complexes" begin
        f=_QSC.fixture("certificate";layers=:one)
        p=_qsc_product(f)
        cells=copy(only(_QSC.blocks(p.source)).nodes);cells[:,2]=reverse(cells[:,2])
        bad=MixedMesh(p.source.coords,[ElementBlock(3,cells)])
        @test_throws r"CAD winding|equal cycles" _QSC.source_complex(bad,f)
        coords=copy(p.volume.coords);coords[3,1]=.125
        @test_throws Exception _QSC.certify(MixedMesh(coords,_QSC.blocks(p.volume)),f,p.source)
        @test _QSC.crosses((0,0),(1,1),(0,1),(1,0))
        @test !_QSC.crosses((0,0),(1,0),(0,1),(1,1))
    end

    @testset "Represented precision and actual emitted centroids" begin
        for laterals in (false,true),height in (-1.,1.)
            f=_QSC.fixture("ulp_plane";height,laterals,layers=:one)
            text=replace(f.source,"Layers{1}"=>"Layers{4}")
            g=_QSC.execute(text;dim=0);offset=2.0^50+2
            for (tag,p) in collect(g.model.points)
                g.model.points[tag]=(p[1],p[2],p[3]+offset)
            end
            actual=merge(f,(;source=text,corners=Tuple((p[1],p[2],p[3]+offset) for p in f.corners),levels=collect(0.:.25:1.)))
            src=mesh_model_surface(g.model,f.surface);before=deepcopy(g.model.curve_params)
            if laterals
                @test_throws r"centroid|degenerate|folded|positivity" mesh_model_volume(g.model,Int(g.lists["sweep"][2]))
                @test g.model.curve_params==before
            else
                volume=mesh_model_volume(g.model,Int(g.lists["sweep"][2]))
                @test _QSC.certify(volume,actual,src).total==1
                @test all(abs.(diff(sort!(unique(volume.coords[3,:])))).==.25)
            end
        end
        for laterals in (false,true)
            f=_QSC.fixture("inplane_ulp";laterals,layers=:one)
            g=_QSC.execute(f.source;dim=0);offset=2.0^51
            for (tag,p) in collect(g.model.points);g.model.points[tag]=(p[1]+offset,p[2],p[3]);end
            actual=merge(f,(;corners=Tuple((p[1]+offset,p[2],p[3]) for p in f.corners)))
            src=mesh_model_surface(g.model,f.surface)
            @test nnodes(src)==6
            if laterals
                before=repr(g.model)
                @test_throws r"centroid|degenerate|folded|positivity" mesh_model_volume(g.model,Int(g.lists["sweep"][2]))
                @test repr(g.model)==before
            else
                @test _QSC.certify(mesh_model_volume(g.model,Int(g.lists["sweep"][2])),actual,src).total==1
            end
        end
    end

    @testset "Unsupported source or transform failures are atomic" begin
        f=_QSC.fixture("guard";layers=:one)
        larger=replace(f.source,"}=3;"=>"}=6;")
        _qsc_atomic(merge(f,(;source=larger)),r"QuadTriNoNewVerts")
        tilted=replace(f.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{0.1,0.0,1.0}")
        @test tilted!=f.source
        _qsc_atomic(merge(f,(;source=tilted)),r"normal|translation")
        cap=f.source*"Transfinite Surface{sweep[0]};"
        _qsc_atomic(merge(f,(;source=cap)),r"constrained-boundary|override")
    end
end

@testset "Native straight transfinite density controls" begin
    # Actual Gmsh 4.15.2 meshGEdge F_Transfinite samples, not the independent
    # analytic grading function. HWall progression retains its established
    # approximate integration bound; all other controls use 2e-11.
    laws=(("Progression",4.,3,[.2000000039104483]),
        ("Progression",.25,3,[.8000000128519222]),
        ("Progression",1e6,3,[1.0074495582442956e-6]),
        ("Progression",1e15,3,[7.450581596923805e-9]),
        ("Progression",1e-6,3,[.9999999925474284]),
        ("Progression",1e-15,3,[.499999996273405]),
        ("Bump",.2,7,[.08333334191164227,.249999996538392,.49999999999795225,.7500000034595244,.9166666580878113]),
        ("Bump",2.,7,[.2113246604743475,.3660253944175726,.4999999999989506,.6339746055809568,.788675339524792]),
        ("Beta",1.1,7,[.06125511290227867,.15546752144623055,.2940833678182301,.4852738556462574,.7267565949164975]),
        ("Progression_HWall",.1,7,[.10000000039746743,.2202793661449751,.3649506245467072,.5389602948057355,.748258021875258]),
        ("Bump_HWall",.01,7,[.014064499421960246,.11176990985389142,.49999999999625205,.8882300901442276,.9859355005778916]),
        ("Beta_HWall",.05,7,[.06037197207705748,.15369595794395216,.2916878498704612,.482881864770212,.7252306565317737]))
    for (law,coef,count,samples) in laws
        text="Point(1)={0,0,0,1};Point(2)={1,0,0,1};Line(1)={1,2};Transfinite Curve{1}=$count Using $law $coef;"
        m=_QSC.execute(text;dim=0).model;s=m.meshing.transfinite_curves[1]
        actual=Tessella.Model._model_curve_transfinite_native_params(m,1,0.,1.,s,"saved native control")
        @test actual[1]==0 && actual[end]==1 && all(diff(actual).>0)
        @test actual[2:end-1]≈samples atol=(law=="Progression_HWall" ? 3e-7 : 2e-11) rtol=0
        m.meshing.transfinite_curves[1]=merge(s,(;num_nodes=2))
        @test Tessella.Model._model_curve_transfinite_native_params(m,1,0.,1.,m.meshing.transfinite_curves[1],"endpoints")==[0.,1.]
    end
    m=_QSC.execute("Point(1)={0,0,0,1};Point(2)={1,0,0,1};Line(1)={1,2};Transfinite Curve{1}=3;";dim=0).model
    m.points[2]=m.points[1]
    @test Tessella.Model._model_curve_transfinite_native_params(m,1,0.,1.,m.meshing.transfinite_curves[1],"zero uniform")==[0.,.5,1.]
    @test isempty(Tessella.Model._model_curve_mesh_parts(m,"zero public curve"))
    val=Tessella.Model._transfinite_val(1,1e-15,1.,3)
    @test val(1.,1.)==0
    @test_throws r"unrepresentable interior" val(1.1,1.)
    # Explicitly seeded exact tiny fractions preserve the independent local
    # endpoint-snap contract, which native extreme grading does not guarantee.
    m.points[2]=(1.,0.,0.);tiny=1/(1+1e15)
    mesh=Mesh([0. tiny 1.;0. 0. 0.;0. 0. 0.];segs=Int32[1 2;2 3])
    entries,tolerance=Tessella.Model._curve_parameter_nodes(m,mesh,1,trues(3),Set([(Int32(1),Int32(2)),(Int32(2),Int32(3))]),1e-12,"seeded exact chain")
    @test first.(entries)==[0.,tiny,1.]
    @test tolerance<tiny/2
end
