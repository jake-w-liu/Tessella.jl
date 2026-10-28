# `plugin("Crack")`/`plugin("HomologyPostProcessing")`/`plugin("VoroMetal")`/
# `plugin("NearToFarField")` — Gmsh 4.15.2 registrations and executions.
# Crack: duplicated interior crack nodes + copied elements on a fresh
# discrete entity + positive-side rebinding. HomologyPostProcessing: chain
# transforms and boundary/trace/project operations. VoroMetal: truthful
# `requires mesh and Voro++` rejection. NearToFarField: complex E/H views
# to a far-field `_NearToFarField` view.
module PluginMeshPostChecks
using Test,Tessella
const A=Tessella.API

# Four-triangle square with a two-segment diagonal crack path.
function _crack_fixture()
    A.model.add("crack")
    A.model.add_discrete_entity(1,1)
    A.model.add_discrete_entity(2,2)
    A.mesh.add_nodes(2,2,[1,2,3,4,5],
        [0.,0.,0., 2.,0.,0., 2.,2.,0., 0.,2.,0., 1.,1.,0.])
    A.mesh.add_elements_by_type(2,2,[11,12,13,14],
        [1,2,5, 2,3,5, 3,4,5, 4,1,5])
    A.mesh.add_nodes(1,1,Int[],Float64[])
    A.mesh.add_elements_by_type(1,1,[21,22],[1,5, 5,3])
    A.model.add_physical_group(1,[1];tag=10)
    A.mesh.reclassify_nodes()
end

@testset "Plugin(Crack) option table and validation" begin
    A.initialize()
    try
        _crack_fixture()
        A.plugin.set_number("Crack","Dimension",1)
        A.plugin.set_number("Crack","PhysicalGroup",10)
        A.plugin.set_number("Crack","OpenBoundaryPhysicalGroup",0)
        A.plugin.set_number("Crack","NormalX",0)
        A.plugin.set_number("Crack","NormalY",0)
        A.plugin.set_number("Crack","NormalZ",1)
        A.plugin.set_number("Crack","NewPhysicalGroup",0)
        A.plugin.set_number("Crack","DebugView",0)
        A.plugin.set_number("Crack","SwapOrientation",0)
        @test_throws ArgumentError A.plugin.set_number(
            "Crack","NoSuchOption",1.)
        # Invalid dimension reports the upstream diagnostic, no mutation.
        A.option("General.AbortOnError",0)
        A.plugin.set_number("Crack","Dimension",3)
        A.plugin.run("Crack")
        @test A.logger.get_last_error()=="Crack dimension should be 1 or 2"
        # Empty physical group reports the upstream diagnostic.
        A.plugin.set_number("Crack","Dimension",1)
        A.plugin.set_number("Crack","PhysicalGroup",99)
        A.plugin.run("Crack")
        @test A.logger.get_last_error()==
            "Physical group 99 (dimension 1) is empty"
        # Under the default AbortOnError=2 the plugin Msg::Error throws and
        # gmsh::plugin::run wraps it as "Unknown plugin or plugin action".
        A.option("General.AbortOnError",2)
        A.plugin.set_number("Crack","Dimension",3)
        @test_throws ErrorException A.plugin.run("Crack")
        @test A.logger.get_last_error()=="Unknown plugin or plugin action"
    finally
        A.option("General.AbortOnError",2)
        A.finalize()
    end
end

@testset "Plugin(Crack) dim-1 surgery" begin
    A.initialize()
    try
        _crack_fixture()
        A.plugin.set_number("Crack","Dimension",1)
        A.plugin.set_number("Crack","PhysicalGroup",10)
        A.plugin.run("Crack")
        # A new discrete curve entity carries copies of both crack lines.
        curves=[t for (d,t) in A.model.get_entities(1)]
        newc=setdiff(curves,[1])
        @test length(newc)==1
        cn=only(newc)
        etypes,etags,enodes=A.mesh.get_elements(1,cn)
        @test etypes==[1]
        copies=only(enodes)
        @test length(copies)==4
        # Interior node 5 is duplicated to a fresh tag; endpoints 1,3 keep
        # the originals: copies are (1,dup) and (dup,3).
        @test copies[1]==1 && copies[4]==3
        @test copies[2]==copies[3] && copies[2]!=5 && copies[2]>5
        dup=copies[2]
        # Positive-side elements (the two tris on the +y side of the
        # crack tangent) rebind node 5 to the duplicate.
        _,_,tnodes=A.mesh.get_elements(2,2)
        tris=only(tnodes)
        @test Set([tris[7],tris[8]])==Set([3,4])
        @test tris[9]==dup
        @test Set([tris[10],tris[11]])==Set([4,1])
        @test tris[12]==dup
        @test tris[3]==5 && tris[6]==5  # negative side unchanged
        # The duplicated node sits at the original coordinates.
        ntags,ncoords=A.mesh.get_nodes(1,cn)
        @test dup in ntags
        # The new entity joins the crack physical group.
        @test cn in A.model.get_entities_for_physical_group(1,10)
    finally
        A.finalize()
    end
end

@testset "Plugin(Crack) NewPhysicalGroup + DebugView + SwapOrientation" begin
    A.initialize()
    try
        _crack_fixture()
        A.plugin.set_number("Crack","Dimension",1)
        A.plugin.set_number("Crack","PhysicalGroup",10)
        A.plugin.set_number("Crack","NewPhysicalGroup",55)
        A.plugin.set_number("Crack","SwapOrientation",1)
        A.plugin.set_number("Crack","DebugView",1)
        A.plugin.run("Crack")
        cn=only(setdiff([t for (d,t) in A.model.get_entities(1)],[1]))
        @test cn in A.model.get_entities_for_physical_group(1,55)
        _,_,enodes=A.mesh.get_elements(1,cn)
        copies=only(enodes)
        dup=copies[1]
        # Reversed copies: (5,1)->(dup,1) and (3,5)->(3,dup).
        @test copies[1]==dup && copies[2]==1
        @test copies[3]==3 && copies[4]==dup
        # DebugView: ElementNodeData flags on the positive-side elements.
        @test !isempty(A.view.get_tags())
    finally
        A.finalize()
    end
end

function _loop_fixture()
    A.model.add("loop")
    A.model.add_discrete_entity(1,1)
    A.model.add_discrete_entity(1,2)
    A.mesh.add_nodes(1,1,[11,12,13],
        [0.,0.,0., 1.,0.,0., 0.,1.,0.])
    A.mesh.add_elements_by_type(1,1,[101,102,103],
        [11,12,12,13,13,11])
    A.mesh.add_nodes(1,2,[14],[2.,0.,0.])
    A.mesh.add_elements_by_type(2,1,[104],[12,14])
    A.model.add_physical_group(1,[1];tag=10)
    A.model.add_physical_group(1,[2];tag=20)
    A.mesh.reclassify_nodes()
end

@testset "Plugin(HomologyPostProcessing) identity transform" begin
    A.initialize()
    try
        _loop_fixture()
        A.plugin.set_string("HomologyPostProcessing",
            "PhysicalGroupsOfOperatedChains","10")
        A.plugin.set_string("HomologyPostProcessing","TransformationMatrix","I")
        A.plugin.run("HomologyPostProcessing")
        groups=A.model.get_physical_groups()
        new=setdiff(groups,[(1,10),(1,20)])
        @test length(new)==1
        (dim,phys)=only(new)
        @test dim==1
        @test A.model.get_physical_name(dim,phys)=="C1 c1"
        ents=A.model.get_entities_for_physical_group(dim,phys)
        @test length(ents)==1
        etypes,_,enodes=A.mesh.get_elements(dim,only(ents))
        @test etypes==[1]
        @test length(only(enodes))==6  # three line cells of the loop chain
    finally
        A.finalize()
    end
end

@testset "Plugin(HomologyPostProcessing) boundary and trace" begin
    A.initialize()
    try
        _loop_fixture()
        # Boundary of the closed loop is the zero 0-chain: not added.
        A.plugin.set_string("HomologyPostProcessing",
            "PhysicalGroupsOfOperatedChains","10")
        A.plugin.set_string("HomologyPostProcessing","TransformationMatrix","I")
        A.plugin.set_number("HomologyPostProcessing",
            "ApplyBoundaryOperatorToResults",1)
        A.plugin.run("HomologyPostProcessing")
        @test isempty(setdiff(A.model.get_physical_groups(),[(1,10),(1,20)]))
        # Boundary of the open path (group 20) = its two endpoints.
        A.plugin.set_string("HomologyPostProcessing",
            "PhysicalGroupsOfOperatedChains","20")
        A.plugin.run("HomologyPostProcessing")
        new=setdiff(A.model.get_physical_groups(),[(1,10),(1,20)])
        @test length(new)==1
        (dim,phys)=only(new)
        @test dim==0
        @test A.model.get_physical_name(dim,phys)=="C0 c1"
        etypes,_,enodes=A.mesh.get_elements(dim,
            only(A.model.get_entities_for_physical_group(dim,phys)))
        @test etypes==[15]
        @test Set(only(enodes))==Set(UInt64[12,14])
    finally
        A.finalize()
    end
end

@testset "Plugin(VoroMetal) registration and truthful rejection" begin
    A.initialize()
    try
        _loop_fixture()
        A.option("General.AbortOnError",0)
        A.plugin.set_number("VoroMetal","ComputeBestSeeds",1)
        A.plugin.set_number("VoroMetal","ComputeMicrostructure",0)
        A.plugin.set_string("VoroMetal","SeedsFile","x.txt")
        @test_throws ArgumentError A.plugin.set_number(
            "VoroMetal","NoSuchOption",1.)
        A.plugin.run("VoroMetal")
        @test A.logger.get_last_error()==
            "Plugin(VoroMetal) requires mesh module and voro++"
    finally
        A.option("General.AbortOnError",2)
        A.finalize()
    end
end

@testset "Plugin(NearToFarField) complex views to far field" begin
    A.initialize()
    try
        _crack_fixture()
        # Two-step vector views on the triangle grid: step 0 = real,
        # step 1 = imaginary; E and H share the grid (same elements).
        ev=A.view.add("E")
        hv=A.view.add("H")
        # VT record: 3 nodes -> 9 coord floats, then 3 comps x 2 steps
        # (real, imag) -> 18 value floats.
        geo=[0.,0.,0., 1.,0.,0., 0.,1.,0.]
        efield=[1.,0.,0., 1.,0.,0., 1.,0.,0.,   # real part
                0.,0.,0., 0.,0.,0., 0.,0.,0.]   # imag part
        hfield=[0.,1.,0., 0.,1.,0., 0.,1.,0.,
                0.,0.,0., 0.,0.,0., 0.,0.,0.]
        A.view.add_list_data(ev,"VT",1,vcat(geo,efield))
        A.view.add_list_data(hv,"VT",1,vcat(geo,hfield))
        A.plugin.set_number("NearToFarField","EView",ev)
        A.plugin.set_number("NearToFarField","HView",hv)
        A.plugin.set_number("NearToFarField","NumPointsPhi",2)
        A.plugin.set_number("NearToFarField","NumPointsTheta",2)
        A.plugin.set_number("NearToFarField","Normalize",0)
        A.plugin.set_number("NearToFarField","dB",0)
        A.plugin.set_string("NearToFarField","MatlabOutputFile","")
        A.plugin.run("NearToFarField")
        tags=A.view.get_tags()
        @test length(tags)==3
        newtag=maximum(tags)
        dt,data=A.view.get_list_data(newtag)
        @test dt==["SQ"]
        flat=only(data)
        @test length(flat)==4*16
        # Oracle 4.15.2 (same fixture, positional views 0/1): farField at
        # theta=pi/2 is 0.001492077591486519, ~0 at the poles; SQ corners
        # carry f[i][j],f[i+1][j],f[i+1][j+1],f[i][j+1] and j iterates
        # innermost, so records alternate the two theta rows.
        for (rec,j) in ((0,0),(1,1),(2,0),(3,1))
            base=rec*16
            @test flat[base+13]≈(j==1 ? 0.001492077591486519 : 0.0) atol=1e-15
            @test flat[base+14]≈(j==1 ? 0.001492077591486519 : 0.0) atol=1e-15
            @test flat[base+15]≈(j==1 ? 0.0 : 0.001492077591486519) atol=1e-15
            @test flat[base+16]≈(j==1 ? 0.0 : 0.001492077591486519) atol=1e-15
        end
        # Incompatible views (mismatched step counts) report upstream's
        # ordering: counts first, complex check second.
        sv=A.view.add("S")
        A.view.add_list_data(sv,"VT",1,vcat(geo,efield[1:9]))
        A.plugin.set_number("NearToFarField","HView",sv)
        A.option("General.AbortOnError",0)
        A.plugin.run("NearToFarField")
        @test A.logger.get_last_error()==
            "Incompatible views for E field and H field"
    finally
        A.option("General.AbortOnError",2)
        A.finalize()
    end
end

end # module
