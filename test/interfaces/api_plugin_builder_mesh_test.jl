module PluginBuilderMeshTests
using Test,Tessella
const A=Tessella.API

function fixture(dir,cached)
    if cached
        path=joinpath(dir,"triangle.msh")
        write(path,raw"""$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
11 0 0 0
21 1 0 0
31 0 1 0
$EndNodes
$Elements
1
301 2 2 0 7 11 21 31
$EndElements
""")
        A.open(path)
    else
        A.model.add("builder")
        A.model.add_discrete_entity(2,7)
        A.mesh.add_nodes(2,7,[11,21,31],[0.,0.,0.,1.,0.,0.,0.,1.,0.])
        A.mesh.add_elements_by_type(7,2,[301],[11,21,31])
    end
    A.model.add_physical_group(2,[7];tag=23)
    return A.model.get_current()
end

function output(name)
    previous=A.view.get_tags()
    A.plugin.run(name)
    return only(setdiff(A.view.get_tags(),previous))
end

function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir(f)
    finally
        A.finalize()
    end
end

@testset "Model-read plugins resolve actual builder meshes, owners and sparse UIDs" begin
    # Native SDK 4.15.2 public controls on this discrete triangle yield an ST
    # for filtered MathEval, all3 nodes for NewView, area.5, the Gauss centroid
    # with UID301, and an ElementData Jacobian ratio1 on UID301.
    for cached in (false,true)
        @testset "cached=$cached" begin
            session() do dir
                model=fixture(dir,cached)
                source=A.view.add("source")
                A.view.add_model_data(source,0,model,"NodeData",[11,21,31],[[0.],[1.],[1.]])
                geometry=(A.mesh.get_nodes(),A.mesh.get_elements())
                cache=A.LAST_MESH[];classification=A.LAST_MESH_CLASS[]
                data=A.view.get_model_data(source,0)

                # The selected view keeps its own model when another model
                # is current; physical ownership must follow that source.
                A.model.add("foreign")
                A.plugin.set_number("MathEval","View",source)
                A.plugin.set_number("MathEval","PhysicalGroup",23)
                A.plugin.set_string("MathEval","Expression0","v0")
                types,records=A.view.get_list_data(output("MathEval"))
                @test types==["ST"]
                @test only(records)==[0.,1.,0.,0.,0.,1.,0.,0.,0.,0.,1.,1.]
                @test A.model.get_current()=="foreign"
                A.model.set_current(model)

                A.plugin.set_number("NewView","Value",7)
                A.plugin.set_number("NewView","NumComp",3)
                constant=output("NewView")
                @test A.view.get_model_data(constant,0)==(["NodeData"],UInt64[11,21,31],fill(7.,9))
                @test A.view.probe(constant,.25,.25,0.)[1]≈fill(7.,3)
                A.plugin.set_number("MeshVolume","Dimension",2)
                A.plugin.set_number("MeshVolume","PhysicalGroup",23)
                types,records=A.view.get_list_data(output("MeshVolume"))
                @test types==["SP"]
                @test only(records)≈[.5,.5,0.,.5]
                A.plugin.set_number("GaussPoints","PhysicalGroup",23)
                types,records=A.view.get_list_data(output("GaussPoints"))
                @test types==["SP"]
                @test only(records)≈[1/3,1/3,0.,301.]
                A.plugin.set_number("AnalyseMeshQuality","CreateView",1)
                A.plugin.set_number("AnalyseMeshQuality","JacobianDeterminant",1)
                quality=output("AnalyseMeshQuality")
                @test A.view.get_model_data(quality,0)==(["ElementData"],UInt64[301],[1.])
                @test A.view.probe(quality,.25,.25,0.)[1]≈[1.]
                @test (A.mesh.get_nodes(),A.mesh.get_elements())==geometry
                @test A.LAST_MESH[]===cache && A.LAST_MESH_CLASS[]===classification
                @test A.view.get_model_data(source,0)==data

                # Caller-owned records can coexist with an older cache on
                # the same entity. Read operators must see their union.
                A.mesh.add_nodes(2,7,[51,61,71],[2.,0.,0.,3.,0.,0.,2.,1.,0.])
                A.mesh.add_elements_by_type(7,2,[901],[51,61,71])
                geometry=(A.mesh.get_nodes(),A.mesh.get_elements())
                constant=output("NewView")
                dtype,ids,values=A.view.get_model_data(constant,0)
                @test dtype==["NodeData"] && Set(ids)==Set(UInt64[11,21,31,51,61,71])
                @test values==fill(7.,18)
                @test A.view.probe(constant,2.25,.25,0.)[1]≈fill(7.,3)
                types,records=A.view.get_list_data(output("MeshVolume"))
                @test types==["SP"]
                @test only(records)[4]≈1.
                quality=output("AnalyseMeshQuality")
                dtype,ids,values=A.view.get_model_data(quality,0)
                @test dtype==["ElementData"] && Set(ids)==Set(UInt64[301,901])
                @test values==[1.,1.]
                @test (A.mesh.get_nodes(),A.mesh.get_elements())==geometry
            end
        end
    end
end

@testset "Model interpolation keeps captured arrays unboxed" begin
    boxed=Symbol[]
    inspected=0
    for name in names(A;all=true)
        replace(String(name),r"^#|#\d+$"=>"")=="_plugin_probe_model" || continue
        isdefined(A,name) || continue
        object=getfield(A,name)
        object isa Function || continue
        for method in methods(object)
            method.module===A || continue
            inspected+=1
            info=Base.uncompressed_ast(method)
            for statement in info.code
                statement isa Expr && statement.head==:(=) || continue
                rhs=statement.args[2]
                (rhs isa Expr && rhs.head==:call && rhs.args[1]==GlobalRef(Core,:Box)) || continue
                push!(boxed,info.slotnames[statement.args[1].id])
            end
        end
    end
    @test inspected>0
    @test isempty(boxed)
end
end
