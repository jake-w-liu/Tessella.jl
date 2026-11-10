module PluginModelSmoothTests
using Test,Tessella
const A=Tessella.API

function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        f()
    finally
        A.finalize()
    end
end

function tri3()
    A.model.add("smooth")
    A.model.add_discrete_entity(2,1)
    A.mesh.add_nodes(2,1,[1,2,3,4],
        [0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.])
    A.mesh.add_elements_by_type(1,2,[11,12],[1,2,3,2,4,3])
end

function model_rows(kind,nc)
    if kind=="ElementData"
        return [[Float64(10i+j) for j in 1:nc] for i in 0:1]
    end
    return [[Float64(10i+j) for j in 1:3nc] for i in 0:1]
end

function expected_rows(kind,nc,rows)
    if kind=="ElementData"
        first=rows[1];second=rows[2]
        shared=[(first[c]+second[c])/2 for c in 1:nc]
        return [first;shared;shared;second]
    end
    first=rows[1];second=rows[2]
    # Second triangle connectivity is (node2,node4,node3), so its row
    # positions 1,2,3 carry nodes 2,4,3 respectively.
    node2=[(first[nc+c]+second[c])/2 for c in 1:nc]
    node3=[(first[2nc+c]+second[2nc+c])/2 for c in 1:nc]
    return [first[1:nc];node2;node3;second[nc+1:2nc]]
end

@testset "Model Smooth averages incident values by UID" begin
    for kind in ("ElementData","ElementNodeData"),nc in (1,2)
        session() do
            tri3()
            rows=model_rows(kind,nc)
            view=A.view.add("v")
            A.view.add_model_data(view,0,"smooth",kind,[11,12],rows,0.0,nc)
            A.plugin.set_number("Smooth","View",view)
            A.plugin.run("Smooth")
            types,tags,values=A.view.get_model_data(view,0)
            @test types==["NodeData"]
            @test tags==UInt64[1,2,3,4]
            @test values≈reduce(vcat,expected_rows(kind,nc,rows)) atol=1e-15
            context=A.SESSION_VIEWS[view].model_contexts[(0,"NodeData")]
            @test context.counts==fill(Int32(nc),4)
            @test context.extent>=104
        end
    end
end

@testset "Model Smooth visits complete higher-order nodes" begin
    session() do
        A.model.add("high order")
        A.model.add_discrete_entity(2,1)
        A.mesh.add_nodes(2,1,[1,2,3,4,5,6],
            [0.,0.,0.,1.,0.,0.,0.,1.,0.,.5,0.,0.,.5,.5,0.,0.,.5,0.])
        A.mesh.add_elements_by_type(1,9,[11],[1,2,3,4,5,6])
        view=A.view.add("v")
        A.view.add_model_data(view,0,"high order","ElementData",[11],[[7.]],0.0,1)
        A.plugin.set_number("Smooth","View",view)
        A.plugin.run("Smooth")
        types,tags,values=A.view.get_model_data(view,0)
        @test types==["NodeData"] && tags==UInt64.(1:6)
        @test values==fill(7.,6)
        A.view.add_model_data(view,1,"high order","ElementNodeData",[11],
            [collect(1.:6.)],0.0,1)
        A.plugin.set_number("Smooth","View",view)
        A.plugin.run("Smooth")
        types,tags,values=A.view.get_model_data(view,1)
        @test types==["NodeData"] && tags==UInt64.(1:6)
        @test values==collect(1.:6.)
    end
end

@testset "Coincident distinct UIDs remain distinct" begin
    session() do
        A.model.add("coincident")
        A.model.add_discrete_entity(2,1)
        A.mesh.add_nodes(2,1,[1,2,3,4,5,6],
            [0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,1.,0.])
        A.mesh.add_elements_by_type(1,2,[11,12],[1,2,3,4,5,6])
        view=A.view.add("v")
        A.view.add_model_data(view,0,"coincident","ElementData",[11,12],
            [[1.],[3.]],0.0,1)
        A.plugin.set_number("Smooth","View",view)
        A.plugin.run("Smooth")
        types,tags,values=A.view.get_model_data(view,0)
        @test types==["NodeData"] && tags==UInt64.(1:6)
        @test values==[1.,1.,1.,3.,3.,3.]
    end
end

@testset "Short ElementNodeData rows do not fabricate zero samples" begin
    session() do
        A.model.add("short")
        A.model.add_discrete_entity(2,1)
        A.mesh.add_nodes(2,1,[1,2,3,4],
            [0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.])
        A.mesh.add_elements_by_type(1,2,[11,12],[1,2,3,4,2,3])
        view=A.view.add("v")
        A.view.add_model_data(view,0,"short","ElementNodeData",[11,12],
            [[10.,20.,30.],[40.]],0.0,1)
        A.plugin.set_number("Smooth","View",view)
        A.plugin.run("Smooth")
        types,tags,values=A.view.get_model_data(view,0)
        @test types==["NodeData"] && tags==UInt64[1,2,3,4]
        @test values==[10.,20.,30.,40.]
    end
end

@testset "Partial ElementNodeData components are rejected atomically" begin
    session() do
        A.model.add("partial")
        A.model.add_discrete_entity(2,1)
        A.mesh.add_nodes(2,1,[1,2,3],
            [0.,0.,0.,1.,0.,0.,0.,1.,0.])
        A.mesh.add_elements_by_type(1,2,[11],[1,2,3])
        view=A.view.add("v")
        A.view.add_model_data(view,0,"partial","ElementNodeData",[11],
            [[1.,2.,3.]],0.0,2)
        before=A.view.get_model_data(view,0)
        A.plugin.set_number("Smooth","View",view)
        @test_throws ArgumentError A.plugin.run("Smooth")
        @test A.view.get_model_data(view,0)==before
    end
end

@testset "Malformed ElementData rows are rejected atomically" begin
    session() do
        A.model.add("partial element")
        A.model.add_discrete_entity(2,1)
        A.mesh.add_nodes(2,1,[1,2,3],
            [0.,0.,0.,1.,0.,0.,0.,1.,0.])
        A.mesh.add_elements_by_type(1,2,[11],[1,2,3])
        view=A.view.add("v")
        A.view.add_model_data(view,0,"partial element","ElementData",[11],
            [[1.,2.,3.]],0.0,2)
        before=A.view.get_model_data(view,0)
        A.plugin.set_number("Smooth","View",view)
        @test_throws ArgumentError A.plugin.run("Smooth")
        @test A.view.get_model_data(view,0)==before
    end
end

@testset "Steps, times, partitions and multiple models survive projection" begin
    session() do
        tri3()
        view=A.view.add("v")
        A.view.add_model_data(view,0,"smooth","ElementData",[11,12],
            [[1.],[3.]],.25,1,2)
        A.view.add_model_data(view,5,"smooth","ElementData",[11,12],
            [[11.],[13.]],1.25,1,7)
        A.model.add("second")
        A.model.add_discrete_entity(2,1)
        A.mesh.add_nodes(2,1,[101,102,103],[0.,0.,0.,2.,0.,0.,0.,2.,0.])
        A.mesh.add_elements_by_type(1,2,[201],[101,102,103])
        A.view.add_model_data(view,9,"second","ElementData",[201],
            [[21.]],2.25,1,3)
        before=(A.mesh.get_nodes(),A.mesh.get_elements(),A.model.get_current())
        A.plugin.set_number("Smooth","View",view)
        A.plugin.run("Smooth")
        @test A.view.get_model_data(view,0)==(["NodeData"],UInt64[1,2,3,4],
            [1.,2.,2.,3.])
        @test A.view.get_model_data(view,5)==(["NodeData"],UInt64[1,2,3,4],
            [11.,12.,12.,13.])
        @test A.view.get_model_data(view,9)==(["NodeData"],UInt64[101,102,103],
            fill(21.,3))
        owner=A.SESSION_VIEWS[view]
        @test owner.model_contexts[(0,"NodeData")].time==.25
        @test owner.model_contexts[(5,"NodeData")].time==1.25
        @test owner.model_contexts[(9,"NodeData")].time==2.25
        @test owner.model_contexts[(0,"NodeData")].partitions==Set([2])
        @test owner.model_contexts[(5,"NodeData")].partitions==Set([7])
        @test owner.model_contexts[(9,"NodeData")].partitions==Set([3])
        @test (A.mesh.get_nodes(),A.mesh.get_elements(),A.model.get_current())==before
    end
end

@testset "Unknown rows, sparse UIDs and repeated runs are bounded" begin
    session() do
        A.model.add("sparse")
        A.model.add_discrete_entity(2,1)
        A.mesh.add_nodes(2,1,[10001,10002,10003,4],
            [0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.])
        A.mesh.add_elements_by_type(1,2,[11,12],[10001,10002,10003,10002,4,10003])
        view=A.view.add("v")
        A.view.add_model_data(view,0,"sparse","ElementData",[11,12,999],
            [[1.],[3.],[99.]],0.0,1)
        A.plugin.set_number("Smooth","View",view)
        A.plugin.run("Smooth")
        first=A.view.get_model_data(view,0)
        @test first==(["NodeData"],UInt64[4,10001,10002,10003],
            [3.,1.,2.,2.])
        context=A.SESSION_VIEWS[view].model_contexts[(0,"NodeData")]
        # Sparse extent grows by the bounded native increment, not to an
        # allocation proportional to the largest UID.
        @test context.extent==UInt64(10101)
        A.plugin.run("Smooth")
        @test A.view.get_model_data(view,0)==first
        @test A.SESSION_VIEWS[view].model_kind=="NodeData"
    end
end

@testset "Smoothed NodeData probes and exports" begin
    session() do
        tri3()
        view=A.view.add("v")
        A.view.add_model_data(view,0,"smooth","ElementData",[11,12],
            [[1.],[3.]],0.0,1)
        A.plugin.set_number("Smooth","View",view)
        A.plugin.run("Smooth")
        @test only(first(A.view.probe(view,1/3,1/3,0.)))≈5/3 rtol=1e-14
        @test only(first(A.view.probe(view,2/3,2/3,0.)))≈7/3 rtol=1e-14
        mktempdir() do directory
            path=joinpath(directory,"smooth.pos")
            A.view.write(view,path)
            @test isfile(path) && filesize(path)>0
            text=read(path,String)
            @test occursin("ST(",text)
            @test occursin("{1,2,2}",text)
            @test occursin("{2,3,2}",text)
        end
    end
end

@testset "Projection preserves aliases and rejects unknown models atomically" begin
    session() do
        tri3()
        view=A.view.add("v")
        A.view.add_model_data(view,0,"smooth","ElementData",[11,12],
            [[1.],[3.]],0.0,1)
        alias=A.view.add_alias(view)
        A.plugin.set_number("Smooth","View",alias)
        A.plugin.run("Smooth")
        expected=(["NodeData"],UInt64[1,2,3,4],[1.,2.,2.,3.])
        @test A.view.get_model_data(view,0)==expected
        @test A.view.get_model_data(alias,0)==expected
        before=A.view.get_model_data(view,0)
        @test_throws ArgumentError A.view.add_model_data(view,1,"missing",
            "ElementData",[11],[[5.]],0.0,1)
        @test A.view.get_model_data(view,0)==before
    end
end
end
