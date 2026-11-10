module PluginDistancePDETests
using Test,Tessella
const A=Tessella.API

function fixture(kind;order=1,scale=1.)
    dim=kind===:line ? 1 : kind===:tet ? 3 : 2
    owner=dim==3 ? 8 : 7
    if kind===:line
        points=[(0.,0.,0.),(.5,0.,0.),(1.,0.,0.)]
        boundary=[[11],[31]];cells=[[11,21],[21,31]]
        boundary_type=15;cell_types=[1,1]
    elseif kind===:tet
        points=[(0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),(.25,.25,.25)]
        boundary=[[11,31,21],[11,21,41],[21,31,41],[31,11,41]]
        cells=[[11,21,31,51],[11,41,21,51],[21,41,31,51],[31,41,11,51]]
        boundary_type=2;cell_types=fill(4,4)
    elseif kind in (:quad,:mixed)
        points=[(i/2,j/2,0.) for j in 0:2 for i in 0:2]
        boundary=[[11,21],[21,31],[31,61],[61,91],[91,81],[81,71],[71,41],[41,11]]
        cells=[[11,21,51,41],[21,31,61,51],[41,51,81,71],[51,61,91,81]]
        boundary_type=1;cell_types=fill(3,4)
        if kind===:mixed
            cells=[cells[1:3];[[51,61,91],[51,91,81]]]
            cell_types=[3,3,3,2,2]
        end
    else
        points=[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),(.5,.5,0.)]
        boundary=[[11,21],[21,31],[31,41],[41,11]]
        cells=[[11,21,51],[21,31,51],[31,41,51],[41,11,51]]
        boundary_type=1;cell_types=fill(2,4)
        if order==2
            append!(points,[(.5,0.,0.),(1.,.5,0.),(.5,1.,0.),(0.,.5,0.),
                (.25,.25,0.),(.75,.25,0.),(.75,.75,0.),(.25,.75,0.)])
            boundary=[[11,21,61],[21,31,71],[31,41,81],[41,11,91]]
            cells=[[11,21,51,61,111,101],[21,31,51,71,121,111],
                [31,41,51,81,131,121],[41,11,51,91,101,131]]
            boundary_type=8;cell_types=fill(9,4)
        end
    end
    A.model.add("diffusion")
    for tag in eachindex(boundary);A.model.add_discrete_entity(dim-1,tag);end
    A.model.add_discrete_entity(dim,owner)
    ids=[10i+1 for i in 1:length(points)]
    A.mesh.add_nodes(dim,owner,ids,scale.*collect(Iterators.flatten(points)))
    for tag in eachindex(boundary)
        A.mesh.add_elements_by_type(tag,boundary_type,[200+tag],boundary[tag])
    end
    for i in eachindex(cells)
        A.mesh.add_elements_by_type(owner,cell_types[i],[300+i],cells[i])
    end
    return ids
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

function output(dtype=1.)
    A.plugin.set_number("Distance","DistanceType",dtype)
    previous=A.view.get_tags();A.plugin.run("Distance")
    return only(setdiff(A.view.get_tags(),previous))
end

function samples(tag)
    types,data=A.view.get_list_data(tag)
    values=Dict{NTuple{3,Float64},Float64}()
    for (type,flat) in zip(types,data)
        n=Dict("SL"=>2,"ST"=>3,"SQ"=>4,"SS"=>4)[type]
        for offset in 0:4n:length(flat)-1,node in 1:n
            point=ntuple(axis->flat[offset+(axis-1)*n+node],3)
            value=flat[offset+3n+node]
            haskey(values,point) && (@test values[point]≈value rtol=1e-11 atol=0.)
            values[point]=value
        end
    end
    return types,values
end

@testset "PDE Distance uses complete nodal FEM on actual model products" begin
    # SDK4.15.2 independent opaque controls; straight fan P1 also has the
    # analytic center s=2/(48mu^2+1), mu=DistanceType in bbox-diagonal units.
    for order in (1,2),reopen in (false,true)
        session() do dir
            ids=fixture(:tri;order)
            if reopen
                path=joinpath(dir,"pde.msh");A.write(path);A.clear();A.open(path)
            end
            before=(A.mesh.get_nodes(),A.mesh.get_elements())
            reference=order==1 ?
                (.02946304896348487,.058934092430339834,.11812507609023656,.24506453586713703,1.3025388268121334) :
                (.022045422176726236,.043783809989140714,.08519397586197644,.15345577644718975,.1851294667245662)
            for (dtype,expected) in zip((2.,1.,.5,.25,.1),reference)
                tag=output(dtype);types,values=samples(tag)
                @test types==["ST"]
                @test values[(.5,.5,0.)]≈expected rtol=2e-11
                @test all(values[point]==0. for point in ((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)))
                @test A.view.get_number(tag,"NbTimeStep")==1.
            end
            @test (A.mesh.get_nodes(),A.mesh.get_elements())==before
            @test length(ids)==(order==1 ? 5 : 13)
        end
    end
end

@testset "Diffusion Distance supports dimensions and mixed families" begin
    # Independent SDK values at DistanceType1 for the single interior node.
    for (kind,point,expected) in ((:line,(.5,0.,0.),.12260232209233239),
            (:quad,(.5,.5,0.),.06647658264152548),
            (:mixed,(.5,.5,0.),.06395794695823259),
            (:tet,(.25,.25,.25),.006016160295109585))
        session() do _
            fixture(kind);before=(A.mesh.get_nodes(),A.mesh.get_elements())
            types,values=samples(output())
            @test !isempty(types)
            @test values[point]≈expected rtol=2e-11
            @test all(value==0. for (p,value) in values if p!=point)
            @test (A.mesh.get_nodes(),A.mesh.get_elements())==before
        end
    end
end

@testset "Curved PDE geometry uses its complete basis and separate integration rules" begin
    session() do _
        ids=fixture(:tri;order=2)
        A.mesh.set_node(51,[.5,.5,.1])
        for (uid,x,y) in ((101,.25,.25),(111,.75,.25),(121,.75,.75),(131,.25,.75))
            A.mesh.set_node(uid,[x,y,.02])
        end
        source=A.view.add("retained curved values")
        A.view.add_model_data(source,0,A.model.get_current(),"NodeData",ids,[[Float64(i)] for i in eachindex(ids)])
        before=(A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(source,0))
        # Independently sampled native SDK4.15.2: all source coordinates were
        # installed before the elements, so no stale geometry cache is involved.
        for (dtype,expected) in zip((2.,1.,.5,.25,.1),
                (.023175447269645176,.04601944282631918,.08947955979389964,
                 .16075569527082229,.18906048602197745))
            _,values=samples(output(dtype))
            @test values[(.5,.5,.1)]≈expected rtol=2e-11
        end
        @test (A.mesh.get_nodes(),A.mesh.get_elements(),A.view.get_model_data(source,0))==before
    end
end

@testset "Selected PDE boundary uses natural data elsewhere and keeps rescaling" begin
    session() do _
        fixture(:tri);A.model.add_physical_group(1,[1];tag=13)
        A.plugin.set_number("Distance","PhysicalLine",13.)
        _,values=samples(output())
        @test values[(0.,0.,0.)]==values[(1.,0.,0.)]==0.
        @test values[(.5,.5,0.)]≈.2860084318937513 rtol=2e-11
        @test values[(0.,1.,0.)]≈.4570229177730143 rtol=2e-11
        @test values[(1.,1.,0.)]≈.4570229177730143 rtol=2e-11
        A.plugin.set_number("Distance","MinScale",1.)
        A.plugin.set_number("Distance","MaxScale",3.)
        _,scaled=samples(output())
        @test collect(extrema(Base.values(scaled)))≈[1.,3.]
        A.plugin.set_number("Distance","PhysicalLine",999.)
        @test isempty(first(A.view.get_list_data(output())))
    end
    for point_count in (1,2)
        session() do _
            fixture(:tri)
            for tag in 1:point_count
                A.model.add_discrete_entity(0,tag)
                A.mesh.add_elements_by_type(tag,15,[100+tag],[10tag+1])
            end
            A.model.add_physical_group(0,collect(1:point_count);tag=17)
            A.plugin.set_number("Distance","PhysicalPoint",17.)
            types,values=samples(output())
            @test types==["ST"]
            if point_count==1
                @test all(iszero,Base.values(values))
            else
                @test values[(.5,.5,0.)]≈.2860084318937513 rtol=2e-11
                @test values[(0.,1.,0.)]≈.4570229177730143 rtol=2e-11
            end
        end
    end
end

@testset "Disconnected natural components retain the constant reaction solution" begin
    session() do _
        fixture(:tri)
        for tag in 21:24;A.model.add_discrete_entity(1,tag);end
        A.model.add_discrete_entity(2,17)
        ids=[1011,1021,1031,1041,1051]
        points=[(3.,0.,0.),(4.,0.,0.),(4.,1.,0.),(3.,1.,0.),(3.5,.5,0.)]
        A.mesh.add_nodes(2,17,ids,collect(Iterators.flatten(points)))
        for (tag,nodes) in zip(21:24,([1011,1021],[1021,1031],[1031,1041],[1041,1011]))
            A.mesh.add_elements_by_type(tag,1,[2000+tag],nodes)
        end
        A.mesh.add_elements_by_type(17,2,[1301,1302,1303,1304],
            [1011,1021,1051,1021,1031,1051,1031,1041,1051,1041,1011,1051])
        A.model.add_physical_group(1,[1,2,3,4];tag=13)
        A.plugin.set_number("Distance","PhysicalLine",13.)
        for dtype in (1.,1e200)
            _,values=samples(output(dtype))
            expected=-log1p(-.9999)*dtype*sqrt(2.)
            @test all(values[point]≈expected for point in points)
            @test values[(0.,0.,0.)]==0.
            @test values[(.5,.5,0.)]>0.
        end
    end
end

@testset "PDE units and tiny solutions retain finite distances" begin
    for scale in (1e-100,1.,1e100)
        session() do _
            fixture(:tri;scale)
            _,values=samples(output())
            @test values[(.5scale,.5scale,0.)]≈.058934092430339834scale rtol=2e-11
        end
    end
    for dtype in (1e-200,1e200)
        session() do _
            fixture(:tri)
            _,values=samples(output(dtype))
            expected=dtype<1 ? -log1p(-.9999)*dtype*sqrt(2.) : sqrt(2.)/(24dtype)
            @test values[(.5,.5,0.)]≈expected rtol=2e-11
        end
    end
    session() do _
        fixture(:tri)
        for (uid,x,y) in ((11,0.,0.),(21,1e100,0.),(31,1e100,1.),
                (41,0.,1.),(51,.5e100,.5))
            A.mesh.set_node(uid,[x,y,0.])
        end
        _,values=samples(output(1e150))
        expected=setprecision(BigFloat,512) do
            width=BigFloat(1e100);mu=BigFloat(1e150)
            diagonal=hypot(width,BigFloat(1))
            s=2/(12mu^2*diagonal^2*(1+inv(width^2))+1)
            Float64(-mu*diagonal*log1p(-s))
        end
        @test values[(.5e100,.5,0.)]≈expected rtol=2e-11
        @test expected>0.
    end
    for (height,dtype) in ((1e-200,1e-100),(1e-200,1e-200),(1e-300,1e-300))
        session() do _
            fixture(:tri)
            for (uid,x,y) in ((11,0.,0.),(21,1.,0.),(31,1.,height),
                    (41,0.,height),(51,.5,.5height))
                A.mesh.set_node(uid,[x,y,0.])
            end
            _,values=samples(output(dtype))
            expected=setprecision(BigFloat,512) do
                h=BigFloat(height);mu=BigFloat(dtype);diagonal=hypot(BigFloat(1),h)
                s=2/(12mu^2*diagonal^2*(1+inv(h^2))+1)
                Float64(-mu*diagonal*log1p(-s))
            end
            @test values[(.5,.5height,0.)]≈expected rtol=2e-11
            @test expected>0.
        end
    end
end
end
