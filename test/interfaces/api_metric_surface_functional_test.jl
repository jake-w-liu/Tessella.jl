module APIMetricSurfaceFunctionalChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const N=Tessella.MeshTypes

function rectangle(;counts=3,generate_lower=true)
    A.model.add("explicit tensor Surface")
    points=[A.model.add_point(x,y,0.;meshSize=.5) for (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))]
    curves=[A.model.add_line(points[i],points[mod1(i+1,4)]) for i in 1:4]
    loop=A.model.add_curve_loop(curves)
    surface=A.model.add_plane_surface([loop])
    A.model.add_physical_group(2,[surface];tag=31,name="material")
    for curve in curves;A.mesh.set_transfinite_curve(curve,counts);end
    for (name,value) in (("General.Terminal",0),("Mesh.Renumber",0),
            ("Mesh.MeshSizeFromPoints",0),("Mesh.MeshSizeFromCurvature",0),
            ("Mesh.MeshSizeExtendFromBoundary",0),("Mesh.Smoothing",0))
        A.option(name,value)
    end
    generate_lower && A.mesh.generate(1)
    return surface,curves
end
function tensor(expressions)
    tag=A.mesh.field.add("MathEvalAniso")
    for (name,expression) in zip(("M11","M22","M33","M12","M13","M23"),expressions)
        A.mesh.field.set_string(tag,name,string(expression))
    end
    A.mesh.field.set_as_background_mesh(tag)
    tag
end
function coordinates()
    tags,xyz,_=A.mesh.get_nodes()
    Dict(tag=>xyz[3i-2:3i] for (i,tag) in enumerate(tags))
end
function triangles(surface)
    families,labels,nodes=A.mesh.get_elements(2,surface)
    families==Int32[2] || error("expected actual Tri3 product")
    labels[1],collect(eachcol(reshape(nodes[1],3,:)))
end
function edges(cells)
    Set(minmax(cell[a],cell[b]) for cell in cells for (a,b) in ((1,2),(2,3),(3,1)))
end
function bounded(cells,xyz,matrix,target)
    all(edges(cells)) do (a,b)
        first,last=xyz[a],xyz[b];displacement=last-first
        all(point->sqrt(dot(displacement,matrix(point)*displacement))<=target+1e-13,
            (first,last,(first+last)/2))
    end
end
function signed_areas(cells,xyz)
    [(xyz[cell[2]][1]-xyz[cell[1]][1])*(xyz[cell[3]][2]-xyz[cell[1]][2])/2-
     (xyz[cell[2]][2]-xyz[cell[1]][2])*(xyz[cell[3]][1]-xyz[cell[1]][1])/2 for cell in cells]
end
state()=(A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements(),
         A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
function session(body)
    A.initialize()
    try;body();finally;A.finalize();end
end

@testset "API tensor front door publishes classified variable-metric cells" begin
    session() do
        surface,curves=rectangle()
        lower=Dict(curve=>A.mesh.get_elements(1,curve) for curve in curves)
        before_xyz=coordinates();before_node_max=A.mesh.get_max_node_tag()
        before_cell_max=A.mesh.get_max_element_tag()
        tensor(("4+12*x",1,1,0,0,0))
        output=A.mesh.generate_metric_surface(surface;max_metric_edge=2.01,max_nodes=100,max_triangles=200)
        labels,cells=triangles(surface);xyz=coordinates()
        @test validate(output).ok && N.ntris(output)==length(cells)>6
        @test bounded(cells,xyz,p->Diagonal([4+12p[1],1.,1.]),2.01)
        @test all(>(0),signed_areas(cells,xyz)) && isapprox(sum(signed_areas(cells,xyz)),1.;atol=1e-12)
        @test all(point->0<=point[1]<=1 && 0<=point[2]<=1 && point[3]==0,values(xyz))
        @test all(curve->A.mesh.get_elements(1,curve)==lower[curve],curves)
        @test all(pair->xyz[first(pair)]==last(pair),before_xyz)
        @test allunique(labels) && minimum(labels)>before_cell_max
        fresh=setdiff(collect(keys(xyz)),collect(keys(before_xyz)))
        @test !isempty(fresh) && all(tag->tag>before_node_max,fresh)
        @test all(tag->A.mesh.get_node(tag)[3:4]==(2,surface),fresh)
        @test all(uid->A.mesh.get_element(uid)[3:4]==(2,surface),labels)
        @test Set(A.mesh.get_nodes_for_physical_group(2,31)[1])==Set(A.mesh.get_nodes(2,surface,true,false)[1])
        @test A.mesh.get_max_node_tag()>=maximum(keys(xyz)) && A.mesh.get_max_element_tag()>=maximum(labels)
        retained=A.mesh.get();output.coords[1,1]=99.
        @test A.mesh.get().coords==retained.coords
        A.mesh.generate_metric_surface(surface;max_metric_edge=2.01,max_nodes=100,max_triangles=200)
        @test all(curve->A.mesh.get_elements(1,curve)==lower[curve],curves)
        @test bounded(last(triangles(surface)),coordinates(),p->Diagonal([4+12p[1],1.,1.]),2.01)
    end
end

@testset "tensor direction reaches the actual API producer" begin
    results=[]
    for diagonal in ((16,1,1),(1,16,1))
        session() do
            surface,_=rectangle();tensor((diagonal...,0,0,0))
            A.option("Mesh.Algorithm",7)
            before=state()
            @test_throws ArgumentError A.mesh.generate(2)
            @test state()==before
            A.mesh.generate_metric_surface(surface;max_metric_edge=2.01)
            _,cells=triangles(surface);xyz=coordinates()
            @test bounded(cells,xyz,p->Diagonal(collect(Float64.(diagonal))),2.01)
            push!(results,(xyz,deepcopy(cells)))
        end
    end
    @test results[1]!=results[2]
end

@testset "selected background factor and scalar cap apply once" begin
    results=[]
    for capped in (false,true)
        session() do
            surface,_=rectangle()
            if capped
                tensor((16,1,1,0,0,0));A.option("Mesh.MeshSizeMax",.5);A.option("Mesh.MeshSizeFactor",2)
            else
                tensor((4,1,1,0,0,0))
            end
            A.mesh.generate_metric_surface(surface;max_metric_edge=1.1)
            _,cells=triangles(surface);xyz=coordinates()
            @test bounded(cells,xyz,p->Diagonal([4.,1.,1.]),1.1)
            push!(results,(xyz,deepcopy(cells)))
        end
    end
    @test results[1]==results[2]
    session() do
        surface,_=rectangle(;generate_lower=false);tensor((4,1,1,0,0,0))
        output=A.mesh.generate_metric_surface(surface;max_metric_edge=1.1)
        @test validate(output).ok && bounded(last(triangles(surface)),coordinates(),p->Diagonal([4.,1.,1.]),1.1)
    end
end

@testset "API metric admission and production failures are atomic" begin
    session() do
        surface,_=rectangle();before=state()
        @test_throws ArgumentError A.mesh.generate_metric_surface(surface)
        @test state()==before
        field=tensor((16,1,1,0,0,0))
        for kwargs in ((max_metric_edge=0.,),(max_metric_edge=true,),
                (max_metric_edge=2.01,max_nodes=8),(max_metric_edge=2.01,max_triangles=6),
                (max_metric_edge=2.01,max_nodes=true))
            @test_throws ArgumentError A.mesh.generate_metric_surface(surface;kwargs...)
            @test state()==before
        end
        A.mesh.field.set_string(field,"M22","-1")
        @test_throws ArgumentError A.mesh.generate_metric_surface(surface;max_metric_edge=2.01)
        @test state()==before
        A.mesh.field.set_string(field,"M22","1")
        A.option("Mesh.ElementOrder",2)
        @test_throws ArgumentError A.mesh.generate_metric_surface(surface)
        @test state()==before
        A.option("Mesh.ElementOrder",1)
        @test_throws ArgumentError A.mesh.generate_metric_surface(true)
        A.mesh.set_order(2);higher=state()
        @test_throws ArgumentError A.mesh.generate_metric_surface(surface;max_metric_edge=2.01)
        @test state()==higher
    end
    session() do
        surface,_=rectangle(;counts=2);tensor((16,1,1,0,0,0));before=state()
        @test_throws ArgumentError A.mesh.generate_metric_surface(surface;max_metric_edge=2.01)
        @test state()==before
    end
end
end
