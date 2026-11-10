# `plugin("BoundaryLayer")` — Gmsh 4.15.2 `BoundaryLayer.cpp`: parses the
# space/comma-separated entity selections (non-numeric words and unknown
# entity tags are silently skipped), computes geometric layer widths, then
# grows layer nodes/quads and untangles the affected surface. These checks
# exercise the actual in-place mesh modification.
module PluginBoundaryLayerChecks
using Test,Tessella
const A=Tessella.API

const FIXTURE=joinpath(@__DIR__,"..","fixtures","bl_unit_square_31.msh")

function _square(name)
    A.model.add(name)
    A.model.add_discrete_entity(1,1)
    A.model.add_discrete_entity(2,2)
    A.mesh.add_nodes(2,2,[1,2,3],
        [0.,0.,0.,  1.,0.,0.,  0.,1.,0.])
    A.mesh.add_elements_by_type(2,2,[10],[1,2,3])
end

@testset "Plugin(BoundaryLayer) exact option table" begin
    A.initialize()
    try
        _square("bl")
        # Upstream numeric options (BoundaryLayer.cpp).
        A.plugin.set_number("BoundaryLayer","Thickness",0.05)
        A.plugin.set_number("BoundaryLayer","Size",0.005)
        A.plugin.set_number("BoundaryLayer","Ratio",1.1)
        A.plugin.set_number("BoundaryLayer","Run",0.0)
        # Upstream string options holding entity selections.
        for opt in ("Volumes","Surfaces","Curves","Points")
            A.plugin.set_string("BoundaryLayer",opt,"")
        end
        A.plugin.set_string("BoundaryLayer","Curves","1, 2,3")
        @test_throws ArgumentError A.plugin.set_number(
            "BoundaryLayer","NoSuchOption",1.)
        @test_throws ArgumentError A.plugin.set_string(
            "BoundaryLayer","NoSuchOption","x")
    finally
        A.finalize()
    end
end

function counts()
    types,tags,nodes=A.mesh.get_elements()
    result=Dict{Int,Int}()
    for (type,tag) in zip(types,tags)
        result[Int(type)]=get(result,Int(type),0)+length(tag)
    end
    return result
end

function conforming_boundary()
    types,_,groups=A.mesh.get_elements()
    edges=Dict{Tuple{UInt64,UInt64},Int}()
    for (type,flat) in zip(types,groups)
        n=type==2 ? 3 : type==3 ? 4 : 0
        n==0 && continue
        for element in 0:length(flat)÷n-1,corner in 1:n
            a=flat[n*element+corner]
            b=flat[n*element+(corner==n ? 1 : corner+1)]
            key=minmax(a,b)
            edges[key]=get(edges,key,0)+1
        end
    end
    lines=Set{Tuple{UInt64,UInt64}}()
    for (type,flat) in zip(types,groups)
        type==1 || continue
        for element in 0:length(flat)÷2-1
            push!(lines,minmax(flat[2element+1],flat[2element+2]))
        end
    end
    return all(pair->pair[2]>1 || pair[1] in lines,edges)
end

function coordinates()
    tags,flat,_=A.mesh.get_nodes()
    d=Dict{UInt64,NTuple{3,Float64}}()
    for (index,tag) in enumerate(tags)
        d[tag]=(flat[3(index-1)+1],flat[3(index-1)+2],flat[3(index-1)+3])
    end
    return d
end

function run_layers(curves;thickness=.05,size=.01,ratio=1.2,surfaces="1")
    A.plugin.set_number("BoundaryLayer","Thickness",thickness)
    A.plugin.set_number("BoundaryLayer","Size",size)
    A.plugin.set_number("BoundaryLayer","Ratio",ratio)
    A.plugin.set_string("BoundaryLayer","Curves",curves)
    A.plugin.set_string("BoundaryLayer","Surfaces",surfaces)
    A.plugin.set_string("BoundaryLayer","Points","")
    return A.plugin.run("BoundaryLayer")
end

function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(FIXTURE)
        f()
    finally
        A.finalize()
    end
end

@testset "BoundaryLayer grows genuine quads into the selected curve" begin
    session() do
        before=coordinates()
        model=A.model.get_current()
        run_layers("1")
        types,tags,_=A.mesh.get_elements()
        now=counts()
        @test now[1]==22 && now[2]==44 && now[3]==12 && now[15]==4
        @test length(A.mesh.get_nodes()[1])==46
        @test [sum(length,A.mesh.get_elements(1,entity)[2]) for entity in 1:4]==
            [4,7,4,7]
        @test conforming_boundary()
        @test A.model.get_current()==model
        # Every original wall node stays fixed.
        for tag in (1,2,5,6,7)
            @test coordinates()[UInt64(tag)]==before[UInt64(tag)]
        end
        # Original triangles and their identities survive retiming.
        tri_tags=reduce(vcat,tags[findall(==(Int32(2)),types)])
        quad_tags=reduce(vcat,tags[findall(==(Int32(3)),types)])
        @test length(tri_tags)==44
        @test length(quad_tags)==12
        tri_quality=A.mesh.get_element_qualities(tri_tags,"minSICN")
        quad_quality=A.mesh.get_element_qualities(quad_tags,"minSICN")
        @test minimum(tri_quality)>0.8
        @test minimum(quad_quality)>0.07
        # The layer front is actually inside the square.
        for tag in 38:46
            point=coordinates()[UInt64(tag)]
            @test 0<point[2]<.04
        end
    end
end

@testset "BoundaryLayer supports two selected curves and layer stopping" begin
    session() do
        run_layers("1,4")
        now=counts()
        @test now[1]==22 && now[2]==44 && now[3]==24 && now[15]==4
        @test length(A.mesh.get_nodes()[1])==58
        for entity in 1:4
            expected=entity in (1,4) ? 3 : 6
            @test length(A.mesh.get_nodes(1,entity,false,false)[1])==expected
        end
        # Larger first layer but a total thickness that admits only two layers.
        A.open(FIXTURE)
        run_layers("1";thickness=.03,size=.01,ratio=1.2)
        now=counts()
        @test now[3]==8
        @test length(A.mesh.get_nodes()[1])==41
    end
    # Adjacent selected curves share one mitered layer column at the corner.
    session() do
        run_layers("1,2")
        now=counts()
        @test now[1]==22 && now[2]==44 && now[3]==24
        @test length(A.mesh.get_nodes()[1])==58
        @test length(A.mesh.get_nodes(1,3,false,false)[1])==6
        @test length(A.mesh.get_nodes(1,4,false,false)[1])==6
        tri_tags=reduce(vcat,A.mesh.get_elements()[2][findall(==(Int32(2)),
            A.mesh.get_elements()[1])])
        @test minimum(A.mesh.get_element_qualities(tri_tags,"minSICN"))>0.
    end
end

@testset "BoundaryLayer no-op and rollback preserve the source mesh" begin
    session() do
        tags,before,_=A.mesh.get_nodes()
        source=copy(before)
        run_layers("1";thickness=.005)
        now=counts()
        @test now[1]==16 && now[2]==44 && get(now,3,0)==0
        @test A.mesh.get_nodes()[1]==tags
        @test A.mesh.get_nodes()[2]==source
        @test_throws ArgumentError run_layers("1";thickness=0.)
        @test_throws ArgumentError run_layers("1";ratio=1e-300)
        empty=counts()
        @test empty[1]==16 && empty[2]==44 && get(empty,3,0)==0
        @test A.mesh.get_nodes()[1]==tags
        @test A.mesh.get_nodes()[2]==source
    end
end

@testset "BoundaryLayer modifies a model-generated surface mesh" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.MeshSizeMin",.25)
        A.option("Mesh.MeshSizeMax",.25)
        A.model.occ.add_rectangle(0.,0.,0.,1.,1.)
        A.model.occ.synchronize()
        A.mesh.generate(2)
        before=counts()
        @test before[2]>0
        run_layers("1";surfaces="1")
        now=counts()
        @test now[1]==before[1]+6
        @test now[2]==before[2]
        @test now[3]==12
        @test length(A.mesh.get_nodes()[1])==length(coordinates())
        tri_tags=reduce(vcat,A.mesh.get_elements()[2][findall(==(Int32(2)),
            A.mesh.get_elements()[1])])
        @test minimum(A.mesh.get_element_qualities(tri_tags,"minSICN"))>0.
    finally
        A.finalize()
    end
end

@testset "BoundaryLayer retimes only the selected surface" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("shared")
        A.model.add_discrete_entity(1,10)
        A.model.add_discrete_entity(2,1)
        A.model.add_discrete_entity(2,2)
        A.mesh.add_nodes(2,1,[1,2,3,4,5],
            [0.,0.,0., 1.,0.,0., 0.,1.,0., 2.,0.,0., 1.,1.,0.])
        A.mesh.add_elements_by_type(10,1,[101],[2,3])
        A.mesh.add_elements_by_type(1,2,[1],[1,2,3])
        A.mesh.add_elements_by_type(2,2,[2],[2,4,5])
        before=A.mesh.get_elements(2,2)
        run_layers("10")
        @test A.mesh.get_elements(2,2)==before
        @test counts()[3]==3
    finally
        A.finalize()
    end
end

@testset "BoundaryLayer orientation works outside the x/y plane" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("xz")
        for tag in 1:4;A.model.add_discrete_entity(1,tag);end
        A.model.add_discrete_entity(2,2)
        A.mesh.add_nodes(2,2,[1,2,3,4],
            [1e8,1e8,0., 1e8+1,1e8,0., 1e8+1,1e8,1., 1e8,1e8,1.])
        for tag in 1:4
            next=tag==4 ? 1 : tag+1
            A.mesh.add_elements_by_type(tag,1,[10+tag],[tag,next])
        end
        A.mesh.add_elements_by_type(2,2,[21,22],[1,2,3,1,3,4])
        run_layers("1")
        now=counts()
        @test now[1]==10 && now[2]==2 && now[3]==3
        @test length(A.mesh.get_nodes()[1])==10
        @test conforming_boundary()
        quad_tags=reduce(vcat,A.mesh.get_elements()[2][findall(==(Int32(3)),
            A.mesh.get_elements()[1])])
        @test minimum(A.mesh.get_element_qualities(quad_tags,"minSICN"))>0.
    finally
        A.finalize()
    end
end

@testset "BoundaryLayer follows a curved closed boundary" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("disk")
        A.model.add_discrete_entity(1,1)
        A.model.add_discrete_entity(2,2)
        count=8
        points=[(cos(2pi*i/count),sin(2pi*i/count),0.) for i in 0:count-1]
        A.mesh.add_nodes(2,2,[1,collect(2:count+1)...],
            collect(Iterators.flatten(vcat([(0.,0.,0.)],points))))
        for i in 1:count
            j=i==count ? 1 : i+1
            A.mesh.add_elements_by_type(1,1,[100+i],[i+1,j+1])
            A.mesh.add_elements_by_type(2,2,[200+i],[1,i+1,j+1])
        end
        run_layers("1";thickness=.1,size=.02)
        now=counts()
        @test now[1]==8 && now[2]==8 && now[3]==24
        @test length(A.mesh.get_nodes()[1])==33
        quad_tags=reduce(vcat,A.mesh.get_elements()[2][findall(==(Int32(3)),
            A.mesh.get_elements()[1])])
        @test minimum(A.mesh.get_element_qualities(quad_tags,"minSICN"))>0.
    finally
        A.finalize()
    end
end

@testset "Surface-only BoundaryLayer is a real no-op" begin
    session() do
        tags,before,_=A.mesh.get_nodes()
        A.plugin.set_string("BoundaryLayer","Curves","")
        A.plugin.set_string("BoundaryLayer","Surfaces","1")
        A.plugin.run("BoundaryLayer")
        @test A.mesh.get_nodes()[1]==tags
        @test A.mesh.get_nodes()[2]==before
    end
end
end
