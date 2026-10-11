using Test
using Tessella

const _CLEAR_ORDER_API=Tessella.API

function _clear_order_records()
    api=_CLEAR_ORDER_API
    return Dict(pair=>(nodes=api.mesh.get_nodes(pair...),
                      elements=api.mesh.get_elements(pair...))
                for pair in api.model.get_entities())
end

function _clear_order_fixture()
    api=_CLEAR_ORDER_API
    api.option("Mesh.Renumber",0)
    for (point,(x,y)) in enumerate(((0.,0.),(1.,0.),(1.,1.),(0.,1.)))
        api.model.add_point(x,y,0.;tag=point)
    end
    for (curve,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        api.model.add_line(a,b;tag=curve)
        api.mesh.set_transfinite_curve(curve,3)
    end
    api.model.add_curve_loop([1,2,3,4];tag=1)
    api.model.add_plane_surface([1];tag=1)
    api.mesh.set_transfinite_surface(1)
    api.mesh.generate(2)
    api.mesh.refine()
end

@testset "Selective clear retains discrete boundary records" begin
    api=_CLEAR_ORDER_API
    for (selection,changed) in (([(0,1)],Tuple{Int,Int}[]),
                               ([(1,1),(0,1)],[(0,1),(1,1)]),
                               ([(0,1),(1,1)],[(1,1)]))
        api.initialize()
        try
            api.model.add_discrete_entity(0,1)
            api.model.add_discrete_entity(0,2)
            api.model.add_discrete_entity(1,1,[1,2])
            api.mesh.add_nodes(0,1,[101],[0.,0.,0.])
            api.mesh.add_nodes(0,2,[102],[1.,0.,0.])
            api.mesh.add_nodes(1,1,[202],[.5,0.,0.],[.5])
            api.mesh.add_elements_by_type(1,15,[401],[101])
            api.mesh.add_elements_by_type(2,15,[402],[102])
            api.mesh.add_elements_by_type(1,1,[501,502],[101,202,202,102])
            before=_clear_order_records()
            api.mesh.clear(selection)
            after=_clear_order_records()
            @test sort!([pair for pair in keys(before) if before[pair]!=after[pair]])==changed
            for pair in keys(before)
                if pair in changed
                    @test after[pair].nodes==(UInt64[],Float64[],Float64[])
                    @test after[pair].elements==(Int32[],Vector{UInt64}[],Vector{UInt64}[])
                else
                    @test after[pair]==before[pair]
                end
            end
            nodes=api.mesh.get_nodes()[1]
            @test allunique(nodes)
            _,_,connectivity=api.mesh.get_elements()
            @test all(node->node in nodes,Iterators.flatten(connectivity))
        finally
            api.finalize()
        end
    end
end

@testset "Selective clear admits volume boundaries after their parent" begin
    api=_CLEAR_ORDER_API
    source="""
Geometry.AutoCoherence=1;Mesh.Renumber=0;
Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={1,1,0,1};Point(4)={0,1,0,1};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
Transfinite Curve{1,2,3,4}=2;Transfinite Surface{1} Left;
sweep[]=Extrude{0,0,1}{Surface{1};Layers{1};Recombine;QuadTriNoNewVerts;};
Physical Surface(7)={1};
"""
    for (selection,changed) in (([(2,1),(1,1),(0,1)],Tuple{Int,Int}[]),
                               ([(3,1),(2,1),(1,1)],[(2,1),(3,1)]),
                               ([(2,1),(3,1)],[(3,1)]),
                               ([(3,1),(2,1)],[(2,1),(3,1)]))
        api.initialize()
        try
            mktemp() do path,io
                write(io,source);close(io)
                api.open_geo!(path)
            end
            api.mesh.generate(3)
            before=_clear_order_records()
            nodes=api.mesh.get_nodes()
            cache=api.LAST_MESH[]
            tets=only(filter(block->block.msh==4,cache.blocks)).nodes
            @test size(tets,2)==6
            volume=sum(eachcol(tets)) do cell
                a,b,c,d=(cache.coords[:,row] for row in cell)
                x=b-a;y=c-a;z=d-a
                signed=(x[1]*(y[2]*z[3]-y[3]*z[2])-
                        x[2]*(y[1]*z[3]-y[3]*z[1])+
                        x[3]*(y[1]*z[2]-y[2]*z[1]))/6
                @test signed>0
                signed
            end
            @test volume≈1 atol=2eps(1.) rtol=0
            api.mesh.clear(selection)
            after=_clear_order_records()
            @test sort!([pair for pair in keys(before) if before[pair]!=after[pair]])==changed
            @test api.mesh.get_nodes()==nodes
            for pair in keys(before)
                if pair in changed
                    @test after[pair].nodes==(UInt64[],Float64[],Float64[])
                    @test after[pair].elements==(Int32[],Vector{UInt64}[],Vector{UInt64}[])
                else
                    @test after[pair]==before[pair]
                end
            end
            @test Tessella.validate(api.LAST_MESH[]).ok
        finally
            api.finalize()
        end
    end
end

@testset "Selective clear admits boundaries in request order" begin
    # Pinned Gmsh 4.15.2: a Curve/Point remains intact while a direct upward
    # entity owns cells. Earlier requests can clear that parent, making a
    # later request eligible. Every surviving record keeps labels and params.
    cases=(
        ([(1,1),(0,1)],Tuple{Int,Int}[],25),
        ([(0,1)],Tuple{Int,Int}[],25),
        ([(2,1)],[(2,1)],16),
        ([(2,1),(1,1),(0,1)],[(1,1),(2,1)],13),
        ([(0,1),(1,1),(2,1)],[(2,1)],16),
        ([(1,1)],Tuple{Int,Int}[],25),
        ([(2,1),(1,1)],[(1,1),(2,1)],13),
        ([(1,1),(0,1),(2,1)],[(2,1)],16),
        ([(2,1),(1,1),(1,2),(1,3),(1,4),(0,1)],
         [(0,1),(1,1),(1,2),(1,3),(1,4),(2,1)],3),
        ([(2,1),(1,1),(0,1),(1,4),(0,1)],
         [(0,1),(1,1),(1,4),(2,1)],9))
    for (selection,changed,count) in cases
        _CLEAR_ORDER_API.initialize()
        try
            _clear_order_fixture()
            before=_clear_order_records()
            all_before=_CLEAR_ORDER_API.mesh.get_nodes()
            _CLEAR_ORDER_API.mesh.clear(selection)
            after=_clear_order_records()
            @test sort!([pair for pair in keys(before) if before[pair]!=after[pair]])==changed
            @test length(_CLEAR_ORDER_API.mesh.get_nodes()[1])==count
            @test _CLEAR_ORDER_API.model.get_entities()==sort!(collect(keys(before)))
            for pair in keys(before)
                if pair in changed
                    @test after[pair].nodes==(UInt64[],Float64[],Float64[])
                    @test after[pair].elements==(Int32[],Vector{UInt64}[],Vector{UInt64}[])
                else
                    @test after[pair]==before[pair]
                end
            end
            isempty(changed) && @test _CLEAR_ORDER_API.mesh.get_nodes()==all_before
            nodes=_CLEAR_ORDER_API.mesh.get_nodes()[1]
            @test allunique(nodes)
            _,_,connectivity=_CLEAR_ORDER_API.mesh.get_elements()
            @test all(node->node in nodes,Iterators.flatten(connectivity))
            cache=_CLEAR_ORDER_API.LAST_MESH[]
            @test Tessella.validate(cache).ok
            if isempty(changed)
                triangles=only(filter(block->block.msh==2,cache.blocks)).nodes
                area=sum(eachcol(triangles)) do cell
                    a,b,c=(cache.coords[:,row] for row in cell)
                    signed=((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]))/2
                    @test signed>0
                    exact_a,exact_b,exact_c=(Rational{BigInt}.(point) for point in (a,b,c))
                    exact_signed=((exact_b[1]-exact_a[1])*(exact_c[2]-exact_a[2])-
                                  (exact_b[2]-exact_a[2])*(exact_c[1]-exact_a[1]))/2
                    @test exact_signed>0
                    exact_signed
                end
                @test area==1
            end
        finally
            _CLEAR_ORDER_API.finalize()
        end
    end
    _CLEAR_ORDER_API.initialize()
    try
        _clear_order_fixture()
        before=_clear_order_records()
        @test_throws ArgumentError _CLEAR_ORDER_API.mesh.clear([(2,1),(1,1),(1,99)])
        @test _clear_order_records()==before
    finally
        _CLEAR_ORDER_API.finalize()
    end
end
