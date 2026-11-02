module APIMetricSurfaceEmbeddedChecks
using Test, Tessella
const A=Tessella.API

function coordinates()
    tags,xyz,_=A.mesh.get_nodes()
    @test allunique(tags)
    Dict(uid=>Tuple(xyz[3i-2:3i]) for (i,uid) in enumerate(tags))
end

function lower_edges(curve)
    kinds,_,nodes=A.mesh.get_elements(1,curve)
    @test kinds==Int32[1]
    Set(minmax(cell...) for cell in eachcol(reshape(only(nodes),2,:)))
end

function check_embedded(kind;recombine=nothing,local_algorithm=false)
    A.initialize()
    try
        for (name,value) in (("General.Terminal",0),("Mesh.Renumber",0),
                ("Mesh.Algorithm",local_algorithm ? 6 : 7),("Mesh.Smoothing",0),
                ("Mesh.MeshSizeFromPoints",0),("Mesh.MeshSizeFromCurvature",0),
                ("Mesh.MeshSizeExtendFromBoundary",0))
            A.option(name,value)
        end
        A.model.add("owned metric embedding")
        corners=[A.model.add_point(x,y,0.;meshSize=.2) for (x,y) in
                 ((0.,0.),(1.,0.),(1.,1.),(0.,1.))]
        boundary=[A.model.add_line(corners[i],corners[mod1(i+1,4)]) for i in 1:4]
        surface=A.model.add_plane_surface([A.model.add_curve_loop(boundary)])
        fullquad=recombine in (2,3)
        for curve in boundary;A.mesh.set_transfinite_curve(curve,fullquad ? 7 : 6);end
        local_algorithm && A.mesh.set_algorithm(2,surface,7)
        recombine===nothing || begin
            A.option("Mesh.RecombinationAlgorithm",recombine)
            A.mesh.set_recombine(2,surface)
        end
        embedded_points=Int[];embedded_curves=Int[]
        if kind in (:curve,:point_on_curve,:ungraded_curve,:coarse_curve)
            first=A.model.add_point(.2,.5,0.;meshSize=.2)
            last=A.model.add_point(.8,.5,0.;meshSize=.2)
            curve=A.model.add_line(first,last)
            push!(embedded_curves,curve)
            kind===:ungraded_curve || A.mesh.set_transfinite_curve(curve,
                kind===:coarse_curve ? 2 : fullquad ? 5 : 4)
            A.mesh.embed(1,[curve],2,surface)
        end
        if !(kind in (:curve,:ungraded_curve,:coarse_curve))
            x,y=kind===:point_on_curve ? (fullquad ? .425 : .5,.5) :
                kind===:point_on_boundary ? (fullquad ? .4375 : .5,0.) : (.37,.43)
            point=A.model.add_point(x,y,0.;meshSize=.2)
            push!(embedded_points,point)
            A.mesh.embed(0,[point],2,surface)
        end
        diagonal=kind===:ungraded_curve ? (1.,16.,1.) : (16.,4.,1.)
        field=A.mesh.field.add("MathEvalAniso")
        for (name,value) in zip(("M11","M22","M33","M12","M13","M23"),
                               (string.(diagonal)...,"0","0","0"))
            A.mesh.field.set_string(field,name,value)
        end
        A.mesh.field.set_as_background_mesh(field)
        A.mesh.generate(1)
        before=coordinates()
        point_uids=UInt64[]
        for point in embedded_points
            tags,_,_=A.mesh.get_nodes(0,point)
            @test length(tags)==1
            append!(point_uids,tags)
        end
        if kind===:coarse_curve
            @test_throws ArgumentError A.mesh.generate(2)
            return nothing
        end
        A.mesh.generate(2)
        after=coordinates()
        @test all(get(after,uid,nothing)==point for (uid,point) in before)
        kinds,_,connectivity=A.mesh.get_elements(2,surface)
        @test all(in((2,3)),kinds) && !isempty(kinds)
        fullquad && @test kinds==Int32[3]
        incidence=Dict{Tuple{UInt64,UInt64},Int}()
        used=Set{UInt64}();areas=Float64[];count=0
        for (family,nodes) in zip(kinds,connectivity)
            width=family==2 ? 3 : 4
            for cell in eachcol(reshape(nodes,width,:))
                count+=1;union!(used,cell)
                signed=0.
                for i in 1:width
                    a,b=cell[i],cell[mod1(i+1,width)]
                    edge=minmax(a,b);incidence[edge]=get(incidence,edge,0)+1
                    p,q=after[a],after[b]
                    signed+=p[1]*q[2]-p[2]*q[1]
                end
                push!(areas,signed/2)
            end
        end
        @test all(>(0),areas) && isapprox(sum(areas),1.;atol=2e-12)
        @test all(in((1,2)),values(incidence))
        @test length(used)-length(incidence)+count==1
        @test all(in(used),point_uids)
        for curve in boundary
            @test all(get(incidence,edge,0)==1 for edge in lower_edges(curve))
        end
        for curve in embedded_curves
            edges=lower_edges(curve)
            @test !isempty(edges) && all(get(incidence,edge,0)==2 for edge in edges)
            # This x-directed span has metric length .6 although the smallest
            # scalar size is .25. The host's Algorithm 7 selector must reach
            # its lower producer before the Surface consumes those cells.
            kind===:ungraded_curve && @test length(edges)==1
            if kind===:point_on_curve
                @test all(any(uid in edge for edge in edges) for uid in point_uids)
                fullquad || @test length(edges)==4
            end
        end
        if recombine===nothing
            @test all(keys(incidence)) do (a,b)
                p,q=after[a],after[b]
                sqrt(sum(diagonal[i]*(p[i]-q[i])^2 for i in 1:3))<=1+1e-12
            end
        end
    finally
        A.finalize()
    end
end

@testset verbose=true "Algorithm 7 owns embedded Point and Curve products" begin
    for kind in (:point,:curve,:point_on_curve,:point_on_boundary)
        @testset "$kind" begin;check_embedded(kind);end
    end
    @testset "embedded Curve uses its host face metric selector" begin
        check_embedded(:ungraded_curve;local_algorithm=true)
    end
    @testset "coarse owned embedded Curve keeps the strict metric certificate" begin
        check_embedded(:coarse_curve)
    end
    for algorithm in (0,1,2,3)
        @testset "embedded Curve recombination $algorithm" begin
            check_embedded(:curve;recombine=algorithm)
        end
    end
    for algorithm in (2,3)
        @testset "Point splits embedded Curve during full-quad $algorithm" begin
            check_embedded(:point_on_curve;recombine=algorithm)
        end
        @testset "Point splits boundary Curve during full-quad $algorithm" begin
            check_embedded(:point_on_boundary;recombine=algorithm)
        end
    end
end
end
