module OrderedMatchingPrimaryTests
using Test, Tessella, LinearAlgebra
const O=Tessella.OrderedMatching
const SYMMETRIC=Float64[1 .5 -.5 -1 -.5 .5 0;0 1 1 0 -1 -1 0;zeros(7)']
const UNIQUE=Float64[1.05 .54 -.48 -1.10 -.49 .56 .02;
    0 .95 .81 -.03 -.90 -.89 .03;zeros(7)']
const TRIANGLES=Int32[7 7 7 7 7 7;1 2 3 4 5 6;2 3 4 5 6 1]
canonical(q)=minimum(ntuple(i->q[mod1(i+shift,length(q))],length(q))
    for shift in 0:length(q)-1)

# This high-precision acos construction is independent of the production
# raw atan/cross arithmetic. Each fixture is convex and has integer cost cuts
# separated by more than .01, so it establishes the literal integer graph.
function angular_cost(coords,q)
    return setprecision(BigFloat,256) do
        p=BigFloat.(coords);deviation=BigFloat(0)
        for index in 1:4
            u=p[:,q[mod1(index-1,4)]]-p[:,q[index]]
            v=p[:,q[mod1(index+1,4)]]-p[:,q[index]]
            angle=acos(dot(u,v)/(norm(u)*norm(v)))*180/big(pi)
            deviation=max(deviation,abs(angle-90))
        end
        value=1000exp(-(1-deviation/90))
        @test min(value-floor(value),ceil(value)-value)>big".01"
        return floor(Int,value),Float64(deviation)
    end
end

function graph(coords,triangles,tags)
    incidence=Dict{Tuple{Int32,Int32},Vector{Int}}()
    for cell in axes(triangles,2), (first,second) in ((1,2),(2,3),(3,1))
        edge=minmax(triangles[first,cell],triangles[second,cell])
        push!(get!(incidence,edge,Int[]),cell)
    end
    pairs=Tuple{Int,Int,Int,Float64}[]
    periodic=Dict{Int32,Tuple{Int,Int}}()
    for edge in sort!(collect(keys(incidence));by=e->minmax(tags[e[1]],tags[e[2]]))
        cells=incidence[edge]
        if length(cells)==2
            a,b=cells
            first=only(setdiff(triangles[:,a],collect(edge)))
            second=only(setdiff(triangles[:,b],collect(edge)))
            cost,priority=angular_cost(coords,(edge[1],first,edge[2],second))
            push!(pairs,(a,b,cost,priority))
        else
            @test length(cells)==1
            cell=only(cells)
            for vertex in edge
                previous=get(periodic,vertex,(0,0))
                if previous[1]==0
                    periodic[vertex]=(cell,0)
                elseif previous[1]!=cell
                    periodic[vertex]=(previous[1],cell)
                else
                    delete!(periodic,vertex)
                end
            end
        end
    end
    # With six internal pairs, the observed native final insertion phase
    # retains the initial public-edge order between equal finite priorities.
    sort!(pairs;by=last,alg=MergeSort)
    edges=[(a,b,cost) for (a,b,cost,_) in pairs]
    for vertex in sort!(collect(keys(periodic));by=v->tags[v])
        first,second=periodic[vertex]
        @test second!=0
        push!(edges,(first,second,100000))
    end
    @test length(edges)==12
    return edges
end

function emitted_cycle(triangles,a,b)
    first=collect(triangles[:,a]);second=collect(triangles[:,b])
    start=only(findall(vertex->!(vertex in second),first))
    other=only(setdiff(second,first))
    return canonical((first[start],first[mod1(start+1,3)],other,first[mod1(start+2,3)]))
end

function area(coords,cycle)
    p=Rational{BigInt}.(coords)
    return sum(p[1,cycle[i]]*p[2,cycle[mod1(i+1,length(cycle))]]-
        p[2,cycle[i]]*p[1,cycle[mod1(i+1,length(cycle))]] for i in eachindex(cycle))/2
end

@testset "Pinned API and executable C6 records, cost and exact partition" begin
    corpus=normpath(joinpath(@__DIR__,"..","..","validation","weighted_matching","corpus.tsv"))
    rows=readlines(corpus)[2:end]
    @test length(rows)==13
    for row in rows
        name,shift_text,reversed_text,records=split(row,'\t')
        coords=name=="unique_minimum" ? UNIQUE : SYMMETRIC
        shift=parse(Int,shift_text)
        triangles=TRIANGLES[:,mod1.(collect(1:6).+shift,6)]
        tags=collect(10003:7:10045)
        reversed_text=="1" && reverse!(tags)
        edges=graph(coords,triangles,tags)
        input=copy(edges)
        result=O.ordered_perfect_matching(6,edges)
        @test result!==nothing
        result===nothing && continue
        expected=[begin
            label,nodes=split(record,':')
            (parse(Int,label),Tuple(parse.(Int,split(nodes,','))))
        end for record in split(records,';')]
        actual=[(30012+index,emitted_cycle(triangles,edges[edge][1],edges[edge][2]))
            for (index,edge) in enumerate(result.selected)]
        @test actual==expected
        @test edges==input
        quads=last.(actual)
        @test all(area(coords,q)>0 for q in quads)
        @test sum(area(coords,q) for q in quads)==sum(area(coords,Tuple(t)) for t in eachcol(triangles))
        incidence=Dict{Tuple{Int,Int},Int}()
        for q in quads,i in 1:4
            edge=minmax(q[i],q[mod1(i+1,4)])
            incidence[edge]=get(incidence,edge,0)+1
        end
        @test all(count in (1,2) for count in values(incidence))
        @test Set(edge for (edge,count) in incidence if count==1)==Set(minmax(i,mod1(i+1,6)) for i in 1:6)
        @test 7-length(incidence)+length(quads)==1
        @test sum(edges[i][3] for i in result.selected)==(name=="unique_minimum" ? 1552 : 1542)
    end
end

# Literal directed input edges and output indices are retained in the pinned
# native corpus. These general graphs exercise certified singleton and
# fractional paths independently of surface geometry or pair sorting.
@testset "Native general graph singleton/fractional searches and odd-cycle handoff" begin
    fixtures=[
        (name="weighted_60",edges=[(1,4,1),(1,5,3),(1,6,5),(2,3,7),
            (2,4,9),(2,5,0),(2,6,2),(3,4,4)],selected=[8,6,3],singleton=true),
        (name="weighted_41",edges=[(1,2,4),(1,5,10),(1,6,1),(2,3,3),
            (2,5,7),(3,4,0)],selected=[6,5,3],singleton=true),
        (name="weighted_43",edges=[(1,2,5),(1,3,7),(1,5,0),(1,6,2),
            (2,4,6),(2,5,8),(3,4,1)],selected=[7,6,4],singleton=false),
        (name="weighted_110",edges=[(1,3,2),(1,4,4),(1,5,6),(2,4,1),
            (2,6,5),(3,4,7),(3,5,9)],selected=[6,3,5],singleton=false),
        (name="weighted_361",edges=[(1,2,10),(1,5,5),(1,6,7),(2,3,9),
            (2,4,0),(2,5,2),(2,6,4),(3,4,6),(3,5,8),(4,5,1)],
            selected=[8,2,7],singleton=false)]
    corpus=readlines(joinpath(@__DIR__,"..","..","validation","native_concorde",
        "native_concorde_integral_literal_v1.txt"))
    for fixture in fixtures
        first_row=only(findall(line->startswith(line,"CASE "*fixture.name*" "),corpus))
        header=split(corpus[first_row]);count,ne,nm=parse.(Int,header[3:5])
        native_edges=[Tuple(parse.(Int,split(corpus[first_row+i])[2:4])) for i in 1:ne]
        native_selected=[parse(Int,last(split(corpus[first_row+ne+i]))) for i in 1:nm]
        @test count==6 && native_edges==fixture.edges
        @test native_selected==fixture.selected
        before=copy(fixture.edges);graph=O._graph(count,fixture.edges)
        @test !O._bipartite(graph)
        singleton=O._bipartite_perfect(graph;reject_odd_cycles=true)
        @test (singleton!==nothing)==fixture.singleton
        if fixture.name in ("weighted_110","weighted_361")
            # These native cases first close a half-cycle, then enter and
            # reopen it from another free root to finish an integral cover.
            @test O._fractional_singleton_perfect(graph)!==nothing
        end
        result=O.ordered_perfect_matching(count,fixture.edges)
        @test result!==nothing && result.selected==native_selected
        @test fixture.edges==before
        @test sort([vertex for i in result.selected for vertex in fixture.edges[i][1:2]])==collect(1:count)
        @test sum(fixture.edges[i][3] for i in result.selected)==parse(Int,last(header))
    end
end
end # module OrderedMatchingPrimaryTests
