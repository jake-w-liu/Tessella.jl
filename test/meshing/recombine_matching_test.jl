module RecombineMatchingRegressionTests
using Test,Tessella
const API=Tessella.API
const R=Tessella.Recombine

# Independent enumeration permits leaving a vertex unmatched or pairing it
# with each available neighbor. It does not use an augmenting-path algorithm.
function exact_matching_size(adj,mask::UInt32)
    mask==0 && return 0
    first=trailing_zeros(mask)+1
    remainder=mask & ~(UInt32(1)<<(first-1))
    best=exact_matching_size(adj,remainder)
    for second in adj[first]
        bit=UInt32(1)<<(second-1)
        remainder & bit==0 && continue
        best=max(best,1+exact_matching_size(adj,remainder & ~bit))
    end
    return best
end

@testset "Edmonds independent exhaustive graphs through six vertices" begin
    for n in 0:6
        edges=[(i,j) for i in 1:n for j in i+1:n]
        for graph in UInt32(0):(UInt32(1)<<length(edges))-UInt32(1)
            adj=[Int[] for _ in 1:n]
            for (index,(a,b)) in enumerate(edges)
                graph & (UInt32(1)<<(index-1))==0 && continue
                push!(adj[a],b);push!(adj[b],a)
            end
            mate=R._edmonds_matching(n,adj)
            @test all(mate[v]==0 ||
                (mate[v]!=v && mate[v] in adj[v] && mate[mate[v]]==v)
                for v in eachindex(mate))
            @test count(!iszero,mate)÷2==exact_matching_size(adj,(UInt32(1)<<n)-UInt32(1))
        end
    end
end

@testset "Edmonds linear temporary allocation" begin
    previous=0
    for n in (256,512,1024)
        # Every two-node connected component has exactly one possible pair.
        adj=[Int[iseven(i) ? i-1 : i+1] for i in 1:n]
        expected=Int[iseven(i) ? i-1 : i+1 for i in 1:n]
        @test R._edmonds_matching(n,adj)==expected
        allocated=minimum(@allocated(R._edmonds_matching(n,adj)) for _ in 1:3)
        @test allocated<=35*n+2048
        previous==0 || @test allocated<=2.15*previous+4096
        previous=allocated
    end
end

function square_islands()
    return Mesh(Float64[0 1 1 0 3 4 3 6 7 6;
                        0 0 1 1 0 0 1 0 0 1;zeros(10)'];
        tris=Int32[1 1 5 8;2 3 6 9;3 4 7 10],tri_tag=fill(Int32(7),4))
end
function star_islands(;islands=false)
    coords=Float64[0 1 .5 .5 1.2 -.2;0 0 1 -.5 .8 .8;zeros(6)']
    tris=Int32[1 2 3 1;2 1 2 3;3 4 5 6]
    if islands
        coords=hcat(coords,Float64[7 8 7 10 11 10;0 0 1 0 0 1;zeros(6)'])
        tris=hcat(tris,Int32[7 10;8 11;9 12])
    end
    return Mesh(coords;tris,tri_tag=fill(Int32(7),size(tris,2)))
end
function rational_area(mesh)
    points=Rational{BigInt}.(mesh.coords)
    blocks=mesh isa Mesh ? (mesh.tris,) :
        (block.nodes for block in mesh.blocks if block.msh in (2,3))
    total=zero(Rational{BigInt})
    for nodes in blocks,column in axes(nodes,2),i in axes(nodes,1)
        first=nodes[i,column];second=nodes[mod1(i+1,size(nodes,1)),column]
        total+=(points[1,first]*points[2,second]-points[1,second]*points[2,first])/2
    end
    return total
end
canonical(nodes)=minimum(ntuple(i->nodes[mod1(i+shift,length(nodes))],length(nodes))
                         for shift in 0:length(nodes)-1)

@testset "Full boundary closure controls even Blossom retry" begin
    mesh=square_islands();snapshot=(copy(mesh.coords),copy(mesh.tris),copy(mesh.tri_tag))
    for angle in (0.,.5,1.)
        result=recombine_triangles(mesh;algorithm=:blossom,recombine_angle=angle)
        @test only(result.blocks).msh==2
        @test only(result.blocks).nodes==mesh.tris
        @test only(result.blocks).tags==mesh.tri_tag
        @test rational_area(result)==2//1
        @test Tessella.Elements.validate(result).ok
        @test_throws ArgumentError recombine_triangles(mesh;algorithm=:blossom,
            recombine_angle=angle,full_quad=true)
    end
    @test (mesh.coords,mesh.tris,mesh.tri_tag)==snapshot
    for angle in (1.01,45.)
        result=recombine_triangles(mesh;algorithm=:blossom,recombine_angle=angle)
        @test [b.msh for b in result.blocks]==[2,3]
        @test result.blocks[1].nodes==mesh.tris[:,3:4]
        @test canonical(vec(result.blocks[2].nodes))==(1,2,3,4)
        @test rational_area(result)==2//1
    end
    # The star has no perfect internal matching, but its three boundary links
    # complete it. Primary Blossom therefore ignores angle zero and forms a quad.
    star=star_islands()
    result=recombine_triangles(star;algorithm=:blossom,recombine_angle=0.)
    @test [b.msh for b in result.blocks]==[2,3]
    @test canonical(vec(result.blocks[2].nodes))==(1,2,5,3)
    @test rational_area(result)==rational_area(star)
    # Adding isolated triangles makes the FULL graph infeasible. Primary retry
    # uses its raw maximum-corner-deviation priority, choosing the other ear.
    star=star_islands(;islands=true)
    result=recombine_triangles(star;algorithm=:blossom,recombine_angle=45.)
    @test canonical(vec(result.blocks[2].nodes))==(1,2,3,6)
    @test result.blocks[1].nodes==star.tris[:,[2,3,5,6]]
    @test rational_area(result)==rational_area(star)
    # No-angle callers retain maximum-cardinality semantics on disjoint islands.
    result=recombine_triangles(mesh;algorithm=:blossom)
    @test canonical(vec(result.blocks[2].nodes))==(1,2,3,4)
end

function seeded_raw(mesh;maximum=90001)
    API.initialize();API.option("Mesh.Renumber",0);API.option("Mesh.RecombinationAlgorithm",1)
    API.model.add_discrete_entity(2,1)
    tags=UInt64[10003+7i for i in 0:size(mesh.coords,2)-1]
    API.mesh.add_nodes(2,1,tags,vec(mesh.coords))
    API.mesh.add_elements_by_type(1,2,UInt64.(30001:30000+size(mesh.tris,2)),tags[vec(mesh.tris)])
    API.model.add_discrete_entity(0,99);API.mesh.add_nodes(0,99,[70001],[20.,20.,20.])
    API.mesh.add_elements_by_type(99,15,[maximum],[70001])
    API.mesh.remove_elements(0,99,[maximum])
    return tags
end
@testset "Even retry reserves both factory passes and preserves public labels" begin
    for (mesh,angle,expected_max,quad,retained) in
        ((square_islands(),.5,90003,(),[1,2,3,4]),
         (square_islands(),45.,90004,(1,2,3,4),[3,4]),
         (star_islands(),0.,90005,(1,2,5,3),[2,4]),
         (star_islands(;islands=true),45.,90008,(1,2,3,6),[2,3,5,6]))
        try
            tags=seeded_raw(mesh);before=API.mesh.get_nodes()
            API.mesh.set_recombine(2,1,angle);API.mesh.recombine()
            types,ids,nodes=API.mesh.get_elements(2,1)
            @test types==(isempty(quad) ? Int32[2] : Int32[2,3])
            @test ids[1]==UInt64.(30000 .+ retained)
            @test nodes[1]==tags[vec(mesh.tris[:,retained])]
            @test API.mesh.get_nodes()==before
            @test API.mesh.get_max_element_tag()==expected_max
            if !isempty(quad)
                @test ids[2]==UInt64[expected_max]
                @test canonical(nodes[2])==Tuple(tags[collect(quad)])
            end
        finally
            API.finalize()
        end
    end
    try
        seeded_raw(square_islands();maximum=typemax(Int32)-1)
        API.mesh.set_recombine(2,1,.5)
        model=API.CURRENT[];payload=(API.mesh.get_nodes(),API.mesh.get_elements())
        counters=(API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])
        @test_throws ArgumentError API.mesh.recombine()
        @test API.CURRENT[]===model
        @test (API.mesh.get_nodes(),API.mesh.get_elements())==payload
        @test (API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])==counters
    finally
        API.finalize()
    end
end
end # module RecombineMatchingRegressionTests
