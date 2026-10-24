# `plugin("SpanningTree")` — Gmsh 4.15.2 `SpanningTree.cpp`: collect every
# element on the entities of the listed `PhysicalCurves`/`PhysicalSurfaces`/
# `PhysicalVolumes`, extract canonical element edges, run sorted-edge Kruskal
# over a shared union-find in dim order (curves, surfaces, volumes), and store
# the tree as MLINEs on a new discrete curve entity (max elementary tag + 1)
# in physical group `OutputPhysical` (or max physical + 1 when -1). On a
# connected mesh the tree reaches every vertex: |E| = |V| - 1.
module PluginSpanningTreeChecks
using Test,Tessella
const A=Tessella.API

function _square()
    A.model.add("sq")
    for (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))
        A.model.occ.add_point(x,y,0.)
    end
    for i in 1:4; A.model.occ.add_line(i,mod(i,4)+1); end
    A.model.occ.add_curve_loop([1,2,3,4])
    A.model.occ.add_plane_surface([1])
    A.model.occ.synchronize()
    A.model.add_physical_group(1,[1,2,3,4];tag=1)
    A.model.add_physical_group(2,[1];tag=2)
    A.mesh.generate(2)
    return nothing
end

function _tree_lines()
    # The highest dim-1 entity is the plugin's tree entity.
    curves=sort([t for (d,t) in A.model.get_entities(1)])
    tree_ent=last(curves)
    types,etags,ntags=A.mesh.get_elements(1,tree_ent)
    li=findfirst(==(Int32(1)),types)
    @test li!==nothing
    return tree_ent,UInt64.(ntags[li])
end

function _find(parent,v)
    while haskey(parent,v) && parent[v]!=v
        v=parent[v]
    end
    return v
end

@testset "SpanningTree emits a deterministic full spanning tree" begin
    A.initialize()
    try
        _square()
        ntags_all,_,_=A.mesh.get_nodes(-1,-1,false,false)
        nverts=length(ntags_all)
        A.plugin.set_string("SpanningTree","PhysicalCurves","1")
        A.plugin.set_string("SpanningTree","PhysicalSurfaces","2")
        A.plugin.run("SpanningTree")
        # New curve entity = max elementary tag + 1 (points 1..4, curves 1..4,
        # surface 1 → tag 5), new physical group = max + 1 = 3.
        @test (1,5) in A.model.get_entities(1)
        @test (1,3) in A.model.get_physical_groups(1)
        @test A.model.get_entities_for_physical_group(1,3)==[5]
        tree_ent,conn=_tree_lines()
        @test tree_ent==5
        nedges=length(conn)÷2
        @test nedges==nverts-1
        # The tree really spans every mesh vertex and carries no cycle:
        # |E| = |V|-1 with full reachability through the line graph.
        parent=Dict{UInt64,UInt64}()
        touched=Set{UInt64}()
        for i in 1:2:length(conn)
            a,b=conn[i],conn[i+1]
            push!(touched,a,b)
            ra,rb=_find(parent,a),_find(parent,b)
            @test ra!=rb          # acyclic: every edge joins two components
            parent[rb]=ra
        end
        @test touched==Set(UInt64.(ntags_all))
        roots=Set(_find(parent,v) for v in touched)
        @test length(roots)==1
        # Determinism: a second identical run appends an identical copy —
        # upstream keeps the old entity, so another entity+group appears.
        A.plugin.run("SpanningTree")
        @test (1,6) in A.model.get_entities(1)
        @test (1,4) in A.model.get_physical_groups(1)
        _,conn2=_tree_lines()
        @test Set(conn)==Set(conn2)
    finally
        A.finalize()
    end

    A.initialize()
    try
        _square()
        # Curves only: the tree spans the four boundary curves' vertices —
        # each connected curve contributes |V_c|-1 edges (they share corner
        # vertices, so the union is one component over all boundary nodes).
        A.plugin.set_string("SpanningTree","PhysicalCurves","1")
        A.plugin.set_number("SpanningTree","OutputPhysical",9.0)
        A.plugin.run("SpanningTree")
        @test (1,9) in A.model.get_physical_groups(1)
        tree_ent,conn=_tree_lines()
        # Boundary vertex count via the curve entities' nodes — boundary
        # vertices (the corner points) are classified on the points, so the
        # query needs `include_boundary`.
        bnodes=Set{UInt64}()
        for c in 1:4
            nt,_,_=A.mesh.get_nodes(1,c,true,false)
            union!(bnodes,UInt64.(nt))
        end
        @test length(conn)÷2==length(bnodes)-1
        touched=Set(conn)
        @test touched==bnodes
    finally
        A.finalize()
    end

    A.initialize()
    try
        _square()
        # No matching physicals: upstream warns "No elements found in the
        # given physcials: abording!" and creates nothing.
        A.plugin.set_string("SpanningTree","PhysicalCurves","77")
        nbefore=length(A.model.get_entities(1))
        A.plugin.run("SpanningTree")
        @test length(A.model.get_entities(1))==nbefore
    finally
        A.finalize()
    end
end

end # module
