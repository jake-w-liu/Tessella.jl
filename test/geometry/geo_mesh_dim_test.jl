using Test
using Tessella
using Tessella.MeshTypes: nnodes, nsegs, ntris, ntets, validate

function _execute_mesh_dim_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _execute_mesh_dim_save(source::AbstractString)
    return mktemp() do path,io
        write(io,source)
        close(io)
        mktemp() do out,_ignore
            execute_geo(path)
            return read(out,String)
        end
    end
end

@testset ".geo Mesh dimension statements" begin
    # `Mesh 0` is an upstream no-op: the model synchronizes but no mesh is
    # generated (`GModel::mesh(0)` early-returns after checking status).
    noop=_execute_mesh_dim_source("""
        Point(1)={0,0,0}; Point(2)={1,0,0};
        Line(1)={1,2};
        Mesh 0;
        """)
    @test noop.mesh===nothing || nnodes(noop.mesh)==0

    # `Mesh 1` grades the curve and emits point+line parts; the unit line at
    # the default characteristic length gives five segments, like pinned
    # Gmsh 4.15.2 (`Info: 6 nodes 7 elements` — five lines + two points).
    one_d=_execute_mesh_dim_source("""
        Point(1)={0,0,0}; Point(2)={1,0,0};
        Line(1)={1,2};
        Mesh 1;
        """)
    @test validate(one_d.mesh).ok
    @test nnodes(one_d.mesh)==6
    @test nsegs(one_d.mesh)==5
    @test one_d.model.curve_params[1]≈[0,0.2,0.4,0.6,0.8,1.0]

    # The merged `Mesh 1` mesh is reusable: a second `Mesh 1` re-grades from
    # the options in effect, matching upstream's unconditional `deMeshGEdge`.
    regraded=_execute_mesh_dim_source("""
        Point(1)={0,0,0}; Point(2)={1,0,0};
        Line(1)={1,2};
        Mesh 1;
        Mesh.CharacteristicLengthFactor = 0.05;
        Mesh 1;
        """)
    @test validate(regraded.mesh).ok
    @test nnodes(regraded.mesh)==91
    @test nsegs(regraded.mesh)==90

    # A coincident-endpoint line is a closed degenerate curve upstream: it
    # warns, stores a single self-loop element, and keeps both point nodes.
    coincident=_execute_mesh_dim_source("""
        Point(1)={0,0,0}; Point(2)={1,0,0};
        Line(1)={1,1};
        Mesh 1;
        """)
    @test validate(coincident.mesh;allow_degenerate_segs=true).ok
    @test nnodes(coincident.mesh)==2
    @test nsegs(coincident.mesh)==1
    @test coincident.mesh.segs[1,1]==coincident.mesh.segs[2,1]

    # `Degenerated Curve` marks the edge tooSmall — skipped entirely, like
    # upstream's `degenerate(1)` gate in `meshGEdge`.
    degenerated=_execute_mesh_dim_source("""
        Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={0,1,0};
        Line(1)={1,2}; Line(2)={1,3};
        Transfinite Curve {1} = 7;
        Degenerated Curve {2};
        Mesh 1;
        """)
    @test validate(degenerated.mesh).ok
    @test nnodes(degenerated.mesh)==8
    @test nsegs(degenerated.mesh)==6
    @test !haskey(degenerated.model.curve_params,2)
    @test degenerated.model.curve_params[1]≈collect(0:6)./6

    # `Mesh 2` reuses the stored curve discretization for the boundary PSLG
    # and writes refinement splits back onto it: every emitted line element
    # is a triangulation boundary edge (upstream `meshGFace` semantics).
    square=_execute_mesh_dim_source("""
        lc = 4 / 5;
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        MeshSize {:} = lc;
        MeshSize {2,4} = lc / 2;
        Characteristic Length {Sqrt(9)} = 3 * lc / 4;
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
        Curve Loop(1)={1:4};
        Plane Surface(1)={1};
        Mesh 2;
        """)
    @test validate(square.mesh).ok
    edges=Set{NTuple{2,Int32}}()
    for cell in axes(square.mesh.tris,2),(u,v) in ((1,2),(2,3),(3,1))
        a=square.mesh.tris[u,cell]; b=square.mesh.tris[v,cell]
        push!(edges,a<b ? (a,b) : (b,a))
    end
    for s in 1:nsegs(square.mesh)
        a=square.mesh.segs[1,s]; b=square.mesh.segs[2,s]
        @test (a<b ? (a,b) : (b,a)) in edges
    end
    for curve in 1:4
        params=square.model.curve_params[curve]
        @test issorted(params)
        @test params[1]==0 && params[end]==1
        @test length(params)>=3
    end

    # Periodic curves keep the master's discretization verbatim — the slave
    # emits the mapped copy like upstream's `copyMesh`.
    periodic=_execute_mesh_dim_source("""
        Point(1)={0,0,0}; Point(2)={1,0,0};
        Point(3)={0,1,0}; Point(4)={1,1,0};
        Line(1)={1,2}; Line(2)={3,4};
        Periodic Curve {2} = {1};
        Mesh 1;
        """)
    @test validate(periodic.mesh).ok
    @test nnodes(periodic.mesh)==18
    @test nsegs(periodic.mesh)==16
    @test periodic.model.curve_params[2]==periodic.model.curve_params[1]
end

@testset "execute_geo mesh_dim with multiple entities" begin
    # The keyword path merges every remaining entity through `_geo_mesh_model`
    # — the same pipeline a `Mesh n` statement runs — so multi-surface and
    # multi-volume models are no longer blocked and both entry points produce
    # bitwise-identical output.
    shared_body="""
        lc = 0.3;
        Point(1)={0,0,0,lc}; Point(2)={1,0,0,lc}; Point(3)={2,0,0,lc};
        Point(4)={0,1,0,lc}; Point(5)={1,1,0,lc}; Point(6)={2,1,0,lc};
        Line(1)={1,2}; Line(2)={2,3};
        Line(3)={4,5}; Line(4)={5,6};
        Line(5)={1,4}; Line(6)={2,5}; Line(7)={3,6};
        Curve Loop(1)={1,6,-3,-5};
        Plane Surface(1)={1};
        Curve Loop(2)={2,7,-4,-6};
        Plane Surface(2)={2};
        """
    merged=_execute_mesh_dim_source(shared_body;mesh_dim=2)
    @test validate(merged.mesh).ok
    # The shared boundary curve 6 emits one copy of its nodes: the merge
    # deduplicates on bitwise coordinates like upstream's `GModel::mesh`.
    seen=Set{NTuple{3,Float64}}()
    @test all(1:nnodes(merged.mesh)) do i
        key=(merged.mesh.coords[1,i],merged.mesh.coords[2,i],
             merged.mesh.coords[3,i])
        key in seen && return false
        push!(seen,key)
        return true
    end
    @test count(i->merged.mesh.coords[1,i]==1.0,1:nnodes(merged.mesh))==
        length(merged.model.curve_params[6])
    # Both surfaces carry triangles on their side of the shared edge.
    for cell in axes(merged.mesh.tris,2)
        centroid_x=sum(merged.mesh.coords[1,merged.mesh.tris[k,cell]]
                       for k in 1:3)/3
        @test 0.0<=centroid_x<=2.0
    end
    @test any(cell->sum(merged.mesh.coords[1,merged.mesh.tris[k,cell]]
                        for k in 1:3)/3<1.0,axes(merged.mesh.tris,2))
    @test any(cell->sum(merged.mesh.coords[1,merged.mesh.tris[k,cell]]
                        for k in 1:3)/3>1.0,axes(merged.mesh.tris,2))

    # `mesh_dim=2` and a mid-file `Mesh 2` run the identical pipeline.
    statement=_execute_mesh_dim_source(shared_body*"Mesh 2;\n")
    @test merged.mesh.coords==statement.mesh.coords
    @test merged.mesh.tris==statement.mesh.tris
    @test merged.mesh.segs==statement.mesh.segs

    # The merged triangle set is the union of the per-surface products.
    parts=[Tessella.mesh_model_surface(merged.model,tag)
           for tag in sort!(collect(keys(merged.model.surfaces)))]
    @test ntris(merged.mesh)==sum(ntris,parts;init=0)

    # Curved-boundary surfaces merge through the same path.
    curved=_execute_mesh_dim_source("""
        lc = 0.3;
        Point(1)={0,0,0,lc}; Point(2)={1,0,0,lc}; Point(3)={-1,0,0,lc};
        Circle(1)={2,1,3}; Circle(2)={3,1,2};
        Curve Loop(1)={1,2}; Plane Surface(1)={1};
        Point(10)={8,0,0,lc}; Point(11)={9,0,0,lc}; Point(12)={8.5,-0.1,0,lc};
        Line(7)={10,11}; Circle(8)={11,12,10};
        Curve Loop(4)={7,8}; Plane Surface(3)={4};
        """;mesh_dim=2)
    @test validate(curved.mesh).ok
    curved_parts=[Tessella.mesh_model_surface(curved.model,tag)
                  for tag in sort!(collect(keys(curved.model.surfaces)))]
    @test ntris(curved.mesh)==sum(ntris,curved_parts;init=0)

    # A periodic curved pair meshes through the merged path; the slave's node
    # coordinates stay a bitwise affine copy of the master surface nodes.
    periodic=_execute_mesh_dim_source("""
        lc = 0.35;
        Point(1)={0,0,0,lc}; Point(2)={1,0,0,lc};
        Point(3)={1,1,0,lc}; Point(4)={0,1,0,lc}; Point(5)={0.5,0.6,0,lc};
        Line(1)={1,2}; Line(2)={2,3}; Circle(3)={3,5,4}; Line(4)={4,1};
        Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
        Point(6)={3,0,0,lc}; Point(7)={4,0,0,lc};
        Point(8)={4,1,0,lc}; Point(9)={3,1,0,lc}; Point(10)={3.5,0.6,0,lc};
        Line(5)={6,7}; Line(6)={7,8}; Circle(7)={8,10,9}; Line(8)={9,6};
        Curve Loop(2)={5,6,7,8}; Plane Surface(2)={2};
        Periodic Surface {2} = {1} Translate {3,0,0};
        """;mesh_dim=2)
    @test validate(periodic.mesh).ok
    slave=Set(NTuple{3,Float64}[
        (periodic.mesh.coords[1,i],periodic.mesh.coords[2,i],
         periodic.mesh.coords[3,i])
        for i in 1:nnodes(periodic.mesh)
        if periodic.mesh.coords[1,i]>2.0])
    # Forward-map the master nodes like the periodic copy does: the affine
    # relation emits `master + 3` verbatim, so the shifted set must match the
    # slave's coordinates bitwise.
    master=Set(NTuple{3,Float64}[
        (periodic.mesh.coords[1,i]+3.0,periodic.mesh.coords[2,i],
         periodic.mesh.coords[3,i])
        for i in 1:nnodes(periodic.mesh)
        if periodic.mesh.coords[1,i]<2.0])
    @test slave==master

    # Disjoint volumes merge for `mesh_dim=3` the same way.
    volumes=_execute_mesh_dim_source("""
        SetFactory("OpenCASCADE");
        Mesh.CharacteristicLengthMin = 0.4;
        Mesh.CharacteristicLengthMax = 0.4;
        Box(1) = {0,0,0, 1,1,1};
        Box(2) = {3,0,0, 1,1,1};
        """;mesh_dim=3)
    @test validate(volumes.mesh).ok
    @test ntets(volumes.mesh)>0
    box1_nodes=count(i->volumes.mesh.coords[1,i]<2.0,1:nnodes(volumes.mesh))
    box2_nodes=count(i->volumes.mesh.coords[1,i]>2.0,1:nnodes(volumes.mesh))
    @test box1_nodes>0 && box2_nodes>0
    @test box1_nodes+box2_nodes==nnodes(volumes.mesh)
end
