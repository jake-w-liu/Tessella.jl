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
