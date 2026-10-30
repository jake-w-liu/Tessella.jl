# `BoundaryLayer Field =` topological surface pass — the field's
# `CurvesList` walls grow vertex columns off the face boundary (upstream's
# `buildAdditionalPoints2D`/`modifyInitialMeshForBoundaryLayers`), the
# certified `_bl_2d_offsets` strip tiles the layer, and the ordinary
# constrained-Delaunay fill covers the shrunken domain. Run-end columns
# ride the adjacent non-wall segments like upstream's
# `addColumnAtTheEndOfTheBL`, so the adjacent edges receive the same
# BL-graded boundary nodes upstream's `addBoundaryLayerPoints` writes.

using Test
using Tessella
using Tessella.MeshTypes: nnodes
using Tessella.Elements: validate, mixed_crc

function _blf_execute(source::AbstractString;mesh_dim=2)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _blf_error(source::AbstractString;mesh_dim=2)
    try
        ex=_blf_execute(source;mesh_dim=mesh_dim)
        !isempty(ex.errors) ?
            ArgumentError(join(ex.errors,"\n")) : nothing
    catch err
        return err
    end
    return nothing
end

const _BLF_SQUARE_HEAD="""
Point(1)={0,0,0,1};
Point(2)={1,0,0,1};
Point(3)={1,1,0,1};
Point(4)={0,1,0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Curve Loop(1)={1,2,3,4};
Plane Surface(1)={1};
"""

function _blf_signed_area(mesh)
    total=0.0
    xs=mesh.coords[1,:];ys=mesh.coords[2,:]
    for c in axes(mesh.tris,2)
        i,j,k=Int.(mesh.tris[:,c])
        total+=0.5*((xs[j]-xs[i])*(ys[k]-ys[i])-
                    (xs[k]-xs[i])*(ys[j]-ys[i]))
    end
    return total
end

# Upstream's non-beta schedule: L = hWall; while L <= thickness
# { push L; L += hWall·ratio^ith }.
function _blf_schedule(hwall::Float64,ratio::Float64,thickness::Float64)
    offsets=Float64[]
    L=hwall;k=0
    while L<=thickness
        push!(offsets,L)
        k+=1
        L+=hwall*ratio^k
    end
    return offsets
end

@testset "BoundaryLayer Field — single wall strip" begin
    hwall,ratio,thickness=0.02,1.3,0.2
    offsets=_blf_schedule(hwall,ratio,thickness)
    ex=_blf_execute(_BLF_SQUARE_HEAD*"""
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].Size=$hwall;
Field[1].Ratio=$ratio;
Field[1].Thickness=$thickness;
BoundaryLayer Field = 1;
Mesh 2;
""")
    mesh=ex.mesh
    @test validate(mesh).ok
    ys=mesh.coords[2,:]
    # The wall's column rows sit at the exact geometric schedule — the
    # rows also ride the two adjacent non-wall edges (upstream's
    # `addColumnAtTheEndOfTheBL` covered spans).
    for off in offsets
        @test any(y->isapprox(y,off;atol=1e-9),ys)
    end
    @test _blf_signed_area(mesh)≈1.0 atol=1e-9
    # The bottom edge stays one boundary segment while edges 2/4 carry
    # the graded covered nodes.
    edges=Dict{Tuple{Int,Int},Int}()
    for c in axes(mesh.tris,2)
        tri=Int.(mesh.tris[:,c])
        for (a,b) in ((tri[1],tri[2]),(tri[2],tri[3]),(tri[3],tri[1]))
            key=a<b ? (a,b) : (b,a)
            edges[key]=get(edges,key,0)+1
        end
    end
    boundary=first.(filter(p->p[2]==1,collect(edges)))
    xs=mesh.coords[1,:]
    right=count(e->abs(xs[e[1]]-1.0)<1e-12 && abs(xs[e[2]]-1.0)<1e-12,
                boundary)
    @test right==length(offsets)+1
end

@testset "BoundaryLayer Field — two disjoint runs and closed ring" begin
    ex=_blf_execute(_BLF_SQUARE_HEAD*"""
Field[1]=BoundaryLayer;
Field[1].CurvesList={1,3};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.25;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test validate(ex.mesh).ok
    @test _blf_signed_area(ex.mesh)≈1.0 atol=1e-9
    ys=sort!(unique(round.(ex.mesh.coords[2,:];digits=9)))
    for off in _blf_schedule(0.05,1.2,0.25)
        @test any(y->isapprox(y,off;atol=1e-7),ys)
        @test any(y->isapprox(y,1.0-off;atol=1e-7),ys)
    end
    ex2=_blf_execute(_BLF_SQUARE_HEAD*"""
Field[1]=BoundaryLayer;
Field[1].CurvesList={1,2,3,4};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.3;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test validate(ex2.mesh).ok
    @test _blf_signed_area(ex2.mesh)≈1.0 atol=1e-9
end

@testset "BoundaryLayer Field — Quads keeps the strip recombined" begin
    ex=_blf_execute(_BLF_SQUARE_HEAD*"""
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.25;
Field[1].Quads=1;
BoundaryLayer Field = 1;
Mesh 2;
""")
    mesh=ex.mesh
    @test validate(mesh).ok
    nquads=sum(block->block.msh==3 ? size(block.nodes,2) : 0,mesh.blocks)
    @test nquads==length(_blf_schedule(0.05,1.2,0.25))
end

@testset "BoundaryLayer Field — shared wall curve rejects like upstream" begin
    err=_blf_error("""
Point(1)={0,0,0,1};
Point(2)={1,0,0,1};
Point(3)={1,1,0,1};
Point(4)={0,1,0,1};
Point(5)={0,0,1,1};
Point(6)={1,0,1,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line(5)={1,5};
Line(6)={5,6};
Line(7)={6,2};
Curve Loop(1)={1,2,3,4};
Plane Surface(1)={1};
Curve Loop(2)={5,6,7,-1};
Plane Surface(2)={2};
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.25;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test err!==nothing
    @test occursin("Only 2D Boundary Layers are supported",sprint(showerror,err))
    # Excluding the sibling face admits the layer on surface 1.
    ex=_blf_execute("""
Point(1)={0,0,0,1};
Point(2)={1,0,0,1};
Point(3)={1,1,0,1};
Point(4)={0,1,0,1};
Point(5)={0,0,1,1};
Point(6)={1,0,1,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line(5)={1,5};
Line(6)={5,6};
Line(7)={6,2};
Curve Loop(1)={1,2,3,4};
Plane Surface(1)={1};
Curve Loop(2)={5,6,7,-1};
Plane Surface(2)={2};
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.25;
Field[1].ExcludedSurfacesList={2};
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test validate(ex.mesh).ok
end

@testset "BoundaryLayer Field — bounded loud rejections" begin
    for (opt,needle) in (
            ("SizesList={0.02,0.03}","SizesList"),
            ("FanPointsSizesList={3}","FanPointsSizesList"),
            ("IntersectMetrics=1","IntersectMetrics"),
            ("AnisoMax=1e6","AnisoMax"))
        err=_blf_error(_BLF_SQUARE_HEAD*"""
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.25;
Field[1].$opt;
BoundaryLayer Field = 1;
Mesh 2;
""")
        @test err!==nothing
        @test occursin(needle,sprint(showerror,err))
    end
    # A wall curve whose layer depth reaches the adjacent boundary
    # vertex — upstream's "Impossible BL Configuration" class.
    err=_blf_error("""
Point(1)={0,0,0,1};
Point(2)={1,0,0,1};
Point(3)={1,0.1,0,1};
Point(4)={0,0.1,0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Curve Loop(1)={1,2,3,4};
Plane Surface(1)={1};
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.3;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test err!==nothing
    @test occursin("adjacent Curve",sprint(showerror,err))
end

@testset "BoundaryLayer Field — transfinite keeps its producer" begin
    ex=_blf_execute(_BLF_SQUARE_HEAD*"""
Transfinite Curve {1,2,3,4} = 5;
Transfinite Surface {1};
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.25;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test validate(ex.mesh).ok
    # The transfinite grid owns the surface — no BL rows appear.
    ys=ex.mesh.coords[2,:]
    @test !any(y->isapprox(y,0.05;atol=1e-9) ||
               isapprox(y,0.11;atol=1e-9),ys)
end
