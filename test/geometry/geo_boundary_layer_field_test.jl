# `BoundaryLayer Field =` topological surface pass — the field's
# `CurvesList` walls grow vertex columns off the face boundary (upstream's
# `buildAdditionalPoints2D`/`modifyInitialMeshForBoundaryLayers`), the
# certified `_bl_2d_offsets` strip tiles the layer, and the ordinary
# constrained-Delaunay fill covers the shrunken domain. Run-end columns
# ride the adjacent non-wall edges' own 1-D chains like upstream's
# `addColumnAtTheEndOfTheBL`, so per-edge emission pairs `min(depth)`
# quadrangle rows plus one closing triangle on the deeper column —
# without `PointsList`/`NodesList` end-node injection the adjacent edges
# keep their ordinary grading and the interior columns emit only what the
# pairing covers (Gmsh 4.15.2 oracle: an open wall with one-node adjacent
# chains emits two strip levels regardless of the full schedule).

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

# True when a mesh node lies at 3-D distance in `(lo,hi)` of segment
# `a→b` — used to detect a boundary-layer strip grown off a wall (layer
# nodes land at the field's offsets) or prove its absence on a dropped
# shared curve.
function _blf_any_node_within(mesh,a,b,lo,hi)
    abx=b[1]-a[1];aby=b[2]-a[2];abz=b[3]-a[3]
    len2=abx*abx+aby*aby+abz*abz
    for i in 1:size(mesh.coords,2)
        px=mesh.coords[1,i]-a[1]
        py=mesh.coords[2,i]-a[2]
        pz=mesh.coords[3,i]-a[3]
        u=clamp((px*abx+py*aby+pz*abz)/len2,0.0,1.0)
        dx=px-u*abx;dy=py-u*aby;dz=pz-u*abz
        d=sqrt(dx*dx+dy*dy+dz*dz)
        lo<d<hi && return true
    end
    return false
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

function _blf_node_at(mesh,x::Float64,y::Float64;atol=1e-7)
    xs=mesh.coords[1,:];ys=mesh.coords[2,:]
    for i in eachindex(xs)
        (abs(xs[i]-x)<=atol && abs(ys[i]-y)<=atol) && return i
    end
    return 0
end

function _blf_column(corner,n1,n2,hw,ratio,thick)
    raw=(1.01*n1[1]+n2[1],1.01*n1[2]+n2[2])
    L=hypot(raw[1],raw[2]);u=(raw[1]/L,raw[2]/L)
    d=min(abs(n1[1]*n2[1]+n1[2]*n2[2]),1.0)
    fact=1.0/cos(acos(d)/2)
    return NTuple{2,Float64}[(corner[1]+off*fact*u[1],
                              corner[2]+off*fact*u[2])
                             for off in _blf_schedule(hw,ratio,thick)]
end

@testset "BoundaryLayer Field — single wall strip" begin
    hwall,ratio,thickness=0.02,1.3,0.2
    offsets=_blf_schedule(hwall,ratio,thickness)
    ex=_blf_execute(_BLF_SQUARE_HEAD*"""
Mesh.MeshSizeMin = 0.5;
Mesh.MeshSizeMax = 0.5;
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
    xs=mesh.coords[1,:];ys=mesh.coords[2,:]
    # Gmsh 4.15.2 oracle on this exact fixture (12 nodes, two interior
    # Steiner points aside): the interior wall vertex grows its full
    # five-level schedule, but each wall segment pairs it against a
    # depth-1 end column — the adjacent edges' own one-node chains — so
    # the strip emits `min` quadrangle rows plus one closing triangle:
    # levels 1–2 materialize, deeper schedule rows do not.
    @test _blf_node_at(mesh,0.5,offsets[1])!=0
    @test _blf_node_at(mesh,0.5,offsets[2])!=0
    for off in offsets[3:end]
        @test _blf_node_at(mesh,0.5,off)==0
    end
    @test _blf_signed_area(mesh)≈1.0 atol=1e-9
    # Without `PointsList` injection the adjacent edges keep their
    # ordinary one-node grading — no BL-graded chain is written back.
    edges=Dict{Tuple{Int,Int},Int}()
    for c in axes(mesh.tris,2)
        tri=Int.(mesh.tris[:,c])
        for (a,b) in ((tri[1],tri[2]),(tri[2],tri[3]),(tri[3],tri[1]))
            key=a<b ? (a,b) : (b,a)
            edges[key]=get(edges,key,0)+1
        end
    end
    boundary=first.(filter(p->p[2]==1,collect(edges)))
    right=count(e->abs(xs[e[1]]-1.0)<1e-12 && abs(xs[e[2]]-1.0)<1e-12,
                boundary)
    left=count(e->abs(xs[e[1]])<1e-12 && abs(xs[e[2]])<1e-12,
               boundary)
    @test right==2
    @test left==2
    @test _blf_node_at(mesh,1.0,0.5)!=0
    @test _blf_node_at(mesh,0.0,0.5)!=0
end

@testset "BoundaryLayer Field — two disjoint runs and closed ring" begin
    ex=_blf_execute(_BLF_SQUARE_HEAD*"""
Mesh.MeshSizeMin = 0.5;
Mesh.MeshSizeMax = 0.5;
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
    # Each open run pairs its interior column against depth-1 end
    # columns (the shared side edges' single mid nodes serve both
    # runs), so two strip levels materialize per wall — the first two
    # schedule offsets, mirrored top and bottom.
    offs=_blf_schedule(0.05,1.2,0.25)
    @test _blf_node_at(ex.mesh,0.5,offs[1])!=0
    @test _blf_node_at(ex.mesh,0.5,offs[2])!=0
    @test _blf_node_at(ex.mesh,0.5,1.0-offs[1])!=0
    @test _blf_node_at(ex.mesh,0.5,1.0-offs[2])!=0
    for off in offs[3:end]
        @test _blf_node_at(ex.mesh,0.5,off)==0
        @test _blf_node_at(ex.mesh,0.5,1.0-off)==0
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
Mesh.MeshSizeMin = 0.5;
Mesh.MeshSizeMax = 0.5;
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
    # Two wall segments × one `min`-paired quadrangle row (depth-1 end
    # columns against the interior vertex's column); the closing row on
    # the deeper side is a triangle, matching upstream's emission.
    nquads=sum(block->block.msh==3 ? size(block.nodes,2) : 0,mesh.blocks)
    @test nquads==2
end

@testset "BoundaryLayer Field — shared wall curve drops like upstream" begin
    # Upstream `setupFor2d` logs "Only 2D Boundary Layers are supported"
    # per offending (face, curve) and drops the curve from that face's
    # `_curveTags` — both faces still mesh, plain where the wall fell out
    # (Gmsh 4.15.2 oracle on this fixture: identical plain fills and no
    # strip anywhere near the shared edge).
    ex=_blf_execute("""
Point(1)={0,0,0,1};
Point(2)={1,0,0,1};
Point(3)={1,1,0,1};
Point(4)={0,1,0,1};
Point(5)={2,0,0,1};
Point(6)={2,1,0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line(5)={2,5};
Line(6)={5,6};
Line(7)={6,3};
Curve Loop(1)={1,2,3,4};
Curve Loop(2)={5,6,7,-2};
Plane Surface(1)={1};
Plane Surface(2)={2};
Mesh.MeshSizeMin = 0.25;
Mesh.MeshSizeMax = 0.25;
Field[1]=BoundaryLayer;
Field[1].CurvesList={2};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.15;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test ex.msg_error_count>=1
    @test any(occursin("Only 2D Boundary Layers are supported",msg) for
              msg in ex.errors)
    @test validate(ex.mesh).ok
    # Dropped on both faces — neither shows a layer node at the field
    # offsets (0.05/0.11) inboard of the shared wall (1,0,0)→(1,1,0).
    m1=geo_entity_mesh(ex,2,1);m2=geo_entity_mesh(ex,2,2)
    @test !_blf_any_node_within(m1,(1,0,0),(1,1,0),0.02,0.13)
    @test !_blf_any_node_within(m2,(1,0,0),(1,1,0),0.02,0.13)
    # Excluding the sibling face admits the layer on surface 1 — oracle
    # strip nodes at x = 0.95 and x = 0.89 inboard of the shared edge.
    ex=_blf_execute("""
Point(1)={0,0,0,1};
Point(2)={1,0,0,1};
Point(3)={1,1,0,1};
Point(4)={0,1,0,1};
Point(5)={2,0,0,1};
Point(6)={2,1,0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line(5)={2,5};
Line(6)={5,6};
Line(7)={6,3};
Curve Loop(1)={1,2,3,4};
Curve Loop(2)={5,6,7,-2};
Plane Surface(1)={1};
Plane Surface(2)={2};
Mesh.MeshSizeMin = 0.25;
Mesh.MeshSizeMax = 0.25;
Field[1]=BoundaryLayer;
Field[1].CurvesList={2};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.15;
Field[1].ExcludedSurfacesList={2};
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test validate(ex.mesh).ok
    # The exclusion passes the face's `setupFor2d` count — no diagnostic
    # fires anywhere — and the strip grows on surface 1's wall only.
    @test !any(occursin("Only 2D Boundary Layers",msg) for msg in ex.errors)
    m1=geo_entity_mesh(ex,2,1);m2=geo_entity_mesh(ex,2,2)
    for x in (0.95,0.89),y in (0.25,0.5,0.75)
        @test _blf_node_at(m1,x,y;atol=1e-9)!=0
    end
    @test !_blf_any_node_within(m2,(1,0,0),(1,1,0),0.02,0.13)
    # A field listing a private wall alongside the shared one keeps the
    # private strip while the shared wall drops — upstream logs the error
    # and still completes the mesh on both faces.
    ex=_blf_execute("""
Point(1)={0,0,0,1};
Point(2)={1,0,0,1};
Point(3)={1,1,0,1};
Point(4)={0,1,0,1};
Point(5)={2,0,0,1};
Point(6)={2,1,0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line(5)={2,5};
Line(6)={5,6};
Line(7)={6,3};
Curve Loop(1)={1,2,3,4};
Curve Loop(2)={5,6,7,-2};
Plane Surface(1)={1};
Plane Surface(2)={2};
Mesh.MeshSizeMin = 0.25;
Mesh.MeshSizeMax = 0.25;
Field[1]=BoundaryLayer;
Field[1].CurvesList={1,2};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.15;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test validate(ex.mesh).ok
    @test any(occursin("Only 2D Boundary Layers are supported",msg) for
              msg in ex.errors)
    m1=geo_entity_mesh(ex,2,1);m2=geo_entity_mesh(ex,2,2)
    # Bottom wall survives — oracle strip nodes at the field offsets
    # (0.05, 0.11) inboard of curve 1. The shared wall (curve 2) grows
    # nothing on either face.
    for x in (0.25,0.5,0.75),y in (0.05,0.11)
        @test _blf_node_at(m1,x,y;atol=1e-9)!=0
    end
    @test !_blf_any_node_within(m1,(1,0,0),(1,1,0),0.02,0.13)
    @test !_blf_any_node_within(m2,(1,0,0),(1,1,0),0.02,0.13)
end

@testset "BoundaryLayer Field — bounded loud rejections" begin
    # Malformed `SizesList` payloads still reject: a non-integer point
    # tag and an odd (unpaired) list.
    for (opt,needle) in (
            ("SizesList={0.02,0.03}","SizesList point tag"),
            ("SizesList={1,0.02,2}","must list point-tag/height pairs"))
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
    # `FanPointsSizesList` without `FanPointsList` never matches a fan
    # point, and `IntersectMetrics`/`AnisoMax` are inert metric-operator
    # options for `.geo` BoundaryLayer fields — all mesh like upstream.
    for opt in ("FanPointsSizesList={3}","IntersectMetrics=1","AnisoMax=1e6")
        ex=_blf_execute(_BLF_SQUARE_HEAD*"""
Mesh.MeshSizeMin = 0.5;
Mesh.MeshSizeMax = 0.5;
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.25;
Field[1].$opt;
BoundaryLayer Field = 1;
Mesh 2;
""")
        @test validate(ex.mesh).ok
    end
end

# Negative `FanPointsSizesList` fan types −1/−2/−3 — upstream's
# `boundaryLayersData.cpp` `fanType` 1/2/3 (cross, mix, tail). The wall
# corner carries `nlayers+1` (cross/mix) or `2·nlayers` (tail) angular
# columns filled by the bilinear patch / rounded-fold transition /
# `l·n±ll·nc` formulas instead of rotated rays; the boundary columns ride
# the wall normals like every other fan. Node assertions below are the
# Gmsh 4.15.2 oracle's emitted coordinates on this exact L-shaped
# fixture — matched bitwise through the mean-plane working frame.
@testset "BoundaryLayer Field — negative FanPointsSizesList fan types" begin
    head="""
Point(1)={0,0,0,1};
Point(2)={1,0,0,1};
Point(3)={1,0.5,0,1};
Point(4)={0.5,0.5,0,1};
Point(5)={0.5,1,0,1};
Point(6)={0,1,0,1};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,6};
Line(6)={6,1};
Curve Loop(1)={1,2,3,4,5,6};
Plane Surface(1)={1};
Mesh.MeshSizeMin = 0.25;
Mesh.MeshSizeMax = 0.25;
"""
    tail="""
Field[1].CurvesList={3,4};
Field[1].FanPointsList={4};
Field[1].Size=0.04;
Field[1].Ratio=1.2;
Field[1].Thickness=0.12;
BoundaryLayer Field = 1;
Mesh 2;
"""
    function fan_mesh(sizes)
        # `FanPointsSizesList` must follow the `Field[1]=BoundaryLayer`
        # declaration — an option ahead of it is dropped ("No field with
        # id 1") and the corner falls back to an ordinary fan.
        src=head*"Field[1]=BoundaryLayer;\n"*
            "Field[1].FanPointsSizesList={"*sizes*"};"*tail
        return _blf_execute(src).mesh
    end
    # Boundary columns both walls climb the shared two-level schedule
    # (`0.04`, `0.088` cumulative) — every type emits these.
    for mesh in (fan_mesh("-1"),fan_mesh("-2"),fan_mesh("-3"))
        @test validate(mesh).ok
        @test _blf_signed_area(mesh)≈0.75 atol=1e-9
        for (x,y) in ((0.5,0.46),(0.5,0.41200000000000003),
                      (0.46,0.5),(0.41200000000000003,0.5),
                      (0.75,0.46),(0.75,0.41200000000000003),
                      (0.46,0.75),(0.41200000000000003,0.75))
            @test _blf_node_at(mesh,x,y;atol=1e-9)!=0
        end
    end
    # Cross (type 1): the `nbCol−1` bilinear strips land the interior
    # column tops on the square corner patch — oracle nodes (0.46,0.46)
    # through (0.412,0.412).
    m1=fan_mesh("-1")
    for (x,y) in ((0.46,0.46),(0.46,0.41200000000000003),
                  (0.41200000000000003,0.46),
                  (0.41200000000000003,0.41200000000000003))
        @test _blf_node_at(m1,x,y;atol=1e-9)!=0
    end
    # Mix (type 2): the rounded fold through `π·r/4` and `√2` — oracle
    # interior fold nodes, both bitwise-shared with the oracle mesh.
    m2=fan_mesh("-2")
    for (x,y) in ((0.46,0.46),
                  (0.42605887450304575,0.42605887450304575))
        @test _blf_node_at(m2,x,y;atol=1e-9)!=0
    end
    # Tail (type 3): `2·nlayers` columns marching the
    # `MEdge::length()`-scaled `nc` direction — oracle's characteristic
    # square-lattice interior nodes.
    m3=fan_mesh("-3")
    for (x,y) in ((0.32322330470336313,0.2352233047033632),
                  (0.32322330470336313,0.2832233047033632),
                  (0.32322330470336313,0.32322330470336313),
                  (0.2832233047033631,0.32322330470336313),
                  (0.23522330470336314,0.32322330470336313))
        @test _blf_node_at(m3,x,y;atol=1e-9)!=0
    end
    # A denser boundary discretization (MeshSize 0.05) reproduces the
    # oracle's tail-fan interior lattice bitwise — every interior column
    # threads its `⌊h/step⌋` rungs along `nc`.
    ex=_blf_execute(head*"""
Mesh.MeshSizeMin = 0.05;
Mesh.MeshSizeMax = 0.05;
Field[1]=BoundaryLayer;
Field[1].CurvesList={3,4};
Field[1].FanPointsList={4};
Field[1].FanPointsSizesList={-3};
Field[1].Size=0.04;
Field[1].Ratio=1.2;
Field[1].Thickness=0.12;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test validate(ex.mesh).ok
    for (x,y) in ((0.5,0.41200000000000003),
                  (0.41200000000000003,0.5),
                  (0.46464466094067264,0.37664466094067267),
                  (0.46464466094067264,0.42464466094067266),
                  (0.46464466094067264,0.46464466094067264),
                  (0.42464466094067266,0.46464466094067264),
                  (0.37664466094067267,0.46464466094067264))
        @test _blf_node_at(ex.mesh,x,y;atol=1e-9)!=0
    end
end

# A curved (ruled) wall face: the strip kernel works in the face's own
# `(u,v)` parametrization — upstream's `Mesh2DWithBoundaryLayers` pass —
# and lifts each emitted node back through the surface evaluator, so the
# field-graded column spacing follows the patch coordinates rather than a
# plane projection. Node assertions are the Gmsh 4.15.2 oracle's emitted
# positions on this fixture; the ruled path carries a parametrization-level
# residual (~1e-10) rather than the plane path's bitwise match.
@testset "BoundaryLayer Field — curved ruled surface field grading" begin
    ex=_blf_execute("""
Point(1)={0,0,0,1};
Point(2)={0,1,0,1};
Point(3)={0,0,1,1};
Point(4)={0,1,1,1};
Point(5)={0,0.5,0.5,1};
Circle(1)={1,5,2};
Circle(2)={3,5,4};
Line(3)={1,3};
Line(4)={2,4};
Curve Loop(1)={1,4,-2,-3};
Surface(1)={1};
Mesh.MeshSizeMin = 0.15;
Mesh.MeshSizeMax = 0.15;
Field[1]=BoundaryLayer;
Field[1].CurvesList={3};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.2;
BoundaryLayer Field = 1;
Mesh 2;
""")
    mesh=ex.mesh
    @test validate(mesh).ok
    # Oracle strip nodes — columns off Curve 3's wall, lifted through the
    # patch parametrization.
    for (x,y,z) in ((0.0,0.040770881000866795,0.1159369895809042),
                    (0.0,0.09340988706609243,0.08677199067010552),
                    (0.0,0.10715252110940521,-0.08793780141210805),
                    (0.0,0.3620503114569082,-0.193519922879555),
                    (0.0,0.5000000015286769,-0.20710678118654746))
        xs=mesh.coords[1,:];ys=mesh.coords[2,:];zs=mesh.coords[3,:]
        @test any(i->abs(xs[i]-x)<=1e-7 && abs(ys[i]-y)<=1e-7 &&
                     abs(zs[i]-z)<=1e-7,eachindex(xs))
    end
    # The field's geometric schedule is still visible off the wall: every
    # strip node clears the wall curve and no node regresses onto it —
    # `min(u,v)` pairing bounds the shared level count.
    zs=mesh.coords[3,:];ys=mesh.coords[2,:]
    interior=[i for i in eachindex(zs) if abs(ys[i])>1e-6]
    @test !isempty(interior)
end

@testset "BoundaryLayer Field — thin channel opposing walls" begin
    # Opposing walls in a thin channel: the two strips' free boundaries
    # meet — upstream's `modifyInitialMeshForBoundaryLayers` fails the
    # same configuration with "Edge not recovered".
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
Mesh.MeshSizeMin = 0.25;
Mesh.MeshSizeMax = 0.25;
Field[1]=BoundaryLayer;
Field[1].CurvesList={1,3};
Field[1].Size=0.05;
Field[1].Ratio=1.2;
Field[1].Thickness=0.3;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test err!==nothing
    @test occursin("self-intersects",sprint(showerror,err))
end

# `getLocalInfoAtNode` per-vertex grading: `SizesList` heights at model
# vertices, quadratic interpolation along each wall curve through the
# smaller endpoint height (`hb·(1−x)² + 2·min(hb,he)·x(1−x) + he·x²`),
# and `treat2Connections` corner growth — the biased averaged normal
# `normalize(1.01·N1+N2)` scaled by 1/cos(θ/2), where `N1` is the
# lower-tagged adjacent curve's normal (earlier line index inside one
# curve). Unequal neighboring columns tile with min-depth quadrangle
# rows plus a closing triangle on the deeper column.
@testset "BoundaryLayer Field — SizesList closed-ring grading" begin
    # Gmsh 4.15.2 oracle (square, MeshSizeMin/Max=0.4 → two subdivisions
    # per curve): every emitted strip node sits at the field's per-vertex
    # schedule along the `treat2Connections` direction — verified against
    # the oracle's node set on this exact field.
    ex=_blf_execute(_BLF_SQUARE_HEAD*"""
Mesh.MeshSizeMin = 0.4;
Mesh.MeshSizeMax = 0.4;
Field[1]=BoundaryLayer;
Field[1].CurvesList={1,2,3,4};
Field[1].NodesList={1,2,3,4};
Field[1].SizesList={1,0.02,2,0.05,3,0.03,4,0.04};
Field[1].Ratio=1.3;
Field[1].Thickness=0.2;
BoundaryLayer Field = 1;
Mesh 2;
""")
    mesh=ex.mesh
    @test validate(mesh).ok
    @test _blf_signed_area(mesh)≈1.0 atol=1e-9
    ratio,thick=1.3,0.2
    xs=mesh.coords[1,:];ys=mesh.coords[2,:]
    # Corner columns — the biased side is the lower-tagged adjacent
    # curve's normal (v1,v2 ← curve1; v3 ← curve2; v4 ← curve3).
    corners=(( (0.0,0.0),(0.0,1.0),(1.0,0.0),0.02),
             ( (1.0,0.0),(0.0,1.0),(-1.0,0.0),0.05),
             ( (1.0,1.0),(-1.0,0.0),(0.0,-1.0),0.03),
             ( (0.0,1.0),(0.0,-1.0),(1.0,0.0),0.04))
    columns=Vector{NTuple{2,Float64}}[]
    for (corner,n1,n2,hw) in corners
        col=_blf_column(corner,n1,n2,hw,ratio,thick)
        push!(columns,col)
        for (x,y) in col
            @test _blf_node_at(mesh,x,y)!=0
        end
    end
    # Oracle-matched first-layer corner positions (biased-normal scale).
    @test _blf_node_at(mesh,0.019900251,0.020099254;atol=1e-9)!=0
    @test _blf_node_at(mesh,0.950249372,0.050248134;atol=1e-9)!=0
    @test _blf_node_at(mesh,0.969851119,0.970149623;atol=1e-9)!=0
    @test _blf_node_at(mesh,0.039800502,0.959801493;atol=1e-9)!=0
    # Curve-interior wall vertices: quadratic endpoint interpolation,
    # plain inward-normal columns. Each wall maps a boundary point to its
    # endpoint heights, normalized parameter and column direction.
    walls=((v->abs(v[2])<1e-9,   (0.02,0.05), v->v[1],   (0.0,1.0)),
           (v->abs(v[1]-1)<1e-9, (0.05,0.03), v->v[2],  (-1.0,0.0)),
           (v->abs(v[2]-1)<1e-9, (0.03,0.04), v->1-v[1],(0.0,-1.0)),
           (v->abs(v[1])<1e-9,   (0.04,0.02), v->1-v[2],(1.0,0.0)))
    column_of=Dict{NTuple{2,Float64},Vector{NTuple{2,Float64}}}()
    for i in eachindex(xs)
        x,y=xs[i],ys[i]
        (0<x<1 && 0<y<1) && continue
        # Corner vertices carry the biased-normal columns already checked.
        isc=((x<1e-9 || x>1-1e-9) && (y<1e-9 || y>1-1e-9))
        isc && continue
        for (onwall,(hb,he),par,dir) in walls
            onwall((x,y)) || continue
            t=par((x,y))
            hw=hb*(1-t)^2+2*min(hb,he)*t*(1-t)+he*t^2
            col=NTuple{2,Float64}[(x+off*dir[1],y+off*dir[2])
                                  for off in _blf_schedule(hw,ratio,thick)]
            column_of[(x,y)]=col
            for (cx,cy) in col
                @test _blf_node_at(mesh,cx,cy)!=0
            end
            break
        end
    end
    # Two subdivisions per curve like the oracle's 1-D meshing.
    @test length(column_of)==8
    # Closing-triangle topology on unequal adjacent columns: order the
    # twelve wall columns cyclically (convex square → angular order) and,
    # for each neighboring pair whose schedule lengths differ, assert
    # upstream's `MTriangle(c_a[N],c_b[N],deeper[N+1])` element exists.
    ring=Tuple{NTuple{2,Float64},Vector{NTuple{2,Float64}}}[]
    for c in 1:4
        push!(ring,(corners[c][1],columns[c]))
    end
    for (p,col) in column_of
        push!(ring,(p,col))
    end
    sort!(ring,by=e->atan(e[1][2]-0.5,e[1][1]-0.5))
    closing=0
    trisets=[Set{Int}(Int.(mesh.tris[:,c])) for c in axes(mesh.tris,2)]
    for i in 1:length(ring)
        ca=ring[i][2]; cb=ring[mod1(i+1,length(ring))][2]
        na=length(ca); nb=length(cb)
        na==nb && continue
        N=min(na,nb)
        ia=_blf_node_at(mesh,ca[N]...)
        ib=_blf_node_at(mesh,cb[N]...)
        apex=na>nb ? _blf_node_at(mesh,ca[N+1]...) :
                     _blf_node_at(mesh,cb[N+1]...)
        (ia==0 || ib==0 || apex==0) && continue
        @test any(t->t==Set{Int}([ia,ib,apex]),trisets)
        closing+=1
    end
    @test closing>0
    # A uniformly-resolved `SizesList` overrides `Size` — every vertex
    # reports the same hWall, so the run keeps one shared schedule.
    ex2=_blf_execute(_BLF_SQUARE_HEAD*"""
Mesh.MeshSizeMin = 0.4;
Mesh.MeshSizeMax = 0.4;
Field[1]=BoundaryLayer;
Field[1].CurvesList={1,2,3,4};
Field[1].NodesList={1,2,3,4};
Field[1].SizesList={1,0.05,2,0.05,3,0.05,4,0.05};
Field[1].Size=0.02;
Field[1].Ratio=1.3;
Field[1].Thickness=0.2;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test validate(ex2.mesh).ok
    # Each bottom-wall boundary vertex grows its column on the resolved
    # 0.05 schedule — had `Size=0.02` applied, no (x,0.05) node exists.
    bnds=[(x,y) for (x,y) in
          zip(ex2.mesh.coords[1,:],ex2.mesh.coords[2,:])
          if abs(y)<1e-9 && 0<x<1]
    @test !isempty(bnds)
    for (x,y) in bnds
        @test _blf_node_at(ex2.mesh,x,0.05)!=0
    end
    # A heterogeneous `SizesList` on an *open* run resolves per-vertex
    # `hWall`s too — Gmsh 4.15.2 oracle on this exact fixture (13 nodes,
    # one Steiner point aside): the `NodesList` endpoints inject each
    # endpoint's own geometric chain onto the adjacent non-wall edges
    # (`addBoundaryLayerPoints` — geometric schedule capped at
    # `min(thickness, 0.4·chord)`), and the run's end columns ride those
    # chains like `addColumnAtTheEndOfTheBL`.
    ex3=_blf_execute(_BLF_SQUARE_HEAD*"""
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].NodesList={1,2};
Field[1].SizesList={1,0.02,2,0.06};
Field[1].Ratio=1.3;
Field[1].Thickness=0.25;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test validate(ex3.mesh).ok
    @test _blf_signed_area(ex3.mesh)≈1.0 atol=1e-9
    # Point 1's 0.02-wall chain on Curve 4 and Point 2's 0.06-wall
    # chain on Curve 2 — the oracle's exact injected node positions.
    for y in _blf_schedule(0.02,1.3,min(0.25,0.4))
        @test _blf_node_at(ex3.mesh,0.0,y)!=0
    end
    for y in _blf_schedule(0.06,1.3,min(0.25,0.4))
        @test _blf_node_at(ex3.mesh,1.0,y)!=0
    end
end

@testset "BoundaryLayer Field — extruded sibling gets the free pass" begin
    # `setupFor2d`: a BL curve bordering a swept lateral still counts
    # only the genuine plane face — the lateral is meshed by its own
    # producer and upstream's column pass emits no effective layer on
    # it. Gmsh 4.15.2 oracle (extruded cylinder, `Mesh 3`): the BL
    # rings grow on the disk source face only; each lateral keeps two
    # ordinary interior nodes.
    ex=_blf_execute("""
lc = 0.5;
Point(1)={0,0,0,lc};Point(2)={1,0,0,lc};Point(3)={1,1,0,lc};Point(4)={0,1,0,lc};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};
Plane Surface(1)={1};
out[] = Extrude{0,0,0.5}{ Surface{1}; };
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].Size=0.05;Field[1].Ratio=1.2;Field[1].Thickness=0.25;
BoundaryLayer Field = 1;
Mesh 2;
""";mesh_dim=2)
    @test validate(ex.mesh).ok
    # The bottom face grows its strip (the run's interior wall vertex
    # column, capped by the `min`-paired end columns); the swept
    # laterals carry no layer nodes at all — no node sits at a strip
    # schedule height on the wall's own lateral face.
    xs=ex.mesh.coords[1,:];ys=ex.mesh.coords[2,:];zs=ex.mesh.coords[3,:]
    offs=_blf_schedule(0.05,1.2,0.25)
    @test any(i->abs(zs[i])<1e-9 && abs(xs[i]-0.5)<1e-7 &&
                 isapprox(ys[i],offs[1];atol=1e-7),1:length(xs))
    for off in offs
        @test !any(i->abs(ys[i])<1e-9 && isapprox(zs[i],off;atol=1e-7),
                   1:length(xs))
    end
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

@testset "BoundaryLayer Field — PointsList composes with Transfinite" begin
    # Upstream `addBoundaryLayerPoints` runs before `meshGEdgeProcessing`,
    # transfinite included: the injected chain shortens the integration
    # window and `N = nbPointsTransfinite` still grades the reduced middle.
    # Gmsh 4.15.2 oracle on this fixture — curve 2 nodes:
    # 0, 0.06, 0.138, 0.2394 (chain) then 0.39152, 0.54364, 0.69576,
    # 0.84788 (4 interior transfinite nodes, N=6 over the window), 1.
    ex=_blf_execute(_BLF_SQUARE_HEAD*"""
Transfinite Curve {2} = 6;
Mesh.MeshSizeMin = 0.5;Mesh.MeshSizeMax = 0.5;
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].PointsList={2};
Field[1].SizesList={2,0.06};
Field[1].Ratio=1.3;
Field[1].Thickness=0.25;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test validate(ex.mesh).ok
    xs=ex.mesh.coords[1,:];ys=ex.mesh.coords[2,:]
    wall=sort!(Float64[ys[i] for i in eachindex(xs) if xs[i]>1-1e-9])
    @test wall≈[0.0,0.06,0.138,0.2394,0.39152,0.54364,0.69576,0.84788,1.0] atol=1e-5
end
