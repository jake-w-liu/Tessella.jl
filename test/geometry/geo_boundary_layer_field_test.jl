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
            ("SizesList={1,0.02,2}","SizesList"),
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

# `getLocalInfoAtNode` per-vertex grading: `SizesList` heights at model
# vertices, quadratic interpolation along each wall curve through the
# smaller endpoint height (`hb·(1−x)² + 2·min(hb,he)·x(1−x) + he·x²`),
# and `treat2Connections` corner growth — the biased averaged normal
# `normalize(1.01·N1+N2)` scaled by 1/cos(θ/2), where `N1` is the
# lower-tagged adjacent curve's normal (earlier line index inside one
# curve). Unequal neighboring columns tile with min-depth quadrangle
# rows plus a closing triangle on the deeper column.
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
    # Uniform `SizesList` on an *open* run works; heterogeneous grading
    # still rejects loudly — upstream's end columns consume the adjacent
    # edge's own chain, which needs the 1-D pass first.
    err=_blf_error(_BLF_SQUARE_HEAD*"""
Field[1]=BoundaryLayer;
Field[1].CurvesList={1};
Field[1].NodesList={1,2};
Field[1].SizesList={1,0.02,2,0.06};
Field[1].Ratio=1.3;
Field[1].Thickness=0.25;
BoundaryLayer Field = 1;
Mesh 2;
""")
    @test err!==nothing
    @test occursin("open wall runs",sprint(showerror,err))
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
