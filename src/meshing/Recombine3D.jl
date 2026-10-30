"""
    Recombine3D

Hex-dominant recombination of a first-order tetrahedral volume complex —
Tessella's own implementation of the tet-to-hex combination algorithm
published in Pellerin et al., "Combining tetrahedra into hexahedra" (IMR
2017), matching the reachable public behaviour of upstream `meshCombine3D`:

  * candidate hexahedra are enumerated by a vertex-star walk over the tet
    complex, pruned by corner sine / scaled-Jacobian thresholds;
  * each candidate is filtered by the trilinear-hex Jacobian Bezier
    certificate (`hexExactValidity`), then its interior/boundary tetrahedra
    are recovered by a flood that stops on the cell's theoretical boundary
    triangulation;
  * a greedy pass ordered by the subdivision-based scaled-Jacobian bound
    selects a compatible subset (edges/diagonals/facet-triangle disjointness);
  * selected cells replace their tetrahedra, boundary quadrangle facets are
    transported onto the owning surface when both covering triangles share
    its colour, and interior hanging faces are made conformal by inserting
    type-140 trihedra.

Upstream gates this path on `Mesh.Algorithm3D=9` (ALGO_3D_RTREE) plus a
successful domhex `Filler3D` point-insertion pass and `Mesh.Recombine3DAll`
or `Recombine Volume{tag}`; the Filler insertion and the post-combine
`RelocateVertices` smoothing are not part of this port. Only the actually
reached cell-type arm is implemented — upstream's `selectCellsGreedy` call
site enables hexes and tets only (`{true,false,false,true}`), so prisms and
pyramids are never emitted upstream either, and `Mesh.Recombine3DLevel` /
`Mesh.Recombine3DConformity` have no consumers in the 4.15.2 sources.
"""
module Recombine3D

import ..MeshTypes: Mesh, nnodes, ntets
import ..Elements: MixedMesh, ElementBlock

export combine_volume_to_hexes, Combine3DResult

# ---------------------------------------------------------------------------
# Reference-cell topology.  Hex vertex numbering follows the MSH type-5
# layout (bottom quad 1-2-3-4, top quad 5-6-7-8, verticals 1-5, 2-6, 3-7,
# 4-8) — the numbering the upstream candidate walk emits.
# ---------------------------------------------------------------------------

# Hex faces in the MSH outward order (`MHexahedron::getFace` / faces_hexa) —
# the emitted boundary quadrangle vertex order.
const _HEX_FACES=((1,4,3,2),(1,2,6,5),(1,5,8,4),(2,3,7,6),(3,4,8,7),(5,6,7,8))
# The six quad-facet vertex sets in upstream `Hex::facetVertex` order — used
# for the theoretical-boundary slot grouping and boundary-tet detection.
const _HEX_FACETS=((1,2,3,4),(5,6,7,8),(1,2,6,5),(2,3,7,6),(4,3,7,8),(1,4,8,5))
const _HEX_EDGES=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),
                  (1,5),(2,6),(3,7),(4,8))
# Per-corner edge neighbours (`Hex::vertexAdjacentVertex`), for the corner
# scaled-Jacobian quality.
const _HEX_CORNER_NEIGHBORS=((2,4,5),(3,1,6),(4,2,7),(1,3,8),
                             (8,6,1),(5,7,2),(6,8,3),(7,5,4))
# The two diagonals of each quad facet (`Hex::quadFacetDiagonalVertex`),
# 12 pairs in facet order.
const _HEX_DIAGONALS=((1,3),(2,4),(5,7),(6,8),(1,6),(2,5),
                      (2,7),(3,6),(4,7),(3,8),(1,8),(4,5))
# The four canonical boundary triangles per quad facet
# (`Hex::quadFacetTriangleVertex`), 24 ordered triples in facet order.
const _HEX_BOUNDARY_TRIS=(
    (1,2,3),(1,3,4),(1,2,4),(2,3,4),
    (5,6,7),(5,7,8),(5,6,8),(6,7,8),
    (1,2,6),(1,6,5),(1,2,5),(2,6,5),
    (2,3,7),(2,7,6),(2,3,6),(3,7,6),
    (4,3,7),(4,7,8),(4,3,8),(3,7,8),
    (1,4,8),(1,8,5),(4,8,5),(4,5,1))
# Ordered cell-vertex triple -> boundary-triangle slot (1..24); every
# permutation maps to its slot, mirroring `vertexToBoundaryTriangle`.
const _HEX_BOUNDARY_SLOTS=let
    slots=Dict{NTuple{3,Int32},Int32}()
    for (slot,tri) in enumerate(_HEX_BOUNDARY_TRIS)
        a,b,c=tri
        for perm in ((a,b,c),(a,c,b),(b,a,c),(b,c,a),(c,a,b),(c,b,a))
            slots[perm]=Int32(slot)
        end
    end
    slots
end
# Facet index (1..6) of a boundary-triangle slot.
@inline _slot_facet(slot)=Int32((slot-1)÷4+1)
# Tet local faces: vertices of the face opposite local corner i.
const _TET_FACES=((2,3,4),(1,3,4),(1,2,4),(1,2,3))

@inline _pt3(coords,row::Integer)=(coords[1,row],coords[2,row],coords[3,row])
@inline _sub3(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _cross3(u,v)=(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
@inline _dot3(u,v)=u[1]*v[1]+u[2]*v[2]+u[3]*v[3]
@inline _norm3(u)=sqrt(_dot3(u,u))
# Signed tet volume x6, `pseudoTetVolume` convention: det[b-a,c-a,d-a].
@inline _orient3(a,b,c,d)=_dot3(_cross3(_sub3(b,a),_sub3(c,a)),_sub3(d,a))

@inline function _sorted_tri(a,b,c)
    a>b && ((a,b)=(b,a))
    b>c && ((b,c)=(c,b))
    a>b && ((a,b)=(b,a))
    return (a,b,c)
end
@inline _sorted_pair(a,b)=a<b ? (a,b) : (b,a)

# `quadCornerApproximateQuality`: sine of the corner angle — the quad-facet
# corner prune.
function _quad_corner_quality(coords,a,b,c)
    p0=_pt3(coords,a)
    u=_sub3(_pt3(coords,b),p0)
    v=_sub3(_pt3(coords,c),p0)
    uv=_norm3(u)*_norm3(v)
    w=_norm3(_cross3(u,v))
    uv<floatmin(Float64) && return w
    return w/uv
end

# `cellCornerQuality<Hex>`: scaled Jacobian at the corner — det of the three
# incident edges over their norm product.
function _hex_corner_quality(coords,v,ci,c1,c2,c3)
    p0=_pt3(coords,v[ci])
    e1=_sub3(_pt3(coords,v[c1]),p0)
    e2=_sub3(_pt3(coords,v[c2]),p0)
    e3=_sub3(_pt3(coords,v[c3]),p0)
    det=_dot3(_cross3(e2,e3),e1)
    det==0.0 && return 0.0
    norm=_norm3(e1)*_norm3(e2)*_norm3(e3)
    return norm>floatmin(Float64) ? det/norm : 0.0
end

# `cellApproximateQuality<Hex>`: minimum corner scaled Jacobian — the
# candidate-stage quality and threshold test.
function _hex_approx_quality(coords,v)
    q=Inf
    for i in 1:8
        n=_HEX_CORNER_NEIGHBORS[i]
        q=min(q,_hex_corner_quality(coords,v,i,n[1],n[2],n[3]))
    end
    return q
end

# ---------------------------------------------------------------------------
# Trilinear-hex Jacobian certificate and the subdivision quality bound —
# inexact (Float64) ports of `computeJacobian`/`computeBezierCoefficients`/
# `hexExactValidity`/`hexSubdivisionBasedQuality`.
# ---------------------------------------------------------------------------

# The 20 oriented pseudo-tet volumes: 8 corner Jacobians + 12 midpoint
# coefficients (tuples mark edge midpoints).
const _HEX_JAC20=((1,2,4,5),(2,3,1,6),(3,4,2,7),(4,1,3,8),
                  (5,8,6,1),(6,5,7,2),(7,6,8,3),(8,7,5,4),
                  (1,2,(3,4),(5,6)),(2,3,(4,1),(6,7)),(3,4,(1,2),(7,8)),(4,1,(2,3),(8,5)),
                  (1,5,(2,6),(4,8)),(2,6,(3,7),(1,5)),(3,7,(4,8),(2,6)),(4,8,(1,5),(3,7)),
                  (5,6,(1,2),(7,8)),(6,7,(2,3),(8,5)),(7,8,(3,4),(5,6)),(8,5,(4,1),(6,7)))

@inline _mid(p,q)=0.5 .* (p .+ q)
function _hex_jac20!(jac,pts)
    for i in 1:20
        a,b,c,d=_HEX_JAC20[i]
        pc=c isa Tuple ? _mid(pts[c[1]],pts[c[2]]) : pts[c]
        pd=d isa Tuple ? _mid(pts[d[1]],pts[d[2]]) : pts[d]
        jac[i]=_orient3(pts[a],pts[b],pc,pd)
        jac[i]<=0 && return false
    end
    return true
end

# The 27 Bezier control coefficients of the tri-quadratic Jacobian
# determinant (`computeBezierCoefficients`, inexact path).
function _hex_bezier27!(bez,jac)
    for i in 1:8; bez[i]=jac[i]; end
    for i in 9:12
        bez[i]=2*jac[i]-.5*(jac[i-8]+jac[(i-8)%4+1])
    end
    for i in 13:16
        bez[i]=2*jac[i]-.5*(jac[i-12]+jac[i-8])
    end
    for i in 17:20
        bez[i]=2*jac[i]-.5*(jac[i-12]+jac[5+(i-16)%4])
    end
    bez[21]=jac[9]+jac[10]+jac[11]+jac[12]-.75*(jac[1]+jac[2]+jac[3]+jac[4])
    for i in 22:25
        bez[i]=jac[i-13]+jac[i-9]+jac[13+(i-21)%4]+jac[i-5]-
            .75*(jac[i-21]+jac[i-17]+jac[(i-21)%4+1]+jac[5+(i-21)%4])
    end
    bez[26]=jac[17]+jac[18]+jac[19]+jac[20]-.75*(jac[5]+jac[6]+jac[7]+jac[8])
    s1=0.0;s2=0.0
    for k in 1:8; s1+=jac[k]; end
    for k in 9:20; s2+=jac[k]; end
    bez[27]=.5*(s2-1.25*s1)
    return bez
end

# Grid layout of the subdivided 5x5x5 coefficient field — the 27 control
# coefficients sit on the even-index lattice points (1-based).
const _BEZ_GRID=((1,1,1),(5,1,1),(5,5,1),(1,5,1),(1,1,5),(5,1,5),(5,5,5),(1,5,5),
                 (3,1,1),(5,3,1),(3,5,1),(1,3,1),(1,1,3),(5,1,3),(5,5,3),(1,5,3),
                 (3,1,5),(5,3,5),(3,5,5),(1,3,5),
                 (3,3,1),(3,1,3),(5,3,3),(3,5,3),(1,3,3),(3,3,5),(3,3,3))

# de Casteljau fill of the 5x5x5 grid from the 27 control points.
function _bezier_fill125!(g,bez)
    fill!(g,0.0)
    for i in 1:27
        x,y,z=_BEZ_GRID[i]
        g[x,y,z]=bez[i]
    end
    for i in 1:2:5, j in 1:2:5
        g[i,j,2]=.5*(g[i,j,1]+g[i,j,3])
        g[i,j,4]=.5*(g[i,j,3]+g[i,j,5])
        g[i,j,3]=.5*(g[i,j,2]+g[i,j,4])
    end
    for i in 1:2:5, k in 1:5
        g[i,2,k]=.5*(g[i,1,k]+g[i,3,k])
        g[i,4,k]=.5*(g[i,3,k]+g[i,5,k])
        g[i,3,k]=.5*(g[i,2,k]+g[i,4,k])
    end
    for j in 1:5, k in 1:5
        g[2,j,k]=.5*(g[1,j,k]+g[3,j,k])
        g[4,j,k]=.5*(g[3,j,k]+g[5,j,k])
        g[3,j,k]=.5*(g[2,j,k]+g[4,j,k])
    end
    return g
end

# Octant anchor offsets within the 5x5x5 grid and the 27 sub-coefficient
# offsets of `copyHexSubBezier`.
const _BEZ_OCTANTS=((0,0,0),(2,0,0),(2,2,0),(0,2,0),(0,0,2),(2,0,2),(2,2,2),(0,2,2))
const _BEZ_SUB27=((0,0,0),(2,0,0),(2,2,0),(0,2,0),(0,0,2),(2,0,2),(2,2,2),(0,2,2),
                  (1,0,0),(2,1,0),(1,2,0),(0,1,0),(0,0,1),(2,0,1),(2,2,1),(0,2,1),
                  (1,0,2),(2,1,2),(1,2,2),(0,1,2),
                  (1,1,0),(1,0,1),(2,1,1),(1,2,1),(0,1,1),(1,1,2),(1,1,1))

function _hex_valid_subdivide(bez,depth)
    depth>4 && return false
    g=_bezier_fill125!(Array{Float64,3}(undef,5,5,5),bez)
    for i in 1:2:5, j in 1:2:5, k in 1:2:5
        g[i,j,k]<=0 && return false
    end
    allpositive=true
    for i in 1:5, j in 1:5, k in 1:5
        g[i,j,k]<=0 && (allpositive=false;break)
    end
    allpositive && return true
    for n in 1:8
        i0,j0,k0=_BEZ_OCTANTS[n]
        sub=Vector{Float64}(undef,27)
        for (idx,(di,dj,dk)) in enumerate(_BEZ_SUB27)
            sub[idx]=g[i0+di+1,j0+dj+1,k0+dk+1]
        end
        _hex_valid_subdivide(sub,depth+1) || return false
    end
    return true
end

# `hexExactValidity` (inexact arm): all 20 pseudo-volumes positive, then the
# Bezier convex-hull certificate with bisection refinement.
function _hex_exact_validity(pts)
    jac=Vector{Float64}(undef,20)
    _hex_jac20!(jac,pts) || return false
    bez=Vector{Float64}(undef,27)
    _hex_bezier27!(bez,jac)
    minimum(bez)>0 && return true
    return _hex_valid_subdivide(bez,1)
end

# Partial derivatives of the trilinear map at the four free-parameter knots
# (`computePartialDerivatives`): dxdxi[component][direction][knot].
function _hex_partial_derivatives(pts)
    dx=Array{Float64,3}(undef,3,3,4)
    for i in 1:3
        dx[i,1,1]=pts[2][i]-pts[1][i]; dx[i,1,2]=pts[3][i]-pts[4][i]
        dx[i,1,3]=pts[7][i]-pts[8][i]; dx[i,1,4]=pts[6][i]-pts[5][i]
        dx[i,2,1]=pts[4][i]-pts[1][i]; dx[i,2,2]=pts[8][i]-pts[5][i]
        dx[i,2,3]=pts[7][i]-pts[6][i]; dx[i,2,4]=pts[3][i]-pts[2][i]
        dx[i,3,1]=pts[5][i]-pts[1][i]; dx[i,3,2]=pts[6][i]-pts[2][i]
        dx[i,3,3]=pts[7][i]-pts[3][i]; dx[i,3,4]=pts[8][i]-pts[4][i]
    end
    return dx
end

# `computeBezierCoeffDenominator` over the per-knot derivative norms Vj.
function _hex_bezier_denominator!(den,dx)
    v=Array{Float64,2}(undef,3,4)
    for i in 1:3, j in 1:4
        v[i,j]=sqrt(dx[1,i,j]^2+dx[2,i,j]^2+dx[3,i,j]^2)
    end
    den[1]=v[1,1]*v[2,1]*v[3,1]; den[2]=v[1,1]*v[2,4]*v[3,2]
    den[3]=v[1,2]*v[2,4]*v[3,3]; den[4]=v[1,2]*v[2,1]*v[3,4]
    den[5]=v[1,4]*v[2,2]*v[3,1]; den[6]=v[1,4]*v[2,3]*v[3,2]
    den[7]=v[1,3]*v[2,3]*v[3,3]; den[8]=v[1,3]*v[2,2]*v[3,4]
    den[9]=.5*(v[1,1]*v[2,1]*v[3,2]+v[1,1]*v[2,4]*v[3,1])
    den[10]=.5*(v[1,1]*v[2,4]*v[3,3]+v[1,2]*v[2,4]*v[3,2])
    den[11]=.5*(v[1,2]*v[2,1]*v[3,3]+v[1,2]*v[2,4]*v[3,4])
    den[12]=.5*(v[1,1]*v[2,1]*v[3,4]+v[1,2]*v[2,1]*v[3,1])
    den[13]=.5*(v[1,1]*v[2,2]*v[3,1]+v[1,4]*v[2,1]*v[3,1])
    den[14]=.5*(v[1,1]*v[2,3]*v[3,2]+v[1,4]*v[2,4]*v[3,2])
    den[15]=.5*(v[1,2]*v[2,3]*v[3,3]+v[1,3]*v[2,4]*v[3,3])
    den[16]=.5*(v[1,2]*v[2,2]*v[3,4]+v[1,3]*v[2,1]*v[3,4])
    den[17]=.5*(v[1,4]*v[2,2]*v[3,2]+v[1,4]*v[2,3]*v[3,1])
    den[18]=.5*(v[1,4]*v[2,3]*v[3,3]+v[1,3]*v[2,3]*v[3,2])
    den[19]=.5*(v[1,3]*v[2,3]*v[3,3]+v[1,3]*v[2,3]*v[3,4])
    den[20]=.5*(v[1,4]*v[2,2]*v[3,4]+v[1,3]*v[2,2]*v[3,1])
    den[21]=.25*(v[1,2]*v[2,4]*v[3,1]+v[1,2]*v[2,1]*v[3,2]+
                 v[1,1]*v[2,1]*v[3,3]+v[1,1]*v[2,4]*v[3,4])
    den[22]=.25*(v[1,4]*v[2,1]*v[3,2]+v[1,1]*v[2,2]*v[3,2]+
                 v[1,1]*v[2,3]*v[3,1]+v[1,4]*v[2,4]*v[3,1])
    den[23]=.25*(v[1,1]*v[2,3]*v[3,3]+v[1,2]*v[2,3]*v[3,2]+
                 v[1,3]*v[2,4]*v[3,2]+v[1,4]*v[2,4]*v[3,3])
    den[24]=.25*(v[1,3]*v[2,1]*v[3,3]+v[1,2]*v[2,2]*v[3,3]+
                 v[1,2]*v[2,3]*v[3,4]+v[1,3]*v[2,4]*v[3,4])
    den[25]=.25*(v[1,1]*v[2,2]*v[3,4]+v[1,2]*v[2,2]*v[3,1]+
                 v[1,3]*v[2,1]*v[3,1]+v[1,4]*v[2,1]*v[3,4])
    den[26]=.25*(v[1,3]*v[2,3]*v[3,1]+v[1,3]*v[2,2]*v[3,2]+
                 v[1,4]*v[2,2]*v[3,3]+v[1,4]*v[2,3]*v[3,4])
    den[27]=.125*(v[1,1]*v[2,2]*v[3,3]+v[1,1]*v[2,3]*v[3,4]+
                  v[1,2]*v[2,2]*v[3,2]+v[1,2]*v[2,3]*v[3,1]+
                  v[1,3]*v[2,1]*v[3,2]+v[1,3]*v[2,4]*v[3,1]+
                  v[1,4]*v[2,1]*v[3,3]+v[1,4]*v[2,4]*v[3,4])
    return den
end

# `computeBoundRational`: certified bound on num/den over the coefficient
# hull — zero on a strict contradiction.
function _bound_rational(num,den)
    ub=1.0;lb=0.0
    for i in 1:27
        if den[i]==0
            num[i]<0 && return 0.0
        elseif den[i]>0
            ub=min(ub,num[i]/den[i])
        else
            lb=max(lb,num[i]/den[i])
        end
    end
    return lb>ub ? 0.0 : ub
end

function _bounds_scaled_jacobian(num,dx)
    den=Vector{Float64}(undef,27)
    _hex_bezier_denominator!(den,dx)
    ub=1.0
    for i in 1:8
        den[i]>0 && (ub=min(ub,num[i]/den[i]))
        den[i]==0 && num[i]<0 && (ub=0.0)
    end
    return _bound_rational(num,den),ub
end

# The 9 de Casteljau knots of each bilinear derivative patch and the octant
# quadrant table of `copyHexSubDxdxi`.
const _DXI_NXI=((1,1,1),(1,4,2),(2,4,3),(2,1,4),(4,2,1),(4,3,2),(3,3,3),(3,2,4))
const _DXI_KNOTS=((1,5,9,8),(5,2,6,9),(9,6,3,7),(8,9,7,4))

function _dxdxi_fill9!(out,dx)
    for i in 1:3, j in 1:3
        out[i,j,1]=dx[i,j,1]; out[i,j,2]=dx[i,j,2]
        out[i,j,3]=dx[i,j,3]; out[i,j,4]=dx[i,j,4]
        out[i,j,5]=.5*(dx[i,j,1]+dx[i,j,2])
        out[i,j,6]=.5*(dx[i,j,2]+dx[i,j,3])
        out[i,j,7]=.5*(dx[i,j,4]+dx[i,j,1])
        out[i,j,8]=.5*(dx[i,j,3]+dx[i,j,4])
        out[i,j,9]=.5*(out[i,j,5]+out[i,j,7])
    end
    return out
end

# `subdivideAndHexQuality`: refine the (lb,ub) pair by octant bisection until
# the gap drops under 1e-3 or a true Jacobian value goes nonpositive.
function _quality_subdivide!(num,dx,depth,lb,ub)
    depth>10 && return lb,ub
    g=_bezier_fill125!(Array{Float64,3}(undef,5,5,5),num)
    for i in 1:2:5, j in 1:2:5, k in 1:2:5
        g[i,j,k]<=0 && return 0.0,0.0
    end
    allknots=_dxdxi_fill9!(Array{Float64,3}(undef,3,3,9),dx)
    lb=ub
    for n in 1:8
        i0,j0,k0=_BEZ_OCTANTS[n]
        sub=Vector{Float64}(undef,27)
        for (idx,(di,dj,dk)) in enumerate(_BEZ_SUB27)
            sub[idx]=g[i0+di+1,j0+dj+1,k0+dk+1]
        end
        subdx=Array{Float64,3}(undef,3,3,4)
        for j in 1:3
            quad=_DXI_KNOTS[_DXI_NXI[n][j]]
            for c in 1:4, i in 1:3
                subdx[i,j,c]=allknots[i,j,quad[c]]
            end
        end
        slb,subb=_bounds_scaled_jacobian(sub,subdx)
        ub=min(ub,subb)
        if ub-slb>1e-3
            slb,ub=_quality_subdivide!(sub,subdx,depth+1,slb,ub)
        end
        lb=min(lb,slb)
    end
    return lb,ub
end

# `hexSubdivisionBasedQuality`: midpoint of the certified scaled-Jacobian
# bound — the greedy ordering metric (zero for invalid cells).
function _hex_subdivision_quality(pts)
    jac=Vector{Float64}(undef,20)
    _hex_jac20!(jac,pts) || return 0.0
    num=Vector{Float64}(undef,27)
    _hex_bezier27!(num,jac)
    dx=_hex_partial_derivatives(pts)
    lb,ub=_bounds_scaled_jacobian(num,dx)
    ub-lb>1e-3 && ((lb,ub)=_quality_subdivide!(num,dx,1,lb,ub))
    return .5*(lb+ub)
end

# ---------------------------------------------------------------------------
# Working mesh context: adjacency, edge neighbourhoods, facet lookup and
# boundary colouring — `TetMeshForCombining` equivalents.
# ---------------------------------------------------------------------------

struct _RTetMesh
    coords::Matrix{Float64}
    tets::Matrix{Int32}
    vertex_tets::Vector{Vector{Int32}}
    vertex_neighbors::Vector{Vector{Int32}}            # descending
    facets::Set{NTuple{3,Int32}}
    adjacent::Matrix{Int32}                            # 4 x ntets, 0 = boundary
    boundary_tris::Set{NTuple{3,Int32}}                # adj==0 faces (old rows)
    face_owner::Dict{NTuple{3,Int32},Int32}            # boundary tri -> surface tag
    forbidden::Set{NTuple{3,Int32}}                    # quad-facet halves that must not pair
end

function _rtet_mesh(mesh::Mesh,boundary_owner,forbidden)
    nn=nnodes(mesh);nt=ntets(mesh);tets=mesh.tets
    vertex_tets=[Int32[] for _ in 1:nn]
    for t in 1:nt, i in 1:4
        push!(vertex_tets[tets[i,t]],Int32(t))
    end
    corner_neighbors=((2,3,4),(1,3,4),(1,2,4),(1,2,3))
    vertex_neighbors=[Int32[] for _ in 1:nn]
    for t in 1:nt, i in 1:4, j in corner_neighbors[i]
        push!(vertex_neighbors[tets[i,t]],tets[j,t])
    end
    for v in 1:nn
        sort!(unique!(vertex_neighbors[v]);rev=true)
    end
    facets=Set{NTuple{3,Int32}}();sizehint!(facets,4nt)
    boundary_tris=Set{NTuple{3,Int32}}();sizehint!(boundary_tris,4nt)
    adjacent=zeros(Int32,4,nt)
    firsthit=Dict{NTuple{3,Int32},Tuple{Int32,Int32}}();sizehint!(firsthit,2nt)
    for t in 1:nt, f in 1:4
        lv=_TET_FACES[f]
        face=_sorted_tri(tets[lv[1],t],tets[lv[2],t],tets[lv[3],t])
        push!(facets,face)
        prev=get(firsthit,face,nothing)
        if prev===nothing
            firsthit[face]=(Int32(t),Int32(f))
        else
            t2,f2=prev
            adjacent[f,t]=t2; adjacent[f2,t2]=Int32(t)
        end
    end
    for (face,(t,f)) in firsthit
        adjacent[f,t]==0 && push!(boundary_tris,face)
    end
    owner=boundary_owner===nothing ? Dict{NTuple{3,Int32},Int32}() : boundary_owner
    return _RTetMesh(mesh.coords,tets,vertex_tets,vertex_neighbors,facets,
                     adjacent,boundary_tris,owner,
                     forbidden===nothing ? Set{NTuple{3,Int32}}() : forbidden)
end

# Common edge-neighbours of the given vertices, in descending order.
function _around(m::_RTetMesh,vs::Tuple)
    out=collect(m.vertex_neighbors[vs[1]])
    for i in 2:length(vs)
        l=m.vertex_neighbors[vs[i]]
        filter!(x->x in l,out)
        isempty(out) && break
    end
    return out
end
@inline _around(m::_RTetMesh,v::Int32)=m.vertex_neighbors[Int(v)]

# Strictly-greater prefix of a descending list (`filterValues`).
@inline function _filter_above(list,floor::Int32)
    k=searchsortedfirst(list,floor;rev=true)-1
    return k<=0 ? Int32[] : view(list,1:k)
end

@inline _has_facet(m::_RTetMesh,a,b,c)=_sorted_tri(a,b,c) in m.facets
@inline function _is_triangulated_quad(m::_RTetMesh,a,b,c,d)
    (_has_facet(m,a,b,c) && _has_facet(m,a,c,d)) ||
    (_has_facet(m,b,c,d) && _has_facet(m,b,d,a))
end

# ---------------------------------------------------------------------------
# Candidate hex enumeration — `computeAllHex` (each candidate emitted once,
# with a the smallest vertex).
# ---------------------------------------------------------------------------

function _candidate_hexes(m::_RTetMesh,minq::Float64)
    coords=m.coords
    out=NTuple{8,Int32}[]
    for a in Int32(1):Int32(size(m.coords,2))
        aroundA=_filter_above(_around(m,a),a)
        for b in aroundA
            candidatesDE=_filter_above(aroundA,b)
            isempty(candidatesDE) && continue
            for d in candidatesDE
                _quad_corner_quality(coords,a,b,d)<minq && continue
                for e in candidatesDE
                    e==d && continue
                    _hex_corner_quality(coords,(a,b,d,e),1,2,3,4)<minq && continue
                    for c in _filter_above(_around(m,(b,d)),a)
                        c==e && continue
                        _quad_corner_quality(coords,b,a,c)<minq && continue
                        _quad_corner_quality(coords,c,b,d)<minq && continue
                        _quad_corner_quality(coords,d,a,c)<minq && continue
                        _is_triangulated_quad(m,a,b,c,d) || continue
                        for f in _filter_above(_around(m,(b,e)),a)
                            (f==c || f==d) && continue
                            _quad_corner_quality(coords,f,b,e)<minq && continue
                            _quad_corner_quality(coords,e,a,f)<minq && continue
                            _hex_corner_quality(coords,(b,c,a,f),1,2,3,4)<minq && continue
                            _is_triangulated_quad(m,a,b,f,e) || continue
                            for h in _filter_above(_around(m,(d,e)),a)
                                (h==b || h==c || h==f) && continue
                                _quad_corner_quality(coords,h,d,e)<minq && continue
                                _hex_corner_quality(coords,(d,a,c,h),1,2,3,4)<minq && continue
                                _hex_corner_quality(coords,(e,a,h,f),1,2,3,4)<minq && continue
                                _is_triangulated_quad(m,a,d,h,e) || continue
                                for g in _filter_above(_around(m,(c,f,h)),a)
                                    (g==b || g==d || g==e) && continue
                                    _hex_corner_quality(coords,(g,f,h,c),1,2,3,4)<minq && continue
                                    _hex_corner_quality(coords,(f,e,g,b),1,2,3,4)<minq && continue
                                    _hex_corner_quality(coords,(h,g,e,d),1,2,3,4)<minq && continue
                                    _hex_corner_quality(coords,(c,b,g,d),1,2,3,4)<minq && continue
                                    _is_triangulated_quad(m,d,c,g,h) || continue
                                    _is_triangulated_quad(m,e,f,g,h) || continue
                                    _is_triangulated_quad(m,b,c,g,f) || continue
                                    push!(out,(a,b,c,d,e,f,g,h))
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Candidate cell: interior/boundary tet recovery + legality checks
# (`computeOneInteriorTet`/`computeInteriorAndBoundaryTets`/
#  `computeRealCellBoundary` + the topology/colour/boundary audits).
# ---------------------------------------------------------------------------

struct _RCell
    vertices::NTuple{8,Int32}
    interior::Vector{Int32}
    boundary::Vector{Int32}
    boundary_facets::Vector{NTuple{2,Int32}}   # 12 packed (tet,localface) slots
end

@inline function _tet_facet(m::_RTetMesh,t,f)
    lv=_TET_FACES[f]
    _sorted_tri(m.tets[lv[1],t],m.tets[lv[2],t],m.tets[lv[3],t])
end
@inline _is_boundary_facet(m::_RTetMesh,t,f)=m.adjacent[f,t]==0
@inline function _cell_index(cellv,v)
    i=findfirst(==(v),cellv)
    return i===nothing ? Int32(0) : Int32(i)
end

# Tet whose four vertices coincide with one quad facet's vertex set —
# a sliver living on the cell surface (`isTetInBoundary`).
function _tet_in_boundary(m::_RTetMesh,t,cellv)
    tv=sort!(Int32[m.tets[i,t] for i in 1:4])
    for facet in _HEX_FACETS
        fv=sort!(Int32[cellv[facet[1]],cellv[facet[2]],cellv[facet[3]],cellv[facet[4]]])
        tv==fv && return true
    end
    return false
end

function _has4cellvertices(m::_RTetMesh,t,cellv)
    for i in 1:4
        m.tets[i,t] in cellv || return false
    end
    return true
end

# Boundary-triangle slot (1..24) if the face is on the theoretical boundary.
@inline function _theoretical_slot(m::_RTetMesh,t,f,cellv)
    lv=_TET_FACES[f]
    c1=_cell_index(cellv,m.tets[lv[1],t]);c1==0 && return Int32(0)
    c2=_cell_index(cellv,m.tets[lv[2],t]);c2==0 && return Int32(0)
    c3=_cell_index(cellv,m.tets[lv[3],t]);c3==0 && return Int32(0)
    return get(_HEX_BOUNDARY_SLOTS,(c1,c2,c3),Int32(0))
end

@inline function _opposite_vertex(m::_RTetMesh,t,f)
    used=_TET_FACES[f]
    for i in 1:4
        i in used || return m.tets[i,t]
    end
    return Int32(0)
end

# One seed tet known to sit inside the cell: either a tet on four cell
# vertices that is not a boundary sliver, or a tet glued to a theoretical-
# boundary triangle with an in-pointing apex (`isVolumePositive` — upstream
# `orient3d` sign convention).
function _one_interior_tet(m::_RTetMesh,cellv)
    for i in 1:8
        for t in m.vertex_tets[cellv[i]]
            if _has4cellvertices(m,t,cellv) && !_tet_in_boundary(m,t,cellv)
                return t
            end
        end
    end
    for tri in _HEX_BOUNDARY_TRIS
        v=(cellv[tri[1]],cellv[tri[2]],cellv[tri[3]])
        key=_sorted_tri(v[1],v[2],v[3])
        for t in m.vertex_tets[v[1]]
            for f in 1:4
                _tet_facet(m,t,f)==key || continue
                # isVolumePositive(v0,v1,v2,v3): the apex sits on the
                # negative side of the canonical triangle orientation.
                _orient3(_pt3(m.coords,v[1]),_pt3(m.coords,v[2]),
                         _pt3(m.coords,v[3]),
                         _pt3(m.coords,_opposite_vertex(m,t,f)))<0 ||
                    continue
                _tet_in_boundary(m,t,cellv) || return t
            end
        end
    end
    return Int32(0)
end

# Flood the cell interior from the seed without crossing theoretical-
# boundary triangles; recover the real boundary facets (24 slots packed per
# quad facet pair) and run the legality audits.
function _candidate_cell(m::_RTetMesh,cellv)
    seed=_one_interior_tet(m,cellv)
    seed==0 && return nothing
    interior=Int32[seed]
    boundary=Int32[]
    incell=falses(size(m.tets,2));incell[seed]=true
    stack=Int32[seed]
    while !isempty(stack)
        t=pop!(stack)
        for f in 1:4
            (length(interior)>15 || length(boundary)>5) && return nothing
            if !_is_boundary_facet(m,t,f)
                adj=m.adjacent[f,t]
                if !incell[adj]
                    if _theoretical_slot(m,t,f,cellv)==0
                        incell[adj]=true
                        push!(interior,adj);push!(stack,adj)
                    elseif _tet_in_boundary(m,adj,cellv)
                        incell[adj]=true
                        push!(boundary,adj)
                    end
                end
            else
                _theoretical_slot(m,t,f,cellv)==0 && return nothing
            end
        end
    end
    sort!(interior);sort!(boundary)
    nbtets=length(interior)+length(boundary)
    (nbtets==0 || !(length(interior)>4 && nbtets<16)) && return nothing
    # `computeRealCellBoundary`: every cell-tet face that is on the mesh
    # boundary or borders a non-cell tet must fill one theoretical slot —
    # packed pairwise per quad facet.
    slots=fill((Int32(0),Int32(0)),12)
    for t in Iterators.flatten((interior,boundary)), f in 1:4
        if !_is_boundary_facet(m,t,f)
            incell[m.adjacent[f,t]] && continue
        end
        slot=_theoretical_slot(m,t,f,cellv)
        slot==0 && return nothing
        base=2*_slot_facet(slot)-1
        if slots[base][1]==0
            slots[base]=(Int32(t),Int32(f))
        elseif slots[base+1][1]==0
            slots[base+1]=(Int32(t),Int32(f))
        else
            return nothing
        end
    end
    any(s->s[1]==0,slots) && return nothing
    # `checkQuadFacetBoundary` + `checkQuadFacetColor` + the forbidden
    # (embedded-sheet) half guard: the two covering tet faces of every quad
    # facet agree on boundary status and surface colour, and neither is a
    # declared interior sheet triangle.
    for f in 1:6
        t0,f0=slots[2*f-1];t1,f1=slots[2*f]
        _is_boundary_facet(m,t0,f0)==_is_boundary_facet(m,t1,f1) || return nothing
        _facet_color(m,t0,f0)==_facet_color(m,t1,f1) || return nothing
        isempty(m.forbidden) && continue
        (_tet_facet(m,t0,f0) in m.forbidden ||
         _tet_facet(m,t1,f1) in m.forbidden) && return nothing
    end
    return _RCell(cellv,interior,boundary,slots)
end

@inline _facet_color(m::_RTetMesh,t,f)=get(m.face_owner,_tet_facet(m,t,f),Int32(0))

# ---------------------------------------------------------------------------
# Greedy compatible selection — `computeCompatibleCellCombination`.
# ---------------------------------------------------------------------------

function _select_cells(cells::Vector{_RCell},m::_RTetMesh)
    nt=size(m.tets,2)
    selected_tets=falses(nt)
    selected_vertices=falses(size(m.coords,2))
    facet_tris=Set{NTuple{3,Int32}}()
    facet_diags=Set{NTuple{2,Int32}}()
    edges=Set{NTuple{2,Int32}}()
    qualities=[_hex_subdivision_quality(_cell_points(m,cell)) for cell in cells]
    order=sortperm(qualities;rev=true)
    selected=falses(length(cells))
    for i in order
        cell=cells[i]
        _cell_compatible(cell,selected_tets,selected_vertices,
                         facet_tris,facet_diags,edges) || continue
        selected[i]=true
        _cell_add_constraints!(cell,selected_tets,selected_vertices,
                               facet_tris,facet_diags,edges)
    end
    return selected
end

@inline _cell_points(m::_RTetMesh,cell)=
    NTuple{8,NTuple{3,Float64}}(_pt3(m.coords,v) for v in cell.vertices)

function _cell_compatible(cell,selected_tets,selected_vertices,
                          facet_tris,facet_diags,edges)
    any(t->selected_tets[t],cell.interior) && return false
    common=count(v->selected_vertices[v],cell.vertices)
    common<2 && return true
    v=cell.vertices
    for e in _HEX_EDGES
        _sorted_pair(v[e[1]],v[e[2]]) in facet_diags && return false
    end
    for d in _HEX_DIAGONALS
        _sorted_pair(v[d[1]],v[d[2]]) in edges && return false
    end
    for f in 1:6
        d0=_sorted_pair(v[_HEX_DIAGONALS[2*f-1][1]],v[_HEX_DIAGONALS[2*f-1][2]]) in facet_diags
        d1=_sorted_pair(v[_HEX_DIAGONALS[2*f][1]],v[_HEX_DIAGONALS[2*f][2]]) in facet_diags
        d0==d1 || return false
    end
    if common>2
        for f in 1:6
            nb=0
            for t in 1:4
                tri=_HEX_BOUNDARY_TRIS[4*(f-1)+t]
                _sorted_tri(v[tri[1]],v[tri[2]],v[tri[3]]) in facet_tris && (nb+=1)
            end
            (nb==0 || nb==4) || return false
        end
    end
    return true
end

function _cell_add_constraints!(cell,selected_tets,selected_vertices,
                                facet_tris,facet_diags,edges)
    v=cell.vertices
    for vert in v; selected_vertices[vert]=true; end
    for t in cell.interior; selected_tets[t]=true; end
    for tri in _HEX_BOUNDARY_TRIS
        push!(facet_tris,_sorted_tri(v[tri[1]],v[tri[2]],v[tri[3]]))
    end
    for d in _HEX_DIAGONALS
        push!(facet_diags,_sorted_pair(v[d[1]],v[d[2]]))
    end
    for e in _HEX_EDGES
        push!(edges,_sorted_pair(v[e[1]],v[e[2]]))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Emission: selected hexes + surviving tets + boundary-quad receipt + the
# trihedron conformalizer.
# ---------------------------------------------------------------------------

"""Result of `combine_volume_to_hexes` — see that function's docstring."""
struct Combine3DResult
    mesh::MixedMesh
    kept_rows::Vector{Int32}
    tet_source::Vector{Int32}
    quads::Vector{NTuple{4,Int32}}
    quad_surfaces::Vector{Int32}
    halves::Vector{NTuple{2,NTuple{3,Int32}}}
    unpaired_boundary_quads::Int
end

"""
    combine_volume_to_hexes(mesh::Mesh, boundary_owner=nothing;
                            forbidden_tris=nothing, min_quality=1e-4)
        -> Union{Nothing,Combine3DResult}

`boundary_owner` maps each sorted boundary-triangle triple to its owning
surface tag (upstream's per-GFace triangle colours; `nothing` leaves every
boundary face uncoloured — candidates still combine but no surface
quadrangle receipt can be emitted). `forbidden_tris` lists triangle triples
(interior sheet facets) a quad facet may never cover.

Returns `nothing` when no hexahedron is selected; otherwise the combined
mixed mesh plus the receipts the caller needs: `kept_rows` (old -> new node
rows), `tet_source` (source tet row per emitted tet column), and the surface
quadrangle substitution list `quads`/`quad_surfaces`/`halves` (outward face
order and the two consumed covering triangles per quad). Boundary quad
facets whose halves do not share a surface stay triangulated (counted in
`unpaired_boundary_quads` — upstream's "wrong mesh topology" warning arm).
"""
function combine_volume_to_hexes(mesh::Mesh,boundary_owner=nothing;
                                 forbidden_tris=nothing,
                                 min_quality::Float64=1e-4)
    m=_rtet_mesh(mesh,boundary_owner,forbidden_tris)
    coords=m.coords
    cells=_RCell[]
    for v in _candidate_hexes(m,min_quality)
        _hex_approx_quality(coords,v)>min_quality || continue
        pts=NTuple{8,NTuple{3,Float64}}(_pt3(coords,x) for x in v)
        _hex_exact_validity(pts) || continue
        cell=_candidate_cell(m,v)
        cell===nothing || push!(cells,cell)
    end
    isempty(cells) && return nothing
    selected=_select_cells(cells,m)
    any(selected) || return nothing
    # --- volume element emission (old index space) --------------------------
    hexes=Matrix{Int32}(undef,8,0)
    surface_quads=NTuple{4,Int32}[]
    quad_surfaces=Int32[]
    surface_halves=NTuple{2,NTuple{3,Int32}}[]
    deleted_tris=Set{NTuple{3,Int32}}()
    used_vertices=falses(size(coords,2))
    consumed=Set{Int32}()
    unpaired=0
    for i in eachindex(cells)
        selected[i] || continue
        v=cells[i].vertices
        hexes=hcat(hexes,collect(v))
        for vert in v; used_vertices[vert]=true; end
        for face in _HEX_FACES
            f=(v[face[1]],v[face[2]],v[face[3]],v[face[4]])
            emitted=false
            # diagonal (w0,w2) split first, then (w1,w3) — upstream order.
            for (h1,h2) in ((_sorted_tri(f[1],f[2],f[3]),_sorted_tri(f[3],f[4],f[1])),
                            (_sorted_tri(f[1],f[2],f[4]),_sorted_tri(f[4],f[2],f[3])))
                s1=get(m.face_owner,h1,Int32(0));s2=get(m.face_owner,h2,Int32(0))
                s1!=0 && s1==s2 || continue
                push!(surface_quads,f);push!(quad_surfaces,s1)
                push!(surface_halves,(h1,h2))
                push!(deleted_tris,h1);push!(deleted_tris,h2)
                emitted=true;break
            end
            if !emitted
                # Quad facet on the domain boundary whose covering halves
                # cannot be paired on one surface — upstream's "wrong mesh
                # topology" warning arm; the triangles stay.
                for h in (_sorted_tri(f[1],f[2],f[3]),_sorted_tri(f[3],f[4],f[1]),
                          _sorted_tri(f[1],f[2],f[4]),_sorted_tri(f[4],f[2],f[3]))
                    h in m.boundary_tris && (unpaired+=1;break)
                end
            end
        end
        union!(consumed,cells[i].interior)
        union!(consumed,cells[i].boundary)
    end
    tets_out=Matrix{Int32}(undef,4,0)
    kept_tets=Int32[]
    for t in 1:size(m.tets,2)
        t in consumed && continue
        push!(kept_tets,Int32(t))
        tets_out=hcat(tets_out,m.tets[:,t])
        for i in 1:4; used_vertices[m.tets[i,t]]=true; end
    end
    # --- odd-incidence faces + trihedron conformalizer (old index space) ----
    face_count=Dict{Tuple{Vararg{Int32}},Tuple{Int,Tuple{Vararg{Int32}}}}()
    function bump(face)
        key=Tuple(sort!(collect(face)))
        entry=get(face_count,key,(0,face))
        face_count[key]=(entry[1]+1,entry[2])
    end
    for col in eachcol(hexes),face in _HEX_FACES
        bump((col[face[1]],col[face[2]],col[face[3]],col[face[4]]))
    end
    for col in eachcol(tets_out),face in _TET_FACES
        bump((col[face[1]],col[face[2]],col[face[3]]))
    end
    substituted=Set{NTuple{4,Int32}}(Tuple(sort!(collect(q))) for q in surface_quads)
    hanging=Tuple{Tuple{Vararg{Int32}},Tuple{Vararg{Int32}}}[]
    hanging_tris=Set{Tuple{Vararg{Int32}}}()
    for (face,(count,witness)) in face_count
        count==1 || continue
        if length(face)==4
            face in substituted && continue
            f=_sorted_tri(face[1],face[2],face[3])
            # Domain-boundary quads are covered by their boundary halves —
            # never a hanging interior face.
            if f in m.boundary_tris || _sorted_tri(face[1],face[2],face[4]) in m.boundary_tris ||
               _sorted_tri(face[1],face[3],face[4]) in m.boundary_tris ||
               _sorted_tri(face[2],face[3],face[4]) in m.boundary_tris
                continue
            end
            push!(hanging,(face,witness))
        else
            face in m.boundary_tris && continue
            face in deleted_tris && continue
            push!(hanging_tris,face)
        end
    end
    trihs=Matrix{Int32}(undef,4,0)
    for (face,witness) in sort!(hanging;by=first)
        w=witness
        length(w)==4 || continue
        t1=_sorted_tri(w[1],w[2],w[3]);t2=_sorted_tri(w[3],w[4],w[1])
        if t1 in hanging_tris && t2 in hanging_tris
            trihs=hcat(trihs,Int32[w[4],w[1],w[2],w[3]])
            continue
        end
        t1=_sorted_tri(w[1],w[2],w[4]);t2=_sorted_tri(w[4],w[2],w[3])
        if t1 in hanging_tris && t2 in hanging_tris
            trihs=hcat(trihs,Int32[w[1],w[2],w[3],w[4]])
        end
    end
    for col in eachcol(trihs)
        for i in 1:4; used_vertices[col[i]]=true; end
    end
    # --- vertex cull + remap --------------------------------------------------
    kept_rows=findall(used_vertices)
    remap=zeros(Int32,size(coords,2))
    for (new,old) in enumerate(kept_rows); remap[old]=Int32(new); end
    remap_col!(A)=begin
        for i in eachindex(A); A[i]=remap[A[i]]; end; A
    end
    blocks=ElementBlock[]
    if size(tets_out,2)>0
        push!(blocks,ElementBlock(4,remap_col!(tets_out),zeros(Int32,size(tets_out,2))))
    end
    if size(hexes,2)>0
        push!(blocks,ElementBlock(5,remap_col!(hexes),zeros(Int32,size(hexes,2))))
    end
    if size(trihs,2)>0
        push!(blocks,ElementBlock(140,remap_col!(trihs),zeros(Int32,size(trihs,2))))
    end
    quads=[ntuple(i->remap[q[i]],4) for q in surface_quads]
    halves=[( _sorted_tri(remap[h[1][1]],remap[h[1][2]],remap[h[1][3]]),
              _sorted_tri(remap[h[2][1]],remap[h[2][2]],remap[h[2][3]]))
            for h in surface_halves]
    return Combine3DResult(MixedMesh(coords[:,kept_rows],blocks),kept_rows,
                         kept_tets,quads,quad_surfaces,halves,unpaired)
end

end # module
