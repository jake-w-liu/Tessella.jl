# ── Order->=3 curved-cell Jacobian certification suite ────────────────────────
#
# Scope        : the publication-time `_ho_cell_jacobian_verdict` /
#                `_ho_constructed_cell_certify` gate in APICellBezier.jl and the
#                `validate(MixedMesh; certify_jacobians=true)` hook in Elements.
#                The certificate samples the Jacobian determinant on the
#                upstream `JacobianBasis::jacobianOrder` Lagrange points,
#                converts to Bernstein coefficients and tightens the control
#                bound by dyadic subdivision; a strictly one-sided bound
#                certifies the curved map for the element's own primary
#                orientation.
# Correctness  : INDEPENDENT oracle — textbook product-form order-3 Lagrange
#                bases (_sfac/_c4 families) evaluated straight from the public
#                `lagrange_nodes` reference coordinates; a dense sample lattice
#                gives an UPPER bound on the true determinant minimum, so
#                `verdict.lower <= sampled_min` is the no-wrong-certificate
#                property. Pyramid maps are rational (no polynomial oracle);
#                there the published `mesh_basis_functions` evaluator supplies
#                the sample bound — it still checks the certificate, not the
#                basis.
# Robustness   : affine, mildly warped, mirrored, folded, degenerate and
#                non-finite inputs across complete tri/qua/tet/hex/pri/pyr and
#                the serendipity qua/pri/pyr families; folded fixtures are
#                verified folded by the independent oracle (sampled min < 0 <
#                sampled max), so a refused certificate is never a vacuous
#                pass.
# Completeness : `validate` rejects a folded order-3 cell under the flag and
#                reports unsupported families explicitly; `set_order(3)` on a
#                real generated surface publishes only certified cells; a
#                mixed first-order block set round-trips through MSH write/read
#                with a stable CRC.
#
# Cost note: each family's determinant space is built once per process. Under
# --compile=min the tet/hex/pri/pyr builds cost roughly 30–100 s; the runner
# keeps one cell per family per case and shares the cached rules between the
# verdict, validation and generation sections.

using Test
using LinearAlgebra: det
using Random: MersenneTwister
import Tessella
import Tessella.Elements as E
const HOAPI = Tessella.API

_vertices(msh)=let ref=E.lagrange_nodes(msh)
    NTuple{3,Float64}[(ref[1,k],ref[2,k],ref[3,k]) for k in 1:size(ref,2)]
end
_verdict(v,msh)=HOAPI._ho_cell_jacobian_verdict(v,msh,"high-order Jacobian certificate test")

# ══ Independent order-3 nodal oracles ═══════════════════════════════════════
# Simplex direction factor phi_m(lambda; p) = prod_{l<m} (p*lambda - l)/(m - l)
# and its derivative; line cardinal factors on {-1,-1/3,1/3,1}. Written from
# the textbook definition, deliberately independent of the certificate's
# Bernstein machinery (and of Tessella's own basis tables).
_sfac(m,p,λ)=prod(l->(p*λ-l)/(m-l),0:m-1;init=one(λ))
function _sfacd(m,p,λ)
    m==0 && return zero(λ)
    # d/dλ prod_l f_l = sum_j f'_j prod_{l!=j} f_l with f'_j = p/(m-j).
    return p*sum(j->prod(l->l==j ? one(λ) : (p*λ-l)/(m-l),0:m-1)/(m-j),0:m-1)
end
const _C4PTS=(-1.0,-1.0/3,1.0/3,1.0)
_c4(i,x)=prod(j->j==i ? one(x) : (x-_C4PTS[j])/(_C4PTS[i]-_C4PTS[j]),1:4)
_c4d(i,x)=sum(k->k==i ? zero(x) :
    prod(j->(j==i||j==k) ? one(x) :
        (x-_C4PTS[j])/(_C4PTS[i]-_C4PTS[j]),1:4)/(_C4PTS[i]-_C4PTS[k]),1:4)
function _c4_index(q)
    return findmin(abs.([_C4PTS...].-q))[2]
end

# For each node index decode its basis slot: tri/tet barycentric triple plus,
# for the prism, the w cardinal index. Decoding from `lagrange_nodes` keeps the
# oracle independent of Tessella's ordering tables.
function _node_slot(ref,k,msh)
    if msh in (21,)
        return _tri_index((ref[1,k],ref[2,k]))
    elseif msh==29
        return _tet_index((ref[1,k],ref[2,k],ref[3,k]))
    elseif msh==36
        return (_c4_index(ref[1,k]),_c4_index(ref[2,k]))
    elseif msh==92
        return (_c4_index(ref[1,k]),_c4_index(ref[2,k]),_c4_index(ref[3,k]))
    elseif msh==90
        return (_tri_index((ref[1,k],ref[2,k])),_c4_index(ref[3,k]))
    end
    throw(ArgumentError("oracle family $msh not implemented"))
end
_tri_index(p)=(round(Int,3*(1-p[1]-p[2])),round(Int,3p[1]),round(Int,3p[2]))
_tet_index(p)=(round(Int,3*(1-p[1]-p[2]-p[3])),round(Int,3p[1]),
               round(Int,3p[2]),round(Int,3p[3]))

# Evaluate the tri/tet basis gradients at a DECODED multi-index slot.
function _tri_slot_grads(idx,r,s)
    λ=(1-r-s,r,s)
    F=ntuple(i->_sfac(idx[i],3,λ[i]),3)
    D=ntuple(i->_sfacd(idx[i],3,λ[i]),3)
    return (-D[1]*F[2]*F[3]+F[1]*D[2]*F[3],
            -D[1]*F[2]*F[3]+F[1]*F[2]*D[3],
            F[1]*F[2]*F[3])
end
function _tet_slot_grads(idx,r,s,t)
    λ=(1-r-s-t,r,s,t)
    F=ntuple(i->_sfac(idx[i],3,λ[i]),4)
    D=ntuple(i->_sfacd(idx[i],3,λ[i]),4)
    return (-D[1]*F[2]*F[3]*F[4]+F[1]*D[2]*F[3]*F[4],
            -D[1]*F[2]*F[3]*F[4]+F[1]*F[2]*D[3]*F[4],
            -D[1]*F[2]*F[3]*F[4]+F[1]*F[2]*F[3]*D[4])
end

function _oracle_jac(vertices,msh,p)
    ref=E.lagrange_nodes(msh);n=size(ref,2)
    J=zeros(3,3)
    if msh==21
        r,s=p
        for k in 1:n
            dr,ds,_=_tri_slot_grads(_node_slot(ref,k,msh),r,s)
            for d in 1:3
                J[d,1]+=dr*vertices[k][d];J[d,2]+=ds*vertices[k][d]
            end
        end
        return J[:,1],J[:,2]
    elseif msh==36
        u,v=p
        for k in 1:n
            i,j=_node_slot(ref,k,msh)
            du,dv=_c4d(i,u)*_c4(j,v),_c4(i,u)*_c4d(j,v)
            for d in 1:3
                J[d,1]+=du*vertices[k][d];J[d,2]+=dv*vertices[k][d]
            end
        end
        return J[:,1],J[:,2]
    elseif msh==29
        r,s,t=p
        for k in 1:n
            dr,ds,dt=_tet_slot_grads(_node_slot(ref,k,msh),r,s,t)
            for d in 1:3
                J[d,1]+=dr*vertices[k][d];J[d,2]+=ds*vertices[k][d];J[d,3]+=dt*vertices[k][d]
            end
        end
        return J[:,1],J[:,2],J[:,3]
    elseif msh==92
        u,v,w=p
        for k in 1:n
            i,j,l=_node_slot(ref,k,msh)
            du=_c4d(i,u)*_c4(j,v)*_c4(l,w)
            dv=_c4(i,u)*_c4d(j,v)*_c4(l,w)
            dw=_c4(i,u)*_c4(j,v)*_c4d(l,w)
            for d in 1:3
                J[d,1]+=du*vertices[k][d];J[d,2]+=dv*vertices[k][d];J[d,3]+=dw*vertices[k][d]
            end
        end
        return J[:,1],J[:,2],J[:,3]
    elseif msh==90
        r,s,w=p
        for k in 1:n
            (a,b,c),wi=_node_slot(ref,k,msh)
            dr,ds,tr=_tri_slot_grads((a,b,c),r,s)
            du,dv,dw=dr*_c4(wi,w),ds*_c4(wi,w),tr*_c4d(wi,w)
            for d in 1:3
                J[d,1]+=du*vertices[k][d];J[d,2]+=dv*vertices[k][d];J[d,3]+=dw*vertices[k][d]
            end
        end
        return J[:,1],J[:,2],J[:,3]
    end
    throw(ArgumentError("oracle family $msh not implemented"))
end

# Sampled range of the (signed) Jacobian determinant — the sampled minimum is
# an upper bound on the true minimum, so the certificate must never exceed it.
# For dim-2 fixtures in the z=0 plane the certificate's measure is the signed
# area J_z.
function _sampled_range(vertices,msh,points)
    lo,hi=Inf,-Inf
    for p in points
        cols=_oracle_jac(vertices,msh,p)
        value=length(cols)==2 ? cols[1][1]*cols[2][2]-cols[1][2]*cols[2][1] :
            det(hcat(cols...))
        lo=min(lo,value);hi=max(hi,value)
    end
    return lo,hi
end

# Published-evaluator sample range for the rational pyramid map: batched
# GradLagrange det evaluation on an interior lattice of the Bergot pyramid
# (base [-1,1]^2 at w=0 tapering to the apex at w=1 — the physical cross
# section at height w is [-w', w']^2 with w' = 1-w).
function _pyr_sampled_range(vertices,msh)
    ref=E.lagrange_nodes(msh);n=size(ref,2)
    flat=Float64[];count=0
    for w in 0.0:1/6:5/6,ur in -1.0:1/3:1.0,vr in -1.0:1/3:1.0
        append!(flat,(ur*(1-w),vr*(1-w),w));count+=1
    end
    _,grads,_=HOAPI.mesh_basis_functions(msh,flat,"GradLagrange";caller="pyr oracle")
    lo,hi=Inf,-Inf
    for sample in 1:count
        J=zeros(3,3)
        for node in 1:n,axis in 1:3,d in 1:3
            J[d,axis]+=grads[3(n*(sample-1)+node-1)+axis]*vertices[node][d]
        end
        value=det(J);lo=min(lo,value);hi=max(hi,value)
    end
    return lo,hi
end

_tri_pts(step)=((r,s) for r in 0.0:step:1.0 for s in 0.0:step:1.0-r)
_qua_pts(step)=((u,v) for u in -1.0:step:1.0 for v in -1.0:step:1.0)
_tet_pts(step)=((r,s,t) for r in 0.0:step:1.0 for s in 0.0:step:1.0-r
               for t in 0.0:step:1.0-r-s)
_hex_pts(step)=((u,v,w) for u in -1.0:step:1.0 for v in -1.0:step:1.0
               for w in -1.0:step:1.0)
_pri_pts(tstep,wstep)=((r,s,w) for (r,s) in _tri_pts(tstep) for w in -1.0:wstep:1.0)

@testset "High-order curved-cell Jacobian certificate" begin

    @testset "verdict classification" begin
        # Non-element inputs resolve before any determinant work.
        @test _verdict(_vertices(8),8).status==:unsupported      # dim-1 seg3
        @test _verdict(_vertices(15),15).status==:unsupported   # order-1 point
        bad=collect(_vertices(21));bad[4]=(NaN,0.0,0.0)
        @test _verdict(bad,21).status==:invalid
        @test occursin("non-finite",_verdict(bad,21).detail)
    end

    @testset "affine order-3 cells certify with det ~ 1" begin
        for msh in (21,36,29,92,90)
            verdict=_verdict(_vertices(msh),msh)
            @test verdict.status==:certified
            @test abs(verdict.lower-1.0)<1e-9
            @test abs(verdict.upper-1.0)<1e-9
        end
        # The pyramid certificate runs the order-8 collapsed-coordinate tensor
        # conversion; its conditioning makes ~1e-9 the honest noise floor even
        # though the represented determinant is exactly 1.
        verdict=_verdict(_vertices(118),118)
        @test verdict.status==:certified
        @test abs(verdict.lower-1.0)<1e-8
        @test abs(verdict.upper-1.0)<1e-8
        for msh in (16,18,19)    # serendipity qua/pri/pyr order-2
            verdict=_verdict(_vertices(msh),msh)
            @test verdict.status==:certified
            @test abs(verdict.lower-1.0)<1e-9
        end
    end

    @testset "curved-but-valid cells certify; bounds bracket the oracle range" begin
        rng=MersenneTwister(20261223)
        # Mild non-corner warps keep det > 0 but non-constant. Surface cells
        # stay in the z=0 plane so the oracle's signed-area measure matches the
        # certificate's primary-normal projection exactly. The tet fixture uses
        # a gentler amplitude: an order-3 tet's single interior slot plus all
        # edge/face slots moving ±0.04 folds the map, which the oracle
        # verifies would be a true fold, not a "valid curved" case.
        for (msh,pts,amp) in ((21,_tri_pts(1/12),0.04),(36,_qua_pts(0.25),0.04),
                              (29,_tet_pts(0.2),0.02),(92,_hex_pts(1/3),0.04),
                              (90,_pri_pts(1/8,1/3),0.04))
            ref=_vertices(msh);warped=collect(ref)
            dim=E.msh_spec(msh).dim
            for k in 5:length(warped)
                dz=dim==3 ? amp*(2rand(rng)-1) : 0.0
                warped[k]=warped[k].+(amp*(2rand(rng)-1),
                    amp*(2rand(rng)-1),dz)
            end
            verdict=_verdict(warped,msh)
            smin,smax=_sampled_range(warped,msh,pts)
            @test smin>0
            @test verdict.status==:certified
            @test verdict.lower<=smin+1e-9*max(1.0,abs(smin))
            @test verdict.upper>=smax-1e-9*max(1.0,abs(smax))
        end
        # Rational pyramid map: mild warp, bound checked against the published
        # evaluator's sampled range.
        ref=_vertices(118);warped=collect(ref)
        for k in 6:length(warped)
            warped[k]=warped[k].+(0.03*(2rand(rng)-1),
                0.03*(2rand(rng)-1),0.03*(2rand(rng)-1))
        end
        verdict=_verdict(warped,118)
        smin,smax=_pyr_sampled_range(warped,118)
        @test smin>0
        @test verdict.status==:certified
        @test verdict.lower<=smin+1e-9*max(1.0,abs(smin))
        @test verdict.upper>=smax-1e-9*max(1.0,abs(smax))
    end

    @testset "folded, mirrored and degenerate cells" begin
        # Each folded fixture is verified genuinely folded by the independent
        # oracle: sampled min < 0 < sampled max — the refused certificate is a
        # true positive, not bound slack.
        tri=collect(_vertices(21));tri[10]=(tri[10][1],tri[10][2]+2.0,tri[10][3])
        smin,smax=_sampled_range(tri,21,_tri_pts(1/10))
        @test smin<0<smax
        @test _verdict(tri,21).status==:invalid

        qua=collect(_vertices(36));qua[13]=(qua[13][1]+3.0,qua[13][2],qua[13][3])
        smin,smax=_sampled_range(qua,36,_qua_pts(0.125))
        @test smin<0<smax
        @test _verdict(qua,36).status==:invalid

        tet=collect(_vertices(29));tet[16]=(tet[16][1],tet[16][2],5.0)
        smin,smax=_sampled_range(tet,29,_tet_pts(0.2))
        @test smin<0<smax
        @test _verdict(tet,29).status==:invalid

        hex=collect(_vertices(92));hex[64]=(hex[64][1]+3.0,hex[64][2],hex[64][3])
        smin,smax=_sampled_range(hex,92,_hex_pts(1/3))
        @test smin<0<smax
        @test _verdict(hex,92).status==:invalid

        pri=collect(_vertices(90));pri[40]=(3.0,3.0,pri[40][3])
        smin,smax=_sampled_range(pri,90,_pri_pts(1/8,1/3))
        @test smin<0<smax
        @test _verdict(pri,90).status==:invalid

        pyr=collect(_vertices(118));pyr[30]=(pyr[30][1],pyr[30][2],pyr[30][3]-2.0)
        smin,smax=_pyr_sampled_range(pyr,118)
        @test smin<0<smax
        @test _verdict(pyr,118).status==:invalid

        # Mirrored affine maps carry orientation -1 and remain publishable.
        for msh in (29,92,90,118)
            mirrored=[(-x,y,z) for (x,y,z) in _vertices(msh)]
            @test _verdict(mirrored,msh).status==:certified
        end
        # A reversed surface cell is likewise a consistent immersion: the
        # z-normal component keeps a strictly negative bound.
        mirrored=[(-x,y,z) for (x,y,z) in _vertices(21)]
        @test _verdict(mirrored,21).status==:certified
        coplanar=[(x,y,0.0) for (x,y,z) in _vertices(29)]
        @test _verdict(coplanar,29).status==:degenerate
        # Serendipity folds: edge node overshoot into the element interior.
        qua8=collect(_vertices(16));qua8[5]=(qua8[5][1],qua8[5][2]+2.5,qua8[5][3])
        @test _verdict(qua8,16).status==:invalid
        pri15=collect(_vertices(18));pri15[7]=(pri15[7][1]+2.0,pri15[7][2]+0.0,
            pri15[7][3]+2.0)
        @test _verdict(pri15,18).status==:invalid
    end

    @testset "validate(MixedMesh; certify_jacobians) integration" begin
        ref_t=E.lagrange_nodes(21);ref_v=E.lagrange_nodes(29)
        coords=zeros(3,30)
        coords[:,1:10]=ref_t;coords[3,1:10].=0.0
        coords[:,11:30]=ref_v
        blocks=[E.ElementBlock(21,reshape(Int32.(collect(1:10)),10,1),Int32[1]),
                E.ElementBlock(29,reshape(Int32.(collect(11:30)),20,1),Int32[2])]
        mesh=E.MixedMesh(coords,blocks)
        @test E.validate(mesh).ok                       # structural only
        @test E.validate(mesh;certify_jacobians=true).ok
        @test_throws ArgumentError E.validate(mesh;certify_jacobians=1)

        folded=copy(coords);folded[3,16]=5.0
        bad=E.MixedMesh(folded,blocks)
        @test E.validate(bad).ok                        # still structurally fine
        result=E.validate(bad;certify_jacobians=true)
        @test !result.ok
        @test any(msg->occursin("Jacobian certificate",msg),result.messages)

        # Reversed-but-valid tet20 passes the flag.
        reversed=copy(coords);reversed[1,11:30].*=-1
        @test E.validate(E.MixedMesh(reversed,blocks);certify_jacobians=true).ok
    end

    @testset "constructed-cell gate throws on uncertified cells" begin
        ref=E.lagrange_nodes(29)
        nodes=[(point=(ref[1,k],ref[2,k],ref[3,k]),) for k in 1:20]
        good_cell=(msh=Int32(29),nodes=Int32.(collect(1:20)))
        @test HOAPI._ho_constructed_cell_certify(nodes,good_cell,"gate test")
        folded=copy(nodes);folded[16]=(point=(folded[16].point[1],folded[16].point[2],5.0),)
        @test_throws ArgumentError HOAPI._ho_constructed_cell_certify(
            folded,good_cell,"gate test")
    end

    @testset "generation publishes certified order-3 cells" begin
        A=Tessella.API
        A.finalize();A.initialize()
        try
            # Discrete surface: two tri3 lift to tri10 through the gate.
            A.model.add("ho cert")
            A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
            A.model.add_discrete_entity(2,7)
            A.mesh.add_nodes(2,7,[1,2,3,4],
                [0.0,0,0, 1,0,0, 0,1,0, 1,1,0])
            A.mesh.add_elements_by_type(7,2,[11,12],[1,2,3, 2,4,3])
            A.mesh.generate(2)
            A.mesh.set_order(3)
            types,ids,nodes=A.mesh.get_elements(2,7)
            @test types==[21]
            node_tags,flat,_=A.mesh.get_nodes(2,7)
            position=Dict(t=>i for (i,t) in enumerate(node_tags))
            for cell_nodes in nodes
                vertices=[Tuple(Float64.(flat[3position[t]-2:3position[t]]))
                          for t in cell_nodes]
                @test _verdict(vertices,21).status==:certified
            end
        finally
            A.finalize()
        end
    end

    @testset "mixed-family MSH transport round-trip" begin
        directory=mktempdir()
        # Unit-cube vertex set: cols 1-8, plus two helper midpoints.
        coordinates=Float64[
            0 1 0 0 1 1 0 1 0.5 0.5 0.5 0.5;
            0 0 1 0 0 1 1 1 0.5 0.5 0.5 0.5;
            0 0 0 1 0 0 1 1 0.0 0.5 0.5 1.0]
        blocks=Any[
            E.ElementBlock(2,Int32[1 2;2 5;3 3],Int32[1,2]),       # tri3
            E.ElementBlock(3,reshape(Int32[1,2,5,3],4,1),Int32[3]),  # qua4
            E.ElementBlock(4,reshape(Int32[1,2,3,4],4,1),Int32[4]),  # tet4
            E.ElementBlock(6,reshape(Int32[1,2,3,4,6,7],6,1),Int32[5]), # pri6
            E.ElementBlock(5,reshape(Int32[1,2,5,3,4,6,8,7],8,1),Int32[6]), # hex8
            E.ElementBlock(7,reshape(Int32[1,2,5,3,12],5,1),Int32[7]),    # pyr5
        ]
        mesh=E.MixedMesh(coordinates,blocks)
        @test E.validate(mesh).ok
        # MSH4 readers attach lossless entity records the source mesh lacks
        # and MSH2 readers attach elementary tags, so the full mixed_crc
        # legitimately changes across a generic round-trip. The transport
        # invariant is the structural digest: project the read-back mesh to
        # coords+blocks+physical_names, exactly as the existing MSH IO suite
        # does with legacy_crc.
        project(m)=E.MixedMesh(m.coords,m.blocks;physical_names=m.physical_names)
        source_crc=E.mixed_crc(project(mesh))
        for version in (2.2,4.1),binary in (false,true)
            path=joinpath(directory,"mixed-$version-$binary.msh")
            E.write_mixed_msh(path,mesh;version,binary)
            back=E.read_mixed_msh(path)
            @test E.validate(back).ok
            @test back.coords==coordinates
            @test E.mixed_crc(project(back)).sha==source_crc.sha
        end
    end
end
