# Certify newly constructed full P2 maps before publishing an API cache. The
# interval path bounds actual derivative Bernstein coefficients over the entire
# reference element; an inconclusive bound uses HighOrder's exact dyadic
# Bernstein coefficients. Reversed volume connectivity keeps its source sign.

const _API_P2_SEG_DERIVATIVES = (Int16[-3 1; -1 3; 4 -4],)
const _API_P2_TRI_DERIVATIVES = (
    Int16[HighOrder._P2TRI_DR[k][q] for k in 1:6, q in 1:3],
    Int16[HighOrder._P2TRI_DS[k][q] for k in 1:6, q in 1:3])
const _API_P2_TET_DERIVATIVES = ntuple(3) do d
    Int16[HighOrder._REF_VERTEX_GRADS[q][k][d] for k in 1:10, q in 1:4]
end
const _API_P2_QUAD_DERIVATIVES = ntuple(2) do d
    weights=zeros(Int16,9,6)
    for k in 1:9, i in 0:(d==1 ? 1 : 2), j in 0:(d==2 ? 1 : 2)
        a,b=HighOrder._P2QUAD_NODES[k]
        weights[k,i*(d==2 ? 2 : 3)+j+1]=d==1 ?
            HighOrder._P2L2_DER[a][i+1]*HighOrder._P2L2_VAL[b][j+1] :
            HighOrder._P2L2_VAL[a][i+1]*HighOrder._P2L2_DER[b][j+1]
    end
    weights
end
const _API_P2_HEX_DERIVATIVES = ntuple(3) do d
    weights=zeros(Int16,27,18)
    for n in 1:27, i in 0:(d==1 ? 1 : 2),
        j in 0:(d==2 ? 1 : 2), k in 0:(d==3 ? 1 : 2)
        a,b,c=HighOrder._P2HEX_NODES[n]
        column=(i*(d==2 ? 2 : 3)+j)*(d==3 ? 2 : 3)+k+1
        first=d==1 ? HighOrder._P2L2_DER[a][i+1] : HighOrder._P2L2_VAL[a][i+1]
        second=d==2 ? HighOrder._P2L2_DER[b][j+1] : HighOrder._P2L2_VAL[b][j+1]
        third=d==3 ? HighOrder._P2L2_DER[c][k+1] : HighOrder._P2L2_VAL[c][k+1]
        weights[n,column]=first*second*third
    end
    weights
end
const _API_P2_PRISM_DERIVATIVES = ntuple(3) do d
    weights=zeros(Int16,18,d==3 ? 12 : 9)
    for n in 1:18
        a,bw=HighOrder._P2PRI_NODES[n]
        if d==3
            for g in 1:6, b in 0:1
                weights[n,(g-1)*2+b+1]=
                    HighOrder._P2TRI_V2[a][g]*HighOrder._P2L2_DER[bw][b+1]
            end
        else
            for v in 1:3, b in 0:2
                coefficient=d==1 ? HighOrder._P2TRI_DR[a][v] : HighOrder._P2TRI_DS[a][v]
                weights[n,(v-1)*3+b+1]=coefficient*HighOrder._P2L2_VAL[bw][b+1]
            end
        end
    end
    weights
end

_api_p2_derivative_weights(::Val{8})=_API_P2_SEG_DERIVATIVES
_api_p2_derivative_weights(::Val{9})=_API_P2_TRI_DERIVATIVES
_api_p2_derivative_weights(::Val{10})=_API_P2_QUAD_DERIVATIVES
_api_p2_derivative_weights(::Val{11})=_API_P2_TET_DERIVATIVES
_api_p2_derivative_weights(::Val{12})=_API_P2_HEX_DERIVATIVES
_api_p2_derivative_weights(::Val{13})=_API_P2_PRISM_DERIVATIVES
_api_p2_node_count(::Val{8})=Val(3)
_api_p2_node_count(::Val{9})=Val(6)
_api_p2_node_count(::Val{10})=Val(9)
_api_p2_node_count(::Val{11})=Val(10)
_api_p2_node_count(::Val{12})=Val(27)
_api_p2_node_count(::Val{13})=Val(18)
_api_p2_node_count(::Val{14})=Val(14)

function _api_p2_derivative_box(differences,weights)
    return ntuple(Val(3)) do axis
        lo=Inf; hi=-Inf
        for coefficient in axes(weights,2)
            value=HighOrder._P2_JAC_ZERO
            for node in axes(weights,1)
                factor=weights[node,coefficient]
                factor==0 && continue
                term=HighOrder._p2_jac_interval_scale(differences[node][axis],Float64(factor))
                value=HighOrder._p2_jac_interval_add(value,term)
            end
            lo=min(lo,value.lo);hi=max(hi,value.hi)
        end
        HighOrder._P2JacInterval(lo,hi)
    end
end

@inline _api_p2_sign(value)=value>0 ? 1 : value<0 ? -1 : 0
_api_p2_primary_orientation(v,::Val{11})=_api_p2_sign(-Model.orient3(v[1],v[2],v[3],v[4]))
_api_p2_primary_orientation(v,::Val{12})=_api_p2_sign(-Model.orient3(v[1],v[2],v[4],v[5]))
_api_p2_primary_orientation(v,::Val{13})=_api_p2_sign(-Model.orient3(v[1],v[2],v[3],v[4]))
_api_p2_primary_orientation(v,::Val{14})=_api_p2_sign(-Model.orient3(v[1],v[2],v[4],v[5]))
_api_p2_primary_orientation(v,::Val)=1

@inline _api_p2_excludes_zero(value)=value.lo>0 || value.hi<0
function _api_p2_interval_positive(v::NTuple{N,NTuple{3,Float64}},kind,orientation) where N
    differences=ntuple(Val(N)) do node
        ntuple(axis->HighOrder._p2_jac_interval_difference(v[node][axis],v[1][axis]),Val(3))
    end
    weights=_api_p2_derivative_weights(kind)
    derivatives=ntuple(i->_api_p2_derivative_box(differences,weights[i]),length(weights))
    if length(derivatives)==1
        return any(_api_p2_excludes_zero,derivatives[1])
    elseif length(derivatives)==2
        a,b=derivatives
        cross=(HighOrder._p2_jac_interval_sub(
                   HighOrder._p2_jac_interval_mul(a[2],b[3]),HighOrder._p2_jac_interval_mul(a[3],b[2])),
               HighOrder._p2_jac_interval_sub(
                   HighOrder._p2_jac_interval_mul(a[3],b[1]),HighOrder._p2_jac_interval_mul(a[1],b[3])),
               HighOrder._p2_jac_interval_sub(
                   HighOrder._p2_jac_interval_mul(a[1],b[2]),HighOrder._p2_jac_interval_mul(a[2],b[1])))
        return any(_api_p2_excludes_zero,cross)
    end
    determinant=HighOrder._p2_jac_interval_det(derivatives...)
    return orientation==1 ? determinant.lo>0 : determinant.hi<0
end

_api_p2_exact_coefficients(coords,cells,column,::Val{8})=HighOrder._p2seg_bernstein_coeffs(coords,cells,column)[1]
_api_p2_exact_coefficients(coords,cells,column,::Val{9})=HighOrder._p2tri_bernstein_coeffs(coords,cells,column)[1]
_api_p2_exact_coefficients(coords,cells,column,::Val{10})=HighOrder._p2quad_bernstein_coeffs(coords,cells,column)[1]
_api_p2_exact_coefficients(coords,cells,column,::Val{11})=HighOrder._p2_bernstein_coeffs(coords,cells,column)[1]
_api_p2_exact_coefficients(coords,cells,column,::Val{12})=HighOrder._p2hex_bernstein_coeffs(coords,cells,column)[1]
_api_p2_exact_coefficients(coords,cells,column,::Val{13})=HighOrder._p2pri_bernstein_coeffs(coords,cells,column)[1]

function _api_p2_cell_certify(v::NTuple{N,NTuple{3,Float64}},kind::Val{MSH},caller;
                              coords=nothing,cells=nothing,column=1) where {N,MSH}
    all(point->all(isfinite,point),v) || throw(ArgumentError(
        "$caller: newly constructed MSH $MSH has non-finite P2 coordinates"))
    MSH==8 && v[1]==v[2] && throw(ArgumentError(
        "$caller: primary line endpoints coincide below Float64 coordinate resolution"))
    orientation=_api_p2_primary_orientation(v,kind)
    orientation!=0 || throw(ArgumentError(
        "$caller: newly constructed MSH $MSH has degenerate primary reference orientation"))
    if MSH==14
        HighOrder._p2_pyramid_jacobian_certify(v,orientation,caller)
        return nothing
    end
    _api_p2_interval_positive(v,kind,orientation) && return nothing
    if coords===nothing
        coords=Matrix{Float64}(undef,3,N)
        cells=Matrix{Int32}(undef,N,1)
        for node in 1:N
            cells[node,1]=Int32(node)
            for axis in 1:3
                coords[axis,node]=v[node][axis]
            end
        end
        column=1
    end
    coefficients=_api_p2_exact_coefficients(coords,cells,column,kind)
    certified=orientation==1 ? all(>(0),coefficients) : all(<(0),coefficients)
    certified || throw(ArgumentError(
        "$caller: newly constructed MSH $MSH lacks a whole-domain P2 Jacobian certificate"))
    return nothing
end

function _api_p2_block_certify(coords,cells,kind::Val{MSH},caller) where MSH
    for column in axes(cells,2)
        vertices=ntuple(node->ntuple(axis->coords[axis,cells[node,column]],Val(3)),_api_p2_node_count(kind))
        _api_p2_cell_certify(vertices,kind,caller;coords,cells,column)
    end
    return nothing
end

function _api_p2_block_certify(coords,cells,msh::Int,caller)
    msh==15 && return nothing
    msh in 8:14 || throw(ArgumentError("$caller: unsupported full P2 MSH type $msh"))
    return _api_p2_block_certify(coords,cells,Val(msh),caller)
end

function _api_p2_constructed_mesh_certify(mesh,caller)
    for block in mesh.blocks
        _api_p2_block_certify(mesh.coords,block.nodes,Int(block.msh),caller)
    end
    return nothing
end

function _api_p2_constructed_cell_certify(nodes,cell,kind,caller)
    vertices=ntuple(i->nodes[cell.nodes[i]].point,_api_p2_node_count(kind))
    return _api_p2_cell_certify(vertices,kind,caller)
end
_api_p2_constructed_cell_certify(nodes,cell,caller)=
    _api_p2_constructed_cell_certify(nodes,cell,Val(cell.msh),caller)

_api_p2_overlay_certify(p::P2TriMesh,caller)=_api_p2_block_certify(p.coords,p.tri6,Val(9),caller)
_api_p2_overlay_certify(p::P2Mesh,caller)=_api_p2_block_certify(p.coords,p.tet10,Val(11),caller)
