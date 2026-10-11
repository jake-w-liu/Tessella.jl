# Native boundary-to-interior transfinite placement from InnerVertexPlacement.cpp.
# Coefficients accumulate in source order, then physical sums visit ascending columns.
function _native_inner_add!(row,indices,key,value)
    column=indices[key]
    row[column]=get(row,column,0.0)+value
    return nothing
end

function _native_inner_boundary_count(family,order,complete)
    family===:tri && return 3order
    family===:qua && return 4order
    family===:tet && return complete-max(0,(order-3)*(order-2)*(order-1)÷6)
    family===:hex && return complete-(order-1)^3
    family===:pri && return complete-(order-1)^2*(order-2)÷2
    family===:pyr && return complete-max(0,(order-2)*(order-1)*(2order-3)÷6)
    throw(ArgumentError("unsupported native interior placement family $family"))
end

function _native_inner_tri_row!(row,indices,coordinate,n)
    u,v=coordinate
    q=n - u - v
    xi=u / n
    eta=v / n
    rho=q / n
    _native_inner_add!(row,indices,(u + v, 0),xi)
    _native_inner_add!(row,indices,(u + q, v),xi)
    _native_inner_add!(row,indices,(0, v + u),eta)
    _native_inner_add!(row,indices,(u, v + q),eta)
    _native_inner_add!(row,indices,(0, v),rho)
    _native_inner_add!(row,indices,(u, 0),rho)
    _native_inner_add!(row,indices,(n, 0),-(xi))
    _native_inner_add!(row,indices,(0, n),-(eta))
    _native_inner_add!(row,indices,(0, 0),-(rho))
    return nothing
end

function _native_inner_qua_row!(row,indices,coordinate,n)
    u,v=coordinate
    xi=u / n
    eta=v / n
    _native_inner_add!(row,indices,(0, v),1 - xi)
    _native_inner_add!(row,indices,(n, v),xi)
    _native_inner_add!(row,indices,(u, 0),1 - eta)
    _native_inner_add!(row,indices,(u, n),eta)
    _native_inner_add!(row,indices,(0, 0),-((1 - xi) * (1 - eta)))
    _native_inner_add!(row,indices,(n, 0),-(xi * (1 - eta)))
    _native_inner_add!(row,indices,(n, n),-(xi * eta))
    _native_inner_add!(row,indices,(0, n),-((1 - xi) * eta))
    return nothing
end

function _native_inner_tet_row!(row,indices,coordinate,n)
    u,v,w=coordinate
    q=n - u - v - w
    xi=u / n
    eta=v / n
    zeta=w / n
    rho=q / n
    _native_inner_add!(row,indices,(u + v, 0, w),xi)
    _native_inner_add!(row,indices,(u + w, v, 0),xi)
    _native_inner_add!(row,indices,(u + q, v, w),xi)
    _native_inner_add!(row,indices,(0, v + u, w),eta)
    _native_inner_add!(row,indices,(u, v + w, 0),eta)
    _native_inner_add!(row,indices,(u, v + q, w),eta)
    _native_inner_add!(row,indices,(0, v, w + u),zeta)
    _native_inner_add!(row,indices,(u, 0, w + v),zeta)
    _native_inner_add!(row,indices,(u, v, w + q),zeta)
    _native_inner_add!(row,indices,(0, v, w),rho)
    _native_inner_add!(row,indices,(u, 0, w),rho)
    _native_inner_add!(row,indices,(u, v, 0),rho)
    _native_inner_add!(row,indices,(n - q, 0, 0),-(xi))
    _native_inner_add!(row,indices,(n - w, 0, w),-(xi))
    _native_inner_add!(row,indices,(n - v, v, 0),-(xi))
    _native_inner_add!(row,indices,(0, n - q, 0),-(eta))
    _native_inner_add!(row,indices,(0, n - w, w),-(eta))
    _native_inner_add!(row,indices,(u, n - u, 0),-(eta))
    _native_inner_add!(row,indices,(0, 0, n - q),-(zeta))
    _native_inner_add!(row,indices,(0, v, n - v),-(zeta))
    _native_inner_add!(row,indices,(u, 0, n - u),-(zeta))
    _native_inner_add!(row,indices,(0, 0, w),-(rho))
    _native_inner_add!(row,indices,(0, v, 0),-(rho))
    _native_inner_add!(row,indices,(u, 0, 0),-(rho))
    _native_inner_add!(row,indices,(n, 0, 0),xi)
    _native_inner_add!(row,indices,(0, n, 0),eta)
    _native_inner_add!(row,indices,(0, 0, n),zeta)
    _native_inner_add!(row,indices,(0, 0, 0),rho)
    return nothing
end

function _native_inner_hex_row!(row,indices,coordinate,n)
    u,v,w=coordinate
    xi=u / n
    eta=v / n
    zeta=w / n
    _native_inner_add!(row,indices,(0, v, w),1 - xi)
    _native_inner_add!(row,indices,(n, v, w),xi)
    _native_inner_add!(row,indices,(u, 0, w),1 - eta)
    _native_inner_add!(row,indices,(u, n, w),eta)
    _native_inner_add!(row,indices,(u, v, 0),1 - zeta)
    _native_inner_add!(row,indices,(u, v, n),zeta)
    _native_inner_add!(row,indices,(0, 0, w),-((1 - xi) * (1 - eta)))
    _native_inner_add!(row,indices,(0, n, w),-((1 - xi) * eta))
    _native_inner_add!(row,indices,(n, 0, w),-(xi * (1 - eta)))
    _native_inner_add!(row,indices,(n, n, w),-(xi * eta))
    _native_inner_add!(row,indices,(u, 0, 0),-((1 - eta) * (1 - zeta)))
    _native_inner_add!(row,indices,(u, 0, n),-((1 - eta) * zeta))
    _native_inner_add!(row,indices,(u, n, 0),-(eta * (1 - zeta)))
    _native_inner_add!(row,indices,(u, n, n),-(eta * zeta))
    _native_inner_add!(row,indices,(0, v, 0),-((1 - zeta) * (1 - xi)))
    _native_inner_add!(row,indices,(n, v, 0),-((1 - zeta) * xi))
    _native_inner_add!(row,indices,(0, v, n),-(zeta * (1 - xi)))
    _native_inner_add!(row,indices,(n, v, n),-(zeta * xi))
    _native_inner_add!(row,indices,(0, 0, 0),(1 - xi) * (1 - eta) * (1 - zeta))
    _native_inner_add!(row,indices,(0, 0, n),(1 - xi) * (1 - eta) * zeta)
    _native_inner_add!(row,indices,(0, n, 0),(1 - xi) * eta * (1 - zeta))
    _native_inner_add!(row,indices,(0, n, n),(1 - xi) * eta * zeta)
    _native_inner_add!(row,indices,(n, 0, 0),xi * (1 - eta) * (1 - zeta))
    _native_inner_add!(row,indices,(n, 0, n),xi * (1 - eta) * zeta)
    _native_inner_add!(row,indices,(n, n, 0),xi * eta * (1 - zeta))
    _native_inner_add!(row,indices,(n, n, n),xi * eta * zeta)
    return nothing
end

function _native_inner_pri_row!(row,indices,coordinate,n)
    u,v,w=coordinate
    q=n - u - v
    xi=u / n
    eta=v / n
    zeta=w / n
    rho=q / n
    _native_inner_add!(row,indices,(u + v, 0, w),xi)
    _native_inner_add!(row,indices,(u + q, v, w),xi)
    _native_inner_add!(row,indices,(0, v + u, w),eta)
    _native_inner_add!(row,indices,(u, v + q, w),eta)
    _native_inner_add!(row,indices,(0, v, w),rho)
    _native_inner_add!(row,indices,(u, 0, w),rho)
    _native_inner_add!(row,indices,(u, v, 0),1 - zeta)
    _native_inner_add!(row,indices,(u, v, n),zeta)
    _native_inner_add!(row,indices,(n, 0, w),-(xi))
    _native_inner_add!(row,indices,(0, n, w),-(eta))
    _native_inner_add!(row,indices,(0, 0, w),-(rho))
    _native_inner_add!(row,indices,(u + v, 0, 0),-(xi * (1 - zeta)))
    _native_inner_add!(row,indices,(u + q, v, 0),-(xi * (1 - zeta)))
    _native_inner_add!(row,indices,(0, v + u, 0),-(eta * (1 - zeta)))
    _native_inner_add!(row,indices,(u, v + q, 0),-(eta * (1 - zeta)))
    _native_inner_add!(row,indices,(0, v, 0),-(rho * (1 - zeta)))
    _native_inner_add!(row,indices,(u, 0, 0),-(rho * (1 - zeta)))
    _native_inner_add!(row,indices,(u + v, 0, n),-(xi * zeta))
    _native_inner_add!(row,indices,(u + q, v, n),-(xi * zeta))
    _native_inner_add!(row,indices,(0, v + u, n),-(eta * zeta))
    _native_inner_add!(row,indices,(u, v + q, n),-(eta * zeta))
    _native_inner_add!(row,indices,(0, v, n),-(rho * zeta))
    _native_inner_add!(row,indices,(u, 0, n),-(rho * zeta))
    _native_inner_add!(row,indices,(n, 0, 0),xi * (1 - zeta))
    _native_inner_add!(row,indices,(0, n, 0),eta * (1 - zeta))
    _native_inner_add!(row,indices,(0, 0, 0),rho * (1 - zeta))
    _native_inner_add!(row,indices,(n, 0, n),xi * zeta)
    _native_inner_add!(row,indices,(0, n, n),eta * zeta)
    _native_inner_add!(row,indices,(0, 0, n),rho * zeta)
    return nothing
end

function _native_inner_pyr_row!(row,indices,coordinate,n)
    u,v,height=coordinate
    w=n-height
    q=n - u - w
    r=n - v - w
    xi=u / n
    eta=v / n
    rho=q / n
    tau=r / n
    xip=u / (n - w)
    etap=v / (n - w)
    _native_inner_add!(row,indices,(0, v, w),1 - xip)
    _native_inner_add!(row,indices,(u, 0, w),1 - etap)
    _native_inner_add!(row,indices,(n - w, v, w),xip)
    _native_inner_add!(row,indices,(u, n - w, w),etap)
    _native_inner_add!(row,indices,(u, v, 0),rho * tau)
    _native_inner_add!(row,indices,(u + w, v, 0),xi * tau)
    _native_inner_add!(row,indices,(u + w, v + w, 0),xi * eta)
    _native_inner_add!(row,indices,(u, v + w, 0),rho * eta)
    _native_inner_add!(row,indices,(0, 0, w),-((1 - xip) * (1 - etap)))
    _native_inner_add!(row,indices,(n - w, 0, w),-(xip * (1 - etap)))
    _native_inner_add!(row,indices,(n - w, n - w, w),-(xip * etap))
    _native_inner_add!(row,indices,(0, n - w, w),-((1 - xip) * etap))
    _native_inner_add!(row,indices,(0, v, 0),-(rho * tau))
    _native_inner_add!(row,indices,(u, 0, 0),-(rho * tau))
    _native_inner_add!(row,indices,(u + w, 0, 0),-(xi * tau))
    _native_inner_add!(row,indices,(n, v, 0),-(xi * tau))
    _native_inner_add!(row,indices,(n, v + w, 0),-(xi * eta))
    _native_inner_add!(row,indices,(u + w, n, 0),-(xi * eta))
    _native_inner_add!(row,indices,(u, n, 0),-(rho * eta))
    _native_inner_add!(row,indices,(0, v + w, 0),-(rho * eta))
    _native_inner_add!(row,indices,(0, 0, 0),rho * tau)
    _native_inner_add!(row,indices,(n, 0, 0),xi * tau)
    _native_inner_add!(row,indices,(n, n, 0),xi * eta)
    _native_inner_add!(row,indices,(0, n, 0),rho * eta)
    return nothing
end

function _native_inner_rows(family,order)
    coordinates=getfield(Elements,Symbol("_",family,"_monomials"))(order,false)
    boundary=_native_inner_boundary_count(family,order,length(coordinates))
    length(coordinates)<=boundary && return Vector{Pair{Int,Float64}}[]
    indices=Dict((family===:pyr ? (q[1],q[2],order-q[3]) : q)=>i
        for (i,q) in enumerate(@view coordinates[1:boundary]))
    build=getfield(@__MODULE__,Symbol("_native_inner_",family,"_row!"))
    rows=Vector{Pair{Int,Float64}}[]
    for coordinate in @view coordinates[boundary+1:end]
        row=Dict{Int,Float64}()
        build(row,indices,coordinate,order)
        push!(rows,sort!(collect(row);by=first))
    end
    return rows
end
