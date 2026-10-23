# Public reference queries follow polynomialBasis.cpp's floating monomial
# interpolation. The stable FEM nodal engine remains independent.
struct _NativeCellPolynomial
    exponents::Vector{NTuple{3,Int}}
    coefficients::Matrix{Float64}
    dimension::Int
end

function _native_cell_edge_patterns(spec)
    primary=Elements._msh_primary_nodes(spec.family)
    edges=_native_order_edges(spec.family)
    p=spec.order
    return Tuple((a+1,b+1,(primary+(index-1)*(p-1)+1:primary+index*(p-1))...)
                 for (index,(a,b)) in enumerate(edges))
end

function _native_cell_face_patterns(spec,face_type,caller)
    spec.family===:pri && spec.serendipity && spec.order>2 && throw(ArgumentError(
        "$caller: Gmsh 4.15.2 full incomplete Prism face closure accesses absent interior nodes"))
    spec.dim==2 && return (Tuple(1:spec.nnodes),)
    primary=Elements._msh_primary_nodes(spec.family)
    edges=_native_order_edges(spec.family)
    edge_width=spec.order-1
    interior=primary+length(edges)*edge_width
    result=Tuple[]
    for face in _native_order_faces(spec.family)
        corners=Int[index+1 for index in face if index>=0]
        face_width=length(corners)
        inner_width=spec.serendipity ? 0 : face_width==3 ?
            (spec.order-1)*(spec.order-2)÷2 : (spec.order-1)^2
        if face_width==face_type
            nodes=copy(corners)
            for index in 1:face_width
                first=corners[index]-1;last=corners[mod1(index+1,face_width)]-1
                edge_index=findfirst(edge->edge==(first,last) || edge==(last,first),edges)
                edge_index===nothing && throw(ArgumentError("$caller: invalid native face incidence"))
                supports=collect(primary+(edge_index-1)*edge_width+1:primary+edge_index*edge_width)
                edges[edge_index]==(first,last) || reverse!(supports)
                append!(nodes,supports)
            end
            append!(nodes,interior+1:interior+inner_width)
            push!(result,Tuple(nodes))
        end
        interior+=inner_width
    end
    return Tuple(result)
end

function _native_cell_reference_nodes(spec)
    points=getfield(Elements,Symbol("_",spec.family,"_monomials"))(spec.order,spec.serendipity)
    result=zeros(Float64,3,length(points))
    p=spec.order
    if spec.family===:pyr
        for (column,point) in enumerate(points)
            w=1.0-Float64(point[3])/p
            result[1,column]=-1.0+w+Float64(point[1])*2.0/p
            result[2,column]=-1.0+w+Float64(point[2])*2.0/p
            result[3,column]=w
        end
        return result
    end
    factor=spec.family in (:tri,:tet) ? 1.0/p : 2.0/p
    for (column,point) in enumerate(points),axis in 1:spec.dim
        result[axis,column]=Float64(point[axis])*factor
        spec.family in (:qua,:hex) && (result[axis,column]+=-1.0)
    end
    return result
end

# Eigen's dynamic PartialPivLU: recursively factor tall panels, apply their
# recorded row interchanges, solve A12 and update the trailing Schur block.
@static if Sys.ARCH in (:x86_64,:i686)
    function _native_cell_cpuid(leaf::UInt32,subleaf::UInt32)
        try
            return Base.llvmcall(raw"""
            %r = call {i32,i32,i32,i32} asm sideeffect "cpuid", "={ax},={bx},={cx},={dx},0,2"(i32 %0, i32 %1)
            %a = extractvalue {i32,i32,i32,i32} %r, 0
            %b = extractvalue {i32,i32,i32,i32} %r, 1
            %c = extractvalue {i32,i32,i32,i32} %r, 2
            %d = extractvalue {i32,i32,i32,i32} %r, 3
            %t0 = insertvalue [4 x i32] undef, i32 %a, 0
            %t1 = insertvalue [4 x i32] %t0, i32 %b, 1
            %t2 = insertvalue [4 x i32] %t1, i32 %c, 2
            %t3 = insertvalue [4 x i32] %t2, i32 %d, 3
            ret [4 x i32] %t3
            """,NTuple{4,UInt32},Tuple{UInt32,UInt32},leaf,subleaf)
        catch err
            err isa ErrorException || rethrow()
            # `llvmcall` needs the compiler — interpreted minimal-JIT runs
            # fall through to the documented default cache sizes.
            return (UInt32(0),UInt32(0),UInt32(0),UInt32(0))
        end
    end
end

function _native_cell_detect_cache_sizes()
    @static if Sys.ARCH in (:x86_64,:i686)
        vendor=_native_cell_cpuid(UInt32(0),UInt32(0))
        l1=0;l2=0;l3=0
        if (vendor[2],vendor[4],vendor[3]) in
            ((0x68747541,0x69746e65,0x444d4163),(0x69444d41,0x74656273,0x21726574))
            maximum=_native_cell_cpuid(0x80000000,UInt32(0))[1]
            if maximum>=0x80000006
                l1=Int(_native_cell_cpuid(0x80000005,UInt32(0))[3]>>24)*1024
                descriptor=_native_cell_cpuid(0x80000006,UInt32(0))
                l2=Int(descriptor[3]>>16)*1024
                l3=Int((descriptor[4]&0x0fffc000)>>18)*512*1024
            end
        elseif vendor[1]>=4
            for index in UInt32(0):UInt32(15)
                a,b,c,_=_native_cell_cpuid(UInt32(4),index)
                kind=a&0x0f
                kind==0 && break
                if kind==1 || kind==3
                    level=(a&0xe0)>>5
                    width=Int((b&0xffc00000)>>22)+1
                    partitions=Int((b&0x003ff000)>>12)+1
                    line=Int(b&0x00000fff)+1
                    capacity=width*partitions*line*(Int(c)+1)
                    level==1 && (l1=capacity)
                    level==2 && (l2=capacity)
                    level==3 && (l3=capacity)
                end
            end
        end
        return (l1>0 ? l1 : 32*1024,l2>0 ? l2 : 256*1024,l3>0 ? l3 : 2*1024*1024)
    else
        return (16*1024,512*1024,512*1024)
    end
end

# Eigen captures this hardware configuration once per process. Re-querying it
# on every element changes floating association after migration between cores
# with different L1 capacities. The snapshot contains three integers only.
const _NATIVE_CELL_CACHE_SIZES=Base.OncePerProcess{NTuple{3,Int}}(_native_cell_detect_cache_sizes)
_native_cell_cache_sizes()=_NATIVE_CELL_CACHE_SIZES()

# SSE2 mr=nr=4, with Eigen's cache heuristic and balanced final blocks.
function _native_cell_blocking(depth,rows,columns,factor,cache)
    max(depth,rows,columns)<48 && return depth,rows,columns
    l1,l2,l3=cache
    max_kc=max(((l1-128)÷(factor*64))&~7,1)
    old_depth=depth
    if depth>max_kc
        remainder=depth%max_kc
        depth=remainder==0 ? max_kc : max_kc-8*((max_kc-1-remainder)÷(8*(depth÷max_kc+1)))
    end
    remaining_l1=l1-128-rows*depth*8
    max_nc=remaining_l1>=32depth ? remaining_l1÷(8depth) : (3*1572864)÷(4max_kc*8)
    nc=min(1572864÷(2depth*8),max_nc)&~3
    if columns>nc
        remainder=columns%nc
        columns=remainder==0 ? nc : nc-4*((nc-remainder)÷(4*(columns÷nc+1)))
    elseif old_depth==depth
        problem_size=depth*columns*8
        actual_lm=1572864;max_mc=rows
        if problem_size<=1024
            actual_lm=l1
        elseif l3!=0 && problem_size<=32768
            actual_lm=l2;max_mc=min(576,max_mc)
        end
        mc=min(actual_lm÷(3depth*8),max_mc)
        mc>4 ? (mc-=mc%4) : mc==0 && return depth,rows,columns
        remainder=rows%mc
        rows=remainder==0 ? mc : mc-4*((mc-remainder)÷(4*(rows÷mc+1)))
    end
    return depth,rows,columns
end

function _native_cell_product!(target,left,right)
    rows,columns=size(target);depth=size(left,2)
    tail=(rows÷4)*4+1;full_columns=(columns÷4)*4;peeled=(depth÷8)*8
    @inbounds for column in 1:columns,row in 1:rows
        value=0.0
        paired=rows%4>=2 && tail<=row<tail+2 && column<=full_columns
        if paired
            even=0.0;odd=0.0
            for k in 1:2:peeled
                even+=left[row,k]*right[k,column]
                odd+=left[row,k+1]*right[k+1,column]
            end
            value=even+odd
        end
        for k in (paired ? peeled+1 : 1):depth
            value+=left[row,k]*right[k,column]
        end
        target[row,column]+=(-1.0)*value
    end
    return target
end

function _native_cell_lower_solve!(tri,rhs,cache=_native_cell_cache_sizes())
    isempty(rhs) && return rhs
    n=size(tri,1);kc,mc,_=_native_cell_blocking(n,n,size(rhs,2),4,cache)
    for block in 1:kc:n
        block_stop=min(n,block+kc-1)
        for start in block:4:block_stop
            stop=min(block_stop,start+3)
        @inbounds for i in start:stop,column in axes(rhs,2)
            value=rhs[i,column]
            for row in i+1:stop
                rhs[row,column]-=value*tri[row,i]
            end
        end
            stop<block_stop && _native_cell_product!(@view(rhs[stop+1:block_stop,:]),
                @view(tri[stop+1:block_stop,start:stop]),@view(rhs[start:stop,:]))
        end
        for first in block_stop+1:mc:n
            last=min(n,first+mc-1)
            _native_cell_product!(@view(rhs[first:last,:]),@view(tri[first:last,block:block_stop]),
                @view(rhs[block:block_stop,:]))
        end
    end
    return rhs
end

function _native_cell_upper_solve!(tri,rhs,cache=_native_cell_cache_sizes())
    isempty(rhs) && return rhs
    n=size(tri,1);kc,mc,_=_native_cell_blocking(n,n,size(rhs,2),4,cache)
    for block_stop in n:-kc:1
        block=max(1,block_stop-kc+1)
        for stop in block_stop:-4:block
            start=max(block,stop-3)
        @inbounds for i in stop:-1:start
            reciprocal=1.0/tri[i,i]
            for column in axes(rhs,2)
                rhs[i,column]*=reciprocal
                value=rhs[i,column]
                for row in start:i-1
                    rhs[row,column]-=value*tri[row,i]
                end
            end
        end
            start>block && _native_cell_product!(@view(rhs[block:start-1,:]),
                @view(tri[block:start-1,start:stop]),@view(rhs[start:stop,:]))
        end
        for first in 1:mc:block-1
            last=min(block-1,first+mc-1)
            _native_cell_product!(@view(rhs[first:last,:]),@view(tri[first:last,block:block_stop]),
                @view(rhs[block:block_stop,:]))
        end
    end
    return rhs
end

function _native_cell_schur!(target,left,right,cache=_native_cell_cache_sizes())
    isempty(target) && return target
    rows,columns=size(target);depth=size(left,2)
    kc,mc,nc=_native_cell_blocking(depth,rows,columns,1,cache)
    for row in 1:mc:rows,k in 1:kc:depth,column in 1:nc:columns
        last_row=min(rows,row+mc-1);last_k=min(depth,k+kc-1);last_column=min(columns,column+nc-1)
        _native_cell_product!(@view(target[row:last_row,column:last_column]),
            @view(left[row:last_row,k:last_k]),@view(right[k:last_k,column:last_column]))
        end
    return target
end

function _native_cell_lu!(lu,permutations,max_block=256,cache=_native_cell_cache_sizes())
    rows,columns=size(lu);n=min(rows,columns)
    if n<=16
        for k in 1:n
            pivot=k;biggest=abs(lu[k,k])
            @inbounds for row in k+1:rows
                score=abs(lu[row,k])
                if score>biggest
                    pivot=row;biggest=score
                end
            end
            biggest!=0.0 || throw(ArgumentError("singular native cell interpolation matrix"))
            permutations[k]=pivot
            @inbounds for column in 1:columns
                lu[k,column],lu[pivot,column]=lu[pivot,column],lu[k,column]
            end
            diagonal=lu[k,k]
            @inbounds for row in k+1:rows
                lu[row,k]/=diagonal
            end
            @inbounds for column in k+1:columns,row in k+1:rows
                lu[row,column]-=lu[row,k]*lu[k,column]
            end
        end
        return lu
    end
    block_size=min(max((n÷8)÷16*16,8),max_block)
    for start in 1:block_size:n
        stop=min(n,start+block_size-1)
        panel=@view lu[start:rows,start:stop]
        local_permutations=@view permutations[start:stop]
        _native_cell_lu!(panel,local_permutations,16,cache)
        @inbounds for i in start:stop
            pivot=(permutations[i]+=start-1)
            for column in 1:start-1
                lu[i,column],lu[pivot,column]=lu[pivot,column],lu[i,column]
            end
            if stop<rows
                for column in stop+1:columns
                    lu[i,column],lu[pivot,column]=lu[pivot,column],lu[i,column]
                end
            end
        end
        if stop<rows
            diagonal=@view lu[start:stop,start:stop]
            right=@view lu[start:stop,stop+1:columns]
            _native_cell_lower_solve!(diagonal,right,cache)
            _native_cell_schur!(@view(lu[stop+1:rows,stop+1:columns]),
                @view(lu[stop+1:rows,start:stop]),right,cache)
        end
    end
    return lu
end

function _native_cell_inverse(matrix)
    n=size(matrix,1)
    size(matrix,2)==n || throw(ArgumentError("native cell interpolation matrix must be square"))
    lu=copy(matrix);permutations=Vector{Int}(undef,n)
    cache=_native_cell_cache_sizes()
    _native_cell_lu!(lu,permutations,256,cache)
    result=zeros(Float64,n,n)
    @inbounds for i in 1:n
        result[i,i]=1.0
    end
    @inbounds for i in 1:n,column in 1:n
        pivot=permutations[i]
        result[i,column],result[pivot,column]=result[pivot,column],result[i,column]
    end
    _native_cell_lower_solve!(lu,result,cache)
    _native_cell_upper_solve!(lu,result,cache)
    return result
end

function _native_cell_polynomial(msh,caller)
    spec=msh_spec(msh)
    spec.family in (:tri,:qua,:tet,:hex) && 1<=spec.order<=10 ||
        throw(ArgumentError("$caller: type $msh has no native polynomial cell evaluator"))
    points=_native_cell_reference_nodes(spec)
    raw=spec.serendipity ? MeshFunctionSpaces._serendipity_exponents(spec.family,spec.order) :
        getfield(Elements,Symbol("_",spec.family,"_monomials"))(spec.order,false)
    exponents=NTuple{3,Int}[ntuple(axis->axis<=length(point) ? point[axis] : 0,3) for point in raw]
    count=length(exponents)
    matrix=Matrix{Float64}(undef,count,count)
    @inbounds for row in 1:count,column in 1:count
        value=1.0
        for axis in 1:spec.dim
            value*=_dim01_native_line_power(points[axis,column],exponents[row][axis])
        end
        matrix[row,column]=value
    end
    coefficients=_native_cell_inverse(matrix)
    all(isfinite,coefficients) || throw(ArgumentError("$caller: native cell coefficients are not finite"))
    return _NativeCellPolynomial(exponents,coefficients,spec.dim)
end

function _native_cell_shape(data::_NativeCellPolynomial,uvw)
    u,v,w=uvw;count=length(data.exponents)
    terms=Vector{Float64}(undef,count)
    powers=ntuple(axis->Float64[_dim01_native_line_power(uvw[axis],degree) for degree in 0:10],3)
    @inbounds for column in 1:count
        exponent=data.exponents[column];value=1.0
        for axis in data.dimension:-1:1
            value*=powers[axis][exponent[axis]+1]
        end
        terms[column]=value
    end
    values=zeros(Float64,count);gradients=zeros(Float64,3,count)
    @inbounds for row in 1:count,column in 1:count
        coefficient=data.coefficients[row,column]
        exponent=data.exponents[column]
        values[row]+=coefficient*terms[column]
        for axis in 1:data.dimension
            degree=exponent[axis]
            degree==0 && continue
            value=coefficient
            for component in 1:data.dimension
                value*=powers[component][exponent[component]+1-(component==axis)]
                component==axis && (value*=degree)
            end
            gradients[axis,row]+=value
        end
    end
    return values,[(gradients[1,row],gradients[2,row],gradients[3,row]) for row in 1:count]
end

function _native_cell_numeric_contract(spec,caller)
    spec.family===:pri && spec.order>2 && throw(ArgumentError(
        "$caller: Gmsh 4.15.2 does not implement Prism nodal full closure above order two"))
    spec.family===:pyr && spec.serendipity && spec.order>2 && throw(ArgumentError(
        "$caller: Gmsh 4.15.2 does not implement incomplete Pyramid bases above order two"))
    return nothing
end

function _native_cell_basis_functions(element_type,local_coord,function_space_type,wanted_orientations)
    caller="API.mesh.get_basis_functions"
    msh=MeshFunctionSpaces._checked_element_type(element_type,caller)
    space=MeshFunctionSpaces._function_space(function_space_type,caller)
    _,_,count,basis_type=MeshFunctionSpaces._basis_element_contract(msh,space,caller)
    spec=msh_spec(basis_type)
    (space.hierarchical || spec.order<=2 || spec.family===:lin) && return nothing
    _native_cell_numeric_contract(spec,caller)
    spec.family in (:tri,:qua,:tet,:hex) || return nothing
    coordinates,point_count=MeshFunctionSpaces._checked_local_coordinates(local_coord,caller)
    orientations=MeshFunctionSpaces._checked_orientation_sequence(wanted_orientations,1,false,caller)
    result=Vector{Float64}(undef,MeshFunctionSpaces._checked_result_length(
        length(orientations),point_count,count,space.components;caller=caller))
    isempty(result) && return Int32(space.components),result,Int32(1)
    data=_native_cell_polynomial(basis_type,caller)
    cursor=0
    for point in 1:point_count
        base=3point-2
        values,gradients=_native_cell_shape(data,(coordinates[base],coordinates[base+1],coordinates[base+2]))
        if space.kind===:grad_lagrange
            for gradient in gradients,component in 1:3
                cursor+=1;result[cursor]=gradient[component]
            end
        else
            for value in values
                cursor+=1;result[cursor]=value
            end
        end
    end
    return Int32(space.components),result,Int32(1)
end

function _native_cell_public_frame(mesh,cells,cell,dimension,shape,caller,tag)
    values,gradients=shape
    mapped=(0.0,0.0,0.0);du=mapped;dv=mapped;dw=mapped
    @inbounds for node in eachindex(values)
        point=_mixed_node(mesh,cells,cell,node);value=values[node];gradient=gradients[node]
        mapped=(mapped[1]+value*point[1],mapped[2]+value*point[2],mapped[3]+value*point[3])
        du=(du[1]+gradient[1]*point[1],du[2]+gradient[1]*point[2],du[3]+gradient[1]*point[3])
        dv=(dv[1]+gradient[2]*point[1],dv[2]+gradient[2]*point[2],dv[3]+gradient[2]*point[3])
        dw=(dw[1]+gradient[3]*point[1],dw[2]+gradient[3]*point[2],dw[3]+gradient[3]*point[3])
    end
    determinant=if dimension==2
        normal=(du[2]*dv[3]-du[3]*dv[2],du[3]*dv[1]-du[1]*dv[3],du[1]*dv[2]-du[2]*dv[1])
        magnitude=sqrt(normal[3]*normal[3]+normal[2]*normal[2]+normal[1]*normal[1])
        reciprocal=magnitude==0.0 ? 1.0 : 1.0/sqrt(normal[1]*normal[1]+normal[2]*normal[2]+normal[3]*normal[3])
        dw=(normal[1]*reciprocal,normal[2]*reciprocal,normal[3]*reciprocal)
        magnitude
    else
        du[1]*dv[2]*dw[3]+du[3]*dv[1]*dw[2]+du[2]*dv[3]*dw[1]-
            du[3]*dv[2]*dw[1]-du[1]*dv[3]*dw[2]-du[2]*dv[1]*dw[3]
    end
    all(isfinite,(mapped...,du...,dv...,dw...,determinant)) || throw(ArgumentError(
        "$caller: reference mapping for element $tag is not Float64-representable"))
    return (du...,dv...,dw...),determinant,mapped
end
