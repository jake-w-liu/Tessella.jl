# Boundary matrices, integral quotients and representative deformation.
_native_homology_matrix(matrices,dim)=0<=dim<length(matrices) ? matrices[dim+1] : nothing

function _native_homology_ker_cod(matrix)
    matrix===nothing && return nothing,nothing
    nf=_native_homology_hermite(matrix,false,true)
    nf===nothing && return nothing,nothing
    rank=0
    while rank<min(size(nf.canonical)...) && nf.canonical[rank+1,rank+1]!=0
        rank+=1
    end
    kernel=rank==size(matrix,2) ? nothing : copy(nf.right[:,rank+1:end])
    image=rank==0 ? nothing : nf.left*nf.canonical[:,1:rank]
    return kernel,image
end

function _native_homology_inclusion(kernel,image)
    (kernel===nothing || image===nothing) && return nothing
    size(image,1)<size(image,2) && return nothing
    size(kernel,1)<size(kernel,2) && return nothing
    nf=_native_homology_smith(kernel,true,true)
    columns=size(kernel,2)
    for i in 1:columns
        nf.canonical[i,i]==0 && return nothing
    end
    coordinates=(nf.left*image)[1:columns,:]
    for i in axes(coordinates,1), j in axes(coordinates,2)
        divisor=nf.canonical[i,i]
        rem(coordinates[i,j],divisor)==0 || return nothing
        coordinates[i,j]=div(coordinates[i,j],divisor)
    end
    return nf.right*coordinates
end

function _native_homology_quotient_integer(inclusion)
    inclusion===nothing && return nothing,BigInt[]
    nf=_native_homology_smith(inclusion,false,false)
    rows,columns=size(inclusion);torsion=BigInt[]
    for i in 1:columns
        nf.canonical[i,i]==0 && return nothing,torsion
        nf.canonical[i,i]>1 && push!(torsion,nf.canonical[i,i])
    end
    rank=columns-length(torsion)
    return rows-rank>0 ? copy(nf.left[:,rank+1:rows]) : nothing,torsion
end

function _native_homology_boundary_matrices(complex)
    indices=[Int[id for id in _native_homology_active(complex,dim)
                 if !complex.cells[id].subdomain] for dim in 0:3]
    matrices=Union{Nothing,Matrix{BigInt}}[nothing for _ in 0:4]
    for dim in 0:3
        upper=indices[dim+1];isempty(upper) && continue
        lower=dim==0 ? Int[] : indices[dim]
        matrix=zeros(BigInt,max(1,length(lower)),length(upper))
        lower_index=Dict{Int,Int}(id=>row for (row,id) in enumerate(lower))
        for (column,id) in enumerate(upper)
            for child in _native_homology_boundary(complex,id)
                row=get(lower_index,child,0)
                row==0 && continue
                matrix[row,column]+=complex.cells[id].boundary[child][1]
            end
        end
        matrices[dim+1]=matrix
    end
    return matrices,indices
end

function _native_homology_integer_basis(complex,dual)
    matrices,indices=_native_homology_boundary_matrices(complex)
    if dual
        for i in eachindex(matrices)
            matrices[i]===nothing || (matrices[i]=permutedims(matrices[i]))
        end
    end
    kernels=Union{Nothing,Matrix{BigInt}}[nothing for _ in 0:4]
    images=Union{Nothing,Matrix{BigInt}}[nothing for _ in 0:4]
    bases=Union{Nothing,Matrix{BigInt}}[nothing for _ in 0:4]
    torsions=[BigInt[] for _ in 0:4]
    for i in -1:3
        low=dual ? complex.dim+1-i : i
        high=dual ? complex.dim-i : i+1
        target=dual ? high : low
        high_matrix=_native_homology_matrix(matrices,high)
        if 0<=high<=3
            kernels[high+1],images[high+1]=_native_homology_ker_cod(high_matrix)
        end
        0<=target<=4 || continue
        low_matrix=_native_homology_matrix(matrices,low)
        if low==0 && !dual && low_matrix!==nothing && high_matrix===nothing
            bases[target+1]=_native_homology_identity(size(low_matrix,2))
        elseif high==0 && dual && high_matrix!==nothing && low_matrix===nothing
            bases[target+1]=_native_homology_identity(size(high_matrix,1))
        elseif _native_homology_matrix(matrices,target)===nothing
            bases[target+1]=nothing
        elseif high_matrix===nothing
            kernel=_native_homology_matrix(kernels,low)
            bases[target+1]=kernel===nothing ? nothing : copy(kernel)
        else
            if low_matrix===nothing && 0<=low<=4
                kernels[low+1]=_native_homology_identity(size(high_matrix,1))
            end
            kernel=_native_homology_matrix(kernels,low)
            image=_native_homology_matrix(images,high)
            inclusion=_native_homology_inclusion(kernel,image)
            quotient,torsion=_native_homology_quotient_integer(inclusion)
            torsions[target+1]=torsion
            if image===nothing
                bases[target+1]=kernel===nothing ? nothing : copy(kernel)
            elseif inclusion===nothing || quotient===nothing
                bases[target+1]=nothing
            else
                bases[target+1]=kernel*quotient
            end
        end
    end
    return bases,torsions,indices
end

function _native_homology_deform!(complex,chain,inside,outside)
    coefficients=Int[];incidences=Int[]
    for id in _native_homology_order(complex,keys(inside))
        cell=complex.cells[id];cell.immune=false
        if !cell.subdomain
            push!(coefficients,get(chain,id,0));push!(incidences,inside[id])
        end
    end
    (isempty(coefficients) || complex.dim==2 && length(coefficients)<2) && return false
    inout=coefficients[1]*incidences[1]
    all(coefficients[i]*incidences[i]==inout for i in eachindex(coefficients)) || return false
    for id in keys(inside)
        haskey(chain,id) && (chain[id]=0)
    end
    for (position,id) in enumerate(_native_homology_order(complex,keys(outside)))
        complex.cells[id].immune=position==2
        coefficient=-inout*outside[id]
        get(chain,id,0)==0 || throw(ArgumentError("compute_homology!: invalid chain deformation"))
        chain[id]=coefficient
    end
    return true
end

function _native_homology_deform_chain!(complex,chain,id,bend)
    cell=complex.cells[id];dim=cell.dim
    # firstCoboundary(true) skips only the leading zero-original entries;
    # later entries remain in the upstream map traversal.
    ordered=_native_homology_order(complex,keys(cell.coboundary))
    first_original=findfirst(parent->cell.coboundary[parent][2]!=0,ordered)
    first_original===nothing && return false
    for parent in ordered[first_original:end]
        inside=Dict{Int,Int}();outside=Dict{Int,Int}()
        for child in _native_homology_boundary(complex,parent;original=true)
            coefficient=complex.cells[parent].boundary[child][2]
            if get(chain,child,0)!=0 || complex.cells[child].subdomain
                inside[child]=coefficient
            else
                outside[child]=coefficient
            end
        end
        skip=false
        for child in _native_homology_order(complex,keys(inside))
            other=complex.cells[child]
            other.immune && (skip=true)
            other.subdomain && (bend=false)
            abs(get(chain,child,0))>1 && (skip=true)
        end
        any(child->complex.cells[child].subdomain,keys(outside)) && (skip=true)
        skip && continue
        ni,no=length(inside),length(outside)
        if (dim==1 && ni==2 && no==1) || (dim==2 && ni==3 && no==1) ||
           (bend && ((dim==1 && ni==1 && no==2) || (dim==2 && ni==2 && no==2))) ||
           (dim==1 && ni==3 && no==0) || (dim==2 && ni==4 && no==0)
            return _native_homology_deform!(complex,chain,inside,outside)
        end
    end
    return false
end

function _native_homology_visit_chain!(complex,chain,dim,bend)
    lastnumber=0
    while true
        following=0;number=typemax(Int)
        for id in keys(chain)
            candidate=complex.cells[id].number
            if lastnumber<candidate<number
                following=id;number=candidate
            end
        end
        following==0 && break
        bend && dim==2 && _native_homology_deform_chain!(complex,chain,following,true)
        _native_homology_deform_chain!(complex,chain,following,false)
        lastnumber=number
    end
    return nothing
end

function _native_homology_smoothen!(complex,chain,dim)
    (!complex.simplicial || length(chain)<2) && return nothing
    useless=0
    for iteration in 1:20
        size=length(chain)
        _native_homology_visit_chain!(complex,chain,dim,true)
        for id in keys(chain)
            complex.cells[id].immune=false
        end
        filter!(pair->pair.second!=0,chain)
        useless=size>=length(chain) ? useless+1 : 0
        useless>5 && break
    end
    for id in keys(chain)
        complex.cells[id].immune=false
    end
    _native_homology_visit_chain!(complex,chain,dim,false)
    filter!(pair->pair.second!=0,chain)
    return nothing
end

function _native_homology_basis_chains!(complex,dimensions;dual::Bool=false)
    _native_homology_restore!(complex)
    _native_homology_reduce!(complex;dual)
    if dual
        for dim in 2:-1:0
            dim in dimensions || _native_homology_combine!(complex,dim+1)
        end
    end
    bases,torsions,indices=_native_homology_integer_basis(complex,dual)
    output=[Dict{Int,Int}[] for _ in 0:3]
    for dim in (dual ? (3:-1:0) : (0:3))
        basis=bases[dim+1];basis===nothing && continue
        for column in axes(basis,2)
            chain=Dict{Int,Int}()
            torsion=column<=length(torsions[dim+1]) ? torsions[dim+1][column] : BigInt(1)
            for (row,id) in enumerate(indices[dim+1])
                coefficient=basis[row,column];coefficient==0 && continue
                for child in _native_homology_order(complex,keys(complex.cells[id].components))
                    value=BigInt(complex.cells[id].components[child])*coefficient*torsion
                    value==0 && continue
                    typemin(Int)<=value<=typemax(Int) || throw(ArgumentError(
                        "compute_homology!: representative coefficient exceeds native capacity"))
                    chain[child]=Int(value)
                end
            end
            !dual && dim in (1,2) && _native_homology_smoothen!(complex,chain,dim)
            isempty(chain) || push!(output[dim+1],chain)
        end
    end
    return output
end
