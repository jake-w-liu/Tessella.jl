# Exact integer Kannan--Bachem normal forms. Transformation order follows
# KBIPACK, including its first nonzero pivot and last indivisible Smith row.
mutable struct _NativeHomologyNormalForm
    canonical::Matrix{BigInt}
    left::Matrix{BigInt}
    right::Matrix{BigInt}
    left_inverted::Bool
    right_inverted::Bool
end

_native_homology_identity(n)=BigInt[i==j ? 1 : 0 for i in 1:n,j in 1:n]
function _native_homology_normal_form(matrix,left_inverted,right_inverted)
    rows,cols=size(matrix)
    rows==0 || cols==0 ? nothing : _NativeHomologyNormalForm(BigInt.(matrix),
        _native_homology_identity(rows),_native_homology_identity(cols),
        left_inverted,right_inverted)
end
function _native_homology_swap!(matrix,i,j,columns::Bool)
    i==j && return nothing
    if columns
        matrix[:,[i,j]]=matrix[:,[j,i]]
    else
        matrix[[i,j],:]=matrix[[j,i],:]
    end
    return nothing
end
function _native_homology_rotate!(matrix,i,j,a,b,c,d,columns::Bool)
    if columns
        for row in axes(matrix,1)
            x,y=matrix[row,i],matrix[row,j]
            matrix[row,i]=a*x+b*y;matrix[row,j]=c*x+d*y
        end
    else
        for col in axes(matrix,2)
            x,y=matrix[i,col],matrix[j,col]
            matrix[i,col]=a*x+b*y;matrix[j,col]=c*x+d*y
        end
    end
    return nothing
end
function _native_homology_add!(matrix,source,target,coefficient,columns::Bool)
    if columns
        for row in axes(matrix,1)
            matrix[row,target]+=coefficient*matrix[row,source]
        end
    else
        for col in axes(matrix,2)
            matrix[target,col]+=coefficient*matrix[source,col]
        end
    end
    return nothing
end
function _native_homology_negate!(matrix,index,columns::Bool)
    if columns
        for row in axes(matrix,1)
            matrix[row,index]=-matrix[row,index]
        end
    else
        for col in axes(matrix,2)
            matrix[index,col]=-matrix[index,col]
        end
    end
    return nothing
end

function _native_homology_hermite_eliminate!(nf,column)
    matrix=nf.canonical
    for index in 1:min(column-1,size(matrix,1))
        element=matrix[index,column]
        element==0 && continue
        pivot=matrix[index,index]
        divisor,a,b=Base.GMP.MPZ.gcdext(pivot,element)
        c=-div(element,divisor);d=div(pivot,divisor)
        _native_homology_rotate!(matrix,index,column,a,b,c,d,true)
        if nf.right_inverted
            _native_homology_rotate!(nf.right,index,column,a,b,c,d,true)
        else
            _native_homology_rotate!(nf.right,index,column,d,-c,-b,a,false)
        end
    end
    return nothing
end
function _native_homology_hermite_reduce!(nf,column)
    matrix=nf.canonical
    column<=size(matrix,1) || return nothing
    for j in column-1:-1:1, i in j+1:column
        pivot=matrix[i,i];element=matrix[i,j]
        if abs(element)>=abs(pivot) || element<0
            coefficient=cld(-element,pivot)
            _native_homology_add!(matrix,i,j,coefficient,true)
            if nf.right_inverted
                _native_homology_add!(nf.right,i,j,coefficient,true)
            else
                _native_homology_add!(nf.right,j,i,-coefficient,false)
            end
        end
    end
    return nothing
end

function _native_homology_make_hermite!(nf)
    matrix=nf.canonical;rows,columns=size(matrix);schur=1
    while schur<=rows && schur<=columns
        schur>1 && _native_homology_hermite_eliminate!(nf,schur)
        pivot=0
        for row in schur:rows
            if matrix[row,schur]!=0
                pivot=row;break
            end
        end
        if pivot==0
            _native_homology_swap!(matrix,schur,columns,true)
            _native_homology_swap!(nf.right,schur,columns,nf.right_inverted)
            columns-=1
            schur>1 && _native_homology_hermite_reduce!(nf,schur-1)
        else
            _native_homology_swap!(matrix,schur,pivot,false)
            _native_homology_swap!(nf.left,schur,pivot,!nf.left_inverted)
            if matrix[schur,schur]<0
                _native_homology_negate!(matrix,schur,true)
                _native_homology_negate!(nf.right,schur,nf.right_inverted)
            end
            _native_homology_hermite_reduce!(nf,schur)
            schur+=1
        end
    end
    for column in schur:columns
        _native_homology_hermite_eliminate!(nf,column)
        _native_homology_hermite_reduce!(nf,schur-1)
    end
    return nf
end

function _native_homology_hermite(matrix,left_inverted=false,right_inverted=true)
    nf=_native_homology_normal_form(matrix,left_inverted,right_inverted)
    nf===nothing && return nothing
    return _native_homology_make_hermite!(nf)
end

function _native_homology_transpose_normal_form!(nf)
    nf.canonical=permutedims(nf.canonical)
    oldleft=nf.left
    nf.left=permutedims(nf.right);nf.right=permutedims(oldleft)
    nf.left_inverted,nf.right_inverted=nf.right_inverted,nf.left_inverted
    return nothing
end

function _native_homology_smith(matrix,left_inverted=false,right_inverted=false)
    nf=_native_homology_normal_form(matrix,left_inverted,right_inverted)
    nf===nothing && return nothing
    rows,columns=size(nf.canonical)
    last_ready=0;first_ready=min(rows,columns)+1
    while first_ready>last_ready+1
        _native_homology_make_hermite!(nf)
        canonical=nf.canonical
        index=1
        while index<first_ready && canonical[index,index]!=0
            index+=1
        end
        first_ready=index
        index=last_ready+1
        while index<first_ready
            if canonical[index,index]==1
                if index!=last_ready+1
                    _native_homology_swap!(canonical,last_ready+1,index,false)
                    _native_homology_swap!(nf.left,last_ready+1,index,!nf.left_inverted)
                    _native_homology_swap!(canonical,last_ready+1,index,true)
                    _native_homology_swap!(nf.right,last_ready+1,index,nf.right_inverted)
                end
                last_ready+=1
            end
            index+=1
        end
        _native_homology_transpose_normal_form!(nf)
        _native_homology_make_hermite!(nf)
        _native_homology_transpose_normal_form!(nf)
        canonical=nf.canonical
        indivisible=0;trailing=0
        while indivisible==0 && last_ready+1<first_ready && trailing==0
            for column in last_ready+2:columns
                if canonical[last_ready+1,column]!=0
                    trailing=column;break
                end
            end
            if trailing==0
                pivot=canonical[last_ready+1,last_ready+1]
                for j in last_ready+2:first_ready-1, i in last_ready+2:j
                    rem(canonical[i,j],pivot)!=0 && (indivisible=i)
                end
                if indivisible!=0
                    _native_homology_add!(canonical,indivisible,last_ready+1,1,false)
                    if nf.left_inverted
                        _native_homology_add!(nf.left,indivisible,last_ready+1,1,false)
                    else
                        _native_homology_add!(nf.left,last_ready+1,indivisible,-1,true)
                    end
                else
                    last_ready+=1
                end
            end
        end
    end
    return nf
end
