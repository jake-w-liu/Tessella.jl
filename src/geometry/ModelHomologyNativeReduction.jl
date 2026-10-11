# Elementary collapse, dual collapse, omission and cellular combinations.
# Numeric Cell order and the iterator advances are part of representative
# selection, rather than just a way to compute the Betti numbers.
function _native_homology_reduction!(complex,dim,omit,omitted;dual::Bool=false)
    1<=dim<=3 || return 0
    count=0;changed=true
    while changed
        changed=false
        id=complex.heads[dual ? dim+1 : dim]
        while id!=0
            cell=complex.cells[id];following=complex.following[id]
            map=dual ? cell.boundary : cell.coboundary
            degree,neighbour=_native_homology_degree(complex,map)
            if degree==1
                other=complex.cells[neighbour]
                if cell.subdomain==other.subdomain && !cell.immune &&
                   !other.immune && abs(map[neighbour][1])<2
                    if dual
                        dim-1==omit && push!(omitted,neighbour)
                        _native_homology_remove_cell!(complex,neighbour)
                        _native_homology_remove_cell!(complex,id)
                    else
                        dim==omit && push!(omitted,neighbour)
                        _native_homology_remove_cell!(complex,neighbour)
                        _native_homology_remove_cell!(complex,id)
                    end
                    count+=1;changed=true
                    # Upstream advances once before deleting the current
                    # cell, then once at the bottom of the traversal.
                    following==0 || (following=complex.following[following])
                end
            end
            isempty(complex.active[dim+1]) && break
            isempty(complex.active[dim]) && break
            id=following
        end
    end
    return count
end

function _native_homology_enqueue!(complex,ids,queue,queued)
    for id in ids
        id in queued && continue
        push!(queued,id);push!(queue,id)
    end
    return nothing
end

function _native_homology_queue_coreduction!(complex,start,omit,omitted)
    queue=Int[start];queued=Set{Int}(queue);head=1;count=0
    while head<=length(queue)
        id=queue[head];head+=1;delete!(queued,id)
        cell=complex.cells[id]
        degree,lower=_native_homology_degree(complex,cell.boundary)
        if degree==1
            other=complex.cells[lower]
            if cell.subdomain==other.subdomain && !cell.immune && !other.immune &&
               abs(cell.boundary[lower][1])<2
                _native_homology_remove_cell!(complex,id)
                _native_homology_enqueue!(complex,
                    _native_homology_coboundary(complex,lower),queue,queued)
                _native_homology_remove_cell!(complex,lower)
                other.dim==omit && push!(omitted,lower)
                count+=1
            end
        elseif degree==0
            _native_homology_enqueue!(complex,
                _native_homology_coboundary(complex,id),queue,queued)
        end
    end
    return count
end

function _native_homology_combined_pair!(complex,first,second,orientation_match,dual)
    a,b=complex.cells[first],complex.cells[second]
    if length(a.components)<length(b.components)
        first,second=second,first;a,b=b,a
    end
    parts=copy(a.components)
    for (id,value) in b.components
        haskey(parts,id) || (parts[id]=orientation_match ? value : -value)
    end
    id=_native_homology_new_cell!(complex,a.dim,UInt64[],a.subdomain;
                                  combined=true,components=parts)
    complex.cells[id].immune=a.immune || b.immune
    for boundary in (true,false)
        amap=boundary ? a.boundary : a.coboundary
        bmap=boundary ? b.boundary : b.coboundary
        for child in _native_homology_adjacency(complex,amap)
            ori=amap[child][1]
            _native_homology_remove_incidence!(complex,child,first,!boundary,false)
            _native_homology_add_incidence!(complex,id,child,ori,boundary,true)
        end
        for child in _native_homology_adjacency(complex,bmap)
            value=bmap[child]
            ori=orientation_match ? value[1] : -value[1]
            bmap[child]=(ori,value[2])
            _native_homology_remove_incidence!(complex,child,second,!boundary,false)
            add_all=boundary ? !dual : dual
            if add_all || get(amap,child,(0,0))[1]==0
                _native_homology_add_incidence!(complex,id,child,ori,boundary,true)
            end
        end
    end
    return id
end

function _native_homology_combine!(complex,dim;dual::Bool=false)
    (dual ? 0<=dim<=2 : 1<=dim<=3) || return 0
    queue=Int[];queued=Set{Int}();head=1;count=0
    id=complex.heads[dim+1]
    while id!=0
        neighbours=dual ? _native_homology_coboundary(complex,id) :
                         _native_homology_boundary(complex,id)
        _native_homology_enqueue!(complex,neighbours,queue,queued)
        while head<=length(queue)
            shared=queue[head];head+=1;s=complex.cells[shared]
            parents=dual ? _native_homology_boundary(complex,shared) :
                           _native_homology_coboundary(complex,shared)
            if length(parents)==2 && (dual || !s.immune)
                first,second=parents;a,b=complex.cells[first],complex.cells[second]
                map=dual ? s.boundary : s.coboundary
                or1,or2=map[first][1],map[second][1]
                if first!=second && abs(or1)==abs(or2) &&
                   s.subdomain==a.subdomain==b.subdomain && a.immune==b.immune
                    _native_homology_remove_cell!(complex,shared)
                    for parent in (first,second)
                        neighbours=dual ? _native_homology_coboundary(complex,parent) :
                                         _native_homology_boundary(complex,parent)
                        _native_homology_enqueue!(complex,neighbours,queue,queued)
                    end
                    joined=_native_homology_combined_pair!(complex,first,second,or1!=or2,dual)
                    _native_homology_remove_cell!(complex,first)
                    _native_homology_remove_cell!(complex,second)
                    _native_homology_activate!(complex,joined)
                    id=complex.heads[dim+1]
                    count+=1
                end
            end
            delete!(queued,shared)
        end
        id=complex.following[id]
    end
    return count
end

function _native_homology_omit!(complex,id,dual)
    _native_homology_remove_cell!(complex,id;other=false)
    omitted=Int[id];dim=complex.cells[id].dim
    if dual
        _native_homology_queue_coreduction!(complex,id,dim,omitted)
        for k in 1:complex.dim
            _native_homology_reduction!(complex,k,dim,omitted;dual=true)
        end
    else
        for k in 3:-1:1
            _native_homology_reduction!(complex,k,dim,omitted)
        end
    end
    first=complex.cells[id]
    combined=_native_homology_new_cell!(complex,dim,UInt64[],first.subdomain;
        combined=true,components=Dict{Int,Int}(part=>1 for part in omitted))
    complex.cells[combined].immune=any(part->complex.cells[part].immune,omitted)
    return combined
end

function _native_homology_remove_subdomain!(complex)
    ids=Int[]
    for dim in 0:3, id in _native_homology_active(complex,dim)
        complex.cells[id].subdomain && push!(ids,id)
    end
    for id in ids
        _native_homology_remove_cell!(complex,id)
    end
    return nothing
end

function _native_homology_coherent!(complex)
    for dim in 0:3, id in _native_homology_active(complex,dim)
        for boundary in (true,false)
            cell=complex.cells[id];map=boundary ? cell.boundary : cell.coboundary
            for other in _native_homology_adjacency(complex,map)
                ori=map[other][1];neighbour=complex.cells[other]
                if !(other in complex.active[neighbour.dim+1])
                    _native_homology_remove_incidence!(complex,id,other,boundary,false)
                end
                reciprocal=boundary ? neighbour.coboundary : neighbour.boundary
                if get(reciprocal,id,(0,0))[1]==0
                    _native_homology_add_incidence!(complex,other,id,ori,!boundary,false)
                end
            end
        end
    end
    return nothing
end

function _native_homology_reduce!(complex;dual::Bool=false,combine::Int=3,
                                 omit::Bool=true,heuristic::Int=1)
    isempty(complex.active[1]) && return nothing
    complex.relative && _native_homology_remove_subdomain!(complex)
    empty=Int[]
    if dual
        for dim in 0:3
            for id in _native_homology_active(complex,dim)
                _native_homology_queue_coreduction!(complex,id,-1,empty)!=0 && break
            end
        end
        for dim in 1:complex.dim
            _native_homology_reduction!(complex,dim,-1,empty;dual=true)
        end
    else
        for dim in 3:-1:1
            _native_homology_reduction!(complex,dim,-1,empty)
        end
    end
    if omit
        newcells=Int[];dim=dual ? 0 : complex.dim
        while !isempty(complex.active[dim+1])
            id=complex.heads[dim+1]
            preferred=heuristic==1 ? complex.biggest : heuristic==-1 ? complex.smallest : 0
            size=heuristic==1 ? complex.biggest_measure : complex.smallest_measure
            if dual && preferred!=0 && size!=0 && preferred in complex.active[1]
                id=preferred
            end
            push!(newcells,_native_homology_omit!(complex,id,dual))
        end
        for id in newcells
            _native_homology_activate!(complex,id)
        end
    end
    for dim in (dual ? (0,1,2) : (3,2,1))
        combine>0 && _native_homology_combine!(complex,dim;dual)
        if combine>2
            for k in (dual ? (1:3) : (3:-1:1))
                _native_homology_reduction!(complex,k,-1,empty;dual)
            end
        elseif combine>1
            _native_homology_reduction!(complex,dual ? dim+1 : dim-1,-1,empty;dual)
        end
    end
    dual && _native_homology_coherent!(complex)
    return nothing
end
