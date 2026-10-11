# Oriented cellular complex used by the native integral homology backend.
# Cells keep their producer orientation; sorted public vertex tags identify
# duplicates and define the initial CellPtrLessThan order.
mutable struct _NativeHomologyCell
    dim::Int
    nodes::Vector{UInt64}
    number::Int
    subdomain::Bool
    immune::Bool
    combined::Bool
    boundary::Dict{Int,NTuple{2,Int}}
    coboundary::Dict{Int,NTuple{2,Int}}
    components::Dict{Int,Int}
end

mutable struct _NativeHomologyComplex
    cells::Vector{_NativeHomologyCell}
    active::Vector{Set{Int}}
    original::Vector{Set{Int}}
    dim::Int
    relative::Bool
    simplicial::Bool
    nextnumber::Int
    smallest::Int
    biggest::Int
    biggest_measure::Float64
    smallest_measure::Float64
    previous::Vector{Int}
    following::Vector{Int}
    heads::Vector{Int}
    tails::Vector{Int}
end

# Cell/ElemChain sort vertex incidences by the full size_t label, then their
# initial cell comparator reads each label through a signed 32-bit int.
# Keep actual vertex objects in UInt64; only this compatibility key narrows.
# Equal keys reuse the first producer cell, as CellComplex::insert does.
@inline _native_homology_vertex_number(tag::UInt64)=
    reinterpret(Int32,UInt32(tag & UInt64(0xffffffff)))
function _native_homology_key(nodes)
    ordered=sort(nodes)
    return (length(ordered),Tuple(_native_homology_vertex_number(tag) for tag in ordered))
end
_native_homology_order(complex,ids)=sort!(collect(ids);by=id->complex.cells[id].number)
function _native_homology_adjacency(complex,map;original::Bool=false)
    slot=original ? 2 : 1
    return _native_homology_order(complex,
        (id for (id,value) in map if value[slot]!=0))
end
_native_homology_boundary(complex,id;original::Bool=false)=
    _native_homology_adjacency(complex,complex.cells[id].boundary;original)
_native_homology_coboundary(complex,id;original::Bool=false)=
    _native_homology_adjacency(complex,complex.cells[id].coboundary;original)
function _native_homology_active(complex,dim)
    ids=Int[];id=complex.heads[dim+1]
    while id!=0
        push!(ids,id);id=complex.following[id]
    end
    return ids
end

function _native_homology_degree(complex,map)
    count=0;first=0;number=typemax(Int)
    for (id,value) in map
        value[1]==0 && continue
        count+=1
        if complex.cells[id].number<number
            first=id;number=complex.cells[id].number
        end
    end
    return count,first
end

function _native_homology_activate!(complex,id)
    dim=complex.cells[id].dim+1
    previous=complex.tails[dim]
    complex.previous[id]=previous;complex.following[id]=0
    previous==0 ? (complex.heads[dim]=id) : (complex.following[previous]=id)
    complex.tails[dim]=id
    push!(complex.active[dim],id)
    return nothing
end

function _native_homology_patterns(dim::Int,n::Int)
    dim==1 && n==2 && return ((1,),(2,))
    dim==2 && n==3 && return ((1,2),(2,3),(3,1))
    dim==2 && n==4 && return ((1,2),(2,3),(3,4),(4,1))
    dim==3 && n==4 && return ((1,3,2),(1,2,4),(1,4,3),(4,2,3))
    dim==3 && n==8 && return ((1,4,3,2),(1,2,6,5),(1,5,8,4),
        (2,3,7,6),(3,4,8,7),(5,6,7,8))
    dim==3 && n==6 && return ((1,3,2),(4,5,6),(1,2,5,4),
        (1,4,6,3),(2,3,6,5))
    dim==3 && n==5 && return ((1,2,5),(4,1,5),(2,3,5),(3,4,5),
        (1,4,3,2))
    dim==0 && n==1 && return ()
    throw(ArgumentError("compute_homology!: unsupported cellular shape ($dim,$n)"))
end

function _native_homology_orientation(nodes,reference)
    n=length(nodes)
    n==length(reference) || return 0
    n==1 && return 1
    if n==2
        nodes==reference && return 1
        nodes[1]==reference[2] && nodes[2]==reference[1] && return -1
        return 0
    end
    for shift in 0:n-1
        all(nodes[i]==reference[mod1(i+shift,n)] for i in 1:n) && return 1
        all(nodes[i]==reference[mod1(shift-i+2,n)] for i in 1:n) && return -1
    end
    return 0
end

function _native_homology_new_cell!(complex,dim,nodes,subdomain;
                                    combined::Bool=false,components=nothing)
    id=length(complex.cells)+1
    parts=components===nothing ? Dict{Int,Int}(id=>1) : components
    number=combined ? (complex.nextnumber+=1) : 0
    push!(complex.cells,_NativeHomologyCell(dim,UInt64[nodes...],number,
        subdomain,false,combined,Dict{Int,NTuple{2,Int}}(),
        Dict{Int,NTuple{2,Int}}(),parts))
    push!(complex.previous,0);push!(complex.following,0)
    return id
end

function _native_homology_add_incidence!(complex,source,target,ori,
                                        boundary::Bool,other::Bool=true)
    ori==0 && return nothing
    cell=complex.cells[source]
    map=boundary ? cell.boundary : cell.coboundary
    current,original=get(map,target,(0,0))
    current=Base.Checked.checked_add(current,ori)
    if current==0
        if original==0
            delete!(map,target)
        else
            map[target]=(0,original)
        end
        _native_homology_remove_incidence!(complex,target,source,!boundary,false)
    else
        map[target]=(current,original)
        other && _native_homology_add_incidence!(complex,target,source,ori,
                                               !boundary,false)
    end
    return nothing
end

function _native_homology_remove_incidence!(complex,source,target,
                                           boundary::Bool,other::Bool=true)
    map=boundary ? complex.cells[source].boundary : complex.cells[source].coboundary
    value=get(map,target,nothing)
    value===nothing && return nothing
    if value[2]==0
        delete!(map,target)
    else
        map[target]=(0,value[2])
    end
    other && _native_homology_remove_incidence!(complex,target,source,!boundary,false)
    return nothing
end

function _native_homology_remove_cell!(complex,id;other::Bool=true)
    cell=complex.cells[id]
    for child in _native_homology_boundary(complex,id)
        _native_homology_remove_incidence!(complex,child,id,false,other)
    end
    for parent in _native_homology_coboundary(complex,id)
        _native_homology_remove_incidence!(complex,parent,id,true,other)
    end
    dim=cell.dim+1
    if id in complex.active[dim]
        previous,following=complex.previous[id],complex.following[id]
        previous==0 ? (complex.heads[dim]=following) : (complex.following[previous]=following)
        following==0 ? (complex.tails[dim]=previous) : (complex.previous[following]=previous)
        complex.previous[id]=0;complex.following[id]=0
        delete!(complex.active[dim],id)
    end
    return nothing
end

function _native_homology_restore!(complex)
    for dim in 0:3
        empty!(complex.active[dim+1]);complex.heads[dim+1]=0;complex.tails[dim+1]=0
        for id in _native_homology_order(complex,complex.original[dim+1])
            _native_homology_activate!(complex,id)
            cell=complex.cells[id]
            for map in (cell.boundary,cell.coboundary)
                for key in collect(keys(map))
                    value=map[key][2]
                    value==0 ? delete!(map,key) : (map[key]=(value,value))
                end
            end
            cell.immune=false
        end
    end
    return complex
end

function _native_homology_primary(msh_type::Int,nodes)
    family=msh_family(msh_type)
    dim,n= family===:pnt ? (0,1) : family===:lin ? (1,2) :
        family===:tri ? (2,3) : family===:qua ? (2,4) :
        family===:tet ? (3,4) : family===:hex ? (3,8) :
        family===:pri ? (3,6) : family===:pyr ? (3,5) : (-1,0)
    dim<0 && throw(ArgumentError("compute_homology!: unsupported element type $msh_type"))
    length(nodes)>=n || throw(ArgumentError("compute_homology!: truncated element connectivity"))
    return dim,UInt64[nodes[i] for i in 1:n]
end

function _native_homology_insert_elements!(complex,elements,subdomain,coordinates;
                                           measure_geometry::Bool=true)
    index=[Dict{Any,Int}() for _ in 0:3]
    for dim in 0:3, id in complex.active[dim+1]
        index[dim+1][_native_homology_key(complex.cells[id].nodes)]=id
    end
    smallest=fill(0,4);biggest=fill(0,4)
    smallsize=fill(-1.0,4);bigsize=fill(-1.0,4)
    top=0
    for (msh_type,nodes) in elements
        dim,primary=_native_homology_primary(Int(msh_type),nodes)
        length(unique(primary))==length(primary) || continue
        (dim==2 && length(primary)==4 || dim==3 && length(primary)!=4) &&
            (complex.simplicial=false)
        top=max(top,dim)
        key=_native_homology_key(primary)
        id=get(index[dim+1],key,0)
        if id==0
            id=_native_homology_new_cell!(complex,dim,primary,subdomain)
            index[dim+1][key]=id;push!(complex.active[dim+1],id)
        elseif subdomain
            complex.cells[id].subdomain=true
        end
        if !subdomain && measure_geometry
            size=_native_homology_measure(Int(msh_type),nodes,coordinates)
            if smallsize[dim+1]<0 || smallsize[dim+1]>size
                smallest[dim+1]=id;smallsize[dim+1]=size
            end
            if bigsize[dim+1]<0 || bigsize[dim+1]<size
                biggest[dim+1]=id;bigsize[dim+1]=size
            end
        end
    end
    smallest_id=smallest[top+1];biggest_id=biggest[top+1]
    for dim in 3:-1:1
        ids=sort!(collect(complex.active[dim+1]);
                  by=id->_native_homology_key(complex.cells[id].nodes))
        for id in ids
            parent=complex.cells[id]
            for (slot,pattern) in enumerate(_native_homology_patterns(dim,length(parent.nodes)))
                nodes=UInt64[parent.nodes[i] for i in pattern]
                length(unique(nodes))==length(nodes) || continue
                key=_native_homology_key(nodes)
                child=get(index[dim],key,0)
                if child==0
                    child=_native_homology_new_cell!(complex,dim-1,nodes,subdomain)
                    index[dim][key]=child;push!(complex.active[dim],child)
                elseif subdomain
                    complex.cells[child].subdomain=true
                end
                if !subdomain
                    # A comparator-equivalent child can retain another full
                    # UID. Native line incidence follows that retained UID,
                    # independently of the boundary slot being visited.
                    child_nodes=complex.cells[child].nodes
                    ori=dim==1 ? (child_nodes[1]==parent.nodes[1] ? -1 :
                        child_nodes[1]==parent.nodes[2] ? 1 : 0) :
                        _native_homology_orientation(child_nodes,nodes)
                    _native_homology_add_incidence!(complex,id,child,ori,true)
                    smallest_id==id && (smallest_id=child)
                    biggest_id==id && (biggest_id=child)
                end
            end
        end
    end
    complex.smallest=smallest_id;complex.biggest=biggest_id
    complex.smallest_measure=smallsize[top+1];complex.biggest_measure=bigsize[top+1]
    return nothing
end

function _native_homology_complex(domain_elements,subdomain_elements,coordinates;
                                  measure_geometry::Bool=true)
    complex=_NativeHomologyComplex(_NativeHomologyCell[],[Set{Int}() for _ in 0:3],
        [Set{Int}() for _ in 0:3],0,false,true,0,0,0,-1.0,-1.0,
        Int[],Int[],zeros(Int,4),zeros(Int,4))
    _native_homology_insert_elements!(complex,subdomain_elements,true,coordinates;measure_geometry)
    complex.relative=!isempty(complex.active[1])
    _native_homology_insert_elements!(complex,domain_elements,false,coordinates;measure_geometry)
    number=0
    for dim in 0:3
        isempty(complex.active[dim+1]) || (complex.dim=dim)
        ids=sort!(collect(complex.active[dim+1]);
                  by=id->_native_homology_key(complex.cells[id].nodes))
        for id in ids
            cell=complex.cells[id];cell.number=(number+=1)
            for map in (cell.boundary,cell.coboundary), (key,value) in map
                map[key]=(value[1],value[1])
            end
        end
        complex.original[dim+1]=copy(complex.active[dim+1])
        empty!(complex.active[dim+1])
        for id in ids
            _native_homology_activate!(complex,id)
        end
    end
    complex.nextnumber=number
    return complex
end

_native_homology_cell_type(cell::_NativeHomologyCell)=
    cell.dim==0 ? 15 : cell.dim==1 ? 1 : cell.dim==2 ?
    (length(cell.nodes)==3 ? 2 : 3) : length(cell.nodes)==4 ? 4 :
    length(cell.nodes)==5 ? 7 : length(cell.nodes)==6 ? 6 : 5
