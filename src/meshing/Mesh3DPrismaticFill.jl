# Bounded rigid-prism recognition uses coordinates only to bind literal source
# rows. Cap/side topology and the final full PLC certificate decide acceptance.
function _prismatic_row_pairing(surface,partner,tolerance)
    coordinates=surface.coords;n=size(coordinates,2)
    vector=ntuple(axis->coordinates[axis,partner]-coordinates[axis,1],3)
    length3=hypot(vector...)
    isfinite(length3) && length3>2tolerance || return nothing
    direction=ntuple(axis->vector[axis]/length3,3)
    bottom=Int32[];top=Int32[]
    for row in 1:n
        height=sum((coordinates[axis,row]-coordinates[axis,1])*direction[axis] for axis in 1:3)
        isfinite(height) || return nothing
        low=abs(height)<=tolerance;high=abs(height-length3)<=tolerance
        low==high && return nothing
        push!(low ? bottom : top,Int32(row))
    end
    length(bottom)==length(top) && length(bottom)>=3 || return nothing
    pairs=Dict{Int32,Int32}();inverse=Dict{Int32,Int32}()
    for row in bottom
        selected=Int32[candidate for candidate in top if hypot(ntuple(axis->
            coordinates[axis,candidate]-(coordinates[axis,row]+vector[axis]),3)...)<=tolerance]
        length(selected)==1 || return nothing
        other=only(selected);haskey(inverse,other) && return nothing
        pairs[row]=other;inverse[other]=row
    end
    length(inverse)==length(top) || return nothing
    return (;bottom,top,pairs,inverse)
end

function _prismatic_source_topology(surface,pairing)
    bottom=Set(pairing.bottom);top=Set(pairing.top)
    caps=NTuple{3,Int32}[];cap_keys=Set{NTuple{3,Int32}}();top_keys=Set{NTuple{3,Int32}}()
    sides=Dict{Tuple{Int32,Int32},Vector{NTuple{3,Int32}}}()
    for column in axes(surface.tris,2)
        triangle=Tuple(surface.tris[:,column]);allunique(triangle) || return nothing
        if all(row->row in bottom,triangle)
            key=_sort3t(triangle...);key in cap_keys && return nothing
            push!(caps,triangle);push!(cap_keys,key)
        elseif all(row->row in top,triangle)
            key=_sort3t((pairing.inverse[row] for row in triangle)...)
            key in top_keys && return nothing
            push!(top_keys,key)
        else
            base=unique(Int32[get(pairing.inverse,row,row) for row in triangle])
            length(base)==2 || return nothing
            key=minmax(base...);push!(get!(()->NTuple{3,Int32}[],sides,key),triangle)
        end
    end
    !isempty(caps) && cap_keys==top_keys || return nothing
    Set(row for triangle in caps for row in triangle)==bottom || return nothing
    edges=Dict{Tuple{Int32,Int32},Int}()
    for triangle in caps,(a,b) in ((1,2),(2,3),(3,1))
        key=minmax(triangle[a],triangle[b]);edges[key]=get(edges,key,0)+1
    end
    all(count->count in (1,2),values(edges)) || return nothing
    border=Set(key for (key,count) in edges if count==1)
    !isempty(border) && Set(keys(sides))==border || return nothing
    successors=Dict(row=>Int32[] for row in pairing.bottom)
    indegree=Dict(row=>0 for row in pairing.bottom)
    for key in sort!(collect(border))
        triangles=sides[key];length(triangles)==2 || return nothing
        a,b=key;ta,tb=pairing.pairs[a],pairing.pairs[b]
        first=Set(triangles[1]);second=Set(triangles[2])
        union(first,second)==Set((a,b,ta,tb)) || return nothing
        diagonal=intersect(first,second);length(diagonal)==2 || return nothing
        precedence=diagonal==Set((ta,b)) ? (a,b) : diagonal==Set((a,tb)) ? (b,a) : nothing
        precedence===nothing && return nothing
        before,after=precedence
        if !(after in successors[before])
            push!(successors[before],after);indegree[after]+=1
        end
    end
    ready=sort!(Int32[row for row in pairing.bottom if indegree[row]==0]);order=Int32[]
    while !isempty(ready)
        row=popfirst!(ready);push!(order,row)
        for other in successors[row]
            indegree[other]-=1
            if indegree[other]==0;push!(ready,other);sort!(ready);end
        end
    end
    length(order)==length(pairing.bottom) || return nothing
    return (;caps,order)
end

function _prismatic_staircase_fill(surface,pairing,topology)
    ranks=Dict(row=>index for (index,row) in enumerate(topology.order))
    cells=NTuple{4,Int32}[]
    for triangle in topology.caps
        a,b,c=sort!(collect(triangle);by=row->ranks[row])
        ta,tb,tc=pairing.pairs[a],pairing.pairs[b],pairing.pairs[c]
        for cell in ((a,b,c,ta),(b,c,ta,tb),(c,ta,tb,tc))
            points=ntuple(slot->ntuple(axis->surface.coords[axis,cell[slot]],3),4)
            orientation=-orient3(points...)
            isfinite(orientation) && orientation!=0 || return nothing
            positive=orientation>0 ? cell : (cell[1],cell[2],cell[4],cell[3])
            push!(cells,positive)
        end
    end
    connectivity=Matrix{Int32}(undef,4,length(cells))
    for (column,cell) in enumerate(cells);connectivity[:,column]=collect(cell);end
    candidate=Mesh(surface.coords;tets=connectivity)
    validate(candidate).ok || return nothing
    # Recognition never alters source coordinates or relaxes PLC conformity.
    accepted,_=_certify_surface_fill(surface,candidate)
    return accepted ? candidate : nothing
end

function _try_prismatic_surface_fill(surface::Mesh)
    n=size(surface.coords,2)
    6<=n<=128 && iseven(n) && size(surface.tris,2)>=8 || return nothing
    extent=ntuple(axis->maximum(surface.coords[axis,:])-minimum(surface.coords[axis,:]),3)
    span=hypot(extent...);isfinite(span) && span>0 || return nothing
    magnitude=maximum(abs,surface.coords)
    tolerance=max(64eps(Float64)*magnitude,1e-10span)
    isfinite(tolerance) && tolerance>0 || return nothing
    for partner in 2:n
        pairing=_prismatic_row_pairing(surface,partner,tolerance);pairing===nothing && continue
        topology=_prismatic_source_topology(surface,pairing);topology===nothing && continue
        candidate=_prismatic_staircase_fill(surface,pairing,topology)
        candidate===nothing || return candidate
    end
    return nothing
end
