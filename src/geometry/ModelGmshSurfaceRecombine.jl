using ..GmshLibm: _gm87_atan2, _gm87_exp

# The native producer follows meshGFaceOptimize.cpp. The fixed-node public
# Recombine kernel retains its separate established contract.
mutable struct _NativeSurfaceRecombine{F}
    world::Matrix{Float64}
    uv::Matrix{Float64}
    owners::Vector{Tuple{Int,Int32}}
    node_tags::Vector{UInt64}
    triangles::Vector{NTuple{3,Int32}}
    triangle_tags::Vector{UInt64}
    quads::Vector{NTuple{4,Int32}}
    quad_tags::Vector{UInt64}
    alive::BitVector
    high_water::UInt64
    lift::F
    caller::String
end

function _native_recombine_state(product,world,frame,caller)
    mesh=product.mesh
    return _NativeSurfaceRecombine(copy(world),copy(mesh.coords[1:2,:]),
        copy(product.node_owners),copy(product.node_tags),
        NTuple{3,Int32}[Tuple(mesh.tris[:,i]) for i in axes(mesh.tris,2)],
        copy(product.element_tags),NTuple{4,Int32}[],UInt64[],trues(nnodes(mesh)),
        product.element_high_water,_GmshNativeSurfaceLift(frame),String(caller))
end

_native_recombine_node_order(s)=sortperm(s.node_tags)

@inline _native_recombine_point(s,node)=
    (s.world[1,node],s.world[2,node],s.world[3,node])
@inline function _native_recombine_unit(a)
    n=Mesh2D._gmsh_norm3(a)
    n==0 && return a
    inverse=1.0/n
    return (a[1]*inverse,a[2]*inverse,a[3]*inverse)
end
@inline function _native_recombine_angle(a,b,c)
    first=Mesh2D._gmsh_sub3(a,b);second=Mesh2D._gmsh_sub3(c,b)
    return _gm87_atan2(Mesh2D._gmsh_norm3(Mesh2D._gmsh_cross3(first,second)),
        Mesh2D._gmsh_dot3(first,second))
end
function _native_recombine_eta(s,q::NTuple{4,Int32})
    p=ntuple(i->_native_recombine_point(s,q[i]),4)
    edges=ntuple(i->Mesh2D._gmsh_sub3(p[mod1(i+1,4)],p[i]),4)
    corners=ntuple(i->Mesh2D._gmsh_cross3(edges[i],edges[mod1(i+1,4)]),4)
    sign=1.0
    for i in 2:4
        Mesh2D._gmsh_dot3(corners[1],corners[i])<0 && (sign=-1.0)
    end
    angles=ntuple(i->min(180.0,180.0*_native_recombine_angle(
        p[i],p[mod1(i+1,4)],p[mod1(i+2,4)])/pi),4)
    deviation=abs(90.0-angles[1])
    for i in 2:4
        deviation=max(abs(90.0-angles[i]),deviation)
    end
    return sign*(1.0-deviation*(1.0/90.0))
end
function _native_recombine_gamma(s,t::NTuple{3,Int32})
    a,b,c=ntuple(i->_native_recombine_point(s,t[i]),3)
    va=_native_recombine_unit(Mesh2D._gmsh_sub3(c,b))
    vb=_native_recombine_unit(Mesh2D._gmsh_sub3(a,c))
    vc=_native_recombine_unit(Mesh2D._gmsh_sub3(b,a))
    sa=Mesh2D._gmsh_norm3(Mesh2D._gmsh_cross3(vb,vc))
    sb=Mesh2D._gmsh_norm3(Mesh2D._gmsh_cross3(vc,va))
    sc=Mesh2D._gmsh_norm3(Mesh2D._gmsh_cross3(va,vb))
    sa==0 && sb==0 && sc==0 && return 0.0
    return 2.0*(2.0*sa*sb*sc/(sa+sb+sc))
end
function _native_recombine_triangle_eta(s,t::NTuple{3,Int32})
    p=ntuple(i->_native_recombine_point(s,t[i]),3)
    angles=ntuple(i->180.0*_native_recombine_angle(
        p[i],p[mod1(i+1,3)],p[mod1(i+2,3)])/pi,3)
    smallest=min(min(angles[1],angles[2]),angles[3])
    return 1.0-abs(60.0-smallest)*(1.0/60.0)
end
function _native_recombine_birth!(s)
    s.high_water<typemax(UInt64) || throw(ArgumentError(
        "$(s.caller): native recombination constructor exhausts UInt64 labels"))
    s.high_water+=UInt64(1)
    return s.high_water
end
struct _NativeRecombinePair
    first::Int
    second::Int
    shared::NTuple{2,Int32}
    opposite::NTuple{2,Int32}
    eta::Float64
    quality::Float64
    cost::Int
end

# The shared independent engines retain each selected directed input edge.
# Graph vertices are triangle indices; this operation reserves no public labels.
function _native_recombine_matching_result(count,pairs,periodic,caller)
    isodd(count) && return nothing
    edges=Tuple{Int,Int,Int}[(p.first,p.second,p.cost) for p in pairs]
    for (a,b) in periodic
        # The native triangle map inserts vertex zero for a null second pointer.
        # Dense Julia triangle indices begin at one. Self-loops remain in the
        # input and are safely ignored for matching coverage by the engines.
        push!(edges,(a,b==0 ? 1 : b,100000))
    end
    result=OrderedMatching.ordered_perfect_matching(count,edges)
    result===nothing && return nothing
    return (edges=edges,selected=result.selected,physical_count=length(pairs))
end

function _native_recombine_matching(count,pairs,periodic,caller)
    result=_native_recombine_matching_result(count,pairs,periodic,caller)
    result===nothing && return nothing
    return Tuple{Int,Int}[(result.edges[index][1],result.edges[index][2])
                          for index in result.selected]
end

# Keep the established private cost-test entrypoint, backed by the single
# shared BSD engine instead of defining a second weighted implementation.
function _native_weighted_perfect_matching(count,input)
    mate=WeightedMatching.minimum_weight_matching(count,input)
    any(iszero,mate) && return nothing
    return Tuple{Int,Int}[(vertex,mate[vertex]) for vertex in 1:count
                          if mate[vertex]>vertex]
end

@inline _native_recombine_pair_less(a,b)=a.quality<b.quality

function _native_recombine_heap_adjust!(pairs,first,hole,count,value)
    top=hole;child=hole
    while child<(count-1)÷2
        child=2*(child+1)
        _native_recombine_pair_less(pairs[first+child],pairs[first+child-1]) && (child-=1)
        pairs[first+hole]=pairs[first+child];hole=child
    end
    if iseven(count) && child==(count-2)÷2
        child=2*(child+1)
        pairs[first+hole]=pairs[first+child-1];hole=child-1
    end
    parent=(hole-1)÷2
    while hole>top && _native_recombine_pair_less(pairs[first+parent],value)
        pairs[first+hole]=pairs[first+parent];hole=parent;parent=(hole-1)÷2
    end
    pairs[first+hole]=value
    return nothing
end

function _native_recombine_heap_sort!(pairs,first,last)
    count=last-first
    count<2 && return nothing
    for parent in (count-2)÷2:-1:0
        _native_recombine_heap_adjust!(pairs,first,parent,count,pairs[first+parent])
    end
    while count>1
        count-=1
        value=pairs[first+count];pairs[first+count]=pairs[first]
        _native_recombine_heap_adjust!(pairs,first,0,count,value)
    end
    return nothing
end

function _native_recombine_pair_partition!(pairs,first,last)
    a=first+1;b=first+(last-first)÷2;c=last-1
    median=if _native_recombine_pair_less(pairs[a],pairs[b])
        _native_recombine_pair_less(pairs[b],pairs[c]) ? b :
            _native_recombine_pair_less(pairs[a],pairs[c]) ? c : a
    elseif _native_recombine_pair_less(pairs[a],pairs[c])
        a
    elseif _native_recombine_pair_less(pairs[b],pairs[c])
        c
    else
        b
    end
    pairs[first],pairs[median]=pairs[median],pairs[first]
    left=first+1;right=last
    while true
        while _native_recombine_pair_less(pairs[left],pairs[first])
            left+=1
        end
        right-=1
        while _native_recombine_pair_less(pairs[first],pairs[right])
            right-=1
        end
        left>=right && return left
        pairs[left],pairs[right]=pairs[right],pairs[left]
        left+=1
    end
end

function _native_recombine_pair_introsort!(pairs,first,last,depth)
    while last-first>16
        if depth==0
            _native_recombine_heap_sort!(pairs,first,last)
            return nothing
        end
        depth-=1
        cut=_native_recombine_pair_partition!(pairs,first,last)
        _native_recombine_pair_introsort!(pairs,cut,last,depth)
        last=cut
    end
    return nothing
end

function _native_recombine_pair_insert!(pairs,index)
    value=pairs[index];previous=index-1
    while _native_recombine_pair_less(value,pairs[previous])
        pairs[index]=pairs[previous];index=previous;previous-=1
    end
    pairs[index]=value
    return nothing
end

function _native_recombine_pair_sort!(pairs)
    count=length(pairs)
    count==0 && return nothing
    # GNU's quality-only std::sort rearranges ties at the median partition,
    # then inserts the first sixteen entries before unguarded insertion.
    # The native DLL whole-mesh controls independently check this policy.
    depth=2*(8*sizeof(Int)-1-leading_zeros(count))
    _native_recombine_pair_introsort!(pairs,1,count+1,depth)
    guarded=min(count,16)
    for index in 2:guarded
        if _native_recombine_pair_less(pairs[index],pairs[1])
            value=pairs[index]
            for previous in index-1:-1:1
                pairs[previous+1]=pairs[previous]
            end
            pairs[1]=value
        else
            _native_recombine_pair_insert!(pairs,index)
        end
    end
    for index in guarded+1:count
        _native_recombine_pair_insert!(pairs,index)
    end
    return nothing
end
@inline function _native_recombine_opposite(triangle,a,b)
    for node in triangle
        node!=a && node!=b && return node
    end
    throw(ArgumentError("native recombination triangle has no opposite vertex"))
end
@inline function _native_recombine_matched_quad(first,second)
    other=Int32(0);start=0
    for node in second
        if !(node in first)
            other=node;break
        end
    end
    for i in 1:3
        if !(first[i] in second)
            start=i;break
        end
    end
    other!=0 && start!=0 || throw(ArgumentError(
        "native recombination triangles have no distinct opposite vertices"))
    return (first[start],first[mod1(start+1,3)],other,first[mod1(start+2,3)])
end
function _native_recombine_pairs!(s,tags,protected)
    adjacency=Dict{NTuple{2,Int32},Vector{Int}}()
    oriented=Dict{NTuple{2,Int32},NTuple{2,Int32}}()
    for (cell,nodes) in enumerate(s.triangles),i in 1:3
        endpoints=(nodes[i],nodes[mod1(i+1,3)])
        edge=minmax(endpoints...)
        # MEdge keeps the first inserted endpoint order; its sorted indices
        # determine only the map comparison (MEdge.h, buildEdgeToElement).
        haskey(oriented,edge) || (oriented[edge]=endpoints)
        push!(get!(adjacency,edge,Int[]),cell)
    end
    edges=sort!(collect(keys(adjacency));by=e->minmax(tags[e[1]],tags[e[2]]))
    pairs=_NativeRecombinePair[]
    periodic=Dict{Int32,Tuple{Int,Int}}()
    for edge in edges
        cells=adjacency[edge]
        protected!==nothing && edge in protected && continue
        if length(cells)==2
            _native_recombine_birth!(s) # RecombineTriangle's temporary quad
            first,second=cells
            a,b=oriented[edge]
            c=_native_recombine_opposite(s.triangles[first],a,b)
            d=_native_recombine_opposite(s.triangles[second],a,b)
            eta=_native_recombine_eta(s,(a,c,b,d))
            p1,p2,p3,p4=ntuple(i->_native_recombine_point(s,(a,d,b,c)[i]),4)
            quality=abs(90.0-180.0*_native_recombine_angle(p1,p2,p3)/pi)
            for (p,q,r) in ((p2,p3,p4),(p3,p4,p1),(p4,p1,p2))
                quality=max(abs(90.0-180.0*_native_recombine_angle(p,q,r)/pi),quality)
            end
            cost=trunc(Int,1000.0*_gm87_exp(-eta))
            boundary=count(node->s.owners[node][1]<2,(a,b,c,d))
            if cost>trunc(Int,1000.0*_gm87_exp(0.1)) && boundary>2
                cost=5000
            elseif cost>=1000 && boundary>2
                cost=10000
            end
            push!(pairs,_NativeRecombinePair(first,second,(a,b),(c,d),eta,quality,cost))
        elseif length(cells)==1
            for node in edge
                old=get(periodic,node,nothing)
                if old===nothing
                    periodic[node]=(cells[1],0)
                elseif old[1]!=cells[1]
                    periodic[node]=(old[1],cells[1])
                else
                    delete!(periodic,node)
                end
            end
        else
            throw(ArgumentError("$(s.caller): native recombination has nonmanifold incidence"))
        end
    end
    _native_recombine_pair_sort!(pairs)
    return pairs,[periodic[node] for node in sort!(collect(keys(periodic));by=n->tags[n])]
end
function _native_recombine_initial!(s,tags,protected,blossom,angle)
    pairs,periodic=_native_recombine_pairs!(s,tags,protected)
    touched=falses(length(s.triangles))
    matches=blossom ? _native_recombine_matching_result(length(touched),pairs,periodic,s.caller) : nothing
    if blossom && iseven(length(touched)) && matches===nothing
        # The native failed perfect-match path re-enters the greedy operation.
        # Rebuilding the pairs reserves its second set of temporary Quad labels.
        return _native_recombine_initial!(s,tags,protected,false,angle)
    end
    if matches!==nothing
        for index in matches.selected
            # A selected closure edge emits no physical Quad, even when its
            # endpoints parallel a real pair. Preserve input edge identity.
            index<=matches.physical_count || continue
            first,second,_=matches.edges[index]
            q=_native_recombine_matched_quad(s.triangles[first],s.triangles[second])
            push!(s.quads,q);push!(s.quad_tags,_native_recombine_birth!(s))
            touched[first]=true;touched[second]=true
        end
    else
        for pair in pairs
            pair.eta<angle || continue
            touched[pair.first] || touched[pair.second] || begin
                a,b=pair.shared;c,d=pair.opposite
                triangle=s.triangles[pair.first]
                slot=findfirst(==(a),triangle)
                forward=triangle[mod1(slot+1,3)]==b
                q=forward ? (a,d,b,c) : (a,c,b,d)
                push!(s.quads,q);push!(s.quad_tags,_native_recombine_birth!(s))
                touched[pair.first]=true;touched[pair.second]=true
            end
        end
    end
    s.triangles=s.triangles[.!touched];s.triangle_tags=s.triangle_tags[.!touched]
    return nothing
end
function _native_recombine_adjacency(s;triangles=true)
    adjacency=[Tuple{Bool,Int}[] for _ in axes(s.world,2)]
    if triangles
        for (cell,nodes) in enumerate(s.triangles),node in nodes
            push!(adjacency[node],(false,cell))
        end
    end
    for (cell,nodes) in enumerate(s.quads),node in nodes
        push!(adjacency[node],(true,cell))
    end
    return adjacency
end
@inline _native_recombine_cell(s,ref)=ref[1] ? s.quads[ref[2]] : s.triangles[ref[2]]

@inline function _native_recombine_uv_triangle_area(s,a,b,c)
    return 0.5*abs((s.uv[1,b]-s.uv[1,a])*(s.uv[2,c]-s.uv[2,a])-
        (s.uv[1,c]-s.uv[1,a])*(s.uv[2,b]-s.uv[2,a]))
end
function _native_recombine_uv_area(s,ref)
    q=_native_recombine_cell(s,ref)
    ref[1] || return _native_recombine_uv_triangle_area(s,q...)
    first=_native_recombine_uv_triangle_area(s,q[1],q[2],q[3])+
        _native_recombine_uv_triangle_area(s,q[3],q[4],q[1])
    second=_native_recombine_uv_triangle_area(s,q[2],q[3],q[4])+
        _native_recombine_uv_triangle_area(s,q[4],q[1],q[2])
    return min(second,first)
end
@inline _native_recombine_cell_eta(s,ref)=ref[1] ?
    _native_recombine_eta(s,s.quads[ref[2]]) :
    _native_recombine_triangle_eta(s,s.triangles[ref[2]])

function _native_recombine_laplace_admit!(s,node,actual,trial,incident)
    area_old=0.0;minimum_old=1.0
    for ref in incident
        area_old+=_native_recombine_uv_area(s,ref)
        minimum_old=min(_native_recombine_cell_eta(s,ref),minimum_old)
    end
    previous=_native_recombine_point(s,node)
    s.uv[:,node].=trial;s.world[:,node].=s.lift(trial...)
    area_new=0.0;minimum_new=1.0
    for ref in incident
        area_new+=_native_recombine_uv_area(s,ref)
        minimum_new=min(_native_recombine_cell_eta(s,ref),minimum_new)
    end
    s.uv[:,node].=actual;s.world[:,node].=previous
    return !((1.0+1e-10)*area_old<area_new || minimum_new<minimum_old)
end

function _native_recombine_laplace!(s,iterations,frame)
    iterations==0 && return nothing
    adjacency=_native_recombine_adjacency(s)
    order=_native_recombine_node_order(s)
    metric=(Mesh2D._gmsh_dot3(frame.first,frame.first),
        Mesh2D._gmsh_dot3(frame.second,frame.first),
        Mesh2D._gmsh_dot3(frame.second,frame.second))
    for _ in 1:iterations,node in order
        s.alive[node] && s.owners[node][1]==2 || continue
        incident=adjacency[node];isempty(incident) && continue
        neighbors=Int32[]
        for ref in incident
            cell=_native_recombine_cell(s,ref)
            slot=findfirst(==(node),cell)
            push!(neighbors,cell[mod1(slot-1,length(cell))],
                cell[mod1(slot+1,length(cell))])
        end
        unique!(neighbors);sort!(neighbors;by=n->s.node_tags[n])
        before=(s.uv[1,node],s.uv[2,node])
        u=0.0;v=0.0;count=0.0
        for neighbor in neighbors
            adj=(s.uv[1,neighbor],s.uv[2,neighbor])
            direction=_native_recombine_unit((adj[1]-before[1],adj[2]-before[2],0.0))
            weight=sqrt(metric[1]*direction[1]*direction[1]+
                2.0*metric[2]*direction[1]*direction[2]+
                metric[3]*direction[2]*direction[2])
            u+=adj[1]*weight;v+=adj[2]*weight;count+=weight
        end
        inverse=1.0/count;after=(u*inverse,v*inverse)
        factor=1.0;actual=before
        for _ in 1:5
            trial=(after[1]*factor+before[1]*(1.0-factor),
                after[2]*factor+before[2]*(1.0-factor))
            if _native_recombine_laplace_admit!(s,node,actual,trial,incident)
                actual=trial
                s.uv[:,node].=trial;s.world[:,node].=s.lift(trial...)
            end
            factor/=1.4
        end
    end
    return nothing
end

function _gmsh_native_surface_smooth(m,tag,product,world,frame,caller)
    state=_native_recombine_state(product,world,frame,caller)
    _native_recombine_laplace!(state,m.meshing.native_smoothing,frame)
    parameters=copy(product.node_parameters)
    owner=(2,Int32(tag))
    for row in eachindex(parameters)
        state.owners[row]==owner && (parameters[row]=collect(state.uv[:,row]))
    end
    coords=copy(product.mesh.coords);coords[1:2,:].=state.uv
    mesh=Mesh(coords;tris=copy(product.mesh.tris))
    result=_ModelSurfaceProduction(mesh,product.node_tags,product.node_owners,
        parameters,product.retained_nodes,product.node_high_water,
        product.element_tags,product.element_high_water)
    return result,state.world
end
function _native_recombine_objective!(s,node,target,original,incident,xi)
    old=_native_recombine_point(s,node)
    u=target[1]*(1.0-xi)+original[1]*xi
    v=target[2]*(1.0-xi)+original[2]*xi
    s.world[:,node].=s.lift(u,v)
    quality=1.0
    for ref in incident
        q=ref[1] ? _native_recombine_eta(s,s.quads[ref[2]]) :
            abs(_native_recombine_gamma(s,s.triangles[ref[2]]))
        quality=min(q,quality)
    end
    s.world[:,node].=old
    return quality
end
function _native_recombine_golden!(s,node,target,original,incident,tolerance)
    lambda=0.5*(sqrt(5.0)-1.0);mu=0.5*(3.0-sqrt(5.0))
    a=0.0;b=1.0
    worst=_native_recombine_objective!(s,node,target,original,incident,0.0)
    worst>0.5 && return 0.0
    x1=b-lambda*(b-a);x2=a+lambda*(b-a)
    f1=_native_recombine_objective!(s,node,target,original,incident,x1)
    f2=_native_recombine_objective!(s,node,target,original,incident,x2)
    tolerance<0 && return f1>f2 ? x1 : x2
    while !(abs(b-a)<tolerance)
        if f1<f2
            a=x1
            abs(b-a)<tolerance && break
            x1=x2;f1=f2;x2=b-mu*(b-a)
            f2=_native_recombine_objective!(s,node,target,original,incident,x2)
        else
            b=x2
            abs(b-a)<tolerance && break
            x2=x1;f2=f1;x1=a+mu*(b-a)
            f1=_native_recombine_objective!(s,node,target,original,incident,x1)
        end
    end
    final=_native_recombine_objective!(s,node,target,original,incident,a)
    return final<worst ? 0.0 : a
end
function _native_recombine_relocate!(s,iterations)
    adjacency=_native_recombine_adjacency(s)
    order=_native_recombine_node_order(s)
    for _ in 1:iterations,node in order
        s.alive[node] && s.owners[node][1]==2 || continue
        incident=adjacency[node];isempty(incident) && continue
        u=0.0;v=0.0;count=0
        for ref in incident,neighbor in _native_recombine_cell(s,ref)
            u+=s.uv[1,neighbor];v+=s.uv[2,neighbor];count+=1
        end
        inverse=1.0/count;target=(u*inverse,v*inverse)
        original=(s.uv[1,node],s.uv[2,node])
        xi=_native_recombine_golden!(s,node,target,original,incident,1e-6)
        next=(target[1]*(1.0-xi)+original[1]*xi,target[2]*(1.0-xi)+original[2]*xi)
        s.uv[:,node].=next;s.world[:,node].=s.lift(next...)
    end
    return nothing
end
function _native_recombine_two_quads!(s,owner)
    total=0
    while true
        adjacency=_native_recombine_adjacency(s)
        touched=falses(length(s.quads));removed=0
        for node in _native_recombine_node_order(s)
            refs=adjacency[node]
            length(refs)==2 && s.owners[node]==owner && all(ref->ref[1],refs) || continue
            first_cell,second_cell=refs[1][2],refs[2][2]
            touched[first_cell] || touched[second_cell] || begin
                q1=s.quads[first_cell];q2=s.quads[second_cell]
                slot=findfirst(==(Int32(node)),q1)
                a=q1[mod1(slot+1,4)];b=q1[mod1(slot+2,4)];c=q1[mod1(slot+3,4)]
                others=[n for n in q2 if n!=a && n!=c && n!=node]
                length(others)==1 || throw(ArgumentError(
                    "$(s.caller): invalid two-quad Surface node incidence"))
                push!(s.quads,(a,b,c,only(others)))
                push!(s.quad_tags,_native_recombine_birth!(s))
                touched[first_cell]=true;touched[second_cell]=true;s.alive[node]=false;removed+=1
            end
        end
        append!(touched,falses(length(s.quads)-length(touched)))
        s.quads=s.quads[.!touched];s.quad_tags=s.quad_tags[.!touched]
        total+=removed;removed==0 && return total
    end
end
function _native_recombine_collapse!(s,adjacency,quad,first,second)
    incidents=vcat(adjacency[first],adjacency[second])
    old1=_native_recombine_point(s,first);old2=_native_recombine_point(s,second)
    uv1=(s.uv[1,first],s.uv[2,first]);uv2=(s.uv[1,second],s.uv[2,second])
    midpoint=(0.5*(uv1[1]+uv2[1]),0.5*(uv1[2]+uv2[2]))
    point=s.lift(midpoint...)
    old_quality=1.0;new_quality=1.0
    for (_,cell) in incidents
        cell==quad && continue
        old_quality=min(old_quality,_native_recombine_eta(s,s.quads[cell]))
        s.world[:,first].=point;s.world[:,second].=point
        s.uv[:,first].=midpoint;s.uv[:,second].=midpoint
        new_quality=min(new_quality,_native_recombine_eta(s,s.quads[cell]))
        s.world[:,first].=old1;s.world[:,second].=old2
        s.uv[:,first].=uv1;s.uv[1,second]=uv2[1]
        # Literal pinned restoration writes v1.v (meshGFaceOptimize.cpp728).
        s.uv[2,first]=uv2[2]
    end
    if new_quality>old_quality
        s.world[:,first].=point;s.uv[:,first].=midpoint
        for (_,cell) in incidents
            cell==quad && continue
            s.quads[cell]=ntuple(i->s.quads[cell][i]==second ? first : s.quads[cell][i],4)
        end
        return true
    end
    return false
end
function _native_recombine_diamonds!(s,owner)
    total=0
    while true
        adjacency=_native_recombine_adjacency(s;triangles=false)
        touched=falses(length(s.alive));removed=falses(length(s.quads));count=0
        for t in s.triangles,node in t
            touched[node]=true
        end
        for cell in eachindex(s.quads)
            q=s.quads[cell]
            any(node->touched[node],q) && continue
            all(node->s.owners[node]==owner,q) || continue
            for (a,b) in ((q[1],q[3]),(q[2],q[4]))
                length(adjacency[a])==3 && length(adjacency[b])==3 || continue
                if _native_recombine_collapse!(s,adjacency,cell,a,b)
                    touched[collect(q)].=true;removed[cell]=true;s.alive[b]=false
                    count+=1;break
                end
            end
        end
        s.quads=s.quads[.!removed];s.quad_tags=s.quad_tags[.!removed]
        total+=count;count==0 && return total
    end
end
function _native_recombine_split_quality!(s,minimum)
    keep=trues(length(s.quads))
    for (cell,q) in enumerate(s.quads)
        _native_recombine_eta(s,q)<minimum+1e-12 || continue
        first=(q[1],q[2],q[3]);second=(q[3],q[4],q[1])
        third=(q[2],q[3],q[4]);fourth=(q[4],q[1],q[2])
        tags=ntuple(_->_native_recombine_birth!(s),4)
        quality1=min(_native_recombine_gamma(s,first),_native_recombine_gamma(s,second))
        quality2=min(_native_recombine_gamma(s,third),_native_recombine_gamma(s,fourth))
        p=ntuple(i->_native_recombine_point(s,q[i]),4)
        n1=Mesh2D._gmsh_cross3(Mesh2D._gmsh_sub3(p[3],p[2]),Mesh2D._gmsh_sub3(p[4],p[2]))
        n2=Mesh2D._gmsh_cross3(Mesh2D._gmsh_sub3(p[1],p[4]),Mesh2D._gmsh_sub3(p[2],p[4]))
        if quality1>quality2 || Mesh2D._gmsh_dot3(n1,n2)<0
            push!(s.triangles,first,second);push!(s.triangle_tags,tags[1],tags[2])
        else
            push!(s.triangles,third,fourth);push!(s.triangle_tags,tags[3],tags[4])
        end
        keep[cell]=false
    end
    s.quads=s.quads[keep];s.quad_tags=s.quad_tags[keep]
    return nothing
end
function _native_recombine_process!(state,m,tag,protected;
                                    minimum_quality=m.meshing.recombine_minimum_quality)
    isempty(state.triangle_tags) && (state.triangle_tags=
        UInt64[_native_recombine_birth!(state) for _ in state.triangles])
    settings=m.meshing
    _native_recombine_initial!(state,state.node_tags,protected,
        settings.recombine_algo in (1,3),get(settings.recombine,(2,tag),45.0))
    !iszero(settings.recombine_node_repositioning) &&
        _native_recombine_relocate!(state,settings.native_smoothing)
    owner=(2,Int32(tag))
    if settings.recombine_topology_passes>0 &&
       all(entity->entity[1]!=2 || entity==owner,state.owners)
        iteration=0
        while true
            two=_native_recombine_two_quads!(state,owner)
            diamonds=_native_recombine_diamonds!(state,owner)
            !iszero(settings.recombine_node_repositioning) &&
                _native_recombine_relocate!(state,settings.native_smoothing)
            iteration+=1
            (iteration>settings.recombine_topology_passes || two+diamonds==0) && break
        end
    end
    minimum_quality>0 && _native_recombine_split_quality!(state,minimum_quality)
    !iszero(settings.recombine_node_repositioning) &&
        _native_recombine_relocate!(state,settings.native_smoothing)
    return nothing
end

function _native_recombine_export(product,state,tag,caller;
                                   retained=product.retained_nodes,
                                   parameters=product.node_parameters,
                                   node_high_water=product.node_high_water)
    rows=findall(state.alive);mapping=zeros(Int32,length(state.alive))
    mapping[rows]=Int32.(1:length(rows))
    blocks=MixedElementBlock[]
    for (msh,cells,width) in ((2,state.triangles,3),(3,state.quads,4))
        isempty(cells) && continue
        nodes=Matrix{Int32}(undef,width,length(cells))
        for (column,cell) in enumerate(cells),i in 1:width
            mapping[cell[i]]>0 || throw(ArgumentError(
                "$caller: native recombination retained a deleted node"))
            nodes[i,column]=mapping[cell[i]]
        end
        push!(blocks,ElementBlock(msh,nodes))
    end
    coords=zeros(Float64,3,length(rows));coords[1:2,:].=state.uv[:,rows]
    result=MixedMesh(coords,blocks)
    parameters=copy(parameters[rows])
    owner=(2,Int32(tag))
    for (i,row) in enumerate(rows)
        state.owners[row]==owner && (parameters[i]=collect(state.uv[:,row]))
    end
    return _ModelSurfaceProduction(result,state.node_tags[rows],state.owners[rows],
        parameters,retained[rows],node_high_water,
        vcat(state.triangle_tags,state.quad_tags),state.high_water)
end

function _gmsh_native_surface_recombine(m,tag,product,world,protected,frame,caller)
    state=_native_recombine_state(product,world,frame,caller)
    _native_recombine_process!(state,m,tag,protected)
    return _native_recombine_export(product,state,tag,caller)
end
