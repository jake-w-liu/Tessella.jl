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

# ── Recombine3D surface fold ────────────────────────────────────────────────
# Upstream `Generator.cpp` `Mesh3D` runs `quadsToTriangles(gf, 1e6)` on every
# boundary GFace of a `recombine3D`-flagged (or `Mesh.Recombine3DAll`) Region
# before `MeshDelaunayVolume`, so the stored face mesh ends with the split
# triangles and the volume fill sees an all-triangle PLC. Folding the same
# dim-2 cells at the product boundary reproduces that contract: the published
# Surface blocks keep the split triangles in quad order and the Region's
# native hybrid admission sees no quadrangles.
#
# The diagonal choice mirrors `meshGFaceOptimize.cpp`'s — the better of the
# two candidate pairs by triangle gamma, with the second diagonal chosen when
# folding the first would revert orientation on a non-convex quad. It is the
# same test `_native_recombine_split_quality!` applies for
# `RecombineMinimumQuality`.
@inline function _recombine3d_tri_gamma(p,q,r)
    va=_native_recombine_unit(Mesh2D._gmsh_sub3(r,q))
    vb=_native_recombine_unit(Mesh2D._gmsh_sub3(p,r))
    vc=_native_recombine_unit(Mesh2D._gmsh_sub3(q,p))
    sa=Mesh2D._gmsh_norm3(Mesh2D._gmsh_cross3(vb,vc))
    sb=Mesh2D._gmsh_norm3(Mesh2D._gmsh_cross3(vc,va))
    sc=Mesh2D._gmsh_norm3(Mesh2D._gmsh_cross3(va,vb))
    sa==0 && sb==0 && sc==0 && return 0.0
    return 2.0*(2.0*sa*sb*sc/(sa+sb+sc))
end

function _recombine3d_fold_pair(coords,nodes::NTuple{4,Int32})
    p=ntuple(i->ntuple(axis->coords[axis,Int(nodes[i])],3),4)
    first=(nodes[1],nodes[2],nodes[3]);second=(nodes[3],nodes[4],nodes[1])
    third=(nodes[2],nodes[3],nodes[4]);fourth=(nodes[4],nodes[1],nodes[2])
    quality1=min(_recombine3d_tri_gamma(p[1],p[2],p[3]),
                 _recombine3d_tri_gamma(p[3],p[4],p[1]))
    quality2=min(_recombine3d_tri_gamma(p[2],p[3],p[4]),
                 _recombine3d_tri_gamma(p[4],p[1],p[2]))
    n1=Mesh2D._gmsh_cross3(Mesh2D._gmsh_sub3(p[3],p[2]),
                           Mesh2D._gmsh_sub3(p[4],p[2]))
    n2=Mesh2D._gmsh_cross3(Mesh2D._gmsh_sub3(p[1],p[4]),
                           Mesh2D._gmsh_sub3(p[2],p[4]))
    return quality1>quality2 || Mesh2D._gmsh_dot3(n1,n2)<0 ?
        ((first,second),(1,2)) : ((third,fourth),(3,4))
end

# Fold every linear quadrangle block into the published triangle set. The
# merged msh-2 block takes the first dim-2 block's position — the same
# triangles-first order upstream's appended `gf->triangles` produce.
function _recombine3d_fold_mesh(mesh::MixedMesh)
    any(block->block isa ElementBlock && block.msh==3 && !isempty(block.nodes),
        mesh.blocks) || return mesh
    triangles=NTuple{3,Int32}[];triangle_tags=Int32[]
    offsets=Vector{NTuple{2,Int}}()
    for block in mesh.blocks
        block isa ElementBlock || continue
        block.msh==2 || continue
        for column in axes(block.nodes,2)
            push!(triangles,Tuple(block.nodes[:,column]))
            push!(triangle_tags,block.tags[column])
        end
    end
    for block in mesh.blocks
        block isa ElementBlock && block.msh==3 || continue
        for column in axes(block.nodes,2)
            pair,offset=_recombine3d_fold_pair(mesh.coords,
                Tuple(block.nodes[:,column]))
            push!(triangles,pair[1],pair[2]);push!(offsets,offset)
            append!(triangle_tags,(block.tags[column],block.tags[column]))
        end
    end
    nodes=Matrix{Int32}(undef,3,length(triangles))
    for column in eachindex(triangles),row in 1:3
        nodes[row,column]=triangles[column][row]
    end
    blocks=ElementBlock[];inserted=false
    for block in mesh.blocks
        if block isa ElementBlock && block.msh in (2,3)
            inserted || (push!(blocks,ElementBlock(2,nodes,triangle_tags));
                         inserted=true)
        else
            push!(blocks,block)
        end
    end
    inserted || push!(blocks,ElementBlock(2,nodes,triangle_tags))
    return MixedMesh(copy(mesh.coords),blocks;
        physical_names=mesh.physical_names,
        elementary_entities=mesh.elementary_entities,
        periodic_links=mesh.periodic_links,
        ancillary_sections=mesh.ancillary_sections,
        data_sections=mesh.data_sections,
        partition_data=mesh.partition_data),
        length(offsets),offsets
end
_recombine3d_fold_mesh(mesh::Mesh)=(mesh,0,NTuple{2,Int}[])

# The fold keeps every node — parameters, owners, retained flags and birth
# rows carry over unchanged. Element births follow the native split-quality
# convention: each folded quad reserves four consecutive labels and emits
# the chosen pair's offsets, so the published tag list matches the native
# recombination receipt when it would have split the same quads.
function _recombine3d_fold_production(m::GeoModel,surface::Int,production)
    _model_surface_recombine3d(m,surface) || return production
    production isa _ModelSurfaceProduction || return production
    folded,nquads,offsets=_recombine3d_fold_mesh(production.mesh)
    folded===production.mesh && return production
    element_tags=production.element_tags
    element_high_water=production.element_high_water
    if !isempty(element_tags)
        tri_cells=sum(size(block.nodes,2) for block in production.mesh.blocks
            if block isa ElementBlock && block.msh==2;init=0)
        length(element_tags)==tri_cells+nquads || throw(ArgumentError(
            "recombine3d fold: production element receipts do not cover its dim-2 cells"))
        element_high_water+UInt64(4)*nquads<=typemax(UInt64) || throw(ArgumentError(
            "recombine3d fold: quad split births exhaust UInt64 labels"))
        base=element_high_water
        folded_tags=UInt64[]
        sizehint!(folded_tags,tri_cells+2*nquads)
        append!(folded_tags,element_tags[1:tri_cells])
        for (index,offset) in enumerate(offsets)
            birth=base+UInt64(4)*(index-1)
            push!(folded_tags,birth+offset[1],birth+offset[2])
        end
        element_tags=folded_tags;element_high_water=base+UInt64(4)*nquads
    end
    constructor_mesh=production.constructor_mesh===nothing ? nothing :
        first(_recombine3d_fold_mesh(production.constructor_mesh))
    return _ModelSurfaceProduction(folded,production.node_tags,
        production.node_owners,production.node_parameters,
        production.retained_nodes,production.node_high_water,
        element_tags,element_high_water,production.birth_rows,constructor_mesh,
        production.source_node_updates,production.native_recombine_events,
        production.curve_chains)
end

function _recombine3d_fold_surface_part(m::GeoModel,surface::Int,part)
    _model_surface_recombine3d(m,surface) || return part
    part isa _ModelSurfaceGeometryProduction && return _ModelSurfaceGeometryProduction(
        _recombine3d_fold_production(m,surface,part.production),
        first(_recombine3d_fold_mesh(part.geometry)))
    part isa _ModelSurfaceProduction && return _recombine3d_fold_production(m,surface,part)
    part isa MixedMesh && return first(_recombine3d_fold_mesh(part))
    return part
end

# ── Recombine3D retained-product fold ───────────────────────────────────────
# Upstream's `Mesh3D` (`Generator.cpp`) runs `quadsToTriangles` on the STORED
# face mesh of every boundary `GFace` of a `recombine3D` Region at volume-fill
# time, so a face meshed before `Recombine Volume`/`Mesh.Recombine3DAll` took
# effect still feeds `MeshDelaunayVolume` an all-triangle PLC — and stays
# split for every later read of that face's mesh.  The retained authority is
# that stored state: fold its actual Surface cells, the shared geometry
# carrier and the Surface constructor receipts identically, minting each
# split child a fresh element label.  `_dim01_ingest!` raises the assembly
# counter over the minted labels like any other retained source, and the
# republished products keep the split for good.
#
# Only the all-triangle Delaunay arm folds; Transfinite and Layers-extruded
# Regions are excluded by `_model_volume_recombine3d` exactly like upstream's
# `delaunay` list.
function _recombine3d_fold_authority(m::GeoModel,authority;caller="recombine3D retained fold",
                                     above::UInt64=UInt64(0))
    authority===nothing && return authority
    mesh=authority.mesh
    mesh isa MixedMesh && mesh.entity_data!==nothing || return authority
    surfaces=Set{Int32}()
    for (volume,shells) in m.volumes
        _model_volume_recombine3d(m,volume) || continue
        for shell in shells
            members=get(m.surface_loops,shell,nothing)
            members===nothing && continue
            for signed in members
                push!(surfaces,Int32(abs(signed)))
            end
        end
    end
    isempty(surfaces) && return authority
    data=mesh.entity_data
    # Mint above every label the assembler can already own — retained cells,
    # the caller's element high-water and discrete-record cells — so ingest
    # never renumbers a folded triangle away from its constructor/carrier
    # sibling labels.
    next=max(above,
        maximum(tags->maximum(tags;init=UInt64(0)),
            data.external_element_tags;init=UInt64(0)),
        maximum(pair->maximum(pair[2].element_tags;init=UInt64(0)),
            _discrete_mesh_records_model(m);init=UInt64(0)))
    # old actual quad label -> (two emitted triangle labels, their node-UID
    # supports in emitted order), so the carrier bindings, the stored records
    # and the constructor receipts adopt identical labels and support.
    tag_map=Dict{UInt64,NTuple{2,UInt64}}()
    uid_map=Dict{UInt64,NTuple{2,NTuple{3,UInt64}}}()
    blocks=ElementBlock[];owners_out=Vector{Int32}[];tags_out=Vector{UInt64}[]
    elementary_out=mesh.elementary_entities===nothing ? nothing :
        Vector{Int32}[]
    folded_any=false
    for (bi,block) in enumerate(mesh.blocks)
        block isa ElementBlock || throw(ArgumentError(
            "$caller: retained authority carries a non-element block"))
        owners=data.block_entities[bi];tags=data.external_element_tags[bi]
        elementary_b=mesh.elementary_entities===nothing ? nothing :
            mesh.elementary_entities[bi]
        if block.msh!=3 || isempty(block.nodes) ||
                !any(owner->owner in surfaces,owners)
            push!(blocks,block);push!(owners_out,owners);push!(tags_out,tags)
            elementary_out===nothing || push!(elementary_out,elementary_b)
            continue
        end
        folded_any=true
        tris=Int32[];tri_phys=Int32[];tri_owners=Int32[];tri_tags=UInt64[]
        quads=Int32[];quad_phys=Int32[];quad_owners=Int32[];quad_tags=UInt64[]
        tri_elementary=Int32[];quad_elementary=Int32[]
        for column in axes(block.nodes,2)
            if owners[column] in surfaces
                pair,_=_recombine3d_fold_pair(mesh.coords,Tuple(block.nodes[:,column]))
                first_tag=next+UInt64(1);second_tag=next+UInt64(2);next=second_tag
                tag_map[tags[column]]=(first_tag,second_tag)
                uid_map[tags[column]]=(Tuple(data.external_node_tags[pair[1][k]]
                                             for k in 1:3),
                                       Tuple(data.external_node_tags[pair[2][k]]
                                             for k in 1:3))
                for (conn,label) in zip(pair,(first_tag,second_tag))
                    append!(tris,conn)
                    push!(tri_phys,block.tags[column])
                    push!(tri_owners,owners[column]);push!(tri_tags,label)
                    elementary_b===nothing ||
                        push!(tri_elementary,elementary_b[column])
                end
            else
                append!(quads,block.nodes[:,column])
                push!(quad_phys,block.tags[column])
                push!(quad_owners,owners[column]);push!(quad_tags,tags[column])
                elementary_b===nothing ||
                    push!(quad_elementary,elementary_b[column])
            end
        end
        isempty(tris) || (push!(blocks,ElementBlock(2,reshape(tris,3,:),tri_phys));
            push!(owners_out,tri_owners);push!(tags_out,tri_tags);
            elementary_out===nothing || push!(elementary_out,tri_elementary))
        isempty(quads) || (push!(blocks,ElementBlock(3,reshape(quads,4,:),quad_phys));
            push!(owners_out,quad_owners);push!(tags_out,quad_tags);
            elementary_out===nothing || push!(elementary_out,quad_elementary))
    end
    folded_any || return authority
    # `Mesh 2` products also live as per-entity stored records, which ingest
    # dedups against the authoritative cells by (type, owner, connectivity).
    # A folded triangle never dedups a stored quad, so the same split must be
    # applied to every flagged Surface record — the second copy of upstream's
    # stored face state.  Elements whose stored support disagrees with the
    # actual quad are left untouched for ingest to arbitrate.
    for (owner,record) in _discrete_mesh_records_model(m)
        owner[1]==2 && Int32(owner[2]) in surfaces || continue
        isempty(record.element_tags) && continue
        any(i->record.element_types[i]==Int32(3) &&
            haskey(tag_map,record.element_tags[i]),
            eachindex(record.element_tags)) || continue
        types=Int32[];tags=UInt64[];nodes=Vector{UInt64}[]
        for i in eachindex(record.element_tags)
            label=record.element_tags[i]
            if record.element_types[i]==Int32(3) && haskey(tag_map,label)
                conns=uid_map[label]
                if Set(Iterators.flatten(conns))==Set(record.element_nodes[i])
                    for (child,support) in zip(tag_map[label],conns)
                        push!(types,Int32(2));push!(tags,child)
                        push!(nodes,collect(support))
                    end
                    continue
                end
            end
            push!(types,record.element_types[i]);push!(tags,label)
            push!(nodes,record.element_nodes[i])
        end
        record.element_types=types;record.element_tags=tags
        record.element_nodes=nodes
    end
    folded_data=Elements._native_mixed_entity_data(data.entities;
        node_entities=data.node_entities,node_parametric=data.node_parametric,
        external_node_tags=data.external_node_tags,
        block_entities=owners_out,external_element_tags=tags_out)
    folded_mesh=MixedMesh(mesh.coords,blocks;
        physical_names=mesh.physical_names,entity_data=folded_data,
        elementary_entities=elementary_out,
        periodic_links=mesh.periodic_links,
        ancillary_sections=mesh.ancillary_sections,
        data_sections=mesh.data_sections,partition_data=mesh.partition_data)
    # Surface constructor receipts mirror their surface's actual cells; each
    # folded quad keeps the two labels the actual fold minted for it.
    folded_inputs=Dict{Tuple{Int,Int32},MixedMesh}()
    for (owner,constructor) in authority.surface_constructor_inputs
        if !(owner[1]==2 && owner[2] in surfaces)
            folded_inputs[owner]=constructor;continue
        end
        cdata=constructor.entity_data
        cdata===nothing && throw(ArgumentError(
            "$caller: retained Surface[$(Int(owner[2]))] constructor receipt lacks element labels"))
        cblocks=ElementBlock[];cowners=Vector{Int32}[];ctags=Vector{UInt64}[]
        celem=constructor.elementary_entities===nothing ? nothing :
            Vector{Int32}[]
        for (bi,block) in enumerate(constructor.blocks)
            owners_b=cdata.block_entities[bi];tags_b=cdata.external_element_tags[bi]
            elem_b=constructor.elementary_entities===nothing ? nothing :
                constructor.elementary_entities[bi]
            if block.msh!=3 || isempty(block.nodes)
                push!(cblocks,block);push!(cowners,owners_b);push!(ctags,tags_b)
                celem===nothing || push!(celem,elem_b)
                continue
            end
            cnodes=Int32[];cphys=Int32[];co=Int32[];ct=UInt64[];ce=Int32[]
            for column in axes(block.nodes,2)
                mapped=get(tag_map,tags_b[column],nothing)
                mapped===nothing && throw(ArgumentError(
                    "$caller: retained Surface[$(Int(owner[2]))] constructor quad lost its actual element label"))
                pair,_=_recombine3d_fold_pair(constructor.coords,
                    Tuple(block.nodes[:,column]))
                for (conn,label) in zip(pair,mapped)
                    append!(cnodes,conn)
                    push!(cphys,block.tags[column])
                    push!(co,owners_b[column]);push!(ct,label)
                    elem_b===nothing || push!(ce,elem_b[column])
                end
            end
            push!(cblocks,ElementBlock(2,reshape(cnodes,3,:),cphys))
            push!(cowners,co);push!(ctags,ct)
            celem===nothing || push!(celem,ce)
        end
        new_data=Elements._native_mixed_entity_data(cdata.entities;
            node_entities=cdata.node_entities,node_parametric=cdata.node_parametric,
            external_node_tags=cdata.external_node_tags,
            block_entities=cowners,external_element_tags=ctags)
        folded_inputs[owner]=MixedMesh(constructor.coords,cblocks;
            physical_names=constructor.physical_names,entity_data=new_data,
            elementary_entities=celem,
            periodic_links=constructor.periodic_links,
            ancillary_sections=constructor.ancillary_sections,
            data_sections=constructor.data_sections,
            partition_data=constructor.partition_data)
    end
    # The geometry carrier folds the same cells through its actual-label
    # bindings; a carrier that is the actual mesh is already folded, and a
    # plain Mesh carrier holds only simplices.
    carrier=authority.carrier
    folded_carrier=carrier;folded_bindings=nothing
    if carrier===mesh
        folded_carrier=folded_mesh
        folded_bindings=tags_out
    elseif carrier isa MixedMesh && !isempty(authority.carrier_cells)
        cdata=carrier.entity_data
        cnext=cdata===nothing ? UInt64(0) : maximum(
            tags->maximum(tags;init=UInt64(0)),cdata.external_element_tags;
            init=UInt64(0))
        cblocks=ElementBlock[];cowners=Vector{Int32}[];ctags=Vector{UInt64}[]
        cbind=Vector{Union{Nothing,UInt64}}[]
        celem=carrier.elementary_entities===nothing ? nothing :
            Vector{Int32}[]
        for (bi,block) in enumerate(carrier.blocks)
            bind_b=authority.carrier_cells[bi]
            owners_b=cdata===nothing ? fill(Int32(0),size(block.nodes,2)) :
                cdata.block_entities[bi]
            tags_b=cdata===nothing ? fill(UInt64(0),size(block.nodes,2)) :
                cdata.external_element_tags[bi]
            elem_b=carrier.elementary_entities===nothing ? nothing :
                carrier.elementary_entities[bi]
            foldable=block.msh==3 && !isempty(block.nodes) &&
                any(column->(bind_b[column]!==nothing &&
                             haskey(tag_map,bind_b[column])) ||
                            (bind_b[column]===nothing &&
                             owners_b[column] in surfaces),
                    axes(block.nodes,2))
            if !foldable
                push!(cblocks,block);push!(cowners,owners_b)
                push!(ctags,tags_b);push!(cbind,bind_b)
                celem===nothing || push!(celem,elem_b)
                continue
            end
            cnodes=Int32[];cphys=Int32[];co=Int32[];ct=UInt64[]
            cb=Union{Nothing,UInt64}[];ce=Int32[]
            for column in axes(block.nodes,2)
                bound=bind_b[column]
                if bound!==nothing && haskey(tag_map,bound)
                    pair,_=_recombine3d_fold_pair(carrier.coords,
                        Tuple(block.nodes[:,column]))
                    mapped=tag_map[bound]
                    for (conn,label) in zip(pair,mapped)
                        append!(cnodes,conn)
                        push!(cphys,block.tags[column]);push!(co,owners_b[column])
                        cnext+=UInt64(1);push!(ct,cnext);push!(cb,label)
                        elem_b===nothing || push!(ce,elem_b[column])
                    end
                elseif bound===nothing && owners_b[column] in surfaces
                    # Stored face quad with no actual backing — upstream
                    # still splits it; the children stay unbound.
                    pair,_=_recombine3d_fold_pair(carrier.coords,
                        Tuple(block.nodes[:,column]))
                    for conn in pair
                        append!(cnodes,conn)
                        push!(cphys,block.tags[column]);push!(co,owners_b[column])
                        cnext+=UInt64(1);push!(ct,cnext);push!(cb,nothing)
                        elem_b===nothing || push!(ce,elem_b[column])
                    end
                else
                    append!(cnodes,block.nodes[:,column])
                    push!(cphys,block.tags[column]);push!(co,owners_b[column])
                    push!(ct,tags_b[column]);push!(cb,bound)
                    elem_b===nothing || push!(ce,elem_b[column])
                end
            end
            push!(cblocks,ElementBlock(2,reshape(cnodes,3,:),cphys))
            push!(cowners,co);push!(ctags,ct);push!(cbind,cb)
            celem===nothing || push!(celem,ce)
        end
        new_data=cdata===nothing ? nothing :
            Elements._native_mixed_entity_data(cdata.entities;
                node_entities=cdata.node_entities,
                node_parametric=cdata.node_parametric,
                external_node_tags=cdata.external_node_tags,
                block_entities=cowners,external_element_tags=ctags)
        folded_carrier=MixedMesh(carrier.coords,cblocks;
            physical_names=carrier.physical_names,entity_data=new_data,
            elementary_entities=celem,
            periodic_links=carrier.periodic_links,
            ancillary_sections=carrier.ancillary_sections,
            data_sections=carrier.data_sections,
            partition_data=carrier.partition_data)
        folded_bindings=cbind
    end
    base=_model_mesh_authority(folded_mesh,authority.done;caller,
        private_rank=authority.private_rank,
        carrier=folded_carrier,carrier_rows=authority.carrier_rows,
        carrier_cells=folded_bindings)
    return _ModelMeshAuthority(base.mesh,base.products,base.done,
        base.private_rank,authority.native_inputs,base.carrier,base.carrier_rows,
        authority.producer_rows,base.carrier_cells,folded_inputs)
end
