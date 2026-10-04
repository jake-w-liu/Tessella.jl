# Actual source data for one native rectangular Quad disk. The planner keeps
# these original primary-column IDs and cyclic cell frames unchanged.
struct _ExtrudeNoNewRectGridSource
    source_cells::Vector{NTuple{4,Int32}}
    source_coordinates::Vector{NTuple{3,Float64}}
    boundary_cycle::Vector{Int32}
    boundary_curves::NTuple{4,Int}
    boundary_chains::NTuple{4,Vector{Int32}}
    boundary_widths::NTuple{4,Int}
    edges::Vector{NTuple{2,Int32}}
    edge_indices::Vector{NTuple{4,Int32}}
    edge_cells::Vector{NTuple{2,Int32}}
    edge_sides::Vector{NTuple{2,UInt8}}
    edge_curve::Vector{UInt8}
    edge_curve_direction::Vector{Int8}
    boundary_vertices::BitVector
    boundary_masks::Vector{UInt8}
    category::Vector{UInt8}
    grid_shape::NTuple{2,Int}
    normal_axis::Int8
    normal_direction::Int8
    source_orientation::Int8
end

@noinline function _extrude_nonew_rect_grid_error(caller,reason)
    throw(ArgumentError("$caller: QuadTriNoNewVerts rectangular quadrangle grid $reason"))
end

@noinline _extrude_nonew_rect_grid_source_error(caller,reason)=
    _extrude_nonew_rect_grid_error(caller,"source $reason")

function _extrude_nonew_rect_grid_sizes(a::Int,b::Int,caller)
    a>=2 && b>=2 || _extrude_nonew_rect_grid_source_error(caller,
        "requires at least two actual quadrangles in both directions")
    counts=try
        vertices=Base.checked_mul(Base.checked_add(a,1),Base.checked_add(b,1))
        cells=Base.checked_mul(a,b)
        boundary=Base.checked_mul(2,Base.checked_add(a,b))
        edges=Base.checked_add(Base.checked_mul(2,cells),Base.checked_add(a,b))
        (vertices,cells,edges,boundary)
    catch err
        err isa OverflowError || rethrow()
        _extrude_nonew_rect_grid_source_error(caller,"size arithmetic overflows Int")
    end
    counts[1]<=_EXTRUDE_NONEW_MAX_NODES && counts[3]<=typemax(Int32) ||
        _extrude_nonew_rect_grid_source_error(caller,"exceeds the checked source node or edge limit")
    return counts
end

# A candidate check precedes source grading/allocation. Bounded categories keep
# their earlier dispatcher priority. The constructor independently rechecks it.
function _extrude_nonew_rect_grid_shape(m::GeoModel,source::Int,sides::Int,caller)
    sides==4 && haskey(m.meshing.transfinite_surfaces,source) &&
        _model_surface_recombined(m,source) || return nothing
    loops=m.surfaces[source]
    length(loops)==1 || return nothing
    signed=m.loops[only(loops)]
    length(signed)==4 && allunique(abs.(signed)) || return nothing
    widths=ntuple(4) do position
        curve=abs(signed[position])
        haskey(m.curves,curve) && _curve_type(m,curve)===:line &&
            _occ_geometry(m,curve)===nothing && !(curve in m.meshing.degenerated) || return 0
        control=get(m.meshing.transfinite_curves,curve,nothing)
        control===nothing && return 0
        _flexible_transfinite_nodes(m,control.num_nodes,curve,caller)
    end
    widths[1]>=3 && widths[2]>=3 && widths[1]==widths[3] && widths[2]==widths[4] ||
        return nothing
    shape=(widths[1]-1,widths[2]-1)
    _extrude_nonew_rect_grid_sizes(shape[1],shape[2],caller)
    return shape
end

function _extrude_nonew_rect_grid_axis(spec,caller)
    spec.type===:translate || _extrude_nonew_rect_grid_source_error(caller,
        "requires axis-normal translation; other transforms require the source-grid planner")
    all(isfinite,spec.T) && count(!iszero,spec.T)==1 ||
        _extrude_nonew_rect_grid_source_error(caller,"requires one finite nonzero translation component")
    axis=spec.T[1]!=0 ? 1 : spec.T[2]!=0 ? 2 : 3
    return axis,spec.T[axis]>0 ? 1 : -1
end

@inline _extrude_nonew_rect_grid_key(p)=
    (p[1]==0 ? 0.0 : p[1],p[2]==0 ? 0.0 : p[2],p[3]==0 ? 0.0 : p[3])

@inline function _extrude_nonew_rect_grid_assign!(positions,node,value,caller)
    previous=positions[node]
    (previous==(-1,-1) || previous==value) || _extrude_nonew_rect_grid_source_error(caller,
        "has inconsistent rectangular node incidence")
    positions[node]=value
    return nothing
end

# Boundary-seeded lattice propagation proves a complete regular rectangular
# disk. Euler counts alone would not exclude pinched/nonmanifold complexes.
# The temporary integer labels never reorder source nodes or cell frames.
function _extrude_nonew_rect_grid_topology(cells,edges,indices,edge_cells,edge_sides,
        chains,signed,a::Int,b::Int,caller)
    vertices=(a+1)*(b+1)
    positions=fill((-1,-1),vertices)
    for side in 1:4
        chain=chains[side];width=length(chain)
        for position in 1:width
            node=chain[signed[side]>0 ? position : width-position+1]
            value=side==1 ? (position-1,0) : side==2 ? (a,position-1) :
                side==3 ? (a-position+1,b) : (0,b-position+1)
            _extrude_nonew_rect_grid_assign!(positions,node,value,caller)
        end
    end
    first_chain=chains[1];last_chain=chains[4]
    corner=first_chain[signed[1]>0 ? 1 : end]
    next_node=first_chain[signed[1]>0 ? 2 : end-1]
    previous_node=last_chain[signed[4]>0 ? end-1 : 2]
    seed=0
    for number in eachindex(cells)
        cell=cells[number]
        if corner in cell && next_node in cell && previous_node in cell
            seed==0 || _extrude_nonew_rect_grid_source_error(caller,"duplicates the corner quadrangle")
            seed=number
        end
    end
    seed>0 || _extrude_nonew_rect_grid_source_error(caller,"has no complete corner quadrangle")
    for node in cells[seed]
        node in (corner,next_node,previous_node) ||
            _extrude_nonew_rect_grid_assign!(positions,node,(1,1),caller)
    end
    seed_cell=cells[seed]
    x1,y1=positions[seed_cell[1]];x2,y2=positions[seed_cell[2]]
    x3,y3=positions[seed_cell[3]]
    sign=(x2-x1)*(y3-y2)-(y2-y1)*(x3-x2)
    abs(sign)==1 || _extrude_nonew_rect_grid_source_error(caller,"has a nonrectangular corner frame")
    queue=Vector{Int32}(undef,length(cells));queue[1]=Int32(seed)
    visited=falses(length(cells));visited[seed]=true
    used=1;head=1
    while head<=used
        number=Int(queue[head]);head+=1
        cell=cells[number]
        for side in 1:4
            first=positions[cell[side]];second=positions[cell[mod1(side+1,4)]]
            dx=second[1]-first[1];dy=second[2]-first[2]
            abs(dx)+abs(dy)==1 || _extrude_nonew_rect_grid_source_error(caller,
                "has a nonunit rectangular cell edge")
            next_position=positions[cell[mod1(side+2,4)]]
            cross=dx*(next_position[2]-second[2])-dy*(next_position[1]-second[1])
            cross==sign || _extrude_nonew_rect_grid_source_error(caller,
                "has an inconsistent rectangular cyclic frame")
            edge=Int(indices[number][side]);owners=edge_cells[edge]
            neighbor=owners[1]==number ? Int(owners[2]) : Int(owners[1])
            neighbor==0 && continue
            other_side=Int(owners[1]==neighbor ? edge_sides[edge][1] : edge_sides[edge][2])
            other=cells[neighbor]
            first_other=positions[other[other_side]]
            second_other=positions[other[mod1(other_side+1,4)]]
            du=second_other[1]-first_other[1];dv=second_other[2]-first_other[2]
            offset=(-sign*dv,sign*du)
            _extrude_nonew_rect_grid_assign!(positions,other[mod1(other_side+2,4)],
                (second_other[1]+offset[1],second_other[2]+offset[2]),caller)
            _extrude_nonew_rect_grid_assign!(positions,other[mod1(other_side+3,4)],
                (first_other[1]+offset[1],first_other[2]+offset[2]),caller)
            if !visited[neighbor]
                used+=1;queue[used]=Int32(neighbor);visited[neighbor]=true
            end
        end
    end
    used==length(cells) || _extrude_nonew_rect_grid_source_error(caller,
        "is disconnected or has an incomplete rectangular cell complex")
    node_slots=zeros(Int32,vertices)
    for node in eachindex(positions)
        x,y=positions[node]
        0<=x<=a && 0<=y<=b || _extrude_nonew_rect_grid_source_error(caller,
            "has an unused or foreign rectangular node")
        slot=1+x+(a+1)*y
        node_slots[slot]==0 || _extrude_nonew_rect_grid_source_error(caller,
            "identifies distinct nodes at the same rectangular slot")
        node_slots[slot]=Int32(node)
    end
    all(!iszero,node_slots) || _extrude_nonew_rect_grid_source_error(caller,
        "omits a rectangular node slot")
    cell_slots=zeros(Int32,length(cells))
    for number in eachindex(cells)
        cell=cells[number]
        x=min(positions[cell[1]][1],positions[cell[2]][1],positions[cell[3]][1],positions[cell[4]][1])
        y=min(positions[cell[1]][2],positions[cell[2]][2],positions[cell[3]][2],positions[cell[4]][2])
        0<=x<a && 0<=y<b || _extrude_nonew_rect_grid_source_error(caller,
            "has a foreign rectangular cell")
        slot=1+x+a*y
        cell_slots[slot]==0 || _extrude_nonew_rect_grid_source_error(caller,
            "duplicates a rectangular cell slot")
        cell_slots[slot]=Int32(number)
    end
    all(!iszero,cell_slots) || _extrude_nonew_rect_grid_source_error(caller,
        "omits a rectangular cell slot")
    return nothing
end

struct _ExtrudeNoNewRectBoundaryNode
    bounds::NTuple{4,Float64}
    first::Int32
    last::Int32
    left::Int32
    right::Int32
end

function _extrude_nonew_rect_grid_boundary_tree!(nodes,planar,boundary,first::Int,last::Int)
    index=length(nodes)+1
    if first==last
        a=planar[boundary[first]];b=planar[boundary[mod1(first+1,length(boundary))]]
        box=(min(a[1],b[1]),min(a[2],b[2]),max(a[1],b[1]),max(a[2],b[2]))
        push!(nodes,_ExtrudeNoNewRectBoundaryNode(box,Int32(first),Int32(last),0,0))
    else
        push!(nodes,_ExtrudeNoNewRectBoundaryNode((0.0,0.0,0.0,0.0),Int32(first),Int32(last),0,0))
        middle=first+(last-first)÷2
        left=_extrude_nonew_rect_grid_boundary_tree!(nodes,planar,boundary,first,middle)
        right=_extrude_nonew_rect_grid_boundary_tree!(nodes,planar,boundary,middle+1,last)
        a=nodes[left].bounds;b=nodes[right].bounds
        box=(min(a[1],b[1]),min(a[2],b[2]),max(a[3],b[3]),max(a[4],b[4]))
        nodes[index]=_ExtrudeNoNewRectBoundaryNode(box,Int32(first),Int32(last),left,right)
    end
    return Int32(index)
end

function _extrude_nonew_rect_grid_boundary_pair(nodes,first::Int32,second::Int32,
        planar,boundary,budget,caller)
    budget[]-=1
    budget[]>=0 || _extrude_nonew_rect_grid_source_error(caller,
        "actual boundary simplicity exceeds the bounded spatial certification budget")
    a=nodes[first];b=nodes[second]
    x=a.bounds;y=b.bounds
    (x[3]<y[1] || y[3]<x[1] || x[4]<y[2] || y[4]<x[2]) && return nothing
    if first==second
        a.left==0 && return nothing
        _extrude_nonew_rect_grid_boundary_pair(nodes,a.left,a.left,planar,boundary,budget,caller)
        _extrude_nonew_rect_grid_boundary_pair(nodes,a.right,a.right,planar,boundary,budget,caller)
        _extrude_nonew_rect_grid_boundary_pair(nodes,a.left,a.right,planar,boundary,budget,caller)
    elseif a.left==0 && b.left==0
        count=length(boundary)
        u=boundary[a.first];v=boundary[mod1(Int(a.first)+1,count)]
        w=boundary[b.first];z=boundary[mod1(Int(b.first)+1,count)]
        (u==w || u==z || v==w || v==z) && return nothing
        !_extrude_nonew_quad_patch_crossing(planar[u],planar[v],planar[w],planar[z]) ||
            _extrude_nonew_rect_grid_source_error(caller,
                "actual boundary has a nonadjacent crossing, touch or overlap")
    elseif b.left==0 || (a.left!=0 && a.last-a.first>=b.last-b.first)
        _extrude_nonew_rect_grid_boundary_pair(nodes,a.left,second,planar,boundary,budget,caller)
        _extrude_nonew_rect_grid_boundary_pair(nodes,a.right,second,planar,boundary,budget,caller)
    else
        _extrude_nonew_rect_grid_boundary_pair(nodes,first,b.left,planar,boundary,budget,caller)
        _extrude_nonew_rect_grid_boundary_pair(nodes,first,b.right,planar,boundary,budget,caller)
    end
    return nothing
end

function _extrude_nonew_rect_grid_boundary_certify(planar,boundary,caller)
    count=length(boundary)
    for position in eachindex(boundary)
        a=planar[boundary[mod1(position-1,count)]]
        b=planar[boundary[position]];c=planar[boundary[mod1(position+1,count)]]
        # Adjacent collinear samples may continue forward, but cannot backtrack
        # into an overlapping segment. Distinctness excludes zero-length edges.
        orient2(a,b,c)!=0 || _extrude_nonew_quad_patch_on_segment(a,c,b) ||
            _extrude_nonew_rect_grid_source_error(caller,"actual boundary backtracks or overlaps")
    end
    nodes=_ExtrudeNoNewRectBoundaryNode[]
    sizehint!(nodes,2count-1)
    root=_extrude_nonew_rect_grid_boundary_tree!(nodes,planar,boundary,1,count)
    depth=1;remaining=count
    while remaining>1
        remaining=(remaining+1)÷2;depth+=1
    end
    budget=Ref(Base.checked_mul(Base.checked_mul(128,count),depth))
    _extrude_nonew_rect_grid_boundary_pair(nodes,root,root,planar,boundary,budget,caller)
    return nothing
end

function _extrude_nonew_rect_grid_source(m::GeoModel,source::Int,source_mesh,spec,caller)
    shape=_extrude_nonew_rect_grid_shape(m,source,4,caller)
    shape!==nothing || _extrude_nonew_rect_grid_source_error(caller,
        "requires four native straight Curve chains with equal opposite effective counts >=3")
    a,b=shape
    vertices,cell_count,edge_count,boundary_count=_extrude_nonew_rect_grid_sizes(a,b,caller)
    source_mesh isa MixedMesh && nnodes(source_mesh)==vertices && length(source_mesh.blocks)==1 ||
        _extrude_nonew_rect_grid_source_error(caller,"requires the complete actual linear Quad grid without extra nodes")
    block=only(source_mesh.blocks)
    block isa ElementBlock && block.msh==3 && size(block.nodes)==(4,cell_count) ||
        _extrude_nonew_rect_grid_source_error(caller,"requires only the complete actual linear Quad cells")
    _surface_type(m,source) in (:plane,:ruled) && !haskey(m.surface_geometry,source) ||
        _extrude_nonew_rect_grid_source_error(caller,"requires native planar filling without auxiliary surface geometry")
    axis,direction=_extrude_nonew_rect_grid_axis(spec,caller)
    coordinates=Vector{NTuple{3,Float64}}(undef,vertices)
    lookup=Dict{NTuple{3,Float64},Int32}();sizehint!(lookup,vertices)
    for node in 1:vertices
        point=(source_mesh.coords[1,node],source_mesh.coords[2,node],source_mesh.coords[3,node])
        all(isfinite,point) || _extrude_nonew_rect_grid_source_error(caller,"contains a nonfinite actual source node")
        key=_extrude_nonew_rect_grid_key(point)
        !haskey(lookup,key) || _extrude_nonew_rect_grid_source_error(caller,
            "contains physically coincident distinct source nodes")
        coordinates[node]=point;lookup[key]=Int32(node)
    end
    height=coordinates[1][axis]
    all(point->point[axis]==height,coordinates) || _extrude_nonew_rect_grid_source_error(caller,
        "is not in an exact axis plane normal to the translation")
    signed=m.loops[only(m.surfaces[source])]
    curves=ntuple(position->abs(signed[position]),4)
    widths=(a+1,b+1,a+1,b+1)
    chains=ntuple(4) do number
        curve=curves[number];first_point,last_point=m.curves[curve]
        first_point!=last_point || _extrude_nonew_rect_grid_source_error(caller,"has a closed boundary Curve")
        actual=_extrude_curve_nodes(m,curve,caller)
        length(actual)==widths[number] &&
            _extrude_nonew_rect_grid_key(actual[1])==_extrude_nonew_rect_grid_key(m.points[first_point]) &&
            _extrude_nonew_rect_grid_key(actual[end])==_extrude_nonew_rect_grid_key(m.points[last_point]) ||
            _extrude_nonew_rect_grid_source_error(caller,"boundary Curve chains do not retain their native nodes")
        chain=Vector{Int32}(undef,widths[number])
        for position in eachindex(chain)
            node=get(lookup,_extrude_nonew_rect_grid_key(actual[position]),Int32(0))
            node>0 || _extrude_nonew_rect_grid_source_error(caller,
                "actual native Curve chain is absent from the retained source mesh")
            chain[position]=node
        end
        allunique(chain) || _extrude_nonew_rect_grid_source_error(caller,"has repeated native boundary samples")
        chain
    end
    boundary=Vector{Int32}(undef,boundary_count)
    curve_edges=Dict{NTuple{2,Int32},Tuple{UInt8,Int8}}();sizehint!(curve_edges,boundary_count)
    position=0
    starts=ntuple(number->chains[number][signed[number]>0 ? 1 : end],4)
    for number in 1:4
        chain=chains[number];width=widths[number]
        chain[signed[number]>0 ? end : 1]==starts[mod1(number+1,4)] ||
            _extrude_nonew_rect_grid_source_error(caller,"has an unclosed four-Curve boundary")
        for local_position in 1:width-1
            u=chain[local_position];v=chain[local_position+1]
            edge=minmax(u,v)
            !haskey(curve_edges,edge) || _extrude_nonew_rect_grid_source_error(caller,"repeats a boundary edge")
            curve_edges[edge]=(UInt8(number),Int8(u<v ? 1 : -1))
            position+=1
            boundary[position]=chain[signed[number]>0 ? local_position : width-local_position+1]
        end
    end
    allunique(boundary) || _extrude_nonew_rect_grid_source_error(caller,"repeats a boundary vertex")
    planar=Vector{NTuple{2,Float64}}(undef,vertices)
    for node in 1:vertices
        planar[node]=_extrude_nonew_quad_patch_plane(coordinates[node],axis)
    end
    corner_orientation=orient2(planar[starts[1]],planar[starts[2]],planar[starts[3]])
    corner_orientation!=0 || _extrude_nonew_rect_grid_source_error(caller,"has a collinear CAD corner turn")
    for corner in 1:4
        orient2(planar[starts[corner]],planar[starts[mod1(corner+1,4)]],
            planar[starts[mod1(corner+2,4)]])==corner_orientation ||
            _extrude_nonew_rect_grid_source_error(caller,"requires strictly convex four-Curve CAD corners")
    end
    cells=Vector{NTuple{4,Int32}}(undef,cell_count)
    indices=Vector{NTuple{4,Int32}}(undef,cell_count)
    edges=NTuple{2,Int32}[];sizehint!(edges,edge_count)
    edge_cells=NTuple{2,Int32}[];sizehint!(edge_cells,edge_count)
    edge_sides=NTuple{2,UInt8}[];sizehint!(edge_sides,edge_count)
    edge_lookup=Dict{NTuple{2,Int32},Int32}();sizehint!(edge_lookup,edge_count)
    orientation=0
    for number in 1:cell_count
        cell=(block.nodes[1,number],block.nodes[2,number],block.nodes[3,number],block.nodes[4,number])
        all(node->1<=node<=vertices,cell) && allunique(cell) ||
            _extrude_nonew_rect_grid_source_error(caller,"has repeated or foreign Quad nodes")
        cells[number]=cell
        for corner in 1:4
            turn=orient2(planar[cell[corner]],planar[cell[mod1(corner+1,4)]],planar[cell[mod1(corner+2,4)]])
            orientation==0 && (orientation=turn)
            turn!=0 && turn==orientation || _extrude_nonew_rect_grid_source_error(caller,
                "actual Quads are not strictly convex with coherent winding")
        end
        local_edges=ntuple(4) do side
            u,v=_extrude_nonew_quad_patch_edge(cell,side);key=minmax(u,v)
            index=get(edge_lookup,key,Int32(0))
            if index==0
                length(edges)<edge_count || _extrude_nonew_rect_grid_source_error(caller,"has too many actual source edges")
                push!(edges,key);index=Int32(length(edges));edge_lookup[key]=index
                push!(edge_cells,(Int32(number),Int32(0)));push!(edge_sides,(UInt8(side),UInt8(0)))
            else
                owners=edge_cells[index];sides=edge_sides[index]
                owners[2]==0 || _extrude_nonew_rect_grid_source_error(caller,"has a nonmanifold source edge")
                first_u,first_v=_extrude_nonew_quad_patch_edge(cells[owners[1]],Int(sides[1]))
                (u,v)==(first_v,first_u) || _extrude_nonew_rect_grid_source_error(caller,
                    "internal source edges require opposite incidences")
                edge_cells[index]=(owners[1],Int32(number));edge_sides[index]=(sides[1],UInt8(side))
            end
            index
        end
        indices[number]=local_edges
    end
    length(edges)==edge_count || _extrude_nonew_rect_grid_source_error(caller,"has an incomplete source-edge graph")
    edge_curve=zeros(UInt8,edge_count);edge_direction=zeros(Int8,edge_count)
    boundary_vertices=falses(vertices);exterior_count=0
    for index in eachindex(edges)
        curve_info=get(curve_edges,edges[index],(UInt8(0),Int8(0)))
        exterior=edge_cells[index][2]==0
        exterior==(curve_info[1]!=0) || _extrude_nonew_rect_grid_source_error(caller,
            "source exterior does not equal the four actual Curve chains")
        if exterior
            exterior_count+=1
            edge_curve[index]=curve_info[1];edge_direction[index]=curve_info[2]
            u,v=edges[index];boundary_vertices[u]=true;boundary_vertices[v]=true
        end
    end
    exterior_count==boundary_count || _extrude_nonew_rect_grid_source_error(caller,"has an incomplete native boundary")
    _extrude_nonew_rect_grid_topology(cells,edges,indices,edge_cells,edge_sides,
        chains,signed,a,b,caller)
    # A regular disk of coherent strictly positive triangle splits with a
    # bijective Jordan boundary is globally embedded (degree/global inversion).
    # Exact boundary contacts, not area equality alone, are indispensable.
    _extrude_nonew_rect_grid_boundary_certify(planar,boundary,caller)
    masks=Vector{UInt8}(undef,cell_count);category=Vector{UInt8}(undef,cell_count)
    for number in eachindex(cells)
        cell=cells[number];mask=UInt8(0)
        for corner in 1:4
            boundary_vertices[cell[corner]] && (mask|=UInt8(1)<<(corner-1))
        end
        count=count_ones(mask)
        count in (0,2,3) || _extrude_nonew_rect_grid_source_error(caller,
            "requires only B0, adjacent-B2 and B3 rectangular source cells")
        if count==2
            adjacent=0
            for corner in 1:4
                adjacent+=boundary_vertices[cell[corner]] && boundary_vertices[cell[mod1(corner+1,4)]]
            end
            adjacent==1 || _extrude_nonew_rect_grid_source_error(caller,"has an opposite-boundary B2 cell")
        end
        masks[number]=mask;category[number]=UInt8(count)
    end
    return _ExtrudeNoNewRectGridSource(cells,coordinates,boundary,curves,chains,widths,
        edges,indices,edge_cells,edge_sides,edge_curve,edge_direction,boundary_vertices,
        masks,category,shape,Int8(axis),Int8(direction),Int8(orientation))
end
