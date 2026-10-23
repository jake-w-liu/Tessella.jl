using ..SizeField: AbstractAnisoField, metric_edge_length,
                   isotropic_metric, BoundedSize, GMSH_MAX_SIZE
import ..SizeField: metric_at
using ..PipelineSupport: pipeline_float, pipeline_limit, PIPELINE_DEFAULT_MAX_NODES
using ..Mesh3D: _EHeap, _hpush!, _hpop!

@inline _metric_surface_edges(cell)=((cell[1],cell[2]),(cell[2],cell[3]),(cell[3],cell[1]))

function _metric_surface_incidence(cells,caller)
    incidence=Dict{NTuple{2,Int32},Vector{Int}}()
    for (column,cell) in enumerate(cells),edge in _metric_surface_edges(cell)
        push!(get!(()->Int[],incidence,minmax(edge...)),column)
    end
    all(hits->length(hits) in (1,2),values(incidence)) || throw(ArgumentError(
        "$caller: initial Surface has nonmanifold edge incidence"))
    return incidence
end

function _metric_surface_birth(high,count,caller,label)
    UInt64(count)<=typemax(UInt64)-high || throw(ArgumentError(
        "$caller: $label births exhaust UInt64 labels"))
    return high+UInt64(count)
end

function _metric_surface_admit(m,tag,caller)
    haskey(m.surfaces,tag) || throw(ArgumentError("$caller: unknown Surface[$tag]"))
    _surface_type(m,tag)===:plane || throw(ArgumentError("$caller: requires a planar Surface"))
    haskey(m.meshing.transfinite_surfaces,tag) && throw(ArgumentError(
        "$caller: transfinite Surface construction is not supported"))
    _extrude_gate(m,2,tag)===nothing || throw(ArgumentError(
        "$caller: extruded Surface construction is not supported"))
    _model_surface_recombined(m,tag) && throw(ArgumentError(
        "$caller: recombination is not supported"))
    isempty(get(m.embeds,(2,tag),NTuple{2,Int}[])) || throw(ArgumentError(
        "$caller: embedded Surface constraints are not supported"))
    get(m.meshing.smoothing,(2,tag),0)==0 || throw(ArgumentError(
        "$caller: Surface smoothing is not supported"))
    curves=_gmsh_native_surface_curves(m,tag)
    all(curve->_curve_type(m,curve)===:line,curves) || throw(ArgumentError(
        "$caller: requires straight boundary Curves"))
    curve_set=Set(curves)
    any(constraint->constraint.dim==2 &&
        (abs(Int(constraint.slave_entity))==tag || abs(Int(constraint.master_entity))==tag) ||
        constraint.dim==1 && (abs(Int(constraint.slave_entity)) in curve_set ||
        abs(Int(constraint.master_entity)) in curve_set),model_periodic_constraints(m)) &&
        throw(ArgumentError("$caller: periodic Surface or boundary Curve construction is not supported"))
    _model_surface_plane(m,tag,caller)
    return nothing
end

# A stored Curve discretization is stale for this producer when any edge it
# would emit fails the same strict certificate the retained boundary is
# checked against below: the world chord of every consecutive parameter pair,
# measured under `field` in this Surface's entity context, must not exceed
# `target`.
function _metric_boundary_curve_stale(m,curve,params,tag,field,target,caller)
    γ=u->_model_curve_point(m,curve,u,caller)
    previous=γ(first(params))
    for u in Iterators.drop(params,1)
        point=γ(u)
        metric_edge_length(field,previous,point;entity=(2,tag))<=target ||
            return true
        previous=point
    end
    return false
end

function _metric_surface_refine(m,tag,field,target,node_limit,triangle_limit,
                                input,options,node_high_water,element_high_water,caller)
    # Boundary cells that are lazily graded while building the lower input
    # must satisfy this producer's strict `metric_edge_length` certificate
    # against `field`. The curve graders only see `options.inner_field` and
    # the `F_Lc_aniso` dispatch on `m.meshing.algorithm_2d`, so transport the
    # certificate requirement into the grading context: bind the resolved
    # field as the background field and select the BAMG integrand for the
    # duration of the lower-dimensional pass. `m` is the working model; the
    # stored flag is restored afterwards so unrelated surfaces keep their
    # caller's resolution.
    options=_mesh1d_options_with_inner_field(options,field)
    saved_algorithm=m.meshing.algorithm_2d
    m.meshing.algorithm_2d=7
    local authority, limit
    try
        # A boundary Curve may already carry parameters graded under a
        # different sizing context — e.g. `mesh_model_surface`'s
        # `Mesh1D`-before-`Mesh2D` front pass (`_attribute_forced_parameters`)
        # grades lazily under the model's scalar options before this
        # producer's dispatch, and `reuse_curve_params` then keeps them
        # verbatim. A retained boundary cell must already satisfy the strict
        # certificate below, so re-grade any stored boundary discretization
        # whose edges violate it; the lower pass then derives it under the
        # transported metric. Periodic dependencies of a stale curve are
        # dropped too — a slave's params are only a copy of its master's.
        input===nothing && for loop_id in m.surfaces[tag],signed in m.loops[loop_id]
            curve=abs(signed)
            params=get(m.curve_params,curve,nothing)
            params===nothing && continue
            _metric_boundary_curve_stale(m,curve,params,tag,field,target,
                                         caller) || continue
            delete!(m.curve_params,curve)
            for dep in _model_curve_mesh_dependencies(m,curve)
                delete!(m.curve_params,dep)
            end
        end
        authority,limit=_gmsh_native_surface_lower_input(m,options,input,
            node_high_water,caller)
    finally
        m.meshing.algorithm_2d=saved_algorithm
    end
    frame=_gmsh_native_surface_frame(m,tag,options,caller)
    source=_gmsh_native_surface_source(m,tag,authority,frame,caller;ctx_lc=options.ctx_lc)
    length(source.tags)<=node_limit || throw(ArgumentError(
        "$caller: retained lower vertices exceed max_nodes"))
    # Without interior constraints, Euler's disk-with-holes relation fixes
    # the initial triangle count before its triangulation is allocated.
    initial_count=length(source.tags)+2(length(m.surfaces[tag])-1)-2
    initial_count<=triangle_limit || throw(ArgumentError(
        "$caller: initial Surface exceeds max_triangles"))
    initial=_gmsh_native_initial_surface(m,tag,options,authority,limit,caller;
        algorithm=3,element_high_water)
    mesh=initial.mesh
    nnodes(mesh)<=node_limit && ntris(mesh)<=triangle_limit || throw(ArgumentError(
        "$caller: initial lower mesh exceeds max_nodes or max_triangles"))
    validate(mesh).ok || throw(ArgumentError("$caller: initial Surface is invalid"))
    coordinates=NTuple{3,Float64}[Tuple(mesh.coords[:,row]) for row in axes(mesh.coords,2)]
    parameters=NTuple{2,Float64}[_gmsh_native_surface_parameters(frame,point) for point in coordinates]
    cells=NTuple{3,Int32}[Tuple(cell) for cell in eachcol(mesh.tris)]
    incidence=_metric_surface_incidence(cells,caller)
    boundary=Set(edge for (edge,hits) in incidence if length(hits)==1)
    original_boundary=Set(minmax(initial.node_tags[edge[1]],initial.node_tags[edge[2]]) for edge in boundary)
    wanted_boundary=Set(minmax(source.tags[edge[1]],source.tags[edge[2]]) for edge in source.boundaries)
    original_boundary==wanted_boundary || throw(ArgumentError(
        "$caller: initial Surface loses its actual lower Curve-cell incidence"))
    score(edge)=metric_edge_length(field,coordinates[edge[1]],coordinates[edge[2]];entity=(2,tag))
    all(edge->score(edge)<=target,boundary) || throw(ArgumentError(
        "$caller: retained boundary exceeds max_metric_edge; provide finer lower Curve cells"))
    labels=copy(initial.element_tags);numbers=copy(initial.node_tags)
    owners=copy(initial.node_owners);node_parameters=deepcopy(initial.node_parameters)
    retained=copy(initial.retained_nodes);birth_rows=copy(initial.birth_rows)
    alive=trues(length(cells));live_count=length(cells)
    node_alive=trues(length(coordinates))
    node_max=initial.node_high_water;cell_max=initial.element_high_water
    heap=_EHeap();queued=Set{NTuple{2,Int32}}()
    function enqueue(edge)
        haskey(incidence,edge) && !(edge in queued) || return
        violation=score(edge)
        isfinite(violation) && violation>=0 || throw(ArgumentError(
            "$caller: edge metric must be finite and nonnegative"))
        if violation>target
            edge in boundary && throw(ArgumentError("$caller: retained boundary metric changed"))
            push!(queued,edge);_hpush!(heap,(violation,edge[1],edge[2]))
        end
    end
    for edge in sort!(collect(keys(incidence)));enqueue(edge);end
    positions=Dict{NTuple{3,Float64},Int32}()
    for (index,point) in enumerate(coordinates);positions[point]=Int32(index);end
    # Dissolved vertex positions stay reserved so a later split cannot birth a
    # vertex back onto the exact razor site the repair just removed.
    retired=Set{NTuple{3,Float64}}()
    # Tangent-frame pullback of the world metric at a surface point. In the
    # affine plane coordinates the quadratic form J'MJ gives exactly the
    # `metric_edge_length` of any in-plane chord, so circumcenter and
    # in-circle computations in (u,v) are the world tests on lifted points.
    function uv_metric(point)
        mm=metric_at(field,point[1],point[2],point[3],(2,tag))
        f=frame.first;s=frame.second
        a=mm.m11*f[1]*f[1]+mm.m22*f[2]*f[2]+mm.m33*f[3]*f[3]+
          2*(mm.m12*f[1]*f[2]+mm.m13*f[1]*f[3]+mm.m23*f[2]*f[3])
        b=mm.m11*f[1]*s[1]+mm.m22*f[2]*s[2]+mm.m33*f[3]*s[3]+
          mm.m12*(f[1]*s[2]+f[2]*s[1])+mm.m13*(f[1]*s[3]+f[3]*s[1])+
          mm.m23*(f[2]*s[3]+f[3]*s[2])
        d=mm.m11*s[1]*s[1]+mm.m22*s[2]*s[2]+mm.m33*s[3]*s[3]+
          2*(mm.m12*s[1]*s[2]+mm.m13*s[1]*s[3]+mm.m23*s[2]*s[3])
        return (a,b,d)
    end
    # Metric circumcircle of a live cell under `metric` (a 2x2 pullback).
    # `_gmsh_sys2x2` hands back (0,0) on a degenerate solve, so the returned
    # center is re-certified by equal metric distance to all three corners.
    function uv_circumcircle(column,metric)
        cell=cells[column]
        pa,pb,pc=parameters[cell[1]],parameters[cell[2]],parameters[cell[3]]
        cc,r2=Mesh2D._gmsh_circumcenter_metric(pa,pb,pc,metric)
        (isfinite(cc[1]) && isfinite(cc[2]) && isfinite(r2) && r2>0) ||
            return nothing
        for corner in (pa,pb,pc)
            other=Mesh2D._gmsh_metric_distance2(cc,corner,metric)
            (isfinite(other) && abs(other-r2)<=1e-9*r2) || return nothing
        end
        return cc,r2
    end
    # Visibility walk to the live cell containing `q` in parameter space,
    # starting from cell `seed`. Every live cell carries the same winding `w`;
    # crossing the edge whose orientation disagrees advances the walk. Returns
    # `nothing` when `q` leaves the domain or the walk cannot converge.
    function uv_locate(seed,q,w)
        column=seed;seen=Set{Int}()
        while true
            column in seen && return nothing
            push!(seen,column)
            alive[column] || return nothing
            cell=cells[column]
            # Cross the edge whose Float64 signed area against `q` is most
            # violated — the steepest-descent visibility walk that resists
            # cycling on non-Delaunay intermediate triangulations. A violated
            # boundary edge only means `q` left the domain through it; another
            # violated edge with a live neighbour still advances the walk.
            best=0.0;next_column=0;crossed_boundary=false
            for index in 1:3
                u=cell[index];v=cell[mod1(index+1,3)]
                x1,y1=parameters[u];x2,y2=parameters[v]
                violation=-w*((x2-x1)*(q[2]-y1)-(y2-y1)*(q[1]-x1))
                (isfinite(violation) && violation>0) || continue
                nb=0
                for candidate in get(incidence,minmax(u,v),Int[])
                    candidate!=column && (nb=candidate)
                end
                nb==0 && (crossed_boundary=true;continue)
                violation>best && (best=violation;next_column=nb)
            end
            next_column!=0 && (column=next_column;continue)
            crossed_boundary && return nothing
            return column
        end
    end
    # Metric-Delaunay Steiner insertion. Conforming midpoint bisection can
    # only shorten edges — it can never reorient connectivity onto the metric
    # ellipse — so under strongly anisotropic fields every direction collapses
    # toward h_min and the fill overshoots the prescribed triangle count by
    # tens of times. The upstream remedy is an anisotropic Delaunay
    # construction; here the worst violating edge's incident cells nominate a
    # Steiner point at their metric circumcenter, inserted through an
    # anisotropic Bowyer-Watson cavity under the constant local metric at the
    # point. The commit keeps every certificate of the split path:
    # the perimeter — including every retained boundary edge — survives
    # verbatim as fan edges, each fan cell carries the cavity winding and a
    # positive physical area, no surviving live vertex lands strictly inside
    # a fan cell, cavity-interior vertices are non-retained and dissolve like
    # the razor repair, and the birth is a distinct representable position.
    # Any failed check falls through to the certified midpoint split.
    function _metric_insert_steiner(sides,hits)
        length(sides)==2 || return false
        (f1,l1,o1),(f2,l2,o2)=sides
        o1!=o2 || return false
        (f2,l2)==(l1,f1) || return false
        w=orient2(parameters[f1],parameters[l1],parameters[o1])
        (isfinite(w) && w!=0 &&
         orient2(parameters[f2],parameters[l2],parameters[o2])==w) ||
            return false
        # Candidate insertion points are the incident cells' metric
        # circumcenters, worst first. An outside center is the usual case for
        # a metric-obtuse violating cell — locate the containing cell by a
        # parameter-space visibility walk and let the cavity flood reach the
        # popped edge's cells when the metric says they belong.
        for column in sort!(collect(hits);by=column->begin
                cell=cells[column]
                pa,pb,pc=parameters[cell[1]],parameters[cell[2]],parameters[cell[3]]
                centroid=((pa[1]+pb[1]+pc[1])/3,(pa[2]+pb[2]+pc[2])/3)
                solved=uv_circumcircle(column,
                    uv_metric(_surface_metric_lift(frame,centroid[1],centroid[2])))
                solved===nothing ? 0.0 : solved[2]
            end,rev=true)
            cell=cells[column]
            pa,pb,pc=parameters[cell[1]],parameters[cell[2]],parameters[cell[3]]
            centroid=((pa[1]+pb[1]+pc[1])/3,(pa[2]+pb[2]+pc[2])/3)
            solved=uv_circumcircle(column,
                uv_metric(_surface_metric_lift(frame,centroid[1],centroid[2])))
            solved===nothing && (steiner_fail[:cc]=get(steiner_fail,:cc,0)+1;continue)
            p_uv,_=solved
            p_world=_surface_metric_lift(frame,p_uv[1],p_uv[2])
            (haskey(positions,p_world) || p_world in retired) && (steiner_fail[:dup]=get(steiner_fail,:dup,0)+1;continue)
            home=uv_locate(column,p_uv,w)
            home===nothing && (steiner_fail[:locate]=get(steiner_fail,:locate,0)+1;continue)
            mp=uv_metric(p_world)
            _metric_steiner_commit(p_uv,p_world,home,mp,w) && return true
        end
        return false
    end
    # Commit one Bowyer-Watson insertion: flood the constant-local-metric
    # cavity from `home` through interior edges while neighbours' metric
    # circumcircles contain p, certify the disk and the retriangulated fan,
    # then apply the shared bookkeeping. Boundary edges carry one live cell,
    # so the flood can never leave the domain.
    function _metric_steiner_commit(p_uv,p_world,home,mp,w)
        cavity=Set{Int}(home);stack=Int[home]
        while !isempty(stack) && length(cavity)<=2048
            column=pop!(stack)
            for cav_edge in _metric_surface_edges(cells[column])
                neighbours=get(incidence,minmax(cav_edge...),Int[])
                length(neighbours)==2 || continue
                for neighbour in neighbours
                    (neighbour in cavity || !alive[neighbour]) && continue
                    solved=uv_circumcircle(neighbour,mp)
                    solved===nothing && continue
                    Mesh2D._gmsh_metric_distance2(solved[1],p_uv,mp)<solved[2] ||
                        continue
                    push!(cavity,neighbour);push!(stack,neighbour)
                end
            end
        end
        (isempty(stack) && length(cavity)<=2048) || return false
        # Cavity perimeter: edges carried by exactly one cavity cell. The
        # cavity must be a topological disk — a degree-two single cycle.
        count_in=Dict{NTuple{2,Int32},Int}()
        for column in cavity
            for cav_edge in _metric_surface_edges(cells[column])
                key=minmax(cav_edge...)
                count_in[key]=get(count_in,key,0)+1
            end
        end
        perimeter=[key for (key,c) in count_in if c==1]
        isempty(perimeter) && (steiner_fail[:disk]=get(steiner_fail,:disk,0)+1;return false)
        adjacency=Dict{Int32,Vector{Int32}}()
        for (first,last) in perimeter
            push!(get!(()->Int32[],adjacency,first),last)
            push!(get!(()->Int32[],adjacency,last),first)
        end
        (length(adjacency)==length(perimeter) &&
            all(choices->length(choices)==2,values(adjacency))) || (steiner_fail[:disk]=get(steiner_fail,:disk,0)+1;return false)
        polygon=Int32[]
        start=perimeter[1][1];cur=start;prev=Int32(0)
        for _ in 1:length(perimeter)+1
            push!(polygon,cur)
            choices=adjacency[cur]
            nxt=choices[1]==prev ? choices[2] : choices[1]
            prev,cur=cur,nxt
            cur==start && break
        end
        (cur==start && length(polygon)==length(perimeter)) || (steiner_fail[:disk]=get(steiner_fail,:disk,0)+1;return false)
        # Orient the walk with the cavity winding, then p must see every
        # perimeter edge from the inside — the fan star-shaped certificate.
        shoelace=sum(parameters[polygon[index]][1]*
                     parameters[polygon[mod1(index+1,length(polygon))]][2]-
                     parameters[polygon[index]][2]*
                     parameters[polygon[mod1(index+1,length(polygon))]][1]
                     for index in 1:length(polygon))
        sign(shoelace)==w || reverse!(polygon)
        inset=Set(polygon)
        interior=Set{Int32}()
        for column in cavity
            for vertex in cells[column]
                vertex in inset || push!(interior,vertex)
            end
        end
        all(vertex->node_alive[vertex] && !retained[vertex],interior) ||
            (steiner_fail[:interior]=get(steiner_fail,:interior,0)+1;return false)
        for vertex in polygon
            Mesh2D._gmsh_metric_distance2(p_uv,parameters[vertex],mp)>
                (0.02*target)^2 || (steiner_fail[:close]=get(steiner_fail,:close,0)+1;return false)
        end
        newcells=NTuple{3,Int32}[]
        for index in 1:length(polygon)
            u=polygon[index];v=polygon[mod1(index+1,length(polygon))]
            orient2(parameters[u],parameters[v],p_uv)==w || (steiner_fail[:star]=get(steiner_fail,:star,0)+1;return false)
            area=triangle_area(coordinates[u],coordinates[v],p_world)
            (isfinite(area) && area>0) || (steiner_fail[:star]=get(steiner_fail,:star,0)+1;return false)
            push!(newcells,(u,v,Int32(-1)))
        end
        # No surviving live vertex may land strictly inside a new cell; the
        # cavity's own interior vertices are dying with it.
        near=Set{Int32}(polygon)
        for column in cavity
            for vertex in cells[column];push!(near,vertex);end
        end
        for column in cavity
            for cav_edge in _metric_surface_edges(cells[column])
                for neighbour in get(incidence,minmax(cav_edge...),Int[])
                    alive[neighbour] && !(neighbour in cavity) || continue
                    for vertex in cells[neighbour];push!(near,vertex);end
                end
            end
        end
        setdiff!(near,interior)
        for child in newcells
            c1,c2=child[1],child[2]
            for wtest in near
                (wtest==c1 || wtest==c2) && continue
                pw=parameters[wtest]
                o1=orient2(parameters[c1],parameters[c2],pw)
                o2=orient2(parameters[c2],p_uv,pw)
                o3=orient2(p_uv,parameters[c1],pw)
                (o1==w && o2==w && o3==w) && (steiner_fail[:inside]=get(steiner_fail,:inside,0)+1;return false)
                # Strictly on a fan edge is no better than strictly inside:
                # the new edge would carry a T-junction vertex it cannot
                # express, so reject rather than emit a nonmanifold cell.
                for (o,(s,e)) in ((o1,(parameters[c1],parameters[c2])),
                                  (o2,(parameters[c2],p_uv)),
                                  (o3,(p_uv,parameters[c1])))
                    o==0 || continue
                    (min(s[1],e[1])<=pw[1]<=max(s[1],e[1]) &&
                     min(s[2],e[2])<=pw[2]<=max(s[2],e[2])) && (steiner_fail[:onedge]=get(steiner_fail,:onedge,0)+1;return false)
                end
            end
        end
        (length(coordinates)-length(interior)<node_limit &&
            live_count-length(cavity)+length(newcells)<=triangle_limit) ||
            (steiner_fail[:limit]=get(steiner_fail,:limit,0)+1;return false)
        next_node=_metric_surface_birth(node_max,1,caller,"node")
        next_cell=_metric_surface_birth(cell_max,length(newcells),
                                        caller,"cell")
        row=Int32(length(coordinates)+1)
        push!(coordinates,p_world);push!(parameters,p_uv)
        positions[p_world]=row
        push!(numbers,next_node);push!(owners,(2,Int32(tag)))
        push!(node_parameters,collect(p_uv))
        push!(retained,false);push!(birth_rows,row);push!(node_alive,true)
        node_max=next_node
        for column in cavity
            alive[column]=false
            for old_edge in _metric_surface_edges(cells[column])
                key=minmax(old_edge...)
                filter!(!=(column),incidence[key])
                isempty(incidence[key]) && delete!(incidence,key)
            end
        end
        new_edges=Set{NTuple{2,Int32}}()
        for child in newcells
            born=(child[1],child[2],row)
            push!(cells,born);push!(alive,true)
            cell_max+=UInt64(1);push!(labels,cell_max)
            for child_edge in _metric_surface_edges(born)
                key=minmax(child_edge...)
                push!(get!(()->Int[],incidence,key),length(cells))
                push!(new_edges,key)
            end
        end
        cell_max==next_cell ||
            error("$caller: internal cell birth count mismatch")
        for vertex in interior
            node_alive[vertex]=false
            push!(retired,coordinates[vertex])
            delete!(positions,coordinates[vertex])
        end
        live_count+=length(newcells)-length(cavity)
        for new_edge in sort!(collect(new_edges));enqueue(new_edge);end
        return true
    end
    steiner_budget=Dict{NTuple{2,Int32},Int}()
    steiner_stats=[0,0]  # temporary diagnostic: attempts, commits
    steiner_fail=Dict{Symbol,Int}()
    function _metric_steiner_attempt(sides,hits)
        steiner_stats[1]+=1
        ok=_metric_insert_steiner(sides,hits)
        ok && (steiner_stats[2]+=1)
        return ok
    end
    deferred=Set{NTuple{2,Int32}}()
    edge_deferrals=Dict{NTuple{2,Int32},Int}()
    deferred_signature=Dict{NTuple{2,Int32},NTuple{2,NTuple{3,Int32}}}()
    retried_stars=Set{Tuple{NTuple{2,Int32},NTuple{2,NTuple{3,Int32}}}}()
    edge_star_retries=Dict{NTuple{2,Int32},Int}()
    while !isempty(heap.d) || !isempty(deferred)
        if isempty(heap.d)
            # Letting neighboring refinement change a deferred edge's star is
            # what the bounded deferral below buys; re-queue every stalled
            # edge in one pass so the loop either progresses or exhausts the
            # per-edge budgets and reports the certificate that remains.
            for stalled in sort!(collect(deferred))
                delete!(deferred,stalled);push!(queued,stalled)
                _hpush!(heap,(score(stalled),stalled[1],stalled[2]))
            end
        end
        _,a,b=_hpop!(heap);edge=(a,b);delete!(queued,edge)
        haskey(incidence,edge) || continue
        score(edge)>target || continue
        hits=copy(incidence[edge]);length(hits)==2 || throw(ArgumentError(
            "$caller: metric refinement would split an owned lower boundary"))
        length(coordinates)<node_limit && live_count+length(hits)<=triangle_limit ||
            throw(ArgumentError("$caller: metric refinement exceeds max_nodes or max_triangles"))
        sides=[begin
            cell=cells[column]
            slot=only(index for index in 1:3 if
                minmax(cell[index],cell[mod1(index+1,3)])==edge)
            (cell[slot],cell[mod1(slot+1,3)],cell[mod1(slot+2,3)])
        end for column in hits]
        # Prefer a metric-Delaunay Steiner insertion over splitting: the
        # cavity commit above carries the same certificates while building
        # the connectivity the metric actually prescribes. The per-edge
        # budget keeps an unresolvable edge on the certified split path.
        if get(steiner_budget,edge,0)<2
            steiner_budget[edge]=get(steiner_budget,edge,0)+1
            if _metric_steiner_attempt(sides,hits)
                # The insertion point can land outside this edge's star — the
                # popped edge then survives its pop and must go back on the
                # queue rather than escape the certificate.
                enqueue(edge)
                continue
            end
        end
        # A split point is admissible when every incident cell keeps both
        # exact-orientation children and a representable physical area — the
        # certificate the canonical midpoint is checked against below.
        function admissible(point,point_uv)
            for (first,last,other) in sides
                winding=orient2(parameters[first],parameters[last],parameters[other])
                first_child=orient2(parameters[first],point_uv,parameters[other])
                second_child=orient2(point_uv,parameters[last],parameters[other])
                isfinite(winding) && isfinite(first_child) && isfinite(second_child) &&
                    winding!=0 && sign(first_child)==sign(winding) &&
                    sign(second_child)==sign(winding) || return false
                child_areas=(triangle_area(coordinates[first],point,coordinates[other]),
                             triangle_area(point,coordinates[last],coordinates[other]))
                all(area->isfinite(area) && area>0,child_areas) || return false
            end
            return true
        end
        midpoint=ntuple(axis->_convex_coordinate(coordinates[a][axis],coordinates[b][axis],.5),3)
        uv=ntuple(axis->_convex_coordinate(parameters[a][axis],parameters[b][axis],.5),2)
        splitpoint=midpoint;split_uv=uv
        if haskey(positions,midpoint) || midpoint in retired ||
               !admissible(midpoint,uv)
            # On a near-degenerate Float64 lattice the canonical midpoint can
            # land on an existing vertex or produce razor children even when
            # the parents are nondegenerate — the same representability
            # failure the volumetric kernel's `fallback_split_point`
            # (refine_to_size, Mesh3D.jl) resolves. Search dyadic positions
            # near 1/2 with widening ULP shifts and keep the most balanced
            # admissible candidate; the certificate itself is unchanged.
            lerp(x,y,t)=signbit(x)==signbit(y) ? muladd(t,y-x,x) : (1-t)*x+t*y
            shiftulp(x,k)=begin
                y=x
                if k<0;for _ in 1:-k;y=prevfloat(y);end
                elseif k>0;for _ in 1:k;y=nextfloat(y);end
                end
                y
            end
            pa=coordinates[a];pb=coordinates[b]
            found=nothing;found_score=nothing
            for radius in 0:2,n in 1:31
                n==16&&continue
                t=n/32
                base=(lerp(pa[1],pb[1],t),lerp(pa[2],pb[2],t),lerp(pa[3],pb[3],t))
                for ox in -radius:radius,oy in -radius:radius,oz in -radius:radius
                    max(abs(ox),abs(oy),abs(oz))==radius||continue
                    ((pa[1]==pb[1]&&ox!=0)||(pa[2]==pb[2]&&oy!=0)||
                     (pa[3]==pb[3]&&oz!=0))&&continue
                    point=(shiftulp(base[1],ox),shiftulp(base[2],oy),
                           shiftulp(base[3],oz))
                    all(i->min(pa[i],pb[i])<=point[i]<=max(pa[i],pb[i]),1:3)||continue
                    (point!=pa && point!=pb)||continue
                    haskey(positions,point)&&continue
                    point in retired && continue
                    # Keep the candidate's parameter consistent with its
                    # world point: for a ULP-perturbed point `lerp` would put
                    # the uv exactly back on the edge's parametric segment, so
                    # razor-chain children would stay exactly flat no matter
                    # how the world coordinate is perturbed.  The affine frame
                    # map carries the perturbation into parameter space and is
                    # the same map used for every born vertex.
                    point_uv=_gmsh_native_surface_parameters(frame,point)
                    admissible(point,point_uv)||continue
                    candidate=(abs(n-16),radius,abs(ox)+abs(oy)+abs(oz),ox,oy,oz)
                    (found_score===nothing||candidate<found_score)&&
                        (found=(point,point_uv);found_score=candidate)
                end
            end
            if found===nothing
                # When no representable split point exists — razor parent
                # cells make every on-edge child flat or inverted — a 2-2
                # diagonal flip is the same edge-star retriangulation the
                # volumetric kernel applies (`retriangulate_edge!`): swap the
                # diagonal (a,b) of quad (a,o1,b,o2) for (o1,o2), removing the
                # unresolvable edge outright. The flip is itself certified:
                # both incident cells must share one consistent winding with
                # opposite edge traversal, the new diagonal must not already
                # exist, and both children must pass the same orientation and
                # physical-area certificates as a split. Perimeter edges —
                # including every boundary edge — are retained verbatim.
                function try_flip()
                    length(sides)==2 || return false
                    (f1,l1,o1),(f2,l2,o2)=sides
                    o1!=o2 || return false
                    (f2,l2)==(l1,f1) || return false
                    haskey(incidence,minmax(o1,o2)) && return false
                    w1=orient2(parameters[f1],parameters[l1],parameters[o1])
                    w2=orient2(parameters[f2],parameters[l2],parameters[o2])
                    (isfinite(w1) && isfinite(w2) && w1!=0 && w1==w2) ||
                        return false
                    orient2(parameters[f1],parameters[o2],parameters[o1])==w1 &&
                        orient2(parameters[o2],parameters[l1],parameters[o1])==w1 ||
                        return false
                    area1=triangle_area(coordinates[f1],coordinates[o2],coordinates[o1])
                    area2=triangle_area(coordinates[o2],coordinates[l1],coordinates[o1])
                    (isfinite(area1) && area1>0 &&
                     isfinite(area2) && area2>0) || return false
                    _metric_surface_birth(cell_max,2,caller,"cell")
                    for column in hits
                        alive[column]=false
                        for old_edge in _metric_surface_edges(cells[column])
                            key=minmax(old_edge...)
                            filter!(!=(column),incidence[key])
                            isempty(incidence[key]) && delete!(incidence,key)
                        end
                    end
                    new_edges=Set{NTuple{2,Int32}}()
                    for child in ((f1,o2,o1),(o2,l1,o1))
                        push!(cells,child);push!(alive,true)
                        cell_max+=UInt64(1);push!(labels,cell_max)
                        for child_edge in _metric_surface_edges(child)
                            key=minmax(child_edge...)
                            push!(get!(()->Int[],incidence,key),length(cells))
                            push!(new_edges,key)
                        end
                    end
                    for new_edge in sort!(collect(new_edges));enqueue(new_edge);end
                    return true
                end
                # The razor-chain signature is an encroached edge: a vertex
                # folded onto the segment interior, so every representable
                # midpoint lands inside an existing flat cell's footprint and
                # no diagonal flip exists. The certified repair is splitting
                # the edge at the encroaching vertex itself — the flat cell
                # folded over the edge is absorbed (its two short edges are
                # exactly the subdivided path), while every other side still
                # contributes its pair of children. The retriangulation is
                # fully simulated before commit: children keep the parent
                # winding and a positive physical area, every affected edge
                # ends deleted, on its one retained boundary cell, or
                # manifold with opposite-side apices, and no neighbouring
                # vertex lands strictly inside a child.
                function try_split_at_vertex()
                    length(sides)==2 || return false
                    (f1,l1,o1),(f2,l2,o2)=sides
                    pa=parameters[a];pb=parameters[b]
                    seglen2=sum((coordinates[b][i]-coordinates[a][i])^2
                                for i in 1:3)
                    for v in (o2,o1)
                        (v==a || v==b) && continue
                        pv=parameters[v]
                        orient2(pa,pb,pv)==0 || triangle_area(
                            coordinates[a],coordinates[b],
                            coordinates[v])<=1e-6seglen2 || continue
                        all(i->min(pa[i],pb[i])<=pv[i]<=max(pa[i],pb[i]),
                            1:2) || continue
                        newcells=NTuple{3,Int32}[];ok=true
                        for (first,last,other) in sides
                            # The apex-on-edge cell is the absorbed flat.
                            other==v && continue
                            w0=orient2(parameters[first],parameters[last],
                                       parameters[other])
                            w0!=0 || (ok=false;break)
                            for child in ((first,v,other),(v,last,other))
                                orient2(parameters[child[1]],
                                        parameters[child[2]],
                                        parameters[child[3]])==w0 ||
                                    (ok=false;break)
                                area=triangle_area(coordinates[child[1]],
                                    coordinates[child[2]],coordinates[child[3]])
                                (isfinite(area) && area>0) || (ok=false;break)
                                push!(newcells,child)
                            end
                            ok || break
                        end
                        (ok && !isempty(newcells)) || continue
                        delta=Dict{NTuple{2,Int32},Int}()
                        for column in hits
                            for old_edge in _metric_surface_edges(cells[column])
                                key=minmax(old_edge...)
                                delta[key]=get(delta,key,0)-1
                            end
                        end
                        for child in newcells
                            for child_edge in _metric_surface_edges(child)
                                key=minmax(child_edge...)
                                delta[key]=get(delta,key,0)+1
                            end
                        end
                        all(keys(delta)) do key
                            final=length(get(incidence,key,Int[]))+delta[key]
                            final==0 ||
                                (key in boundary ? final==1 : final==2)
                        end || continue
                        # A reused path edge pairs each new child with one
                        # pre-existing partner; their apices must sit on
                        # opposite sides of the edge (or on it).
                        manifold=true
                        for key in (minmax(a,v),minmax(b,v))
                            haskey(delta,key) || continue
                            for column in get(incidence,key,Int[])
                                column in hits && continue
                                cell=cells[column]
                                apex=only(c for c in cell if
                                          c!=key[1] && c!=key[2])
                                side_old=orient2(parameters[key[1]],
                                    parameters[key[2]],parameters[apex])
                                for child in newcells
                                    (key[1] in child && key[2] in child) ||
                                        continue
                                    apex_new=only(c for c in child if
                                                  c!=key[1] && c!=key[2])
                                    side_new=orient2(parameters[key[1]],
                                        parameters[key[2]],parameters[apex_new])
                                    (side_old==0 || side_new==0 ||
                                     side_old!=side_new) ||
                                        (manifold=false;break)
                                end
                                manifold || break
                            end
                            manifold || break
                        end
                        manifold || continue
                        # No live vertex may land strictly inside a child;
                        # the union of the stars of a, b and v bounds the
                        # overlap check locally.
                        near=Set{Int32}((a,b,v,o1,o2))
                        for index in eachindex(cells,alive)
                            alive[index] || continue
                            cell=cells[index]
                            (cell[1] in near || cell[2] in near ||
                             cell[3] in near) || continue
                            push!(near,cell[1]);push!(near,cell[2])
                            push!(near,cell[3])
                        end
                        inside=false
                        for child in newcells
                            c1,c2,c3=child
                            wc=orient2(parameters[c1],parameters[c2],
                                       parameters[c3])
                            wc!=0 || (inside=true;break)
                            for wtest in near
                                (wtest==c1 || wtest==c2 || wtest==c3) &&
                                    continue
                                orient2(parameters[c1],parameters[c2],
                                        parameters[wtest])==wc &&
                                orient2(parameters[c2],parameters[c3],
                                        parameters[wtest])==wc &&
                                orient2(parameters[c3],parameters[c1],
                                        parameters[wtest])==wc ||
                                    continue
                                inside=true;break
                            end
                            inside && break
                        end
                        inside && continue
                        _metric_surface_birth(cell_max,length(newcells),
                                              caller,"cell")
                        for column in hits
                            alive[column]=false
                            for old_edge in _metric_surface_edges(cells[column])
                                key=minmax(old_edge...)
                                filter!(!=(column),incidence[key])
                                isempty(incidence[key]) &&
                                    delete!(incidence,key)
                            end
                        end
                        new_edges=Set{NTuple{2,Int32}}()
                        for child in newcells
                            push!(cells,child);push!(alive,true)
                            cell_max+=UInt64(1);push!(labels,cell_max)
                            for child_edge in _metric_surface_edges(child)
                                key=minmax(child_edge...)
                                push!(get!(()->Int[],incidence,key),
                                      length(cells))
                                push!(new_edges,key)
                            end
                        end
                        live_count+=length(newcells)-length(hits)
                        for new_edge in sort!(collect(new_edges))
                            enqueue(new_edge)
                        end
                        return true
                    end
                    return false
                end
                # A razor failure is rarely local to one vertex: the rail is
                # a connected cluster of sliver cells, and removing a single
                # vertex just lets refinement re-birth the same site. The
                # certified repair is a cavity dissolution: flood-fill the
                # razor cluster, take one ring of healthy neighbours so the
                # cavity polygon carries real area, then re-triangulate the
                # polygon from its boundary vertices alone. Polygon vertices
                # — including every retained and boundary vertex — survive
                # verbatim; only razor-interior born vertices die, which is
                # exactly what prevents the removal/rebirth cycle. The whole
                # retriangulation is simulated before commit: the cavity must
                # be a topological disk, interior vertices must be live and
                # non-retained, every ear keeps the cavity winding and a
                # positive physical area, no diagonal may belong to a live
                # cell outside the cavity, and no live vertex may land
                # strictly inside a new cell.
                function try_dissolve_cavity()
                    length(sides)==2 || return false
                    cell_area(column)=begin
                        cell=cells[column]
                        triangle_area(coordinates[cell[1]],coordinates[cell[2]],
                                      coordinates[cell[3]])
                    end
                    cell_span(column)=begin
                        cell=cells[column]
                        maximum(sum((coordinates[cell[mod1(index+1,3)]][axis]-
                                     coordinates[cell[index]][axis])^2
                                    for axis in 1:3) for index in 1:3)
                    end
                    cluster=Set{Int}(hits);stack=collect(hits)
                    while !isempty(stack) && length(cluster)<=96
                        column=pop!(stack)
                        for cavity_edge in _metric_surface_edges(cells[column])
                            key=minmax(cavity_edge...)
                            for neighbour in get(incidence,key,Int[])
                                (neighbour in cluster || !alive[neighbour]) &&
                                    continue
                                cell_area(neighbour)<=1e-10cell_span(neighbour) ||
                                    continue
                                push!(cluster,neighbour);push!(stack,neighbour)
                            end
                        end
                    end
                    isempty(stack) || return false
                    # Grow rings until the perimeter polygon carries real
                    # width — a razor cluster's first ring is often still on
                    # the rail — but never past the cap.
                    cavity=Set{Int}(cluster);ring=Set{Int}(cluster)
                    for _ in 1:4
                        next_ring=Set{Int}()
                        for column in ring
                            for cavity_edge in _metric_surface_edges(
                                    cells[column])
                                key=minmax(cavity_edge...)
                                for neighbour in get(incidence,key,Int[])
                                    (neighbour in cavity ||
                                     !alive[neighbour]) && continue
                                    push!(next_ring,neighbour)
                                    push!(cavity,neighbour)
                                end
                            end
                        end
                        isempty(next_ring) && break
                        length(cavity)>192 &&
                            return false
                        ring=next_ring
                        span=maximum(cell_span(column) for column in cluster)
                        width=maximum(cell_area(column) for column in cavity)
                        width>1e-6span && break
                    end
                    length(cavity)<=192 ||
                        return false
                    # Cavity perimeter: edges carried by exactly one cavity
                    # cell (mesh-boundary perimeter edges included).
                    count_in=Dict{NTuple{2,Int32},Int}()
                    for column in cavity
                        for cavity_edge in _metric_surface_edges(cells[column])
                            key=minmax(cavity_edge...)
                            count_in[key]=get(count_in,key,0)+1
                        end
                    end
                    perimeter=[key for (key,hits_in) in count_in if hits_in==1]
                    # A disk cavity's perimeter has degree two everywhere.
                    adjacency=Dict{Int32,Vector{Int32}}()
                    for (first,last) in perimeter
                        push!(get!(()->Int32[],adjacency,first),last)
                        push!(get!(()->Int32[],adjacency,last),first)
                    end
                    (length(adjacency)==length(perimeter) &&
                        all(choices->length(choices)==2,
                            values(adjacency))) || return false
                    # Walk the single cycle; reject multicavities.
                    polygon=Int32[]
                    start=perimeter[1][1];cur=start;prev=Int32(0)
                    for _ in 1:length(perimeter)+1
                        push!(polygon,cur)
                        choices=adjacency[cur]
                        nxt=choices[1]==prev ? choices[2] : choices[1]
                        prev,cur=cur,nxt
                        cur==start && break
                    end
                    (cur==start && length(polygon)==length(perimeter)) ||
                        return false
                    length(polygon)<=96 || return false
                    inset=Set(polygon)
                    interior=Set{Int32}()
                    for column in cavity
                        for vertex in cells[column]
                            vertex in inset || push!(interior,vertex)
                        end
                    end
                    all(vertex->node_alive[vertex] && !retained[vertex],
                        interior) || return false
                    # Orient the walk with the cavity winding.
                    seed=cells[first(cavity)]
                    w=orient2(parameters[seed[1]],parameters[seed[2]],
                              parameters[seed[3]])
                    w!=0 || return false
                    shoelace=sum(parameters[polygon[index]][1]*
                                 parameters[polygon[mod1(index+1,
                                 length(polygon))]][2]-
                                 parameters[polygon[index]][2]*
                                 parameters[polygon[mod1(index+1,
                                 length(polygon))]][1]
                                 for index in 1:length(polygon))
                    sign(shoelace)==w || reverse!(polygon)
                    newcells=NTuple{3,Int32}[];ear_ok=true
                    remaining=collect(polygon)
                    while length(remaining)>3
                        m=length(remaining);best=0;best_area=-Inf
                        for index in 1:m
                            prev=remaining[mod1(index-1,m)]
                            vtx=remaining[index]
                            nxtv=remaining[mod1(index+1,m)]
                            orient2(parameters[prev],parameters[vtx],
                                    parameters[nxtv])==w || continue
                            area=triangle_area(coordinates[prev],
                                coordinates[vtx],coordinates[nxtv])
                            (isfinite(area) && area>0) || continue
                            any(k->k!=mod1(index-1,m) && k!=index &&
                                    k!=mod1(index+1,m) &&
                                    begin
                                        q=parameters[remaining[k]]
                                        orient2(parameters[prev],
                                                parameters[vtx],q)==w &&
                                        orient2(parameters[vtx],
                                                parameters[nxtv],q)==w &&
                                        orient2(parameters[nxtv],
                                                parameters[prev],q)==w
                                    end,1:m) && continue
                            diagonal=minmax(prev,nxtv)
                            any(column->!(column in cavity),
                                get(incidence,diagonal,Int[])) && continue
                            area>best_area && (best=index;best_area=area)
                        end
                        best==0 && (ear_ok=false;break)
                        index=best
                        prev=remaining[mod1(index-1,m)]
                        vtx=remaining[index]
                        nxtv=remaining[mod1(index+1,m)]
                        push!(newcells,(prev,vtx,nxtv))
                        deleteat!(remaining,index)
                    end
                    if ear_ok
                        prev,vtx,nxtv=remaining
                        area=triangle_area(coordinates[prev],coordinates[vtx],
                                           coordinates[nxtv])
                        ear_ok=orient2(parameters[prev],parameters[vtx],
                                       parameters[nxtv])==w &&
                               isfinite(area) && area>0
                        ear_ok && push!(newcells,(prev,vtx,nxtv))
                    end
                    # When the zigzag perimeter admits no ear, fall back to a
                    # centroid fan: one Steiner vertex on the rail midline
                    # still gives every fan cell the ring's real width.
                    # For a planar surface the affine frame makes the world
                    # centroid and the parameter centroid the same vertex.
                    steiner=nothing
                    if !ear_ok
                        empty!(newcells)
                        nump=length(polygon)
                        cuv=(sum(parameters[v][1] for v in polygon)/nump,
                             sum(parameters[v][2] for v in polygon)/nump)
                        cxyz=(sum(coordinates[v][1] for v in polygon)/nump,
                              sum(coordinates[v][2] for v in polygon)/nump,
                              sum(coordinates[v][3] for v in polygon)/nump)
                        (haskey(positions,cxyz) || cxyz in retired) &&
                            return false
                        fan=NTuple{3,Int32}[];fan_ok=true
                        for index in 1:nump
                            first=polygon[index]
                            last=polygon[mod1(index+1,nump)]
                            orient2(parameters[first],parameters[last],
                                    cuv)==w || (fan_ok=false;break)
                            area=triangle_area(coordinates[first],
                                coordinates[last],cxyz)
                            (isfinite(area) && area>0) || (fan_ok=false;break)
                            push!(fan,(first,last,Int32(-1)))
                        end
                        fan_ok || return false
                        steiner=(cxyz,cuv);newcells=fan
                    end
                    # No surviving live vertex may sit strictly inside a new
                    # cell; the cavity's own interior vertices are dying.
                    near=Set{Int32}(polygon)
                    for column in cavity
                        for vertex in cells[column];push!(near,vertex);end
                    end
                    for column in cavity
                        for cavity_edge in _metric_surface_edges(cells[column])
                            for neighbour in get(incidence,
                                               minmax(cavity_edge...),Int[])
                                alive[neighbour] && !(neighbour in cavity) ||
                                    continue
                                for vertex in cells[neighbour]
                                    push!(near,vertex)
                                end
                            end
                        end
                    end
                    setdiff!(near,interior)
                    vertex_uv=(vertex)->vertex<0 ? steiner[2] : parameters[vertex]
                    for child in newcells
                        c1,c2,c3=child
                        for wtest in near
                            (wtest==c1 || wtest==c2 || wtest==c3) && continue
                            orient2(vertex_uv(c1),vertex_uv(c2),
                                    parameters[wtest])==w &&
                            orient2(vertex_uv(c2),vertex_uv(c3),
                                    parameters[wtest])==w &&
                            orient2(vertex_uv(c3),vertex_uv(c1),
                                    parameters[wtest])==w || continue
                            return false
                        end
                    end
                    # Commit. A fan first births its Steiner vertex, then the
                    # placeholder resolves to the new row.
                    steiner===nothing ||
                        (_metric_surface_birth(node_max,1,caller,"node"))
                    _metric_surface_birth(cell_max,length(newcells),
                                          caller,"cell")
                    if steiner!==nothing
                        cxyz,cuv=steiner
                        row=Int32(length(coordinates)+1)
                        next_node=_metric_surface_birth(node_max,1,caller,
                                                      "node")
                        push!(coordinates,cxyz);push!(parameters,cuv)
                        positions[cxyz]=row
                        push!(numbers,next_node)
                        push!(owners,(2,Int32(tag)))
                        push!(node_parameters,collect(cuv))
                        push!(retained,false);push!(birth_rows,row)
                        push!(node_alive,true);node_max=next_node
                        newcells=NTuple{3,Int32}[
                            (cell[1],cell[2],row) for cell in newcells]
                    end
                    for column in cavity
                        alive[column]=false
                        for old_edge in _metric_surface_edges(cells[column])
                            key=minmax(old_edge...)
                            filter!(!=(column),incidence[key])
                            isempty(incidence[key]) && delete!(incidence,key)
                        end
                    end
                    new_edges=Set{NTuple{2,Int32}}()
                    for child in newcells
                        push!(cells,child);push!(alive,true)
                        cell_max+=UInt64(1);push!(labels,cell_max)
                        for child_edge in _metric_surface_edges(child)
                            key=minmax(child_edge...)
                            push!(get!(()->Int[],incidence,key),length(cells))
                            push!(new_edges,key)
                        end
                    end
                    for vertex in interior
                        node_alive[vertex]=false
                        push!(retired,coordinates[vertex])
                        delete!(positions,coordinates[vertex])
                    end
                    live_count+=length(newcells)-length(cavity)
                    for new_edge in sort!(collect(new_edges))
                        enqueue(new_edge)
                    end
                    return true
                end
                # When the two-cell star admits no valid retriangulation —
                # a razor quadrilateral with a reflex apex — dissolve the
                # whole vertex star instead: remove a non-retained interior
                # vertex and re-triangulate its link polygon. The covered
                # polygon is unchanged, so the area and boundary certificates
                # hold exactly; the new cells only have to pass the same
                # orientation and physical-area certificates. This is the
                # 2-D counterpart of the volumetric kernel's near-coincident
                # cavity dissolution and is what actually releases a locked
                # razor cluster.
                function try_remove_vertex(u)
                    node_alive[u] && !retained[u] || return false
                    starcells=Int32[]
                    for index in eachindex(cells,alive)
                        alive[index] && u in cells[index] && push!(starcells,index)
                    end
                    length(starcells)>=3 || return false
                    link=NTuple{2,Int32}[]
                    for column in starcells
                        cell=cells[column]
                        slot=only(index for index in 1:3 if cell[index]==u)
                        minmax(u,cell[mod1(slot+1,3)]) in boundary &&
                            return false
                        minmax(u,cell[mod1(slot+2,3)]) in boundary &&
                            return false
                        push!(link,(cell[mod1(slot+1,3)],cell[mod1(slot+2,3)]))
                    end
                    length(unique(minmax(e...) for e in link))==length(link) ||
                        return false
                    nxt=Dict{Int32,Int32}(e[1]=>e[2] for e in link)
                    length(nxt)==length(link) || return false
                    polygon=Int32[];start=link[1][1];cur=start
                    for _ in 1:length(link)
                        push!(polygon,cur)
                        cur=get(nxt,cur,Int32(0))
                        cur==0 && return false
                    end
                    cur==start || return false
                    w=orient2(parameters[cells[starcells[1]][1]],
                              parameters[cells[starcells[1]][2]],
                              parameters[cells[starcells[1]][3]])
                    (isfinite(w) && w!=0) || return false
                    # Ear-clip the link polygon in parameter space, enforcing
                    # the star's winding, no vertex strictly inside an ear and
                    # the physical-area certificate on every candidate.
                    remaining=collect(polygon)
                    newcells=NTuple{3,Int32}[]
                    while length(remaining)>3
                        m=length(remaining)
                        best=0;best_area=-Inf
                        for index in 1:m
                            prev=remaining[mod1(index-1,m)]
                            vtx=remaining[index]
                            nxtv=remaining[mod1(index+1,m)]
                            orient2(parameters[prev],parameters[vtx],
                                    parameters[nxtv])==w || continue
                            area=triangle_area(coordinates[prev],coordinates[vtx],
                                               coordinates[nxtv])
                            (isfinite(area) && area>0) || continue
                            any(k->k!=mod1(index-1,m) && k!=index &&
                                    k!=mod1(index+1,m) &&
                                    begin
                                        q=parameters[remaining[k]]
                                        orient2(parameters[prev],parameters[vtx],q)==w &&
                                        orient2(parameters[vtx],parameters[nxtv],q)==w &&
                                        orient2(parameters[nxtv],parameters[prev],q)==w
                                    end,1:m) && continue
                            haskey(incidence,minmax(prev,nxtv)) && continue
                            area>best_area && (best=index;best_area=area)
                        end
                        best==0 && return false
                        index=best
                        prev=remaining[mod1(index-1,m)]
                        vtx=remaining[index]
                        nxtv=remaining[mod1(index+1,m)]
                        push!(newcells,(prev,vtx,nxtv))
                        deleteat!(remaining,index)
                    end
                    prev,vtx,nxtv=remaining
                    orient2(parameters[prev],parameters[vtx],
                            parameters[nxtv])==w || return false
                    area=triangle_area(coordinates[prev],coordinates[vtx],
                                       coordinates[nxtv])
                    (isfinite(area) && area>0) || return false
                    push!(newcells,(prev,vtx,nxtv))
                    _metric_surface_birth(cell_max,length(newcells),caller,"cell")
                    for column in starcells
                        alive[column]=false
                        for old_edge in _metric_surface_edges(cells[column])
                            key=minmax(old_edge...)
                            filter!(!=(column),incidence[key])
                            isempty(incidence[key]) && delete!(incidence,key)
                        end
                    end
                    new_edges=Set{NTuple{2,Int32}}()
                    for child in newcells
                        push!(cells,child);push!(alive,true)
                        cell_max+=UInt64(1);push!(labels,cell_max)
                        for child_edge in _metric_surface_edges(child)
                            key=minmax(child_edge...)
                            push!(get!(()->Int[],incidence,key),length(cells))
                            push!(new_edges,key)
                        end
                    end
                    node_alive[u]=false;push!(retired,coordinates[u])
                    delete!(positions,coordinates[u])
                    live_count-=2
                    for new_edge in sort!(collect(new_edges));enqueue(new_edge);end
                    return true
                end
                star=(edge,Tuple(sort!(NTuple{3,Int32}[
                    Tuple(sort!(collect(s))) for s in sides])))
                retries=get(edge_star_retries,edge,0)
                repaired=false
                if !(star in retried_stars) && retries<32
                    push!(retried_stars,star);edge_star_retries[edge]=retries+1
                    repaired=try_flip() || try_split_at_vertex() ||
                             try_dissolve_cavity()
                    if !repaired && length(sides)==2
                        for u in (sides[1][3],sides[2][3],sides[1][1],
                                  sides[1][2])
                            try_remove_vertex(u) && (repaired=true;break)
                        end
                    end
                end
                repaired && continue
                signature=star[2]
                deferrals=get(edge_deferrals,edge,0)
                if get(deferred_signature,edge,nothing)!=signature && deferrals<8
                    deferred_signature[edge]=signature
                    edge_deferrals[edge]=deferrals+1
                    push!(deferred,edge)
                    continue
                end
                haskey(positions,midpoint) && throw(ArgumentError(
                    "$caller: edge midpoint is not a distinct representable Surface vertex"))
                # Re-run the certificate sequence on the canonical midpoint
                # to report the original failure kind.
                for (first,last,other) in sides
                    winding=orient2(parameters[first],parameters[last],parameters[other])
                    first_child=orient2(parameters[first],uv,parameters[other])
                    second_child=orient2(uv,parameters[last],parameters[other])
                    isfinite(winding) && isfinite(first_child) && isfinite(second_child) &&
                        winding!=0 && sign(first_child)==sign(winding) &&
                        sign(second_child)==sign(winding) || throw(ArgumentError(
                            "$caller: midpoint creates a flat or inverted Surface triangle"))
                    child_areas=(triangle_area(coordinates[first],midpoint,coordinates[other]),
                                 triangle_area(midpoint,coordinates[last],coordinates[other]))
                    all(area->isfinite(area) && area>0,child_areas) || throw(ArgumentError(
                        "$caller: midpoint creates an unrepresentable physical Surface triangle"))
                end
                error("$caller: unreachable split-admissibility state")
            end
            splitpoint,split_uv=found
        end
        next_node=_metric_surface_birth(node_max,1,caller,"node")
        next_cell=_metric_surface_birth(cell_max,2length(hits),caller,"cell")
        row=Int32(length(coordinates)+1)
        children=NTuple{3,Int32}[]
        for (first,last,other) in sides
            push!(children,(first,row,other),(row,last,other))
        end
        push!(coordinates,splitpoint);push!(parameters,split_uv);positions[splitpoint]=row
        push!(numbers,next_node);push!(owners,(2,Int32(tag)));push!(node_parameters,collect(split_uv))
        push!(retained,false);push!(birth_rows,row);push!(node_alive,true);node_max=next_node
        for column in hits
            alive[column]=false
            for old_edge in _metric_surface_edges(cells[column])
                key=minmax(old_edge...);filter!(!=(column),incidence[key])
                isempty(incidence[key]) && delete!(incidence,key)
            end
        end
        new_edges=Set{NTuple{2,Int32}}()
        for child in children
            push!(cells,child);push!(alive,true);cell_max+=UInt64(1);push!(labels,cell_max)
            for child_edge in _metric_surface_edges(child)
                key=minmax(child_edge...);push!(get!(()->Int[],incidence,key),length(cells));push!(new_edges,key)
            end
        end
        cell_max==next_cell || error("$caller: internal cell birth count mismatch")
        live_count+=length(hits)
        for new_edge in sort!(collect(new_edges));enqueue(new_edge);end
    end
    # ── Certified coarsening post-pass ───────────────────────────────────
    # Edge bisection can only shorten edges: it cannot reorient connectivity
    # onto the metric ellipse, so under anisotropic fields the fill packs
    # vertices toward the smallest prescribed direction and overshoots the
    # requested density by orders of magnitude. Two certified sweeps undo
    # that collapse once every edge already satisfies the bound:
    #
    #  * metric diagonal flips trade the longer diagonal of every interior
    #    quad for the shorter one — the same certified 2-2 edge-star
    #    retriangulation the razor repair uses — walking the connectivity
    #    toward the metric-Delaunay triangulation without touching the
    #    boundary or any vertex;
    #  * star dissolution removes a non-retained interior vertex whose link
    #    polygon ear-clips with every new chord at or below `target`, so
    #    removal never recreates a violating edge.
    #
    # Both keep the boundary set, manifold incidence, exact winding and the
    # physical-area certificate; the pair repeats while either makes
    # progress so alignment created by a flip can unblock the next removal.
    vertcells=Dict{Int32,Vector{Int}}()
    for (column,cell) in enumerate(cells)
        alive[column] || continue
        for vertex in cell
            push!(get!(()->Int[],vertcells,vertex),column)
        end
    end
    function starcells_of(vertex)
        star=Int32[]
        for column in get(vertcells,vertex,Int[])
            alive[column] && vertex in cells[column] && push!(star,column)
        end
        return star
    end
    function drop_cell(column)
        alive[column] || return
        alive[column]=false
        for old_edge in _metric_surface_edges(cells[column])
            key=minmax(old_edge...)
            filter!(!=(column),incidence[key])
            isempty(incidence[key]) && delete!(incidence,key)
        end
        for vertex in cells[column]
            hits_v=get(vertcells,vertex,nothing)
            hits_v===nothing || filter!(!=(column),hits_v)
        end
    end
    function add_cell(child)
        push!(cells,child);push!(alive,true)
        cell_max+=UInt64(1);push!(labels,cell_max)
        for child_edge in _metric_surface_edges(child)
            push!(get!(()->Int[],incidence,minmax(child_edge...)),length(cells))
        end
        for vertex in child
            push!(get!(()->Int[],vertcells,vertex),length(cells))
        end
    end
    # Certified 2-2 diagonal flip for the post-pass: the new diagonal must be
    # strictly shorter under the metric (so the sweep cannot cycle), the
    # quad must be strictly convex in parameter space, the diagonal must not
    # already exist, and both children keep the winding and a positive
    # physical area. A convex quad contains no other live vertex, so the
    # children inherit the no-inside-vertex certificate from their parents.
    function metric_flip(edge)
        edge in boundary && return false
        hits=copy(get(incidence,edge,Int[]))
        (length(hits)!=2 ||
         !alive[hits[1]] || !alive[hits[2]]) && return false
        cell1,cell2=cells[hits[1]],cells[hits[2]]
        slot1=only(index for index in 1:3 if
                   minmax(cell1[index],cell1[mod1(index+1,3)])==edge)
        slot2=only(index for index in 1:3 if
                   minmax(cell2[index],cell2[mod1(index+1,3)])==edge)
        f1,l1,o1=cell1[slot1],cell1[mod1(slot1+1,3)],cell1[mod1(slot1+2,3)]
        f2,l2,o2=cell2[slot2],cell2[mod1(slot2+1,3)],cell2[mod1(slot2+2,3)]
        (f2,l2)==(l1,f1) && o1!=o2 || return false
        old=score((f1,l1));new=score((o1,o2))
        (isfinite(new) && new<=target && new<old*(1-1e-9)) || return false
        haskey(incidence,minmax(o1,o2)) && return false
        w1=orient2(parameters[f1],parameters[l1],parameters[o1])
        w2=orient2(parameters[f2],parameters[l2],parameters[o2])
        (isfinite(w1) && isfinite(w2) && w1!=0 && w1==w2) || return false
        # The quad must be strictly convex in parameter space — a reflex
        # corner lets the swapped children cover area outside the two
        # dropped cells and overlap a neighbour.
        orient2(parameters[f1],parameters[o1],parameters[l1])==w1 &&
            orient2(parameters[o1],parameters[l1],parameters[o2])==w1 &&
            orient2(parameters[l1],parameters[o2],parameters[f1])==w1 &&
            orient2(parameters[o2],parameters[f1],parameters[o1])==w1 ||
            return false
        orient2(parameters[f1],parameters[o2],parameters[o1])==w1 &&
            orient2(parameters[o2],parameters[l1],parameters[o1])==w1 ||
            return false
        area1=triangle_area(coordinates[f1],coordinates[o2],coordinates[o1])
        area2=triangle_area(coordinates[o2],coordinates[l1],coordinates[o1])
        (isfinite(area1) && area1>0 && isfinite(area2) && area2>0) ||
            return false
        isapprox(area1+area2,
                 triangle_area(coordinates[f1],coordinates[l1],
                               coordinates[o1])+
                 triangle_area(coordinates[l1],coordinates[f1],
                               coordinates[o2]);rtol=1e-9,atol=0) ||
            return false
        _metric_surface_birth(cell_max,2,caller,"cell")
        drop_cell(hits[1]);drop_cell(hits[2])
        add_cell((f1,o2,o1));add_cell((o2,l1,o1))
        return true
    end
    # Certified vertex removal for the post-pass: ear-clip the interior
    # vertex's link polygon in the star's winding exactly like the razor
    # dissolution, and additionally require every emitted chord to score at
    # or below `target` — coarsening must not hand the loop a violating
    # edge it has already drained.
    function remove_star(vertex)
        node_alive[vertex] && !retained[vertex] || return false
        star=starcells_of(vertex)
        length(star)>=3 || return false
        link=NTuple{2,Int32}[]
        for column in star
            cell=cells[column]
            slot=only(index for index in 1:3 if cell[index]==vertex)
            minmax(vertex,cell[mod1(slot+1,3)]) in boundary && return false
            minmax(vertex,cell[mod1(slot+2,3)]) in boundary && return false
            push!(link,(cell[mod1(slot+1,3)],cell[mod1(slot+2,3)]))
        end
        length(unique(minmax(e...) for e in link))==length(link) ||
            return false
        nxt=Dict{Int32,Int32}(e[1]=>e[2] for e in link)
        length(nxt)==length(link) || return false
        polygon=Int32[];start=link[1][1];cur=start
        for _ in 1:length(link)
            push!(polygon,cur)
            cur=get(nxt,cur,Int32(0));cur==0 && return false
        end
        cur==start || return false
        w=orient2(parameters[cells[star[1]][1]],
                  parameters[cells[star[1]][2]],
                  parameters[cells[star[1]][3]])
        (isfinite(w) && w!=0) || return false
        remaining=collect(polygon);newcells=NTuple{3,Int32}[]
        while length(remaining)>3
            m=length(remaining);best=0;best_area=-Inf
            for index in 1:m
                prev=remaining[mod1(index-1,m)]
                vtx=remaining[index]
                nxtv=remaining[mod1(index+1,m)]
                orient2(parameters[prev],parameters[vtx],
                        parameters[nxtv])==w || continue
                area=triangle_area(coordinates[prev],coordinates[vtx],
                                   coordinates[nxtv])
                (isfinite(area) && area>0) || continue
                any(k->k!=mod1(index-1,m) && k!=index &&
                        k!=mod1(index+1,m) &&
                        begin
                            q=parameters[remaining[k]]
                            orient2(parameters[prev],parameters[vtx],q)==w &&
                            orient2(parameters[vtx],parameters[nxtv],q)==w &&
                            orient2(parameters[nxtv],parameters[prev],q)==w
                        end,1:m) && continue
                haskey(incidence,minmax(prev,nxtv)) && continue
                score((prev,nxtv))<=target || continue
                area>best_area && (best=index;best_area=area)
            end
            best==0 && return false
            index=best
            prev=remaining[mod1(index-1,m)]
            vtx=remaining[index]
            nxtv=remaining[mod1(index+1,m)]
            push!(newcells,(prev,vtx,nxtv))
            deleteat!(remaining,index)
        end
        prev,vtx,nxtv=remaining
        orient2(parameters[prev],parameters[vtx],
                parameters[nxtv])==w || return false
        area=triangle_area(coordinates[prev],coordinates[vtx],
                           coordinates[nxtv])
        (isfinite(area) && area>0) || return false
        score((prev,nxtv))<=target || return false
        push!(newcells,(prev,vtx,nxtv))
        # Overlap certificates. The ear-clip above only inspects polygon
        # corners; a neighbouring vertex could still sit strictly inside or
        # strictly on an emitted ear. Check every vertex of every live cell
        # touching the polygon (the star cells themselves are about to be
        # removed, so only their corners matter as emitters).
        near=Set{Int32}(polygon)
        for column in star
            for corner in cells[column];push!(near,corner);end
        end
        for index in 1:length(polygon)
            u=polygon[index];v=polygon[mod1(index+1,length(polygon))]
            for column in get(incidence,minmax(u,v),Int[])
                column in star && continue
                alive[column] || continue
                for corner in cells[column];push!(near,corner);end
            end
        end
        delete!(near,vertex)
        for child in newcells
            c1,c2,c3=child
            for wtest in near
                (wtest==c1 || wtest==c2 || wtest==c3) && continue
                pw=parameters[wtest]
                o1=orient2(parameters[c1],parameters[c2],pw)
                o2=orient2(parameters[c2],parameters[c3],pw)
                o3=orient2(parameters[c3],parameters[c1],pw)
                (o1==w && o2==w && o3==w) && return false
                for (o,(s,e)) in ((o1,(parameters[c1],parameters[c2])),
                                  (o2,(parameters[c2],parameters[c3])),
                                  (o3,(parameters[c3],parameters[c1])))
                    o==0 || continue
                    (min(s[1],e[1])<=pw[1]<=max(s[1],e[1]) &&
                     min(s[2],e[2])<=pw[2]<=max(s[2],e[2])) && return false
                end
            end
        end
        # The emitted ears must retile exactly the region the star covered:
        # any ear escaping a self-intersecting link polygon, or overlapping
        # a neighbour, unbalances the world-area sum.
        new_area=sum(triangle_area(coordinates[child[1]],
                                   coordinates[child[2]],
                                   coordinates[child[3]])
                     for child in newcells)
        star_area=sum(triangle_area(coordinates[cells[column][1]],
                                    coordinates[cells[column][2]],
                                    coordinates[cells[column][3]])
                      for column in star)
        (isfinite(new_area) && isfinite(star_area) &&
         isapprox(new_area,star_area;rtol=1e-9,
                  atol=1024eps(Float64)*star_area)) || return false
        _metric_surface_birth(cell_max,length(newcells),caller,"cell")
        for column in star;drop_cell(column);end
        for child in newcells;add_cell(child);end
        node_alive[vertex]=false
        push!(retired,coordinates[vertex])
        delete!(positions,coordinates[vertex])
        live_count-=2
        return true
    end
    for _ in 1:64
        moved=false
        for edge in sort!(collect(keys(incidence)));metric_flip(edge) && (moved=true);end
        verts=Int32[vertex for vertex in eachindex(node_alive) if
                    node_alive[vertex] && !retained[vertex]]
        shortest=Dict{Int32,Float64}()
        for (edge,hits) in incidence
            isempty(hits) && continue
            violation=score(edge)
            isfinite(violation) || continue
            for vertex in edge
                node_alive[vertex] || continue
                shortest[vertex]=min(get(shortest,vertex,Inf),violation)
            end
        end
        sort!(verts,by=vertex->get(shortest,vertex,Inf))
        for vertex in verts
            node_alive[vertex] || continue
            remove_star(vertex) && (moved=true)
        end
        moved || break
    end
    get(ENV,"TESS_METRIC_DEBUG","")!="" && println(stderr,
        "metric steiner: attempts=$(steiner_stats[1]) commits=$(steiner_stats[2]) fails=$(sort!(collect(steiner_fail))) tris=$live_count nodes=$(count(node_alive))")
    live=findall(alive)
    # Star dissolutions retire born vertices in place; compact the node rows
    # so the emitted mesh carries only live vertices. Every retained lower row
    # is always live, so the retained prefix keeps its order and coordinates.
    live_nodes=findall(node_alive)
    remap=zeros(Int32,length(coordinates))
    for (index,row) in enumerate(live_nodes);remap[row]=Int32(index);end
    connectivity=Matrix{Int32}(undef,3,length(live))
    for (column,index) in enumerate(live)
        cell=cells[index]
        connectivity[:,column].=(remap[cell[1]],remap[cell[2]],remap[cell[3]])
    end
    output_coordinates=Matrix{Float64}(undef,3,length(live_nodes))
    for (row,point) in enumerate(coordinates[live_nodes]);output_coordinates[:,row].=point;end
    output=Mesh(output_coordinates;tris=connectivity)
    validate(output).ok || throw(ErrorException("$caller: refined Surface is invalid"))
    final_incidence=_metric_surface_incidence(NTuple{3,Int32}[Tuple(cell) for cell in eachcol(output.tris)],caller)
    output_boundary=Set((remap[edge[1]],remap[edge[2]]) for edge in boundary)
    Set(edge for (edge,hits) in final_incidence if length(hits)==1)==output_boundary || throw(ErrorException(
        "$caller: refined Surface changes its retained boundary"))
    output_score(edge)=metric_edge_length(field,
        Tuple(output_coordinates[:,edge[1]]),Tuple(output_coordinates[:,edge[2]]);
        entity=(2,tag))
    if !all(edge->output_score(edge)<=target,keys(final_incidence))
        if get(ENV,"TESS_METRIC_DEBUG","")!=""
            for edge in keys(final_incidence)
                s=output_score(edge)
                s>target && println(stderr,"violating edge $edge s=$s at ",
                    Tuple(output_coordinates[:,edge[1]])," ",Tuple(output_coordinates[:,edge[2]]),
                    " boundary=$(edge in boundary) hits=$(length(final_incidence[edge]))")
            end
        end
        throw(ErrorException("$caller: final edge exceeds max_metric_edge"))
    end
    area(mesh)=sum(triangle_area((Tuple(mesh.coords[:,row]) for row in cell)...) for cell in eachcol(mesh.tris))
    initial_area=area(mesh);final_area=area(output)
    isfinite(initial_area) && initial_area>0 && isfinite(final_area) &&
        isapprox(final_area,initial_area;rtol=256eps(Float64)*max(1,ntris(output)),atol=0) ||
        throw(ErrorException("$caller: refinement changes the Surface area"))
    isequal(output.coords[:,1:nnodes(mesh)],mesh.coords) &&
        numbers[live_nodes][1:nnodes(mesh)]==initial.node_tags ||
        throw(ErrorException("$caller: refinement changes a retained lower vertex"))
    return _ModelSurfaceProduction(output,numbers[live_nodes],owners[live_nodes],
        node_parameters[live_nodes],retained[live_nodes],node_max,
        labels[live],cell_max,Int32[remap[row] for row in birth_rows if node_alive[row]])
end

"""
    mesh_model_surface(model, tag, field::AbstractAnisoField;
                       max_metric_edge=1.0, max_nodes=10_000_000,
                       max_triangles=20_000_000) -> Mesh

Generate an ordinary planar triangular Surface using a spatially varying tensor
field. Interior refinement prefers metric-Delaunay Steiner insertions at
certified circumcenters, falls back to conforming edge bisection that measures
the full world-coordinate metric — equivalently its tangential pullback — at
both endpoints and the midpoint of every edge, and finishes with certified
metric flip and vertex-removal sweeps so the count tracks the prescribed
density. The output satisfies
`metric_edge_length(field,a,b)<=max_metric_edge`.

Actual lower Point/Line products retain their UID, XYZ and boundary incidence.
They must already satisfy the requested bound: coarse boundary edges are rejected
instead of replacing owned Curve cells. Straight planar boundaries, including
holes, are supported. Transfinite/extruded/periodic/recombined Surfaces, embeddings
and Surface smoothing are rejected. New vertices and cells have owned constructor
births; resource and representability failures return no partial mesh or model
mutation. The supplied field, including any explicit wrapper, is used as given.

This explicit tensor overload leaves the existing native scalar generation path
unchanged. It does not claim BAMG ordering, native anisotropic API generation, a
minimum triangle angle, or pointwise bounds between the three metric samples.
"""
function mesh_model_surface(m::GeoModel,tag::Integer,field::AbstractAnisoField;
        max_metric_edge=1.0,max_nodes=PIPELINE_DEFAULT_MAX_NODES,max_triangles=20_000_000,
        _mesh_input=nothing,_mesh1d_options=nothing,_mesh_node_high_water=nothing,
        _mesh_element_high_water=UInt64(0),_return_surface_production=false)
    caller="mesh_model_surface tensor metric"
    target=pipeline_float(max_metric_edge,caller,"max_metric_edge")
    target>0 || throw(ArgumentError("$caller: max_metric_edge must be positive"))
    node_limit=pipeline_limit(max_nodes,caller,"max_nodes",typemax(Int32))
    triangle_limit=pipeline_limit(max_triangles,caller,"max_triangles",typemax(Int32)÷2)
    t=_tag(tag,caller,2);_metric_surface_admit(m,t,caller)
    working=_model_mesh_transaction_copy(m)
    options=_model_algorithm_size_options(working,_mesh1d_options,caller)
    product=_metric_surface_refine(working,t,field,target,node_limit,triangle_limit,
        _mesh_input,options,_mesh_node_high_water,_mesh_element_high_water,caller)
    return _return_surface_production ? product : product.mesh
end

# ── `Mesh.Algorithm=7` (upstream BAMG selector) ──────────────────────────────
#
# The ordinary surface dispatcher resolves algorithm 7 to the conforming
# tensor-bisection producer above — the native metric-driven construction for
# planar Surfaces, not a renamed BAMG. Its metric is the selected tensor-valued
# background field; when the selection is scalar — or no field is selected at
# all — the ordinary size context is wrapped as the isotropic metric I/h²,
# which is BAMG's own scalar fallback. The bisection target is the BAMG unit
# metric length 1.0.

# Scalar sizing context as the isotropic metric I/h²: `_ModelAlgorithmBackgroundSize`
# already evaluates BGM_MeshSize (raw field or the `size_field` argument, the
# curvature and callback terms, the CTX::lc cap, the lcMin/lcMax clamp and
# lcFactor), so wrapping it keeps the fallback on the exact scalar contract.
struct _MetricIsotropicSize{F} <: AbstractAnisoField
    size::F
end

metric_at(field::_MetricIsotropicSize,x,y,z)=isotropic_metric(field.size(x,y,z))

# Resolve the producer's metric field the way `_ModelAlgorithmBackgroundSize`
# resolves its scalar: `options.inner_field` is the raw selected background
# and shadows the `size_field` argument. A tensor selection gets the global
# size policy applied per direction — the cap at CTX::lc, the [lcMin, lcMax]
# clamp and lcFactor act on each directional size exactly as they act on the
# scalar value — while a scalar or absent selection takes the isotropic
# context metric. A `size_field` tensor arriving without `inner_field` (a
# direct Model call) is wrapped by the same policy rather than used raw,
# matching how a scalar `size_field` is treated.
function _metric_surface_algorithm_field(m::GeoModel,tag::Int,options,
        size_field,caller::AbstractString)
    raw=options.inner_field===nothing ? size_field : options.inner_field
    if raw isa AbstractAnisoField
        lc_min=isnan(options.lc_min) ? 0.0 : options.lc_min
        lc_max=isnan(options.lc_max) ? GMSH_MAX_SIZE : options.lc_max
        lo=min(lc_min,lc_max)
        hi=max(lo,min(lc_max,options.ctx_lc))
        return BoundedSize(raw;size_min=lo,size_max=hi,factor=options.lc_factor)
    end
    return _MetricIsotropicSize(
        _ModelAlgorithmBackgroundSize(m,2,tag,options,size_field,caller))
end

# The working model is meshed in place: the producer's lazy boundary grading
# writes `curve_params` exactly like the scalar path, while the public
# tensor endpoint keeps its own transaction copy for its no-mutation
# contract. `forced`/`param_sizes` stay with the CDT refinement path — the
# producer consumes the actual lower Curve products through its source.
function _metric_surface_algorithm_mesh(m::GeoModel,tag::Int,caller::AbstractString;
        size_field=nothing,_mesh_input=nothing,_mesh1d_options=nothing,
        _mesh_node_high_water=nothing,_mesh_element_high_water=UInt64(0))
    _metric_surface_admit(m,tag,caller)
    options=_model_algorithm_size_options(m,_mesh1d_options,caller)
    field=_metric_surface_algorithm_field(m,tag,options,size_field,caller)
    return _metric_surface_refine(m,tag,field,1.0,
        PIPELINE_DEFAULT_MAX_NODES,20_000_000,_mesh_input,options,
        _mesh_node_high_water,_mesh_element_high_water,caller)
end
