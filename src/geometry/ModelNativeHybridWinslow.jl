# Regularized Winslow relocation over actual closed Pyramid-apex stars.
# Patch boundaries stay fixed; typed cells and their event UIDs are unchanged.
module ModelNativeHybridWinslow
const M=parentmodule(@__MODULE__)
const R=M.ModelNativeHybridRelocate
const TET_FACES=((1,3,2),(1,2,4),(1,4,3),(4,2,3))
const PYRAMID_FACES=((1,2,5),(4,1,5),(2,3,5),(3,4,5),(1,4,3,2))
const IDEAL_FACES=((2,4,3),(1,3,4),(4,2,1),(1,2,3))

struct PatchTet
    points::NTuple{4,NTuple{3,Float64}}
    normals::NTuple{4,NTuple{3,Float64}}
    variable::Int
end

volume(points)=R.dot(R.sub(points[1],points[4]),
    R.cross(R.sub(points[2],points[4]),R.sub(points[3],points[4])))/6.0

function energy_gradient(point,patch,epsilon)
    energy=0.0;gradient=(0.0,0.0,0.0);minimum_det=Inf
    for cell in patch
        columns=ntuple(3) do world
            ntuple(3) do axis
                value=0.0
                for slot in 1:4
                    coordinate=slot==cell.variable ? point[world] : cell.points[slot][world]
                    value+=coordinate*cell.normals[slot][axis]
                end
                value
            end
        end
        determinant=R.dot(columns[1],R.cross(columns[2],columns[3]))
        minimum_det=min(minimum_det,determinant)
        root=sqrt(determinant*determinant+epsilon*epsilon)
        # The rationalized negative branch is the same chi, with cancellation
        # avoided for trial points on the inverted side of a patch cell.
        chi=determinant>=0.0 ? 0.5*(determinant+root) :
            epsilon*epsilon/(2.0*(root-determinant))
        chi>0.0 && isfinite(chi) || return (Inf,(Inf,Inf,Inf),minimum_det)
        chip=chi/root
        power=chi^(2.0/3.0)
        shape=(R.dot(columns[1],columns[1])+R.dot(columns[2],columns[2])+
               R.dot(columns[3],columns[3]))/power
        size=(determinant*determinant+1.0)/chi
        energy+=shape+size
        derivative=ntuple(3) do world
            cofactor=R.cross(columns[mod1(world+1,3)],columns[mod1(world+2,3)])
            coefficient=((2.0/3.0)*shape*chip-2.0*determinant+size*chip)/chi
            column=ntuple(axis->2.0/power*columns[world][axis]-coefficient*cofactor[axis],3)
            R.dot(column,cell.normals[cell.variable])
        end
        gradient=ntuple(axis->gradient[axis]+derivative[axis],3)
    end
    return energy,gradient,minimum_det
end

function zoom(point,patch,epsilon,direction,initial,low,high,low_energy)
    value,gradient,_=initial;slope=R.dot(gradient,direction)
    for _ in 1:60
        step=0.5*(low+high)
        trial=ntuple(axis->point[axis]+step*direction[axis],3)
        result=energy_gradient(trial,patch,epsilon);candidate,derivative,_=result
        if !isfinite(candidate) || candidate>value+1.e-4*step*slope || candidate>=low_energy
            high=step
        else
            trial_slope=R.dot(derivative,direction)
            abs(trial_slope)<=-0.9*slope && return trial,result
            trial_slope*(high-low)>=0.0 && (high=low)
            low=step;low_energy=candidate
        end
        abs(high-low)*sqrt(R.dot(direction,direction))<=1.e-12 && break
    end
    return nothing
end

function line_search(point,patch,epsilon,direction,initial,first_step)
    value,gradient,_=initial;slope=R.dot(gradient,direction)
    slope<0.0 || return nothing
    previous=0.0;previous_energy=value;step=first_step
    for iteration in 1:40
        trial=ntuple(axis->point[axis]+step*direction[axis],3)
        result=energy_gradient(trial,patch,epsilon);candidate,derivative,_=result
        if !isfinite(candidate) || candidate>value+1.e-4*step*slope ||
                (iteration>1 && candidate>=previous_energy)
            return zoom(point,patch,epsilon,direction,initial,previous,step,previous_energy)
        end
        trial_slope=R.dot(derivative,direction)
        abs(trial_slope)<=-0.9*slope && return trial,result
        trial_slope>=0.0 && return zoom(point,patch,epsilon,direction,initial,
            step,previous,candidate)
        previous=step;previous_energy=candidate;step*=2.0
    end
    return nothing
end

function minimize(point,patch,epsilon)
    current=energy_gradient(point,patch,epsilon)
    all(isfinite,(current[1],current[2]...)) || throw(ArgumentError(
        "hybrid Winslow initial energy/gradient is not finite"))
    history=Tuple{NTuple{3,Float64},NTuple{3,Float64},Float64}[]
    for iteration in 1:130
        value,gradient,_=current
        magnitude=sqrt(R.dot(gradient,gradient))
        magnitude<=1.e-4 && break
        q=gradient;alphas=Float64[]
        for (step,change,inverse) in Iterators.reverse(history)
            alpha=inverse*R.dot(step,q);push!(alphas,alpha)
            q=ntuple(axis->q[axis]-alpha*change[axis],3)
        end
        scale=isempty(history) ? 1.0 : R.dot(history[end][1],history[end][2])/
            R.dot(history[end][2],history[end][2])
        direction=ntuple(axis->scale*q[axis],3)
        for ((step,change,inverse),alpha) in zip(history,Iterators.reverse(alphas))
            beta=inverse*R.dot(change,direction)
            direction=ntuple(axis->direction[axis]+step[axis]*(alpha-beta),3)
        end
        direction=ntuple(axis->-direction[axis],3)
        if !all(isfinite,direction) || R.dot(direction,gradient)>=0.0
            empty!(history);direction=ntuple(axis->-gradient[axis],3)
        end
        initial_step=iteration==1 ? 1.0/magnitude : 1.0
        product=line_search(point,patch,epsilon,direction,current,initial_step)
        product===nothing && break
        candidate,next=product
        step=R.sub(candidate,point);change=R.sub(next[2],gradient)
        curvature=R.dot(step,change)
        if isfinite(curvature) && curvature>0.0
            push!(history,(step,change,1.0/curvature))
            length(history)>15 && popfirst!(history)
        end
        point=candidate;current=next
        abs(value-next[1])<=1.e-12*max(1.0,abs(value),abs(next[1])) && break
        sqrt(R.dot(step,step))<=1.e-12 && break
    end
    return point,current
end

function pyramid_rows(birth,mesh,caller)
    a,b,c,d=birth.base_rows;apex=birth.apex_row
    side=(a,b,apex);key=Tuple(sort!(collect(side)))
    faces=NTuple{3,Int32}[]
    for cell in eachcol(mesh.tets),face in TET_FACES
        candidate=ntuple(slot->cell[face[slot]],3)
        Tuple(sort!(collect(candidate)))==key && push!(faces,candidate)
    end
    length(faces)==1 || throw(ArgumentError(
        "$caller: Winslow Pyramid side lacks its unique actual Tet carrier"))
    same=any(side==ntuple(slot->only(faces)[mod1(slot+rotation,3)],3) for rotation in 0:2)
    # Native's topological pass swaps corners 1/3 when that first side and
    # its actual carrier have equal winding. No coordinate lookup selects it.
    return same ? (a,d,c,b,apex) : (a,b,c,d,apex)
end

function patch(plan,mesh,birth,node_uids,caller)
    apex=birth.apex_row;pyramid=pyramid_rows(birth,mesh,caller)
    incident=[Tuple(cell) for cell in eachcol(mesh.tets) if apex in cell]
    isempty(incident) && throw(ArgumentError("$caller: Winslow apex has no actual Tet star"))
    counts=Dict{Tuple,Int}();oriented=Dict{Tuple,Tuple}();vertices=Set{Int32}()
    for cell in incident
        union!(vertices,cell)
        for face in TET_FACES
            key=Tuple(sort!(Int32[cell[slot] for slot in face]))
            counts[key]=get(counts,key,0)+1
            oriented[key]=Tuple(cell[slot] for slot in face)
        end
    end
    union!(vertices,pyramid)
    for face in PYRAMID_FACES
        key=Tuple(sort!(Int32[pyramid[slot] for slot in face]))
        counts[key]=get(counts,key,0)+1
        oriented[key]=Tuple(pyramid[slot] for slot in face)
    end
    all(count->count in (1,2),values(counts)) || throw(ArgumentError(
        "$caller: Winslow apex star is non-manifold"))
    boundary=Set{Int32}(row for (face,count) in counts if count==1 for row in face)
    setdiff(vertices,boundary)==Set([apex]) || throw(ArgumentError(
        "$caller: Winslow patch does not bind exactly one interior apex"))
    panels=[oriented[face] for (face,count) in counts if count==1]
    links=Dict{Tuple{Int32,Int32},Vector{Int}}()
    for (index,panel) in enumerate(panels),slot in eachindex(panel)
        edge=minmax(panel[slot],panel[mod1(slot+1,length(panel))])
        push!(get!(()->Int[],links,edge),index)
    end
    all(faces->length(faces)==2,values(links)) &&
        length(boundary)-length(links)+length(panels)==2 || throw(ArgumentError(
        "$caller: Winslow patch boundary is not a closed spherical complex"))
    reached=Set([1]);queue=[1];cursor=1
    while cursor<=length(queue)
        index=queue[cursor];cursor+=1
        for faces in values(links)
            index in faces || continue
            for neighbor in faces
                neighbor in reached || (push!(reached,neighbor);push!(queue,neighbor))
            end
        end
    end
    length(reached)==length(panels) || throw(ArgumentError(
        "$caller: Winslow patch has disconnected boundary components"))
    boundary_rows=sort!(collect(boundary);by=row->node_uids[row])
    center=ntuple(3) do axis
        value=0.0
        for row in boundary_rows;value+=mesh.coords[axis,row];end
        value/length(boundary_rows)
    end
    rows=NTuple{4,Int32}[(cell[1],cell[4],cell[3],cell[2]) for cell in incident]
    a,b,c,d,e=pyramid
    append!(rows,[(a,e,c,b),(c,e,a,d),(b,e,d,c),(d,e,b,a)])
    points=[ntuple(slot->R.point(mesh.coords,cell[slot]),4) for cell in rows]
    ideals=copy(points)
    for index in length(incident)+1:length(rows)
        ideals[index]=ntuple(slot->rows[index][slot]==apex ? center : points[index][slot],4)
    end
    average=sum(volume,points)/length(points)
    isfinite(average) && average>0.0 || throw(ArgumentError(
        "$caller: Winslow patch lacks positive mean current Tet volume"))
    result=PatchTet[]
    for (cell,coordinates,ideal) in zip(rows,points,ideals)
        volume(ideal)>0.0 || throw(ArgumentError(
            "$caller: Winslow patch has a degenerate or inverted actual ideal shape"))
        normals=ntuple(4) do slot
            face=IDEAL_FACES[slot]
            u=R.sub(ideal[face[2]],ideal[face[1]])
            v=R.sub(ideal[face[3]],ideal[face[1]])
            normal=R.cross(v,u)
            ntuple(axis->0.5*normal[axis]*(1.0/(3.0*average)),3)
        end
        variable=only(findall(==(apex),collect(cell)))
        push!(result,PatchTet(coordinates,normals,variable))
    end
    return result
end

function relocate(plan,mesh,caller;birth_order=eachindex(plan.births),node_uids=plan.node_uids)
    length(node_uids)==M.nnodes(mesh) && allunique(node_uids) && all(!iszero,node_uids) &&
        node_uids[1:length(plan.node_uids)]==plan.node_uids &&
        isequal(mesh.coords[:,1:length(plan.source_uids)],plan.source_coordinates) || throw(ArgumentError(
        "$caller: Winslow relocation lost its actual Source and apex birth UID rows"))
    length(birth_order)==length(plan.births) && allunique(birth_order) &&
        all(ordinal->1<=ordinal<=length(plan.births),birth_order) || throw(ArgumentError(
        "$caller: Winslow relocation omits or aliases a Pyramid birth"))
    coordinates=copy(mesh.coords)
    current=M.Mesh(coordinates;tets=copy(mesh.tets),tet_tag=copy(mesh.tet_tag))
    coordinates=current.coords
    # Visit the actual owned birth UIDs deterministically. This does not claim
    # equivalence with a platform allocator's SDK pointer visitation order.
    order=sort!(collect(birth_order);by=ordinal->plan.births[ordinal].apex_uid)
    for ordinal in order
        birth=plan.births[ordinal]
        cells=patch(plan,current,birth,node_uids,caller)
        point=R.point(coordinates,birth.apex_row)
        initial=energy_gradient(point,cells,1.e-6)
        epsilon=sqrt(1.e-12+0.04*min(initial[3],0.0)^2)
        previous=energy_gradient(point,cells,epsilon)[1]
        for _ in 1:2
            point,result=minimize(point,cells,epsilon)
            all(isfinite,(point...,result[1],result[2]...)) || throw(ArgumentError(
                "$caller: Winslow optimizer produces nonfinite coordinates or energy"))
            relative=abs(result[1]-previous)/max(result[1],eps(Float64))
            result[3]>0.0 && relative<1.e-5 && break
            previous=result[1];epsilon=sqrt(1.e-12+0.04*min(result[3],0.0)^2)
        end
        final=energy_gradient(point,cells,epsilon)
        final[3]>0.0 && final[1]<=energy_gradient(R.point(coordinates,birth.apex_row),cells,epsilon)[1] ||
            throw(ArgumentError("$caller: Winslow apex does not have positive, improved patch cells"))
        coordinates[:,birth.apex_row].=point
    end
    all(isfinite,coordinates) && isequal(coordinates[:,1:length(plan.source_uids)],
        plan.source_coordinates) || throw(ArgumentError("$caller: Winslow changes fixed Source coordinates"))
    return M.Mesh(coordinates;tets=copy(mesh.tets),tet_tag=copy(mesh.tet_tag))
end
end
