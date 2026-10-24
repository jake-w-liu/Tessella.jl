# Native old Pyramid relocation uses v2t_cont's immutable-UID comparator.
# Pointer order only serves membership queries in this operation.
module ModelNativeHybridRelocate
const M=parentmodule(@__MODULE__)
point(coords,row)=ntuple(axis->coords[axis,row],3)
sub(a,b)=ntuple(axis->a[axis]-b[axis],3)
dot(a,b)=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]
cross(a,b)=(a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
function unit(a)
    length=sqrt(dot(a,a));length==0.0 && return a
    factor=1.0/length
    return ntuple(axis->a[axis]*factor,3)
end

function condition(a,b,c)
    x1,y1,z1=a;x2,y2,z2=b;x3,y3,z3=c
    determinant=x1*(y2*z3-z2*y3)-x2*(y1*z3-z1*y3)+x3*(y1*z2-z1*y2)
    determinant==0.0 && return 0.0
    i11=y2*z3-y3*z2;i12=x3*z2-x2*z3;i13=x2*y3-x3*y2
    i21=y3*z1-y1*z3;i22=x1*z3-x3*z1;i23=x3*y1-x1*y3
    i31=y1*z2-y2*z1;i32=x2*z1-x1*z2;i33=x1*y2-x2*y1
    normj=x1*x1+x2*x2+x3*x3+y1*y1+y2*y2+y3*y3+z1*z1+z2*z2+z3*z3
    norminv=i11*i11+i12*i12+i13*i13+i21*i21+i22*i22+i23*i23+i31*i31+i32*i32+i33*i33
    return 3.0*determinant/sqrt(normj*norminv)
end

function pyramid_condition(points)
    signs=((-1.0,-1.0),(1.0,-1.0),(1.0,1.0),(-1.0,1.0))
    minimum_value=0.0
    for (sample,(u,v)) in enumerate(((-1.0,-1.0),(1.0,-1.0),(-1.0,1.0),(1.0,1.0)))
        derivatives=ntuple(3) do axis
            ntuple(3) do world
                result=0.0
                for node in 1:5
                    gradient=if node==5
                        axis==3 ? sqrt(2.0) : 0.0
                    else
                        sx,sy=signs[node]
                        axis==1 ? (sx*(1.0+sy*v)*0.25)*2.0 : axis==2 ?
                            (sy*(1.0+sx*u)*0.25)*2.0 :
                            (-0.25+sx*sy*u*v*0.25)*sqrt(2.0)
                    end
                    result+=gradient*points[node][world]
                end
                result
            end
        end
        value=condition(derivatives...)
        sample==1 && (minimum_value=value)
        value<minimum_value && (minimum_value=value)
    end
    return minimum_value
end

function objective(coords,node,incident,target,xi;only_tets=false)
    original=point(coords,node)
    trial=ntuple(axis->(1.0-xi)*original[axis]+xi*target[axis],3)
    minimum_quality=1.0
    for cell in incident
        points=map(row->row==node ? trial : point(coords,row),cell)
        value=if length(cell)==4
            quality,_=M._gopt_gamma(points...)
            M.orient3(points...)>0 ? -quality : quality
        elseif only_tets
            continue
        else
            abs(pyramid_condition(points))*0.2
        end
        minimum_quality=minimum_quality<value ? minimum_quality : value
    end
    return minimum_quality
end

function center_sum(coords,incident;exclude_pyramids=false)
    x=0.0;y=0.0;z=0.0;count=0
    for cell in incident
        exclude_pyramids && length(cell)==5 && continue
        sx=0.0;sy=0.0;sz=0.0
        for row in cell
            sx+=coords[1,row];sy+=coords[2,row];sz+=coords[3,row]
        end
        x+=sx;y+=sy;z+=sz;count+=length(cell)
    end
    count>0 || throw(ArgumentError("hybrid relocation has no current incident Tet constructor"))
    return (x,y,z),count
end

function normal_step!(coords,node,incident,relax)
    sum,count=center_sum(coords,incident;exclude_pyramids=true)
    target=ntuple(axis->sum[axis]/count,3)
    pyramid=last(incident)
    # Native's normal pass visits every movable vertex in the recorded
    # neighborhood; a vertex without a Pyramid is left for the golden pass.
    length(pyramid)==5 || return false
    pyramid[5]==node || throw(ArgumentError(
        "hybrid relocation loses its owned apex Pyramid incidence"))
    rows=(pyramid[1],pyramid[4],pyramid[3],pyramid[2])
    points=map(row->point(coords,row),rows)
    area=M._gopt_area(points[1],points[2],points[3])+
        M._gopt_area(points[1],points[3],points[4])
    unit_normal=unit(unit(cross(sub(points[2],points[1]),sub(points[3],points[1]))))
    center=ntuple(3) do axis
        value=0.0
        for point in points;value+=point[axis];end
        value/4.0
    end
    normal=dot(unit_normal,sub(target,center))<0.0 ?
        (-unit_normal[1],-unit_normal[2],-unit_normal[3]) : unit_normal
    height=0.5*sqrt(abs(area))
    candidate=ntuple(axis->center[axis]+relax*height*normal[axis],3)
    objective(coords,node,incident,candidate,1.0;only_tets=true)>1e-6 || return false
    coords[:,node].=candidate
    return true
end

function golden(coords,node,incident,target,tolerance)
    lambda=0.5*(sqrt(5.0)-1.0);mu=0.5*(3.0-sqrt(5.0))
    a=0.0;b=2.0
    x1=b-lambda*(b-a);x2=a+lambda*(b-a)
    f1=objective(coords,node,incident,target,x1);f2=objective(coords,node,incident,target,x2)
    while !(abs(b-a)<tolerance)
        if f1<f2
            a=x1
            abs(b-a)<tolerance && break
            x1=x2;f1=f2;x2=b-mu*(b-a)
            f2=objective(coords,node,incident,target,x2)
        else
            b=x2
            abs(b-a)<tolerance && break
            x2=x1;f2=f1;x1=a+mu*(b-a)
            f1=objective(coords,node,incident,target,x1)
        end
    end
    return a
end

function golden_step!(coords,node,incident,tolerance)
    sum,count=center_sum(coords,incident)
    target=ntuple(axis->sum[axis]/count,3)
    before=objective(coords,node,incident,target,0.0)
    after=objective(coords,node,incident,target,1.0)
    after>before && (coords[:,node].=target)
    xi=golden(coords,node,incident,target,tolerance)
    original=point(coords,node)
    coords[:,node].=ntuple(axis->(1.0-xi)*original[axis]+xi*sum[axis]/count,3)
    return nothing
end

function relocate(plan,mesh,caller;birth_order=eachindex(plan.births),iterations=3,tolerance=1e-2,
                  node_uids=plan.node_uids)
    length(node_uids)==M.nnodes(mesh) && length(node_uids)>=length(plan.node_uids) &&
        node_uids[1:length(plan.node_uids)]==plan.node_uids && allunique(node_uids) &&
        all(!iszero,node_uids) &&
        isequal(mesh.coords[:,1:length(plan.source_uids)],plan.source_coordinates) &&
        size(mesh.tets,1)==4 && all(row->0<row<=M.nnodes(mesh),mesh.tets) || throw(ArgumentError(
        "$caller: hybrid relocation lacks its exact actual Source and born apex rows"))
    length(birth_order)==length(plan.births) && allunique(birth_order) &&
        all(index->1<=index<=length(plan.births),birth_order) || throw(ArgumentError(
        "$caller: hybrid relocation omits or aliases an owned Pyramid constructor"))
    iterations==3 && tolerance==1e-2 || throw(ArgumentError(
        "$caller: hybrid relocation requires the native generation iteration/tolerance contract"))
    apexes=Set(birth.apex_row for birth in plan.births)
    neighborhood=Set{Int32}()
    for cell in eachcol(mesh.tets)
        any(row->row in apexes,cell) || continue
        union!(neighborhood,cell)
    end
    movable=Set{Int32}()
    for cell in eachcol(mesh.tets)
        any(row->row in neighborhood,cell) || continue
        for row in cell
            # Retained lower nodes keep their actual ownership. Extra rows
            # are the insertion driver's explicitly bound Region births.
            row in apexes || row>length(plan.node_uids) || continue
            push!(movable,row)
        end
    end
    # The relocation objective inspects only `incident` cells, so the native
    # v2t star must be complete: a movable row that also sits in a cell with
    # no neighborhood vertex would otherwise fold that cell unseen.
    adjacency=Dict{Int32,Vector{Tuple}}()
    for cell in eachcol(mesh.tets)
        key=Tuple(cell)
        for row in cell
            row in movable && push!(get!(()->Tuple[],adjacency,row),key)
        end
    end
    for ordinal in birth_order
        birth=plan.births[ordinal]
        haskey(adjacency,birth.apex_row) || throw(ArgumentError(
            "$caller: hybrid apex has no actual current Tet carrier"))
        push!(adjacency[birth.apex_row],(birth.base_rows...,birth.apex_row))
    end
    coords=copy(mesh.coords)
    # meshGFaceOptimize.h v2t_cont explicitly uses MVertexPtrLessThan::getNum.
    order=sort!(collect(keys(adjacency));by=row->node_uids[row])
    for step in 0:9,row in order
        normal_step!(coords,row,adjacency[row],Float64(step)/10.0+1e-6)
    end
    for _ in 1:(iterations+2),row in order
        golden_step!(coords,row,adjacency[row],tolerance)
    end
    all(isfinite,coords) && isequal(coords[:,1:length(plan.source_uids)],plan.source_coordinates) ||
        throw(ArgumentError("$caller: hybrid relocation changes retained Source coordinates"))
    return M.Mesh(coords;tets=copy(mesh.tets),tet_tag=copy(mesh.tet_tag))
end
end
