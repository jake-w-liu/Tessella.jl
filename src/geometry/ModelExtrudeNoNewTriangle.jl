# One isolated source triangle has triangular caps and three lateral faces.
# Its caps have no diagonal state: the retained prism or an existing-corner
# three-tet template can be selected independently on each logical interval.
# This route never introduces a centroid or pads a triangle into a quadrangle.

function _extrude_nonew_prism_pick(preferred,recombined::Bool,caller)
    templates=_extrude_nonew_prism_templates()
    best=0
    best_cost=typemax(Int)
    for index in eachindex(templates)
        faces=templates[index].faces
        if recombined
            all(iszero,faces) || continue
        else
            all(!iszero,faces) || continue
        end
        cost=0
        for side in 1:3
            cost+=faces[side]!=preferred[side]
        end
        if cost<best_cost
            best=index
            best_cost=cost
        end
    end
    best!=0 || throw(ArgumentError(
        "$caller: QuadTriNoNewVerts has no conforming isolated source-triangle template"))
    return templates[best]
end

function _extrude_nonew_prism_diagonals!(edges,v,states)
    for side in 1:3
        state=states[side]
        state==0 && continue
        face=_EXTRUDE_NONEW_PRISM_FACES[side+2]
        a,b=state==1 ? (face[1],face[3]) : (face[2],face[4])
        push!(edges,_extrude_ekey(v[a],v[b]))
    end
    return nothing
end

function _extrude_nonew_triangle_finish(m,t,params,source,source_mesh,
        cell::NTuple{3,Int32},preferred,cols,levels,refs,top,laterals,caller)
    intervals=length(refs)
    template=_extrude_nonew_prism_pick(preferred,params.recomb_laterals,caller)
    faces=zeros(UInt8,5,intervals)
    edges=Set{NTuple{2,NTuple{3,Float64}}}()
    params.recomb_laterals || sizehint!(edges,3intervals)
    for layer in 1:intervals
        v=_extrude_nonew_corners(cols,cell,layer)
        for side in 1:3
            faces[side,layer]=template.faces[side]
        end
        _extrude_nonew_prism_diagonals!(edges,v,template.faces)
    end
    sweep=_ExtrudeVolumeSweep(t,true,cols,
        NTuple{6,NTuple{3,Float64}}[],NTuple{4,NTuple{3,Float64}}[],
        NTuple{8,NTuple{3,Float64}}[],NTuple{6,NTuple{3,Float64}}[],
        NTuple{5,NTuple{3,Float64}}[],NTuple{3,Float64}[],Bool[false])
    # This category emits exactly one prism or three tetrahedra per interval.
    sizehint!(params.recomb_laterals ? sweep.prisms : sweep.tets,
        params.recomb_laterals ? intervals : 3intervals)
    for layer in 1:intervals
        v=_extrude_nonew_corners(cols,cell,layer)
        _extrude_nonew_certify(v,edges,template,caller)
        _extrude_nonew_emit_template!(sweep,v,template)
    end
    catalog=_ExtrudeNoNewCatalog(cell,0x07,levels,refs)
    return (sweep=sweep,edges=edges,source_mesh=source_mesh,source_tag=source,
        top_tag=top,lateral_tags=laterals,catalog=catalog,faces=faces,
        problem_layers=falses(intervals))
end
