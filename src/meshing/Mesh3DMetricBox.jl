# A constant tensor controls the initial structured volume grid directly.
# Existing refinement and surface pipelines already consume field directions.
using ..SizeField: Metric3, _metric_endpoint_displacement,
                   _metric_endpoint_displacement_exact
using ..Mesh1D: _convex_coordinate
using ..PipelineSupport: pipeline_limit

function _metric_box_grid_counts(counts,limit_nodes,limit_tets,caller)
    node_counts=ntuple(axis->_checked_add3(caller,"axis node",counts[axis],1),3)
    nodes=_checked_mul3(caller,"node",node_counts...)
    tets=_checked_mul3(caller,"tetrahedron",6,counts...)
    nodes<=limit_nodes || throw(ArgumentError("$caller: $nodes nodes exceed max_nodes=$limit_nodes"))
    tets<=limit_tets || throw(ArgumentError("$caller: $tets tetrahedra exceed max_tets=$limit_tets"))
    return nodes,tets
end

function _metric_region_counts(interval_counts,caller)
    ntuple(3) do axis
        count=0
        for intervals in interval_counts[axis]
            count=_checked_add3(caller,"axis $axis interval",count,intervals)
        end
        count
    end
end

function _metric_region_axis(planes,counts,total,caller,axis)
    values=Vector{Float64}(undef,total+1)
    plane_indices=Vector{Int}(undef,length(planes))
    values[1]=planes[1];plane_indices[1]=0
    cursor=1;maximum_gap=0.0
    for segment in eachindex(counts)
        first=planes[segment];last=planes[segment+1];count=counts[segment]
        for index in 1:count
            value=_convex_coordinate(first,last,index/count)
            gap=value-values[cursor]
            isfinite(gap) && gap>0 || throw(ArgumentError(
                "$caller: axis $axis intervals are below Float64 coordinate resolution"))
            maximum_gap=max(maximum_gap,nextfloat(gap))
            cursor+=1;values[cursor]=value
        end
        plane_indices[segment+1]=cursor-1
    end
    isfinite(maximum_gap) || throw(ArgumentError(
        "$caller: axis $axis interval magnitude is not representable"))
    return values,plane_indices,maximum_gap
end

@inline _metric_region_parameter(planes,indices,value)=
    Float64(indices[searchsortedfirst(planes,value)])

function _metric_box_axis(first,last,count,caller,axis)
    values=[_convex_coordinate(first,last,index/count) for index in 0:count]
    maximum_gap=0.0
    for index in 2:length(values)
        gap=values[index]-values[index-1]
        isfinite(gap) && gap>0 || throw(ArgumentError(
            "$caller: axis $axis intervals are below Float64 coordinate resolution"))
        # Round upward so this magnitude also bounds the exact dyadic
        # difference of the actual two stored Float64 coordinates.
        maximum_gap=max(maximum_gap,nextfloat(gap))
    end
    isfinite(maximum_gap) || throw(ArgumentError(
        "$caller: axis $axis interval magnitude is not representable"))
    return values,maximum_gap
end

# Every Kuhn edge is one of seven nonempty subsets of the three positive
# axis increments. Absolute cross coefficients bound any rounded interval
# from the actual axis arrays, including metrics with negative cross terms.
function _metric_box_edge_bound(metric,maximum_gaps,caller)
    origin=(0.,0.,0.);bound=0.0
    for mask in 1:7
        displacement=ntuple(axis->(mask & (1<<(axis-1)))==0 ? 0. : maximum_gaps[axis],3)
        bound=max(bound,_metric_endpoint_displacement_exact(metric,origin,displacement,caller))
    end
    return bound
end

"""
    mesh_box(x0, x1, y0, y1, z0, z1, metric::Metric3;
             max_metric_edge=1.0, max_nodes=10_000_000,
             max_tets=60_000_000) -> Mesh

Construct a tetrahedral box mesh for a **constant** symmetric positive-definite
metric. Directional axis lengths choose separate interval counts; offdiagonal
entries participate in a conservative certificate for every resulting edge.
The metric bound is `sqrt((b-a)' * M * (b-a)) ≤ max_metric_edge` for the actual
stored endpoint coordinates. This positional tensor overload is distinct from
the six-argument physical `hmax` contract.

An axis-aligned tensor grid is filled with the existing six-tetrahedron Kuhn
subdivision. With diagonal metrics, the relative axis counts follow the three
directional sizes. Offdiagonal metrics retain their edge norm but use a
conservative axis grid rather than a rotation to principal directions.
The boundary stays on the exact requested box planes, and the output is
validated for positive, nondegenerate tetrahedra and manifold topology.

This is a constant-tensor construction, not a sampler for spatially varying
fields or a claim of native anisotropic mesher parity. Use `refine_to_size`
or the surface pipelines for existing direction-aware field refinement.
Resource and Float64 coordinate-resolution failures throw before a result
can escape. Input metrics are immutable and are not modified.
"""
function mesh_box(x0::Real,x1::Real,y0::Real,y1::Real,z0::Real,z1::Real,metric::Metric3;
                  max_metric_edge::Real=1.0,
                  max_nodes=PIPELINE_DEFAULT_MAX_NODES,max_tets=PIPELINE_DEFAULT_MAX_TETS)
    caller="mesh_box constant metric"
    first=(_finite3(x0,caller,"x0"),_finite3(y0,caller,"y0"),_finite3(z0,caller,"z0"))
    last=(_finite3(x1,caller,"x1"),_finite3(y1,caller,"y1"),_finite3(z1,caller,"z1"))
    spans=ntuple(axis->last[axis]-first[axis],3)
    all(span->isfinite(span) && span>0,spans) || throw(ArgumentError(
        "$caller: box must have a finite positive extent in every axis"))
    target=_finite3(max_metric_edge,caller,"max_metric_edge")
    target>0 || throw(ArgumentError("$caller: max_metric_edge must be positive"))
    node_limit=pipeline_limit(max_nodes,caller,"max_nodes",typemax(Int32))
    tet_limit=pipeline_limit(max_tets,caller,"max_tets",typemax(Int32))
    origin=(0.,0.,0.)
    counts=ntuple(3) do axis
        displacement=ntuple(other->other==axis ? spans[axis] : 0.,3)
        directional=_metric_endpoint_displacement(metric,origin,displacement,caller)
        _ceil_count3((directional/target)*sqrt(3.),caller,"axis $axis interval")
    end
    # In 3D, taking absolute offdiagonals preserves SPD: the leading principal
    # minors are unchanged and its determinant cannot decrease. Nonnegative
    # entries make this upper bound monotone in the interval magnitudes.
    envelope=Metric3(metric.m11,metric.m22,metric.m33,
        abs(metric.m12),abs(metric.m13),abs(metric.m23))
    for attempt in 1:32
        _metric_box_grid_counts(counts,node_limit,tet_limit,caller)
        axes_and_gaps=ntuple(axis->_metric_box_axis(first[axis],last[axis],counts[axis],caller,axis),3)
        grid_axes=ntuple(axis->axes_and_gaps[axis][1],3)
        gaps=ntuple(axis->axes_and_gaps[axis][2],3)
        bound=_metric_box_edge_bound(envelope,gaps,caller)
        if bound<=target
            # Reuse the original box fill on an integer parameter lattice.
            # hmax=sqrt(3) gives exactly unit grid spacing; physical hmax's
            # established implementation and output contract stay unchanged.
            out=mesh_box(0.,Float64(counts[1]),0.,Float64(counts[2]),0.,Float64(counts[3]);hmax=sqrt(3.))
            for node in axes(out.coords,2),axis in 1:3
                index=Int(out.coords[axis,node])+1
                out.coords[axis,node]=grid_axes[axis][index]
            end
            diagnostic=validate(out)
            diagnostic.ok || throw(ArgumentError(
                "$caller: requested physical grid is not representable as a valid mesh — "*
                join(diagnostic.messages,"; ")))
            return out
        end
        factor=bound/target
        counts=ntuple(3) do axis
            proposed=_ceil_count3(counts[axis]*factor,caller,"axis $axis refined interval")
            max(proposed,_checked_add3(caller,"axis refined interval",counts[axis],1))
        end
    end
    throw(ArgumentError("$caller: rounded grid did not meet the metric edge bound within 32 count adjustments"))
end

"""
    mesh_box_regions(regions::AbstractVector{BoxRegion}, metric::Metric3;
                     max_metric_edge=1.0, max_nodes=10_000_000,
                     max_tets=60_000_000) -> Mesh

Construct a conforming assembly of axis-aligned boxes for one **constant**
symmetric positive-definite metric. Every Region face is a plane of the shared
grid. Each gap between planes receives separate directional interval counts;
every actual stored tetrahedron edge satisfies
`sqrt((b-a)' * M * (b-a)) ≤ max_metric_edge`.

The existing Kuhn fill retains the original Region semantics: the last box
containing a cell takes priority, `tag == 0` carves a void, and cells outside all
boxes are dropped. Adjacent materials share identical node and face rows;
`tet_tag` carries their original tags. Input boxes and the metric are unchanged.

Resource limits conservatively bound the entire shared bounding grid before
compaction, including outside and void cells. Sparse, widely separated boxes
can therefore exceed a limit even when their final mesh would be small.
Original face planes are retained exactly, with signed zero treated as one
geometric plane. Coordinate-resolution and empty-domain failures are explicit.

This positional overload uses an axis-aligned constant-tensor grid and a
conservative offdiagonal edge certificate. It does not sample variable fields
or rotate the grid to principal directions. The physical `hmax` overload keeps
its existing contract.
"""
function mesh_box_regions(regions::AbstractVector{BoxRegion},metric::Metric3;
                          max_metric_edge::Real=1.0,
                          max_nodes=PIPELINE_DEFAULT_MAX_NODES,
                          max_tets=PIPELINE_DEFAULT_MAX_TETS)
    caller="mesh_box_regions constant metric"
    isempty(regions) && throw(ArgumentError("$caller: no regions"))
    target=_finite3(max_metric_edge,caller,"max_metric_edge")
    target>0 || throw(ArgumentError("$caller: max_metric_edge must be positive"))
    node_limit=pipeline_limit(max_nodes,caller,"max_nodes",typemax(Int32))
    tet_limit=pipeline_limit(max_tets,caller,"max_tets",typemax(Int32))
    faces=(Float64[],Float64[],Float64[])
    for region in regions
        for (axis,first,last) in ((1,region.x0,region.x1),
                                 (2,region.y0,region.y1),(3,region.z0,region.z1))
            push!(faces[axis],first==0 ? 0. : first,last==0 ? 0. : last)
        end
    end
    planes=ntuple(axis->sort!(unique(faces[axis])),3)
    origin=(0.,0.,0.)
    interval_counts=ntuple(3) do axis
        counts=Int[]
        for segment in 1:length(planes[axis])-1
            span=planes[axis][segment+1]-planes[axis][segment]
            isfinite(span) && span>0 || throw(ArgumentError(
                "$caller: each shared axis gap must have a finite positive extent"))
            displacement=ntuple(other->other==axis ? span : 0.,3)
            directional=_metric_endpoint_displacement(metric,origin,displacement,caller)
            push!(counts,_ceil_count3((directional/target)*sqrt(3.),caller,
                                     "axis $axis interval"))
        end
        counts
    end
    envelope=Metric3(metric.m11,metric.m22,metric.m33,
        abs(metric.m12),abs(metric.m13),abs(metric.m23))
    for attempt in 1:32
        totals=_metric_region_counts(interval_counts,caller)
        _metric_box_grid_counts(totals,node_limit,tet_limit,caller)
        axes_and_gaps=ntuple(axis->_metric_region_axis(planes[axis],
            interval_counts[axis],totals[axis],caller,axis),3)
        grid_axes=ntuple(axis->axes_and_gaps[axis][1],3)
        plane_indices=ntuple(axis->axes_and_gaps[axis][2],3)
        gaps=ntuple(axis->axes_and_gaps[axis][3],3)
        bound=_metric_box_edge_bound(envelope,gaps,caller)
        if bound<=target
            parameter_regions=BoxRegion[BoxRegion(
                _metric_region_parameter(planes[1],plane_indices[1],region.x0),
                _metric_region_parameter(planes[1],plane_indices[1],region.x1),
                _metric_region_parameter(planes[2],plane_indices[2],region.y0),
                _metric_region_parameter(planes[2],plane_indices[2],region.y1),
                _metric_region_parameter(planes[3],plane_indices[3],region.z0),
                _metric_region_parameter(planes[3],plane_indices[3],region.z1),
                region.tag) for region in regions]
            out=mesh_box_regions(parameter_regions;hmax=sqrt(3.))
            for node in axes(out.coords,2),axis in 1:3
                out.coords[axis,node]=grid_axes[axis][Int(out.coords[axis,node])+1]
            end
            diagnostic=validate(out)
            diagnostic.ok || throw(ArgumentError(
                "$caller: requested physical grid is not representable as a valid mesh — "*
                join(diagnostic.messages,"; ")))
            return out
        end
        factor=bound/target
        for counts in interval_counts,segment in eachindex(counts)
            proposed=_ceil_count3(counts[segment]*factor,caller,"refined axis interval")
            counts[segment]=max(proposed,_checked_add3(caller,"refined axis interval",
                                                     counts[segment],1))
        end
    end
    throw(ArgumentError("$caller: rounded grid did not meet the metric edge bound within 32 count adjustments"))
end
