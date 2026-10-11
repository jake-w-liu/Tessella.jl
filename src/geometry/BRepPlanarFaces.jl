# A holed Face is triangulated from its literal planar wires. Its triangles
# then enter the unchanged closed-shell edge/component certificate. Existing
# one-wire Faces keep their original ear-clipping path and emission order.
function _brep_face_outer_wire_order(rings,caller)
    all(ring->length(ring)>=3,rings) || throw(ArgumentError("$caller: Face wire has fewer than three vertices"))
    _check_poly_face(first(rings),0,caller)
    normal=_newell_normal(first(rings));drop=argmax(abs.(normal))
    axes=filter(!=(drop),collect(1:3))
    projected=[[(point[axes[1]],point[axes[2]]) for point in ring] for ring in rings]
    # FACE_BOUND does not designate the outer wire. A unique containment
    # certificate determines its role; the subsequent domain certificate
    # rejects every crossing, touching, nested or nonplanar configuration.
    candidates=[index for index in eachindex(projected) if
        all(other->other==index || _brep_point_in_wire(first(projected[other]),projected[index]),eachindex(projected))]
    length(candidates)==1 || throw(ArgumentError("$caller: Face wires do not identify one unique outer boundary"))
    index=only(candidates)
    return vcat(rings[index:index],rings[1:index-1],rings[index+1:end])
end

function _polyhedral_surface_faces(face_rings,roles,caller)
    all(rings->length(rings)==1,face_rings) &&
        return _polyhedral_surface([only(rings) for rings in face_rings],roles,caller)
    polygons=Vector{NTuple{3,Float64}}[];polygon_roles=Int[]
    for (rings,role) in zip(face_rings,roles)
        if length(rings)==1
            push!(polygons,only(rings));push!(polygon_roles,role)
        else
            triangles=_brep_planar_face_holes(rings,caller)
            append!(polygons,triangles);append!(polygon_roles,fill(role,length(triangles)))
        end
    end
    return _polyhedral_surface(polygons,polygon_roles,caller)
end

function _brep_wire_on_segment(a,b,p)
    orient2(a,b,p)==0 &&
        min(a[1],b[1])<=p[1]<=max(a[1],b[1]) &&
        min(a[2],b[2])<=p[2]<=max(a[2],b[2])
end

function _brep_wires_intersect(a,b,c,d)
    first=orient2(a,b,c)
    second=orient2(a,b,d)
    third=orient2(c,d,a)
    fourth=orient2(c,d,b)
    ((first>0 && second<0 || first<0 && second>0) &&
     (third>0 && fourth<0 || third<0 && fourth>0)) ||
        _brep_wire_on_segment(a,b,c) || _brep_wire_on_segment(a,b,d) ||
        _brep_wire_on_segment(c,d,a) || _brep_wire_on_segment(c,d,b)
end

function _brep_point_in_wire(point,ring)
    inside=false
    for index in eachindex(ring)
        a=ring[index];b=ring[mod1(index+1,length(ring))]
        if (a[2]>point[2])!=(b[2]>point[2])
            orientation=orient2(a,b,point)
            (b[2]>a[2] ? orientation>0 : orientation<0) && (inside=!inside)
        end
    end
    return inside
end

function _brep_planar_face_holes(rings,caller)
    length(rings)>1 || throw(ArgumentError("$caller: holed Face needs an outer and inner wire"))
    outer=first(rings);_check_poly_face(outer,0,caller)
    normal=_brep_vunit(_newell_normal(outer),caller)
    scale=max(_poly_scale(vcat(rings...)),1.0)
    for ring in rings
        length(ring)>=3 && allunique(ring) || throw(ArgumentError(
            "$caller: planar Face wire repeats vertices or has fewer than three vertices"))
        _check_coplanar(ring,normal,first(outer),1e-9*scale,0,caller)
    end
    drop=argmax(abs.(normal));axes=filter(!=(drop),collect(1:3))
    projected=[[(point[axes[1]],point[axes[2]]) for point in ring] for ring in rings]
    points=vcat(projected...)
    allunique(points) || throw(ArgumentError("$caller: planar Face wires touch or repeat a projected vertex"))
    edges=Tuple{Int,Int,NTuple{2,Float64},NTuple{2,Float64}}[]
    for (wire,ring) in enumerate(projected),index in eachindex(ring)
        a=ring[index];b=ring[mod1(index+1,length(ring))]
        a!=b || throw(ArgumentError("$caller: planar Face wire contains a zero-length edge"))
        push!(edges,(wire,index,a,b))
    end
    for first_index in 1:length(edges)-1,second_index in first_index+1:length(edges)
        wire,index,a,b=edges[first_index];other_wire,other_index,c,d=edges[second_index]
        if wire==other_wire && (mod1(index+1,length(projected[wire]))==other_index ||
                              mod1(other_index+1,length(projected[wire]))==index)
            shared=a==c || a==d ? a : b
            left=a==shared ? b : a;right=c==shared ? d : c
            (_brep_wire_on_segment(shared,left,right) || _brep_wire_on_segment(shared,right,left)) &&
                throw(ArgumentError("$caller: adjacent planar Face wire edges overlap"))
            continue
        end
        _brep_wires_intersect(a,b,c,d) && throw(ArgumentError(
            "$caller: planar Face wires cross, overlap, or touch"))
    end
    area=[_signed_area2(first.(ring),last.(ring)) for ring in projected]
    all(value->isfinite(value) && value!=0,area) || throw(ArgumentError(
        "$caller: planar Face wire has zero or non-finite area"))
    for hole in 2:length(projected)
        signbit(area[hole])!=signbit(area[1]) || throw(ArgumentError(
            "$caller: planar Face inner wire winding must oppose its outer wire"))
        _brep_point_in_wire(first(projected[hole]),first(projected)) || throw(ArgumentError(
            "$caller: planar Face inner wire lies outside its declared outer wire"))
        for other in 2:length(projected)
            hole==other && continue
            _brep_point_in_wire(first(projected[hole]),projected[other]) && throw(ArgumentError(
                "$caller: planar Face inner wires cannot nest"))
        end
    end
    expected=abs(first(area))-sum(abs,area[2:end])
    expected>0 || throw(ArgumentError("$caller: planar Face holes remove the entire outer wire"))
    segments=Tuple{Int,Int}[];offset=0
    for ring in projected
        for index in eachindex(ring)
            push!(segments,(offset+index,offset+mod1(index+1,length(ring))))
        end
        offset+=length(ring)
    end
    triangulation=Mesh2D.constrained_delaunay(first.(points),last.(points),segments)
    mesh=Mesh2D.to_mesh(triangulation;interior=Mesh2D.classify_interior(triangulation))
    source=Dict(projected[wire][index]=>rings[wire][index] for wire in eachindex(rings) for index in eachindex(rings[wire]))
    triangles=Vector{NTuple{3,Float64}}[];measured=0.0
    for cell in eachcol(mesh.tris)
        positions=[(mesh.coords[1,row],mesh.coords[2,row]) for row in cell]
        all(point->haskey(source,point),positions) || throw(ArgumentError(
            "$caller: planar Face triangulation adds an unbound source vertex"))
        measured+=abs(_signed_area2(first.(positions),last.(positions)))
        triangle=[source[point] for point in positions]
        signbit(first(area)) && reverse!(triangle)
        push!(triangles,triangle)
    end
    !isempty(triangles) && abs(measured-expected)<=1e-10*max(expected,scale*scale) ||
        throw(ArgumentError("$caller: planar Face triangulation does not cover its exact wire domain"))
    return triangles
end
