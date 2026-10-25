"""
    GeoUnrolled

`.geo_unrolled` serialization — the output of Gmsh 4.15.2's `-0` batch
action. Every model entity emits one `.geo` statement in a fixed order:
`cl__N` mesh-size variables in first-use order, points, then each curve
followed by its `Transfinite Curve`, each surface preceded by its lazily
emitted `Curve Loop`s and followed by `Transfinite Surface`/`Recombine
Surface`/`In` embeds, each volume preceded by its `Surface Loop`s and
followed by `Transfinite Volume`, then physical groups. Periodic
constraints and `Nurbs`/`BSpline Surface`/`Nurbs Surface` records do not
serialize into the unrolled format upstream and are dropped, matching the
observed 4.15.2 output.
"""
module GeoUnrolled

using Printf: @sprintf
using ..Model: GeoModel
using ..Elements: MixedMesh

export write_geo_unrolled

const _UNROLLED_UNSET_SIZE=1e22

_gnum(x)=@sprintf("%.16g",Float64(x))

const _DIM_NAMES=("Point","Curve","Surface","Volume")

function _clvar!(clmap,clnames,size)
    get!(clmap,size) do
        push!(clnames,_gnum(size))
        length(clnames)
    end
end

function _unrolled_curve_statement(tag,kind,curves,control_points)
    kw=kind==:circle ? "Circle" : kind==:ellipse ? "Ellipse" :
       kind==:spline ? "Spline" : kind==:bspline ? "BSpline" :
       kind==:bezier ? "Bezier" : kind==:polyline ? "Spline" : "Line"
    verts=kind==:line ? collect(get(curves,tag,(0,0))) :
          get(control_points,tag,Int[])
    return "$kw($tag) = {$(join(verts,", "))};"
end

# `Physical Dim = {…}` with an auto-assigned tag does not serialize
# upstream — only named groups and explicitly numbered ones survive.
function _unrolled_physical_statement(dim,tag,name)
    label=name===nothing ? "($tag)" : "(\"$(replace(name,"\\"=>"\\\\","\""=>"\\\""))\")"
    return "Physical $(_DIM_NAMES[dim+1])$label = "
end

"""
    write_geo_unrolled(io::IO, m::GeoModel)
    write_geo_unrolled(path::AbstractString, m::GeoModel)

Serialize `m` as a `.geo_unrolled` script matching the Gmsh 4.15.2 `-0`
output documented in the module header.
"""
function write_geo_unrolled(io::IO,m::GeoModel)
    clmap=Dict{Float64,Int}();clnames=String[]
    point_lines=String[]
    for tag in sort!(collect(keys(m.points)))
        x,y,z=m.points[tag]
        size=get(m.point_size,tag,_UNROLLED_UNSET_SIZE)
        n=_clvar!(clmap,clnames,size)
        push!(point_lines,"Point($tag) = {$(_gnum(x)), $(_gnum(y)), " *
            "$(_gnum(z)), cl__$n};")
    end
    for (i,value) in enumerate(clnames)
        println(io,"cl__$i = $value;")
    end
    for line in point_lines
        println(io,line)
    end
    for tag in sort!(collect(keys(m.curves)))
        kind=get(m.curve_types,tag,:line)
        # `GEdge::writeGEO` returns early for `MSH_SEGM_BND_LAYER` — a
        # boundary-layer chapeau has no `.geo` statement.
        kind in (:nurbs,:boundary_layer) && continue
        println(io,_unrolled_curve_statement(
            tag,kind,m.curves,m.curve_control_points))
        tf=get(m.meshing.transfinite_curves,tag,nothing)
        if tf!==nothing
            kindname=tf.kind==:bump ? "Bump" : "Progression"
            println(io,"Transfinite Curve {$tag} = $(tf.num_nodes) " *
                "Using $kindname $(_gnum(tf.coef));")
        end
    end
    emitted_loops=Set{Int}()
    for tag in sort!(collect(keys(m.surfaces)))
        stype=get(m.surface_types,tag,:plane)
        stype in (:nurbs,:bspline_surface,:nurbs_surface) && continue
        for lp in m.surfaces[tag]
            loop=abs(lp)
            loop in emitted_loops && continue
            push!(emitted_loops,loop)
            println(io,"Curve Loop($loop) = {$(join(m.loops[loop],", "))};")
        end
        kw=stype==:plane ? "Plane Surface" : "Surface"
        println(io,"$kw($tag) = {$(join(m.surfaces[tag],", "))};")
        tf=get(m.meshing.transfinite_surfaces,tag,nothing)
        if tf!==nothing
            isempty(tf.corners) ?
                println(io,"Transfinite Surface {$tag};") :
                println(io,"Transfinite Surface {$tag} = {$(join(tf.corners,","))};")
        end
        if haskey(m.meshing.recombine,(2,tag)) ||
           haskey(m.meshing.recombine,(3,tag))
            println(io,"Recombine Surface {$tag};")
        end
        for (edim,etag) in get(m.embeds,(2,tag),NTuple{2,Int}[])
            println(io,"$(_DIM_NAMES[edim+1]) {$etag} In Surface {$tag};")
        end
    end
    for tag in sort!(collect(keys(m.loops)))
        tag in emitted_loops && continue
        println(io,"Curve Loop($tag) = {$(join(m.loops[tag],", "))};")
    end
    emitted_surface_loops=Set{Int}()
    for tag in sort!(collect(keys(m.volumes)))
        for lp in m.volumes[tag]
            loop=abs(lp)
            loop in emitted_surface_loops && continue
            push!(emitted_surface_loops,loop)
            println(io,"Surface Loop($loop) = {$(join(m.surface_loops[loop],", "))};")
        end
        println(io,"Volume($tag) = {$(join(m.volumes[tag],", "))};")
        corners=get(m.meshing.transfinite_volumes,tag,nothing)
        corners===nothing || isempty(corners) ||
            println(io,"Transfinite Volume {$tag} = {$(join(corners,","))};")
        for (edim,etag) in get(m.embeds,(3,tag),NTuple{2,Int}[])
            println(io,"$(_DIM_NAMES[edim+1]) {$etag} In Volume {$tag};")
        end
    end
    for tag in sort!(collect(keys(m.surface_loops)))
        tag in emitted_surface_loops && continue
        println(io,"Surface Loop($tag) = {$(join(m.surface_loops[tag],", "))};")
    end
    for key in sort!(collect(keys(m.physical)))
        dim,tag=key
        entities=m.physical[key]
        isempty(entities) && continue
        name=get(m.physical_names,key,nothing)
        println(io,_unrolled_physical_statement(dim,tag,name) *
            "{$(join(entities,", "))};")
    end
    return nothing
end

function write_geo_unrolled(path::AbstractString,m::GeoModel)
    open(path,"w") do io
        write_geo_unrolled(io,m)
    end
    return path
end

"""
    write_geo_unrolled(io::IO, m::MixedMesh)

Mesh-input `-0`: upstream emits only the model's discrete point records
(without a size field — `cl__N = 1e+22`) and warns that non-built-in
entities are debugging output. Element and higher-dimension entity records
do not serialize; partition/ghost metadata is ignored.
"""
function write_geo_unrolled(io::IO,m::MixedMesh)
    data=m.entity_data
    data===nothing && return nothing
    points=sort!([entity for entity in values(data.entities) if entity.dim==0];
        by=e->Int(e.tag))
    isempty(points) && return nothing
    println(io,"cl__1 = $(_gnum(_UNROLLED_UNSET_SIZE));")
    for entity in points
        b=entity.bbox
        println(io,"Point($(entity.tag)) = {$(_gnum(b[1])), $(_gnum(b[2])), " *
            "$(_gnum(b[3])), cl__1};")
    end
    return nothing
end

function write_geo_unrolled(path::AbstractString,m::MixedMesh)
    open(path,"w") do io
        write_geo_unrolled(io,m)
    end
    return path
end

end # module GeoUnrolled
