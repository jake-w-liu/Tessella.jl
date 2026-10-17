# `occ.importShapes` / `OCC_Internals::importShapes` — ISO-10303-21 STEP
# solids imported as native OCC records (points/lines/wires/plane
# faces/shells/volumes) rather than the triangulated `import_step`
# meshes. The BRep module supplies the parser and the face/bound/edge
# graph walk; this file converts each closed shell's rings into shared
# vertex/edge entities — identical vertices and edges publish once, so
# the shell validator's two-surface-per-curve rule holds by
# construction. Bounded to planar faces (STEP `PLANE` or untyped
# `FACE`); curved faces reject loudly until curved STEP→record
# conversion lands.

# Solid-level dedupe state: vertex positions → point tags, unordered
# vertex-tag pairs → (line id, stored stop tag).
struct _OCCImportDedupe
    ppos::Vector{NTuple{3,Float64}}
    ptag::Vector{Int}
    edge::Dict{NTuple{2,Int},Tuple{UInt64,Int}}
end
_OCCImportDedupe()=_OCCImportDedupe(NTuple{3,Float64}[],Int[],
                                   Dict{NTuple{2,Int},Tuple{UInt64,Int}}())

function _occ_import_point!(m,d:: _OCCImportDedupe,p)
    for i in eachindex(d.ppos)
        _model_occ_vertices_coincide(d.ppos[i],p) && return d.ptag[i]
    end
    t=_model_occ_add_point!(m,p[1],p[2],p[3])
    push!(d.ppos,p); push!(d.ptag,t)
    return t
end

# One ring (ordered positions, no duplicated closing point) → wire tag.
# `flip` winds the loop backwards (hole wires are consumed CW relative
# to the face normal by `add_plane_surface` convention checks — STEP
# bounds already carry the correct orientation so `flip` stays false).
function _occ_import_wire!(m,d,kernel,ring,flip,caller)
    verts=[_occ_import_point!(m,d,p) for p in ring]
    uses=Tuple{UInt64,Bool}[]
    K=length(verts)
    for j in 1:K
        a,b=verts[j],verts[mod1(j+1,K)]
        key=a<b ? (a,b) : (b,a)
        eu=get(d.edge,key,nothing)
        if eu===nothing
            et=_model_occ_add_line!(m,a,b)
            eu=(kernel.bindings[(1,et)],b)
            d.edge[key]=eu
        end
        rev=eu[2]!=b
        push!(uses,flip ? (eu[1],!rev) : (eu[1],rev))
    end
    flip && (uses=reverse(uses))
    return _occ_wire_from_uses!(m,uses,caller)
end

# One STEP face: rings from `_step_face_parts` (outer first). Planar
# surfaces only — a non-plane surface record rejects.
function _occ_import_step_face(m,d,kernel,entities,fe,curvekinds,caller)
    rings,surf,sense,straight=BRep._step_face_parts(entities,fe,
                                                  curvekinds,caller)
    rings===nothing && throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) uses an unsupported loop kind"))
    surf===nothing || surf.kind=="PLANE" || throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) uses $(surf.kind) — only planar " *
        "faces import currently"))
    straight || throw(ArgumentError(
        "$caller: $(fe.kind) #$(fe.id) has non-line edges — only " *
        "straight-edge faces import currently"))
    wires=[_occ_import_wire!(m,d,kernel,r,false,caller) for r in rings]
    return _model_occ_add_plane_surface!(m,wires)
end

# One solid root → [(3, tag)]. Each CLOSED_SHELL becomes a surface
# loop; `BREP_WITH_VOIDS`-style later shells carry role −1 (cavities).
function _occ_import_step_solid(m,brep,entities,caller)
    kernel=m.occ_kernel
    length(brep.args)>=2 || throw(ArgumentError(
        "$caller: $(brep.kind) #$(brep.id) has no shell"))
    d=_OCCImportDedupe(); curvekinds=Set{String}()
    shell_tags=Int[]
    for (si,sarg) in enumerate(brep.args[2:end])
        items=sarg isa Vector ? sarg : Any[sarg]
        for item in items
            se=BRep._deref_entity(entities,item)
            se===nothing && throw(ArgumentError(
                "$caller: $(brep.kind) #$(brep.id) references a missing shell"))
            se.kind=="CLOSED_SHELL" || throw(ArgumentError(
                "$caller: $(brep.kind) #$(brep.id) references $(se.kind) " *
                "#$(se.id) where a CLOSED_SHELL was expected"))
            i=BRep._skip_name(se.args)
            flist=se.args[i]
            flist isa Vector || throw(ArgumentError(
                "$caller: CLOSED_SHELL #$(se.id) face list is malformed"))
            face_tags=Int[]
            for fref in flist
                fe=BRep._deref_entity(entities,fref)
                (fe!==nothing &&
                 fe.kind in BRep._STEP_BREP_FACE_KINDS) || throw(
                    ArgumentError(
                        "$caller: CLOSED_SHELL #$(se.id) member is not a face"))
                push!(face_tags,
                      _occ_import_step_face(m,d,kernel,entities,fe,
                                            curvekinds,caller))
            end
            isempty(face_tags) && throw(ArgumentError(
                "$caller: CLOSED_SHELL #$(se.id) is empty"))
            push!(shell_tags,
                  _model_occ_add_surface_loop!(m,face_tags))
        end
    end
    isempty(shell_tags) && throw(ArgumentError(
        "$caller: $(brep.kind) #$(brep.id) has no shells"))
    return [(3,_model_occ_add_volume!(m,shell_tags))]
end

function _model_occ_import_shapes(m::GeoModel,path;highest_dim_only=true,
                                  caller="OCC.importShapes")
    isfile(path) || throw(ArgumentError("$caller: missing file $path"))
    source=read(path,String)
    occursin("ISO-10303-21",source) || throw(ArgumentError(
        "$caller: $path is not an ISO-10303-21 STEP file"))
    entities=BRep.parse_step_entities(source)
    breps=sort!([e for e in values(entities)
                 if e.kind in BRep._STEP_BREP_SOLID_KINDS]; by=e->e.id)
    isempty(breps) && throw(ArgumentError(
        "$caller: no MANIFOLD_SOLID_BREP/FACETED_BREP/BREP_WITH_VOIDS " *
        "found in $path"))
    out=Tuple{Int,Int}[]
    for brep in breps
        append!(out,_occ_import_step_solid(m,brep,entities,caller))
    end
    return out
end
