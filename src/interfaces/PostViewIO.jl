# Gmsh `.pos` list-format post-processing views (PViewDataList parsed I/O).
# Included inside `module Post` from Post.jl. The parsed format is the text
# `View "name" { TYPE(coords){values}; ... };` layout that Gmsh 4.15.2 writes
# for list-based views and reads for `Merge`/`PostView` background fields.
# Element records whose names carry a suffix (`SL2`, `ST2`, ...) and views with
# `INTERPOLATION_SCHEME` records are the Gmsh high-order/interpolation-matrix
# extension: Gmsh parses them for adapted visualization, and this reader
# resolves them exactly for `PostViewField` background-field evaluation.

using ..SizeField: PostViewField, PostViewAnisoField, AbstractField
using ..Elements: lagrange_nodes
using LinearAlgebra: inv
using Printf: @sprintf

export PosElement, PosScheme, PosText, PosView, read_pos, write_pos,
       postview_field

const _POS_ELEMENT_ARITY=(point=1,line=2,triangle=3,quadrangle=4,
                        tetrahedron=4,hexahedron=8,prism=6,pyramid=5)

# Two-letter record name → (kind, nodes, components), matching the Gmsh
# grammar's SP..TY table (scalar S, vector V, tensor T times each shape).
const _POS_RECORD_TYPES=let table=Dict{String,Tuple{Symbol,Int,Int}}()
    for (shape,kind,arity) in (("P",:point,1),("L",:line,2),("T",:triangle,3),
                              ("Q",:quadrangle,4),("S",:tetrahedron,4),
                              ("H",:hexahedron,8),("I",:prism,6),("Y",:pyramid,5))
        for (head,ncomp) in (("S",1),("V",3),("T",9))
            table[head*shape]=(kind,arity,ncomp)
        end
    end
    table
end
const _POS_RECORD_NAME=let table=Dict{Tuple{Symbol,Int},String}()
    for (name,(kind,arity,ncomp)) in _POS_RECORD_TYPES
        table[(kind,ncomp)]=name
    end
    table
end

const _POS_MAX_NAME_BYTES=1<<20
const _POS_MAX_STRING_BYTES=1<<20

"""
    PosScheme(coefval, expval)
    PosScheme(coefval, expval, coefgeo, expgeo)

An element interpolation scheme as bound by a `.pos`
`INTERPOLATION_SCHEME` record or a high-order `X2` record name. `coefval` is
the `nval × nmonomial` coefficient matrix and `expval` the `nmonomial × d`
(`d ≤ 3`) exponent table, so value basis function `i` evaluates as
`Σⱼ coefval[i,j] · u^a v^b w^c` over exponent row `(a,b,c)` — Gmsh's
`computeShapeFunctions` convention. On `:pyramid` elements every exponent row
instead evaluates Bergot-style as `u^a v^b w^c · (1-w)^(max(a,b)-a-b)`. The
two-matrix form keeps the family's canonical first-order geometry; the
four-matrix form adds `coefgeo`/`expgeo` for a curved geometric map whose node
count is `size(coefgeo,1)`. Inputs are copied.
"""
struct PosScheme
    coefval::Matrix{Float64}
    expval::Matrix{Float64}
    coefgeo::Matrix{Float64}
    expgeo::Matrix{Float64}
    function PosScheme(coefval::AbstractMatrix{<:Real},
                       expval::AbstractMatrix{<:Real},
                       coefgeo::AbstractMatrix{<:Real}=
                           Matrix{Float64}(undef,0,0),
                       expgeo::AbstractMatrix{<:Real}=
                           Matrix{Float64}(undef,0,0))
        caller="PosScheme"
        cv=_pos_scheme_matrix(coefval,caller,"value coefficients")
        ev=_pos_scheme_matrix(expval,caller,"value exponents")
        size(cv,2)==size(ev,1) || throw(ArgumentError(
            "$caller: value coefficient column count $(size(cv,2)) must " *
            "equal the value exponent row count $(size(ev,1))"))
        1<=size(ev,2)<=3 || throw(ArgumentError(
            "$caller: value exponents must have 1 to 3 columns " *
            "(got $(size(ev,2)))"))
        geo_given=!isempty(coefgeo) || !isempty(expgeo)
        if geo_given
            (isempty(coefgeo) || isempty(expgeo)) && throw(ArgumentError(
                "$caller: geometry coefficients and exponents must be " *
                "provided together"))
            cg=_pos_scheme_matrix(coefgeo,caller,"geometry coefficients")
            eg=_pos_scheme_matrix(expgeo,caller,"geometry exponents")
            size(cg,2)==size(eg,1) || throw(ArgumentError(
                "$caller: geometry coefficient column count $(size(cg,2)) " *
                "must equal the geometry exponent row count $(size(eg,1))"))
            1<=size(eg,2)<=3 || throw(ArgumentError(
                "$caller: geometry exponents must have 1 to 3 columns " *
                "(got $(size(eg,2)))"))
        else
            cg=Matrix{Float64}(undef,0,0)
            eg=Matrix{Float64}(undef,0,0)
        end
        return new(cv,ev,cg,eg)
    end
end

Base.:(==)(a::PosScheme,b::PosScheme)=
    a.coefval==b.coefval && a.expval==b.expval &&
    a.coefgeo==b.coefgeo && a.expgeo==b.expgeo
Base.isequal(a::PosScheme,b::PosScheme)=a==b
Base.hash(a::PosScheme,h::UInt)=
    hash(a.expgeo,hash(a.coefgeo,hash(a.expval,hash(a.coefval,h))))

function _pos_scheme_matrix(matrix,caller::AbstractString,
                            what::AbstractString)::Matrix{Float64}
    matrix isa AbstractMatrix{<:Real} || throw(ArgumentError(
        "$caller: $what must be a real matrix"))
    size(matrix,1)>=1 || throw(ArgumentError("$caller: $what is empty"))
    any(v->v isa Bool,matrix) && throw(ArgumentError(
        "$caller: $what must not contain Bool entries"))
    converted=try
        Matrix{Float64}(matrix)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $what must be Float64-representable"))
    end
    all(isfinite,converted) || throw(ArgumentError(
        "$caller: $what must be finite"))
    return converted
end

_pos_scheme_geom_nodes(scheme::PosScheme,kind::Symbol)=
    isempty(scheme.coefgeo) ? _POS_ELEMENT_ARITY[kind] : size(scheme.coefgeo,1)
_pos_scheme_value_nodes(scheme::PosScheme)=size(scheme.coefval,1)

# Order-2 MSH type per record family, matching PViewDataList::setOrder2.
const _POS_ORDER2_MSH=(line=8,triangle=9,quadrangle=10,tetrahedron=11,
                       hexahedron=12,prism=13,pyramid=14)

# Bergot exponent triples (i,j,k), k ≤ 2−max(i,j), for the 14-node order-2
# pyramid — pyramidalBasis enumerates i,j outer / k inner in this order.
const _POS_PYR2_EXPS=NTuple{3,Float64}[
    (i,j,k) for i in 0.0:2.0 for j in 0.0:2.0
            for k in 0.0:(2.0-max(i,j))]

@inline function _pos_monomial(exps::NTuple{3,Float64},u,v,w,pyramid::Bool)
    value=u^exps[1]*v^exps[2]*w^exps[3]
    pyramid || return value
    # pyramidalBasis::VDM/evaluate clamps 1−w away from the apex singularity.
    onemw=max(1e-14,1.0-w)
    return value*onemw^(max(exps[1],exps[2])-exps[1]-exps[2])
end

# Nodal order-2 scheme for a record family: invert the Vandermonde of the
# family monomials on the Gmsh-ordered reference nodes (lagrange_nodes), so
# Nᵢ(nodeⱼ) = δᵢⱼ. Pyramid rows use the Bergot factor like Gmsh.
function _pos_order2_scheme_build(kind::Symbol)::PosScheme
    pyramid=kind===:pyramid
    nodes=lagrange_nodes(_POS_ORDER2_MSH[kind])
    n=size(nodes,2)
    exps=if pyramid
        _POS_PYR2_EXPS
    else
        NTuple{3,Float64}[(Float64(q[1]),Float64(q[2]),Float64(q[3]))
                          for q in _pos_order2_lattice(kind)]
    end
    vandermonde=Matrix{Float64}(undef,n,n)
    @inbounds for point in 1:n
        u,v,w=nodes[1,point],nodes[2,point],nodes[3,point]
        for term in 1:n
            vandermonde[point,term]=_pos_monomial(exps[term],u,v,w,pyramid)
        end
    end
    # Nᵢ(r) = Σₜ coef[i,t]·Mₜ(r) must satisfy Nᵢ(nodeⱼ) = δᵢⱼ, i.e.
    # coef·Vᵀ = I — the coefficient matrix is (Vᵀ)⁻¹, not V⁻¹.
    coef=Matrix(inv(vandermonde'))
    exponent_matrix=Matrix{Float64}(undef,n,3)
    @inbounds for term in 1:n
        exponent_matrix[term,1]=exps[term][1]
        exponent_matrix[term,2]=exps[term][2]
        exponent_matrix[term,3]=exps[term][3]
    end
    return PosScheme(coef,exponent_matrix,coef,exponent_matrix)
end

# Integer exponent lattice per family at order 2. The sets are the family
# monomials (triangle/tetrahedron total-degree, quadrangle/hexahedron tensor,
# prism product, pyramid uses the Bergot triples above); their order only
# affects the row order of the coefficient matrix, not the evaluated basis.
function _pos_order2_lattice(kind::Symbol)
    kind===:line && return NTuple{3,Int}[(0,0,0),(2,0,0),(1,0,0)]
    kind===:triangle && return NTuple{3,Int}[
        (0,0,0),(2,0,0),(0,2,0),(1,1,0),(1,0,0),(0,1,0)]
    kind===:quadrangle && return NTuple{3,Int}[
        (i,j,0) for j in 0:2 for i in 0:2]
    kind===:tetrahedron && return NTuple{3,Int}[
        (i,j,k) for k in 0:2 for j in 0:(2-k) for i in 0:(2-j-k)]
    kind===:hexahedron && return NTuple{3,Int}[
        (i,j,k) for k in 0:2 for j in 0:2 for i in 0:2]
    kind===:prism && return NTuple{3,Int}[
        (i,j,k) for k in 0:2 for j in 0:2 for i in 0:(2-j)]
    throw(ErrorException("PosScheme: unsupported order-2 family $kind"))
end

# Cached shared order-2 schemes — one instance per family so views bind the
# same object and element `===` comparisons stay cheap and deterministic.
const _POS_ORDER2_SCHEMES=Dict{Symbol,PosScheme}(
    kind=>_pos_order2_scheme_build(kind) for
    kind in (:line,:triangle,:quadrangle,:tetrahedron,
             :hexahedron,:prism,:pyramid))

@inline _pos_order2_scheme(kind::Symbol)=_POS_ORDER2_SCHEMES[kind]

@inline function _pos_suffix_valid(suffix::AbstractString)
    return all(c->('a'<=c<='z' || 'A'<=c<='Z' ||
                   '0'<=c<='9' || c=='_'),suffix)
end

"""
    PosElement(kind, coords, values; suffix="", scheme=nothing)

One element record of a `.pos` view. `kind` is `:point`, `:line`, `:triangle`,
`:quadrangle`, `:tetrahedron`, `:hexahedron`, `:prism`, or `:pyramid`;
`coords` is a finite `3×ngeom` matrix and `values` a finite `c×nval` matrix
(single step) or `c×nval×steps` array with `c ∈ (1, 3, 9)` scalar/vector/tensor
components, ordered `[component, node, step]` as in the file.

With `scheme === nothing` the element is a canonical first-order record and
both widths are the family arity (1, 2, 3, 4, 4, 8, 6, 5). With a
[`PosScheme`](@ref), `coords` holds `size(scheme.coefgeo,1)` geometry nodes
(the family arity when the scheme has no geometry matrices) and `values`
holds `size(scheme.coefval,1)` value degrees of freedom. `suffix` is the
record-name suffix (`"2"` for the built-in order-2 elements); a nonempty
suffix without an explicit scheme binds the shared order-2 scheme for the
family (point records keep the suffix verbatim without a scheme). Inputs are
copied.
"""
struct PosElement
    kind::Symbol
    coords::Matrix{Float64}
    values::Array{Float64,3}
    suffix::String
    scheme::Union{Nothing,PosScheme}
    function PosElement(kind::Symbol,coords::AbstractMatrix{<:Real},
                        values::AbstractArray{<:Real};
                        suffix::AbstractString="",
                        scheme::Union{Nothing,PosScheme}=nothing)
        caller="PosElement"
        haskey(_POS_ELEMENT_ARITY,kind) || throw(ArgumentError(
            "$caller: kind $(repr(kind)) is unsupported (expected one of " *
            "point/line/triangle/quadrangle/tetrahedron/hexahedron/prism/pyramid)"))
        (scheme===nothing || kind!==:point) || throw(ArgumentError(
            "$caller: point records cannot carry an interpolation scheme " *
            "(Gmsh defines no high-order point basis)"))
        _pos_suffix_valid(suffix) || throw(ArgumentError(
            "$caller: record-name suffix $(repr(suffix)) may only contain " *
            "letters, digits, and '_'"))
        resolved=scheme
        if resolved===nothing && !isempty(suffix) && kind!==:point
            resolved=_pos_order2_scheme(kind)
        end
        ngeom=resolved===nothing ? _POS_ELEMENT_ARITY[kind] :
            _pos_scheme_geom_nodes(resolved,kind)
        nval=resolved===nothing ? _POS_ELEMENT_ARITY[kind] :
            _pos_scheme_value_nodes(resolved)
        size(coords,1)==3 || throw(ArgumentError(
            "$caller: coords must be a 3×$ngeom matrix"))
        size(coords,2)==ngeom || throw(ArgumentError(
            "$caller: kind $kind requires exactly $ngeom geometry nodes " *
            "(got $(size(coords,2)))"))
        ndims(values) in (2,3) || throw(ArgumentError(
            "$caller: values must be a c×$nval matrix or c×$nval×steps array"))
        size(values,1) in (1,3,9) || throw(ArgumentError(
            "$caller: values must have 1, 3, or 9 components per node"))
        size(values,2)==nval || throw(ArgumentError(
            "$caller: values must have exactly $nval value nodes " *
            "(got $(size(values,2)))"))
        ndims(values)==3 && size(values,3)==0 && throw(ArgumentError(
            "$caller: empty time-step array"))
        any(v->v isa Bool,coords) && throw(ArgumentError(
            "$caller: coordinates must not contain Bool entries"))
        any(v->v isa Bool,values) && throw(ArgumentError(
            "$caller: values must not contain Bool entries"))
        C=try
            Float64.(coords)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "$caller: coordinates must be Float64-representable"))
        end
        V=try
            ndims(values)==3 ? Array{Float64,3}(Float64.(values)) :
                reshape(Float64.(values),size(values,1),nval,1)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$caller: values must be Float64-representable"))
        end
        all(isfinite,C) || throw(ArgumentError(
            "$caller: coordinates must be finite"))
        all(isfinite,V) || throw(ArgumentError(
            "$caller: values must be finite"))
        return new(kind,C,V,String(suffix),resolved)
    end
end

"""
    PosText(dim, coords, style, strings)

A `.pos` text record. `dim` is `2` (`T2`, screen-space x/y plus a style value,
stored with `z = 0`) or `3` (`T3`, world-space x/y/z plus a style value);
`strings` is the nonempty list of displayed strings. String contents may not
contain `"` characters, so the record always re-serializes losslessly.
"""
struct PosText
    dim::Int8
    coords::NTuple{3,Float64}
    style::Float64
    strings::Vector{String}
    function PosText(dim::Integer,coords,style::Real,strings)
        caller="PosText"
        (strings isa AbstractVector || strings isa Tuple) || throw(ArgumentError(
            "$caller: strings must be a vector or tuple of strings"))
        dim isa Bool && throw(ArgumentError("$caller: dim must not be Bool"))
        (dim==2 || dim==3) || throw(ArgumentError(
            "$caller: dim must be 2 or 3 (got $dim)"))
        (coords isa NTuple{3} || (coords isa AbstractVector && length(coords)==3)) ||
            throw(ArgumentError("$caller: coords must be a 3-coordinate point"))
        point=ntuple(i->coords[i],3)
        any(c->c isa Bool,point) && throw(ArgumentError(
            "$caller: coordinates must not contain Bool entries"))
        xyz=try
            (Float64(point[1]),Float64(point[2]),Float64(point[3]))
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$caller: coordinates must be Float64-representable"))
        end
        all(isfinite,xyz) || throw(ArgumentError(
            "$caller: coordinates must be finite"))
        dim==2 && xyz[3]!=0 && throw(ArgumentError(
            "$caller: a 2-D text record requires z == 0"))
        style isa Bool && throw(ArgumentError("$caller: style must not be Bool"))
        s=try
            Float64(style)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$caller: style must be Float64-representable"))
        end
        isfinite(s) || throw(ArgumentError("$caller: style must be finite"))
        isempty(strings) && throw(ArgumentError(
            "$caller: a text record needs at least one string"))
        text=String[]
        for (i,item) in enumerate(strings)
            item isa AbstractString || throw(ArgumentError(
                "$caller: string $i must be a string"))
            item=String(item)
            occursin('"',item) && throw(ArgumentError(
                "$caller: string $i must not contain '\"' characters"))
            push!(text,item)
        end
        return new(Int8(dim),xyz,s,text)
    end
end

# Family precedence of an `INTERPOLATION_SCHEME` record, evaluated against
# the element families present at the record's file position. Gmsh's grammar
# scans line, triangle, quadrangle, tetrahedron, pyramid, prism, hexahedron
# for two-matrix schemes and line, triangle, quadrangle, tetrahedron,
# hexahedron for four-matrix schemes (the curved form cannot bind pyramids or
# prisms).
const _POS_SCHEME_PRECEDENCE2=(:line,:triangle,:quadrangle,:tetrahedron,
                               :pyramid,:prism,:hexahedron)
const _POS_SCHEME_PRECEDENCE4=(:line,:triangle,:quadrangle,:tetrahedron,
                               :hexahedron)

# Resolve the family→scheme binding of a view in file order: `kinds[i]` and
# `suffixes[i]` describe element record i, `positions[s]` counts how many
# element records precede scheme s. First binding wins per family, matching
# `PViewData::setInterpolationMatrices`. Returns `(bound, source)` where
# `source[family]` is the 1-based `interpolations` index that bound the
# family, or `0` when a suffixed/scheme-carrying element bound it.
function _pos_resolve_schemes(kinds::AbstractVector,
                              suffixes::AbstractVector,
                              element_schemes::AbstractVector,
                              interpolations::AbstractVector,
                              positions::AbstractVector,
                              caller::AbstractString)
    bound=Dict{Symbol,PosScheme}()
    source=Dict{Symbol,Int}()
    seen=Dict{Symbol,Bool}()
    nschemes=length(interpolations)
    next_scheme=1
    nelem=length(kinds)
    for pos in 0:nelem
        while next_scheme<=nschemes && positions[next_scheme]==pos
            matrices=interpolations[next_scheme]
            precedence=length(matrices)==4 ?
                _POS_SCHEME_PRECEDENCE4 : _POS_SCHEME_PRECEDENCE2
            family=findfirst(family->get(seen,family,false),precedence)
            if family!==nothing
                target=precedence[family]
                if !haskey(bound,target)
                    bound[target]=length(matrices)==4 ?
                        PosScheme(matrices[1],matrices[2],
                                  matrices[3],matrices[4]) :
                        PosScheme(matrices[1],matrices[2])
                    source[target]=next_scheme
                end
            end
            next_scheme+=1
        end
        pos==nelem && break
        kind=kinds[pos+1]
        seen[kind]=true
        element_scheme=element_schemes[pos+1]
        if !haskey(bound,kind) && kind!==:point &&
           (!isempty(suffixes[pos+1]) || element_scheme!==nothing)
            bound[kind]=element_scheme===nothing ?
                _pos_order2_scheme(kind) : element_scheme
            source[kind]=0
        end
    end
    return bound,source
end

"""
    PosView(name; tag=0, time=Float64[], elements=PosElement[],
            texts=PosText[], interpolations=[], scheme_positions=[])

A parsed `.pos` view: `elements` are [`PosElement`](@ref) records in file
order, `time` lists the `TIME` step values (required to have one entry per
data step when given, padded with `0:n-1` when omitted), `tag` is the view tag
used by `PostView` field resolution (`read_pos` assigns `1:length(views)`).
`interpolations` preserves `INTERPOLATION_SCHEME` matrices as read
(list-of-rows matrices, 2 or 4 per scheme) and `scheme_positions` the number
of element records preceding each scheme (defaults: all schemes last, which
matches `read_pos` output only when every scheme followed all elements).

Each `INTERPOLATION_SCHEME` binds the first element family present at its
file position in Gmsh's precedence order (line, triangle, quadrangle,
tetrahedron, [pyramid, prism for two-matrix schemes], hexahedron); the first
binding wins per family, and a suffixed record name (`SL2`, …) binds its
family to the built-in order-2 nodal scheme at its own position. Elements of
a bound family are retagged with the resolved [`PosScheme`](@ref), which
revalidates their geometry/value widths — the scheme, not the name suffix,
defines the record layout. Elements carrying an explicit scheme that
contradicts the view-level binding raise `ArgumentError`. Every element
record must share the same step count.
"""
struct PosView
    name::String
    tag::Int
    time::Vector{Float64}
    elements::Vector{PosElement}
    texts::Vector{PosText}
    interpolations::Vector{Vector{Matrix{Float64}}}
    scheme_positions::Vector{Int}
    function PosView(name::AbstractString;
                     tag::Integer=0,
                     time::AbstractVector{<:Real}=Float64[],
                     elements::AbstractVector=PosElement[],
                     texts::AbstractVector=PosText[],
                     interpolations::AbstractVector=
                         Vector{Matrix{Float64}}[],
                     scheme_positions::AbstractVector{<:Integer}=Int[])
        caller="PosView"
        view_name=String(name)
        occursin('"',view_name) && throw(ArgumentError(
            "$caller: view name must not contain '\"' characters"))
        tag isa Bool && throw(ArgumentError("$caller: tag must not be Bool"))
        view_tag=try
            Int(tag)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$caller: tag exceeds the platform Int range"))
        end
        0<=view_tag<=typemax(Int32) || throw(ArgumentError(
            "$caller: tag must lie in 0:$(typemax(Int32)) (got $tag)"))
        any(v->v isa Bool,time) && throw(ArgumentError(
            "$caller: time must not contain Bool entries"))
        steps_time=try
            Float64.(collect(time))
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$caller: time values must be Float64-representable"))
        end
        all(isfinite,steps_time) || throw(ArgumentError(
            "$caller: time values must be finite"))
        element_list=PosElement[]
        for (i,element) in enumerate(elements)
            element isa PosElement || throw(ArgumentError(
                "$caller: element $i is not a PosElement"))
            push!(element_list,element)
        end
        text_list=PosText[]
        for (i,text) in enumerate(texts)
            text isa PosText || throw(ArgumentError(
                "$caller: text $i is not a PosText"))
            push!(text_list,text)
        end
        scheme_list=Vector{Matrix{Float64}}[]
        for (i,scheme) in enumerate(interpolations)
            length(scheme) in (2,4) || throw(ArgumentError(
                "$caller: interpolation scheme $i must have 2 or 4 matrices " *
                "(got $(length(scheme)))"))
            matrices=Matrix{Float64}[]
            for (j,matrix) in enumerate(scheme)
                matrix isa AbstractMatrix{<:Real} || throw(ArgumentError(
                    "$caller: interpolation scheme $i matrix $j must be a real matrix"))
                size(matrix,1)>=1 && size(matrix,2)>=1 || throw(ArgumentError(
                    "$caller: interpolation scheme $i matrix $j is empty"))
                any(v->v isa Bool,matrix) && throw(ArgumentError(
                    "$caller: interpolation scheme $i matrix $j must not contain Bool entries"))
                converted=try
                    Float64.(matrix)
                catch err
                    err isa InterruptException && rethrow()
                    throw(ArgumentError(
                        "$caller: interpolation matrix $j of scheme $i must be " *
                        "Float64-representable"))
                end
                all(isfinite,converted) || throw(ArgumentError(
                    "$caller: interpolation scheme $i matrix $j must be finite"))
                push!(matrices,converted)
            end
            push!(scheme_list,matrices)
        end
        positions=Int[]
        if isempty(scheme_positions)
            resize!(positions,length(scheme_list))
            fill!(positions,length(element_list))
        else
            length(scheme_positions)==length(scheme_list) || throw(
                ArgumentError(
                    "$caller: scheme_positions has " *
                    "$(length(scheme_positions)) entries for " *
                    "$(length(scheme_list)) interpolation schemes"))
            previous=0
            for (i,position) in enumerate(scheme_positions)
                position isa Bool && throw(ArgumentError(
                    "$caller: scheme position $i must not be Bool"))
                converted=try
                    Int(position)
                catch err
                    err isa InterruptException && rethrow()
                    throw(ArgumentError(
                        "$caller: scheme position $i is outside Int bounds"))
                end
                0<=converted<=length(element_list) || throw(ArgumentError(
                    "$caller: scheme position $i ($converted) is outside " *
                    "0:$(length(element_list))"))
                converted>=previous || throw(ArgumentError(
                    "$caller: scheme positions must be nondecreasing"))
                previous=converted
                push!(positions,converted)
            end
        end
        # Rebind schemes in file order, then retag every element to its
        # family's resolved scheme. Elements without a scheme inherit the
        # family binding (rebuild revalidates their record widths, like
        # PViewDataList::_stat decoding the flat lists); elements carrying an
        # explicit scheme must agree with the family binding.
        resolved0,_=_pos_resolve_schemes(
            Symbol[e.kind for e in element_list],
            String[e.suffix for e in element_list],
            Union{Nothing,PosScheme}[nothing for _ in element_list],
            scheme_list,positions,caller)
        # An element-level scheme only steers the binding when it differs
        # from what the view-level records already bind for its family —
        # elements retagged by family binding (the read_pos path) carry the
        # bound scheme as data, not as a binding directive.
        resolved,resolved_source=_pos_resolve_schemes(
            Symbol[e.kind for e in element_list],
            String[e.suffix for e in element_list],
            Union{Nothing,PosScheme}[
                let scheme=e.scheme
                    (scheme===nothing || scheme==get(resolved0,e.kind,nothing)) ?
                        nothing : scheme
                end for e in element_list],
            scheme_list,positions,caller)
        first_of_family=Dict{Symbol,Int}()
        for (i,element) in enumerate(element_list)
            get!(first_of_family,element.kind,i)
        end
        for (i,element) in enumerate(element_list)
            family_scheme=get(resolved,element.kind,nothing)
            if element.scheme===nothing && family_scheme!==nothing
                element_list[i]=PosElement(
                    element.kind,element.coords,element.values;
                    suffix=element.suffix,scheme=family_scheme)
            elseif element.scheme!==nothing && family_scheme!==nothing &&
                   element.scheme!=family_scheme
                throw(ArgumentError(
                    "$caller: element $i carries a scheme that contradicts " *
                    "the view-level $(element.kind) binding"))
            end
        end
        # Custom (non-order-2) schemes only serialize through
        # INTERPOLATION_SCHEME records: an element-bound custom scheme would
        # not rebind on re-parse (record names only convey order-2), so views
        # must carry the scheme explicitly via `interpolations` +
        # `scheme_positions`. A nonempty name suffix on a custom-scheme
        # element is safe only when the binding scheme record precedes the
        # element — otherwise the suffix's order-2 binding would shadow it.
        for (i,element) in enumerate(element_list)
            element.scheme===nothing && continue
            element.kind===:point && continue
            element.scheme===_pos_order2_scheme(element.kind) && continue
            source=get(resolved_source,element.kind,0)
            (source==0 || (!isempty(element.suffix) &&
                           i<=positions[source])) && throw(ArgumentError(
                "$caller: element $i carries a custom scheme not bound by " *
                "an INTERPOLATION_SCHEME entry preceding it — add the " *
                "scheme via `interpolations`/`scheme_positions` (suffixed " *
                "names only serialize built-in order-2 elements)"))
        end
        if isempty(element_list)
            length(steps_time)<=typemax(Int32) || throw(ArgumentError(
                "$caller: time-step count exceeds the Int32 limit"))
            normalized_time=steps_time
        else
            steps=size(element_list[1].values,3)
            for (i,element) in enumerate(element_list)
                size(element.values,3)==steps || throw(ArgumentError(
                    "$caller: element $i has $(size(element.values,3)) time " *
                    "steps but the first element has $steps"))
            end
            normalized_time=if isempty(steps_time)
                Float64[Float64(i-1) for i in 1:steps]
            else
                length(steps_time)==steps || throw(ArgumentError(
                    "$caller: TIME lists $(length(steps_time)) steps but the " *
                    "view data has $steps"))
                steps_time
            end
        end
        return new(view_name,view_tag,normalized_time,element_list,text_list,
                   scheme_list,positions)
    end
end

"""
    PosView(view::View; tag=0) -> PosView

Convert an owned scalar node view to a `.pos` view of `SP` point elements, one
per node. Lossless for scalar data; the `.pos` serialization keeps every node
as an independent point record, like Gmsh's `SP` list.
"""
function PosView(view::View;tag::Integer=0)
    n=size(view.coords,2)
    elements=Vector{PosElement}(undef,n)
    @inbounds for i in 1:n
        elements[i]=PosElement(:point,view.coords[:,i:i],
                               reshape([view.values[i]],1,1))
    end
    return PosView(view.name;tag=tag,elements=elements)
end

# ── Reader ──────────────────────────────────────────────────────────────────

@inline function _pos_limit(value,caller::AbstractString,
                            name::AbstractString;ceiling::Int=typemax(Int))
    value isa Integer || throw(ArgumentError(
        "$caller: $name must be an integer"))
    value isa Bool && throw(ArgumentError(
        "$caller: $name must not be Bool"))
    converted=try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$caller: $name is outside Int bounds"))
    end
    0<=converted<=ceiling || throw(ArgumentError(
        "$caller: $name must lie in 0:$ceiling (got $value)"))
    return converted
end

struct _PosReadLimits
    views::Int
    elements::Int
    values::Int
    steps::Int
    name_bytes::Int
    string_bytes::Int
    matrices::Int
end

struct _PosToken
    kind::Symbol
    text::String
    value::Float64
    at::Int
end

mutable struct _PosLexer
    data::Vector{UInt8}
    pos::Int
end

@inline _pos_is_space(c::UInt8)=
    c==0x20 || c==0x09 || c==0x0a || c==0x0d || c==0x0b || c==0x0c
@inline _pos_is_letter(c::UInt8)=
    (0x41<=c<=0x5a) || (0x61<=c<=0x7a) || c==0x5f
@inline _pos_is_digit(c::UInt8)=0x30<=c<=0x39

function _pos_skip_space!(lx::_PosLexer,caller::AbstractString)
    data=lx.data;n=length(data);i=lx.pos
    while i<=n
        c=data[i]
        if _pos_is_space(c)
            i+=1;continue
        elseif c==0x23          # '#'
            while i<=n && data[i]!=0x0a
                i+=1
            end
            continue
        elseif c==0x2f          # '/'
            i+1<=n || throw(ArgumentError(
                "$caller: unexpected '/' at byte $i"))
            if data[i+1]==0x2f
                i+=2
                while i<=n && data[i]!=0x0a
                    i+=1
                end
                continue
            elseif data[i+1]==0x2a
                i+=2;closed=false
                while i+1<=n
                    if data[i]==0x2a && data[i+1]==0x2f
                        i+=2;closed=true;break
                    end
                    i+=1
                end
                closed || throw(ArgumentError(
                    "$caller: unterminated '/*' comment"))
                continue
            end
            throw(ArgumentError("$caller: unexpected '/' at byte $i"))
        else
            break
        end
    end
    lx.pos=i
    return nothing
end

function _pos_scan_number(data::Vector{UInt8},n::Int,i::Int)
    j=i
    if data[j]==0x2b || data[j]==0x2d     # '+' '-'
        j+=1
    end
    digits=0
    while j<=n && _pos_is_digit(data[j])
        j+=1;digits+=1
    end
    if j<=n && data[j]==0x2e              # '.'
        j+=1
        while j<=n && _pos_is_digit(data[j])
            j+=1;digits+=1
        end
    end
    digits==0 && return i,"",false
    if j<=n && (data[j]==0x65 || data[j]==0x45)   # 'e' 'E'
        k=j+1
        if k<=n && (data[k]==0x2b || data[k]==0x2d)
            k+=1
        end
        exponent=0
        while k<=n && _pos_is_digit(data[k])
            k+=1;exponent+=1
        end
        exponent>0 && (j=k)
    end
    return j,String(data[i:j-1]),true
end

function _pos_next!(lx::_PosLexer,caller::AbstractString)::_PosToken
    _pos_skip_space!(lx,caller)
    data=lx.data;n=length(data);i=lx.pos
    i>n && return _PosToken(:eof,"",0.0,n+1)
    c=data[i]
    if c==0x28;lx.pos=i+1;return _PosToken(:lparen,"",0.0,i);end       # (
    if c==0x29;lx.pos=i+1;return _PosToken(:rparen,"",0.0,i);end       # )
    if c==0x7b;lx.pos=i+1;return _PosToken(:lbrace,"",0.0,i);end       # {
    if c==0x7d;lx.pos=i+1;return _PosToken(:rbrace,"",0.0,i);end       # }
    if c==0x2c;lx.pos=i+1;return _PosToken(:comma,"",0.0,i);end        # ,
    if c==0x3b;lx.pos=i+1;return _PosToken(:semi,"",0.0,i);end         # ;
    if c==0x22 || c==0x27             # '"' '\''
        qchar=c;j=i+1
        while j<=n && data[j]!=qchar
            j+=1
        end
        j<=n || throw(ArgumentError(
            "$caller: unterminated string at byte $i"))
        lx.pos=j+1
        return _PosToken(:str,String(data[i+1:j-1]),0.0,i)
    end
    if _pos_is_letter(c)
        j=i+1
        while j<=n && (_pos_is_letter(data[j]) || _pos_is_digit(data[j]))
            j+=1
        end
        lx.pos=j
        return _PosToken(:ident,String(data[i:j-1]),0.0,i)
    end
    if _pos_is_digit(c) || c==0x2b || c==0x2d || c==0x2e
        stop,text,isnumber=_pos_scan_number(data,n,i)
        isnumber || throw(ArgumentError(
            "$caller: unexpected byte 0x$(string(c,base=16,pad=2)) at byte $i"))
        lx.pos=stop
        value=parse(Float64,text)
        isfinite(value) || throw(ArgumentError(
            "$caller: non-finite number '$text' at byte $i"))
        return _PosToken(:num,text,value,i)
    end
    throw(ArgumentError(
        "$caller: unexpected byte 0x$(string(c,base=16,pad=2)) at byte $i"))
end

mutable struct _PosReader
    lx::_PosLexer
    pending::Union{Nothing,_PosToken}
    caller::String
end

@inline function _pos_peek!(r::_PosReader)
    r.pending===nothing && (r.pending=_pos_next!(r.lx,r.caller))
    return r.pending
end
@inline function _pos_take!(r::_PosReader)
    token=_pos_peek!(r)
    r.pending=nothing
    return token
end
function _pos_expect!(r::_PosReader,kind::Symbol)
    t=_pos_take!(r)
    t.kind==kind || throw(ArgumentError(
        "$(r.caller): expected '$kind' at byte $(t.at), got $(_pos_token_desc(t))"))
    return t
end
_pos_token_desc(t::_PosToken)=
    t.kind==:ident ? "'$(t.text)'" :
    t.kind==:num ? "number '$(t.text)'" :
    t.kind==:str ? "string" : "'$(t.kind)'"

function _pos_number_tail!(r::_PosReader,close::Symbol,
                           caller::AbstractString)::Vector{Float64}
    values=Float64[]
    t=_pos_take!(r)
    t.kind==close && return values
    while true
        t.kind==:num || throw(ArgumentError(
            "$caller: expected a number at byte $(t.at), got $(_pos_token_desc(t))"))
        push!(values,t.value)
        t=_pos_take!(r)
        t.kind==:comma || break
        t=_pos_take!(r)
    end
    t.kind==close || throw(ArgumentError(
        "$caller: expected ',' or a closing delimiter at byte $(t.at), got " *
        _pos_token_desc(t)))
    return values
end

function _pos_number_list!(r::_PosReader,open::Symbol,close::Symbol,
                           caller::AbstractString)::Vector{Float64}
    t=_pos_take!(r)
    t.kind==open || throw(ArgumentError(
        "$caller: expected '$open' at byte $(t.at), got $(_pos_token_desc(t))"))
    return _pos_number_tail!(r,close,caller)
end

function _pos_string_list!(r::_PosReader,limits::_PosReadLimits,
                           caller::AbstractString)::Vector{String}
    _pos_expect!(r,:lbrace)
    strings=String[]
    t=_pos_take!(r)
    t.kind==:rbrace && return strings
    while true
        (t.kind==:str || t.kind==:ident) || throw(ArgumentError(
            "$caller: expected a string at byte $(t.at), got $(_pos_token_desc(t))"))
        ncodeunits(t.text)<=limits.string_bytes || throw(ArgumentError(
            "$caller: string at byte $(t.at) exceeds max_string_bytes=" *
            "$(limits.string_bytes)"))
        push!(strings,t.text)
        t=_pos_take!(r)
        t.kind==:comma || break
        t=_pos_take!(r)
    end
    t.kind==:rbrace || throw(ArgumentError(
        "$caller: expected ',' or '}' at byte $(t.at), got $(_pos_token_desc(t))"))
    return strings
end

# One `{ {..}, {..}, .. }` interpolation matrix: rows are the inner lists, a
# bare scalar is a single-entry row (Gmsh's RecursiveListOfListOfDouble).
function _pos_matrix!(r::_PosReader,limits::_PosReadLimits,
                      caller::AbstractString)::Matrix{Float64}
    _pos_expect!(r,:lbrace)
    rows=Vector{Float64}[]
    entries=0
    while true
        t=_pos_take!(r)
        t.kind==:rbrace && break
        if t.kind==:lbrace
            row=_pos_number_tail!(r,:rbrace,caller)
            push!(rows,row)
        elseif t.kind==:num
            push!(rows,[t.value])
        else
            throw(ArgumentError(
                "$caller: expected an interpolation row at byte $(t.at), got " *
                _pos_token_desc(t)))
        end
        entries+=length(rows[end])
        entries<=limits.matrices || throw(ArgumentError(
            "$caller: interpolation entries exceed " *
            "max_interpolation_entries=$(limits.matrices)"))
        t=_pos_take!(r)
        t.kind==:comma && continue
        t.kind==:rbrace && break
        throw(ArgumentError(
            "$caller: expected ',' or '}' at byte $(t.at), got $(_pos_token_desc(t))"))
    end
    isempty(rows) && throw(ArgumentError(
        "$caller: empty interpolation matrix"))
    ncols=length(rows[1])
    for row in rows
        length(row)==ncols || throw(ArgumentError(
            "$caller: interpolation matrix rows must have equal length"))
    end
    matrix=Matrix{Float64}(undef,length(rows),ncols)
    @inbounds for i in eachindex(rows), j in 1:ncols
        matrix[i,j]=rows[i][j]
    end
    return matrix
end

function _pos_unknown_record(word::AbstractString,caller::AbstractString,at::Int)
    return ArgumentError(
        "$caller: unknown view record '$word' at byte $at")
end

# Raw element record collected before the family schemes are resolved: the
# record layout (geometry/value node counts) is only known once the view's
# `INTERPOLATION_SCHEME` records and `X2` name suffixes have been bound.
struct _PosRawElement
    kind::Symbol
    suffix::String
    ncomp::Int
    coords::Vector{Float64}
    values::Vector{Float64}
    at::Int
end

function _pos_parse_view(r::_PosReader,tag::Int,limits::_PosReadLimits,
                         caller::AbstractString)::PosView
    t=_pos_take!(r)
    (t.kind==:str || t.kind==:ident) || throw(ArgumentError(
        "$caller: expected a view name at byte $(t.at), got $(_pos_token_desc(t))"))
    name=t.text
    ncodeunits(name)<=limits.name_bytes || throw(ArgumentError(
        "$caller: view name exceeds max_name_bytes=$(limits.name_bytes)"))
    occursin('"',name) && throw(ArgumentError(
        "$caller: view name contains an unbalanced quote"))
    _pos_expect!(r,:lbrace)
    raw_elements=_PosRawElement[]
    texts=PosText[]
    schemes=Vector{Matrix{Float64}}[]
    scheme_positions=Int[]
    time=Float64[]
    seen_time=false
    while true
        t=_pos_take!(r)
        t.kind==:rbrace && break
        t.kind==:ident || throw(ArgumentError(
            "$caller: expected a view record at byte $(t.at), got " *
            _pos_token_desc(t)))
        word=t.text
        if word=="TIME"
            seen_time && throw(ArgumentError(
                "$caller: duplicate TIME record at byte $(t.at)"))
            seen_time=true
            open_token=_pos_take!(r)
            (open_token.kind==:lbrace || open_token.kind==:lparen) ||
                throw(ArgumentError(
                    "$caller: TIME must be followed by '{' or '(' at byte " *
                    "$(open_token.at)"))
            time=_pos_number_tail!(
                r,open_token.kind==:lbrace ? :rbrace : :rparen,caller)
            isempty(time) && throw(ArgumentError(
                "$caller: empty TIME record at byte $(t.at)"))
            length(time)<=limits.steps || throw(ArgumentError(
                "$caller: TIME step count exceeds max_steps=$(limits.steps)"))
            _pos_expect!(r,:semi)
        elseif word=="T2" || word=="T3"
            nargs=word=="T2" ? 3 : 4
            args=_pos_number_list!(r,:lparen,:rparen,caller)
            length(args)==nargs || throw(ArgumentError(
                "$caller: $word expects $nargs arguments (got $(length(args)))"))
            strings=_pos_string_list!(r,limits,caller)
            push!(texts,PosText(word=="T2" ? Int8(2) : Int8(3),
                                word=="T2" ? (args[1],args[2],0.0) :
                                (args[1],args[2],args[3]),args[end],strings))
            _pos_expect!(r,:semi)
        elseif word=="INTERPOLATION_SCHEME"
            matrices=Matrix{Float64}[]
            while _pos_peek!(r).kind==:lbrace
                push!(matrices,_pos_matrix!(r,limits,caller))
            end
            length(matrices) in (2,4) || throw(ArgumentError(
                "$caller: INTERPOLATION_SCHEME expects 2 or 4 matrices at " *
                "byte $(t.at) (got $(length(matrices)))"))
            push!(schemes,matrices)
            push!(scheme_positions,length(raw_elements))
            _pos_expect!(r,:semi)
        else
            ncodeunits(word)>=2 || throw(_pos_unknown_record(word,caller,t.at))
            spec=get(_POS_RECORD_TYPES,word[1:2],nothing)
            spec===nothing && throw(_pos_unknown_record(word,caller,t.at))
            kind,_,ncomp=spec
            flat_coords=_pos_number_list!(r,:lparen,:rparen,caller)
            length(flat_coords)%3==0 || throw(ArgumentError(
                "$caller: $word record at byte $(t.at) has " *
                "$(length(flat_coords)) coordinates, not a multiple of 3"))
            flat_values=_pos_number_list!(r,:lbrace,:rbrace,caller)
            (length(flat_coords)+length(flat_values))<=limits.values ||
                throw(ArgumentError(
                    "$caller: element data exceeds max_values=$(limits.values)"))
            push!(raw_elements,_PosRawElement(
                kind,word[3:end],ncomp,flat_coords,flat_values,t.at))
            _pos_expect!(r,:semi)
        end
        length(raw_elements)+length(texts)+length(schemes)<=limits.elements ||
            throw(ArgumentError(
                "$caller: view records exceed max_elements=$(limits.elements)"))
    end
    _pos_expect!(r,:semi)
    return _pos_assemble_view(name,tag,time,raw_elements,texts,
                              schemes,scheme_positions,caller)
end

# Decode raw element records after the whole view body is known: resolve the
# family scheme bindings in file order (INTERPOLATION_SCHEME positions plus
# order-2 name suffixes), then validate each record's geometry and value
# widths against its family's bound scheme — PViewDataList::_stat semantics.
function _pos_assemble_view(name::AbstractString,tag::Int,
                            time::Vector{Float64},raw::Vector{_PosRawElement},
                            texts::Vector{PosText},
                            schemes::Vector{Vector{Matrix{Float64}}},
                            scheme_positions::Vector{Int},
                            caller::AbstractString)::PosView
    nelem=length(raw)
    kinds=Symbol[raw[i].kind for i in 1:nelem]
    suffixes=String[raw[i].suffix for i in 1:nelem]
    bound,_=_pos_resolve_schemes(
        kinds,suffixes,Union{Nothing,PosScheme}[nothing for _ in 1:nelem],
        schemes,scheme_positions,caller)
    elements=Vector{PosElement}(undef,nelem)
    nsteps=0
    @inbounds for i in 1:nelem
        record=raw[i]
        scheme=get(bound,record.kind,nothing)
        ngeom=scheme===nothing ? _POS_ELEMENT_ARITY[record.kind] :
            _pos_scheme_geom_nodes(scheme,record.kind)
        nval=scheme===nothing ? _POS_ELEMENT_ARITY[record.kind] :
            _pos_scheme_value_nodes(scheme)
        word=_POS_RECORD_NAME[(record.kind,record.ncomp)]*record.suffix
        length(record.coords)==3ngeom || throw(ArgumentError(
            "$caller: $word record at byte $(record.at) expects $ngeom " *
            "geometry nodes ($(3ngeom) coordinates, got " *
            "$(length(record.coords)))"))
        width=nval*record.ncomp
        (!isempty(record.values) && length(record.values)%width==0) ||
            throw(ArgumentError(
                "$caller: $word record at byte $(record.at) has " *
                "$(length(record.values)) values, not a positive multiple " *
                "of nodes×components=$width"))
        steps=length(record.values)÷width
        nsteps==0 ? (nsteps=steps) : steps==nsteps || throw(ArgumentError(
            "$caller: $word record at byte $(record.at) has $steps time " *
            "steps but the view has $nsteps"))
        coords=Matrix{Float64}(undef,3,ngeom)
        copyto!(coords,record.coords)
        values=Array{Float64,3}(undef,record.ncomp,nval,steps)
        for s in 1:steps, j in 1:nval, c in 1:record.ncomp
            values[c,j,s]=record.values[(s-1)*width+(j-1)*record.ncomp+c]
        end
        elements[i]=PosElement(record.kind,coords,values;
                               suffix=record.suffix,scheme=scheme)
    end
    return PosView(name;tag=tag,time=time,elements=elements,texts=texts,
                   interpolations=schemes,scheme_positions=scheme_positions)
end

function _pos_parse(data::Vector{UInt8},limits::_PosReadLimits,
                    caller::AbstractString)::Vector{PosView}
    r=_PosReader(_PosLexer(data,1),nothing,caller)
    views=PosView[]
    while true
        t=_pos_take!(r)
        t.kind==:eof && break
        t.kind==:ident || throw(ArgumentError(
            "$caller: expected 'View' at byte $(t.at), got $(_pos_token_desc(t))"))
        t.text=="View" || throw(ArgumentError(
            "$caller: expected 'View' at byte $(t.at), got '$(t.text)'"))
        push!(views,_pos_parse_view(r,length(views)+1,limits,caller))
        length(views)<=limits.views || throw(ArgumentError(
            "$caller: view count exceeds max_views=$(limits.views)"))
    end
    return views
end

function _pos_read_bounded(io::IO,file_limit::Int,
                           caller::AbstractString)::Vector{UInt8}
    data=UInt8[]
    while !eof(io)
        room=file_limit-length(data)
        room<=0 && throw(ArgumentError(
            "$caller: file exceeds max_file_bytes=$file_limit"))
        chunk=read(io,min(1<<20,room))
        isempty(chunk) && break
        append!(data,chunk)
    end
    return data
end

function _read_pos(io::IO,limits::_PosReadLimits,file_limit::Int,
                   caller::AbstractString)
    data=_pos_read_bounded(io,file_limit,caller)
    return _pos_parse(data,limits,caller)
end

"""
    read_pos(io; kwargs...) -> Vector{PosView}
    read_pos(path; max_file_bytes=typemax(Int), kwargs...) -> Vector{PosView}

Read Gmsh parsed-format `.pos` views: `View "name" { ... };` blocks with
scalar/vector/tensor element records `SP..TP`, `SL..TL`, `ST..TT`, `SQ..TQ`,
`SS..TS`, `SH..TH`, `SI..TI`, `SY..TY`, `TIME{...}`/`TIME(...)` step values,
`T2`/`T3` text records, `INTERPOLATION_SCHEME` matrices, `//`, `/* */`, and
`#` comments. Returns the views in file order with tags `1:length(views)`.

The reader is strict: malformed syntax, unknown or higher-order (`SL2`, …)
record types, wrong node counts, value counts that are not a multiple of
`nodes×components`, inconsistent or missing time steps, non-finite numbers,
and resource-limit violations all raise `ArgumentError`. Resources are bounded
before record structures are allocated: `max_file_bytes` (path only, checked
against the file size), `max_views`, `max_elements` (records per view),
`max_values` (coordinate+value numbers per record), `max_steps` (TIME entries
per view), `max_name_bytes`, `max_string_bytes`, and
`max_interpolation_entries`.
"""
function read_pos(io::IO;
                  max_views=typemax(Int32),
                  max_elements=typemax(Int32),
                  max_values=typemax(Int32),
                  max_steps=typemax(Int32),
                  max_name_bytes=_POS_MAX_NAME_BYTES,
                  max_string_bytes=_POS_MAX_STRING_BYTES,
                  max_interpolation_entries=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_pos"
    limits=_PosReadLimits(
        _pos_limit(max_views,caller,"max_views";ceiling=Int(typemax(Int32))),
        _pos_limit(max_elements,caller,"max_elements";
                   ceiling=Int(typemax(Int32))),
        _pos_limit(max_values,caller,"max_values";ceiling=Int(typemax(Int32))),
        _pos_limit(max_steps,caller,"max_steps";ceiling=Int(typemax(Int32))),
        _pos_limit(max_name_bytes,caller,"max_name_bytes"),
        _pos_limit(max_string_bytes,caller,"max_string_bytes"),
        _pos_limit(max_interpolation_entries,caller,
                   "max_interpolation_entries";ceiling=Int(typemax(Int32))))
    file_limit=_pos_limit(max_file_bytes,caller,"max_file_bytes")
    return _read_pos(io,limits,file_limit,caller)
end

function read_pos(path::AbstractString;
                  max_views=typemax(Int32),
                  max_elements=typemax(Int32),
                  max_values=typemax(Int32),
                  max_steps=typemax(Int32),
                  max_name_bytes=_POS_MAX_NAME_BYTES,
                  max_string_bytes=_POS_MAX_STRING_BYTES,
                  max_interpolation_entries=typemax(Int32),
                  max_file_bytes=typemax(Int))
    caller="read_pos"
    limits=_PosReadLimits(
        _pos_limit(max_views,caller,"max_views";ceiling=Int(typemax(Int32))),
        _pos_limit(max_elements,caller,"max_elements";
                   ceiling=Int(typemax(Int32))),
        _pos_limit(max_values,caller,"max_values";ceiling=Int(typemax(Int32))),
        _pos_limit(max_steps,caller,"max_steps";ceiling=Int(typemax(Int32))),
        _pos_limit(max_name_bytes,caller,"max_name_bytes"),
        _pos_limit(max_string_bytes,caller,"max_string_bytes"),
        _pos_limit(max_interpolation_entries,caller,
                   "max_interpolation_entries";ceiling=Int(typemax(Int32))))
    file_limit=_pos_limit(max_file_bytes,caller,"max_file_bytes")
    isfile(path) || throw(ArgumentError("$caller: missing regular file $path"))
    filesize(path)<=file_limit || throw(ArgumentError(
        "$caller: file exceeds max_file_bytes=$file_limit"))
    return open(path,"r") do io
        _read_pos(io,limits,file_limit,caller)
    end
end

# ── Writer ──────────────────────────────────────────────────────────────────

@inline _pos_fmt(x::Float64)=@sprintf("%.16g",x)

function _pos_write_element(io::IO,element::PosElement)
    name=_POS_RECORD_NAME[(element.kind,size(element.values,1))]
    suffix=element.suffix
    if isempty(suffix) && element.scheme!==nothing &&
       element.scheme===_pos_order2_scheme(element.kind)
        # Built-in order-2 elements serialize through their suffixed name,
        # which rebinds the shared scheme on re-parse.
        suffix="2"
    end
    print(io,name,suffix,'(')
    for j in axes(element.coords,2)
        j>1 && print(io,',')
        print(io,_pos_fmt(element.coords[1,j]),',',
                _pos_fmt(element.coords[2,j]),',',
                _pos_fmt(element.coords[3,j]))
    end
    print(io,"){")
    ncomp,arity,steps=size(element.values)
    width=arity*ncomp
    @inbounds for s in 1:steps, j in 1:arity, c in 1:ncomp
        index=(s-1)*width+(j-1)*ncomp+c
        index>1 && print(io,',')
        print(io,_pos_fmt(element.values[c,j,s]))
    end
    println(io,"};")
    return nothing
end

function _pos_write_interpolation(io::IO,scheme::AbstractVector)
    print(io,"INTERPOLATION_SCHEME")
    for matrix in scheme
        print(io,'{')
        for i in axes(matrix,1)
            i>first(axes(matrix,1)) && print(io,',')
            print(io,'{',join((_pos_fmt(v) for v in matrix[i,:]),','),'}')
        end
        print(io,'}')
    end
    println(io,';')
    return nothing
end

function _pos_write_view(io::IO,view::PosView)
    nsteps=isempty(view.elements) ? length(view.time) :
        size(view.elements[1].values,3)
    print(io,"View \"",view.name,"\" {\n")
    nsteps>1 && println(io,"TIME{",join((_pos_fmt(t) for t in view.time),','),
                        "};")
    nschemes=length(view.interpolations)
    next_scheme=1
    nelem=length(view.elements)
    # INTERPOLATION_SCHEME records interleave at their recorded element
    # positions so a re-parse binds the same families — PosView validated
    # that element schemes agree with this binding order.
    for index in 0:nelem
        while next_scheme<=nschemes &&
              view.scheme_positions[next_scheme]==index
            _pos_write_interpolation(io,view.interpolations[next_scheme])
            next_scheme+=1
        end
        index==nelem && break
        _pos_write_element(io,view.elements[index+1])
    end
    for text in view.texts
        quoted=join((string('"',s,'"') for s in text.strings),',')
        if text.dim==2
            println(io,"T2(",_pos_fmt(text.coords[1]),",",
                    _pos_fmt(text.coords[2]),",",_pos_fmt(text.style),
                    "){",quoted,"};")
        else
            println(io,"T3(",_pos_fmt(text.coords[1]),",",
                    _pos_fmt(text.coords[2]),",",_pos_fmt(text.coords[3]),",",
                    _pos_fmt(text.style),"){",quoted,"};")
        end
    end
    print(io,"};\n")
    return nothing
end

"""
    write_pos(io, views)
    write_pos(path, views)
    write_pos(io, view::PosView)
    write_pos(io, view::View)

Serialize `PosView` objects in Gmsh's parsed `.pos` text format — the exact
record layout `PViewDataList::writePOS` emits (`%.16g` numbers, one element
record per line, `TIME` only for multi-step views, records in stored order,
`INTERPOLATION_SCHEME` written back as read). Scalar [`View`](@ref) objects are
converted to `SP` records via [`PosView`](@ref). `parse`/`write` round-trips
are stable: `write_pos(read_pos(text))` is a fixed point of `write_pos ∘
read_pos`.
"""
function write_pos(io::IO,views::AbstractVector)
    for (i,view) in enumerate(views)
        view isa PosView || throw(ArgumentError(
            "write_pos: item $i is not a PosView (got $(typeof(view)))"))
        _pos_write_view(io,view)
    end
    return nothing
end
write_pos(io::IO,view::PosView)=write_pos(io,PosView[view])
write_pos(io::IO,view::View)=write_pos(io,PosView[PosView(view)])
function write_pos(path::AbstractString,views::AbstractVector)
    open(path,"w") do io
        write_pos(io,views)
    end
    return nothing
end
write_pos(path::AbstractString,view::PosView)=write_pos(path,PosView[view])
write_pos(path::AbstractString,view::View)=write_pos(path,PosView[PosView(view)])

"""
    postview_field(view::PosView; kwargs...) -> AbstractField

Build a [`PostViewField`](@ref) from a parsed `.pos` view's element records.
The dominant component kind follows Gmsh's `numComponents` precedence —
tensor if any record has 9 components, else vector if any has 3, else scalar —
and only records of that kind participate (matching `searchTensorClosest`
etc.). Tensor views return a [`PostViewAnisoField`](@ref) (usable with
`metric_at`, `size_at`, `refine_to_size`, and `build_geo_size_field`
contexts); scalar/vector views return a `PostViewField`. Keywords are
forwarded: `time`, `crop_negative`, `use_closest`, `reference_tolerance`,
`max_nodes`, `max_elements`.
"""
function postview_field(view::PosView;kwargs...)
    isempty(view.elements) && throw(ArgumentError(
        "postview_field: view has no element records"))
    ncomp=1
    for element in view.elements
        components=size(element.values,1)
        components==9 && (ncomp=9;break)
        components==3 && (ncomp=max(ncomp,3))
    end
    keep=ncomp
    records=[element for element in view.elements
             if size(element.values,1)==keep]
    field=PostViewField(records;kwargs...)
    return keep==9 ? PostViewAnisoField(field) : field
end
