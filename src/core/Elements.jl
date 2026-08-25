"""
    Elements

Gmsh 4.15.2 fixed-node element families and orders in memory, plus a
mixed-element mesh container with MSH v2.2/v4.1 round-trip. Simplex
[`Mesh`](@ref) remains the certified volume kernel. Cut/sub-element records
with ownership links and simplex-decomposed polygon/polyhedron connectivity
use a separate block type instead of pretending to be ordinary nodal elements.
"""
module Elements

using ..MeshTypes: Mesh, MeshDiagnostic, nnodes, nsegs, ntris, ntets
import ..MeshTypes: validate
using SHA
using Printf: @printf, @sprintf

export ElementSpec, MSH_CATALOG, msh_spec, msh_num_nodes, msh_dimension, msh_order, msh_family
export ElementBlock, ElementRef, SpecialElementBlock, MixedEntity, MixedEntityData, MixedMesh
export mixed_crc, simplex_to_mixed, mixed_to_simplex
export write_mixed_msh, read_mixed_msh
export lagrange_nodes, add_block!, validate

@inline function _elements_reject_bool(value,context::AbstractString)
    value isa Bool && throw(ArgumentError("$context must not be Bool"))
    return value
end

function _elements_reject_bool_values(values,context::AbstractString)
    (eltype(values)<:Bool || any(value->value isa Bool,values)) &&
        throw(ArgumentError("$context must not contain Bool values"))
    return values
end

"""One Gmsh `.msh` element type."""
struct ElementSpec
    msh::Int
    family::Symbol
    dim::Int
    order::Int
    nnodes::Int
    serendipity::Bool
end

# Authority: Gmsh 4.15.2 `GmshDefines.h`, `ElementType.cpp`, and
# `pointsGenerators.cpp`. Only element types with a fixed connectivity and a
# well-defined local-node layout are admitted here.
const _MSH_CATALOG_BUILD = Dict{Int,ElementSpec}()

function _expected_num_nodes(family::Symbol, order::Int, serendipity::Bool)
    order >= 0 || return -1
    family === :pnt && return order == 0 && !serendipity ? 1 : -1
    family === :lin && return !serendipity ? order + 1 : -1
    family === :tri && return serendipity ? (order == 0 ? 1 : 3order) :
        (order + 1) * (order + 2) ÷ 2
    family === :qua && return serendipity ? (order == 0 ? 1 : 4order) :
        (order + 1)^2
    family === :tet && return serendipity ? (order == 0 ? 1 : 4 + 6(order - 1)) :
        (order + 1) * (order + 2) * (order + 3) ÷ 6
    family === :hex && return serendipity ? (order == 0 ? 1 : 8 + 12(order - 1)) :
        (order + 1)^3
    family === :pri && return serendipity ? (order == 0 ? 1 : 6 + 9(order - 1)) :
        (order + 1)^2 * (order + 2) ÷ 2
    family === :pyr && return serendipity ? (order == 0 ? 1 : 5 + 8(order - 1)) :
        (order + 1) * (order + 2) * (2order + 3) ÷ 6
    family === :trih && return order == 1 && !serendipity ? 4 : -1
    return -1
end

function _reg(s::ElementSpec)
    1 <= s.msh <= 140 || throw(ErrorException("invalid MSH type $(s.msh)"))
    0 <= s.dim <= 3 || throw(ErrorException("invalid dimension $(s.dim) for MSH type $(s.msh)"))
    expected = _expected_num_nodes(s.family, s.order, s.serendipity)
    expected == s.nnodes || throw(ErrorException(
        "MSH type $(s.msh) has $(s.nnodes) nodes, expected $expected for " *
        "$(s.family) P$(s.order)"))
    haskey(_MSH_CATALOG_BUILD, s.msh) && throw(ErrorException("duplicate MSH type $(s.msh)"))
    _MSH_CATALOG_BUILD[s.msh] = s
    return s
end
for s in (
    ElementSpec(15,:pnt,0,0,1,false),
    ElementSpec(84,:lin,1,0,1,false),
    ElementSpec(1,:lin,1,1,2,false), ElementSpec(8,:lin,1,2,3,false),
    ElementSpec(26,:lin,1,3,4,false), ElementSpec(27,:lin,1,4,5,false),
    ElementSpec(28,:lin,1,5,6,false), ElementSpec(62,:lin,1,6,7,false),
    ElementSpec(63,:lin,1,7,8,false), ElementSpec(64,:lin,1,8,9,false),
    ElementSpec(65,:lin,1,9,10,false), ElementSpec(66,:lin,1,10,11,false),
    ElementSpec(85,:tri,2,0,1,false),
    ElementSpec(2,:tri,2,1,3,false), ElementSpec(9,:tri,2,2,6,false),
    ElementSpec(21,:tri,2,3,10,false), ElementSpec(23,:tri,2,4,15,false),
    ElementSpec(25,:tri,2,5,21,false), ElementSpec(42,:tri,2,6,28,false),
    ElementSpec(43,:tri,2,7,36,false), ElementSpec(44,:tri,2,8,45,false),
    ElementSpec(45,:tri,2,9,55,false), ElementSpec(46,:tri,2,10,66,false),
    ElementSpec(20,:tri,2,3,9,true), ElementSpec(22,:tri,2,4,12,true),
    ElementSpec(24,:tri,2,5,15,true), ElementSpec(52,:tri,2,6,18,true),
    ElementSpec(53,:tri,2,7,21,true), ElementSpec(54,:tri,2,8,24,true),
    ElementSpec(55,:tri,2,9,27,true), ElementSpec(56,:tri,2,10,30,true),
    ElementSpec(86,:qua,2,0,1,false),
    ElementSpec(3,:qua,2,1,4,false), ElementSpec(10,:qua,2,2,9,false),
    ElementSpec(36,:qua,2,3,16,false), ElementSpec(37,:qua,2,4,25,false),
    ElementSpec(38,:qua,2,5,36,false), ElementSpec(47,:qua,2,6,49,false),
    ElementSpec(48,:qua,2,7,64,false), ElementSpec(49,:qua,2,8,81,false),
    ElementSpec(50,:qua,2,9,100,false), ElementSpec(51,:qua,2,10,121,false),
    ElementSpec(16,:qua,2,2,8,true), ElementSpec(39,:qua,2,3,12,true),
    ElementSpec(40,:qua,2,4,16,true), ElementSpec(41,:qua,2,5,20,true),
    ElementSpec(57,:qua,2,6,24,true), ElementSpec(58,:qua,2,7,28,true),
    ElementSpec(59,:qua,2,8,32,true), ElementSpec(60,:qua,2,9,36,true),
    ElementSpec(61,:qua,2,10,40,true),
    ElementSpec(87,:tet,3,0,1,false),
    ElementSpec(4,:tet,3,1,4,false), ElementSpec(11,:tet,3,2,10,false),
    ElementSpec(29,:tet,3,3,20,false), ElementSpec(30,:tet,3,4,35,false),
    ElementSpec(31,:tet,3,5,56,false), ElementSpec(71,:tet,3,6,84,false),
    ElementSpec(72,:tet,3,7,120,false), ElementSpec(73,:tet,3,8,165,false),
    ElementSpec(74,:tet,3,9,220,false), ElementSpec(75,:tet,3,10,286,false),
    ElementSpec(32,:tet,3,4,22,true), ElementSpec(33,:tet,3,5,28,true),
    ElementSpec(79,:tet,3,6,34,true), ElementSpec(80,:tet,3,7,40,true),
    ElementSpec(81,:tet,3,8,46,true), ElementSpec(82,:tet,3,9,52,true),
    ElementSpec(83,:tet,3,10,58,true), ElementSpec(137,:tet,3,3,16,true),
    ElementSpec(88,:hex,3,0,1,false),
    ElementSpec(5,:hex,3,1,8,false), ElementSpec(12,:hex,3,2,27,false),
    ElementSpec(92,:hex,3,3,64,false), ElementSpec(93,:hex,3,4,125,false),
    ElementSpec(94,:hex,3,5,216,false), ElementSpec(95,:hex,3,6,343,false),
    ElementSpec(96,:hex,3,7,512,false), ElementSpec(97,:hex,3,8,729,false),
    ElementSpec(98,:hex,3,9,1000,false),
    ElementSpec(17,:hex,3,2,20,true), ElementSpec(99,:hex,3,3,32,true),
    ElementSpec(100,:hex,3,4,44,true), ElementSpec(101,:hex,3,5,56,true),
    ElementSpec(102,:hex,3,6,68,true), ElementSpec(103,:hex,3,7,80,true),
    ElementSpec(104,:hex,3,8,92,true), ElementSpec(105,:hex,3,9,104,true),
    ElementSpec(89,:pri,3,0,1,false),
    ElementSpec(6,:pri,3,1,6,false), ElementSpec(13,:pri,3,2,18,false),
    ElementSpec(18,:pri,3,2,15,true),
    ElementSpec(90,:pri,3,3,40,false), ElementSpec(91,:pri,3,4,75,false),
    ElementSpec(106,:pri,3,5,126,false), ElementSpec(107,:pri,3,6,196,false),
    ElementSpec(108,:pri,3,7,288,false), ElementSpec(109,:pri,3,8,405,false),
    ElementSpec(110,:pri,3,9,550,false),
    ElementSpec(111,:pri,3,3,24,true), ElementSpec(112,:pri,3,4,33,true),
    ElementSpec(113,:pri,3,5,42,true), ElementSpec(114,:pri,3,6,51,true),
    ElementSpec(115,:pri,3,7,60,true), ElementSpec(116,:pri,3,8,69,true),
    ElementSpec(117,:pri,3,9,78,true),
    ElementSpec(132,:pyr,3,0,1,false),
    ElementSpec(7,:pyr,3,1,5,false), ElementSpec(14,:pyr,3,2,14,false),
    ElementSpec(19,:pyr,3,2,13,true),
    ElementSpec(118,:pyr,3,3,30,false), ElementSpec(119,:pyr,3,4,55,false),
    ElementSpec(120,:pyr,3,5,91,false), ElementSpec(121,:pyr,3,6,140,false),
    ElementSpec(122,:pyr,3,7,204,false), ElementSpec(123,:pyr,3,8,285,false),
    ElementSpec(124,:pyr,3,9,385,false),
    ElementSpec(125,:pyr,3,3,21,true), ElementSpec(126,:pyr,3,4,29,true),
    ElementSpec(127,:pyr,3,5,37,true), ElementSpec(128,:pyr,3,6,45,true),
    ElementSpec(129,:pyr,3,7,53,true), ElementSpec(130,:pyr,3,8,61,true),
    ElementSpec(131,:pyr,3,9,69,true),
    ElementSpec(140,:trih,3,1,4,false),
)
    _reg(s)
end

"""
Fixed-node element specifications keyed by Gmsh 4.15.2 numeric type.

The catalog is immutable so callers cannot corrupt the process-wide element
registry. Individual [`ElementSpec`](@ref) records are immutable as well.
"""
const MSH_CATALOG = Base.ImmutableDict(
    sort!(collect(_MSH_CATALOG_BUILD);by=first)...)
const _MSH_CATALOG_LOOKUP = copy(_MSH_CATALOG_BUILD)
empty!(_MSH_CATALOG_BUILD)

# These are real Gmsh numeric tags, but not ordinary fixed-connectivity nodal
# elements. In Gmsh 4.15.2, MSH 34/35/69 records contain packed triangle/tet
# decompositions, 67/68/69 carry two domain-element links, and 34/35/70/133:136
# can carry a parent-element link. Tags 138/139 only select MINI bases:
# `MElement::getInfoMSH()` and `MElementFactory::create()` have no mesh-record
# cases for them. Keeping the classification here makes every rejection
# deliberate instead of conflating these tags with unknown/reserved IDs.
const MSH_SPECIAL_TYPES = Dict{Int,NamedTuple}(
    34 => (family=:polygon, dim=2, order=1, nnodes=nothing, kind=:decomposed),
    35 => (family=:polyhedron, dim=3, order=1, nnodes=nothing, kind=:decomposed),
    67 => (family=:line_border, dim=1, order=1, nnodes=2, kind=:border),
    68 => (family=:triangle_border, dim=2, order=1, nnodes=3, kind=:border),
    69 => (family=:polygon_border, dim=2, order=1, nnodes=nothing, kind=:border),
    70 => (family=:line_child, dim=1, order=1, nnodes=2, kind=:child),
    133 => (family=:point_xfem, dim=0, order=1, nnodes=1, kind=:subelement),
    134 => (family=:line_xfem, dim=1, order=1, nnodes=2, kind=:subelement),
    135 => (family=:triangle_xfem, dim=2, order=1, nnodes=3, kind=:subelement),
    136 => (family=:tetrahedron_xfem, dim=3, order=1, nnodes=4, kind=:subelement),
    138 => (family=:triangle_mini, dim=2, order=3, nnodes=4, kind=:basis_only),
    139 => (family=:tetrahedron_mini, dim=3, order=3, nnodes=5, kind=:basis_only),
)

# Record widths and link semantics are from `MElementCut.h`, `MSubElement.h`
# and the MSH2 reader/writer in pinned Gmsh 4.15.2. A `unit` greater than zero
# is both the fixed arity and, for decomposed records, the child-simplex width.
const MSH_SPECIAL_RECORDS = Dict{Int,NamedTuple}(
    34 => (unit=3, variable=true, links=:parent),
    35 => (unit=4, variable=true, links=:parent),
    67 => (unit=2, variable=false, links=:domains),
    68 => (unit=3, variable=false, links=:domains),
    69 => (unit=3, variable=true, links=:domains),
    70 => (unit=2, variable=false, links=:parent),
    133 => (unit=1, variable=false, links=:parent),
    134 => (unit=2, variable=false, links=:parent),
    135 => (unit=3, variable=false, links=:parent),
    136 => (unit=4, variable=false, links=:parent),
)

"""Return the fixed-node [`ElementSpec`](@ref) for a Gmsh numeric type."""
function msh_spec(msh::Integer)
    _elements_reject_bool(msh,"Elements: Gmsh element type")
    tag = try
        Int(msh)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("Elements: Gmsh element type is outside Int bounds"))
    end
    s=get(_MSH_CATALOG_LOOKUP, tag, nothing)
    if s === nothing
        haskey(MSH_SPECIAL_TYPES, tag) && throw(ArgumentError(
            "Elements: Gmsh element type $tag is special and is not an ordinary fixed-node element"))
        throw(ArgumentError("Elements: unknown Gmsh element type $msh"))
    end
    return s
end
"""Return the fixed number of local nodes for a Gmsh numeric type."""
msh_num_nodes(msh::Integer)=msh_spec(msh).nnodes
"""Return the topological dimension for a Gmsh numeric type."""
msh_dimension(msh::Integer)=msh_spec(msh).dim
"""Return the polynomial order for a Gmsh numeric type."""
msh_order(msh::Integer)=msh_spec(msh).order
"""Return the element-family symbol for a Gmsh numeric type."""
msh_family(msh::Integer)=msh_spec(msh).family

abstract type AbstractElementBlock end

"""One homogeneous block: `nodes` is `nnodes × ncells` (1-based indices)."""
struct ElementBlock <: AbstractElementBlock
    msh::Int
    nodes::Matrix{Int32}
    tags::Vector{Int32}
    function ElementBlock(msh::Integer, nodes::AbstractMatrix{<:Integer},
                          tags::AbstractVector{<:Integer}=zeros(Int32,size(nodes,2)))
        _elements_reject_bool(msh,"ElementBlock: Gmsh type")
        _elements_reject_bool_values(nodes,"ElementBlock: node indices")
        _elements_reject_bool_values(tags,"ElementBlock: physical tags")
        spec=msh_spec(msh)
        size(nodes,1)==spec.nnodes || throw(ArgumentError(
            "ElementBlock: type $msh expects $(spec.nnodes) nodes per cell, got $(size(nodes,1))"))
        length(tags)==size(nodes,2) || throw(ArgumentError("ElementBlock: tag length mismatch"))
        size(nodes,2) <= typemax(Int32) || throw(ArgumentError(
            "ElementBlock: cell count exceeds the Int32 MSH limit"))
        N = try
            Matrix{Int32}(nodes)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("ElementBlock: node indices must fit Int32: " *
                                sprint(showerror, err)))
        end
        T = try
            Vector{Int32}(tags)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("ElementBlock: physical tags must fit Int32: " *
                                sprint(showerror, err)))
        end
        @inbounds for j in axes(N,2), i in axes(N,1)
            N[i,j]>=1 || throw(ArgumentError("ElementBlock: node index must be positive"))
        end
        @inbounds for t in T
            t>=0 || throw(ArgumentError("ElementBlock: negative physical tag"))
        end
        return new(Int(msh),N,T)
    end
end

"""A 1-based reference to a cell in a [`MixedMesh`](@ref); `(0,0)` means absent."""
struct ElementRef
    block::Int32
    cell::Int32
    function ElementRef(block::Integer,cell::Integer)
        _elements_reject_bool(block,"ElementRef: block index")
        _elements_reject_bool(cell,"ElementRef: cell index")
        b=try Int(block) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("ElementRef: block index is outside Int bounds"))
        end
        c=try Int(cell) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("ElementRef: cell index is outside Int bounds"))
        end
        0<=b<=typemax(Int32) || throw(ArgumentError(
            "ElementRef: block index must be non-negative and fit Int32"))
        0<=c<=typemax(Int32) || throw(ArgumentError(
            "ElementRef: cell index must be non-negative and fit Int32"))
        (b==0)==(c==0) || throw(ArgumentError(
            "ElementRef: block and cell must either both be zero or both be positive"))
        return new(Int32(b),Int32(c))
    end
end
ElementRef()=ElementRef(0,0)

function _element_ref(value,context::AbstractString)
    value isa ElementRef && return value
    value isa Tuple && length(value)==2 || throw(ArgumentError(
        "$context must be an ElementRef or a (block, cell) tuple"))
    return try
        ElementRef(value[1],value[2])
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$context is invalid: "*sprint(showerror,err)))
    end
end

function _special_refs(values,ncells::Int,context::AbstractString)
    values===nothing && return fill(ElementRef(),ncells)
    values isa AbstractVector || throw(ArgumentError("$context must be a vector"))
    length(values)==ncells || throw(ArgumentError("$context length mismatch"))
    out=Vector{ElementRef}(undef,ncells)
    @inbounds for i in 1:ncells
        out[i]=_element_ref(values[i],"$context entry $i")
    end
    return out
end

function _special_domain_refs(values,ncells::Int,context::AbstractString)
    values===nothing && return fill(ElementRef(),2,ncells)
    values isa AbstractMatrix || throw(ArgumentError("$context must be a 2×ncells matrix"))
    size(values)==(2,ncells) || throw(ArgumentError("$context must be a 2×ncells matrix"))
    out=Matrix{ElementRef}(undef,2,ncells)
    @inbounds for j in 1:ncells, i in 1:2
        out[i,j]=_element_ref(values[i,j],"$context entry ($i,$j)")
    end
    return out
end

function _serializable_special_record(msh::Integer)
    _elements_reject_bool(msh,"SpecialElementBlock: Gmsh type")
    tag=try Int(msh) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("SpecialElementBlock: Gmsh type is outside Int bounds"))
    end
    record=get(MSH_SPECIAL_RECORDS,tag,nothing)
    record===nothing && throw(ArgumentError(
        haskey(MSH_SPECIAL_TYPES,tag) ?
        "SpecialElementBlock: Gmsh type $tag is not a serializable mesh record" :
        "SpecialElementBlock: unsupported Gmsh type $tag"))
    return tag,record
end

"""
    SpecialElementBlock(msh, cells, tags=zeros(...);
                        parent_refs=nothing, domain_refs=nothing)

A homogeneous block for Gmsh cut/sub-element records. `cells` contains one
connectivity vector per cell. Types 34, 35 and 69 store packed triangle or
tetrahedron decompositions, so their connectivity lengths must be positive
multiples of 3, 4 and 3 respectively. `parent_refs` and `domain_refs` retain
the MSH2 ownership links as mesh-cell references; an [`ElementRef`](@ref) of
`(0,0)` denotes a missing link.

MINI basis selectors 138/139 cannot be constructed as mesh records.
"""
struct SpecialElementBlock <: AbstractElementBlock
    msh::Int
    connectivity::Vector{Int32}
    offsets::Vector{Int32}
    tags::Vector{Int32}
    parent_refs::Vector{ElementRef}
    domain_refs::Matrix{ElementRef}
    function SpecialElementBlock(msh::Integer,
                                 connectivity::AbstractVector{<:Integer},
                                 offsets::AbstractVector{<:Integer},
                                 tags::AbstractVector{<:Integer};
                                 parent_refs=nothing,domain_refs=nothing)
        _elements_reject_bool_values(
            connectivity,"SpecialElementBlock: node indices")
        _elements_reject_bool_values(offsets,"SpecialElementBlock: offsets")
        _elements_reject_bool_values(tags,"SpecialElementBlock: physical tags")
        tag,record=_serializable_special_record(msh)
        length(tags)<=typemax(Int32) || throw(ArgumentError(
            "SpecialElementBlock: cell count exceeds the Int32 MSH limit"))
        length(connectivity)<typemax(Int32) || throw(ArgumentError(
            "SpecialElementBlock: connectivity count exceeds the Int32 CSR limit"))
        ncells=length(tags)
        length(offsets)==ncells+1 || throw(ArgumentError(
            "SpecialElementBlock: offset length mismatch"))
        offsets[1]==1 || throw(ArgumentError(
            "SpecialElementBlock: offsets must start at 1"))
        offsets[end]==length(connectivity)+1 || throw(ArgumentError(
            "SpecialElementBlock: final offset does not match connectivity length"))
        C=try Vector{Int32}(connectivity) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("SpecialElementBlock: node indices must fit Int32: "*
                                sprint(showerror,err)))
        end
        O=try Vector{Int32}(offsets) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("SpecialElementBlock: offsets must fit Int32: "*
                                sprint(showerror,err)))
        end
        T=try Vector{Int32}(tags) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("SpecialElementBlock: physical tags must fit Int32: "*
                                sprint(showerror,err)))
        end
        parents=_special_refs(parent_refs,ncells,"SpecialElementBlock: parent_refs")
        domains=_special_domain_refs(domain_refs,ncells,"SpecialElementBlock: domain_refs")
        @inbounds for j in 1:ncells
            O[j]<=O[j+1] || throw(ArgumentError(
                "SpecialElementBlock: offsets must be nondecreasing"))
            width=Int(O[j+1])-Int(O[j])
            if record.variable
                width>0 && width%record.unit==0 || throw(ArgumentError(
                    "SpecialElementBlock: type $tag cell $j connectivity must be a positive multiple of $(record.unit)"))
            else
                width==record.unit || throw(ArgumentError(
                    "SpecialElementBlock: type $tag cell $j expects $(record.unit) nodes, got $width"))
            end
            T[j]>=0 || throw(ArgumentError(
                "SpecialElementBlock: negative physical tag"))
            if record.links===:parent
                (domains[1,j]==ElementRef()&&domains[2,j]==ElementRef()) ||
                    throw(ArgumentError(
                        "SpecialElementBlock: type $tag does not carry domain_refs"))
            else
                parents[j]==ElementRef() || throw(ArgumentError(
                    "SpecialElementBlock: type $tag does not carry parent_refs"))
                domains[1,j]==ElementRef() && domains[2,j]!=ElementRef() &&
                    throw(ArgumentError(
                        "SpecialElementBlock: a second domain requires a first domain"))
            end
        end
        @inbounds for node in C
            node>=1 || throw(ArgumentError(
                "SpecialElementBlock: node index must be positive"))
        end
        return new(tag,C,O,T,parents,domains)
    end
end

function SpecialElementBlock(msh::Integer,cells::AbstractVector,
                             tags::AbstractVector{<:Integer}=zeros(Int32,length(cells));
                             parent_refs=nothing,domain_refs=nothing)
    tag,_=_serializable_special_record(msh)
    length(cells)<=typemax(Int32) || throw(ArgumentError(
        "SpecialElementBlock: cell count exceeds the Int32 MSH limit"))
    length(tags)==length(cells) || throw(ArgumentError(
        "SpecialElementBlock: tag length mismatch"))
    offsets=Vector{Int32}(undef,length(cells)+1); offsets[1]=1
    total=0
    for (j,cell) in pairs(cells)
        cell isa AbstractVector || throw(ArgumentError(
            "SpecialElementBlock: cell $j connectivity must be a vector"))
        total=try Base.checked_add(total,length(cell)) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("SpecialElementBlock: connectivity count overflows Int"))
        end
        total<typemax(Int32) || throw(ArgumentError(
            "SpecialElementBlock: connectivity count exceeds the Int32 CSR limit"))
        offsets[j+1]=Int32(total+1)
    end
    connectivity=Vector{Int32}(undef,total); position=1
    for (j,cell) in pairs(cells)
        @inbounds for value in cell
            value isa Integer || throw(ArgumentError(
                "SpecialElementBlock: cell $j node indices must be integers"))
            _elements_reject_bool(
                value,"SpecialElementBlock: cell $j node index")
            connectivity[position]=try Int32(value) catch err
                err isa InterruptException && rethrow()
                throw(ArgumentError(
                    "SpecialElementBlock: cell $j node index must fit Int32"))
            end
            position+=1
        end
    end
    return SpecialElementBlock(tag,connectivity,offsets,tags;
        parent_refs=parent_refs,domain_refs=domain_refs)
end

function SpecialElementBlock(msh::Integer,nodes::AbstractMatrix{<:Integer},
                             tags::AbstractVector{<:Integer}=zeros(Int32,size(nodes,2));
                             parent_refs=nothing,domain_refs=nothing)
    tag,_=_serializable_special_record(msh)
    length(tags)==size(nodes,2) || throw(ArgumentError(
        "SpecialElementBlock: tag length mismatch"))
    cells=[@view nodes[:,j] for j in axes(nodes,2)]
    return SpecialElementBlock(tag,cells,tags;
        parent_refs=parent_refs,domain_refs=domain_refs)
end

const MixedElementBlock=Union{ElementBlock,SpecialElementBlock}

@inline _block_ncells(block::ElementBlock)=size(block.nodes,2)
@inline _block_ncells(block::SpecialElementBlock)=length(block.tags)
@inline _block_dim(block::ElementBlock)=msh_spec(block.msh).dim
@inline _block_dim(block::SpecialElementBlock)=MSH_SPECIAL_TYPES[block.msh].dim
@inline _block_order(block::ElementBlock)=msh_spec(block.msh).order
@inline _block_order(block::SpecialElementBlock)=MSH_SPECIAL_TYPES[block.msh].order
@inline _cell_arity(block::ElementBlock,cell::Int)=size(block.nodes,1)
@inline _cell_arity(block::SpecialElementBlock,cell::Int)=
    Int(block.offsets[cell+1])-Int(block.offsets[cell])
@inline _cell_node(block::ElementBlock,cell::Int,slot::Int)=block.nodes[slot,cell]
@inline _cell_node(block::SpecialElementBlock,cell::Int,slot::Int)=
    block.connectivity[Int(block.offsets[cell])+slot-1]

@inline _copy_mixed_block(block::ElementBlock)=
    ElementBlock(block.msh,block.nodes,block.tags)
@inline _copy_mixed_block(block::SpecialElementBlock)=
    SpecialElementBlock(block.msh,block.connectivity,block.offsets,block.tags;
                        parent_refs=block.parent_refs,
                        domain_refs=block.domain_refs)

function _mixed_positive_int32(value,context::AbstractString)
    value isa Integer || throw(ArgumentError("$context must be an integer"))
    _elements_reject_bool(value,context)
    converted=try Int(value) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$context is outside Int bounds"))
    end
    1<=converted<=typemax(Int32) || throw(ArgumentError(
        "$context must be positive and fit Int32"))
    return Int32(converted)
end

function _mixed_int32_vector(values,context::AbstractString;
                             positive=false,nonzero=false,absolute_fits=false,
                             unique_values=false)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$context values must be a vector or tuple"))
    out=Int32[]; seen=Set{Int32}()
    for value in values
        value isa Integer || throw(ArgumentError("$context must be an integer"))
        _elements_reject_bool(value,context)
        converted=try Int(value) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$context is outside Int bounds"))
        end
        if absolute_fits
            -typemax(Int32)<=converted<=typemax(Int32) || throw(ArgumentError(
                "$context magnitude must fit Int32"))
        else
            typemin(Int32)<=converted<=typemax(Int32) || throw(ArgumentError(
                "$context must fit Int32"))
        end
        positive && converted<=0 && throw(ArgumentError("$context must be positive"))
        nonzero && converted==0 && throw(ArgumentError("$context cannot be zero"))
        item=Int32(converted)
        unique_values && item in seen && throw(ArgumentError(
            "$context values must be unique"))
        push!(seen,item); push!(out,item)
    end
    return out
end

function _mixed_positive_tag_vector(values,context::AbstractString;unique_values=false)
    (values isa AbstractVector || values isa Tuple) || throw(ArgumentError(
        "$context values must be a vector or tuple"))
    out=UInt64[]; seen=Set{UInt64}()
    for value in values
        value isa Integer || throw(ArgumentError("$context must be an integer"))
        _elements_reject_bool(value,context)
        value>0 || throw(ArgumentError("$context must be positive"))
        converted=try UInt64(value) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$context must fit an unsigned 64-bit MSH tag"))
        end
        unique_values && converted in seen && throw(ArgumentError(
            "$context values must be unique"))
        push!(seen,converted); push!(out,converted)
    end
    return out
end

"""One MSH v4 geometric entity, including all physical memberships and signed boundary tags."""
struct MixedEntity
    dim::Int
    tag::Int32
    bbox::NTuple{6,Float64}
    physical_tags::Vector{Int32}
    boundaries::Vector{Int32}
    function MixedEntity(dim::Integer,tag::Integer,bounds;
                         physical_tags=Int32[],boundaries=Int32[])
        _elements_reject_bool(dim,"MixedEntity: dimension")
        d=try Int(dim) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("MixedEntity: dimension is outside Int bounds"))
        end
        0<=d<=3 || throw(ArgumentError("MixedEntity: dimension $d is outside 0:3"))
        t=_mixed_positive_int32(tag,"MixedEntity: entity tag")
        (bounds isa Tuple || bounds isa AbstractVector) || throw(ArgumentError(
            "MixedEntity: bounds must be a tuple or vector"))
        raw=collect(bounds)
        _elements_reject_bool_values(raw,"MixedEntity: bounds")
        if d==0 && length(raw)==3
            raw=vcat(raw,raw)
        end
        length(raw)==6 || throw(ArgumentError(
            "MixedEntity: bounds must contain six values (or three for a point)"))
        converted=ntuple(6) do i
            value=try Float64(raw[i]) catch err
                err isa InterruptException && rethrow()
                throw(ArgumentError("MixedEntity: bounds must be Float64-representable"))
            end
            isfinite(value) || throw(ArgumentError("MixedEntity: non-finite bound"))
            value
        end
        all(converted[i]<=converted[i+3] for i in 1:3) || throw(ArgumentError(
            "MixedEntity: reversed bounding box"))
        if d==0
            all(isequal(converted[i],converted[i+3]) for i in 1:3) || throw(ArgumentError(
                "MixedEntity: point bounds must be degenerate"))
        end
        physical=_mixed_int32_vector(physical_tags,"MixedEntity: physical tag";
                                     positive=true)
        boundary=_mixed_int32_vector(boundaries,"MixedEntity: boundary tag";
                                     nonzero=true,absolute_fits=true)
        d==0 && !isempty(boundary) && throw(ArgumentError(
            "MixedEntity: point entities cannot have boundaries"))
        return new(d,t,converted,physical,boundary)
    end
end

# Suppress the autogenerated positional constructor: all public metadata
# construction must pass through the copying and validation path below.
struct _OwnedMixedEntityData end
const _OWNED_MIXED_ENTITY_DATA = _OwnedMixedEntityData()

"""
    MixedEntityData(entities; ...)

Optional lossless MSH v4 metadata. Node arrays align with mesh coordinates;
`block_entities` and `external_element_tags` align with every block and cell.
`nothing` in `node_parametric` denotes a non-parametric node block, while an
empty vector preserves a parametric dimension-zero block. External node and
element tags use the format's full unsigned 64-bit `size_t` range.
"""
struct MixedEntityData
    entities::Dict{Tuple{Int,Int},MixedEntity}
    node_entities::Vector{Tuple{Int,Int32}}
    node_parametric::Vector{Union{Nothing,Vector{Float64}}}
    external_node_tags::Vector{UInt64}
    block_entities::Vector{Vector{Int32}}
    external_element_tags::Vector{Vector{UInt64}}
    function MixedEntityData(::_OwnedMixedEntityData,
                             entities::Dict{Tuple{Int,Int},MixedEntity},
                             node_entities::Vector{Tuple{Int,Int32}},
                             node_parametric::Vector{Union{Nothing,Vector{Float64}}},
                             external_node_tags::Vector{UInt64},
                             block_entities::Vector{Vector{Int32}},
                             external_element_tags::Vector{Vector{UInt64}})
        return new(entities,node_entities,node_parametric,external_node_tags,
                   block_entities,external_element_tags)
    end
end

function MixedEntityData(entities::AbstractDict=Dict{Tuple{Int,Int},MixedEntity}();
                         node_entities=Tuple{Int,Int32}[],
                         node_parametric=Union{Nothing,Vector{Float64}}[],
                         external_node_tags=UInt64[],
                         block_entities=Vector{Int32}[],
                         external_element_tags=Vector{UInt64}[])
    copied_entities=Dict{Tuple{Int,Int},MixedEntity}()
    for (key,entity) in pairs(entities)
        entity isa MixedEntity || throw(ArgumentError(
            "MixedEntityData: entity values must be MixedEntity objects"))
        key isa Tuple && length(key)==2 || throw(ArgumentError(
            "MixedEntityData: entity keys must be (dimension, tag) tuples"))
        (key[1] isa Integer && key[2] isa Integer) || throw(ArgumentError(
            "MixedEntityData: entity keys must contain integers"))
        _elements_reject_bool(key[1],"MixedEntityData: entity-key dimension")
        _elements_reject_bool(key[2],"MixedEntityData: entity-key tag")
        key_dim=try Int(key[1]) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("MixedEntityData: entity-key dimension is outside Int bounds"))
        end
        key_tag=try Int(key[2]) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("MixedEntityData: entity-key tag is outside Int bounds"))
        end
        (key_dim,key_tag)==(entity.dim,Int(entity.tag)) || throw(ArgumentError(
            "MixedEntityData: entity key does not match its record"))
        copied_entities[(key_dim,key_tag)]=MixedEntity(
            entity.dim,entity.tag,entity.bbox;
            physical_tags=entity.physical_tags,boundaries=entity.boundaries)
    end
    copied_node_entities=Tuple{Int,Int32}[]
    for key in node_entities
        key isa Tuple && length(key)==2 || throw(ArgumentError(
            "MixedEntityData: node classification must be a (dimension, tag) tuple"))
        (key[1] isa Integer && key[2] isa Integer) || throw(ArgumentError(
            "MixedEntityData: node classification must contain integers"))
        _elements_reject_bool(key[1],"MixedEntityData: node dimension")
        dim=try Int(key[1]) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("MixedEntityData: node dimension is outside Int bounds"))
        end
        0<=dim<=3 || throw(ArgumentError(
            "MixedEntityData: node dimension $dim is outside 0:3"))
        tag=_mixed_positive_int32(key[2],"MixedEntityData: node entity tag")
        push!(copied_node_entities,(dim,tag))
    end
    copied_parametric=Vector{Union{Nothing,Vector{Float64}}}()
    for parameters in node_parametric
        if parameters===nothing
            push!(copied_parametric,nothing)
        else
            parameters isa AbstractVector{<:Real} || throw(ArgumentError(
                "MixedEntityData: node parameters must be vectors or nothing"))
            _elements_reject_bool_values(
                parameters,"MixedEntityData: node parameters")
            values=try Vector{Float64}(parameters) catch err
                err isa InterruptException && rethrow()
                throw(ArgumentError("MixedEntityData: node parameters must be Float64-representable"))
            end
            all(isfinite,values) || throw(ArgumentError(
                "MixedEntityData: non-finite node parameter"))
            push!(copied_parametric,values)
        end
    end
    node_tags=_mixed_positive_tag_vector(
        external_node_tags,"MixedEntityData: external node tag";unique_values=true)
    length(copied_node_entities)==length(copied_parametric)==length(node_tags) ||
        throw(ArgumentError(
            "MixedEntityData: node metadata arrays must have equal lengths"))
    for i in eachindex(copied_node_entities)
        parameters=copied_parametric[i]
        parameters===nothing || length(parameters)==copied_node_entities[i][1] ||
            throw(ArgumentError(
                "MixedEntityData: node $i parametric-coordinate length mismatch"))
    end
    copied_block_entities=Vector{Vector{Int32}}()
    for entities_in_block in block_entities
        push!(copied_block_entities,_mixed_int32_vector(
            entities_in_block,"MixedEntityData: cell entity tag";positive=true))
    end
    copied_element_tags=Vector{Vector{UInt64}}()
    seen_element_tags=Set{UInt64}()
    for tags_in_block in external_element_tags
        tags=_mixed_positive_tag_vector(
            tags_in_block,"MixedEntityData: external element tag")
        for tag in tags
            tag in seen_element_tags && throw(ArgumentError(
                "MixedEntityData: duplicate external element tag $tag"))
            push!(seen_element_tags,tag)
        end
        push!(copied_element_tags,tags)
    end
    length(copied_block_entities)==length(copied_element_tags) || throw(ArgumentError(
        "MixedEntityData: cell metadata block arrays must have equal lengths"))
    for i in eachindex(copied_block_entities)
        length(copied_block_entities[i])==length(copied_element_tags[i]) ||
            throw(ArgumentError(
                "MixedEntityData: cell metadata length mismatch in block $i"))
    end
    return MixedEntityData(
        _OWNED_MIXED_ENTITY_DATA,copied_entities,copied_node_entities,
        copied_parametric,node_tags,copied_block_entities,copied_element_tags)
end

# Internal readers may transfer freshly allocated storage without a second
# full copy; public construction always takes the detached path.
struct _OwnedMixedMesh end
const _OWNED_MIXED_MESH = _OwnedMixedMesh()

"""
Mixed-element mesh with owned coordinates, homogeneous MSH blocks, names, and
optional entity metadata. Construction detaches every mutable array supplied by
the caller.
"""
struct MixedMesh
    coords::Matrix{Float64}
    blocks::Vector{MixedElementBlock}
    physical_names::Dict{Tuple{Int,Int},String}
    entity_data::Union{Nothing,MixedEntityData}
    function MixedMesh(::_OwnedMixedMesh,C::Matrix{Float64},
                       B::Vector{MixedElementBlock},
                       names::Dict{Tuple{Int,Int},String},
                       data::Union{Nothing,MixedEntityData})
        mesh=new(C,B,names,data)
        _assert_mixed_structure(mesh,"MixedMesh")
        return mesh
    end
    function MixedMesh(coords::AbstractMatrix{<:Real}, blocks::AbstractVector;
                       physical_names=Dict{Tuple{Int,Int},String}(),entity_data=nothing)
        size(coords,1)==3 || throw(ArgumentError("MixedMesh: coords must be 3×n"))
        _elements_reject_bool_values(coords,"MixedMesh: coordinates")
        size(coords,2) <= typemax(Int32) || throw(ArgumentError(
            "MixedMesh: node count exceeds Int32 indexing"))
        C = try
            Matrix{Float64}(coords)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("MixedMesh: coordinates must be Float64-representable: " *
                                sprint(showerror, err)))
        end
        @inbounds for i in axes(C,2), d in 1:3
            isfinite(C[d,i]) || throw(ArgumentError("MixedMesh: non-finite coordinate"))
        end
        B=MixedElementBlock[]
        for (i,block) in pairs(blocks)
            block isa MixedElementBlock || throw(ArgumentError(
                "MixedMesh: block $i must be an ElementBlock or SpecialElementBlock"))
            push!(B,_copy_mixed_block(block))
        end
        names=_copy_physical_names(physical_names,"MixedMesh")
        data=entity_data===nothing ? nothing :
             entity_data isa MixedEntityData ? _copy_mixed_entity_data(entity_data) :
             throw(ArgumentError("MixedMesh: entity_data must be MixedEntityData or nothing"))
        return MixedMesh(_OWNED_MIXED_MESH,C,B,names,data)
    end
end

function _copy_physical_names(names, context::AbstractString)
    names isa AbstractDict || throw(ArgumentError(
        "$context: physical_names must be a dictionary"))
    out=Dict{Tuple{Int,Int},String}()
    for (key,value) in pairs(names)
        key isa Tuple && length(key)==2 || throw(ArgumentError(
            "$context: physical-name keys must be (dimension, tag) tuples"))
        dim,tag=key
        (dim isa Integer && tag isa Integer) || throw(ArgumentError(
            "$context: physical-name keys must contain integers"))
        _elements_reject_bool(dim,"$context: physical-name dimension")
        _elements_reject_bool(tag,"$context: physical-name tag")
        d = try Int(dim) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$context: physical-name dimension is outside Int bounds"))
        end
        t = try Int(tag) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$context: physical-name tag is outside Int bounds"))
        end
        0<=d<=3 || throw(ArgumentError(
            "$context: physical-name dimension $d is outside 0:3"))
        1<=t<=typemax(Int32) || throw(ArgumentError(
            "$context: physical-name tag $t must be positive and fit Int32"))
        value isa AbstractString || throw(ArgumentError(
            "$context: physical name for ($d,$t) must be a string"))
        s=String(value)
        isvalid(s) || throw(ArgumentError(
            "$context: physical name for ($d,$t) is not valid UTF-8"))
        occursin('\0',s) && throw(ArgumentError(
            "$context: physical name for ($d,$t) contains a NUL byte"))
        haskey(out,(d,t)) && throw(ArgumentError(
            "$context: duplicate physical name for ($d,$t)"))
        out[(d,t)]=s
    end
    return out
end

function _copy_mixed_entity_data(data::MixedEntityData)
    return MixedEntityData(data.entities;
        node_entities=data.node_entities,
        node_parametric=data.node_parametric,
        external_node_tags=data.external_node_tags,
        block_entities=data.block_entities,
        external_element_tags=data.external_element_tags)
end

function _assert_mixed_entity_data(m::MixedMesh,context::AbstractString)
    data=m.entity_data
    data===nothing && return nothing
    for (key,entity) in data.entities
        key==(entity.dim,Int(entity.tag)) || throw(ArgumentError(
            "$context: entity key $key does not match its record"))
        0<=entity.dim<=3 || throw(ArgumentError(
            "$context: entity $key has an invalid dimension"))
        entity.tag>0 || throw(ArgumentError(
            "$context: entity $key has a non-positive tag"))
        all(isfinite,entity.bbox) || throw(ArgumentError(
            "$context: entity $key has non-finite bounds"))
        all(entity.bbox[i]<=entity.bbox[i+3] for i in 1:3) || throw(ArgumentError(
            "$context: entity $key has reversed bounds"))
        if entity.dim==0
            isempty(entity.boundaries) || throw(ArgumentError(
                "$context: point entity $key has boundaries"))
            all(isequal(entity.bbox[i],entity.bbox[i+3]) for i in 1:3) || throw(ArgumentError(
                "$context: point entity $key has non-degenerate bounds"))
        end
        for physical in entity.physical_tags
            physical>0 || throw(ArgumentError(
                "$context: entity $key has a non-positive physical tag"))
        end
        for boundary in entity.boundaries
            boundary!=0 || throw(ArgumentError(
                "$context: entity $key has a zero boundary tag"))
            boundary!=typemin(Int32) || throw(ArgumentError(
                "$context: entity $key has a boundary magnitude outside Int32"))
            boundary_key=(entity.dim-1,abs(Int(boundary)))
            haskey(data.entities,boundary_key) || throw(ArgumentError(
                "$context: entity $key references missing boundary entity $boundary_key"))
        end
    end
    nn=size(m.coords,2)
    length(data.node_entities)==nn || throw(ArgumentError(
        "$context: node entity-classification length mismatch"))
    length(data.node_parametric)==nn || throw(ArgumentError(
        "$context: node parametric-data length mismatch"))
    length(data.external_node_tags)==nn || throw(ArgumentError(
        "$context: external node-tag length mismatch"))
    seen_nodes=Set{UInt64}()
    classified_entities=Set{Tuple{Int,Int}}()
    for i in 1:nn
        dim,tag=data.node_entities[i]
        0<=dim<=3 || throw(ArgumentError(
            "$context: node $i has invalid entity dimension $dim"))
        tag>0 || throw(ArgumentError(
            "$context: node $i has a non-positive entity tag"))
        push!(classified_entities,(dim,Int(tag)))
        parameters=data.node_parametric[i]
        if parameters!==nothing
            length(parameters)==dim || throw(ArgumentError(
                "$context: node $i has $(length(parameters)) parametric coordinates; expected $dim"))
            all(isfinite,parameters) || throw(ArgumentError(
                "$context: node $i has non-finite parametric coordinates"))
        end
        external=data.external_node_tags[i]
        external>0 || throw(ArgumentError(
            "$context: node $i has a non-positive external tag"))
        external in seen_nodes && throw(ArgumentError(
            "$context: duplicate external node tag $external"))
        push!(seen_nodes,external)
    end
    length(data.block_entities)==length(m.blocks) || throw(ArgumentError(
        "$context: block entity-classification length mismatch"))
    length(data.external_element_tags)==length(m.blocks) || throw(ArgumentError(
        "$context: external element-tag block count mismatch"))
    seen_elements=Set{UInt64}()
    for (bi,block) in pairs(m.blocks)
        entities=data.block_entities[bi]
        external_tags=data.external_element_tags[bi]
        ncells=_block_ncells(block)
        length(entities)==ncells || throw(ArgumentError(
            "$context: block $bi entity-classification length mismatch"))
        length(external_tags)==ncells || throw(ArgumentError(
            "$context: block $bi external element-tag length mismatch"))
        dim=_block_dim(block)
        for j in 1:ncells
            entity_tag=entities[j]
            entity_tag>0 || throw(ArgumentError(
                "$context: block $bi cell $j has a non-positive entity tag"))
            key=(dim,Int(entity_tag))
            entity=get(data.entities,key,nothing)
            entity===nothing && !(key in classified_entities) && throw(ArgumentError(
                "$context: block $bi cell $j references an entity not declared or created by nodes: $key"))
            projected=entity===nothing || isempty(entity.physical_tags) ? Int32(0) :
                      first(entity.physical_tags)
            block.tags[j]==projected || throw(ArgumentError(
                "$context: block $bi cell $j legacy physical tag $(block.tags[j]) " *
                "does not match entity $key projection $projected"))
            external=external_tags[j]
            external>0 || throw(ArgumentError(
                "$context: block $bi cell $j has a non-positive external element tag"))
            external in seen_elements && throw(ArgumentError(
                "$context: duplicate external element tag $external"))
            push!(seen_elements,external)
        end
    end
    return nothing
end

"""Append a nonempty, in-range element block to `m` and return `m`."""
function add_block!(m::MixedMesh, block::AbstractElementBlock)
    block isa MixedElementBlock || throw(ArgumentError(
        "add_block!: unsupported element block type"))
    m.entity_data===nothing || throw(ArgumentError(
        "add_block!: cannot append without synchronized v4 entity metadata"))
    owned=_copy_mixed_block(block)
    push!(m.blocks,owned)
    try
        _assert_mixed_structure(m,"add_block!")
    catch
        pop!(m.blocks)
        rethrow()
    end
    return m
end

"""Convert a simplex [`Mesh`](@ref) to linear Gmsh element blocks."""
function simplex_to_mixed(m::Mesh; physical_names=Dict{Tuple{Int,Int},String}())
    blocks=ElementBlock[]
    nsegs(m)>0 && push!(blocks,ElementBlock(1,m.segs,m.seg_tag))
    ntris(m)>0 && push!(blocks,ElementBlock(2,m.tris,m.tri_tag))
    ntets(m)>0 && push!(blocks,ElementBlock(4,m.tets,m.tet_tag))
    return MixedMesh(m.coords,blocks; physical_names=physical_names)
end

"""
    mixed_to_simplex(m) -> Mesh

Fold point/linear-line/linear-triangle/linear-tetrahedron blocks into a simplex
[`Mesh`](@ref). Any nonempty higher-order or non-simplex block is an explicit
blocker.
"""
function mixed_to_simplex(m::MixedMesh)
    _assert_mixed_structure(m,"mixed_to_simplex")
    ns=0; nf=0; nt=0
    for b in m.blocks
        isempty(b.tags) && continue
        b isa SpecialElementBlock && throw(ArgumentError(
            "mixed_to_simplex: cannot fold special type $(b.msh) into Mesh"))
        spec=msh_spec(b.msh)
        if spec.family===:lin && spec.order==1
            ns=Base.checked_add(ns,_block_ncells(b))
        elseif spec.family===:tri && spec.order==1
            nf=Base.checked_add(nf,_block_ncells(b))
        elseif spec.family===:tet && spec.order==1
            nt=Base.checked_add(nt,_block_ncells(b))
        elseif spec.family===:pnt
            continue
        else
            throw(ArgumentError(
                "mixed_to_simplex: cannot fold type $(b.msh) ($(spec.family) P$(spec.order)) into Mesh"))
        end
    end
    segs=Matrix{Int32}(undef,2,ns); st=Vector{Int32}(undef,ns)
    tris=Matrix{Int32}(undef,3,nf); tt=Vector{Int32}(undef,nf)
    tets=Matrix{Int32}(undef,4,nt); qt=Vector{Int32}(undef,nt)
    js=0; jf=0; jt=0
    @inbounds for b in m.blocks
        b isa SpecialElementBlock && continue
        spec=msh_spec(b.msh); n=_block_ncells(b); n==0 && continue
        if spec.family===:lin && spec.order==1
            copyto!(segs,2js+1,b.nodes,1,2n); copyto!(st,js+1,b.tags,1,n); js+=n
        elseif spec.family===:tri && spec.order==1
            copyto!(tris,3jf+1,b.nodes,1,3n); copyto!(tt,jf+1,b.tags,1,n); jf+=n
        elseif spec.family===:tet && spec.order==1
            copyto!(tets,4jt+1,b.nodes,1,4n); copyto!(qt,jt+1,b.tags,1,n); jt+=n
        end
    end
    return Mesh(m.coords; segs=segs, tris=tris, tets=tets, seg_tag=st, tri_tag=tt, tet_tag=qt)
end

struct _MixedCellRef
    block::Int
    cell::Int
end

@inline _missing_ref(ref::ElementRef)=ref.block==0

function _assert_element_ref(m::MixedMesh,ref::ElementRef,context::AbstractString)
    _missing_ref(ref) && return nothing
    block=Int(ref.block); cell=Int(ref.cell)
    1<=block<=length(m.blocks) || throw(ArgumentError(
        "$context references block $block outside 1:$(length(m.blocks))"))
    1<=cell<=_block_ncells(m.blocks[block]) || throw(ArgumentError(
        "$context references cell $cell outside block $block"))
    return nothing
end

@inline function _special_link(block::SpecialElementBlock,cell::Int,slot::Int)
    record=MSH_SPECIAL_RECORDS[block.msh]
    if record.links===:parent
        return slot==1 ? block.parent_refs[cell] : nothing
    end
    return slot==1 ? block.domain_refs[1,cell] :
           slot==2 ? block.domain_refs[2,cell] : nothing
end

function _assert_acyclic_element_links(m::MixedMesh,context::AbstractString)
    states=[zeros(UInt8,_block_ncells(block)) for block in m.blocks]
    stack_blocks=Int[]; stack_cells=Int[]; stack_slots=Int[]
    for (start_block,block) in pairs(m.blocks), start_cell in 1:_block_ncells(block)
        states[start_block][start_cell]==0 || continue
        push!(stack_blocks,start_block); push!(stack_cells,start_cell); push!(stack_slots,1)
        states[start_block][start_cell]=1
        while !isempty(stack_blocks)
            bi=stack_blocks[end]; cell=stack_cells[end]; slot=stack_slots[end]
            current=m.blocks[bi]
            ref=current isa SpecialElementBlock ? _special_link(current,cell,slot) : nothing
            if ref===nothing
                states[bi][cell]=2
                pop!(stack_blocks); pop!(stack_cells); pop!(stack_slots)
                continue
            end
            stack_slots[end]=slot+1
            _missing_ref(ref) && continue
            target_block=Int(ref.block); target_cell=Int(ref.cell)
            state=states[target_block][target_cell]
            state==1 && throw(ArgumentError(
                "$context: element parent/domain links contain a cycle"))
            state==2 && continue
            states[target_block][target_cell]=1
            push!(stack_blocks,target_block); push!(stack_cells,target_cell)
            push!(stack_slots,1)
        end
    end
    return nothing
end

function _assert_mixed_structure(m::MixedMesh, context::AbstractString)
    size(m.coords,1)==3 || throw(ArgumentError("$context: coords must be 3×n"))
    nn=size(m.coords,2)
    nn<=typemax(Int32) || throw(ArgumentError("$context: node count exceeds Int32 indexing"))
    @inbounds for i in axes(m.coords,2),d in 1:3
        isfinite(m.coords[d,i]) || throw(ArgumentError(
            "$context: node $i has a non-finite coordinate"))
    end
    total=0
    for (bi,b) in pairs(m.blocks)
        ncells=_block_ncells(b)
        ncells>0 || throw(ArgumentError("$context: block $bi is empty"))
        length(b.tags)==ncells || throw(ArgumentError(
            "$context: block $bi tag length mismatch"))
        if b isa ElementBlock
            spec=msh_spec(b.msh)
            size(b.nodes,1)==spec.nnodes || throw(ArgumentError(
                "$context: block $bi type $(b.msh) has the wrong node arity"))
        else
            record=get(MSH_SPECIAL_RECORDS,b.msh,nothing)
            record===nothing && throw(ArgumentError(
                "$context: block $bi type $(b.msh) is not a serializable special record"))
            length(b.offsets)==ncells+1 || throw(ArgumentError(
                "$context: block $bi offset length mismatch"))
            b.offsets[1]==1 && b.offsets[end]==length(b.connectivity)+1 ||
                throw(ArgumentError("$context: block $bi has invalid CSR endpoints"))
            length(b.parent_refs)==ncells || throw(ArgumentError(
                "$context: block $bi parent-reference length mismatch"))
            size(b.domain_refs)==(2,ncells) || throw(ArgumentError(
                "$context: block $bi domain-reference shape mismatch"))
            @inbounds for j in 1:ncells
                b.offsets[j]<=b.offsets[j+1] || throw(ArgumentError(
                    "$context: block $bi offsets are not nondecreasing"))
                width=_cell_arity(b,j)
                if record.variable
                    width>0 && width%record.unit==0 || throw(ArgumentError(
                        "$context: block $bi type $(b.msh) cell $j has invalid decomposed connectivity"))
                else
                    width==record.unit || throw(ArgumentError(
                        "$context: block $bi type $(b.msh) cell $j has the wrong node arity"))
                end
                if record.links===:parent
                    _missing_ref(b.domain_refs[1,j]) &&
                        _missing_ref(b.domain_refs[2,j]) || throw(ArgumentError(
                            "$context: block $bi type $(b.msh) cannot carry domain links"))
                else
                    _missing_ref(b.parent_refs[j]) || throw(ArgumentError(
                        "$context: block $bi type $(b.msh) cannot carry a parent link"))
                    _missing_ref(b.domain_refs[1,j]) &&
                        !_missing_ref(b.domain_refs[2,j]) && throw(ArgumentError(
                            "$context: block $bi cell $j has a second domain without a first"))
                end
            end
        end
        total=try Base.checked_add(total,ncells) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("$context: element count overflows Int"))
        end
        total<=typemax(Int32) || throw(ArgumentError(
            "$context: element count exceeds the Int32 MSH limit"))
        @inbounds for j in 1:ncells
            b.tags[j]>=0 || throw(ArgumentError(
                "$context: block $bi cell $j has negative physical tag $(b.tags[j])"))
            for i in 1:_cell_arity(b,j)
                n=_cell_node(b,j,i)
                1<=n<=nn || throw(ArgumentError(
                    "$context: block $bi cell $j references node $n outside 1:$nn"))
            end
            if b isa SpecialElementBlock
                _assert_element_ref(m,b.parent_refs[j],
                    "$context: block $bi cell $j parent")
                _assert_element_ref(m,b.domain_refs[1,j],
                    "$context: block $bi cell $j first domain")
                _assert_element_ref(m,b.domain_refs[2,j],
                    "$context: block $bi cell $j second domain")
            end
        end
    end
    _assert_acyclic_element_links(m,context)
    _copy_physical_names(m.physical_names,context)
    _assert_mixed_entity_data(m,context)
    return total
end

"""
    validate(m::MixedMesh; reject_duplicate_cells=true) -> MeshDiagnostic

Validate finite coordinates and connectivity. The check rejects repeated
indices within each ordinary cell or constituent simplex and, by default,
duplicate cells independent of their local orientation. Lossless v4 metadata is checked for aligned arrays,
valid entity references and topology, legacy-tag consistency, finite
parametrics, and globally unique positive external tags. Geometry/Jacobian
validity for curved high-order cells is not inferred from nodal coordinates by
this structural validator.
"""
function validate(m::MixedMesh; reject_duplicate_cells=true)
    reject_duplicate_cells isa Bool || throw(ArgumentError(
        "validate(MixedMesh): reject_duplicate_cells must be Bool"))
    messages=String[]
    try
        _assert_mixed_structure(m,"validate(MixedMesh)")
    catch err
        err isa InterruptException && rethrow()
        push!(messages,sprint(showerror,err))
        return MeshDiagnostic(false,messages)
    end
    seen=reject_duplicate_cells ? Set{Any}() : nothing
    for (bi,b) in pairs(m.blocks)
        @inbounds for j in 1:_block_ncells(b)
            k=_cell_arity(b,j)
            duplicate_part=false
            if b isa SpecialElementBlock && MSH_SPECIAL_RECORDS[b.msh].variable
                unit=MSH_SPECIAL_RECORDS[b.msh].unit
                parts=Any[]; sizehint!(parts,k÷unit); repeated=false
                for first_slot in 1:unit:k
                    scratch=Vector{Int32}(undef,unit)
                    for i in 1:unit
                        scratch[i]=_cell_node(b,j,first_slot+i-1)
                    end
                    sort!(scratch)
                    for i in 2:unit
                        scratch[i]==scratch[i-1] && (repeated=true)
                    end
                    push!(parts,Tuple(scratch))
                end
                sort!(parts)
                for i in 2:length(parts)
                    if parts[i]==parts[i-1]
                        duplicate_part=true; break
                    end
                end
                key=(b.msh,Tuple(parts))
            else
                scratch=Vector{Int32}(undef,k)
                for i in 1:k
                    scratch[i]=_cell_node(b,j,i)
                end
                sort!(scratch); repeated=false
                for i in 2:k
                    if scratch[i]==scratch[i-1]
                        repeated=true; break
                    end
                end
                key=(b.msh,Tuple(scratch))
            end
            repeated && push!(messages,
                "block $bi type $(b.msh) cell $j repeats a node index within a constituent simplex")
            duplicate_part && push!(messages,
                "block $bi type $(b.msh) cell $j duplicates a constituent simplex")
            if reject_duplicate_cells
                if key in seen
                    push!(messages,
                        "duplicate type $(b.msh) cell at block $bi cell $j")
                else
                    push!(seen,key)
                end
            end
        end
    end
    return MeshDiagnostic(isempty(messages),messages)
end

@inline function _mixed_cell_lt(m::MixedMesh,a::_MixedCellRef,b::_MixedCellRef)
    ba=m.blocks[a.block]; bb=m.blocks[b.block]
    ba.msh!=bb.msh && return ba.msh<bb.msh
    na=_cell_arity(ba,a.cell); nb=_cell_arity(bb,b.cell)
    na!=nb && return na<nb
    @inbounds for i in 1:na
        va=_cell_node(ba,a.cell,i); vb=_cell_node(bb,b.cell,i)
        va!=vb && return va<vb
    end
    ta=ba.tags[a.cell]; tb=bb.tags[b.cell]
    ta!=tb && return ta<tb
    data=m.entity_data
    if data!==nothing
        ea=data.block_entities[a.block][a.cell]
        eb=data.block_entities[b.block][b.cell]
        ea!=eb && return ea<eb
        xa=data.external_element_tags[a.block][a.cell]
        xb=data.external_element_tags[b.block][b.cell]
        xa!=xb && return xa<xb
    end
    return false
end

"""
    mixed_crc(m) -> NamedTuple

Return a deterministic SHA-256 regression record. The digest includes every
coordinate bit pattern, Gmsh type/dimension/order, oriented local connectivity,
special parent/domain link, physical tag, and physical name. When present, lossless v4 entity records,
classification, parametric coordinates, and external tags are also included.
Cell and block iteration order do not affect the digest; reversing a cell does.
Empty meshes have a finite zero bounding box.
"""
function mixed_crc(m::MixedMesh)
    ncells=_assert_mixed_structure(m,"mixed_crc")
    nn=size(m.coords,2)
    if nn==0
        lo=(0.0,0.0,0.0); hi=lo
    else
        l1=l2=l3=Inf; h1=h2=h3=-Inf
        @inbounds for i in axes(m.coords,2)
            x=m.coords[1,i]; y=m.coords[2,i]; z=m.coords[3,i]
            l1=min(l1,x); l2=min(l2,y); l3=min(l3,z)
            h1=max(h1,x); h2=max(h2,y); h3=max(h3,z)
        end
        lo=(l1,l2,l3); hi=(h1,h2,h3)
    end
    refs=Vector{_MixedCellRef}(undef,ncells); q=0
    @inbounds for (bi,b) in pairs(m.blocks),j in 1:_block_ncells(b)
        q+=1; refs[q]=_MixedCellRef(bi,j)
    end
    sort!(refs;lt=(a,b)->_mixed_cell_lt(m,a,b),alg=MergeSort)
    canonical=Dict{Tuple{Int,Int},Int}()
    sizehint!(canonical,ncells)
    @inbounds for (rank,ref) in pairs(refs)
        canonical[(ref.block,ref.cell)]=rank
    end
    ctx=SHA.SHA2_256_CTX()
    data=m.entity_data
    has_special=any(block->block isa SpecialElementBlock,m.blocks)
    prefix=has_special ? (data===nothing ? "Tessella.MixedMesh.CRC.v3\0" :
                                           "Tessella.MixedMesh.CRC.v4\0") :
                         (data===nothing ? "Tessella.MixedMesh.CRC.v1\0" :
                                           "Tessella.MixedMesh.CRC.v2\0")
    SHA.update!(ctx,codeunits(prefix))
    buf8=Vector{UInt8}(undef,8); buf4=Vector{UInt8}(undef,4)
    _sha_u64!(ctx,buf8,UInt64(nn))
    @inbounds for i in axes(m.coords,2),d in 1:3
        _sha_u64!(ctx,buf8,reinterpret(UInt64,m.coords[d,i]))
    end
    names=sort!(collect(m.physical_names);by=first)
    _sha_u64!(ctx,buf8,UInt64(length(names)))
    for ((dim,tag),name) in names
        _sha_i32!(ctx,buf4,Int32(dim)); _sha_i32!(ctx,buf4,Int32(tag))
        bytes=codeunits(name); _sha_u64!(ctx,buf8,UInt64(length(bytes)))
        SHA.update!(ctx,bytes)
    end
    if data!==nothing
        entities=sort!(collect(values(data.entities));by=e->(e.dim,e.tag))
        _sha_u64!(ctx,buf8,UInt64(length(entities)))
        for entity in entities
            _sha_i32!(ctx,buf4,Int32(entity.dim))
            _sha_i32!(ctx,buf4,entity.tag)
            for value in entity.bbox
                _sha_u64!(ctx,buf8,reinterpret(UInt64,value))
            end
            _sha_u64!(ctx,buf8,UInt64(length(entity.physical_tags)))
            for physical in entity.physical_tags
                _sha_i32!(ctx,buf4,physical)
            end
            _sha_u64!(ctx,buf8,UInt64(length(entity.boundaries)))
            for boundary in entity.boundaries
                _sha_i32!(ctx,buf4,boundary)
            end
        end
        @inbounds for i in 1:nn
            dim,entity=data.node_entities[i]
            _sha_u64!(ctx,buf8,UInt64(data.external_node_tags[i]))
            _sha_i32!(ctx,buf4,Int32(dim)); _sha_i32!(ctx,buf4,entity)
            parameters=data.node_parametric[i]
            _sha_i32!(ctx,buf4,Int32(parameters===nothing ? 0 : 1))
            if parameters!==nothing
                _sha_u64!(ctx,buf8,UInt64(length(parameters)))
                for value in parameters
                    _sha_u64!(ctx,buf8,reinterpret(UInt64,value))
                end
            end
        end
    end
    _sha_u64!(ctx,buf8,UInt64(ncells))
    @inbounds for ref in refs
        b=m.blocks[ref.block]
        _sha_i32!(ctx,buf4,Int32(b.msh)); _sha_i32!(ctx,buf4,Int32(_block_dim(b)))
        _sha_i32!(ctx,buf4,Int32(_block_order(b)))
        _sha_i32!(ctx,buf4,Int32(b isa ElementBlock ?
            msh_spec(b.msh).serendipity : -1))
        _sha_i32!(ctx,buf4,b.tags[ref.cell])
        arity=_cell_arity(b,ref.cell)
        _sha_i32!(ctx,buf4,Int32(arity))
        for i in 1:arity
            _sha_i32!(ctx,buf4,_cell_node(b,ref.cell,i))
        end
        if b isa SpecialElementBlock
            record=MSH_SPECIAL_RECORDS[b.msh]
            _sha_i32!(ctx,buf4,Int32(record.variable ? 1 : 0))
            _sha_i32!(ctx,buf4,Int32(record.unit))
            _sha_i32!(ctx,buf4,Int32(record.links===:parent ? 1 : 2))
            parent=b.parent_refs[ref.cell]
            first_domain=b.domain_refs[1,ref.cell]
            second_domain=b.domain_refs[2,ref.cell]
            for link in (parent,first_domain,second_domain)
                rank=_missing_ref(link) ? 0 : canonical[(Int(link.block),Int(link.cell))]
                _sha_u64!(ctx,buf8,UInt64(rank))
            end
        end
        if data!==nothing
            _sha_i32!(ctx,buf4,data.block_entities[ref.block][ref.cell])
            _sha_u64!(ctx,buf8,
                UInt64(data.external_element_tags[ref.block][ref.cell]))
        end
    end
    return (n_nodes=nn,n_blocks=length(m.blocks),n_cells=ncells,
            bbox=(lo,hi),sha=bytes2hex(SHA.digest!(ctx)))
end

@inline function _put_i32!(buf,off,v::Int32)
    u=reinterpret(UInt32,v)
    @inbounds buf[off]=u%UInt8
    @inbounds buf[off+1]=(u>>8)%UInt8
    @inbounds buf[off+2]=(u>>16)%UInt8
    @inbounds buf[off+3]=(u>>24)%UInt8
    return nothing
end
@inline function _put_u64!(buf,off,u::UInt64)
    @inbounds for shift in 0:8:56
        buf[off+(shift>>3)]=(u>>shift)%UInt8
    end
    return nothing
end
@inline function _sha_i32!(ctx,buf,v::Int32)
    _put_i32!(buf,1,v); SHA.update!(ctx,buf); return nothing
end
@inline function _sha_u64!(ctx,buf,v::UInt64)
    _put_u64!(buf,1,v); SHA.update!(ctx,buf); return nothing
end

# ── Gmsh local nodes ────────────────────────────────────────────────────────

# Gmsh's local-node order is topological: vertices, oriented edges, oriented
# face interiors, then volume interiors. These integral monomial lattices are a
# direct translation of `numeric/pointsGenerators.cpp`; scaling happens once at
# the end to avoid accumulated floating-point drift.
const _TRI_EDGES = ((0,1), (1,2), (2,0))
const _QUA_EDGES = ((0,1), (1,2), (2,3), (3,0))
const _TET_EDGES = ((0,1), (1,2), (2,0), (3,0), (3,2), (3,1))
const _TET_FACES = ((0,2,1), (0,1,3), (0,3,2), (3,1,2))
const _PRI_EDGES = ((0,1), (0,2), (0,3), (1,2), (1,4), (2,5),
                    (3,4), (3,5), (4,5))
const _PRI_FACES = ((0,2,1,-1), (3,4,5,-1), (0,1,4,3),
                    (0,3,5,2), (1,2,5,4))
const _HEX_EDGES = ((0,1), (0,3), (0,4), (1,2), (1,5), (2,3),
                    (2,6), (3,7), (4,5), (4,7), (5,6), (6,7))
const _HEX_FACES = ((0,3,2,1), (0,1,5,4), (0,4,7,3),
                    (1,2,6,5), (2,3,7,6), (4,5,6,7))
const _PYR_EDGES = ((0,1), (0,3), (0,4), (1,2),
                    (1,4), (2,3), (2,4), (3,4))
const _PYR_FACES = ((0,1,4,-1), (3,0,4,-1), (1,2,4,-1),
                    (2,3,4,-1), (0,3,2,1))

function _append_edges!(points::Vector{NTuple{D,Int}},
                        vertices::Vector{NTuple{D,Int}}, edges, p::Int) where {D}
    p > 1 || return points
    for (a, b) in edges
        v0, v1 = vertices[a + 1], vertices[b + 1]
        step = ntuple(d -> (v1[d] - v0[d]) ÷ p, D)
        for i in 1:p-1
            push!(points, ntuple(d -> v0[d] + i * step[d], D))
        end
    end
    return points
end

function _append_face!(points::Vector{NTuple{D,Int}},
                       vertices::Vector{NTuple{D,Int}},
                       face, local_points, p::Int) where {D}
    i0, i1 = face[1] + 1, face[2] + 1
    # For quadrilateral faces Gmsh uses vertex 3 as the second axis.
    i2 = (length(face) == 4 && face[4] != -1 ? face[4] : face[3]) + 1
    v0, v1, v2 = vertices[i0], vertices[i1], vertices[i2]
    u = ntuple(d -> (v1[d] - v0[d]) ÷ p, D)
    v = ntuple(d -> (v2[d] - v0[d]) ÷ p, D)
    for q in local_points
        push!(points, ntuple(d -> v0[d] + u[d] * q[1] + v[d] * q[2], D))
    end
    return points
end

function _line_monomials(p::Int)
    points = NTuple{1,Int}[(0,)]
    p == 0 && return points
    push!(points, (p,))
    append!(points, ((i,) for i in 1:p-1))
    return points
end

function _tri_monomials(p::Int, serendipity::Bool=false)
    p == 0 && return NTuple{2,Int}[(0,0)]
    vertices = NTuple{2,Int}[(0,0), (p,0), (0,p)]
    points = copy(vertices)
    _append_edges!(points, vertices, _TRI_EDGES, p)
    if !serendipity && p > 2
        append!(points, ((q[1] + 1, q[2] + 1) for q in _tri_monomials(p - 3)))
    end
    return points
end

function _qua_monomials(p::Int, serendipity::Bool=false)
    p == 0 && return NTuple{2,Int}[(0,0)]
    vertices = NTuple{2,Int}[(0,0), (p,0), (p,p), (0,p)]
    points = copy(vertices)
    _append_edges!(points, vertices, _QUA_EDGES, p)
    if !serendipity && p > 1
        append!(points, ((q[1] + 1, q[2] + 1) for q in _qua_monomials(p - 2)))
    end
    return points
end

function _tet_monomials(p::Int, serendipity::Bool=false)
    p == 0 && return NTuple{3,Int}[(0,0,0)]
    vertices = NTuple{3,Int}[(0,0,0), (p,0,0), (0,p,0), (0,0,p)]
    points = copy(vertices)
    _append_edges!(points, vertices, _TET_EDGES, p)
    if !serendipity && p > 2
        face_points = [(q[1] + 1, q[2] + 1) for q in _tri_monomials(p - 3)]
        for face in _TET_FACES
            _append_face!(points, vertices, face, face_points, p)
        end
        if p > 3
            append!(points, ((q[1] + 1, q[2] + 1, q[3] + 1)
                             for q in _tet_monomials(p - 4)))
        end
    end
    return points
end

function _pri_monomials(p::Int, serendipity::Bool=false)
    p == 0 && return NTuple{3,Int}[(0,0,0)]
    vertices = NTuple{3,Int}[(0,0,0), (p,0,0), (0,p,0),
                             (0,0,p), (p,0,p), (0,p,p)]
    points = copy(vertices)
    _append_edges!(points, vertices, _PRI_EDGES, p)
    if !serendipity && p > 1
        quad_points = [(q[1] + 1, q[2] + 1) for q in _qua_monomials(p - 2)]
        tri_points = p > 2 ?
            [(q[1] + 1, q[2] + 1) for q in _tri_monomials(p - 3)] :
            NTuple{2,Int}[]
        for face in _PRI_FACES
            local_points = face[4] == -1 ? tri_points : quad_points
            isempty(local_points) || _append_face!(points, vertices, face, local_points, p)
        end
        if p > 2
            for q in _tri_monomials(p - 3), r in _line_monomials(p - 2)
                push!(points, (q[1] + 1, q[2] + 1, r[1] + 1))
            end
        end
    end
    return points
end

function _hex_monomials(p::Int, serendipity::Bool=false)
    p == 0 && return NTuple{3,Int}[(0,0,0)]
    vertices = NTuple{3,Int}[(0,0,0), (p,0,0), (p,p,0), (0,p,0),
                             (0,0,p), (p,0,p), (p,p,p), (0,p,p)]
    points = copy(vertices)
    _append_edges!(points, vertices, _HEX_EDGES, p)
    if !serendipity && p > 1
        face_points = [(q[1] + 1, q[2] + 1) for q in _qua_monomials(p - 2)]
        for face in _HEX_FACES
            _append_face!(points, vertices, face, face_points, p)
        end
        append!(points, ((q[1] + 1, q[2] + 1, q[3] + 1)
                         for q in _hex_monomials(p - 2)))
    end
    return points
end

function _pyr_monomials(p::Int, serendipity::Bool=false)
    p == 0 && return NTuple{3,Int}[(0,0,0)]
    vertices = NTuple{3,Int}[(0,0,p), (p,0,p), (p,p,p), (0,p,p), (0,0,0)]
    points = copy(vertices)
    _append_edges!(points, vertices, _PYR_EDGES, p)
    if !serendipity && p > 1
        quad_points = [(q[1] + 1, q[2] + 1) for q in _qua_monomials(p - 2)]
        tri_points = p > 2 ?
            [(q[1] + 1, q[2] + 1) for q in _tri_monomials(p - 3)] :
            NTuple{2,Int}[]
        for face in _PYR_FACES
            local_points = face[4] == -1 ? tri_points : quad_points
            isempty(local_points) || _append_face!(points, vertices, face, local_points, p)
        end
        if p > 2
            append!(points, ((q[1] + 1, q[2] + 1, q[3] + 2)
                             for q in _pyr_monomials(p - 3)))
        end
    end
    return points
end

const _MSH_BY_SHAPE = let by_shape = Dict{Tuple{Symbol,Int,Bool},Int}()
    for (tag, spec) in MSH_CATALOG
        key = (spec.family, spec.order, spec.serendipity)
        haskey(by_shape, key) && throw(ErrorException("duplicate Gmsh element shape $key"))
        by_shape[key] = tag
    end
    by_shape
end

function _scaled_nodes(spec::ElementSpec)
    p, family = spec.order, spec.family
    family === :trih && return Float64[-1 1 1 -1; -1 -1 1 1; 0 0 0 0]

    monomials = if family === :pnt
        NTuple{1,Int}[(0,)]
    elseif family === :lin
        _line_monomials(p)
    elseif family === :tri
        _tri_monomials(p, spec.serendipity)
    elseif family === :qua
        _qua_monomials(p, spec.serendipity)
    elseif family === :tet
        _tet_monomials(p, spec.serendipity)
    elseif family === :hex
        _hex_monomials(p, spec.serendipity)
    elseif family === :pri
        _pri_monomials(p, spec.serendipity)
    elseif family === :pyr
        _pyr_monomials(p, spec.serendipity)
    else
        throw(ErrorException("unsupported catalog family $family"))
    end

    length(monomials) == spec.nnodes || throw(ErrorException(
        "local-node generator produced $(length(monomials)) nodes for Gmsh type " *
        "$(spec.msh), expected $(spec.nnodes)"))
    coordinates = zeros(Float64, 3, length(monomials))
    p == 0 && return coordinates
    @inbounds for (column, q) in enumerate(monomials)
        if family === :lin
            coordinates[1,column] = 2q[1] / p - 1
        elseif family === :tri
            coordinates[1,column] = q[1] / p
            coordinates[2,column] = q[2] / p
        elseif family === :qua
            coordinates[1,column] = 2q[1] / p - 1
            coordinates[2,column] = 2q[2] / p - 1
        elseif family === :tet
            coordinates[1,column] = q[1] / p
            coordinates[2,column] = q[2] / p
            coordinates[3,column] = q[3] / p
        elseif family === :hex
            coordinates[1,column] = 2q[1] / p - 1
            coordinates[2,column] = 2q[2] / p - 1
            coordinates[3,column] = 2q[3] / p - 1
        elseif family === :pri
            coordinates[1,column] = q[1] / p
            coordinates[2,column] = q[2] / p
            coordinates[3,column] = 2q[3] / p - 1
        elseif family === :pyr
            coordinates[3,column] = 1 - q[3] / p
            coordinates[1,column] = (2q[1] - q[3]) / p
            coordinates[2,column] = (2q[2] - q[3]) / p
        end
    end
    return coordinates
end

"""Return Gmsh-ordered local coordinates as a `3 × nnodes` matrix."""
lagrange_nodes(msh::Integer) = _scaled_nodes(msh_spec(msh))

"""
    lagrange_nodes(family, order; serendipity=false)

Return local nodes for the unique catalogued Gmsh type with this shape. The
lookup bounds accepted orders to those defined by Gmsh 4.15.2.
"""
function lagrange_nodes(family::Symbol, order::Integer; serendipity=false)
    _elements_reject_bool(order,"lagrange_nodes: order")
    serendipity isa Bool || throw(ArgumentError(
        "lagrange_nodes: serendipity must be Bool"))
    p = try
        Int(order)
    catch err
        err isa InexactError || rethrow()
        throw(ArgumentError("lagrange_nodes: order is outside Int bounds"))
    end
    p >= 0 || throw(ArgumentError("lagrange_nodes: order must be non-negative"))
    tag = get(_MSH_BY_SHAPE, (family, p, serendipity), nothing)
    tag === nothing && throw(ArgumentError(
        "lagrange_nodes: no Gmsh 4.15.2 element for family $family, order $p, " *
        "serendipity=$serendipity"))
    return lagrange_nodes(tag)
end

# ── MSH I/O for MixedMesh ─────────────────────────────────────────────────────

# Binary layout authority: Gmsh tag `gmsh_4_15_2`,
# `doc/texinfo/gmsh.texi`, `src/geo/GModelIO_MSH2.cpp`,
# `src/geo/GModelIO_MSH4.cpp`, `src/geo/MElementCut.h` and
# `src/geo/MSubElement.h`. MSH2 uses interleaved Int32/Float64 records; MSH4
# separates Int32 headers, UInt64 `size_t` values and Float64 payloads.
const DEFAULT_MAX_MIXED_NAME_BYTES = 1 << 20
const MSH_PHYSICAL_NAME_MAX_BYTES = 128

# `MElementFactory::create()` in pinned Gmsh 4.15.2 has no constructor for most
# of these otherwise-catalogued fixed-node tags. Type 69 parses, but Gmsh 4.15.2
# crashes while destroying the resulting polygon-border element. Keep Tessella
# round-trip support, but make unsafe output an explicit opt-in instead of
# silently producing a file Gmsh cannot consume through its normal lifecycle.
const GMSH_4_15_2_MSH_READER_GAPS_V4 = (
    84,85,86,87,88,             # P0 line/triangle/quad/tet/hex
    100,101,102,103,104,105,    # incomplete hex P4:P9
    125,126,127,128,129,130,131,# incomplete pyramid P3:P9
    132,                         # P0 pyramid
)
const GMSH_4_15_2_MSH_READER_GAPS_V2 = (
    GMSH_4_15_2_MSH_READER_GAPS_V4...,
    69,89,140) # v2 also crashes on polygon-border/P0 prism and rejects trihedron

struct _MixedReadLimits
    max_nodes::Int
    max_elements::Int
    max_connectivity::Int
    max_blocks::Int
    max_entities::Int
    max_physical_names::Int
    max_name_bytes::Int
    max_file_bytes::Int
end

mutable struct _MixedReadBucket
    nodes::Vector{Int32}
    tags::Vector{Int32}
    entities::Vector{Int32}
    external_tags::Vector{UInt64}
    offsets::Vector{Int32}
    parent_tags::Vector{UInt64}
    domain_tags::Vector{NTuple{2,UInt64}}
end

function _MixedReadBucket(etype::Int)
    record=get(MSH_SPECIAL_RECORDS,etype,nothing)
    offsets=record!==nothing && record.variable ? Int32[1] : Int32[]
    return _MixedReadBucket(Int32[],Int32[],Int32[],UInt64[],offsets,
                            UInt64[],NTuple{2,UInt64}[])
end

mutable struct _MixedReadAccum
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    node_map::Dict{UInt64,Int32}
    buckets::Dict{Int,_MixedReadBucket}
    physical_names::Dict{Tuple{Int,Int},String}
    entity_physical::Dict{Tuple{Int,Int},Int32}
    entities::Dict{Tuple{Int,Int},MixedEntity}
    implicit_entities::Set{Tuple{Int,Int}}
    external_node_tags::Vector{UInt64}
    node_entities::Vector{Tuple{Int,Int32}}
    node_parametric::Vector{Union{Nothing,Vector{Float64}}}
    element_tags::Set{UInt64}
    connectivity_entries::Int
    physical_name_records::Int
    node_blocks::Int
    element_blocks::Int
end

_MixedReadAccum() = _MixedReadAccum(
    Float64[],Float64[],Float64[],Dict{UInt64,Int32}(),
    Dict{Int,_MixedReadBucket}(),Dict{Tuple{Int,Int},String}(),
    Dict{Tuple{Int,Int},Int32}(),Dict{Tuple{Int,Int},MixedEntity}(),
    Set{Tuple{Int,Int}}(),UInt64[],Tuple{Int,Int32}[],
    Union{Nothing,Vector{Float64}}[],Set{UInt64}(),0,0,0,0)

function _read_limit(value,name::AbstractString;ceiling=typemax(Int))
    value isa Integer || throw(ArgumentError("read_mixed_msh: $name must be an integer"))
    _elements_reject_bool(value,"read_mixed_msh: $name")
    ceiling isa Integer || throw(ArgumentError("read_mixed_msh: invalid internal limit ceiling"))
    upper=Int(ceiling)
    out = try
        Int(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: $name is outside Int bounds"))
    end
    0<=out<=upper || throw(ArgumentError(
        "read_mixed_msh: $name must be in 0:$upper"))
    return out
end

"""
    read_mixed_msh(path; tessella_extensions=false, resource limits...) -> MixedMesh

Read an ASCII or binary Gmsh MSH v2.2 or v4.1 file. Arbitrary positive node tags are
compacted to `1:N`; oriented connectivity, all fixed-node element types,
serializable cut/sub-element types, physical names, and complete v4 entity
metadata are preserved. Variable-connectivity types 34, 35 and 69 are defined
only by MSH2 ASCII records; MSH2 binary and MSH4 have no per-record width for
them and are rejected explicitly. MSH2 parent/domain element tags become
[`ElementRef`](@ref)s. MSH4 has no such link fields, so linked special records
cannot originate in that format.
[`ElementBlock`](@ref) `tags` remains the legacy projection to the first entity
physical membership (or zero), while [`MixedEntityData`](@ref) retains every
membership, signed boundary, classification, parametric coordinate, and
external node/element tag. Node, element, connectivity, block, entity,
physical-name and file-size limits are checked before bulk allocation. Repeated
physical-name, entity, pre-element node and element sections are merged, with
cumulative resource limits and global external-tag uniqueness. Gmsh treats
backslashes in physical names literally. Set
`tessella_extensions=true` only when reading Tessella's escaped name extension
produced by `write_mixed_msh(...; gmsh_compatible=false)`.

Binary files are byte-swapped when their 32-bit endianness marker requires it.
MSH 2.2 binary payloads use 32-bit integer tags and 64-bit coordinates; MSH
4.1 binary payloads use 32-bit entity/type values, 64-bit `size_t` values and
64-bit coordinates. Unsupported sections in a binary file are rejected
explicitly because their payload cannot be skipped safely without decoding it.
"""
function read_mixed_msh(path::AbstractString;
                        tessella_extensions=false,
                        max_nodes=typemax(Int32),
                        max_elements=typemax(Int32),
                        max_connectivity=typemax(Int32)-1,
                        max_blocks=typemax(Int32),
                        max_entities=typemax(Int32),
                        max_physical_names=typemax(Int32),
                        max_name_bytes=DEFAULT_MAX_MIXED_NAME_BYTES,
                        max_file_bytes=typemax(Int))
    isfile(path) || throw(ArgumentError("read_mixed_msh: missing regular file $path"))
    tessella_extensions isa Bool || throw(ArgumentError(
        "read_mixed_msh: tessella_extensions must be Bool"))
    limits=_MixedReadLimits(
        _read_limit(max_nodes,"max_nodes";ceiling=typemax(Int32)),
        _read_limit(max_elements,"max_elements";ceiling=typemax(Int32)),
        _read_limit(max_connectivity,"max_connectivity";ceiling=typemax(Int32)-1),
        _read_limit(max_blocks,"max_blocks";ceiling=typemax(Int32)),
        _read_limit(max_entities,"max_entities";ceiling=typemax(Int32)),
        _read_limit(max_physical_names,"max_physical_names";ceiling=typemax(Int32)),
        _read_limit(max_name_bytes,"max_name_bytes"),
        _read_limit(max_file_bytes,"max_file_bytes"))
    filesize(path)<=limits.max_file_bytes || throw(ArgumentError(
        "read_mixed_msh: file exceeds max_file_bytes=$(limits.max_file_bytes)"))
    open(path,"r") do io
        try
            return _read_mixed_stream(io,limits,tessella_extensions)
        catch err
            err isa InterruptException && rethrow()
            err isa EOFError && throw(ArgumentError("read_mixed_msh: truncated MSH file"))
            rethrow()
        end
    end
end

@inline function _msh_line(io,context::AbstractString)
    eof(io) && throw(ArgumentError("read_mixed_msh: unexpected end of file in $context"))
    return readline(io)
end

function _msh_int(token::AbstractString,context::AbstractString)
    try
        return parse(Int,token)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: invalid integer in $context"))
    end
end

function _msh_size_t(token::AbstractString,context::AbstractString)
    try
        return parse(UInt64,token)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: invalid unsigned 64-bit value in $context"))
    end
end

function _msh_float(token::AbstractString,context::AbstractString)
    value = try
        parse(Float64,token)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: invalid floating-point value in $context"))
    end
    isfinite(value) || throw(ArgumentError(
        "read_mixed_msh: non-finite floating-point value in $context"))
    return value
end

@inline function _binary_available(io,bytes::Int,context::AbstractString)
    bytes>=0 || throw(ArgumentError(
        "read_mixed_msh: binary byte count overflows Int in $context"))
    remaining=try
        Base.checked_sub(filesize(io),position(io))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "read_mixed_msh: cannot determine remaining bytes in $context"))
    end
    bytes<=remaining || throw(ArgumentError(
        "read_mixed_msh: truncated binary payload in $context"))
    return nothing
end

function _binary_bytes(count::Int,width::Int,context::AbstractString)
    (count>=0&&width>=0) || throw(ArgumentError(
        "read_mixed_msh: negative binary extent in $context"))
    return try
        Base.checked_mul(count,width)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "read_mixed_msh: binary payload size overflows Int in $context"))
    end
end

function _binary_count(value::UInt64,maximum::Int,context::AbstractString)
    value<=UInt64(maximum) || throw(ArgumentError(
        "read_mixed_msh: $context count $value exceeds allowed range 0:$maximum"))
    return Int(value)
end

@inline function _binary_u32(io,swap::Bool,context::AbstractString)
    _binary_available(io,sizeof(UInt32),context)
    return _binary_u32_unchecked(io,swap)
end

@inline function _binary_u32_unchecked(io,swap::Bool)
    value=read(io,UInt32)
    return swap ? bswap(value) : value
end

@inline function _binary_i32(io,swap::Bool,context::AbstractString)
    return reinterpret(Int32,_binary_u32(io,swap,context))
end

@inline function _binary_i32_unchecked(io,swap::Bool)
    return reinterpret(Int32,_binary_u32_unchecked(io,swap))
end

@inline function _binary_u64(io,swap::Bool,context::AbstractString)
    _binary_available(io,sizeof(UInt64),context)
    value=read(io,UInt64)
    return swap ? bswap(value) : value
end

@inline function _binary_f64(io,swap::Bool,context::AbstractString)
    value=reinterpret(Float64,_binary_u64(io,swap,context))
    isfinite(value) || throw(ArgumentError(
        "read_mixed_msh: non-finite floating-point value in $context"))
    return value
end

@inline function _binary_f64_unchecked(io,swap::Bool,context::AbstractString)
    bits=read(io,UInt64); swap && (bits=bswap(bits))
    value=reinterpret(Float64,bits)
    isfinite(value) || throw(ArgumentError(
        "read_mixed_msh: non-finite floating-point value in $context"))
    return value
end

function _binary_i32_vector(io,count::Int,swap::Bool,context::AbstractString)
    bytes=_binary_bytes(count,sizeof(Int32),context)
    _binary_available(io,bytes,context)
    values=Vector{Int32}(undef,count)
    read!(io,values)
    if swap
        @inbounds for i in eachindex(values)
            values[i]=bswap(values[i])
        end
    end
    return values
end

function _binary_u64_vector(io,count::Int,swap::Bool,context::AbstractString)
    bytes=_binary_bytes(count,sizeof(UInt64),context)
    _binary_available(io,bytes,context)
    values=Vector{UInt64}(undef,count)
    read!(io,values)
    if swap
        @inbounds for i in eachindex(values)
            values[i]=bswap(values[i])
        end
    end
    return values
end

function _binary_f64_vector(io,count::Int,swap::Bool,context::AbstractString)
    bytes=_binary_bytes(count,sizeof(Float64),context)
    _binary_available(io,bytes,context)
    values=Vector{Float64}(undef,count)
    read!(io,values)
    if swap
        bits=reinterpret(UInt64,values)
        @inbounds for i in eachindex(bits)
            bits[i]=bswap(bits[i])
        end
    end
    all(isfinite,values) || throw(ArgumentError(
        "read_mixed_msh: non-finite floating-point value in $context"))
    return values
end

function _consume_binary_newline(io,context::AbstractString)
    _binary_available(io,1,context)
    byte=read(io,UInt8)
    if byte==UInt8('\r')
        _binary_available(io,1,context)
        read(io,UInt8)==UInt8('\n') || throw(ArgumentError(
            "read_mixed_msh: expected newline after binary $context"))
    elseif byte!=UInt8('\n')
        throw(ArgumentError(
            "read_mixed_msh: expected newline after binary $context"))
    end
    return nothing
end

_msh_fields(io,context::AbstractString)=split(strip(_msh_line(io,context)))

const _MSH_NUMERIC_TOKEN_MAX_BYTES = 64

mutable struct _MshRecordReader{T<:IO}
    io::T
    ended::Bool
end

_MshRecordReader(io::T) where {T<:IO}=_MshRecordReader{T}(io,false)

@inline function _msh_record_magnitude!(
    reader::_MshRecordReader,context::AbstractString,
    positive_limit::UInt64,negative_limit::UInt64,allow_negative::Bool)
    reader.ended && throw(ArgumentError(
        "read_mixed_msh: missing token in $context"))
    started=false
    signed=false
    negative=false
    digits=0
    token_bytes=0
    value=UInt64(0)
    while true
        if eof(reader.io)
            reader.ended=true
            started || throw(ArgumentError(
                "read_mixed_msh: unexpected end of file in $context"))
            break
        end
        byte=read(reader.io,UInt8)
        if byte==UInt8('\n')
            reader.ended=true
            started || throw(ArgumentError(
                "read_mixed_msh: missing token in $context"))
            break
        elseif byte==UInt8('\r')
            eof(reader.io) && throw(ArgumentError(
                "read_mixed_msh: truncated CRLF in $context"))
            read(reader.io,UInt8)==UInt8('\n') || throw(ArgumentError(
                "read_mixed_msh: expected LF after CR in $context"))
            reader.ended=true
            started || throw(ArgumentError(
                "read_mixed_msh: missing token in $context"))
            break
        elseif byte==UInt8(' ') || byte==UInt8('\t')
            started && break
        else
            started=true
            token_bytes+=1
            token_bytes<=_MSH_NUMERIC_TOKEN_MAX_BYTES || throw(ArgumentError(
                "read_mixed_msh: numeric token in $context exceeds " *
                "$_MSH_NUMERIC_TOKEN_MAX_BYTES bytes"))
            if digits==0 && !signed && (byte==UInt8('+') || byte==UInt8('-'))
                signed=true
                negative=byte==UInt8('-')
                negative && !allow_negative && throw(ArgumentError(
                    "read_mixed_msh: invalid unsigned 64-bit value in $context"))
                continue
            end
            UInt8('0')<=byte<=UInt8('9') || throw(ArgumentError(
                "read_mixed_msh: invalid integer in $context"))
            digits+=1
            digit=UInt64(byte-UInt8('0'))
            limit=negative ? negative_limit : positive_limit
            value<=(limit-digit)÷UInt64(10) || throw(ArgumentError(
                "read_mixed_msh: integer outside supported range in $context"))
            value=UInt64(10)*value+digit
        end
    end
    digits>0 || throw(ArgumentError(
        "read_mixed_msh: invalid integer in $context"))
    return negative,value
end

@inline function _msh_record_int!(reader::_MshRecordReader,
                                  context::AbstractString)
    positive=UInt64(typemax(Int))
    negative=positive+UInt64(1)
    isnegative,magnitude=_msh_record_magnitude!(
        reader,context,positive,negative,true)
    isnegative || return Int(magnitude)
    magnitude==negative && return typemin(Int)
    return -Int(magnitude)
end

@inline function _msh_record_size_t!(reader::_MshRecordReader,
                                     context::AbstractString)
    _,value=_msh_record_magnitude!(
        reader,context,typemax(UInt64),UInt64(0),false)
    return value
end

function _expect_msh_record_end!(reader::_MshRecordReader,context::AbstractString)
    reader.ended && return nothing
    while !eof(reader.io)
        byte=read(reader.io,UInt8)
        if byte==UInt8(' ') || byte==UInt8('\t')
            continue
        elseif byte==UInt8('\n')
            reader.ended=true
            return nothing
        elseif byte==UInt8('\r')
            eof(reader.io) && throw(ArgumentError(
                "read_mixed_msh: truncated CRLF after $context"))
            read(reader.io,UInt8)==UInt8('\n') || throw(ArgumentError(
                "read_mixed_msh: expected LF after CR in $context"))
            reader.ended=true
            return nothing
        end
        throw(ArgumentError(
            "read_mixed_msh: unexpected content after $context"))
    end
    reader.ended=true
    return nothing
end

mutable struct _MshTokenReader
    io::IO
    tokens::Vector{SubString{String}}
    next::Int
end

_MshTokenReader(io::IO)=_MshTokenReader(io,SubString{String}[],1)

function _msh_token!(reader::_MshTokenReader,context::AbstractString)
    while reader.next>length(reader.tokens)
        line=_msh_line(reader.io,context)
        reader.tokens=split(line)
        reader.next=1
    end
    token=reader.tokens[reader.next]
    reader.next+=1
    return token
end

function _expect_msh_token_end!(reader::_MshTokenReader,token::AbstractString)
    actual=_msh_token!(reader,token)
    actual==token || throw(ArgumentError(
        "read_mixed_msh: expected $token, got '$actual'"))
    reader.next>length(reader.tokens) || throw(ArgumentError(
        "read_mixed_msh: unexpected content after $token"))
    return nothing
end

function _expect_msh_end(io,token::AbstractString)
    line=strip(_msh_line(io,token))
    line==token || throw(ArgumentError(
        "read_mixed_msh: expected $token, got '$line'"))
    return nothing
end

function _skip_msh_section(io,endtoken::AbstractString)
    while !eof(io)
        strip(readline(io))==endtoken && return nothing
    end
    throw(ArgumentError("read_mixed_msh: unterminated section; missing $endtoken"))
end

function _section_count(io,context::AbstractString,maximum::Int)
    fields=_msh_fields(io,context)
    length(fields)==1 || throw(ArgumentError(
        "read_mixed_msh: malformed $context count"))
    n=_msh_int(fields[1],context)
    0<=n<=maximum || throw(ArgumentError(
        "read_mixed_msh: $context count $n exceeds allowed range 0:$maximum"))
    return n
end

function _read_mixed_stream(io,limits::_MixedReadLimits,tessella_extensions::Bool)
    acc=_MixedReadAccum()
    version=0.0
    binary=false; swap=false
    seen_format=false; seen_names=false; seen_entities=false
    seen_nodes=false; seen_elements=false
    while !eof(io)
        header=strip(readline(io)); isempty(header) && continue
        if header=="\$MeshFormat"
            seen_format && throw(ArgumentError("read_mixed_msh: duplicate \$MeshFormat section"))
            (seen_names||seen_entities||seen_nodes||seen_elements) && throw(ArgumentError(
                "read_mixed_msh: \$MeshFormat must be the first MSH section"))
            seen_format=true
            fields=_msh_fields(io,"MeshFormat header")
            length(fields)==3 || throw(ArgumentError(
                "read_mixed_msh: malformed MeshFormat header"))
            version=_msh_float(fields[1],"MeshFormat version")
            version in (2.2,4.1) || throw(ArgumentError(
                "read_mixed_msh: unsupported MSH version $version (supported: 2.2, 4.1)"))
            file_type=_msh_int(fields[2],"MeshFormat file type")
            file_type in (0,1) || throw(ArgumentError(
                "read_mixed_msh: MeshFormat file type must be 0 or 1"))
            data_size=_msh_int(fields[3],"MeshFormat data size")
            data_size==8 || throw(ArgumentError(
                "read_mixed_msh: only an 8-byte MeshFormat data size is supported"))
            binary=file_type==1
            if binary
                raw_marker=_binary_u32(io,false,"MeshFormat endianness marker")
                if raw_marker==UInt32(1)
                    swap=false
                elseif bswap(raw_marker)==UInt32(1)
                    swap=true
                else
                    throw(ArgumentError(
                        "read_mixed_msh: invalid binary endianness marker"))
                end
                _consume_binary_newline(io,"MeshFormat endianness marker")
            end
            _expect_msh_end(io,"\$EndMeshFormat")
        elseif header=="\$PhysicalNames"
            seen_format || throw(ArgumentError(
                "read_mixed_msh: \$PhysicalNames appeared before \$MeshFormat"))
            seen_names=true
            _read_mixed_physical_names!(acc,io,limits,tessella_extensions)
        elseif header=="\$Entities"
            seen_format || throw(ArgumentError(
                "read_mixed_msh: \$Entities appeared before \$MeshFormat"))
            version==4.1 || throw(ArgumentError(
                "read_mixed_msh: \$Entities is only valid in MSH v4.1"))
            (seen_nodes||seen_elements) && throw(ArgumentError(
                "read_mixed_msh: \$Entities must precede nodes and elements"))
            seen_entities=true
            binary ? _read_mixed_entities_v4_binary!(acc,io,limits,swap) :
                     _read_mixed_entities_v4!(acc,io,limits)
        elseif header=="\$Nodes"
            seen_format || throw(ArgumentError(
                "read_mixed_msh: \$Nodes appeared before \$MeshFormat"))
            seen_elements && throw(ArgumentError(
                "read_mixed_msh: \$Nodes must precede \$Elements"))
            seen_nodes=true
            if version==2.2
                binary ? _read_mixed_nodes_v2_binary!(acc,io,limits,swap) :
                         _read_mixed_nodes_v2!(acc,io,limits)
            else
                binary ? _read_mixed_nodes_v4_binary!(acc,io,limits,swap) :
                         _read_mixed_nodes_v4!(acc,io,limits)
            end
        elseif header=="\$Elements"
            seen_nodes || throw(ArgumentError(
                "read_mixed_msh: \$Elements appeared before \$Nodes"))
            seen_elements=true
            if version==2.2
                binary ? _read_mixed_elements_v2_binary!(acc,io,limits,swap) :
                         _read_mixed_elements_v2!(acc,io,limits)
            else
                binary ? _read_mixed_elements_v4_binary!(acc,io,limits,swap) :
                         _read_mixed_elements_v4!(acc,io,limits)
            end
        elseif startswith(header,"\$End")
            throw(ArgumentError("read_mixed_msh: unexpected section terminator $header"))
        elseif startswith(header,"\$")
            seen_format || throw(ArgumentError(
                "read_mixed_msh: section $header appeared before \$MeshFormat"))
            binary && throw(ArgumentError(
                "read_mixed_msh: unsupported section $header in binary MSH"))
            _skip_msh_section(io,"\$End"*header[2:end])
        else
            throw(ArgumentError("read_mixed_msh: unexpected content outside a section"))
        end
    end
    seen_format || throw(ArgumentError("read_mixed_msh: missing \$MeshFormat section"))
    (seen_nodes||version==4.1) || throw(ArgumentError(
        "read_mixed_msh: missing \$Nodes section"))
    return _finish_mixed_read(acc,version==4.1)
end

function _unescape_mixed_name(encoded::AbstractString)
    out=IOBuffer(); i=firstindex(encoded)
    while i<=lastindex(encoded)
        c=encoded[i]
        if c=='\\'
            i=nextind(encoded,i)
            i<=lastindex(encoded) || throw(ArgumentError(
                "read_mixed_msh: trailing escape in physical name"))
            e=encoded[i]
            if e=='\\' || e=='\"'
                write(out,e)
            elseif e=='n'
                write(out,'\n')
            elseif e=='r'
                write(out,'\r')
            elseif e=='t'
                write(out,'\t')
            else
                # Gmsh itself treats backslashes literally. Preserve unknown
                # pairs so paths and names from non-Tessella writers survive.
                write(out,'\\'); write(out,e)
            end
        else
            c=='\"' && throw(ArgumentError(
                "read_mixed_msh: unescaped quote in physical name"))
            write(out,c)
        end
        i=nextind(encoded,i)
    end
    return String(take!(out))
end

function _read_mixed_physical_names!(acc,io,limits,tessella_extensions::Bool)
    n=_section_count(io,"physical-name",limits.max_physical_names)
    acc.physical_name_records=try
        Base.checked_add(acc.physical_name_records,n)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: physical-name record count overflows Int"))
    end
    acc.physical_name_records<=limits.max_physical_names || throw(ArgumentError(
        "read_mixed_msh: cumulative physical-name count exceeds max_physical_names"))
    pattern=r"^([+-]?\d+)\s+([+-]?\d+)\s+(.*)$"
    for _ in 1:n
        line=_msh_line(io,"PhysicalNames record")
        isvalid(line) || throw(ArgumentError(
            "read_mixed_msh: physical-name record is not valid UTF-8"))
        line_limit=limits.max_name_bytes>typemax(Int)-128 ? typemax(Int) :
                   limits.max_name_bytes+128
        ncodeunits(line)<=line_limit || throw(ArgumentError(
            "read_mixed_msh: physical-name record exceeds max_name_bytes"))
        record=match(pattern,line)
        record===nothing && throw(ArgumentError(
            "read_mixed_msh: malformed physical-name record"))
        dim=_msh_int(record.captures[1],"physical-name dimension")
        tag=_msh_int(record.captures[2],"physical-name tag")
        0<=dim<=3 || throw(ArgumentError(
            "read_mixed_msh: physical-name dimension $dim is outside 0:3"))
        1<=tag<=typemax(Int32) || throw(ArgumentError(
            "read_mixed_msh: physical-name tag must be positive and fit Int32"))
        key=(dim,tag)
        quoted=strip(record.captures[3])
        ncodeunits(quoted)>=2 && first(quoted)=='\"' && last(quoted)=='\"' ||
            throw(ArgumentError("read_mixed_msh: malformed quoted physical name"))
        encoded=ncodeunits(quoted)==2 ? "" :
            SubString(quoted,nextind(quoted,firstindex(quoted)),
                      prevind(quoted,lastindex(quoted)))
        name=tessella_extensions ? _unescape_mixed_name(encoded) : String(encoded)
        !tessella_extensions && occursin('\"',name) && throw(ArgumentError(
            "read_mixed_msh: unescaped quote in physical name"))
        ncodeunits(name)<=limits.max_name_bytes || throw(ArgumentError(
            "read_mixed_msh: physical name exceeds max_name_bytes"))
        occursin('\0',name) && throw(ArgumentError(
            "read_mixed_msh: physical name contains a NUL byte"))
        if haskey(acc.physical_names,key)
            acc.physical_names[key]==name || throw(ArgumentError(
                "read_mixed_msh: conflicting physical names for ($dim,$tag)"))
        else
            acc.physical_names[key]=name
        end
    end
    _expect_msh_end(io,"\$EndPhysicalNames")
end

function _read_mixed_entities_v4!(acc,io,limits)
    reader=_MshTokenReader(io)
    counts=ntuple(4) do _
        _msh_int(_msh_token!(reader,"Entities header"),"Entities header")
    end
    all(>=(0),counts) || throw(ArgumentError(
        "read_mixed_msh: negative v4 entity count"))
    total=0
    for count in counts
        total=try Base.checked_add(total,count) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("read_mixed_msh: v4 entity count overflows Int"))
        end
    end
    existing=try Base.checked_add(length(acc.entities),length(acc.implicit_entities)) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: cumulative v4 entity count overflows Int"))
    end
    cumulative=try Base.checked_add(existing,total) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: cumulative v4 entity count overflows Int"))
    end
    cumulative<=limits.max_entities || throw(ArgumentError(
        "read_mixed_msh: cumulative v4 entity count $cumulative exceeds " *
        "max_entities=$(limits.max_entities)"))
    for dim in 0:3
        for _ in 1:counts[dim+1]
            tag=_msh_int(_msh_token!(reader,"dimension-$dim entity"),"entity tag")
            1<=tag<=typemax(Int32) || throw(ArgumentError(
                "read_mixed_msh: entity tags must be positive and fit Int32"))
            key=(dim,tag); haskey(acc.entities,key) && throw(ArgumentError(
                "read_mixed_msh: duplicate entity ($dim,$tag)"))
            ncoordinates=dim==0 ? 3 : 6
            box=ntuple(ncoordinates) do _
                _msh_float(_msh_token!(reader,"entity bounds"),"entity bounds")
            end
            if dim>0
                all(box[i]<=box[i+3] for i in 1:3) || throw(ArgumentError(
                    "read_mixed_msh: entity ($dim,$tag) has reversed bounds"))
            end
            nphysical=_msh_int(
                _msh_token!(reader,"entity physical-tag count"),
                "entity physical-tag count")
            nphysical>=0 || throw(ArgumentError(
                "read_mixed_msh: negative physical-tag count on entity ($dim,$tag)"))
            physical_tags=Int32[]
            for _ in 1:nphysical
                value=_msh_int(
                    _msh_token!(reader,"entity physical tag"),"entity physical tag")
                1<=value<=typemax(Int32) || throw(ArgumentError(
                    "read_mixed_msh: entity physical tag must be positive and fit Int32"))
                push!(physical_tags,Int32(value))
            end
            boundaries=Int32[]
            if dim>0
                nboundary=_msh_int(
                    _msh_token!(reader,"entity boundary count"),
                    "entity boundary count")
                nboundary>=0 || throw(ArgumentError(
                    "read_mixed_msh: negative boundary count on entity ($dim,$tag)"))
                for _ in 1:nboundary
                    boundary=_msh_int(
                        _msh_token!(reader,"entity boundary tag"),"entity boundary tag")
                    boundary!=0 || throw(ArgumentError(
                        "read_mixed_msh: entity boundary tag cannot be zero"))
                    boundary_tag=try abs(boundary) catch err
                        err isa InterruptException && rethrow()
                        throw(ArgumentError("read_mixed_msh: entity boundary tag overflows Int"))
                    end
                    boundary_tag<=typemax(Int32) || throw(ArgumentError(
                        "read_mixed_msh: entity boundary tag does not fit Int32"))
                    haskey(acc.entities,(dim-1,boundary_tag)) || throw(ArgumentError(
                        "read_mixed_msh: entity ($dim,$tag) references undeclared boundary entity"))
                    push!(boundaries,Int32(boundary))
                end
            end
            record=MixedEntity(dim,tag,box;
                               physical_tags=physical_tags,boundaries=boundaries)
            acc.entities[key]=record
            acc.entity_physical[key]=isempty(physical_tags) ? Int32(0) : first(physical_tags)
        end
    end
    _expect_msh_token_end!(reader,"\$EndEntities")
end

function _read_mixed_entities_v4_binary!(acc,io,limits,swap::Bool)
    counts=ntuple(4) do dim
        raw=_binary_u64(io,swap,"v4 binary Entities header")
        _binary_count(raw,limits.max_entities,"dimension-$(dim-1) entity")
    end
    total=0
    for count in counts
        total=try Base.checked_add(total,count) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("read_mixed_msh: v4 entity count overflows Int"))
        end
    end
    existing=try Base.checked_add(length(acc.entities),length(acc.implicit_entities)) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: cumulative v4 entity count overflows Int"))
    end
    cumulative=try Base.checked_add(existing,total) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: cumulative v4 entity count overflows Int"))
    end
    cumulative<=limits.max_entities || throw(ArgumentError(
        "read_mixed_msh: cumulative v4 entity count $cumulative exceeds " *
        "max_entities=$(limits.max_entities)"))

    for dim in 0:3
        for _ in 1:counts[dim+1]
            tag=Int(_binary_i32(io,swap,"dimension-$dim entity tag"))
            1<=tag<=typemax(Int32) || throw(ArgumentError(
                "read_mixed_msh: entity tags must be positive and fit Int32"))
            key=(dim,tag)
            haskey(acc.entities,key) && throw(ArgumentError(
                "read_mixed_msh: duplicate entity ($dim,$tag)"))
            ncoordinates=dim==0 ? 3 : 6
            coordinates=_binary_f64_vector(
                io,ncoordinates,swap,"dimension-$dim entity bounds")
            box=Tuple(coordinates)
            if dim>0
                all(box[i]<=box[i+3] for i in 1:3) || throw(ArgumentError(
                    "read_mixed_msh: entity ($dim,$tag) has reversed bounds"))
            end
            nphysical=_binary_count(
                _binary_u64(io,swap,"entity physical-tag count"),
                typemax(Int),"entity physical-tag")
            physical_tags=_binary_i32_vector(
                io,nphysical,swap,"entity physical tags")
            @inbounds for value in physical_tags
                value>=1 || throw(ArgumentError(
                    "read_mixed_msh: entity physical tag must be positive and fit Int32"))
            end
            boundaries=Int32[]
            if dim>0
                nboundary=_binary_count(
                    _binary_u64(io,swap,"entity boundary count"),
                    typemax(Int),"entity boundary")
                boundaries=_binary_i32_vector(
                    io,nboundary,swap,"entity boundary tags")
                @inbounds for boundary in boundaries
                    boundary!=0 || throw(ArgumentError(
                        "read_mixed_msh: entity boundary tag cannot be zero"))
                    boundary_tag=abs(Int(boundary))
                    boundary_tag<=typemax(Int32) || throw(ArgumentError(
                        "read_mixed_msh: entity boundary tag does not fit Int32"))
                    haskey(acc.entities,(dim-1,boundary_tag)) || throw(ArgumentError(
                        "read_mixed_msh: entity ($dim,$tag) references undeclared boundary entity"))
                end
            end
            record=MixedEntity(dim,tag,box;
                               physical_tags=physical_tags,boundaries=boundaries)
            acc.entities[key]=record
            acc.entity_physical[key]=isempty(physical_tags) ? Int32(0) :
                                     first(physical_tags)
        end
    end
    _consume_binary_newline(io,"Entities section")
    _expect_msh_end(io,"\$EndEntities")
    return nothing
end

function _read_mixed_nodes_v2!(acc,io,limits)
    remaining=limits.max_nodes-length(acc.x)
    count=_section_count(io,"v2 node",remaining)
    for _ in 1:count
        fields=_msh_fields(io,"v2 node record")
        length(fields)==4 || throw(ArgumentError(
            "read_mixed_msh: v2 node records require exactly 4 fields"))
        tag=_msh_size_t(fields[1],"v2 node tag")
        tag>0 || throw(ArgumentError("read_mixed_msh: node tags must be positive"))
        haskey(acc.node_map,tag) && throw(ArgumentError(
            "read_mixed_msh: duplicate node tag $tag"))
        push!(acc.x,_msh_float(fields[2],"node coordinate"))
        push!(acc.y,_msh_float(fields[3],"node coordinate"))
        push!(acc.z,_msh_float(fields[4],"node coordinate"))
        acc.node_map[tag]=Int32(length(acc.x))
    end
    _expect_msh_end(io,"\$EndNodes")
end

function _read_mixed_nodes_v2_binary!(acc,io,limits,swap::Bool)
    remaining=limits.max_nodes-length(acc.x)
    count=_section_count(io,"v2 node",remaining)
    record_bytes=sizeof(Int32)+3sizeof(Float64)
    _binary_available(io,_binary_bytes(count,record_bytes,"v2 binary Nodes"),
                      "v2 binary Nodes")
    target=try Base.checked_add(length(acc.x),count) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: cumulative v2 node count overflows Int"))
    end
    sizehint!(acc.x,target); sizehint!(acc.y,target); sizehint!(acc.z,target)
    sizehint!(acc.node_map,target)
    for _ in 1:count
        signed_tag=_binary_i32_unchecked(io,swap)
        signed_tag>0 || throw(ArgumentError(
            "read_mixed_msh: node tags must be positive"))
        tag=UInt64(signed_tag)
        haskey(acc.node_map,tag) && throw(ArgumentError(
            "read_mixed_msh: duplicate node tag $tag"))
        push!(acc.x,_binary_f64_unchecked(io,swap,"node coordinate"))
        push!(acc.y,_binary_f64_unchecked(io,swap,"node coordinate"))
        push!(acc.z,_binary_f64_unchecked(io,swap,"node coordinate"))
        acc.node_map[tag]=Int32(length(acc.x))
    end
    _consume_binary_newline(io,"Nodes section")
    _expect_msh_end(io,"\$EndNodes")
    return nothing
end

function _ensure_mixed_implicit_entity!(acc,dim::Int,tag::Int,limits)
    key=(dim,tag)
    (haskey(acc.entities,key) || key in acc.implicit_entities) && return nothing
    length(acc.entities)+length(acc.implicit_entities)<limits.max_entities || throw(ArgumentError(
        "read_mixed_msh: implicit entity count exceeds max_entities=$(limits.max_entities)"))
    acc.entity_physical[key]=Int32(0)
    push!(acc.implicit_entities,key)
    return nothing
end

function _read_mixed_nodes_v4!(acc,io,limits)
    reader=_MshTokenReader(io)
    nblocks=_msh_int(
        _msh_token!(reader,"v4 Nodes header"),"v4 node-block count")
    count=_msh_int(
        _msh_token!(reader,"v4 Nodes header"),"v4 node count")
    declared_min=_msh_size_t(
        _msh_token!(reader,"v4 Nodes header"),"v4 node tag range")
    declared_max=_msh_size_t(
        _msh_token!(reader,"v4 Nodes header"),"v4 node tag range")
    nblocks>=0 || throw(ArgumentError(
        "read_mixed_msh: negative v4 node-block count"))
    count>=0 || throw(ArgumentError("read_mixed_msh: negative v4 node count"))
    cumulative_blocks=try Base.checked_add(acc.node_blocks,nblocks) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: cumulative v4 node-block count overflows Int"))
    end
    cumulative_blocks<=limits.max_blocks || throw(ArgumentError(
        "read_mixed_msh: cumulative v4 node-block count exceeds max_blocks"))
    count<=limits.max_nodes-length(acc.x) || throw(ArgumentError(
        "read_mixed_msh: cumulative v4 node count exceeds max_nodes"))
    if count==0
        (declared_min==0&&declared_max==0) || throw(ArgumentError(
            "read_mixed_msh: empty v4 Nodes must declare range 0 0"))
    else
        0<declared_min<=declared_max || throw(ArgumentError(
            "read_mixed_msh: invalid v4 node-tag range"))
    end
    actual_min=typemax(UInt64); actual_max=zero(UInt64); nread=0
    for _ in 1:nblocks
        dim=_msh_int(_msh_token!(reader,"v4 node-block header"),
                     "node-block dimension")
        entity=_msh_int(_msh_token!(reader,"v4 node-block header"),
                        "node-block entity")
        parametric=_msh_int(_msh_token!(reader,"v4 node-block header"),
                            "node-block parametric flag")
        nlocal=_msh_int(_msh_token!(reader,"v4 node-block header"),
                        "node-block count")
        0<=dim<=3 || throw(ArgumentError(
            "read_mixed_msh: node-block dimension $dim is outside 0:3"))
        1<=entity<=typemax(Int32) || throw(ArgumentError(
            "read_mixed_msh: node-block entity tags must be positive and fit Int32"))
        # Gmsh creates a discrete entity when a node block refers to one that
        # was not declared in an optional $Entities section.
        _ensure_mixed_implicit_entity!(acc,dim,entity,limits)
        parametric in (0,1) || throw(ArgumentError(
            "read_mixed_msh: node-block parametric flag must be 0 or 1"))
        0<=nlocal<=count-nread || throw(ArgumentError(
            "read_mixed_msh: invalid or excessive v4 node-block size"))
        tags=UInt64[]
        for i in 1:nlocal
            tag=_msh_size_t(_msh_token!(reader,"v4 node tag"),"v4 node tag")
            tag>0 || throw(ArgumentError("read_mixed_msh: node tags must be positive"))
            haskey(acc.node_map,tag) && throw(ArgumentError(
                "read_mixed_msh: duplicate node tag $tag"))
            push!(tags,tag); acc.node_map[tag]=Int32(0)
            actual_min=min(actual_min,tag); actual_max=max(actual_max,tag)
        end
        for i in 1:nlocal
            x=_msh_float(_msh_token!(reader,"v4 node coordinates"),"node coordinate")
            y=_msh_float(_msh_token!(reader,"v4 node coordinates"),"node coordinate")
            z=_msh_float(_msh_token!(reader,"v4 node coordinates"),"node coordinate")
            parameters=parametric==1 ? Float64[] : nothing
            if parameters!==nothing
                for _ in 1:dim
                    push!(parameters,_msh_float(
                        _msh_token!(reader,"parametric node coordinate"),
                        "parametric node coordinate"))
                end
            end
            push!(acc.x,x); push!(acc.y,y); push!(acc.z,z)
            push!(acc.external_node_tags,tags[i])
            push!(acc.node_entities,(dim,Int32(entity)))
            push!(acc.node_parametric,parameters)
            acc.node_map[tags[i]]=Int32(length(acc.x)); nread+=1
        end
    end
    nread==count || throw(ArgumentError(
        "read_mixed_msh: v4 Nodes declared $count nodes but contained $nread"))
    count==0 || (actual_min==declared_min&&actual_max==declared_max) || throw(ArgumentError(
        "read_mixed_msh: v4 node-tag range does not match records"))
    _expect_msh_token_end!(reader,"\$EndNodes")
    acc.node_blocks=cumulative_blocks
end

function _read_mixed_nodes_v4_binary!(acc,io,limits,swap::Bool)
    nblocks=_binary_count(
        _binary_u64(io,swap,"v4 binary Nodes header"),
        limits.max_blocks,"v4 node-block")
    remaining_nodes=limits.max_nodes-length(acc.x)
    count=_binary_count(
        _binary_u64(io,swap,"v4 binary Nodes header"),
        remaining_nodes,"v4 node")
    declared_min=_binary_u64(io,swap,"v4 binary node tag range")
    declared_max=_binary_u64(io,swap,"v4 binary node tag range")
    cumulative_blocks=try Base.checked_add(acc.node_blocks,nblocks) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "read_mixed_msh: cumulative v4 node-block count overflows Int"))
    end
    cumulative_blocks<=limits.max_blocks || throw(ArgumentError(
        "read_mixed_msh: cumulative v4 node-block count exceeds max_blocks"))
    if count==0
        (declared_min==0&&declared_max==0) || throw(ArgumentError(
            "read_mixed_msh: empty v4 Nodes must declare range 0 0"))
    else
        0<declared_min<=declared_max || throw(ArgumentError(
            "read_mixed_msh: invalid v4 node-tag range"))
    end
    target=try Base.checked_add(length(acc.x),count) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("read_mixed_msh: cumulative v4 node count overflows Int"))
    end
    for values in (acc.x,acc.y,acc.z,acc.external_node_tags,
                   acc.node_entities,acc.node_parametric)
        sizehint!(values,target)
    end
    sizehint!(acc.node_map,target)

    actual_min=typemax(UInt64); actual_max=zero(UInt64); nread=0
    for _ in 1:nblocks
        dim=Int(_binary_i32(io,swap,"v4 binary node-block dimension"))
        entity=Int(_binary_i32(io,swap,"v4 binary node-block entity"))
        parametric=Int(_binary_i32(io,swap,"v4 binary node-block parametric flag"))
        nlocal=_binary_count(
            _binary_u64(io,swap,"v4 binary node-block count"),
            count-nread,"v4 node-block")
        0<=dim<=3 || throw(ArgumentError(
            "read_mixed_msh: node-block dimension $dim is outside 0:3"))
        1<=entity<=typemax(Int32) || throw(ArgumentError(
            "read_mixed_msh: node-block entity tags must be positive and fit Int32"))
        _ensure_mixed_implicit_entity!(acc,dim,entity,limits)
        parametric in (0,1) || throw(ArgumentError(
            "read_mixed_msh: node-block parametric flag must be 0 or 1"))
        stride=3+(parametric==1 ? dim : 0)
        coordinate_count=try Base.checked_mul(nlocal,stride) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "read_mixed_msh: v4 node coordinate count overflows Int"))
        end
        # Check the complete block payload before allocating either vector.
        tag_bytes=_binary_bytes(nlocal,sizeof(UInt64),"v4 node tags")
        coordinate_bytes=_binary_bytes(
            coordinate_count,sizeof(Float64),"v4 node coordinates")
        block_bytes=try Base.checked_add(tag_bytes,coordinate_bytes) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "read_mixed_msh: v4 node-block byte count overflows Int"))
        end
        _binary_available(io,block_bytes,"v4 node block")
        tags=_binary_u64_vector(io,nlocal,swap,"v4 node tags")
        @inbounds for tag in tags
            tag>0 || throw(ArgumentError(
                "read_mixed_msh: node tags must be positive"))
            haskey(acc.node_map,tag) && throw(ArgumentError(
                "read_mixed_msh: duplicate node tag $tag"))
            acc.node_map[tag]=Int32(0)
            actual_min=min(actual_min,tag); actual_max=max(actual_max,tag)
        end
        coordinates=_binary_f64_vector(
            io,coordinate_count,swap,"v4 node coordinates")
        offset=1
        @inbounds for i in 1:nlocal
            push!(acc.x,coordinates[offset]);
            push!(acc.y,coordinates[offset+1]);
            push!(acc.z,coordinates[offset+2]);
            offset+=3
            parameters=parametric==1 ? Float64[] : nothing
            if parameters!==nothing
                sizehint!(parameters,dim)
                for _ in 1:dim
                    push!(parameters,coordinates[offset]); offset+=1
                end
            end
            push!(acc.external_node_tags,tags[i])
            push!(acc.node_entities,(dim,Int32(entity)))
            push!(acc.node_parametric,parameters)
            acc.node_map[tags[i]]=Int32(length(acc.x)); nread+=1
        end
    end
    nread==count || throw(ArgumentError(
        "read_mixed_msh: v4 Nodes declared $count nodes but contained $nread"))
    count==0 || (actual_min==declared_min&&actual_max==declared_max) || throw(ArgumentError(
        "read_mixed_msh: v4 node-tag range does not match records"))
    _consume_binary_newline(io,"Nodes section")
    _expect_msh_end(io,"\$EndNodes")
    acc.node_blocks=cumulative_blocks
    return nothing
end

function _mixed_physical_tag(value::Int,context::AbstractString)
    0<=value<=typemax(Int32) || throw(ArgumentError(
        "read_mixed_msh: $context must be non-negative and fit Int32"))
    return Int32(value)
end

function _record_mixed_element_tag!(acc,tag::UInt64)
    tag>0 || throw(ArgumentError("read_mixed_msh: element tags must be positive"))
    tag in acc.element_tags && throw(ArgumentError(
        "read_mixed_msh: duplicate element tag $tag"))
    push!(acc.element_tags,tag); return nothing
end

function _msh_element_layout(etype::Int,context::AbstractString)
    fixed=get(MSH_CATALOG,etype,nothing)
    fixed!==nothing && return (
        dim=fixed.dim,nnodes=fixed.nnodes,variable=false,special=false,
        unit=fixed.nnodes,links=:none)
    record=get(MSH_SPECIAL_RECORDS,etype,nothing)
    if record!==nothing
        special=MSH_SPECIAL_TYPES[etype]
        return (dim=special.dim,
                nnodes=record.variable ? nothing : record.unit,
                variable=record.variable,special=true,
                unit=record.unit,links=record.links)
    end
    if haskey(MSH_SPECIAL_TYPES,etype)
        throw(ArgumentError(
            "read_mixed_msh: Gmsh element type $etype is an internal MINI basis selector, not a serializable $context record"))
    end
    throw(ArgumentError("read_mixed_msh: unknown Gmsh element type $etype"))
end

function _reserve_mixed_connectivity!(acc,limits,count::Int,context::AbstractString)
    count>=0 || throw(ArgumentError(
        "read_mixed_msh: negative connectivity count in $context"))
    total=try Base.checked_add(acc.connectivity_entries,count) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "read_mixed_msh: cumulative connectivity count overflows Int"))
    end
    total<=limits.max_connectivity || throw(ArgumentError(
        "read_mixed_msh: cumulative connectivity count exceeds max_connectivity=$(limits.max_connectivity)"))
    acc.connectivity_entries=total
    return nothing
end

function _append_mixed_cell_metadata!(bucket::_MixedReadBucket,layout,
                                      physical::Int32,entity::Int32,
                                      element_tag::UInt64,parent_tag::UInt64,
                                      domain_tags::NTuple{2,UInt64})
    push!(bucket.tags,physical)
    push!(bucket.entities,entity)
    push!(bucket.external_tags,element_tag)
    if layout.special
        layout.links===:parent ? push!(bucket.parent_tags,parent_tag) :
                                push!(bucket.domain_tags,domain_tags)
        layout.variable && push!(bucket.offsets,Int32(length(bucket.nodes)+1))
    end
    return nothing
end

function _push_mixed_cell!(acc,etype::Int,physical::Int32,entity::Int32,
                           element_tag::UInt64,node_tokens,limits;
                           parent_tag::UInt64=UInt64(0),
                           domain_tags::NTuple{2,UInt64}=(UInt64(0),UInt64(0)))
    layout=_msh_element_layout(etype,"element")
    if layout.variable
        length(node_tokens)>0 && length(node_tokens)%layout.unit==0 ||
            throw(ArgumentError(
                "read_mixed_msh: type $etype connectivity must be a positive multiple of $(layout.unit)"))
    else
        length(node_tokens)==layout.nnodes || throw(ArgumentError(
            "read_mixed_msh: type $etype expects $(layout.nnodes) node tags"))
    end
    _reserve_mixed_connectivity!(acc,limits,length(node_tokens),"element record")
    bucket=get!(acc.buckets,etype) do
        _MixedReadBucket(etype)
    end
    for token in node_tokens
        external=_msh_size_t(token,"element node tag")
        internal=get(acc.node_map,external,Int32(0))
        internal!=0 || throw(ArgumentError(
            "read_mixed_msh: element references unknown node tag $external"))
        push!(bucket.nodes,internal)
    end
    _append_mixed_cell_metadata!(bucket,layout,physical,entity,element_tag,
                                 parent_tag,domain_tags)
    return nothing
end

function _push_mixed_cell_binary!(acc,etype::Int,physical::Int32,entity::Int32,
                                  element_tag::UInt64,external_nodes,limits;
                                  parent_tag::UInt64=UInt64(0),
                                  reserve_connectivity::Bool=true)
    layout=_msh_element_layout(etype,"binary element")
    layout.variable && throw(ArgumentError(
        "read_mixed_msh: binary MSH 2.2 has no record width for variable-connectivity type $etype"))
    length(external_nodes)==layout.nnodes || throw(ArgumentError(
        "read_mixed_msh: type $etype expects $(layout.nnodes) node tags"))
    reserve_connectivity && _reserve_mixed_connectivity!(
        acc,limits,length(external_nodes),"binary element record")
    bucket=get!(acc.buckets,etype) do
        _MixedReadBucket(etype)
    end
    @inbounds for raw in external_nodes
        raw>0 || throw(ArgumentError(
            "read_mixed_msh: element node tags must be positive"))
        external=UInt64(raw)
        internal=get(acc.node_map,external,Int32(0))
        internal!=0 || throw(ArgumentError(
            "read_mixed_msh: element references unknown node tag $external"))
        push!(bucket.nodes,internal)
    end
    _append_mixed_cell_metadata!(bucket,layout,physical,entity,element_tag,
                                 parent_tag,(UInt64(0),UInt64(0)))
    return nothing
end

function _read_mixed_elements_v2!(acc,io,limits)
    remaining=limits.max_elements-length(acc.element_tags)
    count=_section_count(io,"v2 element",remaining)
    for _ in 1:count
        reader=_MshRecordReader(io)
        element_tag=_msh_record_size_t!(reader,"v2 element tag")
        _record_mixed_element_tag!(acc,element_tag)
        etype=_msh_record_int!(reader,"v2 element type")
        layout=_msh_element_layout(etype,"MSH 2.2")
        ntags=_msh_record_int!(reader,"v2 element tag count")
        ntags>=0 || throw(ArgumentError(
            "read_mixed_msh: negative v2 element tag count"))
        physical=Int32(0)
        npart=0
        base_tags=ntags
        first_link=0
        second_link=0
        for i in 1:ntags
            value=_msh_record_int!(reader,"v2 element metadata tag")
            i==1 && (physical=_mixed_physical_tag(value,"physical tag"))
            if i==3 && ntags>3
                npart=value
                npart>0 || throw(ArgumentError(
                    "read_mixed_msh: v2 partition count must be positive when present"))
                base_tags=try Base.checked_add(3,npart) catch err
                    err isa InterruptException && rethrow()
                    throw(ArgumentError(
                        "read_mixed_msh: v2 partition metadata count overflows Int"))
                end
                base_tags<=ntags || throw(ArgumentError(
                    "read_mixed_msh: v2 partition metadata is truncated"))
            end
            i==base_tags+1 && (first_link=value)
            i==base_tags+2 && (second_link=value)
        end
        parent_tag=UInt64(0); domain_tags=(UInt64(0),UInt64(0))
        if layout.special
            ntags==3 && throw(ArgumentError(
                "read_mixed_msh: ambiguous third metadata tag on special type $etype"))
            extra_tags=ntags-base_tags
            allowed=layout.links===:parent ? (extra_tags in (0,1)) :
                                            (extra_tags in (0,2))
            allowed || throw(ArgumentError(
                "read_mixed_msh: unsupported special-element metadata layout for type $etype"))
            if layout.links===:parent && extra_tags==1
                first_link>0 || throw(ArgumentError(
                    "read_mixed_msh: special-element parent tag must be positive"))
                parent_tag=UInt64(first_link)
            elseif layout.links===:domains && extra_tags==2
                first_link>0 || throw(ArgumentError(
                    "read_mixed_msh: first special-element domain tag must be positive"))
                second_link>=0 || throw(ArgumentError(
                    "read_mixed_msh: second special-element domain tag must be non-negative"))
                domain_tags=(UInt64(first_link),UInt64(second_link))
            end
        end
        if layout.variable
            arity=_msh_record_int!(reader,"v2 variable-connectivity count")
            arity>0 && arity%layout.unit==0 || throw(ArgumentError(
                "read_mixed_msh: type $etype connectivity count must be a positive multiple of $(layout.unit)"))
        else
            arity=layout.nnodes
        end
        _reserve_mixed_connectivity!(acc,limits,arity,"v2 element record")
        bucket=get!(acc.buckets,etype) do
            _MixedReadBucket(etype)
        end
        for _ in 1:arity
            external=_msh_record_size_t!(reader,"v2 element node tag")
            internal=get(acc.node_map,external,Int32(0))
            internal!=0 || throw(ArgumentError(
                "read_mixed_msh: element references unknown node tag $external"))
            push!(bucket.nodes,internal)
        end
        _expect_msh_record_end!(reader,"v2 element record")
        _append_mixed_cell_metadata!(bucket,layout,physical,Int32(0),element_tag,
                                     parent_tag,domain_tags)
    end
    _expect_msh_end(io,"\$EndElements")
end

function _read_mixed_elements_v2_binary!(acc,io,limits,swap::Bool)
    remaining=limits.max_elements-length(acc.element_tags)
    count=_section_count(io,"v2 element",remaining)
    sizehint!(acc.element_tags,length(acc.element_tags)+count)
    nread=0
    while nread<count
        etype=Int(_binary_i32(io,swap,"v2 binary element-block type"))
        nlocal=Int(_binary_i32(io,swap,"v2 binary element-block count"))
        ntags=Int(_binary_i32(io,swap,"v2 binary element tag count"))
        layout=_msh_element_layout(etype,"binary MSH 2.2")
        layout.variable && throw(ArgumentError(
            "read_mixed_msh: binary MSH 2.2 has no record width for variable-connectivity type $etype"))
        0<nlocal<=count-nread || throw(ArgumentError(
            "read_mixed_msh: invalid or excessive v2 binary element-block size"))
        ntags>=0 || throw(ArgumentError(
            "read_mixed_msh: negative v2 element tag count"))
        width=try
            Base.checked_add(1,Base.checked_add(ntags,layout.nnodes))
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "read_mixed_msh: v2 binary element record width overflows Int"))
        end
        words=try Base.checked_mul(nlocal,width) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "read_mixed_msh: v2 binary element-block size overflows Int"))
        end
        connectivity_count=try Base.checked_mul(nlocal,layout.nnodes) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "read_mixed_msh: v2 binary connectivity count overflows Int"))
        end
        _reserve_mixed_connectivity!(acc,limits,connectivity_count,
                                     "v2 binary element block")
        data=_binary_i32_vector(io,words,swap,"v2 binary element block")
        offset=1
        @inbounds for _ in 1:nlocal
            signed_element_tag=data[offset]
            signed_element_tag>0 || throw(ArgumentError(
                "read_mixed_msh: element tags must be positive"))
            element_tag=UInt64(signed_element_tag)
            _record_mixed_element_tag!(acc,element_tag)
            physical=Int32(0)
            if ntags>0
                physical=_mixed_physical_tag(
                    Int(data[offset+1]),"physical tag")
            end
            parent_tag=UInt64(0)
            if layout.special
                npart=0; base_tags=min(ntags,2)
                if ntags>3
                    npart=Int(data[offset+3])
                    npart>0 || throw(ArgumentError(
                        "read_mixed_msh: v2 partition count must be positive when present"))
                    base_tags=try Base.checked_add(3,npart) catch err
                        err isa InterruptException && rethrow()
                        throw(ArgumentError(
                            "read_mixed_msh: v2 partition metadata count overflows Int"))
                    end
                    base_tags<=ntags || throw(ArgumentError(
                        "read_mixed_msh: v2 partition metadata is truncated"))
                end
                extra_tags=ntags-base_tags
                allowed=layout.links===:parent ? (extra_tags in (0,1)) :
                                                extra_tags==0
                allowed || throw(ArgumentError(
                    "read_mixed_msh: unsupported binary special-element metadata layout for type $etype"))
                if layout.links===:parent && extra_tags==1
                    raw=data[offset+ntags]
                    raw>0 || throw(ArgumentError(
                        "read_mixed_msh: special-element parent tag must be positive"))
                    parent_tag=UInt64(raw)
                end
            end
            first_node=offset+1+ntags
            last_node=first_node+layout.nnodes-1
            _push_mixed_cell_binary!(
                acc,etype,physical,Int32(0),element_tag,
                @view(data[first_node:last_node]),limits;
                parent_tag=parent_tag,reserve_connectivity=false)
            offset+=width; nread+=1
        end
    end
    nread==count || throw(ArgumentError(
        "read_mixed_msh: v2 Elements declared $count elements but contained $nread"))
    _consume_binary_newline(io,"Elements section")
    _expect_msh_end(io,"\$EndElements")
    return nothing
end

function _read_mixed_elements_v4!(acc,io,limits)
    reader=_MshTokenReader(io)
    nblocks=_msh_int(
        _msh_token!(reader,"v4 Elements header"),"v4 element-block count")
    count=_msh_int(
        _msh_token!(reader,"v4 Elements header"),"v4 element count")
    declared_min=_msh_size_t(
        _msh_token!(reader,"v4 Elements header"),"v4 element tag range")
    declared_max=_msh_size_t(
        _msh_token!(reader,"v4 Elements header"),"v4 element tag range")
    nblocks>=0 || throw(ArgumentError(
        "read_mixed_msh: negative v4 element-block count"))
    count>=0 || throw(ArgumentError("read_mixed_msh: negative v4 element count"))
    cumulative_blocks=try Base.checked_add(acc.element_blocks,nblocks) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "read_mixed_msh: cumulative v4 element-block count overflows Int"))
    end
    cumulative_blocks<=limits.max_blocks || throw(ArgumentError(
        "read_mixed_msh: cumulative v4 element-block count exceeds max_blocks"))
    count<=limits.max_elements-length(acc.element_tags) || throw(ArgumentError(
        "read_mixed_msh: cumulative v4 element count exceeds max_elements"))
    if count==0
        (declared_min==0&&declared_max==0) || throw(ArgumentError(
            "read_mixed_msh: empty v4 Elements must declare range 0 0"))
    else
        0<declared_min<=declared_max || throw(ArgumentError(
            "read_mixed_msh: invalid v4 element-tag range"))
    end
    actual_min=typemax(UInt64); actual_max=zero(UInt64); nread=0
    for _ in 1:nblocks
        dim=_msh_int(
            _msh_token!(reader,"v4 element-block header"),
            "element-block dimension")
        entity=_msh_int(
            _msh_token!(reader,"v4 element-block header"),
            "element-block entity")
        etype=_msh_int(
            _msh_token!(reader,"v4 element-block header"),
            "element-block type")
        nlocal=_msh_int(
            _msh_token!(reader,"v4 element-block header"),
            "element-block count")
        0<=dim<=3 || throw(ArgumentError(
            "read_mixed_msh: element-block dimension $dim is outside 0:3"))
        1<=entity<=typemax(Int32) || throw(ArgumentError(
            "read_mixed_msh: element-block entity tags must be positive and fit Int32"))
        haskey(acc.entity_physical,(dim,entity)) || throw(ArgumentError(
            "read_mixed_msh: element block references undeclared entity ($dim,$entity)"))
        layout=_msh_element_layout(etype,"MSH 4.1")
        layout.variable && throw(ArgumentError(
            "read_mixed_msh: MSH 4.1 has no per-element connectivity count for variable-connectivity type $etype"))
        layout.dim==dim || throw(ArgumentError(
            "read_mixed_msh: type $etype is incompatible with entity dimension $dim"))
        0<=nlocal<=count-nread || throw(ArgumentError(
            "read_mixed_msh: invalid or excessive v4 element-block size"))
        physical=acc.entity_physical[(dim,entity)]
        for _ in 1:nlocal
            tag=_msh_size_t(
                _msh_token!(reader,"v4 element tag"),"v4 element tag")
            _record_mixed_element_tag!(acc,tag)
            actual_min=min(actual_min,tag); actual_max=max(actual_max,tag)
            bucket=get!(acc.buckets,etype) do
                _MixedReadBucket(etype)
            end
            _reserve_mixed_connectivity!(acc,limits,layout.nnodes,
                                         "v4 element record")
            for _ in 1:layout.nnodes
                external=_msh_size_t(
                    _msh_token!(reader,"v4 element node tag"),
                    "v4 element node tag")
                internal=get(acc.node_map,external,Int32(0))
                internal!=0 || throw(ArgumentError(
                    "read_mixed_msh: element references unknown node tag $external"))
                push!(bucket.nodes,internal)
            end
            _append_mixed_cell_metadata!(bucket,layout,physical,Int32(entity),
                                         tag,UInt64(0),(UInt64(0),UInt64(0)))
            nread+=1
        end
    end
    nread==count || throw(ArgumentError(
        "read_mixed_msh: v4 Elements declared $count elements but contained $nread"))
    count==0 || (actual_min==declared_min&&actual_max==declared_max) || throw(ArgumentError(
        "read_mixed_msh: v4 element-tag range does not match records"))
    _expect_msh_token_end!(reader,"\$EndElements")
    acc.element_blocks=cumulative_blocks
end

function _read_mixed_elements_v4_binary!(acc,io,limits,swap::Bool)
    nblocks=_binary_count(
        _binary_u64(io,swap,"v4 binary Elements header"),
        limits.max_blocks,"v4 element-block")
    remaining_elements=limits.max_elements-length(acc.element_tags)
    count=_binary_count(
        _binary_u64(io,swap,"v4 binary Elements header"),
        remaining_elements,"v4 element")
    declared_min=_binary_u64(io,swap,"v4 binary element tag range")
    declared_max=_binary_u64(io,swap,"v4 binary element tag range")
    cumulative_blocks=try Base.checked_add(acc.element_blocks,nblocks) catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError(
            "read_mixed_msh: cumulative v4 element-block count overflows Int"))
    end
    cumulative_blocks<=limits.max_blocks || throw(ArgumentError(
        "read_mixed_msh: cumulative v4 element-block count exceeds max_blocks"))
    if count==0
        (declared_min==0&&declared_max==0) || throw(ArgumentError(
            "read_mixed_msh: empty v4 Elements must declare range 0 0"))
    else
        0<declared_min<=declared_max || throw(ArgumentError(
            "read_mixed_msh: invalid v4 element-tag range"))
    end
    sizehint!(acc.element_tags,length(acc.element_tags)+count)
    actual_min=typemax(UInt64); actual_max=zero(UInt64); nread=0
    for _ in 1:nblocks
        dim=Int(_binary_i32(io,swap,"v4 binary element-block dimension"))
        entity=Int(_binary_i32(io,swap,"v4 binary element-block entity"))
        etype=Int(_binary_i32(io,swap,"v4 binary element-block type"))
        nlocal=_binary_count(
            _binary_u64(io,swap,"v4 binary element-block count"),
            count-nread,"v4 element-block")
        0<=dim<=3 || throw(ArgumentError(
            "read_mixed_msh: element-block dimension $dim is outside 0:3"))
        1<=entity<=typemax(Int32) || throw(ArgumentError(
            "read_mixed_msh: element-block entity tags must be positive and fit Int32"))
        haskey(acc.entity_physical,(dim,entity)) || throw(ArgumentError(
            "read_mixed_msh: element block references undeclared entity ($dim,$entity)"))
        layout=_msh_element_layout(etype,"binary MSH 4.1")
        layout.variable && throw(ArgumentError(
            "read_mixed_msh: binary MSH 4.1 has no record width for variable-connectivity type $etype"))
        layout.dim==dim || throw(ArgumentError(
            "read_mixed_msh: type $etype is incompatible with entity dimension $dim"))
        nlocal==0 && continue
        width=1+layout.nnodes
        words=try Base.checked_mul(nlocal,width) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "read_mixed_msh: v4 element-block size overflows Int"))
        end
        connectivity_count=try Base.checked_mul(nlocal,layout.nnodes) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "read_mixed_msh: v4 connectivity count overflows Int"))
        end
        _reserve_mixed_connectivity!(acc,limits,connectivity_count,
                                     "v4 binary element block")
        data=_binary_u64_vector(io,words,swap,"v4 element block")
        physical=acc.entity_physical[(dim,entity)]
        bucket=get!(acc.buckets,etype) do
            _MixedReadBucket(etype)
        end
        cell_target=try Base.checked_add(length(bucket.tags),nlocal) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "read_mixed_msh: v4 element bucket count overflows Int"))
        end
        node_target=try Base.checked_mul(cell_target,layout.nnodes) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError(
                "read_mixed_msh: v4 connectivity count overflows Int"))
        end
        sizehint!(bucket.nodes,node_target); sizehint!(bucket.tags,cell_target)
        sizehint!(bucket.entities,cell_target); sizehint!(bucket.external_tags,cell_target)
        offset=1
        @inbounds for _ in 1:nlocal
            tag=data[offset]
            _record_mixed_element_tag!(acc,tag)
            actual_min=min(actual_min,tag); actual_max=max(actual_max,tag)
            for i in 1:layout.nnodes
                external=data[offset+i]
                external>0 || throw(ArgumentError(
                    "read_mixed_msh: element node tags must be positive"))
                internal=get(acc.node_map,external,Int32(0))
                internal!=0 || throw(ArgumentError(
                    "read_mixed_msh: element references unknown node tag $external"))
                push!(bucket.nodes,internal)
            end
            _append_mixed_cell_metadata!(bucket,layout,physical,Int32(entity),
                                         tag,UInt64(0),(UInt64(0),UInt64(0)))
            offset+=width; nread+=1
        end
    end
    nread==count || throw(ArgumentError(
        "read_mixed_msh: v4 Elements declared $count elements but contained $nread"))
    count==0 || (actual_min==declared_min&&actual_max==declared_max) || throw(ArgumentError(
        "read_mixed_msh: v4 element-tag range does not match records"))
    _consume_binary_newline(io,"Elements section")
    _expect_msh_end(io,"\$EndElements")
    acc.element_blocks=cumulative_blocks
    return nothing
end

function _finish_mixed_read(acc,is_v4::Bool)
    nn=length(acc.x); coords=Matrix{Float64}(undef,3,nn)
    @inbounds for i in 1:nn
        coords[1,i]=acc.x[i]; coords[2,i]=acc.y[i]; coords[3,i]=acc.z[i]
    end
    blocks=MixedElementBlock[]
    block_entities=Vector{Int32}[]
    external_element_tags=Vector{UInt64}[]
    ordered=sort!(collect(acc.buckets);by=first)
    tag_to_ref=Dict{UInt64,ElementRef}()
    sizehint!(tag_to_ref,length(acc.element_tags))
    for (bi,(_,bucket)) in pairs(ordered)
        for (cell,tag) in pairs(bucket.external_tags)
            tag_to_ref[tag]=ElementRef(bi,cell)
        end
    end
    resolve_ref(tag::UInt64,context::AbstractString)=begin
        tag==0 && return ElementRef()
        ref=get(tag_to_ref,tag,nothing)
        ref===nothing && throw(ArgumentError(
            "read_mixed_msh: $context references unknown element tag $tag"))
        return ref
    end
    for (etype,bucket) in ordered
        layout=_msh_element_layout(etype,"element"); count=length(bucket.tags)
        if layout.special
            if layout.variable
                length(bucket.offsets)==count+1 || throw(ErrorException(
                    "read_mixed_msh: internal variable-connectivity offset mismatch"))
                offsets=bucket.offsets
            else
                expected=try Base.checked_mul(layout.nnodes,count) catch err
                    err isa InterruptException && rethrow()
                    throw(ErrorException(
                        "read_mixed_msh: internal connectivity count overflows Int"))
                end
                length(bucket.nodes)==expected || throw(ErrorException(
                    "read_mixed_msh: internal connectivity accumulation mismatch"))
                offsets=Vector{Int32}(undef,count+1)
                @inbounds for i in 0:count
                    offsets[i+1]=Int32(i*layout.nnodes+1)
                end
            end
            parents=fill(ElementRef(),count)
            domains=fill(ElementRef(),2,count)
            if layout.links===:parent
                length(bucket.parent_tags)==count || throw(ErrorException(
                    "read_mixed_msh: internal parent-link accumulation mismatch"))
                @inbounds for j in 1:count
                    parents[j]=resolve_ref(bucket.parent_tags[j],
                                           "special-element parent")
                end
            else
                length(bucket.domain_tags)==count || throw(ErrorException(
                    "read_mixed_msh: internal domain-link accumulation mismatch"))
                @inbounds for j in 1:count
                    domains[1,j]=resolve_ref(bucket.domain_tags[j][1],
                                             "special-element first domain")
                    domains[2,j]=resolve_ref(bucket.domain_tags[j][2],
                                             "special-element second domain")
                end
            end
            push!(blocks,SpecialElementBlock(
                etype,bucket.nodes,offsets,bucket.tags;
                parent_refs=parents,domain_refs=domains))
        else
            expected=try Base.checked_mul(layout.nnodes,count) catch err
                err isa InterruptException && rethrow()
                throw(ErrorException(
                    "read_mixed_msh: internal connectivity count overflows Int"))
            end
            length(bucket.nodes)==expected || throw(ErrorException(
                "read_mixed_msh: internal connectivity accumulation mismatch"))
            push!(blocks,ElementBlock(
                etype,reshape(bucket.nodes,layout.nnodes,count),bucket.tags))
        end
        push!(block_entities,copy(bucket.entities))
        push!(external_element_tags,copy(bucket.external_tags))
    end
    data=is_v4 ? MixedEntityData(acc.entities;
        node_entities=acc.node_entities,
        node_parametric=acc.node_parametric,
        external_node_tags=acc.external_node_tags,
        block_entities=block_entities,
        external_element_tags=external_element_tags) : nothing
    return MixedMesh(
        _OWNED_MIXED_MESH,coords,blocks,acc.physical_names,data)
end

function _escape_mixed_name(name::AbstractString)
    out=IOBuffer()
    for c in name
        c=='\\' ? write(out,"\\\\") :
        c=='\"' ? write(out,"\\\"") :
        c=='\n' ? write(out,"\\n") :
        c=='\r' ? write(out,"\\r") :
        c=='\t' ? write(out,"\\t") : write(out,c)
    end
    return String(take!(out))
end

function _write_mixed_physical_names(io,names,gmsh_compatible::Bool)
    isempty(names) && return nothing
    println(io,"\$PhysicalNames"); println(io,length(names))
    for ((dim,tag),name) in sort!(collect(names);by=first)
        encoded=gmsh_compatible ? name : _escape_mixed_name(name)
        ncodeunits(encoded)<=MSH_PHYSICAL_NAME_MAX_BYTES || throw(ArgumentError(
            "write_mixed_msh: serialized physical name for ($dim,$tag) exceeds " *
            "$MSH_PHYSICAL_NAME_MAX_BYTES bytes"))
        println(io,dim," ",tag," \"",encoded,"\"")
    end
    println(io,"\$EndPhysicalNames")
    return nothing
end

"""
    write_mixed_msh(path, mesh; version=4.1, binary=false,
                    gmsh_compatible=true) -> path

Validate and atomically write a mixed mesh as ASCII or binary MSH v2.2 or v4.1. V4
entity metadata is preserved exactly when present; legacy meshes receive a
deterministic discrete-entity layout with coordinate-derived bounds. MSH v2.2
ASCII retains supported parent/domain links and variable decompositions; its
binary form retains parent links but has neither variable record widths nor
domain-link fields. MSH4 supports only fixed-width, unlinked special records.
Unsupported combinations are rejected before a temporary file is created.
MSH v2.2 necessarily writes the legacy first-physical-tag projection. By default,
element tags that pinned Gmsh 4.15.2 cannot safely consume through its normal
lifecycle, and quoted or line-breaking names that Gmsh cannot preserve, are
explicit blockers. Backslashes and tabs are standard literal characters. Set
`gmsh_compatible=false` only for Tessella-to-Tessella serialization of those
records. Legacy MSH4 synthesis is likewise blocked for nonzero-physical special
records: Gmsh 4.15.2 can parse such a file, but rewrites an invalid node section
unless compatible node/entity classification metadata is supplied. Read escaped
names with `read_mixed_msh(...; tessella_extensions=true)`.
"""
function write_mixed_msh(path::AbstractString,m::MixedMesh;version=4.1,
                         binary=false,gmsh_compatible=true)
    version isa Real || throw(ArgumentError(
        "write_mixed_msh: version must be real"))
    _elements_reject_bool(version,"write_mixed_msh: version")
    value = try
        Float64(version)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("write_mixed_msh: version must be representable as Float64"))
    end
    value in (2.2,4.1) || throw(ArgumentError(
        "write_mixed_msh: version must be 2.2 or 4.1"))
    gmsh_compatible isa Bool || throw(ArgumentError(
        "write_mixed_msh: gmsh_compatible must be Bool"))
    binary isa Bool || throw(ArgumentError(
        "write_mixed_msh: binary must be Bool"))
    if gmsh_compatible
        reader_gaps=value==2.2 ? GMSH_4_15_2_MSH_READER_GAPS_V2 :
                                GMSH_4_15_2_MSH_READER_GAPS_V4
        gaps=sort!(unique!(Int[b.msh for b in m.blocks
                              if b.msh in reader_gaps]))
        isempty(gaps) || throw(ArgumentError(
            "write_mixed_msh: Gmsh 4.15.2 cannot safely consume element type(s) " *
            join(gaps,",") * "; use gmsh_compatible=false only for a " *
            "Tessella-only round trip"))
        bad_names=Tuple{Int,Int}[]
        for (key,name) in m.physical_names
            any(c->c in ('\"','\n','\r'),name) && push!(bad_names,key)
        end
        isempty(bad_names) || throw(ArgumentError(
            "write_mixed_msh: Gmsh 4.15.2 cannot preserve quoted or " *
            "line-breaking physical " *
            "name(s) " * join(sort!(bad_names),",") *
            "; use gmsh_compatible=false only for a Tessella-only round trip"))
    end
    diagnostic=validate(m)
    diagnostic.ok || throw(ArgumentError(
        "write_mixed_msh: invalid mixed mesh: "*join(diagnostic.messages,"; ")))
    if gmsh_compatible && value==4.1
        data=m.entity_data
        for (bi,block) in pairs(m.blocks)
            block isa SpecialElementBlock || continue
            @inbounds for cell in 1:_block_ncells(block)
                block.tags[cell]==0 && continue
                data===nothing && throw(ArgumentError(
                    "write_mixed_msh: Gmsh 4.15.2 cannot safely rewrite a " *
                    "nonzero-physical special MSH 4.1 record synthesized " *
                    "without node/entity classification metadata; supply " *
                    "compatible MixedEntityData or use gmsh_compatible=false " *
                    "only for a Tessella-only round trip"))
                owner=(_block_dim(block),data.block_entities[bi][cell])
                for slot in 1:_cell_arity(block,cell)
                    node=_cell_node(block,cell,slot)
                    data.node_entities[node]==owner || throw(ArgumentError(
                        "write_mixed_msh: Gmsh 4.15.2 rewrite safety for a " *
                        "nonzero-physical special MSH 4.1 record requires each " *
                        "record node to be classified on its owning entity; " *
                        "supply compatible MixedEntityData or use " *
                        "gmsh_compatible=false only for a Tessella-only round trip"))
                end
            end
        end
    end
    _assert_mixed_msh_format(m,value,binary)
    names=_copy_physical_names(m.physical_names,"write_mixed_msh")
    target=abspath(path); parent=dirname(target)
    isdir(parent) || throw(ArgumentError(
        "write_mixed_msh: parent directory does not exist: $parent"))
    isdir(target) && throw(ArgumentError(
        "write_mixed_msh: destination is a directory: $target"))
    mktemp(parent) do temporary,io
        if value==2.2
            binary ? _write_mixed_v2_binary(io,m,names,gmsh_compatible) :
                     _write_mixed_v2(io,m,names,gmsh_compatible)
        else
            binary ? _write_mixed_v4_binary(io,m,names,gmsh_compatible) :
                     _write_mixed_v4(io,m,names,gmsh_compatible)
        end
        flush(io); close(io)
        mv(temporary,target;force=true)
    end
    return path
end

function _assert_mixed_msh_format(m::MixedMesh,version::Float64,binary::Bool)
    for (bi,block) in pairs(m.blocks)
        block isa SpecialElementBlock || continue
        record=MSH_SPECIAL_RECORDS[block.msh]
        if version==2.2
            binary && record.variable && throw(ArgumentError(
                "write_mixed_msh: binary MSH 2.2 has no record width for variable-connectivity type $(block.msh)"))
            if binary && record.links===:domains
                for j in 1:_block_ncells(block)
                    if !_missing_ref(block.domain_refs[1,j]) ||
                       !_missing_ref(block.domain_refs[2,j])
                        throw(ArgumentError(
                            "write_mixed_msh: binary MSH 2.2 cannot encode domain links for special type $(block.msh) in block $bi"))
                    end
                end
            end
        else
            record.variable && throw(ArgumentError(
                "write_mixed_msh: MSH 4.1 has no sound variable-connectivity record for type $(block.msh) in Gmsh 4.15.2"))
            for j in 1:_block_ncells(block)
                linked=!_missing_ref(block.parent_refs[j]) ||
                       !_missing_ref(block.domain_refs[1,j]) ||
                       !_missing_ref(block.domain_refs[2,j])
                linked && throw(ArgumentError(
                    "write_mixed_msh: MSH 4.1 cannot encode parent/domain links for special type $(block.msh) in block $bi"))
            end
        end
    end
    return nothing
end

function _mixed_v2_order(m::MixedMesh)
    starts=Vector{Int}(undef,length(m.blocks)+1); starts[1]=1
    for (bi,block) in pairs(m.blocks)
        starts[bi+1]=try Base.checked_add(starts[bi],_block_ncells(block)) catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("write_mixed_msh: element count overflows Int"))
        end
    end
    total=starts[end]-1
    refs=Vector{_MixedCellRef}(undef,total)
    indegree=zeros(Int,total)
    dependents=[Int[] for _ in 1:total]
    for (bi,block) in pairs(m.blocks), cell in 1:_block_ncells(block)
        source=starts[bi]+cell-1
        refs[source]=_MixedCellRef(bi,cell)
        block isa SpecialElementBlock || continue
        slot=1
        while true
            link=_special_link(block,cell,slot); link===nothing && break
            slot+=1; _missing_ref(link) && continue
            target=starts[Int(link.block)]+Int(link.cell)-1
            indegree[source]+=1
            push!(dependents[target],source)
        end
    end
    queue=Vector{Int}(undef,total); tail=0
    for i in 1:total
        indegree[i]==0 || continue
        tail+=1; queue[tail]=i
    end
    order=Vector{_MixedCellRef}(undef,total); head=1; written=0
    while head<=tail
        item=queue[head]; head+=1; written+=1; order[written]=refs[item]
        for dependent in dependents[item]
            indegree[dependent]-=1
            if indegree[dependent]==0
                tail+=1; queue[tail]=dependent
            end
        end
    end
    written==total || throw(ArgumentError(
        "write_mixed_msh: element parent/domain links contain a cycle"))
    output_tags=[zeros(Int32,_block_ncells(block)) for block in m.blocks]
    @inbounds for (tag,ref) in pairs(order)
        output_tags[ref.block][ref.cell]=Int32(tag)
    end
    return order,output_tags
end

@inline function _v2_link_tag(output_tags,ref::ElementRef)
    _missing_ref(ref) && return Int32(0)
    return output_tags[Int(ref.block)][Int(ref.cell)]
end

function _v2_entity_ids(m::MixedMesh)
    keys=Set{Tuple{Int,Int32}}()
    @inbounds for b in m.blocks
        dim=_block_dim(b)
        for tag in b.tags
            push!(keys,(dim,tag))
        end
    end
    counters=zeros(Int,4); ids=Dict{Tuple{Int,Int32},Int}()
    for key in sort!(collect(keys))
        dim=key[1]; counters[dim+1]+=1; ids[key]=counters[dim+1]
    end
    return ids
end

function _write_mixed_v2(io,m::MixedMesh,names,gmsh_compatible::Bool)
    println(io,"\$MeshFormat"); println(io,"2.2 0 8"); println(io,"\$EndMeshFormat")
    _write_mixed_physical_names(io,names,gmsh_compatible)
    nn=size(m.coords,2)
    println(io,"\$Nodes"); println(io,nn)
    @inbounds for i in 1:nn
        @printf(io,"%d %.17g %.17g %.17g\n",i,m.coords[1,i],m.coords[2,i],m.coords[3,i])
    end
    println(io,"\$EndNodes")
    nel=_assert_mixed_structure(m,"write_mixed_msh")
    println(io,"\$Elements"); println(io,nel)
    ids=_v2_entity_ids(m); order,output_tags=_mixed_v2_order(m)
    @inbounds for ref in order
        b=m.blocks[ref.block]; j=ref.cell; physical=b.tags[j]
        entity=ids[(_block_dim(b),physical)]
        print(io,output_tags[ref.block][j]," ",b.msh)
        if b isa SpecialElementBlock
            record=MSH_SPECIAL_RECORDS[b.msh]
            if record.links===:parent && !_missing_ref(b.parent_refs[j])
                print(io," 5 ",physical," ",entity," 1 0 ",
                      _v2_link_tag(output_tags,b.parent_refs[j]))
            elseif record.links===:domains && !_missing_ref(b.domain_refs[1,j])
                print(io," 6 ",physical," ",entity," 1 0 ",
                      _v2_link_tag(output_tags,b.domain_refs[1,j])," ",
                      _v2_link_tag(output_tags,b.domain_refs[2,j]))
            else
                print(io," 2 ",physical," ",entity)
            end
            record.variable && print(io," ",_cell_arity(b,j))
        else
            print(io," 2 ",physical," ",entity)
        end
        for i in 1:_cell_arity(b,j)
            print(io," ",_cell_node(b,j,i))
        end
        println(io)
    end
    println(io,"\$EndElements")
    return nothing
end

function _write_mixed_binary_format(io,version::AbstractString)
    println(io,"\$MeshFormat")
    println(io,version," 1 8")
    write(io,Int32(1)); write(io,UInt8('\n'))
    println(io,"\$EndMeshFormat")
    return nothing
end

function _write_mixed_v2_binary(io,m::MixedMesh,names,gmsh_compatible::Bool)
    _write_mixed_binary_format(io,"2.2")
    _write_mixed_physical_names(io,names,gmsh_compatible)
    nn=size(m.coords,2)
    println(io,"\$Nodes"); println(io,nn)
    @inbounds for i in 1:nn
        write(io,Int32(i))
        write(io,m.coords[1,i]); write(io,m.coords[2,i]); write(io,m.coords[3,i])
    end
    write(io,UInt8('\n')); println(io,"\$EndNodes")
    nel=_assert_mixed_structure(m,"write_mixed_msh")
    println(io,"\$Elements"); println(io,nel)
    ids=_v2_entity_ids(m); order,output_tags=_mixed_v2_order(m)
    position=1
    @inbounds while position<=length(order)
        first_ref=order[position]; first_block=m.blocks[first_ref.block]
        parent_link=first_block isa SpecialElementBlock &&
                    MSH_SPECIAL_RECORDS[first_block.msh].links===:parent &&
                    !_missing_ref(first_block.parent_refs[first_ref.cell])
        ntags=parent_link ? 3 : 2
        last=position
        while last<length(order)
            next_ref=order[last+1]; next_block=m.blocks[next_ref.block]
            next_parent=next_block isa SpecialElementBlock &&
                        MSH_SPECIAL_RECORDS[next_block.msh].links===:parent &&
                        !_missing_ref(next_block.parent_refs[next_ref.cell])
            next_block.msh==first_block.msh && (next_parent ? 3 : 2)==ntags || break
            last+=1
        end
        write(io,Int32(first_block.msh)); write(io,Int32(last-position+1))
        write(io,Int32(ntags))
        for index in position:last
            ref=order[index]; block=m.blocks[ref.block]; j=ref.cell
            write(io,output_tags[ref.block][j]); write(io,block.tags[j])
            write(io,Int32(ids[(_block_dim(block),block.tags[j])]))
            if ntags==3
                write(io,_v2_link_tag(output_tags,block.parent_refs[j]))
            end
            for i in 1:_cell_arity(block,j)
                write(io,_cell_node(block,j,i))
            end
        end
        position=last+1
    end
    length(order)==nel || throw(ErrorException(
        "write_mixed_msh: internal v2 binary element count mismatch"))
    write(io,UInt8('\n')); println(io,"\$EndElements")
    return nothing
end

struct _MixedWriteGroup
    dim::Int
    entity::Int
    msh::Int
    cells::Vector{_MixedCellRef}
end

struct _MixedWriteEntity
    dim::Int
    tag::Int
    physical::Int32
    bounds::NTuple{6,Float64}
end

@inline function _grow_bounds!(bounds::Vector{Float64},m::MixedMesh,
                               b::MixedElementBlock,j::Int)
    @inbounds for i in 1:_cell_arity(b,j)
        n=_cell_node(b,j,i); x=m.coords[1,n]; y=m.coords[2,n]; z=m.coords[3,n]
        bounds[1]=min(bounds[1],x); bounds[2]=min(bounds[2],y); bounds[3]=min(bounds[3],z)
        bounds[4]=max(bounds[4],x); bounds[5]=max(bounds[5],y); bounds[6]=max(bounds[6],z)
    end
    return nothing
end

function _all_mixed_bounds(m::MixedMesh)
    size(m.coords,2)==0 && return (0.0,0.0,0.0,0.0,0.0,0.0)
    bounds=Float64[Inf,Inf,Inf,-Inf,-Inf,-Inf]
    @inbounds for i in axes(m.coords,2)
        x=m.coords[1,i]; y=m.coords[2,i]; z=m.coords[3,i]
        bounds[1]=min(bounds[1],x); bounds[2]=min(bounds[2],y); bounds[3]=min(bounds[3],z)
        bounds[4]=max(bounds[4],x); bounds[5]=max(bounds[5],y); bounds[6]=max(bounds[6],z)
    end
    return Tuple(bounds)
end

function _mixed_v4_layout(m::MixedMesh)
    raw_groups=Dict{Tuple{Int,Int32,Int},Vector{_MixedCellRef}}()
    entity_bounds=Dict{Tuple{Int,Int32},Vector{Float64}}()
    point_cells=_MixedCellRef[]
    @inbounds for (bi,b) in pairs(m.blocks)
        dim=_block_dim(b)
        for j in 1:_block_ncells(b)
            ref=_MixedCellRef(bi,j)
            if dim==0
                push!(point_cells,ref)
            else
                key=(dim,b.tags[j],b.msh)
                cells=get!(raw_groups,key) do
                    _MixedCellRef[]
                end
                push!(cells,ref)
                entity_key=(dim,b.tags[j])
                bounds=get!(entity_bounds,entity_key) do
                    Float64[Inf,Inf,Inf,-Inf,-Inf,-Inf]
                end
                _grow_bounds!(bounds,m,b,j)
            end
        end
    end
    sort!(point_cells;lt=(a,b)->_mixed_cell_lt(m,a,b),alg=MergeSort)
    entity_ids=Dict{Tuple{Int,Int32},Int}(); counters=zeros(Int,4)
    for key in sort!(collect(keys(entity_bounds)))
        dim=key[1]; counters[dim+1]+=1; entity_ids[key]=counters[dim+1]
    end
    entities=_MixedWriteEntity[]; groups=_MixedWriteGroup[]
    for ref in point_cells
        counters[1]+=1; entity=counters[1]
        b=m.blocks[ref.block]; n=_cell_node(b,ref.cell,1)
        x=m.coords[1,n]; y=m.coords[2,n]; z=m.coords[3,n]
        push!(entities,_MixedWriteEntity(0,entity,b.tags[ref.cell],(x,y,z,x,y,z)))
        push!(groups,_MixedWriteGroup(0,entity,b.msh,_MixedCellRef[ref]))
    end
    for key in sort!(collect(keys(entity_bounds)))
        dim,physical=key; entity=entity_ids[key]
        push!(entities,_MixedWriteEntity(dim,entity,physical,Tuple(entity_bounds[key])))
    end
    for (key,cells) in sort!(collect(raw_groups);by=first)
        dim,physical,msh=key
        sort!(cells;lt=(a,b)->_mixed_cell_lt(m,a,b),alg=MergeSort)
        push!(groups,_MixedWriteGroup(dim,entity_ids[(dim,physical)],msh,cells))
    end
    node_owner=(3,0)
    if size(m.coords,2)>0
        if isempty(groups)
            counters[4]+=1
            node_owner=(3,counters[4])
            push!(entities,_MixedWriteEntity(
                3,node_owner[2],Int32(0),_all_mixed_bounds(m)))
        else
            highest=maximum(group.dim for group in groups)
            candidates=filter(group->group.dim==highest,groups)
            sort!(candidates;by=group->(group.entity,group.msh))
            owner=first(candidates)
            node_owner=(owner.dim,owner.entity)
        end
    end
    sort!(entities;by=e->(e.dim,e.tag))
    sort!(groups;by=g->(g.dim,g.entity,g.msh))
    return entities,groups,node_owner
end

function _write_mixed_entity(io,entity::_MixedWriteEntity)
    if entity.dim==0
        @printf(io,"%d %.17g %.17g %.17g %d",entity.tag,
                entity.bounds[1],entity.bounds[2],entity.bounds[3],
                entity.physical==0 ? 0 : 1)
        entity.physical==0 || print(io," ",entity.physical)
        println(io)
    else
        @printf(io,"%d %.17g %.17g %.17g %.17g %.17g %.17g %d",
                entity.tag,entity.bounds...,entity.physical==0 ? 0 : 1)
        entity.physical==0 || print(io," ",entity.physical)
        println(io," 0")
    end
    return nothing
end

function _write_mixed_entity(io,entity::MixedEntity)
    nphysical=length(entity.physical_tags)
    if entity.dim==0
        @printf(io,"%d %.17g %.17g %.17g %d",entity.tag,
                entity.bbox[1],entity.bbox[2],entity.bbox[3],nphysical)
        for physical in entity.physical_tags
            print(io," ",physical)
        end
        println(io)
    else
        @printf(io,"%d %.17g %.17g %.17g %.17g %.17g %.17g %d",
                entity.tag,entity.bbox...,nphysical)
        for physical in entity.physical_tags
            print(io," ",physical)
        end
        print(io," ",length(entity.boundaries))
        for boundary in entity.boundaries
            print(io," ",boundary)
        end
        println(io)
    end
    return nothing
end

function _write_mixed_entity_binary(io,entity::_MixedWriteEntity)
    write(io,Int32(entity.tag))
    ncoordinates=entity.dim==0 ? 3 : 6
    @inbounds for i in 1:ncoordinates
        write(io,entity.bounds[i])
    end
    nphysical=entity.physical==0 ? 0 : 1
    write(io,UInt64(nphysical))
    nphysical==0 || write(io,entity.physical)
    entity.dim==0 || write(io,UInt64(0))
    return nothing
end

function _write_mixed_entity_binary(io,entity::MixedEntity)
    write(io,Int32(entity.tag))
    ncoordinates=entity.dim==0 ? 3 : 6
    @inbounds for i in 1:ncoordinates
        write(io,entity.bbox[i])
    end
    write(io,UInt64(length(entity.physical_tags)))
    for physical in entity.physical_tags
        write(io,physical)
    end
    if entity.dim>0
        write(io,UInt64(length(entity.boundaries)))
        for boundary in entity.boundaries
            write(io,boundary)
        end
    end
    return nothing
end

function _mixed_metadata_node_runs(data::MixedEntityData)
    runs=UnitRange{Int}[]
    n=length(data.external_node_tags); first_node=1
    while first_node<=n
        dim,entity=data.node_entities[first_node]
        parametric=data.node_parametric[first_node]!==nothing
        last_node=first_node
        while last_node<n
            next_dim,next_entity=data.node_entities[last_node+1]
            next_parametric=data.node_parametric[last_node+1]!==nothing
            (next_dim==dim && next_entity==entity &&
             next_parametric==parametric) || break
            last_node+=1
        end
        push!(runs,first_node:last_node); first_node=last_node+1
    end
    return runs
end

struct _MixedMetadataElementRun
    block::Int
    cells::UnitRange{Int}
end

function _mixed_metadata_element_runs(m::MixedMesh,data::MixedEntityData)
    runs=_MixedMetadataElementRun[]
    for (bi,block) in pairs(m.blocks)
        entities=data.block_entities[bi]; first_cell=1; ncells=_block_ncells(block)
        while first_cell<=ncells
            entity=entities[first_cell]; last_cell=first_cell
            while last_cell<ncells && entities[last_cell+1]==entity
                last_cell+=1
            end
            push!(runs,_MixedMetadataElementRun(bi,first_cell:last_cell))
            first_cell=last_cell+1
        end
    end
    return runs
end

function _write_mixed_v4_binary_metadata(
    io,m::MixedMesh,names,data::MixedEntityData,gmsh_compatible::Bool)
    _write_mixed_binary_format(io,"4.1")
    _write_mixed_physical_names(io,names,gmsh_compatible)
    entities=sort!(collect(values(data.entities));by=e->(e.dim,e.tag))
    counts=zeros(Int,4)
    for entity in entities
        counts[entity.dim+1]+=1
    end
    println(io,"\$Entities")
    for count in counts
        write(io,UInt64(count))
    end
    for entity in entities
        _write_mixed_entity_binary(io,entity)
    end
    write(io,UInt8('\n')); println(io,"\$EndEntities")

    nn=size(m.coords,2); node_runs=_mixed_metadata_node_runs(data)
    minimum_node=nn==0 ? UInt64(0) : minimum(data.external_node_tags)
    maximum_node=nn==0 ? UInt64(0) : maximum(data.external_node_tags)
    println(io,"\$Nodes")
    write(io,UInt64(length(node_runs))); write(io,UInt64(nn))
    write(io,minimum_node); write(io,maximum_node)
    @inbounds for run in node_runs
        first_node=first(run); dim,entity=data.node_entities[first_node]
        parametric=data.node_parametric[first_node]!==nothing
        write(io,Int32(dim)); write(io,entity); write(io,Int32(parametric ? 1 : 0))
        write(io,UInt64(length(run)))
        for i in run
            write(io,data.external_node_tags[i])
        end
        for i in run
            write(io,m.coords[1,i]); write(io,m.coords[2,i]); write(io,m.coords[3,i])
            parameters=data.node_parametric[i]
            if parameters!==nothing
                for value in parameters
                    write(io,value)
                end
            end
        end
    end
    write(io,UInt8('\n')); println(io,"\$EndNodes")

    nel=_assert_mixed_structure(m,"write_mixed_msh")
    element_runs=_mixed_metadata_element_runs(m,data)
    all_element_tags=Iterators.flatten(data.external_element_tags)
    minimum_element=nel==0 ? UInt64(0) : minimum(all_element_tags)
    maximum_element=nel==0 ? UInt64(0) :
                    maximum(Iterators.flatten(data.external_element_tags))
    println(io,"\$Elements")
    write(io,UInt64(length(element_runs))); write(io,UInt64(nel))
    write(io,minimum_element); write(io,maximum_element)
    @inbounds for run in element_runs
        block=m.blocks[run.block]
        entity=data.block_entities[run.block][first(run.cells)]
        write(io,Int32(_block_dim(block))); write(io,entity); write(io,Int32(block.msh))
        write(io,UInt64(length(run.cells)))
        for j in run.cells
            write(io,data.external_element_tags[run.block][j])
            for i in 1:_cell_arity(block,j)
                internal=_cell_node(block,j,i)
                write(io,data.external_node_tags[internal])
            end
        end
    end
    write(io,UInt8('\n')); println(io,"\$EndElements")
    return nothing
end

function _write_mixed_v4_binary(
    io,m::MixedMesh,names,gmsh_compatible::Bool)
    m.entity_data===nothing || return _write_mixed_v4_binary_metadata(
        io,m,names,m.entity_data,gmsh_compatible)
    _write_mixed_binary_format(io,"4.1")
    _write_mixed_physical_names(io,names,gmsh_compatible)
    entities,groups,node_owner=_mixed_v4_layout(m)
    counts=zeros(Int,4)
    for entity in entities
        counts[entity.dim+1]+=1
    end
    println(io,"\$Entities")
    for count in counts
        write(io,UInt64(count))
    end
    for entity in entities
        _write_mixed_entity_binary(io,entity)
    end
    write(io,UInt8('\n')); println(io,"\$EndEntities")

    nn=size(m.coords,2)
    println(io,"\$Nodes")
    write(io,UInt64(nn==0 ? 0 : 1)); write(io,UInt64(nn))
    write(io,UInt64(nn==0 ? 0 : 1)); write(io,UInt64(nn))
    if nn>0
        write(io,Int32(node_owner[1])); write(io,Int32(node_owner[2]))
        write(io,Int32(0))
        write(io,UInt64(nn))
        for i in 1:nn
            write(io,UInt64(i))
        end
        @inbounds for i in 1:nn
            write(io,m.coords[1,i]); write(io,m.coords[2,i]); write(io,m.coords[3,i])
        end
    end
    write(io,UInt8('\n')); println(io,"\$EndNodes")

    nel=_assert_mixed_structure(m,"write_mixed_msh")
    println(io,"\$Elements")
    write(io,UInt64(length(groups))); write(io,UInt64(nel))
    write(io,UInt64(nel==0 ? 0 : 1)); write(io,UInt64(nel))
    eid=0
    @inbounds for group in groups
        write(io,Int32(group.dim)); write(io,Int32(group.entity))
        write(io,Int32(group.msh)); write(io,UInt64(length(group.cells)))
        for ref in group.cells
            block=m.blocks[ref.block]; eid+=1; write(io,UInt64(eid))
            for i in 1:_cell_arity(block,ref.cell)
                write(io,UInt64(_cell_node(block,ref.cell,i)))
            end
        end
    end
    eid==nel || throw(ErrorException(
        "write_mixed_msh: internal v4 binary element count mismatch"))
    write(io,UInt8('\n')); println(io,"\$EndElements")
    return nothing
end

function _write_mixed_v4_metadata(
    io,m::MixedMesh,names,data::MixedEntityData,gmsh_compatible::Bool)
    println(io,"\$MeshFormat"); println(io,"4.1 0 8"); println(io,"\$EndMeshFormat")
    _write_mixed_physical_names(io,names,gmsh_compatible)
    entities=sort!(collect(values(data.entities));by=e->(e.dim,e.tag))
    counts=zeros(Int,4)
    for entity in entities
        counts[entity.dim+1]+=1
    end
    println(io,"\$Entities")
    println(io,counts[1]," ",counts[2]," ",counts[3]," ",counts[4])
    for entity in entities
        _write_mixed_entity(io,entity)
    end
    println(io,"\$EndEntities")

    nn=size(m.coords,2); node_runs=_mixed_metadata_node_runs(data)
    minimum_node=nn==0 ? 0 : minimum(data.external_node_tags)
    maximum_node=nn==0 ? 0 : maximum(data.external_node_tags)
    println(io,"\$Nodes")
    println(io,length(node_runs)," ",nn," ",minimum_node," ",maximum_node)
    @inbounds for run in node_runs
        first_node=first(run); dim,entity=data.node_entities[first_node]
        parametric=data.node_parametric[first_node]!==nothing
        println(io,dim," ",entity," ",parametric ? 1 : 0," ",length(run))
        for i in run
            println(io,data.external_node_tags[i])
        end
        for i in run
            @printf(io,"%.17g %.17g %.17g",
                    m.coords[1,i],m.coords[2,i],m.coords[3,i])
            parameters=data.node_parametric[i]
            if parameters!==nothing
                for value in parameters
                    @printf(io," %.17g",value)
                end
            end
            println(io)
        end
    end
    println(io,"\$EndNodes")

    nel=_assert_mixed_structure(m,"write_mixed_msh")
    element_runs=_mixed_metadata_element_runs(m,data)
    all_element_tags=Iterators.flatten(data.external_element_tags)
    minimum_element=nel==0 ? 0 : minimum(all_element_tags)
    maximum_element=nel==0 ? 0 : maximum(Iterators.flatten(data.external_element_tags))
    println(io,"\$Elements")
    println(io,length(element_runs)," ",nel," ",minimum_element," ",maximum_element)
    @inbounds for run in element_runs
        block=m.blocks[run.block]
        entity=data.block_entities[run.block][first(run.cells)]
        println(io,_block_dim(block)," ",entity," ",block.msh," ",length(run.cells))
        for j in run.cells
            print(io,data.external_element_tags[run.block][j])
            for i in 1:_cell_arity(block,j)
                internal=_cell_node(block,j,i)
                print(io," ",data.external_node_tags[internal])
            end
            println(io)
        end
    end
    println(io,"\$EndElements")
    return nothing
end

function _write_mixed_v4(io,m::MixedMesh,names,gmsh_compatible::Bool)
    m.entity_data===nothing || return _write_mixed_v4_metadata(
        io,m,names,m.entity_data,gmsh_compatible)
    println(io,"\$MeshFormat"); println(io,"4.1 0 8"); println(io,"\$EndMeshFormat")
    _write_mixed_physical_names(io,names,gmsh_compatible)
    entities,groups,node_owner=_mixed_v4_layout(m)
    counts=zeros(Int,4)
    for entity in entities
        counts[entity.dim+1]+=1
    end
    println(io,"\$Entities")
    println(io,counts[1]," ",counts[2]," ",counts[3]," ",counts[4])
    for entity in entities
        _write_mixed_entity(io,entity)
    end
    println(io,"\$EndEntities")
    nn=size(m.coords,2)
    println(io,"\$Nodes")
    println(io,nn==0 ? 0 : 1," ",nn," ",nn==0 ? 0 : 1," ",nn)
    if nn>0
        println(io,node_owner[1]," ",node_owner[2]," 0 ",nn)
        for i in 1:nn
            println(io,i)
        end
        @inbounds for i in 1:nn
            @printf(io,"%.17g %.17g %.17g\n",m.coords[1,i],m.coords[2,i],m.coords[3,i])
        end
    end
    println(io,"\$EndNodes")
    nel=_assert_mixed_structure(m,"write_mixed_msh")
    println(io,"\$Elements")
    println(io,length(groups)," ",nel," ",nel==0 ? 0 : 1," ",nel)
    eid=0
    @inbounds for group in groups
        println(io,group.dim," ",group.entity," ",group.msh," ",length(group.cells))
        for ref in group.cells
            b=m.blocks[ref.block]; eid+=1; print(io,eid)
            for i in 1:_cell_arity(b,ref.cell)
                print(io," ",_cell_node(b,ref.cell,i))
            end
            println(io)
        end
    end
    println(io,"\$EndElements")
    return nothing
end

end # module
