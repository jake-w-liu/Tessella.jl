# Nondegenerate NoNewVerts hexahedron relation. The factories are the two
# existing-corner constructions in Gmsh 4.15.2 QuadTriExtruded3D.cpp:
# createFullHexElems (4931), with prism subdivisions from QuadToTriPriPyrTet
# (4377). No template introduces a face or body vertex.

struct _ExtrudeNoNewCell
    msh::UInt8
    nodes::NTuple{8,UInt8}
end

struct _ExtrudeNoNewTemplate
    faces::NTuple{6,UInt8}
    ncells::UInt8
    cells::NTuple{6,_ExtrudeNoNewCell}
    kind::UInt8
    choice::UInt8
end

const _EXTRUDE_NONEW_EMPTY_CELL = _ExtrudeNoNewCell(0,ntuple(_->UInt8(0),8))
const _EXTRUDE_NONEW_HEX_FACES = ((1,2,6,5),(2,3,7,6),(3,4,8,7),
    (4,1,5,8),(1,2,3,4),(5,6,7,8))
const _EXTRUDE_NONEW_PRISM_FACES = ((1,3,2),(4,5,6),(1,2,5,4),
    (2,3,6,5),(3,1,4,6))
const _EXTRUDE_NONEW_CUBE = ((0.0,0.0,0.0),(1.0,0.0,0.0),
    (1.0,1.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0),
    (1.0,0.0,1.0),(1.0,1.0,1.0),(0.0,1.0,1.0))

@inline _extrude_nonew_cell(msh,a,b,c,d) =
    _ExtrudeNoNewCell(UInt8(msh),(UInt8(a),UInt8(b),UInt8(c),UInt8(d),0,0,0,0))
@inline _extrude_nonew_cell(msh,a,b,c,d,e) =
    _ExtrudeNoNewCell(UInt8(msh),(UInt8(a),UInt8(b),UInt8(c),UInt8(d),UInt8(e),0,0,0))
@inline _extrude_nonew_cell(msh,a,b,c,d,e,f) =
    _ExtrudeNoNewCell(UInt8(msh),(UInt8(a),UInt8(b),UInt8(c),UInt8(d),UInt8(e),UInt8(f),0,0))
@inline _extrude_nonew_cell(msh,a,b,c,d,e,f,g,h) =
    _ExtrudeNoNewCell(UInt8(msh),(UInt8(a),UInt8(b),UInt8(c),UInt8(d),
                              UInt8(e),UInt8(f),UInt8(g),UInt8(h)))

@inline function _extrude_nonew_swap_node(nodes,k,first_swap,second_swap)
    a,b=first_swap
    k==a && return nodes[b]
    k==b && return nodes[a]
    c,d=second_swap
    k==c && return nodes[d]
    k==d && return nodes[c]
    return nodes[k]
end

# Normalize once on the reference cube. Actual sweep geometry is certified by
# the region planner before emission; orientation normalization cannot repair
# a folded physical cell.
function _extrude_nonew_positive_cell(cell::_ExtrudeNoNewCell)
    total=0.0
    for (a,b,c,d) in _EXTRUDE_CELL_TETS[Int(cell.msh)]
        total+=tet_signed_volume(_EXTRUDE_NONEW_CUBE[cell.nodes[a]],
            _EXTRUDE_NONEW_CUBE[cell.nodes[b]],
            _EXTRUDE_NONEW_CUBE[cell.nodes[c]],
            _EXTRUDE_NONEW_CUBE[cell.nodes[d]])
    end
    total==0 && throw(ErrorException("NoNewVerts: degenerate reference template"))
    total>0 && return cell
    first_swap,second_swap=_EXTRUDE_REVERSE_SWAPS[Int(cell.msh)]
    nodes=ntuple(k->_extrude_nonew_swap_node(cell.nodes,k,first_swap,second_swap),8)
    return _ExtrudeNoNewCell(cell.msh,nodes)
end

@inline function _extrude_nonew_key3(a::UInt8,b::UInt8,c::UInt8)
    a,b=minmax(a,b)
    b,c=minmax(b,c)
    a,b=minmax(a,b)
    return (UInt8(0),a,b,c)
end

@inline function _extrude_nonew_key4(a::UInt8,b::UInt8,c::UInt8,d::UInt8)
    a,b=minmax(a,b)
    c,d=minmax(c,d)
    a,c=minmax(a,c)
    b,d=minmax(b,d)
    b,c=minmax(b,c)
    return (a,b,c,d)
end

function _extrude_nonew_face_counts(cells::Vector{_ExtrudeNoNewCell})
    counts=Dict{NTuple{4,UInt8},UInt8}()
    for cell in cells
        nodes=cell.nodes
        for face in _VOLUME_CELL_FACES[Int(cell.msh)]
            if length(face)==3
                key=_extrude_nonew_key3(nodes[face[1]],nodes[face[2]],nodes[face[3]])
            else
                key=_extrude_nonew_key4(nodes[face[1]],nodes[face[2]],
                                       nodes[face[3]],nodes[face[4]])
            end
            counts[key]=get(counts,key,UInt8(0))+UInt8(1)
        end
    end
    return counts
end

# External quad state: 0 = whole, 1 = AC, 2 = BD, 255 = incompatible complex.
@inline function _extrude_nonew_quad_state(a,b,c,d,counts)
    get(counts,_extrude_nonew_key4(a,b,c,d),UInt8(0))==1 && return UInt8(0)
    if get(counts,_extrude_nonew_key3(a,b,c),UInt8(0))==1 &&
       get(counts,_extrude_nonew_key3(a,c,d),UInt8(0))==1
        return UInt8(1)
    end
    if get(counts,_extrude_nonew_key3(a,b,d),UInt8(0))==1 &&
       get(counts,_extrude_nonew_key3(b,c,d),UInt8(0))==1
        return UInt8(2)
    end
    return typemax(UInt8)
end

# Derive the relation from the emitted typed complex, rather than assuming
# that a construction's requested face choices match its output.
function _extrude_nonew_external_states(cells::Vector{_ExtrudeNoNewCell})
    counts=_extrude_nonew_face_counts(cells)
    all(n->n==1 || n==2,values(counts)) || return nothing
    states=Vector{UInt8}(undef,6)
    expected=0
    for (i,face) in enumerate(_EXTRUDE_NONEW_HEX_FACES)
        a,b,c,d=face
        state=_extrude_nonew_quad_state(UInt8(a),UInt8(b),UInt8(c),UInt8(d),counts)
        state==typemax(UInt8) && return nothing
        states[i]=state
        expected+=state==0 ? 1 : 2
    end
    count(==(UInt8(1)),values(counts))==expected || return nothing
    return (states[1],states[2],states[3],states[4],states[5],states[6])
end

# An apex belongs to three hex faces or three prism faces. Its incident
# quads must split along diagonals through that existing apex. Each remaining
# quad may stay whole or split either way, and each remaining triangle emits
# one tetrahedron. This is also the complete valid full-prism relation.
function _extrude_nonew_corner_cells(faces,apex::Int,choices::Int)
    cells=_ExtrudeNoNewCell[]
    next_choices=choices
    for face in faces
        apex in face && continue
        if length(face)==3
            push!(cells,_extrude_nonew_cell(4,face[1],face[2],face[3],apex))
        else
            state=next_choices%3
            next_choices÷=3
            a,b,c,d=face
            if state==0
                push!(cells,_extrude_nonew_cell(7,a,b,c,d,apex))
            elseif state==1
                push!(cells,_extrude_nonew_cell(4,a,b,c,apex))
                push!(cells,_extrude_nonew_cell(4,a,c,d,apex))
            else
                push!(cells,_extrude_nonew_cell(4,a,b,d,apex))
                push!(cells,_extrude_nonew_cell(4,b,c,d,apex))
            end
        end
    end
    return cells
end

@inline function _extrude_nonew_remap_cell(cell::_ExtrudeNoNewCell,
                                         mapping::NTuple{6,UInt8})
    nodes=ntuple(k->cell.nodes[k]==0 ? UInt8(0) : mapping[Int(cell.nodes[k])],8)
    return _ExtrudeNoNewCell(cell.msh,nodes)
end

@inline _extrude_nonew_shift_mapping(indices::NTuple{6,Int}) =
    ntuple(k->UInt8(indices[k]+1),6)

# Three planes through pairs of opposite face diagonals, with two diagonal
# orientations each. These are prism_v in createFullHexElems, converted to
# one-based local indices. The two wedges share exactly one quadrangle.
function _extrude_nonew_prism_slice(direction::Int,orientation::Int)
    p=direction-1
    if direction<3
        p0=orientation==1 ? p : (p+2)%4
        first=(p0,(p0+1)%4,(p0+1)%4+4,(p0+3)%4,(p0+2)%4,(p0+2)%4+4)
        second=(p0,p0+4,(p0+1)%4+4,(p0+3)%4,(p0+3)%4+4,(p0+2)%4+4)
    else
        p0=orientation-1
        first=(p0,(p0+1)%4,(p0+2)%4,p0+4,(p0+1)%4+4,(p0+2)%4+4)
        second=(p0,(p0+3)%4,(p0+2)%4,p0+4,(p0+3)%4+4,(p0+2)%4+4)
    end
    first_mapping=_extrude_nonew_shift_mapping(first)
    second_mapping=_extrude_nonew_shift_mapping(second)
    return first_mapping,second_mapping
end

function _extrude_nonew_cell_key(cell::_ExtrudeNoNewCell)
    width=msh_num_nodes(Int(cell.msh))
    nodes=sort!(UInt8[cell.nodes[i] for i in 1:width])
    return _ExtrudeNoNewCell(cell.msh,ntuple(k->k<=width ? nodes[k] : UInt8(0),8))
end

function _extrude_nonew_template_key(cells::Vector{_ExtrudeNoNewCell})
    keys=_extrude_nonew_cell_key.(cells)
    sort!(keys;by=cell->(cell.msh,cell.nodes))
    return ntuple(k->k<=length(keys) ? keys[k] : _EXTRUDE_NONEW_EMPTY_CELL,6)
end

function _extrude_nonew_add_template!(templates,seen,cells,kind,choice)
    faces=_extrude_nonew_external_states(cells)
    faces===nothing && return nothing
    length(cells)<=6 || throw(ErrorException("NoNewVerts: reference template exceeds six cells"))
    key=(faces,_extrude_nonew_template_key(cells))
    key in seen && return nothing
    normalized=_extrude_nonew_positive_cell.(cells)
    packed=ntuple(k->k<=length(normalized) ? normalized[k] : _EXTRUDE_NONEW_EMPTY_CELL,6)
    push!(templates,_ExtrudeNoNewTemplate(faces,UInt8(length(cells)),packed,
                                         UInt8(kind),UInt8(choice)))
    push!(seen,key)
    return nothing
end

function _extrude_nonew_build_templates()
    templates=_ExtrudeNoNewTemplate[]
    seen=Set{Tuple{NTuple{6,UInt8},NTuple{6,_ExtrudeNoNewCell}}}()
    _extrude_nonew_add_template!(templates,seen,
        [_extrude_nonew_cell(5,1,2,3,4,5,6,7,8)],0,0)
    for apex in 1:8,choices in 0:26
        cells=_extrude_nonew_corner_cells(_EXTRUDE_NONEW_HEX_FACES,apex,choices)
        _extrude_nonew_add_template!(templates,seen,cells,1,apex)
    end
    prism_cells=Vector{_ExtrudeNoNewCell}[
        [_extrude_nonew_cell(6,1,2,3,4,5,6)]]
    for apex in 1:6,choices in 0:2
        push!(prism_cells,_extrude_nonew_corner_cells(_EXTRUDE_NONEW_PRISM_FACES,
                                                     apex,choices))
    end
    for direction in 1:3,orientation in 1:2
        first_mapping,second_mapping=_extrude_nonew_prism_slice(direction,orientation)
        for first_cells in prism_cells,second_cells in prism_cells
            cells=_ExtrudeNoNewCell[]
            for cell in first_cells
                push!(cells,_extrude_nonew_remap_cell(cell,first_mapping))
            end
            for cell in second_cells
                push!(cells,_extrude_nonew_remap_cell(cell,second_mapping))
            end
            _extrude_nonew_add_template!(templates,seen,cells,2,2*(direction-1)+orientation)
        end
    end
    return templates
end

const _EXTRUDE_NONEW_TEMPLATES = _extrude_nonew_build_templates()

# Private homogeneous table; callers treat it as read-only. Template and cell
# records themselves are immutable and contain no heap-owned connectivity.
@inline _extrude_nonew_templates() = _EXTRUDE_NONEW_TEMPLATES

function _extrude_nonew_emit_template!(sweep::_ExtrudeVolumeSweep,
        v::NTuple{8,NTuple{3,Float64}},template::_ExtrudeNoNewTemplate)
    for i in 1:Int(template.ncells)
        cell=template.cells[i]
        n=cell.nodes
        if cell.msh==4
            push!(sweep.tets,(v[n[1]],v[n[2]],v[n[3]],v[n[4]]))
        elseif cell.msh==5
            push!(sweep.hexes,(v[n[1]],v[n[2]],v[n[3]],v[n[4]],
                                v[n[5]],v[n[6]],v[n[7]],v[n[8]]))
        elseif cell.msh==6
            push!(sweep.prisms,(v[n[1]],v[n[2]],v[n[3]],v[n[4]],v[n[5]],v[n[6]]))
        elseif cell.msh==7
            push!(sweep.pyramids,(v[n[1]],v[n[2]],v[n[3]],v[n[4]],v[n[5]]))
        else
            throw(ErrorException("NoNewVerts: invalid reference-template cell type"))
        end
    end
    return nothing
end
