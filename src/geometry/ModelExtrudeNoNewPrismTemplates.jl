# Existing-corner NoNewVerts relation for a nondegenerate triangular sweep.
# Gmsh 4.15.2 QuadTriExtruded3D.cpp, QuadToTriPriPyrTet (4377): retain the
# prism, or fan its nonincident faces to one of its six existing vertices.
# The reference is the true unit prism; its six IDs are not cube IDs1:6.

struct _ExtrudeNoNewPrismTemplate
    faces::NTuple{3,UInt8}
    ncells::UInt8
    cells::NTuple{3,_ExtrudeNoNewCell}
    kind::UInt8
    choice::UInt8
end

const _EXTRUDE_NONEW_UNIT_PRISM = ((0.0,0.0,0.0),(1.0,0.0,0.0),
    (0.0,1.0,0.0),(0.0,0.0,1.0),(1.0,0.0,1.0),(0.0,1.0,1.0))
const _EXTRUDE_NONEW_PRISM_LATERALS = ((1,2,5,4),(2,3,6,5),(3,1,4,6))
const _EXTRUDE_NONEW_PRISM_CAPS = ((1,3,2),(4,5,6))

function _extrude_nonew_prism_positive_cell(cell::_ExtrudeNoNewCell)
    total=0.0
    for (a,b,c,d) in _EXTRUDE_CELL_TETS[Int(cell.msh)]
        total+=tet_signed_volume(_EXTRUDE_NONEW_UNIT_PRISM[cell.nodes[a]],
            _EXTRUDE_NONEW_UNIT_PRISM[cell.nodes[b]],
            _EXTRUDE_NONEW_UNIT_PRISM[cell.nodes[c]],
            _EXTRUDE_NONEW_UNIT_PRISM[cell.nodes[d]])
    end
    total==0 && throw(ErrorException("NoNewVerts: degenerate unit-prism template"))
    total>0 && return cell
    first_swap,second_swap=_EXTRUDE_REVERSE_SWAPS[Int(cell.msh)]
    nodes=ntuple(k->_extrude_nonew_swap_node(cell.nodes,k,first_swap,second_swap),8)
    return _ExtrudeNoNewCell(cell.msh,nodes)
end

function _extrude_nonew_prism_external_states(cells::Vector{_ExtrudeNoNewCell})
    counts=_extrude_nonew_face_counts(cells)
    all(n->n==1 || n==2,values(counts)) || return nothing
    for cap in _EXTRUDE_NONEW_PRISM_CAPS
        key=_extrude_nonew_key3(UInt8(cap[1]),UInt8(cap[2]),UInt8(cap[3]))
        get(counts,key,UInt8(0))==1 || return nothing
    end
    first,second,third=_EXTRUDE_NONEW_PRISM_LATERALS
    states=ntuple(3) do index
        face=index==1 ? first : index==2 ? second : third
        _extrude_nonew_quad_state(UInt8(face[1]),UInt8(face[2]),
            UInt8(face[3]),UInt8(face[4]),counts)
    end
    any(==(typemax(UInt8)),states) && return nothing
    expected=2+sum(state->state==0 ? 1 : 2,states)
    count(==(UInt8(1)),values(counts))==expected || return nothing
    return states
end

# Exactly one quadrangle and one triangle do not contain a prism apex.
# Its two incident quadrangles receive diagonals through that same apex.
function _extrude_nonew_prism_corner_cells(apex::Int,choice::Int)
    cells=_ExtrudeNoNewCell[]
    for face in (_EXTRUDE_NONEW_PRISM_CAPS...,
                 _EXTRUDE_NONEW_PRISM_LATERALS...)
        apex in face && continue
        if length(face)==3
            push!(cells,_extrude_nonew_cell(4,face[1],face[2],face[3],apex))
        else
            a,b,c,d=face
            if choice==0
                push!(cells,_extrude_nonew_cell(7,a,b,c,d,apex))
            elseif choice==1
                push!(cells,_extrude_nonew_cell(4,a,b,c,apex))
                push!(cells,_extrude_nonew_cell(4,a,c,d,apex))
            elseif choice==2
                push!(cells,_extrude_nonew_cell(4,a,b,d,apex))
                push!(cells,_extrude_nonew_cell(4,b,c,d,apex))
            else
                throw(ErrorException("NoNewVerts: invalid prism corner choice"))
            end
        end
    end
    return cells
end

function _extrude_nonew_prism_template_key(cells::Vector{_ExtrudeNoNewCell})
    keys=_extrude_nonew_cell_key.(cells)
    sort!(keys;by=cell->(cell.msh,cell.nodes))
    return ntuple(k->k<=length(keys) ? keys[k] : _EXTRUDE_NONEW_EMPTY_CELL,3)
end

function _extrude_nonew_prism_add_template!(templates,seen,cells,kind,choice)
    faces=_extrude_nonew_prism_external_states(cells)
    faces===nothing && return nothing
    length(cells)<=3 || throw(ErrorException(
        "NoNewVerts: unit-prism template exceeds three cells"))
    key=(faces,_extrude_nonew_prism_template_key(cells))
    key in seen && return nothing
    normalized=_extrude_nonew_prism_positive_cell.(cells)
    packed=ntuple(k->k<=length(normalized) ? normalized[k] :
        _EXTRUDE_NONEW_EMPTY_CELL,3)
    push!(templates,_ExtrudeNoNewPrismTemplate(faces,UInt8(length(cells)),
        packed,UInt8(kind),UInt8(choice)))
    push!(seen,key)
    return nothing
end

function _extrude_nonew_prism_build_templates()
    templates=_ExtrudeNoNewPrismTemplate[]
    seen=Set{Tuple{NTuple{3,UInt8},NTuple{3,_ExtrudeNoNewCell}}}()
    _extrude_nonew_prism_add_template!(templates,seen,
        [_extrude_nonew_cell(6,1,2,3,4,5,6)],0,0)
    for apex in 1:6,choice in 0:2
        cells=_extrude_nonew_prism_corner_cells(apex,choice)
        _extrude_nonew_prism_add_template!(templates,seen,cells,1,apex)
    end
    length(templates)==13 || throw(ErrorException(
        "NoNewVerts: unit-prism relation must contain thirteen templates"))
    return ntuple(index->templates[index],13)
end

# The complete table is immutable as well as its fixed-size packed records.
const _EXTRUDE_NONEW_PRISM_TEMPLATES = _extrude_nonew_prism_build_templates()
@inline _extrude_nonew_prism_templates() = _EXTRUDE_NONEW_PRISM_TEMPLATES

function _extrude_nonew_emit_template!(sweep::_ExtrudeVolumeSweep,
        v::NTuple{6,NTuple{3,Float64}},template::_ExtrudeNoNewPrismTemplate)
    for index in 1:Int(template.ncells)
        cell=template.cells[index]
        n=cell.nodes
        if cell.msh==4
            push!(sweep.tets,(v[n[1]],v[n[2]],v[n[3]],v[n[4]]))
        elseif cell.msh==6
            push!(sweep.prisms,(v[n[1]],v[n[2]],v[n[3]],v[n[4]],v[n[5]],v[n[6]]))
        elseif cell.msh==7
            push!(sweep.pyramids,(v[n[1]],v[n[2]],v[n[3]],v[n[4]],v[n[5]]))
        else
            throw(ErrorException("NoNewVerts: invalid unit-prism cell type"))
        end
    end
    return nothing
end
