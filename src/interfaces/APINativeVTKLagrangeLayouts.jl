# Integer reference-node ordering, independently checked against public VTK
# GetParametricCoords lattices. Model coordinates, UIDs and values play no role.
function _native_model_vtk_lagrange_lattice(family,p)
    nodes=NTuple{3,Int}[]
    family in (:tri,:tet,:pyr) && p==0 && return NTuple{3,Int}[(0,0,0)]
    if family in (:tri,:tet)
        vertices=NTuple{3,Int}[(0,0,0),(p,0,0),(0,p,0)]
        family===:tet && push!(vertices,(0,0,p))
        append!(nodes,vertices)
        edges=family===:tri ? ((1,2),(2,3),(3,1)) :
            ((1,2),(2,3),(3,1),(1,4),(2,4),(3,4))
        for (a,b) in edges,i in 1:p-1
            push!(nodes,ntuple(axis->((p-i)*vertices[a][axis]+i*vertices[b][axis])÷p,3))
        end
        if family===:tri && p>=3
            append!(nodes,((x+1,y+1,0) for (x,y,_) in _native_model_vtk_lagrange_lattice(:tri,p-3)))
        elseif family===:tet
            if p>=3
                for (a,b,c) in ((1,2,4),(3,4,2),(1,4,3),(1,3,2)),
                    (x,y,_) in _native_model_vtk_lagrange_lattice(:tri,p-3)
                    u=x+1;v=y+1
                    push!(nodes,ntuple(axis->((p-u-v)*vertices[a][axis]+u*vertices[b][axis]+v*vertices[c][axis])÷p,3))
                end
            end
            p>=4 && append!(nodes,((x+1,y+1,z+1) for (x,y,z) in _native_model_vtk_lagrange_lattice(:tet,p-4)))
        end
        return nodes
    elseif family in (:qua,:hex)
        corners=((0,0),(p,0),(p,p),(0,p))
        levels=family===:qua ? (0,) : (0,p)
        append!(nodes,((x,y,z) for z in levels for (x,y) in corners))
        for z in levels
            append!(nodes,((i,0,z) for i in 1:p-1))
            append!(nodes,((p,i,z) for i in 1:p-1))
            append!(nodes,((i,p,z) for i in 1:p-1))
            append!(nodes,((0,i,z) for i in 1:p-1))
        end
        if family===:hex
            append!(nodes,((x,y,z) for (x,y) in corners for z in 1:p-1))
            for x in (0,p);append!(nodes,((x,y,z) for z in 1:p-1 for y in 1:p-1));end
            for y in (0,p);append!(nodes,((x,y,z) for z in 1:p-1 for x in 1:p-1));end
        end
        for z in levels;append!(nodes,((x,y,z) for y in 1:p-1 for x in 1:p-1));end
        family===:hex && append!(nodes,((x,y,z) for z in 1:p-1 for y in 1:p-1 for x in 1:p-1))
        return nodes
    elseif family===:pyr
        # VTK has no vtkLagrangePyramid class, so type 74 has no authored
        # ordering. The convention here follows the Lagrange template on the
        # pyramid topology of vtkPyramid/vtkTriQuadraticPyramid (public VTK
        # tables): five vertices (counterclockwise base, apex last), eight
        # edges (four base edges counterclockwise, then four lateral edges),
        # four lateral triangular faces in VTK face order, the quadrilateral
        # base face, then the interior. The integer reference lattice is
        # {(i,j,k) : 0<=i,j<=k<=p} with the apex at (0,0,0) and the base at
        # k=p; node (i,j,k) sits at ((2i-k)/p,(2j-k)/p,1-k/p).
        vertices=NTuple{3,Int}[(0,0,p),(p,0,p),(p,p,p),(0,p,p),(0,0,0)]
        append!(nodes,vertices)
        for (a,b) in ((1,2),(2,3),(3,4),(4,1),(1,5),(2,5),(3,5),(4,5)),i in 1:p-1
            push!(nodes,ntuple(axis->((p-i)*vertices[a][axis]+i*vertices[b][axis])÷p,3))
        end
        if p>=3
            for (a,b,c) in ((1,2,5),(2,3,5),(3,4,5),(4,1,5)),
                (x,y,_) in _native_model_vtk_lagrange_lattice(:tri,p-3)
                u=x+1;v=y+1
                push!(nodes,ntuple(axis->((p-u-v)*vertices[a][axis]+u*vertices[b][axis]+
                    v*vertices[c][axis])÷p,3))
            end
        end
        p>=2 && append!(nodes,((x,y,p) for y in 1:p-1 for x in 1:p-1))
        p>=3 && append!(nodes,((x+1,y+1,z+2)
            for (x,y,z) in _native_model_vtk_lagrange_lattice(:pyr,p-3)))
        return nodes
    elseif family===:pri
        vertices=((0,0),(p,0),(0,p));edges=((1,2),(2,3),(3,1))
        append!(nodes,((x,y,z) for z in (0,p) for (x,y) in vertices))
        for z in (0,p),(a,b) in edges,i in 1:p-1
            push!(nodes,(((p-i)*vertices[a][1]+i*vertices[b][1])÷p,
                         ((p-i)*vertices[a][2]+i*vertices[b][2])÷p,z))
        end
        append!(nodes,((x,y,z) for (x,y) in vertices for z in 1:p-1))
        # Wedge interiors use flat triangle rows, unlike recursive Tri cells.
        interior=[(x,y) for y in 1:p-2 for x in 1:p-y-1]
        append!(nodes,((x,y,z) for z in (0,p) for (x,y) in interior))
        for (a,b) in edges,z in 1:p-1,i in 1:p-1
            push!(nodes,(((p-i)*vertices[a][1]+i*vertices[b][1])÷p,
                         ((p-i)*vertices[a][2]+i*vertices[b][2])÷p,z))
        end
        append!(nodes,((x,y,z) for z in 1:p-1 for (x,y) in interior))
        return nodes
    end
    throw(ArgumentError("unsupported VTK Lagrange reference family $family"))
end

function _native_model_vtk_lagrange_layout(msh,spec,caller)
    family=spec.family;p=spec.order
    reference=_native_model_vtk_lagrange_lattice(family,p)
    local_nodes=Elements.lagrange_nodes(msh)
    identities=Dict{NTuple{3,Int},Int}()
    for node in axes(local_nodes,2)
        identity=ntuple(3) do axis
            value=local_nodes[axis,node]
            if family===:pyr
                # Integer (i,j,k) on the {(i,j,k):0<=i,j<=k<=p} lattice:
                # k=p(1-z) and x=(2i-k)/p, so i=p(x+1-z)/2.
                value=axis==3 ? 1-value : (value+1-local_nodes[3,node])/2
            elseif family===:hex || (family===:qua && axis<=2) || (family===:pri && axis==3)
                value=(value+1)/2
            end
            scaled=p*value;tolerance=128 * eps(Float64) * p
            isfinite(scaled) && -tolerance<=scaled<=p+tolerance || throw(ArgumentError(
                "$caller: invalid canonical Lagrange coordinate for Gmsh element type $msh"))
            coordinate=round(Int,scaled)
            abs(scaled-coordinate)<=tolerance || throw(ArgumentError(
                "$caller: non-lattice Lagrange coordinate for Gmsh element type $msh"))
            coordinate
        end
        haskey(identities,identity) && throw(ArgumentError(
            "$caller: duplicate canonical Lagrange node for Gmsh element type $msh"))
        identities[identity]=node
    end
    length(identities)==length(reference)==spec.nnodes && length(Set(reference))==spec.nnodes &&
        all(haskey(identities,node) for node in reference) || throw(ArgumentError(
            "$caller: incomplete canonical Lagrange mapping for Gmsh element type $msh"))
    type=family===:tri ? 69 : family===:qua ? 70 : family===:tet ? 71 :
        family===:hex ? 72 : family===:pri ? 73 : 74
    return (type=type,order=[identities[node] for node in reference])
end
