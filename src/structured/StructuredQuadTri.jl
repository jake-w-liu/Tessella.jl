"""Boundary-diagonal transitions for five- and six-face transfinite volumes.

The logical prism/hexahedron is retained unless a boundary face is divided.
A divided cell receives one centroid and cones its faces into tetrahedra and
pyramids, following Gmsh 4.15.2's transfinite QuadTri construction. Boundary
diagonals come from the actual surface meshes, including their arrangements.
"""
module StructuredQuadTri

using ..Elements: ElementBlock, MixedMesh, _VOLUME_CELL_FACES
using ..MeshTypes: nnodes, validate, _throw_simplex_validation, tet_signed_volume
using ..Predicates: orient3
using ..StructuredRecombine: _canon3, _canon4, _node

# Outward cyclic windings for positively oriented MSH cells.
const _HEX_FACES = ((1,4,3,2),(5,6,7,8),(1,2,6,5),
                    (2,3,7,6),(3,4,8,7),(4,1,5,8))
const _PRISM_FACES = ((1,3,2),(4,5,6),(1,2,5,4),
                      (2,3,6,5),(3,1,4,6))
@inline _edge(a::Int32,b::Int32) = minmax(a,b)
@inline _facekey(a::Int32,b::Int32,c::Int32) = (_canon3(a,b,c)...,Int32(0))
@inline _facekey(a::Int32,b::Int32,c::Int32,d::Int32) = _canon4(a,b,c,d)

@inline function _diagonal(v,face,diags)
    length(face)==3 && return UInt8(0)
    a,b,c,d = v[face[1]],v[face[2]],v[face[3]],v[face[4]]
    first = _edge(a,c) in diags
    second = _edge(b,d) in diags
    first && second && throw(ArgumentError(
        "TransfQuadTri: boundary quadrangle has conflicting diagonals"))
    return first ? UInt8(1) : second ? UInt8(2) : UInt8(0)
end

function _surface_data(base,parts,caller)
    lookup=Dict{NTuple{3,Float64},Int32}()
    sizehint!(lookup,nnodes(base))
    for i in 1:nnodes(base)
        get!(lookup,_node(base.coords,Int32(i)),Int32(i))
    end
    diags=Set{NTuple{2,Int32}}()
    tris=NTuple{3,Int32}[]; tri_tags=Int32[]
    quads=NTuple{4,Int32}[]; quad_tags=Int32[]
    for (tag,mesh) in parts
        remap=Vector{Int32}(undef,nnodes(mesh))
        for i in 1:nnodes(mesh)
            id=get(lookup,_node(mesh.coords,Int32(i)),Int32(0))
            id!=0 || throw(ArgumentError(
                "$caller: surface $tag has a vertex outside the transfinite grid"))
            remap[i]=id
        end
        # The triangular wedge band on a recombined triangular surface does
        # not prescribe diagonals. Upstream skips a face with any quads.
        has_quads=mesh isa MixedMesh && any(
            b isa ElementBlock && b.msh==3 && size(b.nodes,2)>0 for b in mesh.blocks)
        if mesh isa MixedMesh
            for block in mesh.blocks
                block isa ElementBlock || continue
                if block.msh==2
                    for i in axes(block.nodes,2)
                        v=ntuple(k->remap[block.nodes[k,i]],3)
                        push!(tris,v);push!(tri_tags,Int32(tag))
                        if !has_quads
                            push!(diags,_edge(v[1],v[2]))
                            push!(diags,_edge(v[2],v[3]))
                            push!(diags,_edge(v[3],v[1]))
                        end
                    end
                elseif block.msh==3
                    for i in axes(block.nodes,2)
                        push!(quads,ntuple(k->remap[block.nodes[k,i]],4))
                        push!(quad_tags,Int32(tag))
                    end
                end
            end
        else
            for i in axes(mesh.tris,2)
                v=ntuple(k->remap[mesh.tris[k,i]],3)
                push!(tris,v);push!(tri_tags,Int32(tag))
                push!(diags,_edge(v[1],v[2]))
                push!(diags,_edge(v[2],v[3]))
                push!(diags,_edge(v[3],v[1]))
            end
        end
    end
    return diags,tris,tri_tags,quads,quad_tags
end

function _centroid(coords,v)
    x=0.0;y=0.0;z=0.0
    for id in v
        x+=coords[1,id];y+=coords[2,id];z+=coords[3,id]
    end
    p=(x/length(v),y/length(v),z/length(v))
    # Preserve the oracle's ordinary summation order; scale only when a
    # finite centroid overflowed during summation at extreme coordinates.
    if !all(isfinite,p)
        x=0.0;y=0.0;z=0.0
        for id in v
            x+=coords[1,id]/length(v)
            y+=coords[2,id]/length(v)
            z+=coords[3,id]/length(v)
        end
        p=(x,y,z)
    end
    return p
end

function _write_tet!(nodes,pos,coords,a,b,c,d,caller)
    # Faces use outward winding: their apex must lie strictly inward. This
    # certifies the centroid fan instead of repairing a folded logical cell.
    orient3(_node(coords,a),_node(coords,b),_node(coords,c),_node(coords,d))>0 ||
        throw(ArgumentError("$caller: centroid fan has a folded or degenerate face"))
    volume=tet_signed_volume(_node(coords,a),_node(coords,c),
                             _node(coords,b),_node(coords,d))
    isfinite(volume) && volume>0 || throw(ArgumentError(
        "$caller: centroid tetrahedron volume is not finite and positive"))
    nodes[1,pos]=a;nodes[2,pos]=c;nodes[3,pos]=b;nodes[4,pos]=d
    return pos+1
end

function _write_pyramid!(nodes,pos,coords,a,b,c,d,apex,caller)
    # Certify both possible base diagonals: neighboring pyramids retain one
    # unsplit quadrangle even if their tetrahedral shadows choose differently.
    for face in ((a,b,c),(a,c,d),(a,b,d),(b,c,d))
        orient3(_node(coords,face[1]),_node(coords,face[2]),
                _node(coords,face[3]),_node(coords,apex))>0 ||
            throw(ArgumentError("$caller: centroid fan has a folded pyramid base"))
        volume=tet_signed_volume(_node(coords,face[1]),_node(coords,face[3]),
                                 _node(coords,face[2]),_node(coords,apex))
        isfinite(volume) && volume>0 || throw(ArgumentError(
            "$caller: centroid pyramid volume is not finite and positive"))
    end
    nodes[1,pos]=c;nodes[2,pos]=b;nodes[3,pos]=a
    nodes[4,pos]=d;nodes[5,pos]=apex
    return pos+1
end

function _certify_boundary(blocks,tris,quads,caller)
    faces=Dict{NTuple{4,Int32},Int}()
    for block in blocks
        local_faces=get(_VOLUME_CELL_FACES,block.msh,nothing)
        local_faces===nothing && continue
        for cell in axes(block.nodes,2),face in local_faces
            key=length(face)==3 ? _facekey(
                block.nodes[face[1],cell],block.nodes[face[2],cell],
                block.nodes[face[3],cell]) : _facekey(
                block.nodes[face[1],cell],block.nodes[face[2],cell],
                block.nodes[face[3],cell],block.nodes[face[4],cell])
            incidence=get(faces,key,0)+1
            incidence<=2 || throw(ArgumentError(
                "$caller: transition volume has a nonmanifold face"))
            faces[key]=incidence
        end
    end
    for cells in (tris,quads),v in cells
        key=_facekey(v...)
        get(faces,key,0)==1 || throw(ArgumentError(
            "$caller: surface cells do not match the transition volume boundary"))
        delete!(faces,key)
    end
    all(==(2),values(faces)) || throw(ArgumentError(
        "$caller: transition volume leaves an uncovered boundary face"))
    return nothing
end

function _matrix(cells::Vector{NTuple{N,Int32}}) where N
    nodes=Matrix{Int32}(undef,N,length(cells))
    for i in eachindex(cells),k in 1:N
        nodes[k,i]=cells[i][k]
    end
    return nodes
end

function _transfinite_quadtri(base::MixedMesh,parts;
        caller="TransfQuadTri",max_nodes=10_000_000,max_tets=60_000_000)
    diags,tris,tri_tags,quads,quad_tags=_surface_data(base,parts,caller)
    ncenter=0;ntet=0;npyr=0;nhex=0;nprism=0
    # Count before allocating output; recompute small face masks on emission
    # instead of retaining a per-cell scratch record for the entire volume.
    for block in base.blocks
        block isa ElementBlock && block.msh in (5,6) || continue
        local_faces=block.msh==5 ? _HEX_FACES : _PRISM_FACES
        for cell in axes(block.nodes,2)
            v=block.msh==5 ? ntuple(k->block.nodes[k,cell],8) :
                             ntuple(k->block.nodes[k,cell],6)
            flags=map(face->_diagonal(v,face,diags),local_faces)
            if all(iszero,flags)
                block.msh==5 ? (nhex+=1) : (nprism+=1)
            else
                ncenter+=1
                for (face,flag) in zip(local_faces,flags)
                    length(face)==3 ? (ntet+=1) :
                        flag==0 ? (npyr+=1) : (ntet+=2)
                end
            end
        end
    end
    total_nodes=BigInt(nnodes(base))+ncenter
    total_nodes<=min(max_nodes,typemax(Int32)) || throw(ArgumentError(
        "$caller: transition node count $total_nodes exceeds max_nodes"))
    shadow_count=BigInt(ntet)+2BigInt(npyr)+6BigInt(nhex)+3BigInt(nprism)
    shadow_count<=max_tets || throw(ArgumentError(
        "$caller: transition tetrahedron count $shadow_count exceeds max_tets"))
    coords=Matrix{Float64}(undef,3,Int(total_nodes))
    copyto!(coords,1,base.coords,1,length(base.coords))
    tetnodes=Matrix{Int32}(undef,4,ntet);tettags=Vector{Int32}(undef,ntet)
    pyrnodes=Matrix{Int32}(undef,5,npyr);pyrtags=Vector{Int32}(undef,npyr)
    hexnodes=Matrix{Int32}(undef,8,nhex);hextags=Vector{Int32}(undef,nhex)
    prinodes=Matrix{Int32}(undef,6,nprism);pritags=Vector{Int32}(undef,nprism)
    ni=nnodes(base);ti=1;pi=1;hi=1;ri=1
    for block in base.blocks
        block isa ElementBlock && block.msh in (5,6) || continue
        local_faces=block.msh==5 ? _HEX_FACES : _PRISM_FACES
        for cell in axes(block.nodes,2)
            v=block.msh==5 ? ntuple(k->block.nodes[k,cell],8) :
                             ntuple(k->block.nodes[k,cell],6)
            flags=map(face->_diagonal(v,face,diags),local_faces)
            tag=block.tags[cell]
            if all(iszero,flags)
                if block.msh==5
                    for k in 1:8;hexnodes[k,hi]=v[k];end
                    hextags[hi]=tag;hi+=1
                else
                    for k in 1:6;prinodes[k,ri]=v[k];end
                    pritags[ri]=tag;ri+=1
                end
                continue
            end
            ni+=1;apex=Int32(ni)
            p=_centroid(coords,v)
            coords[1,ni]=p[1];coords[2,ni]=p[2];coords[3,ni]=p[3]
            for (face,flag) in zip(local_faces,flags)
                a,b,c=v[face[1]],v[face[2]],v[face[3]]
                if length(face)==3
                    tettags[ti]=tag
                    ti=_write_tet!(tetnodes,ti,coords,a,b,c,apex,caller)
                elseif flag==0
                    pyrtags[pi]=tag
                    pi=_write_pyramid!(pyrnodes,pi,coords,a,b,c,v[face[4]],apex,caller)
                else
                    d=v[face[4]]
                    tettags[ti]=tag;tettags[ti+1]=tag
                    if flag==1
                        ti=_write_tet!(tetnodes,ti,coords,a,b,c,apex,caller)
                        ti=_write_tet!(tetnodes,ti,coords,a,c,d,apex,caller)
                    else
                        ti=_write_tet!(tetnodes,ti,coords,a,b,d,apex,caller)
                        ti=_write_tet!(tetnodes,ti,coords,b,c,d,apex,caller)
                    end
                end
            end
        end
    end
    blocks=ElementBlock[]
    !isempty(tris) && push!(blocks,ElementBlock(2,_matrix(tris),tri_tags))
    !isempty(quads) && push!(blocks,ElementBlock(3,_matrix(quads),quad_tags))
    ntet>0 && push!(blocks,ElementBlock(4,tetnodes,tettags))
    nhex>0 && push!(blocks,ElementBlock(5,hexnodes,hextags))
    nprism>0 && push!(blocks,ElementBlock(6,prinodes,pritags))
    npyr>0 && push!(blocks,ElementBlock(7,pyrnodes,pyrtags))
    _certify_boundary(blocks,tris,quads,caller)
    mesh=MixedMesh(coords,blocks)
    diagnostic=validate(mesh)
    diagnostic.ok || _throw_simplex_validation(caller,diagnostic.messages)
    return mesh
end

end
