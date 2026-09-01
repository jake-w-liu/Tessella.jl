# Differential oracle for cached global edge and face topology.
# This uses the locally installed Gmsh 4.15.2 Julia API and never starts the GUI.
using Tessella
using SHA
using Tessella.MeshTypes: mesh_crc

function _gmsh_binding()
    configured=get(ENV,"GMSH_JULIA_API","")
    candidates=String[]
    isempty(configured) || push!(candidates,configured)
    executable=Sys.which("gmsh")
    if executable!==nothing
        prefix=dirname(dirname(realpath(executable)))
        push!(candidates,joinpath(prefix,"lib","gmsh.jl"))
    end
    append!(candidates,["/opt/homebrew/lib/gmsh.jl","/usr/local/lib/gmsh.jl"])
    for candidate in unique(candidates)
        isfile(candidate) && return candidate
    end
    error("mesh-entity-topology differential: Gmsh Julia API not found; " *
          "set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION=="4.15.2" || error(
    "mesh-entity-topology differential requires Gmsh API 4.15.2, found " *
    gmsh.GMSH_API_VERSION)

const _TOPOLOGY_COORDINATES=Float64[
    0 1 0 0 1;
    0 0 1 0 1;
    0 0 0 1 1]
const _TOPOLOGY_SEGMENTS=Int32[1 3;2 1]
const _TOPOLOGY_TRIANGLES=Int32[1 3;2 2;3 4]
const _TOPOLOGY_TETRAHEDRA=Int32[1 1;2 3;3 2;4 5]
const _MANUAL_TOPOLOGY_COORDINATES=hcat(
    _TOPOLOGY_COORDINATES,Float64[2,1,0])

function _topology_mesh()
    return Mesh(
        _TOPOLOGY_COORDINATES;
        segs=_TOPOLOGY_SEGMENTS,tris=_TOPOLOGY_TRIANGLES,
        tets=_TOPOLOGY_TETRAHEDRA)
end

function _manual_topology_mesh()
    return Mesh(
        _MANUAL_TOPOLOGY_COORDINATES;
        segs=_TOPOLOGY_SEGMENTS,tris=_TOPOLOGY_TRIANGLES,
        tets=_TOPOLOGY_TETRAHEDRA)
end

function _add_gmsh_topology_model!(name::AbstractString,coordinates)
    gmsh.model.add(name)
    for (dimension,entity) in ((1,11),(2,22),(3,33))
        gmsh.model.addDiscreteEntity(dimension,entity)
    end
    gmsh.model.mesh.addNodes(
        3,33,UInt64.(1:size(coordinates,2)),collect(vec(coordinates)))
    gmsh.model.mesh.addElementsByType(
        11,1,UInt64[101,102],UInt64.(vec(_TOPOLOGY_SEGMENTS)))
    gmsh.model.mesh.addElementsByType(
        22,2,UInt64[201,202],UInt64.(vec(_TOPOLOGY_TRIANGLES)))
    gmsh.model.mesh.addElementsByType(
        33,4,UInt64[301,302],UInt64.(vec(_TOPOLOGY_TETRAHEDRA)))
    return nothing
end

function _install_topology_mesh!(mesh::Mesh)
    lock(Tessella.API.STATE_LOCK) do
        Tessella.API._replace_mesh_cache_locked!(
            Tessella.API._copy_mesh(mesh))
    end
    return mesh_crc(mesh)
end

function _sorted_entities(tags,nodes,width::Int)
    length(nodes)==width*length(tags) || error(
        "topology differential: tag/node lengths are inconsistent")
    order=sortperm(tags)
    sorted_tags=UInt64[tags[index] for index in order]
    sorted_nodes=UInt64[]
    sizehint!(sorted_nodes,length(nodes))
    for index in order
        offset=width*(index-1)
        append!(sorted_nodes,@view nodes[(offset+1):(offset+width)])
    end
    return sorted_tags,sorted_nodes
end

function _same(label,gmsh_values,tessella_values)
    gmsh_values==tessella_values || error(
        "$label differs: Gmsh=$(repr(gmsh_values)), " *
        "Tessella=$(repr(tessella_values))")
    return nothing
end

function _rejects_argument(f::Function)
    try
        f()
        return false
    catch err
        err isa ArgumentError || rethrow()
        return true
    end
end

function _write_topology_values!(stream,label,values)
    write(stream,codeunits(label));write(stream,UInt8(0))
    write(stream,htol(UInt64(length(values))))
    if eltype(values)==Int32
        foreach(value->write(
            stream,htol(reinterpret(UInt32,value))),values)
    else
        foreach(value->write(stream,htol(UInt64(value))),values)
    end
    return nothing
end

function _topology_array_digest(arrays)
    stream=IOBuffer()
    for values in arrays
        write(stream,htol(UInt64(length(values))))
        if eltype(values)==Int32
            foreach(value->write(
                stream,htol(reinterpret(UInt32,value))),values)
        else
            foreach(value->write(stream,htol(UInt64(value))),values)
        end
    end
    return bytes2hex(SHA.sha256(take!(stream)))
end

function _tessella_rejects_without_change(f::Function,baseline)
    _rejects_argument(f) || error(
        "Tessella accepted an intentionally unsupported Gmsh mutation")
    current=(
        Tessella.API.mesh.get_all_edges(),
        Tessella.API.mesh.get_all_faces(3),
        Tessella.API.mesh.get_all_faces(4))
    current==baseline || error(
        "Tessella changed topology after rejecting an insertion batch")
    return nothing
end

try
    gmsh.initialize(String[],false)
    gmsh.option.setNumber("General.Terminal",0)
    _add_gmsh_topology_model!(
        "mesh-entity-topology",_TOPOLOGY_COORDINATES)

    Tessella.API.initialize()
    digests=try
        fixture=_topology_mesh()
        baseline=_install_topology_mesh!(fixture)

        _same("pre-create edges",gmsh.model.mesh.getAllEdges(),
              Tessella.API.mesh.get_all_edges())
        _same("pre-create triangle faces",gmsh.model.mesh.getAllFaces(3),
              Tessella.API.mesh.get_all_faces(3))

        gmsh.model.mesh.createEdges()
        Tessella.API.mesh.create_edges()
        gmsh_edges=_sorted_entities(gmsh.model.mesh.getAllEdges()...,2)
        tessella_edges=Tessella.API.mesh.get_all_edges()
        _same("global edge tag/node map",gmsh_edges,tessella_edges)

        edge_nodes=UInt64[
            1,2,2,1, 1,3,3,1, 2,4,4,2, 1,5,5,1]
        edge_query=gmsh.model.mesh.getEdges(edge_nodes)
        tessella_edge_query=Tessella.API.mesh.get_edges(edge_nodes)
        _same("edge lookup and orientation",edge_query,tessella_edge_query)

        gmsh.model.mesh.createFaces()
        Tessella.API.mesh.create_faces()
        gmsh_faces=_sorted_entities(gmsh.model.mesh.getAllFaces(3)...,3)
        tessella_faces=Tessella.API.mesh.get_all_faces(3)
        _same("global triangular-face tag/node map",gmsh_faces,tessella_faces)
        _same("global quadrangular faces",gmsh.model.mesh.getAllFaces(4),
              Tessella.API.mesh.get_all_faces(4))

        face_nodes=UInt64[
            1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1,
            1,2,4, 4,2,1]
        face_query=gmsh.model.mesh.getFaces(3,face_nodes)
        tessella_face_query=Tessella.API.mesh.get_faces(3,face_nodes)
        _same("triangular-face lookup",face_query,tessella_face_query)
        all(iszero,last(tessella_face_query)) || error(
            "pinned triangular-face orientations are not zero")

        gmsh.model.mesh.getEdges(UInt64[1,2,3])==
            (UInt64[1],Int32[1]) || error(
                "Gmsh malformed-edge probe changed")
        _rejects_argument(()->Tessella.API.mesh.get_edges(
            UInt64[1,2,3])) || error(
                "Tessella silently truncated an incomplete edge pair")
        gmsh.model.mesh.getFaces(3,UInt64[1,2])==
            (UInt64[],Int32[]) || error(
                "Gmsh malformed-face probe changed")
        _rejects_argument(()->Tessella.API.mesh.get_faces(
            3,UInt64[1,2])) || error(
                "Tessella silently truncated an incomplete face")
        _rejects_argument(()->Tessella.API.mesh.create_edges(
            [(3,33)])) || error(
                "Tessella fabricated entity-selective edge metadata")
        _rejects_argument(()->Tessella.API.mesh.create_faces(
            [(3,33)])) || error(
                "Tessella fabricated entity-selective face metadata")
        _rejects_argument(()->Tessella.API.mesh.get_edges(
            UInt64[4,5])) || error(
                "Tessella accepted an unknown edge")
        _rejects_argument(()->Tessella.API.mesh.get_faces(
            3,UInt64[1,4,5])) || error(
                "Tessella accepted an unknown face")
        mesh_crc(Tessella.API.mesh.get())==baseline || error(
            "topology queries mutated the cached mesh")

        stream=IOBuffer()
        for (label,values) in (
            ("edge-tags",tessella_edges[1]),
            ("edge-nodes",tessella_edges[2]),
            ("edge-query-tags",tessella_edge_query[1]),
            ("edge-query-orientations",tessella_edge_query[2]),
            ("face-tags",tessella_faces[1]),
            ("face-nodes",tessella_faces[2]),
            ("face-query-tags",tessella_face_query[1]),
            ("face-query-orientations",tessella_face_query[2]))
            _write_topology_values!(stream,label,values)
        end
        result=bytes2hex(SHA.sha256(take!(stream)))
        result=="7cb14f79831c1fb7f3420e34759c16513f64e4d8a5cad82e89bd4e71507f57d0" || error(
            "mesh entity topology checksum changed to $result")
        automatic_digest=result

        _add_gmsh_topology_model!(
            "manual-mesh-entity-topology",_MANUAL_TOPOLOGY_COORDINATES)
        manual_fixture=_manual_topology_mesh()
        manual_baseline=_install_topology_mesh!(manual_fixture)

        manual_edge_tags=UInt64[100,200]
        manual_edge_nodes=UInt64[4,5,2,5]
        gmsh.model.mesh.addEdges(manual_edge_tags,manual_edge_nodes)
        Tessella.API.mesh.add_edges(manual_edge_tags,manual_edge_nodes)
        manual_triangle_tags=UInt64[100,300]
        manual_triangle_nodes=UInt64[1,4,5,1,2,5]
        gmsh.model.mesh.addFaces(
            3,manual_triangle_tags,manual_triangle_nodes)
        Tessella.API.mesh.add_faces(
            3,manual_triangle_tags,manual_triangle_nodes)
        gmsh.model.mesh.addFaces(4,UInt64[200],UInt64[1,2,4,5])
        Tessella.API.mesh.add_faces(4,UInt64[200],UInt64[1,2,4,5])

        _same(
            "manual pre-create edges",
            _sorted_entities(gmsh.model.mesh.getAllEdges()...,2),
            Tessella.API.mesh.get_all_edges())
        _same(
            "manual pre-create triangular faces",
            _sorted_entities(gmsh.model.mesh.getAllFaces(3)...,3),
            Tessella.API.mesh.get_all_faces(3))
        _same(
            "manual pre-create quadrangular faces",
            _sorted_entities(gmsh.model.mesh.getAllFaces(4)...,4),
            Tessella.API.mesh.get_all_faces(4))

        gmsh.model.mesh.addEdges(UInt64[100],UInt64[5,4])
        Tessella.API.mesh.add_edges(UInt64[100],UInt64[5,4])
        gmsh.model.mesh.addFaces(3,UInt64[100],UInt64[5,4,1])
        Tessella.API.mesh.add_faces(3,UInt64[100],UInt64[5,4,1])
        gmsh.model.mesh.addFaces(4,UInt64[200],UInt64[5,4,2,1])
        Tessella.API.mesh.add_faces(4,UInt64[200],UInt64[5,4,2,1])

        gmsh.model.mesh.createEdges()
        Tessella.API.mesh.create_edges()
        gmsh.model.mesh.createFaces()
        Tessella.API.mesh.create_faces()
        completed_edges=_sorted_entities(
            gmsh.model.mesh.getAllEdges()...,2)
        completed_triangles=_sorted_entities(
            gmsh.model.mesh.getAllFaces(3)...,3)
        completed_quadrangles=_sorted_entities(
            gmsh.model.mesh.getAllFaces(4)...,4)
        _same(
            "manual create-after-add edge map",completed_edges,
            Tessella.API.mesh.get_all_edges())
        _same(
            "manual create-after-add triangular-face map",
            completed_triangles,Tessella.API.mesh.get_all_faces(3))
        _same(
            "manual create-after-add quadrangular-face map",
            completed_quadrangles,Tessella.API.mesh.get_all_faces(4))
        _same(
            "manual create-after-add edge tags",
            UInt64[3,4,5,6,7,8,9,10,100,200],completed_edges[1])
        _same(
            "manual create-after-add triangle tags",
            UInt64[4,5,6,7,8,9,100,300],completed_triangles[1])

        gmsh.model.mesh.addEdges(UInt64[300],UInt64[5,6])
        Tessella.API.mesh.add_edges(UInt64[300],UInt64[5,6])
        gmsh.model.mesh.addFaces(3,UInt64[400],UInt64[2,5,6])
        Tessella.API.mesh.add_faces(3,UInt64[400],UInt64[2,5,6])

        final_edges=_sorted_entities(gmsh.model.mesh.getAllEdges()...,2)
        final_triangles=_sorted_entities(
            gmsh.model.mesh.getAllFaces(3)...,3)
        final_quadrangles=_sorted_entities(
            gmsh.model.mesh.getAllFaces(4)...,4)
        tessella_final_edges=Tessella.API.mesh.get_all_edges()
        tessella_final_triangles=Tessella.API.mesh.get_all_faces(3)
        tessella_final_quadrangles=Tessella.API.mesh.get_all_faces(4)
        _same("manual final edge map",final_edges,tessella_final_edges)
        _same(
            "manual final triangular-face map",final_triangles,
            tessella_final_triangles)
        _same(
            "manual final quadrangular-face map",final_quadrangles,
            tessella_final_quadrangles)

        manual_edge_query_nodes=UInt64[4,5,5,4,5,6,6,5]
        manual_triangle_query_nodes=UInt64[1,4,5,5,4,1,2,5,6]
        manual_quadrangle_query_nodes=UInt64[1,2,4,5,5,4,2,1]
        manual_edge_query=gmsh.model.mesh.getEdges(
            manual_edge_query_nodes)
        tessella_manual_edge_query=Tessella.API.mesh.get_edges(
            manual_edge_query_nodes)
        _same(
            "manual edge lookup and orientation",manual_edge_query,
            tessella_manual_edge_query)
        manual_triangle_query=gmsh.model.mesh.getFaces(
            3,manual_triangle_query_nodes)
        tessella_manual_triangle_query=Tessella.API.mesh.get_faces(
            3,manual_triangle_query_nodes)
        _same(
            "manual triangular-face lookup",manual_triangle_query,
            tessella_manual_triangle_query)
        manual_quadrangle_query=gmsh.model.mesh.getFaces(
            4,manual_quadrangle_query_nodes)
        tessella_manual_quadrangle_query=Tessella.API.mesh.get_faces(
            4,manual_quadrangle_query_nodes)
        _same(
            "manual quadrangular-face lookup",manual_quadrangle_query,
            tessella_manual_quadrangle_query)
        all(iszero,last(tessella_manual_triangle_query)) || error(
            "manual triangular-face orientations are not zero")
        all(iszero,last(tessella_manual_quadrangle_query)) || error(
            "manual quadrangular-face orientations are not zero")
        mesh_crc(Tessella.API.mesh.get())==manual_baseline || error(
            "manual topology operations mutated the cached mesh")

        manual_digest=_topology_array_digest((
            tessella_final_edges...,
            tessella_manual_edge_query...,
            tessella_final_triangles...,
            tessella_final_quadrangles...,
            tessella_manual_triangle_query...,
            tessella_manual_quadrangle_query...))
        manual_digest==
            "00110d2815d99ced60d72a8344958d0f0797fc5b85c938142dc4061f0abc8b06" ||
            error("manual mesh entity topology checksum changed to " *
                  manual_digest)

        stable_topology=(
            tessella_final_edges,tessella_final_triangles,
            tessella_final_quadrangles)

        gmsh.model.mesh.addEdges(UInt64[999],UInt64[6,5])
        gmsh.model.mesh.getEdges(UInt64[5,6])[1]==UInt64[300] || error(
            "Gmsh duplicate-edge geometry behavior changed")
        _tessella_rejects_without_change(
            ()->Tessella.API.mesh.add_edges(UInt64[999],UInt64[6,5]),
            stable_topology)

        gmsh.model.mesh.addEdges(UInt64[100],UInt64[1,6])
        gmsh.model.mesh.getEdges(UInt64[1,6])[1]==UInt64[100] || error(
            "Gmsh duplicate edge-tag behavior changed")
        _tessella_rejects_without_change(
            ()->Tessella.API.mesh.add_edges(UInt64[100],UInt64[1,6]),
            stable_topology)

        gmsh.model.mesh.addEdges(UInt64[0],UInt64[4,6])
        only(gmsh.model.mesh.getEdges(UInt64[4,6])[1])>0 || error(
            "Gmsh zero edge-tag assignment behavior changed")
        _tessella_rejects_without_change(
            ()->Tessella.API.mesh.add_edges(UInt64[0],UInt64[4,6]),
            stable_topology)

        gmsh.model.mesh.addEdges(UInt64[500],UInt64[6,6])
        gmsh.model.mesh.getEdges(UInt64[6,6])[1]==UInt64[500] || error(
            "Gmsh repeated-node edge behavior changed")
        _tessella_rejects_without_change(
            ()->Tessella.API.mesh.add_edges(UInt64[500],UInt64[6,6]),
            stable_topology)

        partial_failed=false
        try
            gmsh.model.mesh.addEdges(
                UInt64[600,601],UInt64[3,6,3,99])
        catch err
            err isa InterruptException && rethrow()
            (err isa ErrorException && err.msg=="Unknown node 99") ||
                rethrow()
            partial_failed=true
        end
        partial_failed || error(
            "Gmsh invalid-node edge batch unexpectedly succeeded")
        gmsh.model.mesh.getEdges(UInt64[3,6])[1]==UInt64[600] || error(
            "Gmsh partial edge insertion behavior changed")
        _tessella_rejects_without_change(
            ()->Tessella.API.mesh.add_edges(
                UInt64[600,601],UInt64[3,6,3,99]),stable_topology)

        gmsh.model.mesh.addFaces(
            4,UInt64[100],UInt64[1,2,3,6])
        gmsh.model.mesh.getFaces(
            4,UInt64[1,2,3,6])[1]==UInt64[100] || error(
            "Gmsh shared face-tag behavior changed")
        _tessella_rejects_without_change(
            ()->Tessella.API.mesh.add_faces(
                4,UInt64[100],UInt64[1,2,3,6]),stable_topology)

        automatic_digest,manual_digest
    finally
        Tessella.API.finalize()
    end

    println("mesh-entity-topology differential: Gmsh ",
            gmsh.GMSH_API_VERSION,
            " edges=9 triangle_faces=7 edge_queries=8 face_queries=8 " *
            "manual_edges=11 manual_triangle_faces=9 " *
            "manual_quadrangle_faces=1 divergences=6 " *
            "face_orientations=zero automatic_sha=",digests[1],
            " manual_sha=",digests[2])
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
