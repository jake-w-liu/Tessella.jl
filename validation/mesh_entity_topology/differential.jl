# Differential oracle for cached global simplex edge and face topology.
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

function _topology_mesh()
    return Mesh(
        _TOPOLOGY_COORDINATES;
        segs=_TOPOLOGY_SEGMENTS,tris=_TOPOLOGY_TRIANGLES,
        tets=_TOPOLOGY_TETRAHEDRA)
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

try
    gmsh.initialize(String[],false)
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.model.add("mesh-entity-topology")
    for (dimension,entity) in ((1,11),(2,22),(3,33))
        gmsh.model.addDiscreteEntity(dimension,entity)
    end
    gmsh.model.mesh.addNodes(
        3,33,UInt64.(1:size(_TOPOLOGY_COORDINATES,2)),
        collect(vec(_TOPOLOGY_COORDINATES)))
    gmsh.model.mesh.addElementsByType(
        11,1,UInt64[101,102],UInt64.(vec(_TOPOLOGY_SEGMENTS)))
    gmsh.model.mesh.addElementsByType(
        22,2,UInt64[201,202],UInt64.(vec(_TOPOLOGY_TRIANGLES)))
    gmsh.model.mesh.addElementsByType(
        33,4,UInt64[301,302],UInt64.(vec(_TOPOLOGY_TETRAHEDRA)))

    Tessella.API.initialize()
    digest=try
        fixture=_topology_mesh()
        baseline=mesh_crc(fixture)
        lock(Tessella.API.STATE_LOCK) do
            Tessella.API._replace_mesh_cache_locked!(
                Tessella.API._copy_mesh(fixture))
        end

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
        result
    finally
        Tessella.API.finalize()
    end

    println("mesh-entity-topology differential: Gmsh ",
            gmsh.GMSH_API_VERSION,
            " edges=9 triangle_faces=7 edge_queries=8 face_queries=8 " *
            "face_orientations=zero sha=",digest)
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
