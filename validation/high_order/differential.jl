# Differential oracle for a straight 10-node tetrahedron. Gmsh elevates the
# same linear type-4 element; Tessella's in-memory connectivity and its written
# type-11 record must agree with Gmsh in every local node slot.
using Tessella
using Tessella.MeshTypes: nnodes, ntets
import SHA

function _gmsh_binding()
    configured = get(ENV, "GMSH_JULIA_API", "")
    candidates = String[]
    isempty(configured) || push!(candidates, configured)
    executable = Sys.which("gmsh")
    if executable !== nothing
        prefix = dirname(dirname(realpath(executable)))
        push!(candidates, joinpath(prefix, "lib", "gmsh.jl"))
    end
    append!(candidates, ["/opt/homebrew/lib/gmsh.jl", "/usr/local/lib/gmsh.jl"])
    for candidate in unique(candidates)
        isfile(candidate) && return candidate
    end
    error("high-order differential: Gmsh Julia API not found; set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION == "4.15.2" || error(
    "high-order differential requires Gmsh API 4.15.2, found $(gmsh.GMSH_API_VERSION)")

function _type11_coordinates()
    types, element_tags, node_blocks = gmsh.model.mesh.getElements(3)
    index = findfirst(==(Int32(11)), types)
    index === nothing && error("Gmsh returned no type-11 tetrahedron")
    length(element_tags[index]) == 1 || error(
        "Gmsh returned $(length(element_tags[index])) type-11 elements instead of one")
    connectivity = node_blocks[index]
    length(connectivity) == 10 || error(
        "Gmsh returned malformed type-11 connectivity")
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    length(coordinates) == 3length(node_tags) || error(
        "Gmsh returned malformed node coordinates")
    points = Dict{UInt64,NTuple{3,Float64}}()
    for index in eachindex(node_tags)
        points[node_tags[index]] = (coordinates[3index - 2],
                                    coordinates[3index - 1],
                                    coordinates[3index])
    end
    return Tuple(points[tag] for tag in connectivity)
end

function _tessella_fixture()
    return Mesh(Float64[0 1 0 0; 0 0 1 0; 0 0 0 1];
                tets=reshape(Int32[1,2,3,4], 4, 1), tet_tag=Int32[23])
end

try
    gmsh.initialize(String[], false)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("high-order-differential")
    gmsh.model.addDiscreteEntity(3, 301)
    gmsh.model.mesh.addNodes(
        3, 301, UInt64[1, 2, 3, 4],
        Float64[0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1])
    gmsh.model.mesh.addElementsByType(301, 4, UInt64[1], UInt64[1, 2, 3, 4])
    gmsh.model.addPhysicalGroup(3, [301], 23)
    gmsh.model.mesh.setOrder(2)
    gmsh_coordinates = _type11_coordinates()

    source = _tessella_fixture()
    quadratic = p2_tetmesh(source)
    tessella_coordinates = Tuple(
        (quadratic.coords[1, quadratic.tet10[slot, 1]],
         quadratic.coords[2, quadratic.tet10[slot, 1]],
         quadratic.coords[3, quadratic.tet10[slot, 1]]) for slot in 1:10)
    gmsh_coordinates == tessella_coordinates || error(
        "Gmsh and Tessella disagree on type-11 local node coordinates: " *
        "Gmsh=$gmsh_coordinates Tessella=$tessella_coordinates")
    quadratic.tet_tag == Int32[23] || error(
        "Tessella did not preserve the input tetrahedron tag")
    validate(quadratic).ok || error("Tessella quadratic fixture is invalid")

    mktempdir() do directory
        path = joinpath(directory, "tessella-p2.msh")
        write_msh_p2(path, quadratic)
        digest = bytes2hex(SHA.sha256(read(path)))
        digest == "5a83ebe0386bda71c6761148ed3fe2f964f16c2da2f0b66b6951ef558f4927ab" ||
            error("unexpected Tessella P2 file checksum $digest")
        gmsh.clear()
        gmsh.open(path)
        _type11_coordinates() == tessella_coordinates || error(
            "Gmsh changed Tessella's written type-11 local node coordinates")
        (Int32(3), Int32(23)) in gmsh.model.getPhysicalGroups() || error(
            "Gmsh did not recover physical volume 23 from Tessella's P2 file")
    end

    println("HIGH_ORDER_DIFFERENTIAL_OK gmsh=$(gmsh.GMSH_API_VERSION) ",
            "nodes=$(nnodes(quadratic)) tets=$(ntets(quadratic)) ",
            "sha=5a83ebe0386bda71c6761148ed3fe2f964f16c2da2f0b66b6951ef558f4927ab")
finally
    gmsh.isInitialized() != 0 && gmsh.finalize()
end
