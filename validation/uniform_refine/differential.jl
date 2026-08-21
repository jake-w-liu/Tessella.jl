# Differential oracle for one-pass linear simplex refinement. This intentionally
# uses the locally installed Gmsh 4.15.2 Julia API and never starts the GUI.
using Tessella

if !isdefined(Tessella, :Refine)
    Base.include(Tessella, joinpath(@__DIR__, "..", "..", "src", "Refine.jl"))
end

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
    error("uniform-refine differential: Gmsh Julia API not found; set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION == "4.15.2" || error(
    "uniform-refine differential requires Gmsh API 4.15.2, found $(gmsh.GMSH_API_VERSION)")

function _tessella_fixture()
    coords = Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]
    return Tessella.MeshTypes.Mesh(
        coords;
        segs=reshape(Int32[1, 2], 2, 1),
        tris=reshape(Int32[1, 2, 3], 3, 1),
        tets=reshape(Int32[1, 2, 3, 4], 4, 1),
        seg_tag=Int32[11], tri_tag=Int32[22], tet_tag=Int32[33])
end

function _gmsh_cells(dim::Int, entity::Int, element_type::Int, arity::Int,
                     points::Dict{UInt64,NTuple{3,Float64}})
    types, _, node_blocks = gmsh.model.mesh.getElements(dim, entity)
    index = findfirst(==(Int32(element_type)), types)
    index === nothing && error("Gmsh returned no type-$element_type elements on ($dim,$entity)")
    nodes = node_blocks[index]
    length(nodes) % arity == 0 || error("Gmsh returned malformed type-$element_type connectivity")
    return [ntuple(local_node -> points[nodes[arity * (cell - 1) + local_node]], arity)
            for cell in 1:length(nodes) ÷ arity]
end

function _tessella_cells(mesh, cells::Matrix{Int32})
    return [ntuple(local_node -> Tessella.MeshTypes.node(mesh, cells[local_node, cell]),
                   size(cells, 1))
            for cell in axes(cells, 2)]
end

try
    gmsh.initialize(String[], false)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("uniform-refine-differential")
    for (dimension, entity) in ((1, 101), (2, 201), (3, 301))
        gmsh.model.addDiscreteEntity(dimension, entity)
    end
    gmsh.model.mesh.addNodes(
        3, 301, UInt64[1, 2, 3, 4],
        Float64[0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1])
    gmsh.model.mesh.addElementsByType(101, 1, UInt64[1], UInt64[1, 2])
    gmsh.model.mesh.addElementsByType(201, 2, UInt64[2], UInt64[1, 2, 3])
    gmsh.model.mesh.addElementsByType(301, 4, UInt64[3], UInt64[1, 2, 3, 4])
    for (dimension, entity, physical) in ((1, 101, 11), (2, 201, 22), (3, 301, 33))
        gmsh.model.addPhysicalGroup(dimension, [entity], physical)
    end

    gmsh.model.mesh.refine()
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    length(coordinates) == 3length(node_tags) || error("Gmsh returned malformed coordinates")
    points = Dict{UInt64,NTuple{3,Float64}}()
    for index in eachindex(node_tags)
        points[node_tags[index]] = (coordinates[3index - 2], coordinates[3index - 1],
                                    coordinates[3index])
    end

    tessella = Tessella.Refine.refine_uniform(_tessella_fixture())
    length(points) == Tessella.MeshTypes.nnodes(tessella) == 10 ||
        error("Gmsh/Tessella refined node counts differ")
    gmsh_segments = _gmsh_cells(1, 101, 1, 2, points)
    gmsh_triangles = _gmsh_cells(2, 201, 2, 3, points)
    gmsh_tetrahedra = _gmsh_cells(3, 301, 4, 4, points)
    gmsh_segments == _tessella_cells(tessella, tessella.segs) ||
        error("segment child template differs from Gmsh 4.15.2")
    gmsh_triangles == _tessella_cells(tessella, tessella.tris) ||
        error("triangle child template differs from Gmsh 4.15.2")
    gmsh_tetrahedra == _tessella_cells(tessella, tessella.tets) ||
        error("tetrahedron child template differs from Gmsh 4.15.2")

    for (dimension, entity, physical, tags) in
        ((1, 101, 11, tessella.seg_tag),
         (2, 201, 22, tessella.tri_tag),
         (3, 301, 33, tessella.tet_tag))
        gmsh.model.getPhysicalGroupsForEntity(dimension, entity) == Int32[physical] ||
            error("Gmsh physical tag was not preserved on entity ($dimension,$entity)")
        all(==(Int32(physical)), tags) || error("Tessella physical tag was not preserved")
    end

    crc = Tessella.MeshTypes.mesh_crc(tessella)
    crc.sha == "db9a1713d1174be1035ef3e9d6380a01ed419797a91ded9a2b8508d0b038f031" ||
        error("unexpected uniform-refine checksum $(crc.sha)")
    println("uniform-refine differential: Gmsh $(gmsh.GMSH_API_VERSION), ",
            "nodes=$(crc.n_nodes) segs=$(crc.n_segs) tris=$(crc.n_tris) ",
            "tets=$(crc.n_tets) sha=$(crc.sha)")
finally
    gmsh.isInitialized() != 0 && gmsh.finalize()
end
