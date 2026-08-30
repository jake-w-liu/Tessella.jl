# Differential oracle for read-only bulk node and element access. This uses the
# locally installed Gmsh 4.15.2 Julia API and never starts the GUI.
using Tessella

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
    error("mesh-data-query differential: Gmsh Julia API not found; " *
          "set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION=="4.15.2" || error(
    "mesh-data-query differential requires Gmsh API 4.15.2, found " *
    gmsh.GMSH_API_VERSION)

const _QUERY_COORDS=Float64[0 1 0 0;
                            0 0 1 0;
                            0 0 0 1]

function _query_mesh()
    return Mesh(
        _QUERY_COORDS;
        segs=reshape(Int32[1,2],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[1,2,3,4],4,1),
        seg_tag=Int32[11],tri_tag=Int32[22],tet_tag=Int32[33])
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

function _normalized_nodes(values::Vector{UInt64})
    mapping=Dict(UInt64(10)=>UInt64(1),UInt64(20)=>UInt64(2),
                 UInt64(30)=>UInt64(3),UInt64(40)=>UInt64(4))
    all(tag->haskey(mapping,tag),values) || error(
        "Gmsh connectivity contains an unknown fixture node")
    return UInt64[mapping[tag] for tag in values]
end

try
    gmsh.initialize(String[],false)
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.model.add("mesh-data-query-empty")
    isempty(gmsh.model.mesh.getNodes()[1]) || error(
        "Gmsh empty model unexpectedly owns nodes")
    isempty(gmsh.model.mesh.getElements()[1]) || error(
        "Gmsh empty model unexpectedly owns elements")
    gmsh.model.mesh.getMaxNodeTag()==0 || error(
        "Gmsh empty model has a nonzero maximum node tag")
    gmsh.model.mesh.getMaxElementTag()==0 || error(
        "Gmsh empty model has a nonzero maximum element tag")

    gmsh.clear()
    gmsh.model.add("mesh-data-query-fixture")
    for (dimension,entity) in ((1,101),(2,201),(3,301))
        gmsh.model.addDiscreteEntity(dimension,entity)
    end
    gmsh.model.mesh.addNodes(
        3,301,UInt64[10,20,30,40],collect(vec(_QUERY_COORDS)))
    gmsh.model.mesh.addElementsByType(
        101,1,UInt64[100],UInt64[10,20])
    gmsh.model.mesh.addElementsByType(
        201,2,UInt64[300],UInt64[10,20,30])
    gmsh.model.mesh.addElementsByType(
        301,4,UInt64[200],UInt64[10,20,30,40])

    gmsh_node_tags,gmsh_coordinates,gmsh_parameters=
        gmsh.model.mesh.getNodes()
    gmsh_node_tags==UInt64[10,20,30,40] || error(
        "Gmsh did not preserve the explicit node tags")
    gmsh_coordinates==collect(vec(_QUERY_COORDS)) || error(
        "Gmsh changed the fixture coordinates")
    isempty(gmsh_parameters) || error(
        "Gmsh discrete-volume nodes unexpectedly have parametric coordinates")

    gmsh_types,gmsh_element_tags,gmsh_element_nodes=
        gmsh.model.mesh.getElements()
    gmsh_types==Int32[1,2,4] || error(
        "Gmsh returned element types $gmsh_types")
    gmsh_element_tags==[UInt64[100],UInt64[300],UInt64[200]] || error(
        "Gmsh did not preserve the explicit element tags")

    Tessella.API.initialize()
    refined_crc=try
        _rejects_argument(()->Tessella.API.mesh.get_nodes()) || error(
            "Tessella returned bulk nodes without a cached mesh")
        _rejects_argument(()->Tessella.API.mesh.get_elements()) || error(
            "Tessella returned bulk elements without a cached mesh")

        # Exact differential setup: install the same validated simplex fixture in
        # the session cache. The public API intentionally has no add-nodes mutation
        # yet; all operations under test below are public and read-only.
        fixture=_query_mesh()
        lock(Tessella.API.STATE_LOCK) do
            Tessella.API.LAST_MESH[]=Tessella.API._copy_mesh(fixture)
        end

        tessella_node_tags,tessella_coordinates,tessella_parameters=
            Tessella.API.mesh.get_nodes()
        tessella_node_tags==UInt64[1,2,3,4] || error(
            "Tessella dense node tags are not 1:4")
        tessella_coordinates==gmsh_coordinates || error(
            "Tessella/Gmsh bulk coordinates differ")
        isempty(tessella_parameters) || error(
            "Tessella simplex cache returned parametric coordinates")

        tessella_types,tessella_element_tags,tessella_element_nodes=
            Tessella.API.mesh.get_elements()
        tessella_types==gmsh_types || error(
            "Tessella/Gmsh bulk element types differ")
        normalized_gmsh_nodes=_normalized_nodes.(gmsh_element_nodes)
        tessella_element_nodes==normalized_gmsh_nodes || error(
            "Tessella/Gmsh normalized element connectivity differs")
        tessella_element_tags==[UInt64[1],UInt64[2],UInt64[3]] || error(
            "Tessella dense element tags are not globally consecutive")

        for (dimension,expected_type) in ((0,Int32[]),(1,Int32[1]),
                                          (2,Int32[2]),(3,Int32[4]))
            Tessella.API.mesh.get_element_types(dimension)==expected_type || error(
                "Tessella dimension-$dimension type filter differs")
            gmsh.model.mesh.getElementTypes(dimension,-1)==expected_type || error(
                "Gmsh dimension-$dimension type filter differs")
        end
        for element_type in (1,2,4)
            tessella_tags,tessella_nodes=
                Tessella.API.mesh.get_elements_by_type(element_type)
            gmsh_tags,gmsh_nodes=
                gmsh.model.mesh.getElementsByType(element_type)
            tessella_nodes==_normalized_nodes(gmsh_nodes) || error(
                "type-$element_type normalized connectivity differs")
            length(tessella_tags)==length(gmsh_tags)==1 || error(
                "type-$element_type element count differs")
        end
        Tessella.API.mesh.get_elements_by_type(3)==(UInt64[],UInt64[]) || error(
            "Tessella returned absent quadrangles")
        gmsh.model.mesh.getElementsByType(3)==(UInt64[],UInt64[]) || error(
            "Gmsh returned absent quadrangles")
        _rejects_argument(
            ()->Tessella.API.mesh.get_elements_by_type(999)) || error(
                "Tessella accepted an unknown MSH type")
        gmsh_unknown=try
            gmsh.model.mesh.getElementsByType(999)
            false
        catch
            true
        end
        gmsh_unknown || error("Gmsh accepted unknown MSH type 999")

        Tessella.API.mesh.get_max_node_tag()==UInt64(4) || error(
            "Tessella maximum dense node tag is not 4")
        Tessella.API.mesh.get_max_element_tag()==UInt64(3) || error(
            "Tessella maximum dense element tag is not 3")
        gmsh.model.mesh.getMaxNodeTag()==UInt64(40) || error(
            "Gmsh maximum explicit node tag is not 40")
        gmsh.model.mesh.getMaxElementTag()==UInt64(300) || error(
            "Gmsh maximum explicit element tag is not 300")

        _rejects_argument(()->Tessella.API.mesh.get_nodes(3,-1)) || error(
            "Tessella fabricated node classification")
        _rejects_argument(()->Tessella.API.mesh.get_elements(3,301)) || error(
            "Tessella fabricated element classification")
        _rejects_argument(
            ()->Tessella.API.mesh.get_elements_by_type(4,-1,1,2)) || error(
                "Tessella accepted nondefault Julia task partitioning")

        returned_nodes=Tessella.API.mesh.get_nodes()
        returned_elements=Tessella.API.mesh.get_elements()
        returned_nodes[2][1]=99
        returned_elements[2][1][1]=99
        Tessella.API.mesh.get_nodes()[2]==collect(vec(_QUERY_COORDS)) || error(
            "bulk node output aliases the cache")
        Tessella.API.mesh.get_elements()[2][1]==UInt64[1] || error(
            "bulk element output aliases the cache")

        refined=Tessella.API.mesh.refine()
        Tessella.API.mesh.get_max_node_tag()==UInt64(10) || error(
            "refined cache maximum node tag is not 10")
        Tessella.API.mesh.get_max_element_tag()==UInt64(14) || error(
            "refined cache maximum element tag is not 14")
        lengths=length.(last(Tessella.API.mesh.get_elements()))
        lengths==[4,12,32] || error(
            "refined cache connectivity lengths are $lengths")
        crc=Tessella.mesh_crc(refined)
        crc.sha=="db9a1713d1174be1035ef3e9d6380a01ed419797a91ded9a2b8508d0b038f031" ||
            error("refined query fixture checksum changed to $(crc.sha)")
        Tessella.API.mesh.clear()
        _rejects_argument(()->Tessella.API.mesh.get_nodes()) || error(
            "cleared cache retained bulk nodes")
        crc
    finally
        Tessella.API.finalize()
    end

    println("mesh-data-query differential: Gmsh ",gmsh.GMSH_API_VERSION,
            ", types=1,2,4 nodes=4 elements=3 connectivity_entries=9 ",
            "dense_max_tags=4/3 explicit_max_tags=40/300 ",
            "refined_sha=",refined_crc.sha,
            " bounded=no-mesh-and-classification-blockers")
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
