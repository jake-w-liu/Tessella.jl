# Differential oracle for read-only bulk node and element access. This uses the
# locally installed Gmsh 4.15.2 Julia API and never starts the GUI.
using Tessella
using SHA

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

function _query_sha(groups)
    stream=IOBuffer()
    for values in groups
        marker=eltype(values)===UInt64 ? UInt8(1) : UInt8(2)
        write(stream,marker)
        write(stream,htol(UInt64(length(values))))
        for value in values
            bits=value isa Float64 ? reinterpret(UInt64,value) : UInt64(value)
            write(stream,htol(bits))
        end
    end
    return bytes2hex(SHA.sha256(take!(stream)))
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
    refined_crc,derived_sha=try
        _rejects_argument(()->Tessella.API.mesh.get_nodes()) || error(
            "Tessella returned bulk nodes without a cached mesh")
        _rejects_argument(()->Tessella.API.mesh.get_elements()) || error(
            "Tessella returned bulk elements without a cached mesh")
        _rejects_argument(
            ()->Tessella.API.mesh.get_nodes_by_element_type(4)) || error(
            "Tessella returned type nodes without a cached mesh")
        _rejects_argument(
            ()->Tessella.API.mesh.get_barycenters(4,-1,false,false)) || error(
            "Tessella returned barycenters without a cached mesh")

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

        tessella_type_nodes=Dict{Int,Tuple}()
        gmsh_type_nodes=Dict{Int,Tuple}()
        for element_type in (1,2,4)
            tessella_nodes=Tessella.API.mesh.get_nodes_by_element_type(
                element_type)
            gmsh_nodes=gmsh.model.mesh.getNodesByElementType(element_type)
            tessella_nodes[1]==_normalized_nodes(gmsh_nodes[1]) || error(
                "type-$element_type per-element node tags differ")
            tessella_nodes[2]==gmsh_nodes[2] || error(
                "type-$element_type per-element node coordinates differ")
            isempty(tessella_nodes[3]) && isempty(gmsh_nodes[3]) || error(
                "type-$element_type parametric-coordinate behavior differs")
            tessella_type_nodes[element_type]=tessella_nodes
            gmsh_type_nodes[element_type]=gmsh_nodes
        end
        Tessella.API.mesh.get_nodes_by_element_type(3)==
            (UInt64[],Float64[],Float64[]) || error(
            "Tessella returned nodes for absent quadrangles")
        gmsh.model.mesh.getNodesByElementType(3)==
            (UInt64[],Float64[],Float64[]) || error(
            "Gmsh returned nodes for absent quadrangles")

        tessella_barycenters=Dict{Tuple{Int,Bool,Bool},Vector{Float64}}()
        gmsh_barycenters=Dict{Tuple{Int,Bool,Bool},Vector{Float64}}()
        for element_type in (1,2,4),fast in (false,true),
            primary in (false,true)
            key=(element_type,fast,primary)
            tessella_values=Tessella.API.mesh.get_barycenters(
                element_type,-1,fast,primary)
            gmsh_values=gmsh.model.mesh.getBarycenters(
                element_type,-1,fast,primary)
            tessella_values==gmsh_values || error(
                "type-$element_type fast=$fast primary=$primary " *
                "barycenters differ")
            tessella_barycenters[key]=tessella_values
            gmsh_barycenters[key]=gmsh_values
        end

        tessella_edges=Dict{Tuple{Int,Bool},Vector{UInt64}}()
        gmsh_edges=Dict{Tuple{Int,Bool},Vector{UInt64}}()
        for element_type in (1,2,4),primary in (false,true)
            key=(element_type,primary)
            tessella_values=Tessella.API.mesh.get_element_edge_nodes(
                element_type,-1,primary)
            gmsh_values=gmsh.model.mesh.getElementEdgeNodes(
                element_type,-1,primary)
            tessella_values==_normalized_nodes(gmsh_values) || error(
                "type-$element_type primary=$primary edge nodes differ")
            tessella_edges[key]=tessella_values
            gmsh_edges[key]=gmsh_values
        end

        tessella_faces=Dict{Tuple{Int,Int,Bool},Vector{UInt64}}()
        gmsh_faces=Dict{Tuple{Int,Int,Bool},Vector{UInt64}}()
        for element_type in (1,2,4),face_type in (3,4),
            primary in (false,true)
            key=(element_type,face_type,primary)
            tessella_values=Tessella.API.mesh.get_element_face_nodes(
                element_type,face_type,-1,primary)
            gmsh_values=gmsh.model.mesh.getElementFaceNodes(
                element_type,face_type,-1,primary)
            tessella_values==_normalized_nodes(gmsh_values) || error(
                "type-$element_type face-$face_type primary=$primary " *
                "face nodes differ")
            tessella_faces[key]=tessella_values
            gmsh_faces[key]=gmsh_values
        end

        tessella_derived_sha=_query_sha((
            tessella_type_nodes[1][1],tessella_type_nodes[1][2],
            tessella_type_nodes[2][1],tessella_type_nodes[2][2],
            tessella_type_nodes[4][1],tessella_type_nodes[4][2],
            tessella_barycenters[(1,false,false)],
            tessella_barycenters[(2,false,true)],
            tessella_barycenters[(4,false,false)],
            tessella_edges[(1,false)],tessella_edges[(2,true)],
            tessella_edges[(4,false)],tessella_faces[(2,3,false)],
            tessella_faces[(4,3,true)]))
        gmsh_derived_sha=_query_sha((
            _normalized_nodes(gmsh_type_nodes[1][1]),gmsh_type_nodes[1][2],
            _normalized_nodes(gmsh_type_nodes[2][1]),gmsh_type_nodes[2][2],
            _normalized_nodes(gmsh_type_nodes[4][1]),gmsh_type_nodes[4][2],
            gmsh_barycenters[(1,false,false)],
            gmsh_barycenters[(2,false,true)],
            gmsh_barycenters[(4,false,false)],
            _normalized_nodes(gmsh_edges[(1,false)]),
            _normalized_nodes(gmsh_edges[(2,true)]),
            _normalized_nodes(gmsh_edges[(4,false)]),
            _normalized_nodes(gmsh_faces[(2,3,false)]),
            _normalized_nodes(gmsh_faces[(4,3,true)])))
        tessella_derived_sha==gmsh_derived_sha==
            "04e09b72ebf17bdc7ab2f9f96da2927c5a6e892e5c2d98c9ddbb8313dc4cab13" ||
            error("derived query checksum mismatch")

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
            ()->Tessella.API.mesh.get_nodes_by_element_type(4,301)) || error(
            "Tessella fabricated type-node classification")
        _rejects_argument(
            ()->Tessella.API.mesh.get_elements_by_type(4,-1,1,2)) || error(
            "Tessella accepted nondefault Julia task partitioning")
        _rejects_argument(
            ()->Tessella.API.mesh.get_element_edge_nodes(4,-1,false,1,2)) ||
            error("Tessella accepted nondefault edge-node task partitioning")
        _rejects_argument(
            ()->Tessella.API.mesh.get_element_face_nodes(4,2)) || error(
            "Tessella accepted a non-face node count")
        isempty(gmsh.model.mesh.getElementFaceNodes(4,2)) || error(
            "Gmsh returned nodes for face type 2")
        _rejects_argument(
            ()->Tessella.API.mesh.get_nodes_by_element_type(34)) || error(
            "Tessella treated special MSH type 34 as fixed-node data")
        gmsh.model.mesh.getNodesByElementType(34)==
            (UInt64[],Float64[],Float64[]) || error(
            "Gmsh returned nodes for absent special MSH type 34")

        returned_nodes=Tessella.API.mesh.get_nodes()
        returned_elements=Tessella.API.mesh.get_elements()
        returned_edges=Tessella.API.mesh.get_element_edge_nodes(4)
        returned_nodes[2][1]=99
        returned_elements[2][1][1]=99
        returned_edges[1]=99
        Tessella.API.mesh.get_nodes()[2]==collect(vec(_QUERY_COORDS)) || error(
            "bulk node output aliases the cache")
        Tessella.API.mesh.get_elements()[2][1]==UInt64[1] || error(
            "bulk element output aliases the cache")
        Tessella.API.mesh.get_element_edge_nodes(4)[1]==UInt64(1) || error(
            "edge-node output aliases the cache")

        maximum=floatmax(Float64)
        overflow_fixture=Mesh(
            Float64[maximum maximum;0 0;0 0];
            segs=reshape(Int32[1,2],2,1))
        lock(Tessella.API.STATE_LOCK) do
            Tessella.API.LAST_MESH[]=Tessella.API._copy_mesh(overflow_fixture)
        end
        Tessella.API.mesh.get_barycenters(1,-1,false,false)==
            Float64[maximum,0,0] || error(
            "Tessella overflow-safe barycenter changed")
        _rejects_argument(
            ()->Tessella.API.mesh.get_barycenters(1,-1,true,false)) || error(
            "Tessella returned a nonfinite fast barycenter")
        lock(Tessella.API.STATE_LOCK) do
            Tessella.API.LAST_MESH[]=Tessella.API._copy_mesh(fixture)
        end

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
        crc,tessella_derived_sha
    finally
        Tessella.API.finalize()
    end

    gmsh.model.add("mesh-data-query-overflow")
    gmsh.model.addDiscreteEntity(1,901)
    maximum=floatmax(Float64)
    gmsh.model.mesh.addNodes(
        1,901,UInt64[501,502],Float64[maximum,0,0,maximum,0,0])
    gmsh.model.mesh.addElementsByType(
        901,1,UInt64[601],UInt64[501,502])
    gmsh_normal=gmsh.model.mesh.getBarycenters(1,-1,false,false)
    gmsh_fast=gmsh.model.mesh.getBarycenters(1,-1,true,false)
    isinf(gmsh_normal[1]) && isinf(gmsh_fast[1]) || error(
        "Gmsh overflow barycenter behavior changed")

    println("mesh-data-query differential: Gmsh ",gmsh.GMSH_API_VERSION,
            ", types=1,2,4 nodes=4 elements=3 connectivity_entries=9 ",
            "dense_max_tags=4/3 explicit_max_tags=40/300 ",
            "derived_sha=",derived_sha," ",
            "refined_sha=",refined_crc.sha,
            " bounded=no-mesh/classification/task/special-type/face-count ",
            "blockers and finite-barycenter contract")
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
