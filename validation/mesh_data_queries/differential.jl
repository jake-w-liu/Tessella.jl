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
        # Gmsh 4.15.2 returns empty arrays on a model with no mesh; Tessella
        # matches that contract.
        Tessella.API.mesh.get_nodes()==(UInt64[],Float64[],Float64[]) ||
            error("Tessella empty-cache bulk nodes changed")
        Tessella.API.mesh.get_elements()==
            (Int32[],Vector{UInt64}[],Vector{UInt64}[]) || error(
            "Tessella empty-cache bulk elements changed")
        Tessella.API.mesh.get_nodes_by_element_type(4)==
            (UInt64[],Float64[],Float64[]) || error(
            "Tessella empty-cache type nodes changed")
        Tessella.API.mesh.get_barycenters(4,-1,false,false)==Float64[] ||
            error("Tessella empty-cache barycenters changed")

        # Exact differential setup: install the same validated simplex fixture in
        # the session cache; all operations under test below are read-only.
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
        for (task,num_tasks) in ((-1,1),(0,0),(0,-1))
            _rejects_argument(()->Tessella.API.mesh.get_elements_by_type(
                4,-1,task,num_tasks)) || error(
                "Tessella accepted task=$task num_tasks=$num_tasks")
            _rejects_argument(()->Tessella.API.mesh.get_element_edge_nodes(
                4,-1,false,task,num_tasks)) || error(
                "Tessella accepted edge task=$task num_tasks=$num_tasks")
            _rejects_argument(()->Tessella.API.mesh.get_barycenters(
                4,-1,false,false,task,num_tasks)) || error(
                "Tessella accepted barycenter task=$task num_tasks=$num_tasks")
        end
        _rejects_argument(()->Tessella.API.mesh.get_elements_by_type(
            4,-1,true,1)) || error(
            "Tessella accepted a Bool task")
        # With one cached element per type, task 1 of 2 owns the whole
        # 0-based block [0,1) while task 0 of 2 is the silently-empty range.
        Tessella.API.mesh.get_barycenters(4,-1,false,false,1,2)==
            tessella_barycenters[(4,false,false)] || error(
            "single-tet task-1-of-2 barycenters differ")
        gmsh.model.mesh.getBarycenters(4,-1,false,false,1,2)==
            gmsh_barycenters[(4,false,false)] || error(
            "Gmsh single-tet task-1-of-2 barycenters differ")
        isempty(Tessella.API.mesh.get_barycenters(
            4,-1,false,false,0,2)) || error(
            "single-tet task-0-of-2 barycenters are not empty")
        all(iszero,gmsh.model.mesh.getBarycenters(4,-1,false,false,0,2)) ||
            error("Gmsh single-tet task-0-of-2 barycenters are not zero")

        # Five-segment chain: differential cross-check of the contiguous
        # Gmsh block formula begin=(task*count)÷numTasks,
        # end=((task+1)*count)÷numTasks against the pinned implementation.
        # Gmsh returns full-size zero-padded vectors; Tessella returns the
        # computed slice, so the slice must equal the padded block section
        # with zeros exactly outside it (checked on strictly-positive tags).
        chain_model="mesh-data-query-partitions"
        gmsh.model.add(chain_model)
        gmsh.model.addDiscreteEntity(1,701)
        chain_coordinates=Float64[0 1 2 3 4 5;0 0 0 0 0 0;0 0 0 0 0 0]
        gmsh.model.mesh.addNodes(
            1,701,UInt64[11,12,13,14,15,16],collect(vec(chain_coordinates)))
        gmsh.model.mesh.addElementsByType(
            701,1,UInt64[71,72,73,74,75],
            UInt64[11,12,12,13,13,14,14,15,15,16])
        chain=Mesh(chain_coordinates;
                    segs=Int32[1 2 3 4 5;2 3 4 5 6])
        lock(Tessella.API.STATE_LOCK) do
            Tessella.API.LAST_MESH[]=Tessella.API._copy_mesh(chain)
        end
        chain_block(task,num_tasks)=let
            first_element=(task*5)÷num_tasks+1
            last_element=min(((task+1)*5)÷num_tasks,5)
            first_element:last_element
        end
        _chain_nodes(values)=UInt64[value-UInt64(10) for value in values]
        for (task,num_tasks) in
            ((0,1),(0,2),(1,2),(0,3),(1,3),(2,3),(2,2),(0,5),(4,5),(5,5))
            block=chain_block(task,num_tasks)
            tessella_tags,tessella_nodes=
                Tessella.API.mesh.get_elements_by_type(1,-1,task,num_tasks)
            gmsh_tags,gmsh_nodes=
                gmsh.model.mesh.getElementsByType(1,-1,task,num_tasks)
            tessella_bary=Tessella.API.mesh.get_barycenters(
                1,-1,false,false,task,num_tasks)
            gmsh_bary=gmsh.model.mesh.getBarycenters(
                1,-1,false,false,task,num_tasks)
            tessella_edge=Tessella.API.mesh.get_element_edge_nodes(
                1,-1,false,task,num_tasks)
            gmsh_edge=gmsh.model.mesh.getElementEdgeNodes(
                1,-1,false,task,num_tasks)
            length(gmsh_tags)==5 && length(gmsh_nodes)==10 || error(
                "Gmsh task=$task/$num_tasks tag layout changed")
            length(gmsh_bary)==15 && length(gmsh_edge)==10 || error(
                "Gmsh task=$task/$num_tasks value layout changed")
            if isempty(block)
                isempty(tessella_tags) && isempty(tessella_nodes) || error(
                    "task=$task/$num_tasks element slice is not empty")
                isempty(tessella_bary) && isempty(tessella_edge) || error(
                    "task=$task/$num_tasks value slice is not empty")
                all(iszero,gmsh_tags) && all(iszero,gmsh_nodes) || error(
                    "Gmsh task=$task/$num_tasks tags are not zero")
                all(iszero,gmsh_bary) && all(iszero,gmsh_edge) || error(
                    "Gmsh task=$task/$num_tasks values are not zero")
            else
                first_element,last_element=first(block),last(block)
                node_range=(first_element-1)*2+1:last_element*2
                bary_range=(first_element-1)*3+1:last_element*3
                before_tags=first_element>1 ?
                    gmsh_tags[1:first_element-1] : UInt64[]
                after_tags=last_element<5 ?
                    gmsh_tags[last_element+1:5] : UInt64[]
                all(iszero,before_tags) && all(iszero,after_tags) || error(
                    "Gmsh task=$task/$num_tasks tag padding moved")
                gmsh_tags[block] .- UInt64(70)==tessella_tags || error(
                    "task=$task/$num_tasks element tags differ")
                _chain_nodes(gmsh_nodes[node_range])==
                    tessella_nodes || error(
                    "task=$task/$num_tasks element nodes differ")
                gmsh_bary[bary_range]==tessella_bary || error(
                    "task=$task/$num_tasks barycenters differ")
                before_edge=first_element>1 ?
                    gmsh_edge[1:(first_element-1)*2] : UInt64[]
                after_edge=last_element<5 ?
                    gmsh_edge[last_element*2+1:10] : UInt64[]
                all(iszero,before_edge) && all(iszero,after_edge) || error(
                    "Gmsh task=$task/$num_tasks edge padding moved")
                _chain_nodes(gmsh_edge[node_range])==
                    tessella_edge || error(
                    "task=$task/$num_tasks edge nodes differ")
            end
        end
        for num_tasks in (2,3)
            union_tags=UInt64[]
            union_nodes=UInt64[]
            for task in 0:num_tasks-1
                slice_tags,slice_nodes=
                    Tessella.API.mesh.get_elements_by_type(
                        1,-1,task,num_tasks)
                append!(union_tags,slice_tags)
                append!(union_nodes,slice_nodes)
            end
            union_tags==UInt64[1,2,3,4,5] || error(
                "five-segment task union tags differ")
            union_nodes==UInt64[1,2,2,3,3,4,4,5,5,6] || error(
                "five-segment task union nodes differ")
        end
        lock(Tessella.API.STATE_LOCK) do
            Tessella.API.LAST_MESH[]=Tessella.API._copy_mesh(fixture)
        end
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
        Tessella.API.mesh.get_nodes()==(UInt64[],Float64[],Float64[]) ||
            error("cleared cache retained bulk nodes")
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

    # ── Entity-filtered queries on generated meshes ──────────────────────────
    # The two engines triangulate differently, so exact node/element sets cannot
    # be compared. What must agree: entity-existence errors, dim=-1 tag
    # discarding, the transitive include_boundary closure against
    # model.getBoundary, entity-first ordering, per-type entity filters, and
    # task partitioning over the filtered subset.
    gmsh.model.add("mesh-data-query-entity-surface")
    for (point,(x,y)) in enumerate(
            ((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)))
        gmsh.model.geo.addPoint(x,y,0.0,0.35,point)
    end
    for (curve,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        gmsh.model.geo.addLine(a,b,curve)
    end
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.addPoint(0.5,0.5,0.0,0.2,5)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0,[5],2,1)
    gmsh.model.addPhysicalGroup(2,[1],10)
    gmsh.model.addPhysicalGroup(1,[3,1],11)
    gmsh.model.mesh.generate(2)

    Tessella.API.initialize()
    entity_pairs=try
        for (point,(x,y)) in enumerate(
                ((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)))
            Tessella.API.model.add_point(x,y,0.0;tag=point,meshSize=0.35)
        end
        for (curve,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            Tessella.API.model.add_line(a,b;tag=curve)
        end
        Tessella.API.model.add_curve_loop([1,2,3,4];tag=1)
        Tessella.API.model.add_plane_surface([1];tag=1)
        Tessella.API.model.add_point(0.5,0.5,0.0;tag=5,meshSize=0.2)
        Tessella.API.model.embed(0,[5],2,1)
        Tessella.API.model.add_physical_group(2,[1];tag=10)
        Tessella.API.model.add_physical_group(1,[3,1];tag=11)
        Tessella.API.option("Mesh.MeshSizeMax",0.35)
        Tessella.API.mesh.generate(2)

        gmsh_surface_nodes=Set(gmsh.model.mesh.getNodes(2,1,false)[1])
        tessella_surface_nodes=Set(
            Tessella.API.mesh.get_nodes(2,1,false)[1])
        gmsh_surface_with_boundary=gmsh.model.mesh.getNodes(2,1,true)[1]
        tessella_surface_with_boundary=
            Tessella.API.mesh.get_nodes(2,1,true)[1]
        # include_boundary emits the entity's own nodes first.
        gmsh_surface_with_boundary[1:length(gmsh_surface_nodes)] ==
            gmsh.model.mesh.getNodes(2,1,false)[1] || error(
            "Gmsh surface include_boundary does not emit own nodes first")
        tessella_surface_with_boundary[1:length(tessella_surface_nodes)] ==
            Tessella.API.mesh.get_nodes(2,1,false)[1] || error(
            "Tessella surface include_boundary does not emit own nodes first")
        # The closure equals the transitive boundary node union in both;
        # recursive=false walks direct boundaries level by level because Gmsh's
        # recursive=true reports only the lowest-dimension entities.
        function transitive_nodes(get_boundary,get_nodes,entity)
            result=Set(get_nodes(entity...))
            queue=[entity]
            seen=Set{Tuple{Int,Int}}()
            while !isempty(queue)
                current=popfirst!(queue)
                for boundary in get_boundary(current)
                    boundary in seen && continue
                    push!(seen,boundary)
                    union!(result,get_nodes(boundary...))
                    push!(queue,boundary)
                end
            end
            return result
        end
        gmsh_expected=transitive_nodes(
            entity->Tuple{Int,Int}[Tuple{Int,Int}(pair) for pair in
                gmsh.model.getBoundary([entity],false,false,false)],
            (dim,tag)->gmsh.model.mesh.getNodes(dim,tag,false)[1],
            (2,1))
        tessella_expected=transitive_nodes(
            entity->Tuple{Int,Int}[Tuple{Int,Int}(pair) for pair in
                Tessella.API.model.get_boundary([entity],false,false,false)],
            (dim,tag)->Tessella.API.mesh.get_nodes(dim,tag,false)[1],
            (2,1))
        Set(gmsh_surface_with_boundary)==gmsh_expected || error(
            "Gmsh include_boundary is not the transitive boundary closure")
        Set(tessella_surface_with_boundary)==tessella_expected || error(
            "Tessella include_boundary is not the transitive boundary closure")
        # Point entities classify exactly one node on each corner in both.
        for point in 1:4
            length(gmsh.model.mesh.getNodes(0,point,false)[1])==1 || error(
                "Gmsh point $point does not classify one node")
            length(Tessella.API.mesh.get_nodes(0,point,false)[1])==1 || error(
                "Tessella point $point does not classify one node")
        end
        # dim=-1 discards the tag in both engines.
        gmsh.model.mesh.getNodes(-1,999,false)[1]==
            gmsh.model.mesh.getNodes(-1,-1,false)[1] || error(
            "Gmsh dim=-1 did not discard the tag")
        Tessella.API.mesh.get_nodes(-1,999,false)[1]==
            Tessella.API.mesh.get_nodes(-1,-1,false)[1] || error(
            "Tessella dim=-1 did not discard the tag")
        # The single surface owns every cached triangle in both engines.
        gmsh_all=gmsh.model.mesh.getElements(2,-1)
        gmsh_filtered=gmsh.model.mesh.getElements(2,1)
        tessella_all=Tessella.API.mesh.get_elements(2,-1)
        tessella_filtered=Tessella.API.mesh.get_elements(2,1)
        gmsh_all[1]==gmsh_filtered[1]==Int32[2] || error(
            "Gmsh surface element-type filter changed")
        tessella_all[1]==tessella_filtered[1]==Int32[2] || error(
            "Tessella surface element-type filter changed")
        length(gmsh_filtered[2][1])==length(gmsh_all[2][1]) || error(
            "Gmsh surface element filter is not exhaustive")
        length(tessella_filtered[2][1])==length(tessella_all[2][1]) ||
            error("Tessella surface element filter is not exhaustive")
        # Parametric coordinates ride the queried entity's parametrization in
        # both engines: two entries per node on a Plane, one per node on a
        # Line, none for Points, Volumes, or all-dimension queries. Each
        # emitted parameter vector must re-evaluate to the node's physical
        # coordinates through the entity's forward map.
        for (engine,tags,coords,par,evaluate) in (
            ("Gmsh",gmsh.model.mesh.getNodes(2,1,true,true)...,
             par->gmsh.model.getValue(2,1,collect(par))),
            ("Tessella",Tessella.API.mesh.get_nodes(2,1,true,true)...,
             par->Tessella.API.model.get_value(2,1,collect(par))))
            length(par)==2length(tags) || error(
                "$engine surface parametric width mismatch")
            Base.maximum(index->Base.maximum(abs.(
                evaluate(par[2index-1:2index]).-coords[3index-2:3index])),
                eachindex(tags))<1e-9 || error(
                "$engine surface parameters do not re-evaluate to nodes")
        end
        for (engine,tags,coords,par,evaluate) in (
            ("Gmsh",gmsh.model.mesh.getNodes(1,1,true,true)...,
             par->gmsh.model.getValue(1,1,[par])),
            ("Tessella",Tessella.API.mesh.get_nodes(1,1,true,true)...,
             par->Tessella.API.model.get_value(1,1,[par])))
            length(par)==length(tags) || error(
                "$engine curve parametric width mismatch")
            Base.maximum(index->Base.maximum(abs.(
                evaluate(par[index]).-coords[3index-2:3index])),
                eachindex(tags))<1e-9 || error(
                "$engine curve parameters do not re-evaluate to nodes")
        end
        for (engine,empty_par) in (
            ("Gmsh",gmsh.model.mesh.getNodes(0,1,true,true)[3]),
            ("Tessella",Tessella.API.mesh.get_nodes(0,1,true,true)[3]),
            ("Gmsh",gmsh.model.mesh.getNodes(-1,-1,false,true)[3]),
            ("Tessella",Tessella.API.mesh.get_nodes(-1,-1,false,true)[3]),
            ("Gmsh",gmsh.model.mesh.getNodes(2,1,true,false)[3]),
            ("Tessella",Tessella.API.mesh.get_nodes(2,1,true,false)[3]))
            isempty(empty_par) || error(
                "$engine emitted parameters for an unparametrized query")
        end
        # getNodesByElementType packs each repeated node's parameters on its
        # owning entity: surface owners contribute two, curve owners one,
        # point owners zero, in entry order.
        for (engine,tags,coords,par,owners) in (
            ("Gmsh",gmsh.model.mesh.getNodesByElementType(2,-1,true)...,
             let owner=Dict{eltype(gmsh.model.mesh.getNodes(-1,-1)[1]),
                            Tuple{Int,Int}}()
                 for point in 1:5,
                     node in gmsh.model.mesh.getNodes(0,point,false)[1]
                     owner[node]=(0,point)
                 end
                 for curve in 1:4,
                     node in gmsh.model.mesh.getNodes(1,curve,false)[1]
                     owner[node]=(1,curve)
                 end
                 for node in gmsh.model.mesh.getNodes(2,1,false)[1]
                     owner[node]=(2,1)
                 end
                 owner
             end),
            ("Tessella",
             Tessella.API.mesh.get_nodes_by_element_type(2,-1,true)...,
             let owner=Dict{UInt64,Tuple{Int,Int}}()
                 for point in 1:5,
                     node in Tessella.API.mesh.get_nodes(0,point)[1]
                     owner[node]=(0,point)
                 end
                 for curve in 1:4,
                     node in Tessella.API.mesh.get_nodes(1,curve)[1]
                     owner[node]=(1,curve)
                 end
                 for node in Tessella.API.mesh.get_nodes(2,1)[1]
                     owner[node]=(2,1)
                 end
                 owner
             end))
            expected=sum(node->owners[node][1] in (1,2) ?
                             owners[node][1] : 0,tags)
            length(par)==expected || error(
                "$engine by-element-type parametric width mismatch")
            position=1
            worst=0.0
            for (k,node) in enumerate(tags)
                (dim,entity)=owners[node]
                width=dim in (1,2) ? dim : 0
                if width>0
                    evaluated=engine=="Gmsh" ?
                        gmsh.model.getValue(dim,entity,
                            collect(par[position:position+width-1])) :
                        Tessella.API.model.get_value(dim,entity,
                            collect(par[position:position+width-1]))
                    worst=max(worst,Base.maximum(abs.(
                        evaluated.-coords[3k-2:3k])))
                end
                position+=width
            end
            worst<1e-9 || error(
                "$engine by-element-type parameters do not re-evaluate")
        end
        # getNode reports each node's owning entity and reparametrizes the
        # node on it: one `u` for a Line owner, `(u, v)` for a Plane owner,
        # none for Point or Volume owners; the owner must be the entity
        # whose strict (boundary-excluded) node list contains the node.
        for (engine,all_tags,query,listed,evaluate) in (
            ("Gmsh",gmsh.model.mesh.getNodes(-1,-1)[1],
             t->gmsh.model.mesh.getNode(t),
             (d,t)->Set(gmsh.model.mesh.getNodes(d,t,false)[1]),
             (d,e,p)->gmsh.model.getValue(d,e,collect(p))),
            ("Tessella",Tessella.API.mesh.get_nodes(-1,-1)[1],
             t->Tessella.API.mesh.get_node(t),
             (d,t)->Set(Tessella.API.mesh.get_nodes(d,t)[1]),
             (d,e,p)->Tessella.API.model.get_value(d,e,p)))
            for node in all_tags
                coord,par,dim,entity=query(node)
                dim in 0:3 || error(
                    "$engine getNode returned owner dim $dim")
                node in listed(dim,entity) || error(
                    "$engine getNode owner ($dim, $entity) omits the node")
                length(par)==(dim in (1,2) ? dim : 0) || error(
                    "$engine getNode parametric width mismatch")
                if dim in (1,2)
                    Base.maximum(abs.(
                        evaluate(dim,entity,par).-coord))<1e-9 || error(
                        "$engine getNode parameters do not re-evaluate")
                end
            end
        end
        let thrown=false
            try
                Tessella.API.mesh.get_node(
                    length(Tessella.API.mesh.get_nodes(-1,-1)[1])+1)
            catch err
                thrown=err isa ArgumentError
            end
            thrown || error("Tessella getNode accepted an unknown tag")
        end
        # Physical-group queries emit each member entity's nodes plus its
        # transitive boundary and transitively embedded entities' nodes as one
        # sorted unique set; unknown or dimension-mismatched groups are empty.
        function group_expected(get_boundary,get_nodes,get_embedded,members)
            result=Set{UInt64}()
            queue=Tuple{Int,Int}[members...]
            seen=Set{Tuple{Int,Int}}()
            while !isempty(queue)
                current=popfirst!(queue)
                current in seen && continue
                push!(seen,current)
                union!(result,get_nodes(current...))
                append!(queue,get_boundary(current))
                append!(queue,get_embedded(current...))
            end
            return result
        end
        for (engine,query,get_boundary,get_nodes,get_embedded) in (
            ("Gmsh",
             (d,t)->gmsh.model.mesh.getNodesForPhysicalGroup(d,t),
             entity->Tuple{Int,Int}[Tuple{Int,Int}(pair) for pair in
                 gmsh.model.getBoundary([entity],false,false,false)],
             (d,t)->gmsh.model.mesh.getNodes(d,t,false)[1],
             (d,t)->Tuple{Int,Int}[Tuple{Int,Int}(pair) for pair in
                 gmsh.model.mesh.getEmbedded(d,t)]),
            ("Tessella",
             (d,t)->Tessella.API.mesh.get_nodes_for_physical_group(d,t),
             entity->Tuple{Int,Int}[Tuple{Int,Int}(pair) for pair in
                 Tessella.API.model.get_boundary([entity],false,false,false)],
             (d,t)->Tessella.API.mesh.get_nodes(d,t)[1],
             (d,t)->Tuple{Int,Int}[Tuple{Int,Int}(pair) for pair in
                 Tessella.API.mesh.get_embedded(d,t)]))
            for (group,members) in (
                    ((2,10),[(2,1)]),((1,11),[(1,3),(1,1)]))
                tags,coords=query(group...)
                Set(UInt64.(tags))==group_expected(
                    get_boundary,get_nodes,get_embedded,members) || error(
                    "$engine group $group nodes are not the transitive " *
                    "boundary+embedded closure")
                issorted(tags) && allunique(tags) || error(
                    "$engine group $group nodes are not sorted and unique")
                length(coords)==3length(tags) || error(
                    "$engine group $group coordinate width mismatch")
            end
            isempty(query(2,99)[1]) || error(
                "$engine returned nodes for an unknown group")
            isempty(query(3,10)[1]) || error(
                "$engine returned nodes for a dimension-mismatched group")
        end
        # The two models share embeddings and point sizes exactly.
        gmsh.model.mesh.getEmbedded(2,1)==[(0,5)] || error(
            "Gmsh surface embedding record changed")
        Tessella.API.mesh.get_embedded(2,1)==Tuple{Int32,Int32}[(0,5)] ||
            error("Tessella surface embedding record differs")
        (gmsh.model.mesh.getEmbedded(1,1)==Tuple{Int32,Int32}[] &&
         Tessella.API.mesh.get_embedded(1,1)==Tuple{Int32,Int32}[]) ||
            error("empty embedding records differ")
        size_queries=[(0,1),(0,5),(1,1),(2,1),(9,9)]
        gmsh_sizes=gmsh.model.mesh.getSizes(size_queries)
        tessella_sizes=Tessella.API.mesh.get_sizes(size_queries)
        gmsh_sizes==tessella_sizes==[0.35,0.2,0.0,0.0,0.0] || error(
            "getSizes parity broke: gmsh=$gmsh_sizes tessella=$tessella_sizes")
        # Type-level filters resolve the tag in the type's own dimension.
        for (tessella_call,gmsh_call) in (
            (()->Tessella.API.mesh.get_elements_by_type(2,1),
             ()->gmsh.model.mesh.getElementsByType(2,1)),
            (()->Tessella.API.mesh.get_nodes_by_element_type(2,1),
             ()->gmsh.model.mesh.getNodesByElementType(2,1)),
            (()->Tessella.API.mesh.get_barycenters(2,1,false,false),
             ()->gmsh.model.mesh.getBarycenters(2,1,false,false)),
            (()->Tessella.API.mesh.get_element_edge_nodes(2,1,false),
             ()->gmsh.model.mesh.getElementEdgeNodes(2,1,false)),
            (()->Tessella.API.mesh.get_element_face_nodes(2,3,1,false),
             ()->gmsh.model.mesh.getElementFaceNodes(2,3,1,false)),
            (()->Tessella.API.mesh.get_jacobians(2,[0.25,0.25,0.0],1),
             ()->gmsh.model.mesh.getJacobians(2,[0.25,0.25,0.0],1)))
            tessella_value=tessella_call()
            gmsh_value=gmsh_call()
            length(tessella_value[1])>0 && length(gmsh_value[1])>0 || error(
                "an entity-filtered query returned no data")
        end
        tessella_tri_count=length(tessella_all[2][1])
        gmsh_tri_count=length(gmsh_all[2][1])
        length(Tessella.API.mesh.get_jacobians(
            2,[0.25,0.25,0.0],1)[2])==tessella_tri_count || error(
            "Tessella filtered Jacobian count differs")
        all(>(0),Tessella.API.mesh.get_jacobians(
            2,[0.25,0.25,0.0],1)[2]) || error(
            "Tessella filtered Jacobians are not positive")
        Tessella.API.mesh.get_basis_functions_orientation(
            2,"Lagrange",1)==zeros(Int32,tessella_tri_count) || error(
            "Tessella filtered orientations differ")
        _,tessella_key_entities,_=Tessella.API.mesh.get_keys(
            2,"Lagrange",1)
        tessella_surface_nodes=Set(
            Tessella.API.mesh.get_nodes(2,1,true)[1])
        for (edim,etag) in Tessella.API.mesh.get_embedded(2,1)
            union!(tessella_surface_nodes,
                   Tessella.API.mesh.get_nodes(edim,etag,true)[1])
        end
        sort!(unique!(tessella_key_entities))==
            sort!(collect(tessella_surface_nodes)) ||
            error("Tessella filtered Lagrange keys do not cover the entity")
        # Task partitions subdivide the entity-filtered subset.
        union_tags=UInt64[]
        for task in 0:2
            slice,_=Tessella.API.mesh.get_elements_by_type(2,1,task,3)
            append!(union_tags,slice)
        end
        sort!(union_tags)==sort!(tessella_filtered[2][1]) || error(
            "Tessella filtered task union differs")
        # Unknown entities fail explicitly in both engines.
        _rejects_argument(
            ()->Tessella.API.mesh.get_nodes(2,77)) || error(
            "Tessella accepted an unknown surface")
        gmsh_rejected=try
            gmsh.model.mesh.getNodes(2,77)
            false
        catch err
            err isa ErrorException || rethrow()
            true
        end
        gmsh_rejected || error("Gmsh accepted an unknown surface")
        _rejects_argument(
            ()->Tessella.API.mesh.get_elements_by_type(2,77)) || error(
            "Tessella accepted an unknown surface in a type query")
        gmsh_rejected=try
            gmsh.model.mesh.getElementsByType(2,77)
            false
        catch err
            err isa ErrorException || rethrow()
            true
        end
        gmsh_rejected || error(
            "Gmsh accepted an unknown surface in a type query")
        # A surface cache owns no dim-3 entity in either engine.
        _rejects_argument(
            ()->Tessella.API.mesh.get_elements(3,1)) || error(
            "Tessella classified a phantom volume")
        gmsh_rejected=try
            gmsh.model.mesh.getElements(3,1)
            false
        catch err
            err isa ErrorException || rethrow()
            true
        end
        gmsh_rejected || error("Gmsh classified a phantom volume")
        # Element-by-tag classification returns type, nodes, and owner entity.
        for element_tag in tessella_filtered[2][1]
            element_type,element_nodes,entity_dim,entity_tag=
                Tessella.API.mesh.get_element(element_tag)
            (element_type,entity_dim,entity_tag)==(Int32(2),2,1) ||
                error("Tessella element classification mismatch")
            length(element_nodes)==3 || error(
                "Tessella element node count mismatch")
        end
        for element_tag in gmsh_filtered[2][1]
            element_type,element_nodes,entity_dim,entity_tag=
                gmsh.model.mesh.getElement(element_tag)
            (element_type,entity_dim,entity_tag)==(Int32(2),2,1) ||
                error("Gmsh element classification mismatch")
        end
        _rejects_argument(
            ()->Tessella.API.mesh.get_element(0)) || error(
            "Tessella accepted an unknown element tag")
        gmsh_rejected=try
            gmsh.model.mesh.getElement(0)
            false
        catch err
            err isa ErrorException || rethrow()
            true
        end
        gmsh_rejected || error("Gmsh accepted an unknown element tag")

        # ── Entity-selective mutations ──────────────────────────────────────
        # affineTransform(dimTags) moves only nodes classified on the listed
        # entities in both engines; unlisted entities keep their coordinates.
        gmsh_curve=copy(gmsh.model.mesh.getNodes(1,1,false)[2])
        tessella_curve=copy(Tessella.API.mesh.get_nodes(1,1)[2])
        gmsh_surface_coords=copy(gmsh.model.mesh.getNodes(2,1,false)[2])
        tessella_surface_coords=copy(
            Tessella.API.mesh.get_nodes(2,1)[2])
        lift=[1.0,0,0,0, 0,1,0,0, 0,0,1,0.5, 0,0,0,1]
        gmsh.model.mesh.affineTransform(lift,[(1,1)])
        Tessella.API.mesh.affine_transform(lift,[(1,1)])
        gmsh.model.mesh.getNodes(1,1,false)[2][3:3:end]==
            gmsh_curve[3:3:end].+0.5 || error(
            "Gmsh selective transform did not move the listed curve")
        Tessella.API.mesh.get_nodes(1,1)[2][3:3:end]==
            tessella_curve[3:3:end].+0.5 || error(
            "Tessella selective transform did not move the listed curve")
        gmsh.model.mesh.getNodes(2,1,false)[2]==gmsh_surface_coords ||
            error("Gmsh selective transform moved an unlisted surface")
        Tessella.API.mesh.get_nodes(2,1)[2]==tessella_surface_coords ||
            error("Tessella selective transform moved an unlisted surface")
        # Unknown entities fail explicitly in both engines.
        _rejects_argument(
            ()->Tessella.API.mesh.affine_transform(lift,[(1,77)])) || error(
            "Tessella transformed an unknown curve")
        gmsh_rejected=try
            gmsh.model.mesh.affineTransform(lift,[(1,77)])
            false
        catch err
            err isa ErrorException || rethrow()
            true
        end
        gmsh_rejected || error("Gmsh transformed an unknown curve")
        # Selective edge creation tags only the listed entities' cell edges;
        # a later full creation preserves them.
        gmsh.model.mesh.createEdges([(2,1)])
        Tessella.API.mesh.create_edges([(2,1)])
        gmsh_selective=copy(gmsh.model.mesh.getAllEdges()[1])
        tessella_selective=copy(Tessella.API.mesh.get_all_edges()[1])
        isempty(gmsh_selective) && error(
            "Gmsh selective edge creation produced no edges")
        isempty(tessella_selective) && error(
            "Tessella selective edge creation produced no edges")
        Tessella.API.mesh.create_edges([(2,1)])
        Tessella.API.mesh.get_all_edges()[1]==tessella_selective || error(
            "Tessella selective edge creation was not idempotent")
        gmsh.model.mesh.createEdges()
        Tessella.API.mesh.create_edges()
        gmsh_all_edges=gmsh.model.mesh.getAllEdges()[1]
        tessella_all_edges=Tessella.API.mesh.get_all_edges()[1]
        issubset(tessella_selective,tessella_all_edges) || error(
            "Tessella selective edge tags were not preserved")
        issubset(gmsh_selective,gmsh_all_edges) || error(
            "Gmsh selective edge tags were not preserved")
        # removeElements drops only the listed cell on the entity and keeps
        # every node in both engines; dense tags re-index in Tessella, so
        # compare counts and connectivity sets rather than tags.
        # `getElementsByType` reports type-local tags in Gmsh that
        # `removeElements` rejects; entity removal needs `getElements` tags.
        gmsh_tri_tags=gmsh.model.mesh.getElements(2,1)[2][1]
        tessella_tri_tags=
            Tessella.API.mesh.get_elements_by_type(2,1)[1]
        gmsh_node_count=length(gmsh.model.mesh.getNodes()[1])
        tessella_node_count=length(Tessella.API.mesh.get_nodes()[1])
        gmsh.model.mesh.removeElements(2,1,[gmsh_tri_tags[1]])
        Tessella.API.mesh.remove_elements(2,1,[tessella_tri_tags[1]])
        length(gmsh.model.mesh.getElements(2,1)[2][1])==
            length(gmsh_tri_tags)-1 || error(
            "Gmsh removeElements kept the listed triangle")
        length(Tessella.API.mesh.get_elements_by_type(2,1)[1])==
            length(tessella_tri_tags)-1 || error(
            "Tessella remove_elements kept the listed triangle")
        length(gmsh.model.mesh.getNodes()[1])==gmsh_node_count || error(
            "Gmsh removeElements dropped nodes")
        length(Tessella.API.mesh.get_nodes()[1])==tessella_node_count ||
            error("Tessella remove_elements dropped nodes")
        # A tag owned by another block is rejected by both engines.
        _rejects_argument(()->Tessella.API.mesh.remove_elements(
            1,1,Tessella.API.mesh.get_elements_by_type(2,1)[1][1:1])) ||
            error("Tessella removed a triangle through a curve selection")
        gmsh_rejected=try
            gmsh.model.mesh.removeElements(
                1,1,gmsh.model.mesh.getElements(2,1)[2][1][1:1])
            false
        catch err
            err isa ErrorException || rethrow()
            true
        end
        gmsh_rejected || error(
            "Gmsh removed a triangle through a curve selection")
        _rejects_argument(()->Tessella.API.mesh.remove_elements(2,77)) ||
            error("Tessella removed elements on an unknown surface")
        gmsh_rejected=try
            gmsh.model.mesh.removeElements(2,77)
            false
        catch err
            err isa ErrorException || rethrow()
            true
        end
        gmsh_rejected || error("Gmsh removed elements on an unknown surface")
        # Entity-selective reversal uses Gmsh's first-order convention — a
        # triangle (a,b,c) becomes (a,c,b) — in both engines.
        gmsh_before=copy(gmsh.model.mesh.getElementsByType(2,1)[2][1:3])
        tessella_before=copy(
            Tessella.API.mesh.get_elements_by_type(2,1)[2][1:3])
        gmsh.model.mesh.reverse([(2,1)])
        Tessella.API.mesh.reverse([(2,1)])
        gmsh.model.mesh.getElementsByType(2,1)[2][1:3]==
            [gmsh_before[1],gmsh_before[3],gmsh_before[2]] || error(
            "Gmsh reverse did not swap the last two triangle nodes")
        Tessella.API.mesh.get_elements_by_type(2,1)[2][1:3]==
            [tessella_before[1],tessella_before[3],tessella_before[2]] ||
            error("Tessella reverse did not swap the last two triangle nodes")
        gmsh.model.mesh.reverse([(2,1)])
        Tessella.API.mesh.reverse([(2,1)])
        Tessella.API.mesh.get_elements_by_type(2,1)[2][1:3]==
            tessella_before || error("Tessella double reverse differed")
        gmsh.model.mesh.getElementsByType(2,1)[2][1:3]==gmsh_before ||
            error("Gmsh double reverse differed")
        _rejects_argument(()->Tessella.API.mesh.reverse([(2,77)])) || error(
            "Tessella reversed an unknown surface")
        gmsh_rejected=try
            gmsh.model.mesh.reverse([(2,77)])
            false
        catch err
            err isa ErrorException || rethrow()
            true
        end
        gmsh_rejected || error("Gmsh reversed an unknown surface")
        # Tag-level reversal matches in both engines.
        gmsh.model.mesh.reverseElements(
            [gmsh.model.mesh.getElements(2,1)[2][1][1]])
        Tessella.API.mesh.reverse_elements(
            [Tessella.API.mesh.get_elements_by_type(2,1)[1][1]])
        Tessella.API.mesh.get_elements_by_type(2,1)[2][1:3]==
            [tessella_before[1],tessella_before[3],tessella_before[2]] ||
            error("Tessella reverse_elements differed")
        gmsh.model.mesh.getElementsByType(2,1)[2][1:3]==
            [gmsh_before[1],gmsh_before[3],gmsh_before[2]] || error(
            "Gmsh reverseElements differed")
        Tessella.API.mesh.reverse_elements(
            [Tessella.API.mesh.get_elements_by_type(2,1)[1][1]])
        gmsh.model.mesh.reverseElements(
            [gmsh.model.mesh.getElements(2,1)[2][1][1]])
        # reorderElements permutes one entity's type block by 0-based source
        # positions; reversing the surface's triangles must reverse the
        # per-element vertex-coordinate sequence in both engines.
        function element_coord_sequence(engine)
            if engine==:gmsh
                _,_,conn=gmsh.model.mesh.getElements(2,1)
                flat=only(conn)
                coords=Dict{UInt64,NTuple{3,Float64}}(
                    tag=>Tuple(gmsh.model.mesh.getNode(tag)[1]) for tag in
                    unique(flat))
                return [sort!([coords[node] for node in
                        flat[offset:offset+2]])
                        for offset in 1:3:length(flat)]
            else
                flat=Tessella.API.mesh.get_elements_by_type(2,1)[2]
                all_tags,all_coords,_=Tessella.API.mesh.get_nodes()
                nodes=reshape(all_coords,3,:)
                coordinate=Dict{UInt64,NTuple{3,Float64}}(
                    all_tags[i]=>Tuple(nodes[:,i]) for i in
                    eachindex(all_tags))
                return [sort!([coordinate[node] for node in
                        flat[offset:offset+2]])
                        for offset in 1:3:length(flat)]
            end
        end
        gmsh_sequence=element_coord_sequence(:gmsh)
        tessella_sequence=element_coord_sequence(:tessella)
        gmsh.model.mesh.reorderElements(2,1,
            collect(length(gmsh_sequence)-1:-1:0))
        Tessella.API.mesh.reorder_elements(2,1,
            collect(length(tessella_sequence)-1:-1:0))
        element_coord_sequence(:gmsh)==reverse(gmsh_sequence) || error(
            "Gmsh reorderElements did not reverse the entity block")
        element_coord_sequence(:tessella)==reverse(tessella_sequence) ||
            error("Tessella reorder_elements did not reverse the entity block")
        gmsh.model.mesh.reorderElements(2,1,
            collect(0:length(gmsh_sequence)-1))
        Tessella.API.mesh.reorder_elements(2,1,
            collect(0:length(tessella_sequence)-1))
        gmsh_rejected=false
        try
            gmsh.model.mesh.reorderElements(2,1,[0,1])
        catch err
            gmsh_rejected=true
        end
        gmsh_rejected || error("Gmsh accepted a short ordering")
        _rejects_argument(
            ()->Tessella.API.mesh.reorder_elements(2,1,[0,1])) || error(
            "Tessella accepted a short ordering")
        _rejects_argument(
            ()->Tessella.API.mesh.reorder_elements(1,1,[0])) || error(
            "Tessella reordered an empty type block")
        # getDuplicateNodes reports no coincident nodes on either conforming
        # cache, under whole-mesh and entity-filtered scans alike.
        isempty(gmsh.model.mesh.getDuplicateNodes()) || error(
            "Gmsh reported duplicate nodes")
        isempty(Tessella.API.mesh.get_duplicate_nodes()) || error(
            "Tessella reported duplicate nodes")
        for entity in ((2,1),(1,1),(0,1))
            isempty(gmsh.model.mesh.getDuplicateNodes([entity])) || error(
                "Gmsh reported duplicates on $entity")
            isempty(Tessella.API.mesh.get_duplicate_nodes([entity])) ||
                error("Tessella reported duplicates on $entity")
        end
        # Clearing a curve under a meshed surface is a no-op in both engines:
        # the boundary mesh stays part of the surviving surface mesh.
        gmsh_counts=Dict(ent=>length(gmsh.model.mesh.getNodes(
            ent[1],ent[2],false)[1]) for ent in
            ((0,1),(1,1),(1,2),(2,1)))
        tessella_counts=Dict(ent=>length(Tessella.API.mesh.get_nodes(
            ent[1],ent[2])[1]) for ent in ((0,1),(1,1),(1,2),(2,1)))
        gmsh.model.mesh.clear([(1,1),(0,1)])
        Tessella.API.mesh.clear([(1,1),(0,1)])
        for ent in keys(gmsh_counts)
            length(gmsh.model.mesh.getNodes(ent[1],ent[2],false)[1])==
                gmsh_counts[ent] || error(
                "Gmsh boundary clear changed entity $ent")
            length(Tessella.API.mesh.get_nodes(ent[1],ent[2])[1])==
                tessella_counts[ent] || error(
                "Tessella boundary clear changed entity $ent")
        end
        # Clearing the generating surface drops its elements and owned nodes;
        # boundary-owned nodes survive on their entities in both engines.
        gmsh.model.mesh.clear([(2,1)])
        Tessella.API.mesh.clear([(2,1)])
        all(isempty,gmsh.model.mesh.getElements(2,1)[2]) || error(
            "Gmsh retained surface elements after clear")
        all(isempty,Tessella.API.mesh.get_elements(2,1)[2]) || error(
            "Tessella retained surface elements after clear")
        isempty(gmsh.model.mesh.getNodes(2,1,false)[1]) || error(
            "Gmsh retained surface nodes after clear")
        isempty(Tessella.API.mesh.get_nodes(2,1)[1]) || error(
            "Tessella retained surface nodes after clear")
        length(gmsh.model.mesh.getNodes(1,2,false)[1])==
            gmsh_counts[(1,2)] || error(
            "Gmsh dropped surviving curve nodes")
        length(Tessella.API.mesh.get_nodes(1,2)[1])==
            tessella_counts[(1,2)] || error(
            "Tessella dropped surviving curve nodes")
        _rejects_argument(
            ()->Tessella.API.mesh.clear([(2,99)])) || error(
            "Tessella cleared an unknown surface")
        (tessella_tri_count,gmsh_tri_count)
    finally
        Tessella.API.finalize()
    end

    # Positive duplicate-node and duplicate-element parity on handcrafted
    # coincident caches: Gmsh stores them on discrete entities, Tessella on a
    # directly installed cache. Gmsh merges globally by coordinates, keeps the
    # lowest tag, compacts numbering, and remaps connectivity.
    gmsh.model.add("duplicates")
    gmsh.model.addDiscreteEntity(2,1)
    gmsh.model.mesh.addNodes(2,1,[1,2,3,4],
        [0.0,0,0, 1.0,0,0, 0.0,1,0, 0.0,0,0])
    gmsh.model.mesh.addElementsByType(1,2,[10,11],[1,2,3, 4,2,3])
    gmsh.model.addDiscreteEntity(0,5)
    gmsh.model.mesh.addNodes(0,5,[7],[1.0,0,0])
    gmsh.model.mesh.addElementsByType(5,15,[30],[7])
    Tessella.API.initialize()
    try
        # A real model entity set lets the classified-only parity calls below
        # resolve entity (2,1) and gives remove_embedded a genuine record.
        for (point,(x,y)) in enumerate(
                ((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)))
            Tessella.API.model.add_point(x,y,0.0;tag=point,meshSize=0.5)
        end
        for (curve,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            Tessella.API.model.add_line(a,b;tag=curve)
        end
        Tessella.API.model.add_curve_loop([1,2,3,4];tag=1)
        Tessella.API.model.add_plane_surface([1];tag=1)
        Tessella.API.model.add_point(0.5,0.5,0.0;tag=5)
        Tessella.API.model.embed(0,[5],2,1)
        dup_coordinates=Float64[0 1 0 0 1;
                                0 0 1 0 0;
                                0 0 0 0 0]
        lock(Tessella.API.STATE_LOCK) do
            Tessella.API.LAST_MESH[]=Mesh(dup_coordinates;
                tris=Int32[1 4;
                           2 2;
                           3 3])
            Tessella.API.LAST_MESH_CLASS[]=nothing
        end
        gmsh_duplicates=sort!(map(tag->Tuple(
                gmsh.model.mesh.getNode(tag)[1]),
            gmsh.model.mesh.getDuplicateNodes()))
        tessella_tags,tessella_coords,_=Tessella.API.mesh.get_nodes()
        tessella_node_coords=reshape(tessella_coords,3,:)
        tessella_duplicates=sort!([Tuple(tessella_node_coords[:,tag])
            for tag in Tessella.API.mesh.get_duplicate_nodes()])
        tessella_duplicates==gmsh_duplicates || error(
            "duplicate-node reporting differed: $tessella_duplicates vs " *
            "$gmsh_duplicates")
        gmsh.model.mesh.removeDuplicateNodes()
        Tessella.API.mesh.remove_duplicate_nodes()
        gmsh_tags,gmsh_coords,_=gmsh.model.mesh.getNodes()
        sort!(collect(eachcol(reshape(gmsh_coords,3,:))),
            by=Tuple)==sort!(collect(eachcol(reshape(
                Tessella.API.mesh.get_nodes()[2],3,:))),by=Tuple) || error(
            "duplicate-node merge differed")
        # Gmsh renumbers survivors globally in entity order; compare each
        # element's vertex-coordinate set rather than raw tags.
        _,_,gmsh_conn=gmsh.model.mesh.getElements(2,1)
        gmsh_flat=only(gmsh_conn)
        gmsh_element_coords=sort!([sort!([Tuple(
            gmsh.model.mesh.getNode(node)[1])
            for node in gmsh_flat[offset:offset+2]])
            for offset in 1:3:length(gmsh_flat)])
        tessella_conn=reshape(
            Tessella.API.mesh.get_elements_by_type(2)[2],3,:)
        merged_coords=reshape(Tessella.API.mesh.get_nodes()[2],3,:)
        tessella_element_coords=sort!([sort!([Tuple(merged_coords[:,node])
            for node in column]) for column in eachcol(tessella_conn)])
        tessella_element_coords==gmsh_element_coords || error(
            "duplicate-node connectivity remap differed")
        # Both triangles now share connectivity [1,2,3]; dedup keeps one.
        gmsh.model.mesh.removeDuplicateElements()
        _,gmsh_element_tags,_=gmsh.model.mesh.getElements(2,1)
        only(gmsh_element_tags) |> length==1 || error(
            "Gmsh kept duplicate elements")
        lock(Tessella.API.STATE_LOCK) do
            cache=Tessella.API.LAST_MESH[]
            Tessella.API.LAST_MESH_CLASS[]=
                Tessella.API._MeshClassification(cache,(2,Int32(1)),
                    [(2,Int32(1))],fill((2,Int32(1)),3),
                    Dict{Tuple{Int,Int32},Vector{Int32}}(),
                    Int32[],Int32[1,1],Int32[])
        end
        Tessella.API.mesh.remove_duplicate_elements()
        length(Tessella.API.mesh.get_elements_by_type(2)[1])==1 || error(
            "Tessella kept duplicate elements")
        _rejects_argument(
            ()->Tessella.API.mesh.remove_duplicate_elements([(2,99)])) ||
            error("Tessella deduplicated an unknown surface")
        # A tag swap renumbering preserves the coordinate set in both engines;
        # Gmsh's surviving tags follow global entity order while Tessella keeps
        # lowest-tag retention, so only the multiset is comparable.
        gmsh.model.mesh.renumberNodes([1,2],[2,1])
        Tessella.API.mesh.renumber_nodes([1,2],[2,1])
        _,gmsh_renum_coords,_=gmsh.model.mesh.getNodes()
        sort!(collect(eachcol(reshape(gmsh_renum_coords,3,:))),
            by=Tuple)==sort!(collect(eachcol(reshape(
                Tessella.API.mesh.get_nodes()[2],3,:))),by=Tuple) || error(
            "renumber_nodes swap differed")
        # setNode parity: both engines move node tag 1.
        gmsh.model.mesh.setNode(1,[9.0,8.0,7.0],[])
        Tessella.API.mesh.set_node(1,[9.0,8.0,7.0])
        Tuple(gmsh.model.mesh.getNode(1)[1])==Tuple(reshape(
            Tessella.API.mesh.get_nodes()[2],3,:)[:,1]) || error(
            "set_node coordinate update differed")
        _rejects_argument(()->Tessella.API.mesh.set_node(99,[0,0,0.0])) ||
            error("Tessella moved an unknown node")
        gmsh_rejected=false
        try
            gmsh.model.mesh.setNode(99,[0,0,0.0],[])
        catch err
            gmsh_rejected=true
        end
        gmsh_rejected || error("Gmsh moved an unknown node")
        gmsh.model.mesh.renumberNodes()
        Tessella.API.mesh.renumber_nodes()
        gmsh.model.mesh.renumberElements()
        Tessella.API.mesh.renumber_elements()
        _rejects_argument(
            ()->Tessella.API.mesh.renumber_nodes([1],[99])) || error(
            "Tessella accepted a sparse node renumbering")
        # Unpartitioned-cache parity calls return empty or no-op in both.
        isempty(gmsh.model.mesh.getGhostElements(2,1)[1]) || error(
            "Gmsh reported ghost elements")
        Tessella.API.mesh.get_ghost_elements(2,1)==
            (UInt64[],Int32[]) || error(
            "Tessella reported ghost elements")
        isempty(gmsh.model.mesh.getLastEntityError()) || error(
            "Gmsh reported entity errors")
        isempty(Tessella.API.mesh.get_last_entity_error()) || error(
            "Tessella reported entity errors")
        isempty(gmsh.model.mesh.getLastNodeError()) || error(
            "Gmsh reported node errors")
        isempty(Tessella.API.mesh.get_last_node_error()) || error(
            "Tessella reported node errors")
        gmsh.model.mesh.unpartition()
        Tessella.API.mesh.unpartition()
        gmsh.model.mesh.rebuildNodeCache()
        Tessella.API.mesh.rebuild_node_cache()
        gmsh.model.mesh.rebuildElementCache()
        Tessella.API.mesh.rebuild_element_cache()
        gmsh.model.mesh.reclassifyNodes()
        Tessella.API.mesh.reclassify_nodes()
        gmsh.model.mesh.relocateNodes(2,1)
        Tessella.API.mesh.relocate_nodes(2,1)
        # removeConstraints clears per-entity meshing attributes only; neither
        # engine stores cleared attributes here, so both are validated no-ops.
        gmsh.model.mesh.removeConstraints()
        Tessella.API.mesh.remove_constraints()
        gmsh.model.mesh.removeConstraints([(2,1),(0,5)])
        Tessella.API.mesh.remove_constraints([(2,1),(0,5)])
        gmsh_rejected=false
        try
            gmsh.model.mesh.removeConstraints([(1,99)])
        catch err
            gmsh_rejected=true
        end
        gmsh_rejected || error("Gmsh accepted an unknown constraint entity")
        _rejects_argument(
            ()->Tessella.API.mesh.remove_constraints([(1,99)])) || error(
            "Tessella accepted an unknown constraint entity")
        # computeRenumbering returns sorted old node tags and a permutation of
        # 1:n in both engines; the RCM ordering itself is
        # implementation-defined. Restricting to the surviving triangle makes
        # the involved-node counts comparable — Gmsh's point element has no
        # Tessella simplex counterpart.
        for (engine,fn) in (
                ("gmsh",()->gmsh.model.mesh.computeRenumbering()),
                ("tessella",()->Tessella.API.mesh.compute_renumbering()))
            old_tags,new_tags=fn()
            issorted(old_tags) || error("$engine renumbering old tags unsorted")
            sort(Int.(new_tags))==collect(1:length(new_tags)) || error(
                "$engine renumbering new tags are not a dense permutation")
        end
        _,gmsh_triangle_tags,_=gmsh.model.mesh.getElements(2,1)
        gmsh_triangle_tag=Int.(only(gmsh_triangle_tags))[1]
        tessella_triangle_tag=Int(
            Tessella.API.mesh.get_elements_by_type(2)[1][1])
        gmsh_old,_=gmsh.model.mesh.computeRenumbering(
            "RCMK",[gmsh_triangle_tag])
        tessella_old,_=Tessella.API.mesh.compute_renumbering(
            "RCMK",[tessella_triangle_tag])
        length(gmsh_old)==length(tessella_old)==3 || error(
            "restricted renumbering node counts differ: " *
            "Gmsh=$(length(gmsh_old)) Tessella=$(length(tessella_old))")
        gmsh_rejected=false
        try
            gmsh.model.mesh.computeRenumbering("bogus")
        catch err
            gmsh_rejected=true
        end
        gmsh_rejected || error("Gmsh accepted an unknown renumbering method")
        _rejects_argument(
            ()->Tessella.API.mesh.compute_renumbering("bogus")) || error(
            "Tessella accepted an unknown renumbering method")
        # optimize() preserves connectivity and node counts in both engines.
        gmsh_counts_before=length(gmsh.model.mesh.getNodes()[1])
        tessella_counts_before=length(Tessella.API.mesh.get_nodes()[1])
        gmsh.model.mesh.optimize()
        Tessella.API.mesh.optimize()
        length(gmsh.model.mesh.getNodes()[1])==gmsh_counts_before || error(
            "Gmsh optimize changed the node count")
        length(Tessella.API.mesh.get_nodes()[1])==tessella_counts_before ||
            error("Tessella optimize changed the node count")
        _rejects_argument(
            ()->Tessella.API.mesh.optimize("NoSuchOptimizer")) || error(
            "Tessella accepted an unknown optimizer")
        # Per-element visibility is raw display state: default 1, stored values
        # pass through, and unknown tags report 0 in both engines.
        gmsh_triangle_tags,_=gmsh.model.mesh.getElementsByType(2)
        tessella_triangle_tags,_=Tessella.API.mesh.get_elements_by_type(2)
        gmsh.model.mesh.setVisibility(gmsh_triangle_tags[1:1],0)
        Tessella.API.mesh.set_visibility(tessella_triangle_tags[1:1],0)
        Int.(gmsh.model.mesh.getVisibility(gmsh_triangle_tags[1:1]))==
            Int.(Tessella.API.mesh.get_visibility(
                tessella_triangle_tags[1:1]))==[0] || error(
            "element visibility differed")
        Int.(gmsh.model.mesh.getVisibility([987654]))==
            Int.(Tessella.API.mesh.get_visibility([987654]))==[0] || error(
            "unknown-element visibility differed")

        # removeEmbedded lists the parent entities in both engines.
        Tessella.API.mesh.get_embedded(2,1)==Tuple{Int32,Int32}[(0,5)] ||
            error("Tessella lost its embedded-point record")
        gmsh.model.mesh.embed(0,[5],2,1)
        gmsh.model.mesh.removeEmbedded([(2,1)])
        isempty(gmsh.model.mesh.getEmbedded(2,1)) || error(
            "Gmsh kept the embedding")
        Tessella.API.mesh.remove_embedded([(2,1)])
        isempty(Tessella.API.mesh.get_embedded(2,1)) || error(
            "Tessella kept the embedding")
    finally
        Tessella.API.finalize()
    end

    println("mesh-data-query differential: Gmsh ",gmsh.GMSH_API_VERSION,
            ", types=1,2,4 nodes=4 elements=3 connectivity_entries=9 ",
            "dense_max_tags=4/3 explicit_max_tags=40/300 ",
            "derived_sha=",derived_sha," ",
            "refined_sha=",refined_crc.sha,
            " entity_filtered=tris",entity_pairs,
            " selective=transform/edges/remove/reverse/clear" *
            " duplicates=detect/merge/drop" *
            " set_node/renumber/embedded/constraints/rcmk/optimize/visibility" *
            " parametric=entity/owner/node",
            " bounded=no-mesh/classification/special-type/face-count ",
            "blockers and finite-barycenter contract with partitioned slices")
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
