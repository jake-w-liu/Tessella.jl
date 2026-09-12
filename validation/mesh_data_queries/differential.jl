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
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    Tessella.API.initialize()
    entity_pairs=try
        for (point,(x,y)) in enumerate(
                ((0.0,0.0),(1.0,0.0),(1.0,1.0),(0.0,1.0)))
            Tessella.API.model.add_point(x,y,0.0;tag=point)
        end
        for (curve,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            Tessella.API.model.add_line(a,b;tag=curve)
        end
        Tessella.API.model.add_curve_loop([1,2,3,4];tag=1)
        Tessella.API.model.add_plane_surface([1];tag=1)
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
                 for point in 1:4,
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
                 for point in 1:4,
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
        sort!(unique!(tessella_key_entities))==
            Tessella.API.mesh.get_nodes(2,1,true)[1] |> sort |> unique ||
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

    println("mesh-data-query differential: Gmsh ",gmsh.GMSH_API_VERSION,
            ", types=1,2,4 nodes=4 elements=3 connectivity_entries=9 ",
            "dense_max_tags=4/3 explicit_max_tags=40/300 ",
            "derived_sha=",derived_sha," ",
            "refined_sha=",refined_crc.sha,
            " entity_filtered=tris",entity_pairs,
            " selective=transform/edges/clear parametric=entity/owner/node",
            " bounded=no-mesh/classification/special-type/face-count ",
            "blockers and finite-barycenter contract with partitioned slices")
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
