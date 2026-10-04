#!/usr/bin/env julia
# Gmsh supplies an independent, allocation-dependent admissible mesh. Production
# never invokes it, and this gate never pins a random pointer-selected topology.
using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using TOML, SHA
using Tessella.Elements: MixedMesh, ElementBlock, msh_spec
using Tessella.MeshTypes: nnodes

include(joinpath(@__DIR__,"..","..","test","geometry","quadtri_nonew_certificates.jl"))
const QTNN=QuadTriNoNewCertificates
include(joinpath(@__DIR__,"..","..","test","geometry","quadtri_nonew_triangle_certificates.jl"))
const QTNT=QuadTriNoNewTriangleCertificates
include(joinpath(@__DIR__,"..","..","test","geometry","quadtri_nonew_two_tri_certificates.jl"))
const QTT=QuadTriNoNewTwoTriCertificates
include(joinpath(@__DIR__,"..","..","test","geometry","quadtri_nonew_quad_patch_certificates.jl"))
const QTP=QuadTriNoNewQuadPatchCertificates
include(joinpath(@__DIR__,"..","..","test","geometry","quadtri_nonew_quad_strip_certificates.jl"))
const QTS=QuadTriNoNewQuadStripCertificates
include(joinpath(@__DIR__,"..","..","test","geometry","quadtri_nonew_three_quad_strip_certificates.jl"))
const QT3S=QuadTriNoNewThreeQuadStripCertificates
include(joinpath(@__DIR__,"..","..","test","geometry","quadtri_nonew_four_quad_strip_certificates.jl"))
const QT4S=QuadTriNoNewFourQuadStripCertificates
include(joinpath(@__DIR__,"..","..","test","geometry","quadtri_nonew_rect_grid_certificates.jl"))
const QTRG=QuadTriNoNewRectGridCertificates
binding=get(ENV,"GMSH_JULIA_API","")
isfile(binding) || error("set GMSH_JULIA_API to pinned Gmsh4.15.2 gmsh.jl")
include(binding)
gmsh.GMSH_API_VERSION=="4.15.2" || error("Gmsh4.15.2 binding required")

function oracle_volume(entity=1)
    node_tags,xyz,_=gmsh.model.mesh.getNodes()
    coordinates=reshape(xyz,3,:)
    source=Dict(tag=>index for (index,tag) in enumerate(node_tags))
    types,tags,nodes=gmsh.model.mesh.getElements(3,entity)
    used=sort!(unique(vcat(nodes...)))
    remap=Dict(tag=>Int32(index) for (index,tag) in enumerate(used))
    blocks=ElementBlock[]
    for (type,cell_nodes) in zip(types,nodes)
        width=gmsh.model.mesh.getElementProperties(type)[4]
        connectivity=reshape(Int32[remap[tag] for tag in cell_nodes],width,:)
        push!(blocks,ElementBlock(type,connectivity))
    end
    return MixedMesh(coordinates[:,[source[tag] for tag in used]],blocks),vcat(tags...)
end

function native_api(path,action)
    api=Tessella.API
    api.initialize()
    try
        api.open_geo!(path;mesh_dim=0)
        return action(api)
    finally
        api.finalize()
    end
end

# Read the actual public cell product, including the legacy simplex P2 overlay.
# mesh.get() deliberately continues to expose the primary Mesh in that route.
function triangle_api_volume(api,entity=1)
    tags,xyz,_=api.mesh.get_nodes();coordinates=reshape(xyz,3,:)
    positions=Dict(tag=>index for (index,tag) in enumerate(tags))
    types,_,families=api.mesh.get_elements(3,entity)
    used=sort!(unique(vcat(families...)))
    remap=Dict(tag=>Int32(index) for (index,tag) in enumerate(used))
    blocks=ElementBlock[]
    for (type,nodes) in zip(types,families)
        width=msh_spec(type).nnodes
        push!(blocks,ElementBlock(type,reshape(Int32[remap[tag] for tag in nodes],width,:)))
    end
    return MixedMesh(coordinates[:,[positions[tag] for tag in used]],blocks)
end

function quadratic_volume_from_jacobians(mesh,interface;expected=1.)
    total=0.
    for block in mesh.blocks
        if interface===gmsh.model.mesh
            locations,weights=interface.getIntegrationPoints(block.msh,"Gauss4")
            _,determinants,_=interface.getJacobians(block.msh,locations)
        else
            locations,weights=interface.get_integration_points(block.msh,"Gauss4")
            _,determinants,_=interface.get_jacobians(block.msh,locations)
        end
        all(d->isfinite(d)&&d>0,determinants) || error("nonpositive P2 quadrature Jacobian")
        length(determinants)==length(weights)*size(block.nodes,2) || error("P2 quadrature shape")
        total+=sum(determinants .* repeat(weights,size(block.nodes,2)))
    end
    abs(total-expected)<=2e-11 || error("P2 integrated swept volume $total differs from $expected")
    return total
end

function oracle_surface_faces(volume)
    lookup=Dict(QTNN.point(volume,index)=>index for index in 1:nnodes(volume))
    tags,xyz,_=gmsh.model.mesh.getNodes()
    coordinates=reshape(xyz,3,:)
    points=Dict(tag=>Tuple(coordinates[:,index]) for (index,tag) in enumerate(tags))
    result=Set{Tuple}()
    types,_,families=gmsh.model.mesh.getElements(2)
    for (type,nodes) in zip(types,families)
        width=gmsh.model.mesh.getElementProperties(type)[4]
        for cell in eachcol(reshape(nodes,width,:))
            push!(result,QTNN.key(lookup[points[node]] for node in cell))
        end
    end
    return result
end

function same_grid_and_centroid(native,oracle,fixture)
    nc=QTNN.certify_fixture(native,fixture)
    oc=QTNN.certify_fixture(oracle,fixture)
    nnodes(native)==nnodes(oracle) || error("$(fixture.name): node counts differ")
    for i in eachindex(nc.grid_ids)
        p=QTNN.point(native,nc.grid_ids[i]);q=QTNN.point(oracle,oc.grid_ids[i])
        maximum(abs.(p.-q))<=2e-11 || error("$(fixture.name): swept coordinates differ")
    end
    if fixture.laterals
        p=QTNN.point(native,only(nc.extras));q=QTNN.point(oracle,only(oc.extras))
        maximum(abs.(p.-q))<=2e-11 || error("$(fixture.name): final centroid differs")
        nc.family_counts==oc.family_counts || error("$(fixture.name): recombined families differ")
    end
    # Free lateral diagonals and even family counts can legitimately differ.
    QTNN.surface_faces(native_execution[].mesh_parts,native)==Set(keys(nc.boundary)) ||
        error("$(fixture.name): native volume/surface complex differs")
    oracle_surface_faces(oracle)==Set(keys(oc.boundary)) ||
        error("$(fixture.name): oracle volume/surface complex differs")
    return nc,oc
end

const native_execution=Ref{Any}(nothing)
function initialize_oracle()
    gmsh.initialize(["-nopopup"],false,false)
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.option.setNumber("Mesh.ElementOrder",1)
    gmsh.option.setNumber("Mesh.SecondOrderLinear",0)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete",0)
end

# These twelve products were already captured from pinned Gmsh 4.15.2. Replay
# their complete integer topology and provenance without regenerating them.
# Raw oracle tags/templates remain evidence, not native allocation pins.
const TWO_TRI_EDGE_PATTERNS=Dict(
    2=>((1,2),(2,3),(3,1)),3=>((1,2),(2,3),(3,4),(4,1)),
    4=>((1,2),(2,3),(3,1),(1,4),(2,4),(3,4)),
    6=>((1,2),(2,3),(3,1),(4,5),(5,6),(6,4),(1,4),(2,5),(3,6)))
const TWO_TRI_FACE_PATTERNS=Dict(
    4=>((1,3,2),(1,2,4),(2,3,4),(3,1,4)),
    6=>((1,3,2),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6)))

function two_tri_saved_mesh(record,phase,dimension,entity)
    stage=record[phase]
    rows=only(filter(e->e["dim"]==dimension && e["tag"]==entity,stage["entities"]))["cells"]
    used=sort!(unique(vcat((row["nodes"] for row in rows)...)))
    index=Dict(tag=>Int32(i) for (i,tag) in enumerate(used))
    nodes=Dict(row["tag"]=>row for row in stage["nodes"])
    coordinates=hcat((Float64.(nodes[tag]["coordinates"]) for tag in used)...)
    types=unique(row["type"] for row in rows)
    blocks=ElementBlock[]
    for msh in types
        chosen=filter(row->row["type"]==msh,rows)
        cells=hcat((Int32[index[tag] for tag in row["nodes"]] for row in chosen)...)
        push!(blocks,ElementBlock(msh,cells))
    end
    mesh=dimension==2 && types==[2] ? Tessella.MeshTypes.Mesh(coordinates;tris=only(blocks).nodes) :
                                      MixedMesh(coordinates,blocks)
    return mesh,index
end

function two_tri_add_carrier!(catalog,support,owner)
    push!(get!(catalog,Tuple(sort!(collect(support))),Set{Tuple{Int,Int}}()),owner)
end

function two_tri_actual_carriers(rows)
    catalog=Dict{Tuple,Set{Tuple{Int,Int}}}()
    for (msh,cell,owner) in rows
        if msh==1
            two_tri_add_carrier!(catalog,cell,owner)
        elseif haskey(TWO_TRI_EDGE_PATTERNS,msh)
            for (a,b) in TWO_TRI_EDGE_PATTERNS[msh]
                two_tri_add_carrier!(catalog,(cell[a],cell[b]),owner)
            end
            if owner[1]==2
                two_tri_add_carrier!(catalog,cell,owner)
            elseif haskey(TWO_TRI_FACE_PATTERNS,msh)
                for face in TWO_TRI_FACE_PATTERNS[msh]
                    two_tri_add_carrier!(catalog,Tuple(cell[i] for i in face),owner)
                end
            end
        end
    end
    result=Dict{Tuple,Tuple{Int,Int}}()
    for (key,owners) in catalog
        dimension=minimum(first,owners)
        actual=filter(owner->owner[1]==dimension,collect(owners))
        length(actual)==1 || error("two-Tri actual support has ambiguous lower carrier")
        result[key]=only(actual)
    end
    return result
end

function two_tri_saved_carriers(record)
    return two_tri_actual_carriers((Int(row["type"]),Int.(row["nodes"]),
                (Int(entity["dim"]),Int(entity["tag"])))
        for entity in record["p1"]["entities"] for row in entity["cells"])
end

function two_tri_saved_boundary(record,index)
    return Set(QTNN.key(index[tag] for tag in row["nodes"])
        for entity in record["p1"]["entities"] if entity["dim"]==2 for row in entity["cells"])
end

function two_tri_source_signature(source)
    return Set(minimum(Tuple(QTT.point(source,cell[mod1(i+k,3)]) for i in 1:3)
                       for k in 0:2) for cell in eachcol(source.tris))
end

function two_tri_saved_integrity(record,tolerance)
    bytes2hex(sha256(record["input_geo"]))==record["input_sha256"] || error("two-Tri exact input hash")
    # Windows input-file bytes have their own preserved hash; do not substitute
    # the Unicode text hash or invent raw-label equality with native output.
    all(length(record[key])==64 for key in
        ("input_file_sha256","saved_json_sha256")) || error("two-Tri saved provenance hashes")
    for phase in ("p1","p2"),node in record[phase]["nodes"],i in 1:3
        string(reinterpret(UInt64,Float64(node["coordinates"][i]));base=16,pad=16)==
            node["coordinate_bits"][i] || error("two-Tri saved Float64 coordinate bits")
    end
    old=Dict(node["tag"]=>node for node in record["p1"]["nodes"])
    new=Dict(node["tag"]=>node for node in record["p2"]["nodes"])
    pairs=record["primary_tag_remap"]
    length(pairs)==length(old) && length(unique(first.(pairs)))==length(old) &&
        length(unique(last.(pairs)))==length(old) || error("two-Tri primary identity remap")
    inverse=Dict(pair[2]=>pair[1] for pair in pairs)
    for pair in pairs
        old[pair[1]]["coordinate_bits"]==new[pair[2]]["coordinate_bits"] &&
            old[pair[1]]["owner"]==new[pair[2]]["owner"] || error("two-Tri primary identity changed")
    end
    carriers=two_tri_saved_carriers(record)
    for support in record["support_carriers"]
        ids=support["primary_support"]
        key=Tuple(sort!([inverse[id] for id in ids]))
        expected=carriers[key]
        Tuple(support["carrier"])==expected && Tuple(new[support["node"]]["owner"])==expected ||
            error("two-Tri saved support ownership is not actual classified incidence")
        mean=[sum(Float64(new[id]["coordinates"][i]) for id in ids)/length(ids) for i in 1:3]
        maximum(abs.(Float64.(new[support["node"]]["coordinates"]).-mean))<=tolerance ||
            error("two-Tri saved support is not affine P1 interpolation")
    end
    for inverse_record in record["computed_surface_parameters"]
        node=new[inverse_record["node"]]
        length(inverse_record["computed_uv"])==2 &&
            maximum(abs.(Float64.(inverse_record["evaluated_coordinates"]).-
                         Float64.(node["coordinates"])))<=tolerance || error("saved computed UV roundtrip")
    end
    return new
end

function two_tri_native_supports(api,linear,linear_data,projected,msh,volume,tolerance)
    types,tags,connections=linear_data
    actual_types,actual_tags,actual_connections=api.mesh.get_elements(3,volume)
    length(types)==length(actual_types)==1 && length(only(tags))==length(only(actual_tags)) ||
        error("two-Tri native P2 changed primary cells")
    primary=msh_spec(only(types)).nnodes
    before=reshape(only(connections),primary,:)
    after=reshape(only(actual_connections),msh_spec(msh).nnodes,:)
    inverse=Dict{UInt64,Int}()
    for column in axes(before,2),row in 1:primary
        current=after[row,column];original=Int(before[row,column])
        get(inverse,current,original)==original || error("two-Tri native primary support identity")
        inverse[current]=original
    end
    length(inverse)==nnodes(linear) && length(unique(values(inverse)))==nnodes(linear) ||
        error("two-Tri native primary identity is not bijective")
    projected.coords==linear.coords || error("two-Tri actual projection changed primary IDs")
    data=projected.entity_data
    carriers=two_tri_actual_carriers((Int(block.msh),Tuple(cell),
                (Tessella.Elements.msh_dimension(block.msh),Int(data.block_entities[b][column])))
        for (b,block) in enumerate(projected.blocks) for (column,cell) in enumerate(eachcol(block.nodes)))
    supports=Dict{UInt64,Tuple}()
    for edge in eachcol(reshape(api.mesh.get_element_edge_nodes(msh,volume,false),3,:))
        key=Tuple(sort!([inverse[id] for id in edge[1:2]]))
        get(supports,edge[3],key)==key || error("two-Tri native shared P2 edge identity")
        supports[edge[3]]=key
    end
    if msh==13
        for face in eachcol(reshape(api.mesh.get_element_face_nodes(msh,4,volume,false),9,:))
            key=Tuple(sort!([inverse[id] for id in face[1:4]]))
            get(supports,face[9],key)==key || error("two-Tri native shared P2 face identity")
            supports[face[9]]=key
        end
    end
    for (node,support) in supports
        p,_,dim,entity=api.mesh.get_node(node)
        (dim,entity)==carriers[support] || error("two-Tri native P2 carrier differs from its actual P1 boundary")
        expected=[sum(linear.coords[i,id] for id in support)/length(support) for i in 1:3]
        maximum(abs.(p.-expected))<=tolerance || error("two-Tri native P2 support geometry")
    end
    return length(supports)
end

function two_tri_saved_case(record;oracle_only=false)
    tolerance=2e-11;n=record["intervals"];laterals=record["recombine_laterals"]
    levels=record["intended_levels"]
    layers=record["profile"]=="L1" ? :one : record["profile"]=="L3" ? :three : :graded
    f=QTT.fixture(record["name"];height=record["direction"],laterals,layers,pins=(1,2,3,4))
    f=merge(f,(;source=record["input_geo"]))
    saved_nodes=two_tri_saved_integrity(record,tolerance)
    owner_counts=[count(node->node["owner"][1]==dim,values(saved_nodes)) for dim in 0:3]
    owner_counts==[8,8n+4,8n-2,2n-1] &&
        owner_counts==[record["p2_owner_counts"][string(dim)] for dim in 0:3] ||
        error("two-Tri saved owner dimension counts")
    source,_=two_tri_saved_mesh(record,"p1",2,record["source_tag"])
    linear,index=two_tri_saved_mesh(record,"p1",3,record["volume_tag"])
    quadratic,_=two_tri_saved_mesh(record,"p2",3,record["volume_tag"])
    saved=QTT.certify(linear,f,source)
    two_tri_saved_boundary(record,index)==Set(keys(saved.boundary)) || error("two-Tri saved lower/volume boundary mismatch")
    QTNN.certify_quadratic(quadratic,linear)
    nnodes(quadratic)==18n+9 || error("two-Tri saved P2 node count")
    empty_parameters=sum(length(row["empty_owned_uv_nodes"]) for row in record["laterals"])
    empty_parameters==(laterals ? 4n : 0) || error("two-Tri saved stored UV provenance")
    for lateral in record["laterals"]
        lateral["owned"]==2n-1 && lateral["closure"]==6n+3 || error("two-Tri saved lateral ownership")
    end
    if !oracle_only
        execution=QTT.execute(f.source)
        volume=geo_entity_mesh(execution,3,record["volume_tag"])
        actual_source=geo_entity_mesh(execution,2,record["source_tag"])
        two_tri_source_signature(actual_source)==two_tri_source_signature(source) || error("two-Tri oriented source differs from pinned source")
        native=QTT.certify(volume,f,actual_source)
        QTNN.surface_faces(execution.mesh_parts,volume)==Set(keys(native.boundary)) || error("two-Tri native finalized lower/volume boundary mismatch")
        QTNN.counts(volume)==QTNN.counts(linear) || error("two-Tri saved/native volume family counts")
        nnodes(volume)==nnodes(linear) || error("two-Tri saved/native primary node counts")
        for p in eachcol(volume.coords)
            any(q->maximum(abs.(p.-q))<=tolerance,eachcol(linear.coords)) || error("two-Tri swept grid geometry")
        end
        mktempdir() do directory
            path=joinpath(directory,"two_tri_saved.geo");write(path,f.source)
            native_api(path,api->begin
                primary_cache=api.mesh.generate(3)
                first=QTT.public_volume(api,record["volume_tag"])
                primary_data=api.mesh.get_elements(3,record["volume_tag"])
                QTNN.typed_signature(first)==QTNN.typed_signature(volume) || error("two-Tri API/GEO primary product")
                isempty(api.mesh.get_elements(2)[1]) || error("two-Tri API3 lower-cell contract changed")
                projected=model_to_mixed(api.CURRENT[],primary_cache,3,record["volume_tag"])
                api.mesh.set_order(2)
                actual=QTT.public_volume(api,record["volume_tag"])
                QTNN.certify_quadratic(actual,first)
                quadratic_volume_from_jacobians(actual,api.mesh;expected=1.)
                two_tri_native_supports(api,primary_cache,primary_data,projected,laterals ? 13 : 11,
                    record["volume_tag"],tolerance)==14n+5 || error("two-Tri native interpolation support count")
                owners=zeros(Int,4)
                for node in api.mesh.get_nodes()[1]
                    p,uv,dim,entity=api.mesh.get_node(node);owners[dim+1]+=1
                    if dim==2
                        length(uv)==2 && maximum(abs.(api.model.get_value(2,entity,uv).-p))<=tolerance ||
                            error("two-Tri native computed surface UV")
                    end
                end
                owners==[8,8n+4,8n-2,2n-1] || error("two-Tri native owner dimension counts")
                for lateral in record["lateral_tags"]
                    own,_,uv=api.mesh.get_nodes(2,lateral,false,true)
                    length(own)==2n-1 && length(uv)==2length(own) || error("two-Tri native lateral owned query")
                    closure,coords,uv=api.mesh.get_nodes(2,lateral,true,true)
                    length(closure)==6n+3 && length(uv)==2length(closure) &&
                        maximum(abs.(api.model.get_value(2,lateral,uv).-coords))<=tolerance ||
                            error("two-Tri native lateral closure/computed UV")
                end
                api.mesh.set_order(1)
                QTNN.typed_signature(QTT.public_volume(api,record["volume_tag"]))==QTNN.typed_signature(first) ||
                    error("two-Tri actual P2/P1 cell roundtrip")
            end)
        end
    end
    return empty_parameters
end
# Saved patch cases extend the gate without regenerating any saved oracle.
# These catalogs are independent literal MSH primary support patterns.
const QUAD_PATCH_EDGE_PATTERNS=merge(TWO_TRI_EDGE_PATTERNS,Dict(
    5=>((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8)),
    7=>((1,2),(2,3),(3,4),(4,1),(1,5),(2,5),(3,5),(4,5))))
const QUAD_PATCH_FACE_PATTERNS=merge(TWO_TRI_FACE_PATTERNS,Dict(
    5=>((1,4,3,2),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8)),
    7=>((1,4,3,2),(1,2,5),(2,3,5),(3,4,5),(4,1,5))))

function quad_patch_actual_carriers(rows)
    catalog=Dict{Tuple,Set{Tuple{Int,Int}}}()
    for (msh,cell,owner) in rows
        if msh==1
            two_tri_add_carrier!(catalog,cell,owner)
        elseif haskey(QUAD_PATCH_EDGE_PATTERNS,msh)
            for edge in QUAD_PATCH_EDGE_PATTERNS[msh]
                two_tri_add_carrier!(catalog,Tuple(cell[i] for i in edge),owner)
            end
            if owner[1]==2
                two_tri_add_carrier!(catalog,cell,owner)
            elseif haskey(QUAD_PATCH_FACE_PATTERNS,msh)
                for face in QUAD_PATCH_FACE_PATTERNS[msh]
                    two_tri_add_carrier!(catalog,Tuple(cell[i] for i in face),owner)
                end
                msh==5 && two_tri_add_carrier!(catalog,cell,owner)
            end
        end
    end
    result=Dict{Tuple,Tuple{Int,Int}}()
    for (key,owners) in catalog
        dimension=minimum(first,owners)
        actual=filter(owner->owner[1]==dimension,collect(owners))
        length(actual)==1 || error("quad-patch actual support has ambiguous lower carrier")
        result[key]=only(actual)
    end
    return result
end

function quad_patch_saved_integrity(record,tolerance)
    bytes2hex(sha256(record["input_geo"]))==record["input_sha256"] || error("quad-patch exact input hash")
    length(record["saved_json_sha256"])==64 || error("quad-patch saved provenance hash")
    for phase in ("p1","p2"),node in record[phase]["nodes"]
        for (values,bits) in ((node["coordinates"],node["coordinate_bits"]),
                              (node["stored_parameters"],node["stored_parameter_bits"]))
            length(values)==length(bits) || error("quad-patch packed node payload")
            for i in eachindex(values)
                string(reinterpret(UInt64,Float64(values[i]));base=16,pad=16)==bits[i] ||
                    error("quad-patch saved Float64 payload bits")
            end
        end
    end
    old=Dict(node["tag"]=>node for node in record["p1"]["nodes"])
    new=Dict(node["tag"]=>node for node in record["p2"]["nodes"])
    pairs=record["primary_remap"]
    length(pairs)==length(old) && length(unique(first.(pairs)))==length(old) &&
        length(unique(last.(pairs)))==length(old) || error("quad-patch primary identity remap")
    inverse=Dict(pair[2]=>pair[1] for pair in pairs)
    for pair in pairs
        old[pair[1]]["coordinate_bits"]==new[pair[2]]["coordinate_bits"] &&
            old[pair[1]]["owner"]==new[pair[2]]["owner"] || error("quad-patch primary identity changed")
    end
    carriers=quad_patch_actual_carriers((Int(row["type"]),Int.(row["nodes"]),
        (Int(entity["dim"]),Int(entity["tag"])))
        for entity in record["p1"]["entities"] for row in entity["cells"])
    supports=record["supports"]
    Set(row["node"] for row in supports)==setdiff(Set(keys(new)),Set(keys(inverse))) ||
        error("quad-patch saved interpolation support coverage")
    for support in supports
        ids=support["primary"];key=Tuple(sort!([inverse[id] for id in ids]))
        expected=carriers[key]
        Tuple(support["owner"])==expected && Tuple(new[support["node"]]["owner"])==expected ||
            error("quad-patch saved support owner differs from actual incidence")
        mean=[sum(Float64(new[id]["coordinates"][i]) for id in ids)/length(ids) for i in 1:3]
        maximum(abs.(Float64.(new[support["node"]]["coordinates"]).-mean))<=tolerance ||
            error("quad-patch saved support geometry")
    end
    for surface in record["surface_queries"]
        queries=haskey(surface,"inverse") ? surface["inverse"] : surface["computed"]
        for query in queries
            uv=haskey(query,"computed_uv") ? query["computed_uv"] : query["uv"]
            length(uv)==2 && maximum(abs.(Float64.(query["evaluated"]).-
                Float64.(new[query["node"]]["coordinates"])))<=tolerance || error("quad-patch saved computed UV")
        end
    end
    return new
end

function quad_patch_native_supports(api,linear,linear_data,projected,volume,tolerance)
    types,_,connections=linear_data
    actual_types,_,actual_connections=api.mesh.get_elements(3,volume)
    targets=Dict(4=>11,5=>12,7=>14);inverse=Dict{UInt64,Int}()
    length(types)==length(actual_types) || error("quad-patch native P2 family change")
    for (msh,nodes) in zip(types,connections)
        target=targets[Int(msh)];slot=findfirst(==(target),actual_types)
        slot===nothing && error("quad-patch native P2 family missing")
        before=reshape(nodes,msh_spec(msh).nnodes,:)
        after=reshape(actual_connections[slot],msh_spec(target).nnodes,:)
        size(before,2)==size(after,2) || error("quad-patch native P2 cell count change")
        for column in axes(before,2),row in axes(before,1)
            current=after[row,column];original=Int(before[row,column])
            get(inverse,current,original)==original || error("quad-patch native primary identity conflict")
            inverse[current]=original
        end
    end
    length(inverse)==nnodes(linear) && Set(values(inverse))==Set(1:nnodes(linear)) ||
        error("quad-patch native primary identity is not bijective")
    projected.coords==linear.coords || error("quad-patch actual projection changed primary IDs")
    data=projected.entity_data
    carriers=quad_patch_actual_carriers((Int(block.msh),Tuple(cell),
        (Tessella.Elements.msh_dimension(block.msh),Int(data.block_entities[b][column])))
        for (b,block) in enumerate(projected.blocks) for (column,cell) in enumerate(eachcol(block.nodes)))
    supports=Dict{UInt64,Tuple}()
    function add_support(node,primary)
        key=Tuple(sort!([inverse[id] for id in primary]))
        get(supports,node,key)==key || error("quad-patch native shared support identity")
        supports[node]=key
    end
    for msh in actual_types
        for edge in eachcol(reshape(api.mesh.get_element_edge_nodes(msh,volume,false),3,:))
            add_support(edge[3],edge[1:2])
        end
        if msh in (12,14)
            for face in eachcol(reshape(api.mesh.get_element_face_nodes(msh,4,volume,false),9,:))
                add_support(face[9],face[1:4])
            end
        end
        if msh==12
            slot=only(findall(==(msh),actual_types))
            for cell in eachcol(reshape(actual_connections[slot],27,:));add_support(cell[27],cell[1:8]);end
        end
    end
    Set(keys(supports))==setdiff(Set(api.mesh.get_nodes()[1]),Set(keys(inverse))) ||
        error("quad-patch native interpolation support coverage")
    for (node,support) in supports
        p,_,dim,entity=api.mesh.get_node(node)
        (dim,entity)==carriers[support] || error("quad-patch native P2 owner differs from actual P1 carriers")
        expected=[sum(linear.coords[i,id] for id in support)/length(support) for i in 1:3]
        maximum(abs.(p.-expected))<=tolerance || error("quad-patch native P2 support geometry")
    end
    return length(supports)
end

function quad_patch_source_match(native,saved,tolerance)
    ids=Dict{Int,Int}()
    for node in 1:nnodes(native)
        matches=findall(i->maximum(abs.(native.coords[:,node].-saved.coords[:,i]))<=tolerance,1:nnodes(saved))
        length(matches)==1 || error("quad-patch source geometry is not uniquely matched")
        ids[node]=only(matches)
    end
    Set(values(ids))==Set(1:nnodes(saved)) || error("quad-patch source identity match is not bijective")
    actual=Set(QTP.cycle(Tuple(ids[n] for n in cell)) for block in QTP.blocks(native) for cell in eachcol(block.nodes))
    expected=Set(QTP.cycle(Tuple(cell)) for block in QTP.blocks(saved) for cell in eachcol(block.nodes))
    actual==expected || error("quad-patch oriented source complex differs from saved oracle")
    return nothing
end

function quad_patch_saved_case(record;oracle_only=false)
    tolerance=2e-11;n=record["intervals"];laterals=record["recombine_laterals"]
    layers=record["profile"]=="L1" ? :one : record["profile"]=="L3" ? :three : :graded
    shape=record["geometry_shape"]=="unit_square" ? :unit : :rounded
    f=QTP.fixture(record["name"];shape,height=record["direction"],laterals,layers,pins=(1,2,3,4))
    f=merge(f,(;source=record["input_geo"]))
    saved_nodes=quad_patch_saved_integrity(record,tolerance)
    owners=[count(node->node["owner"][1]==dim,values(saved_nodes)) for dim in 0:3]
    owners==[8,8n+20,24n+6,18n-9] || error("quad-patch saved owner dimension counts")
    source,_=two_tri_saved_mesh(record,"p1",2,record["source_tag"])
    linear,index=two_tri_saved_mesh(record,"p1",3,record["volume_tag"])
    quadratic,_=two_tri_saved_mesh(record,"p2",3,record["volume_tag"])
    saved=QTP.certify(linear,f,source)
    two_tri_saved_boundary(record,index)==Set(keys(saved.boundary)) || error("quad-patch saved typed exterior mismatch")
    QTNN.certify_quadratic(quadratic,linear)
    nnodes(quadratic)==50n+25 || error("quad-patch saved P2 node count")
    empty_parameters=count(node->node["owner"][1]==2 && isempty(node["stored_parameters"]),values(saved_nodes))
    empty_parameters==(laterals ? 8n : 0)==record["empty_stored_surface_uv_nodes"] || error("quad-patch stored UV provenance count")
    if !oracle_only
        execution=QTP.execute(f.source)
        volume=geo_entity_mesh(execution,3,record["volume_tag"])
        actual_source=geo_entity_mesh(execution,2,record["source_tag"])
        quad_patch_source_match(actual_source,source,tolerance)
        native=QTP.certify(volume,f,actual_source)
        QTNN.surface_faces(execution.mesh_parts,volume)==Set(keys(native.boundary)) || error("quad-patch native typed exterior mismatch")
        QTNN.counts(volume)==QTNN.counts(linear) && nnodes(volume)==nnodes(linear) || error("quad-patch native primary product counts")
        for p in eachcol(volume.coords)
            count(q->maximum(abs.(p.-q))<=tolerance,eachcol(linear.coords))==1 || error("quad-patch actual column geometry")
        end
        mktempdir() do directory
            path=joinpath(directory,"quad_patch_saved.geo");write(path,f.source)
            native_api(path,api->begin
                primary_cache=api.mesh.generate(3)
                first=QTP.public_volume(api,record["volume_tag"])
                primary_data=api.mesh.get_elements(3,record["volume_tag"])
                QTNN.typed_signature(first)==QTNN.typed_signature(volume) || error("quad-patch API/GEO primary cells")
                isempty(api.mesh.get_elements(2)[1]) || error("quad-patch API3 lower-cell contract changed")
                projected=model_to_mixed(api.CURRENT[],primary_cache,3,record["volume_tag"])
                api.mesh.set_order(2)
                actual=QTP.public_volume(api,record["volume_tag"])
                QTNN.certify_quadratic(actual,first)
                quadratic_volume_from_jacobians(actual,api.mesh;expected=Float64(native.total))
                quad_patch_native_supports(api,primary_cache,primary_data,projected,record["volume_tag"],tolerance)==41n+16 ||
                    error("quad-patch native interpolation support count")
                native_owners=zeros(Int,4)
                for node in api.mesh.get_nodes()[1]
                    p,uv,dim,entity=api.mesh.get_node(node);native_owners[dim+1]+=1
                    if dim==2
                        length(uv)==2 && maximum(abs.(api.model.get_value(2,entity,uv).-p))<=tolerance ||
                            error("quad-patch native computed surface UV")
                    end
                end
                native_owners==owners || error("quad-patch native owner dimension counts")
                for lateral in record["lateral_tags"]
                    own,_,uv=api.mesh.get_nodes(2,lateral,false,true)
                    length(own)==6n-3 && length(uv)==2length(own) || error("quad-patch lateral owned query")
                    closure,coords,uv=api.mesh.get_nodes(2,lateral,true,true)
                    length(closure)==10n+5 && length(uv)==2length(closure) &&
                        maximum(abs.(api.model.get_value(2,lateral,uv).-coords))<=tolerance || error("quad-patch lateral closure/computed UV")
                end
                for cap in (record["source_tag"],record["top_tag"])
                    own,_,uv=api.mesh.get_nodes(2,cap,false,true)
                    length(own)==9 && length(uv)==18 || error("quad-patch cap owned query")
                    closure,coords,uv=api.mesh.get_nodes(2,cap,true,true)
                    length(closure)==25 && length(uv)==50 &&
                        maximum(abs.(api.model.get_value(2,cap,uv).-coords))<=tolerance || error("quad-patch cap closure/computed UV")
                end
                api.mesh.set_order(1)
                QTNN.typed_signature(QTP.public_volume(api,record["volume_tag"]))==QTNN.typed_signature(first) ||
                    error("quad-patch actual P2/P1 cell roundtrip")
            end)
        end
    end
    return empty_parameters
end

# Strip-specific exact integration uses the actual reference map. A warped
# bilinear Prism6 face cannot be replaced by an arbitrary flat tetra split.
const STRIP_Q=Rational{BigInt}
const STRIP_HEX_SIGNS=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
                       (0,0,1),(1,0,1),(1,1,1),(0,1,1))
strip_det(a,b,c)=a[1]*(b[2]*c[3]-b[3]*c[2])-a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1])
strip_sub(a,b)=ntuple(i->a[i]-b[i],3)
strip_tet(a,b,c,d)=strip_det(strip_sub(b,a),strip_sub(c,a),strip_sub(d,a))

function strip_reference_volume(msh,points)
    v=Tuple(ntuple(i->STRIP_Q(p[i]),3) for p in points)
    msh==4 && return strip_tet(v...)/6
    if msh==7
        return sum(strip_tet(v[k],v[mod1(k+1,4)],v[mod1(k-1,4)],v[5]) for k in 1:4)/12
    elseif msh==6
        a=strip_sub(v[2],v[1]);b=strip_sub(v[3],v[1])
        A=strip_sub(strip_sub(v[5],v[4]),a);B=strip_sub(strip_sub(v[6],v[4]),b)
        total=zero(STRIP_Q)
        for k in 1:3
            c=strip_sub(v[k+3],v[k])
            q0=strip_det(a,b,c)
            q1=strip_det(A,b,c)+strip_det(a,B,c)
            q2=strip_det(A,B,c)
            q0>0 && q0+q1+q2>0 || error("strip nonpositive Prism6 endpoint map")
            if q2>0 && -2q2<q1<0
                4q0*q2-q1*q1>0 || error("strip nonpositive Prism6 interior map")
            end
            total+=(q0+q1/2+q2/3)/6
        end
        return total
    elseif msh==5
        # det(J) has degree at most2 per unit-cube axis: exact tensor Simpson.
        samples=((STRIP_Q(0),STRIP_Q(1//6)),(STRIP_Q(1//2),STRIP_Q(2//3)),(STRIP_Q(1),STRIP_Q(1//6)))
        total=zero(STRIP_Q)
        for (u,wu) in samples,(s,ws) in samples,(w,ww) in samples
            x=(u,s,w)
            jac=ntuple(3) do axis
                ntuple(3) do d
                    sum((v[k][d]-v[1][d])*(STRIP_HEX_SIGNS[k][axis]==0 ? -1 : 1)*
                        prod(STRIP_HEX_SIGNS[k][j]==0 ? 1-x[j] : x[j] for j in 1:3 if j!=axis) for k in 1:8)
                end
            end
            total+=wu*ws*ww*strip_det(jac...)
        end
        return total
    end
    error("strip unsupported P1 reference volume family $msh")
end

function strip_integrated_volume(mesh)
    return sum(strip_reference_volume(block.msh,Tuple(QTNN.point(mesh,n) for n in cell))
               for block in QTS.blocks(mesh) if Tessella.Elements.msh_dimension(block.msh)==3
               for cell in eachcol(block.nodes))
end

function strip_saved_fixture(record)
    phase=record["p1"];source=only(filter(e->e["dim"]==2 && e["tag"]==record["source_tag"],phase["entities"]))
    nodes=Dict(p["tag"]=>p for p in phase["nodes"])
    successors=Dict{Int,Int}()
    for (dimension,signed) in source["oriented_boundary"]
        dimension==1 || error("strip source CAD boundary dimension")
        curve=only(filter(e->e["dim"]==1 && e["tag"]==abs(signed),phase["entities"]))
        cells=curve["cells"];a=first(cells)["nodes"][1];b=last(cells)["nodes"][2]
        signed<0 && ((a,b)=(b,a))
        haskey(successors,a) && error("strip repeated CAD corner successor")
        successors[a]=b
    end
    length(successors)==4 && Set(keys(successors))==Set(values(successors)) || error("strip CAD corner cycle")
    first_corner=argmin(n->nodes[n]["owner"][2],collect(keys(successors)))
    corner_ids=ntuple(4) do k
        node=first_corner
        for _ in 2:k;node=successors[node];end
        node
    end
    successors[last(corner_ids)]==first_corner && length(unique(corner_ids))==4 || error("strip CAD corner cycle closure")
    corners=Tuple(Tuple(Float64.(nodes[n]["coordinates"])) for n in corner_ids)
    axis=Int(record["normal_axis"]);other=Tuple(i for i in 1:3 if i!=axis);axes=(other...,axis)
    winding=sign((corners[2][other[1]]-corners[1][other[1]])*(corners[3][other[2]]-corners[1][other[2]])-
                 (corners[2][other[2]]-corners[1][other[2]])*(corners[3][other[1]]-corners[1][other[1]]))
    plane=record["plane"]=="XZ" ? :ZX : Symbol(record["plane"])
    f=QTS.fixture(record["name"];plane,shape=record["geometry_shape"]=="unit_square" ? :unit : :rounded,
                  height=record["translation_height"],laterals=record["recombine_laterals"],layers=:one)
    return merge(f,(;source=record["input_geo"],axes,corners,winding=Int(winding),
                    levels=Float64.(record["intended_levels"]),height=Float64(record["translation_height"]),
                    surface=Int(record["source_tag"]),curve_tags=Tuple(Int.(record["source_curve_tags"])),
                    point_tags=Tuple(Int(nodes[n]["owner"][2]) for n in corner_ids)))
end

function strip_saved_reference_integrity(record)
    for phase in ("p1","p2")
        certificates=record["reference_map_certificates"][phase]
        cells=[(e["dim"],e["tag"],c["tag"],c["type"]) for e in record[phase]["entities"] if e["dim"]>0 for c in e["cells"]]
        Set((r["entity"][1],r["entity"][2],r["cell"],r["type"]) for r in certificates)==Set(cells) ||
            error("strip saved whole-reference certificate coverage")
        for r in certificates
            pieces=split(r["minimum_exact"],'/')
            q=length(pieces)==1 ? STRIP_Q(parse(BigInt,only(pieces))) : parse(BigInt,pieces[1])//parse(BigInt,pieces[2])
            q>0 && Float64(q)==r["minimum"] || error("strip saved whole-reference bound")
        end
    end
    for entry in record["p2_macro_reference_integration"]
        entry["exact_actual_p2_map_volume"]==entry["exact_actual_p2_boundary_flux"] ||
            error("strip actual P2 macro map/boundary integration")
    end
    return nothing
end

function strip_fraction(text)
    parts=split(text,'/')
    return length(parts)==1 ? STRIP_Q(parse(BigInt,only(parts))) :
        parse(BigInt,parts[1])//parse(BigInt,parts[2])
end

function strip_native_supports(api,linear,linear_data,projected,entity,tolerance)
    types,_,connections=linear_data
    actual_types,_,actual_connections=api.mesh.get_elements(3,entity)
    targets=Dict(4=>11,5=>12,6=>13,7=>14);inverse=Dict{UInt64,Int}()
    length(types)==length(actual_types) || error("strip native P2 family change")
    for (msh,nodes) in zip(types,connections)
        target=targets[Int(msh)];slot=findfirst(==(target),actual_types)
        slot===nothing && error("strip native P2 family missing")
        before=reshape(nodes,msh_spec(msh).nnodes,:)
        after=reshape(actual_connections[slot],msh_spec(target).nnodes,:)
        size(before,2)==size(after,2) || error("strip native P2 cell count change")
        for column in axes(before,2),row in axes(before,1)
            current=after[row,column];original=Int(before[row,column])
            get(inverse,current,original)==original || error("strip native primary identity conflict")
            inverse[current]=original
        end
    end
    length(inverse)==nnodes(linear) && Set(values(inverse))==Set(1:nnodes(linear)) ||
        error("strip native primary identity is not bijective")
    projected.coords==linear.coords || error("strip actual projection changed primary IDs")
    data=projected.entity_data
    carriers=quad_patch_actual_carriers((Int(block.msh),Tuple(cell),
        (Tessella.Elements.msh_dimension(block.msh),Int(data.block_entities[b][column])))
        for (b,block) in enumerate(projected.blocks) for (column,cell) in enumerate(eachcol(block.nodes)))
    supports=Dict{UInt64,Tuple}()
    function add_support(node,primary)
        support=Tuple(sort!([inverse[id] for id in primary]))
        get(supports,node,support)==support || error("strip native shared support identity")
        supports[node]=support
    end
    for msh in actual_types
        for edge in eachcol(reshape(api.mesh.get_element_edge_nodes(msh,entity,false),3,:))
            add_support(edge[3],edge[1:2])
        end
        if msh in (12,13,14)
            for face in eachcol(reshape(api.mesh.get_element_face_nodes(msh,4,entity,false),9,:))
                add_support(face[9],face[1:4])
            end
        end
        if msh==12
            slot=only(findall(==(msh),actual_types))
            for cell in eachcol(reshape(actual_connections[slot],27,:))
                add_support(cell[27],cell[1:8])
            end
        end
    end
    Set(keys(supports))==setdiff(Set(api.mesh.get_nodes()[1]),Set(keys(inverse))) ||
        error("strip native interpolation support coverage")
    for (node,support) in supports
        p,_,dim,owner=api.mesh.get_node(node)
        (dim,owner)==carriers[support] || error("strip native P2 owner differs from actual P1 carrier")
        expected=[sum(linear.coords[i,id] for id in support)/length(support) for i in 1:3]
        maximum(abs.(p.-expected))<=tolerance || error("strip native P2 support geometry")
    end
    return length(supports)
end

function strip_saved_case(record;oracle_only=false)
    tolerance=2e-11;n=Int(record["intervals"]);laterals=record["recombine_laterals"]
    f=strip_saved_fixture(record)
    saved_nodes=quad_patch_saved_integrity(record,tolerance)
    strip_saved_reference_integrity(record)
    owners=[count(node->node["owner"][1]==dim,values(saved_nodes)) for dim in 0:3]
    owners==[8,8n+12,16n-2,6n+(laterals ? 13 : -3)] || error("strip saved owner dimension counts")
    source,_=two_tri_saved_mesh(record,"p1",2,record["source_tag"])
    linear,index=two_tri_saved_mesh(record,"p1",3,record["volume_tag"])
    quadratic,_=two_tri_saved_mesh(record,"p2",3,record["volume_tag"])
    saved=QTS.certify(linear,f,source;oracle=true)
    two_tri_saved_boundary(record,index)==Set(keys(saved.boundary)) || error("strip saved actual typed exterior")
    exact_saved_volume=strip_integrated_volume(linear)
    exact_saved_volume==saved.total && exact_saved_volume==
        sum(strip_fraction(row["exact_partition_volume"]) for row in record["strip_macro_certificate"]["macro_cells"]) ||
        error("strip saved actual reference integration differs from macro partition")
    QTNN.certify_quadratic(quadratic,linear)
    nnodes(quadratic)==30n+(laterals ? 31 : 15) || error("strip saved P2 node count")
    empty_parameters=count(node->node["owner"][1]==2 && isempty(node["stored_parameters"]),values(saved_nodes))
    empty_parameters==(laterals ? 6n : 0)==record["empty_stored_surface_uv_nodes"] ||
        error("strip saved stored UV provenance count")
    for (lateral,width) in zip(record["lateral_tags"],record["lateral_chain_widths"])
        entity=only(e for e in record["p2"]["entities"] if e["dim"]==2 && e["tag"]==lateral)
        length(entity["owned"]["tags"])==(width==3 ? 6n-3 : 2n-1) &&
            length(entity["closure"]["tags"])==(width==3 ? 10n+5 : 6n+3) ||
            error("strip saved lateral owner/closure counts")
    end
    if !oracle_only
        execution=QTS.execute(f.source)
        volume=geo_entity_mesh(execution,3,record["volume_tag"])
        actual_source=geo_entity_mesh(execution,2,record["source_tag"])
        quad_patch_source_match(actual_source,source,tolerance)
        native=QTS.certify(volume,f,actual_source)
        strip_integrated_volume(volume)==native.total || error("strip native actual reference integration")
        QTNN.surface_faces(execution.mesh_parts,volume)==Set(keys(native.boundary)) || error("strip native finalized typed exterior")
        nnodes(volume)==nnodes(linear) || error("strip native primary node count")
        # Free finalized masks can choose different valid existing-node factories.
        # Every primary column and actual terminal centroid must still agree.
        for p in eachcol(volume.coords)
            count(q->maximum(abs.(p.-q))<=tolerance,eachcol(linear.coords))==1 || error("strip actual primary/centroid geometry")
        end
        mktempdir() do directory
            path=joinpath(directory,"quad_strip_saved.geo");write(path,f.source)
            native_api(path,api->begin
                cache=api.mesh.generate(3)
                first=QTS.public_volume(api,record["volume_tag"])
                primary_data=api.mesh.get_elements(3,record["volume_tag"])
                QTNN.typed_signature(first)==QTNN.typed_signature(volume) || error("strip API/GEO actual primary cells")
                isempty(api.mesh.get_elements(2)[1]) || error("strip API3 lower-cell contract changed")
                projected=model_to_mixed(api.CURRENT[],cache,3,record["volume_tag"])
                api.mesh.set_order(2)
                actual=QTS.public_volume(api,record["volume_tag"])
                QTNN.certify_quadratic(actual,first)
                quadratic_volume_from_jacobians(actual,api.mesh;expected=Float64(native.total))
                strip_native_supports(api,cache,primary_data,projected,record["volume_tag"],tolerance)==
                    24n+(laterals ? 23 : 9) || error("strip native P2 support count")
                native_owners=zeros(Int,4)
                for node in api.mesh.get_nodes()[1]
                    p,uv,dim,owner=api.mesh.get_node(node);native_owners[dim+1]+=1
                    if dim==2
                        length(uv)==2 && maximum(abs.(api.model.get_value(2,owner,uv).-p))<=tolerance ||
                            error("strip native computed surface UV")
                    end
                end
                native_owners==owners || error("strip native owner dimension counts")
                for (lateral,width) in zip(record["lateral_tags"],record["lateral_chain_widths"])
                    own,_,uv=api.mesh.get_nodes(2,lateral,false,true)
                    length(own)==(width==3 ? 6n-3 : 2n-1) && length(uv)==2length(own) ||
                        error("strip native lateral owned query")
                    closure,coords,uv=api.mesh.get_nodes(2,lateral,true,true)
                    length(closure)==(width==3 ? 10n+5 : 6n+3) && length(uv)==2length(closure) &&
                        maximum(abs.(api.model.get_value(2,lateral,uv).-coords))<=tolerance ||
                        error("strip native lateral closure/computed UV")
                end
                for cap in (record["source_tag"],record["top_tag"])
                    own,_,uv=api.mesh.get_nodes(2,cap,false,true)
                    length(own)==3 && length(uv)==6 || error("strip native cap owned query")
                    closure,coords,uv=api.mesh.get_nodes(2,cap,true,true)
                    length(closure)==15 && length(uv)==30 &&
                        maximum(abs.(api.model.get_value(2,cap,uv).-coords))<=tolerance ||
                        error("strip native cap closure/computed UV")
                end
                api.mesh.set_order(1)
                QTNN.typed_signature(QTS.public_volume(api,record["volume_tag"]))==QTNN.typed_signature(first) ||
                    error("strip actual P2/P1 roundtrip")
            end)
        end
    end
    return empty_parameters
end

# Three-Quad saved products retain their own counters and parameter provenance.
function three_strip_saved_fixture(record)
    phase=record["p1"];source=only(filter(e->e["dim"]==2 && e["tag"]==record["source_tag"],phase["entities"]))
    nodes=Dict(p["tag"]=>p for p in phase["nodes"])
    successors=Dict{Int,Int}()
    for (dimension,signed) in source["oriented_boundary"]
        dimension==1 || error("strip source CAD boundary dimension")
        curve=only(filter(e->e["dim"]==1 && e["tag"]==abs(signed),phase["entities"]))
        cells=curve["cells"];a=first(cells)["nodes"][1];b=last(cells)["nodes"][2]
        signed<0 && ((a,b)=(b,a))
        haskey(successors,a) && error("strip repeated CAD corner successor")
        successors[a]=b
    end
    length(successors)==4 && Set(keys(successors))==Set(values(successors)) || error("strip CAD corner cycle")
    first_corner=argmin(n->nodes[n]["owner"][2],collect(keys(successors)))
    corner_ids=ntuple(4) do k
        node=first_corner
        for _ in 2:k;node=successors[node];end
        node
    end
    successors[last(corner_ids)]==first_corner && length(unique(corner_ids))==4 || error("strip CAD corner cycle closure")
    corners=Tuple(Tuple(Float64.(nodes[n]["coordinates"])) for n in corner_ids)
    axis=Int(record["normal_axis"]);other=Tuple(i for i in 1:3 if i!=axis);axes=(other...,axis)
    winding=sign((corners[2][other[1]]-corners[1][other[1]])*(corners[3][other[2]]-corners[1][other[2]])-
                 (corners[2][other[2]]-corners[1][other[2]])*(corners[3][other[1]]-corners[1][other[1]]))
    plane=record["plane"]=="XZ" ? :ZX : Symbol(record["plane"])
    f=QT3S.fixture(record["name"];plane,shape=record["geometry_shape"]=="unit_square" ? :unit : :rounded,
                  height=record["translation_height"],laterals=record["recombine_laterals"],layers=:one)
    return merge(f,(;source=record["input_geo"],axes,corners,winding=Int(winding),
                    levels=Float64.(record["intended_levels"]),height=Float64(record["translation_height"]),
                    surface=Int(record["source_tag"]),curve_tags=Tuple(Int.(record["source_curve_tags"])),
                    point_tags=Tuple(Int(nodes[n]["owner"][2]) for n in corner_ids),
                    intervals=Int(record["intervals"]),translation=ntuple(k->k==axis ? Float64(record["translation_height"]) : 0.,3)))
end


function three_strip_saved_case(record;oracle_only=false)
    tolerance=2e-11;n=Int(record["intervals"]);laterals=record["recombine_laterals"]
    f=three_strip_saved_fixture(record)
    saved_nodes=quad_patch_saved_integrity(record,tolerance)
    strip_saved_reference_integrity(record)
    owners=[count(node->node["owner"][1]==dim,values(saved_nodes)) for dim in 0:3]
    owners==[8,8n+20,24n-2,10n-5+(laterals ? 24 : 0)] || error("three-strip saved owner dimension counts")
    source,_=two_tri_saved_mesh(record,"p1",2,record["source_tag"])
    linear,index=two_tri_saved_mesh(record,"p1",3,record["volume_tag"])
    quadratic,_=two_tri_saved_mesh(record,"p2",3,record["volume_tag"])
    saved=QT3S.certify(linear,f,source;oracle=true)
    two_tri_saved_boundary(record,index)==Set(keys(saved.boundary)) || error("three-strip saved actual typed exterior")
    exact_saved_volume=strip_integrated_volume(linear)
    exact_saved_volume==saved.total && exact_saved_volume==
        sum(strip_fraction(row["exact_partition_volume"]) for row in record["strip_macro_certificate"]["macro_cells"]) ||
        error("three-strip saved actual reference integration differs from macro partition")
    QTNN.certify_quadratic(quadratic,linear)
    nnodes(quadratic)==42n+21+(laterals ? 24 : 0) || error("three-strip saved P2 node count")
    empty_parameters=count(node->node["owner"][1]==2 && isempty(node["stored_parameters"]),values(saved_nodes))
    empty_parameters==(laterals ? 8n : 0)==record["empty_stored_surface_uv_nodes"] ||
        error("three-strip saved stored UV provenance count")
    for (lateral,width) in zip(record["lateral_tags"],record["lateral_chain_widths"])
        entity=only(e for e in record["p2"]["entities"] if e["dim"]==2 && e["tag"]==lateral)
        length(entity["owned"]["tags"])==(width==4 ? 10n-5 : 2n-1) &&
            length(entity["closure"]["tags"])==(width==4 ? 14n+7 : 6n+3) ||
            error("three-strip saved lateral owner/closure counts")
    end
    if !oracle_only
        execution=QT3S.execute(f.source)
        volume=geo_entity_mesh(execution,3,record["volume_tag"])
        actual_source=geo_entity_mesh(execution,2,record["source_tag"])
        quad_patch_source_match(actual_source,source,tolerance)
        native=QT3S.certify(volume,f,actual_source)
        strip_integrated_volume(volume)==native.total || error("three-strip native actual reference integration")
        QTNN.surface_faces(execution.mesh_parts,volume)==Set(keys(native.boundary)) || error("three-strip native finalized typed exterior")
        nnodes(volume)==nnodes(linear) || error("three-strip native primary node count")
        # Free finalized masks can choose different valid existing-node factories.
        # Every primary column and actual terminal centroid must still agree.
        for p in eachcol(volume.coords)
            count(q->maximum(abs.(p.-q))<=tolerance,eachcol(linear.coords))==1 || error("three-strip actual primary/centroid geometry")
        end
        mktempdir() do directory
            path=joinpath(directory,"three_quad_strip_saved.geo");write(path,f.source)
            native_api(path,api->begin
                cache=api.mesh.generate(3)
                first=QT3S.public_volume(api,record["volume_tag"])
                primary_data=api.mesh.get_elements(3,record["volume_tag"])
                QTNN.typed_signature(first)==QTNN.typed_signature(volume) || error("three-strip API/GEO actual primary cells")
                isempty(api.mesh.get_elements(2)[1]) || error("three-strip API3 lower-cell contract changed")
                projected=model_to_mixed(api.CURRENT[],cache,3,record["volume_tag"])
                api.mesh.set_order(2)
                actual=QT3S.public_volume(api,record["volume_tag"])
                QTNN.certify_quadratic(actual,first)
                quadratic_volume_from_jacobians(actual,api.mesh;expected=Float64(native.total))
                strip_native_supports(api,cache,primary_data,projected,record["volume_tag"],tolerance)==
                    34n+13+(laterals ? 21 : 0) || error("three-strip native P2 support count")
                native_owners=zeros(Int,4)
                for node in api.mesh.get_nodes()[1]
                    p,uv,dim,owner=api.mesh.get_node(node);native_owners[dim+1]+=1
                    if dim==2
                        length(uv)==2 && maximum(abs.(api.model.get_value(2,owner,uv).-p))<=tolerance ||
                            error("three-strip native computed surface UV")
                    end
                end
                native_owners==owners || error("three-strip native owner dimension counts")
                for (lateral,width) in zip(record["lateral_tags"],record["lateral_chain_widths"])
                    own,_,uv=api.mesh.get_nodes(2,lateral,false,true)
                    length(own)==(width==4 ? 10n-5 : 2n-1) && length(uv)==2length(own) ||
                        error("three-strip native lateral owned query")
                    closure,coords,uv=api.mesh.get_nodes(2,lateral,true,true)
                    length(closure)==(width==4 ? 14n+7 : 6n+3) && length(uv)==2length(closure) &&
                        maximum(abs.(api.model.get_value(2,lateral,uv).-coords))<=tolerance ||
                        error("three-strip native lateral closure/computed UV")
                end
                for cap in (record["source_tag"],record["top_tag"])
                    own,_,uv=api.mesh.get_nodes(2,cap,false,true)
                    length(own)==5 && length(uv)==10 || error("three-strip native cap owned query")
                    closure,coords,uv=api.mesh.get_nodes(2,cap,true,true)
                    length(closure)==21 && length(uv)==42 &&
                        maximum(abs.(api.model.get_value(2,cap,uv).-coords))<=tolerance ||
                        error("three-strip native cap closure/computed UV")
                end
                api.mesh.set_order(1)
                QTNN.typed_signature(QT3S.public_volume(api,record["volume_tag"]))==QTNN.typed_signature(first) ||
                    error("three-strip actual P2/P1 roundtrip")
            end)
        end
    end
    return empty_parameters
end


# Four-Quad saved products retain their own counters and parameter provenance.
function four_strip_saved_fixture(record)
    phase=record["p1"];source=only(filter(e->e["dim"]==2 && e["tag"]==record["source_tag"],phase["entities"]))
    nodes=Dict(p["tag"]=>p for p in phase["nodes"])
    successors=Dict{Int,Int}()
    for (dimension,signed) in source["oriented_boundary"]
        dimension==1 || error("strip source CAD boundary dimension")
        curve=only(filter(e->e["dim"]==1 && e["tag"]==abs(signed),phase["entities"]))
        cells=curve["cells"];a=first(cells)["nodes"][1];b=last(cells)["nodes"][2]
        signed<0 && ((a,b)=(b,a))
        haskey(successors,a) && error("strip repeated CAD corner successor")
        successors[a]=b
    end
    length(successors)==4 && Set(keys(successors))==Set(values(successors)) || error("strip CAD corner cycle")
    first_corner=argmin(n->nodes[n]["owner"][2],collect(keys(successors)))
    corner_ids=ntuple(4) do k
        node=first_corner
        for _ in 2:k;node=successors[node];end
        node
    end
    successors[last(corner_ids)]==first_corner && length(unique(corner_ids))==4 || error("strip CAD corner cycle closure")
    corners=Tuple(Tuple(Float64.(nodes[n]["coordinates"])) for n in corner_ids)
    axis=Int(record["normal_axis"]);other=Tuple(i for i in 1:3 if i!=axis);axes=(other...,axis)
    winding=sign((corners[2][other[1]]-corners[1][other[1]])*(corners[3][other[2]]-corners[1][other[2]])-
                 (corners[2][other[2]]-corners[1][other[2]])*(corners[3][other[1]]-corners[1][other[1]]))
    plane=record["plane"]=="XZ" ? :ZX : Symbol(record["plane"])
    f=QT4S.fixture(record["name"];plane,shape=record["geometry_shape"]=="unit_square" ? :unit : :rounded,
                  height=record["translation_height"],laterals=record["recombine_laterals"],layers=:one)
    return merge(f,(;source=record["input_geo"],axes,corners,winding=Int(winding),
                    levels=Float64.(record["intended_levels"]),height=Float64(record["translation_height"]),
                    surface=Int(record["source_tag"]),curve_tags=Tuple(Int.(record["source_curve_tags"])),
                    point_tags=Tuple(Int(nodes[n]["owner"][2]) for n in corner_ids),
                    intervals=Int(record["intervals"]),translation=ntuple(k->k==axis ? Float64(record["translation_height"]) : 0.,3)))
end


function four_strip_saved_case(record;oracle_only=false)
    tolerance=2e-11;n=Int(record["intervals"]);laterals=record["recombine_laterals"]
    f=four_strip_saved_fixture(record)
    saved_nodes=quad_patch_saved_integrity(record,tolerance)
    strip_saved_reference_integrity(record)
    owners=[count(node->node["owner"][1]==dim,values(saved_nodes)) for dim in 0:3]
    owners==[8,8n+28,32n-2,14n-7+(laterals ? 32 : 0)] || error("four-strip saved owner dimension counts")
    source,_=two_tri_saved_mesh(record,"p1",2,record["source_tag"])
    linear,index=two_tri_saved_mesh(record,"p1",3,record["volume_tag"])
    quadratic,_=two_tri_saved_mesh(record,"p2",3,record["volume_tag"])
    saved=QT4S.certify(linear,f,source;oracle=true)
    two_tri_saved_boundary(record,index)==Set(keys(saved.boundary)) || error("four-strip saved actual typed exterior")
    exact_saved_volume=strip_integrated_volume(linear)
    exact_saved_volume==saved.total && exact_saved_volume==
        sum(strip_fraction(row["exact_partition_volume"]) for row in record["strip_macro_certificate"]["macro_cells"]) ||
        error("four-strip saved actual reference integration differs from macro partition")
    QTNN.certify_quadratic(quadratic,linear)
    nnodes(quadratic)==54n+27+(laterals ? 32 : 0) || error("four-strip saved P2 node count")
    empty_parameters=count(node->node["owner"][1]==2 && isempty(node["stored_parameters"]),values(saved_nodes))
    empty_parameters==(laterals ? 10n : 0)==record["empty_stored_surface_uv_nodes"] ||
        error("four-strip saved stored UV provenance count")
    for (lateral,width) in zip(record["lateral_tags"],record["lateral_chain_widths"])
        entity=only(e for e in record["p2"]["entities"] if e["dim"]==2 && e["tag"]==lateral)
        length(entity["owned"]["tags"])==(width==5 ? 14n-7 : 2n-1) &&
            length(entity["closure"]["tags"])==(width==5 ? 18n+9 : 6n+3) ||
            error("four-strip saved lateral owner/closure counts")
    end
    if !oracle_only
        execution=QT4S.execute(f.source)
        volume=geo_entity_mesh(execution,3,record["volume_tag"])
        actual_source=geo_entity_mesh(execution,2,record["source_tag"])
        quad_patch_source_match(actual_source,source,tolerance)
        native=QT4S.certify(volume,f,actual_source)
        strip_integrated_volume(volume)==native.total || error("four-strip native actual reference integration")
        QTNN.surface_faces(execution.mesh_parts,volume)==Set(keys(native.boundary)) || error("four-strip native finalized typed exterior")
        nnodes(volume)==nnodes(linear) || error("four-strip native primary node count")
        # Free finalized masks can choose different valid existing-node factories.
        # Every primary column and actual terminal centroid must still agree.
        for p in eachcol(volume.coords)
            count(q->maximum(abs.(p.-q))<=tolerance,eachcol(linear.coords))==1 || error("four-strip actual primary/centroid geometry")
        end
        mktempdir() do directory
            path=joinpath(directory,"four_quad_strip_saved.geo");write(path,f.source)
            native_api(path,api->begin
                cache=api.mesh.generate(3)
                first=QT4S.public_volume(api,record["volume_tag"])
                primary_data=api.mesh.get_elements(3,record["volume_tag"])
                QTNN.typed_signature(first)==QTNN.typed_signature(volume) || error("four-strip API/GEO actual primary cells")
                isempty(api.mesh.get_elements(2)[1]) || error("four-strip API3 lower-cell contract changed")
                projected=model_to_mixed(api.CURRENT[],cache,3,record["volume_tag"])
                api.mesh.set_order(2)
                actual=QT4S.public_volume(api,record["volume_tag"])
                QTNN.certify_quadratic(actual,first)
                quadratic_volume_from_jacobians(actual,api.mesh;expected=Float64(native.total))
                strip_native_supports(api,cache,primary_data,projected,record["volume_tag"],tolerance)==
                    44n+17+(laterals ? 28 : 0) || error("four-strip native P2 support count")
                native_owners=zeros(Int,4)
                for node in api.mesh.get_nodes()[1]
                    p,uv,dim,owner=api.mesh.get_node(node);native_owners[dim+1]+=1
                    if dim==2
                        length(uv)==2 && maximum(abs.(api.model.get_value(2,owner,uv).-p))<=tolerance ||
                            error("four-strip native computed surface UV")
                    end
                end
                native_owners==owners || error("four-strip native owner dimension counts")
                for (lateral,width) in zip(record["lateral_tags"],record["lateral_chain_widths"])
                    own,_,uv=api.mesh.get_nodes(2,lateral,false,true)
                    length(own)==(width==5 ? 14n-7 : 2n-1) && length(uv)==2length(own) ||
                        error("four-strip native lateral owned query")
                    closure,coords,uv=api.mesh.get_nodes(2,lateral,true,true)
                    length(closure)==(width==5 ? 18n+9 : 6n+3) && length(uv)==2length(closure) &&
                        maximum(abs.(api.model.get_value(2,lateral,uv).-coords))<=tolerance ||
                        error("four-strip native lateral closure/computed UV")
                end
                for cap in (record["source_tag"],record["top_tag"])
                    own,_,uv=api.mesh.get_nodes(2,cap,false,true)
                    length(own)==7 && length(uv)==14 || error("four-strip native cap owned query")
                    closure,coords,uv=api.mesh.get_nodes(2,cap,true,true)
                    length(closure)==27 && length(uv)==54 &&
                        maximum(abs.(api.model.get_value(2,cap,uv).-coords))<=tolerance ||
                        error("four-strip native cap closure/computed UV")
                end
                api.mesh.set_order(1)
                QTNN.typed_signature(QT4S.public_volume(api,record["volume_tag"]))==QTNN.typed_signature(first) ||
                    error("four-strip actual P2/P1 roundtrip")
            end)
        end
    end
    return empty_parameters
end


# Definitions to append to the existing saved-case driver. No prior gate or
# pointer-sensitive family signature is changed.
function rect_grid_saved_fixture(record)
    audit=record["audit"];raw=record["payload"]
    f=QTRG.fixture(raw["name"];grid=Tuple(Int.(audit["grid"])),
        layers=audit["intervals"]==1 ? :one : :three,
        height=Float64(audit["direction"]),laterals=Bool(audit["recombine_laterals"]))
    if haskey(raw,"variant")
        variant=raw["variant"];axes=Tuple(Int.(variant["axes"]).+1)
        corners=Tuple(ntuple(k->k==axes[1] ? Float64(xy[1]) : k==axes[2] ? Float64(xy[2]) : 0.,3)
            for xy in variant["corners"])
        return merge(f,(;source=raw["input_geo"],axes,corners,
            height=Float64(variant["height"]),levels=Float64.(raw["requested_levels"])))
    end
    return merge(f,(;source=raw["input_geo"]))
end

function rect_grid_saved_integrity(record)
    raw=record["payload"];audit=record["audit"];tolerance=2e-11
    raw["status"]=="certified" && isempty(raw["warnings_or_errors"]) || error("rect-grid primary certification status")
    literal_hash=bytes2hex(sha256(raw["input_geo"]))
    if haskey(audit,"input_literal_sha256")
        literal_hash==audit["input_literal_sha256"] || error("rect-grid exact input text hash")
        raw["input_sha256"]==audit["capture_input_sha256"] || error("rect-grid preserved input-file hash")
        record["capture_sha256"]==audit["payload_sha256"] || error("rect-grid saved payload provenance")
    else
        # The new geometric captures record the actual file hash directly.
        # Keep Windows CRLF bytes distinct from the decoded recipe text hash.
        raw["input_sha256"] in (literal_hash,bytes2hex(sha256(replace(raw["input_geo"],"\n"=>"\r\n")))) || error("rect-grid geometric input-file hash")
        length(record["capture_sha256"])==64 || error("rect-grid geometric saved provenance hash")
    end
    for phase in ("p1","p2"),node in raw[phase]["nodes"]
        for (values,bits) in ((node["coordinates"],node["coordinate_bits"]),(node["stored_parameters"],node["stored_parameter_bits"]))
            length(values)==length(bits) || error("rect-grid packed payload width")
            all(string(reinterpret(UInt64,Float64(value));base=16,pad=16)==bit for (value,bit) in zip(values,bits)) || error("rect-grid saved Float64 payload bits")
        end
    end
    original=Dict(Int(node["tag"])=>node for node in raw["p1"]["nodes"])
    quadratic=Dict(Int(node["tag"])=>node for node in raw["p2"]["nodes"])
    remap=Dict(parse(Int,node)=>Int(mapped) for (node,mapped) in raw["primary_remap"])
    Set(keys(remap))==Set(keys(original)) && length(Set(values(remap)))==length(original) || error("rect-grid primary label bijection")
    for (old,new) in remap
        original[old]["coordinate_bits"]==quadratic[new]["coordinate_bits"] && original[old]["owner"]==quadratic[new]["owner"] || error("rect-grid saved primary identity changed")
    end
    inverse=Dict(new=>old for (old,new) in remap)
    supports=Dict(parse(Int,node)=>Tuple(Int.(support)) for (node,support) in raw["supports"])
    Set(keys(supports))==setdiff(Set(keys(quadratic)),Set(values(remap))) || error("rect-grid saved interpolation support coverage")
    length(Set(values(supports)))==length(supports) || error("rect-grid saved distinct support identity")
    carriers=quad_patch_actual_carriers((Int(cell["type"]),Int.(cell["nodes"]),(Int(entity["dim"]),Int(entity["tag"])))
        for entity in raw["p1"]["entities"] for cell in entity["cells"])
    for (node,support) in supports
        carrier=carriers[Tuple(sort!([inverse[id] for id in support]))]
        Tuple(Int.(quadratic[node]["owner"]))==carrier || error("rect-grid saved interpolation carrier")
        mean=[sum(Float64(quadratic[id]["coordinates"][axis]) for id in support)/length(support) for axis in 1:3]
        maximum(abs.(Float64.(quadratic[node]["coordinates"]).-mean))<=tolerance || error("rect-grid saved interpolation placement")
    end
    for phase in ("p1","p2")
        certificates=raw[phase*"_whole_maps"]
        cells=Set((Int(entity["dim"]),Int(entity["tag"]),Int(cell["tag"]),Int(cell["type"])) for entity in raw[phase]["entities"] if entity["dim"]>0 for cell in entity["cells"])
        Set((Int(row["entity"][1]),Int(row["entity"][2]),Int(row["cell"]),Int(row["type"])) for row in certificates)==cells || error("rect-grid saved whole-map coverage")
        all(strip_fraction(row["minimum_exact"])>0 && Float64(strip_fraction(row["minimum_exact"]))==row["minimum"] for row in certificates) || error("rect-grid saved positive exact whole-map bounds")
    end
    recorded=Dict((Tuple(Int.(row["entity"])),Int(row["cell"]))=>row for row in raw["p2_whole_maps"])
    for entity in raw["p2"]["entities"]
        entity["dim"]==3 || continue
        for cell in entity["cells"]
            map=QTRG.quadratic_map_certificate(Int(cell["type"]),Tuple(Tuple(Float64.(quadratic[Int(node)]["coordinates"])) for node in cell["nodes"]))
            stored=recorded[((3,Int(entity["tag"])),Int(cell["tag"]))]
            map.minimum==strip_fraction(stored["minimum_exact"]) && map.coefficients==stored["coefficients"] && collect(map.degrees)==stored["degrees"] || error("rect-grid actual P2 polynomial differs from independent exact primary proof")
        end
    end
    for query in raw["computed_uv"]
        length(query["computed_uv"])==2 && all(isfinite,query["computed_uv"]) || error("rect-grid saved computed surface UV width")
        # The original six captures retain inverse UVs without a getValue
        # result. Later captures additionally preserve their reevaluation.
        # Native reevaluation is checked independently for every case below.
        if haskey(query,"evaluated")
            maximum(abs.(Float64.(query["evaluated"]).-Float64.(quadratic[Int(query["node"])]["coordinates"])))<=tolerance || error("rect-grid saved computed surface UV")
        end
    end
    # Compare full face-node identities through actual supports, not merely
    # the common primary coordinates or a quadrature sample.
    basis=Dict{Int,Tuple}(node=>(node,) for node in values(remap));merge!(basis,supports)
    face_rows=Dict{Tuple,Vector{Tuple}}()
    for entity in raw["p2"]["entities"]
        entity["dim"]==3 || continue
        for cell in entity["cells"],pattern in QTRG.FACES[Int(cell["type"])-7]
            nodes=Int.(cell["nodes"]);primary=QTRG.key(nodes[k] for k in pattern)
            actual=QTRG.key(node for node in nodes if Set(basis[node])⊆Set(primary))
            length(actual)==(length(primary)==3 ? 6 : 9) || error("rect-grid saved quadratic face width")
            push!(get!(face_rows,primary,Tuple[]),actual)
        end
    end
    all(length(rows) in (1,2) && (length(rows)==1 || rows[1]==rows[2]) for rows in values(face_rows)) || error("rect-grid saved complete P2 face conformity")
    boundary=Dict(face=>only(rows) for (face,rows) in face_rows if length(rows)==1)
    lower=Dict(QTRG.key(cell["nodes"][1:(Int(cell["type"])==9 ? 3 : 4)])=>QTRG.key(cell["nodes"])
        for entity in raw["p2"]["entities"] if entity["dim"]==2 for cell in entity["cells"])
    boundary==lower || error("rect-grid saved actual P2 typed exterior")
    return quadratic
end

function rect_grid_saved_surface_queries(record)
    haskey(record["audit"],"surface_queries") && return record["audit"]["surface_queries"]
    raw=record["payload"]
    return [Dict("tag"=>Int(entity["tag"]),"owned"=>length(entity["owned"]["tags"]),
        "closure"=>length(entity["closure"]["tags"]),"empty_stored_uv"=>count(node->node["owner"]==[2,entity["tag"]] && isempty(node["stored_parameters"]),raw["p2"]["nodes"]))
        for entity in raw["p2"]["entities"] if entity["dim"]==2]
end

function rect_grid_check_factories(certificate,relation)
    for (macro_id,state) in certificate.states
        families=Tuple(sort!([cell.msh for cell in certificate.domains[macro_id]]))
        (state,families) in relation || error("rect-grid actual macro differs from immutable315 state/family relation")
    end
end

function rect_grid_saved_case(record,masks,relation;oracle_only=false)
    raw=record["payload"];audit=record["audit"];f=rect_grid_saved_fixture(record);n=f.intervals;a,b=f.grid
    saved_nodes=rect_grid_saved_integrity(record)
    source,_=two_tri_saved_mesh(raw,"p1",2,1);linear,index=two_tri_saved_mesh(raw,"p1",3,1)
    saved=QTRG.certify(linear,f,source;oracle=true)
    all(state in masks for state in values(saved.states)) || error("rect-grid primary factory relation")
    rect_grid_check_factories(saved,relation)
    two_tri_saved_boundary(raw,index)==Set(keys(saved.boundary)) || error("rect-grid primary typed exterior")
    isempty(raw["logical"]["extra_nodes"]) && nnodes(linear)==(a+1)*(b+1)*(n+1) || error("rect-grid primary no-center product")
    nnodes(two_tri_saved_mesh(raw,"p2",3,1)[1])==(2a+1)*(2b+1)*(2n+1) || error("rect-grid primary P2 lattice count")
    owners=[count(node->node["owner"][1]==dim,values(saved_nodes)) for dim in 0:3]
    empty_parameters=count(node->node["owner"][1]==2 && isempty(node["stored_parameters"]),values(saved_nodes))
    surface_queries=rect_grid_saved_surface_queries(record)
    empty_parameters==sum(Int(row["empty_stored_uv"]) for row in surface_queries) || error("rect-grid stored-UV gap provenance")
    if !oracle_only
        execution=QTRG.execute(f.source);volume=geo_entity_mesh(execution,3,1);actual_source=geo_entity_mesh(execution,2,1)
        quad_patch_source_match(actual_source,source,2e-11)
        actual=QTRG.certify(volume,f,actual_source)
        all(state in masks for state in values(actual.states)) || error("rect-grid native factory relation")
        rect_grid_check_factories(actual,relation)
        QTNN.surface_faces(execution.mesh_parts,volume)==Set(keys(actual.boundary)) || error("rect-grid native complete typed exterior")
        projection=model_to_mixed(execution.model,volume,3,1)
        nodes=Dict(Int(node["tag"])=>node for node in raw["p1"]["nodes"])
        for node in 1:nnodes(volume)
            matches=[tag for tag in keys(index) if maximum(abs.(QTRG.point(volume,node).-Tuple(Float64.(nodes[tag]["coordinates"]))))<=2e-11]
            length(matches)==1 || error("rect-grid native actual primary geometry bijection")
            Tuple(Int.(nodes[only(matches)]["owner"]))==projection.entity_data.node_entities[node] || error("rect-grid native primary owner")
        end
        chains=haskey(raw,"source_chains") ? raw["source_chains"] : audit["source_chains"]
        for chain in chains
            params=execution.model.curve_params[Int(chain["curve"])]
            expected=Float64.(chain["stored_owned_parameters"])
            length(params)==length(expected)+2 && maximum(abs.(sort(params[2:end-1]).-sort(expected)))<=2e-11 || error("rect-grid actual sampled curve parameters")
        end
        mktempdir() do directory
            path=joinpath(directory,"rect_grid_saved.geo");write(path,f.source)
            native_api(path,api->begin
                cache=api.mesh.generate(3);first=QTRG.public_volume(api,1)
                QTNN.typed_signature(first)==QTNN.typed_signature(volume) || error("rect-grid API/GEO primary cells")
                data=api.mesh.get_elements(3,1);projected=model_to_mixed(api.CURRENT[],cache,3,1)
                isempty(api.mesh.get_elements(2)[1]) || error("rect-grid API3 lower-cell contract")
                api.mesh.set_order(2);p2=QTRG.public_volume(api,1)
                certificate=QTRG.certify_quadratic(p2,f,actual_source)
                abs(Float64(certificate.total)-Float64(actual.total))<=2e-11 || error("rect-grid actual P2 reference integral")
                strip_native_supports(api,cache,data,projected,1,2e-11)==nnodes(p2)-nnodes(first) || error("rect-grid complete actual P2 supports")
                native_owners=zeros(Int,4)
                for node in api.mesh.get_nodes()[1]
                    p,uv,dim,owner=api.mesh.get_node(node);native_owners[dim+1]+=1
                    if dim in (1,2)
                        length(uv)==dim && maximum(abs.(api.model.get_value(dim,owner,uv).-p))<=2e-11 || error("rect-grid native computed carrier parameters")
                    end
                end
                native_owners==owners || error("rect-grid native P2 owner counts")
                for row in surface_queries,include_boundary in (false,true)
                    ids,xyz,uv=api.mesh.get_nodes(2,Int(row["tag"]),include_boundary,true)
                    length(ids)==row[include_boundary ? "closure" : "owned"] && length(uv)==2length(ids) || error("rect-grid native surface owned/closure query")
                    maximum(abs.(api.model.get_value(2,Int(row["tag"]),uv).-xyz))<=2e-11 || error("rect-grid native boundary computed UV")
                end
                api.mesh.set_order(1)
                QTNN.typed_signature(QTRG.public_volume(api,1))==QTNN.typed_signature(first) || error("rect-grid P2/P1 actual cell roundtrip")
            end)
        end
    end
    return empty_parameters
end


initialize_oracle()
oracle_only=get(ENV,"QUADTRI_NONEW_ORACLE_ONLY","")=="1"
completed=Ref(0);oracle_samples=Ref(0);quadratic_cases=Ref(0)
height_errors=Ref(0);isolated_cases=Ref(0)
native_helical_cases=Ref(0);helical_oracle_gaps=Ref(0)
triangle_cases=Ref(0);triangle_quadratic_cases=Ref(0);triangle_isolated_cases=Ref(0)
triangle_parameter_provenance_gaps=Ref(0)
two_tri_cases=Ref(0);two_tri_quadratic_cases=Ref(0);two_tri_parameter_provenance_gaps=Ref(0)
two_tri_saved_empty_parameters=Ref(0)
quad_patch_cases=Ref(0);quad_patch_quadratic_cases=Ref(0);quad_patch_trapezoid_cases=Ref(0)
quad_patch_parameter_provenance_gaps=Ref(0);quad_patch_empty_parameters=Ref(0)
quad_strip_cases=Ref(0);quad_strip_quadratic_cases=Ref(0);quad_strip_variant_cases=Ref(0)
quad_strip_parameter_provenance_gaps=Ref(0);quad_strip_empty_parameters=Ref(0)
three_strip_cases=Ref(0);three_strip_quadratic_cases=Ref(0);three_strip_variant_cases=Ref(0)
three_strip_parameter_provenance_gaps=Ref(0);three_strip_empty_parameters=Ref(0)
four_strip_cases=Ref(0);four_strip_quadratic_cases=Ref(0);four_strip_variant_cases=Ref(0)
four_strip_parameter_provenance_gaps=Ref(0);four_strip_empty_parameters=Ref(0)
rect_grid_cases=Ref(0);rect_grid_quadratic_cases=Ref(0);rect_grid_variant_cases=Ref(0)
rect_grid_parameter_provenance_gaps=Ref(0);rect_grid_empty_parameters=Ref(0)
try
    startswith(gmsh.option.getString("General.Version"),"4.15.2") ||
        error("Gmsh4.15.2 runtime required")
    selected=get(ENV,"QUADTRI_NONEW_CASE","")
    for fixture in QTNN.fixtures()
        isempty(selected) || occursin(selected,fixture.name) || continue
        mktempdir() do directory
            path=joinpath(directory,"$(fixture.name).geo")
            write(path,fixture.source*"Mesh 3;\n")
            # Repeats expose allocator choices; validity does not depend on
            # seeing every possible outcome during this bounded gate.
            repeat_count=startswith(fixture.name,"uniform") ? 3 : 1
            native=oracle_only ? nothing : QTNN.execute(fixture.source)
            native_execution[]=native
            nvolume=oracle_only ? nothing : geo_entity_mesh(native,3,1)
            outcomes=Set{Tuple}()
            for repeat in 1:repeat_count
                gmsh.clear();gmsh.open(path)
                oracle,element_tags=oracle_volume()
                oc=QTNN.certify_fixture(oracle,fixture)
                oracle_surface_faces(oracle)==Set(keys(oc.boundary)) ||
                    error("$(fixture.name): oracle surface boundary mismatch")
                qualities=gmsh.model.mesh.getElementQualities(element_tags,"minDetJac")
                all(q->isfinite(q)&&q>0,qualities) || error("$(fixture.name): invalid oracle Jacobian")
                oracle_only || same_grid_and_centroid(nvolume,oracle,fixture)
                push!(outcomes,Tuple(sort!(collect(oc.family_counts))))
                oracle_samples[]+=1
            end
            if !oracle_only
                validate(native.mesh).ok || error("$(fixture.name): invalid native global mesh")
                again=QTNN.execute(fixture.source)
                mixed_crc(again.mesh)==mixed_crc(native.mesh) ||
                    error("$(fixture.name): nondeterministic native CRC")
                native_api(path,api->begin
                    actual=api.mesh.generate(3)
                    QTNN.certify_fixture(actual,fixture)
                    QTNN.typed_signature(actual)==QTNN.typed_signature(nvolume) ||
                        error("$(fixture.name): API/GEO volume differs")
                    isempty(api.mesh.get_elements(2)[1]) ||
                        error("legacy API3 cache lower-cell contract changed")
                    length(api.mesh.get_nodes(3,1,true)[1])==nnodes(actual) ||
                        error("$(fixture.name): API volume closure nodes differ")
                    sum(length,api.mesh.get_elements(3,1)[2])==
                        QTNN.certify_complex(actual).ncell || error("API3 cell query differs")
                end)
            end
            completed[]+=1
            println("QUADTRI_NONEW_CASE_OK name=$(fixture.name) oracle_outcomes=$(collect(outcomes))")
        end
    end
    if isempty(selected)
        for fixture in filter(f->startswith(f.name,"uniform"),QTNN.fixtures())
            mktempdir() do directory
                path=joinpath(directory,"quadratic.geo");write(path,fixture.source)
                gmsh.clear();gmsh.open(path);gmsh.model.mesh.generate(3)
                linear,_=oracle_volume()
                gmsh.model.mesh.setOrder(2)
                elevated,tags=oracle_volume()
                QTNN.certify_quadratic(elevated,linear)
                quadratic_volume_from_jacobians(elevated,gmsh.model.mesh)
                all(>(0),gmsh.model.mesh.getElementQualities(tags,"minDetJac")) ||
                    error("oracle P2 quality")
                if !oracle_only
                    native_api(path,api->begin
                        first=api.mesh.generate(3)
                        api.mesh.set_order(2);quadratic=api.mesh.get()
                        QTNN.certify_quadratic(quadratic,first)
                        quadratic_volume_from_jacobians(quadratic,api.mesh)
                        api.mesh.set_order(1)
                        QTNN.typed_signature(api.mesh.get())==QTNN.typed_signature(first) ||
                            error("P2/P1 roundtrip cell set changed")
                    end)
                end
                quadratic_cases[]+=1;oracle_samples[]+=1
                println("QUADTRI_NONEW_QUADRATIC_OK name=$(fixture.name)")
            end
        end
        for coincident in (false,true)
            mktempdir() do directory
                path=joinpath(directory,"isolated.geo");write(path,QTNN.paired_source())
                gmsh.clear();gmsh.open(path)
                if coincident
                    for (_,tag) in gmsh.model.getEntities(0)
                        x,y,z=gmsh.model.getValue(0,tag,Float64[])
                        x>=3 || continue
                        gmsh.model.setCoordinates(tag,x-3,y,z)
                    end
                end
                gmsh.model.mesh.generate(3)
                entities=last.(gmsh.model.getEntities(3))
                length(entities)==2 || error("isolated oracle region count")
                volumes=[oracle_volume(tag)[1] for tag in entities]
                fixtures=(QTNN.fixture("uniform3_false","Layers{3}",[0.,1/3,2/3,1.],false),
                          QTNN.fixture("uniform3_true","Layers{3}",[0.,1/3,2/3,1.],true))
                for (i,volume) in enumerate(volumes)
                    shifted=i==2 && !coincident ? MixedMesh(volume.coords.-[3.,0,0],volume.blocks) : volume
                    QTNN.certify_fixture(shifted,fixtures[i])
                end
                tags=gmsh.model.mesh.getNodes()[1]
                length(tags)==33 || error("isolated oracle nodes $(length(tags)) differ from 33")
                a=gmsh.model.mesh.getNodes(3,entities[1],true)[1]
                b=gmsh.model.mesh.getNodes(3,entities[2],true)[1]
                isempty(intersect(a,b)) || error("independent oracle topology welded")
                if !oracle_only
                    native_api(path,api->begin
                        if coincident
                            for (tag,p) in collect(api.CURRENT[].points)
                                p[1]>=3 || continue
                                api.model.set_coordinates(tag,p[1]-3,p[2],p[3])
                            end
                        end
                        mesh=api.mesh.generate(3)
                        nnodes(mesh)==33 || error("isolated native nodes differ")
                        a=api.mesh.get_nodes(3,entities[1],true)[1]
                        b=api.mesh.get_nodes(3,entities[2],true)[1]
                        length(a)==16 && length(b)==17 && isempty(intersect(a,b)) ||
                            error("isolated native ownership differs")
                        validate(mesh).ok || error("isolated native mesh invalid")
                    end)
                end
                isolated_cases[]+=1;oracle_samples[]+=1
                println("QUADTRI_NONEW_ISOLATED_OK coincident=$coincident")
            end
        end
        # Large-angle Gmsh extrusion takes a different route from this bounded
        # native planner. Record actual outcomes without claiming cell parity.
        for (pitch,laterals) in ((-3.,false),(-3.,true),(3.,false),(3.,true),
                                 (-.1,true),(.1,true))
            fixture=QTNN.helical_fixture("helical_oracle_gap",pitch,laterals)
            mktempdir() do directory
                path=joinpath(directory,"helix.geo");write(path,fixture.source)
                # After a fatal FindDiagonalEdgeIndices error, clear() alone
                # leaves upstream parser state unusable for later fixtures.
                gmsh.finalize();initialize_oracle();gmsh.open(path)
                failure=try gmsh.model.mesh.generate(3);nothing catch e;e end
                source_types,source_tags,_=gmsh.model.mesh.getElements(2,1)
                source_types==[3] && length(only(source_tags))==1 ||
                    error("helical oracle source quad changed")
                types,tags,_=gmsh.model.mesh.getElements(3,1)
                if !laterals
                    failure!==nothing && occursin("FindDiagonalEdgeIndices",sprint(showerror,failure)) ||
                        error("unexpected free helical oracle outcome")
                    isempty(types) || error("failed free helical oracle emitted volume cells")
                    println("QUADTRI_NONEW_HELIX_ORACLE_GAP pitch=$pitch laterals=$laterals kind=generation_error volume_cells=0")
                else
                    failure===nothing || error("unexpected recombined helical oracle error: $failure")
                    oracle,element_tags=oracle_volume()
                    all(>(0),gmsh.model.mesh.getElementQualities(element_tags,"minDetJac")) ||
                        error("large-angle oracle local quality changed")
                    # A positive sampled quality alone is not a complete cell,
                    # face-complex or global disjointness certificate.
                    local_result=try
                        QTNN.certify_complex(oracle);"certified"
                    catch e
                        sprint(showerror,e)
                    end
                    parity_failure=try QTNN.certify_fixture(oracle,fixture);nothing catch e;e end
                    parity_failure!==nothing || error("large-angle oracle gap disappeared; reassess parity")
                    families=Tuple((Int(type),length(family)) for (type,family) in zip(types,tags))
                    println("QUADTRI_NONEW_HELIX_ORACLE_GAP pitch=$pitch laterals=$laterals kind=different_cell_route volume_nodes=$(nnodes(oracle)) families=$families local_complex=$(repr(local_result))")
                end
                if !oracle_only
                    if abs(pitch)==3
                        native=QTNN.execute(fixture.source)
                        nvolume=geo_entity_mesh(native,3,1)
                        QTNN.certify_fixture(nvolume,fixture)
                        native_api(path,api->begin
                            actual=api.mesh.generate(3)
                            QTNN.certify_fixture(actual,fixture)
                            QTNN.typed_signature(actual)==QTNN.typed_signature(nvolume) ||
                                error("certified native helical API/GEO cells differ")
                        end)
                        native_helical_cases[]+=1
                    else
                        failure=try QTNN.execute(fixture.source);nothing catch e;e end
                        failure isa ArgumentError && occursin("global P1 cell-hull disjointness is not certified",sprint(showerror,failure)) ||
                            error("native overlapping helix was not rejected precisely")
                    end
                end
                helical_oracle_gaps[]+=1;oracle_samples[]+=1
            end
        end
        # Literal layer heights do not move the CAD cap, which stays at u=1.
        # These upstream errors are recorded separately from valid fixtures.
        for height in (.5,2.),laterals in (false,true)
            fixture=QTNN.fixture("terminal_height","Layers{{1},{$height}}",[0.,height],laterals)
            mktempdir() do directory
                path=joinpath(directory,"invalid_height.geo");write(path,fixture.source)
                gmsh.finalize();initialize_oracle();gmsh.open(path)
                failure=try gmsh.model.mesh.generate(3);nothing catch e;e end
                failure!==nothing && occursin("Could not find extruded node",sprint(showerror,failure)) ||
                    error("unexpected non-unit height oracle behavior")
                isempty(gmsh.model.mesh.getElements(3,1)[1]) || error("invalid height created volume cells")
                if !oracle_only
                    failure=try QTNN.execute(fixture.source);nothing catch e;e end
                    failure isa ArgumentError && occursin("normalized final layer height of 1.0",sprint(showerror,failure)) ||
                        error("native non-unit height was not rejected precisely")
                end
                height_errors[]+=1;oracle_samples[]+=1
            end
        end
    end
    for fixture in QTNT.fixtures()
        isempty(selected) || occursin(selected,fixture.name) || continue
        mktempdir() do directory
            path=joinpath(directory,"$(fixture.name).geo")
            write(path,fixture.source*"Mesh 3;\n")
            gmsh.clear();gmsh.open(path)
            oracle,element_tags=oracle_volume()
            oc=QTNT.certify(oracle,fixture)
            source_types,_,source_nodes=gmsh.model.mesh.getElements(2,1)
            source_types==[2] && length(only(source_nodes))==3 || error("triangle oracle source cell changed")
            source_points=Tuple(Tuple(gmsh.model.mesh.getNode(node)[1]) for node in only(source_nodes))
            q=map(QTNT.exact_point,source_points)
            sign(QTNT.determinant(q[2].-q[1],q[3].-q[1],(0,0,1)))==fixture.winding ||
                error("triangle oracle source winding changed")
            oracle_surface_faces(oracle)==Set(keys(oc.boundary)) ||
                error("$(fixture.name): triangle oracle surface boundary mismatch")
            all(>(0),gmsh.model.mesh.getElementQualities(element_tags,"minDetJac")) ||
                error("$(fixture.name): triangle oracle Jacobian quality")
            if !oracle_only
                execution=QTNN.execute(fixture.source)
                native=geo_entity_mesh(execution,3,1);nc=QTNT.certify(native,fixture)
                QTNT.certify_cap_winding(execution.mesh_parts,native,fixture,nc)
                QTNN.surface_faces(execution.mesh_parts,native)==Set(keys(nc.boundary)) ||
                    error("$(fixture.name): native triangle surface boundary mismatch")
                nnodes(native)==nnodes(oracle) || error("triangle volume node counts differ")
                for i in eachindex(nc.grid_ids)
                    p=QTNN.point(native,nc.grid_ids[i]);q=QTNN.point(oracle,oc.grid_ids[i])
                    maximum(abs.(p.-q))<=2e-11 || error("triangle swept coordinates differ")
                end
                nc.family_counts==oc.family_counts || error("triangle volume family counts differ")
                QTNT.crc(QTNN.execute(fixture.source).mesh)==QTNT.crc(execution.mesh) ||
                    error("nondeterministic native triangle CRC")
                native_api(path,api->begin
                    actual=api.mesh.generate(3)
                    QTNT.certify(actual,fixture)
                    QTNN.typed_signature(actual)==QTNN.typed_signature(native) ||
                        error("triangle API/GEO volume differs")
                    isempty(api.mesh.get_elements(2)[1]) || error("legacy API3 lower-cell contract changed")
                    length(api.mesh.get_nodes(3,1,true)[1])==nnodes(actual) ||
                        error("triangle API volume closure differs")
                    isempty(api.mesh.get_nodes(3,1,false)[1]) || error("triangle API invented body nodes")
                end)
            end
            triangle_cases[]+=1;oracle_samples[]+=1
            println("QUADTRI_NONEW_TRIANGLE_OK name=$(fixture.name) families=$(oc.family_counts)")
        end
    end
    if isempty(selected) || selected=="triangle"
        for laterals in (false,true),winding in (-1,1)
            fixture=QTNT.fixture("triangle_p2","Layers{3}",[0.,1/3,2/3,1.],laterals;winding)
            mktempdir() do directory
                path=joinpath(directory,"triangle_p2.geo");write(path,fixture.source)
                gmsh.clear();gmsh.open(path);gmsh.model.mesh.generate(3)
                linear,_=oracle_volume();QTNT.certify(linear,fixture)
                gmsh.model.mesh.setOrder(2);quadratic,tags=oracle_volume()
                QTNN.certify_quadratic(quadratic,linear)
                quadratic_volume_from_jacobians(quadratic,gmsh.model.mesh;expected=.5)
                all(>(0),gmsh.model.mesh.getElementQualities(tags,"minDetJac")) ||
                    error("triangle oracle P2 quality")
                empty_face_parameters=0
                oracle_lateral_count=0
                for (_,surface) in gmsh.model.getEntities(2)
                    owned,_,parameters=gmsh.model.mesh.getNodes(2,surface,false,true)
                    isempty(owned) && continue
                    oracle_lateral_count+=1
                    length(owned)==5 && length(parameters)==(laterals ? 4 : 10) ||
                        error("oracle lateral ownership/parameter dimensions changed")
                    boundary_nodes,_,boundary_parameters=gmsh.model.mesh.getNodes(2,surface,true,true)
                    length(boundary_nodes)==21 && length(boundary_parameters)==(laterals ? 36 : 42) ||
                        error("oracle partial boundary parameter provenance changed")
                    for node in owned
                        p,uv,dim,owner=gmsh.model.mesh.getNode(node)
                        dim==2 && owner==surface || error("oracle P2 lateral ownership")
                        if isempty(uv)
                            empty_face_parameters+=1
                            # The stored absence is distinct from the complete
                            # geometric inverse available on the same surface.
                            computed=gmsh.model.getParametrization(2,surface,p)
                            maximum(abs.(gmsh.model.getValue(2,surface,computed).-p))<=2e-11 ||
                                error("oracle geometric inverse roundtrip")
                        end
                    end
                end
                oracle_lateral_count==3 || error("oracle P2 lateral ownership missing")
                empty_face_parameters==(laterals ? 9 : 0) || error("oracle face-center parameter provenance")
                if laterals
                    triangle_parameter_provenance_gaps[]+=1
                    println("QUADTRI_NONEW_TRIANGLE_PARAMETER_PROVENANCE_GAP winding=$winding oracle_empty_face_nodes=$empty_face_parameters native_contract=computed_uv")
                end
                if !oracle_only
                    native_api(path,api->begin
                        api.mesh.generate(3)
                        first=triangle_api_volume(api)
                        api.mesh.set_order(2);actual=triangle_api_volume(api)
                        QTNN.certify_quadratic(actual,first)
                        quadratic_volume_from_jacobians(actual,api.mesh;expected=.5)
                        node_tags=api.mesh.get_nodes()[1]
                        for node in node_tags
                            p,uv,dim,entity=api.mesh.get_node(node)
                            dim==2 || continue
                            length(uv)==2 || error("native computed lateral UV missing")
                            maximum(abs.(api.model.get_value(2,entity,uv).-p))<=2e-11 ||
                                error("native computed lateral UV does not reproduce actual node")
                        end
                        lateral_count=0
                        for (_,surface) in api.model.get_entities(2)
                            # The source and cap have no interior P2 nodes.
                            # All three actual laterals own five supports.
                            owned,_,parameters=api.mesh.get_nodes(2,surface,false,true)
                            isempty(owned) && continue
                            lateral_count+=1
                            length(owned)==5 && length(parameters)==10 ||
                                error("native computed lateral parameter dimensions")
                            boundary_nodes,coordinates,uv=api.mesh.get_nodes(2,surface,true,true)
                            length(boundary_nodes)==21 && length(uv)==42 ||
                                error("native computed boundary parameter dimensions")
                            maximum(abs.(api.model.get_value(2,surface,uv).-coordinates))<=2e-11 ||
                                error("native computed boundary UV roundtrip")
                        end
                        lateral_count==3 || error("native P2 lateral ownership missing")
                        api.mesh.set_order(1)
                        QTNN.typed_signature(triangle_api_volume(api))==QTNN.typed_signature(first) ||
                            error("triangle P2/P1 actual cell roundtrip differs")
                    end)
                end
                triangle_quadratic_cases[]+=1;oracle_samples[]+=1
                println("QUADTRI_NONEW_TRIANGLE_P2_OK laterals=$laterals winding=$winding")
            end
        end
        for coincident in (false,true)
            mktempdir() do directory
                path=joinpath(directory,"triangle_isolated.geo");write(path,QTNT.paired_source())
                gmsh.clear();gmsh.open(path)
                if coincident
                    for (_,tag) in gmsh.model.getEntities(0)
                        x,y,z=gmsh.model.getValue(0,tag,Float64[])
                        x>=3 || continue
                        gmsh.model.setCoordinates(tag,x-3,y,z)
                    end
                end
                gmsh.model.mesh.generate(3)
                entities=last.(gmsh.model.getEntities(3));length(entities)==2 || error("triangle region count")
                for (index,entity) in enumerate(entities)
                    volume=oracle_volume(entity)[1]
                    if index==2 && !coincident
                        volume=MixedMesh(volume.coords.-[3.,0,0],volume.blocks)
                    end
                    f=QTNT.fixture("triangle_isolated","Layers{3}",[0.,1/3,2/3,1.],index==2)
                    QTNT.certify(volume,f)
                end
                a=gmsh.model.mesh.getNodes(3,entities[1],true)[1]
                b=gmsh.model.mesh.getNodes(3,entities[2],true)[1]
                length(a)==length(b)==12 && isempty(intersect(a,b)) || error("triangle oracle identities welded")
                length(gmsh.model.mesh.getNodes()[1])==24 || error("triangle oracle global node count")
                if !oracle_only
                    native_api(path,api->begin
                        if coincident
                            for (tag,p) in collect(api.CURRENT[].points)
                                p[1]>=3 || continue
                                api.model.set_coordinates(tag,p[1]-3,p[2],p[3])
                            end
                        end
                        generated=api.mesh.generate(3)
                        a=api.mesh.get_nodes(3,entities[1],true)[1]
                        b=api.mesh.get_nodes(3,entities[2],true)[1]
                        length(a)==length(b)==12 && isempty(intersect(a,b)) || error("triangle native identities welded")
                        nnodes(generated)==24 || error("triangle native global node count")
                    end)
                end
                triangle_isolated_cases[]+=1;oracle_samples[]+=1
                println("QUADTRI_NONEW_TRIANGLE_ISOLATED_OK coincident=$coincident")
            end
        end
    end
    catalog=TOML.parsefile(joinpath(@__DIR__,"..","..","test","artifacts","quadtri_nonew_two_tri_oracle.toml"))
    catalog["gmsh_version"]=="4.15.2" && catalog["coordinate_abs_tolerance"]==2e-11 &&
        length(catalog["fixtures"])==12 || error("two-Tri saved oracle catalog contract")
    for record in catalog["fixtures"]
        isempty(selected) || occursin(selected,record["name"]) || continue
        empty_parameters=two_tri_saved_case(record;oracle_only)
        two_tri_cases[]+=1;two_tri_quadratic_cases[]+=1
        two_tri_saved_empty_parameters[]+=empty_parameters
        empty_parameters>0 && (two_tri_parameter_provenance_gaps[]+=1)
        println("QUADTRI_NONEW_TWO_TRI_SAVED_OK name=$(record["name"]) p1_nodes=$(length(record["p1"]["nodes"])) p2_nodes=$(length(record["p2"]["nodes"])) oracle_empty_stored_uv=$empty_parameters native_contract=computed_uv")
    end
    patch_catalog=TOML.parsefile(joinpath(@__DIR__,"..","..","test","artifacts","quadtri_nonew_quad_patch_oracle.toml"))
    patch_catalog["gmsh_version"]=="4.15.2" && patch_catalog["coordinate_abs_tolerance"]==2e-11 &&
        length(patch_catalog["fixtures"])==16 || error("quad-patch saved oracle catalog contract")
    for record in patch_catalog["fixtures"]
        isempty(selected) || occursin(selected,record["name"]) || continue
        empty_parameters=quad_patch_saved_case(record;oracle_only)
        quad_patch_cases[]+=1;quad_patch_quadratic_cases[]+=1
        quad_patch_empty_parameters[]+=empty_parameters
        empty_parameters>0 && (quad_patch_parameter_provenance_gaps[]+=1)
        record["geometry_shape"]=="trapezoid" && (quad_patch_trapezoid_cases[]+=1)
        println("QUADTRI_NONEW_QUAD_PATCH_SAVED_OK name=$(record["name"]) p1_nodes=$(length(record["p1"]["nodes"])) p2_nodes=$(length(record["p2"]["nodes"])) oracle_empty_stored_uv=$empty_parameters native_contract=computed_uv")
    end
    strip_catalog=TOML.parsefile(joinpath(@__DIR__,"..","..","test","artifacts","quadtri_nonew_quad_strip_oracle.toml"))
    strip_catalog["gmsh_version"]=="4.15.2" && strip_catalog["coordinate_abs_tolerance"]==2e-11 &&
        length(strip_catalog["fixtures"])==16 || error("quad-strip saved oracle catalog contract")
    for record in strip_catalog["fixtures"]
        isempty(selected) || occursin(selected,record["name"]) || continue
        empty_parameters=strip_saved_case(record;oracle_only)
        quad_strip_cases[]+=1;quad_strip_quadratic_cases[]+=1
        quad_strip_empty_parameters[]+=empty_parameters
        empty_parameters>0 && (quad_strip_parameter_provenance_gaps[]+=1)
        record["capture_kind"]=="independent_variant" && (quad_strip_variant_cases[]+=1)
        println("QUADTRI_NONEW_QUAD_STRIP_SAVED_OK name=$(record["name"]) p1_nodes=$(length(record["p1"]["nodes"])) p2_nodes=$(length(record["p2"]["nodes"])) oracle_empty_stored_uv=$empty_parameters native_contract=computed_uv")
    end
    three_strip_catalog=TOML.parsefile(joinpath(@__DIR__,"..","..","test","artifacts","quadtri_nonew_three_quad_strip_oracle.toml"))
    three_strip_catalog["gmsh_version"]=="4.15.2" && three_strip_catalog["coordinate_abs_tolerance"]==2e-11 &&
        length(three_strip_catalog["fixtures"])==16 || error("three-strip saved oracle catalog contract")
    for record in three_strip_catalog["fixtures"]
        isempty(selected) || occursin(selected,record["name"]) || continue
        empty_parameters=three_strip_saved_case(record;oracle_only)
        three_strip_cases[]+=1;three_strip_quadratic_cases[]+=1
        three_strip_empty_parameters[]+=empty_parameters
        empty_parameters>0 && (three_strip_parameter_provenance_gaps[]+=1)
        record["capture_kind"]=="independent_variant" && (three_strip_variant_cases[]+=1)
        println("QUADTRI_NONEW_THREE_QUAD_STRIP_SAVED_OK name=$(record["name"]) p1_nodes=$(length(record["p1"]["nodes"])) p2_nodes=$(length(record["p2"]["nodes"])) oracle_empty_stored_uv=$empty_parameters native_contract=computed_uv")
    end
    four_strip_catalog=TOML.parsefile(joinpath(@__DIR__,"..","..","test","artifacts","quadtri_nonew_four_quad_strip_oracle.toml"))
    four_strip_catalog["gmsh_version"]=="4.15.2" && four_strip_catalog["coordinate_abs_tolerance"]==2e-11 &&
        length(four_strip_catalog["fixtures"])==24 && four_strip_catalog["empty_stored_surface_uv_nodes"]==340 ||
        error("four-strip saved oracle catalog contract")
    for record in four_strip_catalog["fixtures"]
        isempty(selected) || occursin(selected,record["name"]) || continue
        empty_parameters=four_strip_saved_case(record;oracle_only)
        four_strip_cases[]+=1;four_strip_quadratic_cases[]+=1
        four_strip_empty_parameters[]+=empty_parameters
        empty_parameters>0 && (four_strip_parameter_provenance_gaps[]+=1)
        record["capture_kind"]=="independent_variant" && (four_strip_variant_cases[]+=1)
        println("QUADTRI_NONEW_FOUR_QUAD_STRIP_SAVED_OK name=$(record["name"]) p1_nodes=$(length(record["p1"]["nodes"])) p2_nodes=$(length(record["p2"]["nodes"])) oracle_empty_stored_uv=$empty_parameters native_contract=computed_uv")
    end
    rect_catalog=TOML.parsefile(joinpath(@__DIR__,"..","..","test","artifacts","quadtri_nonew_rect_grid_oracle.toml"))
    rect_variant_path=joinpath(@__DIR__,"..","..","test","artifacts","quadtri_nonew_rect_grid_variants_oracle.toml")
    rect_variants=TOML.parsefile(rect_variant_path)
    # This new artifact was independently checked against all eight original
    # JSON/file hashes and losslessly decoded. Bind that entire payload as
    # well as its primary manifest; canonical LF permits platform checkout.
    bytes2hex(sha256(replace(read(rect_variant_path,String),"\r\n"=>"\n")))==
        "133cff82ad4215232d8c2c5f1bda0c9f0cdda788048e2702da3b4c564d3e02f3" &&
        rect_variants["primary_manifest_sha256"]=="e609a845958f15ab29b55b95aa3bc19183ff9fdc48b93c517fbac7d56843c377" ||
        error("rect-grid independent geometric capture provenance")
    rect_catalog["gmsh_version"]==rect_variants["gmsh_version"]=="4.15.2" &&
        rect_catalog["coordinate_abs_tolerance"]==rect_variants["coordinate_abs_tolerance"]==2e-11 &&
        rect_catalog["fixture_count"]==16 && rect_variants["fixture_count"]==8 || error("rect-grid saved oracle catalog contract")
    rect_masks=Set(Tuple(Int.(row["states"])) for row in rect_catalog["factory_masks"])
    rect_relation=Set((Tuple(Int.(row["states"])),Tuple(sort!(Int.(row["types"])))) for row in rect_catalog["factory_masks"])
    length(rect_catalog["factory_masks"])==315 && length(rect_masks)==253 || error("rect-grid immutable factory relation")
    for record in vcat(rect_catalog["fixtures"],rect_variants["fixtures"])
        name=record["payload"]["name"]
        isempty(selected) || selected=="rect_grid" || occursin(selected,"rect_grid_"*name) || continue
        empty_parameters=rect_grid_saved_case(record,rect_masks,rect_relation;oracle_only)
        rect_grid_cases[]+=1;rect_grid_quadratic_cases[]+=1
        rect_grid_empty_parameters[]+=empty_parameters
        empty_parameters>0 && (rect_grid_parameter_provenance_gaps[]+=1)
        haskey(record["payload"],"variant") && (rect_grid_variant_cases[]+=1)
        println("QUADTRI_NONEW_RECT_GRID_SAVED_OK name=$name C=0 p1_nodes=$(length(record["payload"]["p1"]["nodes"])) p2_nodes=$(length(record["payload"]["p2"]["nodes"])) oracle_empty_stored_uv=$empty_parameters native_contract=computed_uv")
    end
    println("QUADTRI_NONEW_RECT_GRID_DIFFERENTIAL_OK saved_cases=$(rect_grid_cases[]) p2_cases=$(rect_grid_quadratic_cases[]) independent_variants=$(rect_grid_variant_cases[]) parameter_provenance_gaps=$(rect_grid_parameter_provenance_gaps[]) empty_stored_uv=$(rect_grid_empty_parameters[])")
    println("QUADTRI_NONEW_FOUR_QUAD_STRIP_DIFFERENTIAL_OK saved_cases=$(four_strip_cases[]) p2_cases=$(four_strip_quadratic_cases[]) independent_variants=$(four_strip_variant_cases[]) parameter_provenance_gaps=$(four_strip_parameter_provenance_gaps[]) empty_stored_uv=$(four_strip_empty_parameters[])")
    println("QUADTRI_NONEW_THREE_QUAD_STRIP_DIFFERENTIAL_OK saved_cases=$(three_strip_cases[]) p2_cases=$(three_strip_quadratic_cases[]) independent_variants=$(three_strip_variant_cases[]) parameter_provenance_gaps=$(three_strip_parameter_provenance_gaps[]) empty_stored_uv=$(three_strip_empty_parameters[])")
    println("QUADTRI_NONEW_QUAD_STRIP_DIFFERENTIAL_OK saved_cases=$(quad_strip_cases[]) p2_cases=$(quad_strip_quadratic_cases[]) independent_variants=$(quad_strip_variant_cases[]) parameter_provenance_gaps=$(quad_strip_parameter_provenance_gaps[]) empty_stored_uv=$(quad_strip_empty_parameters[])")
    println("QUADTRI_NONEW_DIFFERENTIAL_OK gmsh=4.15.2 p1_cases=$(completed[]) p2_cases=$(quadratic_cases[]) native_helical_cases=$(native_helical_cases[]) helical_oracle_gaps=$(helical_oracle_gaps[]) height_errors=$(height_errors[]) isolated_cases=$(isolated_cases[]) triangle_cases=$(triangle_cases[]) triangle_p2_cases=$(triangle_quadratic_cases[]) triangle_isolated_cases=$(triangle_isolated_cases[]) triangle_parameter_provenance_gaps=$(triangle_parameter_provenance_gaps[]) oracle_samples=$(oracle_samples[]) two_tri_saved_cases=$(two_tri_cases[]) two_tri_p2_cases=$(two_tri_quadratic_cases[]) two_tri_parameter_provenance_gaps=$(two_tri_parameter_provenance_gaps[]) two_tri_empty_stored_uv=$(two_tri_saved_empty_parameters[]) quad_patch_saved_cases=$(quad_patch_cases[]) quad_patch_p2_cases=$(quad_patch_quadratic_cases[]) quad_patch_trapezoid_cases=$(quad_patch_trapezoid_cases[]) quad_patch_parameter_provenance_gaps=$(quad_patch_parameter_provenance_gaps[]) quad_patch_empty_stored_uv=$(quad_patch_empty_parameters[]) oracle_only=$oracle_only")
finally
    gmsh.finalize()
end
