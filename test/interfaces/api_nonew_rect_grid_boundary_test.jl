if !isdefined(@__MODULE__,:QuadTriNoNewRectGridCertificates)
    include("../geometry/quadtri_nonew_rect_grid_certificates.jl")
end
if !isdefined(@__MODULE__,:QuadTriNoNewTriangleCertificates)
    include("../geometry/quadtri_nonew_triangle_certificates.jl")
end
module NoNewRectGridBoundaryTests
using Test,Tessella,SHA
using Tessella.Elements:msh_dimension
using ..QuadTriNoNewRectGridCertificates
using ..QuadTriNoNewTriangleCertificates
const API=Tessella.API
const Cert=QuadTriNoNewRectGridCertificates
const PROFILES=(:one,:three,:graded)
const GRIDS=((2,3),(3,3))
function with_model(f,action)
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"rect_grid_boundary.geo");write(path,f.source)
            opened=API.open_geo!(path;mesh_dim=0);out=Int.(opened.lists["sweep"])
            return action((;source=f.surface,top=out[1],volume=out[2],laterals=out[3:end]))
        end
    finally
        API.finalize()
    end
end

function snapshot()
    model=API.CURRENT[]
    return (;model,value=repr(model),params=deepcopy(model.curve_params),
        cache=API.LAST_MESH[],class=API.LAST_MESH_CLASS[],
        overlay=API.LAST_MESH_HIGH_ORDER[],edges=API.LAST_MESH_EDGES[],faces=API.LAST_MESH_FACES[],
        payload=(API.mesh.get_nodes(),API.mesh.get_elements()),
        history=(API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[]))
end

function unchanged(before)
    @test API.CURRENT[]===before.model && repr(before.model)==before.value
    @test before.model.curve_params==before.params
    @test API.LAST_MESH[]===before.cache && API.LAST_MESH_CLASS[]===before.class
    @test API.LAST_MESH_HIGH_ORDER[]===before.overlay
    @test API.LAST_MESH_EDGES[]===before.edges && API.LAST_MESH_FACES[]===before.faces
    @test (API.mesh.get_nodes(),API.mesh.get_elements())==before.payload
    @test (API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])==before.history
end

function actual_carriers(carriers)
    model=API.CURRENT[]
    source=Tessella.Model.mesh_model_surface(model,carriers.source)
    volume=Tessella.Model.mesh_model_volume(model,carriers.volume)
    projected=Tessella.Model.model_to_mixed(model,volume,3,carriers.volume)
    tags,xyz,_=API.mesh.get_nodes(3,carriers.volume,true)
    points=reshape(xyz,3,:)
    positions=Dict(Tuple(p)=>tag for (tag,p) in zip(tags,eachcol(points)))
    @test length(positions)==length(tags)
    remap=[positions[Tuple(projected.coords[:,i])] for i in axes(projected.coords,2)]
    edges=Dict{Tuple,Tuple{Int,Int}}();quads=Dict{Tuple,Tuple{Int,Int}}()
    for (index,block) in enumerate(projected.blocks)
        dim=msh_dimension(block.msh);dim in (1,2) || continue
        owners=projected.elementary_entities[index]
        for (column,cell) in enumerate(eachcol(block.nodes))
            carrier=(dim,Int(owners[column]))
            if block.msh==1
                edges[Cert.key(remap[v] for v in cell)]=carrier
            else
                for i in eachindex(cell)
                    support=Cert.key((remap[cell[i]],remap[cell[mod1(i+1,length(cell))]]))
                    haskey(edges,support) && edges[support][1]<dim && continue
                    @test !haskey(edges,support) || edges[support]==carrier
                    edges[support]=carrier
                end
                block.msh==3 && (quads[Cert.key(remap[v] for v in cell)]=carrier)
            end
        end
    end
    return (;edges,quads,volume=carriers.volume,source)
end

function support_nodes(volume)
    types,tags,connectivity=API.mesh.get_elements(3,volume)
    result=Dict{UInt64,Tuple}()
    for (msh,ids,nodes) in zip(types,tags,connectivity)
        for edge in eachcol(reshape(API.mesh.get_element_edge_nodes(msh,volume,false),3,:))
            support=Cert.key(edge[1:2]);node=edge[3]
            @test !haskey(result,node) || result[node]==support
            result[node]=support
        end
        if msh in (12,13,14)
            for face in eachcol(reshape(API.mesh.get_element_face_nodes(msh,4,volume,false),9,:))
                support=Cert.key(face[1:4]);node=face[9]
                @test !haskey(result,node) || result[node]==support
                result[node]=support
            end
        end
        if msh==12
            for cell in eachcol(reshape(nodes,27,:))
                result[cell[27]]=Cert.key(cell[1:8])
            end
        end
    end
    return result
end

function check_parameters(dim,entity,boundary,count)
    tags,xyz,uv=API.mesh.get_nodes(dim,entity,boundary,true)
    @test length(tags)==count && allunique(tags)
    @test length(xyz)==3count && length(uv)==dim*count
    @test maximum(abs.(API.model.get_value(dim,entity,uv).-xyz))<=2e-11
    for node in tags
        p,parameters,d,owner=API.mesh.get_node(node)
        @test length(parameters)==d
        @test maximum(abs.(API.model.get_value(d,owner,parameters).-p))<=2e-11
    end
    return Set(tags)
end

function node_geometry(volume=nothing)
    tags=volume===nothing ? API.mesh.get_nodes()[1] : API.mesh.get_nodes(3,volume,true)[1]
    Set((Tuple(p),dim,owner,Tuple(uv)) for node in tags
        for (p,uv,dim,owner) in (API.mesh.get_node(node),))
end

function cell_geometry(volume)
    types,tags,connectivity=API.mesh.get_elements(3,volume)
    Set((Int(msh),Tuple(sort!([Tuple(API.mesh.get_node(v)[1]) for v in cell])))
        for (msh,ids,nodes) in zip(types,tags,connectivity)
        for cell in eachcol(reshape(nodes,Tessella.Elements.msh_spec(msh).nnodes,:)))
end

function file_cells(mesh;owners=true)
    return sort!([(block.msh,Tuple(cell),block.tags[column],
        !owners ? Int32(0) : mesh.entity_data===nothing ? mesh.elementary_entities[index][column] :
            mesh.entity_data.block_entities[index][column])
        for (index,block) in enumerate(mesh.blocks)
        for (column,cell) in enumerate(eachcol(block.nodes))])
end

function check_files(mesh,directory,stem)
    classified=mesh.entity_data!==nothing || mesh.elementary_entities!==nothing
    for version in (2.2,4.1),binary in (false,true)
        path=joinpath(directory,"$(stem)_$(version)_$(binary).msh")
        Tessella.Elements.write_mixed_msh(path,mesh;version,binary)
        back=Tessella.Elements.read_mixed_msh(path)
        @test back.physical_names==mesh.physical_names && back.coords==mesh.coords
        @test file_cells(back;owners=classified)==file_cells(mesh;owners=classified)
        if version==4.1 && mesh.entity_data!==nothing
            @test back.entity_data.node_entities==mesh.entity_data.node_entities
            @test back.entity_data.node_parametric==mesh.entity_data.node_parametric
        end
    end
end

function artifact_record(name,order)
    base=QuadTriNoNewTriangleCertificates.api_record(name,API)
    buffer=IOBuffer();println(buffer,"rect_grid_support_query_v1");println(buffer,base.sha)
    class=API.LAST_MESH_CLASS[]
    for (label,catalog) in (("edge",class.edge_entities),("tri",class.face_entities),("quad",class.quad_entities))
        println(buffer,label);write(buffer,UInt64(length(catalog)))
        for (support,owner) in sort!(collect(catalog);by=first)
            write(buffer,UInt64(length(support)))
            for node in support;write(buffer,UInt64(node));end
            write(buffer,Int64(owner[1]),Int64(owner[2]))
        end
    end
    for (_,surface) in sort(API.model.get_entities(2)),boundary in (false,true)
        tags,xyz,uv=API.mesh.get_nodes(2,surface,boundary,true)
        write(buffer,UInt64(surface),UInt8(boundary),UInt64(length(tags)))
        for node in tags;write(buffer,UInt64(node));end
        for array in (xyz,uv)
            write(buffer,UInt64(length(array)))
            for value in array;write(buffer,reinterpret(UInt64,Float64(value)));end
        end
    end
    return merge(base,(;kind="api_p$(order)_supports",sha=bytes2hex(sha256(take!(buffer)))))
end


function check_quadratic(carriers,f,expected)
    a,b=f.grid;n=f.intervals;boundary=2(a+b);inside=(2a-1)*(2b-1)
    primary=(a+1)*(b+1)*(n+1)
    total=(2a+1)*(2b+1)*(2n+1)
    volume_nodes=API.mesh.get_nodes(3,carriers.volume,true)[1]
    @test length(volume_nodes)==total
    counts=zeros(Int,4)
    for node in volume_nodes
        p,uv,dim,owner=API.mesh.get_node(node);counts[dim+1]+=1
        @test dim in (1,2) ? length(uv)==dim : isempty(uv)
        dim in (1,2) && (@test maximum(abs.(API.model.get_value(dim,owner,uv).-p))<=2e-11)
    end
    @test counts==[8,8n+4boundary-12,2inside+(2boundary-4)*(2n-1),inside*(2n-1)]
    supports=support_nodes(carriers.volume)
    @test length(supports)==total-primary
    signature=Dict{Tuple,Tuple}()
    for (node,support) in supports
        p,uv,dim,owner=API.mesh.get_node(node)
        carrier=length(support)==2 ? get(expected.edges,support,(3,carriers.volume)) :
            length(support)==4 ? get(expected.quads,support,(3,carriers.volume)) : (3,carriers.volume)
        @test (dim,owner)==carrier
        corners=[API.mesh.get_node(v)[1] for v in support]
        mean=[sum(c[d] for c in corners)/length(corners) for d in 1:3]
        @test maximum(abs.(p.-mean))<=2e-11
        signature[Tuple(sort!(Tuple.(corners)))]=(Tuple(p),carrier)
    end
    for (position,surface) in pairs(carriers.laterals)
        # Resolve the actual generatrix identity; output lateral order can
        # differ from the native source Curve ordering.
        link=API.CURRENT[].meshing.extrude_sources[(2,surface)]
        curve=abs(link[2]);slot=findfirst(==(curve),f.curve_tags)
        @test slot!==nothing
        segments=slot in (1,3) ? a : b
        own=check_parameters(2,surface,false,(2segments-1)*(2n-1))
        closure=check_parameters(2,surface,true,(2segments+1)*(2n+1))
        @test own⊆closure
    end
    for surface in (carriers.source,carriers.top)
        own=check_parameters(2,surface,false,inside)
        closure=check_parameters(2,surface,true,(2a+1)*(2b+1))
        @test own⊆closure
    end
    for (slot,curve) in pairs(f.curve_tags)
        segments=slot in (1,3) ? a : b
        check_parameters(1,curve,false,2segments-1)
        check_parameters(1,curve,true,2segments+1)
    end
    for msh in API.mesh.get_element_types(3,carriers.volume)
        q,_=API.mesh.get_integration_points(msh,"Gauss3")
        jac,dets,xyz=API.mesh.get_jacobians(msh,q,carriers.volume)
        @test all(isfinite,jac) && all(>(0),dets) && all(isfinite,xyz)
        @test length(jac)==9length(dets) && length(xyz)==3length(dets)
    end
    return signature
end

artifact_name(grid,profile,laterals,direction,order)=
    "rect_grid_$(grid[1])x$(grid[2])_$(profile)_R$(laterals)_D$(direction)_P$(order)"

function artifact_records()
    records=NamedTuple[]
    for grid in GRIDS,profile in PROFILES,laterals in (false,true),direction in (-1,1)
        f=Cert.fixture("artifact";grid,layers=profile,laterals,height=direction,pins=(1,2,3,4))
        with_model(f,carriers->begin
            API.mesh.generate(3)
            push!(records,artifact_record(artifact_name(grid,profile,laterals,direction,1),1))
            API.mesh.set_order(2)
            push!(records,artifact_record(artifact_name(grid,profile,laterals,direction,2),2))
        end)
    end
    return records
end

function artifact_pins()
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_rect_grid_crc.txt")
    result=Dict{String,NamedTuple}()
    for line in readlines(path)
        isempty(line) || startswith(line,"#") || begin
            fields=split(line,'\t');length(fields)==6 || error("rect grid artifact row width")
            name,kind,n_nodes,n_cells,families,sha=String.(fields)
            haskey(result,name) && error("duplicate rect grid artifact name")
            result[name]=(;name,kind,n_nodes=parse(Int,n_nodes),n_cells=parse(Int,n_cells),families,sha)
        end
    end
    length(result)==48 || error("rect grid artifact must contain48 products")
    return result
end

function run_tests()
@testset "Rectangular-grid actual P2 carriers and public lifecycle" begin
    @testset "Twenty-four P1/P2 native regression products" begin
        pins=artifact_pins()
        for grid in GRIDS,profile in PROFILES,laterals in (false,true),direction in (-1,1)
            f=Cert.fixture("golden";grid,layers=profile,laterals,height=direction,pins=(1,2,3,4))
            with_model(f,carriers->begin
                API.mesh.generate(3)
                @test length(API.mesh.get_nodes()[1])==prod(grid.+1)*length(f.levels)
                @test artifact_record(artifact_name(grid,profile,laterals,direction,1),1)==pins[artifact_name(grid,profile,laterals,direction,1)]
                expected=actual_carriers(carriers)
                @test Cert.certify(Cert.public_volume(API,carriers.volume),f,expected.source).total==1
                API.mesh.set_order(2)
                @test artifact_record(artifact_name(grid,profile,laterals,direction,2),2)==pins[artifact_name(grid,profile,laterals,direction,2)]
                initial=check_quadratic(carriers,f,expected)
                @test Cert.certify_quadratic(Cert.public_volume(API,carriers.volume),f,expected.source).total==1
                before=snapshot();API.mesh.set_order(2)
                @test (API.mesh.get_nodes(),API.mesh.get_elements())==before.payload
                API.mesh.set_order(1)
                @test length(API.mesh.get_nodes()[1])==prod(grid.+1)*length(f.levels)
                API.mesh.set_order(2)
                @test check_quadratic(carriers,f,expected)==initial
            end)
        end
    end
    @testset "Opposite winding, axes and sparse CAD tags" begin
        for plane in (:XY,:YZ,:ZX),laterals in (false,true)
            f=Cert.fixture("owners";grid=(2,3),plane,laterals,height=-1.,layers=:graded,
                shape=:rounded,winding=-1,pins=:opposite,curve_reverse=5,tags=:sparse)
            with_model(f,carriers->begin
                API.mesh.generate(3);expected=actual_carriers(carriers)
                API.mesh.set_order(2);check_quadratic(carriers,f,expected)
            end)
        end
    end
    @testset "Actual supports survive lifecycle operations" begin
        for laterals in (false,true),direction in (-1,1)
            f=Cert.fixture("lifecycle";grid=(3,3),laterals,height=direction,layers=:graded,pins=(1,2,3,4))
            with_model(f,carriers->begin
                API.mesh.generate(3);expected=actual_carriers(carriers)
                API.mesh.set_order(2);initial=check_quadratic(carriers,f,expected)
                nodes=node_geometry();cells=cell_geometry(carriers.volume);before=snapshot()
                @test_throws ArgumentError API.mesh.renumber_nodes([1,2],[1,1]);unchanged(before)
                API.mesh.reverse([(3,carriers.volume)])
                @test node_geometry()==nodes && cell_geometry(carriers.volume)==cells
                API.mesh.reverse([(3,carriers.volume)])
                primary_count=length(API.LAST_MESH_CLASS[].node_entities)
                API.mesh.renumber_nodes(collect(1:primary_count),collect(primary_count:-1:1))
                expected=actual_carriers(carriers)
                @test check_quadratic(carriers,f,expected)==initial
                for msh in API.mesh.get_element_types(3,carriers.volume)
                    tags=API.mesh.get_elements_by_type(msh,carriers.volume)[1]
                    API.mesh.reorder_elements(msh,carriers.volume,collect(length(tags)-1:-1:0))
                end
                @test node_geometry()==nodes && cell_geometry(carriers.volume)==cells
                types,tags,_=API.mesh.get_elements(3,carriers.volume)
                tag=first(first(tags));msh,connectivity,_,_=API.mesh.get_element(tag)
                removed=(Int(msh),Tuple(sort!([Tuple(API.mesh.get_node(v)[1]) for v in connectivity])))
                API.mesh.remove_elements(3,carriers.volume,[tag])
                @test cell_geometry(carriers.volume)==setdiff(cells,Set([removed]))
                API.mesh.clear();@test isempty(API.mesh.get_elements()[1])
            end)
        end
    end
    @testset "Coincident independent CAD carriers remain distinct" begin
        left=Cert.fixture("left";grid=(2,3),laterals=false,layers=:three)
        right=Cert.fixture("right";grid=(2,3),laterals=true,layers=:three,tags=:sparse,offset=(3.,0.,0.))
        pair=replace(left.source,"sweep[]="=>"left[]=")*replace(right.source,"sweep[]="=>"right[]=")
        API.initialize()
        try
            mktempdir() do directory
                path=joinpath(directory,"independent_grids.geo");write(path,pair)
                opened=API.open_geo!(path;mesh_dim=0)
                for (tag,p) in collect(API.CURRENT[].points)
                    p[1]>=3 || continue
                    API.model.set_coordinates(tag,p[1]-3,p[2],p[3])
                end
                outputs=[Int.(opened.lists[name]) for name in ("left","right")]
                carriers=[(;source=f.surface,top=out[1],volume=out[2],laterals=out[3:end])
                    for (f,out) in zip((left,right),outputs)]
                API.mesh.generate(3);expected=actual_carriers.(carriers);API.mesh.set_order(2)
                first_nodes=API.mesh.get_nodes(3,carriers[1].volume,true)[1]
                second_nodes=API.mesh.get_nodes(3,carriers[2].volume,true)[1]
                @test isempty(intersect(first_nodes,second_nodes))
                @test length(API.mesh.get_nodes()[1])==490
                for (c,f,e) in zip(carriers,(left,right),expected);check_quadratic(c,f,e);end
                before=(node_geometry(carriers[2].volume),cell_geometry(carriers[2].volume))
                API.mesh.clear([(3,carriers[1].volume)])
                @test isempty(API.mesh.get_elements(3,carriers[1].volume)[1])
                @test (node_geometry(carriers[2].volume),cell_geometry(carriers[2].volume))==before
            end
        finally
            API.finalize()
        end
    end
    @testset "Refined child supports retain complete CAD carriers" begin
        for laterals in (false,true),initial_order in (1,2)
            f=Cert.fixture("refine";grid=(2,3),laterals,layers=:one,pins=(1,2,3,4))
            with_model(f,carriers->begin
                API.mesh.generate(3);initial_order==2 && API.mesh.set_order(2)
                API.mesh.refine();types,tags,_=API.mesh.get_elements(3,carriers.volume)
                @test Dict(Int(msh)=>length(ids) for (msh,ids) in zip(types,tags))==
                    (laterals ? Dict(4=>160,7=>64) : Dict(4=>240,7=>24))
                API.mesh.set_order(2);supports=support_nodes(carriers.volume)
                cad=[(dim,entity,API.model.get_bounding_box(dim,entity)) for dim in (1,2)
                    for (_,entity) in API.model.get_entities(dim)]
                for (node,support) in supports
                    points=[API.mesh.get_node(v)[1] for v in support];owner=(3,carriers.volume)
                    for (dim,entity,bounds) in cad
                        all(all(bounds[d]<=p[d]<=bounds[d+3] for d in 1:3) for p in points) || continue
                        owner=(dim,entity);break
                    end
                    @test API.mesh.get_node(node)[3:4]==owner
                end
                for surface in carriers.laterals
                    curve=abs(API.CURRENT[].meshing.extrude_sources[(2,surface)][2])
                    slot=findfirst(==(curve),f.curve_tags);segments=slot in (1,3) ? 2 : 3
                    own=check_parameters(2,surface,false,3(4segments-1))
                    closure=check_parameters(2,surface,true,5(4segments+1));@test own⊆closure
                end
                for surface in (carriers.source,carriers.top)
                    check_parameters(2,surface,false,77);check_parameters(2,surface,true,117)
                end
            end)
        end
    end
    @testset "Actual primary and quadratic file payloads round trip" begin
        for laterals in (false,true)
            f=Cert.fixture("files";grid=(3,3),laterals,layers=:three,tags=:sparse)
            with_model(f,carriers->begin
                API.mesh.generate(3)
                linear=Tessella.Model.mesh_model_volume(API.CURRENT[],carriers.volume)
                projected=Tessella.Model.model_to_mixed(API.CURRENT[],linear,3,carriers.volume)
                API.mesh.set_order(2)
                mktempdir() do directory
                    check_files(projected,directory,"classified_grid")
                    check_files(API.mesh.get(),directory,"quadratic_grid")
                end
            end)
        end
    end
    @testset "Unrepresentable new P2 nodes reject atomically" begin
        for laterals in (false,true)
            f=Cert.fixture("precision";grid=(2,3),laterals,layers=:four,offset=(0.,0.,2.0^50+2.))
            with_model(f,carriers->begin
                API.mesh.generate(3);before=snapshot()
                @test_throws ArgumentError API.mesh.set_order(2);unchanged(before)
            end)
        end
    end
    @testset "Invalid new constraints preserve operation state" begin
        for laterals in (false,true)
            f=Cert.fixture("atomic";grid=(2,3),laterals,layers=:three)
            with_model(f,carriers->begin
                API.mesh.generate(3);API.mesh.set_order(2)
                API.mesh.set_transfinite_surface(carriers.top)
                before=snapshot();@test_throws r"QuadTriNoNewVerts" API.mesh.generate(3);unchanged(before)
            end)
        end
    end
end
end
get(ENV,"TESSELLA_RECT_GRID_API_SKIP_TESTS","0")=="1" || run_tests()
end
