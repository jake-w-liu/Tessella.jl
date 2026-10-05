if !isdefined(@__MODULE__,:QuadTriNoNewTriangleCertificates)
    include("../geometry/quadtri_nonew_triangle_certificates.jl")
end
if !isdefined(@__MODULE__,:QuadTriNoNewB4FreeStripCertificates)
    include("../geometry/quadtri_nonew_b4_free_strip_certificates.jl")
end
module NoNewB4FreeStripBoundaryTests
using Test,Tessella,SHA,TOML
using ..QuadTriNoNewB4FreeStripCertificates
using ..QuadTriNoNewTriangleCertificates
const API=Tessella.API
const Cert=QuadTriNoNewB4FreeStripCertificates

function with_model(f,action)
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"free_b4_boundary.geo");write(path,f.source)
            opened=API.open_geo!(path;mesh_dim=0);out=Int.(opened.lists["sweep"])
            action((;source=f.surface,top=out[1],volume=out[2],laterals=out[3:end]))
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
        history=(API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[]),
        options=copy(API.OPTIONS),callback=model.meshing.size_callback,
        visibility=copy(API.ELEMENT_VISIBILITY[]),locator=API.LAST_MESH_LOCATOR[],
        mixed_locator=API.LAST_MIXED_MESH_LOCATOR[])
end

function unchanged(before)
    @test API.CURRENT[]===before.model && repr(before.model)==before.value
    @test before.model.curve_params==before.params
    @test API.LAST_MESH[]===before.cache && API.LAST_MESH_CLASS[]===before.class
    @test API.LAST_MESH_HIGH_ORDER[]===before.overlay
    @test API.LAST_MESH_EDGES[]===before.edges && API.LAST_MESH_FACES[]===before.faces
    @test (API.mesh.get_nodes(),API.mesh.get_elements())==before.payload
    @test (API.NODE_TAG_MAX[],API.ELEMENT_TAG_MAX[])==before.history
    @test API.OPTIONS==before.options && before.model.meshing.size_callback===before.callback
    @test API.ELEMENT_VISIBILITY[]==before.visibility && API.LAST_MESH_LOCATOR[]===before.locator
    @test API.LAST_MIXED_MESH_LOCATOR[]===before.mixed_locator
end

function node_geometry(volume=nothing)
    tags=volume===nothing ? API.mesh.get_nodes()[1] : API.mesh.get_nodes(3,volume,true)[1]
    return Set((Tuple(p),dim,owner,Tuple(uv)) for node in tags
        for (p,uv,dim,owner) in (API.mesh.get_node(node),))
end

function cell_geometry(volume=nothing)
    types,_,connectivity=volume===nothing ? API.mesh.get_elements() : API.mesh.get_elements(3,volume)
    return Set((Int(msh),Tuple(sort!([Tuple(API.mesh.get_node(node)[1]) for node in cell])))
        for (msh,nodes) in zip(types,connectivity)
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
    buffer=IOBuffer();println(buffer,"b4_free_strip_support_query_v1");println(buffer,base.sha)
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

include("api_nonew_b4_free_strip_actual_supports.jl")

saved_records()=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","quadtri_nonew_b4_free_strip_oracle.toml"))["fixtures"]
artifact_name(record,order)="free_b4_$(record["payload"]["variant"]["name"])_P$(order)"

function artifact_pins()
    result=Dict{String,NamedTuple}()
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_b4_free_strip_crc.txt")
    for line in readlines(path)
        isempty(line) || startswith(line,"#") || begin
            fields=split(line,'\t');length(fields)==6 || error("free B4 artifact row width")
            name,kind,n_nodes,n_cells,families,sha=String.(fields)
            haskey(result,name) && error("duplicate free B4 artifact name")
            result[name]=(;name,kind,n_nodes=parse(Int,n_nodes),n_cells=parse(Int,n_cells),families,sha)
        end
    end
    length(result)==24 || error("free B4 artifact must contain 24 products")
    return result
end

function check_source_bits(carriers,source)
    points=Dict(Tuple(p)=>Tuple(reinterpret(UInt64,p)) for p in eachcol(source.coords))
    tags,xyz,_=API.mesh.get_nodes(3,carriers.volume,true)
    actual=Dict(Tuple(p)=>Tuple(reinterpret(UInt64,p)) for p in eachcol(reshape(xyz,3,:)))
    @test all(haskey(actual,p) && actual[p]==bits for (p,bits) in points)
end

function artifact_records()
    records=NamedTuple[]
    for record in saved_records()
        f=Cert.saved_fixture(record)
        with_model(f,carriers->begin
            API.mesh.generate(3);expected=free_actual_carriers(carriers)
            Cert.certify(Cert.public_volume(API,carriers.volume),f,expected.source)
            check_source_bits(carriers,expected.source)
            push!(records,artifact_record(artifact_name(record,1),1))
            API.mesh.set_order(2);free_check_actual_p2(carriers,f,expected)
            push!(records,artifact_record(artifact_name(record,2),2))
        end)
    end
    return records
end

function native_construction_cases()
    pins=artifact_pins()
    for record in saved_records()
        f=Cert.saved_fixture(record)
        with_model(f,carriers->begin
            API.mesh.generate(3);expected=free_actual_carriers(carriers)
            linear=Cert.certify(Cert.public_volume(API,carriers.volume),f,expected.source)
            @test linear.total==linear.source.area*abs(Cert.Q(last(linear.heights))-Cert.Q(first(linear.heights)))
            check_source_bits(carriers,expected.source)
            @test isempty(API.mesh.get_elements(2)[1])
            @test artifact_record(artifact_name(record,1),1)==pins[artifact_name(record,1)]
            API.mesh.set_order(2);initial=free_check_actual_p2(carriers,f,expected)
            @test artifact_record(artifact_name(record,2),2)==pins[artifact_name(record,2)]
            before=snapshot();API.mesh.set_order(2)
            @test API.LAST_MESH[]===before.cache && API.LAST_MESH_CLASS[]===before.class
            @test API.LAST_MESH_HIGH_ORDER[]===before.overlay
            @test (API.mesh.get_nodes(),API.mesh.get_elements())==before.payload
            API.mesh.set_order(1)
            @test Set(API.mesh.get_nodes(3,carriers.volume,true)[1])==expected.primary
            check_source_bits(carriers,expected.source)
            API.mesh.set_order(2)
            @test free_check_actual_p2(carriers,f,expected).supports==initial.supports
        end)
    end
end

function lower_primary_tags()
    primary=Set{UInt64}()
    types,_,connectivity=API.mesh.get_elements()
    for (msh,nodes) in zip(types,connectivity)
        msh_dimension(msh)<=1 || continue
        width=msh_spec(msh).nnodes;count=msh==15 ? 1 : 2
        for cell in eachcol(reshape(nodes,width,:)),node in cell[1:count];push!(primary,node);end
    end
    return primary
end

function check_lower_regeneration(entity,old,current,order,old_primary,new_primary)
    oldtags,oldxyz,oldparams=old;tags,xyz,params=current
    oldpoints=reshape(oldxyz,3,:);points=reshape(xyz,3,:)
    olduv=reshape(oldparams,entity[1],length(oldtags));uv=reshape(params,entity[1],length(tags))
    @test length(tags)==length(oldtags) && allunique(tags)
    @test Set(Tuple(p) for p in eachcol(uv))==Set(Tuple(p) for p in eachcol(olduv))
    signature(a,b,ids,selected)=Set((Tuple(reinterpret(UInt64,a[:,i])),Tuple(b[:,i]))
        for i in eachindex(ids) if ids[i] in selected)
    @test signature(points,uv,tags,new_primary)==signature(oldpoints,olduv,oldtags,old_primary)
    if order==2
        # The generation prepass retains primary sampling, then recreates
        # requested quadratic Curve supports by evaluating native CAD.
        @test xyz==API.model.get_value(entity[1],entity[2],params)
    else
        @test signature(points,uv,tags,Set(tags))==signature(oldpoints,olduv,oldtags,Set(oldtags))
    end
end

# Supported public attached records are an independent API elevation route.
# These seeded products do not assert generated-source cache reuse.
function seed_classified_primary(projected)
    public=UInt64[10000+3i for i in axes(projected.coords,2)]
    owners=projected.entity_data.node_entities
    for owner in sort!(unique(owners))
        dim,entity=Int.(owner)
        selected=findall(==(owner),owners)
        xyz=vec(projected.coords[:,selected])
        uv=dim in (1,2) ? API.model.get_parametrization(dim,entity,xyz) : Float64[]
        API.mesh.add_nodes(dim,entity,public[selected],xyz,uv)
    end
    next_element=50000
    for (index,block) in enumerate(projected.blocks)
        dim=msh_dimension(block.msh);owners=projected.entity_data.block_entities[index]
        for owner in sort!(unique(owners))
            selected=findall(==(owner),owners)
            elements=UInt64[next_element+7i for i in eachindex(selected)]
            next_element+=7length(selected)
            API.mesh.add_elements_by_type(Int(owner),block.msh,elements,vec(public[block.nodes[:,selected]]))
        end
    end
    @test Set(API.mesh.get_nodes()[1])==Set(public)
    return public
end

function native_record_cases()
    f=Cert.fixture("raw_closures";strip_length=5,layers=:one)
    with_model(f,carriers->begin
        positions=Dict{UInt64,Vector{Float64}}()
        point_nodes=Dict{Int,UInt64}()
        for (i,(_,entity)) in enumerate(API.model.get_entities(0))
            node=UInt64(50000+3i);point=API.model.get_value(0,entity,Float64[])
            API.mesh.add_nodes(0,entity,[node],point);positions[node]=point
            point_nodes[entity]=node
        end
        @test length(positions)==8 && Set(API.mesh.get_nodes()[1])==Set(keys(positions))
        for dim in (1,2,3)
            emitted=UInt64[]
            for (_,entity) in API.model.get_entities(dim)
                bounds=API.model.get_bounding_box(dim,entity)
                expected=Set(node for (node,p) in positions if all(bounds[k]<=p[k]<=bounds[k+3] for k in 1:3))
                tags,xyz,uv=API.mesh.get_nodes(dim,entity,true,true)
                @test Set(tags)==expected && allunique(tags)
                if dim in (2,3)
                    # GFace/GRegion vertices() uses GEntityPtrLessThan (tag
                    # order), independently confirmed by the raw primary query.
                    @test tags==[point_nodes[vertex] for vertex in sort!(collect(keys(point_nodes)))
                        if point_nodes[vertex] in expected]
                end
                @test isempty(API.mesh.get_nodes(dim,entity,false)[1])
                @test length(tags)==(dim==1 ? 2 : dim==2 ? 4 : 8)
                @test all(positions[node]==p for (node,p) in zip(tags,eachcol(reshape(xyz,3,:))))
                # Independent raw native queries have inverse-map residuals up
                # to 4.44e-16; retain the established CAD roundtrip tolerance.
                @test dim==3 ? isempty(uv) : length(uv)==dim*length(tags) &&
                    maximum(abs.(API.model.get_value(dim,entity,uv).-xyz))<=2e-11
                append!(emitted,tags)
            end
            aggregate=API.mesh.get_nodes(dim,-1,true,true)
            @test aggregate[1]==emitted
            @test length(aggregate[3])==(dim==3 ? 0 : dim*length(emitted))
        end
        global_boundary=API.mesh.get_nodes(-1,-1,true,true)
        @test length(global_boundary[1])==64 && isempty(global_boundary[3])
        @test Set(global_boundary[1])==Set(keys(positions))
    end)
    API.initialize()
    try
        text="Point(1)={.4,.3,.5,1};\nPoint(2)={.5,.4,.5,1};\nPoint(3)={.6,.3,.5,1};\nPoint(4)={.5,.2,.5,1};\nSpline(10)={1,2,3,4,1};\n"
        mktemp() do path,io
            write(io,text);close(io);API.open_geo!(path;mesh_dim=0)
        end
        API.mesh.add_nodes(0,1,[50001],[.4,.3,.5])
        @test API.model.get_boundary([(1,10)],false,false,false)==[(0,1),(0,1)]
        @test API.mesh.get_nodes(1,10,true,true)==
            (UInt64[50001,50001],[.4,.3,.5,.4,.3,.5],[0.,0.])
        @test API.mesh.get_nodes(-1,-1,true,false)==
            (UInt64[50001,50001,50001],[.4,.3,.5,.4,.3,.5,.4,.3,.5],Float64[])
    finally
        API.finalize()
    end
    extra=merge(f,(;source=f.source*"\nPoint(40000)={10,10,10,1};\n"))
    with_model(extra,carriers->begin
        API.option("Mesh.Renumber",0)
        API.model.add_discrete_entity(0,41000)
        for (entity,node,element,point) in ((40000,50001,60001,[10.,10.,10.]),(41000,50003,60003,[12.,13.,14.]))
            API.mesh.add_nodes(0,entity,[node],point)
            API.mesh.add_elements_by_type(entity,15,[element],[node])
        end
        API.mesh.generate(1)
        payload=Dict(entity=>(API.mesh.get_nodes(0,entity),API.mesh.get_elements(0,entity)) for entity in (40000,41000))
        for order in (1,2)
            API.mesh.generate(3);order==2 && API.mesh.set_order(2)
            @test allunique(API.mesh.get_nodes()[1])
            for (entity,before) in payload
                @test (API.mesh.get_nodes(0,entity),API.mesh.get_elements(0,entity))==before
            end
            @test !isempty(API.mesh.get_elements(3,carriers.volume)[1])
        end
    end)
    for renumber in (0,1)
    with_model(f,carriers->begin
        API.option("Mesh.Renumber",renumber)
        native_nodes=UInt64[50001,50002,50003,50004]
        native_coords=[.2,.2,.2,.4,.2,.2,.2,.4,.2,.2,.2,.4]
        API.mesh.add_nodes(3,carriers.volume,native_nodes,native_coords)
        API.mesh.add_elements_by_type(carriers.volume,4,[60001],native_nodes)
        API.model.add_discrete_entity(3,42000)
        discrete_nodes=UInt64[70001,70002,70003,70004]
        discrete_coords=[10.,0.,0.,11.,0.,0.,10.,1.,0.,10.,0.,1.]
        API.mesh.add_nodes(3,42000,discrete_nodes,discrete_coords)
        API.mesh.add_elements_by_type(42000,4,[80001],discrete_nodes)
        before=(API.mesh.get_nodes(3,42000),API.mesh.get_elements(3,42000))
        API.mesh.generate(2)
        @test isempty(API.mesh.get_nodes(3,carriers.volume)[1]) && isempty(API.mesh.get_elements(3,carriers.volume)[1])
        if renumber==0
            @test (API.mesh.get_nodes(3,42000),API.mesh.get_elements(3,42000))==before
        else
            # Pinned primary renumbers uncovered records in the same global
            # pass. The discrete Volume follows every native lower entity.
            @test API.mesh.get_nodes(3,42000)==
                (UInt64[25,26,27,28],discrete_coords,Float64[])
            all_elements=vcat(API.mesh.get_elements()[2]...)
            @test sort(all_elements)==UInt64.(1:length(all_elements))
            @test API.mesh.get_elements(3,42000)==
                (Int32[4],[UInt64[length(all_elements)]],[UInt64[25,26,27,28]])
        end
        @test isempty(intersect(Set(native_nodes),Set(API.mesh.get_nodes()[1])))
        @test allunique(API.mesh.get_nodes()[1])
    end)
    end
    # Native allocation keeps every surviving public identity and allocates
    # independently above its history; raw Point labels can coincide with the
    # former dense cache indices without aliasing actual primary vertices.
    for (node,element) in ((1,60001),(50001,1))
        with_model(extra,carriers->begin
            API.option("Mesh.Renumber",0)
            API.mesh.add_nodes(0,40000,[node],[10.,10.,10.])
            API.mesh.add_elements_by_type(40000,15,[element],[node])
            API.mesh.generate(0)
            retained=(API.mesh.get_nodes(0,40000),API.mesh.get_elements(0,40000))
            API.mesh.generate(3)
            @test (API.mesh.get_nodes(0,40000),API.mesh.get_elements(0,40000))==retained
            @test API.mesh.get_node(node)==([10.,10.,10.],Float64[],0,40000)
            @test API.mesh.get_element(element)==(15,UInt64[node],0,40000)
            expected=free_actual_carriers(carriers)
            @test Cert.certify(Cert.public_volume(API,carriers.volume),extra,expected.source).total==1
            @test !(UInt64(node) in expected.primary)
            tags=API.mesh.get_nodes()[1];elements=vcat(API.mesh.get_elements()[2]...)
            @test allunique(tags) && allunique(elements)
            @test all(n in Set(tags) for nodes in API.mesh.get_elements()[3] for n in nodes)
            API.mesh.set_order(2);free_check_actual_p2(carriers,extra,expected)
            @test (API.mesh.get_nodes(0,40000),API.mesh.get_elements(0,40000))==retained
        end)
    end
# Native Curve meshing calls deMeshGEdge; an uncovered discrete Face cannot
# silently retain a reference to that discarded Curve-owned node identity.
with_model(f,carriers->begin
    API.option("Mesh.Renumber",0)
    curve=f.curve_tags[first(f.long_pair)]
    point=API.model.get_value(1,curve,[.5])
    API.mesh.add_nodes(1,curve,[100001],point,[.5])
    API.model.add_discrete_entity(2,41001,[curve])
    API.mesh.add_nodes(2,41001,[100003,100005],[12.,13.,14.,13.,13.,14.])
    API.mesh.add_elements_by_type(41001,2,[120001],[100001,100003,100005])
    API.mesh.generate(0);before=snapshot()
    @test_throws r"surviving record.*100001.*no retained-node identity" API.mesh.generate(3)
    unchanged(before)
end)

    for policy in (:free,:rec)
        recipe=policy==:free ? extra : merge(extra,(;laterals=true,
            source=replace(extra.source,"QuadTriNoNewVerts"=>"QuadTriNoNewVerts RecombLaterals")))
        with_model(recipe,carriers->begin
            API.option("Mesh.Renumber",0)
            API.mesh.add_nodes(0,40000,[80],[10.,10.,10.])
            API.mesh.add_elements_by_type(40000,15,[60001],[80])
            API.mesh.generate(3)
            retained=(API.mesh.get_nodes(0,40000),API.mesh.get_elements(0,40000))
            expected=free_actual_carriers(carriers)
            @test !(UInt64(80) in expected.primary)
            @test allunique(API.mesh.get_nodes()[1])
            API.mesh.set_order(2)
            result=free_check_actual_p2(carriers,recipe,expected;certify_geometry=policy==:free)
            if policy==:rec
                # The actual support, ownership and typed traces above are
                # shared; the recombined shell has its own independent proof.
                Cert.Rec.certify_quadratic(Cert.public_volume(API,carriers.volume),recipe,expected.source)
            end
            @test !(UInt64(80) in keys(result.supports))
            @test (API.mesh.get_nodes(0,40000),API.mesh.get_elements(0,40000))==retained
            @test API.mesh.get_node(80)==([10.,10.,10.],Float64[],0,40000)
            @test allunique(API.mesh.get_nodes()[1]) && allunique(vcat(API.mesh.get_elements()[2]...))
            API.mesh.set_order(1)
            @test (API.mesh.get_nodes(0,40000),API.mesh.get_elements(0,40000))==retained
            certificate=policy==:free ? Cert.certify : Cert.Rec.certify
            @test certificate(Cert.public_volume(API,carriers.volume),recipe,expected.source).total==1
        end)
    end
    API.initialize()
    try
        API.option("Mesh.Renumber",0)
        API.model.add_box(0,0,0,1,1,1;tag=1)
        API.model.add_discrete_entity(0,41000)
        generated=API.mesh.generate(3)
        @test generated isa Tessella.Mesh && API.LAST_MESH[] isa Tessella.Mesh
        node=size(generated.coords,2)+1
        API.mesh.add_nodes(0,41000,[node],[10.,10.,10.])
        retained=API.mesh.get_nodes(0,41000)
        tags=Set(API.mesh.get_nodes(3,1,true)[1])
        @test !(UInt64(node) in tags) && allunique(API.mesh.get_nodes()[1])
        cells=API.mesh.get_elements(3,1)[3]
        p1total=zero(Cert.Q)
        for connectivity in cells,cell in eachcol(reshape(connectivity,4,:))
            p1total+=Cert.cell_volume(Tuple(Cert.exact(API.mesh.get_node(n)[1]) for n in cell),4)
        end
        @test p1total==1
        API.mesh.set_order(2)
        @test API.LAST_MESH[] isa Tessella.Mesh
        @test API.mesh.get_nodes(0,41000)==retained && API.mesh.get_node(node)==([10.,10.,10.],Float64[],0,41000)
        @test API.mesh.get_element_types(3,1)==Int32[11]
        supports=free_public_support_nodes(1)
        @test !(UInt64(node) in keys(supports))
        for (midpoint,edge) in supports
            xyz=API.mesh.get_node(midpoint)[1]
            mean=sum(API.mesh.get_node(n)[1] for n in edge)/2
            @test maximum(abs.(xyz.-mean))<=2e-11
        end
        total=zero(Cert.Q)
        for cell in eachcol(reshape(API.mesh.get_elements_by_type(11,1)[2],10,:))
            certificate=Cert.quadratic_map_certificate(11,Tuple(API.mesh.get_node(n)[1] for n in cell))
            @test certificate.volume>0
            total+=certificate.volume
        end
        @test total==p1total && allunique(API.mesh.get_nodes()[1])
        API.mesh.set_order(1)
        @test API.mesh.get_nodes(0,41000)==retained && API.LAST_MESH[] isa Tessella.Mesh
    finally
        API.finalize()
    end

end

# Pinned Gmsh 4.15.2 retains Point identities through these six operations.
function native_point_identity_cases()
    for renumber in (0,1)
        API.initialize()
        try
            API.option("Mesh.Renumber",renumber)
            API.option("Mesh.ElementOrder",1)
            API.model.add_box(0,0,0,1,1,1;tag=1)
            point=API.model.get_value(0,1,Float64[])
            @test point==[0.,0.,1.]
            API.mesh.add_nodes(0,1,[111],point)
            API.mesh.add_elements_by_type(1,15,[777],[111])
            for (_,curve) in API.model.get_entities(1)
                API.mesh.set_transfinite_curve(curve,2)
            end
            for (_,surface) in API.model.get_entities(2)
                API.mesh.set_transfinite_surface(surface)
                API.mesh.set_recombine(2,surface)
            end
            API.mesh.set_transfinite_volume(1)
            node=UInt64(renumber==0 ? 111 : 1)
            element=UInt64(renumber==0 ? 777 : 1)
            for (step,volume_type,expected_nodes) in ((:generate3,5,8),(:order2,12,27),
                    (:order1,5,8),(:generate3,5,8),(:generate2,0,8),(:generate3,5,8))
                if step===:order2
                    API.mesh.set_order(2)
                elseif step===:order1
                    API.mesh.set_order(1)
                else
                    API.mesh.generate(step===:generate2 ? 2 : 3)
                end
                @test API.mesh.get_nodes(0,1,true,true)==([node],point,Float64[])
                @test API.mesh.get_elements(0,1)==(Int32[15],[[element]],[[node]])
                @test API.mesh.get_node(node)==(point,Float64[],0,1)
                tags,xyz,params=API.mesh.get_nodes()
                @test length(tags)==expected_nodes && allunique(tags) && count(==(node),tags)==1
                @test isempty(params) && length(xyz)==3expected_nodes
                types,ids,connectivity=API.mesh.get_elements()
                @test allunique(vcat(ids...)) && all(n in Set(tags) for nodes in connectivity for n in nodes)
                if volume_type==0
                    @test isempty(API.mesh.get_elements(3,1)[1])
                else
                    @test API.mesh.get_element_types(3,1)==Int32[volume_type]
                    @test length(API.mesh.get_elements_by_type(volume_type,1)[1])==1
                end
            end
            # Higher-to-lower generation retains the Point but remeshes Curves.
            for (_,curve) in API.model.get_entities(1)
                API.mesh.set_transfinite_curve(curve,3)
            end
            API.mesh.generate(1)
            @test API.mesh.get_nodes(0,1)==([node],point,Float64[])
            @test API.mesh.get_elements(0,1)==(Int32[15],[[element]],[[node]])
            @test isempty(API.mesh.get_elements(3,1)[1])
            @test length(API.mesh.get_nodes()[1])==20 && allunique(API.mesh.get_nodes()[1])
        finally
            API.finalize()
        end
    end
end

function legacy_unclassified_aggregate_cases()
    # Preserve the supported detached-cache contract pinned in
    # api_mesh_data_test.jl, including unrelated CAD metadata. The fixture
    # installs real Mesh data through the existing cache utility; queries use
    # the public API without replacing any implementation method.
    coordinates=Float64[0 1 0 0;0 0 1 0;0 0 0 1]
    for cad in (false,true)
        API.initialize()
        try
            cad && API.model.add_point(10.,10.,10.;tag=40000)
            mesh=Tessella.Mesh(coordinates;tets=reshape(Int32[1,2,3,4],4,1))
            lock(API.STATE_LOCK) do
                API._replace_mesh_cache_locked!(API._copy_mesh(mesh))
            end
            expected=(UInt64[1,2,3,4],vec(coordinates),Float64[])
            @test API.mesh.get_nodes()==expected
            @test API.mesh.get_nodes(-1,-2,true,false)==expected
            @test API.mesh.get_nodes(-1,40000,true,true)==expected
        finally
            API.finalize()
        end
    end
end


function run_tests()
@testset "Dynamic free B4 actual API contracts" begin
    @testset "Twelve original native public P1/P2 recipes" begin
        native_construction_cases()
    end
    @testset "Sparse CAD identities, opposite chains and signed frames" begin
        for (m,plane,long_pair,winding) in ((5,:XY,(1,3),1),(7,:YZ,(2,4),-1),(9,:XZ,(1,3),-1)),direction in (-1,1)
            f=Cert.fixture("sparse";strip_length=m,plane,height=.75direction,layers=:nonbinary,
                shape=:rounded,tags=:sparse,pins=(3,4,1,2),curve_reverse=5,winding,
                long_pair,law="Progression",coefficient=4.)
            with_model(f,carriers->begin
                API.mesh.generate(3);expected=free_actual_carriers(carriers)
                Cert.certify(Cert.public_volume(API,carriers.volume),f,expected.source)
                check_source_bits(carriers,expected.source)
                API.mesh.set_order(2);free_check_actual_p2(carriers,f,expected)
            end)
        end
    end
    @testset "Actual attached CAD closures and independent records" begin
        native_record_cases()
        native_point_identity_cases()
        legacy_unclassified_aggregate_cases()
    end
    @testset "Actual retained fans through supported public seeded P2" begin
        for product in Cert.retained_products()
            f=product.f
            with_model(f,carriers->begin
                API.option("Mesh.Renumber",0)
                public=seed_classified_primary(product.projected)
                expected=free_actual_carriers(carriers;source=product.source,
                    volume=product.volume,projected=product.projected,public_tags=public)
                @test allunique(API.mesh.get_nodes()[1])
                for (entity,closure) in expected.closures
                    entity[1] in (1,2) || continue
                    owned=Set(node for (node,owner) in expected.owners if owner==entity)
                    free_check_parameters(entity,closure,owned)
                end
                linear=Cert.certify(Cert.public_volume(API,carriers.volume),f,product.source)
                @test length(linear.centers)==product.expected_centers
                @test product.projected.coords==product.volume.coords
                centers=public[collect(product.certificate.centers)]
                @test length(centers)==product.expected_centers
                @test all(API.mesh.get_node(node)[3:4]==(3,carriers.volume) for node in centers)
                if product.nonterminal
                    @test any(min(first(linear.heights),linear.heights[end-1])<API.mesh.get_node(node)[1][f.axes[3]]<
                        max(first(linear.heights),linear.heights[end-1]) for node in centers)
                end
                primary_geometry=node_geometry()
                API.mesh.set_order(2);quadratic=free_check_actual_p2(carriers,f,expected)
                for center in centers
                    radial=[node for (node,support) in quadratic.supports if center in support]
                    @test length(radial)==8
                    @test all(API.mesh.get_node(node)[3:4]==(3,carriers.volume) for node in radial)
                    @test API.mesh.get_node(center)[3:4]==(3,carriers.volume)
                end
                before=snapshot();API.mesh.set_order(2)
                @test (API.mesh.get_nodes(),API.mesh.get_elements())==before.payload
                API.mesh.set_order(1);@test node_geometry()==primary_geometry
                API.mesh.set_order(2)
                @test Set(values(free_check_actual_p2(carriers,f,expected).supports))==Set(values(quadratic.supports))
                geometry=node_geometry();cells=cell_geometry()
                tags=API.mesh.get_nodes()[1];labels=UInt64[20000+5i for i in eachindex(tags)]
                API.mesh.renumber_nodes(tags,labels)
                @test Set(API.mesh.get_nodes()[1])==Set(labels)
                @test node_geometry()==geometry && cell_geometry()==cells
                types,ids,_=API.mesh.get_elements(3,carriers.volume)
                for group in ids;API.mesh.renumber_elements(group,reverse(group));end
                selected=UInt64.(vcat(API.mesh.get_elements(3,carriers.volume)[2]...))
                API.mesh.set_visibility(selected,0);@test all(==(0),API.mesh.get_visibility(selected))
                API.mesh.set_visibility(selected,1);@test all(==(1),API.mesh.get_visibility(selected))
                before=API.mesh.get_elements(3,carriers.volume)
                API.mesh.reverse([(3,carriers.volume)]);API.mesh.reverse([(3,carriers.volume)])
                @test API.mesh.get_elements(3,carriers.volume)==before
                @test node_geometry()==geometry && cell_geometry()==cells
                if product.nonterminal && f.height>0
                    mktempdir() do directory;check_files(API.mesh.get(),directory,"seeded_nonterminal_fan");end
                    API.mesh.set_order(1)
                    types,ids,_=API.mesh.get_elements(3,carriers.volume)
                    parents=Dict(Int(msh)=>length(tags) for (msh,tags) in zip(types,ids))
                    API.mesh.refine();types,ids,_=API.mesh.get_elements(3,carriers.volume)
                    children=Dict(4=>8*(get(parents,4,0)+get(parents,7,0)),
                        5=>8get(parents,5,0),6=>8get(parents,6,0),7=>4get(parents,7,0))
                    filter!(pair->last(pair)>0,children)
                    @test Dict(Int(msh)=>length(tags) for (msh,tags) in zip(types,ids))==children
                    API.mesh.set_order(2);refined=Cert.public_volume(API,carriers.volume)
                    total=zero(Cert.Q)
                    for block in refined.blocks,cell in eachcol(block.nodes)
                        map=Cert.quadratic_map_certificate(Int(block.msh),Tuple(Cert.point(refined,node) for node in cell))
                        @test map.volume>0;total+=map.volume
                    end
                    @test abs(total-linear.total)<=Cert.Q(2e-11)
                    for point in (Tuple(product.volume.coords[:,node]) for node in product.certificate.centers)
                        preserved=[node for node in API.mesh.get_nodes(3,carriers.volume,true)[1]
                            if Tuple(API.mesh.get_node(node)[1])==point]
                        @test length(preserved)==1 && API.mesh.get_node(only(preserved))[3:4]==(3,carriers.volume)
                    end
                end
                lower=Dict((dim,tag)=>API.mesh.get_nodes(dim,tag,false,true)
                    for dim in (0,1,2) for (_,tag) in API.model.get_entities(dim))
                history=(API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag())
                API.mesh.clear([(3,carriers.volume)])
                @test isempty(API.mesh.get_elements(3,carriers.volume)[1])
                @test (API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag())==history
                for (entity,payload) in lower;@test API.mesh.get_nodes(entity...,false,true)==payload;end
            end)
        end
    end
    @testset "Real source and volume generation across all dimensions" begin
        f=Cert.fixture("dimensions";strip_length=7,layers=:nonbinary,plane=:YZ,
            shape=:skew,height=-.75,law="Progression",coefficient=4.)
        with_model(f,carriers->begin
            for dim in (0,1,2)
                API.mesh.generate(dim)
                types,_,_=API.mesh.get_elements()
                # Pinned Gmsh generate(0) is an empty pass; Mesh0D is run by
                # generate(1). Positive dimensions must contain actual cells.
                @test dim==0 ? isempty(types) && isempty(API.mesh.get_nodes()[1]) :
                    !isempty(types) && maximum(msh_dimension.(types))==dim
                @test isempty(API.mesh.get_elements(3,carriers.volume)[1])
                @test allunique(API.mesh.get_nodes()[1])
                for dimension in 0:dim,(_,entity) in API.model.get_entities(dimension)
                    tags,xyz,_=API.mesh.get_nodes(dimension,entity,false,true)
                    for (node,p) in zip(tags,eachcol(reshape(xyz,3,:)))
                        actual,_,owner_dim,owner=API.mesh.get_node(node)
                        @test (owner_dim,owner)==(dimension,entity) && actual==p
                    end
                end
            end
            types,ids,_=API.mesh.get_elements(2,carriers.source)
            @test Dict(Int(msh)=>length(tags) for (msh,tags) in zip(types,ids))==Dict(3=>7)
            source_nodes=Set(Tuple(reinterpret(UInt64,API.mesh.get_node(node)[1]))
                for node in API.mesh.get_nodes(2,carriers.source,true)[1])
            @test length(source_nodes)==16
            API.mesh.generate(3);expected=free_actual_carriers(carriers)
            @test Set(Tuple(reinterpret(UInt64,p)) for p in eachcol(expected.source.coords))==source_nodes
            Cert.certify(Cert.public_volume(API,carriers.volume),f,expected.source)
            API.mesh.set_order(2);free_check_actual_p2(carriers,f,expected)
        end)
    end
    @testset "OnlyEmpty lower regeneration retains native primary sampling" begin
        for order in (1,2)
            f=Cert.fixture("bridge";strip_length=5,layers=:nonbinary)
            with_model(f,carriers->begin
                API.mesh.generate(3);primary=Set(API.mesh.get_nodes()[1])
                order==2 && API.mesh.set_order(2)
                before=Dict((dim,tag)=>API.mesh.get_nodes(dim,tag,false,true)
                    for dim in (0,1) for (_,tag) in API.model.get_entities(dim))
                history=(API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag())
                API.option("Mesh.ElementOrder",order);API.option("Mesh.MeshOnlyEmpty",1)
                API.mesh.generate(1);new_primary=lower_primary_tags()
                for (entity,payload) in before
                    check_lower_regeneration(entity,payload,
                        API.mesh.get_nodes(entity...,false,true),order,primary,new_primary)
                end
                @test API.mesh.get_max_node_tag()>=history[1] && API.mesh.get_max_element_tag()>=history[2]
                geometry=node_geometry();cells=cell_geometry()
                tags=API.mesh.get_nodes()[1];labels=UInt64[10000+3i for i in eachindex(tags)]
                API.mesh.renumber_nodes(tags,labels)
                @test Set(API.mesh.get_nodes()[1])==Set(labels)
                @test node_geometry()==geometry && cell_geometry()==cells
            end)
        end
    end
    @testset "Actual geometry, cache and selections survive lifecycle operations" begin
        f=Cert.fixture("lifecycle";strip_length=5,layers=:nonbinary,tags=:sparse)
        with_model(f,carriers->begin
            API.mesh.generate(3);API.mesh.set_order(2)
            geometry=node_geometry();cells=cell_geometry()
            nodes=API.mesh.get_nodes()[1];API.mesh.renumber_nodes(nodes,reverse(nodes))
            @test node_geometry()==geometry && cell_geometry()==cells
            groups=API.mesh.get_elements()[2];elements=UInt64.(vcat(groups...))
            before_elements=Dict(element=>API.mesh.get_element(element) for element in elements)
            API.mesh.renumber_elements(elements,reverse(elements))
            @test all(API.mesh.get_element(new)==before_elements[old]
                for (old,new) in zip(elements,reverse(elements)))
            @test node_geometry()==geometry && cell_geometry()==cells
            before=snapshot()
            @test_throws r"duplicate element tag" API.mesh.renumber_elements(elements[1:2],fill(elements[1],2))
            unchanged(before)
            for group in API.mesh.get_elements()[2];API.mesh.renumber_elements(group,reverse(group));end
            @test node_geometry()==geometry && cell_geometry()==cells
            selected=first(API.mesh.get_elements(3,carriers.volume)[2])
            API.mesh.set_visibility(selected,0)
            @test all(==(0),API.mesh.get_visibility(selected))
            @test node_geometry()==geometry && cell_geometry()==cells
            API.mesh.set_visibility(selected,1)
            @test all(==(1),API.mesh.get_visibility(selected))
            API.mesh.reverse([(3,carriers.volume)]);API.mesh.reverse([(3,carriers.volume)])
            @test node_geometry()==geometry && cell_geometry()==cells
            for msh in API.mesh.get_element_types(3,carriers.volume)
                count=length(API.mesh.get_elements_by_type(msh,carriers.volume)[1])
                API.mesh.reorder_elements(msh,carriers.volume,collect(count-1:-1:0))
            end
            @test node_geometry()==geometry && cell_geometry()==cells
            before=snapshot()
            @test_throws ArgumentError API.mesh.clear([(3,carriers.volume),(3,90001)])
            unchanged(before)
            history=(API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag())
            API.mesh.clear([(3,carriers.volume)])
            @test isempty(API.mesh.get_elements(3,carriers.volume)[1])
            @test (API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag())==history
        end)
    end
    @testset "Independent coincident regions preserve entity identity" begin
        left=Cert.fixture("left";strip_length=5,layers=:one)
        right=Cert.fixture("right";strip_length=5,layers=:one,tags=:sparse,offset=(3.,0.,0.))
        source=replace(left.source,"sweep[]"=>"left[]")*replace(right.source,"sweep[]"=>"right[]")
        API.initialize()
        try
            mktempdir() do directory
                path=joinpath(directory,"independent_free_b4.geo");write(path,source)
                opened=API.open_geo!(path;mesh_dim=0)
                for (_,tag) in API.model.get_entities(0)
                    p=API.model.get_value(0,tag,Float64[]);p[1]>=3 || continue
                    API.model.set_coordinates(tag,p[1]-3,p[2],p[3])
                end
                output=[Int.(opened.lists[name]) for name in ("left","right")]
                carriers=[(;source=f.surface,top=out[1],volume=out[2],laterals=out[3:end])
                    for (f,out) in zip((left,right),output)]
                API.mesh.generate(3)
                @test allunique(API.mesh.get_nodes()[1])
                expected=free_actual_carriers.(carriers)
                API.mesh.set_order(2)
                sets=[Set(API.mesh.get_nodes(3,c.volume,true)[1]) for c in carriers]
                @test isempty(intersect(sets...)) && Set(API.mesh.get_nodes()[1])==union(sets...)
                @test Set(first.(node_geometry(carriers[1].volume)))==Set(first.(node_geometry(carriers[2].volume)))
                shifted_right=merge(right,(;corners=Tuple((p[1]-3,p[2],p[3]) for p in right.corners),offset=(0.,0.,0.)))
                for (point,position) in zip(right.point_tags,shifted_right.corners)
                    @test Tuple(API.model.get_value(0,point,Float64[]))==position
                end
                for (c,f,e) in zip(carriers,(left,shifted_right),expected);free_check_actual_p2(c,f,e);end
                before=(node_geometry(carriers[2].volume),cell_geometry(carriers[2].volume))
                API.mesh.clear([(3,carriers[1].volume)])
                @test isempty(API.mesh.get_elements(3,carriers[1].volume)[1])
                @test (node_geometry(carriers[2].volume),cell_geometry(carriers[2].volume))==before
            end
        finally
            API.finalize()
        end
    end
    @testset "Refined actual child maps and support carriers" begin
        for m in (5,7),initial_order in (1,2)
            f=Cert.fixture("refine";strip_length=m,layers=:one,pins=(1,2,3,4))
            with_model(f,carriers->begin
                API.mesh.generate(3);initial_order==2 && API.mesh.set_order(2)
                types,ids,_=API.mesh.get_elements(3,carriers.volume)
                parent_counts=Dict(Int(msh)=>length(tags) for (msh,tags) in zip(types,ids))
                parent(msh)=get(parent_counts,initial_order==2 ? Dict(4=>11,5=>12,6=>13,7=>14)[msh] : msh,0)
                children=Dict(4=>8*(parent(4)+parent(7)),5=>8parent(5),6=>8parent(6),7=>4parent(7))
                filter!(pair->last(pair)>0,children)
                API.mesh.refine();types,ids,_=API.mesh.get_elements(3,carriers.volume)
                @test Dict(Int(msh)=>length(tags) for (msh,tags) in zip(types,ids))==children
                primary=Cert.public_volume(API,carriers.volume);total=zero(Cert.Q)
                for block in primary.blocks,cell in eachcol(block.nodes)
                    points=Tuple(Cert.exact(Cert.point(primary,node)) for node in cell)
                    @test Cert.map_positive(points,Int(block.msh))
                    total+=Cert.cell_volume(points,Int(block.msh))
                end
                @test total==1
                API.mesh.set_order(2);supports=free_public_support_nodes(carriers.volume)
                cad=[(dim,entity,API.model.get_bounding_box(dim,entity)) for dim in (1,2)
                    for (_,entity) in API.model.get_entities(dim)]
                for (node,support) in supports
                    corners=[API.mesh.get_node(v)[1] for v in support];owner=(3,carriers.volume)
                    for (dim,entity,bounds) in cad
                        all(all(bounds[d]<=p[d]<=bounds[d+3] for d in 1:3) for p in corners) || continue
                        owner=(dim,entity);break
                    end
                    @test API.mesh.get_node(node)[3:4]==owner
                end
                quadratic=Cert.public_volume(API,carriers.volume);qtotal=zero(Cert.Q)
                for block in quadratic.blocks,cell in eachcol(block.nodes)
                    certificate=Cert.quadratic_map_certificate(Int(block.msh),Tuple(Cert.point(quadratic,node) for node in cell))
                    @test certificate.volume>0
                    qtotal+=certificate.volume
                end
                @test abs(qtotal-total)<=Cert.Q(2e-11)
            end)
        end
    end
    @testset "Physical metadata and actual P2 file payloads" begin
        f=Cert.fixture("files";strip_length=7,layers=:nonbinary,tags=:sparse)
        text=f.source*"\nPhysical Curve(\"source lines\",501)={$(join(f.curve_tags,','))};\n"*
            "Physical Surface(\"caps\",502)={$(f.surface),sweep[0]};\n"*
            "Physical Surface(\"walls\",504)={sweep[2],sweep[3],sweep[4],sweep[5]};\n"*
            "Physical Volume(\"body\",503)={sweep[1]};\n"
        with_model(merge(f,(;source=text)),carriers->begin
            names=Dict((1,501)=>"source lines",(2,502)=>"caps",(2,504)=>"walls",(3,503)=>"body")
            @test Set(API.model.get_physical_groups())==Set(keys(names))
            @test API.model.get_entities_for_physical_group(3,503)==[carriers.volume]
            @test Set(API.model.get_entities_for_physical_group(2,502))==Set((carriers.source,carriers.top))
            for (entity,name) in names;@test API.model.get_physical_name(entity...)==name;end
            API.mesh.generate(3);expected=free_actual_carriers(carriers)
            @test expected.projected.physical_names==names
            @test all(all(==(503),block.tags) for block in expected.projected.blocks if msh_dimension(block.msh)==3)
            API.mesh.set_order(2);free_check_actual_p2(carriers,f,expected)
            mktempdir() do directory
                check_files(expected.projected,directory,"classified_free_b4")
                check_files(API.mesh.get(),directory,"actual_quadratic_free_b4")
            end
        end)
    end
    @testset "Checked public resource and source controls fail atomically" begin
        f=Cert.fixture("node_reserve";strip_length=5,layers=:one)
        # Columns alone fit the ten-million-node bound; reserved retained
        # centers push the checked preflight beyond it, before level allocation.
        @test 12*(588235+1)<=10_000_000<12*(588235+1)+5*588235
        oversized=merge(f,(;source=replace(f.source,"Layers{1}"=>"Layers{588235}")))
        with_model(oversized,carriers->begin
            API.mesh.generate(0);API.mesh.set_order(2)
            before=snapshot()
            @test_throws r"QuadTriNoNewVerts.*(resource|node|bound|large)" API.mesh.generate(3)
            unchanged(before)
        end)
        with_model(f,carriers->begin
            API.mesh.generate(0);API.mesh.set_order(2)
            for slot in f.long_pair;API.mesh.set_transfinite_curve(f.curve_tags[slot],typemax(Int32));end
            before=snapshot()
            @test_throws r"(QuadTriNoNewVerts|transfinite).*(resource|node|bound|large)" API.mesh.generate(3)
            unchanged(before)
        end)
        with_model(f,carriers->begin
            API.mesh.generate(3);API.mesh.set_order(2)
            API.mesh.create_edges();API.mesh.create_faces()
            msh=first(API.mesh.get_element_types(3,carriers.volume))
            conn=API.mesh.get_elements_by_type(msh,carriers.volume)[2]
            primary=msh_spec(msh_type(msh_spec(msh).family,1)).nnodes
            point=sum(API.mesh.get_node(node)[1] for node in conn[1:primary])/primary
            @test API.mesh.get_element_by_coordinates(point...,3,true)[1]>0
            @test API.LAST_MESH_EDGES[]!==nothing && API.LAST_MESH_FACES[]!==nothing && API.LAST_MIXED_MESH_LOCATOR[]!==nothing
            API.mesh.set_transfinite_curve(f.curve_tags[last(f.long_pair)],7)
            before=snapshot();@test_throws r"QuadTriNoNewVerts" API.mesh.generate(3);unchanged(before)
        end)
        collapsed=Cert.fixture("collapsed";strip_length=5,layers=:three,height=.25,offset=(0.,0.,2.0^50+2.))
        with_model(collapsed,carriers->begin
            API.mesh.generate(0);API.mesh.set_order(2)
            before=snapshot();@test_throws r"QuadTriNoNewVerts" API.mesh.generate(3);unchanged(before)
        end)
        base=Cert.fixture("p2_precision";strip_length=5,layers=:one,offset=(0.,0.,2.0^50+2.))
        precision=merge(base,(;source=replace(base.source,"Layers{1}"=>"Layers{4}"),
            levels=[0.,.25,.5,.75,1.],intervals=4))
        with_model(precision,carriers->begin
            API.mesh.generate(3);expected=free_actual_carriers(carriers)
            @test Cert.certify(Cert.public_volume(API,carriers.volume),precision,expected.source).total==1
            before=snapshot();@test_throws ArgumentError API.mesh.set_order(2);unchanged(before)
        end)
    end
end
end
get(ENV,"TESSELLA_B4_FREE_API_DEFINE_ONLY","0")=="1" || run_tests()
end
