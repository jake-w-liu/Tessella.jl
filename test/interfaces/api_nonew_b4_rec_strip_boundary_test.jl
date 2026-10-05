if !isdefined(@__MODULE__,:QuadTriNoNewTriangleCertificates)
    include("../geometry/quadtri_nonew_triangle_certificates.jl")
end
if !isdefined(@__MODULE__,:QuadTriNoNewB4RecStripCertificates)
    include("../geometry/quadtri_nonew_b4_rec_strip_certificates.jl")
end
module NoNewB4RecStripBoundaryTests
using Test,Tessella,SHA
using ..QuadTriNoNewB4RecStripCertificates
using ..QuadTriNoNewTriangleCertificates
const API=Tessella.API
const Cert=QuadTriNoNewB4RecStripCertificates
const LENGTHS=(5,7,9)
const PROFILES=(:one,:nonbinary)
expected_primary(m,n)=2*(m+1)*(n+1)+m
expected_quadratic(m,n)=(12m+6)*n+14m+3
expected_owners(m,n,order)=order==1 ? [8,4n+4m-8,(2m-2)*(n-1),m] :
    [8,8n+8m-4,8m*n-2,(2m-1)*(2n-1)+8m]
expected_cap(m)=(2m-1,6m+3)
expected_long(m,n)=((2m-1)*(2n-1),(4m+2)*n+2m+1)
expected_short(n)=(2n-1,6n+3)

function with_model(f,action)
    API.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"b4_boundary.geo");write(path,f.source)
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
        visibility=copy(API.ELEMENT_VISIBILITY[]),locator=API.LAST_MESH_LOCATOR[])
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
end

function node_geometry(volume=nothing)
    tags=volume===nothing ? API.mesh.get_nodes()[1] : API.mesh.get_nodes(3,volume,true)[1]
    Set((Tuple(p),dim,owner,Tuple(uv)) for node in tags
        for (p,uv,dim,owner) in (API.mesh.get_node(node),))
end

function cell_geometry(volume=nothing)
    types,tags,connectivity=volume===nothing ? API.mesh.get_elements() : API.mesh.get_elements(3,volume)
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
    buffer=IOBuffer();println(buffer,"b4_rec_strip_support_query_v1");println(buffer,base.sha)
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


include("api_nonew_b4_rec_strip_actual_supports.jl")

artifact_name(m,profile,direction,order)="b4_M$(m)_$(profile)_D$(direction)_P$(order)"

function artifact_records()
    records=NamedTuple[]
    for m in LENGTHS,profile in PROFILES,direction in (-1,1)
        f=Cert.fixture("artifact";strip_length=m,layers=profile,height=direction,pins=(1,2,3,4))
        with_model(f,carriers->begin
            API.mesh.generate(3);check_owner_totals(carriers.volume,m,f.intervals,1)
            expected=b4_actual_carriers(carriers)
            linear=Cert.certify(Cert.public_volume(API,carriers.volume),f,expected.source)
            @test linear.total==1 && length(linear.centers)==m
            push!(records,artifact_record(artifact_name(m,profile,direction,1),1))
            API.mesh.set_order(2);b4_check_actual_p2(carriers,f,expected)
            push!(records,artifact_record(artifact_name(m,profile,direction,2),2))
        end)
    end
    return records
end

function artifact_pins()
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_b4_rec_strip_crc.txt")
    result=Dict{String,NamedTuple}()
    for line in readlines(path)
        isempty(line) || startswith(line,"#") || begin
            fields=split(line,'\t');length(fields)==6 || error("B4 artifact row width")
            name,kind,n_nodes,n_cells,families,sha=String.(fields)
            haskey(result,name) && error("duplicate B4 artifact name")
            result[name]=(;name,kind,n_nodes=parse(Int,n_nodes),n_cells=parse(Int,n_cells),families,sha)
        end
    end
    length(result)==24 || error("B4 artifact must contain 24 products")
    return result
end

function lower_primary_tags()
    primary=Set{UInt64}()
    types,tags,connectivity=API.mesh.get_elements()
    for (msh,nodes) in zip(types,connectivity)
        msh_dimension(msh)<=1 || continue
        width=Tessella.Elements.msh_spec(msh).nnodes
        count=msh==15 ? 1 : 2
        for cell in eachcol(reshape(nodes,width,:)),node in cell[1:count]
            push!(primary,node)
        end
    end
    return primary
end

function check_lower_regeneration(entity,old,current,order,old_primary,new_primary)
    oldtags,oldxyz,oldparams=old
    tags,xyz,params=current
    oldpoints=reshape(oldxyz,3,:);points=reshape(xyz,3,:)
    olduv=reshape(oldparams,entity[1],length(oldtags));uv=reshape(params,entity[1],length(tags))
    @test length(tags)==length(oldtags) && allunique(tags)
    @test Set(Tuple(p) for p in eachcol(uv))==Set(Tuple(p) for p in eachcol(olduv))
    signature(a,b,ids,selected)=Set((Tuple(reinterpret(UInt64,a[:,i])),Tuple(b[:,i])) for i in eachindex(ids) if ids[i] in selected)
    @test signature(points,uv,tags,new_primary)==signature(oldpoints,olduv,oldtags,old_primary)
    if order==2
        # Gmsh Generator.cpp runs SetOrder1 before MeshOnlyEmpty and creates
        # requested P2 supports on CAD afterward. Primary sampling is retained;
        # old linear-average P2 bits need not equal new CAD evaluation bits.
        # A pinned edited-midpoint control independently proves this rebuild.
        @test xyz==API.model.get_value(entity[1],entity[2],params)
    else
        @test signature(points,uv,tags,Set(tags))==signature(oldpoints,olduv,oldtags,Set(oldtags))
    end
end

function check_owner_totals(volume,m,n,order)
    own=zeros(Int,4)
    for node in API.mesh.get_nodes(3,volume,true)[1]
        own[API.mesh.get_node(node)[3]+1]+=1
    end
    @test own==expected_owners(m,n,order)
    @test sum(own)==(order==1 ? expected_primary(m,n) : expected_quadratic(m,n))
end

function run_tests()
@testset "Dynamic recombined B4 actual API contracts" begin
    @testset "Twelve healthy public P1/P2 constructions" begin
        b4_native_construction_cases()
    end
    @testset "Sparse CAD identities and actual captured recipes" begin
        for m in LENGTHS,plane in (:XY,:YZ,:ZX),direction in (-1,1)
            f=Cert.fixture("sparse";strip_length=m,plane,height=direction,layers=:nonbinary,
                shape=:rounded,tags=:sparse,pins=:opposite,curve_reverse=5,law="Progression",coefficient=4.)
            with_model(f,carriers->begin
                API.mesh.generate(3);expected=b4_actual_carriers(carriers)
                @test length(Cert.certify(Cert.public_volume(API,carriers.volume),f,expected.source).centers)==m
                API.mesh.set_order(2);b4_check_actual_p2(carriers,f,expected)
            end)
        end
    end
    @testset "Actual geometry survives permutations and exact selections" begin
        f=Cert.fixture("permutations";strip_length=5,layers=:nonbinary,tags=:sparse)
        with_model(f,carriers->begin
            API.mesh.generate(3);API.mesh.set_order(2)
            geometry=node_geometry();cells=cell_geometry()
            old=API.mesh.get_nodes()[1];API.mesh.renumber_nodes(old,reverse(old))
            @test node_geometry()==geometry && cell_geometry()==cells
            groups=API.mesh.get_elements()[2]
            element=UInt64.(vcat(groups...))
            before_elements=Dict(tag=>API.mesh.get_element(tag) for tag in element)
            source=Tessella.Model.mesh_model_surface(API.CURRENT[],carriers.source)
            # Native public labels rename cells simultaneously across families;
            # ordered connectivity and actual entity identity stay with the cell.
            # Independent Gmsh 4.15.2 cross-family label evidence is recorded in
            # b4_cross_family_renumber_primary_v1.json (D46F930B...ACB87D29).
            @test any(before_elements[old][1]!=before_elements[new][1]
                for (old,new) in zip(element,reverse(element)))
            API.mesh.renumber_elements(element,reverse(element))
            @test all(API.mesh.get_element(new)==before_elements[old]
                for (old,new) in zip(element,reverse(element)))
            public_elements=UInt64.(vcat(API.mesh.get_elements()[2]...))
            @test allunique(public_elements) && Set(public_elements)==Set(element)
            @test node_geometry()==geometry && cell_geometry()==cells
            certified=Cert.certify_quadratic(Cert.public_volume(API,carriers.volume),f,source)
            @test certified.p1_total==1 && length(certified.centers)==5
            before=snapshot()
            @test_throws r"duplicate element tag" API.mesh.renumber_elements(
                element[1:2],fill(element[1],2))
            unchanged(before)
            # Query the current family memberships after the global permutation.
            for group in API.mesh.get_elements()[2]
                API.mesh.renumber_elements(group,reverse(group))
            end
            @test node_geometry()==geometry && cell_geometry()==cells
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
    @testset "OnlyEmpty native lower regeneration and sparse public labels" begin
        for order in (1,2)
            f=Cert.fixture("bridge";strip_length=5,layers=:nonbinary)
            with_model(f,carriers->begin
                API.mesh.generate(3)
                primary=Set(API.mesh.get_nodes()[1]);order==2 && API.mesh.set_order(2)
                before=Dict((dim,tag)=>API.mesh.get_nodes(dim,tag,false,true)
                    for dim in (0,1) for (_,tag) in API.model.get_entities(dim))
                history=(API.mesh.get_max_node_tag(),API.mesh.get_max_element_tag())
                API.option("Mesh.ElementOrder",order);API.option("Mesh.MeshOnlyEmpty",1)
                API.mesh.generate(1)
                new_primary=lower_primary_tags()
                # The native generation prepass removes higher-dimensional
                # cells and rebuilds P2 Curve supports on CAD. Compare retained
                # primary bits by actual carrier and parameter, allowing the
                # generation renumber option to change public labels.
                for (entity,payload) in before
                    check_lower_regeneration(entity,payload,
                        API.mesh.get_nodes(entity...,false,true),order,primary,new_primary)
                end
                @test API.mesh.get_max_node_tag()>=history[1] && API.mesh.get_max_element_tag()>=history[2]
                before_geometry=node_geometry();before_cells=cell_geometry()
                tags=API.mesh.get_nodes()[1];labels=UInt64[10000+3i for i in eachindex(tags)]
                API.mesh.renumber_nodes(tags,labels)
                @test Set(API.mesh.get_nodes()[1])==Set(labels)
                @test node_geometry()==before_geometry && cell_geometry()==before_cells
            end)
        end
    end
    @testset "Real refined child support carriers" begin
        for m in (5,7),initial_order in (1,2)
            f=Cert.fixture("refine";strip_length=m,layers=:one,pins=(1,2,3,4))
            with_model(f,carriers->begin
                API.mesh.generate(3);initial_order==2 && API.mesh.set_order(2)
                API.mesh.refine();types,tags,_=API.mesh.get_elements(3,carriers.volume)
                # Actual fixed rec parents:2M Tet and5M Pyr. Native family
                # subdivision is8Tet/Tet and8Tet+4Pyr/Pyr, independently pinned.
                @test Dict(Int(msh)=>length(ids) for (msh,ids) in zip(types,tags))==Dict(4=>56m,7=>20m)
                API.mesh.set_order(2);supports=b4_public_support_nodes(carriers.volume)
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
                    slot=findfirst(==(curve),f.curve_tags);segments=slot in f.long_pair ? m : 1
                    own=b4_check_parameters(2,surface,false,3(4segments-1))
                    closure=b4_check_parameters(2,surface,true,5(4segments+1));@test own⊆closure
                end
                for surface in (carriers.source,carriers.top)
                    b4_check_parameters(2,surface,false,3(4m-1));b4_check_parameters(2,surface,true,5(4m+1))
                end
            end)
        end
    end
    @testset "Classified and actual quadratic file encodings" begin
        f=Cert.fixture("files";strip_length=7,layers=:nonbinary,tags=:sparse)
        text=f.source*"\nPhysical Curve(\"source lines\",501)={$(join(f.curve_tags,','))};\n"*
            "Physical Surface(\"caps\",502)={$(f.surface),sweep[0]};\n"*
            "Physical Surface(\"walls\",504)={sweep[2],sweep[3],sweep[4],sweep[5]};\n"*
            "Physical Volume(\"body\",503)={sweep[1]};\n"
        f=merge(f,(;source=text))
        with_model(f,carriers->begin
            @test Set(API.model.get_physical_groups())==Set([(1,501),(2,502),(2,504),(3,503)])
            @test API.model.get_entities_for_physical_group(3,503)==[carriers.volume]
            @test Set(API.model.get_entities_for_physical_group(2,502))==Set((carriers.source,carriers.top))
            for (dim,tag,name) in ((1,501,"source lines"),(2,502,"caps"),(2,504,"walls"),(3,503,"body"))
                @test API.model.get_physical_name(dim,tag)==name
            end
            API.mesh.generate(3)
            linear=Tessella.Model.mesh_model_volume(API.CURRENT[],carriers.volume)
            projected=Tessella.Model.model_to_mixed(API.CURRENT[],linear,3,carriers.volume)
            @test projected.physical_names==Dict((1,501)=>"source lines",(2,502)=>"caps",(2,504)=>"walls",(3,503)=>"body")
            @test all(all(==(503),block.tags) for block in projected.blocks if Tessella.Elements.msh_dimension(block.msh)==3)
            API.mesh.set_order(2)
            mktempdir() do directory
                check_files(projected,directory,"classified_b4")
                check_files(API.mesh.get(),directory,"actual_quadratic_b4")
            end
        end)
    end
    @testset "Independent coincident region identity and selective clear" begin
        left=Cert.fixture("left";strip_length=5,layers=:one)
        right=Cert.fixture("right";strip_length=5,layers=:one,tags=:sparse,offset=(3.,0.,0.))
        source=replace(left.source,"sweep[]"=>"left[]")*replace(right.source,"sweep[]"=>"right[]")
        API.initialize()
        try
            mktempdir() do directory
                path=joinpath(directory,"independent_b4.geo");write(path,source)
                opened=API.open_geo!(path;mesh_dim=0)
                for (_,tag) in API.model.get_entities(0)
                    p=API.model.get_value(0,tag,Float64[]);p[1]>=3 || continue
                    API.model.set_coordinates(tag,p[1]-3,p[2],p[3])
                end
                output=[Int.(opened.lists[name]) for name in ("left","right")]
                carriers=[(;source=f.surface,top=out[1],volume=out[2],laterals=out[3:end])
                    for (f,out) in zip((left,right),output)]
                API.mesh.generate(3);API.mesh.set_order(2)
                sets=[Set(API.mesh.get_nodes(3,c.volume,true)[1]) for c in carriers]
                @test isempty(intersect(sets...)) && length(API.mesh.get_nodes()[1])==2expected_quadratic(5,1)
                before=(node_geometry(carriers[2].volume),cell_geometry(carriers[2].volume))
                API.mesh.clear([(3,carriers[1].volume)])
                @test isempty(API.mesh.get_elements(3,carriers[1].volume)[1])
                @test (node_geometry(carriers[2].volume),cell_geometry(carriers[2].volume))==before
            end
        finally
            API.finalize()
        end
    end
    @testset "Unsupported and precision requests are atomic" begin
        for m in (5,7)
            f=Cert.fixture("atomic";strip_length=m,layers=:one)
            long_tags=Tuple(f.curve_tags[i] for i in f.long_pair)
            malformed=replace(replace(f.source," RecombLaterals"=>""),
                "Transfinite Curve{$(join(long_tags,','))}=$(m+1);"=>
                "Transfinite Curve{$(first(long_tags))}=$(m+1);Transfinite Curve{$(last(long_tags))}=$(m+2);")
            invalid=(malformed,
                replace(f.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{0.125,0.0,1.0}"),
                replace(f.source,"Layers{1}"=>"Layers{{1},{0.5}}"))
            for source in invalid
                @test source!=f.source
                with_model(merge(f,(;source)),carriers->begin
                    before=snapshot();@test_throws r"QuadTriNoNewVerts" API.mesh.generate(3);unchanged(before)
                end)
            end
        end
        f=Cert.fixture("collapsed";strip_length=5,layers=:three,height=.25,offset=(0.,0.,2.0^50+2.))
        with_model(f,carriers->begin
            before=snapshot();@test_throws r"QuadTriNoNewVerts" API.mesh.generate(3);unchanged(before)
        end)
        # Each last fan mean is represented for N2, but the requested elevated
        # edge nodes can collapse at a large normal offset. Certify P1 first.
        f=Cert.fixture("p2_precision";strip_length=5,layers=:one,height=1.,offset=(0.,0.,2.0^50+2.))
        f=merge(f,(;source=replace(f.source,"Layers{1}"=>"Layers{2}"),levels=[0.,.5,1.],intervals=2))
        with_model(f,carriers->begin
            API.mesh.generate(3);expected=b4_actual_carriers(carriers)
            @test Cert.certify(Cert.public_volume(API,carriers.volume),f,expected.source).total==1
            before=snapshot();@test_throws ArgumentError API.mesh.set_order(2);unchanged(before)
        end)
    end
end
end
get(ENV,"TESSELLA_B4_API_DEFINE_ONLY","0")=="1" || run_tests()
end
