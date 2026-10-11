module APINativeRecombineAdapterTests
using Test,Tessella,SHA,Serialization,TOML
if Sys.iswindows() && Sys.ARCH===:x86_64
const A=Tessella.API
const G=Tessella.GmshLibm
call(function_value,args...;keywords...)=Base.invokelatest(function_value,args...;keywords...)
fingerprint(value)=begin io=IOBuffer();serialize(io,value);take!(io);end
rawbits(value)=reinterpret(UInt64,value)
parse_values(text,T)=parse.(T,split(text,','))
function start_plane(reposition,topology;stored=true)
    call(A.initialize)
    # Each pinned native fixture starts with a clear process Msg state.
    call(A.mesh.generate,0)
    call(A.option,"Mesh.Renumber",0)
    call(A.option,"Mesh.RecombineNodeRepositioning",reposition)
    call(A.option,"Mesh.RecombineOptimizeTopology",topology)
    corners=((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.))
    for (tag,point) in enumerate(corners)
        call(A.model.add_point,point...;tag,meshSize=.1)
    end
    for tag in 1:4;call(A.model.add_line,tag,mod1(tag+1,4);tag);end
    call(A.model.add_curve_loop,[1,2,3,4];tag=1)
    call(A.model.add_plane_surface,[1];tag=1)
    for (tag,point) in enumerate(corners)
        call(A.mesh.add_nodes,0,tag,[10000+tag],collect(point))
    end
    midpoint=((.5,0.,0.),(1.,.5,0.),(.5,1.,0.),(0.,.5,0.))
    for (curve,point) in enumerate(midpoint)
        parameter=call(A.model.get_parametrization,1,curve,collect(point))
        call(A.mesh.add_nodes,1,curve,[10004+curve],collect(point),parameter)
        n1=10000+curve;nm=10004+curve;n2=10000+mod1(curve+1,4)
        call(A.mesh.add_elements_by_type,curve,1,[20000+2curve-1,20000+2curve],[n1,nm,nm,n2])
    end
    parameter=stored ? call(A.model.get_parametrization,2,1,[.25,.3,0.]) : Float64[]
    call(A.mesh.add_nodes,2,1,[10009],[.25,.3,0.],parameter)
    cells=[10001,10005,10009,10008,10005,10002,10006,10009,
           10009,10006,10003,10007,10008,10009,10007,10004]
    call(A.mesh.add_elements_by_type,1,3,[30001,30002,30003,30004],cells)
    return cells
end

@testset "Actual public defaults, Plane XYZ/UV, labels and boundaries" begin
    @test bytes2hex(sha256(read(joinpath(@__DIR__,"../artifacts/native_recombine_reposition_primary.json"))))==
        "d6c4e81a9c58314a0e67d5d21280f43c6de762609d2138579327738d4558480b"
    for reposition in (0,1),topology in (0,5)
        cells=start_plane(reposition,topology)
        @test call(A.option,"Mesh.Smoothing")==1
        @test call(A.option,"Mesh.RecombineMinimumQuality")==.01
        before=[call(A.mesh.get_node,node) for node in 10001:10008]
        control=G._win_x87_control_word()
        call(A.mesh.recombine)
        expected=reposition==0 ? [.25,.3,0.] : topology==0 ? [.484375,.4875,0.] : [.49609375,.496875,0.]
        xyz,uv,dimension,entity=call(A.mesh.get_node,10009)
        @test rawbits.(xyz)==rawbits.(expected)
        @test rawbits.(uv)==rawbits.(expected[[2,1]])
        @test (dimension,entity)==(2,1)
        @test [call(A.mesh.get_node,node) for node in 10001:10008]==before
        types,tags,nodes=call(A.mesh.get_elements,2,1)
        @test types==[3] && only(tags)==[30001,30002,30003,30004]
        @test only(nodes)==cells
        @test call(A.mesh.get_max_element_tag)==30004
        @test call(A.mesh.get_max_node_tag)==10009
        @test G._win_x87_control_word()==control
        call(A.finalize)
    end
end

@testset "Actual public degree-two topology and deleted-node high-water" begin
    path=joinpath(@__DIR__,"../artifacts/native_recombine_degree_two_primary.tsv")
    @test bytes2hex(sha256(read(path)))=="46dc14373b2d83dc691513130bf3e051d7d7848ab05963a608b0ab839d2a552d"
    for line in readlines(path)[2:end]
        f=split(line,'\t');labels=parse_values(f[4],UInt64)
        coords=reshape(parse_values(f[2],Float64),3,length(labels))
        owners=[(Int(pair[1]),Int32(pair[2])) for pair in
            (parse_values(text,Int) for text in split(f[3],';'))]
        cells=hcat((parse_values(text,Int32) for text in split(f[5],';'))...)
        element_tags=parse_values(f[6],UInt64)
        call(A.initialize)
        # Each pinned native fixture starts with a clear process Msg state.
        call(A.mesh.generate,0)
        call(A.option,"Mesh.Renumber",0)
        @test call(A.option,"Mesh.RecombineNodeRepositioning")==1
        @test call(A.option,"Mesh.RecombineOptimizeTopology")==5
        @test call(A.option,"Mesh.Smoothing")==1
        call(A.model.add_discrete_entity,2,1)
        for (row,(dimension,entity)) in enumerate(owners)
            dimension==0 && call(A.model.add_discrete_entity,0,entity)
            call(A.mesh.add_nodes,dimension,entity,[labels[row]],coords[:,row])
        end
        call(A.mesh.add_elements_by_type,1,3,element_tags,vec(labels[cells]))
        before=Dict(labels[row]=>call(A.mesh.get_node,labels[row]) for row in eachindex(labels) if owners[row][1]==0)
        call(A.mesh.recombine)
        types,tags,nodes=call(A.mesh.get_elements,2,1)
        @test types==[3]
        raw=reshape(only(nodes),4,:)
        actual=sort!([(tag,Tuple(raw[:,column])) for (column,tag) in enumerate(only(tags))])
        expected=sort!([(entry[1],Tuple(entry[2:5])) for entry in
            (parse_values(text,UInt64) for text in split(f[7],';'))])
        @test actual==expected
        @test sort(first(call(A.mesh.get_nodes)))==sort(collect(keys(before)))
        @test all(call(A.mesh.get_node,node)==payload for (node,payload) in before)
        @test call(A.mesh.get_max_element_tag)==parse(UInt64,f[8])
        @test call(A.mesh.get_max_node_tag)==maximum(labels)
        for row in findall(owner->owner[1]==2,owners)
            @test_throws ErrorException call(A.mesh.get_node,labels[row])
        end
        call(A.finalize)
    end
end

@testset "Exact original Plane without raw Face UV relocates natively" begin
    raw=joinpath(@__DIR__,"../artifacts/native_recombine_adapter_missing_face_primary_v1.json")
    path=joinpath(@__DIR__,"../artifacts/native_recombine_adapter_missing_face_primary_v1.toml")
    @test bytes2hex(sha256(read(raw)))=="3f8739c17b0ef6c27361c5529b5ca459680ce9934935a3aac6873bfaebaa9e81"
    @test bytes2hex(sha256(read(path)))=="c1f2f65bb3c8673c3e53abad4bcb168c24b5dd1ba048ebb71c7b7fdf49293435"
    native=TOML.parsefile(path)
    @test native["gmsh_version"]=="4.15.2"
    @test native["primary_json_sha256"]==bytes2hex(sha256(read(raw)))
    cells=start_plane(1,5;stored=false)
    try
        function check_native_nodes(expected)
            for (label,record) in expected
                xyz,uv,dimension,entity=call(A.mesh.get_node,parse(UInt64,label))
                @test rawbits.(xyz)==rawbits.(Float64.(record["xyz"]))
                @test rawbits.(uv)==rawbits.(Float64.(record["uv"]))
                @test (dimension,entity)==(record["dimension"],record["entity"])
            end
        end
        check_native_nodes(native["before"]["nodes"])
        boundary=[call(A.mesh.get_node,node) for node in 10001:10008]
        elements=deepcopy(call(A.mesh.get_elements))
        old=(A.CURRENT[],A.LAST_MESH[],A.LAST_MESH_CLASS[])
        original=fingerprint(old)
        clocks=(A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[],A.NODE_COUNTER_WRAPPED[],A.ELEMENT_COUNTER_WRAPPED[])
        last_error=call(A.logger.get_last_error);errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
        word=G._win_x87_control_word()
        call(A.mesh.recombine)
        check_native_nodes(native["after"]["nodes"])
        @test isempty(call(A.mesh.get_node,10009)[2])
        @test [call(A.mesh.get_node,node) for node in 10001:10008]==boundary
        @test call(A.mesh.get_elements)==elements
        @test call(A.mesh.get_elements,2,1)==([3],[[30001,30002,30003,30004]],[cells])
        @test call(A.mesh.get_max_node_tag)==native["after"]["max_node"]==10009
        @test call(A.mesh.get_max_element_tag)==native["after"]["max_element"]==30004
        @test clocks==(A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[],A.NODE_COUNTER_WRAPPED[],A.ELEMENT_COUNTER_WRAPPED[])
        @test fingerprint(old)==original
        @test call(A.logger.get_last_error)==last_error
        @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==errors
        @test G._win_x87_control_word()==word
    finally
        call(A.finalize)
    end
end

@testset "Finite malformed stored Face UV rejects atomically" begin
    start_plane(1,5;stored=true)
    try
        # Deliberately corrupt retained runtime storage: a finite nonempty
        # one-vector is distinct from the newly supported absent parameter.
        call(A.mesh.get_nodes)
        A.CURRENT[].meshing.attached[(2,1)].node_native_params[UInt64(10009)]=[.3]
        @test call(A.mesh.get_node,10009)[2]==[.3]
        model=A.CURRENT[];cache=A.LAST_MESH[];class=A.LAST_MESH_CLASS[]
        before=fingerprint((model,cache,class,A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[],
            A.NODE_COUNTER_WRAPPED[],A.ELEMENT_COUNTER_WRAPPED[],
            call(A.mesh.get_nodes),call(A.mesh.get_elements)))
        last_error=call(A.logger.get_last_error);errors=A.NATIVE_LOOKUP_QUERY_ERRORS[]
        word=G._win_x87_control_word()
        caught=try call(A.mesh.recombine);nothing catch error;error end
        @test caught isa ArgumentError
        @test caught isa ArgumentError && occursin("two finite stored Surface parameters",sprint(showerror,caught))
        @test A.CURRENT[]===model && A.LAST_MESH[]===cache && A.LAST_MESH_CLASS[]===class
        @test fingerprint((model,cache,class,A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[],
            A.NODE_COUNTER_WRAPPED[],A.ELEMENT_COUNTER_WRAPPED[],
            call(A.mesh.get_nodes),call(A.mesh.get_elements)))==before
        @test call(A.logger.get_last_error)==last_error
        @test A.NATIVE_LOOKUP_QUERY_ERRORS[]==errors
        @test G._win_x87_control_word()==word
    finally
        call(A.finalize)
    end
end

function public_snapshot()
    types,tags,nodes=call(A.mesh.get_elements)
    node_tags=first(call(A.mesh.get_nodes))
    return (types=Int.(types),elements=Dict(Int(kind)=>(UInt64.(labels),UInt64.(connectivity))
        for (kind,labels,connectivity) in zip(types,tags,nodes)),
        nodes=Dict(UInt64(tag)=>call(A.mesh.get_node,tag) for tag in node_tags),
        max_node=call(A.mesh.get_max_node_tag),max_element=call(A.mesh.get_max_element_tag))
end

function canonical_cycle(nodes)
    best=Tuple(nodes)
    for start in 2:length(nodes)
        rotated=ntuple(offset->nodes[mod1(start+offset-1,length(nodes))],length(nodes))
        isless(rotated,best) && (best=rotated)
    end
    return best
end

@testset "Public weighted C6 matches held ordered native controls" begin
    path=joinpath(@__DIR__,"../artifacts/native_recombine_weighted_c6.tsv")
    @test bytes2hex(sha256(read(path)))=="87de72f99007c3018eb192a2d1cb1138fbd1649c48eb9f2d5ce5c419f22d867b"
    for line in readlines(path)[2:end]
        name,shift_text,reversed_text,records=split(line,'\t')
        xyz=name=="unique_minimum" ?
            Float64[1.05 .54 -.48 -1.10 -.49 .56 .02;0 .95 .81 -.03 -.90 -.89 .03;0 0 0 0 0 0 0] :
            Float64[1 .5 -.5 -1 -.5 .5 0;0 1 1 0 -1 -1 0;0 0 0 0 0 0 0]
        labels=UInt64.(10003:7:10045)
        reversed_text=="1" && reverse!(labels)
        triangles=Int32[(row==1 ? 7 : row==2 ? column : mod1(column+1,6)) for row in 1:3,column in 1:6]
        shift=parse(Int,shift_text)
        triangles=triangles[:,vcat(shift+1:6,1:shift)]
        call(A.initialize)
        # Each pinned native fixture starts with a clear process Msg state.
        call(A.mesh.generate,0)
        try
            call(A.option,"Mesh.Renumber",0)
            call(A.option,"Mesh.RecombinationAlgorithm",1)
            call(A.option,"Mesh.RecombineNodeRepositioning",0)
            call(A.option,"Mesh.RecombineOptimizeTopology",0)
            call(A.option,"Mesh.RecombineMinimumQuality",0)
            call(A.model.add_discrete_entity,2,1)
            call(A.mesh.add_nodes,2,1,labels,vec(xyz))
            call(A.mesh.add_elements_by_type,1,2,30001:30006,vec(labels[triangles]))
            call(A.mesh.set_recombine,2,1,0)
            before=call(A.mesh.get_nodes)
            word=G._win_x87_control_word()
            call(A.mesh.recombine)
            @test call(A.mesh.get_nodes)==before
            types,tags,nodes=call(A.mesh.get_elements,2,1)
            @test types==[3]
            inverse=Dict(tag=>row for (row,tag) in enumerate(labels))
            cells=reshape(only(nodes),4,:)
            actual=sort!([(tag,canonical_cycle([inverse[node] for node in cells[:,column]]))
                for (column,tag) in enumerate(only(tags))])
            expected=sort!([(parse(UInt64,first(split(record,':'))),
                Tuple(parse.(Int,split(last(split(record,':')),',')))) for record in split(records,';')])
            @test actual==expected
            @test call(A.mesh.get_max_element_tag)==30015
            @test G._win_x87_control_word()==word
        finally
            call(A.finalize)
        end
    end
end
function check_snapshot(actual,expected)
    @test actual.types==expected["types"]
    @test actual.max_node==parse(UInt64,expected["max_node"])
    @test actual.max_element==parse(UInt64,expected["max_element"])
    @test Set(keys(actual.nodes))==Set(parse(UInt64,row["tag"]) for row in expected["nodes"])
    for row in expected["elements"]
        @test actual.elements[row["msh"]]==(parse.(UInt64,row["tags"]),parse.(UInt64,row["nodes"]))
    end
    for row in expected["nodes"]
        xyz,parameters,dimension,entity=actual.nodes[parse(UInt64,row["tag"])]
        @test rawbits.(xyz)==rawbits.(Float64.(row["xyz"]))
        @test rawbits.(parameters)==rawbits.(Float64.(row["parameters"]))
        @test (dimension,entity)==(row["dimension"],row["entity"])
    end
end

@testset "Actual generated classified Surface uses retained native parameters" begin
    path=joinpath(@__DIR__,"../artifacts/native_recombine_generated_primary.toml")
    cases=TOML.parsefile(path)
    @test bytes2hex(sha256(read(joinpath(@__DIR__,"../artifacts/native_recombine_generated_primary.json"))))==
        cases["primary_json_sha256"]
    for case in cases["cases"]
        call(A.initialize)
        # Each pinned native fixture starts with a clear process Msg state.
        call(A.mesh.generate,0)
        try
            mktempdir() do directory
                source=joinpath(directory,"square.geo")
                write(source,raw"""
Point(1)={0,0,0,0.1};Point(2)={1,0,0,0.1};Point(3)={1,1,0,0.1};Point(4)={0,1,0,0.1};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
Transfinite Curve{1,2,3,4}=3;Transfinite Surface{1};Recombine Surface{1};
""")
                call(A.open_geo!,source;mesh_dim=0)
            end
            call(A.option,"Mesh.Renumber",0)
            call(A.option,"Mesh.RecombinationAlgorithm",1)
            call(A.option,"Mesh.RecombineNodeRepositioning",case["reposition"])
            call(A.option,"Mesh.RecombineOptimizeTopology",case["topology"])
            call(A.option,"Mesh.Smoothing",1)
            call(A.option,"Mesh.RecombineMinimumQuality",.01)
            call(A.mesh.generate,2)
            @test A.LAST_MESH[] isa Tessella.Elements.MixedMesh
            check_snapshot(public_snapshot(),case["initial"])
            center=UInt64(case["center"])
            uv=call(A.model.get_parametrization,2,1,[.25,.3,0.])
            call(A.mesh.set_node,center,[.25,.3,0.],uv)
            check_snapshot(public_snapshot(),case["before"])
            old_model=A.CURRENT[];old_mesh=A.LAST_MESH[]
            before=fingerprint((old_model,old_mesh))
            word=G._win_x87_control_word()
            call(A.mesh.recombine)
            check_snapshot(public_snapshot(),case["after"])
            @test fingerprint((old_model,old_mesh))==before
            @test G._win_x87_control_word()==word
        finally
            call(A.finalize)
        end
    end
end

function classified_payload_fixture(;consumed=false,extra_curve=false,extra_surface=false)
    E=Tessella.Elements
    call(A.initialize)
    # Each pinned native fixture starts with a clear process Msg state.
    call(A.mesh.generate,0)
    call(A.option,"Mesh.Renumber",0)
    for entity in ((0,99),(1,1),(1,2),(2,1),(3,1))
        call(A.model.add_discrete_entity,entity...)
    end
    given_xyz=[2.,3.,0.,3.,3.,0.]
    given_nodes=UInt64[91001,91002];given_cells=UInt64[60001]
    if extra_curve
        call(A.model.add_discrete_entity,1,3)
        call(A.mesh.add_nodes,1,3,given_nodes,given_xyz,[.2,.8])
        call(A.mesh.add_elements_by_type,3,1,given_cells,given_nodes)
    end
    surface_xyz=[5.,0.,0.,6.,0.,0.,6.,1.,0.,5.,1.,0.]
    surface_nodes=UInt64[92001,92002,92003,92004]
    if extra_surface
        call(A.model.add_discrete_entity,2,2)
        call(A.mesh.add_nodes,2,2,surface_nodes,surface_xyz)
        call(A.mesh.add_elements_by_type,2,3,[61001],surface_nodes)
    end
    xyz=Float64[0 1 1 0 10 12 13 12 12 0 1 0 1;
                0 0 1 1 10 0 0 1 0 3 3 4 4;
                0 0 0 0 10 0 0 0 1 0 0 0 0]
    entities=Dict((2,1)=>E.MixedEntity(2,1,(0.,0.,0.,1.,1.,0.);physical_tags=[7]),
        (0,99)=>E.MixedEntity(0,99,(10.,10.,10.)),
        (3,1)=>E.MixedEntity(3,1,(12.,0.,0.,13.,1.,1.);physical_tags=[91]),
        (1,1)=>E.MixedEntity(1,1,(0.,3.,0.,1.,3.,0.);physical_tags=[42]),
        (1,2)=>E.MixedEntity(1,2,(0.,4.,0.,1.,4.,0.);physical_tags=[43]))
    blocks=[E.ElementBlock(2,Int32[1 1;2 3;3 4],Int32[7,7]),
        E.ElementBlock(4,reshape(Int32[6,7,8,9],4,1),Int32[91]),
        E.ElementBlock(1,Int32[10 12;11 13],Int32[42,43])]
    owners=[Int32[1,1],Int32[1],Int32[1,2]]
    tags=[UInt64[30001,30002],UInt64[50001],UInt64[20001,20002]]
    node_tags=UInt64[10003,10010,10017,10024,70001,80001,80002,80003,80004,90001,90002,90003,90004]
    node_owners=vcat(fill((2,Int32(1)),4),[(0,Int32(99))],fill((3,Int32(1)),4),
        [(1,Int32(1)),(1,Int32(1)),(1,Int32(2)),(1,Int32(2))])
    projected=Union{Nothing,Vector{Float64}}[nothing for _ in node_tags]
    projected[5]=Float64[]
    for row in 10:13;projected[row]=[isodd(row) ? .9 : .1];end
    data=E.MixedEntityData(entities;node_entities=node_owners,node_parametric=projected,
        external_node_tags=node_tags,block_entities=owners,external_element_tags=tags)
    sections=[E.MshDataSection("NodeData",["temperature"],[0.],Int64[0,1,1],
        Int32[],Int32[6],Int32[],false,[.125],Int32(4)),
        E.MshDataSection("ElementData",["volume data"],[0.],Int64[0,1,1],
        Int32[3],Int32[],Int32[],false,[.375],Int32(4)),
        E.MshDataSection("ElementNodeData",["volume supports"],[0.],Int64[0,1,1],
        Int32[3],Int32[6,7,8,9],Int32[4],false,[1.,2.,3.,4.],Int32(4))]
    consumed && push!(sections,E.MshDataSection("ElementData",["consumed Surface"],[0.],Int64[0,1,1],
        Int32[1],Int32[],Int32[],false,[.5],Int32(4)))
    partition=E.MixedPartitionData(2;ghost_elements=[E.MixedGhostElement(3,1,[2])])
    affine=vec(permutedims(Float64[1 0 0 0;0 1 0 1;0 0 1 0;0 0 0 1]))
    link=E.MixedPeriodicLink(1,2,1,Int32[12,13],Int32[10,11];affine)
    mesh=E.MixedMesh(xyz,blocks;entity_data=data,elementary_entities=owners,
        physical_names=Dict((2,7)=>"material",(3,91)=>"solid"),periodic_links=[link],
        ancillary_sections=[E.MshAncillarySection("Comments",false,collect(codeunits("retained bytes\n")),Int32(4))],
        data_sections=sections,partition_data=partition)
    class=call(A._higher_classification_from_data,mesh,"native recombine payload fixture")
    runtime=deepcopy(projected)
    for row in 10:13;runtime[row]=[.01row,.02row];end
    class=call(A._mixed_rebind_class,class,mesh;node_parametric=runtime)
    lock(A.STATE_LOCK) do
        A._replace_mesh_cache_locked!(mesh,class)
    end
    return (;mesh,given_xyz,given_nodes,given_cells,surface_xyz,surface_nodes)
end

@testset "Public classified cache preserves runtime, projected and MSH payloads" begin
    fixture=classified_payload_fixture(;extra_curve=true,extra_surface=true)
    source=fixture.mesh
    try
        # The cache does not cover the already published Curve3 record.
        given_xyz=fixture.given_xyz;given_nodes=fixture.given_nodes;given_cells=fixture.given_cells
        @test length(A.LAST_MESH[].data_sections)==3
        @test !((1,Int32(3)) in A.LAST_MESH_CLASS[].public_tags.authority)
        @test !((2,Int32(2)) in A.LAST_MESH_CLASS[].public_tags.authority)
        old_model=A.CURRENT[];old_cache=A.LAST_MESH[]
        before=fingerprint((old_model,old_cache,source,given_xyz,given_nodes,given_cells,
            fixture.surface_xyz,fixture.surface_nodes))
        lower=Dict(tag=>call(A.mesh.get_node,tag) for tag in source.entity_data.external_node_tags[5:end])
        other_surface=Dict(tag=>call(A.mesh.get_node,tag) for tag in fixture.surface_nodes)
        word=G._win_x87_control_word()
        call(A.mesh.recombine)
        result=A.LAST_MESH[];data=result.entity_data
        @test fingerprint((old_model,old_cache,source,given_xyz,given_nodes,given_cells,
            fixture.surface_xyz,fixture.surface_nodes))==before
        @test all(call(A.mesh.get_node,tag)==value for (tag,value) in lower)
        @test all(call(A.mesh.get_node,tag)==value for (tag,value) in other_surface)
        @test G._win_x87_control_word()==word
        @test call(A.mesh.get_max_element_tag)==61003
        @test call(A.mesh.get_elements,2,1)[1]==[3]
        @test only(call(A.mesh.get_elements,2,1)[2])==UInt64[61003]
        @test sort(only(call(A.mesh.get_elements,2,1)[3]))==sort(source.entity_data.external_node_tags[1:4])
        @test call(A.mesh.get_elements,1,3)==([1],[UInt64[60001]],[UInt64[91001,91002]])
        @test call(A.mesh.get_elements,2,2)==([3],[UInt64[61001]],[UInt64[92001,92002,92003,92004]])
        @test source.coords!==result.coords
        @test result.data_sections[1].values==[.125]
        @test data.external_node_tags[result.data_sections[1].nodes]==UInt64[80001]
        volume=only(findall(block->block.msh==4,result.blocks))
        volume_flat=sum(size(result.blocks[index].nodes,2) for index in 1:volume-1;init=0)+1
        @test result.data_sections[2].elements==Int32[volume_flat]
        @test result.data_sections[3].elements==Int32[volume_flat]
        @test data.external_node_tags[result.data_sections[3].nodes]==UInt64[80001,80002,80003,80004]
        @test result.data_sections[3].values==[1.,2.,3.,4.]
        @test only(result.partition_data.ghost_elements).element==volume_flat
        @test only(result.partition_data.ghost_elements).ghost_partitions==Int32[2]
        @test result.ancillary_sections[1].payload==collect(codeunits("retained bytes\n"))
        @test result.physical_names==source.physical_names
        @test result.blocks[volume].tags==Int32[91]
        for uid in UInt64[20001,20002]
            positions=[(bi,column) for (bi,labels) in enumerate(data.external_element_tags)
                for (column,tag) in enumerate(labels) if tag==uid]
            bi,column=only(positions)
            @test result.blocks[bi].tags[column]==(uid==20001 ? 42 : 43)
        end
        for tag in UInt64[90001,90002,90003,90004]
            old=findfirst(==(tag),source.entity_data.external_node_tags)
            row=findfirst(==(tag),data.external_node_tags)
            @test data.node_parametric[row]==source.entity_data.node_parametric[old]
            @test call(A.mesh.get_node,tag)[2]==[.01old,.02old]
        end
        link=only(result.periodic_links)
        @test data.external_node_tags[link.slave_nodes]==UInt64[90003,90004]
        @test data.external_node_tags[link.master_nodes]==UInt64[90001,90002]
        @test link.affine==only(source.periodic_links).affine
        result.data_sections[1].values[1]=9.
        result.ancillary_sections[1].payload[1]=0x00
        data.entities[(2,1)].physical_tags[1]=8
        @test source.data_sections[1].values==[.125]
        @test source.ancillary_sections[1].payload==collect(codeunits("retained bytes\n"))
        @test source.entity_data.entities[(2,1)].physical_tags==Int32[7]
        given_xyz[1]=99.;given_nodes[1]=99999;given_cells[1]=99999
        @test call(A.mesh.get_node,91001)[1]==[2.,3.,0.]
        @test only(call(A.mesh.get_elements,1,3)[2])==UInt64[60001]
    finally
        call(A.finalize)
    end
end

@testset "Blocked cell payload publication is atomic" begin
    source=classified_payload_fixture(;consumed=true).mesh
    try
        model=A.CURRENT[];cache=A.LAST_MESH[];class=A.LAST_MESH_CLASS[]
        before=fingerprint((model,cache,class,public_snapshot(),A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))
        failure=try call(A.mesh.recombine);nothing catch exception;exception end
        @test failure isa ArgumentError
        @test occursin("view-data transfer",sprint(showerror,failure))
        @test A.CURRENT[]===model && A.LAST_MESH[]===cache && A.LAST_MESH_CLASS[]===class
        @test fingerprint((model,cache,class,public_snapshot(),A.NODE_TAG_MAX[],A.ELEMENT_TAG_MAX[]))==before
        @test source===cache
    finally
        call(A.finalize)
    end
end

@testset "Empty parameter storage shapes remain empty" begin
    for matrix_columns in (0,1)
        start_plane(0,0)
        try
            record=A.CURRENT[].meshing.attached[(0,1)]
            record.node_params=zeros(0,matrix_columns)
            before=fingerprint(record);node=call(A.mesh.get_node,10001)
            call(A.mesh.recombine)
            @test size(A.CURRENT[].meshing.attached[(0,1)].node_params)==(0,matrix_columns)
            @test fingerprint(record)==before
            @test call(A.mesh.get_node,10001)==node
        finally
            call(A.finalize)
        end
    end
    for matrix_columns in (0,4)
        call(A.initialize)
        # Each pinned native fixture starts with a clear process Msg state.
        call(A.mesh.generate,0)
        call(A.option,"Mesh.Renumber",0)
        try
            call(A.model.add_discrete_entity,2,1)
            given_xyz=[0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.]
            given_nodes=UInt64[101,102,103,104]
            call(A.mesh.add_nodes,2,1,given_nodes,given_xyz)
            call(A.mesh.add_elements_by_type,1,3,[301],given_nodes)
            record=A.CURRENT[].discrete[(2,1)]
            record.node_params=zeros(0,matrix_columns)
            before=fingerprint(record);snapshot=public_snapshot()
            call(A.mesh.recombine)
            @test size(A.CURRENT[].discrete[(2,1)].node_params)==(0,matrix_columns)
            @test fingerprint(record)==before
            @test public_snapshot()==snapshot
            @test given_nodes==UInt64[101,102,103,104]
            @test given_xyz==[0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.]
        finally
            call(A.finalize)
        end
    end
end
end # Windows x86_64 SDK context
end # module
