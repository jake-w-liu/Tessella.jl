using Test
using Tessella

const _API01=Tessella.API

function _api01_begin(;renumber=true,order=1,only_empty=false)
    _API01.initialize()
    for (option,value) in (("Mesh.Renumber",renumber ? 1 : 0),
                           ("Mesh.ElementOrder",order),
                           ("Mesh.MeshOnlyEmpty",only_empty ? 1 : 0))
        _API01.option(option,value)
    end
end

function _api01_line(;transfinite=3)
    _API01.model.add_point(0.,0,0;tag=1,meshSize=0.2)
    _API01.model.add_point(1.,0,0;tag=2,meshSize=0.2)
    _API01.model.add_line(1,2;tag=1)
    transfinite===nothing || _API01.mesh.set_transfinite_curve(1,transfinite)
end

_api01_snapshot()=(_API01.mesh.get_nodes(),_API01.mesh.get_elements())

function _api01_check_references(;allow_duplicate_cells=false)
    tags,flat,_=_API01.mesh.get_nodes()
    @test allunique(tags)
    @test length(flat)==3length(tags)
    coordinates=Dict(tag=>flat[3i-2:3i] for (i,tag) in enumerate(tags))
    types,element_tags,connectivity=_API01.mesh.get_elements()
    @test issorted(types)
    @test allunique(vcat(element_tags...))
    for (msh,els,nodes) in zip(types,element_tags,connectivity)
        width=Tessella.Elements.msh_num_nodes(msh)
        @test length(nodes)==width*length(els)
        @test Set(nodes)⊆Set(tags)
        for (i,element) in enumerate(els)
            actual,cell,dim,owner=_API01.mesh.get_element(element)
            @test actual==msh
            @test cell==nodes[width*(i-1)+1:width*i]
            @test dim==Tessella.Elements.msh_dimension(msh)
            @test element in vcat(_API01.mesh.get_elements(dim,owner)[2]...)
        end
    end
    for tag in tags
        coord,_,_,_=_API01.mesh.get_node(tag)
        @test coord==coordinates[tag]
    end
    cache=_API01.mesh.get()
    @test validate(cache;reject_duplicate_cells=!allow_duplicate_cells).ok
end

function _api01_components()
    types,_,blocks=_API01.mesh.get_elements()
    adjacency=Dict{UInt64,Set{UInt64}}()
    for (msh,flat) in zip(types,blocks)
        width=Tessella.Elements.msh_num_nodes(msh)
        for offset in 1:width:length(flat)
            nodes=flat[offset:offset+width-1]
            for node in nodes
                neighbors=get!(adjacency,node,Set{UInt64}())
                union!(neighbors,nodes)
            end
        end
    end
    seen=Set{UInt64}();count=0
    for seed in keys(adjacency)
        seed in seen && continue
        count+=1;pending=UInt64[seed]
        while !isempty(pending)
            node=pop!(pending)
            node in seen && continue
            push!(seen,node)
            append!(pending,setdiff(adjacency[node],seen))
        end
    end
    return count
end

function _api01_raw_line(;native=true,quadratic=false)
    if native
        _api01_line()
    else
        _API01.model.add_discrete_entity(0,1)
        _API01.model.add_discrete_entity(0,2)
        _API01.model.add_discrete_entity(1,1,[1,2])
    end
    _API01.mesh.add_nodes(0,1,[11],[0.,0,0])
    _API01.mesh.add_nodes(0,2,[22],[1.,0,0])
    _API01.mesh.add_nodes(1,1,[33],[0.5,quadratic ? 0.2 : 0.,0])
    _API01.mesh.add_elements_by_type(1,15,[101],[11])
    _API01.mesh.add_elements_by_type(2,15,[102],[22])
    if quadratic
        _API01.mesh.add_elements_by_type(1,8,[201],[11,22,33])
    else
        _API01.mesh.add_elements_by_type(1,1,[201,202],[11,33,33,22])
    end
end

@testset "Generation zero preserves fresh geometry; generation one owns point and curve cells" begin
    _api01_begin()
    try
        for (tag,x) in ((1,0.),(2,1.),(3,0.),(4,1.),(5,0.5))
            _API01.model.add_point(x,0,0;tag=tag)
        end
        for (curve,a,b) in ((1,1,2),(2,3,4))
            _API01.model.add_line(a,b;tag=curve)
            _API01.mesh.set_transfinite_curve(curve,3)
        end
        _API01.model.add_physical_group(1,[1];tag=7)
        _API01.model.add_physical_group(0,[5];tag=8)
        _API01.mesh.generate(0)
        @test isempty(_API01.mesh.get_nodes()[1])
        @test isempty(_API01.mesh.get_elements()[1])
        _API01.mesh.generate(1)
        @test length(_API01.mesh.get_nodes()[1])==7
        types,tags,nodes=_API01.mesh.get_elements()
        @test types==Int32[1,15]
        @test length.(tags)==[4,5]
        @test tags==[UInt64[2,3,8,9],UInt64[4,5,6,7,1]]
        @test isempty(intersect(_API01.mesh.get_nodes(1,1,true)[1],
                                _API01.mesh.get_nodes(1,2,true)[1]))
        @test length(_API01.mesh.get_nodes_for_physical_group(1,7)[1])==3
        @test length(_API01.mesh.get_nodes_for_physical_group(0,8)[1])==1
        @test _api01_components()==3
        before=_api01_snapshot()
        _API01.mesh.generate(0)
        @test _api01_snapshot()==before
        for bad in (true,-1,4)
            @test_throws ArgumentError _API01.mesh.generate(bad)
            @test _api01_snapshot()==before
        end
        # Detached query payloads cannot change the authoritative cache.
        types[1]=15;tags[1][1]=999;nodes[1][1]=999
        @test _api01_snapshot()==before
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "Point15 generation uses actual stored cells and last-vertex synthesis" begin
    for native in (true,false),explicit in (:none,:first,:both),renumber in (true,false)
        _api01_begin(;renumber)
        try
            native ? _API01.model.add_point(0.,0,0;tag=1) :
                     _API01.model.add_discrete_entity(0,1)
            _API01.mesh.add_nodes(0,1,[101,205],[0.,0,0,2.,0,0])
            explicit===:first && _API01.mesh.add_elements_by_type(1,15,[501],[101])
            explicit===:both && _API01.mesh.add_elements_by_type(1,15,[501,502],[101,205])
            _API01.mesh.generate(0)
            @test length(_API01.mesh.get_nodes()[1])==2
            @test length(_API01.mesh.get_elements_by_type(15)[1])==
                (explicit===:none ? 0 : explicit===:first ? 1 : 2)
            _API01.mesh.generate(1)
            point_tags,coords,_=_API01.mesh.get_nodes(0,1)
            @test point_tags==(renumber ? UInt64[1,2] : UInt64[101,205])
            @test coords==[0.,0,0,2.,0,0]
            cells,conn=_API01.mesh.get_elements_by_type(15,1)
            @test length(cells)==(explicit===:both ? 2 : 1)
            expected=explicit===:none ? [point_tags[2]] :
                explicit===:first ? [point_tags[1]] : point_tags
            @test conn==expected
            # Renumbering changes live tags, while the allocator high-water remains.
            @test _API01.mesh.get_max_node_tag()==205
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

@testset "Source tags preserve coincident raw point vertices and disconnected chains" begin
    for dim in (0,1),renumber in (true,false)
        _api01_begin(;renumber)
        try
            _API01.model.add_discrete_entity(0,1)
            _API01.mesh.add_nodes(0,1,[901,902],zeros(6))
            _API01.mesh.generate(dim)
            @test length(_API01.mesh.get_nodes()[1])==2
            @test _API01.mesh.get_nodes()[2]==zeros(6)
            @test _API01.mesh.get_elements_by_type(15)[2]==
                (dim==0 ? UInt64[] : renumber ? UInt64[2] : UInt64[902])
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
    for only_empty in (false,true),renumber in (true,false)
        _api01_begin(;only_empty,renumber)
        try
            _API01.model.add_discrete_entity(1,1)
            raw_tags=[11,12,13,21,22,23]
            raw_conn=[11,12,12,13,21,22,22,23]
            _API01.mesh.add_nodes(1,1,raw_tags,[0.,0,0,.5,0,0,1,0,0,
                                               0,0,0,.5,0,0,1,0,0])
            _API01.mesh.add_elements_by_type(1,1,[101,102,103,104],raw_conn)
            _API01.mesh.generate(1)
            @test length(_API01.mesh.get_nodes()[1])==6
            @test _api01_components()==2
            tags,nodes=_API01.mesh.get_elements_by_type(1,1)
            @test tags==UInt64.(renumber ? (1:4) : (101:104))
            @test nodes==UInt64.(renumber ? [1,2,2,3,4,5,5,6] : raw_conn)
            @test length(_API01.mesh.get_nodes(1,1)[1])==6
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

@testset "Raw elements reference globally classified nodes atomically" begin
    _api01_begin(;renumber=false,only_empty=true)
    try
        _API01.model.add_discrete_entity(0,1)
        _API01.model.add_discrete_entity(0,2)
        _API01.model.add_discrete_entity(1,1,[1,2])
        _API01.mesh.add_nodes(0,1,[11],[0.,0,0])
        _API01.mesh.add_nodes(0,2,[22],[1.,0,0])
        _API01.mesh.add_nodes(1,1,[33],[0.5,0,0])
        _API01.mesh.add_elements_by_type(1,15,[101],[11])
        before=_api01_snapshot()
        @test_throws r"node 999" _API01.mesh.add_elements_by_type(1,1,[201],[11,999])
        @test _api01_snapshot()==before
        @test _API01.mesh.get_max_element_tag()==101
        # The rejected tag201 is reusable. Shared node11 remains Point-owned.
        _API01.mesh.add_elements_by_type(1,1,[201,202],[11,33,33,22])
        _API01.mesh.generate(1)
        @test _API01.mesh.get_elements_by_type(1,1)==
            (UInt64[201,202],UInt64[11,33,33,22])
        @test _API01.mesh.get_node(11)[3:4]==(0,1)
        @test _API01.mesh.get_node(22)[3:4]==(0,2)
        @test _API01.mesh.get_node(33)[3:4]==(1,1)
        @test _API01.mesh.get_nodes(1,1,true)[1]==UInt64[33,11,22]
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "MeshOnlyEmpty retains raw native curves without duplicate query records" begin
    for only_empty in (false,true),renumber in (true,false)
        _api01_begin(;only_empty,renumber)
        try
            _api01_raw_line()
            _API01.mesh.generate(1)
            @test length(_API01.mesh.get_nodes()[1])==3
            @test _API01.mesh.get_elements()[1]==Int32[1,15]
            @test length.(_API01.mesh.get_elements()[2])==[2,2]
            curve_nodes=_API01.mesh.get_nodes(1,1)[1]
            @test length(curve_nodes)==1
            if !renumber
                @test curve_nodes==UInt64[only_empty ? 33 : 34]
                @test _API01.mesh.get_elements_by_type(1)[1]==
                    UInt64.(only_empty ? (203:204) : (205:206))
            end
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

@testset "Discrete P2 geometry survives lower-dimensional generation prepasses" begin
    for native in (false,true),dim in (0,1),order in (1,2)
        _api01_begin(;renumber=false,only_empty=true,order)
        try
            _api01_raw_line(;native,quadratic=true)
            before=_api01_snapshot()
            _API01.mesh.generate(dim)
            if !native
                @test _api01_snapshot()==before
                @test _API01.mesh.get_node(33)[1]==[0.5,0.2,0]
            else
                @test _API01.mesh.get_elements()[1]==Int32[order==1 ? 1 : 8,15]
                @test !(UInt64(33) in _API01.mesh.get_nodes()[1])
                if order==2
                    @test _API01.mesh.get_node(34)[1]≈[0.5,0,0] atol=3e-12
                end
            end
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

@testset "Generation one keeps discrete higher-dimensional cells and ownership" begin
    for dim in (0,1)
        _api01_begin(;renumber=false)
        try
            for (entity_dim,tag) in ((1,11),(2,22),(3,33))
                _API01.model.add_discrete_entity(entity_dim,tag)
            end
            _API01.mesh.add_nodes(3,33,[101,205,309,407],
                                  [0.,0,0,1,0,0,0,1,0,0,0,1])
            _API01.mesh.add_elements_by_type(11,1,[501],[101,205])
            _API01.mesh.add_elements_by_type(22,2,[601],[101,205,309])
            _API01.mesh.add_elements_by_type(33,4,[701],[101,205,309,407])
            before=_api01_snapshot()
            _API01.mesh.generate(dim)
            @test _api01_snapshot()==before
            @test _API01.mesh.get_nodes(3,33)[1]==UInt64[101,205,309,407]
            @test isempty(_API01.mesh.get_nodes(1,11)[1])
            @test _API01.mesh.get_element(501)[3:4]==(1,11)
            @test _API01.mesh.get_element(601)[3:4]==(2,22)
            @test _API01.mesh.get_element(701)[3:4]==(3,33)
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

@testset "Generation uses Mesh.ElementOrder independently of set_order" begin
    _api01_begin()
    try
        _api01_line()
        _API01.mesh.set_order(2)
        @test _API01.option("Mesh.ElementOrder")==1
        _API01.mesh.generate(1)
        @test _API01.mesh.get_elements()[1]==Int32[1,15]
        @test length(_API01.mesh.get_nodes()[1])==3
        _API01.option("Mesh.ElementOrder",2)
        _API01.mesh.generate(0)
        @test _API01.mesh.get_elements()[1]==Int32[8,15]
        @test length(_API01.mesh.get_nodes()[1])==5
        _API01.mesh.set_order(1)
        @test _API01.option("Mesh.ElementOrder")==2
        _API01.mesh.generate(0)
        @test _API01.mesh.get_elements()[1]==Int32[8,15]
        _API01.option("Mesh.ElementOrder",1)
        _API01.mesh.generate(0)
        @test _API01.mesh.get_elements()[1]==Int32[1,15]
        @test length(_API01.mesh.get_nodes()[1])==3
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "Failed 1D callbacks preserve the current cache and curve parameters" begin
    _api01_begin()
    try
        _api01_line(;transfinite=nothing)
        _API01.mesh.generate(1)
        before=_api01_snapshot()
        params=deepcopy(_API01.CURRENT[].curve_params)
        _API01.mesh.set_size_callback((dim,tag,x,y,z,lc)->error("api01 callback failure"))
        @test_throws r"api01 callback failure" _API01.mesh.generate(1)
        @test _api01_snapshot()==before
        @test _API01.CURRENT[].curve_params==params
        _API01.mesh.remove_size_callback()
        _API01.mesh.generate(1)
        @test _api01_snapshot()==before
    finally
        _API01.finalize()
    end
end

@testset "Generation one replaces native higher cells while respecting OnlyEmpty" begin
    for only_empty in (false,true),order in (1,2)
        _api01_begin()
        try
            _API01.model.add_box(0,0,0,1,1,1;tag=1)
            for (_,curve) in _API01.model.get_entities(1)
                _API01.mesh.set_transfinite_curve(curve,2)
            end
            for (_,surface) in _API01.model.get_entities(2)
                _API01.mesh.set_transfinite_surface(surface)
                _API01.mesh.set_recombine(2,surface)
            end
            _API01.mesh.set_transfinite_volume(1)
            _API01.mesh.generate(3)
            @test _API01.mesh.get_element_types(3,1)==Int32[5]
            for (_,curve) in _API01.model.get_entities(1)
                _API01.mesh.set_transfinite_curve(curve,3)
            end
            _API01.option("Mesh.MeshOnlyEmpty",only_empty ? 1 : 0)
            _API01.option("Mesh.ElementOrder",order)
            _API01.mesh.generate(1)
            @test isempty(_API01.mesh.get_elements(3,1)[1])
            @test isempty(_API01.mesh.get_elements(2,-1)[1])
            @test _API01.mesh.get_elements()[1]==Int32[order==1 ? 1 : 8,15]
            @test length.(_API01.mesh.get_elements()[2])==[only_empty ? 12 : 24,8]
            expected=only_empty ? (order==1 ? 8 : 20) : (order==1 ? 20 : 44)
            @test length(_API01.mesh.get_nodes()[1])==expected
            @test isempty(_API01.mesh.get_nodes(3,1)[1])
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

@testset "Native 1D fields apply clipping and factor once" begin
    _api01_begin()
    try
        _api01_line(;transfinite=nothing)
        _API01.option("Mesh.MeshSizeFromPoints",0)
        _API01.option("Mesh.MeshSizeMin",0.2)
        _API01.option("Mesh.MeshSizeMax",0.2)
        _API01.option("Mesh.MeshSizeFactor",2)
        field=_API01.mesh.field.add("MathEval")
        _API01.mesh.field.set_string(field,"F","0.1")
        _API01.mesh.field.set_as_background_mesh(field)
        observed=Float64[]
        _API01.mesh.set_size_callback((dim,tag,x,y,z,lc)->begin
            push!(observed,lc)
            return lc
        end)
        _API01.mesh.generate(1)
        # Upstream clips the raw field 0.1 to0.2, then multiplies by2.
        # Callback receives the pre-clipped field. Uniform length1/h0.4 rounds
        # to3 lines; wrapping the field first would apply the factor twice.
        @test !isempty(observed)
        @test all(==(0.1),observed)
        @test length(_API01.mesh.get_elements_by_type(1)[1])==3
        @test length(_API01.mesh.get_nodes()[1])==4
        _api01_check_references()
        _API01.mesh.remove_size_callback()
        _API01.option("Mesh.CharacteristicLengthFactor",1)
        _API01.mesh.generate(1)
        @test length(_API01.mesh.get_elements_by_type(1)[1])==5
        @test _API01.option("Mesh.MeshSizeFactor")==1
    finally
        _API01.finalize()
    end
end

@testset "Periodic grading copies the master and preserves ownership" begin
    _api01_begin()
    try
        for (tag,x,y) in ((1,0.,0.),(2,1.,0.),(3,0.,2.),(4,1.,2.))
            _API01.model.add_point(x,y,0;tag)
        end
        _API01.model.add_line(1,2;tag=1)
        _API01.model.add_line(3,4;tag=2)
        _API01.mesh.set_transfinite_curve(1,5,"Progression",2.)
        _API01.mesh.set_transfinite_curve(2,8)
        transform=(1.,0,0,0,0,1.,0,2.,0,0,1.,0,0,0,0,1.)
        _API01.mesh.set_periodic(1,[2],[1],transform)
        _API01.mesh.generate(1)
        @test length(_API01.mesh.get_elements_by_type(1,1)[1])==4
        @test length(_API01.mesh.get_elements_by_type(1,2)[1])==4
        master,slave_tags,master_tags,affine=_API01.mesh.get_periodic_nodes(1,2)
        @test master==1
        @test length(slave_tags)==length(master_tags)==5
        @test collect(affine)==collect(transform)
        for (slave,leader) in zip(slave_tags,master_tags)
            @test _API01.mesh.get_node(slave)[1]==_API01.mesh.get_node(leader)[1]+[0.,2.,0]
        end
        _,master_xyz,_=_API01.mesh.get_nodes(1,1,true)
        # Saved Gmsh 4.15.2 master samples from this periodic native Line.
        @test sort(master_xyz[1:3:end])≈[0.,.06666666781522945,
            .19999999863373244,.46666666554782726,1.] atol=3e-12
        before=_api01_snapshot()
        _API01.mesh.generate(0)
        @test _api01_snapshot()==before
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "MeshOnlyVisible skips hidden native curves and retains visible point cells" begin
    _api01_begin()
    try
        _api01_line()
        _API01.model.set_visibility([(1,1)],0)
        _API01.option("Mesh.MeshOnlyVisible",1)
        _API01.mesh.generate(1)
        @test _API01.mesh.get_elements()[1]==Int32[15]
        @test length(_API01.mesh.get_nodes()[1])==2
        _API01.model.set_visibility([(1,1)],1)
        _API01.mesh.generate(1)
        @test length(_API01.mesh.get_nodes()[1])==3
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "Detached generation checks capacity and unsupported order before commit" begin
    _api01_begin()
    try
        _api01_line()
        _API01.mesh.generate(1)
        before=_api01_snapshot()
        source=_API01.CURRENT[]
        params=deepcopy(source.curve_params)
        options=_API01._session_mesh1d_options(source,"api01 capacity test")
        for kwargs in ((max_nodes=2,), (max_cells=3,),
                       (element_order=2,max_nodes=4,),
                       (max_nodes=true,), (max_cells=-1,))
            @test_throws ArgumentError _API01._generate_dim01_plan(
                source,_API01.LAST_MESH[],_API01.LAST_MESH_CLASS[],1,options,
                "api01 capacity test";kwargs...)
            @test _API01.CURRENT[]===source
            @test source.curve_params==params
            @test _api01_snapshot()==before
        end
        _API01.option("Mesh.ElementOrder",3)
        @test_throws r"order 3" _API01.mesh.generate(1)
        @test _api01_snapshot()==before
        @test source.curve_params==params
        _API01.option("Mesh.ElementOrder",1)
        _API01.mesh.generate(1)
        @test _api01_snapshot()==before
    finally
        _API01.finalize()
    end
end

@testset "1D option aliases and minimum counts reach generation" begin
    _api01_begin()
    try
        _api01_line(;transfinite=nothing)
        _API01.option("Mesh.CharacteristicLengthFromPoints",0)
        _API01.option("Mesh.CharacteristicLengthMin",1)
        _API01.option("Mesh.CharacteristicLengthMax",1)
        _API01.option("Mesh.MinimumLineNodes",5)
        @test _API01.option("Mesh.MeshSizeFromPoints")==0
        @test _API01.option("Mesh.MeshSizeMin")==1
        @test _API01.option("Mesh.MinLineNodes")==5
        _API01.mesh.generate(1)
        @test length(_API01.mesh.get_elements_by_type(1)[1])==4
        @test length(_API01.mesh.get_nodes()[1])==5
        _API01.option("Mesh.MinimumLineNodes",0)
        @test _API01.option("Mesh.MinLineNodes")==2
        _API01.mesh.generate(1)
        @test length(_API01.mesh.get_elements_by_type(1)[1])==1
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "Native Line grading retains Gmsh density-integration samples" begin
    # Gmsh 4.15.2 samples native Lines by integrating F_Transfinite and
    # inverting its sampled density primitive. These saved public getNodes
    # coordinates differ from the separate analytic parameter helper's laws.
    cases=(("Progression",2.,[0.,.03225806585200275,.09677419010543135,
                              .2258064432066724,.4838709572937826,1.]),
           ("Bump",2.,[0.,.2113246604743475,.3660253944175726,
                       .4999999999989506,.6339746055809568,.788675339524792,1.]),
           ("Beta",1.2,[0.,.0865292929320693,.2036267615379391,
                        .3559899524133273,.544421476246288,.7633519289437464,1.]))
    for (law,coefficient,expected) in cases
        _api01_begin()
        try
            _api01_line(;transfinite=nothing)
            _API01.mesh.set_transfinite_curve(1,length(expected),law,coefficient)
            _API01.mesh.generate(1)
            _,xyz,_=_API01.mesh.get_nodes(1,1,true)
            @test sort(xyz[1:3:end])≈expected atol=2e-11 rtol=0
        finally
            _API01.finalize()
        end
    end
end

@testset "Sparse source tags reach native queries and mirrored mutations" begin
    _api01_begin(;renumber=false)
    try
        _API01.model.add_discrete_entity(1,1)
        _API01.mesh.add_nodes(1,1,[500,100],[0.,0,0,1,0,0])
        _API01.mesh.add_elements_by_type(1,1,[700],[500,100])
        _API01.mesh.generate(1)
        @test isempty(_API01.mesh.get_duplicate_nodes())
        @test _API01.mesh.get_nodes_by_element_type(1)==
            (UInt64[500,100],[0.,0,0,1,0,0],Float64[])
        @test _API01.mesh.get_element_edge_nodes(1)==UInt64[500,100]
        keys=_API01.mesh.get_keys(1,"Lagrange")
        @test keys==(Int32[0,0],UInt64[500,100],[0.,0,0,1,0,0])
        @test _API01.mesh.get_keys_for_element(700,"Lagrange")==keys
        @test _API01.mesh.get_basis_functions_orientation(1,"H1Legendre1")==Int32[1]
        @test _API01.mesh.get_basis_functions_orientation_for_element(700,"H1Legendre1")==1
        located=_API01.mesh.get_element_by_coordinates(.5,0,0,1,true)
        @test located[1:3]==(UInt64(700),Int32(1),UInt64[500,100])
        @test _API01.mesh.get_elements_by_coordinates(.5,0,0,1,true)==UInt64[700]
        @test collect(_API01.mesh.get_local_coordinates_in_element(700,.25,0,0))≈[-.5,0.,0.]
        @test _API01.mesh.get_element_qualities([700],"minEdge")==[1.]
        @test _API01.mesh.get_jacobian(700,[0.,0,0])[2]==[.5]
        _API01.mesh.create_edges()
        @test _API01.mesh.get_all_edges()==(UInt64[1],UInt64[500,100])
        @test _API01.mesh.get_edges([500,100,100,500])==(UInt64[1,1],Int32[-1,1])
        _API01.mesh.reverse_elements([700])
        @test _API01.mesh.get_elements_by_type(1)==(UInt64[700],UInt64[100,500])
        @test _API01.mesh.get_basis_functions_orientation(1,"H1Legendre1")==Int32[0]
        _API01.mesh.generate(0)
        @test _API01.mesh.get_elements_by_type(1)==(UInt64[700],UInt64[100,500])
        _API01.mesh.set_node(100,[2.,0,0])
        @test _API01.mesh.get_node(100)[1]==[2.,0,0]
        @test _API01.mesh.get_element_qualities([700],"minEdge")==[2.]
        _API01.mesh.generate(0)
        @test _API01.mesh.get_node(100)[1]==[2.,0,0]
        _API01.mesh.remove_elements(1,1,[700])
        @test isempty(_API01.mesh.get_elements()[1])
        @test _API01.mesh.get_nodes()[1]==UInt64[500,100]
        @test _API01.mesh.get_max_element_tag()==700
        _API01.mesh.clear()
        @test isempty(_API01.mesh.get_nodes()[1])
        @test _API01.mesh.get_max_node_tag()==500
        @test _API01.mesh.get_max_element_tag()==700
        _API01.mesh.generate(0)
        @test isempty(_API01.mesh.get_nodes()[1])
        @test _API01.mesh.get_max_node_tag()==500
        @test _API01.mesh.get_max_element_tag()==700
    finally
        _API01.finalize()
    end
end

@testset "Sparse refinement preserves upstream support pruning and tag allocation" begin
    for quadratic in (false,true),renumber in (false,true)
        _api01_begin(;renumber=false)
        try
            _API01.model.add_discrete_entity(1,1)
            tags=quadratic ? [500,100,200] : [500,100]
            xyz=quadratic ? [0.,0,0,1,0,0,.5,.2,0] : [0.,0,0,1,0,0]
            _API01.mesh.add_nodes(1,1,tags,xyz)
            _API01.mesh.add_elements_by_type(1,quadratic ? 8 : 1,[700],tags)
            _API01.mesh.generate(0)
            before=_api01_snapshot()
            @test_throws ArgumentError _API01.mesh.refine(max_nodes=2)
            @test _api01_snapshot()==before
            _API01.option("Mesh.Renumber",renumber ? 1 : 0)
            _API01.mesh.refine()
            @test _API01.mesh.get_nodes()[1]==(renumber ? UInt64[1,2,3] : UInt64[100,500,501])
            @test _API01.mesh.get_nodes()[2]==[1.,0,0,0,0,0,.5,0,0]
            @test _API01.mesh.get_elements_by_type(1)==
                (renumber ? UInt64[1,2] : UInt64[703,704],
                 renumber ? UInt64[2,3,3,1] : UInt64[500,501,501,100])
            @test _API01.mesh.get_max_node_tag()==501
            @test _API01.mesh.get_max_element_tag()==704
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

@testset "Sparse duplicate removal uses entity and source order" begin
    for reverse_tags in (false,true)
        _api01_begin(;renumber=false)
        try
            _API01.model.add_discrete_entity(1,2)
            _API01.model.add_discrete_entity(1,1)
            first,second=reverse_tags ? (100,500) : (500,100)
            _API01.mesh.add_nodes(1,1,[first,first+1],[0.,0,0,1,0,0])
            _API01.mesh.add_nodes(1,2,[second,second+1],[0.,0,0,2,0,0])
            _API01.mesh.add_elements_by_type(1,1,[701],[first,first+1])
            _API01.mesh.add_elements_by_type(2,1,[702],[second,second+1])
            _API01.mesh.generate(1)
            @test _API01.mesh.get_duplicate_nodes()==UInt64[100,500]
            _API01.mesh.remove_duplicate_nodes()
            @test _API01.mesh.get_nodes()[1]==UInt64[first,first+1,second+1]
            @test _API01.mesh.get_elements_by_type(1,2)[2]==UInt64[first,second+1]
            @test _API01.mesh.get_node(first)[3:4]==(1,1)
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
    for reverse_source in (false,true)
        _api01_begin(;renumber=false)
        try
            _API01.model.add_discrete_entity(1,1)
            _API01.mesh.add_nodes(1,1,[500,100],[0.,0,0,1,0,0])
            tags=reverse_source ? [300,700] : [700,300]
            _API01.mesh.add_elements_by_type(1,1,tags,[500,100,500,100])
            _API01.mesh.generate(1)
            _API01.mesh.remove_duplicate_elements()
            @test _API01.mesh.get_elements_by_type(1)==(UInt64[tags[1]],UInt64[500,100])
            _API01.mesh.generate(0)
            @test _API01.mesh.get_elements_by_type(1)[1]==UInt64[tags[1]]
        finally
            _API01.finalize()
        end
    end
end

@testset "Sparse renumbering maps labels simultaneously and assigns omitted labels" begin
    cases=((:unknown,[999],[800],UInt64[801,802]),
           (:unknown_and_matched,[999,500],[800,10],UInt64[10,801]),
           (:empty,Int[],Int[],UInt64[1,2]),
           (:partial,[500],[800],UInt64[800,801]),
           (:collision_unmapped,[500],[100],UInt64[100,101]),
           (:swap,[500,100],[100,500],UInt64[100,500]),
           (:duplicate_old,[500,500],[10,20],UInt64[20,21]),
           (:duplicate_old_high_previous,[500,500],[800,20],UInt64[20,21]))
    for kind in (:node,:element),(name,old,new,expected) in cases
        _api01_begin(;renumber=false)
        try
            _API01.model.add_discrete_entity(1,1)
            _API01.mesh.add_nodes(1,1,[500,100],[0.,0,0,1,0,0])
            _API01.mesh.add_elements_by_type(1,1,[500,100],[500,100,100,500])
            _API01.mesh.generate(1)
            if kind===:node
                _API01.mesh.renumber_nodes(old,new)
                @test _API01.mesh.get_nodes()[1]==expected
                @test _API01.mesh.get_elements_by_type(1)[2]==expected[[1,2,2,1]]
                @test _API01.mesh.get_max_node_tag()==max(UInt64(500),maximum(expected))
            else
                _API01.mesh.renumber_elements(old,new)
                @test _API01.mesh.get_elements_by_type(1)[1]==expected
                @test _API01.mesh.get_nodes()[1]==UInt64[500,100]
                @test _API01.mesh.get_max_element_tag()==max(UInt64(500),maximum(expected))
            end
            before=_api01_snapshot()
            _API01.mesh.generate(0)
            @test _api01_snapshot()==before
            _api01_check_references(;allow_duplicate_cells=true)
        finally
            _API01.finalize()
        end
    end
    for kind in (:node,:element),targets in ([0,10],[10,10])
        _api01_begin(;renumber=false)
        try
            _API01.model.add_discrete_entity(1,1)
            _API01.mesh.add_nodes(1,1,[500,100],[0.,0,0,1,0,0])
            _API01.mesh.add_elements_by_type(1,1,[500,100],[500,100,100,500])
            _API01.mesh.generate(1)
            before=_api01_snapshot()
            mutate=kind===:node ? _API01.mesh.renumber_nodes : _API01.mesh.renumber_elements
            # Upstream accepts zero/duplicate final labels, which make its
            # tag-indexed mesh ambiguous. The native cache rejects these atomically.
            @test_throws ArgumentError mutate([500,100],targets)
            @test _api01_snapshot()==before
        finally
            _API01.finalize()
        end
    end
end

@testset "Mixed discrete face cells renumber in GFace family order" begin
    for quad_first in (false,true)
        _api01_begin()
        try
            _API01.model.add_discrete_entity(2,1)
            _API01.mesh.add_nodes(2,1,[101,205,309,407],[0.,0,0,1,0,0,1,1,0,0,1,0])
            for msh in (quad_first ? (3,2) : (2,3))
                _API01.mesh.add_elements_by_type(1,msh,[msh==2 ? 700 : 300],
                    msh==2 ? [101,205,309] : [101,205,309,407])
            end
            _API01.mesh.generate(0)
            @test _API01.mesh.get_elements()[1]==Int32[2,3]
            @test _API01.mesh.get_elements()[2]==[UInt64[1],UInt64[2]]
            @test _API01.mesh.get_element(1)==(Int32(2),UInt64[1,2,3],2,1)
            @test _API01.mesh.get_element(2)==(Int32(3),UInt64[1,2,3,4],2,1)
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

@testset "Automatic additions after clear retain allocator state per model" begin
    _api01_begin(;renumber=false)
    try
        _API01.model.add("allocator_a")
        _API01.model.add_discrete_entity(1,1)
        _API01.mesh.add_nodes(1,1,[500,100],[0.,0,0,1,0,0])
        _API01.mesh.add_elements_by_type(1,1,[700],[500,100])
        _API01.mesh.generate(0)
        _API01.model.add("allocator_b")
        _API01.model.add_discrete_entity(1,1)
        _API01.mesh.add_nodes(1,1,[10,20],[0.,1,0,1,1,0])
        _API01.mesh.add_elements_by_type(1,1,[30],[10,20])
        _API01.mesh.generate(0)
        @test (_API01.mesh.get_max_node_tag(),_API01.mesh.get_max_element_tag())==(20,30)
        _API01.model.set_current("allocator_a")
        @test (_API01.mesh.get_max_node_tag(),_API01.mesh.get_max_element_tag())==(500,700)
        _API01.mesh.clear()
        _API01.mesh.add_nodes(1,1,Int[],[0.,0,0,1,0,0])
        @test _API01.mesh.get_nodes()[1]==UInt64[501,502]
        _API01.mesh.add_elements_by_type(1,1,Int[],[501,502])
        @test _API01.mesh.get_elements_by_type(1)==(UInt64[701],UInt64[501,502])
        _API01.mesh.generate(0)
        @test (_API01.mesh.get_max_node_tag(),_API01.mesh.get_max_element_tag())==(502,701)
        _API01.model.set_current("allocator_b")
        @test _API01.mesh.get_nodes()[1]==UInt64[10,20]
        @test (_API01.mesh.get_max_node_tag(),_API01.mesh.get_max_element_tag())==(20,30)
        _API01.model.set_current("allocator_a")
        @test _API01.mesh.get_nodes()[1]==UInt64[501,502]
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "Selective clear removes owned curve vertices and retains point cells" begin
    _api01_begin(;renumber=false)
    try
        for point in 1:3;_API01.model.add_discrete_entity(0,point);end
        for curve in 1:2;_API01.model.add_discrete_entity(1,curve,[curve,curve+1]);end
        for (tag,node,x) in ((1,11,0.),(2,22,1.),(3,33,2.))
            _API01.mesh.add_nodes(0,tag,[node],[x,0,0])
            _API01.mesh.add_elements_by_type(tag,15,[100+tag],[node])
        end
        _API01.mesh.add_nodes(1,1,[44],[.5,0,0])
        _API01.mesh.add_nodes(1,2,[55],[1.5,0,0])
        _API01.mesh.add_elements_by_type(1,1,[201,202],[11,44,44,22])
        _API01.mesh.add_elements_by_type(2,1,[301,302],[22,55,55,33])
        _API01.mesh.generate(1)
        _API01.mesh.clear([(1,1)])
        @test _API01.mesh.get_nodes()[1]==UInt64[11,22,33,55]
        @test _API01.mesh.get_elements()[2]==[UInt64[301,302],UInt64[101,102,103]]
        @test isempty(_API01.mesh.get_nodes(1,1)[1])
        @test _API01.mesh.get_nodes(1,2,true)[1]==UInt64[55,22,33]
        @test (_API01.mesh.get_max_node_tag(),_API01.mesh.get_max_element_tag())==(55,302)
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "Periodic Lagrange keys include full quadratic curve nodes" begin
    _api01_begin(;order=2)
    try
        for (tag,x,y) in ((1,0.,0.),(2,1.,0.),(3,0.,2.),(4,1.,2.))
            _API01.model.add_point(x,y,0;tag)
        end
        _API01.model.add_line(1,2;tag=1);_API01.model.add_line(3,4;tag=2)
        for curve in 1:2;_API01.mesh.set_transfinite_curve(curve,3);end
        _API01.mesh.set_periodic(1,[2],[1],[1.,0,0,0,0,1.,0,2.,0,0,1.,0,0,0,0,1.])
        _API01.mesh.generate(1)
        master,types,master_types,keys,master_keys,coords,master_coords=
            _API01.mesh.get_periodic_keys(8,"Lagrange",2)
        @test master==1
        @test types==master_types==zeros(Int32,6)
        @test keys==UInt64[3,8,9,8,4,10]
        @test master_keys==UInt64[1,5,6,5,2,7]
        # Gmsh's getPeriodicKeys protocol returns master coordinates first.
        @test coords[2:3:end]==zeros(6)
        @test master_coords[2:3:end]==fill(2.,6)
        @test sort(unique(coords[1:3:end]))==[0.,.25,.5,.75,1.]
        primary=_API01.mesh.get_periodic_nodes(1,2)
        full=_API01.mesh.get_periodic_nodes(1,2,true)
        @test Dict(zip(primary.slave_nodes,primary.master_nodes))==Dict(3=>1,4=>2,8=>5)
        @test Dict(zip(full.slave_nodes,full.master_nodes))==Dict(3=>1,4=>2,8=>5,9=>6,10=>7)
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "Sparse element display state and homology use public identities" begin
    _api01_begin(;renumber=false)
    try
        _API01.model.add_discrete_entity(1,1)
        _API01.mesh.add_nodes(1,1,[500,100,300],[0.,0,0,.5,0,0,1.,0,0])
        _API01.mesh.add_elements_by_type(1,1,[700,300],[500,100,100,300])
        _API01.mesh.generate(0)
        @test _API01.mesh.get_visibility([700,300,1])==Int32[1,1,0]
        _API01.mesh.set_visibility([700,1],5)
        @test _API01.mesh.get_visibility([700,300,1])==Int32[5,1,0]
        _API01.mesh.set_visibility_per_window(700,0,2)
        @test _API01.MODEL_SLOTS[_API01._current_slot_index_locked()].window_element_visibility[2][700]==0
        _API01.mesh.reorder_elements(1,1,[1,0])
        @test _API01.mesh.get_elements_by_type(1)==
            (UInt64[300,700],UInt64[100,300,500,100])
        @test _API01.mesh.get_visibility([700,300])==Int32[5,1]
        _API01.mesh.add_homology_request("Homology",Int[],Int[],[0])
        groups=_API01.mesh.compute_homology()
        @test length(groups)==1
        @test groups[1][1]==0
        @test Set(_API01.mesh.get_elements_by_type(15)[2])⊆Set(UInt64[500,100,300])
    finally
        _API01.finalize()
    end
end

@testset "Coordinate lookup prioritizes curve cells over endpoint Point15 cells" begin
    _api01_begin()
    try
        _api01_line()
        _API01.mesh.generate(1)
        for (x,line_tag,point_tag,coordinate) in ((0.,3,1,-1.),(1.,4,2,1.)),
            dimension in (-1,0,1),strict in (false,true)
            single=_API01.mesh.get_element_by_coordinates(x,0,0,dimension,strict)
            many=_API01.mesh.get_elements_by_coordinates(x,0,0,dimension,strict)
            @test single[1]==(dimension==0 ? point_tag : line_tag)
            @test single[2]==(dimension==0 ? 15 : 1)
            @test single[4]==(dimension==0 ? 0. : coordinate)
            @test many==UInt64.(dimension<0 ? [line_tag,point_tag] :
                               dimension==0 ? [point_tag] : [line_tag])
        end
    finally
        _API01.finalize()
    end
end

@testset "Tagged identity-preserving edits retain element display state" begin
    for action in (:set_node,:reverse,:reverse_elements,:affine,:add_node,:add_element,
                   :renumber_node,:renumber_element,:remove_other_element,:remove_element)
        _api01_begin(;renumber=false)
        try
            _API01.model.add_discrete_entity(1,1)
            _API01.mesh.add_nodes(1,1,[500,100,300],[0.,0,0,.5,0,0,1.,0,0])
            _API01.mesh.add_elements_by_type(1,1,[700,300],[500,100,100,300])
            _API01.mesh.generate(0)
            _API01.mesh.set_visibility([700],5)
            _API01.mesh.set_visibility_per_window(700,0,2)
            tag=700
            if action===:set_node
                _API01.mesh.set_node(100,[.6,0,0])
            elseif action===:reverse
                _API01.mesh.reverse([(1,1)])
            elseif action===:reverse_elements
                _API01.mesh.reverse_elements([700])
            elseif action===:affine
                _API01.mesh.affine_transform([1.,0,0,0,0,1.,0,1.,0,0,1.,0,0,0,0,1.])
            elseif action===:add_node
                _API01.mesh.add_nodes(1,1,[800],[2.,0,0])
                @test isempty(_API01.mesh.get_node(800)[2])
            elseif action===:add_element
                _API01.mesh.add_elements_by_type(1,1,[900],[500,300])
            elseif action===:renumber_node
                _API01.mesh.renumber_nodes([500],[800])
            elseif action===:renumber_element
                _API01.mesh.renumber_elements([700],[800]);tag=800
            elseif action===:remove_other_element
                _API01.mesh.remove_elements(1,1,[300])
            else
                _API01.mesh.remove_elements(1,1,[700])
            end
            @test _API01.mesh.get_visibility([tag])==Int32[action===:remove_element ? 0 : 5]
            window=_API01.MODEL_SLOTS[_API01._current_slot_index_locked()].window_element_visibility
            @test get(get(window,2,Dict{Int,Int32}()),tag,Int32(-1))==
                (action===:remove_element ? -1 : 0)
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

include("api_generate01_artifacts.jl")
@testset "Generation 0/1 typed artifacts retain exact metadata across runtimes" begin
    expected=Dict{String,Vector{String}}()
    for line in eachline(joinpath(@__DIR__,"..","artifacts","api_generate01_crc.txt"))
        (isempty(line) || startswith(line,'#')) && continue
        fields=split(line,'\t')
        @test length(fields)==7
        expected[fields[1]]=fields
    end
    records=APIGenerate01Artifacts.records()
    @test length(records)==length(expected)==9
    for record in records
        fields=expected[record.name]
        @test [record.n_nodes,record.n_blocks,record.n_cells,
               Int(record.max_node),Int(record.max_element)]==parse.(Int,fields[2:6])
        @test record.sha==fields[7]
    end
end

@testset "Legacy native higher caches reject untracked attachment allocations atomically" begin
    _api01_begin(;renumber=false)
    try
        _API01.model.add_box(0,0,0,1,1,1;tag=1)
        xyz=collect(_API01.model.get_bounding_box(0,1)[1:3])
        _API01.mesh.add_nodes(0,1,[111],xyz)
        _API01.mesh.add_elements_by_type(1,15,[777],[111])
        for (_,curve) in _API01.model.get_entities(1)
            _API01.mesh.set_transfinite_curve(curve,2)
        end
        for (_,surface) in _API01.model.get_entities(2)
            _API01.mesh.set_transfinite_surface(surface)
            _API01.mesh.set_recombine(2,surface)
        end
        _API01.mesh.set_transfinite_volume(1)
        _API01.mesh.generate(3)
        for (_,curve) in _API01.model.get_entities(1)
            _API01.mesh.set_transfinite_curve(curve,3)
        end
        model=_API01.CURRENT[];cache=_API01.LAST_MESH[]
        class=_API01.LAST_MESH_CLASS[];before=_api01_snapshot()
        point_payload=_API01.mesh.get_nodes(0,1)
        parameters=deepcopy(model.curve_params)
        attributes=deepcopy(model.meshing.transfinite_curves)
        counters=(_API01.mesh.get_max_node_tag(),_API01.mesh.get_max_element_tag())
        @test_throws r"untracked source tag allocations" _API01.mesh.generate(1)
        @test _API01.CURRENT[]===model
        @test _API01.LAST_MESH[]===cache
        @test _API01.LAST_MESH_CLASS[]===class
        @test _api01_snapshot()==before
        @test model.curve_params==parameters
        @test model.meshing.transfinite_curves==attributes
        @test (_API01.mesh.get_max_node_tag(),_API01.mesh.get_max_element_tag())==counters
        @test _API01.mesh.get_nodes(0,1)==point_payload
    finally
        _API01.finalize()
    end
end

@testset "Unchanged quadratic cells preserve visibility through order operations" begin
    _api01_begin(;renumber=false,order=2)
    try
        _api01_line()
        _API01.mesh.generate(1)
        tag=first(_API01.mesh.get_elements_by_type(8)[1])
        _API01.mesh.set_visibility([tag],5)
        _API01.mesh.set_order(2)
        @test _API01.mesh.get_visibility([tag])==Int32[5]
        _API01.mesh.set_order(1)
        @test _API01.mesh.get_visibility([tag])==Int32[0]
        @test _API01.mesh.get_visibility(_API01.mesh.get_elements_by_type(1)[1])==Int32[1,1]
        _api01_check_references()
    finally
        _API01.finalize()
    end
    _api01_begin(;renumber=false)
    try
        _api01_raw_line(;native=false,quadratic=true)
        _API01.mesh.generate(0)
        _API01.mesh.set_visibility([201],5)
        _API01.mesh.generate(0)
        @test _API01.mesh.get_visibility([201])==Int32[5]
        @test _API01.mesh.get_element_types()==Int32[8,15]
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "Immediate order changes ingest raw point nodes and curve cells" begin
    for order in (1,2),explicit in (false,true),renumber in (false,true)
        _api01_begin(;renumber)
        try
            _API01.model.add_discrete_entity(0,1)
            _API01.mesh.add_nodes(0,1,[11,13],[0.,0,0,2.,0,0])
            explicit && _API01.mesh.add_elements_by_type(1,15,[101],[11])
            @test _API01.LAST_MESH[]===nothing
            _API01.mesh.set_order(order)
            retained=order==1 ? [11,13] : explicit ? [11] : Int[]
            @test _API01.mesh.get_nodes()[1]==UInt64.(renumber ? (1:length(retained)) : retained)
            @test length(_API01.mesh.get_elements_by_type(15)[1])==(explicit ? 1 : 0)
            @test _API01.option("Mesh.ElementOrder")==1
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
    for quadratic in (false,true)
        _api01_begin(;renumber=false)
        try
            _api01_raw_line(;native=false,quadratic)
            @test _API01.LAST_MESH[]===nothing
            _API01.mesh.set_order(quadratic ? 1 : 2)
            @test _API01.mesh.get_element_types(1,1)==Int32[quadratic ? 1 : 8]
            @test length(_API01.mesh.get_nodes()[1])==(quadratic ? 2 : 5)
            @test _API01.mesh.get_elements_by_type(15)==(UInt64[101,102],UInt64[11,22])
            quadratic && @test _API01.mesh.get_elements_by_type(1)==(UInt64[202],UInt64[11,22])
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

function _api01_triangle_square(offset,surface)
    for (slot,p) in enumerate(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)))
        _API01.model.add_point(p...;tag=offset+slot,meshSize=10)
    end
    for edge in 1:4
        _API01.model.add_line(offset+edge,offset+mod1(edge+1,4);tag=offset+edge)
        _API01.mesh.set_transfinite_curve(offset+edge,2)
    end
    _API01.model.add_curve_loop(offset .+ collect(1:4);tag=surface)
    _API01.model.add_plane_surface([surface];tag=surface)
    _API01.mesh.set_transfinite_surface(surface)
end

function _api01_tetra_geometry()
    for (tag,p) in enumerate(((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.)))
        _API01.model.add_point(p...;tag,meshSize=10)
    end
    for (tag,(a,b)) in enumerate(((1,2),(2,3),(3,1),(1,4),(2,4),(3,4)))
        _API01.model.add_line(a,b;tag);_API01.mesh.set_transfinite_curve(tag,2)
    end
    for (tag,curves) in enumerate(([1,2,3],[1,5,-4],[2,6,-5],[3,4,-6]))
        _API01.model.add_curve_loop(curves;tag);_API01.model.add_plane_surface([tag];tag)
        _API01.mesh.set_transfinite_surface(tag)
    end
    _API01.model.add_surface_loop([-1,2,3,4];tag=1)
    _API01.model.add_volume([1];tag=1)
end

@testset "Simplex P2 overlays complete their actual native boundary source" begin
    for geometry in (:square,:coincident_squares,:tetra),order in (1,2)
        _api01_begin(;order)
        try
            if geometry===:tetra
                _api01_tetra_geometry()
            else
                _api01_triangle_square(0,1)
                geometry===:coincident_squares && _api01_triangle_square(4,2)
            end
            _API01.model.add_physical_group(geometry===:tetra ? 3 : 2,[1];tag=7,name="source")
            _API01.mesh.generate(geometry===:tetra ? 3 : 2)
            original=_API01.LAST_MESH[];original_class=_API01.LAST_MESH_CLASS[]
            @test original isa Tessella.Mesh
            actual,class=_API01._dim01_actual_cache(original,original_class;
                physical_names=_API01.CURRENT[].physical_names)
            source,source_class=_API01._dim01_complete_native_cache(
                _API01.CURRENT[],actual,class,"api01 simplex completion regression")
            @test source.coords==actual.coords
            @test _API01.LAST_MESH[]===original && _API01.LAST_MESH_CLASS[]===original_class
            @test source.physical_names[(geometry===:tetra ? 3 : 2,7)]=="source"
            counts=Dict(block.msh=>size(block.nodes,2) for block in source.blocks)
            factor=geometry===:coincident_squares ? 2 : 1
            expected=geometry===:tetra ? Dict(15=>4,(order==1 ? 1 : 8)=>6,
                (order==1 ? 2 : 9)=>4,(order==1 ? 4 : 11)=>4) :
                Dict(15=>4factor,(order==1 ? 1 : 8)=>4factor,(order==1 ? 2 : 9)=>2factor)
            @test counts==expected
            @test size(source.coords,2)==(geometry===:tetra ? (order==1 ? 5 : 15) :
                                        factor*(order==1 ? 4 : 9))
            if factor==2
                cells=only(block for block in source.blocks if Tessella.Elements.msh_dimension(block.msh)==2)
                block_index=findfirst(block->block===cells,source.blocks)
                owners=source.entity_data.block_entities[block_index]
                @test isempty(intersect(Set(vec(cells.nodes[:,owners.==1])),Set(vec(cells.nodes[:,owners.==2]))))
            end
            @test all(owner[1]!=1 || source.entity_data.node_parametric[node]!==nothing
                for (node,owner) in enumerate(source_class.node_entities))
            _API01.option("Mesh.MeshOnlyEmpty",1)
            _API01.mesh.generate(1)
            @test isempty(_API01.mesh.get_elements(2,-1)[1])
            @test isempty(_API01.mesh.get_elements(3,-1)[1])
            _api01_check_references()
        finally
            _API01.finalize()
        end
    end
end

@testset "Periodic compaction and selective clear publish only live pairs atomically" begin
    _api01_begin()
    try
        for (tag,x,y) in ((1,0.,0.),(2,1.,0.),(3,0.,2.),(4,1.,2.))
            _API01.model.add_point(x,y,0;tag)
        end
        _API01.model.add_line(1,2;tag=1);_API01.model.add_line(3,4;tag=2)
        for curve in 1:2;_API01.mesh.set_transfinite_curve(curve,3);end
        affine=[1.,0,0,0,0,1.,0,2.,0,0,1.,0,0,0,0,1.]
        _API01.mesh.set_periodic(1,[2],[1],affine)
        _API01.mesh.generate(1)
        _API01.mesh.set_order(2)
        @test length(_API01.mesh.get_periodic_nodes(1,2,true).slave_nodes)==5
        _API01.mesh.set_order(1)
        @test length(_API01.mesh.get_periodic_nodes(1,2,true).slave_nodes)==3
        model=_API01.CURRENT[];cache=_API01.LAST_MESH[];class=_API01.LAST_MESH_CLASS[]
        before=_api01_snapshot();params=deepcopy(model.curve_params)
        raw=deepcopy(model.meshing.attached[(1,2)].element_nodes)
        @test_throws ArgumentError _API01.mesh.clear([(1,2),(1,999)])
        @test _API01.CURRENT[]===model && _API01.LAST_MESH[]===cache && _API01.LAST_MESH_CLASS[]===class
        @test _api01_snapshot()==before
        @test model.curve_params==params
        @test model.meshing.attached[(1,2)].element_nodes==raw
        _API01.mesh.clear([(1,2)])
        mapping=_API01.mesh.get_periodic_nodes(1,2,true)
        @test isempty(mapping.slave_nodes) && isempty(mapping.master_nodes)
        mesh=_API01.mesh.get()
        @test length(mesh.periodic_links)==1
        @test isempty(mesh.periodic_links[1].slave_nodes) && isempty(mesh.periodic_links[1].master_nodes)
        @test validate(mesh).ok
        mktempdir() do directory
            path=joinpath(directory,"periodic_partial_clear.msh")
            Tessella.Elements.write_mixed_msh(path,mesh;version=4.1)
            reread=Tessella.Elements.read_mixed_msh(path)
            @test length(reread.periodic_links)==1
            @test isempty(reread.periodic_links[1].slave_nodes) && isempty(reread.periodic_links[1].master_nodes)
        end
        _api01_check_references()
    finally
        _API01.finalize()
    end
end

@testset "Unmeshed hidden curved entities do not block straight curve refinement" begin
    _api01_begin()
    try
        mktempdir() do directory
            path=joinpath(directory,"hidden_refine.geo")
            write(path,"Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={2,0,0,1};" *
                "Point(4)={2,1,0,1};Point(5)={3,1,0,1};Line(1)={1,2};" *
                "Circle(2)={3,4,5};Transfinite Curve{1}=3;")
            _API01.open_geo!(path)
            _API01.model.set_visibility([(1,2)],0)
            _API01.option("Mesh.MeshOnlyVisible",1)
            _API01.mesh.generate(1)
            @test length(_API01.mesh.get_elements_by_type(1,1)[1])==2
            @test isempty(_API01.mesh.get_elements(1,2)[1])
            _API01.mesh.refine()
            @test length(_API01.mesh.get_elements_by_type(1,1)[1])==4
            @test isempty(_API01.mesh.get_elements(1,2)[1])
            @test sort(_API01.mesh.get_nodes(1,1)[2][1:3:end])==[.25,.5,.75]
            _api01_check_references()
        end
    finally
        _API01.finalize()
    end
end
