using Test,Tessella
using Tessella.Elements: ElementBlock,MixedMesh
using SHA
include(joinpath(@__DIR__,"..","artifacts","native_public_square_quality_v1.jl"))
const _PQ=Tessella.API

function _pq_square(quadratic)
    _PQ.initialize();_PQ.option("Mesh.Renumber",0)
    for (x,y) in ((0,0),(1,0),(1,1),(0,1));_PQ.model.add_point(x,y,0);end
    for (a,b) in ((1,2),(2,3),(3,4),(4,1));_PQ.model.add_line(a,b);end
    _PQ.model.add_curve_loop([1,2,3,4]);_PQ.model.add_plane_surface([1])
    _PQ.model.add_physical_group(2,[1];tag=99)
    for curve in 1:4;_PQ.mesh.set_transfinite_curve(curve,3);end
    _PQ.mesh.set_transfinite_surface(1);_PQ.mesh.generate(2)
    quadratic && _PQ.mesh.set_order(2)
end

@testset "Native ordinary Point and P2 Line metrics retain public payload" begin
    for quadratic in (false,true)
        _pq_square(quadratic)
        try
            native=quadratic ? _PQ_NATIVE_SQUARE.p2 : _PQ_NATIVE_SQUARE.p1
            @test _pq_public_square_payload()==native.payload
            tags=vcat(native.payload.elements[2]...)
            @test eltype(_PQ.mesh.get_nodes()[1])===UInt64
            @test eltype(_PQ.mesh.get_elements()[1])===Int32
            @test all(block->eltype(block)===UInt64,_PQ.mesh.get_elements()[2])
            for metric in ("volume","minEdge","maxEdge")
                expected=getproperty(native,Symbol(metric))
                @test _PQ.mesh.get_element_qualities(tags,metric)==expected
                for (position,tag) in enumerate(tags)
                    @test only(_PQ.mesh.get_element_qualities([tag],metric))==expected[position]
                end
            end
            @test _pq_public_square_payload()==native.payload
            pointtags,_=_PQ.mesh.get_elements_by_type(15)
            @test _PQ.mesh.get_element_qualities(pointtags,"volume")==fill(1.0,4)
            @test _PQ.mesh.get_element_qualities(pointtags,"minEdge")==fill(1e25,4)
            @test _PQ.mesh.get_element_qualities(pointtags,"maxEdge")==zeros(4)
            lines,_=_PQ.mesh.get_elements_by_type(quadratic ? 8 : 1)
            selected=findall(tag->tag in Set(lines),tags)
            @test _PQ.mesh.get_element_qualities(lines,"volume")==native.volume[selected]
            for absent in (quadratic ? (1,2) : (8,9))
                emptytags,_=_PQ.mesh.get_elements_by_type(absent)
                @test isempty(emptytags)
                @test _PQ.mesh.get_element_qualities(emptytags,"volume")==Float64[]
            end
            records=Dict(tag=>_PQ.mesh.get_element(tag) for tag in tags)
            before_nodes=_PQ.mesh.get_nodes();before_elements=_PQ.mesh.get_elements()
            new=reverse(tags).+UInt64(1000)
            _PQ.mesh.renumber_elements(tags,new)
            @test _PQ.mesh.get_element_qualities(new,"volume")==native.volume
            @test _PQ.mesh.get_nodes()==before_nodes
            @test _PQ.mesh.get_elements()[1]==before_elements[1]
            @test _PQ.mesh.get_elements()[3]==before_elements[3]
            for (old,tag) in zip(tags,new);@test _PQ.mesh.get_element(tag)==records[old];end
            @test _PQ.mesh.get_max_element_tag()==native.payload.maxima[2]+UInt64(1000)
        finally
            _PQ.finalize()
        end
    end
end

function _pq_scale_fixture(count)
    coords=zeros(Float64,3,4count);points=Matrix{Int32}(undef,1,count);lines=Matrix{Int32}(undef,3,count)
    for cell in 1:count
        row=4cell-3;points[1,cell]=Int32(row)
        lines[:,cell].=Int32[row+1,row+2,row+3]
        coords[1,row]=3cell;coords[1,row+1]=3cell;coords[1,row+2]=3cell+1;coords[1,row+3]=3cell+.5
    end
    return MixedMesh(coords,[ElementBlock(15,points),ElementBlock(8,lines)])
end

function _pq_has_box(value)
    value===Core.Box && return true
    value isa GlobalRef && value.mod===Core && value.name===:Box && return true
    value isa Expr && return any(_pq_has_box,value.args)
    value isa Core.CodeInfo && return any(_pq_has_box,value.code)
    value isa AbstractArray && return any(_pq_has_box,value)
    return false
end

function _pq_box_positive_control()
    captured=0
    increment=()->(captured+=1)
    return increment
end

@testset "Actual Point and Line quality method and linear allocations" begin
    @test _pq_has_box(Base.uncompressed_ast(which(_pq_box_positive_control,Tuple{})))
    allocations=Int[];records=Any[]
    for count in (1000,2000,4000)
        mesh=_pq_scale_fixture(count);tags=UInt64.(1:2count)
        method=which(_PQ.mesh_element_qualities,(typeof(mesh),typeof(tags),String))
        @test basename(String(method.file))=="APIMixedQueries.jl"
        @test normpath(String(method.file))==normpath(joinpath(pkgdir(Tessella),"src","interfaces","APIMixedQueries.jl"))
        lowered=Base.uncompressed_ast(method)
        @test !_pq_has_box(lowered)
        @test bytes2hex(sha256(read(String(method.file))))=="b453eb5f7c2751efc0eebc38ce8cef14e3d834494963947dbce4036a483011a4"
        result=_PQ.mesh_element_qualities(mesh,tags,"volume")
        @test result==ones(2count)
        @test length(result)==2count
        @test Base.summarysize(result)<=64+8length(result)
        _PQ.mesh_element_qualities(mesh,tags,"volume");GC.gc()
        push!(allocations,minimum(@allocated(_PQ.mesh_element_qualities(mesh,tags,"volume")) for _ in 1:3))
        push!(records,(method=string(method),source_sha256=bytes2hex(sha256(read(String(method.file)))),
            lowered_sha256=bytes2hex(sha256(repr(lowered))),Core_Box=_pq_has_box(lowered)))
    end
    for index in 2:3;@test allocations[index]<=2.5allocations[index-1]+262144;end
    @info "Point and Line quality original native measures" counts=(1000,2000,4000) allocations=Tuple(allocations) actual_methods=records
end

const _PL_RENUMBER_SOURCE_SHA=Dict("API.jl"=>"3e97a24e6b32bff0cee4f884e36bc74449dfa90f157996cf159be1c3953dc259","APINativeMeshTags.jl"=>"c79e2fa26113bf3abed468dedc9962bdf7ebcf5403529356687df1856ac5c154","APITaggedMutations.jl"=>"079714c1cfe489a85a135db631f012c69e34de7340f781a5f5d721d0cd875af7")

include(joinpath(@__DIR__,"..","artifacts","native_public_labels_lifecycle_v3.jl"))

function _pl_square_geometry(;physical=true,save_all=0)
    _PQ.initialize();_PQ.option("Mesh.Renumber",0);_PQ.option("Mesh.SaveAll",save_all)
    for (x,y) in ((0,0),(1,0),(1,1),(0,1));_PQ.model.add_point(x,y,0);end
    for (a,b)in ((1,2),(2,3),(3,4),(4,1));_PQ.model.add_line(a,b);end
    _PQ.model.add_curve_loop([1,2,3,4]);_PQ.model.add_plane_surface([1])
    physical && _PQ.model.add_physical_group(2,[1];tag=99)
end

@testset "Native empty element labels prioritize physical entities across cached and raw routes" begin
    profiles=((:physical_empty_default,true,0),(:physical_empty_saveall0,true,0),
        (:physical_empty_saveall1,true,1),(:no_physical_empty_default,false,0))
    for quadratic in (false,true),(profile,physical,save_all)in profiles,repeat in 1:3
        _pl_square_geometry(;physical,save_all)
        try
            for curve in 1:4;_PQ.mesh.set_transfinite_curve(curve,3);end
            _PQ.mesh.set_transfinite_surface(1);_PQ.mesh.generate(2)
            quadratic && _PQ.mesh.set_order(2)
            @test _pq_public_square_payload()==_pl_payload(quadratic,profile,:before_operation)
            old=vcat(_PQ.mesh.get_elements()[2]...)
            _PQ.mesh.renumber_elements(old,reverse(old).+UInt64(1000))
            @test _pq_public_square_payload()==_pl_payload(quadratic,profile,:before_empty_renumber)
            _PQ.mesh.renumber_elements()
            @test _pq_public_square_payload()==_pl_payload(quadratic,profile,:after_operation)
        finally
            _PQ.finalize()
        end
    end
    rawprofiles=((:raw_physical_empty_default,true,0),(:raw_physical_empty_saveall1,true,1),
        (:no_physical_raw_empty_default,false,0))
    for quadratic in (false,true),(profile,physical,save_all)in rawprofiles,repeat in 1:3
        _pl_square_geometry(;physical,save_all)
        try
            seed=quadratic ? _PQ_NATIVE_SQUARE.p2.payload : _PQ_NATIVE_SQUARE.p1.payload
            for row in seed.owners
                d,t=row.owner
                isempty(row.nodes[1]) || _PQ.mesh.add_nodes(d,t,row.nodes[1],row.nodes[2],row.nodes[3])
            end
            for row in seed.owners
                d,t=row.owner
                for (kind,tags,nodes)in zip(row.elements...)
                    _PQ.mesh.add_elements_by_type(t,kind,tags,nodes)
                end
            end
            @test _PQ.LAST_MESH[]===nothing
            @test _pq_public_square_payload()==_pl_payload(quadratic,profile,:before_operation)
            old=vcat(_PQ.mesh.get_elements()[2]...)
            _PQ.mesh.renumber_elements(old,reverse(old).+UInt64(1000))
            @test _pq_public_square_payload()==_pl_payload(quadratic,profile,:before_empty_renumber)
            _PQ.mesh.renumber_elements()
            @test _pq_public_square_payload()==_pl_payload(quadratic,profile,:after_operation)
            @test _PQ.LAST_MESH[]===nothing
        finally
            _PQ.finalize()
        end
    end
end

@testset "Native accepted zero unknown and signed-low32 element labels retain complete bulk payload" begin
    # This entire testset uses bulk/perentity getters only. Scalar/cache/quality
    # calls are deliberately absent after zero or upper-half public labels.
    profiles=(:old_unknown_zero_full,:wide_2pow31,:int32_max,:zero_singleton,
        :empty_both,:old_unknown_singleton,:duplicate_old_last_wins,:full_cycle,:same_singleton)
    for quadratic in (false,true),profile in profiles,repeat in 1:3
        _pq_square(quadratic)
        try
            native=_PL_NATIVE_BOUNDARY[(quadratic ? 2 : 1,profile)]
            old=vcat(_PQ.mesh.get_elements()[2]...)
            _PQ.mesh.renumber_elements(old,reverse(old).+UInt64(1000))
            @test _pl_public_bulk_payload()==native.before
            _PQ.mesh.renumber_elements(native.old,native.new)
            @test _pl_public_bulk_payload()==native.after
            @test eltype(_PQ.mesh.get_nodes()[1])===UInt64
            @test eltype(_PQ.mesh.get_elements()[1])===Int32
            @test all(tags->eltype(tags)===UInt64,_PQ.mesh.get_elements()[2])
        finally
            _PQ.finalize()
        end
    end
end

@testset "Actual empty-map source methods preserve checked routing bodies" begin
    _pq_square(false)
    try
        model=_PQ.CURRENT[];mesh=_PQ.LAST_MESH[];class=_PQ.LAST_MESH_CLASS[]
        bindings=((_PQ._native_empty_element_physical_owners,(typeof(model),),"APITaggedMutations.jl"),
            (_PQ._tagged_renumber_order,(typeof(model),typeof(mesh),typeof(class),Bool,String),"APITaggedMutations.jl"),
            (_PQ._tagged_renumber_plan,(typeof(model),typeof(mesh),typeof(class),Vector{UInt64},Vector{UInt64},Symbol,String),"APITaggedMutations.jl"),
            (_PQ._renumber_elements_native_lookup_live,(Vector{UInt64},Vector{UInt64}),"API.jl"))
        @test _pq_has_box(Base.uncompressed_ast(which(_pq_box_positive_control,Tuple{})))
        rows=Any[]
        for (fn,args,file)in bindings
            method=which(fn,args);actual=normpath(String(method.file))
            @test actual==normpath(joinpath(pkgdir(Tessella),"src","interfaces",file))
            @test bytes2hex(sha256(read(actual)))==_PL_RENUMBER_SOURCE_SHA[file]
            lowered=Base.uncompressed_ast(method);@test !_pq_has_box(lowered)
            push!(rows,(method=string(method),source_sha256=bytes2hex(sha256(read(actual))),
                lowered_sha256=bytes2hex(sha256(repr(lowered))),Core_Box=_pq_has_box(lowered)))
        end
        # Bind the actual raw-record keyword body as well as its selected API caller.
        records=Tessella.Model.DiscreteEntity[]
        kwargs=(model=model,)
        method=which(Core.kwcall,(typeof(kwargs),typeof(_PQ._native_record_renumber_mapping),
            typeof(records),Vector{UInt64},Vector{UInt64},Bool,String))
        actual=normpath(String(method.file));@test actual==normpath(joinpath(pkgdir(Tessella),"src","interfaces","APINativeMeshTags.jl"))
        @test bytes2hex(sha256(read(actual)))==_PL_RENUMBER_SOURCE_SHA["APINativeMeshTags.jl"]
        lowered=Base.uncompressed_ast(method);@test !_pq_has_box(lowered)
        push!(rows,(method=string(method),source_sha256=bytes2hex(sha256(read(actual))),
            lowered_sha256=bytes2hex(sha256(repr(lowered))),Core_Box=_pq_has_box(lowered)))
        @info "Actual empty element renumber methods" records=rows
    finally
        _PQ.finalize()
    end
end

@testset "Physical-first raw renumber records keep original linear growth caps" begin
    allocations=Int[];retained=Int[]
    for count in (1000,2000,4000)
        _PQ.initialize()
        try
            _PQ.model.add_discrete_entity(2,1)
            _PQ.model.add_physical_group(2,[1];tag=99)
            coords=zeros(Float64,3,3count)
            for cell in 1:count
                coords[:,3cell-2:3cell]=Float64[3cell 3cell+1 3cell;0 0 1;0 0 0]
            end
            _PQ.mesh.add_nodes(2,1,UInt64.(1:3count),vec(coords))
            _PQ.mesh.add_elements_by_type(1,2,UInt64.(1:count),UInt64.(1:3count))
            @test _PQ.LAST_MESH[]===nothing
            _PQ.mesh.renumber_elements(UInt64.(1:count),UInt64.(10001:10000+count))
            _PQ.mesh.renumber_elements()
            @test _PQ.mesh.get_elements_by_type(2)[1]==UInt64.(1:count)
            @test _PQ.mesh.get_max_element_tag()==10000+count
            _PQ.mesh.renumber_elements();GC.gc()
            push!(allocations,minimum(@allocated(_PQ.mesh.renumber_elements()) for _ in 1:3))
            push!(retained,Base.summarysize(_PQ.CURRENT[]))
            @test _PQ.LAST_MESH[]===nothing
        finally
            _PQ.finalize()
        end
    end
    for index in 2:3
        @test allocations[index]<=2.5allocations[index-1]+262144
        @test retained[index]<=2.5retained[index-1]+262144
    end
    @info "Physical-first raw element record storage" cells=(1000,2000,4000) allocations=Tuple(allocations) retained=Tuple(retained)
end
