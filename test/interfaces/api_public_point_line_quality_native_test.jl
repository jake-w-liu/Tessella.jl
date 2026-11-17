using Test,Tessella
using Tessella.Elements: ElementBlock,MixedMesh
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
            before=(_PQ.mesh.get_nodes(),_PQ.mesh.get_elements())
            types,labels,_=before[2]
            @test types==Int32[quadratic ? 8 : 1,quadratic ? 9 : 2,15]
            @test length.(labels)==[8,8,4]
            tags=vcat(labels...)
            kinds=Dict(tag=>kind for (kind,block) in zip(types,labels) for tag in block)
            maxima=(_PQ.mesh.get_max_node_tag(),_PQ.mesh.get_max_element_tag())
            expected_volume=[kinds[tag]==15 ? 1.0 : kinds[tag] in (1,8) ? .5 : .125 for tag in tags]
            volumes=_PQ.mesh.get_element_qualities(tags,"volume")
            @test eltype(_PQ.mesh.get_nodes()[1])===UInt64
            @test eltype(_PQ.mesh.get_elements()[1])===Int32
            @test all(block->eltype(block)===UInt64,_PQ.mesh.get_elements()[2])
            for metric in ("volume","minEdge","maxEdge")
                expected=metric=="volume" ? expected_volume :
                    [kinds[tag]==15 ? (metric=="minEdge" ? 1e25 : 0.0) :
                     kinds[tag] in (1,8) || metric=="minEdge" ? .5 : sqrt(.5) for tag in tags]
                @test _PQ.mesh.get_element_qualities(tags,metric)≈expected atol=2e-11 rtol=0
                for (position,tag) in enumerate(tags)
                    @test only(_PQ.mesh.get_element_qualities([tag],metric))≈expected[position] atol=2e-11 rtol=0
                end
            end
            @test (_PQ.mesh.get_nodes(),_PQ.mesh.get_elements())==before
            pointtags,_=_PQ.mesh.get_elements_by_type(15)
            @test _PQ.mesh.get_element_qualities(pointtags,"volume")==fill(1.0,4)
            @test _PQ.mesh.get_element_qualities(pointtags,"minEdge")==fill(1e25,4)
            @test _PQ.mesh.get_element_qualities(pointtags,"maxEdge")==zeros(4)
            lines,_=_PQ.mesh.get_elements_by_type(quadratic ? 8 : 1)
            @test _PQ.mesh.get_element_qualities(lines,"volume")≈fill(.5,length(lines)) atol=2e-11 rtol=0
            for absent in (quadratic ? (1,2) : (8,9))
                emptytags,_=_PQ.mesh.get_elements_by_type(absent)
                @test isempty(emptytags)
                @test _PQ.mesh.get_element_qualities(emptytags,"volume")==Float64[]
            end
            records=Dict(tag=>_PQ.mesh.get_element(tag) for tag in tags)
            before_nodes=_PQ.mesh.get_nodes();before_elements=_PQ.mesh.get_elements()
            new=reverse(tags).+UInt64(1000)
            _PQ.mesh.renumber_elements(tags,new)
            @test _PQ.mesh.get_element_qualities(new,"volume")==volumes
            @test _PQ.mesh.get_nodes()==before_nodes
            @test _PQ.mesh.get_elements()[1]==before_elements[1]
            @test _PQ.mesh.get_elements()[3]==before_elements[3]
            for (old,tag) in zip(tags,new);@test _PQ.mesh.get_element(tag)==records[old];end
            @test _PQ.mesh.get_max_element_tag()==maxima[2]+UInt64(1000)
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

@testset "Actual Point and Line quality method and linear allocations" begin
    allocations=Int[]
    for count in (1000,2000,4000)
        mesh=_pq_scale_fixture(count);tags=UInt64.(1:2count)
        result=_PQ.mesh_element_qualities(mesh,tags,"volume")
        @test result==ones(2count)
        @test length(result)==2count
        @test Base.summarysize(result)<=64+8length(result)
        _PQ.mesh_element_qualities(mesh,tags,"volume");GC.gc()
        push!(allocations,minimum(@allocated(_PQ.mesh_element_qualities(mesh,tags,"volume")) for _ in 1:3))
    end
    for index in 2:3;@test allocations[index]<=2.5allocations[index-1]+262144;end
    @info "Point and Line quality original native measures" counts=(1000,2000,4000) allocations=Tuple(allocations)
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
