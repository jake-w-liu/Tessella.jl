module PartitionParallelDispatchChecks
using Test,Tessella
const A=Tessella.API
const E=Tessella.Elements

struct WatchedIndices <: AbstractVector{Int}
    values::Vector{Int}
    active::Threads.Atomic{Int}
    peak::Threads.Atomic{Int}
end
WatchedIndices(values)=WatchedIndices(values,Threads.Atomic{Int}(0),Threads.Atomic{Int}(0))
Base.size(indices::WatchedIndices)=size(indices.values)
function Base.getindex(indices::WatchedIndices,index::Int)
    current=Threads.atomic_add!(indices.active,1)+1
    while true
        previous=indices.peak[]
        current<=previous && break
        Threads.atomic_cas!(indices.peak,previous,current)==previous && break
    end
    try
        sleep(.001)
        return indices.values[index]
    finally
        Threads.atomic_sub!(indices.active,1)
    end
end

function callback_probe(workers;count=8)
    active=Threads.Atomic{Int}(0)
    peak=Threads.Atomic{Int}(0)
    seen=zeros(Int,count)
    owners=Vector{Task}(undef,count)
    A._partition_parallel_for(7:count+6,workers) do index
        current=Threads.atomic_add!(active,1)+1
        while true
            previous=peak[]
            current<=previous && break
            Threads.atomic_cas!(peak,previous,current)==previous && break
        end
        try
            slot=index-6
            seen[slot]+=1
            owners[slot]=current_task()
            # Yield with the callback active so overlapping scopes are counted.
            # There is no elapsed-time or physical-thread identity assertion.
            sleep(.01)
        finally
            Threads.atomic_sub!(active,1)
        end
    end
    expected=min(max(workers,1),Threads.threadpoolsize(:default),count)
    @test active[]==0
    @test all(==(1),seen)
    @test 1<=peak[]<=expected
    @test length(unique(owners))==expected
    return nothing
end

@testset "partition callbacks honor the requested worker cap and join" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        pool=Threads.threadpoolsize(:default)
        for requested in (0,1,2,pool+1)
            A.option("General.NumThreads",requested)
            workers=A._partition_workers()
            @test workers==(requested==0 ? pool : requested)
            callback_probe(workers)
        end
        callback_probe(pool+1;count=1)
        callback_probe(pool+1;count=3)
    finally
        A.finalize()
    end
    called=Ref(false)
    A._partition_parallel_for(1:0,2) do _
        called[]=true
    end
    @test !called[]
    active=Threads.Atomic{Int}(0)
    failure=try
        A._partition_parallel_for(1:8,2) do index
            Threads.atomic_add!(active,1)
            try
                index==1 && error("partition callback probe")
                sleep(.01)
            finally
                Threads.atomic_sub!(active,1)
            end
        end
        nothing
    catch error
        error
    end
    @test failure!==nothing
    @test occursin("partition callback probe",sprint(showerror,failure))
    @test active[]==0
end

@testset "mesh structure and quality scans honor the requested worker cap" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("General.NumThreads",2)
        workers=Tessella.MeshTypes._parallel_workers()
        @test workers==2
        expected=min(workers,Threads.threadpoolsize(:default))
        indices=WatchedIndices(collect(1:32))
        flags=zeros(Int8,8)
        Tessella.MeshTypes._scan_out_of_range!(flags,reshape(indices,4,8),32,workers)
        @test all(iszero,flags)
        @test indices.active[]==0
        @test 1<=indices.peak[]<=expected
        mesh=Tessella.Mesh([0. 1. 0. 0.;0. 0. 1. 0.;0. 0. 0. 1.];
            tets=reshape(Int32[1,2,3,4],4,1))
        # Repeated requests are part of the public quality-query contract.
        tags=WatchedIndices(ones(Int,8))
        qualities=zeros(8)
        Tessella.MeshElementQuality._scan_qualities!(qualities,mesh,tags,
            "volume",0,0,workers)
        @test qualities==fill(1/6,8)
        @test tags.active[]==0
        @test 1<=tags.peak[]<=expected
    finally
        A.finalize()
    end
end

function product(dimension,workers)
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.Renumber",0)
        A.model.add("partition dispatch")
        A.model.add_discrete_entity(dimension,11)
        A.model.add_discrete_entity(dimension,22)
        if dimension==3
            coordinates=[0. 1. 0. 0. 0.;0. 0. 1. 0. 0.;0. 0. 0. 1. -1.]
            A.mesh.add_nodes(3,11,11:15,vec(coordinates))
            A.mesh.add_elements_by_type(11,4,[501],[11,12,13,14])
            A.mesh.add_elements_by_type(22,4,[502],[11,13,12,15])
        else
            A.mesh.add_nodes(2,11,11:13,[0.,0.,0.,1.,0.,0.,1.,1.,0.])
            A.mesh.add_nodes(2,22,[14],[0.,1.,0.])
            A.mesh.add_elements_by_type(11,2,[501],[11,12,13])
            A.mesh.add_elements_by_type(22,2,[502],[11,13,14])
        end
        A.mesh.generate(dimension)
        A.option("Mesh.PartitionCreateTopology",1)
        A.option("Mesh.PartitionCreateGhostCells",1)
        A.option("General.NumThreads",workers)
        A.mesh.partition(2,[501,502],[1,2])
        mesh=A.mesh.get()
        @test Tessella.validate(mesh).ok
        @test !isempty(mesh.partition_data.ghost_elements)
        @test count(entity->entity.dim==dimension && entity.parent[1]==dimension &&
            entity.parent[2] in (11,22),values(mesh.partition_data.entities))==2
        return mktempdir() do directory
            path=joinpath(directory,"partition.msh")
            E.write_mixed_msh(path,mesh;version=4.1)
            read(path)
        end
    finally
        A.finalize()
    end
end

@testset "parallel partition preserves complete volume and surface products" begin
    # With -t4 these exercise actual concurrent callbacks. On a single-thread
    # runtime they still check the requested-parallel dispatch and its output.
    for dimension in (3,2)
        serial=product(dimension,1)
        parallel=product(dimension,3)
        @test parallel==serial
    end
end
end
