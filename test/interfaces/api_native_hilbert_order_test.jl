module NativeHilbertIntervalTests
using Tessella,Test
const K=Tessella.Model.ModelNativeInitialTet
const H=K.HilbertIntervals

# A mathematical index oracle, independent of the partition engine. Its
# state update is the geometric Gray-code entry/exit construction in
# Hamilton's Compact Hilbert Indices, not an unstable partition algorithm.
function index_for_lattice(point,entry,exit_axis,bits)
    index=0
    for level in bits-1:-1:0
        cell=sum(((point[axis]>>level)&1)<<(axis-1) for axis in 1:3)
        shift=mod(exit_axis+1,3)
        reflected=((cell⊻entry)>>shift)|((cell⊻entry)<<(3-shift))
        reflected&=7
        rank=reflected⊻(reflected>>1)⊻(reflected>>2)
        index=(index<<3)|rank
        entry_word=rank==0 ? 0 : (2*((rank-1)÷2))
        entry_word=entry_word⊻(entry_word>>1)
        entry=entry⊻(((entry_word<<shift)|(entry_word>>(3-shift)))&7)
        turn=rank==0 ? 0 : mod(trailing_ones(UInt(isodd(rank) ? rank : rank-1)),3)
        exit_axis=mod(exit_axis+turn+1,3)
    end
    return index
end
function check_lattice(bits,entry,exit_axis,permutation)
    side=1<<bits
    lattice=NTuple{3,Int}[(x,y,z) for z in 0:side-1 for y in 0:side-1 for x in 0:side-1]
    points=NTuple{3,Float64}[ntuple(axis->(lattice[row][axis]+0.5)/side,3) for row in permutation]
    before=copy(points)
    expected=sortperm([index_for_lattice(lattice[row],entry,exit_axis,bits) for row in permutation])
    order=Int32.(1:length(points))
    H.sort_cell!(order,points,1,length(order),H.Frame(entry,exit_axis),(0.,0.,0.),(1.,1.,1.),0)
    @test order==expected
    @test isequal(points,before)
    @test sort(order)==Int32.(1:length(points))
end
@testset "Hilbert partitions follow independent geometric lattice indices" begin
    for entry in 0:7,axis in 0:2
        check_lattice(1,entry,axis,[8,3,5,1,7,2,6,4])
        check_lattice(2,entry,axis,[mod(17i+11,64)+1 for i in 0:63])
    end
    # Equal coordinates remain separately owned rows, and the final leaf
    # threshold keeps a two-row sequence in its supplied order.
    duplicate=[(0.,-0.,0.),(0.,-0.,0.)]
    before=copy(duplicate)
    order=Int32[2,1]
    H.sort_cell!(order,duplicate,1,2,H.Frame(0,0),(-1.,-1.,-1.),(1.,1.,1.),0)
    @test order==Int32[2,1]
    @test isequal(duplicate,before)
    result=K.hilbert(duplicate,2)
    result[1]=2
    @test isequal(duplicate,before)
end
function resource_points(n)
    return NTuple{3,Float64}[(Float64(mod(13i,1009))/1009,
        Float64(mod(29i,1013))/1013,Float64(mod(43i,1019))/1019) for i in 1:n]
end
function partition_allocation(order,points)
    H.sort_cell!(order,points,1,length(order),H.Frame(0,0),(0.,0.,0.),(1.,1.,1.),0)
    return @allocated H.sort_cell!(order,points,1,length(order),H.Frame(0,0),
        (0.,0.,0.),(1.,1.,1.),0)
end
@testset "Hilbert ordering owns only its retained row vector" begin
    points=resource_points(128);order=Int32.(1:128)
    @test partition_allocation(order,points)==0
    previous=0;measurements=NamedTuple[]
    # Fixed before measurement. A fresh output row vector is required; the
    # allowance rejects a per-partition collection while tolerating wrappers.
    for n in (1000,2000,4000)
        points=resource_points(n);before=copy(points)
        result=K.hilbert(points,n)
        @test sort(result)==Int32.(1:n)
        @test isequal(points,before)
        GC.gc();bytes=@allocated K.hilbert(points,n)
        @test bytes<=8n+65_536
        previous==0 || @test bytes<=2.15previous+65_536
        push!(measurements,(points=n,bytes=bytes,retained_bytes=sizeof(result)))
        previous=bytes
    end
    println("HILBERT_INTERVAL_ALLOCATIONS ",measurements)
end
end
