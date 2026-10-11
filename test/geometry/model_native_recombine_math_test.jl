module ModelNativeRecombineMathTests
using Test,Tessella,SHA
const M=Tessella.Model
const G=Tessella.GmshLibm
const CORPUS=joinpath(@__DIR__,"..","artifacts","model_native_recombine_cost_primary.tsv")
floatbits(text)=reinterpret(Float64,parse(UInt64,text;base=16))
function state(shear)
    world=[0.0 1.0 (1.0+shear) shear;0.0 0.0 1.0 1.0;0.0 0.0 0.0 0.0]
    return M._NativeSurfaceRecombine(world,world[1:2,:],fill((2,Int32(1)),4),
        UInt64[10003,10010,10017,10024],NTuple{3,Int32}[(1,2,3),(1,3,4)],
        UInt64[30001,30002],NTuple{4,Int32}[],UInt64[],trues(4),UInt64(30002),identity,
        "pinned native producer integer pair-cost control")
end
function check_pair(shear,eta,cost)
    source=state(shear);before=copy(source.world)
    pairs,periodic=M._native_recombine_pairs!(source,source.node_tags,nothing)
    pair=only(pairs)
    @test reinterpret(UInt64,pair.eta)==reinterpret(UInt64,eta)
    @test pair.cost==cost
    @test source.high_water==30003
    @test periodic==[(1,2),(1,2)]
    @test source.world==before && source.triangles==[(1,2,3),(1,3,4)]
    return nothing
end
if Sys.iswindows() && Sys.ARCH===:x86_64
@testset "Native Surface producer angular and integer weighted costs" begin
    @test bytes2hex(sha256(read(CORPUS)))=="dd5438d9842239c5d81db724285600001469d1de25a64f704305fd22cd57bc61"
    saved=G._win_x87_control_word();selection=G._win_cli_precision()
    rows=0;controls=Tuple{UInt16,Float64,Float64,Int}[]
    try
        task_local_storage(G._WIN_CLI_PRECISION_KEY,false) do
            current=saved
            for line in Iterators.drop(eachline(CORPUS),1)
                mode,index,shear_text,eta_text,exp_text,cost_text=split(line,'\t')
                requested=UInt16((saved&~0x0f00)|parse(UInt16,mode;base=16))
                if current!=requested;G._win_x87_set_control_word!(requested);current=requested;end
                shear=floatbits(shear_text);eta=floatbits(eta_text);cost=parse(Int,cost_text)
                check_pair(shear,eta,cost)
                @test reinterpret(UInt64,G._gm87_exp(-eta))==parse(UInt64,exp_text;base=16)
                @test G._win_x87_control_word()==requested
                index=="145" && push!(controls,(requested,shear,eta,cost))
                rows+=1
            end
        end
        @test rows==37920
        @test length(controls)==12
        # The direct API inherits hardware precision. CLI math selects PC64
        # for each native leaf and restores every caller control-word bit.
        for (requested,shear,_,_) in controls
            G._win_x87_set_control_word!(requested)
            effective64=UInt16((requested&~0x0300)|0x0300)
            expected=only(control for control in controls if control[1]==effective64)
            G._with_win_cli_precision() do
                check_pair(shear,expected[3],expected[4])
                @test G._win_x87_control_word()==requested
                G._with_win_cli_precision() do
                    check_pair(shear,expected[3],expected[4])
                    @test G._win_x87_control_word()==requested
                end
            end
            @test G._win_cli_precision()==selection
        end
    finally
        G._win_x87_set_control_word!(saved)
    end
    @test G._win_x87_control_word()==saved
    @test G._win_cli_precision()==selection
end
end
end
