module AggregateCaseLoopRegressionTests
using Test

# Exercise the actual aggregate case loop with controlled process/mesh results.
# This verifies failure propagation without invoking another mesher from Pkg.test.
const SOURCE=joinpath(@__DIR__,"..","..","validation","run_all.jl")
const TEXT=replace(read(SOURCE,String),"\r\n"=>"\n")
const START=first(findfirst("for c in cases\n",TEXT))
const STOP=first(findnext("\n# ── enclosure coax junction:",TEXT,START))
const BODY=TEXT[START:prevind(TEXT,STOP)]
const HERE="unused-fault-injection-directory"
const haveg=true
const rows=String[]
const cases=NamedTuple[]
const REFERENCE_OK=Ref(true)
const REFERENCE=Ref((nnodes=4,ntets=1,volume=1.0))
const TESSELLA=Ref((nnodes=4,ntets=1,volume=1.0))

surface_stub()=nothing
relerr(value,reference)=abs(value-reference)/abs(reference)
tessella_fill(surface;smooth=false)=(TESSELLA[],0.0)
run_gmsh(input,output)=(ok=REFERENCE_OK[],seconds=0.0)
gmsh_metrics(output)=REFERENCE[]
fmtrow(case,tool,mesh,seconds,reference)="$case/$tool/$(mesh.volume)"

function actual_loop(;ok=true,reference=(nnodes=4,ntets=1,volume=1.0),
                     tessella=(nnodes=4,ntets=1,volume=1.0),model=1.0,truevolume=1.0)
    REFERENCE_OK[]=ok;REFERENCE[]=reference;TESSELLA[]=tessella
    empty!(rows);empty!(cases)
    push!(cases,(dir="injected_flat",geo="fixture.geo",surf=surface_stub,
                 model_vol=model,true_vol=truevolume,note="controlled volume oracle"))
    try
        redirect_stdout(devnull) do
            include_string(@__MODULE__,BODY,"aggregate-case-loop-regression.jl")
        end
        return nothing
    catch error
        return error isa LoadError ? error.error : error
    end
end

@testset "Aggregate reference failures propagate" begin
    @test actual_loop()===nothing
    @test length(rows)==2
    failure=actual_loop(ok=false)
    @test failure isa AssertionError
    @test occursin("required reference meshing failed",sprint(showerror,failure))
    @test length(rows)==1
    for volume in (0.0,-1.0,NaN,Inf)
        @test actual_loop(reference=(nnodes=4,ntets=1,volume=volume)) isa AssertionError
        @test actual_loop(tessella=(nnodes=4,ntets=1,volume=volume)) isa AssertionError
    end
    for (nodes,tets) in ((0,1),(4,0))
        @test actual_loop(reference=(nnodes=nodes,ntets=tets,volume=1.0)) isa AssertionError
        @test actual_loop(tessella=(nnodes=nodes,ntets=tets,volume=1.0)) isa AssertionError
    end
end

@testset "Exact flat-volume oracle and curved reporting remain distinct" begin
    @test actual_loop(reference=(nnodes=4,ntets=1,volume=17.0)) isa AssertionError
    @test actual_loop(reference=(nnodes=4,ntets=1,volume=1.0+5e-7))===nothing
    @test actual_loop(reference=(nnodes=4,ntets=1,volume=1.0+2e-6)) isa AssertionError
    @test actual_loop(tessella=(nnodes=4,ntets=1,volume=1.0+2e-6)) isa AssertionError
    # The inscribed Tessella model and true curved solid have different volumes.
    @test actual_loop(model=0.9,truevolume=1.0,
        tessella=(nnodes=4,ntets=1,volume=0.9),
        reference=(nnodes=4,ntets=1,volume=0.8))===nothing
    @test length(rows)==2
    # A sphere-like case has no independent polyhedron volume value.
    @test actual_loop(model=NaN,truevolume=1.0,
        tessella=(nnodes=4,ntets=1,volume=0.9),
        reference=(nnodes=4,ntets=1,volume=0.8))===nothing
    @test actual_loop(model=NaN,truevolume=1.0,
        tessella=(nnodes=4,ntets=0,volume=0.0)) isa AssertionError
end
end # module AggregateCaseLoopRegressionTests
