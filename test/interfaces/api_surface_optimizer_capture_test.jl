module SurfaceOptimizerCaptureChecks
using Test,Tessella
const A=Tessella.API
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)

function _boxed_slots(f)
    slots=Symbol[]
    for method in methods(f)
        method.module===A || continue
        info=Base.uncompressed_ast(method)
        for statement in info.code
            statement isa Expr && statement.head==:(=) || continue
            rhs=statement.args[2]
            (rhs isa Expr && rhs.head==:call && rhs.args[1]==GlobalRef(Core,:Box)) || continue
            push!(slots,info.slotnames[statement.args[1].id])
        end
    end
    return slots
end

function _fixture()
    call(A.initialize)
    call(A.option,"General.Terminal",0)
    for (tag,point) in enumerate(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)))
        call(A.model.add_point,point...;tag)
    end
    for tag in 1:4;call(A.model.add_line,tag,mod1(tag+1,4);tag);end
    call(A.model.add_curve_loop,collect(1:4);tag=1)
    call(A.model.add_plane_surface,[1];tag=1)
    coords=[0. 1. 1. 0. .35;0. 0. 1. 1. .4;0. 0. 0. 0. 0.]
    flat=call(A.model.get_parametrization,2,1,vec(coords))
    parameters=[collect(reshape(flat,2,:)[:,i]) for i in 1:5]
    cells=[(:tri,Int32[1,2,5]),(:tri,Int32[2,3,5]),
           (:tri,Int32[3,4,5]),(:tri,Int32[4,1,5])]
    return coords,parameters,cells,A.CURRENT[]
end

@testset "Surface optimizer closures retain typed local state" begin
    @test isempty(_boxed_slots(A._surface_optimize_good_move))
    @test isempty(_boxed_slots(A._surface_relocate_vertex!))
    coords,parameters,cells,model=_fixture()
    try
        original=copy(coords)
        uv_for=node->parameters[node]
        trial=call(A.model.get_parametrization,2,1,[.5,.5,0.])
        result=call(A._surface_optimize_good_move,coords,uv_for,5,collect(1:4),cells,model,1,trial)
        @test result isa Bool
        @test coords==original
        for stored in (true,false)
            coords.=original
            lookup=stored ? uv_for : node->node==5 ? nothing : parameters[node]
            result=call(A._surface_relocate_vertex!,coords,lookup,5,collect(1:4),cells,model,1)
            @test stored ? result isa Vector{Float64} : result===true
            @test coords[:,1:4]==original[:,1:4]
            @test all(isfinite,coords[:,5])
            @test coords[3,5]==0.
            @test 0. <= coords[1,5] <= 1. && 0. <= coords[2,5] <= 1.
            if stored
                @test call(A.model.get_value,2,1,result)==coords[:,5]
            end
        end
    finally
        call(A.finalize)
    end
end
end
