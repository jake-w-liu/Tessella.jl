module HybridWinslowAllocationChecks
using Test,Tessella
const W=Tessella.Model.ModelNativeHybridWinslow

function fixture(count)
    points=((1.,0.,0.),(0.,1.,0.),(0.,0.,1.),(0.,0.,0.))
    normals=((1.,0.,0.),(0.,1.,0.),(0.,0.,1.),(-1.,-1.,-1.))
    fill(W.PatchTet(points,normals,4),count)
end

function independent_energy(point,count,epsilon)
    determinant=1.0-sum(point)
    root=hypot(determinant,epsilon)
    chi=determinant>=0 ? (determinant+root)/2 :
        epsilon^2/(2*(root-determinant))
    count*((3-2sum(point)+3sum(abs2,point))/chi^(2/3)+
        (determinant^2+1)/chi)
end

measure(f)=(f();f();@allocated f())

@testset "Winslow tuple state and warmed allocations" begin
    for count in (1,12),epsilon in (1.e-3,1.e-6)
        patch=fixture(count)
        for point in ((0.08,0.1,-0.05),(0.8,0.5,0.4))
            energy,gradient,determinant=W.energy_gradient(point,patch,epsilon)
            @test energy≈independent_energy(point,count,epsilon) rtol=2e-14
            @test determinant≈1.0-sum(point) atol=3eps(Float64)
            for axis in 1:3
                delta=1.e-6
                plus=Base.setindex(point,point[axis]+delta,axis)
                minus=Base.setindex(point,point[axis]-delta,axis)
                difference=(independent_energy(plus,count,epsilon)-
                    independent_energy(minus,count,epsilon))/(2delta)
                @test gradient[axis]≈difference rtol=2e-7
            end
            if Base.JLOptions().compile_enabled!=3
                @test measure(()->W.energy_gradient(point,patch,epsilon))==0
            end
        end
        point=(0.08,0.1,-0.05)
        initial=W.energy_gradient(point,patch,epsilon)
        direction=map(-,initial[2])
        trial,result=W.line_search(point,patch,epsilon,direction,initial,0.01)
        @test result[1]<initial[1]
        @test result[3]>0.0
        optimal,minimum=W.minimize(point,patch,epsilon)
        @test all(value->abs(value)<1.e-4,optimal)
        @test minimum[1]≈5count rtol=3epsilon^2 atol=count*1.e-8
        @test sqrt(sum(abs2,minimum[2]))<1.e-4
        if Base.JLOptions().compile_enabled!=3
            @test measure(()->W.line_search(point,patch,epsilon,direction,initial,0.01))==0
            # L-BFGS owns a bounded history and alpha buffers; tuple updates add no boxes.
            @test measure(()->W.minimize(point,patch,epsilon))<=4096
        end
    end
end
end
