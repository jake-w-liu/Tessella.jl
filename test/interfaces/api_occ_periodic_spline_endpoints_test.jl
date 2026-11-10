module OCCPeriodicSplineEndpointChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

function circle_chart(parameter)
    width=2π/3
    index=clamp(floor(Int,parameter/width),0,2)
    lower=index*width
    q=tan(width/4)
    z=(2(parameter-lower)/width-1)*q
    angle=lower+width/2+2atan(z)
    first=4q/(width*(1+z*z))
    second=-first*first*z
    point=[cos(angle),sin(angle),0.]
    tangent=[-sin(angle),cos(angle),0.]
    return point,first.*tangent,second.*tangent.-(first*first).*point
end

@testset "Periodic splines evaluate their repeated active upper knot" begin
    A.initialize()
    try
        fixture=joinpath(@__DIR__,"..","fixtures","brep_native_periodic_circle.brep")
        @test A.model.occ.import_shapes(fixture)==[(1,1)]
        A.model.occ.synchronize()
        lo,hi=A.model.get_parametrization_bounds(1,1)
        @test lo==[0.] && hi==[6.28318530717959]
        @test hi[1]>2π # Preserve the writer-rounded BREP edge range.
        geometry=A.CURRENT[].curve_geometry[1].curve
        for parameter in (lo[1],.21,2π/3,4π/3,hi[1]-1e-7,2π,hi[1])
            expected=circle_chart(parameter)
            actual=M._model_occ_curve_d2(geometry,parameter)
            @test A.model.get_value(1,1,[parameter])≈expected[1] atol=8e-15
            @test A.model.get_derivative(1,1,[parameter])≈expected[2] atol=8e-15
            @test all(isfinite,A.model.get_second_derivative(1,1,[parameter]))
            for order in 1:3
                @test collect(actual[order])≈expected[order] atol=3e-14
            end
        end
        @test A.model.occ.get_mass(1,1)≈2π rtol=1e-11
        @test collect(A.model.occ.get_center_of_mass(1,1))≈zeros(3) atol=2e-14
        @test A.model.occ.get_matrix_of_inertia(1,1)≈
            [π,0.,0.,0.,π,0.,0.,0.,2π] atol=8e-14
        A.option("Mesh.MeshSizeMin",.4);A.option("Mesh.MeshSizeMax",.4)
        A.mesh.generate(1)
        tags,xyz,parameters=A.mesh.get_nodes(1,1)
        @test !isempty(tags) && length(tags)==length(parameters)
        @test A.model.get_value(1,1,parameters)≈xyz atol=2e-13
        families,labels,cells=A.mesh.get_elements(1,1)
        @test families==Int32[1] && !isempty(only(labels))
        incidence=Dict{UInt64,Int}()
        for node in only(cells)
            incidence[node]=get(incidence,node,0)+1
        end
        @test all(==(2),values(incidence))
        @test Tessella.validate(A.LAST_MESH[]).ok
        A.mesh.clear()
        A.model.occ.dilate([(1,1)],0.,0.,0.,1e-12,1e-12,1e-12)
        A.model.occ.synchronize()
        @test A.model.occ.get_mass(1,1)≈2π*1e-12 rtol=1e-11
        @test collect(A.model.occ.get_center_of_mass(1,1))≈zeros(3) atol=2e-26
        @test A.model.occ.get_matrix_of_inertia(1,1)≈
            1e-36.*[π,0.,0.,0.,π,0.,0.,0.,2π] rtol=1e-11 atol=1e-47
    finally
        A.finalize()
    end
end
end
