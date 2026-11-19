module OCCBooleanFaceMeasureChecks
using Test,Tessella
const A=Tessella.API
const M=Tessella.Model

# `occ.get_mass(2, …)` used the carrier's whole chart for every
# sphere/cylinder/cone face. A face whose wire bounds only part of the chart —
# a Boolean kept piece — must integrate over its wire instead; the trimmed
# cases are checked on real Boolean results in the sphere-pole and
# clipped-cylinder tests. Here: the test that tells the two apart, and the
# closed forms it must keep for whole primitives.
@testset "OCC quadric faces: chart-window wires keep their closed forms" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        occ.add_sphere(0,0,0,1.5)
        occ.add_cylinder(3,0,0,0,0,2,0.5)
        occ.add_cone(6,0,0,0,0,1,0.5,0.2)
        kernel=A.CURRENT[].occ_kernel
        expected=Dict(:sphere=>4π*1.5^2,:cylinder=>2π*0.5*2,
                      :cone=>π*0.7*hypot(1,0.3))
        seen=Symbol[]
        for ((d,tag),id) in kernel.bindings
            d==2 || continue
            shape=kernel.shapes[id];rec=shape.geometry.geometry
            kind=get(rec,:occ,nothing)
            haskey(expected,kind) || continue
            push!(seen,kind)
            window=M._occ_quadric_window(rec)
            @test M._occ_quadric_window_wire(kernel,shape,rec,window...)
            @test occ.get_mass(2,tag)≈expected[kind] rtol=1e-14
            # a narrower window is not the wire's: the face reads as trimmed
            u0,u1,v0,v1=window
            @test !M._occ_quadric_window_wire(kernel,shape,rec,
                                              u0,u1,v0,(v0+v1)/2)
        end
        @test sort!(seen)==[:cone,:cylinder,:sphere]
    finally
        A.finalize()
    end
end
end
