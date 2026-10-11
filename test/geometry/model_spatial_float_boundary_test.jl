module SpatialFloatBoundaryControlsV1
using Tessella,Test
const M=Tessella.Model
function square(z;offset=0.)
    m=GeoModel()
    for (tag,x,y) in ((1,0.,0.),(2,1.,0.),(3,1.,1.),(4,0.,1.))
        add_point!(m,offset+x,offset+y,z;tag)
    end
    for (tag,a,b) in ((1,1,2),(2,2,3),(3,3,4),(4,4,1))
        add_line!(m,a,b;tag)
    end
    add_curve_loop!(m,[1,2,3,4];tag=1)
    add_plane_surface!(m,[1];tag=1)
    return m
end
geometry_state(m)=(points=copy(m.points),curves=copy(m.curves),loops=deepcopy(m.loops),surfaces=deepcopy(m.surfaces))
@testset "Float64 distance projector after exact native Plane lift" begin
    for base in (0.,2.0^40,2.0^50)
        m=square(base;offset=base)
        add_point!(m,base+.5,base+.5,base+2.;tag=9)
        add_point!(m,base+3.,base+.5,base;tag=10)
        before=geometry_state(m)
        for (point,expected) in ((9,2.),(10,2.))
            d,pa,pb=M.model_distance(m,0,point,2,1)
            @test d isa Float64
            @test pa isa NTuple{3,Float64} && pb isa NTuple{3,Float64}
            @test d≈expected atol=1e-9
            reverse_d,reverse_a,reverse_b=M.model_distance(m,2,1,0,point)
            @test (reverse_d,reverse_a,reverse_b)==(d,pb,pa)
            @test hypot((pa.-pb)...)≈d atol=1e-9
        end
        @test geometry_state(m)==before
    end
end
end
