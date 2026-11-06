module PostViewTranslatedLineChecks
using Test,Tessella
const F=Tessella.SizeField
function line_field(coordinates,values)
    F.PostViewField(coordinates,values;lines=reshape(Int32[1,2],2,1),
        crop_negative=false,use_closest=false)
end
@testset "plain PostView lines preserve physical query offsets" begin
    coordinates=[-1e16 0.;0. 0.;0. 0.]
    values=[1e16,0.]
    before_coordinates=copy(coordinates);before_values=copy(values)
    for order in ([1,2],[2,1])
        field=line_field(coordinates[:,order],values[order])
        stored_coordinates=copy(field.coords);stored_values=copy(field.values)
        for query in (-.25,-1.,-5e15,-1e16+2)
            @test F.field_value(field,query,0.,0.)==-query
        end
        @test isequal(field.coords,stored_coordinates)
        @test isequal(field.values,stored_values)
    end
    @test coordinates==before_coordinates
    @test values==before_values
    vectors=line_field(coordinates,[1e16 0.;2e16 0.;3e16 0.])
    @test F.field_value(vectors,-.25,0.,0.)≈sqrt(14)/4
    # Both endpoint differences lose the small midpoint shift; the existing
    # exact interpolation must receive the original physical query.
    wide=line_field([-1e308 1e308;0. 0.;0. 0.],[1e308,-1e308])
    @test F.field_value(wide,1.,0.,0.)==-1.
    # Even the normalized endpoint weight underflows for this finite result.
    tiny=line_field([-1e308 0.;0. 0.;0. 0.],[1e308,0.])
    @test F.field_value(tiny,-nextfloat(0.),0.,0.)==nextfloat(0.)
end
end
