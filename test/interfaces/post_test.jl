using Test
using SHA
using Tessella
using Tessella.Post: View,view_value,add_plugin!,apply_plugin

function _post_view_sha(view::View)
    coordinate_bits=join(reinterpret(UInt64,vec(view.coords)),',')
    value_bits=join(reinterpret(UInt64,view.values),',')
    payload="$(ncodeunits(view.name)):$(view.name)|$(size(view.coords,1))×"*
            "$(size(view.coords,2))|$coordinate_bits|$value_bits"
    return bytes2hex(SHA.sha256(codeunits(payload)))
end

@testset "owned scalar post views and plugins" begin
    coords=Float64[0 1;0 0;0 0]
    values=[-1.0,2.0]
    view=View("samples",coords,values)
    @test view.name=="samples" && view.coords==coords && view.values==values
    @test view.coords!==coords && view.values!==values
    coords[1,1]=9.0;values[1]=9.0
    @test view.coords[1,1]==0.0 && view_value(view,1)==-1.0

    mesh=Mesh(Float64[0 1;0 0;0 0])
    mesh_view=View(SubString("xmeshx",2,5),mesh,[3,4])
    @test mesh_view.name=="mesh" && mesh_view.values==[3.0,4.0]
    @test _post_view_sha(view)==
          "56e766682618b76029640b47caa692205eb967a97438ee383eb057fc47cd96cd"

    absolute=apply_plugin("Abs",view)
    @test absolute isa View
    @test absolute.name=="samples_abs" && absolute.values==[1.0,2.0]
    @test absolute.coords==view.coords && absolute.coords!==view.coords
    @test view.values==[-1.0,2.0]
    zero_view=View("zero",zeros(3,2),[0.0,-0.0])
    @test apply_plugin("IsosurfaceZeroCount",zero_view)==2

    plugin_name=SubString("xPostTestDoublex",2,15)
    @test add_plugin!(plugin_name,v->View(v.name*"_double",v.coords,2 .* v.values))==
          "PostTestDouble"
    doubled=apply_plugin("PostTestDouble",view)
    @test doubled.values==[-2.0,4.0] && view.values==[-1.0,2.0]
    add_plugin!("PostTestDouble",v->sum(v.values))
    @test apply_plugin("PostTestDouble",view)==1.0
    add_plugin!("PostTestError",_ -> error("post plugin failure"))
    @test_throws ErrorException apply_plugin("PostTestError",view)
    @test_throws ArgumentError add_plugin!("",identity)
    @test_throws ArgumentError apply_plugin("NoSuchPostPlugin",view)

    @test_throws ArgumentError View("bad",zeros(2,1),[1.0])
    @test_throws ArgumentError View("bad",zeros(3,2),[1.0])
    @test_throws ArgumentError View("bad",reshape([NaN,0.0,0.0],3,1),[1.0])
    @test_throws ArgumentError View("bad",zeros(3,1),[Inf])
    empty_view=View("empty",zeros(3,0),Float64[])
    @test isempty(empty_view.values)
    @test_throws ArgumentError view_value(view,false)
    @test_throws ArgumentError view_value(view,0)
    @test_throws ArgumentError view_value(view,3)
    @test_throws ArgumentError view_value(view,big(typemax(Int))+1)
    @test isempty(Docs.undocumented_names(Tessella.Post;private=false))
end
