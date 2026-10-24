# `plugin("Distance")` — Gmsh 4.15.2 `Distance.cpp` `DistanceType <= 0` arm:
# every mesh vertex's unsigned distance to the boundary entities' LINE and
# TRIANGLE elements, emitted as one scalar record per element on the
# max-dimension entities in a list view named `distance`. All Physical*
# options 0 targets dim-(maxDim-1) entities; a `Physical*` id targets the
# entities of that physical group in the matching dimension. On a unit
# square the boundary segments tile the four edges exactly, so every
# vertex's field value is analytically `min(x, 1-x, y, 1-y)`; restricting
# to the bottom+right physical curves gives `min(y, 1-x)`.
module PluginDistanceChecks
using Test,Tessella
const A=Tessella.API

function _square(;physical=false)
    A.model.add("sq")
    for (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))
        A.model.occ.add_point(x,y,0.)
    end
    for i in 1:4; A.model.occ.add_line(i,mod(i,4)+1); end
    A.model.occ.add_curve_loop([1,2,3,4])
    A.model.occ.add_plane_surface([1])
    A.model.occ.synchronize()
    physical && A.model.add_physical_group(1,[1,2];tag=1)
    A.mesh.generate(2)
    return nothing
end

function _distance_view()
    A.plugin.run("Distance")
    tags=A.view.get_tags()
    @test !isempty(tags)
    types,data=A.view.get_list_data(last(tags))
    st=findfirst(==("ST"),types)
    @test st!==nothing
    d=data[st]
    # Each ST record: 3×x, 3×y, 3×z, then 3 nodal distances.
    return reshape(d,12,length(d)÷12)
end

@testset "Distance emits the exact analytic boundary-distance field" begin
    A.initialize()
    try
        _square()
        block=_distance_view()
        @test size(block,2)>0
        for e in axes(block,2)
            for node in 1:3
                x,y,z=block[node,e],block[3+node,e],block[6+node,e]
                d=block[9+node,e]
                @test z==0.0
                @test d≈min(x,1-x,y,1-y) atol=1e-10
            end
        end
    finally
        A.finalize()
    end

    A.initialize()
    try
        _square(physical=true)
        # `PhysicalLine=1` restricts the boundary set to physical curves 1,2
        # — the bottom and right edges — so the field reads min(y, 1-x).
        A.plugin.set_number("Distance","PhysicalLine",1.0)
        block=_distance_view()
        @test size(block,2)>0
        for e in axes(block,2)
            for node in 1:3
                x,y=block[node,e],block[3+node,e]
                @test block[9+node,e]≈min(y,1-x) atol=1e-10
            end
        end
    finally
        A.finalize()
    end

    A.initialize()
    try
        _square()
        # MinScale/MaxScale rescale to [min,max] when both are positive.
        A.plugin.set_number("Distance","MinScale",1.0)
        A.plugin.set_number("Distance","MaxScale",3.0)
        block=_distance_view()
        vals=vec(block[10:12,:])
        @test minimum(vals)≈1.0 atol=1e-12
        @test maximum(vals)≈3.0 atol=1e-12
        # MinScale alone shifts the raw field.
        A.plugin.set_number("Distance","MaxScale",0.0)
        block=_distance_view()
        vals=vec(block[10:12,:])
        @test minimum(vals)≈1.0 atol=1e-12
        for e in axes(block,2),node in 1:3
            x,y=block[node,e],block[3+node,e]
            @test block[9+node,e]≈1.0+min(x,1-x,y,1-y) atol=1e-10
        end
    finally
        A.finalize()
    end
end

end # module
