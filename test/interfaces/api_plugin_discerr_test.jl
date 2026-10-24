# `plugin("DiscretizationError")` — Gmsh 4.15.2 `DiscretizationError.cpp`:
# every surface element is supersampled (quadrangles on an n×n grid over
# [-1+5tol,~1-5tol]^2, triangles on the (n+1)n/2-point simplicial grid),
# each sample evaluates the element's `pnt` and is projected onto the
# owning face; the (n-1)^2 sub-elements are emitted as SQ/ST records
# carrying the unsigned projection distance. On a unit sphere the
# projection is radial, so every value is exactly `R - |p|`.
module PluginDiscretizationErrorChecks
using Test,Tessella
const A=Tessella.API

@testset "DiscretizationError supersamples surface elements onto the CAD face" begin
    A.initialize()
    try
        A.model.add("sph")
        A.model.occ.add_sphere(0.,0.,0.,1.0)
        A.model.occ.synchronize()
        A.mesh.generate(2)
        A.plugin.set_number("DiscretizationError","SuperSamplingNodes",10.0)
        A.plugin.run("DiscretizationError")
        tags=A.view.get_tags()
        @test !isempty(tags)
        types,data=A.view.get_list_data(last(tags))
        st=findfirst(==("ST"),types)
        @test st!==nothing
        d=data[st]
        ne=length(d)÷12
        @test ne>0
        # Every record: 3x,3y,3z then 3 distances — on a unit sphere the
        # projection is radial, so v == 1 - |p| within float noise.
        for e in 1:ne
            base=12(e-1)
            for k in 0:2
                x,y,z=d[base+1+k],d[base+4+k],d[base+7+k]
                v=d[base+10+k]
                @test v≈1-hypot(x,y,z) atol=1e-6
            end
        end
        # Interior samples sag below the sphere — nonzero error that stays
        # below the element sag bound; corner samples are ~0.
        vals=d[10:12:end]
        @test maximum(vals)>1e-4
        @test maximum(vals)<0.1
        @test minimum(vals)<1e-6
        # With SuperSamplingNodes=3 each element emits (3-1)^2 = 4 records.
        A.plugin.set_number("DiscretizationError","SuperSamplingNodes",3.0)
        A.plugin.run("DiscretizationError")
        types2,data2=A.view.get_list_data(last(A.view.get_tags()))
        st2=findfirst(==("ST"),types2)
        @test st2!==nothing
        @test length(data2[st2])÷12==ne÷(10-1)^2*(3-1)^2
    finally
        A.finalize()
    end

    A.initialize()
    try
        # A planar face projects to itself — every value ~0.
        A.model.add("sq")
        for (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))
            A.model.occ.add_point(x,y,0.)
        end
        for i in 1:4; A.model.occ.add_line(i,mod(i,4)+1); end
        A.model.occ.add_curve_loop([1,2,3,4])
        A.model.occ.add_plane_surface([1])
        A.model.occ.synchronize()
        A.mesh.generate(2)
        A.plugin.run("DiscretizationError")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        st=findfirst(==("ST"),types)
        @test st!==nothing
        vals=data[st][10:12:end]
        @test all(<(1e-6),vals)
    finally
        A.finalize()
    end
end

end # module
