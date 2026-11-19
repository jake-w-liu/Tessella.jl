module APIMeshAdaptScalarGradingChecks
using Test, Tessella
const A=Tessella.API
const M=Tessella.Model
const F=Tessella.SizeField

# `Mesh.Algorithm=1` (MeshAdapt) certifies the boundary Curves of its faces
# against the face's own size. When an emitted Curve cell failed that
# certificate the Curve was regraded under `_MetricScaledAniso(inner, scale)`
# with `inner` the selected background field — which is `nothing` without a
# field and a scalar field otherwise, so every such model threw
# `MethodError: no method matching _MetricScaledAniso(::Nothing, ::Float64)`
# (a unit OCC rectangle at size 0.25 was enough: its side cells sit on the
# certificate bound, so two of its Curves are regraded).
@testset "MeshAdapt regrades Curves without a tensor background" begin
    @testset "tightened grading field" begin
        certificates=Tuple{Int,Any}[
            (1,M._MetricIsotropicSize((x,y,z)->0.5)),
            (2,M._MetricIsotropicSize((x,y,z)->0.2+x))]
        # No selection, or a scalar one: the certificates are the base, the
        # tightest of them at each point, divided by sqrt(scale).
        for inner in (nothing,F.ConstantSize(3.0))
            grading=M._model_curve_tightened_grading(inner,certificates,4.0)
            @test grading isa F.AbstractAnisoField
            @test F.field_value(grading,0.0,0.0,0.0)≈0.1
            @test F.field_value(grading,1.0,0.0,0.0,(1,7))≈0.25
        end
        # A tensor selection is tightened itself, as before.
        tensor=M._MetricIsotropicSize((x,y,z)->1.0)
        grading=M._model_curve_tightened_grading(tensor,certificates,4.0)
        @test grading isa M._MetricScaledAniso
        @test F.field_value(grading,0.0,0.0,0.0)≈0.5
    end

    @testset "OCC rectangle whose Curves are regraded" begin
        A.initialize()
        try
            A.model.occ.add_rectangle(0,0,0,1,1);A.model.occ.synchronize()
            A.option("Mesh.MeshSizeMin",0.0)
            A.option("Mesh.MeshSizeMax",0.25);A.option("Mesh.MeshSizeMin",0.25)
            A.option("Mesh.Algorithm",1)
            A.mesh.generate(2)
            tags,xyz,_=A.mesh.get_nodes()
            row=Dict(tag=>index for (index,tag) in enumerate(tags))
            point(tag)=(xyz[3row[tag]-2],xyz[3row[tag]-1],xyz[3row[tag]])
            types,_,nodes=A.mesh.get_elements(2)
            @test types==[2]
            area=0.0
            for k in 1:3:length(nodes[1])
                a,b,c=point(nodes[1][k]),point(nodes[1][k+1]),point(nodes[1][k+2])
                area+=((b[1]-a[1])*(c[2]-a[2])-(c[1]-a[1])*(b[2]-a[2]))/2
            end
            @test abs(area)≈1.0 atol=1e-12
            @test all(tag->point(tag)[3]==0.0,tags)
        finally
            A.finalize()
        end
    end
end
end
