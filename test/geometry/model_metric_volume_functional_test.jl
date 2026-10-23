module ModelMetricVolumeFunctionalChecks
using Test,Tessella,LinearAlgebra
const M=Tessella.Model
const F=Tessella.SizeField
const A=Tessella.API

# Independent metric-edge audit of a tetrahedral fill. `_volume_metric_edge_field`
# scales the resolved metric by 1/_GMSH_TET_EDGE_OVER_RADIUS² (=1/4), so the
# fill criterion `metric_edge_length ≤ 1` on the wrapped field reads
# √(dᵀMd) ≤ 2 in raw input-metric units — the anisotropic restatement of the
# scalar `edge ≤ 2·lc` bound. Samples both endpoints and the midpoint of every
# tet edge exactly like `metric_edge_length`; `metric(p)` returns the raw
# (unscaled) 3×3 input metric at p. The fill runs `best_effort=true` — the
# same flag the scalar arm uses — so an edge whose star never becomes
# splittable is left unsplit and exempt from the strict postcondition; the
# audit therefore reports the violating-edge count alongside the maximum and
# the tests require both a tiny violating fraction and a bounded residual.
function metric_edge_audit(mesh,metric;bound=2.0)
    volumes=Float64[];edges=Set{Tuple{Int,Int}}()
    for vertices in eachcol(mesh.tets)
        a,b,c,d=(mesh.coords[:,vertex] for vertex in vertices)
        push!(volumes,det(hcat(b-a,c-a,d-a))/6)
        for i in 1:3,j in i+1:4
            push!(edges,minmax(Int(vertices[i]),Int(vertices[j])))
        end
    end
    maximum_metric=0.0;violating=0
    for (a,b) in edges
        first,last=mesh.coords[:,a],mesh.coords[:,b]
        displacement=last-first
        L=maximum(sqrt(dot(displacement,metric(point)*displacement))
                  for point in (first,last,(first+last)/2))
        maximum_metric=max(maximum_metric,L)
        L>bound+1e-9 && (violating+=1)
    end
    return (;volumes,maximum_metric,violating,total_edges=length(edges))
end

# Best-effort fill contract: nearly all edges satisfy the metric bound and
# every remnant stays within a bounded factor of it.
function check_metric_bound(checked;fraction=0.01,cap=2.5)
    @test checked.violating<=ceil(Int,fraction*checked.total_edges)
    @test checked.maximum_metric<=cap
end

function box_model(;extend=0)
    model=GeoModel()
    M.add_box!(model,0.,0.,0.,1.,1.,1.;tag=1)
    model.meshing.lc_extend_from_boundary=extend
    return model
end

@testset "tensor size_field drives the factory volume edge fill" begin
    # diag(64,1,1): strong direction x with h_x=1/8, weak directions y,z with
    # h=1. The scalar arm would collapse to h=1/8 everywhere (edges ≤ 1/4 in
    # every direction); the tensor arm bounds √(dᵀMd)≤2, i.e. x-edges ≤1/4
    # while y/z edges may span the whole unit box.
    model=box_model()
    mesh=mesh_model_volume(model,1;size_field=MathEvalAnisoField(m11="64",m22="1",m33="1"))
    checked=metric_edge_audit(mesh,point->Diagonal([64.,1.,1.]))
    @test validate(mesh).ok && size(mesh.tets,2)>0
    @test all(>(0),checked.volumes) && isapprox(sum(checked.volumes),1.;atol=1e-12)
    check_metric_bound(checked)

    # The isotropic collapse the tensor arm replaces: same unit box meshed by
    # the scalar field equal to the metric's smallest directional size.
    scalar=box_model()
    smesh=mesh_model_volume(scalar,1;size_field=ConstantSize(0.125))
    @test validate(smesh).ok
    @test size(mesh.tets,2)<size(smesh.tets,2)

    # A different strong direction produces a different fill, and satisfies
    # its own metric bound — the anisotropy actually transports.
    ymodel=box_model()
    ymesh=mesh_model_volume(ymodel,1;
        size_field=MathEvalAnisoField(m11="1",m22="64",m33="1"))
    ychecked=metric_edge_audit(ymesh,point->Diagonal([1.,64.,1.]))
    @test validate(ymesh).ok
    check_metric_bound(ychecked)
    @test !isequal(mesh.coords,ymesh.coords) || mesh.tets!=ymesh.tets
end

@testset "boundary sizes still cap a tensor volume fill" begin
    # With the default boundary-size extension the PostView boundary field is
    # composed into the tensor through the Loewner minimum (MinAnisoField):
    # the box boundary's isotropic sizes bound every direction. The fill stays
    # valid and satisfies the raw tensor bound (the composed metric is tighter
    # than the raw tensor, so the raw bound is implied).
    model=box_model(;extend=1)
    mesh=mesh_model_volume(model,1;size_field=MathEvalAnisoField(m11="64",m22="1",m33="1"))
    checked=metric_edge_audit(mesh,point->Diagonal([64.,1.,1.]))
    @test validate(mesh).ok && size(mesh.tets,2)>0
    @test all(>(0),checked.volumes) && isapprox(sum(checked.volumes),1.;atol=1e-12)
    check_metric_bound(checked)
end

@testset "spatially varying tensor grades the volume fill" begin
    # m11=(8x+1)²: h_x shrinks from 1 to 1/9 across the box while y,z stay at
    # h=1 — refinement must concentrate toward x=1.
    model=box_model()
    field=MathEvalAnisoField(m11="(8*x+1)^2",m22="1",m33="1")
    mesh=mesh_model_volume(model,1;size_field=field)
    checked=metric_edge_audit(mesh,point->Diagonal([(8*point[1]+1)^2,1.,1.]))
    @test validate(mesh).ok && size(mesh.tets,2)>0
    @test all(>(0),checked.volumes) && isapprox(sum(checked.volumes),1.;atol=1e-12)
    check_metric_bound(checked)
    right=count(vertices->sum(mesh.coords[1,vertices])/4>0.5,eachcol(mesh.tets))
    left=size(mesh.tets,2)-right
    @test right>left
end

@testset "isotropic-varying tensor keeps the local edge bound" begin
    # h(x)=1/(4x+1) in every direction on the unit box: the wrapped metric
    # certificate is |edge| <= 2*h(midpoint), so every edge's raw
    # sqrt(d'Md) must stay <= 2 — assert that per x-bin, with zero
    # violating edges for this smooth field, plus the density ordering the
    # fill actually produces (the right end may only coarsen to its bound;
    # a monotone count decrease across x is not required — the criterion is
    # an upper bound over an inherited initial-Delaunay floor).
    model=box_model()
    field=MathEvalAnisoField(m11="(4*x+1)^2",m22="(4*x+1)^2",m33="(4*x+1)^2")
    mesh=mesh_model_volume(model,1;size_field=field)
    @test validate(mesh).ok && size(mesh.tets,2)>0
    checked=metric_edge_audit(mesh,point->Diagonal(fill((4*point[1]+1)^2,3)))
    @test checked.violating==0 && checked.maximum_metric<=2.0
    # h halves from x=0 (h=1) to x=1 (h=1/5): the fine end must not carry
    # fewer cells than the coarse end.
    right=count(vertices->sum(mesh.coords[1,vertices])/4>0.5,eachcol(mesh.tets))
    left=size(mesh.tets,2)-right
    @test right>=left
    # The isotropic tensor arm and the equivalent scalar FunctionSize run
    # mathematically the same edge criterion; their fills must track the
    # same density — best-effort ordering can shift the count a few
    # percent, so compare with a generous tolerance.
    smodel=box_model()
    smesh=mesh_model_volume(smodel,1;
        size_field=FunctionSize((x,y,z)->1.0/(4x+1)))
    @test isapprox(size(mesh.tets,2),size(smesh.tets,2);rtol=0.15)
end

@testset "scalar size_callback keeps the scalar arm" begin
    # An installed size_callback is a scalar lc contract: the tensor field
    # must be scalarized through _volume_size_field (callback is invoked) and
    # the output remains a valid mesh.
    model=box_model()
    calls=Ref(0)
    model.meshing.size_callback=(dim,tag,x,y,z,lc)->(calls[]+=1;lc)
    mesh=mesh_model_volume(model,1;size_field=MathEvalAnisoField(m11="64",m22="1",m33="1"))
    @test calls[]>0
    @test validate(mesh).ok && size(mesh.tets,2)>0
end

@testset "API background MathEvalAniso reaches generate(3)" begin
    A.initialize()
    try
        A.model.add("metric volume session")
        A.model.add_box(0.,0.,0.,1.,1.,1.;tag=1)
        for (name,value) in (("General.Terminal",0),("Mesh.Renumber",0),
                ("Mesh.MeshSizeFromPoints",0),("Mesh.MeshSizeFromCurvature",0),
                ("Mesh.MeshSizeExtendFromBoundary",0),("Mesh.Smoothing",0))
            A.option(name,value)
        end
        field=A.mesh.field.add("MathEvalAniso")
        for (name,expression) in zip(("M11","M22","M33","M12","M13","M23"),
                                     ("64","1","1","0","0","0"))
            A.mesh.field.set_string(field,name,expression)
        end
        A.mesh.field.set_as_background_mesh(field)
        A.mesh.generate(3)
        @test A.mesh.get_element_types(3,1)==Int32[4]
        tags,xyz,_=A.mesh.get_nodes()
        coordinates=Dict(tag=>xyz[3i-2:3i] for (i,tag) in enumerate(tags))
        families,labels,nodes=A.mesh.get_elements(3,1)
        @test families==Int32[4] && !isempty(labels)
        cells=collect(eachcol(reshape(nodes[1],4,:)))
        volumes=[det(hcat(coordinates[cell[2]]-coordinates[cell[1]],
                          coordinates[cell[3]]-coordinates[cell[1]],
                          coordinates[cell[4]]-coordinates[cell[1]]))/6
                 for cell in cells]
        @test all(>(0),volumes) && isapprox(sum(volumes),1.;atol=1e-12)
        lengths=[begin
            first,last=coordinates[a],coordinates[b]
            displacement=last-first
            maximum(sqrt(dot(displacement,
                             Diagonal([64.,1.,1.])*displacement))
                    for point in (first,last,(first+last)/2))
        end for (a,b) in Set(minmax(cell[i],cell[j]) for cell in cells
                            for i in 1:3 for j in i+1:4)]
        # The session's _GmshBackgroundAniso wrapper intersects the tensor with
        # an isotropic model-lc cap, and the API boundary shell differs from
        # the direct-fill one, so a somewhat larger best-effort remnant share
        # is expected here than on the model-level arm.
        @test count(>(2.0+1e-9),lengths)<=ceil(Int,0.05*length(lengths))
        @test maximum(lengths)<=2.5
    finally
        A.finalize()
    end
end

@testset ".geo MathEvalAniso background field reaches Box meshing" begin
    source="""
    Mesh.MeshSizeExtendFromBoundary = 0;
    Mesh.MeshSizeFromPoints = 0;
    Mesh.MeshSizeFromCurvature = 0;
    Box(1) = {0,0,0, 1,1,1};
    Field[1] = MathEvalAniso;
    Field[1].M11 = "64";
    Field[1].M22 = "1";
    Field[1].M33 = "1";
    Field[1].M12 = "0";
    Field[1].M13 = "0";
    Field[1].M23 = "0";
    Background Field = 1;
    """
    mktempdir() do directory
        path=joinpath(directory,"metric_volume.geo")
        write(path,source)
        result=Tessella.execute_geo(path;mesh_dim=3)
        @test result.msg_error_count==0
        mesh=result.mesh
        @test mesh isa Mesh && size(mesh.tets,2)>0
        checked=metric_edge_audit(mesh,point->Diagonal([64.,1.,1.]))
        @test validate(mesh).ok && size(mesh.tets,2)>0
        @test all(>(0),checked.volumes) &&
              isapprox(sum(checked.volumes),1.;atol=1e-12)
        check_metric_bound(checked;fraction=0.05)
    end
end

end
