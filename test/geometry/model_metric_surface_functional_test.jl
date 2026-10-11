module ModelMetricSurfaceFunctionalChecks
using Test,Tessella,LinearAlgebra
const K=Tessella.Model
const N=Tessella.MeshTypes

function rectangle(;oblique=false,counts=3,hole=false)
    model=GeoModel()
    origin=oblique ? [2.,-1.,3.] : zeros(3)
    basis=oblique ? [1. 0.;0. .6;0. .8] : [1. 0.;0. 1.;0. 0.]
    loops=Int[]
    polygons=hole ? (((0.,0.),(1.,0.),(1.,1.),(0.,1.)),
                     ((.25,.25),(.25,.75),(.75,.75),(.75,.25))) :
                    (((0.,0.),(1.,0.),(1.,1.),(0.,1.)),)
    for polygon in polygons
        points=[add_point!(model,(origin+basis*collect(point))...;mesh_size=.5) for point in polygon]
        curves=[add_line!(model,points[index],points[mod1(index+1,4)]) for index in 1:4]
        for curve in curves;K.set_transfinite_curve!(model,curve,counts);end
        push!(loops,add_curve_loop!(model,curves))
    end
    tag=add_plane_surface!(model,loops)
    return model,tag,origin,basis
end

function edge_incidence(mesh)
    incidence=Dict{Tuple{Int32,Int32},Int}()
    for cell in eachcol(mesh.tris),edge in ((cell[1],cell[2]),(cell[2],cell[3]),(cell[3],cell[1]))
        key=minmax(edge...);incidence[key]=get(incidence,key,0)+1
    end
    incidence
end
function area(mesh)
    sum(norm(cross(mesh.coords[:,cell[2]]-mesh.coords[:,cell[1]],
                   mesh.coords[:,cell[3]]-mesh.coords[:,cell[1]]))/2 for cell in eachcol(mesh.tris))
end
function metric_bound(mesh,matrix,target)
    all(keys(edge_incidence(mesh))) do (a,b)
        first=mesh.coords[:,a];last=mesh.coords[:,b];displacement=last-first
        all(point->sqrt(dot(displacement,matrix(point)*displacement))<=target+1e-13,
            (first,last,(first+last)/2))
    end
end
function initial_product(model,tag)
    working=K._model_mesh_transaction_copy(model)
    options=K._model_default_mesh1d_options(working,"metric functional reference")
    K._gmsh_native_initial_surface(working,tag,options,nothing,nothing,
        "metric functional reference";algorithm=3)
end
function unchanged_model(model,before)
    model.points==before.points && model.curves==before.curves && model.loops==before.loops &&
        model.surfaces==before.surfaces && model.curve_params==before.curve_params &&
        model.meshing.transfinite_curves==before.meshing.transfinite_curves
end

@testset "ordinary Surface consumes a spatially varying tensor" begin
    model,tag,origin,basis=rectangle();before=deepcopy(model);initial=initial_product(model,tag)
    field=MathEvalAnisoField(m11="4+12*x",m22="1",m33="1")
    matrix(point)=Diagonal([4+12point[1],1.,1.])
    mesh=mesh_model_surface(model,tag,field;max_metric_edge=2.01,max_nodes=100,max_triangles=200)
    @test validate(mesh).ok && N.nnodes(mesh)>N.nnodes(initial.mesh)
    @test metric_bound(mesh,matrix,2.01)
    @test isapprox(area(mesh),1.;atol=1e-12) && all(cell->N.triangle_area(
        (Tuple(mesh.coords[:,row]) for row in cell)...)>0,eachcol(mesh.tris))
    @test all(point->0<=point[1]<=1 && 0<=point[2]<=1 && point[3]==0,eachcol(mesh.coords))
    @test unchanged_model(model,before)
    # Exercise the same public front door's detached constructor metadata.
    product=mesh_model_surface(model,tag,field;max_metric_edge=2.01,_return_surface_production=true)
    @test isequal(product.mesh.coords,mesh.coords) && product.mesh.tris==mesh.tris
    retained=findall(product.retained_nodes)
    @test product.node_tags[retained]==initial.node_tags &&
        isequal(product.mesh.coords[:,retained],initial.mesh.coords)
    @test product.node_owners[retained]==initial.node_owners &&
        product.node_parameters[retained]==initial.node_parameters
    @test all(owner->owner==(2,Int32(tag)),product.node_owners[.!product.retained_nodes])
    @test allunique(product.node_tags) && allunique(product.element_tags) &&
        length(product.element_tags)==N.ntris(product.mesh)
    @test product.node_high_water==maximum(product.node_tags) &&
        product.element_high_water>=maximum(product.element_tags)
    original_boundary=Set(minmax(initial.node_tags[a],initial.node_tags[b]) for
        ((a,b),count) in edge_incidence(initial.mesh) if count==1)
    actual_boundary=Set(minmax(product.node_tags[a],product.node_tags[b]) for
        ((a,b),count) in edge_incidence(product.mesh) if count==1)
    @test actual_boundary==original_boundary
    product.mesh.coords[1,1]=99.
    @test unchanged_model(model,before) && mesh.coords[1,1]!=99.
end

@testset "tensor directions and tangential cross terms affect actual cells" begin
    model,tag,_,_=rectangle()
    xfield=MathEvalAnisoField(m11="16",m22="1",m33="1")
    yfield=MathEvalAnisoField(m11="1",m22="16",m33="1")
    xmesh=mesh_model_surface(model,tag,xfield;max_metric_edge=2.01)
    ymesh=mesh_model_surface(model,tag,yfield;max_metric_edge=2.01)
    @test metric_bound(xmesh,point->Diagonal([16.,1.,1.]),2.01)
    @test metric_bound(ymesh,point->Diagonal([1.,16.,1.]),2.01)
    @test !isequal(xmesh.coords,ymesh.coords) || xmesh.tris!=ymesh.tris
    @test isapprox(area(xmesh),1.;atol=1e-12) && isapprox(area(ymesh),1.;atol=1e-12)
    model,tag,origin,basis=rectangle(;oblique=true)
    tensor=[4. .5 .25;.5 2. .4;.25 .4 3.]
    field=MathEvalAnisoField(m11="4",m22="2",m33="3",m12=".5",m13=".25",m23=".4")
    mesh=mesh_model_surface(model,tag,field;max_metric_edge=1.5)
    @test validate(mesh).ok && metric_bound(mesh,point->tensor,1.5)
    @test isapprox(area(mesh),1.;atol=1e-12)
    @test all(point->isapprox(point[3]-origin[3],(4/3)*(point[2]-origin[2]);atol=1e-12),eachcol(mesh.coords))
    @test all(keys(edge_incidence(mesh))) do (a,b)
        displacement=mesh.coords[:,b]-mesh.coords[:,a];local_displacement=transpose(basis)*displacement
        isapprox(dot(displacement,tensor*displacement),
            dot(local_displacement,(transpose(basis)*tensor*basis)*local_displacement);atol=1e-12)
    end
end

@testset "literal hole boundaries remain owned" begin
    model,tag,_,_=rectangle(;hole=true)
    field=MathEvalAnisoField(m11="4+12*x",m22="1",m33="1")
    mesh=mesh_model_surface(model,tag,field;max_metric_edge=2.01)
    @test validate(mesh).ok && isapprox(area(mesh),.75;atol=1e-12)
    @test metric_bound(mesh,point->Diagonal([4+12point[1],1.,1.]),2.01)
    @test all(cell->begin
        center=sum(mesh.coords[:,cell];dims=2)/3
        !(.25<center[1]<.75 && .25<center[2]<.75)
    end,eachcol(mesh.tris))
end

@testset "metric and ownership failures leave the model unchanged" begin
    model,tag,_,_=rectangle();before=deepcopy(model)
    field=MathEvalAnisoField(m11="16",m22="1",m33="1")
    @test_throws ArgumentError mesh_model_surface(model,tag,field;max_metric_edge=0)
    @test_throws ArgumentError mesh_model_surface(model,tag,field;max_metric_edge=true)
    @test_throws ArgumentError mesh_model_surface(model,tag,field;max_metric_edge=2.01,max_nodes=8)
    @test_throws ArgumentError mesh_model_surface(model,tag,field;max_metric_edge=2.01,max_triangles=6)
    @test_throws ArgumentError mesh_model_surface(model,tag,field;max_metric_edge=2.01,_mesh_element_high_water=typemax(UInt64))
    bad=MathEvalAnisoField(m11="1",m22="-1",m33="1")
    @test_throws ArgumentError mesh_model_surface(model,tag,bad;max_metric_edge=2.01)
    coarse,coarse_tag,_,_=rectangle(;counts=2)
    @test_throws ArgumentError mesh_model_surface(coarse,coarse_tag,field;max_metric_edge=2.01)
    @test isempty(coarse.curve_params)
    @test unchanged_model(model,before)
end
end
