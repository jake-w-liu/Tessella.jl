module ModelMetricAlgorithm7FullQuadChecks
using Test,Tessella,LinearAlgebra
const K=Tessella.Model
const N=Tessella.MeshTypes

# Unit square with straight boundary Lines; transfinite `counts` pins an even
# per-curve cell count so upstream's "1D mesh cannot be divided by 2" contract
# does not fire.
function rectangle(;counts=7)
    model=GeoModel()
    points=[add_point!(model,x,y,0.;mesh_size=.5) for (x,y) in
            ((0.,0.),(1.,0.),(1.,1.),(0.,1.))]
    curves=[add_line!(model,points[i],points[mod1(i+1,4)]) for i in 1:4]
    counts===nothing || foreach(c->K.set_transfinite_curve!(model,c,counts),curves)
    loop=add_curve_loop!(model,curves)
    tag=add_plane_surface!(model,[loop])
    return model,tag
end

# Skew ruled patch — genuinely curved, no planar frame exists.
function ruled(;counts=5)
    model=GeoModel()
    p1=add_point!(model,0.,0.,0.;mesh_size=.4);p2=add_point!(model,1.,0.,0.;mesh_size=.4)
    p3=add_point!(model,0.,1.,1.;mesh_size=.4);p4=add_point!(model,1.,1.,1.;mesh_size=.4)
    curves=[add_line!(model,p1,p2),add_line!(model,p2,p4),
            add_line!(model,p4,p3),add_line!(model,p3,p1)]
    counts===nothing || foreach(c->K.set_transfinite_curve!(model,c,counts),curves)
    tag=add_ruled_surface!(model,[add_curve_loop!(model,curves)])
    return model,tag
end

cell_counts(mesh)=(
    tris=sum(size(b.nodes,2) for b in mesh.blocks
        if b isa Tessella.ElementBlock && b.msh==2;init=0),
    quads=sum(size(b.nodes,2) for b in mesh.blocks
        if b isa Tessella.ElementBlock && b.msh==3;init=0))

function edge_incidence(mesh)
    incidence=Dict{Tuple{Int32,Int32},Int}()
    for b in mesh.blocks
        b isa Tessella.ElementBlock && b.msh in (2,3) || continue
        width=size(b.nodes,1)
        for cell in eachcol(b.nodes),index in 1:width
            key=minmax(cell[index],cell[mod1(index+1,width)])
            incidence[key]=get(incidence,key,0)+1
        end
    end
    return incidence
end

function quad_area(mesh)
    total=0.0
    for b in mesh.blocks
        b isa Tessella.ElementBlock && b.msh==3 || continue
        for cell in eachcol(b.nodes)
            pa,pb,pc,pd=(mesh.coords[:,cell[i]] for i in 1:4)
            total+=norm(cross(pb-pa,pc-pa))/2+norm(cross(pd-pc,pa-pc))/2
        end
    end
    return total
end

@testset "planar algorithm-7 full-quad recombination" begin
    # Oracle (gmsh 4.15.2, same fixture): Mesh.Algorithm=7 +
    # RecombinationAlgorithm=2/3 with a tensor field produces an
    # all-quadrilateral Surface — Simple emits 60 quads/73 nodes, Blossom 56
    # quads/69 nodes; both keep the graded 1-D boundary restored. The certified
    # producer's topology family matches: quads only, valid, closed boundary.
    for algo in (2,3)
        model,tag=rectangle()
        model.meshing.algorithm_2d=7
        model.meshing.recombine[(2,tag)]=45.0
        model.meshing.recombine_algo=algo
        mesh=mesh_model_surface(model,tag;
            size_field=MathEvalAnisoField(m11="16",m22="1",m33="1"))
        @test mesh isa Tessella.MixedMesh
        @test validate(mesh).ok
        tris,quads=cell_counts(mesh)
        @test tris==0 && quads>0
        # Transfinite N=7 keeps six boundary cells per curve — the full-quad
        # stage halves them then restores the stored Curve middles, so the
        # published boundary still carries all 24 retained segments.
        incidence=edge_incidence(mesh)
        boundary=Set(key for (key,count) in incidence if count==1)
        @test length(boundary)==4*(7-1)
        @test all(boundary) do (a,b)
            pa,pb=mesh.coords[:,a],mesh.coords[:,b]
            (pa[1]==pb[1] && pa[1] in (0.,1.)) || (pa[2]==pb[2] && pa[2] in (0.,1.))
        end
        @test all(point->0<=point[1]<=1 && 0<=point[2]<=1 && point[3]==0,
                  eachcol(mesh.coords))
        @test isapprox(quad_area(mesh),1.;atol=1e-12)
    end

    # Determinism: the certified producer and pairing are bitwise stable.
    model,tag=rectangle()
    model.meshing.algorithm_2d=7
    model.meshing.recombine[(2,tag)]=45.0
    model.meshing.recombine_algo=3
    mesh=mesh_model_surface(model,tag;
        size_field=MathEvalAnisoField(m11="16",m22="1",m33="1"))
    model2,tag2=rectangle()
    model2.meshing.algorithm_2d=7
    model2.meshing.recombine[(2,tag2)]=45.0
    model2.meshing.recombine_algo=3
    mesh2=mesh_model_surface(model2,tag2;
        size_field=MathEvalAnisoField(m11="16",m22="1",m33="1"))
    @test isequal(mesh.coords,mesh2.coords)
    @test [b.nodes for b in mesh.blocks if b isa Tessella.ElementBlock] ==
        [b.nodes for b in mesh2.blocks if b isa Tessella.ElementBlock]
end

@testset "ruled algorithm-7 full-quad recombination stays on the patch" begin
    # Oracle (gmsh 4.15.2 on the same warped RuledSurface + tensor field):
    # Blossom emits 44 quads, 0 triangles — the all-quad topology family.
    # The certified ruled chart subdivides through its own (u,v) domain and
    # lifts every born vertex back onto the blend.
    lc_hint=nothing
    for algo in (2,3)
        model,tag=ruled()
        model.meshing.algorithm_2d=7
        model.meshing.recombine[(2,tag)]=45.0
        model.meshing.recombine_algo=algo
        mesh=mesh_model_surface(model,tag;
            size_field=MathEvalAnisoField(m11="16",m22="1",m33="1"))
        @test mesh isa Tessella.MixedMesh
        @test validate(mesh).ok
        tris,quads=cell_counts(mesh)
        @test tris==0 && quads>0
        lc=K._model_lc(model);silent=Returns(nothing)
        @test all(eachcol(mesh.coords)) do point
            u,v=K._ruled_xyz_to_uv(model,tag,Tuple(point),1.0,lc,
                true,silent,silent,"test",false)
            isfinite(u) && isfinite(v) && -1e-7<=u<=1+1e-7 &&
                -1e-7<=v<=1+1e-7
        end
    end
end

@testset "explicit tensor overload and OCC plane full-quad" begin
    # The public `mesh_model_surface(model,tag,field)` overload routes the
    # same certified front — recombined planar tensor surfaces reach the
    # shared full-quad lifecycle.
    model,tag=rectangle()
    model.meshing.recombine[(2,tag)]=45.0
    model.meshing.recombine_algo=2
    mesh=mesh_model_surface(model,tag,
        MathEvalAnisoField(m11="16",m22="1",m33="1"))
    @test mesh isa Tessella.MixedMesh && validate(mesh).ok
    tris,quads=cell_counts(mesh)
    @test tris==0 && quads>0

    # OCC box faces are Plane records: a factory-encoded Surface meshes the
    # same planar metric front, middles mapped through
    # `_occ_edge_reparam_on_face`.
    model=GeoModel()
    add_box!(model,0.,0.,0.,1.,1.,1.;tag=1)
    for c in keys(model.curves);K.set_transfinite_curve!(model,c,5);end
    face=sort!(collect(keys(model.surfaces)))[1]
    model.meshing.algorithm_2d=7
    model.meshing.recombine[(2,face)]=45.0
    model.meshing.recombine_algo=2
    mesh=mesh_model_surface(model,face;
        size_field=MathEvalAnisoField(m11="16",m22="1",m33="1"))
    @test mesh isa Tessella.MixedMesh && validate(mesh).ok
    tris,quads=cell_counts(mesh)
    @test tris==0 && quads>0
end

@testset "unsupported recombination contexts still reject loudly" begin
    # Bipartite recombination (algorithm 4) stays a loud not-implemented
    # rejection on the metric producer.
    model,tag=rectangle()
    model.meshing.algorithm_2d=7
    model.meshing.recombine[(2,tag)]=45.0
    model.meshing.recombine_algo=4
    @test_throws ArgumentError mesh_model_surface(model,tag)
    # Transfinite Surface construction keeps rejecting in the certified
    # producer's admission — the explicit tensor overload consults
    # `_metric_surface_admit` before its refine stage runs.
    model,tag=rectangle()
    K.set_transfinite_surface!(model,tag)
    model.meshing.recombine[(2,tag)]=45.0
    model.meshing.recombine_algo=3
    @test_throws ArgumentError mesh_model_surface(model,tag,
        MathEvalAnisoField(m11="16",m22="1",m33="1"))
end

end
