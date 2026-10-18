module ModelMetricAlgorithm7FunctionalChecks
using Test,Tessella,LinearAlgebra
const K=Tessella.Model
const N=Tessella.MeshTypes
const A=Tessella.API

# Unit square with straight boundary Lines. `counts` selects a transfinite
# boundary subdivision so retained Curve cells are deterministic; `nothing`
# leaves the ordinary scalar curve grading in place.
function rectangle(;counts=6)
    model=GeoModel()
    points=[add_point!(model,x,y,0.;mesh_size=.5) for (x,y) in
            ((0.,0.),(1.,0.),(1.,1.),(0.,1.))]
    curves=[add_line!(model,points[i],points[mod1(i+1,4)]) for i in 1:4]
    for curve in curves
        counts===nothing || K.set_transfinite_curve!(model,curve,counts)
    end
    loop=add_curve_loop!(model,curves)
    tag=add_plane_surface!(model,[loop])
    return model,tag
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
        all(point->sqrt(dot(displacement,matrix(point)*displacement))<=target+1e-12,
            (first,last,(first+last)/2))
    end
end
function edges(cells)
    Set(minmax(cell[a],cell[b]) for cell in cells for (a,b) in ((1,2),(2,3),(3,1)))
end

@testset "Mesh.Algorithm=7 routes the ordinary selector to the metric producer" begin
    # Scalar fallback: no tensor field, so the ordinary size context is the
    # isotropic metric I/h². Algorithm 7 must produce a valid triangulation
    # rather than reject or alias another algorithm.
    model,tag=rectangle()
    model.meshing.algorithm_2d=7
    mesh=mesh_model_surface(model,tag)
    @test validate(mesh).ok && N.ntris(mesh)>0 && N.nnodes(mesh)>=9
    @test isapprox(area(mesh),1.;atol=1e-12)
    @test all(cell->N.triangle_area(
        (Tuple(mesh.coords[:,row]) for row in cell)...)>0,eachcol(mesh.tris))
    incidence=edge_incidence(mesh)
    @test all(in((1,2)),values(incidence))
    boundary=Set(key for (key,count) in incidence if count==1)
    # Transfinite N=6 places six nodes per curve — five boundary edges each,
    # twenty total across the four sides.
    @test length(boundary)==4*(6-1)
    @test all(boundary) do (a,b)
        pa,pb=mesh.coords[:,a],mesh.coords[:,b]
        (pa[1]==pb[1] && pa[1] in (0.,1.)) || (pa[2]==pb[2] && pa[2] in (0.,1.))
    end
    @test all(point->0<=point[1]<=1 && 0<=point[2]<=1 && point[3]==0,eachcol(mesh.coords))

    # Per-surface selection: algorithm 7 through `set_algorithm!` (the API
    # `setAlgorithm`/`.geo MeshAlgorithm` store) with a tensor `size_field`
    # must respect the anisotropy. Boundary cells at 0.2 satisfy metric 1.0
    # under diag(16,1,1) (0.2*4=0.8), and every edge is metric-bounded by the
    # BAMG unit target.
    model,tag=rectangle()
    K.set_algorithm!(model,2,tag,7)
    xfield=MathEvalAnisoField(m11="16",m22="1",m33="1")
    yfield=MathEvalAnisoField(m11="1",m22="16",m33="1")
    xmesh=mesh_model_surface(model,tag;size_field=xfield)
    ymesh=mesh_model_surface(model,tag;size_field=yfield)
    @test validate(xmesh).ok && validate(ymesh).ok
    @test metric_bound(xmesh,point->Diagonal([16.,1.,1.]),1.0)
    @test metric_bound(ymesh,point->Diagonal([1.,16.,1.]),1.0)
    @test isapprox(area(xmesh),1.;atol=1e-12) && isapprox(area(ymesh),1.;atol=1e-12)
    @test !isequal(xmesh.coords,ymesh.coords) || xmesh.tris!=ymesh.tris
    # The other algorithms are unchanged on the same fixture.
    model.meshing.algorithm[(2,tag)]=6
    @test validate(mesh_model_surface(model,tag)).ok
end

@testset "algorithm 7 keeps unsupported surfaces loud" begin
    # A circle-bounded face is planar but its boundary Curve is not a Line —
    # the metric producer contract rejects it rather than substituting the
    # scalar mesher.
    model=GeoModel()
    center=add_point!(model,.5,.5,0.;mesh_size=.5)
    rim=[add_point!(model,x,y,0.;mesh_size=.5) for (x,y) in
         ((1.,.5),(.5,1.),(0.,.5),(.5,0.))]
    arcs=[add_circle_arc!(model,rim[i],center,rim[mod1(i+1,4)]) for i in 1:4]
    loop=add_curve_loop!(model,arcs)
    tag=add_plane_surface!(model,[loop])
    model.meshing.algorithm_2d=7
    @test_throws ArgumentError mesh_model_surface(model,tag)

    # A curved surface kind stays rejected: a ruled filling between two
    # skew edges is not a plane Surface.
    model=GeoModel()
    p1=add_point!(model,0.,0.,0.;mesh_size=.5);p2=add_point!(model,1.,0.,0.;mesh_size=.5)
    p3=add_point!(model,0.,1.,1.;mesh_size=.5);p4=add_point!(model,1.,1.,1.;mesh_size=.5)
    curves=[add_line!(model,p1,p2),add_line!(model,p2,p4),
            add_line!(model,p4,p3),add_line!(model,p3,p1)]
    tag=add_ruled_surface!(model,[add_curve_loop!(model,curves)])
    model.meshing.algorithm_2d=7
    @test_throws ArgumentError mesh_model_surface(model,tag)

    # Unknown algorithms still reject through the same admission.
    model,tag=rectangle()
    model.meshing.algorithm_2d=9
    @test_throws ArgumentError mesh_model_surface(model,tag)
end

@testset "API Mesh.Algorithm=7 generates through the tensor path" begin
    A.initialize()
    try
        A.model.add("algorithm 7 session")
        points=[A.model.add_point(x,y,0.;meshSize=.5) for (x,y) in
                ((0.,0.),(1.,0.),(1.,1.),(0.,1.))]
        curves=[A.model.add_line(points[i],points[mod1(i+1,4)]) for i in 1:4]
        loop=A.model.add_curve_loop(curves)
        surface=A.model.add_plane_surface([loop])
        for curve in curves;A.mesh.set_transfinite_curve(curve,6);end
        for (name,value) in (("General.Terminal",0),("Mesh.Renumber",0),
                ("Mesh.MeshSizeFromPoints",0),("Mesh.MeshSizeFromCurvature",0),
                ("Mesh.MeshSizeExtendFromBoundary",0),("Mesh.Smoothing",0))
            A.option(name,value)
        end
        field=A.mesh.field.add("MathEvalAniso")
        for (name,expression) in zip(("M11","M22","M33","M12","M13","M23"),
                                     ("16","1","1","0","0","0"))
            A.mesh.field.set_string(field,name,expression)
        end
        A.mesh.field.set_as_background_mesh(field)
        A.option("Mesh.Algorithm",7)
        A.mesh.generate(2)
        tags,xyz,_=A.mesh.get_nodes()
        coordinates=Dict(tag=>xyz[3i-2:3i] for (i,tag) in enumerate(tags))
        families,labels,nodes=A.mesh.get_elements(2,surface)
        @test families==Int32[2] && !isempty(labels)
        cells=collect(eachcol(reshape(nodes[1],3,:)))
        @test all(edges(cells)) do (a,b)
            first,last=coordinates[a],coordinates[b];displacement=last-first
            all(point->sqrt(dot(displacement,Diagonal([16.,1.,1.])*displacement))<=1.0+1e-12,
                (first,last,(first+last)/2))
        end
        signed=[(coordinates[cell[2]][1]-coordinates[cell[1]][1])*
                (coordinates[cell[3]][2]-coordinates[cell[1]][2])/2-
                (coordinates[cell[2]][2]-coordinates[cell[1]][2])*
                (coordinates[cell[3]][1]-coordinates[cell[1]][1])/2 for cell in cells]
        @test all(>(0),signed) && isapprox(sum(signed),1.;atol=1e-12)
    finally
        A.finalize()
    end
end

@testset ".geo Mesh.Algorithm=7 reaches the same producer" begin
    path=tempname()*".geo"
    write(path,"""
        Mesh.Algorithm = 7;
        Point(1) = {0,0,0,0.5}; Point(2) = {1,0,0,0.5};
        Point(3) = {1,1,0,0.5}; Point(4) = {0,1,0,0.5};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Transfinite Curve {1,2,3,4} = 6;
        Curve Loop(1) = {1,2,3,4};
        Plane Surface(1) = {1};
        Mesh 2;
        """)
    try
        execution=Tessella.execute_geo(path)
        @test execution.msg_error_count==0
        mesh=execution.mesh
        @test validate(mesh).ok && N.ntris(mesh)>0
        @test isapprox(area(mesh),1.;atol=1e-12)
    finally
        rm(path;force=true)
    end
end

end
