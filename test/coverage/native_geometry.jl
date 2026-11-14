using Test,LinearAlgebra,Tessella
import Tessella.MeshTypes as MT

function boundary_sides(corners,nu,nv)
    a,b,c,d=corners
    segment(p,q,n)=[ntuple(axis->(1-i/n)*p[axis]+(i/n)*q[axis],3) for i in 0:n]
    (segment(a,b,nu),segment(b,c,nv),segment(c,d,nu),segment(d,a,nv))
end

@testset "Structured native geometry and topology" begin
    for (corners,warped) in (
        (((0.,0.,0.),(3.,0.,0.),(3.,2.,0.),(0.,2.,0.)),false),
        (((2.,-1.,3.),(5.,-1.,3.),(6.,1.,4.),(3.,1.,4.)),false),
        (((0.,0.,0.),(1.,0.,0.),(1.,1.,0.5),(0.,1.,0.)),true)),
        (nu,nv) in ((1,1),(3,2),(4,4)),
        arrangement in (:left,:right,:alternate_left,:alternate_right)
        sides=boundary_sides(corners,nu,nv)
        saved=deepcopy(sides)
        mesh=mesh_transfinite_patch(sides...;allow_warped=warped,arrangement,
            face_tag=7,side_tags=(1,2,3,4))
        @test validate(mesh).ok
        @test (MT.nnodes(mesh),MT.ntris(mesh),MT.nsegs(mesh))==
            ((nu+1)*(nv+1),2nu*nv,2(nu+nv))
        @test sides==saved
        @test all(==(7),mesh.tri_tag)
        @test Set(mesh.seg_tag)==Set(Int32[1,2,3,4])
        edge_set,incidence=MT.boundary_edges(mesh.tris)
        @test incidence==2
        @test length(edge_set)==2(nu+nv)
        a,b,c,d=corners
        for j in 0:nv,i in 0:nu
            u=i/nu;v=j/nv
            expected=ntuple(axis->(1-u)*(1-v)*a[axis]+u*(1-v)*b[axis]+
                u*v*c[axis]+(1-u)*v*d[axis],3)
            @test any(column->isapprox(mesh.coords[:,column],collect(expected);atol=2e-14,rtol=2e-14),axes(mesh.coords,2))
        end
        if !warped
            expected_area=norm(cross(collect(b).-collect(a),collect(d).-collect(a)))
            area=sum(MT.triangle_area((MT.node(mesh,row) for row in cell)...)
                for cell in eachcol(mesh.tris))
            @test area≈expected_area rtol=2e-14
        end
    end
    sides=boundary_sides(((0.,0.,0.),(3.,0.,0.),(3.,2.,0.),(0.,2.,0.)),3,2)
    @test_throws ArgumentError mesh_transfinite_patch(sides...;max_nodes=11)
    @test_throws ArgumentError mesh_transfinite_patch(sides...;max_triangles=11)
    @test_throws ArgumentError mesh_transfinite_patch(sides...;arrangement=:invalid)
    @test_throws ArgumentError mesh_transfinite_patch(sides...;face_tag=true)
    invalid=copy(sides[1]);invalid[2]=(NaN,0.,0.)
    @test_throws ArgumentError mesh_transfinite_patch(invalid,sides[2],sides[3],sides[4])
    mismatched=copy(sides[2]);push!(mismatched,mismatched[end])
    @test_throws ArgumentError mesh_transfinite_patch(sides[1],mismatched,sides[3],sides[4])
    bowtie=([(0.,0.,0.),(1.,1.,0.)],[(1.,1.,0.),(0.,1.,0.)],
        [(0.,1.,0.),(1.,0.,0.)],[(1.,0.,0.),(0.,0.,0.)])
    @test_throws ArgumentError mesh_transfinite_patch(bowtie...)
    warped=boundary_sides(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.5),(0.,1.,0.)),2,2)
    @test_throws ArgumentError mesh_transfinite_patch(warped...)
end

@testset "Direct 3-D Delaunay fill and boundary conservation" begin
    for (origin,scale) in (((0.,0.,0.),1.),((2.,-3.,4.),2.)),interior in (false,true)
        points=[(0.,0.,0.),(2.,0.,0.),(0.,3.,0.),(0.,0.,4.)]
        interior && push!(points,(0.2,0.3,0.4))
        transformed=[ntuple(axis->origin[axis]+scale*point[axis],3) for point in points]
        x=first.(transformed);y=getindex.(transformed,2);z=last.(transformed)
        source=(copy(x),copy(y),copy(z))
        triangulation=Tessella.Mesh3D.delaunay3d(x,y,z;perturb=false)
        mesh=Tessella.Mesh3D.to_mesh3(triangulation)
        @test validate(mesh).ok
        @test Tessella.Mesh3D.check_consistency3(triangulation)[1]
        @test MT.ntets(mesh)==(interior ? 4 : 1)
        @test (x,y,z)==source
        @test Set(Tuple(column) for column in eachcol(mesh.coords))==Set(transformed)
        volumes=[MT.tet_signed_volume((MT.node(mesh,row) for row in cell)...)
            for cell in eachcol(mesh.tets)]
        @test all(>(0),volumes)
        @test sum(volumes)≈4scale^3 rtol=2e-14
        faces,incidence=MT.boundary_faces(mesh.tets)
        @test incidence==(interior ? 2 : 1) && length(faces)==4
    end
end
