module TopologyNativeGenerationCertificate
using Test, Tessella
const A=Tessella.API
function _add_api_topology_tetrahedron()
    for (tag,x,y,z) in ((10,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (7,0.0,1.0,0.0),(5,0.0,0.0,1.0))
        A.model.add_point(x,y,z;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((8,10,2),(3,2,7),(11,7,10),(6,10,5),(12,2,5),(4,7,5))
        A.model.add_line(first_point,last_point;tag=tag)
    end
    for (tag,curves) in ((21,[8,3,11]),(22,[8,12,-6]),
                         (23,[3,4,-12]),(24,[11,6,-4]))
        A.model.add_curve_loop(curves;tag=tag)
        A.model.add_plane_surface([tag];tag=tag)
    end
    A.model.add_surface_loop([21,-22,23,-24];tag=30)
    A.model.add_volume([30];tag=40)
    return nothing
end
@testset "Gmsh4152 native h1 tetrahedron records and exact volume" begin
    A.finalize()
    A.initialize()
    try
        _add_api_topology_tetrahedron()
        generated=A.mesh.generate(3)
        @test Tessella.validate(generated).ok
        nodes,flat,_=A.mesh.get_nodes();xyz=Dict(tag=>Rational{BigInt}.(flat[3i-2:3i]) for (i,tag) in enumerate(nodes))
        labels,cells=A.mesh.get_elements_by_type(4,40)
        volume=zero(Rational{BigInt});faces=Dict{NTuple{3,UInt64},Int}()
        patterns=((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        for cell in eachcol(reshape(cells,4,:))
            p,q,r,s=[xyz[tag] for tag in cell];u=q-p;v=r-p;w=s-p
            determinant=u[1]*(v[2]*w[3]-v[3]*w[2])-u[2]*(v[1]*w[3]-v[3]*w[1])+u[3]*(v[1]*w[2]-v[2]*w[1])
            @test determinant>0;volume+=determinant/6
            for pattern in patterns
                key=Tuple(sort!(collect(cell[collect(pattern)])));faces[key]=get(faces,key,0)+1
            end
        end
        @test volume==1//6
        _,triangles=A.mesh.get_elements_by_type(2)
        @test Set(face for (face,count) in faces if count==1)==Set(Tuple(sort!(collect(cell))) for cell in eachcol(reshape(triangles,3,:)))
        @test all(value->value in (1,2),values(faces))

        crc=Tessella.mesh_crc(generated)
        @test crc.n_nodes==7
        @test crc.n_tets==4
        @test crc.n_boundary_faces==10
        @test Tessella.mesh_crc(A.mesh.get())==crc
        @test crc.sha=="357e66dce7c727fc08cba1f1d0e9e6970a18ffd0e9f78aee3326b227c250a094"
        # The returned facade uses a separate storage order. Certify its
        # represented coordinates independently of the public entity catalog.
        carrier_volume=zero(Rational{BigInt})
        for tet in eachcol(generated.tets)
            p,q,r,s=[Rational{BigInt}.(generated.coords[:,row]) for row in tet]
            u=q-p;v=r-p;w=s-p
            determinant=u[1]*(v[2]*w[3]-v[3]*w[2])-u[2]*(v[1]*w[3]-v[3]*w[1])+u[3]*(v[1]*w[2]-v[2]*w[1])
            @test determinant>0
            carrier_volume+=determinant/6
        end
        @test carrier_volume==1//6
    finally
        A.finalize()
    end
end
end
