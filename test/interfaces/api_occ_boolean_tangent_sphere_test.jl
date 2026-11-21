module OCCBooleanTangentSphereChecks
using Test,Tessella
const A=Tessella.API

# A sphere touching a box face in one point was refused ("tangent
# plane×sphere contact is unsupported"). There is no section curve: the
# operands stay whole, and what the operation returns is Gmsh 4.15.2's. The
# sphere's volume is 4π·0.3³/3.
function touching(run,z)
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        ball=occ.add_sphere(0,0,z,0.3)
        return run(occ,box,ball)
    finally
        A.finalize()
    end
end
entity_counts()=Tuple(length(A.model.get_entities(d)) for d in 0:3)
faces_of(volume)=sort!([abs(t) for (_,t) in A.model.get_boundary([(3,volume)],false,true)])
const BALL=4π*0.3^3/3

@testset "OCC Booleans of a sphere resting on a box face" begin
    # The sphere's pole is the contact point. OpenCASCADE makes it an
    # internal vertex of the box face, so fuse and fragments return the box
    # as a modified shape (volume 3 and faces 8 to 13 from the fuse, volume
    # 1 with faces 8 to 13 from the fragments) and the sphere as it was.
    touching(1.3) do occ,box,ball
        out,map=occ.fuse([(3,box)],[(3,ball)])
        occ.synchronize()
        @test out==[(3,3),(3,2)]
        @test map==[Tuple{Int,Int}[],[(3,2)]]
        @test entity_counts()==(10,15,7,2)
        @test faces_of(3)==collect(8:13)
        @test faces_of(2)==[7]
        @test occ.get_mass(3,3)≈8 rtol=1e-12
        @test occ.get_mass(3,2)≈BALL rtol=1e-12
    end
    touching(1.3) do occ,box,ball
        out,map=occ.fragment([(3,box)],[(3,ball)])
        occ.synchronize()
        @test out==[(3,1),(3,2)]
        @test map==[[(3,1)],[(3,2)]]
        @test faces_of(1)==collect(8:13)
        @test faces_of(2)==[7]
    end
    touching(1.3) do occ,box,ball
        out,map=occ.cut([(3,box)],[(3,ball)])
        occ.synchronize()
        @test out==[(3,1)]
        @test map==[[(3,1)],Tuple{Int,Int}[]]
        @test faces_of(1)==collect(1:6)
        @test occ.get_mass(3,1)≈8 rtol=1e-12
    end
    touching(1.3) do occ,box,ball
        out,_=occ.intersect([(3,box)],[(3,ball)])
        occ.synchronize()
        @test isempty(out)
        @test entity_counts()==(0,0,0,0)
    end
end

@testset "OCC Booleans of a sphere touching a box face from inside" begin
    # The cut is the box with a cavity. The cavity's shell was taken for a
    # solid of its own, with volume -0.113: its containment in the box was
    # asked at the sphere's pole, which is the contact point.
    touching(0.7) do occ,box,ball
        out,_=occ.cut([(3,box)],[(3,ball)])
        occ.synchronize()
        @test out==[(3,1)]
        @test entity_counts()==(10,15,7,1)
        @test length(A.model.get_boundary([(3,1)],false,true))==7
        @test occ.get_mass(3,1)≈8-BALL rtol=1e-12
        # the contact point, a vertex of the cavity, is an internal vertex
        # of the box face it touches: both surface meshes get the node
        touched=[t for (_,t) in A.model.get_entities(2)
                 if !isempty(A.mesh.get_embedded(2,t))]
        @test length(touched)==1
        point=only(A.mesh.get_embedded(2,only(touched)))
        @test point[1]==0
        @test collect(A.model.get_value(0,point[2],Float64[]))≈[0,0,1] atol=1e-12
    end
    touching(0.7) do occ,box,ball
        out,map=occ.fragment([(3,box)],[(3,ball)])
        occ.synchronize()
        @test out==[(3,3),(3,2)]
        @test map==[[(3,3),(3,2)],[(3,2)]]
        @test occ.get_mass(3,3)≈8-BALL rtol=1e-12
        @test occ.get_mass(3,2)≈BALL rtol=1e-12
    end
    touching(0.7) do occ,box,ball
        out,_=occ.intersect([(3,box)],[(3,ball)])
        occ.synchronize()
        @test out==[(3,2)]
        @test entity_counts()==(2,3,1,1)
    end
    touching(0.7) do occ,box,ball
        out,_=occ.fuse([(3,box)],[(3,ball)])
        occ.synchronize()
        @test out==[(3,1)]
        @test occ.get_mass(3,1)≈8 rtol=1e-12
    end
end
end
