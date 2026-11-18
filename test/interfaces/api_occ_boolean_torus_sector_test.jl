module OCCBooleanTorusSectorChecks
using Test,Tessella
const A=Tessella.API

# A plane through a torus' axis leaves a genuine half-turn sector whose
# carrier record still says `angle=2π`. The face area, the face box and the
# solid box used the carrier's sweep instead of the wire's window, so the
# half torus reported the whole torus' area and a box reaching across the
# cutting plane.
@testset "OCC torus sector cut out by a Boolean" begin
    R=0.2925;r=0.1575
    for (name,volume) in (("intersect",π^2*R*r^2),("cut",8-π^2*R*r^2))
        @testset "$name" begin
            A.initialize()
            try
                A.option("General.Terminal",0)
                occ=A.model.occ
                box=occ.add_box(-1,-1,-1,2,2,2)
                torus=occ.add_torus(1,0.05,0.1,R,r)
                out,_=getfield(occ,Symbol(name))([(3,box)],[(3,torus)])
                occ.synchronize()
                tag=out[1][2]
                @test occ.get_mass(3,tag)≈volume rtol=1e-8
                torus_faces=[t for (_,t) in A.model.get_entities(2)
                             if A.model.get_type(2,t)=="Torus"]
                @test length(torus_faces)==1
                @test occ.get_mass(2,only(torus_faces))≈2π^2*R*r rtol=1e-10
                half=[1-R-r,0.05-R-r,0.1-r,1.0,0.05+R+r,0.1+r]
                @test collect(occ.get_bounding_box(2,only(torus_faces)))≈half atol=1e-9
                expected=name=="intersect" ? half : [-1.,-1.,-1.,1.,1.,1.]
                @test collect(A.model.get_bounding_box(3,tag))≈expected atol=1e-9
                @test collect(occ.get_bounding_box(3,tag))≈expected atol=1e-9
            finally
                A.finalize()
            end
        end
    end
end
end
