module OCCBooleanClippedCylinderChecks
using Test,Tessella
const A=Tessella.API

# A cylinder wall clipped between two planes keeps the carrier's seam
# rectangle but starts at v≠0. The wire-orientation sign of its divergence
# integral was read from an un-unwrapped plane shoelace, which is negative
# for such a rectangle — the wall's contribution entered the solid's volume
# with the wrong sign (box minus a through-hole reported MORE than the box).
@testset "OCC volume of solids bounded by a clipped cylinder wall" begin
    r=0.245
    hole=π*r^2*2
    for (name,volume,center) in (
            ("cut",8-hole,-hole.*[0.03,0.02,0.0]./(8-hole)),
            ("intersect",hole,[0.03,0.02,0.0]))
        @testset "$name" begin
            A.initialize()
            try
                A.option("General.Terminal",0)
                occ=A.model.occ
                box=occ.add_box(-1,-1,-1,2,2,2)
                cylinder=occ.add_cylinder(0.03,0.02,-1.6,0,0,3.2,r)
                out,_=getfield(occ,Symbol(name))([(3,box)],[(3,cylinder)])
                @test length(out)==1
                @test occ.get_mass(3,out[1][2])≈volume rtol=1e-10
                @test collect(occ.get_center_of_mass(3,out[1][2]))≈center atol=1e-10
                inertia=occ.get_matrix_of_inertia(3,out[1][2])
                @test all(>(0),(inertia[1],inertia[5],inertia[9]))
            finally
                A.finalize()
            end
        end
    end
    A.initialize()
    try
        A.option("General.Terminal",0)
        occ=A.model.occ
        box=occ.add_box(-1,-1,-1,2,2,2)
        cylinder=occ.add_cylinder(0.03,0.02,-1.6,0,0,3.2,r)
        out,_=occ.fragment([(3,box)],[(3,cylinder)])
        masses=sort!(Float64[occ.get_mass(3,t) for (_,t) in out])
        @test masses≈[π*r^2*0.6,π*r^2*0.6,hole,8-hole] rtol=1e-10
    finally
        A.finalize()
    end
end
end
