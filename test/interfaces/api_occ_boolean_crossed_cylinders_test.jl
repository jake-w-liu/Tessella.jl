module OCCBooleanCrossedCylindersChecks
using Test,Tessella
const A=Tessella.API

# Equal radii (boolean.geo's cylinders, r = 0.98, length 5.6): the section
# is two ellipses crossing where the cylinders touch. The fuse once returned
# one cylinder unchanged and the fragment both untouched. Exact volumes: the
# common part 16r³/3 (Steinmetz), the three-cylinder common part 8(2−√2)r³;
# Gmsh 4.15.2 gives 28.7727878844, pieces 5.93827430/5.01969066 and
# 40.0403459044. Gmsh also has 8/14/8/1 and 8/14/12/5 entities (it keeps
# the cylinders' seam vertices); Tessella has 6/12/8/1 and 6/12/12/5.
# Interpreted these take over half an hour (the fitted ellipse pcurves and
# the sphere): compiled runs only.
@testset "OCC Booleans of crossed cylinders of equal radius" begin
  if Base.JLOptions().compile_enabled==1
    r=1.4*0.7;L=5.6
    cyl=π*r^2*L;common=16r^3/3
    function run(f,n)
        A.initialize()
        try
            A.option("General.Terminal",0)
            occ=A.model.occ
            occ.add_cylinder(-L/2,0,0,L,0,0,r;tag=4)
            occ.add_cylinder(0,-L/2,0,0,L,0,r;tag=5)
            n==3 && occ.add_cylinder(0,0,-L/2,0,0,L,r;tag=6)
            occ.synchronize()
            f(occ)
        finally
            A.finalize()
        end
    end
    run(2) do occ
        out,_=occ.fuse([(3,4)],[(3,5)];tag=7)
        @test out==[(3,7)]
        @test occ.get_mass(3,7)≈2cyl-common rtol=1e-8
        @test_broken Tuple(length(occ.get_entities(d)) for d in 0:3)==(8,14,8,1)
    end
    run(2) do occ
        out,images=occ.fragment([(3,4)],[(3,5)])
        @test length(out)==5
        ms=sort([occ.get_mass(3,t) for (_,t) in out])
        @test ms≈sort([fill((cyl-common)/2,4);common]) rtol=1e-8
        @test images==[[(3,1),(3,2),(3,3)],[(3,4),(3,2),(3,5)]] broken=true
    end
    run(3) do occ
        out,_=occ.fuse([(3,4)],[(3,5),(3,6)];tag=7)
        @test out==[(3,7)]
        @test occ.get_mass(3,7)≈3cyl-3common+8*(2-sqrt(2))*r^3 rtol=1e-8
        # boolean.geo's last step: (box ∩ sphere) less the cross. The box
        # face around each cylinder is an annulus whose interior probe sat
        # in the hole ("cannot orient boolean result face"). Gmsh 4.15.2:
        # 3.98182885452
        occ.add_box(-1.4,-1.4,-1.4,2.8,2.8,2.8;tag=1)
        occ.add_sphere(0,0,0,1.75;tag=2)
        occ.intersect([(3,1)],[(3,2)];tag=3)
        out,_=occ.cut([(3,3)],[(3,7)];tag=8)
        @test out==[(3,8)]
        @test occ.get_mass(3,8)≈3.98182885452 rtol=1e-6
    end
  end
end

# Two cylinders crossing at right angles, radii 0.5 and 0.3: the fuse once
# measured 0.63. The common part of two perpendicular cylinders through one
# axis point is 8 ∫₀^r √(r²−x²) √(R²−x²) dx = 0.2693703275 for these radii,
# so fuse = π/2 + 0.18π − 0.2693703275 and cut = π/2 − 0.2693703275. Gmsh
# 4.15.2 gives 1.8669200576 and 1.3014335187, off in the sixth digit. The
# section is a traced curve, minutes when interpreted: compiled runs only.
@testset "OCC Booleans of two crossed cylinders" begin
    if Base.JLOptions().compile_enabled==1
        common=0.2693703275
        for (operation,exact) in ((:fuse,π/2+0.18π-common),
                (:cut,π/2-common),(:intersect,common))
            A.initialize()
            try
                A.option("General.Terminal",0)
                occ=A.model.occ
                upright=occ.add_cylinder(0,0,0,0,0,2,0.5)
                across=occ.add_cylinder(-1,0,1,2,0,0,0.3)
                out,_=getproperty(occ,operation)([(3,upright)],[(3,across)])
                @test length(out)==1
                @test occ.get_mass(3,out[1][2])≈exact atol=2e-6
            finally
                A.finalize()
            end
        end
    end
end
end
