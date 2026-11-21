module GeoExactExtrusionParametrizationTest

using Test,Tessella

# A surface swept by a rigid `Extrude` evaluates its source curve moved by the
# sweep (`InterpolateExtrudedSurface`, `Geometry.ExactExtrusion=1`), not the
# blend of its four generatrices: the twisted laterals' spline generatrices
# only approximate the sweep, so the blend missed the swept mesh nodes by
# about 1e-3 and their inversion ran every relaxation level before returning
# a wrong (u,v). Values below are Gmsh 4.15.2's.
const _GEO="""
Point(1)={0,0,0,0.1};Point(2)={0.3,0,0,0.1};Point(3)={0.3,0.3,0,0.1};Point(4)={0,0.3,0,0.1};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
Transfinite Curve{1:4}=3;Transfinite Surface{1};Recombine Surface{1};
out[]=Extrude{{0,0,0.2},{0,0,1},{0.15,0.15,0},Pi/2}{Surface{1};Layers{3};Recombine;};
c[]=Extrude{{0,1,0},{0.6,0,0},Pi/2}{Curve{2};Layers{2};};
"""

@testset "exact extrusion evaluates swept surfaces" begin
    A=Tessella.API
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            path=joinpath(directory,"twist.geo");write(path,_GEO)
            A.open(path)
            # Twisted lateral (translate + rotate) and a surface of revolution.
            @test A.model.get_value(2,13,[0.25,0.5])==
                [0.20303300858899104,-0.00909902576697319,0.1]
            @test A.model.get_value(2,30,[0.25,0.75])==
                [0.48519497029047304,0.075,0.27716385975338603]
            A.mesh.generate(2)
            tags,coords,params=A.mesh.get_nodes(2,13,false,true)
            @test tags==[49,50]
            @test coords==[0.22499999999967005,0.020096189432143696,0.06666666666666667,
                           0.27990381056747526,0.07499999999967005,0.13333333333333333]
            @test params==[0.49999999999873,0.3333333333333333,
                           0.49999999999873007,0.6666666666666666]
            @test A.mesh.get_nodes(2,21,false,true)[3]==
                [0.4999999999987167,0.33333333333333337,
                 0.4999999999987168,0.6666666666666666]
            @test A.mesh.get_nodes(2,30,false,true)[3]==[0.5,0.5]
            @test A.model.get_parametrization(2,13,coords[1:3])==params[1:2]
            if Base.JLOptions().compile_enabled==1
                A.model.get_parametrization(2,13,coords[1:3])
                # One inversion: a few surface evaluations, no per-step arrays.
                @test (@allocated A.model.get_parametrization(2,13,coords[1:3]))<
                    16*1024
            end
        end
    finally
        A.finalize()
    end
end

end
