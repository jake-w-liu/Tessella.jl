# `-0` unrolled output details pinned against the Gmsh 4.15.2 binary on the
# same script: a point declared without a size writes `1e+22`, a reversed
# `Transfinite Curve {-1}` keeps its sign, and embedded curves are spelled
# `Line {..} In Surface {..}`.
using Test
using Tessella

@testset "CLI -0 unrolled details" begin
    dir=mktempdir()
    geo=joinpath(dir,"details.geo")
    write(geo,"""
        Point(1) = {0,0,0,0.5}; Point(2) = {1,0,0,0.5}; Point(3) = {1,1,0}; Point(4) = {0,1,0,0.5};
        Point(5) = {0.2,0.5,0,0.5}; Point(6) = {0.8,0.5,0,0.5};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Line(5) = {5,6};
        Curve Loop(1) = {1,2,3,4};
        Plane Surface(1) = {1};
        Line{5} In Surface{1};
        Transfinite Curve{-1, 3} = 5 Using Progression 1.2;
        """)
    out=joinpath(dir,"details.geo_unrolled")
    Tessella.CLI.main([geo,"-0","-o",out])
    expected="""
        cl__1 = 0.5;
        cl__2 = 1e+22;
        Point(1) = {0, 0, 0, cl__1};
        Point(2) = {1, 0, 0, cl__1};
        Point(3) = {1, 1, 0, cl__2};
        Point(4) = {0, 1, 0, cl__1};
        Point(5) = {0.2, 0.5, 0, cl__1};
        Point(6) = {0.8, 0.5, 0, cl__1};
        Line(1) = {1, 2};
        Transfinite Curve {-1} = 5 Using Progression 1.2;
        Line(2) = {2, 3};
        Line(3) = {3, 4};
        Transfinite Curve {3} = 5 Using Progression 1.2;
        Line(4) = {4, 1};
        Line(5) = {5, 6};
        Curve Loop(1) = {1, 2, 3, 4};
        Plane Surface(1) = {1};
        Line {5} In Surface {1};
        """
    @test replace(read(out,String),"\r\n"=>"\n")==expected
end
