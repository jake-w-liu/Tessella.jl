module GeoT3SweptParametrizationTest

using Test,Tessella

# Gmsh's tutorial t3: a translated, a rotated and a twisted (translate and
# rotate) structured extrusion. Every swept-surface node's parametric
# coordinates are Gmsh's inversion of the exactly swept surface; a few are
# pinned here to Gmsh 4.15.2's values.
const _T3="""
lc = 1e-2;
Point(1) = {0, 0, 0, lc};
Point(2) = {.1, 0,  0, lc};
Point(3) = {.1, .3, 0, lc};
Point(4) = {0,  .3, 0, lc};
Line(1) = {1, 2};
Line(2) = {3, 2};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {4, 1, -2, 3};
Plane Surface(1) = {1};
Physical Curve(5) = {1, 2, 4};
Physical Surface("My surface") = {1};
h = 0.1;
Extrude {0,0,h} { Surface{1}; Layers{ {8,2}, {0.5,1} }; }
Extrude { {0,1,0} , {-0.1,0,0.1} , -Pi/2 } { Surface{28}; Layers{7}; Recombine; }
out[] = Extrude { {-2*h,0,0}, {1,0,0} , {0,0.15,0.25} , Pi/2 } {
  Surface{50}; Layers{10}; Recombine;
};
"""

# Inverted parameters come out of a Newton iteration. They are Gmsh's bit for
# bit on x86-64, where these values were taken; another processor's arithmetic
# (fused multiply-add on ARM) moves their last digits.
_same_parameters(a,b)=Sys.ARCH===:x86_64 ? a==b : isapprox(a,b;atol=1e-12,rtol=0)

_node(tags,coords,params,tag)=begin
    i=findfirst(==(tag),tags)
    (coords[3i-2:3i],params[2i-1:2i])
end

@testset "t3 swept surfaces carry Gmsh's parametric coordinates" begin
    A=Tessella.API
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            path=joinpath(directory,"t3.geo");write(path,_T3)
            A.open(path)
            A.mesh.generate(2)
            expected=Dict(
                # translated lateral
                15=>(759,[0.0,0.29000000000001513,0.00625],
                     [0.03333333333328277,0.0625]),
                # rotated lateral
                37=>(1768,[-0.009903113209758091,0.29000000000001513,0.14338837391175582],
                     [0.03333333333328277,0.28571428571428614]),
                # twisted laterals
                59=>(2547,[-0.12000000000000001,0.2960980909353458,0.2225164080758778],
                     [0.03333333333328294,0.10000000000000003]),
                67=>(2890,[-0.14,0.0014012379999578939,0.2542904466022738],
                     [0.0333333333334272,0.19999999999999993]))
            for (surface,(tag,xyz,uv)) in expected
                tags,coords,params=A.mesh.get_nodes(2,surface,false,true)
                @test length(tags)==(surface==37 ? 174 : 261)
                node_xyz,node_uv=_node(tags,coords,params,tag)
                @test node_xyz==xyz
                @test _same_parameters(node_uv,uv)
            end
        end
    finally
        A.finalize()
    end
end

end
