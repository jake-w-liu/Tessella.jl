module TransfiniteOrdinaryCurvesTest
# A Transfinite Surface or Volume bounded by curves without a transfinite
# record takes the nodes those curves carry, as Gmsh 4.15.2's
# MeshTransfiniteSurface / MeshTransfiniteVolume do (examples/api/square.geo).
# Counts are Gmsh's: nodes, then (element type, count).
using Test, Tessella

const _A=Tessella.API

function _mesh(source,dim)
    _A.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"case.geo")
            write(path,source)
            _A.open(path)
            _A.mesh.generate(dim)
            types,tags,_=_A.mesh.get_elements(dim)
            return length(_A.mesh.get_nodes()[1]),
                [(Int(type),length(cells)) for (type,cells) in zip(types,tags)]
        end
    finally
        _A.finalize()
    end
end

const _BOX="""Point(1)={0,0,0,LC};Point(2)={1,0,0,LC};Point(3)={1,1,0,LC};Point(4)={0,1,0,LC};
Point(5)={0,0,1,LC};Point(6)={1,0,1,LC};Point(7)={1,1,1,LC};Point(8)={0,1,1,LC};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={1,5};Line(10)={2,6};Line(11)={3,7};Line(12)={4,8};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};Curve Loop(2)={5,6,7,8};Plane Surface(2)={2};
Curve Loop(3)={1,10,-5,-9};Plane Surface(3)={3};Curve Loop(4)={2,11,-6,-10};Plane Surface(4)={4};
Curve Loop(5)={3,12,-7,-11};Plane Surface(5)={5};Curve Loop(6)={4,9,-8,-12};Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6};Volume(1)={1};
Transfinite Curve{1,2,3,4,5,6,7,8,9,10,11}=N;Transfinite Surface{:};
Transfinite Volume{1}={1,2,3,4,5,6,7,8};
"""

@testset "transfinite on ordinarily meshed curves" begin
    @test _mesh("""SetFactory("OpenCASCADE");Rectangle(1)={0,0,0,1,1};
        Mesh.MeshSizeMin=1;Mesh.MeshSizeMax=1;Transfinite Surface{1};""",2)==(4,[(2,2)])
    @test _mesh("""SetFactory("OpenCASCADE");Rectangle(1)={0,0,0,2,1};
        Mesh.MeshSizeMax=0.3;Transfinite Surface{1};Recombine Surface{1};""",2)==(77,[(3,60)])
    # Curve 12 gets 4 nodes from the point sizes, matching the others.
    @test _mesh(replace(_BOX,"LC"=>"0.4","N"=>"4"),3)==(64,[(4,162)])
    @test _mesh(replace(_BOX,"LC"=>"0.5","N"=>"3"),3)==(27,[(4,48)])
    # With 3 and 4 nodes the opposite edges differ: Gmsh stops too
    # ("Incompatible surface 5 in transfinite volume 1").
    @test_throws ArgumentError _mesh(replace(_BOX,"LC"=>"0.4","N"=>"3"),3)
end

end
