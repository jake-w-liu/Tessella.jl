module APISubdivisionAlgorithmTest
# Mesh.SubdivisionAlgorithm 1 / 2 after generate and in refine(), against
# Gmsh 4.15.2 (gmsh.open, option.setNumber, mesh.generate / mesh.refine):
# node counts and element counts by type.
using Test, Tessella

const _A=Tessella.API
const _SQUARE="""Point(1)={0,0,0,0.5};Point(2)={1,0,0,0.5};Point(3)={1,1,0,0.5};Point(4)={0,1,0,0.5};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
"""
const _GRID=_SQUARE*"Transfinite Curve{:}=3;Transfinite Surface{1};Recombine Surface{1};\n"

function _counts()
    out=Dict{Int,Int}()
    for dim in 0:3
        types,tags,_=_A.mesh.get_elements(dim)
        for (type,cells) in zip(types,tags)
            out[Int(type)]=get(out,Int(type),0)+length(cells)
        end
    end
    return length(_A.mesh.get_nodes()[1]),out
end

function _run(source,dim,algorithm;refine=false)
    _A.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"case.geo")
            write(path,source)
            _A.open_geo!(path;mesh_dim=0)
            refine || _A.option("Mesh.SubdivisionAlgorithm",algorithm)
            _A.mesh.generate(dim)
            if refine
                _A.option("Mesh.SubdivisionAlgorithm",algorithm)
                _A.mesh.refine()
            end
            @test validate(_A.mesh.get()).ok
            return _counts()
        end
    finally
        _A.finalize()
    end
end

@testset "Mesh.SubdivisionAlgorithm as in Gmsh" begin
    # Triangles into three quadrangles after generate(2).
    @test _run(_SQUARE,2,1)==(51,Dict(1=>16,3=>42,15=>4))
    # Algorithm 1 leaves a volume mesh alone, 2 a surface mesh.
    @test _run(_GRID*"Extrude{0,0,1}{Surface{1};Layers{2};Recombine;}\n",3,1)==
        (27,Dict(1=>24,3=>24,5=>8,15=>8))
    @test _run(_SQUARE,2,2)==(12,Dict(1=>8,2=>14,15=>4))
    # refine() splits by the option too, here into quadrangles.
    @test _run(_SQUARE,2,2;refine=true)==(51,Dict(1=>16,3=>42,15=>4))
    # Prisms into six hexahedra, the lateral quadrangles into four.
    @test _run(_SQUARE*"Extrude{0,0,1}{Surface{1};Layers{2};}\n",3,2)==
        (467,Dict(1=>48,3=>180,5=>336,15=>8))
    # Hexahedra into eight.
    @test _run(_GRID*"Extrude{0,0,1}{Surface{1};Layers{2};Recombine;}\n",3,2)==
        (125,Dict(1=>48,3=>96,5=>64,15=>8))
end

end
