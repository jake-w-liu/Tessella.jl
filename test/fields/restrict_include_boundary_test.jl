module RestrictIncludeBoundaryTest

using Test
using Tessella

# Restrict/Constant fields default to IncludeBoundary=1: the bounding points of
# a listed curve take the restricted size too, so curves meeting it start
# finer at the shared endpoint (upstream `getRestrictEntities`). Gmsh 4.15.2
# node counts per curve: 51, 14, 35, 15.
const _RESTRICT_GEO = """
Point(1)={0,0,0,0.1}; Point(2)={1,0,0,0.1}; Point(3)={1,1,0,0.1}; Point(4)={0,1,0,0.1};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Mesh.MeshSizeExtendFromBoundary=0; Mesh.MeshSizeFromPoints=0; Mesh.MeshSizeFromCurvature=0;
Field[1]=MathEval; Field[1].F="0.02"; Field[2]=Restrict; Field[2].InField=1; Field[2].CurvesList={1};
Field[3]=MathEval; Field[3].F="0.1"; Field[4]=Min; Field[4].FieldsList={2,3};
Field[5]=Constant; Field[5].VIn=0.03; Field[5].VOut=0.1; Field[5].CurvesList={3};
Field[6]=Min; Field[6].FieldsList={4,5}; Background Field=6;
"""

@testset "Restrict/Constant IncludeBoundary reaches curve endpoints" begin
    path = joinpath(mktempdir(), "restrict.geo")
    write(path, _RESTRICT_GEO)
    A = Tessella.API
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.open(path)
        A.mesh.generate(1)
        counts = [length(A.mesh.get_nodes(1, c, true)[1]) for c in 1:4]
        @test counts == [51, 14, 35, 15]
    finally
        A.finalize()
    end
    # The `.geo` executor's own field builder follows the same rule.
    write(path, _RESTRICT_GEO * "Mesh 1;\n")
    execution = execute_geo(path)
    @test Tessella.MeshTypes.nnodes(execution.mesh) == 51 + 14 + 35 + 15 - 4
end

end # module
