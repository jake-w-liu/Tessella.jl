using Test
using Tessella

const _XA = Tessella.API

function _xgeo_mesh(source::AbstractString, dim::Int)
    dir = mktempdir()
    path = joinpath(dir, "case.geo")
    write(path, source)
    _XA.initialize()
    try
        _XA.option("General.Terminal", 0)
        _XA.open(path)
        _XA.mesh.generate(dim)
        counts = Dict{Int,Int}()
        for d in 0:dim
            types, tags, _ = _XA.mesh.get_elements(d)
            for (t, tg) in zip(types, tags)
                counts[t] = get(counts, t, 0) + length(tg)
            end
        end
        nt, _, _ = _XA.mesh.get_nodes()
        return (nodes=length(nt), counts=counts)
    finally
        _XA.finalize()
    end
end

const _XSQUARE = """
Point(1)={0,0,0,0.2}; Point(2)={1,0,0,0.2}; Point(3)={1,1,0,0.2}; Point(4)={0,1,0,0.2};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
"""

@testset "extrusion Gmsh 4.15.2 parity regressions" begin
    # Boundary-layer extrusion of a recombined (quadrangle) source face used to
    # throw a MethodError in the quad normal helper. Gmsh 4.15.2: 144 nodes,
    # 110 quadrangles, 75 hexahedra.
    bl = _xgeo_mesh(_XSQUARE * """
        Recombine Surface{1};
        Transfinite Curve{1:4}=6; Transfinite Surface{1};
        e[] = Extrude { Surface{1}; Layers{{1,1,1},{0.01,0.03,0.07}}; Recombine; };
        """, 3)
    @test bl.nodes == 144
    @test bl.counts[3] == 110
    @test bl.counts[5] == 75

    # OCC-factory `.geo` translation extrusion with Layers used to throw a
    # MethodError (`_affine_translation` missing its caller argument).
    # Gmsh 4.15.2: 150 nodes, 112 quadrangles, 80 hexahedra.
    oc = _xgeo_mesh("""
        SetFactory("OpenCASCADE");
        Rectangle(1)={0,0,0,1,0.5};
        Transfinite Curve{:}=5; Transfinite Surface{1}; Recombine Surface{1};
        Extrude {0,0,0.7} { Surface{1}; Layers{{2,3},{0.4,1}}; Recombine; }
        """, 3)
    @test oc.nodes == 150
    @test oc.counts[3] == 112
    @test oc.counts[5] == 80

end
