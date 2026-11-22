module IOWindowsScanfTest

using Test, Tessella

const API = Tessella.API
const LIBM = Tessella.GmshLibm

# Gmsh 4.15.2's Windows build reads `%lf` through the 80-bit long double and
# then rounds to double, so "-0.529113" reads as -0.5291129999999999. Its
# MSH 2/4, MESH, UNV, VTK and STL readers use `%lf`; `.geo` and free-field
# BDF numbers round correctly. These values were read with the Gmsh 4.15.2
# Python wheel on Windows; elsewhere Gmsh rounds correctly everywhere.
const TWICE = Sys.iswindows() ? -0.5291129999999999 : -0.529113

@testset "%lf double rounding" begin
    @test LIBM._gm_scanf_double("-0.529113") === TWICE
    @test LIBM._gm_scanf_double("0.5") === 0.5
    @test LIBM._gm_scanf_double("1e300") === 1e300
    @test LIBM._gm_scanf_double("-0x1.0ee7e62dc6e2bp-1") === -0.529113
    exact(s) = Sys.iswindows() ?
        Float64(setprecision(() -> parse(BigFloat, s), BigFloat, 64)) : parse(Float64, s)
    for s in ("-0.529113", "0.1", "123.456e-7", "7.0000000000000001", "1e23",
              "4.35", "0.3", "2.675", "1234567890123456789e-30", "9007199254740993")
        @test LIBM._gm_scanf_double(s) === exact(s)
    end
end

const V = "-0.529113"
const FILES = [
    ("a.msh", "\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Nodes\n1 2 1 2\n1 1 0 2\n1\n2\n$V 0 0\n1 0 0\n" *
              "\$EndNodes\n\$Elements\n1 1 1 1\n1 1 1 1\n1 1 2\n\$EndElements\n", TWICE),
    ("b.msh", "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 $V 0 0\n2 1 0 0\n\$EndNodes\n" *
              "\$Elements\n1\n1 1 2 0 1 1 2\n\$EndElements\n", TWICE),
    ("e.mesh", " MeshVersionFormatted 2\n Dimension\n 3\n Vertices\n 2\n $V 0 0 0\n 1 0 0 0\n" *
               " Edges\n 1\n 1 2 1\n End\n", TWICE),
    ("f.unv", "    -1\n  2411\n         1         1         1        11\n   $V   0.0   0.0\n" *
              "         2         1         1        11\n   1.0   0.0   0.0\n    -1\n    -1\n  2412\n" *
              "         1        21         1         1         7         2\n         0         0         0\n" *
              "         1         2\n    -1\n", TWICE),
    ("g.vtk", "# vtk DataFile Version 2.0\nx\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 2 double\n$V 0 0\n" *
              "1 0 0\nCELLS 1 3\n2 0 1\nCELL_TYPES 1\n3\n", TWICE),
    ("h.bdf", "GRID,1,0,$V,0.,0.\nGRID,2,0,1.,0.,0.\nCBAR,1,1,1,2\nENDDATA\n", -0.529113),
]

@testset "ASCII readers read numbers as Gmsh does" begin
    mktempdir() do dir
        for (name, text, expected) in FILES
            path = joinpath(dir, name); write(path, text)
            API.initialize()
            try
                API.option("General.Terminal", 0)
                API.open(path)
                @test API.mesh.get_nodes()[2][1] === expected
            finally
                API.finalize()
            end
        end
        path = joinpath(dir, "c.geo"); write(path, "Point(1) = {$V, 0, 0};\n")
        API.initialize()
        try
            API.option("General.Terminal", 0)
            API.open(path)
            @test API.model.get_value(0, 1, Float64[])[1] === -0.529113
        finally
            API.finalize()
        end
    end
end

end
