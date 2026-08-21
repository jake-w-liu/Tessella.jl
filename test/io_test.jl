# ── Stage-0 CRC suite: .msh v2/v4 round-trip, STL ingest, .geo scan ─────────────
#
# Correctness  : hand-written gmsh-format samples (an independent transcription of
#                the documented format) are parsed and checked node-by-node /
#                element-by-element; round-trips are checked by mesh_crc equality.
# Robustness   : non-contiguous node tags, physical groups, empty mesh, both
#                format versions, ASCII+binary STL, welding of coincident verts.
# Completeness : cross-format (v2→v4→v2) preserves connectivity CRC exactly.

using Test
import Tessella
using Tessella.MeshTypes
using Tessella.IO

# Positively-oriented 6-tet Kuhn cube (same helper as meshtypes_test, redefined
# locally so this file runs standalone).
function _cube()
    coords = Float64[0 1 1 0 0 1 1 0; 0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1]
    raw = [ (1,2,3,7), (1,2,6,7), (1,4,3,7), (1,4,8,7), (1,5,6,7), (1,5,8,7) ]
    tets = Matrix{Int32}(undef, 4, length(raw))
    tags = Int32[]
    for (k, t) in enumerate(raw)
        f(i) = (coords[1,i],coords[2,i],coords[3,i])
        v = tet_signed_volume(f(t[1]),f(t[2]),f(t[3]),f(t[4]))
        tt = v < 0 ? (t[1],t[2],t[4],t[3]) : t
        tets[:,k] = Int32[tt...]
        push!(tags, Int32(4))     # physical tag 4 ("case")
    end
    return Mesh(coords; tets=tets, tet_tag=tags)
end

@testset "IO (Stage 0)" begin
    dir = mktempdir()

    @testset "v2 round-trip preserves connectivity CRC" begin
        m = _cube()
        p = joinpath(dir, "cube_v2.msh")
        write_msh(p, m; version=2.2, physical_names=Dict((3,4)=>"case"))
        f = read_msh(p)
        @test mesh_crc(f.mesh).sha == mesh_crc(m).sha
        @test ntets(f.mesh) == 6
        @test all(f.mesh.tet_tag .== 4)
        @test f.physical_names[(3,4)] == "case"
    end

    @testset "physical name with interior whitespace round-trips verbatim" begin
        # regression: the reader must NOT collapse runs of spaces/tabs inside a quoted name
        # (split/join-on-whitespace did; the quoted-substring parse preserves it exactly).
        m = _cube(); p = joinpath(dir, "cube_names.msh")
        for (v, nm) in ((2.2, "air  region"), (4.1, "co ax\tpin"))
            write_msh(p, m; version=v, physical_names=Dict((3,4)=>nm))
            @test read_msh(p).physical_names[(3,4)] == nm
        end
        nm = "quoted \"name\" \\ path\nline two"
        write_msh(p, m; version=2.2, physical_names=Dict((3,4)=>nm))
        @test read_msh(p).physical_names[(3,4)] == nm
    end

    @testset "MSH format errors are explicit and writes are atomic" begin
        m=_cube(); p=joinpath(dir,"atomic.msh")
        @test_throws ArgumentError write_msh(p,m;version=3.0)
        write(p,"sentinel")
        @test_throws ArgumentError write_msh(p,m;physical_names=Dict((3,4)=>17))
        @test read(p,String)=="sentinel"                 # failed write did not truncate target

        missing=joinpath(dir,"missing_format.msh")
        write(missing,"\$Nodes\n0\n\$EndNodes\n")
        @test_throws ArgumentError read_msh(missing)
        unsupported=joinpath(dir,"unsupported.msh")
        write(unsupported,"\$MeshFormat\n3.0 0 8\n\$EndMeshFormat\n")
        @test_throws ArgumentError read_msh(unsupported)
        duplicate=joinpath(dir,"duplicate_node.msh")
        write(duplicate,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n1 1 0 0\n\$EndNodes\n")
        @test_throws ArgumentError read_msh(duplicate)
        unknown=joinpath(dir,"unknown_node.msh")
        write(unknown,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n1\n1 0 0 0\n\$EndNodes\n\$Elements\n1\n1 1 0 1 2\n\$EndElements\n")
        @test_throws ArgumentError read_msh(unknown)

        badentity=joinpath(dir,"bad_entity_header.msh")
        write(badentity,"\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Entities\n0 0\n\$EndEntities\n")
        @test_throws ArgumentError read_msh(badentity)
        negentity=joinpath(dir,"negative_entity.msh")
        write(negentity,"\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Entities\n-1 0 0 0\n\$EndEntities\n")
        @test_throws ArgumentError read_msh(negentity)
        zeroentity=joinpath(dir,"zero_entity.msh")
        write(zeroentity,"\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Entities\n0 0 0 1\n0 0 0 0 1 1 1 0 0\n\$EndEntities\n")
        @test_throws ArgumentError read_msh(zeroentity)
        missingbound=joinpath(dir,"missing_boundary_entity.msh")
        write(missingbound,"\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Entities\n0 1 0 0\n1 0 0 0 1 1 1 0 1 9\n\$EndEntities\n")
        @test_throws ArgumentError read_msh(missingbound)
        dupv4=joinpath(dir,"duplicate_v4_node.msh")
        write(dupv4,"\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Nodes\n1 2 1 2\n3 1 0 2\n1\n1\n0 0 0\n1 0 0\n\$EndNodes\n")
        @test_throws ArgumentError read_msh(dupv4)
        dupelem=joinpath(dir,"duplicate_element.msh")
        write(dupelem,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n2\n1 1 0 1 2\n1 1 0 1 2\n\$EndElements\n")
        @test_throws ArgumentError read_msh(dupelem)
        zeroelem=joinpath(dir,"zero_element.msh")
        write(zeroelem,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n1\n0 1 0 1 2\n\$EndElements\n")
        @test_throws ArgumentError read_msh(zeroelem)
        bigphys=joinpath(dir,"big_physical.msh")
        write(bigphys,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n1\n1 1 1 999999999999 1 2\n\$EndElements\n")
        @test_throws ArgumentError read_msh(bigphys)
        unsupported_elem=joinpath(dir,"unsupported_element.msh")
        write(unsupported_elem,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 1 1 0\n4 0 1 0\n\$EndNodes\n\$Elements\n1\n1 3 0 1 2 3 4\n\$EndElements\n")
        @test_throws ArgumentError read_msh(unsupported_elem)
        multiphys=joinpath(dir,"multiple_physical_groups.msh")
        write(multiphys,"\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Entities\n0 0 0 1\n1 0 0 0 1 1 1 2 11 12 0\n\$EndEntities\n\$Nodes\n0 0 0 0\n\$EndNodes\n")
        @test_throws ArgumentError read_msh(multiphys)

        dirty=_cube(); dirty.coords[1,1]=NaN; write(p,"sentinel")
        @test_throws ArgumentError write_msh(p,dirty)
        @test read(p,String)=="sentinel"
    end

    @testset "v4.1 round-trip preserves connectivity CRC" begin
        m = _cube()
        p = joinpath(dir, "cube_v4.msh")
        write_msh(p, m; version=4.1, physical_names=Dict((3,4)=>"case"))
        f = read_msh(p)
        @test mesh_crc(f.mesh).sha == mesh_crc(m).sha
        @test ntets(f.mesh) == 6
        @test all(f.mesh.tet_tag .== 4)
        @test f.physical_names[(3,4)] == "case"
    end

    @testset "cross-format v2 → v4 → v2 keeps CRC" begin
        m = _cube()
        p2 = joinpath(dir, "x_v2.msh"); p4 = joinpath(dir, "x_v4.msh"); p2b = joinpath(dir, "x_v2b.msh")
        write_msh(p2, m; version=2.2)
        m2 = read_msh(p2).mesh
        write_msh(p4, m2; version=4.1)
        m4 = read_msh(p4).mesh
        write_msh(p2b, m4; version=2.2)
        m2b = read_msh(p2b).mesh
        @test mesh_crc(m).sha == mesh_crc(m2b).sha
    end

    @testset "mixed-dim mesh (segs+tris+tets) round-trip" begin
        coords = Float64[0 1 0 0; 0 0 1 0; 0 0 0 1]
        m = Mesh(coords;
                 segs=reshape(Int32[1,2],2,1), seg_tag=Int32[6],
                 tris=reshape(Int32[1,2,3],3,1), tri_tag=Int32[3],
                 tets=reshape(Int32[1,2,3,4],4,1), tet_tag=Int32[1])
        for ver in (2.2, 4.1)
            p = joinpath(dir, "mixed_$(ver).msh")
            write_msh(p, m; version=ver)
            f = read_msh(p)
            @test mesh_crc(f.mesh).sha == mesh_crc(m).sha
            @test nsegs(f.mesh)==1 && ntris(f.mesh)==1 && ntets(f.mesh)==1
            @test f.mesh.seg_tag[1]==6 && f.mesh.tri_tag[1]==3 && f.mesh.tet_tag[1]==1
        end
    end

    @testset "parse hand-written gmsh v2 sample (format oracle)" begin
        sample = """
        \$MeshFormat
        2.2 0 8
        \$EndMeshFormat
        \$PhysicalNames
        1
        2 5 "surf"
        \$EndPhysicalNames
        \$Nodes
        3
        10 0.0 0.0 0.0
        20 2.0 0.0 0.0
        30 0.0 3.0 0.0
        \$EndNodes
        \$Elements
        1
        1 2 2 5 1 10 20 30
        \$EndElements
        """
        p = joinpath(dir, "sample_v2.msh"); write(p, sample)
        f = read_msh(p)
        @test nnodes(f.mesh) == 3
        @test ntris(f.mesh) == 1
        @test f.mesh.tri_tag[1] == 5
        @test f.physical_names[(2,5)] == "surf"
        # non-contiguous tags (10,20,30) relabelled but geometry intact:
        @test triangle_area(node(f.mesh,f.mesh.tris[1,1]), node(f.mesh,f.mesh.tris[2,1]), node(f.mesh,f.mesh.tris[3,1])) ≈ 3.0
    end

    @testset "parse hand-written gmsh v4.1 sample (format oracle)" begin
        sample = """
        \$MeshFormat
        4.1 0 8
        \$EndMeshFormat
        \$PhysicalNames
        1
        3 1 "vol"
        \$EndPhysicalNames
        \$Entities
        0 0 0 1
        1 0 0 0 1 1 1 1 1 0
        \$EndEntities
        \$Nodes
        1 4 1 4
        3 1 0 4
        1
        2
        3
        4
        0 0 0
        1 0 0
        0 1 0
        0 0 1
        \$EndNodes
        \$Elements
        1 1 1 1
        3 1 4 1
        1 1 2 3 4
        \$EndElements
        """
        p = joinpath(dir, "sample_v4.msh"); write(p, sample)
        f = read_msh(p)
        @test nnodes(f.mesh) == 4
        @test ntets(f.mesh) == 1
        @test f.mesh.tet_tag[1] == 1
        @test f.physical_names[(3,1)] == "vol"
        @test tet_volume(node(f.mesh,1),node(f.mesh,2),node(f.mesh,3),node(f.mesh,4)) ≈ 1/6
        @test validate(f.mesh).ok
    end

    @testset "STL ingest welds coincident vertices (ASCII)" begin
        # tetrahedron surface: 4 triangles, 4 unique vertices (each shared by 3 faces)
        A=(0.0,0.0,0.0); B=(1.0,0.0,0.0); C=(0.0,1.0,0.0); D=(0.0,0.0,1.0)
        faces = [(A,B,C),(A,B,D),(A,C,D),(B,C,D)]
        buf = IOBuffer()
        println(buf, "solid tet")
        for (p,q,r) in faces
            println(buf, "facet normal 0 0 0"); println(buf, "  outer loop")
            for v in (p,q,r); println(buf, "    vertex $(v[1]) $(v[2]) $(v[3])"); end
            println(buf, "  endloop"); println(buf, "endfacet")
        end
        println(buf, "endsolid tet")
        pstl = joinpath(dir, "tet.stl"); write(pstl, String(take!(buf)))
        m = read_stl(pstl)
        @test nnodes(m) == 4          # welded from 12 raw vertices
        @test ntris(m) == 4
        be, maxinc = boundary_edges(m.tris)
        @test isempty(be)             # closed surface: no boundary edges
        @test maxinc == 2             # manifold
    end

    @testset "STL ingest (binary)" begin
        A=(0.0f0,0.0f0,0.0f0); B=(1.0f0,0.0f0,0.0f0); C=(0.0f0,1.0f0,0.0f0); D=(0.0f0,0.0f0,1.0f0)
        faces = [(A,B,C),(A,B,D),(A,C,D),(B,C,D)]
        buf = IOBuffer()
        write(buf, zeros(UInt8, 80))            # header
        write(buf, UInt32(length(faces)))       # triangle count
        for (p,q,r) in faces
            write(buf, 0.0f0, 0.0f0, 0.0f0)     # normal
            for v in (p,q,r); write(buf, v[1], v[2], v[3]); end
            write(buf, UInt16(0))               # attribute
        end
        pstl = joinpath(dir, "tet_bin.stl"); write(pstl, take!(buf))
        m = read_stl(pstl)
        @test nnodes(m) == 4
        @test ntris(m) == 4
    end

    @testset "STL robustness: far-from-origin welding + ASCII detection (BOM / UTF-8 name)" begin
        write_ascii_tri(path, name, tri) = open(path, "w") do io
            println(io, "solid ", name); println(io, "facet normal 0 0 0"); println(io, "  outer loop")
            for v in tri; println(io, "    vertex $(v[1]) $(v[2]) $(v[3])"); end
            println(io, "  endloop"); println(io, "endfacet"); println(io, "endsolid ", name)
        end
        # far from the origin: welding must quantize relative to the bbox, not the
        # absolute coordinate (else round(Int, coord/tol) overflows Int64 → crash).
        pfar = joinpath(dir, "far.stl")
        write_ascii_tri(pfar, "far", ((1e11,0.0,0.0),(1e11+1,0.0,0.0),(1e11,1.0,0.0)))
        mf = read_stl(pfar)                                # must not throw InexactError
        @test ntris(mf) == 1 && nnodes(mf) == 3

        # ASCII STL prefixed with a UTF-8 BOM must not be misdetected as binary.
        bom = joinpath(dir, "bom.stl")
        open(bom, "w") do io
            write(io, UInt8[0xEF,0xBB,0xBF])
            println(io, "solid tet"); println(io, "facet normal 0 0 0"); println(io, "  outer loop")
            println(io, "    vertex 0 0 0"); println(io, "    vertex 1 0 0"); println(io, "    vertex 0 1 0")
            println(io, "  endloop"); println(io, "endfacet"); println(io, "endsolid tet")
        end
        mb = read_stl(bom); @test ntris(mb) == 1 && nnodes(mb) == 3

        # ASCII STL whose solid name carries a non-ASCII (UTF-8) char, likewise.
        pacc = joinpath(dir, "acc.stl")
        write_ascii_tri(pacc, "pièce_de_test", ((0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0)))
        ma = read_stl(pacc); @test ntris(ma) == 1 && nnodes(ma) == 3

        # a genuine BINARY STL whose 80-byte header text begins with "solid" is still
        # detected as binary via the size check (not misread as ASCII).
        buf = IOBuffer(); hdr = zeros(UInt8, 80); hdr[1:5] = collect(codeunits("solid")); write(buf, hdr)
        write(buf, UInt32(1)); write(buf, 0f0,0f0,0f0, 0f0,0f0,0f0, 1f0,0f0,0f0, 0f0,1f0,0f0, UInt16(0))
        pbs = joinpath(dir, "binsolid.stl"); write(pbs, take!(buf))
        @test ntris(read_stl(pbs)) == 1
    end

    @testset "STL parser and geometric welding contracts" begin
        @test_throws ArgumentError read_stl(joinpath(dir,"tet.stl");merge_tol=-1)
        malformed=joinpath(dir,"malformed.stl")
        write(malformed,"solid x\nvertex 0 0 0\nvertex 1 0 0\nendsolid x\n")
        @test_throws ArgumentError read_stl(malformed)
        nonfinite=joinpath(dir,"nonfinite.stl")
        write(nonfinite,"solid x\nvertex NaN 0 0\nvertex 1 0 0\nvertex 0 1 0\nendsolid x\n")
        @test_throws ArgumentError read_stl(nonfinite)
        soup=joinpath(dir,"unframed_soup.stl")
        write(soup,"solid x\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendsolid x\n")
        @test_throws ArgumentError read_stl(soup)
        four=joinpath(dir,"four_vertex_facet.stl")
        write(four,"solid x\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nvertex 1 1 0\nendloop\nendfacet\nendsolid x\n")
        @test_throws ArgumentError read_stl(four)

        # Two vertices within tolerance but on opposite hash-cell sides must weld.
        within = NTuple{9,Float64}[
            (0.099,0.1,0.1, 10.,0.,0., 10.,1.,0.),
            (0.101,0.1,0.1, 0.,10.,0., 1.,10.,0.)]
        mw = Tessella.IO._weld_triangles(within,0.0071) # abs tol ≈0.100; gap=0.002
        @test nnodes(mw)==5 && ntris(mw)==2

        # Sharing a hash neighbourhood is insufficient: Euclidean distance is the
        # contract, so a 3-D diagonal farther than tolerance stays distinct.
        apart = NTuple{9,Float64}[
            (0.10,0.10,0.10, 10.,0.,0., 10.,1.,0.),
            (0.18,0.18,0.18, 0.,10.,0., 1.,10.,0.)]
        ma = Tessella.IO._weld_triangles(apart,0.0071) # distance≈0.139 > tol≈0.100
        @test nnodes(ma)==6 && ntris(ma)==2
    end

    @testset "read_geo_params on the enclosure fixture" begin
        gp = read_geo_params(joinpath(@__DIR__, "fixtures", "enclosure_coax_junction.geo"))
        @test gp.mesh_size_min ≈ 0.00026669999999999933 rtol=1e-9
        @test gp.mesh_size_max == 0.012
        @test gp.mesh_size_factor == 1.0
        @test gp.random_seed == 1
        # physical groups declared in the .geo (dim, tag) => name
        @test gp.physical_groups[(3,1)] == "air"
        @test gp.physical_groups[(3,2)] == "coax_pin"
        @test gp.physical_groups[(3,4)] == "case"
        @test gp.physical_groups[(2,3)] == "resistor"
        @test gp.physical_groups[(2,7)] == "radiation"
        @test gp.physical_groups[(1,6)] == "p1_line"
        # the three volumes the mesher must fill are all present
        vols = sort([name for ((d,t),name) in gp.physical_groups if d==3])
        @test vols == ["air", "case", "coax_pin"]

        @test gp.background_field == 8
        @test length(gp.fields) == 8
        @test [gp.fields[i].kind for i in 1:8] ==
              ["Distance","Threshold","Distance","Threshold",
               "Distance","Threshold","Box","Min"]
        @test gp.fields[1].options["SurfacesList"] == "{sm_coax_pin_pec[]}"
        @test gp.fields[2].options["InField"] == "1"
        @test gp.fields[8].options["FieldsList"] == "{2, 4, 6, 7}"
        @test isempty(gp.boundary_layer_fields)
        @test isnan(gp.geometry_tolerance)

        # The original four-argument constructor remains source-compatible.
        old=GeoParams(0.1,1.0,7,Dict{Tuple{Int,Int},String}())
        @test old.mesh_size_factor==1.0 && isempty(old.fields) && old.background_field==0
        @test isempty(old.boundary_layer_fields)
        @test isnan(old.geometry_tolerance)
        @test old.mesh_boundary_layer_fan_elements==5

        construction_defaults=joinpath(dir,"field_construction_defaults.geo")
        write(construction_defaults,
              "Mesh.MeshSizeFromCurvature = 0; Field[1] = AutomaticMeshSizeField;\n" *
              "Mesh.MinimumElementsPerTwoPi = 37.9; Field[2] = AutomaticMeshSizeField;\n" *
              "Mesh.BoundaryLayerFanElements = 30.8;\n" *
              "Mesh.BoundaryLayerFanPoints = 12.9; Field[3] = BoundaryLayer;\n")
        construction_params=read_geo_params(construction_defaults)
        @test construction_params.fields[1].creation_mesh_size_from_curvature==0
        @test construction_params.fields[2].creation_mesh_size_from_curvature==37
        @test construction_params.mesh_boundary_layer_fan_elements==12

        boundary=joinpath(dir,"boundary_layer_fields.geo")
        write(boundary,"Field[1] = BoundaryLayer; Field[2] = BoundaryLayer;\n" *
                       "BoundaryLayer Field = {1, 2, 1}; Background Field = {2};\n")
        boundary_params=read_geo_params(boundary)
        @test boundary_params.boundary_layer_fields==[1,2]
        @test boundary_params.background_field==2
        missing_boundary=joinpath(dir,"missing_boundary_layer_field.geo")
        write(missing_boundary,"BoundaryLayer Field = 3;\n")
        @test_throws ArgumentError read_geo_params(missing_boundary)
        malformed_boundary=joinpath(dir,"malformed_boundary_layer_field.geo")
        write(malformed_boundary,"Field[1] = BoundaryLayer; BoundaryLayer Field = {1, x};\n")
        @test_throws ArgumentError read_geo_params(malformed_boundary)
        quoted=joinpath(dir,"quoted_field_options.geo")
        write(quoted,raw"""Field[1] = Structured;
                           Field[1].FileName = "http://host/C:\new\tab\grid.bin";
                           Field[2] = ExternalProcess;
                           Field[2].CommandLine = "printf a;b"; // outside comment
                           Background Field = 1;
                           """)
        quoted_params=read_geo_params(quoted)
        @test quoted_params.fields[1].options["FileName"]==
              raw"\"http://host/C:\new\tab\grid.bin\""
        @test quoted_params.fields[2].options["CommandLine"]==raw"\"printf a;b\""
        @test quoted_params.fields[1].option_order==["FileName"]
        unterminated_quote=joinpath(dir,"unterminated_quote.geo")
        write(unterminated_quote,"Field[1] = Structured; Field[1].FileName = \"bad;\n")
        @test_throws ArgumentError read_geo_params(unterminated_quote)

        duplicate=joinpath(dir,"duplicate_field.geo")
        write(duplicate,"Field[1] = Box; Field[1] = Min; Background Field = 1;\n")
        @test_throws ArgumentError read_geo_params(duplicate)
        undeclared=joinpath(dir,"undeclared_field.geo")
        write(undeclared,"Field[2].VIn = 1;\n")
        @test_throws ArgumentError read_geo_params(undeclared)
        malformed=joinpath(dir,"malformed_size.geo")
        write(malformed,"Mesh.MeshSizeMin = nope;\n")
        @test_throws ArgumentError read_geo_params(malformed)
        invalid=joinpath(dir,"unsupported_size_expression.geo")
        write(invalid,"Mesh.MeshSizeMin = 1oops;\n")
        @test_throws ArgumentError read_geo_params(invalid)
        valid_numeric=joinpath(dir,"numeric_mesh_options.geo")
        write(valid_numeric,"Mesh.MeshSizeMin = 1 + 2 * 3; Mesh.RandomSeed = 17;\n")
        numeric_params=read_geo_params(valid_numeric)
        @test numeric_params.mesh_size_min==7.0
        @test numeric_params.random_seed==17
        geometry_option=joinpath(dir,"geometry_tolerance.geo")
        write(geometry_option,"Geometry.Tolerance = 2e-5;\n")
        @test read_geo_params(geometry_option).geometry_tolerance==2e-5
        invalid=joinpath(dir,"unsupported_seed_expression.geo")
        write(invalid,"Mesh.RandomSeed = 1oops;\n")
        @test_throws ArgumentError read_geo_params(invalid)

        @testset "finite Gmsh constant expressions" begin
            expressions=joinpath(dir,"constant_expressions.geo")
            write(expressions,raw"""
                base = 2;
                half = base^-1;
                angle = Pi / 2;
                Physical Point("expression point", Min(10, base + 1.9)) = {1};
                Mesh.MeshSizeMin = half * Sin(angle);
                Mesh.MeshSizeMax = 2^3^2;
                Mesh.MeshSizeFactor = Sqrt(9) / 3;
                Mesh.RandomSeed = 3.9;
                Geometry.Tolerance = Min(1e-4, Hypot(3, 4) * 1e-5);
                curvature = 4.9;
                Mesh.MeshSizeFromCurvature = curvature;
                field_tag = 1.9;
                Field[field_tag] = AutomaticMeshSizeField;
                curvature = 7.9;
                Mesh.MinimumElementsPerTwoPi = curvature;
                field_tag = 2.9;
                Field[field_tag] = AutomaticMeshSizeField;
                field_tag = 3.9;
                Field[field_tag] = BoundaryLayer;
                Field[field_tag].Size = half;
                Field[field_tag].hwall_n = base / 8;
                Field[field_tag].SizesList = {half, Atan2(1, 1), Round(-1.5)};
                Field[field_tag].SurfacesList = {surface_group[], base + 2};
                Field[Sin[Pi / 2] + 3] = Min;
                Field[Min[8, Sin[Pi / 2] + 3]].FieldsList = {1, Min(4, 2) + 2};
                Field[Round(4.5)] = MathEval;
                Field[Round(4.5)].F = "Sin(x) + base; // source, not a scanner expression";
                BoundaryLayer Field = {base + 1};
                Background Field = base + 2;
                """)
            params=read_geo_params(expressions)
            @test params.mesh_size_min==0.5
            @test params.mesh_size_max==512.0
            @test params.mesh_size_factor==1.0
            @test params.random_seed==3
            @test params.geometry_tolerance==5e-5
            @test params.physical_groups[(0,3)]=="expression point"
            @test params.fields[1].creation_mesh_size_from_curvature==4
            @test params.fields[2].creation_mesh_size_from_curvature==7
            @test params.fields[3].options["Size"]=="0.5"
            @test params.fields[3].options["hwall_n"]=="0.25"
            @test parse.(Float64,Tessella.SizeField._geo_list(
                params.fields[3],"SizesList")) ≈
                [0.5,0.7853981633974483,-1.0]
            @test params.fields[3].options["SurfacesList"]==
                  "{surface_group[], 4}"
            @test params.fields[4].options["FieldsList"]=="{1, 4}"
            @test params.fields[3].option_order==
                  ["Size","hwall_n","SizesList","SurfacesList"]
            @test params.fields[5].options["F"]==
                  raw"\"Sin(x) + base; // source, not a scanner expression\""
            @test params.boundary_layer_fields==[3]
            @test params.background_field==4

            # Gmsh 4.15.2 truncates then clamps RandomSeed to UInt32's range.
            for (source,expected) in (("-2",0),("0.9",0),("3.9",3),
                                      ("1e300",Int(typemax(UInt32))))
                write(expressions,"Mesh.RandomSeed = $source;\n")
                @test read_geo_params(expressions).random_seed==expected
            end
        end

        @testset "Gmsh 4.15.2 numeric-function oracle" begin
            # Expected values were independently obtained with the local 4.15.2
            # parser and Printf("%.17g", expression).
            oracles=(
                "Acos(1)"=>0.0,
                "Asin(1)"=>1.5707963267948966,
                "Atan(1)"=>0.78539816339744828,
                "Atan2(1, 1)"=>0.78539816339744828,
                "Ceil(1.1)"=>2.0,
                "Cos(0)"=>1.0,
                "Cosh(1)"=>1.5430806348152437,
                "Exp(1)"=>2.7182818284590451,
                "Fabs(-2)"=>2.0,
                "Abs(-2)"=>2.0,
                "Floor(1.9)"=>1.0,
                "Hypot(3, 4)"=>5.0,
                "Log(Exp(1))"=>1.0,
                "Log10(1000)"=>3.0,
                "Max(3, 2)"=>3.0,
                "Min(3, 2)"=>2.0,
                "Modulo(7, 4)"=>3.0,
                "Fmod(-7, 4)"=>-3.0,
                "Round(-1.5)"=>-1.0,
                "Sqrt(9)"=>3.0,
                "Sin[Pi / 2]"=>1.0,
                "Sinh(1)"=>1.1752011936438014,
                "Step(-0.)"=>1.0,
                "Tan(Pi / 4)"=>0.99999999999999989,
                "Tanh(1)"=>0.76159415595576485)
            oracle_file=joinpath(dir,"numeric_function_oracle.geo")
            for (expression,expected) in oracles
                write(oracle_file,
                    "Field[1] = Box; Field[1].VIn = $expression; Background Field = 1;\n")
                parsed=read_geo_params(oracle_file)
                actual=parse(Float64,parsed.fields[1].options["VIn"])
                @test actual ≈ expected rtol=4eps(Float64) atol=4eps(Float64)
            end
            precedence=("1 + 2 * 3"=>7.0,"(1 + 2) * 3"=>9.0,
                        "2^3^2"=>512.0,"-2^2"=>-4.0,
                        "2^-2^2"=>0.0625,"- -2"=>2.0,"+ +2"=>2.0,
                        "7%(-4)"=>3.0,"1.e-2"=>0.01)
            for (expression,expected) in precedence
                write(oracle_file,
                    "Mesh.MeshSizeMin = $expression;\n")
                @test read_geo_params(oracle_file).mesh_size_min==expected
            end
        end

        @testset "finite Gmsh list ranges" begin
            ranges=joinpath(dir,"constant_list_ranges.geo")
            write(ranges,raw"""
                first = 1;
                last = 5;
                increment = 2;
                Physical Point("ranged membership", 11) = {first:last:increment};
                Field[1] = Box;
                Field[2] = Box;
                Field[3] = Box;
                Field[4] = Box;
                Field[5] = Box;
                Field[6] = Min;
                Field[6].FieldsList = {first:last:increment, 5:-2:1, 4};
                Field[7] = BoundaryLayer;
                Field[7].SizesList = {Min(0, 1):3 / 4:1 / 4};
                Field[7].PointsList = {point_group[], 5:1:-2};
                Field[8] = BoundaryLayer;
                Field[9] = Min;
                Field[9].FieldsList = {5:-2:1};
                BoundaryLayer Field = {7:8, 8:7};
                Background Field = {5:-2:1, 6};
                """)
            params=read_geo_params(ranges)
            @test params.physical_groups[(0,11)]=="ranged membership"
            @test params.fields[6].options["FieldsList"]=="{1, 3, 5, 4}"
            @test parse.(Float64,Tessella.SizeField._geo_list(
                params.fields[7],"SizesList"))==[0.0,0.25,0.5,0.75]
            @test params.fields[7].options["PointsList"]==
                  "{point_group[], 5, 3, 1}"
            @test params.fields[7].option_order==["SizesList","PointsList"]
            @test params.fields[9].options["FieldsList"]=="{}"
            @test params.boundary_layer_fields==[7,8]
            @test params.background_field==6

            multiplied=joinpath(dir,"multiplied_field_range.geo")
            write(multiplied,"""
                Field[2] = Box;
                Field[4] = BoundaryLayer;
                Field[6] = BoundaryLayer;
                Background Field = 2 * {1:1};
                BoundaryLayer Field = 2 * {1:3};
                """)
            multiplied_params=read_geo_params(multiplied)
            @test multiplied_params.background_field==2
            @test multiplied_params.boundary_layer_fields==[2,4,6]

            wrapped=joinpath(dir,"wrapped_field_ranges.geo")
            write(wrapped,"""
                Field[1] = BoundaryLayer;
                Field[1].PointsList = -{-1:-3};
                Field[1].SizesList = -{-0.1:-0.3:-0.1};
                Field[3] = Box;
                Background Field = (1 + 2) * {1:1};
                """)
            wrapped_params=read_geo_params(wrapped)
            @test wrapped_params.fields[1].options["PointsList"]=="{1, 2, 3}"
            @test parse.(Float64,Tessella.SizeField._geo_list(
                wrapped_params.fields[1],"SizesList"))==[0.1,0.2]
            @test wrapped_params.background_field==3

            modulo=joinpath(dir,"integer_modulo_ranges.geo")
            write(modulo,"""
                Field[1] = BoundaryLayer;
                Field[1].SizesList = {(5.5 % 2.2):(5.5 % 2.2)};
                Field[11] = Box;
                Background Field = (5.5 % 2.2) * {11};
                """)
            modulo_params=read_geo_params(modulo)
            @test modulo_params.fields[1].options["SizesList"]=="{1.0}"
            @test modulo_params.background_field==11

            opaque_physical=joinpath(dir,"opaque_physical_ranges.geo")
            write(opaque_physical,"""
                Physical Surface("all", 11) = {Surface{:}};
                Physical Volume("mixed", 12) = {
                    Volume In BoundingBox{0,0,0,1,1,1}, 1:2};
                Physical Point("ternary", 13) = {1 ? 1 : 2};
                """)
            opaque_params=read_geo_params(opaque_physical)
            @test opaque_params.physical_groups[(2,11)]=="all"
            @test opaque_params.physical_groups[(3,12)]=="mixed"
            @test opaque_params.physical_groups[(0,13)]=="ternary"

            function range_error(source)
                path=joinpath(dir,"list_range_error.geo")
                write(path,source)
                try
                    read_geo_params(path)
                    return nothing
                catch err
                    return err
                end
            end
            invalid_ranges=(
                "Field[1] = Min; Field[1].FieldsList = {1:5:0};\n"=>"nonzero",
                "Field[1] = Min; Field[1].FieldsList = {:5};\n"=>"must not be empty",
                "Field[1] = Min; Field[1].FieldsList = {1:};\n"=>"must not be empty",
                "Field[1] = Min; Field[1].FieldsList = {1:2:3:4};\n"=>"at most two",
                "Field[1] = Min; Field[1].FieldsList = {missing:5};\n"=>"unknown scalar",
                "Physical Point(\"dynamic\", 1) = {points[]:5};\n"=>"unknown numeric function",
                "Field[1] = Distance; Field[1].PointsList = {1e20:1e20 + 1};\n"=>
                    "does not advance",
                "Field[1] = Box; Background Field = 0.4 * {3:4};\n"=>
                    "requires exactly one",
                "Field[1] = Box; Field[3] = Box; " *
                    "Background Field = 1 + 2 * {1:1};\n"=>"must be parenthesized",
                "Field[1] = BoundaryLayer; " *
                    "Field[1].SizesList = {1/* comment */2:12};\n"=>"unexpected token",
                "Field[1] = BoundaryLayer; BoundaryLayer Field = -{};\n"=>
                    "must not be empty",
                "Field[1] = BoundaryLayer; BoundaryLayer Field = 2 * {};\n"=>
                    "must not be empty",
                "Field[1] = Distance; Field[1].PointsList = -{};\n"=>
                    "must not be empty",
                "Field[1] = BoundaryLayer; BoundaryLayer Field = ;\n"=>
                    "must not be empty",
                "Field[1] = Distance; " *
                    "Field[1].PointsList = {2147483647:2147483648};\n"=>
                    "signed 32-bit",
                "Field[1] = Distance; Field[1].Sampling = 2147483648;\n"=>
                    "signed 32-bit")
            for (source,message) in invalid_ranges
                err=range_error(source)
                @test err isa ArgumentError
                @test occursin(message,sprint(showerror,err))
            end
            for source in (
                "Field[1] = Distance; Field[1].PointsList = {1:65537};\n",
                "Field[1] = Distance; Field[1].PointsList = {point_group[], 1:65536};\n",
                "Physical Point(\"too large\", 1) = {1:65537};\n")
                err=range_error(source)
                @test err isa ArgumentError
                @test occursin("65536 entries",sprint(showerror,err))
            end
            huge_range=
                "Field[1] = Distance; Field[1].PointsList = {1:1000000000};\n"
            range_error(huge_range)
            GC.gc()
            @test (@allocated range_error(huge_range))<=100_000

            gmsh=Sys.which("gmsh")
            @test gmsh!==nothing
            if gmsh!==nothing
                gmsh_version=strip(read(`$gmsh --version`,String))
                @test gmsh_version=="4.15.2" ||
                      startswith(gmsh_version,"4.15.2-")
                function gmsh_range_output(code)
                    output=IOBuffer()
                    command=ignorestatus(`$gmsh /dev/null -parse_and_exit -string $code`)
                    run(pipeline(command;stdout=output,stderr=output))
                    return String(take!(output))
                end
                values=gmsh_range_output(raw"""
                    a[] = {1:5}; b[] = {5:1:-2};
                    c[] = {5:-2:1}; d[] = {5:-2};
                    Printf("R %g %g %g %g %g %g %g %g %g",
                           #a[], a[0], a[4], #b[], b[0], b[2], #c[], #d[], d[7]);
                    """)
                @test occursin("R 5 1 5 3 5 1 0 8 -2",values)
                accepted=gmsh_range_output(raw"""
                    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
                    Point(3)={2,0,0,1}; Point(4)={3,0,0,1};
                    Point(5)={4,0,0,1};
                    Physical Point("ranged", 11)={1:5};
                    Field[1]=Box; Field[2]=Box; Field[3]=Box;
                    Field[4]=Box; Field[5]=Box; Field[6]=Min;
                    Field[6].FieldsList={1:5};
                    Field[7]=BoundaryLayer; Field[8]=BoundaryLayer;
                    BoundaryLayer Field={7:8};
                    """)
                @test !occursin("Error",accepted)
                wrapped_accepted=gmsh_range_output(raw"""
                    Field[1]=BoundaryLayer;
                    Field[1].PointsList=-{-1:-3};
                    Field[1].SizesList=-{-0.1:-0.3:-0.1};
                    Field[3]=Box; Background Field=(1+2)*{1:1};
                    """)
                @test !occursin("Error",wrapped_accepted)
                modulo_values=gmsh_range_output(raw"""
                    a[]={(5.5%2.2):(5.5%2.2)};
                    Printf("MOD %g %g", #a[], a[0]);
                    """)
                @test occursin("MOD 1 1",modulo_values)
                precedence_rejected=gmsh_range_output(
                    "Field[1]=Box; Field[3]=Box; " *
                    "Background Field=1+2*{1:1};")
                @test occursin("syntax error",precedence_rejected)
                cardinality_rejected=gmsh_range_output(
                    "Field[1]=Box; Background Field=0.4*{3:4};")
                @test occursin("Only 1 field",cardinality_rejected)
                comment_rejected=gmsh_range_output(
                    "Field[1]=BoundaryLayer; " *
                    "Field[1].SizesList={1/* comment */2:12};")
                @test occursin("syntax error",comment_rejected)
                for source in (
                    "Field[1]=BoundaryLayer; BoundaryLayer Field=-{};",
                    "Field[1]=BoundaryLayer; BoundaryLayer Field=2*{};",
                    "Field[1]=Distance; Field[1].PointsList=-{};",
                    "Field[1]=BoundaryLayer; BoundaryLayer Field=;")
                    @test occursin("syntax error",gmsh_range_output(source))
                end
                rejected=gmsh_range_output("a[] = {1:5:0};")
                @test occursin("Wrong increment in '1:5:0'",rejected)
            end

            allocation_path=joinpath(dir,"list_range_allocation.geo")
            function range_allocation(count)
                write(allocation_path,
                    "Field[1] = Distance; Field[1].PointsList = {1:$count};\n")
                read_geo_params(allocation_path)
                GC.gc()
                return @allocated read_geo_params(allocation_path)
            end
            allocation_2k=range_allocation(2_000)
            allocation_4k=range_allocation(4_000)
            @test allocation_2k<=1_000_000
            @test allocation_4k<=1_800_000
            @test allocation_4k<=2.15*allocation_2k+75_000
        end

        @testset "expression safety and resource limits" begin
            function geo_error(source)
                path=joinpath(dir,"expression_error.geo")
                write(path,source)
                try
                    read_geo_params(path)
                    return nothing
                catch err
                    return err
                end
            end
            invalid_expressions=(
                "missing"=>"unknown scalar identifier",
                "sin(1)"=>"unknown numeric function",
                "Rand(1)"=>"non-constant or externally stateful",
                "newp"=>"side-effecting Gmsh symbol",
                "1 < 2"=>"outside the supported arithmetic subset",
                "1 / 0"=>"non-finite",
                "1e309"=>"must be finite",
                "Hypot(1e154, 1e154)"=>"non-finite",
                "Sqrt(-1)"=>"finite real domain",
                "2++"=>"increment and decrement",
                "(1 + 2"=>"closing parenthesis")
            for (expression,message) in invalid_expressions
                err=geo_error("Mesh.MeshSizeMin = $expression;\n")
                @test err isa ArgumentError
                @test occursin(message,sprint(showerror,err))
            end

            err=geo_error("a = Rand(1); Mesh.MeshSizeMin = a;\n")
            @test err isa ArgumentError
            @test occursin("scalar variable a is unavailable",sprint(showerror,err))
            err=geo_error("a = 1; a += 1; Mesh.MeshSizeMin = a;\n")
            @test err isa ArgumentError
            @test occursin("compound assignment",sprint(showerror,err))
            err=geo_error("a = 1; a[] = {2}; Mesh.MeshSizeMin = a;\n")
            @test err isa ArgumentError
            @test occursin("list assignment",sprint(showerror,err))
            err=geo_error("a = 1; For i In {1:2}\n a += 1; EndFor\n" *
                          "junk = newp; Mesh.MeshSizeMin = a;\n")
            @test err isa ArgumentError
            @test occursin("unsupported loop",sprint(showerror,err))
            err=geo_error("For i In {1:2}\nMesh.MeshSizeMin = i; EndFor\n")
            @test err isa ArgumentError
            @test occursin("unsupported control-flow",sprint(showerror,err))
            err=geo_error("If (0)\n dummy = 1;\n Mesh.MeshSizeMin = 5;\n" *
                          "EndIf\n Mesh.MeshSizeMax = 9;\n")
            @test err isa ArgumentError
            @test occursin("unsupported control-flow",sprint(showerror,err))
            control_ignored=joinpath(dir,"control_ignored.geo")
            write(control_ignored,"If (0)\n dummy = 1;\n EndIf\n" *
                                  "Mesh.MeshSizeMax = 9;\n")
            @test read_geo_params(control_ignored).mesh_size_max==9.0
            err=geo_error("a = 1; If (0)\n dummy = 1; a = 5; EndIf\n" *
                          "Mesh.MeshSizeMin = a;\n")
            @test err isa ArgumentError
            @test occursin("scalar variable a is unavailable",sprint(showerror,err))

            @test geo_error("Field[1] = Box; Background Field = 1e100;\n") isa ArgumentError
            @test geo_error("Field[1] = Threshold; Field[1].InField = 1e100;\n") isa ArgumentError
            for source in (
                "Physical Point(\"bad\", missing) = {1};\n",
                "Physical Point(\"bad\", -1) = {1};\n",
                "Physical Point(\"bad\", 0.9) = {1};\n",
                "Physical Point(\"bad\", 2147483648) = {1};\n",
                "Physical Point(\"bad\", 1e100) = {1};\n",
                "Field[0.9] = Box;\n",
                "Field[newf] = Box;\n",
                "Field[2147483648] = Box;\n",
                "Field[1e100] = Box;\n",
                "Field[1] = Box; Field[missing].VIn = 1;\n",
                "Field[1] = Box; Background Field = 2147483648;\n",
                "Field[1] = Min; Field[1].FieldsList = {2147483648};\n")
                @test geo_error(source) isa ArgumentError
            end
            @test geo_error("Field[1 + 1] = Box; Field[2.9] = Min;\n") isa ArgumentError

            too_many_tokens=join(fill("1",div(Tessella.IO._MAX_GEO_EXPRESSION_TOKENS,2)+1),"+")
            err=geo_error("Mesh.MeshSizeMin = $too_many_tokens;\n")
            @test err isa ArgumentError
            @test occursin("tokens",sprint(showerror,err))
            err=geo_error("Physical Point(\"too many\", $too_many_tokens) = {1};\n")
            @test err isa ArgumentError
            @test occursin("tokens",sprint(showerror,err))
            too_deep=repeat("(",Tessella.IO._MAX_GEO_EXPRESSION_DEPTH+1)*"1"*
                     repeat(")",Tessella.IO._MAX_GEO_EXPRESSION_DEPTH+1)
            err=geo_error("Mesh.MeshSizeMin = $too_deep;\n")
            @test err isa ArgumentError
            @test occursin("nesting",sprint(showerror,err))
            too_long="1"*repeat(" ",Tessella.IO._MAX_GEO_EXPRESSION_BYTES)*"+0"
            err=geo_error("Mesh.MeshSizeMin = $too_long;\n")
            @test err isa ArgumentError
            @test occursin("bytes",sprint(showerror,err))
            too_many_items=join(fill("1",Tessella.IO._MAX_GEO_LIST_ITEMS+1),",")
            err=geo_error("Field[1] = Min; Field[1].FieldsList = {$too_many_items};\n")
            @test err isa ArgumentError
            @test occursin("entries",sprint(showerror,err))
        end

        @testset "comments and quoted strings cannot inject assignments" begin
            safe=joinpath(dir,"safe_quoted_expressions.geo")
            write(safe,raw"""
                base = 2 /* Mesh.MeshSizeMin = 999; */ + 1;
                Mesh.MeshSizeMin = base; // Field[999] = Box;
                Field[1] = ExternalProcess;
                Field[1].CommandLine = "printf '; Mesh.MeshSizeMax = 999; // still text'";
                Field[2] = MathEval;
                Field[2].F = "Sin(x); Background Field = 999;";
                Background Field = 2;
                """)
            params=read_geo_params(safe)
            @test params.mesh_size_min==3.0
            @test isnan(params.mesh_size_max)
            @test params.background_field==2
            @test params.fields[1].options["CommandLine"]==
                  raw"\"printf '; Mesh.MeshSizeMax = 999; // still text'\""
            @test params.fields[2].options["F"]==
                  raw"\"Sin(x); Background Field = 999;\""
        end
    end

    @testset "empty mesh round-trips" begin
        e = Mesh(Matrix{Float64}(undef, 3, 0))
        for ver in (2.2, 4.1)
            p = joinpath(dir, "empty_$(ver).msh")
            write_msh(p, e; version=ver)
            f = read_msh(p)
            @test nnodes(f.mesh) == 0 && ntets(f.mesh) == 0
        end
    end
end
