using Test
using Tessella
using SHA
using Tessella.MeshTypes: ntris, nnodes, nsegs, ntets, validate
using Tessella.Elements: ElementBlock, MixedMesh
using Tessella.Model: model_to_mixed, model_physical_groups,
                      model_entities_for_physical_group,
                      model_entities_for_physical_name,
                      model_physical_groups_for_entity

const _GEO_SQUARE = raw"""
Point(1) = {0,0,0,0.2};
Point(2) = {1,0,0,0.2};
Point(3) = {1,1,0,0.2};
Point(4) = {0,1,0,0.2};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
"""

const _GEO_BOX = raw"""
Point(1) = {0,0,0,0.4};
Point(2) = {1,0,0,0.4};
Point(3) = {1,1,0,0.4};
Point(4) = {0,1,0,0.4};
Point(5) = {0,0,1,0.4};
Point(6) = {1,0,1,0.4};
Point(7) = {1,1,1,0.4};
Point(8) = {0,1,1,0.4};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
Line(9) = {1,5}; Line(10) = {2,6}; Line(11) = {3,7}; Line(12) = {4,8};
Curve Loop(1) = {1,2,3,4}; Curve Loop(2) = {5,6,7,8};
Curve Loop(3) = {1,10,-5,-9}; Curve Loop(4) = {2,11,-6,-10};
Curve Loop(5) = {3,12,-7,-11}; Curve Loop(6) = {4,9,-8,-12};
Plane Surface(1) = {1}; Plane Surface(2) = {2};
Plane Surface(3) = {3}; Plane Surface(4) = {4};
Plane Surface(5) = {5}; Plane Surface(6) = {6};
Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};
"""

# Two coplanar squares sharing edge Line(2) — surface 2 reuses curve 2 with a
# negative sign so curve 2 is owned by both surfaces.
const _GEO_TWO_SQUARES = raw"""
Point(1) = {0,0,0,0.2};
Point(2) = {1,0,0,0.2};
Point(3) = {1,1,0,0.2};
Point(4) = {0,1,0,0.2};
Point(5) = {2,0,0,0.2};
Point(6) = {2,1,0,0.2};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {6,3};
Line(7) = {2,5};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Curve Loop(2) = {7,5,6,-2};
Plane Surface(2) = {2};
"""

function _execute_constraint_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _constraint_error(source::AbstractString;mesh_dim=0)
    try
        _execute_constraint_source(source;mesh_dim=mesh_dim)
        return nothing
    catch err
        return err
    end
end

@testset ".geo Transfinite Curve" begin
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Curve{1,3} = 6 Using Progression 2;
        Transfinite Curve{-2} = 7 Using "Bump" 1.5;
        Transfinite Curve{4,99} = 8;
        """)
    meshing=execution.model.meshing
    @test meshing.transfinite_curves[1]==
          (num_nodes=6,kind=:progression,coef=2.0,reversed=false)
    @test meshing.transfinite_curves[3]==
          (num_nodes=6,kind=:progression,coef=2.0,reversed=false)
    # A negative list tag negates the stored type — `reversed` on the record.
    @test meshing.transfinite_curves[2]==
          (num_nodes=7,kind=:bump,coef=1.5,reversed=true)
    # Missing tags skip silently (Gmsh `FindCurve` miss).
    @test meshing.transfinite_curves[4]==
          (num_nodes=8,kind=:progression,coef=1.0,reversed=false)
    @test !haskey(meshing.transfinite_curves,99)

    # `{:}` wildcard applies to every curve that exists at execution time.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Curve{:} = 9;
        Line(5) = {1,3};
        """)
    meshing=execution.model.meshing
    for curve in 1:4
        @test meshing.transfinite_curves[curve].num_nodes==9
    end
    @test !haskey(meshing.transfinite_curves,5)

    # An entry truncating to 0 is the `setTransfiniteLine(0, ...)` wildcard;
    # its sign selects the stored type orientation.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Curve{1, -0.5} = 4;
        """)
    meshing=execution.model.meshing
    for curve in 1:4
        @test meshing.transfinite_curves[curve].reversed
    end
    @test meshing.transfinite_curves[1].num_nodes==4

    # An exact-zero entry stores `type * gmsh_sign(0)` = type 0 — Gmsh's
    # `F_Transfinite` hits the unknown-type fallback and emits a uniform
    # distribution rather than the requested law.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Curve{0} = 6 Using Progression 2;
        """)
    meshing=execution.model.meshing
    for curve in 1:4
        @test meshing.transfinite_curves[curve]==
              (num_nodes=6,kind=:uniform,coef=2.0,reversed=false)
    end

    # `(int)` truncation on tags: 1.9 targets curve 1.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Curve{1.9} = 5;
        """)
    @test haskey(execution.model.meshing.transfinite_curves,1)
    @test !haskey(execution.model.meshing.transfinite_curves,2)

    # `n < 2` clamps to two nodes like Gmsh's `nbPointsTransfinite`.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Curve{1} = -3;
        """)
    @test execution.model.meshing.transfinite_curves[1].num_nodes==2

    # Every `.geo` `Using` law spelling maps to its record kind; the
    # grammar-only Beta_Symmetrical kinds store distinct kinds that mesh as
    # the unknown-type uniform fallback.
    for (law,kind) in ("Progression"=>:progression,"Power"=>:progression,
                       "Bump"=>:bump,"Beta"=>:beta,
                       "Progression_HWall"=>:progression_hwall,
                       "Bump_HWall"=>:bump_hwall,"Beta_HWall"=>:beta_hwall,
                       "Beta_Symmetrical"=>:beta_symmetrical,
                       "Beta_Symmetrical_HWall"=>:beta_symmetrical_hwall)
        execution=_execute_constraint_source(_GEO_SQUARE * """
            Transfinite Curve{1} = 6 Using $law 0.05;
            """)
        @test execution.model.meshing.transfinite_curves[1].kind==kind
    end

    err=_constraint_error(_GEO_SQUARE * raw"""
        Transfinite Curve{1} = 6 Using Nonsense 1;
        """)
    @test err isa ArgumentError
    @test occursin("unknown transfinite mesh type",string(err))
end

@testset ".geo Transfinite Surface/Volume + TransfQuadTri" begin
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Surface{1} = {1,2,3,4} Right;
        """)
    @test execution.model.meshing.transfinite_surfaces[1]==
          (arrangement=:right,corners=[1,2,3,4])

    # Arrangement words: Left/Right/AlternateLeft/AlternateRight, bare
    # `Alternate` and unknown words -> AlternateRight, none -> Left.
    for (word,arrangement) in ("Left"=>:left,"Right"=>:right,
                               "AlternateLeft"=>:alternate_left,
                               "AlternateRight"=>:alternate_right,
                               "Alternate"=>:alternate_right,
                               "AnythingElse"=>:alternate_right)
        execution=_execute_constraint_source(_GEO_SQUARE * """
            Transfinite Surface{1} $word;
            """)
        @test execution.model.meshing.transfinite_surfaces[1].arrangement==
              arrangement
    end
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Surface{1};
        """)
    @test execution.model.meshing.transfinite_surfaces[1].arrangement==:left

    # Negative tags skip silently; a zero entry wildcards with empty corners.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Plane Surface(2) = {1};
        Transfinite Surface{-1,2};
        Delete{Surface{2};}
        """)
    @test !haskey(execution.model.meshing.transfinite_surfaces,1)
    @test !haskey(execution.model.meshing.transfinite_surfaces,2)

    # A corner list with a bad count errors on a live surface (Gmsh's
    # internals setter reports it), while a missing surface skips silently.
    err=_constraint_error(_GEO_SQUARE * raw"""
        Transfinite Surface{1} = {1,2};
        """)
    @test err isa ArgumentError
    @test occursin("3 or 4 corner",string(err))
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Surface{99} = {1,2};
        """)
    @test !haskey(execution.model.meshing.transfinite_surfaces,99)

    # An unknown corner point on a live surface is a Gmsh error.
    err=_constraint_error(_GEO_SQUARE * raw"""
        Transfinite Surface{1} = {1,2,3,99};
        """)
    @test err isa ArgumentError
    @test occursin("unknown corner Point[99]",string(err))

    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Transfinite Volume{1} = {1,2,3,4,5,6,7,8};
        TransfQuadTri{1};
        """)
    @test execution.model.meshing.transfinite_volumes[1]==collect(1:8)
    @test 1 in execution.model.meshing.quad_tri

    # 5-corner (and other non-6/8) lists are silently dropped.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Transfinite Volume{1} = {1,2,3,4,5};
        """)
    @test execution.model.meshing.transfinite_volumes[1]==Int[]

    # TransfQuadTri wildcard flags every current volume; the mesher consumes
    # the flag when a transfinite volume is generated.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        TransfQuadTri{:};
        """)
    @test 1 in execution.model.meshing.quad_tri
end

@testset ".geo transfinite volume meshing" begin
    using Tessella.MeshTypes: nnodes, ntets, tet_signed_volume

    constraints = raw"""
        Transfinite Curve{:} = 4;
        Transfinite Surface{:};
        Transfinite Volume{1} = {1,2,3,4,5,6,7,8};
        """
    execution=_execute_constraint_source(_GEO_BOX * constraints;mesh_dim=3)
    @test validate(execution.mesh).ok
    @test (nnodes(execution.mesh),ntets(execution.mesh))==(64,162)

    # Every boundary surface must be transfinite (Gmsh's "Incompatible
    # surface" blocker): dropping one declaration leaves the volume unmeshed.
    err=_constraint_error(_GEO_BOX * raw"""
        Transfinite Curve{:} = 4;
        Transfinite Surface{1,2,3,4,5};
        Transfinite Volume{1} = {1,2,3,4,5,6,7,8};
        """;mesh_dim=3)
    @test err isa ArgumentError
    @test occursin("Incompatible surface 6",string(err))

    # Every boundary curve must be transfinite.
    err=_constraint_error(_GEO_BOX * raw"""
        Transfinite Curve{1,2,3,4,5,6,7,8,9,10,11} = 4;
        Transfinite Surface{:};
        Transfinite Volume{1} = {1,2,3,4,5,6,7,8};
        """;mesh_dim=3)
    @test err isa ArgumentError
    @test occursin("requires boundary Curve[12] to be transfinite",string(err))

    # Opposite edges in one direction must carry equal node counts; the
    # boundary surfaces mesh before the volume, so the patch kernel's
    # side-mismatch diagnostic fires first (the `Mesh 3` statement path
    # reports the same error).
    err=_constraint_error(_GEO_BOX * raw"""
        Transfinite Curve{1} = 5;
        Transfinite Curve{2,3,4,5,6,7,8,9,10,11,12} = 4;
        Transfinite Surface{:};
        Transfinite Volume{1} = {1,2,3,4,5,6,7,8};
        """;mesh_dim=3)
    @test err isa ArgumentError
    @test occursin("non-matching node counts",string(err))

    # A genuinely non-affine block: Point(7) moves off the corner
    # parallelepiped, so its two incident ruled faces mesh as warped
    # transfinite patches and the interior follows the face Coons
    # interpolation instead of the trilinear corner map.
    warped_box=raw"""
        Point(1) = {0,0,0,0.5}; Point(2) = {1,0,0,0.5};
        Point(3) = {1,1,0,0.5}; Point(4) = {0,1,0,0.5};
        Point(5) = {0,0,1,0.5}; Point(6) = {1,0,1,0.5};
        Point(7) = {1.2,1,1,0.5}; Point(8) = {0,1,1,0.5};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
        Line(9) = {1,5}; Line(10) = {2,6}; Line(11) = {3,7}; Line(12) = {4,8};
        Curve Loop(1) = {1,2,3,4}; Curve Loop(2) = {5,6,7,8};
        Curve Loop(3) = {1,10,-5,-9}; Curve Loop(4) = {2,11,-6,-10};
        Curve Loop(5) = {3,12,-7,-11}; Curve Loop(6) = {4,9,-8,-12};
        Plane Surface(1) = {1}; Plane Surface(2) = {2};
        Plane Surface(3) = {3}; Surface(4) = {4};
        Plane Surface(5) = {5}; Plane Surface(6) = {6};
        Surface Loop(1) = {1,2,3,4,5,6};
        Volume(1) = {1};
        """ * constraints
    execution=_execute_constraint_source(warped_box;mesh_dim=3)
    @test validate(execution.mesh).ok
    @test (nnodes(execution.mesh),ntets(execution.mesh))==(64,162)
    # The shifted vertex and its subdivided incident edges are in the mesh —
    # x>1 nodes cannot come from the affine corner parallelepiped.
    @test maximum(execution.mesh.coords[1,:])==1.2
    @test count(>(1.0),execution.mesh.coords[1,:])>4
    # A v-direction edge of the warped w-max face runs (1,0,1)→(1.2,1,1):
    # its transfinite subdivision lands at x=1+0.2t.
    @test any(isapprox(1.0+0.2/3;atol=1e-12),
              execution.mesh.coords[1,:])

    # A curved boundary edge can leave its `Plane Surface`'s declared plane:
    # upstream's `planeSurface` mean plane comes from on-curve boundary
    # samples — the arc's control point never enters it — so the front
    # face's interior interpolates exactly on y=0 while the arc keeps its
    # true positions. The arc bulging OUT of the box stays a valid
    # transfinite volume, and so does a shallow inward bulge.
    curved_volume(center)=raw"""
        Point(1) = {0,0,0,0.4}; Point(2) = {1,0,0,0.4};
        Point(3) = {1,1,0,0.4}; Point(4) = {0,1,0,0.4};
        Point(5) = {0,0,1,0.4}; Point(6) = {1,0,1,0.4};
        Point(7) = {1,1,1,0.4}; Point(8) = {0,1,1,0.4};
        Point(9) = """ * center * raw""";
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Line(5) = {2,6}; Circle(6) = {6,9,5}; Line(7) = {5,1};
        Line(8) = {3,7}; Line(9) = {7,8}; Line(10) = {8,4};
        Line(11) = {6,7}; Line(12) = {8,5};
        Curve Loop(1) = {1,2,3,4}; Curve Loop(2) = {1,5,6,7};
        Curve Loop(3) = {8,9,10,-3}; Curve Loop(4) = {-4,-10,12,7};
        Curve Loop(5) = {2,8,-11,-5}; Curve Loop(6) = {-6,11,9,12};
        Plane Surface(1) = {1}; Plane Surface(2) = {2};
        Plane Surface(3) = {3}; Plane Surface(4) = {4};
        Plane Surface(5) = {5}; Surface(6) = {6};
        Surface Loop(1) = {1,2,3,4,5,6};
        Volume(1) = {1};
        Transfinite Curve{:} = 5;
        Transfinite Surface{:};
        Transfinite Volume{1};
        """
    for center in ("{0.5,0.3,1,0.4}","{0.5,-0.25,1,0.4}")
        execution=_execute_constraint_source(
            curved_volume(center);mesh_dim=3)
        @test validate(execution.mesh).ok
        @test ntets(execution.mesh)==384
        # Every tetrahedron is strictly positive — the volume did not fold.
        @test all(t->tet_signed_volume(
                      execution.mesh.coords[:,execution.mesh.tets[1,t]],
                      execution.mesh.coords[:,execution.mesh.tets[2,t]],
                      execution.mesh.coords[:,execution.mesh.tets[3,t]],
                      execution.mesh.coords[:,execution.mesh.tets[4,t]])>0,
                  axes(execution.mesh.tets,2))
        # The arc's interior nodes keep their true off-plane coordinates.
        @test count(i->abs(execution.mesh.coords[2,i])>1e-9,
                    axes(execution.mesh.coords,2))>=3
    end

    # The strongly inward-bulging mirror (an antipodal arc dipping deep into
    # the box) is REJECTED by the boundary-fold audit: Gmsh 4.15.2 silently
    # emits a self-intersecting mesh for it — an interior tet edge pierces
    # the emitted boundary band — and Tessella refuses to emit the defect.
    err=_constraint_error(curved_volume("{0.5,0,1,0.4}");mesh_dim=3)
    @test err isa ArgumentError
    @test occursin("folds",string(err))

    # A five-face transfinite volume: Gmsh's legacy `Mesh.TransfiniteTri=0`
    # path meshes the prism as a degenerate hexahedron (s3≡s0, s7≡s4), with
    # collapsed-triangle grids on the two triangular faces. Gmsh 4.15.2 emits
    # 52 nodes / 246 dim≥1 elements (27 segments, 84 triangles, 135
    # tetrahedra) for this fixture.
    prism=raw"""
        Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={0,1,0};
        Point(4)={0,0,1}; Point(5)={1,0,1}; Point(6)={0,1,1};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1};
        Line(4)={1,4}; Line(5)={2,5}; Line(6)={3,6};
        Line(7)={4,5}; Line(8)={5,6}; Line(9)={6,4};
        Curve Loop(1)={1,2,3}; Curve Loop(2)={7,8,9};
        Curve Loop(3)={1,5,-7,-4}; Curve Loop(4)={2,6,-8,-5};
        Curve Loop(5)={3,4,-9,-6};
        Surface(1)={1}; Surface(2)={2};
        Surface(3)={3}; Surface(4)={4}; Surface(5)={5};
        Surface Loop(1)={1,2,3,4,5};
        Volume(1)={1};
        Transfinite Curve{:}=4;
        Transfinite Surface{:};
        """
    execution=_execute_constraint_source(prism*raw"""
        Transfinite Volume{1}={1,2,3,4,5,6};
        Mesh 3;
        """)
    using Tessella.MeshTypes: nsegs
    @test validate(execution.mesh).ok
    @test (nnodes(execution.mesh),nsegs(execution.mesh),
           ntris(execution.mesh),ntets(execution.mesh))==(52,27,84,135)
    # The volume-only part carries no boundary triangles.
    part=geo_entity_mesh(execution,3,1)
    @test (nnodes(part),ntris(part),ntets(part))==(52,0,135)
    # Six corner points are the tetrahedral degree-6 vertices of the unit
    # prism and every tet is positively oriented.
    @test all(t->tet_signed_volume(
                  execution.mesh.coords[:,execution.mesh.tets[1,t]],
                  execution.mesh.coords[:,execution.mesh.tets[2,t]],
                  execution.mesh.coords[:,execution.mesh.tets[3,t]],
                  execution.mesh.coords[:,execution.mesh.tets[4,t]])>0,
              axes(execution.mesh.tets,2))

    # Corner detection without an explicit list seeds the apex from a
    # triangular face — the same five-face mesh results.
    execution=_execute_constraint_source(prism*raw"""
        Transfinite Volume{1};
        Mesh 3;
        """)
    @test validate(execution.mesh).ok
    @test (nnodes(execution.mesh),ntets(execution.mesh))==(52,135)

    # The compact `Mesh.TransfiniteTri=1` algorithm keeps the triangular
    # faces on equal-sided compact lattices (Gmsh's `transfinite3` branch):
    # the volume's slot grid still reads the face lattice with the j>i
    # slots aliased to the diagonal, and the cell subdivision emits
    # SIM_10–SIM_12 on diagonal cells and SIM_7–SIM_12 on strictly-lower
    # cells. Gmsh 4.15.2 emits 46 nodes / 189 dim≥1 elements (27 segments,
    # 72 triangles, 81 tetrahedra) for this fixture — the two interior tab
    # slots behind the compact diagonal are unreferenced orphan nodes with
    # their own evaluated coordinates, exactly like upstream.
    execution=_execute_constraint_source(prism*raw"""
        Mesh.TransfiniteTri=1;
        Transfinite Volume{1}={1,2,3,4,5,6};
        Mesh 3;
        """)
    @test validate(execution.mesh).ok
    @test (nnodes(execution.mesh),nsegs(execution.mesh),
           ntris(execution.mesh),ntets(execution.mesh))==(46,27,72,81)
    # The raw volume part carries all 46 slot columns — the two
    # behind-diagonal orphans hold distinct evaluated coordinates like
    # Gmsh's written node entries.
    part=geo_entity_mesh(execution,3,1)
    @test (nnodes(part),ntris(part),ntets(part))==(46,0,81)
    @test all(t->tet_signed_volume(
                  execution.mesh.coords[:,execution.mesh.tets[1,t]],
                  execution.mesh.coords[:,execution.mesh.tets[2,t]],
                  execution.mesh.coords[:,execution.mesh.tets[3,t]],
                  execution.mesh.coords[:,execution.mesh.tets[4,t]])>0,
              axes(execution.mesh.tets,2))
    # The triangular face mesh and the volume's canonical grid come from
    # the same canonicalized boundary chain — every node welds bitwise.
    surface_part=geo_entity_mesh(execution,2,2)
    volume_nodes=Set(Tuple(part.coords[:,c]) for c in axes(part.coords,2))
    @test all(Tuple(surface_part.coords[:,c]) in volume_nodes
              for c in axes(surface_part.coords,2))

    # Sign-reversed loops yield the same canonicalized boundary chain — the
    # mesh is identical to the forward declaration.
    execution=_execute_constraint_source(replace(
        prism,"Curve Loop(2)={7,8,9}"=>"Curve Loop(2)={-7,-9,-8}")*raw"""
        Mesh.TransfiniteTri=1;
        Transfinite Volume{1}={1,2,3,4,5,6};
        Mesh 3;
        """)
    @test validate(execution.mesh).ok
    @test (nnodes(execution.mesh),ntets(execution.mesh))==(46,81)

    # A triangular face whose canonical corner order does not begin on the
    # prism apex cannot occupy the degenerate slot — Gmsh reports the same
    # "Incompatible surface" rejection.
    err=_constraint_error(replace(
        prism,"Curve Loop(2)={7,8,9}"=>"Curve Loop(2)={9,7,8}")*raw"""
        Mesh.TransfiniteTri=1;
        Transfinite Volume{1}={1,2,3,4,5,6};
        Mesh 3;
        """)
    @test err isa ArgumentError
    @test occursin("Incompatible surface 2",string(err))

    # Compact triangles require equal side divisions — a mismatched edge
    # count is rejected (the quad whose opposite sides desynchronize fires
    # first, the same "non-matching" audit upstream performs).
    err=_constraint_error(prism*raw"""
        Mesh.TransfiniteTri=1;
        Transfinite Curve{7,8,9}=3;
        Transfinite Volume{1}={1,2,3,4,5,6};
        Mesh 3;
        """)
    @test err isa ArgumentError
    @test occursin("non-matching node counts",string(err))

    # A five-face volume whose boundary is not a prism (a square pyramid:
    # one quadrilateral + four triangular faces) is rejected by the
    # boundary-topology audit rather than silently hexahedral.
    err=_constraint_error(raw"""
        Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={1,1,0};
        Point(4)={0,1,0}; Point(5)={0.5,0.5,1};
        Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
        Line(5)={1,5}; Line(6)={2,5}; Line(7)={3,5}; Line(8)={4,5};
        Curve Loop(1)={1,2,3,4}; Curve Loop(2)={1,6,-5};
        Curve Loop(3)={2,7,-6}; Curve Loop(4)={3,8,-7};
        Curve Loop(5)={4,5,-8};
        Surface(1)={1}; Surface(2)={2}; Surface(3)={3};
        Surface(4)={4}; Surface(5)={5};
        Surface Loop(1)={1,2,3,4,5};
        Volume(1)={1};
        Transfinite Curve{:}=4;
        Transfinite Surface{:};
        Transfinite Volume{1};
        Mesh 3;
        """)
    @test err isa ArgumentError
    @test occursin("prismatic boundary topology",string(err))
end

@testset ".geo transfinite curves on curved edges" begin
    using Tessella.MeshTypes: nnodes, ntets, ntris

    # A quarter-circle boundary at 5 nodes subdivides by angle fraction —
    # the nodes land on the analytic 22.5° marks of the arc.
    circle_patch=raw"""
        Point(1) = {1,0,0,1}; Point(2) = {0,1,0,1}; Point(9) = {0,0,0,1};
        Point(3) = {0.2,1.7,0,1}; Point(4) = {1,1.5,0,1};
        Circle(1) = {1,9,2};
        Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Curve Loop(1) = {1,2,3,4};
        Plane Surface(1) = {1};
        Transfinite Curve{:} = 4;
        Transfinite Curve{1,3} = 5;
        Transfinite Surface{1};
        """
    execution=_execute_constraint_source(circle_patch*"Mesh 2;\n";mesh_dim=0)
    @test validate(execution.mesh).ok
    @test (nnodes(execution.mesh),ntris(execution.mesh))==(21,24)
    for degrees in (22.5,45.0,67.5)
        expected=[cosd(degrees),sind(degrees),0.0]
        @test any(i->isapprox(execution.mesh.coords[:,i],expected;
                              atol=1e-9),
                  axes(execution.mesh.coords,2))
    end

    # A spline boundary subdivides its normalized parameter by arc length
    # for the uniform law — matching Gmsh's `F_Transfinite` default arm.
    spline_head=raw"""
        Point(1) = {0,0,0,1}; Point(2) = {1,0.1,0,1};
        Point(5) = {0.3,0.45,0,1}; Point(6) = {0.7,0.38,0,1};
        Point(3) = {1,1.5,0,1}; Point(4) = {0,1.4,0,1};
        Spline(1) = {1,5,6,2};
        Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
        Curve Loop(1) = {1,2,3,4};
        Plane Surface(1) = {1};
        Transfinite Curve{:} = 4;
        Transfinite Curve{1,3} = 5;
        Transfinite Surface{1};
        """
    execution=_execute_constraint_source(spline_head*"Mesh 2;\n";mesh_dim=0)
    @test validate(execution.mesh).ok
    @test (nnodes(execution.mesh),ntris(execution.mesh))==(22,24)
    # Nonuniform laws are DENSITIES over the native parameter upstream — the
    # node parameters are NOT the law positions arc-length-inverted. The
    # values below are the Gmsh 4.15.2 stored parameters for this case
    # (observed via its mesh node parameters); a length-fraction mapping
    # would give ~0.1108/0.2722/0.5719 instead.
    execution=_execute_constraint_source(spline_head*raw"""
        Transfinite Curve{1} = 5 Using Progression 1.4;
        Mesh 2;
        """;mesh_dim=0)
    @test validate(execution.mesh).ok
    @test (nnodes(execution.mesh),ntris(execution.mesh))==(22,24)
    params=execution.model.curve_params[1]
    @test isapprox(params,[0.0,0.1163858,0.2770245,0.5685709,1.0];atol=2e-6)
    # The reversed declaration mirrors the distribution.
    execution=_execute_constraint_source(spline_head*raw"""
        Transfinite Curve{-1} = 5 Using Progression 1.4;
        Mesh 2;
        """;mesh_dim=0)
    params=execution.model.curve_params[1]
    @test isapprox(params,[0.0,0.3084843,0.6201999,0.8477194,1.0];atol=2e-6)
    # `Beta_Symmetrical` has no `F_Transfinite` case upstream — the curve
    # falls to the unknown-type `val = 1` arm, which on a curved
    # parametrization is parameter-uniform, not length-uniform.
    execution=_execute_constraint_source(spline_head*raw"""
        Transfinite Curve{1} = 5 Using Beta_Symmetrical 2.5;
        Mesh 2;
        """;mesh_dim=0)
    params=execution.model.curve_params[1]
    @test params==[0.0,0.25,0.5,0.75,1.0]
    # The transfinite side chain and the stored curve discretization share
    # the same native parameters bitwise — sibling parts cannot duplicate
    # boundary nodes.
    spec=execution.model.meshing.transfinite_curves[1]
    t0,t1=Tessella.Model._model_curve_param_bounds(
        execution.model,1,"test")
    @test Tessella.Model._model_curve_transfinite_native_params(
        execution.model,1,t0,t1,spec,"test")==params
    # An annular-sector volume: two ruled cylindrical faces plus four planar
    # faces, all transfinite. The circle boundaries feed the volume's
    # canonical face grids and interior — Gmsh emits the identical counts
    # (the two extra nodes are the circle-center point entities).
    annulus=raw"""
        Point(1) = {0.5,0,0,1}; Point(2) = {1,0,0,1};
        Point(3) = {0,0.5,0,1}; Point(4) = {0,1,0,1};
        Point(5) = {0.5,0,1,1}; Point(6) = {1,0,1,1};
        Point(7) = {0,0.5,1,1}; Point(8) = {0,1,1,1};
        Point(9) = {0,0,0,1}; Point(10) = {0,0,1,1};
        Circle(1) = {1,9,3}; Circle(2) = {2,9,4};
        Circle(3) = {5,10,7}; Circle(4) = {6,10,8};
        Line(5) = {1,2}; Line(6) = {3,4};
        Line(7) = {5,6}; Line(8) = {7,8};
        Line(9) = {1,5}; Line(10) = {2,6};
        Line(11) = {3,7}; Line(12) = {4,8};
        Curve Loop(1) = {1,6,-2,-5}; Curve Loop(2) = {3,8,-4,-7};
        Curve Loop(3) = {1,11,-3,-9}; Curve Loop(4) = {2,12,-4,-10};
        Curve Loop(5) = {5,10,-7,-9}; Curve Loop(6) = {6,12,-8,-11};
        Plane Surface(1) = {1}; Plane Surface(2) = {2};
        Surface(3) = {3}; Surface(4) = {4};
        Plane Surface(5) = {5}; Plane Surface(6) = {6};
        Surface Loop(1) = {1,2,3,4,5,6};
        Volume(1) = {1};
        Transfinite Curve{:} = 4;
        Transfinite Surface{:};
        Transfinite Volume{1};
        """
    execution=_execute_constraint_source(annulus*"Mesh 3;\n";mesh_dim=0)
    @test validate(execution.mesh).ok
    @test (nnodes(execution.mesh),ntets(execution.mesh),
           ntris(execution.mesh))==(66,162,108)
end

@testset ".geo Recombine/Smoother/Algorithm/SizeFromBoundary" begin
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Recombine Surface{1} = 30;
        Smoother Surface{1} = 4;
        MeshAlgorithm Surface{1} = 6;
        MeshSizeFromBoundary Surface{1} = 1;
        """)
    meshing=execution.model.meshing
    @test meshing.recombine[(2,1)]==30.0
    @test meshing.smoothing[(2,1)]==4
    @test meshing.algorithm[(2,1)]==6
    @test meshing.size_from_boundary[(2,1)]==1

    # Default recombine angle is 45; `(int)` truncation applies to the angle.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Recombine Surface{1};
        """)
    @test execution.model.meshing.recombine[(2,1)]==45.0
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Recombine Surface{1} = 44.7;
        """)
    @test execution.model.meshing.recombine[(2,1)]==44.0

    # `Recombine Volume` carries no angle.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Recombine Volume{1};
        """)
    @test execution.model.meshing.recombine[(3,1)]==0.0

    # Wildcard forms expand over entities alive at execution time.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Plane Surface(2) = {1};
        Recombine Surface{:} = 20;
        Smoother Surface{:} = 2;
        """)
    meshing=execution.model.meshing
    @test meshing.recombine[(2,1)]==20.0
    @test meshing.recombine[(2,2)]==20.0
    @test meshing.smoothing[(2,1)]==2
    @test meshing.smoothing[(2,2)]==2

    # `MeshSizeFromBoundary = 0` stores the verbatim 0 — upstream
    # `getMeshSizeFromBoundary` resolves it to "do not extend", so the record
    # carries the disable, not a deletion back to the enabled default.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        MeshSizeFromBoundary Surface{1} = 1;
        MeshSizeFromBoundary Surface{1} = 0;
        """)
    @test execution.model.meshing.size_from_boundary[(2,1)]==0

    # MeshAlgorithm/MeshSizeFromBoundary only exist for surfaces in the .geo
    # grammar; other dimensions are syntax errors upstream.
    err=_constraint_error(_GEO_SQUARE * raw"""
        MeshAlgorithm Curve{1} = 5;
        """)
    @test err isa ArgumentError
    err=_constraint_error(_GEO_BOX * raw"""
        MeshSizeFromBoundary Volume{1} = 1;
        """)
    @test err isa ArgumentError

    # Missing surfaces are skipped silently (signed FindSurface miss).
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        MeshAlgorithm Surface{99,-1} = 5;
        MeshSizeFromBoundary Surface{99} = 1;
        """)
    @test isempty(execution.model.meshing.algorithm)
    @test isempty(execution.model.meshing.size_from_boundary)
end

@testset ".geo Reverse/Relocate/Reorient/Degenerated/Compound" begin
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Reverse Curve{1,2};
        ReverseMesh Surface{1};
        """)
    meshing=execution.model.meshing
    @test meshing.reverse[(1,1)] && meshing.reverse[(1,2)]
    @test meshing.reverse[(2,1)]

    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Degenerated Curve{3,99,-4};
        """)
    @test execution.model.meshing.degenerated==Set([3])

    # RelocateMesh/ReorientMesh parse and validate their lists; neither has
    # any effect without an existing mesh.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        RelocateMesh Point{1,2};
        RelocateMesh Curve{:};
        ReorientMesh Volume{1};
        """)
    @test isempty(execution.model.meshing.outward_orientation)

    # Compound records the raw member list per dimension; `MeshAlgorithm n`
    # appends `-(int)n`.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Compound Curve{1,2,-3};
        Compound Surface{1} MeshAlgorithm 5;
        """)
    @test execution.model.meshing.compounds==
          [1=>[1,2,-3],2=>[1,-5]]

    err=_constraint_error(_GEO_SQUARE * raw"""
        ReorientMesh Volume{99};
        """)
    @test err===nothing
end

@testset ".geo OptimizeMesh" begin
    # Without a mesh every method name validates but nothing moves, like
    # `GModel::optimizeMesh` iterating an entity list with no cells.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        OptimizeMesh "Gmsh";
        OptimizeMesh "Laplace2D";
        """)
    @test execution.mesh===nothing

    # Upstream validates the method name before any entity iteration, so an
    # unknown optimizer errors even on an empty model.
    err=_constraint_error(raw"""
        OptimizeMesh "Bogus";
        """)
    @test err isa ArgumentError
    @test occursin("optimization method", sprint(showerror,err))

    # A 2-D cache has no regions: the 3-D methods are silent no-ops while the
    # 2-D methods run Laplacian smoothing and keep the mesh.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        OptimizeMesh "Gmsh";
        OptimizeMesh "Optimize";
        OptimizeMesh "Relocate3D";
        OptimizeMesh "Laplace2D";
        OptimizeMesh "Relocate2D";
        """; mesh_dim=2)
    @test execution.mesh!==nothing
    @test ntris(execution.mesh)>0
    @test validate(execution.mesh).ok

    err=_constraint_error(_GEO_SQUARE * raw"""
        OptimizeMesh "Netgen";
        """; mesh_dim=2)
    @test err isa ArgumentError
    @test occursin("optimization method", sprint(showerror,err))
end

@testset ".geo Delete per-entity semantics" begin
    # Ownership refusal: a curve owned by a surviving surface cannot be
    # deleted; once the surface goes, it can.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Delete{Curve{1};}
        """)
    @test haskey(execution.model.curves,1)
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Delete{Surface{1};}
        Delete{Curve{1};}
        """)
    @test !haskey(execution.model.surfaces,1)
    @test !haskey(execution.model.curves,1)

    # Point deletion is refused while a surviving curve references it (as an
    # endpoint or a control point).
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Delete{Point{1};}
        """)
    @test haskey(execution.model.points,1)

    # Negative surface/volume tags are silent no-ops (signed lookup), while
    # point and curve deletion tries both signs.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Delete{Surface{-1};}
        """)
    @test haskey(execution.model.surfaces,1)
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Delete{Surface{1};}
        Delete{Curve{-1};}
        """)
    @test !haskey(execution.model.curves,1)

    # Recursive Delete removes the pre-collected boundary closure flat;
    # entities still owned by survivors are refused.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Recursive Delete{Volume{1};}
        """)
    model=execution.model
    @test isempty(model.volumes)
    @test isempty(model.surfaces)
    @test isempty(model.curves)
    @test isempty(model.points)
    # The dangling loop records survive, like Gmsh's internals.
    @test haskey(model.surface_loops,1)
    @test haskey(model.loops,1)

    # Recursive Delete does not cross ownership: a surface owned by a
    # surviving volume is refused outright.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Recursive Delete{Surface{1};}
        """)
    @test haskey(execution.model.surfaces,1)

    # Boundary entities shared with surviving surfaces are refused even under
    # Recursive Delete — deleting square 1 leaves shared curve 2 (owned by
    # square 2) plus every entity still referenced by a survivor.
    execution=_execute_constraint_source(_GEO_TWO_SQUARES * raw"""
        Recursive Delete{Surface{1};}
        """)
    model=execution.model
    @test !haskey(model.surfaces,1)
    @test haskey(model.surfaces,2)
    @test sort!(collect(keys(model.curves)))==[2,5,6,7]
    @test sort!(collect(keys(model.points)))==[2,3,5,6]
    # The dangling loop record survives, like Gmsh's internals.
    @test haskey(model.loops,1)

    # Deleting a point owned by removed curves in the same statement works
    # when the curve deletions come first (list order); a point still owned
    # by a surviving curve is refused.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Delete{Surface{1};}
        Delete{Curve{3,4}; Point{4};}
        """)
    @test !haskey(execution.model.points,4)
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Delete{Surface{1};}
        Delete{Curve{4}; Point{4};}
        """)
    @test haskey(execution.model.points,4)

    # Deleting an embedded source drops the embedding; deleting a target does
    # not drop unrelated embeds.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Point(9) = {0.5,0.5,0,0.1};
        Point{9} In Surface{1};
        Delete{Point{9};}
        """)
    @test isempty(execution.model.embeds)

    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Delete Embedded{Surface{1};}
        """)
    @test isempty(execution.model.embeds)

    # Physical selectors inside `Delete` expand live members only (the curves
    # must be unowned first — they are refused while the surface lives).
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Physical Curve(5) = {1,2};
        Delete{Surface{1};}
        Delete{Physical Curve{5};}
        """)
    model=execution.model
    @test !haskey(model.curves,1) && !haskey(model.curves,2)
    @test haskey(model.curves,3)

    # Wildcard `{:}` inside Delete removes everything of that dimension that
    # is not owned by a survivor.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Delete{Volume{1};}
        Delete{Surface{:};}
        """)
    @test isempty(execution.model.surfaces)
    @test !isempty(execution.model.curves)
end

@testset ".geo Delete physical-group staleness" begin
    # `.geo` deletion keeps stale physical member integers in the raw parser
    # registry and keeps names; the synchronized view filters to live members.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Physical Curve("border",5) = {1,2};
        Delete{Surface{1};}
        Delete{Curve{1,2,3,4};}
        """)
    model=execution.model
    # Stale raw member records persist in the parser registry, but the
    # synchronized view drops unresolvable members — the group vanishes.
    @test !haskey(model.physical,(1,5))
    @test model_physical_groups(model)==Tuple{Int,Int}[]
    # `getEntitiesForPhysicalName` errors when the name resolves to nothing.
    @test_throws ArgumentError model_entities_for_physical_name(model,"border")
    @test model.physical_names[(1,5)]=="border" # the name survives

    # Recreating a deleted tag resurrects its physical membership.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Physical Point(5) = {1};
        Delete{Surface{1};}
        Delete{Curve{1,2,3,4};}
        Delete{Point{1,2,3,4};}
        Point(1) = {0,0,0,0.1};
        """)
    model=execution.model
    @test model_entities_for_physical_group(model,0,5)==[1]
    @test model_physical_groups_for_entity(model,0,1)==[5]

    # `Delete Physicals` drops the records but keeps names and the counter.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Physical Curve("border",5) = {1,2};
        Delete Physicals;
        Physical Point("later") = {1};
        """)
    model=execution.model
    @test !haskey(model.physical,(1,5))        # the dropped group is gone
    @test model.physical[(0,6)]==[1]           # the new declaration stands
    @test model.physical_names[(1,5)]=="border"
    @test haskey(model.physical_names,(0,6))   # counter kept: auto tag is 6
end

@testset ".geo Delete Model/All/Variables/Options" begin
    # `Delete Model` clears geometry, physical records and counters while the
    # name table and user variables survive.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        x = 42;
        Physical Point("kept",7) = {1};
        Delete Model;
        Point(3) = {0,0,0,1};
        out[] = {newp, x};
        """)
    model=execution.model
    @test sort!(collect(keys(model.points)))==[3]
    @test model.next_tag[1]==3
    @test isempty(model.physical)
    @test model.physical_names[(0,7)]=="kept"
    @test execution.lists["out"]==[4.0,42.0]

    # `Delete All` additionally clears variables and the physical-name table.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        x = 42;
        Physical Point("gone",7) = {1};
        Delete All;
        Point(3) = {0,0,0,1};
        out[] = {newp};
        """)
    model=execution.model
    @test sort!(collect(keys(model.points)))==[3]
    @test isempty(model.physical_names)
    @test execution.lists["out"]==[4.0]

    err=_constraint_error(raw"""
        x = 1;
        Delete Variables;
        y = x;
        """)
    @test err isa ArgumentError

    # `Delete Meshes` is a no-op mid-execution; `Delete Options` resets the
    # tracked option state.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Mesh.TransfiniteTri = 1;
        Delete Options;
        Delete Meshes;
        """)
    @test execution.model.meshing.transfinite_tri==0
    @test execution.transfinite_tri==0

    # Named variable deletion; unknown names error like Gmsh.
    execution=_execute_constraint_source(raw"""
        x = 1;
        y[] = {1,2};
        Delete x;
        Delete y;
        z = 9;
        """)
    @test !haskey(execution.lists,"y")

    err=_constraint_error(raw"""
        Delete nothing_here;
        """)
    @test err isa ArgumentError
    @test occursin("Unknown object or expression to delete",string(err))
end

@testset ".geo SetMaxTag/SetTag" begin
    # `SetMaxTag` assigns the counter unconditionally (not max()).
    execution=_execute_constraint_source(raw"""
        Point(10) = {0,0,0,1};
        SetMaxTag Point(3);
        Point(newp) = {1,0,0,1};
        """)
    @test sort!(collect(keys(execution.model.points)))==[4,10]
    @test execution.model.next_tag[1]==4

    # The `GeoEntity{dim}` spelling selects the same counters.
    execution=_execute_constraint_source(raw"""
        Point(1) = {0,0,0,1};
        SetMaxTag GeoEntity{0}(30);
        Point(newp) = {1,0,0,1};
        """)
    @test sort!(collect(keys(execution.model.points)))==[1,31]

    # `SetTag` can never resolve mid-parse — Gmsh reports an unknown model
    # entity, so it is an explicit error here.
    err=_constraint_error(_GEO_SQUARE * raw"""
        SetTag Surface(1) = 9;
        """)
    @test err isa ArgumentError
    @test occursin("SetTag",string(err))
end

@testset ".geo constraints drive meshing" begin
    # Transfinite + Recombine produce Gmsh's structured patch: 36 nodes on a
    # 6-by-6 boundary grid.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Curve{1,3} = 6;
        Transfinite Curve{2,4} = 6;
        Transfinite Surface{1} = {1,2,3,4};
        """;mesh_dim=2)
    @test size(execution.mesh.coords,2)==36

    # A degenerated boundary curve collapses out of the transfinite side
    # chain, turning the surface into a three-sided transfinite patch.
    execution=_execute_constraint_source(raw"""
        Point(1) = {0,0,0,0.2};
        Point(2) = {1,0,0,0.2};
        Point(3) = {0.5,1,0,0.2};
        Line(1) = {1,2};
        Line(2) = {2,3};
        Line(3) = {3,1};
        Curve Loop(1) = {1,2,3};
        Plane Surface(1) = {1};
        Transfinite Curve{:} = 5;
        Transfinite Surface{1};
        """;mesh_dim=2)
    @test size(execution.mesh.coords,2)>0
end

@testset ".geo Mesh.FlexibleTransfinite" begin
    # `Mesh.FlexibleTransfinite` divides transfinite curve node counts by
    # `Mesh.CharacteristicLengthFactor` (upstream `meshGEdge`: `N /= lcFactor`,
    # truncating). Counts verified against Gmsh 4.15.2: 11 -> 11/5/3.
    flexible_square(options::AbstractString;count=11)=
        _execute_constraint_source(_GEO_SQUARE * """
            Transfinite Curve{:} = $count;
            Transfinite Surface{1};
            $options
            """;mesh_dim=2)
    @test size(flexible_square("Mesh.FlexibleTransfinite = 0;").mesh.coords,2)==121
    @test size(flexible_square("""
        Mesh.FlexibleTransfinite = 1;
        Mesh.CharacteristicLengthFactor = 2;
        """).mesh.coords,2)==25
    @test size(flexible_square("""
        Mesh.FlexibleTransfinite = 1;
        Mesh.CharacteristicLengthFactor = 3;
        """).mesh.coords,2)==9
    # `Mesh.MeshSizeFactor` is the same upstream option (`lcFactor`) under a
    # second name — writes to either alias drive the division.
    @test size(flexible_square("""
        Mesh.FlexibleTransfinite = 1;
        Mesh.MeshSizeFactor = 2;
        """).mesh.coords,2)==25
    # Non-positive `lcFactor` writes are ignored upstream (`if(val > 0)` in
    # `opt_mesh_lc_factor`), so the divisor stays 1.
    @test size(flexible_square("""
        Mesh.FlexibleTransfinite = 1;
        Mesh.CharacteristicLengthFactor = 0;
        """).mesh.coords,2)==121
    # Fractional factors truncate (`N /= lcFactor` on an int).
    @test size(flexible_square("""
        Mesh.FlexibleTransfinite = 1;
        Mesh.CharacteristicLengthFactor = 2.9;
        """).mesh.coords,2)==9
    # `Recombine` on an adjacent face forces the reduced count odd so blossom
    # can pair the boundary (12 -> 6 -> 7 per side); without the option the
    # declared count is used unchanged.
    @test size(flexible_square("""
        Recombine Surface{1};
        Mesh.FlexibleTransfinite = 1;
        Mesh.CharacteristicLengthFactor = 2;
        """;count=12).mesh.coords,2)==49
    @test size(flexible_square("""
        Recombine Surface{1};
        Mesh.CharacteristicLengthFactor = 2;
        """;count=12).mesh.coords,2)==144
    # `Mesh.RecombinationAlgorithm = 0` disables the odd-count forcing.
    @test size(flexible_square("""
        Recombine Surface{1};
        Mesh.FlexibleTransfinite = 1;
        Mesh.CharacteristicLengthFactor = 2;
        Mesh.RecombinationAlgorithm = 0;
        """;count=12).mesh.coords,2)==36
    # `Mesh.RecombineAll` forces odd counts without a per-surface flag.
    @test size(flexible_square("""
        Mesh.FlexibleTransfinite = 1;
        Mesh.CharacteristicLengthFactor = 2;
        Mesh.RecombineAll = 1;
        """;count=12).mesh.coords,2)==49
    # A factor that truncates the count below two clamps to the endpoint-only
    # curve, matching Gmsh's degenerate corner grid (3 -> 1 -> 4 nodes).
    @test size(flexible_square("""
        Mesh.FlexibleTransfinite = 1;
        Mesh.CharacteristicLengthFactor = 2;
        """;count=3).mesh.coords,2)==4
    # The option state mirrors onto the model for `mesh_dim` generation.
    execution=flexible_square("""
        Mesh.FlexibleTransfinite = 1;
        Mesh.CharacteristicLengthFactor = 2;
        """)
    @test execution.model.meshing.flexible_transfinite==true
    @test execution.model.meshing.lc_factor==2.0
    @test execution.model.meshing.recombine_algo==1
end

@testset ".geo flexible transfinite density and recombination" begin
    # Gmsh 4.15.2 meshGEdge uses the declared count (divided by lcFactor)
    # in F_Transfinite, then adjusts the output count only when its numerical
    # density primitive exceeds 0.75. These controls were captured from
    # actual native Line meshes, including both neighbors of the threshold.
    function flexible_line(law,coefficient,count,factor,algorithm,all_faces)
        flags=all_faces ? "Mesh.RecombineAll=1;" : "Recombine Surface{1};"
        return _execute_constraint_source(_GEO_SQUARE * """
            Mesh.FlexibleTransfinite=1;
            Mesh.CharacteristicLengthFactor=$factor;
            Mesh.RecombinationAlgorithm=$algorithm;
            $flags
            Transfinite Curve{1}=$count Using $law $coefficient;
            Transfinite Curve{2,3,4}=3;
            Mesh 1;
            """)
    end
    for (coefficient,blossom_count) in
            ((0.5,6),(prevfloat(0.75),6),(0.75,6),(nextfloat(0.75),6),(0.9,7)),
        algorithm in (0,1),all_faces in (false,true)
        execution=flexible_line("Beta",coefficient,12,2.0,algorithm,all_faces)
        part=geo_entity_mesh(execution,1,1)
        count=algorithm==0 ? 6 : blossom_count
        @test nnodes(part)==count
        @test nsegs(part)==count-1
        @test validate(part).ok
        @test part.coords[1,:]≈collect(range(0.0,1.0;length=count)) atol=2e-11 rtol=0
        @test execution.model.curve_params[1]≈part.coords[1,:] atol=0 rtol=0
        @test all(diff(execution.model.curve_params[1]).>0)
        @test Tessella.Model._flexible_transfinite_nodes(
            execution.model,12,1,"density control")==count
    end

    # The density still uses six points when recombination produces seven;
    # rebuilding the law with seven would give different actual positions.
    progression_samples=(
        (0,Float64[0,0.0029325535916704293,0.014662747401729634,
                    0.061583559284854296,0.24926683194892127,1]),
        (1,Float64[0,0.00244379307574998,0.010752676685410667,
                    0.03812314575766151,0.12414462647064331,
                    0.3743889670179139,1]))
    for (algorithm,expected) in progression_samples
        execution=flexible_line("Progression",4.0,3,0.5,algorithm,false)
        part=geo_entity_mesh(execution,1,1)
        @test part.coords[1,:]≈expected atol=2e-11 rtol=0
        @test nsegs(part)==length(expected)-1
        @test execution.model.meshing.transfinite_curves[1].num_nodes==3
        @test execution.model.curve_params[1]≈part.coords[1,:] atol=0 rtol=0
    end
    for (law,coefficient) in (("Bump",2.0),("Beta",1.1),
                             ("Progression_HWall",0.1),
                             ("Bump_HWall",0.01),("Beta_HWall",0.05))
        execution=flexible_line(law,coefficient,7,0.5,1,false)
        part=geo_entity_mesh(execution,1,1)
        @test nnodes(part)==15
        @test nsegs(part)==14
        @test validate(part).ok
        @test all(diff(execution.model.curve_params[1]).>0)
    end
    # Positive HWall types transform signed/zero wall heights before the
    # mass gate. A negative coefficient denotes the opposite wall, whereas
    # a negative type (reversed curve declaration) follows another arm.
    for (law,coefficient) in (("Progression_HWall",-0.1),("Bump_HWall",-0.01),
                             ("Beta_HWall",-0.05),("Progression_HWall",0.0),
                             ("Bump_HWall",0.0),("Beta_HWall",0.0))
        execution=flexible_line(law,coefficient,12,2.0,1,false)
        part=geo_entity_mesh(execution,1,1)
        @test nnodes(part)==7
        @test nsegs(part)==6
        @test validate(part).ok
        @test all(diff(execution.model.curve_params[1]).>0)
    end
    # The truncated law count of one needs no interior primitive marks.
    execution=flexible_line("Progression",4.0,3,2.0,1,false)
    @test geo_entity_mesh(execution,1,1).coords[1,:]==[0.0,1.0]

    # Flexible Circle controls retain the declared-count HWall transform
    # and the divided density-law count, independently of the final odd N.
    # These native OCC Circle samples were captured from Gmsh 4.15.2; the
    # established curved-law comparison bound is 3e-7.
    for (law,coefficient,expected) in (
            ("Progression",2.0,Float64[0,0.1689028348257925,0.47292790152443426,
                1.0134169226425884,1.9592727403206374,3.5807399265602187]),
            ("Progression_HWall",-0.1,Float64[0,1.6898564402531275,
                3.0526184812969497,4.1476425356450335,5.024710690883849,
                5.725188895849936]))
        execution=_execute_constraint_source("""
            SetFactory("OpenCASCADE");Cylinder(1)={0,0,0,0,0,1,1};
            Mesh.FlexibleTransfinite=1;Mesh.MeshSizeFactor=2;
            Mesh.RecombineAll=1;Mesh.RecombinationAlgorithm=1;
            Transfinite Curve{1}=12 Using $law $coefficient;Mesh 1;
            """)
        part=geo_entity_mesh(execution,1,1)
        stored_params=execution.model.curve_params[1]
        params=stored_params[1:end-1]
        @test nnodes(part)==nsegs(part)==length(expected)==6
        @test stored_params[end]==2pi
        @test maximum(abs.(params.-expected))<=3e-7
        @test all(diff(params).>0)
        @test part.segs[:,end]==Int32[6,1]
        @test validate(part).ok
        for node in axes(part.coords,2)
            @test maximum(abs.(part.coords[:,node].-
                [cos(expected[node]),sin(expected[node]),1.0]))<=3e-7
        end
        @test execution.model.meshing.transfinite_curves[1].num_nodes==12
    end
    # The actual Circle primitive has a different numerical threshold from
    # the unit native Line above. Saved OCC controls keep five unique nodes
    # for Beta0.5 and six for all three 0.75 neighbors and Beta0.9.
    for (coefficient,unique_count) in ((0.5,5),(prevfloat(0.75),6),
                                      (0.75,6),(nextfloat(0.75),6),(0.9,6))
        execution=_execute_constraint_source("""
            SetFactory("OpenCASCADE");Cylinder(1)={0,0,0,0,0,1,1};
            Mesh.FlexibleTransfinite=1;Mesh.MeshSizeFactor=2;
            Mesh.RecombineAll=1;Mesh.RecombinationAlgorithm=1;
            Transfinite Curve{1}=12 Using Beta $coefficient;Mesh 1;
            """)
        part=geo_entity_mesh(execution,1,1)
        params=execution.model.curve_params[1]
        @test nnodes(part)==nsegs(part)==unique_count
        @test length(params)==unique_count+1
        @test params[end]==2pi
        expected=Float64[2pi*k/unique_count for k in 0:unique_count]
        @test maximum(abs.(params.-expected))<=2e-11
        @test Tessella.Model._flexible_transfinite_nodes(
            execution.model,12,1,"Circle mass control")==unique_count+1
    end
end

@testset ".geo ordinary curve recombination counts" begin
    # Independent Gmsh 4.15.2 controls: meshGEdge's increaseN is an identity,
    # and the odd-N adjustment is disabled below primitive mass 0.75.
    # The digest covers actual directed Line1 coordinates/connectivity,
    # elementary node/cell carriers and stored native parameters in curve
    # order, keeping legacy global dense tag numbering a separate contract.
    function curve_digest(model,part)
        chain=Int[part.segs[1,1]]
        for cell in axes(part.segs,2)
            @test part.segs[1,cell]==chain[end]
            push!(chain,Int(part.segs[2,cell]))
        end
        @test allunique(chain)
        owners=Tessella.GeoExec._geo_mesh_part_node_entities(
            model,[(1,1,part)],"ordinary curve control")[1]
        bytes=IOBuffer()
        write(bytes,htol(UInt32(length(chain))),htol(UInt32(nsegs(part))),UInt8(1))
        for node in chain
            for x in part.coords[:,node]
                write(bytes,htol(reinterpret(UInt64,x)))
            end
            for carrier in owners[node]
                write(bytes,htol(Int32(carrier)))
            end
        end
        for cell in axes(part.segs,2)
            write(bytes,htol(Int32(cell)),htol(Int32(cell+1)),
                  htol(Int32(1)),htol(Int32(1)))
        end
        for parameter in model.curve_params[1]
            write(bytes,htol(reinterpret(UInt64,parameter)))
        end
        return bytes2hex(sha256(take!(bytes))),chain,owners
    end
    for size in (0.18,2.0),algorithm in (0,1,2,4)
        execution=_execute_constraint_source("""
            Mesh.RecombineAll=1;Mesh.RecombinationAlgorithm=$algorithm;
            Mesh.MeshSizeMin=$size;Mesh.MeshSizeMax=$size;Mesh.MeshSizeFromPoints=0;
            Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};Mesh 1;
            """)
        part=geo_entity_mesh(execution,1,1)
        digest,chain,owners=curve_digest(execution.model,part)
        expected_count=size==0.18 ? 7 : 2
        expected_digest=size==0.18 ?
            "664ada8df794669435356a46cb39bfeac867612c74e5a56c4f4c0e83538c8968" :
            "f7f5b024dbb79850422e5e5dd4b816b64f4d23df516ebd1e4890106362018247"
        @test length(chain)==expected_count
        @test nsegs(part)==expected_count-1
        @test [owners[node] for node in chain]==
            [(0,1);fill((1,1),expected_count-2);(0,2)]
        @test execution.model.curve_params[1]==[part.coords[1,node] for node in chain]
        @test execution.model.curve_params[1]≈
            collect(range(0.0,1.0;length=expected_count)) atol=2e-11 rtol=0
        @test validate(part).ok
        @test digest==expected_digest
    end
    # A positive primitive with an initially even N exercises the actual
    # adjustment: native size 0.2 gives N6 without blossom and N7 with it.
    for algorithm in (0,1,2,4)
        execution=_execute_constraint_source("""
            Mesh.RecombineAll=1;Mesh.RecombinationAlgorithm=$algorithm;
            Mesh.MeshSizeMin=0.2;Mesh.MeshSizeMax=0.2;Mesh.MeshSizeFromPoints=0;
            Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};Mesh 1;
            """)
        part=geo_entity_mesh(execution,1,1)
        expected_count=algorithm==0 ? 6 : 7
        @test nnodes(part)==expected_count
        @test nsegs(part)==expected_count-1
        @test execution.model.curve_params[1]≈
            collect(range(0.0,1.0;length=expected_count)) atol=2e-11 rtol=0
    end
    # Gmsh's native Line derivative uses a 1e-5 stencil. At these three
    # independently captured sizes its integrated mass is just below 0.75,
    # even when the placement primitive rounds above that threshold.
    for size in (prevfloat(4/3),4/3,nextfloat(4/3))
        execution=_execute_constraint_source("""
            Mesh.RecombineAll=1;Mesh.RecombinationAlgorithm=1;
            Mesh.MeshSizeMin=$size;Mesh.MeshSizeMax=$size;Mesh.MeshSizeFromPoints=0;
            Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};Mesh 1;
            """)
        part=geo_entity_mesh(execution,1,1)
        digest,chain,owners=curve_digest(execution.model,part)
        @test nnodes(part)==2 && nsegs(part)==1
        @test part.coords[:,chain]==Float64[0 1;0 0;0 0]
        @test owners[chain]==[(0,1),(0,2)]
        @test execution.model.curve_params[1]==[0.0,1.0]
        @test digest=="f7f5b024dbb79850422e5e5dd4b816b64f4d23df516ebd1e4890106362018247"
    end
    # Count-only native integration reuses the already smoothed sizes.
    # Fixing the resulting edge count explicitly must reproduce the entire
    # field-call trace and placement, including filterPoints callbacks.
    function sampled_curve(;exact_edges=nothing)
        trace=Float64[]
        count_calls=Ref(0)
        smoothed=Ref(false)
        law=x->x<=0.5 ? 0.2 : 0.3
        field=Tessella.SizeField.FunctionSize((x,y,z)->(push!(trace,x);law(x)))
        primitive=function(samples)
            count_calls[]+=1
            smoothed[]=any(p->p.h!=law(p.t),samples)
            calls_before=length(trace)
            value=Tessella.Model._model_native_line_recombination_mass(
                (0.0,0.0,0.0),(1.0,0.0,0.0),samples)
            @test length(trace)==calls_before
            return value
        end
        points,parameters=Tessella.Mesh1D.mesh_curve(t->(t,0.0,0.0),field;
            derivative=t->(1.0,0.0,0.0),nsample=32,smooth_ratio=1.01,
            force_odd=true,_recombine_odd_nodes=true,
            _recombination_primitive=primitive,exact_edges=exact_edges)
        return points,parameters,trace,count_calls[],smoothed[]
    end
    counted=sampled_curve()
    fixed=sampled_curve(;exact_edges=length(counted[2])-1)
    @test counted[4]==1
    @test counted[5]
    @test fixed[4]==0
    @test counted[1:3]==fixed[1:3]
    # Upstream N includes both endpoints even for a closed curve; the
    # published mesh omits its repeated endpoint. A Cylinder's seam circles
    # therefore contain seven unique nodes for algorithm0 and eight for
    # algorithms1/2/4, as independently captured from Gmsh 4.15.2.
    for algorithm in (0,1,2,4)
        execution=_execute_constraint_source("""
            SetFactory("OpenCASCADE");Cylinder(1)={0,0,0,0,0,1,1};
            Mesh.RecombineAll=1;Mesh.RecombinationAlgorithm=$algorithm;
            Mesh.MeshSizeMin=1;Mesh.MeshSizeMax=1;Mesh.MeshSizeFromPoints=0;
            Mesh 1;
            """)
        expected_count=algorithm==0 ? 7 : 8
        for curve in (1,3)
            part=geo_entity_mesh(execution,1,curve)
            model=execution.model
            params=model.curve_params[curve]
            owners=Tessella.GeoExec._geo_mesh_part_node_entities(
                model,[(1,curve,part)],"closed curve control")[1]
            @test nnodes(part)==nsegs(part)==expected_count
            @test part.segs[:,end]==Int32[expected_count,1]
            @test model.curves[curve][1]==model.curves[curve][2]
            @test owners[1]==(0,model.curves[curve][1])
            @test owners[2:end]==fill((1,curve),expected_count-1)
            @test all(diff(params).>0)
            expected_params=Float64[2pi*k/expected_count for k in 0:expected_count-1]
            @test maximum(abs.(params.-expected_params))<=2e-11
            @test all(node->isapprox(hypot(part.coords[1,node],part.coords[2,node]),
                                    1.0;atol=4eps(Float64),rtol=0),
                      axes(part.coords,2))
            @test validate(part).ok
        end
    end
end

@testset ".geo model-level recombination" begin
    # Unstructured `Recombine Surface` emits a MixedMesh with quadrangle and
    # leftover-triangle blocks per entity. The default (blossom) algorithm
    # pairs a boundary regraded to the odd-node count (24 segments here).
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Recombine Surface{1};
        """;mesh_dim=2)
    @test execution.mesh isa MixedMesh
    nquad=sum(size(b.nodes,2) for b in execution.mesh.blocks if b.msh==3)
    ntri=sum(size(b.nodes,2) for b in execution.mesh.blocks if b.msh==2)
    @test (nquad,ntri)==(47,4)
    @test sum(size(b.nodes,2) for b in execution.mesh.blocks if b.msh==1)==24

    # `Mesh.RecombinationAlgorithm = 0` keeps the unforced boundary grading,
    # so every quadrangle replaces exactly two triangles of the simplex
    # product generated on the same grading.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Recombine Surface{1};
        Mesh.RecombinationAlgorithm = 0;
        """;mesh_dim=2)
    @test execution.mesh isa MixedMesh
    nquad=sum(size(b.nodes,2) for b in execution.mesh.blocks if b.msh==3)
    ntri=sum(size(b.nodes,2) for b in execution.mesh.blocks if b.msh==2)
    @test nquad>0
    execution_simplex=_execute_constraint_source(_GEO_SQUARE;mesh_dim=2)
    @test 2*nquad+ntri==ntris(execution_simplex.mesh)

    # The full-quad arms (2/3) run upstream's triangle-subdivision pipeline,
    # which the pairing kernel does not implement — the rejection is explicit.
    err=_constraint_error(_GEO_SQUARE * raw"""
        Recombine Surface{1};
        Mesh.RecombinationAlgorithm = 3;
        """;mesh_dim=2)
    @test err isa ArgumentError
    @test occursin("RecombinationAlgorithm",err.msg)

    # Transfinite + Recombine gives an all-quadrangle MixedMesh — Gmsh's
    # structured patch, not the simplex path.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        Transfinite Curve{:} = 5;
        Transfinite Surface{1};
        Recombine Surface{1};
        """;mesh_dim=2)
    @test execution.mesh isa MixedMesh
    quads=only([b for b in execution.mesh.blocks if b.msh==3])
    @test size(quads.nodes,2)==16
    @test !any(b.msh==2 for b in execution.mesh.blocks)

    # A fully recombined transfinite volume emits hexahedra plus quadrangle
    # boundary sheets through the MixedMesh merge.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Transfinite Curve{:} = 4;
        Transfinite Surface{:};
        Recombine Surface{:};
        Transfinite Volume{1};
        """;mesh_dim=3)
    @test execution.mesh isa MixedMesh
    hexes=only([b for b in execution.mesh.blocks if b.msh==5])
    @test size(hexes.nodes,2)==27
    @test sum(size(b.nodes,2) for b in execution.mesh.blocks if b.msh==3)==54
    @test !any(b.msh==4 for b in execution.mesh.blocks)
    @test !any(b.msh==2 for b in execution.mesh.blocks)

    # `Mesh.RecombineAll` drives the same volume mask without per-surface
    # flags.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Transfinite Curve{:} = 4;
        Transfinite Surface{:};
        Transfinite Volume{1};
        Mesh.RecombineAll = 1;
        """;mesh_dim=3)
    @test execution.mesh isa MixedMesh
    @test any(b.msh==5 && size(b.nodes,2)==27 for b in execution.mesh.blocks)

    # An unrecombined opposite face pair emits prisms through the same mask
    # path (x=0/x=1 faces left simplex).
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Transfinite Curve{:} = 4;
        Transfinite Surface{:};
        Recombine Surface{1,2,3,5};
        Transfinite Volume{1};
        """;mesh_dim=3)
    @test execution.mesh isa MixedMesh
    @test any(b.msh==6 && size(b.nodes,2)==54 for b in execution.mesh.blocks)

    # `Mesh 3` on a recombined boundary folds boundary quadrangles into the
    # PLC and still produces a valid all-tet interior — the merged product
    # is mixed because the surface parts keep their quadrangle blocks.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Recombine Surface{:};
        """;mesh_dim=3)
    @test execution.mesh isa MixedMesh
    ntets=sum(size(b.nodes,2) for b in execution.mesh.blocks if b.msh==4)
    @test ntets>0
    @test validate(execution.mesh).ok

    # Compact transfinite prism volumes prune behind-diagonal slots exactly
    # like Gmsh's writer — every emitted node is referenced by a cell.
    prism_source=raw"""
        Point(1) = {0,0,0,0.3};
        Point(2) = {1,0,0,0.3};
        Point(3) = {0,1,0,0.3};
        Point(4) = {0,0,1,0.3};
        Point(5) = {1,0,1,0.3};
        Point(6) = {0,1,1,0.3};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,1};
        Line(4) = {4,5}; Line(5) = {5,6}; Line(6) = {6,4};
        Line(7) = {1,4}; Line(8) = {2,5}; Line(9) = {3,6};
        Curve Loop(1) = {1,2,3};
        Curve Loop(2) = {4,5,6};
        Curve Loop(3) = {1,8,-4,-7};
        Curve Loop(4) = {2,9,-5,-8};
        Curve Loop(5) = {3,7,-6,-9};
        Plane Surface(1) = {1}; Plane Surface(2) = {2};
        Plane Surface(3) = {3}; Plane Surface(4) = {4};
        Plane Surface(5) = {5};
        Surface Loop(1) = {1,2,3,4,5};
        Volume(1) = {1};
        Mesh.TransfiniteTri = 1;
        Transfinite Curve{:} = 4;
        Transfinite Surface{:};
        Recombine Surface{:};
        Transfinite Volume{1};
        """
    execution=_execute_constraint_source(prism_source;mesh_dim=3)
    @test execution.mesh isa MixedMesh
    # 46 merged nodes like Gmsh's written set — two compact-tab interior
    # slots behind the diagonal are unreferenced orphan nodes with their
    # own evaluated coordinates.
    @test size(execution.mesh.coords,2)==46
    used=falses(size(execution.mesh.coords,2))
    for b in execution.mesh.blocks
        used[b.nodes[:]].=true
    end
    @test count(used)==44

    # Simplex-only mesh statements reject a mixed product explicitly.
    err=_constraint_error(_GEO_SQUARE * raw"""
        Transfinite Curve{:} = 5;
        Transfinite Surface{1};
        Recombine Surface{1};
        Mesh 2;
        RefineMesh;
        """)
    @test err isa ArgumentError
end

@testset ".geo classified projection of recombined volumes" begin
    function _block_counts(mesh::MixedMesh)
        counts=Dict{Int,Int}()
        for block in mesh.blocks
            block isa ElementBlock || continue
            counts[Int(block.msh)]=get(counts,Int(block.msh),0)+
                size(block.nodes,2)
        end
        return counts
    end

    # All-hexahedral transfinite box: the classified projection keeps the
    # native hex block, emits quad boundary faces, and classifies the full
    # entity hierarchy.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Transfinite Curve{:} = 4;
        Transfinite Surface{:};
        Recombine Surface{:};
        Transfinite Volume{1} = {1,2,3,4,5,6,7,8};
        Mesh 3;
        """;mesh_dim=3)
    part=geo_entity_mesh(execution,3,1)
    @test part isa MixedMesh
    projected=model_to_mixed(execution.model,part,3,1)
    @test projected isa MixedMesh
    @test validate(projected).ok
    counts=_block_counts(projected)
    @test counts[15]==8 && counts[1]==36
    @test counts[3]==54 && counts[5]==27
    @test !haskey(counts,2) && !haskey(counts,4)
    entities=projected.entity_data.entities
    @test count(e->e.dim==0,values(entities))==8
    @test count(e->e.dim==1,values(entities))==12
    @test count(e->e.dim==2,values(entities))==6
    @test count(e->e.dim==3,values(entities))==1
    @test entities[(3,1)].boundaries==Int32[1,2,3,4,5,6]

    # Recombined five-face prism: quads and triangle caps coexist, and the
    # folded-quad diagonal on a transfinite side face does not mint a
    # phantom surface crease during fill certification.
    prism_source=raw"""
        Point(1) = {0,0,0,0.4}; Point(2) = {1,0,0,0.4};
        Point(3) = {0,1,0,0.4}; Point(4) = {0,0,1,0.4};
        Point(5) = {1,0,1,0.4}; Point(6) = {0,1,1,0.4};
        Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,1};
        Line(4) = {4,5}; Line(5) = {5,6}; Line(6) = {6,4};
        Line(7) = {1,4}; Line(8) = {2,5}; Line(9) = {3,6};
        Curve Loop(1) = {1,2,3}; Curve Loop(2) = {4,5,6};
        Curve Loop(3) = {1,8,-4,-7}; Curve Loop(4) = {2,9,-5,-8};
        Curve Loop(5) = {3,7,-6,-9};
        Plane Surface(1) = {1}; Plane Surface(2) = {2};
        Plane Surface(3) = {3}; Plane Surface(4) = {4};
        Plane Surface(5) = {5};
        Surface Loop(1) = {1,2,3,4,5};
        Volume(1) = {1};
        Transfinite Curve{:} = 4;
        Transfinite Surface{:};
        Recombine Surface{3,4,5};
        Transfinite Volume{1} = {1,2,3,4,5,6};
        Mesh 3;
        """
    execution=_execute_constraint_source(prism_source;mesh_dim=3)
    part=geo_entity_mesh(execution,3,1)
    @test part isa MixedMesh
    projected=model_to_mixed(execution.model,part,3,1)
    @test projected isa MixedMesh
    @test validate(projected).ok
    counts=_block_counts(projected)
    @test counts[15]==6 && counts[1]==27
    @test counts[3]==27 && counts[2]==30 && counts[6]==45
    entities=projected.entity_data.entities
    @test count(e->e.dim==0,values(entities))==6
    @test count(e->e.dim==1,values(entities))==9
    @test count(e->e.dim==2,values(entities))==5
    @test entities[(3,1)].boundaries==Int32[1,2,3,4,5]

    # Periodic boundary surfaces still emit the affine node links through
    # quad faces on a recombined volume.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Periodic Surface 2 {5,6,7,8} = 1 {1,2,3,4};
        Transfinite Curve{:} = 4;
        Transfinite Surface{:};
        Recombine Surface{:};
        Transfinite Volume{1} = {1,2,3,4,5,6,7,8};
        Mesh 3;
        """;mesh_dim=3)
    part=geo_entity_mesh(execution,3,1)
    projected=model_to_mixed(execution.model,part,3,1)
    @test validate(projected).ok
    surface_links=[l for l in projected.periodic_links if l.dim==2]
    @test length(surface_links)==1
    @test only(surface_links).slave_entity==2
    @test only(surface_links).master_entity==1
    @test length(only(surface_links).slave_nodes)==16

    # Physical groups ride through quad surface blocks and hex cells.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        Physical Surface(10) = {1,3,5};
        Physical Volume(20) = {1};
        Transfinite Curve{:} = 4;
        Transfinite Surface{:};
        Recombine Surface{:};
        Transfinite Volume{1} = {1,2,3,4,5,6,7,8};
        Mesh 3;
        """;mesh_dim=3)
    part=geo_entity_mesh(execution,3,1)
    projected=model_to_mixed(execution.model,part,3,1)
    quad_tags=Int32[]
    hex_tags=Int32[]
    for block in projected.blocks
        block isa ElementBlock || continue
        block.msh==3 && append!(quad_tags,block.tags)
        block.msh==5 && append!(hex_tags,block.tags)
    end
    @test 10 in quad_tags
    @test all(==(20),hex_tags)

    # A MixedMesh carrying non-volume blocks is not a volume part.
    bad=MixedMesh(part.coords,
                  [ElementBlock(2,ones(Int32,3,1),Int32[0]),
                   first(b for b in part.blocks if b isa ElementBlock)])
    @test_throws ArgumentError model_to_mixed(execution.model,bad,3,1)

    # Simplex volumes still project through the shared path.
    execution=_execute_constraint_source(_GEO_BOX * "Mesh 3;\n";mesh_dim=3)
    part=geo_entity_mesh(execution,3,1)
    @test part isa Mesh
    projected=model_to_mixed(execution.model,part,3,1)
    @test validate(projected).ok
    counts=_block_counts(projected)
    @test counts[4]>0 && counts[2]>0 && !haskey(counts,3) && !haskey(counts,5)
end

# `.geo` `Extrude ... Layers` — the structured sweep upstream's
# `meshGRegionExtruded`/`meshGFaceExtruded`/`SubdivideExtrudedMesh`
# produces: lateral quad strips or constrained triangle pairs, verbatim
# top copies, prisms/hexahedra under `Recombine`, and the global
# diagonal-compatible three-tetrahedron subdivision otherwise. Element
# counts below are the Gmsh 4.15.2 outputs for the same files.
const _GEO_EXTRUDE_TRI = raw"""
Point(1) = {0,0,0,0.5};
Point(2) = {1,0,0,0.5};
Point(3) = {0,1,0,0.5};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,1};
Curve Loop(1) = {1,2,3};
Plane Surface(1) = {1};
"""

function _extrude_block_counts(mesh)
    counts=Dict{Int,Int}()
    mesh isa MixedMesh || return counts
    for block in mesh.blocks
        block isa ElementBlock || continue
        counts[block.msh]=get(counts,block.msh,0)+size(block.nodes,2)
    end
    return counts
end

@testset ".geo Extrude Layers sweep" begin
    @testset "triangle source, Recombine -> prisms" begin
        execution=_execute_constraint_source(_GEO_EXTRUDE_TRI * raw"""
            Extrude{0,0,2}{Surface{1};Layers{3};Recombine;}
            """;mesh_dim=3)
        counts=_extrude_block_counts(execution.mesh)
        # Gmsh: 21 quads (lateral strips), 14 tris (bottom+top copies),
        # 21 prisms, 23 segs.
        @test counts[3]==21
        @test counts[2]==14
        @test counts[6]==21
        @test counts[1]==23
        @test nnodes(execution.mesh)==32
        volume=geo_entity_mesh(execution,3,1)
        @test volume isa MixedMesh
        @test only(size(b.nodes,2) for b in volume.blocks
                   if b isa ElementBlock && b.msh==6)==21
        # The classified projection accepts the swept volume part.
        projected=model_to_mixed(execution.model,volume,3,1)
        @test projected isa MixedMesh
        pcounts=_extrude_block_counts(projected)
        @test pcounts[6]==21
        @test pcounts[3]==21
    end
    @testset "triangle source, no Recombine -> tets" begin
        execution=_execute_constraint_source(_GEO_EXTRUDE_TRI * raw"""
            Extrude{0,0,2}{Surface{1};Layers{3};}
            """;mesh_dim=3)
        # Gmsh: 23 segs, 56 tris (14 copies + 42 lateral), 63 tets
        # (7 triangles x 3 layers x 3).
        @test execution.mesh isa Mesh
        @test nsegs(execution.mesh)==23
        @test ntris(execution.mesh)==56
        @test ntets(execution.mesh)==63
        @test nnodes(execution.mesh)==32
        @test validate(execution.mesh).ok
        # Determinism — the sweep must reproduce bitwise.
        again=_execute_constraint_source(_GEO_EXTRUDE_TRI * raw"""
            Extrude{0,0,2}{Surface{1};Layers{3};}
            """;mesh_dim=3)
        @test again.mesh.coords==execution.mesh.coords
        @test again.mesh.tets==execution.mesh.tets
        @test again.mesh.tris==execution.mesh.tris
    end
    @testset "quadrangle source, Recombine -> hexahedra" begin
        execution=_execute_constraint_source(_GEO_SQUARE * raw"""
            Recombine Surface{1};
            Extrude{0,0,2}{Surface{1};Layers{3};Recombine;}
            """;mesh_dim=3)
        counts=_extrude_block_counts(execution.mesh)
        # The swept hexahedra outnumber any recombination-leftover
        # prisms and the laterals are pure quadrangle strips.
        @test counts[5]>0
        volume=geo_entity_mesh(execution,3,1)
        @test volume isa MixedMesh
        btypes=sort!(collect(Set(b.msh for b in volume.blocks
                                 if b isa ElementBlock)))
        @test 5 in btypes
    end
    @testset "layer groups and heights" begin
        execution=_execute_constraint_source(_GEO_EXTRUDE_TRI * raw"""
            Extrude{0,0,2}{Surface{1};Layers{{2,3},{0.4,1.0}};Recombine;}
            """;mesh_dim=3)
        counts=_extrude_block_counts(execution.mesh)
        # Gmsh: 35 prisms (7 triangles x (2+3) levels), 35 lateral quads,
        # 14 tris, 29 segs, 48 nodes.
        @test counts[6]==35
        @test counts[3]==35
        @test counts[2]==14
        @test counts[1]==29
        @test nnodes(execution.mesh)==48
        # The intermediate level lands at the normalized group height —
        # z = 2*(0.4 + k/3*0.6) for k=1,2.
        zs=sort!(unique(round.(execution.mesh.coords[3,:];digits=9)))
        @test 2*0.4+2*1.2/3*1 in zs
        @test 2*(0.4+2/3*0.6) in zs
    end
    @testset "rotation extrusion" begin
        execution=_execute_constraint_source(raw"""
            Point(1)={1,0,0,0.4};
            Point(2)={2,0,0,0.4};
            Point(3)={2,1,0,0.4};
            Point(4)={1,1,0,0.4};
            Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
            Curve Loop(1)={1,2,3,4};
            Plane Surface(1)={1};
            Extrude{{0,1,0},{0,0,0},Pi/2}{Surface{1};Layers{4};}
            """;mesh_dim=3)
        # A quarter-turn twist sweep — interior nodes must follow the arc,
        # not the chord.
        @test validate(execution.mesh).ok
        zs=execution.mesh.coords[3,:]
        xs=execution.mesh.coords[1,:]
        @test any(>(0.5),abs.(xs))
        # Arc points are off the extrusion plane (y stays, x/z rotate).
        @test all(>=(0.0),round.(execution.mesh.coords[2,:];digits=9))
    end
    @testset "QuadTriAddVerts transition and NoNewVerts source-grid guard" begin
        execution=_execute_constraint_source(_GEO_EXTRUDE_TRI * raw"""
            Extrude{0,0,2}{Surface{1};Layers{3};Recombine;QuadTriAddVerts;}
            """;mesh_dim=3)
        counts=_extrude_block_counts(execution.mesh)
        @test counts[4]>0
        @test counts[7]>0
        @test validate(execution.mesh).ok
        err=_constraint_error(_GEO_EXTRUDE_TRI * raw"""
            Extrude{0,0,2}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;}
            """;mesh_dim=3)
        @test err !== nothing
        @test occursin("QuadTriNoNewVerts",sprint(showerror,err))
    end
    @testset "no Layers falls back to unstructured filling" begin
        execution=_execute_constraint_source(_GEO_EXTRUDE_TRI * raw"""
            Extrude{0,0,2}{Surface{1};}
            """;mesh_dim=3)
        @test execution.mesh isa Mesh
        @test ntets(execution.mesh)>0
        @test nsegs(execution.mesh)>0
    end
end
