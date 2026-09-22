using Test
using Tessella
using Tessella.Model: model_physical_groups, model_entities_for_physical_group,
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

    # TransfQuadTri wildcard flags every current volume; meshing a flagged
    # volume is the native-kernel QuadTri blocker.
    execution=_execute_constraint_source(_GEO_BOX * raw"""
        TransfQuadTri{:};
        """)
    @test 1 in execution.model.meshing.quad_tri
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
    @test meshing.size_from_boundary[(2,1)]

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

    # `MeshSizeFromBoundary = 0` clears the flag.
    execution=_execute_constraint_source(_GEO_SQUARE * raw"""
        MeshSizeFromBoundary Surface{1} = 1;
        MeshSizeFromBoundary Surface{1} = 0;
        """)
    @test !haskey(execution.model.meshing.size_from_boundary,(2,1))

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
