using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, validate, nnodes, ntris

function _periodic_geo_square(periodic_statement::AbstractString;
                              mesh_size=0.5)
    return """
        Point(1) = {0, 0, 0, $mesh_size};
        Point(2) = {1, 0, 0, $mesh_size};
        Point(3) = {1, 1, 0, $mesh_size};
        Point(4) = {0, 1, 0, $mesh_size};
        Line(1) = {1, 2};
        Line(2) = {2, 3};
        Line(3) = {3, 4};
        Line(4) = {4, 1};
        Curve Loop(1) = {1, 2, 3, 4};
        Plane Surface(1) = {1};
        $periodic_statement
        """
end

function _execute_geo_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _geo_periodic_affine_point(affine,coordinates)
    x,y,z=coordinates
    return (
        affine[4]+muladd(affine[3],z,
                         muladd(affine[2],y,affine[1]*x)),
        affine[8]+muladd(affine[7],z,
                         muladd(affine[6],y,affine[5]*x)),
        affine[12]+muladd(affine[11],z,
                          muladd(affine[10],y,affine[9]*x)),
    )
end

function _periodic_geo_rotation()
    return """
        Point(1) = {2, 1, 0, 0.5};
        Point(2) = {3, 1, 0, 0.5};
        Point(3) = {1, 2, 0, 0.5};
        Point(4) = {1, 3, 0, 0.5};
        Line(1) = {1, 2};
        Line(2) = {2, 4};
        Line(3) = {3, 4};
        Line(4) = {3, 1};
        Curve Loop(1) = {1, 2, -3, 4};
        Plane Surface(1) = {1};
        axisZ = 2 * 0.5;
        centerX = Cos(0);
        centerY = Sqrt(1);
        quarterTurn = Pi / 2;
        Periodic Curve {3} = {1}
          Rotate {{Atan2(0, 1), Sin(0), axisZ},
                  {centerX, centerY, 0}, quarterTurn};
        """
end

function _periodic_geo_embedded_curves()
    return """
        Point(1) = {0, 0, 0, 0.5};
        Point(2) = {1, 0, 0, 0.5};
        Point(3) = {1, 1, 0, 0.5};
        Point(4) = {0, 1, 0, 0.5};
        Line(1) = {1, 2};
        Line(2) = {2, 3};
        Line(3) = {3, 4};
        Line(4) = {4, 1};
        Curve Loop(1) = {1, 2, 3, 4};
        Plane Surface(1) = {1};
        Point(5) = {0.25, 0.25, 0, 0.5};
        Point(6) = {0.75, 0.25, 0, 0.5};
        Point(7) = {0.25, 0.75, 0, 0.5};
        Point(8) = {0.75, 0.75, 0, 0.5};
        Point(9) = {0.5, 0.25, 0, 0.5};
        Line(5) = {5, 6};
        Line(6) = {7, 8};
        Point{9} In Surface{1};
        Line{5, 6} In Surface{1};
        Periodic Curve {6} = {5} Translate {0, 0.5, 0};
        """
end

function _periodic_geo_curve_graph(mode::Symbol;cycle::Bool=false)
    leaf_master=mode==:branch ? 30 : mode in (:chain,:expressions) ? 20 : throw(
        ArgumentError("unknown periodic graph mode $mode"))
    leaf_offset=mode==:branch ? 0.6 : 0.3
    cycle_statement=cycle ?
        "Periodic Curve {30} = {10} Translate {0, -0.6, 0};" : ""
    periodic_statements=if mode==:expressions
        """
        tagFraction = 9 / 10;
        slaveBegin = Sqrt(100) + tagFraction;
        slaveEnd = 4 * 5 + tagFraction;
        masterBegin = slaveEnd;
        masterEnd = 3 * 10 + tagFraction;
        curveStep = 20 / 2;
        verticalShift = 3 / 10;
        exactZero = Atan2(0, 1);
        Periodic Curve {slaveBegin:slaveEnd:curveStep} =
          {masterBegin:masterEnd:curveStep}
          Translate {exactZero, verticalShift, Sin(exactZero)};
        """
    else
        """
        Periodic Curve {20} = {30} Translate {0, 0.3, 0};
        Periodic Curve {10} = {$leaf_master} Translate {0, $leaf_offset, 0};
        """
    end
    return """
        Point(1) = {0, 0, 0, 0.4};
        Point(2) = {1, 0, 0, 0.4};
        Point(3) = {1, 1, 0, 0.4};
        Point(4) = {0, 1, 0, 0.4};
        Line(1) = {1, 2};
        Line(2) = {2, 3};
        Line(3) = {3, 4};
        Line(4) = {4, 1};
        Curve Loop(1) = {1, 2, 3, 4};
        Plane Surface(1) = {1};
        Point(101) = {0.2, 0.2, 0, 0.4};
        Point(102) = {0.8, 0.2, 0, 0.4};
        Point(103) = {0.2, 0.5, 0, 0.4};
        Point(104) = {0.8, 0.5, 0, 0.4};
        Point(105) = {0.2, 0.8, 0, 0.4};
        Point(106) = {0.8, 0.8, 0, 0.4};
        Point(107) = {0.425, 0.8, 0, 0.4};
        Line(30) = {101, 102};
        Line(20) = {103, 104};
        Line(10) = {106, 105};
        Point{107} In Surface{1};
        Line{30, 20, 10} In Surface{1};
        $periodic_statements
        $cycle_statement
        """
end

const _PERIODIC_SURFACE_VOLUME_GEO=normpath(joinpath(
    @__DIR__,"..","fixtures","periodic_surface_volume.geo"))

# Two unit squares stacked in z: surface 1 at z=0 with curves 1–4, surface 2
# at z=1 with curves 5–8 — the `Periodic Surface slave {curves} = master
# {curves}` edge-map fixture.
function _periodic_geo_two_squares(periodic_statement::AbstractString;
                                   mesh_size=0.5)
    return """
        Point(1) = {0, 0, 0, $mesh_size}; Point(2) = {1, 0, 0, $mesh_size};
        Point(3) = {1, 1, 0, $mesh_size}; Point(4) = {0, 1, 0, $mesh_size};
        Point(5) = {0, 0, 1, $mesh_size}; Point(6) = {1, 0, 1, $mesh_size};
        Point(7) = {1, 1, 1, $mesh_size}; Point(8) = {0, 1, 1, $mesh_size};
        Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
        Line(5) = {5, 6}; Line(6) = {6, 7}; Line(7) = {7, 8}; Line(8) = {8, 5};
        Curve Loop(1) = {1, 2, 3, 4};
        Curve Loop(2) = {5, 6, 7, 8};
        Plane Surface(1) = {1};
        Plane Surface(2) = {2};
        $periodic_statement
        """
end

@testset "bounded .geo transform-free periodic curve execution" begin
    # `Periodic Curve {slave} = {master}` — upstream's orientation-only
    # relation: no affine is stored, signed tags carry the reversal and no
    # geometric endpoint check runs (`GEdge::setMeshMaster(source, ori)`).
    for (statement,reversed) in (
            ("Periodic Curve {2} = {4};",false),
            ("Periodic Curve {2} = {-4};",true),
            ("Periodic Curve {-2} = {4};",true),
            ("Periodic Curve {-2} = {-4};",false),
            ("Periodic Line {2} = {4};",false))
        built=_execute_geo_source(_periodic_geo_square(statement))
        constraint=only(model_periodic_constraints(built.model))
        @test constraint.dim==1
        @test (constraint.slave_entity,constraint.master_entity)==(2,4)
        @test constraint.affine===nothing
        @test constraint.reversed==reversed
    end

    # Orientation-only ignores geometry entirely — a longer slave curve is
    # accepted where an affine relation would drop on its endpoint check.
    mismatched=_execute_geo_source("""
        Point(1) = {0, 0, 0, 0.4};
        Point(2) = {1, 0, 0, 0.4};
        Point(3) = {3, 0, 1, 0.4};
        Point(4) = {0, 0, 1, 0.4};
        Line(1) = {1, 2};
        Line(2) = {4, 3};
        Periodic Curve {2} = {1};
        """)
    @test only(model_periodic_constraints(mismatched.model)).affine===nothing
    dropped=_execute_geo_source("""
        Point(1) = {0, 0, 0, 0.4};
        Point(2) = {1, 0, 0, 0.4};
        Point(3) = {3, 0, 1, 0.4};
        Point(4) = {0, 0, 1, 0.4};
        Line(1) = {1, 2};
        Line(2) = {4, 3};
        Periodic Curve {2} = {1} Translate {0, 0, 1};
        """)
    @test dropped.msg_error_count==0
    @test isempty(model_periodic_constraints(dropped.model))

    # Meshing: the slave inherits the master's parameter distribution and the
    # pairing is parameter-based — both transforms and orientation signs map
    # onto the same synchronized nodes.
    meshed=_execute_geo_source(
        _periodic_geo_square("Periodic Curve {2} = {4};");mesh_dim=2)
    @test validate(meshed.mesh).ok
    mapping=model_periodic_nodes(meshed.model,meshed.mesh,1,2)
    @test mapping.master_entity==4
    @test mapping.affine===nothing
    @test !isempty(mapping.slave_nodes)
    slave_coordinates=[Tuple(meshed.mesh.coords[:,node])
                       for node in mapping.slave_nodes]
    master_coordinates=[Tuple(meshed.mesh.coords[:,node])
                        for node in mapping.master_nodes]
    # Forward pairing: equal parameters — curve 2 ascends (2→3) while curve 4
    # descends (4→1), so t corresponds to opposite endpoints.
    for (slave,master) in zip(slave_coordinates,master_coordinates)
        @test slave[1]≈1.0 && master[1]≈0.0
        @test slave[2]+master[2]≈1.0 atol=1e-12
    end
    reversed=_execute_geo_source(
        _periodic_geo_square("Periodic Curve {2} = {-4};");mesh_dim=2)
    reversed_mapping=model_periodic_nodes(
        reversed.model,reversed.mesh,1,2)
    for (slave,master) in zip(
            (Tuple(reversed.mesh.coords[:,node])
             for node in reversed_mapping.slave_nodes),
            (Tuple(reversed.mesh.coords[:,node])
             for node in reversed_mapping.master_nodes))
        @test slave[2]≈master[2] atol=1e-12
    end

    # The zero-affine record round-trips through MSH4 exactly like Gmsh's own
    # transform-free output — vertex links carry no affine either.
    mixed=model_to_mixed(meshed.model,meshed.mesh,1)
    link_kinds=sort!([(link.dim,link.affine===nothing)
                      for link in mixed.periodic_links])
    @test link_kinds==[(0,true),(0,true),(1,true)]
    path=tempname()*".msh"
    Tessella.Elements.write_mixed_msh(path,mixed)
    back=Tessella.Elements.read_mixed_msh(path)
    @test sort!([(link.dim,link.affine===nothing,length(link.slave_nodes))
                 for link in back.periodic_links])==
          sort!([(link.dim,link.affine===nothing,length(link.slave_nodes))
                 for link in mixed.periodic_links])
    rm(path;force=true)
end

@testset "bounded .geo periodic surface edge-mapping execution" begin
    # `Periodic Surface j {slave curves} = k {master curves}` derives the
    # transform from the mapped boundary vertices (a translation or a
    # rotation about the mean-plane intersection) and then takes the
    # ordinary affine surface path.
    translated=_execute_geo_source(_periodic_geo_two_squares(
        "Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3, 4};"))
    constraint=only(model_periodic_constraints(translated.model))
    @test constraint.dim==2
    @test (constraint.slave_entity,constraint.master_entity)==(2,1)
    @test constraint.affine==
          (1.0,0.0,0.0,0.0,
           0.0,1.0,0.0,0.0,
           0.0,0.0,1.0,1.0,
           0.0,0.0,0.0,1.0)

    # Signed curve entries flip the endpoint correspondence but still derive
    # the same translation.
    flipped=_execute_geo_source(_periodic_geo_two_squares(
        "Periodic Surface 2 {5, 6, 7, -8} = 1 {1, 2, 3, -4};"))
    @test only(model_periodic_constraints(flipped.model)).affine==
          constraint.affine

    # A rotation: the slave square is the master square turned 90° about the
    # x axis — the mean planes meet on that axis.
    rotated=_execute_geo_source("""
        Point(1) = {0, 0, 0, 0.5}; Point(2) = {1, 0, 0, 0.5};
        Point(3) = {1, 1, 0, 0.5}; Point(4) = {0, 1, 0, 0.5};
        Point(5) = {0, 0, 0, 0.5}; Point(6) = {1, 0, 0, 0.5};
        Point(7) = {1, 0, 1, 0.5}; Point(8) = {0, 0, 1, 0.5};
        Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
        Line(5) = {5, 6}; Line(6) = {6, 7}; Line(7) = {7, 8}; Line(8) = {8, 5};
        Curve Loop(1) = {1, 2, 3, 4};
        Curve Loop(2) = {5, 6, 7, 8};
        Plane Surface(1) = {1};
        Plane Surface(2) = {2};
        Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3, 4};
        """)
    affine=only(model_periodic_constraints(rotated.model)).affine
    @test affine!==nothing
    # The derived map must send every master vertex onto its slave image.
    model=rotated.model
    for (master_point,slave_point) in ((1,5),(2,6),(3,7),(4,8))
        expected=model.points[slave_point]
        @test all(isapprox.(collect(_geo_periodic_affine_point(
                affine,model.points[master_point]).-expected),0.0;atol=1e-12))
    end

    # A count mismatch is a `yymsg(0)` diagnostic — accumulated and thrown at
    # end of parse like upstream's nonzero exit.
    thrown=try
        _execute_geo_source(_periodic_geo_two_squares(
            "Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3};"));nothing
    catch err
        err
    end
    @test thrown isa ArgumentError
    @test occursin("Wrong number of surface curves in periodicity " *
                   "constraint",sprint(showerror,thrown))
    # `Msg::Error` diagnostics mark the run failed (msg_error_count) but the
    # parse completes and returns a model without the relation — upstream's
    # nonfatal `addPeriodicFace` contract. A matched-count declaration that
    # still leaves a boundary curve unmapped fails the same way.
    for (statement,needle) in (
            ("Periodic Surface 2 {5, 6, 7} = 1 {1, 2, 3};",
             "Could not find curve counterpart 8 in slave surface 1"),
            ("Periodic Surface 2 {5, 6, 7, 9} = 1 {1, 2, 3, 4};",
             "Could not find curve counterpart 8 in slave surface 1"),
            ("Periodic Surface 9 {5, 6, 7, 8} = 1 {1, 2, 3, 4};",
             "Could not find surface 9 or 1 for periodic copy"),
            ("Periodic Surface 2 {5, 6, 7, 8} = 9 {1, 2, 3, 4};",
             "Could not find surface 2 or 9 for periodic copy"))
        built=_execute_geo_source(_periodic_geo_two_squares(statement))
        @test built.msg_error_count>0
        @test isempty(model_periodic_constraints(built.model))
    end
    # A master counterpart that is not a curve remains a hard diagnostic —
    # upstream dereferences the missing edge outright here.
    @test_throws ArgumentError _execute_geo_source(_periodic_geo_two_squares(
        "Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3, 9};"))

    # `Periodic Surface 2 = 1` is a parse-level `syntax error (=)` upstream.
    bare=try
        _execute_geo_source(_periodic_geo_two_squares(
            "Periodic Surface 2 = 1;"));nothing
    catch err
        err
    end
    @test bare isa ArgumentError
    @test occursin("syntax error (=)",sprint(showerror,bare))
end

@testset "bounded .geo periodic affine arity" begin
    affine12="{1,0,0,0, 0,1,0,0, 0,0,1,1}"
    affine16="{1,0,0,0, 0,1,0,0, 0,0,1,1, 0,0,0,1}"
    curve_geo="""
        Point(1) = {0,0,0,0.2}; Point(2) = {1,0,0,0.2};
        Point(3) = {0,0,1,0.2}; Point(4) = {1,0,1,0.2};
        Line(1) = {1,2}; Line(2) = {3,4};
        """
    # `Affine{16}` on matching endpoints stores the affine relation.
    stored=_execute_geo_source(
        curve_geo*"Periodic Curve {2} = {1} Affine $affine16;")
    @test only(model_periodic_constraints(stored.model)).affine!==nothing
    # `Affine{}` records the orientation-only relation, like an empty
    # `PeriodicTransform`.
    empty_list=_execute_geo_source(
        curve_geo*"Periodic Curve {2} = {1} Affine {};")
    @test only(model_periodic_constraints(empty_list.model)).affine===nothing
    # Curves with 12–15 or >16 entries: the endpoint check runs on the first
    # twelve entries; matching geometry then reports `GEntity::setMeshMaster`'s
    # exact-16 `Msg::Error` for the edge and each endpoint vertex — the
    # relation drops, the parse still completes.
    for statement in (
            "Periodic Curve {2} = {1} Affine $affine12;",
            "Periodic Curve {2} = {1} Affine " *
                "{1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,9};")
        built=_execute_geo_source(curve_geo*statement)
        @test built.msg_error_count==3
        @test isempty(model_periodic_constraints(built.model))
    end
    # `Affine` with 1–11 entries is a `yymsg(0)` arity error — accumulated
    # and thrown at end of parse (upstream still runs the orientation path,
    # invisible behind the thrown diagnostic).
    thrown=try
        _execute_geo_source(
            curve_geo*"Periodic Curve {2} = {1} Affine {1,0,0,0};");nothing
    catch err
        err
    end
    @test thrown isa ArgumentError
    @test occursin("Affine transformation requires at least 12 entries " *
                   "(4 provided)",sprint(showerror,thrown))
    # A bare `Affine` keyword is a `syntax error (;)` upstream.
    thrown=try
        _execute_geo_source(
            curve_geo*"Periodic Curve {2} = {1} Affine;");nothing
    catch err
        err
    end
    @test thrown isa ArgumentError
    @test occursin("syntax error (;)",sprint(showerror,thrown))
    # Mismatched curve geometry is a silent `Msg::Info` drop upstream: the
    # parse succeeds, no relation records, nothing is marked failed.
    dropped=_execute_geo_source("""
        Point(1) = {0,0,0,0.2}; Point(2) = {1,0,0,0.2};
        Point(3) = {0,0,1,0.2}; Point(4) = {5,0,1,0.2};
        Line(1) = {1,2}; Line(2) = {3,4};
        Periodic Curve {2} = {1} Translate {0,0,1};
        """)
    @test dropped.msg_error_count==0
    @test isempty(model_periodic_constraints(dropped.model))
    # The strict public API still rejects the same geometry — `set_periodic!`
    # is a validator, the `.geo` path mirrors upstream's tolerant drop.
    strict_model=_execute_geo_source("""
        Point(1) = {0,0,0,0.2}; Point(2) = {1,0,0,0.2};
        Point(3) = {0,0,1,0.2}; Point(4) = {5,0,1,0.2};
        Line(1) = {1,2}; Line(2) = {3,4};
        """).model
    @test_throws ArgumentError set_periodic!(
        strict_model,1,[2],[1],
        (1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,1.0, 0.0,0.0,0.0,1.0))
    # Surfaces need at least 12 entries and zero-pad to 16 — only the first
    # twelve are ever applied (`SPoint3::transform` reads a 3×4 map).
    for statement in (
            "Periodic Surface {2} = {1} Affine $affine12;",
            "Periodic Surface {2} = {1} Affine $affine16;",
            "Periodic Surface {2} = {1} Affine " *
                "{1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,5,5,5,5};")
        built=_execute_geo_source(_periodic_geo_two_squares(statement))
        @test built.msg_error_count==0
        @test only(model_periodic_constraints(built.model)).affine!==nothing
    end
    # The parser's `transfo(16, 0)` copy stores the fourth row verbatim — a
    # 12-entry surface `Affine` keeps `(0, 0, 0, 0)`, matching the bytes
    # upstream writes to `$Periodic`.
    padded=_execute_geo_source(_periodic_geo_two_squares(
        "Periodic Surface {2} = {1} Affine $affine12;"))
    @test only(model_periodic_constraints(padded.model)).affine[13:16]==
          (0.0,0.0,0.0,0.0)
    for (statement,needle) in (
            ("Periodic Surface {2} = {1} Affine {1,0,0,0};","(4 provided)"),
            ("Periodic Surface {2} = {1};","(0 provided)"))
        thrown=try
            _execute_geo_source(_periodic_geo_two_squares(statement));nothing
        catch err
            err
        end
        @test thrown isa ArgumentError
        @test occursin("Affine transformation requires at least 12 entries " *
                       needle,sprint(showerror,thrown))
    end
end

@testset "bounded .geo periodic planar-surface execution" begin
    built=execute_geo(_PERIODIC_SURFACE_VOLUME_GEO)
    @test built.mesh===nothing
    @test [(constraint.dim,Int(constraint.slave_entity),
            Int(constraint.master_entity))
           for constraint in model_periodic_constraints(built.model)]==
          [(2,4,6),(2,5,3)]

    meshed=execute_geo(_PERIODIC_SURFACE_VOLUME_GEO;mesh_dim=3)
    @test meshed.mesh!==nothing
    @test validate(meshed.mesh).ok
    @test mesh_crc(meshed.mesh).sha==
          "a58374071a4c485a339e1c5b48b8b0f3e69bf362ff0e41f57ca1a665139e81df"
    @test length(model_periodic_nodes(
        meshed.model,meshed.mesh,2,4).slave_nodes)==25
    @test length(model_periodic_nodes(
        meshed.model,meshed.mesh,2,5).slave_nodes)==25

    source=read(_PERIODIC_SURFACE_VOLUME_GEO,String)
    rotated=_execute_geo_source(replace(
        source,
        "Periodic Surface {4} = {6} Translate {1, 0, 0};"=>
        "Periodic Surface {4} = {6} " *
        "Rotate {{0, 0, 1}, {0.5, 0.5, 0}, Pi};");mesh_dim=3)
    @test validate(rotated.mesh).ok
    @test mesh_crc(rotated.mesh).sha==mesh_crc(meshed.mesh).sha
    rotated_mapping=model_periodic_nodes(
        rotated.model,rotated.mesh,2,4)
    @test length(rotated_mapping.slave_nodes)==25
    for (slave,master) in zip(rotated_mapping.slave_nodes,
                              rotated_mapping.master_nodes)
        actual=Tuple(rotated.mesh.coords[:,slave])
        expected=_geo_periodic_affine_point(
            rotated_mapping.affine,Tuple(rotated.mesh.coords[:,master]))
        @test hypot((actual.-expected)...)<=1e-15
    end

    affine_source=replace(
        source,
        "Periodic Surface {4} = {6} Translate {1, 0, 0};"=>
        "Periodic Surface {4} = {6} Affine " *
        "{1,0,0,1, 0,1,0,0, 0,0,1,0};")
    affine_built=_execute_geo_source(affine_source)
    # The applied 3×4 map is identical — but upstream pads a 12-entry
    # `Affine` to row4 `(0,0,0,0)` while `Translate` derives `(0,0,0,1)`, so
    # the stored records legitimately differ in the fourth row.
    affine_constraint=first(model_periodic_constraints(affine_built.model))
    translate_constraint=first(model_periodic_constraints(meshed.model))
    @test affine_constraint.affine[1:12]==translate_constraint.affine[1:12]
    @test affine_constraint.affine[13:16]==(0.0,0.0,0.0,0.0)

    @test_throws ArgumentError _execute_geo_source(replace(
        source,"Periodic Surface {4}"=>"Periodic Volume {4}"))
    @test_throws ArgumentError _execute_geo_source(replace(
        source,"Periodic Surface {4}"=>"Periodic Point {4}"))
end

@testset "bounded .geo periodic straight-curve execution" begin
    translated=_execute_geo_source(
        _periodic_geo_square(
            "Periodic Line {2} = {4} Translate {1, 0, 0};");
        mesh_dim=2)
    @test translated.mesh!==nothing
    @test validate(translated.mesh).ok
    @test mesh_crc(translated.mesh).sha==
          "3511d556ca0894daa79152eaf56abc6961024a72fa4f7e94f3357a7aa3cf0ff5"
    translation_constraint=only(
        model_periodic_constraints(translated.model))
    @test translation_constraint.affine==
          (1.0,0.0,0.0,1.0,
           0.0,1.0,0.0,0.0,
           0.0,0.0,1.0,0.0,
           0.0,0.0,0.0,1.0)
    translation_mapping=model_periodic_nodes(
        translated.model,translated.mesh,1,2)
    @test length(translation_mapping.slave_nodes)==5
    for (slave,master) in zip(translation_mapping.slave_nodes,
                              translation_mapping.master_nodes)
        @test Tuple(translated.mesh.coords[:,slave])==
              (translated.mesh.coords[1,master]+1,
               translated.mesh.coords[2,master],
               translated.mesh.coords[3,master])
    end

    affine=_execute_geo_source(_periodic_geo_square("""
        PeriodicAffineOne = Cos(0);
        affineZero = Atan2(0, 1);
        affineDx = Sqrt(4) / 2;
        Periodic Curve {2} = {4} Affine
          {PeriodicAffineOne,affineZero,affineZero,affineDx,
           affineZero,PeriodicAffineOne,affineZero,affineZero,
           affineZero,affineZero,PeriodicAffineOne,affineZero,
           affineZero,affineZero,affineZero,PeriodicAffineOne};
        """))
    @test only(model_periodic_constraints(affine.model)).affine==
          translation_constraint.affine

    rotated=_execute_geo_source(_periodic_geo_rotation();mesh_dim=2)
    @test validate(rotated.mesh).ok
    @test mesh_crc(rotated.mesh).sha==
          "f6ad616e56d52d7e10a598a4079db2de9b3d5f2a777f492f5a2366946d8ea990"
    rotation_constraint=only(model_periodic_constraints(rotated.model))
    @test !rotation_constraint.reversed
    @test rotation_constraint.affine[2]≈-1.0 atol=1e-15
    @test rotation_constraint.affine[4]≈2.0 atol=1e-15
    @test rotation_constraint.affine[5]≈1.0 atol=1e-15
    rotation_mapping=model_periodic_nodes(rotated.model,rotated.mesh,1,3)
    @test length(rotation_mapping.slave_nodes)==3
    for (slave,master) in zip(rotation_mapping.slave_nodes,
                              rotation_mapping.master_nodes)
        @test Tuple(rotated.mesh.coords[:,slave])==
              _geo_periodic_affine_point(
                  rotation_constraint.affine,
                  Tuple(rotated.mesh.coords[:,master]))
    end

    embedded=_execute_geo_source(
        _periodic_geo_embedded_curves();mesh_dim=2)
    @test validate(embedded.mesh).ok
    @test mesh_crc(embedded.mesh).sha==
          "9794a65ea5402683d0d50612522c2f71f7c98ec2a9f6b9e6b49a61e62cd85cf2"
    embedded_mapping=model_periodic_nodes(
        embedded.model,embedded.mesh,1,6)
    @test embedded_mapping.master_entity==5
    @test length(embedded_mapping.slave_nodes)==3
    for (slave,master) in zip(embedded_mapping.slave_nodes,
                              embedded_mapping.master_nodes)
        @test Tuple(embedded.mesh.coords[:,slave])==
              (embedded.mesh.coords[1,master],
               embedded.mesh.coords[2,master]+0.5,
               embedded.mesh.coords[3,master])
    end

    fractional=_execute_geo_source(_periodic_geo_square(
        "Periodic Curve {29 / 10} = {41 / 10} " *
        "Translate {Cos(0), Sin(0), Atan2(0, 1)};");mesh_dim=2)
    fractional_constraint=only(model_periodic_constraints(fractional.model))
    @test (fractional_constraint.slave_entity,
           fractional_constraint.master_entity)==(2,4)
    @test mesh_crc(fractional.mesh).sha==
          "3511d556ca0894daa79152eaf56abc6961024a72fa4f7e94f3357a7aa3cf0ff5"

    for (mode,masters,offsets) in (
            (:branch,Dict(10=>30,20=>30),Dict(10=>0.6,20=>0.3)),
            (:chain,Dict(10=>20,20=>30),Dict(10=>0.3,20=>0.3)),
            (:expressions,Dict(10=>20,20=>30),Dict(10=>0.3,20=>0.3)))
        graph=_execute_geo_source(
            _periodic_geo_curve_graph(mode);mesh_dim=2)
        @test validate(graph.mesh).ok
        @test mesh_crc(graph.mesh).sha==
              "dad04f30f3b17630127c3f1b4f5b5a4776ae5ff20d3c89afa6c674fac24d5338"
        graph_constraints=model_periodic_constraints(graph.model)
        @test Int.(getproperty.(graph_constraints,:slave_entity))==[10,20]
        @test Int.(getproperty.(graph_constraints,:master_entity))==
              [masters[10],masters[20]]
        @test graph_constraints[1].reversed
        for slave_entity in (10,20)
            graph_mapping=model_periodic_nodes(
                graph.model,graph.mesh,1,slave_entity)
            @test graph_mapping.master_entity==masters[slave_entity]
            @test length(graph_mapping.slave_nodes)==9
            for (slave,master) in zip(graph_mapping.slave_nodes,
                                      graph_mapping.master_nodes)
                @test Tuple(graph.mesh.coords[:,slave])==
                      (graph.mesh.coords[1,master],
                       graph.mesh.coords[2,master]+offsets[slave_entity],
                       graph.mesh.coords[3,master])
            end
        end
    end
    # A consistent cyclic declaration stores like upstream — `setMeshMaster`
    # has no cycle check — and all three relations serialize. Upstream's
    # deferred mesh copy starves the cycle members (empty link node lists);
    # Tessella's parameter sync converges on the cycle instead, so the links
    # carry real node pairs — a strict superset of upstream's output.
    cyclic=_execute_geo_source(
        _periodic_geo_curve_graph(:chain;cycle=true);mesh_dim=2)
    @test cyclic.msg_error_count==0
    @test validate(cyclic.mesh).ok
    @test mesh_crc(cyclic.mesh).sha==
          "dad04f30f3b17630127c3f1b4f5b5a4776ae5ff20d3c89afa6c674fac24d5338"
    cyclic_constraints=model_periodic_constraints(cyclic.model)
    @test Int.(getproperty.(cyclic_constraints,:slave_entity))==[10,20,30]
    @test Int.(getproperty.(cyclic_constraints,:master_entity))==[20,30,10]
    for (slave_entity,master_entity,offset) in
            ((10,20,0.3),(20,30,0.3),(30,10,-0.6))
        cyclic_mapping=model_periodic_nodes(
            cyclic.model,cyclic.mesh,1,slave_entity)
        @test cyclic_mapping.master_entity==master_entity
        @test length(cyclic_mapping.slave_nodes)==9
        for (slave,master) in zip(cyclic_mapping.slave_nodes,
                                  cyclic_mapping.master_nodes)
            @test Tuple(cyclic.mesh.coords[:,slave])==
                  (cyclic.mesh.coords[1,master],
                   cyclic.mesh.coords[2,master]+offset,
                   cyclic.mesh.coords[3,master])
        end
    end
    # All three stored relations serialize — including the cycle-closing edge —
    # plus their six endpoint links, matching upstream's nine `$Periodic`
    # records (upstream's starved curve links carry empty node lists; the
    # converged sync here populates them).
    cyclic_mixed=model_to_mixed(cyclic.model,cyclic.mesh,1)
    @test length(cyclic_mixed.periodic_links)==9
    cyclic_counts=Dict{Int,Int}()
    for link in cyclic_mixed.periodic_links
        cyclic_counts[link.dim]=get(cyclic_counts,link.dim,0)+1
    end
    @test cyclic_counts==Dict(0=>6,1=>3)
    cyclic_curve_masters=Dict{Int,Int}(
        Int(link.slave_entity)=>Int(link.master_entity)
        for link in cyclic_mixed.periodic_links if link.dim==1)
    @test cyclic_curve_masters==Dict(10=>20,20=>30,30=>10)
    cyclic_point_masters=Dict{Int,Int}(
        Int(link.slave_entity)=>Int(link.master_entity)
        for link in cyclic_mixed.periodic_links if link.dim==0)
    @test cyclic_point_masters==
          Dict(106=>104,105=>103,103=>101,104=>102,102=>106,101=>105)

    # Hard failures — syntax aborts and accumulated `yymsg(0)` errors are
    # still thrown as `ArgumentError` at end of parse.
    invalid_statements=(
        "Periodic Curve {2} = {4} Translate {1,0};",
        "Periodic Curve {2} = {4} Affine {1,0,0,1};",
        "Periodic Curve {2} = {4} Rotate {{0,0,1},{0,0,0},missingAngle};",
        "Periodic Curve {2} = {4} Translate {NaN,0,0};",
        "Periodic Curve {2} = {4} Translate {1 / 0,0,0};",
        "Periodic Curve {unknownTag} = {4} Translate {1,0,0};",
        "Periodic Curve {2:4:0} = {4} Translate {1,0,0};",
        "Periodic Curve {1:65537} = {4} Translate {1,0,0};",
        "Periodic Curve {2} = {4} Translate {Atan2(0),0,0};",
        "Periodic Curve {2} = {4} Mirror {1,0,0};",
        "Periodic Curve {2,} = {4} Translate {1,0,0};",
    )
    for statement in invalid_statements
        @test_throws ArgumentError _execute_geo_source(
            _periodic_geo_square(statement))
    end
    # `Msg::Error` diagnostics mark the run failed but return the model —
    # upstream's nonfatal `addPeriodicEdge`/`addPeriodicFace` contract.
    for statement in (
            "Periodic Surface {1} = {1} Translate {1,0,0};",
            "dynamicTag = newl; Periodic Curve {dynamicTag} = {4} " *
                "Translate {1,0,0};")
        built=_execute_geo_source(_periodic_geo_square(statement))
        @test built.msg_error_count>0
        @test isempty(model_periodic_constraints(built.model))
    end
    # A 16-entry `Affine` with a junk homogeneous row is stored verbatim like
    # upstream's `GEntity::setMeshMaster` — only the first twelve entries are
    # applied, so the stored record and its `$Periodic` bytes match Gmsh's.
    junk_row=_execute_geo_source(_periodic_geo_square(
        "Periodic Curve {2} = {4} Affine " *
        "{1,0,0,1, 0,1,0,0, 0,0,1,0, 0,0,0,2};"))
    junk_constraint=only(model_periodic_constraints(junk_row.model))
    @test junk_constraint.affine[13:16]==(0.0,0.0,0.0,2.0)
    @test junk_constraint.reversed
    # A degenerate transform whose endpoints never match is a silent
    # `Msg::Info` drop upstream — no relation, no error.
    dropped=_execute_geo_source(_periodic_geo_square(
        "Periodic Curve {2} = {4} Rotate {{0,0,0},{0,0,0},1};"))
    @test dropped.msg_error_count==0
    @test isempty(model_periodic_constraints(dropped.model))
    @test_throws ArgumentError _execute_geo_source(_periodic_geo_square(
        "newl = 2; Periodic Curve {newl} = {4} Translate {1,0,0};"))
    pi_error=try
        _execute_geo_source(_periodic_geo_square(
            "Pi = 3; Periodic Curve {2} = {4} Translate {1,0,0};"))
        nothing
    catch err
        err
    end
    @test pi_error isa ArgumentError
    # `Pi` is `tPi` upstream — an affectation head rejects it at parse time.
    @test occursin("syntax error (Pi)",sprint(showerror,pi_error))
    @test isempty(Docs.undocumented_names(Tessella.GeoExec;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.GeoExec;recursive=true))
end

function _periodic_geo_cube(offset::Int,shift::Float64)
    return """
        Point($(offset+1)) = {$shift, 0, 0, 1}; Point($(offset+2)) = {$(shift+1), 0, 0, 1};
        Point($(offset+3)) = {$(shift+1), 1, 0, 1}; Point($(offset+4)) = {$shift, 1, 0, 1};
        Point($(offset+5)) = {$shift, 0, 1, 1}; Point($(offset+6)) = {$(shift+1), 0, 1, 1};
        Point($(offset+7)) = {$(shift+1), 1, 1, 1}; Point($(offset+8)) = {$shift, 1, 1, 1};
        Line($(offset+1)) = {$(offset+1), $(offset+2)}; Line($(offset+2)) = {$(offset+2), $(offset+3)};
        Line($(offset+3)) = {$(offset+3), $(offset+4)}; Line($(offset+4)) = {$(offset+4), $(offset+1)};
        Line($(offset+5)) = {$(offset+5), $(offset+6)}; Line($(offset+6)) = {$(offset+6), $(offset+7)};
        Line($(offset+7)) = {$(offset+7), $(offset+8)}; Line($(offset+8)) = {$(offset+8), $(offset+5)};
        Line($(offset+9)) = {$(offset+1), $(offset+5)}; Line($(offset+10)) = {$(offset+2), $(offset+6)};
        Line($(offset+11)) = {$(offset+3), $(offset+7)}; Line($(offset+12)) = {$(offset+4), $(offset+8)};
        Curve Loop($(offset+1)) = {$(offset+1), $(offset+2), $(offset+3), $(offset+4)};
        Curve Loop($(offset+2)) = {$(offset+5), $(offset+6), $(offset+7), $(offset+8)};
        Curve Loop($(offset+3)) = {$(offset+1), $(offset+10), -$(offset+5), -$(offset+9)};
        Curve Loop($(offset+4)) = {$(offset+2), $(offset+11), -$(offset+6), -$(offset+10)};
        Curve Loop($(offset+5)) = {$(offset+3), $(offset+12), -$(offset+7), -$(offset+11)};
        Curve Loop($(offset+6)) = {$(offset+4), $(offset+9), -$(offset+8), -$(offset+12)};
        Plane Surface($(offset+1)) = {$(offset+1)}; Plane Surface($(offset+2)) = {$(offset+2)};
        Plane Surface($(offset+3)) = {$(offset+3)}; Plane Surface($(offset+4)) = {$(offset+4)};
        Plane Surface($(offset+5)) = {$(offset+5)}; Plane Surface($(offset+6)) = {$(offset+6)};
        Surface Loop($(offset+7)) = {$(offset+1), $(offset+2), $(offset+3),
                               $(offset+4), $(offset+5), $(offset+6)};
        """
end

@testset "bounded .geo periodic volume storage" begin
    # Gmsh 4.15.2 has no `Periodic Volume` `.geo` production — the statement
    # is a parse-level `syntax error (Volume)`. Tessella retains volume
    # periodicity as a model-level extension through `set_periodic!`.
    source=_periodic_geo_cube(0,0.0)*"\n"*_periodic_geo_cube(100,2.0)*"""
        Volume(1) = {7}; Volume(2) = {107};
        """
    built=_execute_geo_source(source)
    @test isempty(model_periodic_constraints(built.model))
    set_periodic!(built.model,3,[2],[1],
        (1.0,0.0,0.0,2.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0))
    constraint=only(model_periodic_constraints(built.model))
    @test constraint.dim==3
    @test constraint.slave_entity==2 && constraint.master_entity==1
    @test constraint.affine[4]==2.0

    for statement in (
            "Periodic Volume {2} = {1} Translate {2, 0, 0};",
            "Periodic Volume {2} = {1} " *
                "Affine {1,0,0,2, 0,1,0,0, 0,0,1,0, 0,0,0,1};",
            "Periodic Volume {2} = {9} Translate {2, 0, 0};",
            "Periodic Volume {9} = {1} Translate {2, 0, 0};",
            "Periodic Volume {2} = {2} Translate {2, 0, 0};")
        error=try
            _execute_geo_source(source*statement)
            nothing
        catch err
            err
        end
        @test error isa ArgumentError
        @test occursin("syntax error (Volume)",sprint(showerror,error))
    end
    @test_throws ArgumentError _execute_geo_source(
        _periodic_geo_square("Periodic Volume {1} = {1} Translate {1,0,0};"))

    model_only=_execute_geo_source(source).model
    @test_throws ArgumentError set_periodic!(
        model_only,3,[2],[9],
        (1.0,0.0,0.0,2.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0))
    @test_throws ArgumentError set_periodic!(
        model_only,3,[9],[1],
        (1.0,0.0,0.0,2.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0))
    @test_throws ArgumentError set_periodic!(
        model_only,3,[2],[2],
        (1.0,0.0,0.0,2.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0))
    set_periodic!(model_only,3,[2],[1],
        (1.0,0.0,0.0,2.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0))
    # A cyclic relation stores like upstream — `setMeshMaster` has no cycle
    # check — and volume relations are mesh-inert, so nothing downstream
    # needs the acyclic guarantee.
    set_periodic!(model_only,3,[1],[2],
        (1.0,0.0,0.0,-2.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0))
    @test sort!(Int.(getproperty.(
        model_periodic_constraints(model_only),:slave_entity)))==[1,2]
    @test_throws ArgumentError set_periodic!(
        model_only,3,[1],[1],
        (1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0))
end

@testset "cyclic periodic surface dependency" begin
    # Upstream stores the cycle — `GFace::setMeshMaster` has no cycle check —
    # and its deferred mesh copy starves both members forever. Tessella's
    # slave-copy path is recursive, so the cycle is detected in the ancestry
    # walk and reported explicitly instead of recursing.
    cyclic=_execute_geo_source(_periodic_geo_two_squares("""
        Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3, 4};
        Periodic Surface 1 {1, 2, 3, 4} = 2 {5, 6, 7, 8};
        """))
    @test cyclic.msg_error_count==0
    surface_masters=Dict{Int,Int}(
        Int(constraint.slave_entity)=>Int(constraint.master_entity)
        for constraint in model_periodic_constraints(cyclic.model)
        if constraint.dim==2)
    @test surface_masters==Dict(1=>2,2=>1)
    cycle_error=try
        mesh_model_surface(cyclic.model,1)
        nothing
    catch err
        err
    end
    @test cycle_error isa ArgumentError
    @test occursin(
        "cyclic periodic dependency Surface[1] -> Surface[2] -> Surface[1]",
        sprint(showerror,cycle_error))
    other_error=try
        mesh_model_surface(cyclic.model,2)
        nothing
    catch err
        err
    end
    @test other_error isa ArgumentError
    @test occursin(
        "cyclic periodic dependency Surface[2] -> Surface[1] -> Surface[2]",
        sprint(showerror,other_error))
end

@testset "periodic surface induced-edge resolution" begin
    # `GFace::setMeshMaster(master, tfo)` resolves every slave boundary and
    # embedded edge to a master counterpart at declaration, disambiguating
    # shared endpoint signatures through the transformed curve midpoint —
    # the pairs surface as the induced `$Periodic` curve links.
    built=_execute_geo_source(_periodic_geo_two_squares(
        "Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3, 4};"))
    @test built.msg_error_count==0
    constraint=only(model_periodic_constraints(built.model))
    @test constraint.dim==2
    _,pairs=Tessella.Model._model_periodic_surface_edge_map(
        built.model,2,1,constraint.affine,constraint.atol,"test")
    @test Dict(pairs)==Dict(5=>1,6=>2,7=>3,8=>4)

    # A slave surface meshes as a verbatim copy of its master — upstream's
    # meshGFace path meshes the master on demand, so the copy is exact even
    # when the master was never meshed directly.
    master_mesh=mesh_model_surface(built.model,1)
    slave_mesh=mesh_model_surface(built.model,2)
    @test slave_mesh.tris==master_mesh.tris
    @test slave_mesh.coords==master_mesh.coords .+ [0.0,0.0,1.0]

    # Embedded edges join the induced map: a matching embedded pair resolves
    # like a boundary edge, and an unbalanced embedding aborts the whole
    # relation with upstream's `Msg::Error` (the run is marked failed but
    # returns a model).
    matched=_execute_geo_source(_periodic_geo_two_squares("""
        Point(9) = {0.25, 0.5, 0, 0.5}; Point(10) = {0.75, 0.5, 0, 0.5};
        Point(11) = {0.25, 0.5, 1, 0.5}; Point(12) = {0.75, 0.5, 1, 0.5};
        Line(9) = {9, 10}; Line(10) = {11, 12};
        Line{9} In Surface{1}; Line{10} In Surface{2};
        Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3, 4};
        """))
    @test matched.msg_error_count==0
    constraint=only(model_periodic_constraints(matched.model))
    _,pairs=Tessella.Model._model_periodic_surface_edge_map(
        matched.model,2,1,constraint.affine,constraint.atol,"test")
    @test Dict(pairs)[10]==9

    unbalanced=_execute_geo_source(_periodic_geo_two_squares("""
        Point(9) = {0.25, 0.5, 1, 0.5}; Point(10) = {0.75, 0.5, 1, 0.5};
        Line(9) = {9, 10};
        Line{9} In Surface{2};
        Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3, 4};
        """))
    @test unbalanced.msg_error_count>0
    @test isempty(model_periodic_constraints(unbalanced.model))

    # Two embedded curves sharing one directed endpoint pair are ambiguous —
    # the transformed candidate midpoint decides (a straight chord and an
    # arc bulging away from the chord line resolve to their own kinds;
    # upstream's sampled bounding box fallback cannot separate them, so the
    # midpoint must hit first).
    ambiguous=_execute_geo_source("""
        Point(1) = {0, 0, 0, 0.5}; Point(2) = {1, 0, 0, 0.5};
        Point(3) = {1, 1, 0, 0.5}; Point(4) = {0, 1, 0, 0.5};
        Point(5) = {0, 0, 1, 0.5}; Point(6) = {1, 0, 1, 0.5};
        Point(7) = {1, 1, 1, 0.5}; Point(8) = {0, 1, 1, 0.5};
        Point(9) = {0, 0.5, 0, 0.5}; Point(10) = {1, 0.5, 0, 0.5};
        Point(11) = {0, 0.5, 1, 0.5}; Point(12) = {1, 0.5, 1, 0.5};
        Point(13) = {0.5, 0.7, 0, 0.5}; Point(14) = {0.5, 0.7, 1, 0.5};
        Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4};
        Line(4) = {4, 1};
        Line(5) = {5, 6}; Line(6) = {6, 7}; Line(7) = {7, 8};
        Line(8) = {8, 5};
        Line(9) = {9, 10}; Circle(10) = {9, 13, 10};
        Line(11) = {11, 12}; Circle(12) = {11, 14, 12};
        Curve Loop(1) = {1, 2, 3, 4}; Curve Loop(2) = {5, 6, 7, 8};
        Plane Surface(1) = {1}; Plane Surface(2) = {2};
        Line{9, 10} In Surface{1}; Line{11, 12} In Surface{2};
        Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3, 4};
        """)
    @test ambiguous.msg_error_count==0
    constraint=only(model_periodic_constraints(ambiguous.model))
    _,pairs=Tessella.Model._model_periodic_surface_edge_map(
        ambiguous.model,2,1,constraint.affine,constraint.atol,"test")
    map=Dict(pairs)
    @test map[11]==9 && map[12]==10

    # A rotation-derived map resolves the same way: surface 2 is surface 1
    # turned a quarter turn about the x axis, so the arc midpoint still
    # identifies its counterpart.
    rotated=_execute_geo_source("""
        Point(1) = {0, 0, 0, 0.5}; Point(2) = {1, 0, 0, 0.5};
        Point(3) = {1, 1, 0, 0.5}; Point(4) = {0, 1, 0, 0.5};
        Point(5) = {0, 0, 0, 0.5}; Point(6) = {1, 0, 0, 0.5};
        Point(7) = {1, 0, 1, 0.5}; Point(8) = {0, 0, 1, 0.5};
        Point(9) = {0, 0.5, 0, 0.5}; Point(10) = {1, 0.5, 0, 0.5};
        Point(11) = {0, 0, 0.5, 0.5}; Point(12) = {1, 0, 0.5, 0.5};
        Point(13) = {0.5, 0.7, 0, 0.5}; Point(14) = {0.5, 0, 0.7, 0.5};
        Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4};
        Line(4) = {4, 1};
        Line(5) = {5, 6}; Line(6) = {6, 7}; Line(7) = {7, 8};
        Line(8) = {8, 5};
        Line(9) = {9, 10}; Circle(10) = {9, 13, 10};
        Line(11) = {11, 12}; Circle(12) = {11, 14, 12};
        Curve Loop(1) = {1, 2, 3, 4}; Curve Loop(2) = {5, 6, 7, 8};
        Plane Surface(1) = {1}; Plane Surface(2) = {2};
        Line{9, 10} In Surface{1}; Line{11, 12} In Surface{2};
        Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3, 4};
        """)
    @test rotated.msg_error_count==0
    constraint=only(model_periodic_constraints(rotated.model))
    _,pairs=Tessella.Model._model_periodic_surface_edge_map(
        rotated.model,2,1,constraint.affine,constraint.atol,"test")
    map=Dict(pairs)
    @test map[11]==9 && map[12]==10
end

@testset "multi-surface periodic projection" begin
    # Two disjoint unit squares linked by a `Periodic Surface` edge map —
    # upstream writes nine links for this fixture (4 induced point pairs,
    # 4 induced curve pairs, and the surface map).
    function disjoint_pair(extra::AbstractString="")
        return """
            Point(1) = {0, 0, 0, 0.5};
            Point(2) = {1, 0, 0, 0.5};
            Point(3) = {1, 1, 0, 0.5};
            Point(4) = {0, 1, 0, 0.5};
            Line(1) = {1, 2}; Line(2) = {2, 3};
            Line(3) = {3, 4}; Line(4) = {4, 1};
            Curve Loop(1) = {1, 2, 3, 4};
            Plane Surface(1) = {1};
            Point(5) = {2, 0, 0, 0.5};
            Point(6) = {3, 0, 0, 0.5};
            Point(7) = {3, 1, 0, 0.5};
            Point(8) = {2, 1, 0, 0.5};
            Line(5) = {5, 6}; Line(6) = {6, 7};
            Line(7) = {7, 8}; Line(8) = {8, 5};
            Curve Loop(2) = {5, 6, 7, 8};
            Plane Surface(2) = {2};
            Periodic Surface 2 {5, 6, 7, 8} = 1 {1, 2, 3, 4};
            $extra
            """
    end
    execution=_execute_geo_source(disjoint_pair())
    @test execution.msg_error_count==0
    m=execution.model
    parts=[(tag,mesh_model_surface(m,tag)) for tag in (1,2)]
    mixed=model_to_mixed(m,parts)
    @test size(mixed.coords,2)==nnodes(parts[1][2])+nnodes(parts[2][2])
    @test length(mixed.periodic_links)==9
    counts=Dict{Int,Int}()
    for link in mixed.periodic_links
        counts[link.dim]=get(counts,link.dim,0)+1
        @test link.affine!==nothing && length(link.affine)==16
    end
    @test counts==Dict(0=>4,1=>4,2=>1)
    surface_link=only(link for link in mixed.periodic_links if link.dim==2)
    @test Int(surface_link.slave_entity)==2 &&
          Int(surface_link.master_entity)==1
    @test length(surface_link.slave_nodes)==nnodes(parts[2][2])
    @test length(Set(surface_link.slave_nodes))==nnodes(parts[2][2])
    curve_masters=Dict{Int,Int}(Int(link.slave_entity)=>Int(link.master_entity)
        for link in mixed.periodic_links if link.dim==1)
    @test curve_masters==Dict(5=>1,6=>2,7=>3,8=>4)
    point_masters=Dict{Int,Int}(Int(link.slave_entity)=>Int(link.master_entity)
        for link in mixed.periodic_links if link.dim==0)
    @test point_masters==Dict(5=>1,6=>2,7=>3,8=>4)

    # The serialized MSH4 reloads with the same nine links.
    roundtrip=mktemp() do path,io
        close(io)
        Tessella.Elements.write_mixed_msh(path,mixed)
        Tessella.Elements.read_mixed_msh(path)
    end
    @test length(roundtrip.periodic_links)==9
    @test Set((link.dim,Int(link.slave_entity),Int(link.master_entity))
              for link in roundtrip.periodic_links) ==
          Set((link.dim,Int(link.slave_entity),Int(link.master_entity))
              for link in mixed.periodic_links)

    # An explicit curve relation duplicating an induced pair folds into the
    # same link — one master per slave is preserved.
    duplicated=_execute_geo_source(disjoint_pair(
        "Periodic Curve {5} = {1} Translate {2, 0, 0};"))
    @test duplicated.msg_error_count==0
    dup_parts=[(tag,mesh_model_surface(duplicated.model,tag)) for tag in (1,2)]
    dup_mixed=model_to_mixed(duplicated.model,dup_parts)
    @test length(dup_mixed.periodic_links)==9

    # Rejections: partial coverage, duplicates, unknown or mis-shaped entries.
    @test_throws ArgumentError model_to_mixed(m,Tuple{Int,Mesh}[])
    @test_throws ArgumentError model_to_mixed(m,[parts[1]])
    @test_throws ArgumentError model_to_mixed(m,[parts[1],parts[1]])
    @test_throws ArgumentError model_to_mixed(m,[(99,parts[1][2])])
    @test_throws ArgumentError model_to_mixed(
        m,[(1,parts[1][2]),(2,parts[1][2])])
    @test_throws ArgumentError model_to_mixed(
        m,[(1,"not a mesh"),(2,parts[2][2])])

    # Adjacent coplanar squares share a boundary curve — the merge unifies
    # the shared nodes and emits each shared cell once.
    adjacent=_execute_geo_source("""
        Point(1) = {0, 0, 0, 0.5};
        Point(2) = {1, 0, 0, 0.5};
        Point(3) = {2, 0, 0, 0.5};
        Point(4) = {2, 1, 0, 0.5};
        Point(5) = {1, 1, 0, 0.5};
        Point(6) = {0, 1, 0, 0.5};
        Line(1) = {1, 2}; Line(2) = {2, 5}; Line(3) = {5, 6};
        Line(4) = {6, 1}; Line(5) = {2, 3}; Line(6) = {3, 4};
        Line(7) = {4, 5};
        Curve Loop(1) = {1, 2, 3, 4};
        Curve Loop(2) = {5, 6, 7, -2};
        Plane Surface(1) = {1};
        Plane Surface(2) = {2};
        """)
    am=adjacent.model
    aparts=[(tag,mesh_model_surface(am,tag)) for tag in (1,2)]
    amixed=model_to_mixed(am,aparts)
    @test isempty(amixed.periodic_links)
    @test size(amixed.coords,2)<sum(nnodes(p[2]) for p in aparts)
    tri_count=sum(block.msh==2 ? size(block.nodes,2) : 0
                  for block in amixed.blocks)
    @test tri_count==sum(ntris(p[2]) for p in aparts)
    point_count=sum(block.msh==15 ? size(block.nodes,2) : 0
                    for block in amixed.blocks)
    @test point_count==6
    seen_cells=Set{Tuple{Int,NTuple{4,Int32}}}()
    for block in amixed.blocks,column in axes(block.nodes,2)
        nodes=sort!(vec(Int32.(block.nodes[:,column])))
        key=(block.msh,ntuple(
            slot->slot<=length(nodes) ? nodes[slot] : Int32(0),4))
        @test key ∉ seen_cells
        push!(seen_cells,key)
    end
    # Shared curve 2 keeps curve-level ownership of its interior nodes in the
    # merged classification.
    data=amixed.entity_data
    @test any(owner->owner==(1,Int32(2)),data.node_entities)
end
