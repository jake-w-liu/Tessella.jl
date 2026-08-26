using Test
using Tessella
using Tessella.MeshTypes: mesh_crc, validate

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
        Periodic Curve {3} = {1}
          Rotate {{0, 0, 2}, {1, 1, 0}, 1.5707963267948966};
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

    affine=_execute_geo_source(_periodic_geo_square(
        "Periodic Curve {2} = {4} Affine " *
        "{1,0,0,1, 0,1,0,0, 0,0,1,0};"))
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

    invalid_statements=(
        "Periodic Surface {1} = {1} Translate {1,0,0};",
        "Periodic Curve {2} = {4} Translate {1,0};",
        "Periodic Curve {2} = {4} Affine {1,0,0,1};",
        "Periodic Curve {2} = {4} Rotate {{0,0,0},{0,0,0},1};",
        "Periodic Curve {2} = {4} Rotate {{0,0,1},{0,0,0},Pi/2};",
        "Periodic Curve {2} = {4} Translate {NaN,0,0};",
        "Periodic Curve {2} = {4} Mirror {1,0,0};",
        "Periodic Curve {2,} = {4} Translate {1,0,0};",
    )
    for statement in invalid_statements
        @test_throws ArgumentError _execute_geo_source(
            _periodic_geo_square(statement))
    end
    @test isempty(Docs.undocumented_names(Tessella.GeoExec;private=false))
    @test isempty(Test.detect_ambiguities(Tessella.GeoExec;recursive=true))
end
