# A STEP product assembly imported through the OpenCASCADE route places one
# solid per usage path, as Gmsh 4.15.2's STEP reader does. Tessella imported
# each part once at its own origin: examples/boolean/as1-tu-203.stp
# (step_assembly.geo) gave 5 volumes where Gmsh has 18 (84/126/53/5 entities
# against 236/354/160/18), so its -3 mesh had 8071 nodes instead of ~10700.
# Here the box fixture becomes a part used twice: translated by (5,0,0), and
# turned a quarter about z then moved to (0,3,1). Gmsh: 16/24/12/2 entities,
# bounding boxes [5,0,0,6,1,1] and [-1,3,1,0,4,2], names
# "Shapes/asm/BOX::1/<part>" and "Shapes/asm/BOX::2/<part>".
module ApiOccStepAssemblyTest
using Test, Tessella

const A = Tessella.API

@testset "OCC STEP assembly instances (Gmsh 4.15.2)" begin
    box = read(joinpath(@__DIR__, "..", "fixtures", "brep_occ_box.step"), String)
    extra = """
        #9001 = PRODUCT('asm','asm','',(#8));
        #9002 = PRODUCT_DEFINITION_FORMATION('','',#9001);
        #9003 = PRODUCT_DEFINITION('design','',#9002,#9);
        #9010 = CARTESIAN_POINT('',(0.,0.,0.));
        #9011 = DIRECTION('',(0.,0.,1.));
        #9012 = DIRECTION('',(1.,0.,0.));
        #9013 = AXIS2_PLACEMENT_3D('',#9010,#9011,#9012);
        #9014 = CARTESIAN_POINT('',(5.,0.,0.));
        #9015 = AXIS2_PLACEMENT_3D('',#9014,#9011,#9012);
        #9016 = CARTESIAN_POINT('',(0.,3.,1.));
        #9018 = DIRECTION('',(0.,1.,0.));
        #9017 = AXIS2_PLACEMENT_3D('',#9016,#9011,#9018);
        #9020 = SHAPE_REPRESENTATION('',(#9013,#9015,#9017),#345);
        #9030 = SHAPE_DEFINITION_REPRESENTATION(#9031,#9020);
        #9031 = PRODUCT_DEFINITION_SHAPE('','',#9003);
        #9021 = NEXT_ASSEMBLY_USAGE_OCCURRENCE('BOX::1','','BOX::1',#9003,#5,\$);
        #9022 = PRODUCT_DEFINITION_SHAPE('','',#9021);
        #9023 = ITEM_DEFINED_TRANSFORMATION('','',#11,#9015);
        #9024 = ( REPRESENTATION_RELATIONSHIP('','',#10,#9020) REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#9023) SHAPE_REPRESENTATION_RELATIONSHIP() );
        #9025 = CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#9024,#9022);
        #9041 = NEXT_ASSEMBLY_USAGE_OCCURRENCE('BOX::2','','BOX::2',#9003,#5,\$);
        #9042 = PRODUCT_DEFINITION_SHAPE('','',#9041);
        #9043 = ITEM_DEFINED_TRANSFORMATION('','',#11,#9017);
        #9044 = ( REPRESENTATION_RELATIONSHIP('','',#10,#9020) REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#9043) SHAPE_REPRESENTATION_RELATIONSHIP() );
        #9045 = CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#9044,#9042);
        """
    i = findlast("ENDSEC;", box)
    path = tempname() * ".step"
    write(path, box[1:first(i)-1] * extra * box[first(i):end])
    A.initialize()
    try
        A.option("General.Terminal", 0)
        A.merge(path)
        m = A.model
        @test [length(m.get_entities(d)) for d in 0:3] == [16, 24, 12, 2]
        @test collect(m.get_bounding_box(3, 1)) ≈ [5, 0, 0, 6, 1, 1] atol=1e-6
        @test collect(m.get_bounding_box(3, 2)) ≈ [-1, 3, 1, 0, 4, 2] atol=1e-6
        @test startswith(m.get_entity_name(3, 1), "Shapes/asm/BOX::1/")
        @test startswith(m.get_entity_name(3, 2), "Shapes/asm/BOX::2/")
        @test m.occ.get_mass(3, 2) ≈ 1 rtol=1e-12
    finally
        A.finalize(); rm(path; force=true)
    end
end

end
