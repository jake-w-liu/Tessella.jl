@testset "STEP mapped occurrences retain real placements" begin
    mesh=import_step(joinpath(@__DIR__,"..","fixtures","brep_mapped_boxes.step");fill=false)
    @test mesh.entity_data.entities[(3,1)].bbox==(3.,0.,0.,4.,1.,1.) &&
        mesh.entity_data.entities[(3,2)].bbox==(5.,0.,0.,6.,1.,1.)
    @test validate(mesh).ok && size(mesh.coords,2)==16
    # A nonzero Cartesian-operator scale is admitted: the mapped unit box
    # becomes the 2-edge box at the operator origin. Zero scale rejects.
    source=read(joinpath(@__DIR__,"..","fixtures","brep_mapped_boxes.step"),String)
    target="#905 = AXIS2_PLACEMENT_3D('',#902,#13,#904);"
    mktempdir() do directory
        path=joinpath(directory,"scaled.step")
        write(path,replace(source,target=>
            "#905 = CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',\$,\$,\$,#902,2.,\$);"))
        scaled=import_step(path;fill=false)
        @test validate(scaled).ok &&
            scaled.entity_data.entities[(3,2)].bbox==(6.,0.,0.,8.,2.,2.)
        write(path,replace(source,target=>
            "#905 = CARTESIAN_TRANSFORMATION_OPERATOR_3D('','',\$,\$,\$,#902,0.,\$);"))
        @test_throws ArgumentError import_step(path;fill=false)
    end
end
