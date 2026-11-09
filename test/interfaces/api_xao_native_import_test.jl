module XAONativeImportChecks
using Test,Tessella
const A=Tessella.API
const FIXTURES=joinpath(@__DIR__,"..","fixtures")
const XAO=joinpath(FIXTURES,"group_box.xao")

function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        f()
    finally
        A.finalize()
    end
end

@testset "XAO native topology, names and group accumulation" begin
    session() do
        @test A.open(XAO)===nothing
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(8,12,6,1)
        @test A.model.occ.get_mass(3,1)≈1.
        @test A.model.get_entities_for_physical_group(2,1)==[6]
        @test A.model.get_entities_for_physical_group(3,2)==[1]
        @test A.model.get_physical_name(2,1)=="top"
        @test A.merge(XAO)===nothing
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(16,24,12,2)
        @test A.model.get_entities_for_physical_group(2,1)==[6,12]
        @test A.model.get_entities_for_physical_group(3,2)==[1,2]
    end
    session() do
        A.open(XAO)
        expected=A.model.occ.get_bounding_box(2,1)
        text=read(XAO,String)
        text=replace(text,"face index=\"0\""=>"face index=\"5\"",
                          "face index=\"5\""=>"face index=\"0\"")
        text=replace(text,"face index=\"5\" name=\"\""=>"face index=\"5\" name=\"left face\"")
        mktempdir() do directory
            path=joinpath(directory,"permuted.xao")
            write(path,text)
            A.open(path)
            @test A.model.get_entity_name(2,6)=="left face"
            @test collect(A.model.occ.get_bounding_box(2,6))≈collect(expected)
            @test A.model.get_entities_for_physical_group(2,1)==[6]
            @test A.model.occ.get_mass(3,1)≈1.
        end
    end
end

@testset "XAO external shape uses its relative filename" begin
    session() do
        text=read(XAO,String)
        shape=String(match(r"(?s)<!\[CDATA\[(.*?)\]\]>",text).captures[1])
        mktempdir() do directory
            write(joinpath(directory,"source body.brep"),shape)
            external=replace(text,r"(?s)<shape\b.*?</shape>"=>
                "<shape format=\"BREP\" file=\"source body.brep\"/>")
            path=joinpath(directory,"external.xao")
            write(path,external)
            metadata=Tessella.IO.read_xao(path)
            @test metadata.format=="BREP" && metadata.brep===nothing
            @test A.open(path)===nothing
            @test A.model.occ.get_mass(3,1)≈1.
            @test A.model.get_entities_for_physical_group(2,1)==[6]
        end
    end
end

@testset "XAO merge retains actual mesh, fields, views and failed-import state" begin
    session() do
        A.model.add("retained")
        A.model.occ.add_box(3.,0.,0.,1.,1.,1.)
        A.model.occ.synchronize()
        A.option("Mesh.MeshSizeMin",2.)
        A.option("Mesh.MeshSizeMax",2.)
        A.mesh.generate(1)
        tags,coordinates,_=A.mesh.get_nodes()
        elements=A.mesh.get_elements()
        slot=A.MODEL_SLOTS[A._current_slot_index_locked()]
        field=A.mesh.field.add("MathEval")
        A.mesh.field.set_string(field,"F","2")
        A.mesh.field.set_as_background_mesh(field)
        view=A.view.add("held")
        A.view.add_model_data(view,0,A.model.get_current(),"NodeData",tags,
            [Float64[1.] for _ in tags])
        data=A.view.get_model_data(view,0)
        @test A.merge(XAO)===nothing
        @test A.MODEL_SLOTS[A._current_slot_index_locked()]===slot
        @test A.model.get_current()=="retained"
        @test A.model.get_file_name()==""
        @test A.mesh.get_nodes()[1:2]==(tags,coordinates)
        @test A.mesh.get_elements()==elements
        @test A.view.get_model_data(view,0)==data
        @test A.BACKGROUND_FIELD[]==field
        @test A.mesh.field.get_string(field,"F")=="2"
        @test A.model.get_entities_for_physical_group(3,2)==[2]
        before=(entities=A.model.get_entities(),nodes=A.mesh.get_nodes(),
            elements=A.mesh.get_elements(),groups=A.model.get_physical_groups(),
            views=A.view.get_tags(),data=A.view.get_model_data(view,0),
            field=A.mesh.field.get_string(field,"F"),name=A.model.get_current())
        mktempdir() do directory
            path=joinpath(directory,"invalid.xao")
            source=read(XAO,String)
            for broken in (replace(source,"reference=\"7\""=>"reference=\"999999\""),
                           replace(source,"<element index=\"5\""=>"<element index=\"999999\""))
                write(path,broken)
                @test_throws ArgumentError A.merge(path)
                @test A.MODEL_SLOTS[A._current_slot_index_locked()]===slot
                @test (entities=A.model.get_entities(),nodes=A.mesh.get_nodes(),
                    elements=A.mesh.get_elements(),groups=A.model.get_physical_groups(),
                    views=A.view.get_tags(),data=A.view.get_model_data(view,0),
                    field=A.mesh.field.get_string(field,"F"),name=A.model.get_current())==before
            end
        end
    end
end
end
