module CADFileImportChecks
using Test, Tessella
const A=Tessella.API
const M=Tessella.Model
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function native_model(path)
    model=GeoModel()
    M._model_occ_import_shapes(model,path)
    M._model_occ_synchronize!(model)
    return model
end

function state()
    return (models=A.model.list(),name=A.model.get_current(),file=A.model.get_file_name(),
        entities=A.model.get_entities(),nodes=A.mesh.get_nodes(),cells=A.mesh.get_elements(),
        views=A.view.get_tags(),fields=Dict(tag=>(spec.kind,copy(spec.options),copy(spec.order))
            for (tag,spec) in A.FIELD_SPECS[]),background=A.BACKGROUND_FIELD[],
        node_max=A.NODE_TAG_MAX[],cell_max=A.ELEMENT_TAG_MAX[],
        context=A.NATIVE_CONTEXT_LC[],forced=A.NATIVE_CONTEXT_FORCED_BBOX[])
end

function malformed_step(directory)
    source=read(joinpath(FIXTURES,"brep_occ_box.step"),String)
    position=findlast("ENDSEC;",source).start
    path=joinpath(directory,"partially valid.step")
    write(path,source[1:prevind(source,position)]*
        "\n#999999=MANIFOLD_SOLID_BREP('broken',#888888);\n"*source[position:end])
    return path
end

@testset "STEP and IGES open publish actual native CAD models" begin
    for name in ("brep_occ_box.step","brep_occ_box.igs")
        path=joinpath(FIXTURES,name)
        expected=native_model(path)
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.model.get_current()
            initial=A.MODEL_SLOTS[A._current_slot_index_locked()]
            @test A.open(path)===nothing
            @test A.MODEL_SLOTS[A._current_slot_index_locked()]===initial
            @test length(A.model.list())==1
            @test A.model.get_current()==splitext(name)[1]
            @test A.model.get_file_name()==path
            @test A.model.get_entities()==M.model_entities(expected)
            @test collect(A.model.get_bounding_box(-1,-1))≈collect(M.model_bounding_box(expected,-1,-1))
            @test isempty(first(A.mesh.get_nodes()))
            @test all(A.model.get_type(d,t)!="Discrete surface" for (d,t) in A.model.get_entities(2))
            if endswith(name,".step")
                A.option("Mesh.Algorithm",6)
                A.option("Mesh.MeshSizeMin",.75);A.option("Mesh.MeshSizeMax",.75)
                A.mesh.generate(1)
                @test !isempty(first(A.mesh.get_nodes()))
                @test !isempty(first(A.mesh.get_elements(1)))
            end
            curve=first(A.model.get_entities(1))[2]
            lower,upper=A.model.get_parametrization_bounds(1,curve)
            parameter=(only(lower)+only(upper))/2
            value=A.model.get_value(1,curve,[parameter])
            bounds=collect(A.model.get_bounding_box(-1,-1))
            roots=A.model.get_entities(3)
            isempty(roots) && (roots=A.model.get_entities(2))
            A.model.occ.translate(roots,2.,3.,4.)
            A.model.occ.synchronize()
            @test A.model.get_value(1,curve,[parameter])≈value+[2.,3.,4.]
            @test collect(A.model.get_bounding_box(-1,-1))≈bounds+[2.,3.,4.,2.,3.,4.] atol=1e-6
        finally
            A.finalize()
        end
    end
end

@testset "CAD merge retains the current mesh, fields and view owners" begin
    mktempdir() do directory
        broken=malformed_step(directory)
        for name in ("brep_occ_box.step","brep_occ_box.igs")
            A.initialize()
            try
                A.option("General.Terminal",0)
                A.model.add("kept CAD session");A.model.set_file_name("kept-associated.msh")
                A.model.add_discrete_entity(1,99)
                A.mesh.add_nodes(1,99,[501,502],[5.,0.,0.,6.,0.,0.])
                A.mesh.add_elements_by_type(99,1,[901],[501,502])
                held=A.view.add("retained source")
                A.view.add_model_data(held,0,"kept CAD session","NodeData",[501],[[7.25]])
                view=A.SESSION_VIEWS[held]
                slot=A.MODEL_SLOTS[A._current_slot_index_locked()]
                field=A.mesh.field.add("MathEval",701)
                A.mesh.field.set_string(field,"F","0.25")
                A.mesh.field.set_as_background_mesh(field)
                before=state();data=A.view.get_model_data(held,0)
                A.NATIVE_CONTEXT_FORCED_BBOX[]=true;A.NATIVE_CONTEXT_LC[]=123.
                @test A.merge(joinpath(FIXTURES,name))===nothing
                @test A.MODEL_SLOTS[A._current_slot_index_locked()]===slot
                @test A.model.get_current()==before.name && A.model.get_file_name()==before.file
                @test A.mesh.get_nodes()==before.nodes && A.mesh.get_elements()==before.cells
                @test A.SESSION_VIEWS[held]===view && A.view.get_model_data(held,0)==data
                @test all(context.model_slot===slot for context in values(view.model_contexts))
                @test A.mesh.field.get_string(field,"F")=="0.25" && A.BACKGROUND_FIELD[]==field
                @test !isempty(A.model.get_entities(2)) && (1,99) in A.model.get_entities(1)
                @test A.NATIVE_CONTEXT_FORCED_BBOX[] && A.NATIVE_CONTEXT_LC[]==123.
                before=state();model=A.CURRENT[]
                for endpoint in (A.open,A.merge)
                    @test_throws Exception endpoint(broken)
                    @test isequal(state(),before)
                    @test A.CURRENT[]===model && A.MODEL_SLOTS[A._current_slot_index_locked()]===slot
                    @test A.SESSION_VIEWS[held]===view && A.view.get_model_data(held,0)==data
                end
                @test A.open(joinpath(FIXTURES,name))===nothing
                @test length(A.model.list())==length(before.models)+1
                @test A.SESSION_VIEWS[held]===view && A.view.get_model_data(held,0)==data
                A.model.set_current(before.name)
                @test A.MODEL_SLOTS[A._current_slot_index_locked()]===slot
                @test A.mesh.get_nodes()==before.nodes && A.mesh.get_elements()==before.cells
            finally
                A.finalize()
            end
        end
    end
end

@testset "CAD suffix aliases and filenames are data" begin
    mktempdir() do directory
        for (name,suffix) in (("brep_occ_box.step",".STP"),("brep_occ_box.igs",".IGES"))
            path=joinpath(directory,"native shape 台灣"*suffix)
            cp(joinpath(FIXTURES,name),path)
            A.initialize()
            try
                A.option("General.Terminal",0)
                @test A.open(path)===nothing
                @test A.model.get_current()==splitext(basename(path))[1]
                @test !isempty(A.model.get_entities(2))
            finally
                A.finalize()
            end
        end
    end
end
end
