module OnelabTypedChecks
using Test,Tessella
const A=Tessella.API

@testset "independent ONELAB number and string namespaces" begin
    A.initialize();A.option("General.Terminal",0)
    try
        A.onelab.set_number("same",[1.,2.])
        A.onelab.set_string("same",["a","b"])
        @test A.onelab.get_number("same")==[1.,2.]
        @test A.onelab.get_string("same")==["a","b"]
        @test A.onelab.get_names()==["same","same"]
        records=A._onelab_json_parse(A.onelab.get(),"test")
        @test [record["type"] for record in records]==["number","string"]
        @test only(A._onelab_json_parse(A.onelab.get("same"),"test"))["type"]=="number"

        # Register consumers through the production wire projection. Changes
        # to one typed record leave the other record's consumers unchanged.
        A._onelab_store_wire(A._ONELAB_SESSION_STORE,A.SESSION_ONELAB[(:number,"same")],"number reader")
        A._onelab_store_wire(A._ONELAB_SESSION_STORE,A.SESSION_ONELAB[(:string,"same")],"string reader")
        A.onelab.clear_changed()
        A.onelab.set_string("same",["changed"])
        @test !A.onelab.get_changed("number reader")
        @test A.onelab.get_changed("string reader")
        @test A.onelab.get_number("same")==[1.,2.]

        mktempdir() do directory
            path=joinpath(directory,"typed.geo")
            write(path,"n=GetNumber(\"same\"); s=GetString(\"same\");\n"*
                "SetNumber(\"same\",3); SetString(\"same\",\"geo\");\n")
            result=Base.invokelatest(A.open_geo!,path)
            @test result.values["n"]==1. && result.strings["s"]==["changed"]
            @test A.onelab.get_number("same")==[3.]
            @test A.onelab.get_string("same")==["geo"]
        end
        # Public Gmsh control: an untyped named clear selects number first;
        # the next clear selects the independent string with that name.
        A.onelab.clear("same")
        @test A.onelab.get_string("same")==["geo"]
        @test only(A._onelab_json_parse(A.onelab.get("same"),"test"))["type"]=="string"
        A.onelab.clear("same")
        @test isempty(A.onelab.get_names())
        A.onelab.set_number("same",1.);A.onelab.set_string("same","a")
        A.onelab.clear()
        @test isempty(A.onelab.get_names()) && isempty(A._ONELAB_SESSION_STORE.typed_changed)
    finally
        A.finalize()
    end
end

@testset "template expansion checks its remaining budget before copying" begin
    A.initialize()
    try
        A.onelab.set_string("expansion",[repeat("x",1024*1024)])
        state=A._OnelabTemplateState(A._ONELAB_SESSION_STORE,"budget",false,0,0,Set{String}(),message->nothing)
        A._onelab_template_substitute(state,"OL.get(expansion)")
        state.output_bytes=64*1024*1024-32
        @test_throws ArgumentError A._onelab_template_substitute(state,"prefix OL.get(expansion)")
        allocation=@allocated try
            A._onelab_template_substitute(state,"prefix OL.get(expansion)")
        catch exception
            exception isa ArgumentError || rethrow()
        end
        @test allocation<65536
        @test A._onelab_template_substitute(state,"λ OL.get(expansion)";limit=2*1024*1024)==
            "λ "*repeat("x",1024*1024)
    finally
        A.finalize()
    end
end
end
