module OnelabMetadataChecks
using Test,Tessella
const A=Tessella.API

function parameter(kind,name)
    return only(record for record in A._onelab_json_parse(A.onelab.get(),"metadata test")
        if record["type"]==String(kind) && record["name"]==name)
end

@testset "GEO declarations retain full ONELAB metadata" begin
    A.initialize();A.option("General.Terminal",0)
    try
        mktempdir() do directory
            path=joinpath(directory,"metadata.geo")
            write(path,"counter=0; n=DefineNumber(2,Name \"Controls/N\",Min 1,Max 10,"*
                "Step 0.5,Choices{counter++=\"zero\",counter++=\"one\"},"*
                "Label \"Count\",Help \"Number help\",ReadOnly 1,Visible 0,"*
                "Closed 1,Highlight \"Yellow\",Loop \"2\",Graph \"3\",AutoCheck 0);\n"*
                "s=DefineString(\"a\",Name \"Controls/S\",Label \"Mode\",Help \"String help\","*
                "Choices{\"a\",\"b\"},Kind \"file\",MultipleSelection \"10\",Macro \"Action\");\n")
            @test A.parser.parse(path)===nothing
            n=parameter(:number,"Controls/N");s=parameter(:string,"Controls/S")
            @test n["values"]==[2.]
            @test (n["min"],n["max"],n["step"])==(1.,10.,.5)
            @test n["choices"]==[0.,1.]
            @test n["valueLabels"]==Dict("0"=>"zero","1"=>"one")
            @test n["label"]=="Count" && n["shortHelp"]=="Number help"
            @test n["readOnly"] && !n["visible"] && n["closed"]
            @test n["attributes"]==Dict("Choices"=>"zero","Closed"=>"1",
                "Highlight"=>"Yellow","Loop"=>"2","Graph"=>"3","AutoCheck"=>"0")
            @test A.parser.get_number("counter")==[2.]
            @test s["label"]=="Mode" && s["shortHelp"]=="String help"
            @test s["choices"]==["a","b"] && s["kind"]=="file"
            @test s["attributes"]==Dict("Choices"=>"a","Kind"=>"file",
                "MultipleSelection"=>"10","Macro"=>"Action")
            wire=A._onelab_store_wire(A._ONELAB_SESSION_STORE,
                A.SESSION_ONELAB[(:number,"Controls/N")],"metadata reader")
            decoded=A._onelab_wire_decode(A._onelab_wire_encode(wire),"metadata test")
            @test (decoded.min,decoded.max,decoded.step)==(1.,10.,.5)
            @test decoded.value_labels==Dict(0.0=>"zero",1.0=>"one")
            @test decoded.read_only && !decoded.visible && decoded.help=="Number help"

            # Full record metadata and client ownership publish together only
            # after the containing file succeeds.
            snapshot=A.onelab.get();changed=A.onelab.get_changed("metadata reader")
            write(path,"n=DefineNumber(99,Name \"Controls/N\",Label \"bad\",Min -99);"*
                "s=DefineString(\"bad\",Name \"Controls/S\",ReadOnly 1); Exit 7;\n")
            @test_throws ArgumentError A.parser.parse(path)
            @test A.onelab.get()==snapshot
            @test A.onelab.get_changed("metadata reader")==changed
            @test A.parser.get_number("n")==[2.]
        end
    finally
        A.finalize()
    end
end

@testset "GEO metadata updates retain server values and independent typed records" begin
    A.initialize();A.option("General.Terminal",0)
    try
        A.onelab.set("""[
          {"type":"number","name":"same","values":[5,6],"min":4,"max":20,
           "step":2,"choices":[4,8],"label":"Old","shortHelp":"Old help",
           "visible":false,"readOnly":true,"attributes":{"Custom":"Old"}},
          {"type":"string","name":"same","values":["old","tail"],"label":"String old"}
        ]""")
        mktempdir() do directory
            path=joinpath(directory,"update.geo")
            write(path,"n=DefineNumber(2,Name \"same\"); s=DefineString(\"new\",Name \"same\");\n")
            @test A.parser.parse(path)===nothing
            n=parameter(:number,"same");s=parameter(:string,"same")
            @test n["values"]==[2.] && n["readOnly"]
            @test (n["min"],n["max"],n["step"],n["choices"])==(4.,20.,2.,[4.,8.])
            @test n["label"]=="Old" && n["shortHelp"]=="Old help" && !n["visible"]
            @test n["attributes"]==Dict("Custom"=>"Old")
            @test s["values"]==["old","tail"] && s["label"]=="String old"
            @test A.parser.get_number("n")==[2.] && A.parser.get_string("s")==["old"]

            write(path,"n=DefineNumber(3,Name \"same\",ReadOnly 0,Min 1,Max 10,Step 0.5,"*
                "Choices{1,2},Range{100,200,10},Label \"New\");"*
                "s=DefineString(\"new\",Name \"same\",ReadOnly 1,Label \"String new\");\n")
            @test A.parser.parse(path)===nothing
            n=parameter(:number,"same");s=parameter(:string,"same")
            @test n["values"]==[2.] && !get(n,"readOnly",false)
            @test (n["min"],n["max"],n["step"],n["choices"])==(1.,10.,.5,[1.,2.])
            @test !haskey(n["attributes"],"Range") && n["attributes"]["Custom"]=="Old"
            @test n["label"]=="New" && n["shortHelp"]=="Old help"
            @test s["values"]==["new"] && s["label"]=="String new" && s["readOnly"]
            @test A.onelab.get_names()==["same","same"]
        end
    finally
        A.finalize()
    end
end

@testset "standalone GEO declarations reach a native client with metadata" begin
    mktempdir() do directory
        path=joinpath(directory,"standalone.geo");output=joinpath(directory,"metadata.txt")
        quote_arg(value)=Sys.iswindows() ? A._onelab_windows_quote(String(value)) :
            "\""*replace(String(value),'\\'=>"\\\\",'\"'=>"\\\"")*"\""
        fixture=joinpath(@__DIR__,"..","fixtures","onelab_native_client.jl")
        command=join(quote_arg.([String(Base.julia_cmd().exec[1]),"--startup-file=no",
            "--compile=min","-O0",fixture,"metadata",output])," ")
        write(path,"n=DefineNumber(2,Name \"Controls/N\",Min 1,Max 10,Step 0.5,"*
            "Label \"Count\",Help \"Number help\");"*
            "s=DefineString(\"a\",Name \"Controls/S\",Choices{\"a\",\"b\"},Kind \"file\");"*
            "OnelabRun(\"metadata\",'"*command*"');\n")
        result=Tessella.GeoExec.execute_geo(path)
        @test result.msg_error_count==0
        records=readlines(output)
        @test length(records)==2
        n=A._onelab_wire_decode(records[1],"metadata test")
        s=A._onelab_wire_decode(records[2],"metadata test")
        @test n.values_f==[2.] && (n.min,n.max,n.step)==(1.,10.,.5)
        @test n.label=="Count" && n.help=="Number help"
        @test s.values_s==["a"] && s.choices_s==["a","b"] && s.kind_name=="file"
        @test isempty(A._ONELAB_ACTIVE_ENGINES)
    end
end
end
