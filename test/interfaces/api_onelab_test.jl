# `onelab` module UI01: session parameter records, JSON exchange, glob name
# queries, changed flags, clear, honest `run` rejection, and the `.geo`
# DefineNumber/SetNumber/GetNumber exchange through the session file frontends.
module OnelabChecks
using Test, Tessella
const A=Tessella.API

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "typed setters/getters and name queries" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # Scalar and vector number writes; vectors round-trip.
        A.onelab.set_number("Group/a",1.5)
        A.onelab.set_number("Group/b",[2.0,3.0])
        @test A.onelab.get_number("Group/a")==[1.5]
        @test A.onelab.get_number("Group/b")==[2.0,3.0]
        A.onelab.set_string("Group/s","hello")
        @test A.onelab.get_string("Group/s")==["hello"]
        # Name uniqueness across types: set_string replaces the number record.
        A.onelab.set_string("Group/a","now a string")
        @test A.onelab.get_string("Group/a")==["now a string"]
        @test occursin("no ONELAB number parameter",
            _errtext(()->A.onelab.get_number("Group/a")))
        A.onelab.set_number("Group/a",4.0)
        # Glob name queries.
        @test A.onelab.get_names()==["Group/a","Group/b","Group/s"]
        @test A.onelab.get_names("Group/*")==["Group/a","Group/b","Group/s"]
        @test A.onelab.get_names("Group/?")==["Group/a","Group/b","Group/s"]
        @test A.onelab.get_names("*/b")==["Group/b"]
        @test A.onelab.get_names("nothing*")==String[]
        # Missing/type-mismatched getters fail loudly.
        @test occursin("no ONELAB number parameter",
            _errtext(()->A.onelab.get_number("missing")))
        @test occursin("no ONELAB string parameter",
            _errtext(()->A.onelab.get_string("Group/b")))
        # Input validation mirrors the API's numeric/string conventions.
        @test_throws ArgumentError A.onelab.set_number("x",true)
        @test_throws ArgumentError A.onelab.set_number("x",Inf)
        @test_throws ArgumentError A.onelab.set_number("x","s")
        @test_throws ArgumentError A.onelab.set_number("",1.0)
        @test_throws ArgumentError A.onelab.set_number(3,1.0)
        @test_throws ArgumentError A.onelab.set_string("x",1.0)
        @test_throws ArgumentError A.onelab.get_number(3)
    finally
        A.finalize()
    end
end

@testset "JSON exchange" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # A full-record JSON set: values plus optional members.
        A.onelab.set("""[
            {"type":"number","name":"Geo/lc","values":[0.1,0.2],
             "min":0.0,"max":1.0,"step":0.05,"choices":[0.1,0.2],
             "valueLabels":{"0.1":"coarse","0.2":"fine"},
             "attributes":{"Highlight":"Pink"},
             "label":"Geo","shortHelp":"mesh size","readOnly":true},
            {"type":"string","name":"Geo/name","values":["square"],
             "choices":["square","disk"],"visible":false,"closed":true,
             "kind":"attribute"}
        ]""")
        json=A.onelab.get()
        @test occursin("\"name\": \"Geo/lc\"",json)
        @test occursin("\"values\": [0.1, 0.2]",json)
        @test occursin("\"min\": 0",json)
        @test occursin("\"max\": 1",json)
        @test occursin("\"step\": 0.05",json)
        @test occursin("\"choices\": [0.1, 0.2]",json)
        @test occursin("\"valueLabels\": {\"0.1\": \"coarse\", \"0.2\": \"fine\"}",json)
        @test occursin("\"attributes\": {\"Highlight\": \"Pink\"}",json)
        @test occursin("\"label\": \"Geo\"",json)
        @test occursin("\"shortHelp\": \"mesh size\"",json)
        @test occursin("\"readOnly\": true",json)
        @test occursin("\"visible\": false",json)
        @test occursin("\"closed\": true",json)
        @test occursin("\"kind\": \"attribute\"",json)
        @test occursin("\"values\": [\"square\"]",json)
        @test occursin("\"choices\": [\"square\", \"disk\"]",json)
        # Round-trip: the emitted JSON re-sets into an identical store.
        snapshot=json
        A.onelab.clear()
        @test A.onelab.get()=="[]"
        A.onelab.set(snapshot)
        @test A.onelab.get()==snapshot
        @test A.onelab.get_number("Geo/lc")==[0.1,0.2]
        @test A.onelab.get_string("Geo/name")==["square"]
        # Name-filtered get returns a one-element array; unknown names give [].
        one=A.onelab.get("Geo/lc")
        @test occursin("\"Geo/lc\"",one) && !occursin("Geo/name",one)
        @test A.onelab.get("absent")=="[]"
        # Typed set replaces the record's optional members (upstream assign).
        A.onelab.set_number("Geo/lc",9.0)
        @test !occursin("\"min\"",A.onelab.get("Geo/lc"))
        # set validates the payload and the format.
        @test_throws ArgumentError A.onelab.set("not json")
        @test_throws ArgumentError A.onelab.set("""[{"type":"number"}]""")
        @test_throws ArgumentError A.onelab.set("""[{"type":"sphere","name":"x"}]""")
        @test_throws ArgumentError A.onelab.set("[]","xml")
        @test_throws ArgumentError A.onelab.get("","xml")
        @test_throws ArgumentError A.onelab.set(3)
        # Unicode escapes and UTF-8 names decode on input.
        A.onelab.set("""[{"type":"string","name":"u\\u00e9","values":["a\\nb"]}]""")
        @test A.onelab.get_string("u\u00e9")==["a\nb"]
        @test occursin("u\u00e9",A.onelab.get("u\u00e9"))
    finally
        A.finalize()
    end
end

@testset "changed flags, clear, run" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        @test !A.onelab.get_changed("a")
        A.onelab.set_number("a",1.0)
        @test A.onelab.get_changed("a")
        @test !A.onelab.get_changed("b")
        A.onelab.clear_changed()
        @test !A.onelab.get_changed("a")
        A.onelab.set_changed("b")
        @test A.onelab.get_changed("b")
        # clear(name) removes matching records and their flags.
        A.onelab.set_number("keep",1.0);A.onelab.set_number("drop1",2.0)
        A.onelab.set_number("drop2",3.0)
        A.onelab.clear("drop*")
        @test A.onelab.get_names()==["a","keep"]
        @test !A.onelab.get_changed("drop1")
        @test A.onelab.get_changed("b")
        A.onelab.clear()
        @test isempty(A.onelab.get_names())
        @test !A.onelab.get_changed("b")
        # run(): zero registered clients is the upstream no-op; a real client
        # name or command line rejects loudly (solver clients unimplemented).
        @test A.onelab.run()===nothing
        @test occursin("solver clients are not implemented",
            _errtext(()->A.onelab.run("GetDP")))
        @test occursin("solver clients are not implemented",
            _errtext(()->A.onelab.run("","getdp solver.pro")))
        @test_throws ArgumentError A.onelab.run(3)
        @test_throws ArgumentError A.onelab.get_names(3)
        # The store needs an initialized session.
        A.finalize()
        @test occursin("no current model",
            _errtext(()->A.onelab.get_names()))
        A.initialize()
        # clear() preserves the ONELAB store (model/view-only reset upstream).
        A.option("General.Terminal",0)
        A.onelab.set_number("persist",8.0)
        A.clear()
        @test A.onelab.get_number("persist")==[8.0]
        # initialize() resets the store.
        A.initialize()
        A.option("General.Terminal",0)
        @test isempty(A.onelab.get_names())
    finally
        A.finalize()
    end
end

@testset ".geo parameter exchange" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # The session value wins the script's DefineNumber default and the
        # script's own writes/adoptions land in the session store.
        A.onelab.set_number("pv",3.0)
        A.onelab.set_string("gw","preset")
        path=tempname()*".geo"
        write(path,"""
pv = DefineNumber(5, Name "pv");
rx = GetNumber("pv");
SetNumber("echo_pv", rx);
SetNumber("qv", 7.5);
SetString("gw", "kept");
SetString("gs", "fresh");
ds = DefineString("df = 9");
UndefineConstant[pv];
Point(1) = {rx, 0, 0, 1};
""")
        A.open_geo!(path)
        # pv's server value (3) drove the script; the script mirrored it back
        # into `echo_pv`, and `UndefineConstant` removed `pv` itself.
        @test A.onelab.get_number("echo_pv")==[3.0]
        @test occursin("no ONELAB number parameter",
            _errtext(()->A.onelab.get_number("pv")))
        @test A.onelab.get_number("qv")==[7.5]
        # A pre-existing string parameter is not overwritten by SetString...
        @test A.onelab.get_string("gw")==["kept"]
        @test A.onelab.get_string("gs")==["fresh"]
        # DefineString stores the full spec string under the first-word name.
        @test A.onelab.get_string("df")==["df = 9"]
        # Script-set parameters arrive changed; the Point used rx=3.
        @test A.onelab.get_changed("qv")
        xyz=A.model.get_value(0,1,Float64[])
        @test xyz[1:3]≈[3.0,0.0,0.0]
        # A failing script leaves the ONELAB store untouched (atomic).
        bad=tempname()*".geo"
        write(bad,"SetNumber(\"lost\", 1);\nExit 3;\n")
        before=A.onelab.get()
        @test_throws ArgumentError A.open_geo!(bad)
        @test A.onelab.get()==before
        @test !("lost" in A.onelab.get_names())
    finally
        A.finalize()
    end
end

end # module OnelabChecks
