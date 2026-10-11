module NativeCurveContextEventsTests
using Tessella,Test
const A=Tessella.API
include("../../validation/native_curve_events/native_curve_context_events_fixture_v1.jl")
float_bits(value)=string(reinterpret(UInt64,Float64(value));base=16,pad=16)
function open_text(text)
    mktemp() do path,stream
        write(stream,text);close(stream)
        A.open_geo!(path)
    end
end
function check_event(name;nodes=false)
    expected=NATIVE_CURVE_CONTEXT_EVENTS["api_events"][name]
    @test float_bits(A.option("General.BoundingBoxSize"))==expected["lc_bits"]
    if nodes
        for row in expected["curve_nodes"]
            tags,coordinates,parameters=A.mesh.get_nodes(1,row["owner"])
            @test float_bits.(parameters)==row["parameter_bits"]
            @test float_bits.(coordinates)==row["coordinate_bits"]
        end
    end
end
@testset "Native global Curve context synchronization events" begin
    A.initialize()
    try
        check_event("initialized")
        A.option("Mesh.Renumber",0)
        A.model.add("small")
        open_text("Point(1)={0,0,0};Point(2)={.3,.8,.2};Point(3)={.7,-.2,.5};Point(4)={1,.3,0};Spline(11)={1,2,3,4};Transfinite Curve{11}=4;")
        check_event("small:geo-synchronize")
        A.mesh.generate(1);check_event("small:generate-first";nodes=true)
        A.mesh.generate(1);check_event("small:generate-repeat";nodes=true)
        A.model.add("large");check_event("large:add")
        A.model.add_point(-100.,-200.,-300.;tag=1)
        A.model.add_point(100.,200.,300.;tag=2)
        A.model.add_line(1,2;tag=11)
        A.model.geo.synchronize();check_event("large:geo-synchronize")
        A.model.set_current("small");check_event("small:set-current-after-large-sync")
        A.mesh.generate(1);check_event("small:generate-after-large-sync";nodes=true)
        A.mesh.clear();check_event("small:clear-mesh-after-large-sync")
        A.mesh.generate(1);check_event("small:generate-after-clear";nodes=true)
        A.model.set_tag(1,11,91);check_event("small:retag-curve")
        A.mesh.generate(1);check_event("small:generate-retagged";nodes=true)
        A.model.set_coordinates(2,.3,.9,.2);check_event("small:set-point-coordinates")
        A.model.synchronize()
        check_event("small:forced-geo-synchronize")
        A.model.set_current("large");check_event("large:set-current")
        A.model.remove();check_event("large:remove-current")
    finally
        A.finalize()
    end
    @testset "GEO Point, Mesh, SyncModel and BoundingBox events" begin
        for row in NATIVE_CURVE_CONTEXT_EVENTS["geo_rows"]
            A.initialize()
            try
                result=open_text(row["geo"])
                for (name,expected) in row["variable_bits"]
                    @test float_bits(result.values[name])==only(expected)
                end
                @test float_bits(result.native_context_lc)==row["final_lc_bits"]
                @test float_bits(A.option("General.BoundingBoxSize"))==row["final_lc_bits"]
            finally
                A.finalize()
            end
        end
    end
end
end
