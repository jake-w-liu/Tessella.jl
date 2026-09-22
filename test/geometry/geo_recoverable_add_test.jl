using Test
using Tessella

const _GE=Tessella.GeoExec
const _GI=Tessella.IO

# `execute_geo` throws the accumulated `yymsg(0)` diagnostics at end of parse;
# to inspect the model state those recoverable statements left behind, the
# internal driver runs the same sequence without the final throw.
function _exec_model(source::AbstractString)
    return mktemp() do path,io
        write(io,source)
        close(io)
        params=_GI.read_geo_params(path)
        model=Tessella.Model.GeoModel()
        context=_GI._GeoNumericContext()
        context.file_name=String(path)
        context.fields=empty!(params.fields)
        context.soft_unknown_reads=true
        context.stderr_diagnostics=false
        context.entity_name_lookup=(dim,tag,kind)->
            kind===:physical ? get(model.physical_names,(dim,tag),"") :
                get(model.entity_names,(dim,tag),"")
        allocator_state=_GE._GeoTagAllocatorState()
        context.exec_hook=src->_GE._geo_exec_value_term(
            model,src,context,allocator_state)
        statements=_GE._geo_exec_statements(path)
        _GE._exec_geo_statements!(
            model,statements,firstindex(statements),lastindex(statements),
            context,allocator_state,Ref(0))
        _GE._geo_sync_physical_view!(model,context)
        return model,context
    end
end

function _execute_source(source::AbstractString;mesh_dim=0)
    return mktemp() do path,io
        write(io,source)
        close(io)
        execute_geo(path;mesh_dim=mesh_dim)
    end
end

function _exec_error(source::AbstractString)
    try
        _execute_source(source)
        return nothing
    catch err
        return err
    end
end

@testset "recoverable duplicate entity definitions" begin
    # Gmsh 4.15.2: `GEO <type> with tag N already exists` (Msg::Error) +
    # `Could not add <name>` (yymsg 0); the stream continues and the original
    # entity is kept.
    for (stmt,name) in (("Point(1)={0,0,0};","point"),
                        ("Line(1)={1,2};","line"),
                        ("Circle(1)={1,2,3};","circle"),
                        ("Spline(1)={1,2,3};","spline"),
                        ("BSpline(1)={1,2,3};","BSpline"),
                        ("Bezier(1)={1,2,3};","Bezier"),
                        ("Curve Loop(1)={1};","curve loop"),
                        ("Plane Surface(1)={1};","plane surface"),
                        ("Surface Loop(1)={1};","surface loop"),
                        ("Volume(1)={1};","volume"))
        err=_exec_error(stmt*stmt)
        @test err!==nothing
        @test occursin("Could not add $name",err.msg)
    end
    # The first definition wins.
    m,_=_exec_model("Point(1)={0,0,0}; Point(1)={9,9,9};")
    @test m.points[1]==(0.0,0.0,0.0)
    m,_=_exec_model("Point(1)={0,0,0}; Point(2)={1,0,0};" *
                    "Line(1)={1,2}; Line(1)={2,1};")
    @test m.curves[1]==(1,2)
    # OCC-gated solids use the `OpenCASCADE` wording.
    err=_exec_error(
        "SetFactory(\"OpenCASCADE\");" *
        "Sphere(9)={0,0,0,1}; Sphere(9)={0,0,0,2};")
    @test err!==nothing
    @test occursin("Could not add sphere",err.msg)
    # `Curve(n)` is a line alias: same duplicate semantics.
    err=_exec_error("Point(1)={0,0,0}; Point(2)={1,0,0};" *
                    "Curve(1)={1,2}; Curve(1)={2,1};")
    @test occursin("Could not add line",err.msg)
end

@testset "partial entity storage after failed creation" begin
    # Unknown control points are dropped; the curve is still stored with the
    # found subset and the tag is consumed.
    m,c=_exec_model("Point(9)={0,0,0}; Line(1)={9,999};")
    @test haskey(m.curves,1)
    @test m.curves[1]==(9,9)
    @test any(e->occursin("Could not add line",e),c.exec_errors)
    # A second attempt hits the consumed tag.
    err=_exec_error("Point(9)={0,0,0}; Line(1)={9,999}; Line(1)={9,9};")
    @test count("Could not add line",err.msg)==2
    # Unknown curves in a curve loop: stored blindly, tag consumed.
    m,c=_exec_model("Curve Loop(1)={999};")
    @test m.loops[1]==[999]
    @test isempty(c.exec_errors)
    # Empty lists store an empty loop/surface loop/volume; empty surfaces and
    # lines fail without storing.
    m,_=_exec_model("Curve Loop(1)={}; Surface Loop(2)={}; Volume(3)={};")
    @test m.loops[1]==Int[] && m.surface_loops[2]==Int[] && m.volumes[3]==Int[]
    m,_=_exec_model("Plane Surface(1)={}; Surface(2)={}; Line(3)={};")
    @test isempty(m.surfaces) && isempty(m.curves)
    # Failed volume construction still consumes the tag — the second attempt
    # reports the duplicate and emits `Could not add volume` again.
    err=_exec_error("Volume(1)={999}; Volume(1)={888};")
    @test err!==nothing
    @test count("Could not add volume",err.msg)==2
    m,_=_exec_model("Volume(1)={999};")
    @test m.volumes[1]==[999]
    m,_=_exec_model("Plane Surface(1)={999};")
    @test m.surfaces[1]==[999]
    # Unknown-member ruled surface is NOT stored.
    m,_=_exec_model("Surface(2)={999};")
    @test isempty(m.surfaces)
    # Surface loop members are stored without resolution.
    m,_=_exec_model("Surface Loop(1)={999,888};")
    @test m.surface_loops[1]==[999,888]
end

@testset "curve loop partial sorting" begin
    # A dead-end while chaining stores the sorted prefix and reports
    # `Curve loop N is wrong` through Msg::Error; the tag is consumed.
    m,c=_exec_model("Point(1)={0,0,0};Point(2)={1,0,0};Point(3)={0,1,0};" *
                    "Line(1)={1,2};Line(8)={3,1};Curve Loop(1)={1,8};")
    @test m.loops[1]==[1]
    @test any(e->occursin("Could not add curve loop",e),c.exec_errors)
    err=_exec_error("Point(1)={0,0,0};Point(2)={1,0,0};Point(3)={0,1,0};" *
                    "Line(1)={1,2};Line(8)={3,1};" *
                    "Curve Loop(1)={1,8};Curve Loop(1)={1,8};")
    @test count("Could not add curve loop",err.msg)==2
    # A closed chain sorts cleanly.
    m,c=_exec_model("Point(1)={0,0,0};Point(2)={1,0,0};Point(3)={0,1,0};" *
                    "Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};" *
                    "Curve Loop(1)={1,2,3};")
    @test m.loops[1]==[1,2,3]
    @test isempty(c.exec_errors)
end

@testset "circle/ellipse EndCurve diagnostics" begin
    # Non-cocircular control points emit forward AND reversed diagnostics.
    _,c=_exec_model("Point(1)={2,0,0};Point(2)={0,0,0};Point(3)={0,1,0};" *
                    "Circle(1)={1,2,3};")
    @test c.msg_error_count==2
    @test any(e->occursin("Could not add circle",e),c.exec_errors)
    # `Plane{n}` normal runs the check twice per orientation.
    _,c=_exec_model("Point(1)={2,0,0};Point(2)={0,0,0};Point(3)={0,1,0};" *
                    "Circle(1)={1,2,3} Plane{0,0,1};")
    @test c.msg_error_count==4
    # A valid circle is silent.
    _,c=_exec_model("Point(1)={1,0,0};Point(2)={0,0,0};Point(3)={0,1,0};" *
                    "Circle(1)={1,2,3};")
    @test c.msg_error_count==0 && isempty(c.exec_errors)
    # Malformed short ellipse: upstream's reversed pass reads past the list.
    _,c=_exec_model("Point(1)={1,0,0};Point(2)={0,0,0};Point(3)={0,1,0};" *
                    "Ellipse(2)={1,2,999};")
    @test any(e->occursin("Could not add ellipse",e),c.exec_errors)
    @test c.msg_error_count>=2
end

@testset "degenerate line warning" begin
    # Start==end within geometrical tolerance: stored with a warning.
    m,c=_exec_model("Point(9)={1,1,0}; Line(1)={9,9};")
    @test m.curves[1]==(9,9)
    @test any(w->occursin("closer than the geometrical tolerance",w),
              c.exec_warnings)
    @test isempty(c.exec_errors)
    # Multi-point line is a polyline.
    m,_=_exec_model("Point(9)={0,0,0};Point(10)={1,0,0};Point(11)={2,0,0};" *
                    "Line(1)={9,10,11};")
    @test m.curves[1]==(9,11)
end

@testset "parametric sphere registry" begin
    # `Sphere(n)={center,point}` registers a parametric surface — no volume.
    m,c=_exec_model("Point(1)={0,0,0};Point(2)={1,0,0};Sphere(9)={1,2};")
    @test isempty(m.volumes)
    @test haskey(c.parametric_surfaces,9)
    @test c.parametric_surface!==nothing
    # A duplicate tag reports through Msg::Error and replaces the entry.
    m,c=_exec_model("Point(1)={0,0,0};Point(2)={1,0,0};" *
                    "Sphere(9)={1,2};Sphere(9)={1,2};")
    @test c.msg_error_count==1
    @test isempty(c.exec_errors)
    # `Coordinates Surface n` selects; a missing tag reports and clears.
    _,c=_exec_model("Coordinates Surface 77;")
    @test c.msg_error_count==1
    @test c.parametric_surface===nothing
    _,c=_exec_model("Point(1)={0,0,0};Point(2)={1,0,0};" *
                    "Sphere(9)={1,2};Euclidian Coordinates;")
    @test c.parametric_surface===nothing
    @test haskey(c.parametric_surfaces,9)
    # `Delete Model` clears the registry: no `already exists` on re-register.
    _,c=_exec_model("Point(1)={0,0,0};Point(2)={1,0,0};Sphere(9)={1,2};" *
                    "Delete Model;Sphere(9)={1,2};")
    @test !any(e->occursin("already exists",e),c.exec_errors)
    # `PolarSphere` registers identically.
    _,c=_exec_model("Point(1)={0,0,0};Point(2)={1,0,0};" *
                    "PolarSphere(9)={1,2};Coordinates Surface 9;")
    @test isempty(c.exec_errors) && c.msg_error_count==0
end

@testset "named and alternate grammar forms" begin
    # `Curve "name"(n)` / `Surface "name"(n)` are syntax errors upstream.
    for bad in ("Curve \"x\"(1)={1};","Surface \"y\"(2)={1};")
        err=_exec_error(bad)
        @test err!==nothing
        @test occursin("syntax error (\")",err.msg)
    end
    # `Surface Loop(n)={...} Using Sewing` is accepted and ignored.
    m,c=_exec_model("Point(1)={0,0,0};Point(2)={1,0,0};Point(3)={0,1,0};" *
                    "Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};" *
                    "Curve Loop(1)={1,2,3};Plane Surface(1)={1};" *
                    "Surface Loop(1)={1} Using Sewing;")
    @test m.surface_loops[1]==[1]
    @test isempty(c.exec_errors)
end

@testset "exec_errors vs msg_error_count separation" begin
    # Inner Msg::Error diagnostics count but do not throw.
    parsed=_execute_source("Point(1)={0,0,0};Point(2)={1,0,0};" *
                           "Sphere(9)={1,2};Sphere(9)={1,2};")
    @test parsed.msg_error_count==1
    # `Could not add` yymsg diagnostics accumulate and throw.
    err=_exec_error("Point(1)={0,0,0};Point(1)={0,0,0};")
    @test err!==nothing
    @test occursin("Could not add point",err.msg)
end
