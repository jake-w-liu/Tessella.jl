# Differential oracle for Point/Line/Triangle/Tetrahedron reference quadrature.
# This uses the locally installed Gmsh 4.15.2 Julia API and never starts the GUI.
using Tessella
using SHA
using Tessella.Elements: MSH_CATALOG

function _gmsh_binding()
    configured=get(ENV,"GMSH_JULIA_API","")
    candidates=String[]
    isempty(configured) || push!(candidates,configured)
    executable=Sys.which("gmsh")
    if executable!==nothing
        prefix=dirname(dirname(realpath(executable)))
        push!(candidates,joinpath(prefix,"lib","gmsh.jl"))
    end
    append!(candidates,["/opt/homebrew/lib/gmsh.jl","/usr/local/lib/gmsh.jl"])
    for candidate in unique(candidates)
        isfile(candidate) && return candidate
    end
    error("mesh-quadrature differential: Gmsh Julia API not found; " *
          "set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION=="4.15.2" || error(
    "mesh-quadrature differential requires Gmsh API 4.15.2, found " *
    gmsh.GMSH_API_VERSION)

function _compare_values(label,gmsh_values,tessella_values)
    length(gmsh_values)==length(tessella_values) || error(
        "$label lengths differ: Gmsh=$(length(gmsh_values)), " *
        "Tessella=$(length(tessella_values))")
    for index in eachindex(gmsh_values,tessella_values)
        isapprox(tessella_values[index],gmsh_values[index];
                 atol=2e-14,rtol=2e-14) || error(
            "$label differs at $index: Gmsh=$(gmsh_values[index]), " *
            "Tessella=$(tessella_values[index])")
    end
    return nothing
end

function _write_case!(stream,element_type,rule,coordinates,weights)
    write(stream,htol(Int64(element_type)))
    write(stream,htol(UInt64(ncodeunits(rule))))
    write(stream,codeunits(rule))
    for values in (coordinates,weights)
        write(stream,htol(UInt64(length(values))))
        foreach(value->write(
            stream,htol(reinterpret(UInt64,Float64(value)))),values)
    end
    return nothing
end

function _rejects_argument(f::Function)
    try
        f()
        return false
    catch err
        err isa ArgumentError || rethrow()
        return true
    end
end

# Mutant analysis:
# - Replacing Duffy rules with Cartesian tensor points is rejected by every
#   triangle/tetrahedron coordinate and weight comparison.
# - Using the line point-count law for a simplex is rejected by exact array
#   lengths at CompositeGauss0, CompositeGauss5, and CompositeGauss29.
# - Transposing the nested integration loops is rejected by sequential flattened
#   coordinate comparisons at asymmetric points.
# - Dispatching by interpolation order instead of parent family is rejected by
#   every higher-order and serendipity catalog entry.
# - Silently treating malformed suffixes as order zero is rejected locally; Gmsh
#   is deliberately not called with those unsafe inputs.

try
    gmsh.initialize(String[],false)
    gmsh.option.setNumber("General.Terminal",0)
    stream=IOBuffer()
    comparison_count=0

    function compare_case(element_type::Int,rule::String)
        gmsh_coordinates,gmsh_weights=
            gmsh.model.mesh.getIntegrationPoints(element_type,rule)
        tessella_coordinates,tessella_weights=
            Tessella.API.mesh.get_integration_points(element_type,rule)
        label="type-$element_type $rule"
        _compare_values(
            "$label coordinates",gmsh_coordinates,tessella_coordinates)
        _compare_values("$label weights",gmsh_weights,tessella_weights)
        _write_case!(
            stream,element_type,rule,tessella_coordinates,tessella_weights)
        comparison_count+=1
        return nothing
    end

    supported=sort!([
        (Int(element_type),spec.family)
        for (element_type,spec) in MSH_CATALOG
        if spec.family in (:pnt,:lin,:tri,:tet)];by=first)
    for (element_type,_) in supported
        compare_case(element_type,"Gauss2")
        compare_case(element_type,"CompositeGauss5")
    end

    for rule in ("Gauss","Gauss0","Gauss1","Gauss2","Gauss3",
                 "Gauss4","Gauss5","Gauss0005",
                 "CompositeGauss","CompositeGauss0","CompositeGauss1",
                 "CompositeGauss2","CompositeGauss5",
                 "CompositeGauss12","CompositeGauss29")
        for element_type in (15,1,2,4)
            compare_case(element_type,rule)
        end
    end

    length(Tessella.API.mesh.get_integration_points(
        1,"CompositeGauss255")[2])==128 || error(
        "bounded 128-point line rule is unavailable")
    for (element_type,rule) in (
        (2,"Gauss6"),(4,"Gauss6"),(1,"CompositeGauss256"),
        (2,"CompositeGauss255"),(4,"CompositeGauss198"),
        (3,"Gauss2"),(34,"Gauss2"),(999,"Gauss2"),
        (2,"Gauss-1"),(2,"Gauss 2"),(2,"gauss2"),
        (2,"Gauss999999999999999999999999999999"))
        _rejects_argument(()->Tessella.API.mesh.get_integration_points(
            element_type,rule)) || error(
            "Tessella accepted invalid quadrature request " *
            "($element_type, $(repr(rule)))")
    end

    digest=bytes2hex(SHA.sha256(take!(stream)))
    digest=="a7e167c24bde2a871b1c9e4e4e5ae6c6bbbec0675ca2c56eda0602760e2888c2" || error(
        "mesh quadrature checksum changed to $digest")
    println("mesh-quadrature differential: Gmsh ",
            gmsh.GMSH_API_VERSION," fixed_types=",length(supported),
            " comparisons=",comparison_count," sha=",digest)
finally
    gmsh.finalize()
end
