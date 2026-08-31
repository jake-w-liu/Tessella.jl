# Differential oracle for cached linear-simplex reference maps and Jacobians.
# This uses the locally installed Gmsh 4.15.2 Julia API and never starts the GUI.
using Tessella
using SHA
using Tessella.MeshTypes: mesh_crc

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
    error("mesh-Jacobian differential: Gmsh Julia API not found; " *
          "set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION=="4.15.2" || error(
    "mesh-Jacobian differential requires Gmsh API 4.15.2, found " *
    gmsh.GMSH_API_VERSION)

const _JACOBIAN_COORDINATES=Float64[
    1 4  -3 -1  -2 1 0.5   2 3 -0.5   0 2 0 0   1 0 0 0;
   -2 2   1  4   0 0.2 3  -1 2  1     0 0 3 0   0 0 1 0;
    0.5 -1 2 3.5 0 1 -0.4  0.5 1 4    0 0 0 4   0 1 0 4]
const _JACOBIAN_SEGMENTS=Int32[1 3;2 4]
const _JACOBIAN_TRIANGLES=Int32[5 8;6 9;7 10]
const _JACOBIAN_TETRAHEDRA=Int32[11 15;12 16;13 17;14 18]
const _GM_TO_DENSE_JACOBIAN=Dict(
    UInt64(101)=>UInt64(1),UInt64(102)=>UInt64(2),
    UInt64(201)=>UInt64(3),UInt64(202)=>UInt64(4),
    UInt64(301)=>UInt64(5),UInt64(302)=>UInt64(6))

function _jacobian_mesh()
    return Mesh(
        _JACOBIAN_COORDINATES;
        segs=_JACOBIAN_SEGMENTS,tris=_JACOBIAN_TRIANGLES,
        tets=_JACOBIAN_TETRAHEDRA)
end

function _compare_values(label,gmsh_values,tessella_values)
    length(gmsh_values)==length(tessella_values) || error(
        "$label lengths differ: Gmsh=$(length(gmsh_values)), " *
        "Tessella=$(length(tessella_values))")
    for index in eachindex(gmsh_values,tessella_values)
        isapprox(tessella_values[index],gmsh_values[index];
                 atol=2.0e-13,rtol=2.0e-13) || error(
            "$label differs at index $index: Tessella=" *
            "$(tessella_values[index]), Gmsh=$(gmsh_values[index])")
    end
    return nothing
end

function _write_values!(stream,label,values)
    write(stream,codeunits(label));write(stream,UInt8(0))
    write(stream,htol(UInt64(length(values))))
    foreach(value->write(
        stream,htol(reinterpret(UInt64,Float64(value)))),values)
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

try
    gmsh.initialize(String[],false)
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.model.add("mesh-jacobians")
    for (dimension,entity) in ((1,101),(2,201),(3,301))
        gmsh.model.addDiscreteEntity(dimension,entity)
    end
    gmsh.model.mesh.addNodes(
        3,301,UInt64.(1:size(_JACOBIAN_COORDINATES,2)),
        collect(vec(_JACOBIAN_COORDINATES)))
    gmsh.model.mesh.addElementsByType(
        101,1,UInt64[101,102],UInt64.(vec(_JACOBIAN_SEGMENTS)))
    gmsh.model.mesh.addElementsByType(
        201,2,UInt64[201,202],UInt64.(vec(_JACOBIAN_TRIANGLES)))
    gmsh.model.mesh.addElementsByType(
        301,4,UInt64[301,302],UInt64.(vec(_JACOBIAN_TETRAHEDRA)))

    Tessella.API.initialize()
    digest=try
        fixture=_jacobian_mesh()
        baseline=mesh_crc(fixture)
        lock(Tessella.API.STATE_LOCK) do
            Tessella.API._replace_mesh_cache_locked!(
                Tessella.API._copy_mesh(fixture))
        end

        stream=IOBuffer()
        cases=(
            (1,101,Float64[-1,0,0, 0,0,0, 1,0,0]),
            (2,201,Float64[0,0,0, 0.2,0.3,0, 1,0,0]),
            (4,301,Float64[0,0,0, 0.2,0.3,0.1, 1,0,0]))
        for (element_type,entity,local_coordinates) in cases
            gmsh_result=gmsh.model.mesh.getJacobians(
                element_type,local_coordinates,entity)
            tessella_result=Tessella.API.mesh.get_jacobians(
                element_type,local_coordinates)
            for (name,gmsh_values,tessella_values) in zip(
                ("jacobians","determinants","coordinates"),
                gmsh_result,tessella_result)
                label="type-$element_type:$name"
                _compare_values(label,gmsh_values,tessella_values)
                _write_values!(stream,label,tessella_values)
            end
        end

        single_local=Float64[0.2,0.3,0.1, 0.1,0.2,0.3]
        for gmsh_tag in UInt64[101,202,301,302]
            dense_tag=_GM_TO_DENSE_JACOBIAN[gmsh_tag]
            gmsh_result=gmsh.model.mesh.getJacobian(
                gmsh_tag,single_local)
            tessella_result=Tessella.API.mesh.get_jacobian(
                dense_tag,single_local)
            for (name,gmsh_values,tessella_values) in zip(
                ("jacobians","determinants","coordinates"),
                gmsh_result,tessella_result)
                label="element-$gmsh_tag:$name"
                _compare_values(label,gmsh_values,tessella_values)
                _write_values!(stream,label,tessella_values)
            end
        end

        Tessella.API.mesh.get_jacobian(
            UInt64(6),Float64[0,0,0])[2][1]<0 || error(
            "inverted tetrahedron lost its signed determinant")
        Tessella.API.mesh.get_jacobians(
            3,Float64[0,0,0])==(Float64[],Float64[],Float64[]) || error(
            "known absent type did not return empty arrays")
        for malformed in (
            Float64[0],Float64[0,0],Float64[0,0,0,0],
            Float64[0,0,0,0,0])
            _rejects_argument(()->Tessella.API.mesh.get_jacobians(
                1,malformed)) || error(
                "Tessella accepted malformed local-coordinate length " *
                "$(length(malformed))")
        end
        for invalid in (Float64[NaN,0,0],Float64[Inf,0,0])
            _rejects_argument(()->Tessella.API.mesh.get_jacobians(
                1,invalid)) || error(
                "Tessella accepted non-finite local coordinates")
        end
        _rejects_argument(()->Tessella.API.mesh.get_jacobians(
            1,Float64[0,0,0],101)) || error(
            "Tessella accepted entity filtering without classification")
        _rejects_argument(()->Tessella.API.mesh.get_jacobians(
            1,Float64[0,0,0],-1,1,2)) || error(
            "Tessella accepted detached-array task partitioning")
        mesh_crc(Tessella.API.mesh.get())==baseline || error(
            "Jacobian queries mutated the cached mesh")

        result=bytes2hex(SHA.sha256(take!(stream)))
        result=="85803b9ce40eed41e887957e73a58a1a3a1c1251a72ef3c1e5315b165772d631" || error(
            "mesh Jacobian checksum changed to $result")
        result
    finally
        Tessella.API.finalize()
    end

    println("mesh-Jacobian differential: Gmsh ",
            gmsh.GMSH_API_VERSION,
            " elements=6 evaluation_points=3 bulk_types=3 " *
            "single_elements=4 inverted_tet=true sha=",digest)
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
