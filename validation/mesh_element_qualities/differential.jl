# Differential oracle for cached linear-simplex element qualities. This uses the
# locally installed Gmsh 4.15.2 Julia API and never starts the GUI.
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
    error("mesh-element-quality differential: Gmsh Julia API not found; " *
          "set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION=="4.15.2" || error(
    "mesh-element-quality differential requires Gmsh API 4.15.2, found " *
    gmsh.GMSH_API_VERSION)

const _QUALITY_COORDINATES=Float64[
    1 4   0 2 0.3   -1 1.4 -0.2   0.2 2.2 0.6 0.1   1 1 3 1;
   -2 2   0 0.2 3    0 0.7  2.1  -0.1 0.4 2.7 0.2  -1 2 -1 -1;
    0.5 -1 0 1 -0.4  0.5 1.2 0.8   0.3 0.1 0.8 3.4   0 0 0 2]
const _QUALITY_SEGMENTS=reshape(Int32[1,2],2,1)
const _QUALITY_TRIANGLES=Int32[3 6;4 7;5 8]
const _QUALITY_TETRAHEDRA=Int32[9 13;10 14;11 15;12 16]
const _GM_TO_DENSE=Dict(
    UInt64(101)=>UInt64(1),UInt64(201)=>UInt64(2),UInt64(202)=>UInt64(3),
    UInt64(301)=>UInt64(4),UInt64(302)=>UInt64(5))
const _ALL_QUALITY_NAMES=(
    "minDetJac","maxDetJac","minSJ","minSICN","minSIGE","gamma",
    "innerRadius","outerRadius","minIsotropy","angleShape","minEdge",
    "maxEdge","volume")
const _SEGMENT_QUALITY_NAMES=(
    "minSJ","minSICN","gamma","innerRadius","outerRadius","angleShape",
    "minEdge","maxEdge","volume")

function _quality_mesh()
    return Mesh(
        _QUALITY_COORDINATES;
        segs=_QUALITY_SEGMENTS,tris=_QUALITY_TRIANGLES,
        tets=_QUALITY_TETRAHEDRA)
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

function _write_quality_result!(stream,name,values)
    write(stream,codeunits(name))
    write(stream,UInt8(0))
    write(stream,htol(UInt64(length(values))))
    for value in values
        write(stream,htol(reinterpret(UInt64,Float64(value))))
    end
end

try
    gmsh.initialize(String[],false)
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.model.add("mesh-element-qualities")
    for (dimension,entity) in ((1,101),(2,201),(3,301))
        gmsh.model.addDiscreteEntity(dimension,entity)
    end
    gmsh.model.mesh.addNodes(
        3,301,UInt64.(1:size(_QUALITY_COORDINATES,2)),
        collect(vec(_QUALITY_COORDINATES)))
    gmsh.model.mesh.addElementsByType(
        101,1,UInt64[101],UInt64.(vec(_QUALITY_SEGMENTS)))
    gmsh.model.mesh.addElementsByType(
        201,2,UInt64[201,202],UInt64.(vec(_QUALITY_TRIANGLES)))
    gmsh.model.mesh.addElementsByType(
        301,4,UInt64[301,302],UInt64.(vec(_QUALITY_TETRAHEDRA)))

    Tessella.API.initialize()
    digest=try
        fixture=_quality_mesh()
        baseline=mesh_crc(fixture)
        lock(Tessella.API.STATE_LOCK) do
            Tessella.API._replace_mesh_cache_locked!(
                Tessella.API._copy_mesh(fixture))
        end

        stream=IOBuffer()
        gmsh_tags=UInt64[302,201,301,202,301]
        dense_tags=UInt64[_GM_TO_DENSE[tag] for tag in gmsh_tags]
        for quality in _ALL_QUALITY_NAMES
            gmsh_values=gmsh.model.mesh.getElementQualities(
                gmsh_tags,quality)
            tessella_values=Tessella.API.mesh.get_element_qualities(
                dense_tags,quality)
            length(gmsh_values)==length(tessella_values) || error(
                "$quality result lengths differ")
            for index in eachindex(gmsh_values,tessella_values)
                isapprox(tessella_values[index],gmsh_values[index];
                         atol=2.0e-13,rtol=2.0e-13) || error(
                    "$quality differs at request index $index: Tessella=" *
                    "$(tessella_values[index]), Gmsh=$(gmsh_values[index])")
            end
            _write_quality_result!(stream,quality,tessella_values)
        end

        for quality in _SEGMENT_QUALITY_NAMES
            gmsh_values=gmsh.model.mesh.getElementQualities(
                UInt64[101,101],quality)
            tessella_values=Tessella.API.mesh.get_element_qualities(
                UInt64[1,1],quality)
            tessella_values==gmsh_values || error(
                "segment $quality differs: Tessella=$tessella_values, " *
                "Gmsh=$gmsh_values")
            _write_quality_result!(stream,"segment:"*quality,tessella_values)
        end

        Tessella.API.mesh.get_element_qualities(
            UInt64[5],"volume")[1]<0 || error(
            "the inverted-tetrahedron fixture is not inverted")

        for quality in ("minDetJac","maxDetJac","minSIGE","minIsotropy")
            _rejects_argument(()->
                Tessella.API.mesh.get_element_qualities(
                    UInt64[2,1],quality)) || error(
                "Tessella accepted undefined segment quality $quality")
        end
        mesh_crc(Tessella.API.mesh.get())==baseline || error(
            "quality queries mutated the cached mesh")

        result=bytes2hex(SHA.sha256(take!(stream)))
        result=="4b32e86fee56ef55e7ad571c640155d88f40ceaff857186768554acb4250f3e0" ||
            error("mesh element-quality checksum changed to $result")
        result
    finally
        Tessella.API.finalize()
    end

    println("mesh-element-quality differential: Gmsh ",
            gmsh.GMSH_API_VERSION,
            " simplex_cases=5 mixed_order=true inverted_tet=true ",
            "qualities=13 segment_blockers=4 sha=",digest)
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
