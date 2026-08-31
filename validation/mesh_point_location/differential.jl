# Differential oracle for cached linear-simplex point location and reference
# coordinates. This uses the locally installed Gmsh 4.15.2 Julia API and never
# starts the GUI.
using Tessella
using SHA

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
    error("mesh-point-location differential: Gmsh Julia API not found; " *
          "set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION=="4.15.2" || error(
    "mesh-point-location differential requires Gmsh API 4.15.2, found " *
    gmsh.GMSH_API_VERSION)

const _LOCATION_COORDINATES=Float64[0 1 0 0;
                                    0 0 1 0;
                                    0 0 0 1]
const _GM_TO_DENSE_ELEMENT=Dict(
    UInt64(100)=>UInt64(1),UInt64(300)=>UInt64(2),UInt64(200)=>UInt64(3))
const _GM_TO_DENSE_NODE=Dict(
    UInt64(10)=>UInt64(1),UInt64(20)=>UInt64(2),
    UInt64(30)=>UInt64(3),UInt64(40)=>UInt64(4))

function _location_mesh()
    return Mesh(
        _LOCATION_COORDINATES;
        segs=reshape(Int32[1,2],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[1,2,3,4],4,1))
end

function _normalized_element_tags(tags)
    all(tag->haskey(_GM_TO_DENSE_ELEMENT,tag),tags) || error(
        "Gmsh point location returned an unknown fixture element tag")
    return UInt64[_GM_TO_DENSE_ELEMENT[tag] for tag in tags]
end

function _normalized_node_tags(tags)
    all(tag->haskey(_GM_TO_DENSE_NODE,tag),tags) || error(
        "Gmsh point location returned an unknown fixture node tag")
    return UInt64[_GM_TO_DENSE_NODE[tag] for tag in tags]
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

function _gmsh_rejects(f::Function)
    try
        f()
        return false
    catch err
        err isa InterruptException && rethrow()
        return true
    end
end

function _write_location_result!(stream,tags,coordinates)
    write(stream,htol(UInt64(length(tags))))
    foreach(tag->write(stream,htol(UInt64(tag))),tags)
    for value in coordinates
        write(stream,htol(reinterpret(UInt64,Float64(value))))
    end
end

try
    gmsh.initialize(String[],false)
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.model.add("mesh-point-location")
    for (dimension,entity) in ((1,101),(2,201),(3,301))
        gmsh.model.addDiscreteEntity(dimension,entity)
    end
    gmsh.model.mesh.addNodes(
        3,301,UInt64[10,20,30,40],collect(vec(_LOCATION_COORDINATES)))
    gmsh.model.mesh.addElementsByType(
        101,1,UInt64[100],UInt64[10,20])
    gmsh.model.mesh.addElementsByType(
        201,2,UInt64[300],UInt64[10,20,30])
    gmsh.model.mesh.addElementsByType(
        301,4,UInt64[200],UInt64[10,20,30,40])

    Tessella.API.initialize()
    digest=try
        fixture=_location_mesh()
        lock(Tessella.API.STATE_LOCK) do
            Tessella.API._replace_mesh_cache_locked!(
                Tessella.API._copy_mesh(fixture))
        end

        stream=IOBuffer()
        cases=(
            (0.25,0.0,0.0,-1,true),
            (0.2,0.3,0.0,-1,true),
            (0.1,0.2,0.3,-1,true),
            (0.25,0.0,0.0,1,true),
            (0.25,0.0,0.0,2,true),
            (0.25,0.0,0.0,3,true),
        )
        for (x,y,z,dimension,strict) in cases
            gmsh_tags=gmsh.model.mesh.getElementsByCoordinates(
                x,y,z,dimension,strict)
            normalized_tags=_normalized_element_tags(gmsh_tags)
            tessella_tags=Tessella.API.mesh.get_elements_by_coordinates(
                x,y,z,dimension,strict)
            tessella_tags==normalized_tags || error(
                "all-element location differs at ($x, $y, $z), dim=$dimension")

            gmsh_first=gmsh.model.mesh.getElementByCoordinates(
                x,y,z,dimension,strict)
            tessella_first=Tessella.API.mesh.get_element_by_coordinates(
                x,y,z,dimension,strict)
            tessella_first[1]==_GM_TO_DENSE_ELEMENT[gmsh_first[1]] || error(
                "first element tag differs at ($x, $y, $z), dim=$dimension")
            tessella_first[2]==gmsh_first[2] || error(
                "element type differs at ($x, $y, $z), dim=$dimension")
            tessella_first[3]==_normalized_node_tags(gmsh_first[3]) || error(
                "element nodes differ at ($x, $y, $z), dim=$dimension")
            isapprox(collect(tessella_first[4:6]),collect(gmsh_first[4:6]);
                     atol=32eps(Float64),rtol=32eps(Float64)) || error(
                "reference coordinates differ at ($x, $y, $z), dim=$dimension")
            _write_location_result!(
                stream,tessella_tags,tessella_first[4:6])
        end

        local_cases=(
            (UInt64(100),UInt64(1),0.25,0.0,0.0),
            (UInt64(300),UInt64(2),0.2,0.3,0.0),
            (UInt64(200),UInt64(3),0.1,0.2,0.3),
            (UInt64(200),UInt64(3),-0.25,0.5,1.25),
        )
        for (gmsh_tag,tessella_tag,x,y,z) in local_cases
            gmsh_coordinates=gmsh.model.mesh.getLocalCoordinatesInElement(
                gmsh_tag,x,y,z)
            tessella_coordinates=
                Tessella.API.mesh.get_local_coordinates_in_element(
                    tessella_tag,x,y,z)
            isapprox(collect(tessella_coordinates),collect(gmsh_coordinates);
                     atol=32eps(Float64),rtol=32eps(Float64)) || error(
                "local coordinates differ for dense element $tessella_tag")
            _write_location_result!(
                stream,UInt64[tessella_tag],tessella_coordinates)
        end

        outside=(-5.0e-6,0.2,0.2)
        _gmsh_rejects(()->gmsh.model.mesh.getElementsByCoordinates(
            outside...,3,true)) || error(
                "Gmsh strict tetrahedron tolerance changed")
        _rejects_argument(
            ()->Tessella.API.mesh.get_elements_by_coordinates(
                outside...,3,true)) || error(
                    "Tessella strict tetrahedron search accepted relaxed point")
        _normalized_element_tags(
            gmsh.model.mesh.getElementsByCoordinates(outside...,3,false))==
            UInt64[3] || error("Gmsh relaxed tetrahedron tolerance changed")
        Tessella.API.mesh.get_elements_by_coordinates(outside...,3,false)==
            UInt64[3] || error(
                "Tessella relaxed tetrahedron search did not match Gmsh")

        far=(2.0,2.0,2.0)
        _gmsh_rejects(()->gmsh.model.mesh.getElementsByCoordinates(
            far...,-1,false)) || error("Gmsh no-match behavior changed")
        _rejects_argument(
            ()->Tessella.API.mesh.get_elements_by_coordinates(
                far...,-1,false)) || error(
                    "Tessella returned an element for a no-match point")

        # Gmsh's line Newton inversion and triangle dominant-plane inversion can
        # return off-span artifacts. Tessella deliberately returns stable orthogonal
        # projections with unused coordinates fixed at zero.
        Tessella.API.mesh.get_local_coordinates_in_element(
            1,0.25,2.0,3.0)==(-0.5,0.0,0.0) || error(
                "Tessella off-span line projection changed")
        Tessella.API.mesh.get_local_coordinates_in_element(
            2,0.2,0.3,5.0)==(0.2,0.3,0.0) || error(
                "Tessella off-span triangle projection changed")

        result=bytes2hex(SHA.sha256(take!(stream)))
        result=="afdef11844fca534f32a62914a8c21432c282788c6c0504c74cfb2a2205dbdb7" || error(
            "mesh point-location checksum changed to $result")
        result
    finally
        Tessella.API.finalize()
    end

    println("mesh-point-location differential: Gmsh ",
            gmsh.GMSH_API_VERSION,
            " exact_cases=6 local_cases=4 strict_tolerance=1e-6 ",
            "relaxed_tolerance=true deterministic_dimension_order=true ",
            "orthogonal_off_span=true sha=",digest)
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
