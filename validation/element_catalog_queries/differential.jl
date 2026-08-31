# Differential oracle for fixed element-type and property lookup. This uses the
# locally installed Gmsh 4.15.2 Julia API and never starts the GUI.
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
    error("element-catalog-query differential: Gmsh Julia API not found; " *
          "set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION=="4.15.2" || error(
    "element-catalog-query differential requires Gmsh API 4.15.2, found " *
    gmsh.GMSH_API_VERSION)

const _CATALOG_FAMILY_NAMES=Dict(
    :pnt=>"Point",:lin=>"Line",:tri=>"Triangle",:qua=>"Quadrangle",
    :tet=>"Tetrahedron",:hex=>"Hexahedron",:pri=>"Prism",
    :pyr=>"Pyramid",:trih=>"Trihedron")
const _CATALOG_PRIMARY_NODES=Dict(
    :pnt=>1,:lin=>2,:tri=>3,:qua=>4,:tet=>4,
    :hex=>8,:pri=>6,:pyr=>5,:trih=>4)
const _CATALOG_PROPERTY_GAPS=Set([
    90,91,106,107,108,109,110,111,112,113,114,115,116,117,140])
const _CATALOG_SPECIAL_PROPERTIES=Dict(
    34=>("Polygon",2),35=>("Polyhedron",3),
    69=>("Polygon Border",2),133=>("Point Xfem",0),
    134=>("Line Xfem",1),135=>("Triangle Xfem",2),
    136=>("Tetrahedron Xfem",3))
const _CATALOG_SPECIAL_PROPERTY_ERRORS=(67,68,70,138,139)
const _CATALOG_SHA=
    "4c8e4086febda2bc3f482087d589eb157dbe57394f30781449a17a64e1b2a65c"

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
    catch
        return true
    end
end

function _catalog_write_i32!(stream,value)
    write(stream,htol(reinterpret(UInt32,Int32(value))))
end

function _catalog_sha()
    stream=IOBuffer()
    tags=vcat(sort!(collect(keys(MSH_CATALOG))),
              sort!(collect(keys(_CATALOG_SPECIAL_PROPERTIES))))
    for tag in tags
        name,dimension,order,num_nodes,coordinates,num_primary_nodes=
            Tessella.API.mesh.get_element_properties(tag)
        write(stream,htol(UInt64(tag)))
        bytes=codeunits(name)
        write(stream,htol(UInt64(length(bytes))))
        write(stream,bytes)
        for value in (dimension,order,num_nodes,num_primary_nodes)
            _catalog_write_i32!(stream,value)
        end
        write(stream,htol(UInt64(length(coordinates))))
        for value in coordinates
            write(stream,htol(reinterpret(UInt64,value)))
        end
    end
    return bytes2hex(SHA.sha256(take!(stream)))
end

try
    gmsh.initialize(String[],false)
    gmsh.option.setNumber("General.Terminal",0)
    Tessella.API.finalize()

    fixed_tags=sort!(collect(keys(MSH_CATALOG)))
    for tag in fixed_tags
        spec=msh_spec(tag)
        family_name=_CATALOG_FAMILY_NAMES[spec.family]
        gmsh_type=gmsh.model.mesh.getElementType(
            family_name,spec.order,spec.serendipity)
        gmsh_type==tag || error(
            "Gmsh type lookup for $family_name order $(spec.order) " *
            "serendipity=$(spec.serendipity) returned $gmsh_type, expected $tag")
        Tessella.API.mesh.get_element_type(
            family_name,spec.order,spec.serendipity)==Int32(tag) || error(
            "Tessella type lookup differs for Gmsh type $tag")

        tessella_properties=Tessella.API.mesh.get_element_properties(tag)
        if tag in _CATALOG_PROPERTY_GAPS
            _gmsh_rejects(()->gmsh.model.mesh.getElementProperties(tag)) ||
                error("Gmsh unexpectedly returned properties for gap type $tag")
            coordinate_dimension=max(spec.dim,1)
            expected_coordinates=vec(
                lagrange_nodes(tag)[1:coordinate_dimension,:])
            tessella_properties[2:4]==
                (Int32(spec.dim),Int32(spec.order),Int32(spec.nnodes)) || error(
                "Tessella property dimensions differ for gap type $tag")
            tessella_properties[5]==expected_coordinates || error(
                "Tessella local nodes differ for gap type $tag")
            tessella_properties[6]==
                Int32(_CATALOG_PRIMARY_NODES[spec.family]) || error(
                "Tessella primary-node count differs for gap type $tag")
        else
            gmsh_properties=gmsh.model.mesh.getElementProperties(tag)
            tessella_properties[1:4]==gmsh_properties[1:4] || error(
                "Tessella/Gmsh scalar properties differ for type $tag")
            tessella_properties[6]==gmsh_properties[6] || error(
                "Tessella/Gmsh primary-node count differs for type $tag")
            isapprox(tessella_properties[5],gmsh_properties[5];
                     atol=8eps(Float64),rtol=8eps(Float64)) || error(
                "Tessella/Gmsh local coordinates differ for type $tag")
        end
    end

    for (tag,(name,dimension)) in _CATALOG_SPECIAL_PROPERTIES
        gmsh_properties=gmsh.model.mesh.getElementProperties(tag)
        expected=(name,Int32(dimension),Int32(1),Int32(0),
                  Float64[],Int32(0))
        gmsh_properties==expected || error(
            "Gmsh special properties changed for type $tag")
        Tessella.API.mesh.get_element_properties(tag)==expected || error(
            "Tessella/Gmsh special properties differ for type $tag")
    end
    for tag in _CATALOG_SPECIAL_PROPERTY_ERRORS
        _gmsh_rejects(()->gmsh.model.mesh.getElementProperties(tag)) ||
            error("Gmsh unexpectedly returned special properties for type $tag")
        _rejects_argument(
            ()->Tessella.API.mesh.get_element_properties(tag)) || error(
            "Tessella did not reject unsupported special type $tag")
    end

    for (family,order,incomplete,expected) in (
        ("Triangle",1,true,2),("Line",10,true,66),
        ("Pyramid",1,true,7),("Point",-1,true,15),
        ("Trihedron",typemax(Int32),false,140))
        gmsh.model.mesh.getElementType(family,order,incomplete)==expected ||
            error("Gmsh fallback lookup changed for $family order $order")
        Tessella.API.mesh.get_element_type(family,order,incomplete)==expected ||
            error("Tessella fallback lookup differs for $family order $order")
    end

    gmsh.model.mesh.getElementType("Unknown",1,false)==0 || error(
        "Gmsh unknown-family sentinel changed")
    _rejects_argument(
        ()->Tessella.API.mesh.get_element_type("Unknown",1,false)) || error(
        "Tessella did not reject an unknown family")

    digest=_catalog_sha()
    digest==_CATALOG_SHA || error(
        "element catalog checksum changed to $digest")
    detached=Tessella.API.mesh.get_element_properties(4)
    detached[5][1]=99
    Tessella.API.mesh.get_element_properties(4)[5][1]==0 || error(
        "element-property coordinates alias catalog state")

    println("element-catalog-query differential: Gmsh ",
            gmsh.GMSH_API_VERSION," fixed_types=",length(fixed_tags),
            " direct_properties=",length(fixed_tags)-length(_CATALOG_PROPERTY_GAPS),
            " explicit_property_gaps=",length(_CATALOG_PROPERTY_GAPS),
            " special_properties=",length(_CATALOG_SPECIAL_PROPERTIES),
            " sha=",digest,
            " strict_unknown_family=true detached_coordinates=true")
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
