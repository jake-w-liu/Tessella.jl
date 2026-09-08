# Differential oracle for fixed-family first-order nodal bases plus simplex
# hierarchical bases, orientations, and keys.
# This uses the locally installed Gmsh 4.15.2 Julia API and never starts the GUI.
using Tessella
using Random
using SHA
using Tessella.Elements: MSH_CATALOG
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
    error("mesh-function-space differential: Gmsh Julia API not found; " *
          "set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION=="4.15.2" || error(
    "mesh-function-space differential requires Gmsh API 4.15.2, found " *
    gmsh.GMSH_API_VERSION)

const _FUNCTION_COORDINATES=Float64[
    0 2 0 0;
    0 0 3 0;
    0 0 0 4]
const _FUNCTION_SEGMENTS=reshape(Int32[2,1],2,1)
const _FUNCTION_TRIANGLES=reshape(Int32[3,1,2],3,1)
const _FUNCTION_TETRAHEDRA=reshape(Int32[4,2,1,3],4,1)
const _FUNCTION_SPACES=(
    "Lagrange","IsoParametric","Lagrange1",
    "GradLagrange","GradIsoParametric","GradLagrange1",
    "H1Legendre1","GradH1Legendre1",
    "HcurlLegendre0","CurlHcurlLegendre0")
const _HIERARCHICAL_SPACES=Set((
    "H1Legendre1","GradH1Legendre1",
    "HcurlLegendre0","CurlHcurlLegendre0"))
const _FIXED_NODAL_SPACES=("Lagrange1","GradLagrange1")
const _FIXED_NODAL_TYPES=sort!([
    Int(element_type)
    for (element_type,spec) in MSH_CATALOG
    if spec.family in (:pnt,:lin,:tri,:qua,:tet,:hex,:pri,:pyr)])
const _LINEAR_FAMILY_TYPES=(15,3,5,6,7)
const _LINEAR_FAMILY_SPACES=(
    "Lagrange","IsoParametric","GradLagrange","GradIsoParametric")
const _MAX_ABS_DIFFERENCE=Ref(0.0)

function _function_mesh()
    return Mesh(
        _FUNCTION_COORDINATES;
        segs=_FUNCTION_SEGMENTS,tris=_FUNCTION_TRIANGLES,
        tets=_FUNCTION_TETRAHEDRA)
end

function _install_function_mesh!(mesh)
    lock(Tessella.API.STATE_LOCK) do
        Tessella.API._replace_mesh_cache_locked!(
            Tessella.API._copy_mesh(mesh))
    end
    return nothing
end

function _build_gmsh_function_model(name)
    gmsh.model.add(name)
    for (dimension,tag) in ((1,1),(2,2),(3,3))
        gmsh.model.addDiscreteEntity(dimension,tag)
    end
    gmsh.model.mesh.addNodes(
        3,3,UInt64[1,2,3,4],collect(vec(_FUNCTION_COORDINATES)))
    gmsh.model.mesh.addElementsByType(
        1,1,UInt64[101],UInt64.(vec(_FUNCTION_SEGMENTS)))
    gmsh.model.mesh.addElementsByType(
        2,2,UInt64[102],UInt64.(vec(_FUNCTION_TRIANGLES)))
    gmsh.model.mesh.addElementsByType(
        3,4,UInt64[103],UInt64.(vec(_FUNCTION_TETRAHEDRA)))
    return nothing
end

function _lexicographic_permutations(count)
    result=Vector{Vector{Int32}}()
    function visit(prefix,remaining)
        if isempty(remaining)
            push!(result,Int32.(prefix))
            return
        end
        for index in eachindex(remaining)
            visit([prefix;remaining[index]],
                  remaining[[i for i in eachindex(remaining) if i!=index]])
        end
    end
    visit(Int[],collect(1:count))
    return result
end

function _build_orientation_models()
    segment_cells=hcat(_lexicographic_permutations(2)...)
    triangle_cells=hcat(_lexicographic_permutations(3)...)
    tetrahedron_cells=hcat(_lexicographic_permutations(4)...)
    gmsh.model.add("mesh-function-orientations")
    for (dimension,tag) in ((1,1),(2,2),(3,3))
        gmsh.model.addDiscreteEntity(dimension,tag)
    end
    gmsh.model.mesh.addNodes(
        3,3,UInt64[1,2,3,4],collect(vec(_FUNCTION_COORDINATES)))
    gmsh.model.mesh.addElementsByType(
        1,1,UInt64.(101:102),UInt64.(vec(segment_cells)))
    gmsh.model.mesh.addElementsByType(
        2,2,UInt64.(201:206),UInt64.(vec(triangle_cells)))
    gmsh.model.mesh.addElementsByType(
        3,4,UInt64.(301:324),UInt64.(vec(tetrahedron_cells)))
    tessella=Mesh(
        _FUNCTION_COORDINATES;
        segs=segment_cells,tris=triangle_cells,tets=tetrahedron_cells)
    return tessella
end

function _compare_float(label,gmsh_values,tessella_values)
    length(gmsh_values)==length(tessella_values) || error(
        "$label lengths differ: Gmsh=$(length(gmsh_values)), " *
        "Tessella=$(length(tessella_values))")
    for index in eachindex(gmsh_values,tessella_values)
        difference=abs(tessella_values[index]-gmsh_values[index])
        isapprox(tessella_values[index],gmsh_values[index];
                 atol=8e-14,rtol=8e-14) || error(
            "$label differs at $index: Gmsh=$(gmsh_values[index]), " *
            "Tessella=$(tessella_values[index])")
        _MAX_ABS_DIFFERENCE[]=max(_MAX_ABS_DIFFERENCE[],difference)
    end
    return nothing
end

function _compare_exact(label,gmsh_values,tessella_values)
    gmsh_values==tessella_values || error(
        "$label differs: Gmsh=$(repr(gmsh_values)), " *
        "Tessella=$(repr(tessella_values))")
    return nothing
end

function _write_ints!(stream,label,values)
    write(stream,codeunits(label));write(stream,UInt8(0))
    write(stream,htol(UInt64(length(values))))
    foreach(value->write(stream,htol(Int64(value))),values)
    return nothing
end

function _write_floats!(stream,label,values)
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

# Mutant analysis:
# - Reordering tetrahedron edges is rejected by the exact sequential key arrays.
# - Using parity instead of lexicographic permutation rank is rejected by the
#   reversed segment/triangle/tetrahedron orientation values 1, 4, and 20.
# - Flattening point-major before orientation-major is rejected by multi-point,
#   all-orientation full-array comparison against Gmsh.
# - Completing the whole edge cache during a type-specific key query is rejected
#   by the intermediate one-edge and four-edge catalog checks.
# - Dispatching explicit order-one nodal functions by interpolation order is
#   rejected across every higher-order, serendipity, and order-zero fixed type.
# - Permuting quadrangle, hexahedron, prism, or pyramid vertices is rejected by
#   sequential full-array comparisons at asymmetric evaluation points.
# - Replacing the pyramid's rational basis by a polynomial tensor basis is
#   rejected at interior, near-apex, and exact-apex points and by its gradients.
# - Rounding the pyramid height before a cancelling or overflowing factor is
#   rejected by the exact-rational extreme-coordinate oracle.
# - Omitting a prism line factor or a tensor derivative factor is rejected by
#   every non-simplex gradient comparison.

try
    gmsh.initialize(String[],false)
    gmsh.option.setNumber("General.Terminal",0)
    _build_gmsh_function_model("mesh-function-spaces")
    Tessella.API.initialize()
    digest,nodal_point_count=try
        fixture=_function_mesh()
        baseline=mesh_crc(fixture)
        _install_function_mesh!(fixture)
        stream=IOBuffer()

        rng=Xoshiro(0x4d65736846756e63)
        local_coordinates=Float64[
            -0.8,0.2,0.7,
             0.0,0.0,0.0,
             0.2,0.3,0.1,
             1.0,0.0,0.0]
        append!(local_coordinates,2rand(rng,Float64,36).-1)
        for element_type in (1,2,4),space in _FUNCTION_SPACES
            gmsh_basis=gmsh.model.mesh.getBasisFunctions(
                element_type,local_coordinates,space)
            tessella_basis=Tessella.API.mesh.get_basis_functions(
                element_type,local_coordinates,space)
            _compare_exact(
                "type-$element_type $space components",
                gmsh_basis[1],tessella_basis[1])
            _compare_float(
                "type-$element_type $space basis",
                gmsh_basis[2],tessella_basis[2])
            _compare_exact(
                "type-$element_type $space orientation count",
                gmsh_basis[3],tessella_basis[3])
            _compare_exact(
                "type-$element_type $space number of orientations",
                gmsh.model.mesh.getNumberOfOrientations(element_type,space),
                Tessella.API.mesh.get_number_of_orientations(
                    element_type,space))
            _compare_exact(
                "type-$element_type $space number of keys",
                gmsh.model.mesh.getNumberOfKeys(element_type,space),
                Tessella.API.mesh.get_number_of_keys(element_type,space))

            orientation_count=Int(gmsh_basis[3])
            wanted=if space in _HIERARCHICAL_SPACES
                orientation_count==2 ? Int32[1,0] :
                    Int32[orientation_count-1,0,orientation_count÷2]
            else
                Int32[0]
            end
            gmsh_selected=gmsh.model.mesh.getBasisFunctions(
                element_type,local_coordinates,space,wanted)
            tessella_selected=Tessella.API.mesh.get_basis_functions(
                element_type,local_coordinates,space,wanted)
            _compare_float(
                "type-$element_type $space selected basis",
                gmsh_selected[2],tessella_selected[2])
            _write_ints!(
                stream,"type-$element_type:$space:meta",
                (tessella_basis[1],tessella_basis[3]))
            _write_floats!(
                stream,"type-$element_type:$space:basis",tessella_basis[2])
            _write_floats!(
                stream,"type-$element_type:$space:selected",
                tessella_selected[2])
        end

        nodal_coordinates=Float64[
            0.13,-0.21,0.37,
            -0.4,0.7,-0.2,
            0.0,0.0,1.0,
            0.0,0.0,prevfloat(1.0)]
        append!(nodal_coordinates,2rand(rng,Float64,48).-1)
        for element_type in _FIXED_NODAL_TYPES,
            space in _FIXED_NODAL_SPACES
            gmsh_basis=gmsh.model.mesh.getBasisFunctions(
                element_type,nodal_coordinates,space)
            tessella_basis=Tessella.API.mesh.get_basis_functions(
                element_type,nodal_coordinates,space)
            label="type-$element_type $space fixed-family"
            _compare_exact(
                "$label components",gmsh_basis[1],tessella_basis[1])
            _compare_float(
                "$label basis",gmsh_basis[2],tessella_basis[2])
            _compare_exact(
                "$label orientation count",gmsh_basis[3],tessella_basis[3])
            _compare_float(
                "$label selected basis",
                gmsh.model.mesh.getBasisFunctions(
                    element_type,nodal_coordinates,space,Int32[0])[2],
                Tessella.API.mesh.get_basis_functions(
                    element_type,nodal_coordinates,space,Int32[0])[2])
            gmsh_orientations=
                gmsh.model.mesh.getNumberOfOrientations(element_type,space)
            tessella_orientations=
                Tessella.API.mesh.get_number_of_orientations(
                    element_type,space)
            _compare_exact(
                "$label number of orientations",
                gmsh_orientations,tessella_orientations)
            gmsh_key_count=gmsh.model.mesh.getNumberOfKeys(
                element_type,space)
            tessella_key_count=Tessella.API.mesh.get_number_of_keys(
                element_type,space)
            _compare_exact(
                "$label number of keys",gmsh_key_count,tessella_key_count)
            type_keys=zeros(Int32,gmsh_key_count)
            entity_keys=UInt64.(1:gmsh_key_count)
            _compare_exact(
                "$label key information",
                gmsh.model.mesh.getKeysInformation(
                    type_keys,entity_keys,element_type,space),
                Tessella.API.mesh.get_keys_information(
                    type_keys,entity_keys,element_type,space))
            _write_ints!(
                stream,"$label:meta",
                (tessella_basis[1],tessella_basis[3],
                 tessella_orientations,tessella_key_count))
            _write_floats!(stream,"$label:basis",tessella_basis[2])
        end

        for element_type in _LINEAR_FAMILY_TYPES,
            space in _LINEAR_FAMILY_SPACES
            gmsh_basis=gmsh.model.mesh.getBasisFunctions(
                element_type,nodal_coordinates,space)
            tessella_basis=Tessella.API.mesh.get_basis_functions(
                element_type,nodal_coordinates,space)
            label="type-$element_type $space linear-family"
            _compare_exact(
                "$label components",gmsh_basis[1],tessella_basis[1])
            _compare_float(
                "$label basis",gmsh_basis[2],tessella_basis[2])
            _compare_exact(
                "$label orientation count",gmsh_basis[3],tessella_basis[3])
            gmsh_orientations=
                gmsh.model.mesh.getNumberOfOrientations(element_type,space)
            tessella_orientations=
                Tessella.API.mesh.get_number_of_orientations(
                    element_type,space)
            _compare_exact(
                "$label number of orientations",
                gmsh_orientations,tessella_orientations)
            gmsh_key_count=gmsh.model.mesh.getNumberOfKeys(
                element_type,space)
            tessella_key_count=Tessella.API.mesh.get_number_of_keys(
                element_type,space)
            _compare_exact(
                "$label number of keys",gmsh_key_count,tessella_key_count)
            type_keys=zeros(Int32,gmsh_key_count)
            entity_keys=UInt64.(1:gmsh_key_count)
            _compare_exact(
                "$label key information",
                gmsh.model.mesh.getKeysInformation(
                    type_keys,entity_keys,element_type,space),
                Tessella.API.mesh.get_keys_information(
                    type_keys,entity_keys,element_type,space))
            _write_ints!(
                stream,"$label:meta",
                (tessella_basis[1],tessella_basis[3],
                 tessella_orientations,tessella_key_count))
            _write_floats!(stream,"$label:basis",tessella_basis[2])
        end

        gmsh_element_tags=Dict(1=>UInt64(101),2=>UInt64(102),4=>UInt64(103))
        dense_element_tags=Dict(1=>UInt64(1),2=>UInt64(2),4=>UInt64(3))
        for element_type in (1,2,4),space in _FUNCTION_SPACES
            gmsh_orientation=gmsh.model.mesh.getBasisFunctionsOrientation(
                element_type,space)
            tessella_orientation=
                Tessella.API.mesh.get_basis_functions_orientation(
                    element_type,space)
            _compare_exact(
                "type-$element_type $space orientations",
                gmsh_orientation,tessella_orientation)
            gmsh_single_orientation=
                gmsh.model.mesh.getBasisFunctionsOrientationForElement(
                    gmsh_element_tags[element_type],space)
            tessella_single_orientation=
                Tessella.API.mesh.get_basis_functions_orientation_for_element(
                    dense_element_tags[element_type],space)
            _compare_exact(
                "type-$element_type $space single orientation",
                gmsh_single_orientation,tessella_single_orientation)

            gmsh_keys=gmsh.model.mesh.getKeys(
                element_type,space,-1,true)
            tessella_keys=Tessella.API.mesh.get_keys(
                element_type,space,-1,true)
            _compare_exact(
                "type-$element_type $space type keys",
                gmsh_keys[1],tessella_keys[1])
            _compare_exact(
                "type-$element_type $space entity keys",
                gmsh_keys[2],tessella_keys[2])
            _compare_float(
                "type-$element_type $space key coordinates",
                gmsh_keys[3],tessella_keys[3])
            gmsh_info=gmsh.model.mesh.getKeysInformation(
                gmsh_keys[1],gmsh_keys[2],element_type,space)
            tessella_info=Tessella.API.mesh.get_keys_information(
                tessella_keys[1],tessella_keys[2],element_type,space)
            _compare_exact(
                "type-$element_type $space key information",
                gmsh_info,tessella_info)
            _write_ints!(
                stream,"type-$element_type:$space:type-keys",
                tessella_keys[1])
            _write_ints!(
                stream,"type-$element_type:$space:entity-keys",
                tessella_keys[2])
            _write_floats!(
                stream,"type-$element_type:$space:key-coordinates",
                tessella_keys[3])
        end

        all_edges=Tessella.API.mesh.get_all_edges()
        _compare_exact("completed edge tags",UInt64[1,2,3,4,5,6],all_edges[1])
        _compare_exact(
            "completed edge nodes",
            UInt64[2,1, 3,1, 2,3, 4,2, 1,4, 3,4],all_edges[2])
        mesh_crc(Tessella.API.mesh.get())==baseline || error(
            "basis or key queries mutated the cached mesh")

        gmsh.clear()
        orientation_mesh=_build_orientation_models()
        _install_function_mesh!(orientation_mesh)
        for (element_type,count) in ((1,2),(2,6),(4,24)),space in
            ("H1Legendre1","HcurlLegendre0")
            expected=Int32.(0:count-1)
            gmsh_orientations=
                gmsh.model.mesh.getBasisFunctionsOrientation(
                    element_type,space)
            tessella_orientations=
                Tessella.API.mesh.get_basis_functions_orientation(
                    element_type,space)
            _compare_exact(
                "type-$element_type $space exhaustive orientations",
                expected,gmsh_orientations)
            _compare_exact(
                "type-$element_type $space exhaustive Tessella orientations",
                expected,tessella_orientations)
            _write_ints!(
                stream,"type-$element_type:$space:orientations",
                tessella_orientations)
        end

        gmsh.clear()
        _build_gmsh_function_model("mesh-function-single")
        _install_function_mesh!(fixture)
        gmsh_single=gmsh.model.mesh.getKeysForElement(
            UInt64(103),"HcurlLegendre0",true)
        tessella_single=Tessella.API.mesh.get_keys_for_element(
            UInt64(3),"HcurlLegendre0",true)
        _compare_exact("single tetrahedron type keys",
                       gmsh_single[1],tessella_single[1])
        _compare_exact("single tetrahedron entity keys",
                       gmsh_single[2],tessella_single[2])
        _compare_float("single tetrahedron key coordinates",
                       gmsh_single[3],tessella_single[3])
        _compare_exact(
            "single tetrahedron edge catalog",
            (UInt64[1,2,3,4,5,6],
             UInt64[4,2, 2,1, 1,4, 3,4, 3,1, 3,2]),
            Tessella.API.mesh.get_all_edges())
        gmsh_triangle=gmsh.model.mesh.getKeys(
            2,"HcurlLegendre0",-1,false)
        tessella_triangle=Tessella.API.mesh.get_keys(
            2,"HcurlLegendre0",-1,false)
        _compare_exact("post-single triangle keys",
                       gmsh_triangle,tessella_triangle)

        before=Tessella.API.mesh.get_all_edges()
        for invalid in (
            ()->Tessella.API.mesh.get_basis_functions(
                4,Float64[0,0],"Lagrange"),
            ()->Tessella.API.mesh.get_basis_functions(
                4,Float64[0,0,0],"HcurlLegendre1"),
            ()->Tessella.API.mesh.get_basis_functions(
                4,Float64[0,0,0],"HcurlLegendre0",Int32[24]),
            ()->Tessella.API.mesh.get_keys(
                4,"HcurlLegendre0",-1,1),
            ()->Tessella.API.mesh.get_keys(
                4,"HcurlLegendre0",3,true),
            ()->Tessella.API.mesh.get_basis_functions_orientation(
                4,"HcurlLegendre0",-1,1,2),
            ()->Tessella.API.mesh.get_keys_information(
                Int32[1],UInt64[1],4,"HcurlLegendre0"),
            ()->Tessella.API.mesh.get_basis_functions(
                10,Float64[0,0,0],"Lagrange"),
            ()->Tessella.API.mesh.get_basis_functions(
                3,Float64[0,0,0],"H1Legendre1"),
            ()->Tessella.API.mesh.get_basis_functions(
                140,Float64[0,0,0],"Lagrange1"),
        )
            _rejects_argument(invalid) || error(
                "Tessella accepted an invalid function-space query")
            Tessella.API.mesh.get_all_edges()==before || error(
                "a rejected function-space query changed the edge catalog")
            mesh_crc(Tessella.API.mesh.get())==baseline || error(
                "a rejected function-space query changed the mesh")
        end

        result=bytes2hex(SHA.sha256(take!(stream)))
        result=="26d28400c434fa81836b6c0e83cf8afeabc8f812c86f6176d464ac7520816612" ||
            error("mesh function-space checksum changed to $result")
        result,length(nodal_coordinates)÷3
    finally
        Tessella.API.finalize()
    end
    println("mesh-function-space differential: Gmsh ",
            gmsh.GMSH_API_VERSION,
            " basis_cases=30 orientation_cases=30 key_cases=30 " *
            "fixed_nodal_types=",length(_FIXED_NODAL_TYPES),
            " fixed_nodal_cases=",
            length(_FIXED_NODAL_TYPES)*length(_FIXED_NODAL_SPACES),
            " fixed_nodal_points=",nodal_point_count,
            " linear_family_alias_cases=",
            length(_LINEAR_FAMILY_TYPES)*length(_LINEAR_FAMILY_SPACES),
            " all_orientations=32 max_abs_difference=",
            _MAX_ABS_DIFFERENCE[],
            " lazy_edge_catalog=true sha=",digest)
finally
    gmsh.finalize()
end
