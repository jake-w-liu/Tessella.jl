# Differential oracle for whole-mesh affine coordinate transformation. This uses
# the locally installed Gmsh 4.15.2 Julia API and never starts the GUI.
using Tessella

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
    error("mesh-affine-transform differential: Gmsh Julia API not found; " *
          "set GMSH_JULIA_API")
end

include(_gmsh_binding())
gmsh.GMSH_API_VERSION=="4.15.2" || error(
    "mesh-affine-transform differential requires Gmsh API 4.15.2, found " *
    gmsh.GMSH_API_VERSION)

const _SOURCE_COORDS=Float64[0 1 0 0;
                             0 0 1 0;
                             0 0 0 1]
const _SOURCE_TETS=reshape(Int32[1,2,3,4],4,1)
const _POSITIVE12=Float64[2,1,0,3,
                          0,1,0,-2,
                          0,0,1,4]
const _POSITIVE16=vcat(_POSITIVE12,[0.0,0.0,0.0,1.0])
const _REFLECTION12=Float64[-1,0,0,2,
                            0,1,0,3,
                            0,0,1,4]
const _SINGULAR12=Float64[0,0,0,0,
                          0,1,0,0,
                          0,0,1,0]

function _source_mesh()
    Tessella.MeshTypes.Mesh(_SOURCE_COORDS;tets=_SOURCE_TETS)
end

function _add_gmsh_fixture(name::String)
    gmsh.clear()
    gmsh.model.add(name)
    gmsh.model.addDiscreteEntity(3,301)
    gmsh.model.mesh.addNodes(
        3,301,UInt64[1,2,3,4],collect(vec(_SOURCE_COORDS)))
    gmsh.model.mesh.addElementsByType(
        301,4,UInt64[1],UInt64[1,2,3,4])
    return nothing
end

function _gmsh_coordinates()
    tags,flat,_=gmsh.model.mesh.getNodes()
    length(tags)==4 || error("Gmsh returned $(length(tags)) nodes, expected 4")
    length(flat)==12 || error("Gmsh returned malformed node coordinates")
    result=Matrix{Float64}(undef,3,4)
    columns=Dict(Int(tag)=>index for (index,tag) in enumerate(tags))
    Set(keys(columns))==Set(1:4) || error("Gmsh changed the fixture node tags")
    for tag in 1:4
        source=columns[tag]
        result[:,tag].=flat[3source-2:3source]
    end
    return result
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

    source=_source_mesh()
    positive_matrix=[2.0 1.0 0.0;
                     0.0 1.0 0.0;
                     0.0 0.0 1.0]
    positive_expected=Tessella.affine_transform(
        source,positive_matrix;translation=(3.0,-2.0,4.0))

    _add_gmsh_fixture("mesh-affine-positive-12")
    gmsh.model.mesh.affineTransform(_POSITIVE12)
    gmsh_positive12=_gmsh_coordinates()
    positive_error=maximum(abs.(gmsh_positive12.-positive_expected.coords))
    positive_error<=8eps(Float64)*maximum(abs,positive_expected.coords) || error(
        "Gmsh/Tessella positive affine coordinate difference is $positive_error")
    _,positive_connectivity=gmsh.model.mesh.getElementsByType(4,301)
    positive_connectivity==UInt64[1,2,3,4] || error(
        "Gmsh changed positive-affine tetrahedron connectivity")

    _add_gmsh_fixture("mesh-affine-positive-16")
    gmsh.model.mesh.affineTransform(_POSITIVE16)
    gmsh_positive16=_gmsh_coordinates()
    gmsh_positive16==gmsh_positive12 || error(
        "Gmsh 12- and 16-entry affine transforms differ")

    _add_gmsh_fixture("mesh-affine-reflection")
    gmsh.model.mesh.affineTransform(_REFLECTION12)
    gmsh_reflected=_gmsh_coordinates()
    reflected_expected=Tessella.affine_transform(
        source,[-1.0 0.0 0.0;0.0 1.0 0.0;0.0 0.0 1.0];
        translation=(2.0,3.0,4.0))
    reflection_error=maximum(abs.(gmsh_reflected.-reflected_expected.coords))
    reflection_error<=8eps(Float64)*maximum(abs,reflected_expected.coords) || error(
        "Gmsh/Tessella reflected coordinate difference is $reflection_error")
    _,gmsh_reflected_tet=gmsh.model.mesh.getElementsByType(4,301)
    gmsh_reflected_tet==UInt64[1,2,3,4] || error(
        "Gmsh unexpectedly rewound reflected tetrahedron connectivity")
    reflected_expected.tets[:,1]==Int32[2,1,3,4] || error(
        "Tessella did not rewind reflected tetrahedron connectivity")
    Tessella.validate(reflected_expected).ok || error(
        "Tessella reflected mesh is invalid")

    _add_gmsh_fixture("mesh-affine-singular")
    gmsh.model.mesh.affineTransform(_SINGULAR12)
    all(iszero,_gmsh_coordinates()[1,:]) || error(
        "Gmsh did not apply the measured singular transform")
    _rejects_argument(()->Tessella.affine_transform(
        source,[0.0 0.0 0.0;0.0 1.0 0.0;0.0 0.0 1.0])) || error(
        "Tessella accepted a singular affine transform")

    _add_gmsh_fixture("mesh-affine-extra-entry")
    gmsh.model.mesh.affineTransform(vcat(_POSITIVE12,[99.0]))
    _gmsh_coordinates()==gmsh_positive12 || error(
        "Gmsh no longer ignores an affine entry after its first 12")

    Tessella.API.initialize()
    api_crc=try
        Tessella.API.model.add_box(0,0,0,1,1,1;tag=1)
        generated=Tessella.API.mesh.generate(3)
        expected=Tessella.affine_transform(
            generated,positive_matrix;translation=(3.0,-2.0,4.0))
        transformed=Tessella.API.mesh.affine_transform(_POSITIVE12)
        transformed.coords==expected.coords || error(
            "session affine coordinates differ from the canonical kernel")
        transformed.tets==expected.tets || error(
            "session affine connectivity differs from the canonical kernel")
        transformed_crc=Tessella.mesh_crc(transformed)
        transformed_crc==Tessella.mesh_crc(expected) || error(
            "session affine checksum differs from the canonical kernel")
        snapshot=copy(transformed.coords)
        transformed.coords[1,1]+=1
        Tessella.API.mesh.get().coords==snapshot || error(
            "session affine result aliases the cached mesh")
        for rejected in (_SINGULAR12,vcat(_POSITIVE12,[99.0]))
            _rejects_argument(
                ()->Tessella.API.mesh.affine_transform(rejected)) || error(
                    "session affine transform accepted a bounded divergence")
            Tessella.API.mesh.get().coords==snapshot || error(
                "rejected session affine transform changed the cache")
        end
        Tessella.API.model.get_bounding_box(-1,-1)==
            (0.0,0.0,0.0,1.0,1.0,1.0) || error(
                "session affine transform changed model geometry")
        Tessella.API.mesh.clear()
        Tessella.mesh_crc(Tessella.API.mesh.generate(3))==
            Tessella.mesh_crc(generated) || error(
                "clear/regenerate did not restore the geometry-derived mesh")
        transformed_crc
    finally
        Tessella.API.finalize()
    end

    println("mesh-affine-transform differential: Gmsh ",
            gmsh.GMSH_API_VERSION,
            ", max_positive_error=",positive_error,
            " max_reflection_error=",reflection_error,
            " orientation=Tessella-rewound api_bbox=",api_crc.bbox,
            " api_sha=",api_crc.sha)
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
