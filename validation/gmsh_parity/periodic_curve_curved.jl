#!/usr/bin/env julia
# P4/P6: Tessella curved periodic curve pairs vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: validate
using Tessella.Model: model_periodic_nodes

function find_gmsh_api()
    explicit=get(ENV,"GMSH_JULIA_API","")
    !isempty(explicit) && isfile(explicit) && return explicit
    executable=Sys.which("gmsh")
    executable===nothing && error("gmsh is not on PATH")
    prefix=dirname(dirname(realpath(executable)))
    for path in (joinpath(prefix,"lib","gmsh.jl"),
                 "/opt/homebrew/opt/gmsh/lib/gmsh.jl")
        isfile(path) && return path
    end
    error("could not locate gmsh.jl")
end

include(find_gmsh_api())
const GEO=joinpath(@__DIR__,"periodic_curve_curved.geo")
const EXPECTED_AFFINE=[1.0,0.0,0.0,0.0,
                       0.0,1.0,0.0,3.0,
                       0.0,0.0,1.0,0.0,
                       0.0,0.0,0.0,1.0]
const PAIRS=(
    (slave=2,master=1,count=9,name="circle"),
    (slave=12,master=11,count=7,name="spline"))

gmsh.initialize(["gmsh","-v","0"])
try
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "curved periodic differential requires Gmsh 4.15.2, got $(gmsh.GMSH_API_VERSION)")
    gmsh.open(GEO)
    gmsh.model.mesh.generate(2)

    node_tags,node_coordinates,_=gmsh.model.mesh.getNodes()
    coordinates=Dict{Int,NTuple{3,Float64}}()
    for (i,tag) in enumerate(node_tags)
        coordinates[Int(tag)]=(node_coordinates[3i-2],node_coordinates[3i-1],
                               node_coordinates[3i])
    end

    gmsh_pair_coordinates=Dict{Int,Vector{NTuple{2,NTuple{3,Float64}}}}()
    for case in PAIRS
        master_entity,slave_tags,master_tags,affine=
            gmsh.model.mesh.getPeriodicNodes(1,case.slave)
        master_entity==case.master || error(
            "$(case.name): Gmsh periodic master curve is $master_entity, " *
            "expected $(case.master)")
        length(slave_tags)==length(master_tags)==case.count || error(
            "$(case.name): Gmsh periodic pair count is $(length(slave_tags)), " *
            "expected $(case.count)")
        length(affine)==16 || error(
            "$(case.name): Gmsh periodic affine is not 4×4")
        maximum(abs.(affine.-EXPECTED_AFFINE))<=1e-14 || error(
            "$(case.name): Gmsh periodic affine is $affine")
        pairs=NTuple{2,NTuple{3,Float64}}[]
        max_gmsh_error=0.0
        for (master_tag,slave_tag) in zip(master_tags,slave_tags)
            master=coordinates[Int(master_tag)]
            slave=coordinates[Int(slave_tag)]
            # Gmsh re-evaluates the slave curve at the inherited parameter
            # rather than transforming the master coordinate, so its pairs are
            # only affine-identical up to its own evaluator precision.
            max_gmsh_error=max(max_gmsh_error,
                hypot(slave[1]-master[1],slave[2]-master[2]-3.0,
                      slave[3]-master[3]))
            push!(pairs,(master,slave))
        end
        max_gmsh_error<=1e-7 || error(
            "$(case.name): Gmsh periodic correspondence error is $max_gmsh_error")
        gmsh_pair_coordinates[case.slave]=pairs
    end

    execution=execute_geo(GEO)
    mesh=execution.mesh
    mesh===nothing && error("Tessella curved periodic `.geo` did not mesh")
    validate(mesh).ok || error(
        "Tessella curved periodic `.geo` mesh is invalid")
    for case in PAIRS
        mapping=model_periodic_nodes(execution.model,mesh,1,case.slave)
        mapping.master_entity==case.master || error(
            "$(case.name): Tessella periodic master curve is " *
            "$(mapping.master_entity), expected $(case.master)")
        mapping.affine==Tuple(EXPECTED_AFFINE) || error(
            "$(case.name): Tessella changed the periodic affine transform")
        length(mapping.slave_nodes)==length(mapping.master_nodes)==
            case.count || error(
            "$(case.name): Tessella periodic pair count is " *
            "$(length(mapping.slave_nodes)), expected $(case.count)")
        gmsh_pairs=gmsh_pair_coordinates[case.slave]
        max_difference=0.0
        for (master_node,slave_node) in
                zip(mapping.master_nodes,mapping.slave_nodes)
            master=Tuple(mesh.coords[:,master_node])
            slave=Tuple(mesh.coords[:,slave_node])
            slave==(master[1],master[2]+3.0,master[3]) || error(
                "$(case.name): Tessella periodic pair was not transformed exactly")
            best=argmin(gmsh_pairs) do (gmsh_master,_)
                hypot((master.-gmsh_master)...)
            end
            max_difference=max(max_difference,
                hypot((master.-best[1])...),hypot((slave.-best[2])...))
        end
        max_difference<=1e-7 || error(
            "$(case.name): Tessella/Gmsh periodic pair difference is " *
            "$max_difference")
    end
    println("CURVED_PERIODIC_DIFFERENTIAL_OK gmsh=4.15.2 " *
            "cases=$(length(PAIRS)) pairs=$(sum(case.count for case in PAIRS))")
finally
    gmsh.finalize()
end
