#!/usr/bin/env julia
# P6: periodic dependency graphs and expression transforms vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: mesh_crc
using Tessella.Elements: mixed_crc,read_mixed_msh,write_mixed_msh

const AFFINE_Y03=(
    1.0,0.0,0.0,0.0,
    0.0,1.0,0.0,0.3,
    0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,1.0)
const AFFINE_Y06=(
    1.0,0.0,0.0,0.0,
    0.0,1.0,0.0,0.6,
    0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,1.0)
const AFFINE_YN06=(
    1.0,0.0,0.0,0.0,
    0.0,1.0,0.0,-0.6,
    0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,1.0)
const MESH_CRC=
    "dad04f30f3b17630127c3f1b4f5b5a4776ae5ff20d3c89afa6c674fac24d5338"
const CASES=(
    (name=:branch,
     path=joinpath(@__DIR__,"periodic_curve_branch.geo"),
     masters=Dict(10=>30,20=>30),offsets=Dict(10=>0.6,20=>0.3),
     point_links=[(103,101),(104,102),(105,101),(106,102)],
     projected_crc="6eb5020b186e9abdd8472312791bf375e1e9512066e7cc67b1fd2c1990a6d82f",
     msh2_crc="bdc360ec06e87180d06b5bd0bb2091c6b3fa0d4ea27d8122d39556824d97014c"),
    (name=:chain,
     path=joinpath(@__DIR__,"periodic_curve_chain.geo"),
     masters=Dict(10=>20,20=>30),offsets=Dict(10=>0.3,20=>0.3),
     point_links=[(103,101),(104,102),(105,103),(106,104)],
     projected_crc="3f98267cc70f9326ebe490c854cb59a9987c638e6aaabcba326d086bfb887ab1",
     msh2_crc="7fb21b1038e4d7f21b98ca52cb85951b8a6ed776c445f8151653fb1950878c88"),
    (name=:expressions,
     path=joinpath(@__DIR__,"periodic_curve_expressions.geo"),
     masters=Dict(10=>20,20=>30),offsets=Dict(10=>0.3,20=>0.3),
     point_links=[(103,101),(104,102),(105,103),(106,104)],
     projected_crc="3f98267cc70f9326ebe490c854cb59a9987c638e6aaabcba326d086bfb887ab1",
     msh2_crc="7fb21b1038e4d7f21b98ca52cb85951b8a6ed776c445f8151653fb1950878c88"),
)
const EXPRESSION_TRANSFORMS=(
    (name=:affine,
     path=joinpath(@__DIR__,"periodic_curve_affine_expressions.geo"),
     slave=2,master=4,tessella_pairs=5,gmsh_pairs=3,
     mesh_crc="3511d556ca0894daa79152eaf56abc6961024a72fa4f7e94f3357a7aa3cf0ff5"),
    (name=:rotate,
     path=joinpath(@__DIR__,"periodic_curve_rotate_expressions.geo"),
     slave=3,master=1,tessella_pairs=3,gmsh_pairs=3,
     mesh_crc="f6ad616e56d52d7e10a598a4079db2de9b3d5f2a777f492f5a2366946d8ea990"),
)

function expected_affine(offset)
    offset==0.3 && return AFFINE_Y03
    offset==0.6 && return AFFINE_Y06
    error("unsupported periodic graph offset $offset")
end

function expression_transform_point(name::Symbol,point)
    name==:affine && return (point[1]+1,point[2],point[3])
    name==:rotate && return (2-point[2],point[1],point[3])
    error("unsupported periodic expression transform $name")
end

function validate_graph_links(mixed,case,pair_count)
    validate(mixed).ok || error("$(case.name) periodic graph projection is invalid")
    mixed.physical_names==Dict(
        (1,41)=>"periodic traces",(2,42)=>"domain") || error(
        "$(case.name) periodic graph lost physical names")
    point_links=sort([
        (Int(link.slave_entity),Int(link.master_entity))
        for link in mixed.periodic_links if link.dim==0])
    point_links==case.point_links || error(
        "$(case.name) periodic graph changed endpoint links: $point_links")
    curve_links=Dict(
        Int(link.slave_entity)=>link
        for link in mixed.periodic_links if link.dim==1)
    Set(keys(curve_links))==Set((10,20)) || error(
        "$(case.name) periodic graph lost a curve link")
    for slave in (10,20)
        link=curve_links[slave]
        link.master_entity==case.masters[slave] || error(
            "$(case.name) Curve[$slave] has master $(link.master_entity)")
        length(link.slave_nodes)==length(link.master_nodes)==pair_count || error(
            "$(case.name) Curve[$slave] pair count is not $pair_count")
        link.affine==expected_affine(case.offsets[slave]) || error(
            "$(case.name) Curve[$slave] affine changed")
    end
    return nothing
end

projected=Dict{Symbol,MixedMesh}()
for case in CASES
    execution=execute_geo(case.path;mesh_dim=2)
    mesh=execution.mesh
    mesh===nothing && error("$(case.name) periodic graph produced no mesh")
    validate(mesh).ok || error("$(case.name) periodic graph mesh is invalid")
    mesh_crc(mesh).sha==MESH_CRC || error(
        "$(case.name) periodic graph mesh CRC changed: $(mesh_crc(mesh).sha)")
    constraints=model_periodic_constraints(execution.model)
    Int.(getproperty.(constraints,:slave_entity))==[10,20] || error(
        "$(case.name) periodic constraints lost slave order")
    Int.(getproperty.(constraints,:master_entity))==
        [case.masters[10],case.masters[20]] || error(
        "$(case.name) periodic constraints changed masters")
    constraints[1].reversed || error(
        "$(case.name) reversed Curve[10] orientation was not retained")
    for slave in (10,20)
        mapping=model_periodic_nodes(execution.model,mesh,1,slave)
        mapping.master_entity==case.masters[slave] || error(
            "$(case.name) Curve[$slave] mapped to the wrong master")
        length(mapping.slave_nodes)==length(mapping.master_nodes)==9 || error(
            "$(case.name) Curve[$slave] compact pair count is not 9")
        offset=case.offsets[slave]
        for (slave_node,master_node) in
                zip(mapping.slave_nodes,mapping.master_nodes)
            master=Tuple(mesh.coords[:,master_node])
            slave_coordinate=Tuple(mesh.coords[:,slave_node])
            slave_coordinate==(master[1],master[2]+offset,master[3]) || error(
                "$(case.name) Curve[$slave] pair is not exactly translated")
        end
    end
    mixed=model_to_mixed(execution.model,mesh,1)
    validate_graph_links(mixed,case,9)
    mixed.entity_data.entities[(2,1)].embedded_curves==Int32[10,20,30] ||
        error("$(case.name) projection lost embedded curves")
    mixed_crc(mixed).sha==case.projected_crc || error(
        "$(case.name) projected CRC changed: $(mixed_crc(mixed).sha)")
    projected[case.name]=mixed
end

expression_affines,max_tessella_transform_error=let
    affines=Dict{Symbol,NTuple{16,Float64}}()
    max_error=0.0
    for case in EXPRESSION_TRANSFORMS
        execution=execute_geo(case.path;mesh_dim=2)
        mesh=execution.mesh
        mesh===nothing && error(
            "$(case.name) expression transform produced no mesh")
        validate(mesh).ok || error(
            "$(case.name) expression transform mesh is invalid")
        mesh_crc(mesh).sha==case.mesh_crc || error(
            "$(case.name) expression transform CRC changed: $(mesh_crc(mesh).sha)")
        constraint=only(model_periodic_constraints(execution.model))
        constraint.slave_entity==case.slave &&
            constraint.master_entity==case.master || error(
            "$(case.name) expression transform changed its curve entities")
        mapping=model_periodic_nodes(execution.model,mesh,1,case.slave)
        length(mapping.slave_nodes)==length(mapping.master_nodes)==
            case.tessella_pairs || error(
            "$(case.name) Tessella expression pair count changed")
        for (slave_node,master_node) in
                zip(mapping.slave_nodes,mapping.master_nodes)
            expected=expression_transform_point(
                case.name,Tuple(mesh.coords[:,master_node]))
            got=Tuple(mesh.coords[:,slave_node])
            max_error=max(max_error,hypot((got.-expected)...))
        end
        affines[case.name]=constraint.affine
    end
    affines,max_error
end
max_tessella_transform_error<=1e-15 || error(
    "Tessella expression transform error is $max_tessella_transform_error")

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

function gmsh_coordinates()
    node_tags,node_values,_=gmsh.model.mesh.getNodes()
    coordinates=Dict{UInt64,NTuple{3,Float64}}()
    for (index,tag) in pairs(node_tags)
        coordinates[UInt64(tag)]=(
            node_values[3index-2],node_values[3index-1],node_values[3index])
    end
    return coordinates
end

gmsh.initialize(["gmsh","-v","0"])
try
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "periodic graph differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    max_gmsh_error=0.0
    for case in CASES
        gmsh.clear()
        gmsh.open(case.path)
        Set(gmsh.model.mesh.getEmbedded(2,1))==Set(
            Tuple{Int32,Int32}[(0,107),(1,10),(1,20),(1,30)]) || error(
            "Gmsh $(case.name) source lost embedded entities")
        gmsh.model.mesh.generate(2)
        _,surface_elements,_=gmsh.model.mesh.getElements(2,1)
        sum(length,surface_elements;init=0)==44 || error(
            "Gmsh $(case.name) source triangle count changed")
        coordinates=gmsh_coordinates()
        for slave in (10,20)
            master,slave_nodes,master_nodes,affine=
                gmsh.model.mesh.getPeriodicNodes(1,slave)
            master==case.masters[slave] || error(
                "Gmsh $(case.name) Curve[$slave] has master $master")
            length(slave_nodes)==length(master_nodes)==3 || error(
                "Gmsh $(case.name) Curve[$slave] pair count is not 3")
            affine==collect(expected_affine(case.offsets[slave])) || error(
                "Gmsh $(case.name) Curve[$slave] affine changed")
            offset=case.offsets[slave]
            for (slave_node,master_node) in zip(slave_nodes,master_nodes)
                slave_coordinate=coordinates[UInt64(slave_node)]
                master_coordinate=coordinates[UInt64(master_node)]
                max_gmsh_error=max(max_gmsh_error,
                    hypot(slave_coordinate[1]-master_coordinate[1],
                          slave_coordinate[2]-master_coordinate[2]-offset,
                          slave_coordinate[3]-master_coordinate[3]))
            end
        end
        for (slave,master) in case.point_links
            gmsh.model.mesh.getPeriodicNodes(0,slave)[1]==master || error(
                "Gmsh $(case.name) lost Point[$slave] -> Point[$master]")
        end

        mktempdir() do directory
            for version in (2.2,4.1),binary in (false,true)
                path=joinpath(
                    directory,"$(case.name)-$version-$binary.msh")
                write_mixed_msh(
                    path,projected[case.name];version=version,binary=binary)
                reread=read_mixed_msh(path)
                validate_graph_links(reread,case,9)
                expected_crc=version==4.1 ?
                    case.projected_crc : case.msh2_crc
                mixed_crc(reread).sha==expected_crc || error(
                    "$(case.name) MSH$version binary=$binary CRC changed")
                if version==4.1
                    reread.entity_data.entities[(2,1)].embedded_curves==
                        Int32[10,20,30] || error(
                        "$(case.name) MSH4 lost embedded curves")
                else
                    reread.entity_data===nothing || error(
                        "$(case.name) MSH2 unexpectedly retained MSH4 entities")
                end

                gmsh.clear()
                gmsh.option.setNumber("Mesh.IgnorePeriodicity",0)
                gmsh.open(path)
                for slave in (10,20)
                    master,slave_nodes,master_nodes,affine=
                        gmsh.model.mesh.getPeriodicNodes(1,slave)
                    master==case.masters[slave] &&
                        length(slave_nodes)==length(master_nodes)==9 &&
                        affine==collect(expected_affine(case.offsets[slave])) ||
                        error("Gmsh lost $(case.name) MSH$version Curve[$slave]")
                end
                for (slave,master) in case.point_links
                    gmsh.model.mesh.getPeriodicNodes(0,slave)[1]==master ||
                        error("Gmsh lost $(case.name) MSH$version Point[$slave]")
                end
                expected_embedded=version==4.1 && !binary ?
                    Tuple{Int32,Int32}[(1,10),(1,20),(1,30)] :
                    Tuple{Int32,Int32}[]
                gmsh.model.mesh.getEmbedded(2,1)==expected_embedded || error(
                    "Gmsh embedding mismatch for $(case.name) MSH$version " *
                    "binary=$binary")
            end
        end
    end
    max_gmsh_error<=1e-12 || error(
        "Gmsh periodic graph coordinate error is $max_gmsh_error")

    max_gmsh_transform_error=0.0
    max_affine_difference=0.0
    for case in EXPRESSION_TRANSFORMS
        gmsh.clear()
        gmsh.open(case.path)
        gmsh.model.mesh.generate(2)
        coordinates=gmsh_coordinates()
        master,slave_nodes,master_nodes,affine=
            gmsh.model.mesh.getPeriodicNodes(1,case.slave)
        master==case.master || error(
            "Gmsh $(case.name) expression master is $master")
        length(slave_nodes)==length(master_nodes)==case.gmsh_pairs || error(
            "Gmsh $(case.name) expression pair count changed")
        length(affine)==16 || error(
            "Gmsh $(case.name) expression affine is not 4×4")
        max_affine_difference=max(max_affine_difference,
            maximum(abs.(affine.-collect(expression_affines[case.name]))))
        for (slave_node,master_node) in zip(slave_nodes,master_nodes)
            expected=expression_transform_point(
                case.name,coordinates[UInt64(master_node)])
            got=coordinates[UInt64(slave_node)]
            max_gmsh_transform_error=max(
                max_gmsh_transform_error,hypot((got.-expected)...))
        end
    end
    max_gmsh_transform_error<=3e-12 || error(
        "Gmsh expression transform error is $max_gmsh_transform_error")
    max_affine_difference<=1e-15 || error(
        "Tessella/Gmsh expression affine difference is $max_affine_difference")

    gmsh.clear()
    gmsh.open(CASES[2].path)
    gmsh.model.mesh.setPeriodic(1,[30],[10],collect(AFFINE_YN06))
    gmsh.model.mesh.generate(2)
    all(slave->isempty(gmsh.model.mesh.getPeriodicNodes(1,slave)[2]),
        (10,20,30)) || error(
        "Gmsh 4.15.2 cyclic periodic-curve probe no longer returns empty maps")

    println("GMSH_PARITY_PERIODIC_GRAPH_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) modes=branch,chain,expressions " *
            "tessella_nodes=$(size(projected[:expressions].coords,2)) " *
            "tessella_pairs=9 gmsh_tris=44 gmsh_pairs=3 " *
            "gmsh_max_error=$max_gmsh_error transforms=translate,rotate,affine " *
            "transform_max_error=$(max(max_tessella_transform_error,max_gmsh_transform_error)) " *
            "affine_difference=$max_affine_difference cycle_maps=empty " *
            "msh2_msh4_ascii_binary=ok")
finally
    gmsh.finalize()
end
