#!/usr/bin/env julia
# P6: embedded straight-curve periodic lifecycle vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, ntris
using Tessella.Elements: mixed_crc, read_mixed_msh, write_mixed_msh

const GEO=joinpath(@__DIR__,"periodic_embedded_curve.geo")
const EXPECTED_AFFINE=(
    1.0,0.0,0.0,0.0,
    0.0,1.0,0.0,0.5,
    0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,1.0)
const MESH_CRC=
    "9794a65ea5402683d0d50612522c2f71f7c98ec2a9f6b9e6b49a61e62cd85cf2"
const PROJECTED_CRC=
    "e32e8317842c099bc4a91cdd94d02d0f816884f0e091d7194bac56e95bbfeade"
const MSH2_CRC=
    "d6da1835be0a570f81b99ccec03acd47bd46722ed69316007d3fc4fa020b2445"

execution=execute_geo(GEO;mesh_dim=2)
mesh=execution.mesh
mesh===nothing && error("Tessella embedded periodic curves produced no mesh")
validate(mesh).ok || error("Tessella embedded periodic-curve mesh is invalid")
mesh_crc(mesh).sha==MESH_CRC || error(
    "Tessella embedded periodic-curve mesh CRC changed: $(mesh_crc(mesh).sha)")
mapping=model_periodic_nodes(execution.model,mesh,1,102)
mapping.master_entity==101 || error(
    "Tessella periodic master Curve is $(mapping.master_entity), expected 101")
mapping.affine==EXPECTED_AFFINE ||
    error("Tessella changed the embedded periodic affine transform")
length(mapping.slave_nodes)==length(mapping.master_nodes)==3 || error(
    "Tessella embedded periodic pair count is not 3")
for (slave_node,master_node) in zip(mapping.slave_nodes,mapping.master_nodes)
    master=Tuple(mesh.coords[:,master_node])
    slave=Tuple(mesh.coords[:,slave_node])
    slave==(master[1],master[2]+0.5,master[3]) || error(
        "Tessella embedded periodic node pair is not exactly translated")
end

projected=model_to_mixed(execution.model,mesh,1)
validate(projected).ok ||
    error("Tessella embedded periodic projection is invalid")
projected.entity_data.entities[(2,1)].embedded_curves==Int32[101,102] ||
    error("Tessella projection lost embedded Curve[101] or Curve[102]")
projected_crc=mixed_crc(projected)
projected_crc.sha==PROJECTED_CRC || error(
    "Tessella embedded periodic projection CRC changed: $(projected_crc.sha)")
sort([(Int(link.slave_entity),Int(link.master_entity))
      for link in projected.periodic_links if link.dim==0])==
    [(103,101),(104,102)] ||
    error("Tessella projection changed periodic endpoint links")
curve_link=only(filter(link->link.dim==1,projected.periodic_links))
curve_link.slave_entity==102 && curve_link.master_entity==101 &&
    curve_link.slave_nodes==mapping.slave_nodes &&
    curve_link.master_nodes==mapping.master_nodes ||
    error("Tessella projection changed the embedded periodic curve link")

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
gmsh.initialize(["gmsh","-v","0"])
try
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "embedded periodic differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh.open(GEO)
    Set(gmsh.model.mesh.getEmbedded(2,1))==Set(
        Tuple{Int32,Int32}[(0,105),(1,101),(1,102)]) ||
        error("Gmsh source model lost embedded curve constraints")
    gmsh.model.mesh.generate(2)
    gmsh_master,gmsh_slaves,gmsh_masters,gmsh_affine=
        gmsh.model.mesh.getPeriodicNodes(1,102)
    gmsh_master==101 ||
        error("Gmsh embedded periodic master is $gmsh_master, expected 101")
    length(gmsh_slaves)==length(gmsh_masters)==2 || error(
        "Gmsh embedded periodic pair count is $(length(gmsh_slaves)), expected 2")
    gmsh_affine==collect(EXPECTED_AFFINE) ||
        error("Gmsh changed the embedded periodic affine transform")
    gmsh.model.mesh.getPeriodicNodes(0,103)[1]==101 ||
        error("Gmsh lost Point[103] -> Point[101] periodicity")
    gmsh.model.mesh.getPeriodicNodes(0,104)[1]==102 ||
        error("Gmsh lost Point[104] -> Point[102] periodicity")
    _,gmsh_surface_elements,_=gmsh.model.mesh.getElements(2,1)
    gmsh_surface_triangles=sum(length,gmsh_surface_elements;init=0)
    gmsh_surface_triangles==16 || error(
        "Gmsh embedded periodic surface has $gmsh_surface_triangles triangles")

    gmsh_node_tags,gmsh_node_coordinates,_=gmsh.model.mesh.getNodes()
    gmsh_coordinates=Dict{UInt64,NTuple{3,Float64}}()
    for (index,tag) in pairs(gmsh_node_tags)
        gmsh_coordinates[UInt64(tag)]=(
            gmsh_node_coordinates[3index-2],
            gmsh_node_coordinates[3index-1],
            gmsh_node_coordinates[3index])
    end
    max_gmsh_error=0.0
    for (slave_tag,master_tag) in zip(gmsh_slaves,gmsh_masters)
        slave=gmsh_coordinates[UInt64(slave_tag)]
        master=gmsh_coordinates[UInt64(master_tag)]
        max_gmsh_error=max(max_gmsh_error,
            hypot(slave[1]-master[1],slave[2]-master[2]-0.5,
                  slave[3]-master[3]))
    end
    max_gmsh_error<=1e-12 || error(
        "Gmsh embedded periodic coordinate error is $max_gmsh_error")

    mktempdir() do directory
        gmsh.option.setNumber("Mesh.SaveAll",1)
        gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        msh2_crcs=Set{String}()
        for binary in (false,true)
            gmsh.option.setNumber("Mesh.Binary",binary ? 1 : 0)
            native_path=joinpath(
                directory,"gmsh-periodic-embedded-$binary.msh")
            gmsh.write(native_path)
            native=read_mixed_msh(native_path)
            validate(native).ok ||
                error("Tessella rejected Gmsh's embedded periodic mesh")
            native.entity_data.entities[(2,1)].embedded_curves==
                Int32[101,102] ||
                error("Tessella lost Gmsh's embedded curves")
            native_curve=only(filter(
                link->link.dim==1,native.periodic_links))
            native_curve.slave_entity==102 &&
                native_curve.master_entity==101 &&
                length(native_curve.slave_nodes)==2 ||
                error("Tessella lost Gmsh's embedded periodic curve link")

            for version in (2.2,4.1)
                path=joinpath(
                    directory,"tessella-periodic-embedded-$version-$binary.msh")
                write_mixed_msh(path,projected;version=version,binary=binary)
                reread=read_mixed_msh(path)
                validate(reread).ok || error(
                    "Tessella embedded periodic MSH$version binary=$binary is invalid")
                reread_curve=only(filter(
                    link->link.dim==1,reread.periodic_links))
                length(reread_curve.slave_nodes)==3 || error(
                    "Tessella MSH$version lost embedded periodic node pairs")
                if version==4.1
                    mixed_crc(reread)==projected_crc || error(
                        "Tessella embedded periodic MSH4 CRC mismatch")
                    reread.entity_data.entities[(2,1)].embedded_curves==
                        Int32[101,102] ||
                        error("Tessella MSH4 lost embedded curves")
                else
                    reread.entity_data===nothing ||
                        error("Tessella MSH2 unexpectedly retained MSH4 entities")
                    push!(msh2_crcs,mixed_crc(reread).sha)
                end

                gmsh.clear()
                gmsh.option.setNumber("Mesh.IgnorePeriodicity",0)
                gmsh.open(path)
                projected_master,projected_slaves,projected_masters,
                projected_affine=gmsh.model.mesh.getPeriodicNodes(1,102)
                projected_master==101 &&
                    length(projected_slaves)==length(projected_masters)==3 &&
                    projected_affine==collect(EXPECTED_AFFINE) || error(
                        "Gmsh lost Tessella's MSH$version embedded periodic curve")
                gmsh.model.mesh.getPeriodicNodes(0,103)[1]==101 &&
                    gmsh.model.mesh.getPeriodicNodes(0,104)[1]==102 || error(
                        "Gmsh lost Tessella's periodic endpoint links")
                expected_embedded=version==4.1 && !binary ?
                    Tuple{Int32,Int32}[(1,101),(1,102)] :
                    Tuple{Int32,Int32}[]
                gmsh.model.mesh.getEmbedded(2,1)==expected_embedded || error(
                    "Gmsh embedding relation mismatch for MSH$version " *
                    "binary=$binary")
                for (dim,tags) in ((0,[1,2,3,4,101,102,103,104,105]),
                                   (1,[1,2,3,4,101,102]),(2,[1])),tag in tags
                    _,elements,_=gmsh.model.mesh.getElements(dim,tag)
                    sum(length,elements;init=0)>0 || error(
                        "Gmsh lost entity ($dim,$tag) elements for " *
                        "MSH$version binary=$binary")
                end
            end

            gmsh.clear()
            gmsh.open(GEO)
            gmsh.model.mesh.generate(2)
            gmsh.option.setNumber("Mesh.SaveAll",1)
            gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        end
        msh2_crcs==Set([MSH2_CRC]) || error(
            "Tessella embedded periodic MSH2 CRC mismatch: $msh2_crcs")
    end
    println("GMSH_PARITY_PERIODIC_EMBEDDED_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) tessella_nodes=$(nnodes(mesh)) " *
            "tessella_tris=$(ntris(mesh)) tessella_pairs=3 " *
            "gmsh_tris=$gmsh_surface_triangles gmsh_pairs=2 " *
            "gmsh_max_error=$max_gmsh_error " *
            "projected_crc=$(projected_crc.sha) " *
            "msh2_msh4_ascii_binary=ok")
finally
    gmsh.finalize()
end
