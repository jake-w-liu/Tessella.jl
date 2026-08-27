#!/usr/bin/env julia
# P6: planar periodic boundaries of an explicit volume vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, node, ntets, tet_volume
using Tessella.Elements: mixed_crc, read_mixed_msh, write_mixed_msh

const GEO=normpath(joinpath(
    @__DIR__,"..","..","test","fixtures","periodic_surface_volume.geo"))
const AFFINE_X=(
    1.0,0.0,0.0,1.0,
    0.0,1.0,0.0,0.0,
    0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,1.0)
const AFFINE_Y=(
    1.0,0.0,0.0,0.0,
    0.0,1.0,0.0,1.0,
    0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,1.0)
const RELATIONS=(
    (dim=0,slave=2,master=1,tessella_pairs=1,gmsh_pairs=1,affine=AFFINE_X),
    (dim=0,slave=3,master=2,tessella_pairs=1,gmsh_pairs=1,affine=AFFINE_Y),
    (dim=0,slave=4,master=1,tessella_pairs=1,gmsh_pairs=1,affine=AFFINE_Y),
    (dim=0,slave=6,master=5,tessella_pairs=1,gmsh_pairs=1,affine=AFFINE_X),
    (dim=0,slave=7,master=6,tessella_pairs=1,gmsh_pairs=1,affine=AFFINE_Y),
    (dim=0,slave=8,master=5,tessella_pairs=1,gmsh_pairs=1,affine=AFFINE_Y),
    (dim=1,slave=2,master=4,tessella_pairs=2,gmsh_pairs=4,affine=AFFINE_X),
    (dim=1,slave=3,master=1,tessella_pairs=2,gmsh_pairs=4,affine=AFFINE_Y),
    (dim=1,slave=6,master=8,tessella_pairs=2,gmsh_pairs=4,affine=AFFINE_X),
    (dim=1,slave=7,master=5,tessella_pairs=2,gmsh_pairs=4,affine=AFFINE_Y),
    (dim=1,slave=10,master=9,tessella_pairs=2,gmsh_pairs=4,affine=AFFINE_X),
    (dim=1,slave=11,master=10,tessella_pairs=2,gmsh_pairs=4,affine=AFFINE_Y),
    (dim=1,slave=12,master=9,tessella_pairs=2,gmsh_pairs=4,affine=AFFINE_Y),
    (dim=2,slave=4,master=6,tessella_pairs=5,gmsh_pairs=21,affine=AFFINE_X),
    (dim=2,slave=5,master=3,tessella_pairs=4,gmsh_pairs=20,affine=AFFINE_Y),
)
const EXPECTED_RELATIONS=Set(
    (relation.dim,relation.slave,relation.master) for relation in RELATIONS)
const PHYSICAL_NAMES=Dict(
    (0,61)=>"corners",(0,65)=>"face probes",(1,62)=>"edges",
    (2,63)=>"boundary",(3,64)=>"domain")

function affine_point(affine,point)
    x,y,z=point
    return (
        muladd(affine[3],z,muladd(affine[2],y,
            muladd(affine[1],x,affine[4]))),
        muladd(affine[7],z,muladd(affine[6],y,
            muladd(affine[5],x,affine[8]))),
        muladd(affine[11],z,muladd(affine[10],y,
            muladd(affine[9],x,affine[12]))),
    )
end

function relation_error(coordinates,slave_nodes,master_nodes,affine)
    maximum(zip(slave_nodes,master_nodes);init=0.0) do (slave,master)
        expected=affine_point(affine,coordinates[Int(master)])
        hypot((coordinates[Int(slave)].-expected)...)
    end
end

function mixed_coordinates(mesh)
    Dict(node=>Tuple(mesh.coords[:,node]) for node in axes(mesh.coords,2))
end

function validate_mixed_relations(mesh,pair_field::Symbol)
    validate(mesh).ok || error("periodic volume projection is invalid")
    mesh.physical_names==PHYSICAL_NAMES || error(
        "periodic volume projection changed physical names")
    links=Dict(
        (link.dim,Int(link.slave_entity),Int(link.master_entity))=>link
        for link in mesh.periodic_links)
    Set(keys(links))==EXPECTED_RELATIONS || error(
        "periodic volume projection changed its relation forest")
    coordinates=mixed_coordinates(mesh)
    max_error=0.0
    for relation in RELATIONS
        link=links[(relation.dim,relation.slave,relation.master)]
        expected_pairs=getproperty(relation,pair_field)
        length(link.slave_nodes)==length(link.master_nodes)==expected_pairs ||
            error("projected periodic entity ($(relation.dim)," *
                  "$(relation.slave)) has the wrong pair count")
        link.affine==relation.affine || error(
            "projected periodic entity ($(relation.dim)," *
            "$(relation.slave)) changed its affine transform")
        max_error=max(max_error,relation_error(
            coordinates,link.slave_nodes,link.master_nodes,link.affine))
    end
    return max_error
end

execution=execute_geo(GEO;mesh_dim=3)
mesh=execution.mesh
mesh===nothing && error("Tessella periodic volume produced no mesh")
validate(mesh).ok || error("Tessella periodic volume mesh is invalid")
nnodes(mesh)==11 && ntets(mesh)==16 || error(
    "Tessella periodic volume size changed to $(nnodes(mesh)) nodes and " *
    "$(ntets(mesh)) tetrahedra")
mesh_crc(mesh).sha==
    "2fc8151cb4a8176a9a81e02c9c3e56ca66f9f9a46baf0d14f25f751a977ad808" ||
    error("Tessella periodic volume mesh CRC changed")
tessella_volume=sum(tet_volume(
    node(mesh,mesh.tets[1,cell]),node(mesh,mesh.tets[2,cell]),
    node(mesh,mesh.tets[3,cell]),node(mesh,mesh.tets[4,cell]))
    for cell in 1:ntets(mesh))
abs(tessella_volume-1)<=1e-12 || error(
    "Tessella periodic volume is $tessella_volume, expected 1")
constraints=model_periodic_constraints(execution.model)
[(constraint.dim,Int(constraint.slave_entity),
  Int(constraint.master_entity)) for constraint in constraints]==
    [(2,4,6),(2,5,3)] || error(
        "Tessella changed the stored periodic surface relations")
for relation in filter(relation->relation.dim==2,RELATIONS)
    mapping=model_periodic_nodes(
        execution.model,mesh,relation.dim,relation.slave)
    mapping.master_entity==relation.master || error(
        "Tessella changed Surface[$(relation.slave)] master")
    length(mapping.slave_nodes)==relation.tessella_pairs || error(
        "Tessella changed Surface[$(relation.slave)] pair count")
end

projected=model_to_mixed(execution.model,mesh,3,1)
max_tessella_error=validate_mixed_relations(projected,:tessella_pairs)
max_tessella_error==0 || error(
    "Tessella periodic volume coordinate error is $max_tessella_error")
projected_crc=mixed_crc(projected)
projected_crc.sha==
    "27417f652cf93e0d6aad41c2f1b6c65af3751dfb3cb3166432d2e798f25a6493" ||
    error("Tessella periodic volume projection CRC changed")

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
    tags,values,_=gmsh.model.mesh.getNodes()
    return Dict(Int(tag)=>(values[3index-2],values[3index-1],values[3index])
                for (index,tag) in pairs(tags))
end

function gmsh_relation_keys()
    keys=Set{Tuple{Int,Int,Int}}()
    for dim in 0:2,(_,tag) in gmsh.model.getEntities(dim)
        master,slaves,_,_=gmsh.model.mesh.getPeriodicNodes(dim,tag)
        isempty(slaves) || push!(keys,(dim,Int(tag),Int(master)))
    end
    return keys
end

function validate_gmsh_relations(pair_field::Symbol)
    gmsh_relation_keys()==EXPECTED_RELATIONS || error(
        "Gmsh periodic relation forest changed")
    coordinates=gmsh_coordinates()
    max_error=0.0
    for relation in RELATIONS
        master,slaves,masters,affine=
            gmsh.model.mesh.getPeriodicNodes(relation.dim,relation.slave)
        master==relation.master || error(
            "Gmsh changed periodic entity ($(relation.dim)," *
            "$(relation.slave)) master")
        expected_pairs=getproperty(relation,pair_field)
        length(slaves)==length(masters)==expected_pairs || error(
            "Gmsh periodic entity ($(relation.dim),$(relation.slave)) " *
            "has $(length(slaves)) pairs; expected $expected_pairs")
        affine==collect(relation.affine) || error(
            "Gmsh changed periodic entity ($(relation.dim)," *
            "$(relation.slave)) affine transform")
        max_error=max(max_error,relation_error(
            coordinates,slaves,masters,relation.affine))
    end
    return max_error
end

gmsh.initialize(["gmsh","-v","0"])
try
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "periodic-volume differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh.open(GEO)
    gmsh.model.mesh.generate(3)
    native_coordinates=gmsh_coordinates()
    length(native_coordinates)==83 || error(
        "Gmsh periodic volume node count changed to " *
        "$(length(native_coordinates))")
    types,element_tags,_=gmsh.model.mesh.getElements(3,1)
    types==Int32[4] && length(only(element_tags))==188 || error(
        "Gmsh periodic volume tetrahedron count changed")
    max_gmsh_error=validate_gmsh_relations(:gmsh_pairs)
    max_gmsh_error<=1e-12 || error(
        "Gmsh periodic volume coordinate error is $max_gmsh_error")

    master_center=only(first(gmsh.model.mesh.getNodes(0,101)))
    slave_center=only(first(gmsh.model.mesh.getNodes(0,102)))
    _,surface_slaves,surface_masters,_=
        gmsh.model.mesh.getPeriodicNodes(2,4)
    any(pair->pair==(slave_center,master_center),
        zip(surface_slaves,surface_masters)) || error(
            "Gmsh periodic Surface[4] map omitted its embedded face-center pair")
    isempty(gmsh.model.mesh.getPeriodicNodes(0,102)[2]) || error(
        "Gmsh unexpectedly emitted a separate embedded Point[102] relation")

    projected_crcs=Dict(2.2=>Set{String}(),4.1=>Set{String}())
    max_roundtrip_error=0.0
    mktempdir() do directory
        for version in (2.2,4.1),binary in (false,true)
            path=joinpath(
                directory,"periodic-volume-$version-$binary.msh")
            write_mixed_msh(
                path,projected;version=version,binary=binary)
            reread=read_mixed_msh(path)
            reread_error=validate_mixed_relations(
                reread,:tessella_pairs)
            reread_error==0 || error(
                "Tessella MSH$version periodic coordinate error is " *
                "$reread_error")
            push!(projected_crcs[version],mixed_crc(reread).sha)
            if version==4.1
                reread.entity_data.entities[(3,1)].boundaries==
                    Int32.(1:6) || error(
                        "Tessella MSH4 lost periodic volume boundaries")
            else
                reread.entity_data===nothing || error(
                    "Tessella MSH2 unexpectedly retained MSH4 entities")
            end

            gmsh.clear()
            gmsh.option.setNumber("Mesh.IgnorePeriodicity",0)
            gmsh.open(path)
            max_roundtrip_error=max(
                max_roundtrip_error,
                validate_gmsh_relations(:tessella_pairs))
            reopened_boundary=gmsh.model.getBoundary(
                [(3,1)],false,true,false)
            expected_boundary=version==4.1 ?
                [(Int32(2),Int32(surface)) for surface in 1:6] :
                Tuple{Int32,Int32}[]
            reopened_boundary==expected_boundary || error(
                "Gmsh MSH$version periodic volume boundary changed")
        end
    end
    projected_crcs[4.1]==Set([projected_crc.sha]) || error(
        "periodic volume MSH4 CRC depends on file mode")
    projected_crcs[2.2]==Set([
        "9cc65eb95bbcca5508016ff7cc1340a6d1a7311d0482c2759444c16ce4120502"]) ||
        error("periodic volume MSH2 CRC changed or depends on file mode")
    max_roundtrip_error<=1e-12 || error(
        "Gmsh periodic volume round-trip error is $max_roundtrip_error")

    println("GMSH_PARITY_PERIODIC_SURFACE_VOLUME_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) tessella_nodes=$(nnodes(mesh)) " *
            "tessella_tets=$(ntets(mesh)) gmsh_nodes=$(length(native_coordinates)) " *
            "gmsh_tets=$(length(only(element_tags))) relations=$(length(RELATIONS)) " *
            "embedded_surface_pair=ok max_error=" *
            "$(max(max_gmsh_error,max_roundtrip_error)) " *
            "msh2_msh4_ascii_binary=ok")
finally
    gmsh.finalize()
end
