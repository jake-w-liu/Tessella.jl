#!/usr/bin/env julia
# P6: holed Surface-In-Volume lifecycle vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: nnodes, ntets, node, tet_volume, triangle_area
using Tessella.Elements: mixed_crc, read_mixed_msh, write_mixed_msh

const GEO=joinpath(@__DIR__,"embed_sheet_hole.geo")
const SURFACE_BOUNDARIES=Int32[101,102,103,104,-108,-107,-106,-105]
const PROJECTED_CRC=
    "0e92af2702054065d564461a691f4035ab5bace358bd5f37821c9dcb5f54730d"
const MSH2_CRC=
    "250f6627ef3712e881a363b0e6d8999a6e77ea263503a0d813c6a0d42169a400"

function triangle_stats(coordinate,connectivity)
    area=0.0
    hole_centroid_hits=0
    for cell in axes(connectivity,2)
        points=ntuple(slot->coordinate(connectivity[slot,cell]),3)
        area+=triangle_area(points...)
        centroid=ntuple(axis->sum(point[axis] for point in points)/3,3)
        0.4<centroid[1]<0.6 && 0.4<centroid[2]<0.6 &&
            (hole_centroid_hits+=1)
    end
    return area,hole_centroid_hits
end

result=execute_geo(GEO;mesh_dim=3)
mesh=result.mesh
mesh===nothing && error("Tessella holed sheet produced no mesh")
validate(mesh).ok || error("Tessella holed-sheet volume mesh is invalid")
ntets(mesh)>0 || error("Tessella holed sheet has no tetrahedra")
volume=sum(tet_volume(
    node(mesh,mesh.tets[1,cell]),node(mesh,mesh.tets[2,cell]),
    node(mesh,mesh.tets[3,cell]),node(mesh,mesh.tets[4,cell]))
    for cell in 1:ntets(mesh))
abs(volume-1)<=1e-12 || error("Tessella holed-sheet volume $volume != 1")

projected=model_to_mixed(result.model,mesh,3,1)
validate(projected).ok || error("Tessella holed-sheet projection is invalid")
surface_index=only(findall(block->block.msh==2,projected.blocks))
surface_block=projected.blocks[surface_index]
all(==(Int32(101)),projected.entity_data.block_entities[surface_index]) ||
    error("Tessella holed-sheet triangles lost Surface[101] ownership")
coordinate(node_tag)=(projected.coords[1,node_tag],
                      projected.coords[2,node_tag],
                      projected.coords[3,node_tag])
projected_area,projected_hole_hits=
    triangle_stats(coordinate,surface_block.nodes)
abs(projected_area-0.45)<=1e-12 ||
    error("Tessella holed-sheet projected area $projected_area != 0.45")
projected_hole_hits==0 || error(
    "Tessella classified $projected_hole_hits triangle centroids inside the hole")
surface_entity=projected.entity_data.entities[(2,101)]
surface_entity.boundaries==SURFACE_BOUNDARIES || error(
    "Tessella holed-sheet boundaries changed: $(surface_entity.boundaries)")
surface_entity.embedded_curves==Int32[109] ||
    error("Tessella holed sheet lost nested Curve[109]")
projected_crc=mixed_crc(projected)
projected_crc.sha==PROJECTED_CRC || error(
    "Tessella holed-sheet projection CRC changed: $(projected_crc.sha)")

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
        "holed-sheet differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh.open(GEO)
    gmsh.model.getBoundary([(2,101)],false,true,false)==
        [(Int32(1),tag) for tag in SURFACE_BOUNDARIES] ||
        error("Gmsh source Surface[101] boundary mismatch")
    Set(gmsh.model.mesh.getEmbedded(2,101))==Set(
        Tuple{Int32,Int32}[(0,111),(1,109)]) ||
        error("Gmsh source model lost nested sheet constraints")
    gmsh.model.mesh.getEmbedded(3,1)==Tuple{Int32,Int32}[(2,101)] ||
        error("Gmsh source model lost Surface[101] In Volume[1]")
    gmsh.option.setNumber("Mesh.MeshSizeMax",0.5)
    gmsh.model.mesh.generate(3)
    node_tags,node_coordinates,_=gmsh.model.mesh.getNodes()
    gmsh_coordinates=Dict{UInt64,NTuple{3,Float64}}()
    for (index,tag) in pairs(node_tags)
        offset=3index-2
        gmsh_coordinates[UInt64(tag)]=(
            node_coordinates[offset],node_coordinates[offset+1],
            node_coordinates[offset+2])
    end
    types,surface_tags,surface_nodes=gmsh.model.mesh.getElements(2,101)
    types==Int32[2] || error("Gmsh holed sheet emitted element types $types")
    connectivity=reshape(Int32.(only(surface_nodes)),3,:)
    gmsh_coordinate(tag)=gmsh_coordinates[UInt64(tag)]
    gmsh_area,gmsh_hole_hits=triangle_stats(gmsh_coordinate,connectivity)
    abs(gmsh_area-0.45)<=1e-12 ||
        error("Gmsh holed-sheet area $gmsh_area != 0.45")
    gmsh_hole_hits==0 ||
        error("Gmsh emitted $gmsh_hole_hits triangle centroids inside the hole")
    volume_types,volume_elements,_=gmsh.model.mesh.getElements(3,1)
    volume_types==Int32[4] ||
        error("Gmsh holed sheet emitted volume element types $volume_types")
    gmsh_tets=sum(length,volume_elements;init=0)
    gmsh_tets>0 || error("Gmsh holed sheet has no volume elements")

    mktempdir() do directory
        gmsh.option.setNumber("Mesh.SaveAll",1)
        gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        msh2_crcs=Set{String}()
        for binary in (false,true)
            gmsh.option.setNumber("Mesh.Binary",binary ? 1 : 0)
            native_path=joinpath(directory,"gmsh-holed-sheet-$binary.msh")
            gmsh.write(native_path)
            native=read_mixed_msh(native_path)
            validate(native).ok ||
                error("Tessella did not validate Gmsh's holed-sheet mesh")
            native.entity_data.entities[(2,101)].boundaries==
                SURFACE_BOUNDARIES ||
                error("Tessella did not read Gmsh's holed-sheet boundaries")
            native.entity_data.entities[(2,101)].embedded_curves==Int32[109] ||
                error("Tessella did not read Gmsh's nested Curve[109]")

            for version in (2.2,4.1)
                path=joinpath(
                    directory,"tessella-holed-sheet-$version-$binary.msh")
                write_mixed_msh(path,projected;version=version,binary=binary)
                reread=read_mixed_msh(path)
                validate(reread).ok || error(
                    "Tessella holed-sheet MSH$version binary=$binary is invalid")
                if version==4.1
                    mixed_crc(reread)==projected_crc || error(
                        "Tessella holed-sheet MSH4 binary=$binary CRC mismatch")
                    reread.entity_data.entities[(2,101)].boundaries==
                        SURFACE_BOUNDARIES ||
                        error("Tessella MSH4 lost holed-sheet boundaries")
                else
                    reread.entity_data===nothing ||
                        error("Tessella MSH2 unexpectedly retained MSH4 entities")
                    push!(msh2_crcs,mixed_crc(reread).sha)
                end

                gmsh.clear()
                gmsh.open(path)
                isempty(gmsh.model.mesh.getEmbedded(3,1)) || error(
                    "Gmsh reconstructed a nonexistent volume-embedding record")
                expected_relation=version==4.1 && !binary ?
                    Tuple{Int32,Int32}[(1,109)] : Tuple{Int32,Int32}[]
                gmsh.model.mesh.getEmbedded(2,101)==expected_relation || error(
                    "Gmsh nested relation mismatch for MSH$version " *
                    "binary=$binary")
                expected_boundary=version==4.1 ?
                    [(Int32(1),tag) for tag in SURFACE_BOUNDARIES] :
                    Tuple{Int32,Int32}[]
                gmsh.model.getBoundary([(2,101)],false,true,false)==
                    expected_boundary || error(
                    "Gmsh boundary mismatch for MSH$version binary=$binary")
                for (dim,tags) in ((0,101:111),(1,101:109),
                                   (2,101:101),(3,1:1)),tag in tags
                    _,elements,_=gmsh.model.mesh.getElements(dim,tag)
                    sum(length,elements;init=0)>0 || error(
                        "Gmsh lost entity ($dim,$tag) elements for " *
                        "MSH$version binary=$binary")
                end
            end

            gmsh.clear()
            gmsh.open(GEO)
            gmsh.option.setNumber("Mesh.MeshSizeMax",0.5)
            gmsh.model.mesh.generate(3)
            gmsh.option.setNumber("Mesh.SaveAll",1)
            gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        end
        msh2_crcs==Set([MSH2_CRC]) || error(
            "Tessella holed-sheet MSH2 CRC mismatch: $msh2_crcs")
    end
    println("GMSH_PARITY_HOLED_SHEET_OK gmsh=$(gmsh.GMSH_API_VERSION) " *
            "tessella_volume=$volume tessella_tets=$(ntets(mesh)) " *
            "tessella_nodes=$(nnodes(mesh)) gmsh_tets=$gmsh_tets " *
            "gmsh_surface_tris=$(length(only(surface_tags))) " *
            "sheet_area=$projected_area projected_crc=$(projected_crc.sha) " *
            "msh2_msh4_ascii_binary=ok hole_centroid_hits=0")
finally
    gmsh.finalize()
end
