#!/usr/bin/env julia
# P6: explicit native model-topology queries vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: mesh_crc, nnodes, ntets, validate

const EXPECTED_MESH_CRC=
    "71ab10cf31fa64d469e1bc3985bd8c50bb240d1cdefaebbc17101bce22e7008b"

function add_tessella_topology!()
    for (tag,x,y,z) in ((10,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (7,0.0,1.0,0.0),(5,0.0,0.0,1.0))
        Tessella.API.model.add_point(x,y,z;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((8,10,2),(3,2,7),(11,7,10),(6,10,5),(12,2,5),(4,7,5))
        Tessella.API.model.add_line(first_point,last_point;tag=tag)
    end
    for (tag,curves) in ((21,[8,3,11]),(22,[8,12,-6]),
                         (23,[3,4,-12]),(24,[11,6,-4]))
        Tessella.API.model.add_curve_loop(curves;tag=tag)
        Tessella.API.model.add_plane_surface([tag];tag=tag)
    end
    Tessella.API.model.add_surface_loop([21,-22,23,-24];tag=30)
    Tessella.API.model.add_volume([30];tag=40)
    return nothing
end

function add_tessella_tetra_shell!(offset,scale)
    coordinates=((0.0,0.0,0.0),(scale,0.0,0.0),
                 (0.0,scale,0.0),(0.0,0.0,scale))
    for (index,point) in pairs(coordinates)
        Tessella.API.model.add_point(point...;tag=offset+index)
    end
    for (index,(first_point,last_point)) in
            pairs(((1,2),(2,3),(3,1),(1,4),(2,4),(3,4)))
        Tessella.API.model.add_line(
            offset+first_point,offset+last_point;tag=offset+index)
    end
    for (index,curves) in pairs(((1,2,3),(1,5,-4),(2,6,-5),(3,4,-6)))
        signed_curves=Int[sign(curve)*(offset+abs(curve)) for curve in curves]
        Tessella.API.model.add_curve_loop(signed_curves;tag=offset+index)
        Tessella.API.model.add_plane_surface([offset+index];tag=offset+index)
    end
    Tessella.API.model.add_surface_loop(
        [offset+1,-(offset+2),offset+3,-(offset+4)];tag=offset+1)
    return offset+1
end

function add_tessella_heldout_topology!()
    for (tag,x,y) in
            ((201,0.0,0.0),(202,2.0,0.0),(203,2.0,2.0),(204,0.0,2.0),
             (205,0.5,0.5),(206,1.5,0.5),(207,1.5,1.5),(208,0.5,1.5))
        Tessella.API.model.add_point(x,y,0.0;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((201,201,202),(202,202,203),(203,203,204),(204,204,201),
             (205,205,206),(206,206,207),(207,207,208),(208,208,205))
        Tessella.API.model.add_line(first_point,last_point;tag=tag)
    end
    Tessella.API.model.add_curve_loop([201,202,203,204];tag=201)
    Tessella.API.model.add_curve_loop([205,206,207,208];tag=202)
    Tessella.API.model.add_plane_surface([201,202];tag=201)

    outer=add_tessella_tetra_shell!(300,2.0)
    inner=add_tessella_tetra_shell!(400,1.0)
    Tessella.API.model.add_volume([outer,inner];tag=300)

    for (tag,x,y) in ((501,0.2,0.2),(502,0.4,0.2),(503,0.2,0.4))
        Tessella.API.model.add_point(x,y,0.0;tag=tag)
    end
    Tessella.API.model.add_line(501,502;tag=501)
    Tessella.API.model.embed(0,[503],2,21)
    Tessella.API.model.embed(1,[501],2,21)
    return nothing
end

function tessella_heldout_results()
    return (
        hole_oriented=Tessella.API.model.get_boundary(
            [(2,201)],false,true,false),
        cavity_oriented=Tessella.API.model.get_boundary(
            [(3,300)],false,true,false),
        negative_volume=Tessella.API.model.get_boundary(
            [(3,-300)],false,true,false),
        hole_adjacencies=Tessella.API.model.get_adjacencies(2,201),
        cavity_adjacencies=Tessella.API.model.get_adjacencies(3,300),
        embedded_endpoint_adjacencies=
            Tessella.API.model.get_adjacencies(0,501),
        embedded_point_adjacencies=
            Tessella.API.model.get_adjacencies(0,503),
        embedded_curve_adjacencies=
            Tessella.API.model.get_adjacencies(1,501),
        embedded_surface_adjacencies=
            Tessella.API.model.get_adjacencies(2,21),
    )
end

function tessella_topology_results()
    Tessella.API.finalize()
    return try
        Tessella.API.initialize()
        add_tessella_topology!()
        results=(
            entities=Tessella.API.model.get_entities(),
            surfaces=Tessella.API.model.get_entities(2),
            dimension=Tessella.API.model.get_dimension(),
            point_direct=Tessella.API.model.get_boundary(
                [(0,10)],false,false,false),
            point_recursive=Tessella.API.model.get_boundary(
                [(0,10)],false,false,true),
            line_direct=Tessella.API.model.get_boundary(
                [(1,8)],false,false,false),
            negative_line=Tessella.API.model.get_boundary(
                [(1,-8)],false,false,false),
            combined_line=Tessella.API.model.get_boundary(
                [(1,8)],true,false,false),
            surface_direct=Tessella.API.model.get_boundary(
                [(2,22)],false,false,false),
            surface_oriented=Tessella.API.model.get_boundary(
                [(2,22)],false,true,false),
            surface_combined=Tessella.API.model.get_boundary(
                [(2,22)],true,true,false),
            surface_recursive=Tessella.API.model.get_boundary(
                [(2,22)],false,false,true),
            volume_oriented=Tessella.API.model.get_boundary(
                [(3,40)],false,true,false),
            two_direct=Tessella.API.model.get_boundary(
                [(2,21),(2,22)],false,true,false),
            two_combined=Tessella.API.model.get_boundary(
                [(2,21),(2,22)],true,false,false),
            two_combined_oriented=Tessella.API.model.get_boundary(
                [(2,21),(2,22)],true,true,false),
            two_recursive=Tessella.API.model.get_boundary(
                [(2,21),(2,22)],true,false,true),
            duplicate_combined=Tessella.API.model.get_boundary(
                [(2,21),(2,21)],true,true,false),
            triple_combined=Tessella.API.model.get_boundary(
                [(2,21),(2,21),(2,21)],true,true,false),
            mixed_recursive=Tessella.API.model.get_boundary(
                [(1,8),(2,21)],true,false,true),
            point_adjacencies=Tessella.API.model.get_adjacencies(0,10),
            curve_adjacencies=Tessella.API.model.get_adjacencies(1,8),
            surface_adjacencies=Tessella.API.model.get_adjacencies(2,21),
            volume_adjacencies=Tessella.API.model.get_adjacencies(3,40),
        )
        mesh=Tessella.API.mesh.generate(3)
        validate(mesh).ok || error("Tessella topology-query mesh is invalid")
        nnodes(mesh)==5 && ntets(mesh)==4 || error(
            "Tessella topology-query mesh size changed")
        crc=mesh_crc(mesh).sha
        crc==EXPECTED_MESH_CRC || error(
            "Tessella topology-query mesh CRC changed: $crc")
        cached=Tessella.API.LAST_MESH[]
        Tessella.API.model.get_boundary(
            [(2,21),(2,22)],true,true,false)==results.two_combined_oriented ||
            error("Tessella topology query changed after meshing")
        Tessella.API.LAST_MESH[]===cached || error(
            "Tessella topology query invalidated the mesh cache")
        mesh_crc(Tessella.API.mesh.get()).sha==crc || error(
            "Tessella topology query changed the cached mesh")
        add_tessella_heldout_topology!()
        heldout=tessella_heldout_results()
        (;results,heldout,crc)
    finally
        Tessella.API.finalize()
    end
end

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

gmsh_dim_tags(values)=Tuple{Int,Int}[
    (Int(dimension),Int(tag)) for (dimension,tag) in values]
gmsh_tags(values)=Int.(values)

function gmsh_adjacencies(dimension,tag)
    upward,downward=gmsh.model.getAdjacencies(dimension,tag)
    return gmsh_tags(upward),gmsh_tags(downward)
end

function gmsh_boundary(dim_tags,combined,oriented,recursive)
    return gmsh_dim_tags(gmsh.model.getBoundary(
        dim_tags,combined,oriented,recursive))
end

function gmsh_topology_results()
    gmsh.clear()
    gmsh.model.add("model_topology_queries")
    for (tag,x,y,z) in ((10,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (7,0.0,1.0,0.0),(5,0.0,0.0,1.0))
        gmsh.model.geo.addPoint(x,y,z,1.0,tag)
    end
    for (tag,first_point,last_point) in
            ((8,10,2),(3,2,7),(11,7,10),(6,10,5),(12,2,5),(4,7,5))
        gmsh.model.geo.addLine(first_point,last_point,tag)
    end
    for (tag,curves) in ((21,[8,3,11]),(22,[8,12,-6]),
                         (23,[3,4,-12]),(24,[11,6,-4]))
        gmsh.model.geo.addCurveLoop(curves,tag)
        gmsh.model.geo.addPlaneSurface([tag],tag)
    end
    gmsh.model.geo.addSurfaceLoop([21,-22,23,-24],30)
    gmsh.model.geo.addVolume([30],40)
    gmsh.model.geo.synchronize()
    return (
        entities=gmsh_dim_tags(gmsh.model.getEntities()),
        surfaces=gmsh_dim_tags(gmsh.model.getEntities(2)),
        dimension=Int(gmsh.model.getDimension()),
        point_direct=gmsh_boundary([(0,10)],false,false,false),
        point_recursive=gmsh_boundary([(0,10)],false,false,true),
        line_direct=gmsh_boundary([(1,8)],false,false,false),
        negative_line=gmsh_boundary([(1,-8)],false,false,false),
        combined_line=gmsh_boundary([(1,8)],true,false,false),
        surface_direct=gmsh_boundary([(2,22)],false,false,false),
        surface_oriented=gmsh_boundary([(2,22)],false,true,false),
        surface_combined=gmsh_boundary([(2,22)],true,true,false),
        surface_recursive=gmsh_boundary([(2,22)],false,false,true),
        volume_oriented=gmsh_boundary([(3,40)],false,true,false),
        two_direct=gmsh_boundary([(2,21),(2,22)],false,true,false),
        two_combined=gmsh_boundary([(2,21),(2,22)],true,false,false),
        two_combined_oriented=gmsh_boundary(
            [(2,21),(2,22)],true,true,false),
        two_recursive=gmsh_boundary([(2,21),(2,22)],true,false,true),
        duplicate_combined=gmsh_boundary(
            [(2,21),(2,21)],true,true,false),
        triple_combined=gmsh_boundary(
            [(2,21),(2,21),(2,21)],true,true,false),
        mixed_recursive=gmsh_boundary([(1,8),(2,21)],true,false,true),
        point_adjacencies=gmsh_adjacencies(0,10),
        curve_adjacencies=gmsh_adjacencies(1,8),
        surface_adjacencies=gmsh_adjacencies(2,21),
        volume_adjacencies=gmsh_adjacencies(3,40),
    )
end

function add_gmsh_tetra_shell!(offset,scale)
    coordinates=((0.0,0.0,0.0),(scale,0.0,0.0),
                 (0.0,scale,0.0),(0.0,0.0,scale))
    for (index,point) in pairs(coordinates)
        gmsh.model.geo.addPoint(point...,1.0,offset+index)
    end
    for (index,(first_point,last_point)) in
            pairs(((1,2),(2,3),(3,1),(1,4),(2,4),(3,4)))
        gmsh.model.geo.addLine(
            offset+first_point,offset+last_point,offset+index)
    end
    for (index,curves) in pairs(((1,2,3),(1,5,-4),(2,6,-5),(3,4,-6)))
        signed_curves=Int[sign(curve)*(offset+abs(curve)) for curve in curves]
        gmsh.model.geo.addCurveLoop(signed_curves,offset+index)
        gmsh.model.geo.addPlaneSurface([offset+index],offset+index)
    end
    gmsh.model.geo.addSurfaceLoop(
        [offset+1,-(offset+2),offset+3,-(offset+4)],offset+1)
    return offset+1
end

function add_gmsh_heldout_topology!()
    for (tag,x,y) in
            ((201,0.0,0.0),(202,2.0,0.0),(203,2.0,2.0),(204,0.0,2.0),
             (205,0.5,0.5),(206,1.5,0.5),(207,1.5,1.5),(208,0.5,1.5))
        gmsh.model.geo.addPoint(x,y,0.0,1.0,tag)
    end
    for (tag,first_point,last_point) in
            ((201,201,202),(202,202,203),(203,203,204),(204,204,201),
             (205,205,206),(206,206,207),(207,207,208),(208,208,205))
        gmsh.model.geo.addLine(first_point,last_point,tag)
    end
    gmsh.model.geo.addCurveLoop([201,202,203,204],201)
    gmsh.model.geo.addCurveLoop([205,206,207,208],202)
    gmsh.model.geo.addPlaneSurface([201,202],201)

    outer=add_gmsh_tetra_shell!(300,2.0)
    inner=add_gmsh_tetra_shell!(400,1.0)
    gmsh.model.geo.addVolume([outer,inner],300)

    for (tag,x,y) in ((501,0.2,0.2),(502,0.4,0.2),(503,0.2,0.4))
        gmsh.model.geo.addPoint(x,y,0.0,1.0,tag)
    end
    gmsh.model.geo.addLine(501,502,501)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0,[503],2,21)
    gmsh.model.mesh.embed(1,[501],2,21)
    return nothing
end

function gmsh_heldout_results()
    return (
        hole_oriented=gmsh_boundary([(2,201)],false,true,false),
        cavity_oriented=gmsh_boundary([(3,300)],false,true,false),
        negative_volume=gmsh_boundary([(3,-300)],false,true,false),
        hole_adjacencies=gmsh_adjacencies(2,201),
        cavity_adjacencies=gmsh_adjacencies(3,300),
        embedded_endpoint_adjacencies=gmsh_adjacencies(0,501),
        embedded_point_adjacencies=gmsh_adjacencies(0,503),
        embedded_curve_adjacencies=gmsh_adjacencies(1,501),
        embedded_surface_adjacencies=gmsh_adjacencies(2,21),
    )
end

tessella=tessella_topology_results()
gmsh.initialize(["gmsh","-v","0"])
try
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "model-topology differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    expected=gmsh_topology_results()
    add_gmsh_heldout_topology!()
    heldout=gmsh_heldout_results()
    tessella.results==expected || error(
        "model-topology queries differ: Tessella=$(tessella.results) " *
        "Gmsh=$expected")
    tessella.heldout==heldout || error(
        "held-out model-topology queries differ: " *
        "Tessella=$(tessella.heldout) Gmsh=$heldout")
    println("GMSH_PARITY_MODEL_TOPOLOGY_OK " *
            "gmsh=$(gmsh.GMSH_API_VERSION) entities=$(length(expected.entities)) " *
            "boundary_cases=20 adjacency_cases=10 mesh_crc=$(tessella.crc)")
finally
    gmsh.finalize()
end
