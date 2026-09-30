#!/usr/bin/env julia
# P4: curved-boundary planar surface meshing vs Gmsh 4.15.2 — disk, annulus,
# circular segment, and a periodic surface pair with a Circle-arc edge.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: nnodes, ntris
using Tessella.Model: mesh_model_surface

const GEO=joinpath(@__DIR__,"curved_surface.geo")
const PERIODIC_GEO=joinpath(@__DIR__,"curved_periodic_surface.geo")

execution=execute_geo(GEO;mesh_dim=0)
m=execution.model
meshes=Dict{Int,typeof(mesh_model_surface(m,1))}()
for tag in (1,2,3)
    mesh=mesh_model_surface(m,tag)
    validate(mesh).ok || error("Tessella curved surface $tag mesh is invalid")
    meshes[tag]=mesh
end

# Boundary nodes must follow the evaluated arcs, not endpoint chords: the disk
# rim carries interior nodes at radius 1, the annulus hole rim at radius 0.5,
# and the segment's arc bulges above its chord.
disk=meshes[1]
rim=count(i->abs(hypot(disk.coords[1,i],disk.coords[2,i])-1.0)<1e-9,
          1:nnodes(disk))
rim>=4 || error("Tessella disk boundary carries only $rim rim nodes")
all(i->hypot(disk.coords[1,i],disk.coords[2,i])<=1.0+1e-9,1:nnodes(disk)) ||
    error("Tessella disk mesh escaped the circular boundary")
annulus=meshes[2]
hole_rim=count(i->abs(hypot(annulus.coords[1,i]-4.0,
                            annulus.coords[2,i])-0.5)<1e-9,1:nnodes(annulus))
hole_rim>=3 || error("Tessella annulus hole carries only $hole_rim rim nodes")
all(i->(r=hypot(annulus.coords[1,i]-4.0,annulus.coords[2,i]);
        r>=0.5-1e-9 && r<=1.0+1e-9),1:nnodes(annulus)) ||
    error("Tessella annulus mesh escaped its circular boundaries")
segment=meshes[3]
any(i->segment.coords[2,i]>0.3,1:nnodes(segment)) ||
    error("Tessella segment lost the arc bulge above its chord")

# Periodic curved pair: the slave is a bitwise affine copy of the master and
# the derived curve masters cover the Circle arcs.
pexecution=execute_geo(PERIODIC_GEO;mesh_dim=0)
pm=pexecution.model
curve_masters=Dict{Int,Int}(
    Int(constraint.slave_entity)=>Int(constraint.master_entity)
    for ((dim,_),constraint) in pm.periodic if dim==1)
curve_masters==Dict(5=>1,6=>2,7=>3,8=>4) ||
    error("Tessella derived curve masters changed: $curve_masters")
master=mesh_model_surface(pm,1)
slave=mesh_model_surface(pm,2)
validate(master).ok && validate(slave).ok ||
    error("Tessella periodic curved surface mesh is invalid")
nnodes(master)==nnodes(slave) && ntris(master)==ntris(slave) ||
    error("Tessella periodic curved surfaces diverged")
affine=pm.periodic[(2,2)].affine
master_set=Set(NTuple{3,Float64}[
    (affine[4]+muladd(affine[3],master.coords[3,i],
                      muladd(affine[2],master.coords[2,i],
                             affine[1]*master.coords[1,i])),
     affine[8]+muladd(affine[7],master.coords[3,i],
                      muladd(affine[6],master.coords[2,i],
                             affine[5]*master.coords[1,i])),
     affine[12]+muladd(affine[11],master.coords[3,i],
                       muladd(affine[10],master.coords[2,i],
                              affine[9]*master.coords[1,i])))
    for i in 1:nnodes(master)])
slave_set=Set(NTuple{3,Float64}[
    (slave.coords[1,i],slave.coords[2,i],slave.coords[3,i])
    for i in 1:nnodes(slave)])
slave_set==master_set ||
    error("Tessella periodic curved slave is not a bitwise affine copy")

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
        "curved-surface differential requires Gmsh 4.15.2, got " *
        gmsh.GMSH_API_VERSION)
    gmsh.open(GEO)
    gmsh.model.mesh.generate(2)
    gmsh_counts=Dict{Int,Int}()
    for tag in (1,2,3)
        _,elements,_=gmsh.model.mesh.getElements(2,tag)
        triangles=sum(length,elements;init=0)
        triangles>0 || error("Gmsh produced no triangles on surface $tag")
        gmsh_counts[tag]=triangles
    end
    gmsh.clear()
    gmsh.open(PERIODIC_GEO)
    gmsh.model.mesh.generate(2)
    gmsh_master,gmsh_slaves,gmsh_masters,gmsh_affine=
        gmsh.model.mesh.getPeriodicNodes(2,2)
    gmsh_master==1 ||
        error("Gmsh periodic curved master is $gmsh_master, expected 1")
    length(gmsh_slaves)==length(gmsh_masters)>0 ||
        error("Gmsh produced no periodic surface node pairs")
    gmsh_node_tags,gmsh_node_coordinates,_=gmsh.model.mesh.getNodes()
    gmsh_coordinates=Dict{UInt64,NTuple{3,Float64}}()
    for (index,tag) in pairs(gmsh_node_tags)
        gmsh_coordinates[UInt64(tag)]=(
            gmsh_node_coordinates[3index-2],
            gmsh_node_coordinates[3index-1],
            gmsh_node_coordinates[3index])
    end
    max_gmsh_error=maximum(zip(gmsh_slaves,gmsh_masters);init=0.0) do (s,t)
        slave=gmsh_coordinates[UInt64(s)]
        master=gmsh_coordinates[UInt64(t)]
        hypot(slave[1]-master[1]-3.0,slave[2]-master[2],slave[3]-master[3])
    end
    # Gmsh copies the master's parameter list and re-evaluates on the slave
    # curve — its own OCC/builtin arc evaluators leave ~1e-9 residue.
    max_gmsh_error<=1e-8 || error(
        "Gmsh periodic curved node error is $max_gmsh_error")
    _,gmsh_slave_elements,_=gmsh.model.mesh.getElements(2,2)
    gmsh_slave_triangles=sum(length,gmsh_slave_elements;init=0)
    gmsh_slave_triangles>0 ||
        error("Gmsh produced no triangles on the periodic slave surface")
    println("GMSH_PARITY_CURVED_SURFACE_OK gmsh=$(gmsh.GMSH_API_VERSION) " *
            "tessella_disk=$(ntris(disk)) tessella_annulus=$(ntris(annulus)) " *
            "tessella_segment=$(ntris(segment)) " *
            "tessella_periodic=$(ntris(master))+$(ntris(slave)) " *
            "rim=$rim hole_rim=$hole_rim " *
            "gmsh_surfs=$(gmsh_counts[1])+$(gmsh_counts[2])+$(gmsh_counts[3]) " *
            "gmsh_periodic_pairs=$(length(gmsh_slaves)) " *
            "gmsh_slave_tris=$gmsh_slave_triangles " *
            "gmsh_max_error=$max_gmsh_error")
finally
    gmsh.finalize()
end
