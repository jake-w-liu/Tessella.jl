#!/usr/bin/env julia
# P6: Tessella translation-periodic node correspondence vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: Mesh,validate,mesh_crc

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
const GEO=joinpath(@__DIR__,"periodic_translation.geo")
gmsh.initialize(["gmsh","-v","0"])
try
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "periodic differential requires Gmsh 4.15.2, got $(gmsh.GMSH_API_VERSION)")
    gmsh.open(GEO)
    gmsh.model.mesh.generate(2)
    master_entity,slave_tags,master_tags,affine=
        gmsh.model.mesh.getPeriodicNodes(1,2)
    master_entity==4 || error("Gmsh periodic master curve is $master_entity, expected 4")
    length(slave_tags)==length(master_tags)==5 || error(
        "Gmsh periodic curve pair count is $(length(slave_tags)), expected 5")
    expected_affine=[1.0,0.0,0.0,1.0,
                     0.0,1.0,0.0,0.0,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0]
    length(affine)==16 || error("Gmsh periodic affine transform is not 4×4")
    maximum(abs.(affine.-expected_affine))<=1e-14 || error(
        "Gmsh periodic affine transform is $affine")

    node_tags,node_coordinates,_=gmsh.model.mesh.getNodes()
    coordinates=Dict{Int,NTuple{3,Float64}}()
    for (i,tag) in enumerate(node_tags)
        coordinates[Int(tag)]=(node_coordinates[3i-2],node_coordinates[3i-1],
                               node_coordinates[3i])
    end
    all(tag->haskey(coordinates,Int(tag)),master_tags) || error(
        "Gmsh omitted a periodic master node from getNodes")
    all(tag->haskey(coordinates,Int(tag)),slave_tags) || error(
        "Gmsh omitted a periodic slave node from getNodes")
    order=sortperm(eachindex(master_tags);by=i->coordinates[Int(master_tags[i])][2])
    n=length(order)
    tessella_coordinates=Matrix{Float64}(undef,3,2n)
    max_gmsh_error=0.0
    for (column,pair_index) in enumerate(order)
        master=coordinates[Int(master_tags[pair_index])]
        slave=coordinates[Int(slave_tags[pair_index])]
        expected=(master[1]+1.0,master[2],master[3])
        max_gmsh_error=max(max_gmsh_error,
            hypot(slave[1]-expected[1],slave[2]-expected[2],slave[3]-expected[3]))
        tessella_coordinates[:,column].=master
        tessella_coordinates[:,n+column].=slave
    end
    max_gmsh_error<=1e-11 || error(
        "Gmsh periodic node correspondence error is $max_gmsh_error")

    segments=Matrix{Int32}(undef,2,2(n-1))
    for i in 1:n-1
        segments[:,i].=(Int32(i),Int32(i+1))
        segments[:,n-1+i].=(Int32(n+i),Int32(n+i+1))
    end
    input=Mesh(tessella_coordinates;segs=segments)
    output=periodic_identify(input,(1.0,0.0,0.0),collect(1:n),collect(n+1:2n);
                             atol=1e-11)
    validate(output).ok || error("Tessella periodic output is invalid")
    output.segs==input.segs || error("Tessella periodic operation changed connectivity")
    for i in 1:n
        output.coords[:,n+i]==output.coords[:,i].+[1.0,0.0,0.0] || error(
            "Tessella periodic pair $i was not snapped exactly")
    end
    println("GMSH_PARITY_PERIODIC_OK gmsh=$(gmsh.GMSH_API_VERSION) pairs=$n "*
            "gmsh_max_error=$max_gmsh_error sha=$(mesh_crc(output).sha)")
finally
    gmsh.finalize()
end
