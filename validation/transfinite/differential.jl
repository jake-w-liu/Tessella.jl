#!/usr/bin/env julia
# In-memory differential for four-sided planar transfinite patches against the
# pinned Gmsh 4.15.2 public API. No mesh or geometry file is created.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: nsegs,nnodes,ntris,node,validate

if !isdefined(Tessella,:Transfinite)
    Base.include(Tessella, joinpath(
        @__DIR__, "..", "..", "src", "structured", "Transfinite.jl"))
end
using Tessella.Transfinite: mesh_transfinite_patch

const TARGET_GMSH_VERSION="4.15.2"

function find_gmsh_executable()
    explicit=get(ENV,"GMSH_EXECUTABLE","")
    if !isempty(explicit)
        isfile(explicit) || error("GMSH_EXECUTABLE does not name a file: $explicit")
        return realpath(explicit)
    end
    executable=Sys.which("gmsh")
    executable!==nothing && return realpath(executable)
    fallback="/opt/homebrew/bin/gmsh"
    isfile(fallback) && return realpath(fallback)
    error("Gmsh 4.15.2 is required; install it or set GMSH_EXECUTABLE")
end

function find_gmsh_api(executable)
    explicit=get(ENV,"GMSH_JULIA_API","")
    if !isempty(explicit)
        isfile(explicit) || error("GMSH_JULIA_API does not name a file: $explicit")
        return realpath(explicit)
    end
    prefix=dirname(dirname(executable))
    candidates=(joinpath(prefix,"lib","gmsh.jl"),
                joinpath(prefix,"lib64","gmsh.jl"),
                "/opt/homebrew/lib/gmsh.jl",
                "/opt/homebrew/opt/gmsh/lib/gmsh.jl",
                "/usr/local/opt/gmsh/lib/gmsh.jl")
    for candidate in candidates
        isfile(candidate) && return realpath(candidate)
    end
    error("could not locate gmsh.jl for $executable; set GMSH_JULIA_API")
end

const GMSH_EXECUTABLE=find_gmsh_executable()
const GMSH_CLI_VERSION=strip(read(`$GMSH_EXECUTABLE --version`,String))
(GMSH_CLI_VERSION==TARGET_GMSH_VERSION ||
 startswith(GMSH_CLI_VERSION,TARGET_GMSH_VERSION*"-")) ||
    error("expected Gmsh $TARGET_GMSH_VERSION, got $GMSH_CLI_VERSION")
const GMSH_API_FILE=find_gmsh_api(GMSH_EXECUTABLE)
include(GMSH_API_FILE)

@inline _distance(a,b)=hypot(a[1]-b[1],a[2]-b[2],a[3]-b[3])

function curve_points(curve)
    tags,coordinates,parameters=gmsh.model.mesh.getNodes(1,curve,true,true)
    length(tags)==length(parameters) || error(
        "Gmsh curve $curve did not return one parameter per node")
    order=sortperm(parameters)
    return [(coordinates[3index-2],coordinates[3index-1],coordinates[3index])
            for index in order]
end

function canonical_triangles(connectivity)
    length(connectivity)%3==0 || error("triangle connectivity is not divisible by 3")
    result=NTuple{3,Int32}[]
    for index in 1:3:length(connectivity)
        values=sort(Int32[connectivity[index],connectivity[index+1],
                          connectivity[index+2]])
        push!(result,(values[1],values[2],values[3]))
    end
    sort!(result)
end

function canonical_segments(connectivity)
    length(connectivity)%2==0 || error("segment connectivity is not divisible by 2")
    result=NTuple{2,Int32}[]
    for index in 1:2:length(connectivity)
        a=Int32(connectivity[index]);b=Int32(connectivity[index+1])
        push!(result,a<b ? (a,b) : (b,a))
    end
    sort!(result)
end

function add_curved_patch(arrangement)
    gmsh.clear();gmsh.model.add("transfinite_"*arrangement)
    p1=gmsh.model.geo.addPoint(0.,0.,0.)
    pb1=gmsh.model.geo.addPoint(1.,-.4,0.)
    pb2=gmsh.model.geo.addPoint(2.5,-.7,0.)
    pb3=gmsh.model.geo.addPoint(3.2,-.2,0.)
    p2=gmsh.model.geo.addPoint(4.,0.,0.)
    pr1=gmsh.model.geo.addPoint(4.4,1.,0.)
    pr2=gmsh.model.geo.addPoint(4.2,2.,0.)
    p3=gmsh.model.geo.addPoint(4.,3.,0.)
    pt1=gmsh.model.geo.addPoint(3.,3.6,0.)
    pt2=gmsh.model.geo.addPoint(1.3,3.2,0.)
    p4=gmsh.model.geo.addPoint(0.,3.,0.)
    pl1=gmsh.model.geo.addPoint(-.3,2.1,0.)
    pl2=gmsh.model.geo.addPoint(-.6,.7,0.)
    curves=(gmsh.model.geo.addSpline([p1,pb1,pb2,pb3,p2]),
            gmsh.model.geo.addSpline([p2,pr1,pr2,p3]),
            gmsh.model.geo.addSpline([p3,pt1,pt2,p4]),
            gmsh.model.geo.addSpline([p4,pl1,pl2,p1]))
    loop=gmsh.model.geo.addCurveLoop(collect(curves))
    surface=gmsh.model.geo.addPlaneSurface([loop])
    for (curve,count) in zip(curves,(5,4,5,4))
        gmsh.model.geo.mesh.setTransfiniteCurve(curve,count)
    end
    gmsh.model.geo.mesh.setTransfiniteSurface(
        surface,arrangement,[p1,p2,p3,p4])
    gmsh.model.geo.synchronize();gmsh.model.mesh.generate(2)
    return curves,surface
end

function gmsh_to_tessella_node_map(mesh,surface)
    types,_,element_nodes=gmsh.model.mesh.getElements(2,surface)
    triangle_position=findfirst(==(Int32(2)),types)
    triangle_position===nothing && error("Gmsh emitted no type-2 triangles")
    tags=sort!(unique(element_nodes[triangle_position]))
    length(tags)==nnodes(mesh) || error(
        "node-count mismatch: Gmsh $(length(tags)), Tessella $(nnodes(mesh))")
    used=falses(nnodes(mesh));mapping=Dict{UInt64,Int32}();maximum_error=0.0
    coordinates=Dict(tag=>gmsh.model.mesh.getNode(tag)[1] for tag in tags)
    scale=maximum((maximum(abs,point;init=0.0) for point in values(coordinates));init=1.0)
    tolerance=256eps(Float64)*max(scale,1.0)
    for source in eachindex(tags)
        raw=coordinates[tags[source]];point=(raw[1],raw[2],raw[3])
        best=0;best_error=Inf
        for destination in 1:nnodes(mesh)
            used[destination] && continue
            error=_distance(point,node(mesh,destination))
            if error<best_error
                best=destination;best_error=error
            end
        end
        best!=0&&best_error<=tolerance || error(
            "no Tessella node matches Gmsh node $(tags[source]); nearest error=$best_error, tolerance=$tolerance")
        mapping[tags[source]]=Int32(best);used[best]=true
        maximum_error=max(maximum_error,best_error)
    end
    all(used) || error("some Tessella nodes were not matched to Gmsh nodes")
    return mapping,maximum_error
end

function check_arrangement(arrangement,symbol)
    curves,surface=add_curved_patch(arrangement)
    sides=map(curve_points,curves)
    mesh=mesh_transfinite_patch(sides...;arrangement=symbol,
                                face_tag=21,side_tags=(11,12,13,14))
    validate(mesh).ok || error("Tessella $arrangement patch did not validate")
    (nnodes(mesh),nsegs(mesh),ntris(mesh))==(20,14,24) || error(
        "unexpected Tessella $arrangement counts")
    mapping,maximum_error=gmsh_to_tessella_node_map(mesh,surface)

    types,_,element_nodes=gmsh.model.mesh.getElements(2,surface)
    triangle_position=findfirst(==(Int32(2)),types)
    triangle_position===nothing && error("Gmsh $arrangement emitted no type-2 triangles")
    mapped_triangles=Int32[mapping[tag] for tag in element_nodes[triangle_position]]
    canonical_triangles(mapped_triangles)==canonical_triangles(vec(mesh.tris)) ||
        error("$arrangement triangle connectivity differs from Gmsh")

    mapped_segments=Int32[]
    for curve in curves
        line_types,_,line_nodes=gmsh.model.mesh.getElements(1,curve)
        line_position=findfirst(==(Int32(1)),line_types)
        line_position===nothing && error("Gmsh curve $curve emitted no type-1 lines")
        append!(mapped_segments,(mapping[tag] for tag in line_nodes[line_position]))
    end
    canonical_segments(mapped_segments)==canonical_segments(vec(mesh.segs)) ||
        error("$arrangement boundary connectivity differs from Gmsh")
    return maximum_error
end

function add_rectangle(counts;hole=false)
    outer=[gmsh.model.geo.addPoint(x,y,0.)
           for (x,y) in ((0.,0.),(3.,0.),(3.,2.),(0.,2.))]
    curves=[gmsh.model.geo.addLine(outer[i],outer[mod1(i+1,4)]) for i in 1:4]
    loops=[gmsh.model.geo.addCurveLoop(curves)]
    if hole
        inner=[gmsh.model.geo.addPoint(x,y,0.)
               for (x,y) in ((1.,.6),(2.,.6),(2.,1.4),(1.,1.4))]
        inner_curves=[gmsh.model.geo.addLine(inner[i],inner[mod1(i+1,4)]) for i in 1:4]
        push!(loops,gmsh.model.geo.addCurveLoop(inner_curves))
        for curve in inner_curves
            gmsh.model.geo.mesh.setTransfiniteCurve(curve,3)
        end
    end
    surface=gmsh.model.geo.addPlaneSurface(loops)
    for (curve,count) in zip(curves,counts)
        gmsh.model.geo.mesh.setTransfiniteCurve(curve,count)
    end
    gmsh.model.geo.mesh.setTransfiniteSurface(surface,"Left",outer)
    gmsh.model.geo.synchronize()
    return surface
end

function captured_generate()
    gmsh.logger.start();thrown=nothing;messages=String[]
    try
        try
            gmsh.model.mesh.generate(2)
        catch err
            thrown=sprint(showerror,err)
        end
        messages=gmsh.logger.get()
    finally
        gmsh.logger.stop()
    end
    return thrown,messages
end

function expect_argument_error(f,fragment)
    try
        f()
    catch err
        err isa ArgumentError || error("expected ArgumentError, got $(typeof(err)): $(sprint(showerror,err))")
        occursin(fragment,sprint(showerror,err)) || error(
            "error did not contain '$fragment': $(sprint(showerror,err))")
        return nothing
    end
    error("expected an ArgumentError containing '$fragment'")
end

function check_failure_behavior()
    gmsh.clear();gmsh.model.add("mismatched")
    surface=add_rectangle((4,3,5,3))
    thrown,messages=captured_generate()
    thrown===nothing || error("Gmsh mismatch unexpectedly threw: $thrown")
    any(occursin("non-matching number of nodes on opposite sides",message)
        for message in messages) || error("Gmsh mismatch warning was absent")
    types,tags,_=gmsh.model.mesh.getElements(2,surface)
    types==Int32[2] || error("Gmsh mismatch did not emit only first-order triangles")
    sum(length,tags;init=0)==23 || error("Gmsh mismatch emitted an unexpected element count")
    length(gmsh.model.mesh.getNodes()[1])==18 || error(
        "Gmsh mismatch emitted an unexpected node count")

    bottom=[(0.,0.,0.),(1.,0.,0.),(2.,0.,0.),(3.,0.,0.)]
    right=[(3.,0.,0.),(3.,1.,0.),(3.,2.,0.)]
    top=[(3.,2.,0.),(2.25,2.,0.),(1.5,2.,0.),(.75,2.,0.),(0.,2.,0.)]
    left=[(0.,2.,0.),(0.,1.,0.),(0.,0.,0.)]
    expect_argument_error(
        ()->mesh_transfinite_patch(bottom,right,top,left),"non-matching node counts")

    gmsh.clear();gmsh.model.add("hole")
    surface=add_rectangle((4,3,4,3);hole=true)
    thrown,messages=captured_generate()
    thrown!==nothing&&occursin("is transfinite but has 1 hole",thrown) || error(
        "Gmsh hole error differed: $thrown")
    isempty(gmsh.model.mesh.getElements(2,surface)[1]) || error(
        "Gmsh emitted surface elements after the transfinite-hole error")
    return nothing
end

gmsh.initialize([GMSH_EXECUTABLE,"-nopopup"],false,false)
try
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.option.setNumber("Mesh.Smoothing",1)
    gmsh.option.setNumber("Mesh.QuasiTransfinite",0)
    gmsh.option.setNumber("Mesh.RecombineAll",0)
    api_version=gmsh.option.getString("General.Version")
    (api_version==TARGET_GMSH_VERSION ||
     startswith(api_version,TARGET_GMSH_VERSION*"-")) || error(
        "expected Gmsh API $TARGET_GMSH_VERSION, got $api_version")
    cases=(("Left",:left),("Right",:right),
           ("AlternateLeft",:alternate_left),
           ("AlternateRight",:alternate_right))
    errors=Float64[]
    for (arrangement,symbol) in cases
        push!(errors,check_arrangement(arrangement,symbol))
    end
    check_failure_behavior()
    println("TRANSFINITE_DIFFERENTIAL_OK gmsh=$api_version arrangements=$(length(cases)) " *
            "coordinate_samples=$(20length(cases)) max_node_error=$(maximum(errors)) " *
            "mismatch_nodes=18 mismatch_triangles=23 hole_elements=0")
finally
    gmsh.finalize()
end
