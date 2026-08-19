# Differential validation of Tessella's Gmsh-compatible Distance→Threshold,
# Box, Ball, Cylinder, and Frustum size fields against an installed Gmsh Julia API.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.Mesh1D: mesh_segment
using Statistics: median

function find_gmsh_api()
    explicit=get(ENV,"GMSH_JULIA_API","")
    !isempty(explicit) && isfile(explicit) && return explicit
    executable=Sys.which("gmsh")
    executable===nothing && error("gmsh is not on PATH; set GMSH_JULIA_API to gmsh.jl")
    prefix=dirname(dirname(realpath(executable)))
    candidates=(joinpath(prefix,"lib","gmsh.jl"),
                joinpath(prefix,"lib64","gmsh.jl"),
                "/opt/homebrew/opt/gmsh/lib/gmsh.jl",
                "/usr/local/opt/gmsh/lib/gmsh.jl")
    for path in candidates
        isfile(path) && return path
    end
    error("could not locate gmsh.jl next to $(realpath(executable)); set GMSH_JULIA_API")
end

include(find_gmsh_api())

function common_gmsh_model(name)
    gmsh.clear()
    gmsh.model.add(name)
    gmsh.model.geo.addPoint(-2.0,0.0,0.0,1.0,1)
    gmsh.model.geo.addPoint( 2.0,0.0,0.0,1.0,2)
    gmsh.model.geo.addLine(1,2,1)
end

function finish_gmsh_line()
    gmsh.model.geo.synchronize()
    for (option,value) in (("Mesh.MeshSizeExtendFromBoundary",0.0),
                           ("Mesh.MeshSizeFromPoints",0.0),
                           ("Mesh.MeshSizeFromCurvature",0.0),
                           ("Mesh.MeshSizeMin",0.01),
                           ("Mesh.MeshSizeMax",1.0))
        gmsh.option.setNumber(option,value)
    end
    gmsh.model.mesh.generate(1)
    _,coordinates,_=gmsh.model.mesh.getNodes(1,1,true,false)
    return sort!(unique(coordinates[1:3:end]))
end

function gmsh_threshold_line()
    common_gmsh_model("threshold_line")
    gmsh.model.geo.addPoint(0.0,0.0,0.0,1.0,3)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.field.add("Distance",1)
    gmsh.model.mesh.field.setNumbers(1,"PointsList",[3])
    gmsh.model.mesh.field.add("Threshold",2)
    for (option,value) in (("InField",1.0),("DistMin",0.0),("DistMax",1.0),
                           ("SizeMin",0.1),("SizeMax",0.5))
        gmsh.model.mesh.field.setNumber(2,option,value)
    end
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    return finish_gmsh_line()
end

function gmsh_box_line()
    common_gmsh_model("box_line")
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.field.add("Box",1)
    for (option,value) in (("VIn",0.1),("VOut",0.5),
                           ("XMin",-0.5),("XMax",0.5),
                           ("YMin",-1.0),("YMax",1.0),
                           ("ZMin",-1.0),("ZMax",1.0),("Thickness",0.5))
        gmsh.model.mesh.field.setNumber(1,option,value)
    end
    gmsh.model.mesh.field.setAsBackgroundMesh(1)
    return finish_gmsh_line()
end

function gmsh_ball_line()
    common_gmsh_model("ball_line")
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.field.add("Ball",1)
    for (option,value) in (("VIn",0.1),("VOut",0.5),
                           ("XCenter",0.0),("YCenter",0.0),("ZCenter",0.0),
                           ("Radius",0.5),("Thickness",0.5))
        gmsh.model.mesh.field.setNumber(1,option,value)
    end
    gmsh.model.mesh.field.setAsBackgroundMesh(1)
    return finish_gmsh_line()
end

function gmsh_cylinder_line()
    common_gmsh_model("cylinder_line")
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.field.add("Cylinder",1)
    for (option,value) in (("VIn",0.1),("VOut",0.5),
                           ("XCenter",0.0),("YCenter",0.0),("ZCenter",0.0),
                           ("XAxis",0.5),("YAxis",0.0),("ZAxis",0.0),
                           ("Radius",1.0))
        gmsh.model.mesh.field.setNumber(1,option,value)
    end
    gmsh.model.mesh.field.setAsBackgroundMesh(1)
    return finish_gmsh_line()
end

function gmsh_frustum_line()
    common_gmsh_model("frustum_line")
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.field.add("Frustum",1)
    for (option,value) in (("X1",0.0),("Y1",0.0),("Z1",-1.0),
                           ("X2",0.0),("Y2",0.0),("Z2",1.0),
                           ("InnerR1",0.0),("OuterR1",1.0),
                           ("InnerR2",0.0),("OuterR2",1.0),
                           ("InnerV1",0.1),("OuterV1",0.5),
                           ("InnerV2",0.1),("OuterV2",0.5))
        gmsh.model.mesh.field.setNumber(1,option,value)
    end
    gmsh.model.mesh.field.add("Box",2)
    gmsh.model.mesh.field.setNumber(2,"VIn",0.5)
    gmsh.model.mesh.field.setNumber(2,"VOut",0.5)
    gmsh.model.mesh.field.add("Min",3)
    gmsh.model.mesh.field.setNumbers(3,"FieldsList",[1,2])
    gmsh.model.mesh.field.setAsBackgroundMesh(3)
    return finish_gmsh_line()
end

function tessella_line(field)
    points,_=mesh_segment((-2.0,0.0,0.0),(2.0,0.0,0.0),field)
    return [point[1] for point in points]
end

function spacing_metrics(nodes)
    gaps=diff(nodes)
    midpoints=(nodes[1:end-1]+nodes[2:end])./2
    return (count=length(nodes), minimum=minimum(gaps), maximum=maximum(gaps),
            near=median(gaps[abs.(midpoints).<0.4]),
            far=median(gaps[abs.(midpoints).>1.2]))
end

gmsh.initialize(["gmsh","-v","0"])
try
    gt=gmsh_threshold_line()
    distance=DistanceField(points=[(0.0,0.0,0.0)])
    threshold=ThresholdField(distance;dist_min=0.0,dist_max=1.0,
                             size_min=0.1,size_max=0.5)
    tt=tessella_line(threshold)
    gm=spacing_metrics(gt);tm=spacing_metrics(tt)
    @assert abs(gm.count-tm.count)<=1 "Threshold node-count mismatch: gmsh=$gm Tessella=$tm"
    @assert 0.75<=tm.minimum/gm.minimum<=1.25 "Threshold fine-spacing mismatch: gmsh=$gm Tessella=$tm"
    @assert 0.90<=tm.maximum/gm.maximum<=1.10 "Threshold coarse-spacing mismatch: gmsh=$gm Tessella=$tm"

    gb=gmsh_box_line()
    box=BoxField(-0.5,0.5,-1.0,1.0,-1.0,1.0;
                 vin=0.1,vout=0.5,thickness=0.5)
    tb=tessella_line(box)
    gbm=spacing_metrics(gb);tbm=spacing_metrics(tb)
    @assert abs(gbm.count-tbm.count)<=4 "Box node-count mismatch: gmsh=$gbm Tessella=$tbm"
    @assert 0.95<=tbm.near/gbm.near<=1.05 "Box interior-spacing mismatch: gmsh=$gbm Tessella=$tbm"
    @assert 0.80<=tbm.far/gbm.far<=1.25 "Box exterior-spacing mismatch: gmsh=$gbm Tessella=$tbm"
    @assert tbm.near<tbm.far && gbm.near<gbm.far "Box field did not refine its interior"

    gba=gmsh_ball_line()
    ball=BallField((0.0,0.0,0.0),0.5;vin=0.1,vout=0.5,thickness=0.5)
    tba=tessella_line(ball)
    gbam=spacing_metrics(gba);tbam=spacing_metrics(tba)
    @assert abs(gbam.count-tbam.count)<=4 "Ball node-count mismatch: gmsh=$gbam Tessella=$tbam"
    @assert 0.90<=tbam.near/gbam.near<=1.10 "Ball interior-spacing mismatch: gmsh=$gbam Tessella=$tbam"
    @assert 0.80<=tbam.far/gbam.far<=1.25 "Ball exterior-spacing mismatch: gmsh=$gbam Tessella=$tbam"
    @assert tbam.near<tbam.far && gbam.near<gbam.far "Ball field did not refine its interior"

    gc=gmsh_cylinder_line()
    cylinder=CylinderField((0.0,0.0,0.0),(0.5,0.0,0.0),1.0;
                           vin=0.1,vout=0.5)
    tc=tessella_line(cylinder)
    gcm=spacing_metrics(gc);tcm=spacing_metrics(tc)
    @assert abs(gcm.count-tcm.count)<=4 "Cylinder node-count mismatch: gmsh=$gcm Tessella=$tcm"
    @assert 0.75<=tcm.near/gcm.near<=1.25 "Cylinder interior-spacing mismatch: gmsh=$gcm Tessella=$tcm"
    @assert 0.75<=tcm.far/gcm.far<=1.25 "Cylinder exterior-spacing mismatch: gmsh=$gcm Tessella=$tcm"
    @assert tcm.near<tcm.far && gcm.near<gcm.far "Cylinder field did not refine its interior"

    gf=gmsh_frustum_line()
    frustum=MinSize((FrustumField((0.0,0.0,-1.0),(0.0,0.0,1.0);
        inner_r1=0.0,outer_r1=1.0,inner_r2=0.0,outer_r2=1.0,
        inner_v1=0.1,outer_v1=0.5,inner_v2=0.1,outer_v2=0.5),ConstantSize(0.5)))
    tf=tessella_line(frustum)
    gfm=spacing_metrics(gf);tfm=spacing_metrics(tf)
    @assert abs(gfm.count-tfm.count)<=4 "Frustum node-count mismatch: gmsh=$gfm Tessella=$tfm"
    @assert 0.75<=tfm.near/gfm.near<=1.25 "Frustum inner-spacing mismatch: gmsh=$gfm Tessella=$tfm"
    @assert 0.80<=tfm.far/gfm.far<=1.25 "Frustum outer-spacing mismatch: gmsh=$gfm Tessella=$tfm"
    @assert tfm.near<tfm.far && gfm.near<gfm.far "Frustum field did not grade radially"

    println("SIZE_FIELD_DIFFERENTIAL_OK gmsh=$(gmsh.GMSH_API_VERSION)")
    println("  Threshold: gmsh=$gm Tessella=$tm")
    println("  Box:       gmsh=$gbm Tessella=$tbm")
    println("  Ball:      gmsh=$gbam Tessella=$tbam")
    println("  Cylinder:  gmsh=$gcm Tessella=$tcm")
    println("  Frustum:   gmsh=$gfm Tessella=$tfm")
finally
    gmsh.finalize()
end
