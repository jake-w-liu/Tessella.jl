#!/usr/bin/env julia
# P6: Tessella IGES-128 bilinear patch vs analytic centre and vs Gmsh OCC BSpline surface.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.NURBS: NURBSCurve, NURBSSurface, nurbs_eval

patch=NURBSSurface(1,1,[0,0,1,1],[0,0,1,1],
                   [(0.0,0.0,0.0) (0.0,1.0,0.0); (1.0,0.0,0.0) (1.0,1.0,0.0)])
export_curve=NURBSCurve(2,[0,0,0,1,1,1],
                        [(0.0,0.0,2.0),(0.5,0.5,2.0),(1.0,0.0,2.0)])
imported=mktemp() do path,io
    close(io)
    export_iges_nurbs(path,[patch])
    only(import_nurbs_iges(path))
end
q=nurbs_eval(imported,0.5,0.5)
hypot(q[1]-0.5,q[2]-0.5,q[3]-0.0)<=1e-14 || error("Tessella NURBS centre $q != (0.5,0.5,0)")
for (u,v,target) in ((0.0,0.0,(0.0,0.0,0.0)),(1.0,0.0,(1.0,0.0,0.0)),
                     (0.0,1.0,(0.0,1.0,0.0)),(1.0,1.0,(1.0,1.0,0.0)))
    p=nurbs_eval(imported,u,v)
    hypot(p[1]-target[1],p[2]-target[2],p[3]-target[3])<=1e-14 ||
        error("Tessella NURBS corner ($u,$v) = $p != $target")
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
gmsh.initialize(["gmsh","-v","0"])
try
    p11=gmsh.model.occ.addPoint(0,0,0)
    p21=gmsh.model.occ.addPoint(1,0,0)
    p12=gmsh.model.occ.addPoint(0,1,0)
    p22=gmsh.model.occ.addPoint(1,1,0)
    gmsh.model.occ.addBSplineSurface(Int32[p11,p21,p12,p22], 2, -1, 1, 1,
                                    Float64[], Float64[0.0,1.0], Float64[0.0,1.0],
                                    Int32[2,2], Int32[2,2])
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
    gmsh.model.mesh.generate(2)
    function surface_metrics(label)
        types,tags,conn=gmsh.model.mesh.getElements(2)
        isempty(types) && error("$label has no surface elements")
        nodes,coord,_=gmsh.model.mesh.getNodes()
        lookup=Dict{Int,NTuple{2,Float64}}()
        for (k,tag) in enumerate(nodes)
            lookup[Int(tag)]=(coord[3k-2],coord[3k-1])
        end
        acc=0.0; ntri=0
        for (t,tg,cn) in zip(types,tags,conn)
            nper=Int(length(cn)÷length(tg))
            t==2 && (ntri+=length(tg))
            for i in 1:length(tg)
                ids=cn[(i-1)*nper+1:i*nper]
                xs=Float64[]; ys=Float64[]
                for id in ids
                    xy=lookup[Int(id)]
                    push!(xs,xy[1]); push!(ys,xy[2])
                end
                a=0.0
                n=length(xs)
                for k in 1:n
                    k2=k==n ? 1 : k+1
                    a+=xs[k]*ys[k2]-ys[k]*xs[k2]
                end
                acc+=abs(a)/2
            end
        end
        return acc,ntri,length(nodes)
    end
    area,gmsh_tris,_=surface_metrics("Gmsh OCC NURBS surface")
    abs(area-1)<=1e-8 || error("Gmsh NURBS surface area $area != 1")
    gmsh_iges_center=mktempdir() do directory
        path=joinpath(directory,"gmsh_nurbs.iges")
        gmsh.write(path)
        objects=import_nurbs_iges(path)
        surfaces=filter(object->object isa NURBSSurface,objects)
        length(surfaces)==1 || error(
            "Tessella did not recover one IGES 128 surface from Gmsh")
        nurbs_eval(only(surfaces),0.5,0.5)
    end
    hypot(gmsh_iges_center[1]-0.5,gmsh_iges_center[2]-0.5,gmsh_iges_center[3])<=1e-12 ||
        error("Tessella imported Gmsh's IGES centre as $gmsh_iges_center")

    exported_area,exported_tris,exported_nodes=mktempdir() do directory
        path=joinpath(directory,"tessella_nurbs.iges")
        export_iges_nurbs(path,[export_curve,patch])
        all(ncodeunits(line)==80 for line in readlines(path)) ||
            error("Tessella IGES export contains a non-80-column record")
        gmsh.clear()
        gmsh.open(path)
        imported_entities=gmsh.model.getEntities()
        any(entity->entity[1]==2,imported_entities) ||
            error("Gmsh did not import a surface from Tessella's IGES 128/144 export")
        count(entity->entity[1]==1,imported_entities)>=5 ||
            error("Gmsh did not import Tessella's standalone IGES 126 curve")
        gmsh.option.setNumber("Mesh.MeshSizeMax",0.5)
        gmsh.model.mesh.generate(2)
        surface_metrics("Tessella-exported IGES surface")
    end
    abs(exported_area-1)<=1e-8 || error(
        "Tessella-exported IGES surface area $exported_area != 1")
    println("GMSH_PARITY_NURBS_OK gmsh=$(gmsh.GMSH_API_VERSION) tessella_centre=0.5 " *
            "gmsh_area=1 gmsh_tris=$gmsh_tris gmsh_iges_centre=0.5 export_area=1 " *
            "export_nodes=$exported_nodes export_tris=$exported_tris")
finally
    gmsh.finalize()
end
