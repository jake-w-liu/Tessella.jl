#!/usr/bin/env julia
# P6: Tessella IGES-128 bilinear patch vs analytic centre and vs Gmsh OCC BSpline surface.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."); io=devnull)
using Tessella
using Tessella.NURBS: NURBSSurface, nurbs_eval

patch=NURBSSurface(1,1,[0,0,1,1],[0,0,1,1],
                   [(0.0,0.0,0.0) (0.0,1.0,0.0); (1.0,0.0,0.0) (1.0,1.0,0.0)])
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
    types, tags, conn = gmsh.model.mesh.getElements(2)
    isempty(types) && error("Gmsh NURBS surface has no elements")
    nodes, coord, _ = gmsh.model.mesh.getNodes()
    lookup=Dict{Int,NTuple{2,Float64}}()
    for (k,tag) in enumerate(nodes)
        lookup[Int(tag)]=(coord[3k-2],coord[3k-1])
    end
    area, gmsh_tris = let acc=0.0, ntri=0
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
        acc, ntri
    end
    abs(area-1)<=1e-8 || error("Gmsh NURBS surface area $area != 1")
    println("GMSH_PARITY_NURBS_OK gmsh=$(gmsh.GMSH_API_VERSION) tessella_centre=0.5 gmsh_area=1 gmsh_tris=$gmsh_tris")
finally
    gmsh.finalize()
end
