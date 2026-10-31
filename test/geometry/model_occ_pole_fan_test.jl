# Pole-fan meshing for degenerate-edge OCC faces.
#
# A degenerate OCC edge is a single 3-D vertex (a sphere pole or a cone
# apex) whose stored pcurve still bounds the face's (u,v) domain — on the
# sphere it is a full collapsed row at constant latitude, on a cone apex the
# zero-radius row. `_model_occ_curved_surface_mesh` emits that row as an
# ordinary uv boundary chain whose samples all lift through the canonical
# pole vertex: the bitwise 3-D weld folds the row into the fan apex, CDT
# cells touching two row vertices collapse out of the lifted mesh, and the
# surviving cells fan the apex against the first non-collapsed ring. The
# collapsed row is deliberately not u-periodic unwrapped — its two endpoint
# representatives legitimately differ by the full 2π domain period — and a
# degenerate pcurve that does not lift to its recorded pole vertex is a
# corrupt record, rejected loudly rather than meshing.

using Test
using Tessella
const M=Tessella.Model
const MT=Tessella.MeshTypes

_error(f)=try;f();nothing catch err;err end

# Number of triangles incident to one welded node — the fan degree at a pole.
function _fan_degree(mesh,node)
    degree=0
    for c in axes(mesh.tris,2)
        if mesh.tris[1,c]==node || mesh.tris[2,c]==node ||
           mesh.tris[3,c]==node
            degree+=1
        end
    end
    return degree
end

# The node indices sitting bitwise on a canonical model point.
function _nodes_at(mesh,p)
    return [j for j in axes(mesh.coords,2) if
            isequal((mesh.coords[1,j],mesh.coords[2,j],mesh.coords[3,j]),p)]
end

@testset "full sphere face meshes as a two-pole fan" begin
    m=GeoModel()
    M._model_occ_add_sphere!(m,0.0,0.0,0.0,1.5;tag=1)
    M._model_occ_synchronize!(m)
    sface=findfirst(==(:sphere),m.surface_types)
    degs=sort!([c for (c,k) in m.curve_types if k===:degenerate])
    @test length(degs)==2          # north + south collapsed pole edges
    mesh=M.mesh_model_surface(m,sface)
    @test MT.validate(mesh).ok
    @test MT.ntris(mesh)>0
    # Every lifted node lies on the analytic sphere, including the fanned
    # apex: off-surface lift or a stray seam copy would show up here.
    for j in axes(mesh.coords,2)
        @test abs(sqrt(mesh.coords[1,j]^2+mesh.coords[2,j]^2+
                       mesh.coords[3,j]^2)-1.5)<=1e-9
    end
    # Each collapsed uv row welds to ONE apex node — bitwise on the model
    # vertex — with a genuine fan of triangles around it.
    for deg in degs
        pole=m.points[m.curves[deg][1]]
        nodes=_nodes_at(mesh,pole)
        @test length(nodes)==1
        @test _fan_degree(mesh,nodes[1])>=4
    end
    # uv-CCW cells lifted through the right-handed parametrization face
    # outward from the sphere centre.
    c=(0.0,0.0,0.0)
    for t in axes(mesh.tris,2)
        a=mesh.coords[:,mesh.tris[1,t]];b=mesh.coords[:,mesh.tris[2,t]]
        d=mesh.coords[:,mesh.tris[3,t]]
        nrm=(b[2]-a[2])*(d[3]-a[3])-(b[3]-a[3])*(d[2]-a[2]),
            (b[3]-a[3])*(d[1]-a[1])-(b[1]-a[1])*(d[3]-a[3]),
            (b[1]-a[1])*(d[2]-a[2])-(b[2]-a[2])*(d[1]-a[1])
        g=(a.+b.+d)./3 .-c
        @test nrm[1]*g[1]+nrm[2]*g[2]+nrm[3]*g[3]>0.0
    end
end

@testset "sphere sector face keeps its pole fan" begin
    # angle3=π: a half-sphere sector — the sphere face still carries both
    # collapsed pole rows (partial row spans stay supported) beside two
    # trimmed meridians.
    m=GeoModel()
    M._model_occ_add_sphere!(m,0.0,0.0,0.0,1.0;tag=1,angle3=pi)
    M._model_occ_synchronize!(m)
    sface=findfirst(==(:sphere),m.surface_types)
    mesh=M.mesh_model_surface(m,sface)
    @test MT.validate(mesh).ok
    @test MT.ntris(mesh)>0
    for j in axes(mesh.coords,2)
        @test abs(sqrt(mesh.coords[1,j]^2+mesh.coords[2,j]^2+
                       mesh.coords[3,j]^2)-1.0)<=1e-9
    end
    for deg in [c for (c,k) in m.curve_types if k===:degenerate]
        nodes=_nodes_at(mesh,m.points[m.curves[deg][1]])
        @test length(nodes)==1
        @test _fan_degree(mesh,nodes[1])>=3
    end
end

@testset "cone apex face meshes as a fan" begin
    # r2=0: the top rim collapses to a degenerate apex edge whose pcurve is
    # the v=len row of the cone lateral.
    m=GeoModel()
    M._model_occ_add_cone!(m,0.0,0.0,0.0,0.0,0.0,2.0,1.0,0.0;tag=1)
    M._model_occ_synchronize!(m)
    deg=first(c for (c,k) in m.curve_types if k===:degenerate)
    apex=m.points[m.curves[deg][1]]
    lat=findfirst(==(:cone),m.surface_types)
    mesh=M.mesh_model_surface(m,lat)
    @test MT.validate(mesh).ok
    @test MT.ntris(mesh)>0
    # radius(z) = r1 + (r2-r1)·z/h on the analytic cone.
    for j in axes(mesh.coords,2)
        @test abs(hypot(mesh.coords[1,j],mesh.coords[2,j])-
                  (1.0-mesh.coords[3,j]/2.0))<=1e-8
    end
    nodes=_nodes_at(mesh,apex)
    @test length(nodes)==1
    @test _fan_degree(mesh,nodes[1])>=3
end

@testset "corrupt degenerate pcurves still reject" begin
    # A nonzero-span degenerate pcurve that does not collapse onto its edge
    # vertex is a corrupt record, not a fan boundary — the kernel rejects
    # loudly instead of meshing a self-contradictory domain.
    m=GeoModel()
    M._model_occ_add_sphere!(m,0.0,0.0,0.0,1.5;tag=1)
    M._model_occ_synchronize!(m)
    sface=findfirst(==(:sphere),m.surface_types)
    deg=first(c for (c,k) in m.curve_types if k===:degenerate)
    sg=m.surface_geometry[sface]
    pcs=Dict{Int,NamedTuple}(sg.pcurves)
    # Equator row instead of the pole row: every lift lands on the equator,
    # 1.5 away from the pole vertex.
    pcs[deg]=(fwd=(lin2d=((0.0,0.0),(1.0,0.0)),),rev=nothing)
    m.surface_geometry[sface]=merge(sg,(pcurves=pcs,))
    err=_error(()->M.mesh_model_surface(m,sface))
    @test err isa ArgumentError
    @test occursin("does not collapse onto its pole vertex",
                   sprint(showerror,err))
end

@testset "explicit sphere shell volume meshes through the pole fan" begin
    # The compact `m.spheres` encoding drives a dedicated geodesic mesher;
    # dropping it (the state `_reconcile_curved_encodings!` leaves after an
    # independent sub-entity move) routes the shell through the generic
    # `_model_volume_surface_mesh` path — the same kernel `mesh_model_surface`
    # exercises — so the pole fan carries a sphere volume end to end.
    m=GeoModel()
    M._model_occ_add_sphere!(m,0.0,0.0,0.0,1.0;tag=1)
    M._model_occ_synchronize!(m)
    delete!(m.spheres,1)
    mesh=M.mesh_model_volume(m,1)
    @test MT.ntets(mesh)>0
    @test MT.validate(mesh).ok
    # Boundary nodes sit on the sphere; the single interior Steiner vertex
    # is the only node off the radius.
    on=sum(abs(sqrt(mesh.coords[1,j]^2+mesh.coords[2,j]^2+
                    mesh.coords[3,j]^2)-1.0)<=1e-9 for j in
           axes(mesh.coords,2))
    @test on>=MT.nnodes(mesh)-1
end

@testset "degenerate pole edges emit no Line cells" begin
    # Upstream marks an OCC `degenerate(0)` edge DONE at entry — the 4.15.2
    # oracle emits zero Line elements on a `Sphere` face's pole edges, so the
    # vertex's Point product is the only lower cell a collapsed row leaves.
    # The earlier `N=1` reading stored a single endpoint-edge discretization
    # whose (vertex,vertex) cell reached the merged public mesh as a
    # zero-length Line `validate` rejects; the collapsed row's incidence
    # inside the surface product is the welded `(uid,uid)` boundary pair
    # keyed by the Point product, never a curve cell.
    m=GeoModel()
    M._model_occ_add_sphere!(m,0.0,0.0,0.0,1.0;tag=1)
    M._model_occ_synchronize!(m)
    degs=sort!([c for (c,k) in m.curve_types if k===:degenerate])
    @test length(degs)==2
    sface=findfirst(==(:sphere),m.surface_types)
    mesh=M.mesh_model_surface(m,sface)
    @test MT.validate(mesh).ok && MT.ntris(mesh)>0
    # The collapsed edges grade to `:keep` — no stored discretization —
    # while the shared meridian seam keeps its ordinary grading.
    @test !any(deg->haskey(m.curve_params,deg),degs)
    seam=only([c for (c,k) in m.curve_types if k!==:degenerate])
    @test haskey(m.curve_params,seam)
    # The merged `.geo` route agrees on both algorithm arms: the public mesh
    # carries the meridian's Line cells and no zero-length self-loops.
    for tail in ("","""
Mesh.Algorithm = 7;
Field[1] = MathEvalAniso;
Field[1].M11 = 1;
Field[1].M22 = 1;
Field[1].M33 = 1;
Field[1].M12 = 0;
Field[1].M13 = 0;
Field[1].M23 = 0;
Background Field = 1;
""")
        src="SetFactory(\"OpenCASCADE\");\nSphere(1)={0,0,0,1};\n"*tail
        result=mktemp() do path,io
            write(io,src);close(io)
            Tessella.execute_geo(path;mesh_dim=2)
        end
        @test result.msg_error_count==0
        @test MT.validate(result.mesh).ok
        @test MT.ntris(result.mesh)>0
        @test size(result.mesh.segs,2)>0
        @test !any(c->result.mesh.segs[1,c]==result.mesh.segs[2,c],
                   axes(result.mesh.segs,2))
    end
end

@testset "initial-only boundary mesh stays a loud rejection" begin
    # `Mesh.Algorithm=3` suppresses the interior uv lattice; the
    # boundary-only PSLG is then just the two seam columns plus the two
    # collapsed pole rows, and the CDT triangulates it entirely from
    # same-height seam pairs — every cell weld-collapses. The kernel reports
    # the empty lift loudly rather than emitting a sliver "bead" sphere.
    m=GeoModel()
    M._model_occ_add_sphere!(m,0.0,0.0,0.0,1.0;tag=1)
    M._model_occ_synchronize!(m)
    sface=findfirst(==(:sphere),m.surface_types)
    M.set_algorithm!(m,2,sface,3)
    err=_error(()->M.mesh_model_surface(m,sface))
    @test err isa ArgumentError
    @test occursin("collapsed entirely",sprint(showerror,err))
end
