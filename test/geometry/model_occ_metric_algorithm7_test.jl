# `Mesh.Algorithm=7` on OCC-curved faces — the certified parametric metric
# producer.
#
# Upstream `meshGFaceBamg` refines a curved GFace in its own (u,v) domain under
# the background metric's tangent-plane pullback. The native producer mirrors
# that contract on the OCC pcurve sheet: the same constrained-Delaunay initial
# triangulation the Delaunay-family kernel lifts, refined by the shared
# certified core under J'MJ evaluated through `_occ_surface_d1`, and folded
# back onto the closed face by the bitwise weld the lifted path already
# certifies. A scalar or absent field has no surface-parameter meaning on a
# non-affine chart, so the selector keeps a loud requirement rejection instead
# of falling back to the planar isotropic context.

using Test
using Tessella
const M=Tessella.Model
const MT=Tessella.MeshTypes
const SF=Tessella.SizeField

_error(f)=try;f();nothing catch err;err end

function _geo_exec(source::AbstractString;mesh_dim=2)
    return mktemp() do path,io
        write(io,source);close(io)
        Tessella.execute_geo(path;mesh_dim=mesh_dim)
    end
end

# A full 2π OCC cylinder synchronized into explicit shell topology: surface 1
# is the closed :cylinder lateral — the periodic-seam weld case.
function _full_lateral()
    m=GeoModel()
    M._model_occ_add_cylinder!(m,0.0,0.0,0.0,0.0,0.0,2.0,1.5;tag=1)
    M._model_occ_synchronize!(m)
    return m,findfirst(==(:cylinder),m.surface_types)
end

# An explicit π-sector shell: the lateral is an open :cylinder face — the
# ordinary boundary-only case.
function _sector_lateral()
    m=GeoModel()
    M._model_occ_add_cylinder!(m,0.0,0.0,0.0,0.0,0.0,2.0,1.5;tag=1,angle=pi)
    M._model_occ_synchronize!(m)
    return m,findfirst(==(:cylinder),m.surface_types)
end

# A cylinder cut by an axis-interior box: the surviving lateral band is a
# *trimmed* :cylinder face — its unwrapped sheet covers only part of the
# analytic (u,v) support rectangle, the case the boundary-chain domain
# certificate exists for.
function _trimmed_lateral()
    m=GeoModel()
    M.add_cylinder!(m,0.0,0.0,-1.0,0.0,0.0,4.0,1.5;tag=1)
    M.add_box!(m,-5.0,-5.0,0.0,10.0,10.0,4.5;tag=2)
    M.boolean_volumes_multi!(m,:difference,[1],[2])
    lateral=first(t for t in sort!(collect(keys(m.surface_geometry))) if
                  get(m.surface_geometry[t],:occ,nothing)===:cylinder)
    m.meshing.algorithm_2d=7
    return m,lateral
end

function _edge_incidence(mesh)
    incidence=Dict{NTuple{2,Int32},Int}()
    for c in axes(mesh.tris,2),(a,b) in ((1,2),(2,3),(3,1))
        key=minmax(mesh.tris[a,c],mesh.tris[b,c])
        incidence[key]=get(incidence,key,0)+1
    end
    return incidence
end
# Coordinate-pair boundary chain — the bitwise weld certificate sibling faces
# share through their common Curve products.
function _boundary_coords(mesh)
    incidence=_edge_incidence(mesh)
    return Set(minmax(Tuple(mesh.coords[:,a]),Tuple(mesh.coords[:,b]))
               for (a,b) in keys(incidence) if incidence[(a,b)]==1)
end

function _assert_metric_cert(mesh,field,tag)
    bounded=SF.BoundedSize(field;size_min=0.0,size_max=SF.GMSH_MAX_SIZE,
        factor=1.0)
    return all(keys(_edge_incidence(mesh))) do (a,b)
        SF.metric_edge_length(bounded,Tuple(mesh.coords[:,a]),
            Tuple(mesh.coords[:,b]);entity=(2,tag))<=1.0+1e-12
    end
end

@testset "closed OCC cylinder lateral: algorithm 7 welds the seam" begin
    m,lateral=_full_lateral()
    m.meshing.algorithm_2d=7
    field=MathEvalAnisoField(m11="1",m22="1",m33="1")
    production=M.mesh_model_surface(m,lateral;size_field=field,
        _return_surface_production=true)
    @test production isa M._ModelSurfaceProduction
    mesh=production.mesh
    @test MT.validate(mesh).ok
    @test MT.ntris(mesh)>0 && MT.nnodes(mesh)>0
    # Every vertex lifts onto the analytic lateral — the weld folded the two
    # seam sheets onto shared bitwise rows and the split/interior vertices
    # evaluate `_occ_surface_point`, so nothing sits off the cylinder.
    @test all(k->abs(mesh.coords[1,k]^2+mesh.coords[2,k]^2-1.5^2)<1e-9,
              1:MT.nnodes(mesh))
    # Face-owned rows carry the (u,v) provenance their coordinates evaluate.
    sg=m.surface_geometry[lateral]
    @test all(eachindex(production.node_tags)) do i
        production.node_owners[i][1]!=2 ||
            (p=production.node_parameters[i]; p!==nothing &&
             M._occ_surface_point(sg,p[1],p[2])==Tuple(mesh.coords[:,i]))
    end
    # Retained rows are exactly the lower-dimensional product identities;
    # births remain a complete output permutation.
    @test count(production.retained_nodes)>0
    @test all(eachindex(production.node_tags)) do i
        production.retained_nodes[i]==(production.node_owners[i][1]<2)
    end
    @test sort(production.birth_rows)==collect(1:length(production.node_tags))
    # The unit-metric certificate holds on lifted world chords.
    @test _assert_metric_cert(mesh,field,lateral)
    # Boundary ownership weld: every incidence-1 perimeter vertex is a
    # retained lower-dimensional product — the curve cells the metric path
    # may not split or move — and a same-field rerun reproduces the exact
    # coordinate chain (deterministic bitwise weld).
    incidence=_edge_incidence(mesh)
    perimeter=[edge for (edge,hits) in incidence if hits==1]
    @test !isempty(perimeter)
    @test all(edge->production.retained_nodes[edge[1]] &&
                  production.retained_nodes[edge[2]],perimeter)
    m2,_=_full_lateral()
    m2.meshing.algorithm_2d=7
    rerun=M.mesh_model_surface(m2,lateral;size_field=field)
    @test _boundary_coords(rerun)==_boundary_coords(mesh)
end

@testset "open OCC sector lateral: algorithm 7 surface + volume" begin
    m,lateral=_sector_lateral()
    m.meshing.algorithm_2d=7
    field=MathEvalAnisoField(m11="1",m22="1",m33="1")
    mesh=M.mesh_model_surface(m,lateral;size_field=field)
    @test MT.validate(mesh).ok && MT.ntris(mesh)>0
    @test all(k->abs(mesh.coords[1,k]^2+mesh.coords[2,k]^2-1.5^2)<1e-9,
              1:MT.nnodes(mesh))
    @test _assert_metric_cert(mesh,field,lateral)
    # Per-face selection inside a volume: the caps keep their own algorithm
    # while the lateral refines under the metric.
    mv,lat=_sector_lateral()
    mv.meshing.algorithm_2d=6
    M.set_algorithm!(mv,2,lat,7)
    volume=M.mesh_model_volume(mv,1;size_field=field)
    @test MT.validate(volume).ok && MT.ntets(volume)>0
end

@testset "algorithm 7 keeps the field contract loud on curved faces" begin
    # No field, a scalar field and the ordinary scalar context all reject —
    # the producer requires a tensor background it can pull back through J'MJ.
    for field in (nothing,MathEvalField("0.5"))
        m,lateral=_full_lateral()
        m.meshing.algorithm_2d=7
        err=_error(()->M.mesh_model_surface(m,lateral;
            size_field=field))
        @test err isa ArgumentError
        @test occursin("anisotropic",sprint(showerror,err))
        # The rejection is transactional: the ordinary Delaunay-family path
        # still meshes the same face afterwards.
        m.meshing.algorithm_2d=5
        @test MT.validate(M.mesh_model_surface(m,lateral)).ok
    end
    # Global algorithm 7 through the volume path likewise names the
    # requirement on the curved face.
    m,lateral=_sector_lateral()
    m.meshing.algorithm_2d=7
    err=_error(()->M.mesh_model_volume(m,1))
    @test err isa ArgumentError
    @test occursin("anisotropic",sprint(showerror,err))
end

@testset "trimmed OCC cylinder face: algorithm 7 covers the cut domain" begin
    m,lateral=_trimmed_lateral()
    field=MathEvalAnisoField(m11="1",m22="1",m33="1")
    production=M.mesh_model_surface(m,lateral;size_field=field,
        _return_surface_production=true)
    @test production isa M._ModelSurfaceProduction
    mesh=production.mesh
    @test MT.validate(mesh).ok && MT.ntris(mesh)>0
    # The surviving band is z∈[-1,0] — every vertex lifts onto the analytic
    # lateral, and no vertex crosses the trimming box's section plane.
    @test all(k->abs(mesh.coords[1,k]^2+mesh.coords[2,k]^2-1.5^2)<1e-9,
              1:MT.nnodes(mesh))
    @test all(k->-1.0-1e-12<=mesh.coords[3,k]<=1e-12,1:MT.nnodes(mesh))
    # The face-owned rows evaluate their own (u,v) provenance.
    sg=m.surface_geometry[lateral]
    @test all(eachindex(production.node_tags)) do i
        production.node_owners[i][1]!=2 ||
            (p=production.node_parameters[i]; p!==nothing &&
             M._occ_surface_point(sg,p[1],p[2])==Tuple(mesh.coords[:,i]))
    end
    # The perimeter is exactly the two boundary circles — the section curve
    # at z=0 and the rim at z=-1 — every edge a retained lower product.
    incidence=_edge_incidence(mesh)
    perimeter=[edge for (edge,hits) in incidence if hits==1]
    @test !isempty(perimeter)
    @test all(edge->production.retained_nodes[edge[1]] &&
                  production.retained_nodes[edge[2]],perimeter)
    @test sort!(collect(Set(round(mesh.coords[3,e[1]];digits=9)
                            for e in perimeter)))==[-1.0,0.0]
    @test _assert_metric_cert(mesh,field,lateral)
    # Deterministic bitwise weld on a same-field rerun.
    m2,lateral2=_trimmed_lateral()
    rerun=M.mesh_model_surface(m2,lateral2;size_field=field)
    @test _boundary_coords(rerun)==_boundary_coords(mesh)
end

@testset "OCC plane face with ellipse boundary: algorithm 7" begin
    # `OCC.addDisk` binds a `:plane` face whose single boundary edge is an
    # ellipse — the upstream metric selector reaches curved boundaries on
    # plane faces through exactly this path.
    m=GeoModel()
    disk=M._model_occ_add_disk!(m,0.0,0.0,0.0,1.0,1.0)
    M._model_occ_synchronize!(m)
    m.meshing.algorithm_2d=7
    field=MathEvalAnisoField(m11="1",m22="1",m33="1")
    mesh=M.mesh_model_surface(m,disk;size_field=field)
    @test MT.validate(mesh).ok && MT.ntris(mesh)>0
    @test all(k->abs(mesh.coords[3,k])<1e-12,1:MT.nnodes(mesh))
    @test all(k->mesh.coords[1,k]^2+mesh.coords[2,k]^2<=1.0+1e-9,
              1:MT.nnodes(mesh))
    @test _assert_metric_cert(mesh,field,disk)
end

@testset "builtin circle-bounded plane: algorithm 7" begin
    m=GeoModel()
    M.add_point!(m,1.0,0.0,0.0;tag=1)
    M.add_point!(m,-1.0,0.0,0.0;tag=2)
    M.add_point!(m,0.0,0.0,0.0;tag=3)
    M.add_circle_arc!(m,1,3,2;tag=1)
    M.add_circle_arc!(m,2,3,1;tag=2)
    loop=M.add_curve_loop!(m,[1,2];tag=1)
    disk=M.add_plane_surface!(m,[loop];tag=1)
    m.meshing.algorithm_2d=7
    field=MathEvalAnisoField(m11="1",m22="1",m33="1")
    mesh=M.mesh_model_surface(m,disk;size_field=field)
    @test MT.validate(mesh).ok && MT.ntris(mesh)>0
    @test all(k->abs(mesh.coords[3,k])<1e-12,1:MT.nnodes(mesh))
    @test all(k->mesh.coords[1,k]^2+mesh.coords[2,k]^2<=1.0+1e-9,
              1:MT.nnodes(mesh))
    @test _assert_metric_cert(mesh,field,disk)
end

@testset "lofted tensor face: algorithm 7 refines the pcurve domain" begin
    m=GeoModel()
    c1=M._model_occ_add_circle!(m,0.0,0.0,0.0,1.0)
    c2=M._model_occ_add_circle!(m,0.0,0.0,2.0,0.5)
    w1=M._model_occ_add_wire!(m,[c1])
    w2=M._model_occ_add_wire!(m,[c2])
    M._model_occ_thru_sections!(m,-1,[w1,w2];make_solid=false)
    M._model_occ_synchronize!(m)
    lateral=only([t for t in keys(m.surface_types)
                  if m.surface_types[t]===:tensor_bspline])
    m.meshing.algorithm_2d=7
    field=MathEvalAnisoField(m11="1",m22="1",m33="1")
    mesh=M.mesh_model_surface(m,lateral;size_field=field)
    @test MT.validate(mesh).ok && MT.ntris(mesh)>0
    # The loft is the conical frustum r(z)=1-z/4 between the two wires.
    @test all(k->abs(sqrt(mesh.coords[1,k]^2+mesh.coords[2,k]^2)-
                     (1.0-0.25*mesh.coords[3,k]))<1e-7,1:MT.nnodes(mesh))
    @test _assert_metric_cert(mesh,field,lateral)
end

@testset "scalar pass then algorithm 7: closed Curve params stay complete" begin
    # Regression for the stored-parameter contract on closed curves:
    # `mesh_curve(closed=true)` returns the nedge station list — the
    # closing t1 is implicit — while every `curve_params` consumer holds
    # the open-style `[t0,…,t1]` contract. A first scalar surface pass
    # used to store the truncated list verbatim; the algorithm-7 reuse
    # path then re-emitted one Curve cell short of the owned lower
    # product and the incidence certificate rejected. The grading layer
    # now pins the endpoint before storage.
    m=GeoModel()
    M.add_cylinder!(m,0.0,0.0,-1.0,0.0,0.0,4.0,1.5;tag=1)
    M.add_box!(m,-5.0,-5.0,0.0,10.0,10.0,4.5;tag=2)
    M.boolean_volumes_multi!(m,:difference,[1],[2])
    lateral=first(t for t in sort!(collect(keys(m.surface_geometry))) if
                  get(m.surface_geometry[t],:occ,nothing)===:cylinder)
    m.meshing.algorithm_2d=5
    mesh5=M.mesh_model_surface(m,lateral)
    @test MT.validate(mesh5).ok && MT.ntris(mesh5)>0
    # Every stored closed-curve list carries both pinned endpoints.
    stored_closed=[(c,p) for (c,p) in m.curve_params
                   if get(m.curves,c,(0,1))[1]==get(m.curves,c,(0,1))[2]]
    @test !isempty(stored_closed)
    for (c,p) in stored_closed
        t0,t1=M._model_curve_param_bounds(m,c,"test")
        @test p[1]==t0 && p[end]==t1
    end
    # The metric reuse path re-emits the same closed discretization —
    # no incidence loss against the owned lower products.
    m.meshing.algorithm_2d=7
    field=MathEvalAnisoField(m11="1",m22="1",m33="1")
    mesh7=M.mesh_model_surface(m,lateral;size_field=field)
    @test MT.validate(mesh7).ok && MT.ntris(mesh7)>0
    @test all(k->abs(mesh7.coords[1,k]^2+mesh7.coords[2,k]^2-1.5^2)<1e-9,
              1:MT.nnodes(mesh7))
    @test _assert_metric_cert(mesh7,field,lateral)
end

@testset ".geo higher route: anisotropic fields reach OCC arc grading" begin
    # `generate(dim>1)` emits the boundary Curve products before the Surface
    # phase (`_boundary_regrade=false`), so each boundary curve must already
    # resolve the algorithm-7 certificate field at grade time. The OCC arc
    # grader used to size circles with the scalarized `speed/h` law even
    # while its incident face certified under the tensor — under an M=16
    # transverse metric the rim circles came back under-resolved and the
    # retained-boundary certificate rejected the whole assembly.
    src="""
SetFactory("OpenCASCADE");
Cylinder(1) = {0,0,0, 0,0,2, 1.5};
Mesh.Algorithm = 7;
Field[1] = MathEvalAniso;
Field[1].M11 = 16;
Field[1].M22 = 16;
Field[1].M33 = 1;
Field[1].M12 = 0;
Field[1].M13 = 0;
Field[1].M23 = 0;
Background Field = 1;
"""
    result=_geo_exec(src)
    @test result.msg_error_count==0
    mesh=result.mesh
    @test MT.validate(mesh).ok && MT.ntris(mesh)>0
    @test all(k->(abs(hypot(mesh.coords[1,k],mesh.coords[2,k])-1.5)<1e-9 ||
                  abs(mesh.coords[3,k])<1e-9 || abs(mesh.coords[3,k]-2)<1e-9),
              1:MT.nnodes(mesh))
end
