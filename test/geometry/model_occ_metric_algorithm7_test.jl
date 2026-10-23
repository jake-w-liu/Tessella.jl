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
