# Shared-edge discretization stability for OCC volume shells.
#
# A native `Mesh.Algorithm=3` plane face inside an OCC shell consumes its
# boundary Curve products through a detached dim-1 generation. That plan used
# to regrade every Curve from an `options.ctx_lc` derived from a bounding box
# that samples the stored `curve_params` nodes themselves: the first surface
# mesh populated `curve_params`, lifting the sector's degenerate planar bbox
# off its padding branch, so every later face regraded a denser — then
# ulp-drifting — subdivision. Sibling faces welded exact coordinates over
# different parameter lists and `_model_oriented_shell_volume` reported
# incidence-1 edges.
#
# Upstream `meshGFace` reuses each `GEdge`'s stored `mesh_vertices` verbatim;
# the detached plan must emit stored `curve_params`, grading only curves that
# still lack a discretization.

using Test
using Tessella
const M=Tessella.Model
const MT=Tessella.MeshTypes

function _sector_model(algorithm)
    m=GeoModel()
    M._model_occ_add_cylinder!(m,0.0,0.0,0.0,0.0,0.0,2.0,1.0;tag=1,angle=pi)
    M._model_occ_synchronize!(m)
    m.meshing.algorithm_2d=algorithm
    return m
end

_coord_bits(mesh,j)=ntuple(i->reinterpret(UInt64,mesh.coords[i,j]),3)

function _coord_set(mesh)
    out=Set{NTuple{3,UInt64}}()
    for j in axes(mesh.coords,2)
        push!(out,_coord_bits(mesh,j))
    end
    return out
end

@testset "OCC sector: stored curve parameters do not ratchet across faces" begin
    m=_sector_model(3)
    # The first plane face lazily grades and stores every boundary curve.
    M.mesh_model_surface(m,2)
    stored=Dict(k=>copy(v) for (k,v) in m.curve_params)
    @test !isempty(stored)
    for s in (3,4,5)
        M.mesh_model_surface(m,s)
        @test m.curve_params==stored
    end
end

@testset "OCC sector: native plane boundary shares curved-kernel nodes bitwise" begin
    m=_sector_model(3)
    curved=M._model_occ_curved_surface_mesh(m,1,"test";initial_only=true)
    cap=M.mesh_model_surface(m,2)
    shared=intersect(_coord_set(curved),_coord_set(cap))
    # Six arc parameters bound both faces at the top rim.
    @test length(shared)==length(m.curve_params[1])
    @test length(shared)>=6
end

@testset "OCC sector: volume shell welds exactly under Mesh.Algorithm=3" begin
    m=_sector_model(3)
    mesh=M.mesh_model_volume(m,1)
    @test MT.nnodes(mesh)>0
end
