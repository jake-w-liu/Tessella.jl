"""
    validation/common.jl

Shared helpers for the Tessella-vs-external-tools validation suite. The convention
(see README.md): for each geometry, mesh it with BOTH Tessella and a reference tool
(gmsh) on the same domain, then cross-check element quality, the meshed volume
against the analytic value, and wall-clock time. The reference `.geo` scripts are
retained in each case folder.
"""

using Tessella
using Tessella.MeshTypes
using Tessella.IO: read_msh
using Printf

gmsh_available() = try success(pipeline(`gmsh --version`; stdout=devnull, stderr=devnull)); catch; false; end

"gmsh version string (or \"n/a\")."
function gmsh_version()
    try
        return strip(read(pipeline(`gmsh --version`; stderr=`cat`), String))
    catch
        return "n/a"
    end
end

"""
    run_gmsh(geo, out; algo=1) -> (ok::Bool, seconds::Float64)

Mesh `geo` to `out` (MSH 2.2) with gmsh's 3-D mesher, timing the wall clock.
`algo` selects `Mesh.Algorithm3D` (1 = Delaunay, 10 = HXT).
"""
function run_gmsh(geo::AbstractString, out::AbstractString; algo::Integer=1)
    isfile(out) && rm(out)
    cmd = `gmsh -3 $geo -o $out -format msh2 -algo del3d -clscale 1.0 -v 2 -nt 1`
    ok = false
    seconds = @elapsed (ok = success(pipeline(cmd; stdout=devnull, stderr=devnull)))
    return (ok = ok && isfile(out), seconds = seconds)
end

"Quality/size metrics of a linear tet [`Mesh`](@ref)."
function tet_metrics(m::Mesh)
    nt = ntets(m)
    vol = 0.0
    @inbounds for t in 1:nt
        vol += tet_volume(node(m,m.tets[1,t]),node(m,m.tets[2,t]),node(m,m.tets[3,t]),node(m,m.tets[4,t]))
    end
    q = mesh_quality(m)
    return (nnodes = nnodes(m), ntets = nt, volume = vol,
            min_dih = q.min_dihedral_deg, mean_dih = q.mean_dihedral_deg,
            slivers = q.n_slivers)
end

"Read a gmsh .msh and return its tet metrics (tets only)."
gmsh_metrics(msh::AbstractString) = tet_metrics(read_msh(msh).mesh)

"Mesh a Tessella surface into a volume, timed; returns (metrics, seconds)."
function tessella_fill(surface::Mesh; kwargs...)
    m = mesh_volume(surface; kwargs...)     # warm up / correctness path
    seconds = @elapsed (m = mesh_volume(surface; kwargs...))
    return (tet_metrics(m), seconds)
end

"Pretty relative error vs an analytic reference (percent)."
relerr(v, ref) = ref == 0 ? abs(v) : abs(v - ref) / abs(ref)
pct(x) = @sprintf("%.4g%%", 100x)
