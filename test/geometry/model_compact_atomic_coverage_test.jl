module ModelCompactAtomicCoverageTests
using Test
const Tessella=getfield(parentmodule(@__MODULE__),:Tessella)
const M=Tessella.Model

@testset "Compact shifted Quad retains oriented atomic coverage" begin
    # Six producer vertices on the lattice 0 <= j <= i <= 2. The shifted
    # Quad joins atoms from two logical panels; these exact integer points
    # also make its area independently checkable.
    points=((0,0),(1,0),(1,1),(2,0),(2,1),(2,2))
    ids=UInt64[1,2,4,1,3,5,1,3,6]
    panels=M._model_transfinite_grid_panels(ids,(2,2),"compact test")
    atoms=M._model_transfinite_grid_compact_atoms(ids,
        Int32[1 2 2 5;2 3 6 6;5 6 5 9],panels,"compact test")
    @test atoms==NTuple{3,UInt64}[(1,2,3),(2,4,5),(2,5,3),(3,5,6)]
    cells=Tuple[(1,2,3),(2,4,5),(2,5,6,3)]
    area(cell)=abs(sum(points[cell[i]][1]*points[cell[mod1(i+1,length(cell))]][2]-
                       points[cell[i]][2]*points[cell[mod1(i+1,length(cell))]][1]
                       for i in eachindex(cell)))//2
    @test sum(area,cells)==2
    @test area(cells[3])==1
    @test_throws ArgumentError M._model_transfinite_grid_partition(panels,cells,"old panel test")
    @test M._model_transfinite_grid_atomic_partition(atoms,cells,"compact test")===nothing
    @test M._model_transfinite_grid_atomic_partition(atoms,reverse.(cells),"reversed sheet")===nothing
    @test M._model_transfinite_grid_partition(panels,atoms,"global evidence")!==nothing

    @test_throws ArgumentError M._model_transfinite_grid_atomic_partition(atoms,cells[1:2],"dropped cell")
    @test_throws ArgumentError M._model_transfinite_grid_atomic_partition(atoms,vcat(cells,[cells[1]]),"duplicate cell")
    @test_throws ArgumentError M._model_transfinite_grid_atomic_partition(atoms,[(1,4,6,3)],"cross-grid cell")
    @test_throws ArgumentError M._model_transfinite_grid_atomic_partition(atoms,[cells[1],cells[2],(2,6,5,3)],"bow-tie cell")
    @test_throws ArgumentError M._model_transfinite_grid_atomic_partition(atoms,[reverse(cells[1]),cells[2],cells[3]],"partial reversal")
    bad_atoms=copy(atoms);bad_atoms[3]=reverse(bad_atoms[3])
    @test_throws ArgumentError M._model_transfinite_grid_atomic_partition(bad_atoms,cells,"inconsistent constructor")
    @test_throws ArgumentError M._model_transfinite_grid_partition(vcat(panels,panels),vcat(atoms,atoms),"overlapping Surface owners")
    @test_throws ArgumentError M._model_transfinite_grid_compact_atoms(ids,reshape(Int32[1,2,99],3,1),panels,"absent row")
end
end
