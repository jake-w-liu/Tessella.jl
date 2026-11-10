module FEMLinearSolverTests
using Test, LinearAlgebra
using Tessella
const solve=Tessella.FEMLinearSolvers._fem_solve_spd
@testset "Sparse FEM independent manufactured solutions" begin
 for n in (1,7,100, Base.JLOptions().compile_enabled==3 ? 1000 : 10000)
  rows=Int[]; cols=Int[]; vals=Float64[]
  for i in 1:n
   push!(rows,i);push!(cols,i);push!(vals,2.0)
   if i<n
    append!(rows,(i,i+1));append!(cols,(i+1,i));append!(vals,(-0.5,-0.5))
   end
  end
  exact=[sin(i/19)+0.3 for i in 1:n]; b=zeros(n)
  for k in eachindex(vals);b[rows[k]]+=vals[k]*exact[cols[k]];end
  for scale in (1e-200,1.0,1e200)
   x=solve(n,rows,cols,vals.*scale,b.*scale)
   @test norm(x-exact)<=1e-9*norm(exact)
  end
 end
 @test solve(0,Int[],Int[],Float64[],Float64[])==[]
 @test solve(1,[1,1],[1,1],[1.,2.],[6.])≈[2.]
 @test_throws ArgumentError solve(2,[1,2,1],[1,2,2],[2.,2.,1.],[1.,1.])
 @test_throws ArgumentError solve(1,[1],[1],[-1.],[1.])
 @test_throws ArgumentError solve(1,[1],[1],[1.],[1.];maxiter=0)
 @test_throws ArgumentError solve(1,[1],[1],[1.],[Inf])
end

end
