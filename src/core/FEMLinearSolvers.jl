"""Internal sparse finite-element linear solvers, with no external solver runtime."""
module FEMLinearSolvers

using LinearAlgebra: norm, dot

# Assemble duplicate COO entries deterministically, then equilibrate symmetrically.
# Dirichlet unknowns must already have been eliminated by the caller.
function _fem_solve_spd(n::Int, rows::Vector{Int}, columns::Vector{Int},
        values::Vector{Float64}, rhs::Vector{Float64};
        rtol::Real=1e-11, maxiter::Int=max(200, 4n))
    n >= 0 || throw(ArgumentError("FEM system size must be nonnegative"))
    length(rhs) == n || throw(DimensionMismatch("FEM right-hand side"))
    length(rows) == length(columns) == length(values) ||
        throw(DimensionMismatch("FEM sparse triplets"))
    isfinite(rtol) && 0 < rtol < 1 || throw(ArgumentError("invalid FEM tolerance"))
    maxiter >= 0 || throw(ArgumentError("negative FEM iteration limit"))
    all(isfinite, rhs) || throw(ArgumentError("nonfinite FEM right-hand side"))
    entries = Dict{Tuple{Int,Int},Float64}()
    for k in eachindex(values)
        i, j, v = rows[k], columns[k], values[k]
        1 <= i <= n && 1 <= j <= n || throw(ArgumentError("FEM sparse index out of range"))
        isfinite(v) || throw(ArgumentError("nonfinite FEM matrix entry"))
        key = (i,j)
        total = get(entries, key, 0.0) + v
        isfinite(total) || throw(ArgumentError("FEM matrix accumulation overflow"))
        entries[key] = total
    end
    n == 0 && return Float64[]
    diagonal = [get(entries, (i,i), 0.0) for i in 1:n]
    all(v -> isfinite(v) && v > 0, diagonal) ||
        throw(ArgumentError("FEM matrix needs positive diagonal entries"))
    # An asymmetric assembly is an error, rather than input for a different solver.
    for ((i,j),v) in entries
        w = get(entries, (j,i), 0.0)
        abs(v-w) <= 32eps(Float64)*max(abs(v),abs(w),floatmin(Float64)) ||
            throw(ArgumentError("FEM matrix is not symmetric"))
    end
    roots = sqrt.(diagonal)
    keys_sorted = sort!(collect(keys(entries)))
    rowptr = ones(Int,n+1)
    for (i,_) in keys_sorted
        rowptr[i+1] += 1
    end
    for i in 2:n+1
        rowptr[i] += rowptr[i-1]-1
    end
    cols = Vector{Int}(undef,length(keys_sorted))
    coefficients = Vector{Float64}(undef,length(keys_sorted))
    for (k,(i,j)) in enumerate(keys_sorted)
        cols[k] = j
        coefficients[k] = i == j ? 1.0 : (entries[(i,j)]/roots[i])/roots[j]
        isfinite(coefficients[k]) || throw(ArgumentError("FEM equilibration overflow"))
    end
    b = rhs ./ roots
    all(isfinite,b) || throw(ArgumentError("FEM right-hand-side equilibration overflow"))
    scale = maximum(abs,b)
    scale == 0 && return zeros(n)
    b ./= scale
    target = Float64(rtol)*norm(b)
    x = zeros(n)
    residual = copy(b)
    direction = copy(b)
    product = zeros(n)
    rr = dot(residual,residual)
    function multiply!(out, input)
        for i in 1:n
            s = 0.0
            for k in rowptr[i]:rowptr[i+1]-1
                s = muladd(coefficients[k],input[cols[k]],s)
            end
            out[i] = s
        end
    end
    converged = false
    for iteration in 1:maxiter
        multiply!(product,direction)
        curvature = dot(direction,product)
        isfinite(curvature) && curvature > 0 ||
            throw(ArgumentError("FEM matrix is not positive definite"))
        alpha = rr/curvature
        for i in 1:n
            x[i] += alpha*direction[i]
            residual[i] -= alpha*product[i]
        end
        next_rr = dot(residual,residual)
        if sqrt(next_rr) <= target
            # Recurrence residuals can drift; acceptance uses the actual operator.
            multiply!(product,x)
            for i in 1:n
                residual[i] = b[i]-product[i]
            end
            next_rr = dot(residual,residual)
            if sqrt(next_rr) <= target
                converged = true
                break
            end
            copyto!(direction,residual)
        else
            beta = next_rr/rr
            for i in 1:n
                direction[i] = residual[i]+beta*direction[i]
            end
        end
        isfinite(next_rr) || throw(ArgumentError("nonfinite FEM solver residual"))
        rr = next_rr
    end
    converged || throw(ArgumentError("FEM conjugate-gradient solve did not converge"))
    for i in 1:n
        x[i] = (x[i]*scale)/roots[i]
    end
    all(isfinite,x) || throw(ArgumentError("FEM solution is not representable"))
    return x
end

end # module
