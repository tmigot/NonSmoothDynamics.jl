"""
    boyle_dykstra!(sol::Vector{Float64}, x0::Vector{Float64}, projections::Vector{Function}; tol::Float64=1e-6, max_iter::Int=1000)

Implements the Boyle-Dykstra algorithm to project a point onto the intersection of multiple convex sets.

# Arguments
- `sol::Vector{Float64}`: The solution vector to be updated.
- `x0::Vector{Float64}`: The initial point to be projected.
- `projections::Vector{Function}`: A vector of projection functions, where each function computes the projection onto a specific convex set.
- `tol::Float64=1e-6`: Convergence tolerance. The algorithm stops if the change in `x` between iterations is less than this value.
- `max_iter::Int=1000`: Maximum number of iterations to run the algorithm.

# Returns
- `x::Vector{Float64}`: The projection of `x0` onto the intersection of the convex sets.
- `num_iter::Int`: The number of iterations performed.
- `converged::Bool`: Whether the algorithm converged within the maximum number of iterations.

# Example
```julia
# Define projection functions
project_C1(sol, x) = sol .= clamp.(x, 0.0, 1.0)    # Projection onto the box [0, 1]^n
project_C2(sol, x) = sol .= x ./ sum(x)            # Projection onto the simplex

# Initial point
x0 = [1.5, 2.0, -0.5]

# Run the Boyle-Dykstra algorithm
x, num_iter, converged = boyle_dykstra(x0, [project_C1, project_C2])

println("Projection: ", x)
println("Converged: ", converged, " in ", num_iter, " iterations.")
```
"""
function boyle_dykstra!(
  sol::S,
  x0::S,
  projections::Vector{Function};
  tol::Float64 = sqrt(eps(eltype(S))),
  max_iter::Int = 1000,
  verbose::Int = 0,
  kwargs...,
) where {S}
  T = eltype(S)
  m = length(projections) # Number of projection functions
  sol .= copy(x0) # Current point
  y_k = similar(x0)
  corrections = [zeros(length(x0)) for _ = 1:m] # Correction terms for each projection

  verbose > 0 && @info log_header([:iter, :nx, :proj_success], [Int, T, Int])

  for iter = 1:max_iter
    x_old = copy(sol)

    # Cycle through each projection function
    for k = 1:m
      # Apply projection with correction
      proj_success = projections[k](y_k, sol + corrections[k])

      if !proj_success
        @error "$k-th Projection failed"
        return false
      end

      # Update the correction term
      corrections[k] .= sol .+ corrections[k] .- y_k

      # Update the current point
      sol .= y_k
    end

    nx = norm(sol - x_old)
    verbose > 0 && @info log_row(Any[iter, nx, Int(proj_success)])
    # Check for convergence
    if nx < tol
      return true
    end
  end

  return false
end
