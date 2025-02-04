"""
    boyle_dykstra(x0::Vector{Float64}, projections::Vector{Function}; tol::Float64=1e-6, max_iter::Int=1000)

Implements the Boyle-Dykstra algorithm to project a point onto the intersection of multiple convex sets.

# Arguments
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
project_C1(x) = clamp.(x, 0.0, 1.0)    # Projection onto the box [0, 1]^n
project_C2(x) = x ./ sum(x)            # Projection onto the simplex

# Initial point
x0 = [1.5, 2.0, -0.5]

# Run the Boyle-Dykstra algorithm
x, num_iter, converged = boyle_dykstra(x0, [project_C1, project_C2])

println("Projection: ", x)
println("Converged: ", converged, " in ", num_iter, " iterations.")
```
"""
function boyle_dykstra(
  x0::Vector{Float64},
  projections::Vector{Function};
  tol::Float64 = 1e-6,
  max_iter::Int = 1000,
  kwargs...,
)
  m = length(projections) # Number of projection functions
  x = copy(x0) # Current point
  corrections = [zeros(length(x0)) for _ = 1:m] # Correction terms for each projection
  for iter = 1:max_iter
    x_old = copy(x)

    # Cycle through each projection function
    for k = 1:m
      # Apply projection with correction
      y_k = projections[k](x + corrections[k])

      # Update the correction term
      corrections[k] = x + corrections[k] - y_k

      # Update the current point
      x = y_k
    end

    # Check for convergence
    if norm(x - x_old) < tol
      return x
    end
  end

  return x
end
