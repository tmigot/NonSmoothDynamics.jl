"""
    numerical_projection(x::Vector{Float64}, Aeq::AbstractMatrix, beq::AbstractVector,
                         Ain::AbstractMatrix, bin::AbstractVector,
                         l::Vector{Float64}, u::Vector{Float64}) -> Vector{Float64}

Computes the numerical projection of a vector `x` onto a constrained set defined by
equality and inequality constraints, as well as variable bounds.

# Arguments
- `x::Vector{Float64}`: The input vector to be projected.
- `Aeq::AbstractMatrix`: Coefficient matrix for equality constraints \( A_{eq} x = b_{eq} \).
- `beq::AbstractVector`: Right-hand side of the equality constraints.
- `Ain::AbstractMatrix`: Coefficient matrix for inequality constraints \( A_{in} x \leq b_{in} \).
- `bin::AbstractVector`: Right-hand side of the inequality constraints.
- `l::Vector{Float64}`: Lower bounds for the variables.
- `u::Vector{Float64}`: Upper bounds for the variables.

# Returns
- `Vector{Float64}`: The projection of `x` that satisfies the constraints.

# Example
```julia
using ADNLSModel, NLPModelsIpopt

# Example data
x = [0.5, -1.0]
Aeq = [1.0 1.0]          # Equality constraint: x₁ + x₂ = 1
beq = [1.0]
Ain = [1.0 0.0; 0.0 1.0] # Inequality constraints: x₁ ≤ 0.8, x₂ ≤ 0.6
bin = [0.8, 0.6]
l = [-Inf, -Inf]         # No lower bounds
u = [Inf, Inf]           # No upper bounds

# Compute the projection
proj_x = numerical_projection(x, Aeq, beq, Ain, bin, l, u)
println("Projected x: ", proj_x)
```
"""
function numerical_projection(
  x::Vector{Float64},
  Aeq,
  beq,
  Ain,
  bin,
  l::Vector{Float64},
  u::Vector{Float64},
)
  n, m = length(x), length(bin)
  lcon = vcat(beq, -Inf * ones(m))
  ucon = vcat(beq, bin)

  # Combine equality and inequality constraints into bounds
  lcon = vcat(beq, -Inf * ones(m))  # Lower bounds for constraints
  ucon = vcat(beq, bin)             # Upper bounds for constraints

  # Define the non-linear least squares model
  nls = ADNLSModel(
    x,                                      # Initial guess
    d -> norm(x - d)^2,                    # Objective: minimize distance to x
    n,                                     # Number of variables
    vcat(Aeq, Ain),                        # Combined constraint matrix
    lcon,                                  # Combined lower constraint bounds
    ucon,                                  # Combined upper constraint bounds
    l,                                     # Variable lower bounds
    u,                                      # Variable upper bounds
  )

  # Solve the problem using IPOPT
  stats = ipopt(nls)

  # Return the projected solution
  return nls.solution
end
