"""
    numerical_projection!(sol::Vector{Float64},
                          x::Vector{Float64}, Px::AbstractMatrix,
                          Aeq::AbstractMatrix, beq::AbstractVector,
                          Ain::AbstractMatrix, bin::AbstractVector,
                          l::Vector{Float64}, u::Vector{Float64}) -> Vector{Float64}

Computes the numerical projection of a vector `x` onto a constrained set defined by
equality and inequality constraints, as well as variable bounds.

# Arguments
- `sol::Vector{Float64}`: The solution vector to be updated.
- `x::AbstractVector`: The input vector to be projected.
- `Px::AbstractMatrix`: Projection matrix.
- `Aeq::AbstractMatrix`: Coefficient matrix for equality constraints Aeq x = beq.
- `beq::AbstractVector`: Right-hand side of the equality constraints.
- `Ain::AbstractMatrix`: Coefficient matrix for inequality constraints Ain x ≤ bin.
- `bin::AbstractVector`: Right-hand side of the inequality constraints.
- `l::AbstractVector`: Lower bounds for the variables.
- `u::AbstractVector`: Upper bounds for the variables.

# Returns
- `AbstractVector`: The projection of `x` that satisfies the constraints.

# Example
```julia
using ADNLSModel, NLPModelsIpopt

# Example data
x   = [0.5, -1.0]
Px  = I
Aeq = [1.0 1.0]          # Equality constraint: x₁ + x₂ = 1
beq = [1.0]
Ain = [1.0 0.0; 0.0 1.0] # Inequality constraints: x₁ ≤ 0.8, x₂ ≤ 0.6
bin = [0.8, 0.6]
l = [-Inf, -Inf]         # No lower bounds
u = [Inf, Inf]           # No upper bounds

# Compute the projection
proj_x = numerical_projection(x, Px, Aeq, beq, Ain, bin, l, u)
println("Projected x: ", proj_x)
```
"""
function numerical_projection!(
  sol::S,
  x::S,
  Px = I,
  Aeq = spzeros(0, length(x)),
  beq = ones(0),
  Ain = spzeros(0, length(x)),
  bin = ones(0),
  l = -Inf * ones(length(x)),
  u = Inf * ones(length(x));
  atol = sqrt(eps(eltype(x))),
  rtol = sqrt(eps(eltype(x))),
  kwargs...,
) where {S}
  n, m = length(x), length(bin)
  lcon = vcat(beq, -Inf * ones(m))
  ucon = vcat(beq, bin)

  # Combine equality and inequality constraints into bounds
  lcon = vcat(beq, -Inf * ones(m))  # Lower bounds for constraints
  ucon = vcat(beq, bin)             # Upper bounds for constraints

  # Define the non-linear least squares model
  nls = ADNLSModel(
    d -> Px * x - d,                       # Objective: minimize distance to x
    x,                                     # Initial guess
    n,                                     # Number of variables
    l,                                     # Variable lower bounds
    u,                                     # Variable upper bounds
    vcat(Aeq, Ain),                        # Combined constraint matrix
    x -> eltype(S)[],                      # No equality constraints
    lcon,                                  # Combined lower constraint bounds
    ucon,                                  # Combined upper constraint bounds
  )

  # Solve the problem using IPOPT
  stats = percival(nls; atol = atol, rtol = rtol, kwargs...)

  # Return the projected solution
  sol .= stats.solution
  OK = if stats.primal_feas == Inf
    stats.status == :first_order
  else
    (norm(stats.primal_feas) < atol + norm(x) * rtol)
  end
  return OK
end
