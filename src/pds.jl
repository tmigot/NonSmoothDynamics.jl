"""
    projected_dynamical_system(x0::Vector{Float64}, F::Function, project_C::Function;
                               step_size::Float64=0.01, max_iter::Int=1000,
                               atol::Float64=1e-6, rtol::Float64=1e-6)

Simulates the discretized dynamics of a Projected Dynamical System (PDS):
    x(t+1) = P_{C}(x(t) - step_size * F(x(t)))

# Arguments
- `x0::Vector{Float64}`: Initial state (starting point within the set C).
- `F::Function`: The vector field defining the dynamics (e.g., gradient, payoff vector, etc.).
- `project_C::Function`: A function to compute the projection onto the feasible set C.
- `step_size::Float64=0.01`: The time step for discretization.
- `max_iter::Int=1000`: Maximum number of iterations.
- `atol::Float64=1e-6`: Convergence tolerance. The algorithm stops if the norm of the change is less than this value.
- `rtol::Float64=1e-6`: Convergence tolerance.

# Returns
- `x_vals::Matrix{Float64}`: Matrix where each column is the state of the system at a time step.
- `num_iter::Int`: Number of iterations performed.
- `converged::Bool`: Whether the system reached convergence.

# Example
```julia
# Define the vector field F(x) = x - 1 (gradient of f(x) = 0.5 * ||x - 1||^2)
F(x) = x - 1.0

# Define projection onto the box [0, 1]^2
project_C(x) = clamp.(x, 0.0, 1.0)

# Initial point
x0 = [2.0, -1.0]

# Simulate the PDS
x_vals, num_iter, converged = projected_dynamical_system(x0, F, project_C)

println("Final State: ", x_vals[:, end])
println("Converged: ", converged, " in ", num_iter, " iterations.")
```
"""
function projected_dynamical_system(
  x0::Vector{Float64},
  F::Function,
  project_C::Function;
  step_size::Float64 = 0.01,
  max_iter::Int = 1000,
  atol::Float64 = 1e-6,
  rtol::Float64 = 1e-6,
)
  # Initialize variables
  x = copy(x0)
  x_vals = [x] # Store the trajectory
  converged = false
  iter = 0
  for iter = 1:max_iter
    # Compute the update step:
    x_new = project_C(
      x - step_size * F(x);
      step_size = step_size,
      x = x,
      atol = atol,
      rtol = rtol,
    )
    push!(x_vals, x_new)

    # Check for convergence
    if norm(x_new - x) < atol
      converged = true
      break
    end

    # Update for the next iteration
    x = x_new
  end

  # Return trajectory, iterations, and convergence status
  return hcat(x_vals...), iter, converged
end
