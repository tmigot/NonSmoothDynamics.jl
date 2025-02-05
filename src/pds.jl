@doc (@doc projected_dynamical_system!) function projected_dynamical_system(
  x0::Vector{Float64},
  F::Function,
  project_fun!::Function,
  args...;
  kwargs...,
)
  vals = InterpolationStruct(length(x0), args...)
  converged = projected_dynamical_system!(vals, x0, F, project_fun!; kwargs...)
  return vals.x_vals, vals.t_vals, converged
end

"""
    projected_dynamical_system!(vals::InterpolationStruct,
                                x0::Vector{Float64}, F::Function, project_fun!::Function;
                                kwargs...)
    projected_dynamical_system(x0::Vector{Float64}, F::Function, project_fun!::Function,
                               args...; kwargs...)

The `args...` are passed to InterpolationStruct constructor.

Simulates the discretized dynamics of a Projected Dynamical System (PDS):
    x(t+1) = P_{C}(x(t) - h * F(x(t)))

# Arguments
- `vals::InterpolationStruct`: Interpolation struct containing the time and state values.
- `x0::Vector{Float64}`: Initial state (starting point within the set C).
- `F::Function`: The vector field defining the dynamics (e.g., gradient, payoff vector, etc.).
- `project_fun!::Function`: A function to compute the projection onto the feasible set C.
- `verbose::Int=0`: Verbosity level.
- `proj_verbose::Int=0`: Verbosity level for the projection step.
- `project_x0::Bool = true`:
- `atol::Float64=1e-6`: Convergence tolerance. The algorithm stops if the norm of the change is less than this value.
- `rtol::Float64=1e-6`: Convergence tolerance.
- `stop_first_stable::Bool=true`: Stop the simulation if the first stable point is found.

# Returns
- `stats::Bool`: Whether the system computation was successful.

`vals.x_vals::Matrix{Float64}`: Matrix where each column is the state of the system at a time step.

# Example
```julia
# Define the vector field F(x) = x - 1 (gradient of f(x) = 0.5 * ||x - 1||^2)
F(x) = x - 1.0

# Define projection onto the box [0, 1]^2
project_fun!(x) = clamp.(x, 0.0, 1.0)

# Initial point
x0 = [2.0, -1.0]

# Simulate the PDS
x_vals, t_vals, converged = projected_dynamical_system(x0, F, project_fun!)
```
"""
function projected_dynamical_system!(
  vals::InterpolationStruct{T1,T2},
  x0::Vector{Float64},
  F::Function,
  project_fun!::Function;
  verbose::Int = 0,
  proj_verbose::Int = 0,
  project_x0::Bool = true,
  atol::Float64 = sqrt(eps(T1)),
  rtol::Float64 = sqrt(eps(T1)),
  stop_first_stable::Bool = false,
) where {T1,T2}
  t_vals, x_vals = vals.t_vals, vals.x_vals

  # Initialize variables
  n0 = norm(x0)
  Fx = F(x0)
  stats = true
  tempx = similar(x0)

  verbose > 0 && @info log_header(
    [:iter, :step_size, :diff, :nFx, :proj_success],
    [Int, T1, T2, T2, Int],
    hdr_override = Dict(:step_size => "h", :diff => "‖x⁺-x‖", :nFx => "‖F(x)‖"),
  )

  # Main loop
  for iter in eachindex(t_vals)
    if (iter == 1)
      if project_x0
        proj_success = project_fun!(
          tempx,
          x0;
          step_size = 0.0,
          x = x0,
          atol = atol,
          rtol = rtol,
          verbose = proj_verbose,
        )
        x_vals[:, iter] .= tempx
        if !proj_success
          @error "Projection failed"
          return false
        end
        if norm(x_vals[:, iter] - x0) > atol + n0 * rtol
          @warn "x0 is not feasible, project on C $(proj_success) and continue."
        end
      else
        x_vals[:, iter] .= x0
      end
      verbose > 0 && @info log_row(Any[iter, NaN, NaN, norm(Fx), Int(proj_success)])
      continue
    end

    # Compute the update step:
    step_size = t_vals[iter] - t_vals[iter-1]
    proj_success = project_fun!(
      tempx,
      x_vals[:, iter-1] + step_size * Fx;
      step_size = step_size,
      x = x_vals[:, iter-1],
      atol = atol,
      rtol = rtol,
      verbose = proj_verbose,
    )
    x_vals[:, iter] .= tempx

    if !proj_success
      # @error "Projection failed"
      @warn "Projection failed"
      #return false
    end

    # Check for stability
    xdiff = norm(x_vals[:, iter] - x_vals[:, iter-1])
    if stop_first_stable && (xdiff < atol + n0 * rtol)
      break
    end

    Fx = F(x_vals[:, iter])
    verbose > 0 && @info log_row(Any[iter, step_size, xdiff, norm(Fx), Int(proj_success)])
  end

  return stats
end
