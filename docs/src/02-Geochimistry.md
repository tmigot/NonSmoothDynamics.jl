```@example ex1
using LinearAlgebra, SparseArrays

x0 = [1; 0.5]
S = [-1 1]' # (nxm) matrix
Q = [1 1]' # (nx(n-m)) matrix
tau(x; kr = 1, kp = 2) = kr * x[1] - kp # size m

Stau(x) = (S * tau(x))[:, 1]
```

```@example ex1
using NonSmoothDynamics

project_C0(y; kwargs...) = NonSmoothDynamics.numerical_projection(
  y, I, 
  sparse(Q'), Q' * x0,
  spzeros(0, 2), ones(0),
  zeros(2), Inf * ones(2),
)
function project_rate(y; step_size = step_size, x = x, kwargs...)
   λ = NonSmoothDynamics.numerical_projection(y, diagm(0 => Stau(x)))
   return x + step_size * diagm(0 => Stau(x)) * λ
end
function project_intersection(y; step_size = step_size, x = x, kwargs...)
  return NonSmoothDynamics.boyle_dykstra(
    y,
    [project_C0, y -> project_rate(y, step_size = step_size, x = x)],
)
end
```

# Simulate the PDS

```@example ex1
x_vals, num_iter, converged = NonSmoothDynamics.projected_dynamical_system(x0, x -> Stau(x), project_intersection)

println("Final State: ", x_vals[:, end])
println("Converged: ", converged, " in ", num_iter, " iterations.")
```

## Visualization

To plot the trajectory:

```@example ex1
using Plots

# Extract x1 and x2 trajectories
x1_vals = x_vals[1, :]
x2_vals = x_vals[2, :]

# Plot the trajectory
scatter!(x1_vals, x2_vals, label="Trajectory", color=:red, markersize=4)
scatter!([x0[1]], [x0[2]], label="Initial Point", color=:orange, marker=:star, markersize=8)
scatter!([x_vals[1, end]], [x_vals[2, end]], label="Final Point", color=:purple, markersize=8)
```

```@example ex1
t_interp = collect(0:0.001:(203 * 0.1))  # Fine-grained time vector
x_interp = hcat([NonSmoothDynamics.interpolated_solution(t) for t in t_interp]...)

step_size = 0.01
num_iters = size(x_vals, 2) - 1
t_interp = collect(0:step_size:(num_iters * step_size))  # Fine-grained time vector
f = NonSmoothDynamics.get_interpolation(x0, x_vals, step_size, num_iters)

plot(t_interp, x_vals[1, :], label="x₁ function of time", color=:red, markersize=4)
plot!(t_interp, x_vals[2, :], label="x₂ function of time", color=:green, lw=2)
```
