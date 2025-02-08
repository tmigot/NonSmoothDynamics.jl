# Projected dynamical system for Optimization

## Problem definition

Let’s minimize the function:
```math
\begin{aligned}
f(x) = \frac{1}{2} \|x - 1\|^2_2,
\end{aligned}
```
subject to $x \in [0, 1]^n$.

- **Gradient:** $F(x) = \nabla f(x) = x - 1$.
- **Projection:** $P_C(x) = \text{clamp}(x, 0, 1)$.

```@example ex1
# Define the vector field F(x)
F(x) = -(x .- 1.0)

# Define projection onto the box [0, 1]^2
function project_C(sol, x; kwargs...)
    sol .= clamp.(x, 0.0, 1.0)
    return true
end

# Initial state
x0 = [2.0, -1.0]
```

## PDS

```@example ex1
# Initial state
x0 = Vector{Float64}(undef, 2)
project_C(x0, [2.0, -1.0])

using NonSmoothDynamics
# Simulate the PDS
x_vals, t_vals, converged = NonSmoothDynamics.projected_dynamical_system(x0, F, project_C, 0.0, 1.0)

# Print results
println("Final State: ", x_vals[:, end])
```

## Visualization

To plot the trajectory:

```@example ex1
using Plots

# Extract x1 and x2 trajectories
x1_vals = x_vals[1, :]
x2_vals = x_vals[2, :]

# Create the box [0, 1]^2 for visualization
box_x = [0.0, 1.0, 1.0, 0.0, 0.0]
box_y = [0.0, 0.0, 1.0, 1.0, 0.0]

# Plot the trajectory
plot(box_x, box_y, label="Feasible Set [0, 1]^2", color=:blue)
scatter!(x1_vals, x2_vals, label="Trajectory", color=:red, markersize=4)
scatter!([x0[1]], [x0[2]], label="Initial Point", color=:orange, marker=:star, markersize=8)
scatter!([x_vals[1, end]], [x_vals[2, end]], label="Final Point", color=:purple, markersize=8)
plot!(xlabel="x₁", ylabel="x₂", title="Projected Dynamical System Trajectory", aspect_ratio=:equal)
```
