using NonSmoothDynamics

# Define the vector field F(x)
F(x) = x .- 1.0

# Define projection onto the box [0, 1]^2
project_C(x; kwargs...) = clamp.(x, 0.0, 1.0)

# Initial state
x0 = project_C([2.0, -1.0])

# Simulate the PDS
x_vals, t_vals, converged =
  NonSmoothDynamics.projected_dynamical_system(x0, F, project_C, 0.0, 1.0)
