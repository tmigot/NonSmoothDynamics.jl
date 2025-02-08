using LinearAlgebra, SparseArrays
using NonSmoothDynamics

# Define the vector field F(x)
F(x) = -(x .- 1.0)

# Define projection onto the box [0, 1]^2
project_C(x; kwargs...) = clamp.(x, 0.0, 1.0)

# Initial state
x0 = project_C([2.0, -1.0])

# Simulate the PDS
x_vals, t_vals, converged =
  NonSmoothDynamics.projected_dynamical_system(x0, F, project_C, 0.0, 1.0, verbose = 1)

# Example 2
F(x) = -[
  2 * x[1] + 8 / 3 * x[2] - 34
  2 * x[2] + 5 / 4 * x[1] - 24.25
]

function project_moving_set!(sol, y; x = x, kwargs...)
  λ = similar(sol)
  Ain = [
    1 0.9
    0.9 1
  ]
  bin = [14.4; 14.1]
  proj_success = NonSmoothDynamics.numerical_projection!(
    λ,
    y,
    I,
    spzeros(0, 2),
    ones(0),
    Ain,
    bin,
    zeros(2),
    Inf * ones(2),
  )
  sol .= λ
  return proj_success
end

x0 = zeros(2)
t0, tf = 0.0, 30.0
x_vals, t_vals, converged =
  NonSmoothDynamics.projected_dynamical_system(x0, F, project_moving_set!, t0, tf, 300)

@test x_vals[:, end] ≈ [5, 9]

# Example 3

x0 = [1; 0.5]
n = length(x0) # number of species
m = 2 # number of reactions

# Stoichiometry matrix
S = [-1 1]' # (nxm) matrix
# Matrix of conservation
Q = [1 1]' # (nx(n-m)) matrix

# Reaction rate
tau(x; kr = 1, kp = 2) = kr * x[1] - kp # size m
Stau(x) = (S*tau(x))[:, 1]

project_C0!(sol, y; kwargs...) = NonSmoothDynamics.numerical_projection!(
  sol,
  y,
  I,
  sparse(Q'),
  Q' * x0, # satisfy conservation equation
  spzeros(0, n),
  ones(0),
  zeros(n),
  Inf * ones(n), # non-negative
)
function project_rate!(sol, y; step_size = step_size, x = x, kwargs...)
  λ = similar(sol)
  proj_success = NonSmoothDynamics.numerical_projection!(
    λ,
    y,
    diagm(0 => Stau(x)),
    spzeros(0, 2),
    ones(0),
    spzeros(0, 2),
    ones(0),
    zeros(2),
    ones(2),
  )
  sol .= x + step_size * diagm(0 => Stau(x)) * λ
  return proj_success
end
function project_intersection!(sol, y; step_size = step_size, x = x, kwargs...)
  return NonSmoothDynamics.boyle_dykstra!(
    sol,
    y,
    [
      project_C0!,
      (sol, y; kwargs...) -> project_rate!(sol, y, step_size = step_size, x = x, kwargs...),
    ],
  )
end

x_vals, t_vals, converged = NonSmoothDynamics.projected_dynamical_system(
  x0,
  x -> Stau(x),
  project_intersection!,
  0.0,
  2.0,
  100; # run the reaction with 100 discretization point between 0 and 2.
  project_x0 = false,
)
