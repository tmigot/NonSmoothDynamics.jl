@testset "Projection Functions" begin
  # Test projection onto a hyperplane
  x = [2.0, 3.0]
  a = [1.0, 1.0]
  b = 4.0
  proj_x = similar(x)
  NonSmoothDynamics.project_hyperplane!(proj_x, x, a, b)
  @test isapprox(proj_x, [1.5, 2.5], atol = 1e-6)
  # NonSmoothDynamics.plot_hyperplane_projection(x, a, b)

  # Test projection onto a ball
  x = [2.0, 2.0]
  c = [0.0, 0.0]
  r = 1.0
  proj_x = similar(x)
  NonSmoothDynamics.project_ball!(proj_x, x, c, r)
  @test isapprox(proj_x, [0.707107, 0.707107], atol = 1e-6)
  # NonSmoothDynamics.plot_ball_projection([2.0, 2.0], c, r)
  NonSmoothDynamics.project_ball!(proj_x, [0.5, 0.5], c, r)
  @test isapprox(proj_x, [0.5, 0.5], atol = 1e-6)
  # NonSmoothDynamics.plot_ball_projection([0.5, 0.5], c, r)

  # Test projection onto a box
  x = [1.5, -0.5]
  l = [0.0, 0.0]
  u = [1.0, 1.0]
  proj_x = similar(x)
  NonSmoothDynamics.project_box!(proj_x, x, l, u)
  @test isapprox(proj_x, [1.0, 0.0], atol = 1e-6)
  # NonSmoothDynamics.plot_box_projection([1.5, -0.5], l, u)

  # Test projection onto a simplex
  x = [0.5, 0.5, 0.5]
  proj_x = similar(x)
  NonSmoothDynamics.project_simplex!(proj_x, x)
  @test isapprox(proj_x, [1 / 3, 1 / 3, 1 / 3], atol = 1e-6)
  proj_x = similar(x)
  NonSmoothDynamics.project_simplex!(proj_x, [1.0, 0.5, 0.0])
  @test isapprox(proj_x, [0.75, 0.25, 0.0], atol = 1e-6)
  # NonSmoothDynamics.plot_simplex_projection([0.5, 0.5, 0.5])

  # Test projection onto an L2 norm ball
  x = [3.0, 4.0]
  proj_x = similar(x)
  NonSmoothDynamics.project_l2_ball!(proj_x, x, 5.0)
  @test isapprox(proj_x, [3.0, 4.0], atol = 1e-6)
  proj_x = similar(x)
  NonSmoothDynamics.project_l2_ball!(proj_x, [6.0, 8.0], 5.0)
  @test isapprox(proj_x, [3.0, 4.0], atol = 1e-6)

  # Test projection onto the positive orthant
  x = [-1.0, 2.0]
  proj_x = similar(x)
  NonSmoothDynamics.project_positive_orthant!(proj_x, x)
  @test isapprox(proj_x, [0.0, 2.0], atol = 1e-6)

  # Test projection onto a half-space
  x = [2.0, 2.0]
  a = [1.0, 1.0]
  b = 4.0
  proj_x = similar(x)
  NonSmoothDynamics.project_halfspace!(proj_x, x, a, b)
  @test isapprox(proj_x, [2.0, 2.0], atol = 1e-6)
  x = [3.0, 3.0]
  NonSmoothDynamics.project_halfspace!(proj_x, x, a, b)
  @test isapprox(proj_x, [2.0, 2.0], atol = 1e-6)
end

# Define projection functions
function project_C1(sol, x; kwargs...)
  sol .= clamp.(x, 0.0, 1.0)    # Projection onto the box [0, 1]^n
  return true
end
function project_C2(sol, x; kwargs...)
  sol .= x ./ sum(x)            # Projection onto the simplex
  return true
end

# Initial point
x0 = [1.5, 2.0, -0.5]
x = similar(x0)

# Run the Boyle-Dykstra algorithm
NonSmoothDynamics.boyle_dykstra!(x, x0, [project_C1, project_C2])

println("Projection: ", x)

# Plot the box and simplex
box_x = [0.0, 1.0, 1.0, 0.0, 0.0]
box_y = [0.0, 0.0, 1.0, 1.0, 0.0]

triangle_x = [0.0, 1.0, 0.0]
triangle_y = [0.0, 0.0, 1.0]

using Plots

plot(box_x, box_y, label = "Box [0, 1]^n", color = :blue, lw = 2)
plot!(
  [triangle_x..., triangle_x[1]],
  [triangle_y..., triangle_y[1]],
  label = "Simplex",
  color = :green,
  lw = 2,
)
scatter!(
  [x0[1]],
  [x0[2]],
  label = "Initial Point",
  color = :orange,
  marker = :star,
  markersize = 8,
)
scatter!([x[1]], [x[2]], label = "Final Projection", color = :purple, markersize = 8)
plot!(
  aspect_ratio = :equal,
  xlabel = "x₁",
  ylabel = "x₂",
  title = "Projection onto Box ∩ Simplex using Boyle-Dykstra",
)
