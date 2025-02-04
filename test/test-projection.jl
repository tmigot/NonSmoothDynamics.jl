@testset "Projection Functions" begin
  # Test projection onto a hyperplane
  x = [2.0, 3.0]
  a = [1.0, 1.0]
  b = 4.0
  @test isapprox(NonSmoothDynamics.project_hyperplane(x, a, b), [1.5, 2.5], atol = 1e-6)
  # NonSmoothDynamics.plot_hyperplane_projection(x, a, b)

  # Test projection onto a ball
  c = [0.0, 0.0]
  r = 1.0
  @test isapprox(
    NonSmoothDynamics.project_ball([2.0, 2.0], c, r),
    [0.707107, 0.707107],
    atol = 1e-6,
  )
  # NonSmoothDynamics.plot_ball_projection([2.0, 2.0], c, r)
  @test isapprox(NonSmoothDynamics.project_ball([0.5, 0.5], c, r), [0.5, 0.5], atol = 1e-6)
  # NonSmoothDynamics.plot_ball_projection([0.5, 0.5], c, r)

  # Test projection onto a box
  l = [0.0, 0.0]
  u = [1.0, 1.0]
  @test isapprox(NonSmoothDynamics.project_box([1.5, -0.5], l, u), [1.0, 0.0], atol = 1e-6)
  # NonSmoothDynamics.plot_box_projection([1.5, -0.5], l, u)

  # Test projection onto a simplex
  @test isapprox(
    NonSmoothDynamics.project_simplex([0.5, 0.5, 0.5]),
    [1 / 3, 1 / 3, 1 / 3],
    atol = 1e-6,
  )
  @test isapprox(
    NonSmoothDynamics.project_simplex([1.0, 0.5, 0.0]),
    [0.75, 0.25, 0.0],
    atol = 1e-6,
  )
  # NonSmoothDynamics.plot_simplex_projection([0.5, 0.5, 0.5])

  # Test projection onto an L2 norm ball
  @test isapprox(
    NonSmoothDynamics.project_l2_ball([3.0, 4.0], 5.0),
    [3.0, 4.0],
    atol = 1e-6,
  )
  @test isapprox(
    NonSmoothDynamics.project_l2_ball([6.0, 8.0], 5.0),
    [3.0, 4.0],
    atol = 1e-6,
  )

  # Test projection onto the positive orthant
  @test isapprox(
    NonSmoothDynamics.project_positive_orthant([-1.0, 2.0]),
    [0.0, 2.0],
    atol = 1e-6,
  )

  # Test projection onto a half-space
  x = [2.0, 2.0]
  a = [1.0, 1.0]
  b = 4.0
  @test isapprox(NonSmoothDynamics.project_halfspace(x, a, b), [2.0, 2.0], atol = 1e-6)
  x = [3.0, 3.0]
  @test isapprox(NonSmoothDynamics.project_halfspace(x, a, b), [2.0, 2.0], atol = 1e-6)
end

# Define projection functions
project_C1(x; kwargs...) = clamp.(x, 0.0, 1.0)    # Projection onto the box [0, 1]^n
project_C2(x; kwargs...) = x ./ sum(x)            # Projection onto the simplex

# Initial point
x0 = [1.5, 2.0, -0.5]

# Run the Boyle-Dykstra algorithm
x = NonSmoothDynamics.boyle_dykstra(x0, [project_C1, project_C2])

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
