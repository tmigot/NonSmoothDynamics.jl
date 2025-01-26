using Plots

# Function to plot projection onto a line or hyperplane
function plot_hyperplane_projection(x::Vector{Float64}, a::Vector{Float64}, b::Float64)
    proj_x = project_hyperplane(x, a, b)

    xs = range(-1, stop=5, length=100)
    ys = (b .- a[1] .* xs) ./ a[2]

    plot(xs, ys, label="Hyperplane", color=:blue)
    scatter!([x[1]], [x[2]], label="Original Vector", color=:red)
    scatter!([proj_x[1]], [proj_x[2]], label="Projection", color=:green)
    plot!([x[1], proj_x[1]], [x[2], proj_x[2]], label="Projection Path", linestyle=:dash)
end

# Function to plot projection onto a ball
function plot_ball_projection(x::Vector{Float64}, c::Vector{Float64}, r::Float64)
    proj_x = project_ball(x, c, r)

    θ = range(0, stop=2π, length=100)
    circle_x = r .* cos.(θ) .+ c[1]
    circle_y = r .* sin.(θ) .+ c[2]

    plot(circle_x, circle_y, label="Ball Boundary", color=:blue, aspect_ratio=1)
    scatter!([x[1]], [x[2]], label="Original Vector", color=:red)
    scatter!([proj_x[1]], [proj_x[2]], label="Projection", color=:green)
    plot!([x[1], proj_x[1]], [x[2], proj_x[2]], label="Projection Path", linestyle=:dash)
end

# Function to plot projection onto a box
function plot_box_projection(x::Vector{Float64}, l::Vector{Float64}, u::Vector{Float64})
    proj_x = project_box(x, l, u)

    box_x = [l[1], u[1], u[1], l[1], l[1]]
    box_y = [l[2], l[2], u[2], u[2], l[2]]

    plot(box_x, box_y, label="Box", color=:blue)
    scatter!([x[1]], [x[2]], label="Original Vector", color=:red)
    scatter!([proj_x[1]], [proj_x[2]], label="Projection", color=:green)
    plot!([x[1], proj_x[1]], [x[2], proj_x[2]], label="Projection Path", linestyle=:dash)
end

# Function to plot projection onto a simplex
function plot_simplex_projection(x::Vector{Float64})
    proj_x = project_simplex(x)

    simplex_x = [1.0, 0.0, 0.0, 1.0]
    simplex_y = [0.0, 1.0, 0.0, 0.0]

    plot(simplex_x, simplex_y, label="Simplex", lw=2, color=:blue)
    scatter!([x[1]], [x[2]], label="Original Vector", color=:red)
    scatter!([proj_x[1]], [proj_x[2]], label="Projection", color=:green)
    plot!([x[1], proj_x[1]], [x[2], proj_x[2]], label="Projection Path", linestyle=:dash)
end
