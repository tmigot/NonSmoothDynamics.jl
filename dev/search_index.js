var documenterSearchIndex = {"docs":
[{"location":"91-developer/#dev_docs","page":"Developer documentation","title":"Developer documentation","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"note: Contributing guidelines\nIf you haven't, please read the Contributing guidelines first.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If you want to make contributions to this package that involves code, then this guide is for you.","category":"page"},{"location":"91-developer/#First-time-clone","page":"Developer documentation","title":"First time clone","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: If you have writing rights\nIf you have writing rights, you don't have to fork. Instead, simply clone and skip ahead. Whenever upstream is mentioned, use origin instead.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If this is the first time you work with this repository, follow the instructions below to clone the repository.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fork this repo\nClone your repo (this will create a git remote called origin)\nAdd this repo as a remote:\ngit remote add upstream https://github.com/tmigot/NonSmoothDynamics.jl","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"This will ensure that you have two remotes in your git: origin and upstream. You will create branches and push to origin, and you will fetch and update your local main branch from upstream.","category":"page"},{"location":"91-developer/#Linting-and-formatting","page":"Developer documentation","title":"Linting and formatting","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Install a plugin on your editor to use EditorConfig. This will ensure that your editor is configured with important formatting settings.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We use https://pre-commit.com to run the linters and formatters. In particular, the Julia code is formatted using JuliaFormatter.jl, so please install it globally first:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # Press ]\npkg> activate\npkg> add JuliaFormatter","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To install pre-commit, we recommend using pipx as follows:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"# Install pipx following the link\npipx install pre-commit","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"With pre-commit installed, activate it as a pre-commit hook:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit install","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To run the linting and formatting manually, enter the command below:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit run -a","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Now, you can only commit if all the pre-commit tests pass.","category":"page"},{"location":"91-developer/#Testing","page":"Developer documentation","title":"Testing","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"As with most Julia packages, you can just open Julia in the repository folder, activate the environment, and run test:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # press ]\npkg> activate .\npkg> test","category":"page"},{"location":"91-developer/#Working-on-a-new-issue","page":"Developer documentation","title":"Working on a new issue","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We try to keep a linear history in this repo, so it is important to keep your branches up-to-date.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fetch from the remote and fast-forward your local main\ngit fetch upstream\ngit switch main\ngit merge --ff-only upstream/main\nBranch from main to address the issue (see below for naming)\ngit switch -c 42-add-answer-universe\nPush the new local branch to your personal remote repository\ngit push -u origin 42-add-answer-universe\nCreate a pull request to merge your remote branch into the org main.","category":"page"},{"location":"91-developer/#Branch-naming","page":"Developer documentation","title":"Branch naming","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If there is an associated issue, add the issue number.\nIf there is no associated issue, and the changes are small, add a prefix such as \"typo\", \"hotfix\", \"small-refactor\", according to the type of update.\nIf the changes are not small and there is no associated issue, then create the issue first, so we can properly discuss the changes.\nUse dash separated imperative wording related to the issue (e.g., 14-add-tests, 15-fix-model, 16-remove-obsolete-files).","category":"page"},{"location":"91-developer/#Commit-message","page":"Developer documentation","title":"Commit message","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Use imperative or present tense, for instance: Add feature or Fix bug.\nHave informative titles.\nWhen necessary, add a body with details.\nIf there are breaking changes, add the information to the commit message.","category":"page"},{"location":"91-developer/#Before-creating-a-pull-request","page":"Developer documentation","title":"Before creating a pull request","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: Atomic git commits\nTry to create \"atomic git commits\" (recommended reading: The Utopic Git History).","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Make sure the tests pass.\nMake sure the pre-commit tests pass.\nFetch any main updates from upstream and rebase your branch, if necessary:\ngit fetch upstream\ngit rebase upstream/main BRANCH_NAME\nThen you can open a pull request and work with the reviewer to address any issues.","category":"page"},{"location":"91-developer/#Building-and-viewing-the-documentation-locally","page":"Developer documentation","title":"Building and viewing the documentation locally","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Following the latest suggestions, we recommend using LiveServer to build the documentation. Here is how you do it:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Run julia --project=docs to open Julia in the environment of the docs.\nIf this is the first time building the docs\nPress ] to enter pkg mode\nRun pkg> dev . to use the development version of your package\nPress backspace to leave pkg mode\nRun julia> using LiveServer\nRun julia> servedocs()","category":"page"},{"location":"91-developer/#Making-a-new-release","page":"Developer documentation","title":"Making a new release","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To create a new release, you can follow these simple steps:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Create a branch release-x.y.z\nUpdate version in Project.toml\nUpdate the CHANGELOG.md:\nRename the section \"Unreleased\" to \"[x.y.z] - yyyy-mm-dd\" (i.e., version under brackets, dash, and date in ISO format)\nAdd a new section on top of it named \"Unreleased\"\nAdd a new link in the bottom for version \"x.y.z\"\nChange the \"[unreleased]\" link to use the latest version - end of line, vx.y.z ... HEAD.\nCreate a commit \"Release vx.y.z\", push, create a PR, wait for it to pass, merge the PR.\nGo back to main screen and click on the latest commit (link: https://github.com/tmigot/NonSmoothDynamics.jl/commit/main)\nAt the bottom, write @JuliaRegistrator register","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"After that, you only need to wait and verify:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Wait for the bot to comment (should take < 1m) with a link to a PR to the registry\nFollow the link and wait for a comment on the auto-merge\nThe comment should said all is well and auto-merge should occur shortly\nAfter the merge happens, TagBot will trigger and create a new GitHub tag. Check on https://github.com/tmigot/NonSmoothDynamics.jl/releases\nAfter the release is create, a \"docs\" GitHub action will start for the tag.\nAfter it passes, a deploy action will run.\nAfter that runs, the stable docs should be updated. Check them and look for the version number.","category":"page"},{"location":"95-reference/#reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"95-reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Modules = [NonSmoothDynamics]","category":"page"},{"location":"95-reference/#NonSmoothDynamics.boyle_dykstra!-Union{Tuple{S}, Tuple{S, S, Vector{Function}}} where S","page":"Reference","title":"NonSmoothDynamics.boyle_dykstra!","text":"boyle_dykstra!(sol::Vector{Float64}, x0::Vector{Float64}, projections::Vector{Function}; tol::Float64=1e-6, max_iter::Int=1000)\n\nImplements the Boyle-Dykstra algorithm to project a point onto the intersection of multiple convex sets.\n\nArguments\n\nsol::Vector{Float64}: The solution vector to be updated.\nx0::Vector{Float64}: The initial point to be projected.\nprojections::Vector{Function}: A vector of projection functions, where each function computes the projection onto a specific convex set.\ntol::Float64=1e-6: Convergence tolerance. The algorithm stops if the change in x between iterations is less than this value.\nmax_iter::Int=1000: Maximum number of iterations to run the algorithm.\n\nReturns\n\nx::Vector{Float64}: The projection of x0 onto the intersection of the convex sets.\nnum_iter::Int: The number of iterations performed.\nconverged::Bool: Whether the algorithm converged within the maximum number of iterations.\n\nExample\n\n# Define projection functions\nproject_C1(sol, x) = sol .= clamp.(x, 0.0, 1.0)    # Projection onto the box [0, 1]^n\nproject_C2(sol, x) = sol .= x ./ sum(x)            # Projection onto the simplex\n\n# Initial point\nx0 = [1.5, 2.0, -0.5]\n\n# Run the Boyle-Dykstra algorithm\nx, num_iter, converged = boyle_dykstra(x0, [project_C1, project_C2])\n\nprintln(\"Projection: \", x)\nprintln(\"Converged: \", converged, \" in \", num_iter, \" iterations.\")\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#NonSmoothDynamics.numerical_projection!-Union{Tuple{S}, Tuple{S, S}, Tuple{S, S, Any}, Tuple{S, S, Any, Any}, Tuple{S, S, Any, Any, Any}, Tuple{S, S, Vararg{Any, 4}}, Tuple{S, S, Vararg{Any, 5}}, Tuple{S, S, Vararg{Any, 6}}, Tuple{S, S, Vararg{Any, 7}}} where S","page":"Reference","title":"NonSmoothDynamics.numerical_projection!","text":"numerical_projection!(sol::Vector{Float64},\n                      x::Vector{Float64}, Px::AbstractMatrix,\n                      Aeq::AbstractMatrix, beq::AbstractVector,\n                      Ain::AbstractMatrix, bin::AbstractVector,\n                      l::Vector{Float64}, u::Vector{Float64}) -> Vector{Float64}\n\nComputes the numerical projection of a vector x onto a constrained set defined by equality and inequality constraints, as well as variable bounds.\n\nArguments\n\nsol::Vector{Float64}: The solution vector to be updated.\nx::AbstractVector: The input vector to be projected.\nPx::AbstractMatrix: Projection matrix.\nAeq::AbstractMatrix: Coefficient matrix for equality constraints Aeq x = beq.\nbeq::AbstractVector: Right-hand side of the equality constraints.\nAin::AbstractMatrix: Coefficient matrix for inequality constraints Ain x ≤ bin.\nbin::AbstractVector: Right-hand side of the inequality constraints.\nl::AbstractVector: Lower bounds for the variables.\nu::AbstractVector: Upper bounds for the variables.\n\nReturns\n\nAbstractVector: The projection of x that satisfies the constraints.\n\nExample\n\n# Example data\nx   = [0.5, -1.0]\nPx  = I\nAeq = [1.0 1.0]          # Equality constraint: x₁ + x₂ = 1\nbeq = [1.0]\nAin = [1.0 0.0; 0.0 1.0] # Inequality constraints: x₁ ≤ 0.8, x₂ ≤ 0.6\nbin = [0.8, 0.6]\nl = [-Inf, -Inf]         # No lower bounds\nu = [Inf, Inf]           # No upper bounds\n\n# Compute the projection\nproj_x = numerical_projection(x, Px, Aeq, beq, Ain, bin, l, u)\nprintln(\"Projected x: \", proj_x)\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#NonSmoothDynamics.project_ball!-Union{Tuple{T}, Tuple{S}, Tuple{S, S, S, T}} where {S, T}","page":"Reference","title":"NonSmoothDynamics.project_ball!","text":"project_ball!(sol::Vector{Float64}, x::Vector{Float64}, c::Vector{Float64}, r::Float64) -> Vector{Float64}\n\nProjects a vector x onto the ball centered at c with radius r.\n\nArguments\n\nsol::Vector{Float64}: The solution vector to be updated.\nx::Vector{Float64}: The vector to project.\nc::Vector{Float64}: The center of the ball.\nr::Float64: The radius of the ball.\n\nReturns\n\ntrue if the projection is successful.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#NonSmoothDynamics.project_box!-Union{Tuple{S}, NTuple{4, S}} where S","page":"Reference","title":"NonSmoothDynamics.project_box!","text":"project_box!(sol::Vector{Float64}, x::Vector{Float64}, l::Vector{Float64}, u::Vector{Float64}) -> Vector{Float64}\n\nProjects a vector x onto the box defined by lower bounds l and upper bounds u.\n\nArguments\n\nsol::Vector{Float64}: The solution vector to be updated.\nx::Vector{Float64}: The vector to project.\nl::Vector{Float64}: The lower bounds of the box.\nu::Vector{Float64}: The upper bounds of the box.\n\nReturns\n\ntrue if the projection is successful.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#NonSmoothDynamics.project_halfspace!-Union{Tuple{T}, Tuple{S}, Tuple{S, S, S, T}} where {S, T}","page":"Reference","title":"NonSmoothDynamics.project_halfspace!","text":"project_halfspace!(sol::Vector{Float64}, x::Vector{Float64}, a::Vector{Float64}, b::Float64) -> Vector{Float64}\n\nProjects a vector x onto the half-space defined by a^T * x <= b.\n\nArguments\n\nsol::Vector{Float64}: The solution vector to be updated.\nx::Vector{Float64}: The vector to project.\na::Vector{Float64}: The normal vector defining the half-space.\nb::Float64: The offset of the half-space.\n\nReturns\n\ntrue if the projection is successful.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#NonSmoothDynamics.project_hyperplane!-Union{Tuple{T}, Tuple{S}, Tuple{S, S, S, T}} where {S, T}","page":"Reference","title":"NonSmoothDynamics.project_hyperplane!","text":"project_hyperplane!(sol::Vector{Float64}, x::Vector{Float64}, a::Vector{Float64}, b::Float64) -> Vector{Float64}\n\nProjects a vector x onto the hyperplane defined by a^T * x = b.\n\nArguments\n\nsol::Vector{Float64}: The solution vector to be updated.\nx::Vector{Float64}: The vector to project.\na::Vector{Float64}: The normal vector defining the hyperplane.\nb::Float64: The offset of the hyperplane.\n\nReturns\n\ntrue if the projection is successful.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#NonSmoothDynamics.project_l2_ball!-Union{Tuple{T}, Tuple{S}, Tuple{S, S, T}} where {S, T}","page":"Reference","title":"NonSmoothDynamics.project_l2_ball!","text":"project_l2_ball!(sol::Vector{Float64}, x::Vector{Float64}, r::Float64) -> Vector{Float64}\n\nProjects a vector x onto the L2 norm ball of radius r.\n\nArguments\n\nsol::Vector{Float64}: The solution vector to be updated.\nx::Vector{Float64}: The vector to project.\nr::Float64: The radius of the L2 norm ball.\n\nReturns\n\ntrue if the projection is successful.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#NonSmoothDynamics.project_positive_orthant!-Union{Tuple{S}, Tuple{S, S}} where S","page":"Reference","title":"NonSmoothDynamics.project_positive_orthant!","text":"project_positive_orthant!(sol::Vector{Float64}, x::Vector{Float64}) -> Vector{Float64}\n\nProjects a vector x onto the positive orthant (all elements >= 0).\n\nArguments\n\nsol::Vector{Float64}: The solution vector to be updated.\nx::Vector{Float64}: The vector to project.\n\nReturns\n\ntrue if the projection is successful.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#NonSmoothDynamics.project_simplex!-Union{Tuple{S}, Tuple{S, S}} where S","page":"Reference","title":"NonSmoothDynamics.project_simplex!","text":"project_simplex!(sol::Vector{Float64}, x::Vector{Float64}) -> Vector{Float64}\n\nProjects a vector x onto the probability simplex, defined as:     {y ∈ R^n : yi ≥ 0, ∑ yi = 1}\n\nThis implementation avoids sorting by iteratively adjusting the elements of x to meet the simplex constraints.\n\nArguments\n\nsol::Vector{Float64}: The solution vector to be updated.\nx::Vector{Float64}: The input vector to be projected.\n\nReturns\n\ntrue if the projection is successful.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#NonSmoothDynamics.projected_dynamical_system!-Union{Tuple{T2}, Tuple{T1}, Tuple{NonSmoothDynamics.InterpolationStruct{T1, T2}, Vector{Float64}, Function, Function}} where {T1, T2}","page":"Reference","title":"NonSmoothDynamics.projected_dynamical_system!","text":"projected_dynamical_system!(vals::InterpolationStruct,\n                            x0::Vector{Float64}, F::Function, project_fun!::Function;\n                            kwargs...)\nprojected_dynamical_system(x0::Vector{Float64}, F::Function, project_fun!::Function,\n                           args...; kwargs...)\n\nThe args... are passed to InterpolationStruct constructor.\n\nSimulates the discretized dynamics of a Projected Dynamical System (PDS):     x(t+1) = P_{C}(x(t) - h * F(x(t)))\n\nArguments\n\nvals::InterpolationStruct: Interpolation struct containing the time and state values.\nx0::Vector{Float64}: Initial state (starting point within the set C).\nF::Function: The vector field defining the dynamics (e.g., gradient, payoff vector, etc.).\nproject_fun!::Function: A function to compute the projection onto the feasible set C.\nverbose::Int=0: Verbosity level.\nproj_verbose::Int=0: Verbosity level for the projection step.\nproject_x0::Bool = true:\natol::Float64=1e-6: Convergence tolerance. The algorithm stops if the norm of the change is less than this value.\nrtol::Float64=1e-6: Convergence tolerance.\nstop_first_stable::Bool=true: Stop the simulation if the first stable point is found.\n\nReturns\n\nstats::Bool: Whether the system computation was successful.\n\nvals.x_vals::Matrix{Float64}: Matrix where each column is the state of the system at a time step.\n\nExample\n\n# Define the vector field F(x) = x - 1 (gradient of f(x) = 0.5 * ||x - 1||^2)\nF(x) = x - 1.0\n\n# Define projection onto the box [0, 1]^2\nproject_fun!(x) = clamp.(x, 0.0, 1.0)\n\n# Initial point\nx0 = [2.0, -1.0]\n\n# Simulate the PDS\nx_vals, t_vals, converged = projected_dynamical_system(x0, F, project_fun!)\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#NonSmoothDynamics.projected_dynamical_system-Tuple{Vector{Float64}, Function, Function, Vararg{Any}}","page":"Reference","title":"NonSmoothDynamics.projected_dynamical_system","text":"No documentation found.\n\nBinding NonSmoothDynamics.projected_dynamical_system! does not exist.\n\n\n\n\n\n","category":"method"},{"location":"01-PDS/#Projected-dynamical-system-for-Optimization","page":"Projected dynamical system for Optimization","title":"Projected dynamical system for Optimization","text":"","category":"section"},{"location":"01-PDS/#Problem-definition","page":"Projected dynamical system for Optimization","title":"Problem definition","text":"","category":"section"},{"location":"01-PDS/","page":"Projected dynamical system for Optimization","title":"Projected dynamical system for Optimization","text":"Let’s minimize the function: [ f(x) = \\frac{1}{2} \\|x - 1\\|^2, ] subject to (x \\in [0, 1]^n).","category":"page"},{"location":"01-PDS/","page":"Projected dynamical system for Optimization","title":"Projected dynamical system for Optimization","text":"Gradient: (F(x) = \\nabla f(x) = x - 1).\nProjection: (P_C(x) = \\text{clamp}(x, 0, 1)).","category":"page"},{"location":"01-PDS/","page":"Projected dynamical system for Optimization","title":"Projected dynamical system for Optimization","text":"# Define the vector field F(x)\nF(x) = -(x .- 1.0)\n\n# Define projection onto the box [0, 1]^2\nfunction project_C(sol, x; kwargs...)\n    sol .= clamp.(x, 0.0, 1.0)\n    return true\nend\n\n# Initial state\nx0 = [2.0, -1.0]","category":"page"},{"location":"01-PDS/#PDS","page":"Projected dynamical system for Optimization","title":"PDS","text":"","category":"section"},{"location":"01-PDS/","page":"Projected dynamical system for Optimization","title":"Projected dynamical system for Optimization","text":"# Initial state\nx0 = Vector{Float64}(undef, 2)\nproject_C(x0, [2.0, -1.0])\n\nusing NonSmoothDynamics\n# Simulate the PDS\nx_vals, t_vals, converged = NonSmoothDynamics.projected_dynamical_system(x0, F, project_C, 0.0, 1.0)\n\n# Print results\nprintln(\"Final State: \", x_vals[:, end])","category":"page"},{"location":"01-PDS/#Visualization","page":"Projected dynamical system for Optimization","title":"Visualization","text":"","category":"section"},{"location":"01-PDS/","page":"Projected dynamical system for Optimization","title":"Projected dynamical system for Optimization","text":"To plot the trajectory:","category":"page"},{"location":"01-PDS/","page":"Projected dynamical system for Optimization","title":"Projected dynamical system for Optimization","text":"using Plots\n\n# Extract x1 and x2 trajectories\nx1_vals = x_vals[1, :]\nx2_vals = x_vals[2, :]\n\n# Create the box [0, 1]^2 for visualization\nbox_x = [0.0, 1.0, 1.0, 0.0, 0.0]\nbox_y = [0.0, 0.0, 1.0, 1.0, 0.0]\n\n# Plot the trajectory\nplot(box_x, box_y, label=\"Feasible Set [0, 1]^2\", color=:blue)\nscatter!(x1_vals, x2_vals, label=\"Trajectory\", color=:red, markersize=4)\nscatter!([x0[1]], [x0[2]], label=\"Initial Point\", color=:orange, marker=:star, markersize=8)\nscatter!([x_vals[1, end]], [x_vals[2, end]], label=\"Final Point\", color=:purple, markersize=8)\nplot!(xlabel=\"x₁\", ylabel=\"x₂\", title=\"Projected Dynamical System Trajectory\", aspect_ratio=:equal)","category":"page"},{"location":"03-GeochimistryWM/#Projected-dynamical-system-with-constraints-on-derivative","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"","category":"section"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"In this tutorial, we present another application of PDS in geochemical system with precipitation-dissolution reactions. Ordinary differential equations are very often used to model aqueous reactions. However, non smoothness is induced by full dissolution. Moreover, one reaction can include several minerals and one mineral can participate in several reactions.","category":"page"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"The article Erhel, J., Hamlat, B., & Migot, T. (2024). A projected dynamical system approach to mineral precipitation-dissolution reactions in geochemistry. introduces a PDS model for this type of reaction involving a constraints on the reaction rate.","category":"page"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"(Image: \"PDS in geochemistry\")","category":"page"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"using LinearAlgebra, SparseArrays\nusing NonSmoothDynamics\nusing Plots","category":"page"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"The first example is a single reaction with an aqueous species W and a mineral species M","category":"page"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"(Image: \"Equation W->M\")","category":"page"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"First, define constants relative to this example.","category":"page"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"x0 = [1; 0.5]\nn = length(x0) # number of species\nm = 2 # number of reactions\n\n# Stoichiometry matrix\nS = [-1 1]' # (nxm) matrix\n# Matrix of conservation\nQ = [1 1]' # (nx(n-m)) matrix\n\n# Reaction rate\ntau(x; kr = 1, kp = 2) = kr * x[1] - kp # size m\nStau(x) = (S * tau(x))[:, 1]","category":"page"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"Then, define the projection operator.","category":"page"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"project_C0!(sol, y; kwargs...) = NonSmoothDynamics.numerical_projection!(\n  sol, y, I,\n  sparse(Q'), Q' * x0, # satisfy conservation equation\n  spzeros(0, n), ones(0),\n  zeros(n), Inf * ones(n), # non-negative\n)\nfunction project_rate!(sol, y; step_size = step_size, x = x, kwargs...)\n   λ = similar(sol)\n   proj_success = NonSmoothDynamics.numerical_projection!(\n   λ, y, diagm(0 => Stau(x)),\n   spzeros(0, 2), ones(0),\n   spzeros(0, 2), ones(0),\n   zeros(2), ones(2))\n   sol .= x + step_size * diagm(0 => Stau(x)) * λ\n   return proj_success\nend\nfunction project_intersection!(sol, y; step_size = step_size, x = x, kwargs...)\n  return NonSmoothDynamics.boyle_dykstra!(\n    sol, y,\n    [project_C0!, (sol, y; kwargs...) -> project_rate!(sol, y, step_size = step_size, x = x, kwargs...)],\n)\nend","category":"page"},{"location":"03-GeochimistryWM/#Simulate-the-PDS","page":"Projected dynamical system with constraints on derivative","title":"Simulate the PDS","text":"","category":"section"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"x_vals, t_vals, converged = NonSmoothDynamics.projected_dynamical_system(\n  x0, x -> Stau(x), project_intersection!,\n  0.0, 2.0, 100; # run the reaction with 100 discretization point between 0 and 2.\n  project_x0 = false\n)\n\nprintln(\"Final State: \", x_vals[:, end])","category":"page"},{"location":"03-GeochimistryWM/#Visualization","page":"Projected dynamical system with constraints on derivative","title":"Visualization","text":"","category":"section"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"To plot the trajectory:","category":"page"},{"location":"03-GeochimistryWM/","page":"Projected dynamical system with constraints on derivative","title":"Projected dynamical system with constraints on derivative","text":"using Plots\n\nplot(t_vals, x_vals[1, :], label=\"x₁ function of time\", color=:red, markersize=4)\nplot!(t_vals, x_vals[2, :], label=\"x₂ function of time\", color=:green, lw=2)","category":"page"},{"location":"02-MovingSet/#Projected-dynamical-system-with-moving-set","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"","category":"section"},{"location":"02-MovingSet/","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"It is possible to apply the PDS framework to constraint set that depend on the state x. One example is the study of nonsmooth dynamics in game theory.","category":"page"},{"location":"02-MovingSet/","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"In Migot, T., & Cojocaru, M. G. (2020). Nonsmooth dynamics of generalized Nash games. J. Nonlinear Var. Anal, 1(4), 27-44., the authors study a nonsmooth dynamics whose stable points are generalized Nash equilibrium for noncooperative games.","category":"page"},{"location":"02-MovingSet/","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"(Image: \"Example 6.1\")","category":"page"},{"location":"02-MovingSet/","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"using LinearAlgebra, SparseArrays\nusing NonSmoothDynamics\nusing Plots","category":"page"},{"location":"02-MovingSet/","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"First, define constants relative to this example.","category":"page"},{"location":"02-MovingSet/","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"F(x) = -[\n    2 * x[1] + 8/3 * x[2] - 34;\n    2 * x[2] + 5/4 * x[1] - 24.25\n]\n\nfunction project_moving_set!(sol, y; x = x, kwargs...)\n   λ = similar(sol)\n   Ain = [\n    1 0.9\n    0.9 1\n   ]\n   bin = [14.4; 14.1]\n   proj_success = NonSmoothDynamics.numerical_projection!(\n   λ, y, I,\n   spzeros(0, 2), ones(0),\n   Ain, bin,\n   zeros(2), Inf * ones(2))\n   sol .= λ\n   return proj_success\nend","category":"page"},{"location":"02-MovingSet/#Simulate-the-PDS","page":"Projected dynamical system with moving set","title":"Simulate the PDS","text":"","category":"section"},{"location":"02-MovingSet/","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"Run the reaction with 300 discretization point between 0 and 30.","category":"page"},{"location":"02-MovingSet/","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"x0 = zeros(2)\nt0, tf = 0.0, 30.0\nx_vals, t_vals, converged = NonSmoothDynamics.projected_dynamical_system(\n  x0, F, project_moving_set!,\n  t0, tf, 300\n)\n\nprintln(\"Final State: \", x_vals[:, end], \"is close to the known Nash equilibrium [5, 9]\")","category":"page"},{"location":"02-MovingSet/#Visualization","page":"Projected dynamical system with moving set","title":"Visualization","text":"","category":"section"},{"location":"02-MovingSet/","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"To plot the trajectory:","category":"page"},{"location":"02-MovingSet/","page":"Projected dynamical system with moving set","title":"Projected dynamical system with moving set","text":"using Plots\n\nplot(t_vals, x_vals[1, :], label=\"x₁ function of time\", color=:red, markersize=4)\nplot!(t_vals, x_vals[2, :], label=\"x₂ function of time\", color=:green, markersize=4)","category":"page"},{"location":"90-contributing/#contributing","page":"Contributing guidelines","title":"Contributing guidelines","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"First of all, thanks for the interest!","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"We welcome all kinds of contribution, including, but not limited to code, documentation, examples, configuration, issue creating, etc.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"Be polite and respectful, and follow the code of conduct.","category":"page"},{"location":"90-contributing/#Bug-reports-and-discussions","page":"Contributing guidelines","title":"Bug reports and discussions","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you think you found a bug, feel free to open an issue. Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.","category":"page"},{"location":"90-contributing/#Working-on-an-issue","page":"Contributing guidelines","title":"Working on an issue","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you found an issue that interests you, comment on that issue what your plans are. If the solution to the issue is clear, you can immediately create a pull request (see below). Otherwise, say what your proposed solution is and wait for a discussion around it.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"tip: Tip\nFeel free to ping us after a few days if there are no responses.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If your solution involves code (or something that requires running the package locally), check the developer documentation. Otherwise, you can use the GitHub interface directly to create your pull request.","category":"page"},{"location":"","page":"NonSmoothDynamics","title":"NonSmoothDynamics","text":"CurrentModule = NonSmoothDynamics","category":"page"},{"location":"#NonSmoothDynamics","page":"NonSmoothDynamics","title":"NonSmoothDynamics","text":"","category":"section"},{"location":"","page":"NonSmoothDynamics","title":"NonSmoothDynamics","text":"Documentation for NonSmoothDynamics.","category":"page"},{"location":"#Contributors","page":"NonSmoothDynamics","title":"Contributors","text":"","category":"section"},{"location":"","page":"NonSmoothDynamics","title":"NonSmoothDynamics","text":"<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->\n<!-- prettier-ignore-start -->\n<!-- markdownlint-disable -->\n<table>\n  <tbody>\n    <tr>\n      <td align=\"center\" valign=\"top\" width=\"14.28%\"><a href=\"http://tmigot.github.io\"><img src=\"https://avatars.githubusercontent.com/u/25304288?v=4?s=100\" width=\"100px;\" alt=\"Tangi Migot\"/><br /><sub><b>Tangi Migot</b></sub></a><br /><a href=\"#infra-tmigot\" title=\"Infrastructure (Hosting, Build-Tools, etc)\">🚇</a> <a href=\"#test-tmigot\" title=\"Tests\">⚠️</a> <a href=\"#code-tmigot\" title=\"Code\">💻</a></td>\n    </tr>\n  </tbody>\n</table>\n\n<!-- markdownlint-restore -->\n<!-- prettier-ignore-end -->\n\n<!-- ALL-CONTRIBUTORS-LIST:END -->","category":"page"}]
}
