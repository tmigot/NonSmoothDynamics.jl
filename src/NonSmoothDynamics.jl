module NonSmoothDynamics

using LinearAlgebra, SparseArrays

include("projection/project.jl")

using ADNLPModels, JSOSuite, SolverCore

include("projection/numerical_projection.jl")
include("projection/projection_plot.jl")
include("projection/projection_intersection.jl")

using Interpolations

include("interpolation.jl")

include("utils.jl")

using Logging

include("pds.jl")

end
