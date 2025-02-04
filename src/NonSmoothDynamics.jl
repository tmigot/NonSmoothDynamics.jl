module NonSmoothDynamics

using LinearAlgebra, SparseArrays

include("projection/project.jl")

using ADNLPModels, NLPModelsIpopt, Percival

include("projection/numerical_projection.jl")
include("projection/projection_plot.jl")
include("projection/projection_intersection.jl")

using Interpolations

include("interpolation.jl")

include("utils.jl")

using Logging, SolverCore

include("pds.jl")

end
