module NonSmoothDynamics

using LinearAlgebra

include("projection/project.jl")

using ADNLSModel, NLPModelsIpopt

include("projection/numerical_projection.jl")
include("projection/projection_plot.jl")
include("projection/projection_intersection.jl")

using Interpolations

include("interpolation.jl")

include("pds.jl")

end
