module NonSmoothDynamics

using LinearAlgebra

include("projection/project.jl")
include("projection/projection_plot.jl")
include("projection/projection_intersection.jl")

using Interpolations

include("interpolation.jl")

include("pds.jl")

end
