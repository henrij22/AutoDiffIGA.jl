module AutoDiffIGA

using ForwardDiff
using NURBS
using StatsBase
using FastGaussQuadrature
using StaticArrays

include("geometry/geometryhelper1d.jl")
include("geometry/geometryhelper2d.jl")

export makeGeometry1d, makeGeometry2d

include("basis/bspline.jl")
include("basis/tensorbspline.jl")

export BSplineBasis, BSplineTensorBasis2d

include("autodiff/bernoulli.jl")
include("autodiff/linearElastic.jl")
include("autodiff/timoshenko.jl")

include("autodiff/quadrature.jl")

export LinearElastic, TimoshenkoBeam, BernoulliBeam
export computeK, computeF

end
