module PiecewiseHilbert

using HypergeometricFunctions: _₂F₁
using Piecewise
using PolyLog: li2

import Base.show, Base.print, Base.broadcastable

export
    # Types
    HilbertTransform,

    # Methods
    hilbert_transform

include("types.jl")
include("private-methods.jl")
include("public-methods.jl")

end
