module PiecewiseHilbert

using PolyLog: li2
using HypergeometricFunctions: _₂F₁

using Piecewise

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
