# Ultimately, this should provide the L^m transform for all Formula and all m.
# For the moment, only m=1,2,3 are provided, and only for POLY, LOG, PLS, and XLOG.

module PiecewiseLorentz

using PolyLog: li2
using HypergeometricFunctions: _₂F₁

using Piecewise

import Base.show, Base.print, Base.broadcastable

export
    # Types
    LorentzTransform,

    # Methods
    lorentz_transform

include("types.jl")
include("private-methods.jl")
include("public-methods.jl")

end
