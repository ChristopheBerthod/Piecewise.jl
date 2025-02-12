module Piecewise

using Printf
using HypergeometricFunctions: _₂F₁
using LsqFit: curve_fit
using Base.Threads: @threads, @spawn, fetch

import Base.show, Base.print, Printf.format
import Base.broadcastable, Base.+, Base.-, Base.*, Base./, Base.sum

export

    # Types
    Formula,
    Piece,
    PiecewiseFunction,

    # Constants
    POLY,
    TAIL,
    LOG,
    ISRS,
    PLS,
    XLOG,
    XISRS,

    # Methods
    domains,
    intervals,
    support,
    singularities,
    formulas,
    integraltransform,
    moment,
    piecewisefit,
    format

include("types.jl")
include("formulas.jl")
include("private-methods.jl")
include("public-methods.jl")

end
