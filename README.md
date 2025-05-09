# Piecewise

[![CI](https://github.com/ChristopheBerthod/Piecewise.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ChristopheBerthod/Piecewise.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ChristopheBerthod/Piecewise.jl/graph/badge.svg?token=cXaZZi9hdM)](https://codecov.io/gh/ChristopheBerthod/Piecewise.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ChristopheBerthod.github.io/Piecewise.jl/dev)

This repository contains three modules written in [Julia](https://julialang.org/):

- [Piecewise](https://christopheberthod.github.io/Piecewise.jl/dev/index.html) – Tools for defining piecewise functions made of user-defined elementary functions with parameters, called *formulas*. Properly configured formulas enable fast **integral transforms** of the piecewise function. The module defines seven formulas and a method for fitting a piecewise function with arbitrary formulas to a real function of a real variable.

- [PiecewiseHilbert](https://christopheberthod.github.io/Piecewise.jl/dev/hilbert.html) – Add methods to the formulas defined in [Piecewise](https://christopheberthod.github.io/Piecewise.jl/dev/index.html), enabling fast Hilbert transform of the piecewise functions that use these formulas.

- [PiecewiseLorentz](https://christopheberthod.github.io/Piecewise.jl/dev/lorentz.html) – Add methods to some of the formulas defined in [Piecewise](https://christopheberthod.github.io/Piecewise.jl/dev/index.html), enabling what we call a Lorentz transform of the piecewise functions.

### Dependencies

[HypergeometricFunctions](https://github.com/JuliaMath/HypergeometricFunctions.jl)&nbsp;&nbsp;|&nbsp;&nbsp;[LsqFit](https://github.com/JuliaNLSolvers/LsqFit.jl)&nbsp;&nbsp;|&nbsp;&nbsp;[PolyLog](https://github.com/Expander/PolyLog.jl)

### Installation

```julia
using Pkg; Pkg.add.(["Piecewise", "PiecewiseHilbert", "PiecewiseLorentz"]);
```

### Examples

- [Tutorial 1](https://github.com/ChristopheBerthod/Piecewise.jl/blob/main/notebooks/Tutorial-1.ipynb): Constructing approximations with [`Piecewise.piecewisefit`](https://christopheberthod.github.io/Piecewise.jl/dev/index.html#Piecewise.piecewisefit)
- [Tutorial 2](https://github.com/ChristopheBerthod/Piecewise.jl/blob/main/notebooks/Tutorial-2.ipynb): Solving an implicit equation using [`PiecewiseHilbert`](https://christopheberthod.github.io/Piecewise.jl/dev/hilbert.html)
- Use case: [MagnetoTransport.jl](https://christopheberthod.github.io/MagnetoTransport.jl/dev/index.html#MagnetoTransport.jl)

### Contributing

Contributions and suggestions are welcome. An expert advice on the following would be especially valuable:
- The evaluation of `PiecewiseFunction` objects seems slow and may benefit from optimizations.
- The parallelism of `piecewisefit` needs improvement.
