# Piecewise

[![CI](https://github.com/ChristopheBerthod/Piecewise.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ChristopheBerthod/Piecewise.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ChristopheBerthod/Piecewise.jl/graph/badge.svg?token=cXaZZi9hdM)](https://codecov.io/gh/ChristopheBerthod/Piecewise.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ChristopheBerthod.github.io/Piecewise.jl/dev)
[![status](https://joss.theoj.org/papers/c54fbe7bd14b9ef757580b362235db46/status.svg)](https://joss.theoj.org/papers/c54fbe7bd14b9ef757580b362235db46)

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
- Use case: [MagnetoTransport.jl](https://github.com/ChristopheBerthod/MagnetoTransport.jl)

### Contributing

Contributions and suggestions are welcome. If you're interested in contributing, please feel free to reach out via [email](mailto:christophe.berthod@unige.ch). Bug reports can be submitted through the [issue tracker](https://github.com/ChristopheBerthod/Piecewise.jl/issues), but email is preferred for quicker responses.

This policy may evolve if additional contributors become actively involved.

Expert input would be especially valuable in the following areas:

* **Evaluation performance**: The evaluation of `PiecewiseFunction` objects appears to be slower than desirable. More efficient strategies to identify the relevant domain could yield significant improvements.

* **Parallel execution**: The parallelism in `piecewisefit` could benefit from the scrutiny of those with a deeper understanding of Julia's parallel execution.