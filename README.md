# Piecewise <a href='https://ChristopheBerthod.github.io/Piecewise.jl/dev'><img src="docs/src/assets/logo.png" align="right" height="138.5" /></a>

[![CI](https://github.com/ChristopheBerthod/Piecewise.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ChristopheBerthod/Piecewise.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ChristopheBerthod/Piecewise.jl/graph/badge.svg?token=cXaZZi9hdM)](https://codecov.io/gh/ChristopheBerthod/Piecewise.jl)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/ChristopheBerthod/Piecewise.jl/LICENSE)
[![Documentation](https://github.com/ChristopheBerthod/Piecewise.jl/actions/workflows/Documenter.yml/badge.svg)](https://github.com/ChristopheBerthod/Piecewise.jl/actions/workflows/Documenter.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ChristopheBerthod.github.io/Piecewise.jl/dev)
<!--[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ChristopheBerthod.github.io/Piecewise.jl/stable)-->

The [Julia](https://julialang.org/) module **Piecewise** provides tools for defining a piecewise function made of arbitrary user-defined elementary functions with parameters. Properly configured elementary functions enable fast **integral transforms** of the piecewise function. The module defines seven elementary functions and a method for fitting a piecewise function to a real function of a real variable.

The module **PiecewiseHilbert** adds methods to the elementary functions defined in **Piecewise**, enabling a fast Hilbert transform of the piecewise functions that use these elementary functions.

The module **PiecewiseLorentz** adds methods to some of the elementary functions defined in **Piecewise**, enabling what we call a Lorentz transform of the piecewise functions.
