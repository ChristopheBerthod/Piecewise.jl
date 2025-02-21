# PiecewiseHilbert

## Installation

The package won't work without [Piecewise](https://github.com/ChristopheBerthod/Piecewise.jl). Both can be installed with `Pkg.add`.
```julia
using Pkg; Pkg.add("Piecewise"); Pkg.add("PiecewiseHilbert")
```



## Introduction

The Hilbert transform of a function ``f(x)`` of a real variable ``x`` is defined as
```math
(H\circ f)(z) = \int_{-\infty}^{\infty}dx\,\frac{f(x)}{z-x},
```
where ``z\in\mathbb{C}\setminus\mathbb{R}`` is a complex number with non-zero imaginary part. The transform is well-defined if ``f(x)`` has no non-integrable singularity and if it vanishes at infinity faster than ``|x|^{-\epsilon}`` with ``\epsilon>0`` [^1]. In the limiting case ``\epsilon=0``, the transform is well-defined only if ``f(+\infty)=f(-\infty)``.

[^1]: We have ``\int^x |x|^{-\epsilon}/(z-x)=(\mathrm{sign}(x)/z)^{\epsilon}B_{x/z}(1-\epsilon,0)``, where ``B_z(a, b)`` is the [incomplete beta function](https://en.wikipedia.org/wiki/Beta_function), which approaches ``-\pi (\mp1/z)^{\epsilon}/\sin(\pi\epsilon)`` for ``x\to\pm\infty``.

The module [PiecewiseHilbert](https://github.com/ChristopheBerthod/Piecewise.jl) adds to the [Formulas](@ref) defined in the module [Piecewise](@ref) the primitive functions corresponding to the Hilbert kernel ``1/(z-x)``, enabling the fast Hilbert transform of piecewise functions using these formulas.

The method [`hilbert_transform`](@ref) returns the Hilbert transform of a [`PiecewiseFunction`](@ref) object, as calculated using the definition above. It also works with user-defined [`Formula`](@ref) objects. Depending on the formulas used, however, the result can be numerically unstable at large ``|z|``. This is circumvented by means of a [moment expansion](@ref momentexpansion-H). Both the definition and the moment expansion are implemented in a [`HilbertTransform`](@ref) object created as:
```
H = HilbertTransform(f::PiecewiseFunction[, radius::Real])
```
The moment expansion is used if ``|z|>`` `radius` (set automatically by default). Once initialized, the Hilbert transform can be evaluated as `H(z)`.



## [Moment expansion](@id momentexpansion-H)

The large-``|z|`` expansion of the Hilbert transform is
```math
(H\circ f)(z) = \sum_{n=0}^{\infty}\frac{(M\circ f)(n)}{z^{n+1}},\qquad
(M\circ f)(n) = \int_{-\infty}^{\infty}dx\,f(x)x^n.
```
The structure [`HilbertTransform`](@ref) stores the piecewise function and its moments, as provided by the method [`moment`](@ref) of [Piecewise](@ref), such that a numerically stable result is produced using the definition of the Hilbert transform within a given disk in the complex plane and the moment expansion beyond that disk.



## Primitives

Here, we describe the primitives added by [PiecewiseHilbert](https://github.com/ChristopheBerthod/Piecewise.jl) to the [Formulas](@ref) provided by [Piecewise](@ref).

[`POLY`](@ref POLY-H) | 
[`TAIL`](@ref TAIL-H) | 
[`LOG`](@ref LOG-H) |
[`ISRS`](@ref ISRS-H) |
[`PLS`](@ref PLS-H) |
[`XLOG`](@ref XLOG-H) |
[`XISRS`](@ref XISRS-H)

---

##### [`POLY`](@id POLY-H)

The function `POLY.value(x, a)` is
```math
F(x,\mathbf{a}) = \sum_{i=1}^k a_i x^{i-1}.
```
A primitive function for the Hilbert kernel is
```math
\mathcal{F}(x,\mathbf{a},z) = -\ln(x-z)\sum_{i=1}^k a_i z^{i-1}
-\sum_{i=0}^{k-2}z^i\sum_{j=i}^{k-2}a_{j+2}\frac{x^{j-i+1}}{j-i+1}.
```
This primitive is interesting, because it involves only elementary functions. Unfortunately, it suffers from numerical instability at large ``z``. We use another primitive that is stable at large ``z``, although slightly slower as it requires the [hypergeometric function](https://en.wikipedia.org/wiki/Hypergeometric_function) ``{_2F_1}(a,b,c,z)``:
```math
\mathcal{F}(x,\mathbf{a},z) = \frac{1}{z}\sum_{i=1}^k\frac{a_ix^i}{i}
{_2F_1}\left(1,i,i+1,\frac{x}{z}\right).
```
This primitive is continuous for ``x\in\mathbb{R}``.

---

##### [`TAIL`](@id TAIL-H)

The function `TAIL.value(x, a)` is
```math
F(x,\mathbf{a}) = \frac{a_1+a_2x}{a_3+a_4x+a_5x^2}.
```
The primitive function for the Hilbert kernel is
```math
\begin{align*}
\mathcal{F}(x,\mathbf{a},z) &= -\frac{1}{a_3+a_4z+a_5z^2}\left\{
(a_1+a_2z)\left[\ln(z-x)-\ln\sqrt{a_3+a_4x+a_5x^2}\right]
\phantom{\frac{\frac{1_1^1}{1_1^1}}{1_1^1}}\right.
\\&\quad\left. -\left[2a_2a_3-a_1a_4+(a_2a_4-2a_1a_5)z\right]
\frac{\tanh^{-1}\left(\frac{a_4+2a_5x}{\Delta}\right)}{\Delta}\right\},
\end{align*}
```
where ``\Delta=\sqrt{a_4^2-4a_3a_5}``. This primitive is obviously continuous for ``x\in\mathbb{R}`` if ``a_4^2-4a_3a_5<0``, because in that case the zeros ``x_{\pm}`` of ``a_3+a_4x+a_5x^2`` are not on the real axis, such that ``a_3+a_4x+a_5x^2`` does not change sign and, therefore, has the sign of ``a_3``. If ``a_4^2-4a_3a_5>0``, both ``\ln\sqrt{\cdots}`` and ``\tanh^{-1}(\cdots)`` have a discontinuous imaginary part at ``x=x_{\pm}``. However, since none of these zeros lies inside the domain, the function is continuous in the domain.

---

##### [`LOG`](@id LOG-H)

The function `LOG.value(x, a)` is
```math
F(x,\mathbf{a}) = a_2\ln|x-a_1|.
```
The primitive function for the Hilbert kernel is
```math
\mathcal{F}(x,\mathbf{a},z) = -a_2\left[\ln|x-a_1|
\ln\left(\frac{z-x}{z-a_1}\right)
+\mathrm{Li}_2\left(\frac{x-a_1}{z-a_1}\right)\right],
```
where ``\mathrm{Li}_2(z)`` is the [polylogarithm](https://en.wikipedia.org/wiki/Polylogarithm). Both functions ``\ln(z)`` and ``\mathrm{Li}_2(z)`` are continuous for ``z\in\mathbb{C}``, except on the real axis, where the imaginary part has a branch cut. Since both ``(z-x)/(z-a_1)`` and ``(x-a_1)/(z-a_1)`` cross the real axis at ``x=a_1``, the primitive is continuous in the domains where it can be used.

---

##### [`ISRS`](@id ISRS-H)

The function `ISRS.value(x, a)` is
```math
F(x,\mathbf{a}) = \frac{a_2}{\sqrt{|x^2-a_1^2|}}.
```
The primitive function for the Hilbert kernel is
```math
\mathcal{F}(x,\mathbf{a},z) = -\frac{a_2}{\sqrt{|x^2-a_1^2|}}
\frac{\sqrt{x^2-a_1^2}}{\sqrt{z^2-a_1^2}}\tanh^{-1}
\left(\frac{a_1^2-xz}{\sqrt{x^2-a_1^2}\sqrt{z^2-a_1^2}}\right).
```
This primitive is discontinuous at ``x=\pm a_1`` but continuous everywhere else, so it is continuous inside the domains where it can be used.

---

##### [`PLS`](@id PLS-H)

The function `PLS.value(x, a)` is
```math
F(x,\mathbf{a}) = a_3|x-a_1|^{a_2}.
```
The primitive function for the Hilbert kernel is
```math
\mathcal{F}(x,\mathbf{a},z) = \frac{a_3|x-a_1|^{a_2}}{1+a_2}
\frac{x-a_1}{z-a_1}{_2F_1}\left(1,a_2+1,a_2+2,\frac{x-a_1}{z-a_1}\right),
```
where ``{_2F_1}(a,b,c,z)`` is the [hypergeometric function](https://en.wikipedia.org/wiki/Hypergeometric_function). ``{_2F_1}(a,b,c,z)`` is continuous for ``z\in\mathbb{C}``, except on the real axis, where the imaginary part has a branch cut. Since ``(x-a_1)/(z-a_1)`` crosses the real axis at ``x=a_1``, the primitive is continuous in the domain, which must exclude ``a_1``.

---

##### [`XLOG`](@id XLOG-H)

The function `XLOG.value(x, a)` is
```math
F(x,\mathbf{a}) = a_2x\ln|x-a_1|.
```
The primitive function for the Hilbert kernel is
```math
\mathcal{F}(x,\mathbf{a},z) = a_2\left\{x-\ln|x-a_1|\left[x-a_1
+z\ln\left(\frac{z-x}{z-a_1}\right)\right]
-z\,\mathrm{Li}_2\left(\frac{x-a_1}{z-a_1}\right)\right\},
```
where ``\mathrm{Li}_2(z)`` is the [polylogarithm](https://en.wikipedia.org/wiki/Polylogarithm). Both functions ``\ln(z)`` and ``\mathrm{Li}_2(z)`` are continuous for ``z\in\mathbb{C}``, except on the real axis, where the imaginary part has a branch cut. Since both ``(z-x)/(z-a_1)`` and ``(x-a_1)/(z-a_1)`` cross the real axis at ``x=a_1``, the primitive is continuous in the domains where it can be used.

---
 
##### [`XISRS`](@id XISRS-H)

The function `XISRS.value(x, a)` is
```math
F(x,\mathbf{a}) = \frac{a_2x}{\sqrt{|x^2-a_1^2|}}.
```
The primitive function for the Hilbert kernel is
```math
\mathcal{F}(x,\mathbf{a},z) = -\frac{a_2\sqrt{x^2-a_1^2}}{\sqrt{|x^2-a_1^2|}}
\left[\tanh^{-1}\left(\frac{x}{\sqrt{x^2-a_1^2}}\right)
+\frac{z}{\sqrt{z^2-a_1^2}}\tanh^{-1}\left(\frac{a_1^2-xz}
{\sqrt{x^2-a_1^2}\sqrt{z^2-a_1^2}}\right)\right].
```
This primitive is discontinuous at ``x=\pm a_1`` but continuous everywhere else, so it is continuous inside the domains where it can be used.



## Example

The following example produces an image similar to the [Piecewise](@ref) logo.

```@example
ENV["GKSwstype"] = "100" # hide
using ..Piecewise, ..PiecewiseHilbert # hide
f = PiecewiseFunction([
    Piece((-1, -6/10), (true, false), LOG, [-6/10, -1]),
    Piece((-6/10, -3/10), (false, true), [POLY, LOG], [[2, 10/3], [-6/10, -1]]),
    Piece((-3/10, 0), (false, true), [POLY, LOG, XLOG],
        [[1+log(10/3)], [-3/10, 1+1/log(10/3)], [-3/10, 10/3*(1+1/log(10/3))]]),
    Piece((2/10, 5/10), (false, false), [XISRS, LOG], [[2/10, 1], [5/10, -1]]),
    Piece((5/10, 1), (false, false), [POLY, LOG, PLS],
        [[log(1/2)], [5/10, -1], [1, 1/4, 1]])
])

H = HilbertTransform(f)

using Plots
z1 = -2:0.001:1.5
plot(xlim=(-1.1, 1.1), ylim=(0, 6), axis=false, grid=false, legend=:none)
for z2 in reverse(0.01:0.015:0.6)
    plot!(z1 .+ z2*4/3, -imag.(H.(z1 .+ im*z2^2))/π .+ z2*20/3,
        color=:grey, fill=(0, :white))
end
plot!(z1, f.(z1), linewidth = 5, color = :red, fill=(0, :pink))
savefig("logo.svg"); nothing # hide
```

![](logo.svg)



## Public interface

### Index

```@index
Pages   = ["hilbert.md"]
Modules = [PiecewiseHilbert]
Order   = [:type, :function]
```


### Type

```@docs
HilbertTransform
```


### Methods

```@docs
hilbert_transform
```
