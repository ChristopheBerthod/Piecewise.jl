# PiecewiseLorentz

## Installation

The package won't work without [Piecewise](https://github.com/ChristopheBerthod/Piecewise.jl). Both can be installed with `Pkg.add`.
```julia
using Pkg; Pkg.add("Piecewise"); Pkg.add("PiecewiseLorentz")
```



## Introduction

The order-``m`` Lorentz transform of a function ``f(x)`` of a real variable ``x`` is defined as
```math
\left(L^m\circ f\right)(y, z) = \int_{-\infty}^{\infty}dx\,f(x)
\left[\frac{-\mathrm{Im}\,z/\pi}{(y-\mathrm{Re}\,z-x)^2+(\mathrm{Im}\,z)^2}\right]^m,
```
where ``m`` is a positive integer, ``y`` is a real number, and ``z`` is a complex number with strictly negative imaginary part. The name *Lorentz transform* was chosen, because the function in square brackets is a [Lorentzian](https://en.wikipedia.org/wiki/Cauchy_distribution) of half-width ``\mathrm{Im}\,z`` centered at ``x=y-\mathrm{Re}\,z``.

The module [PiecewiseLorentz](https://github.com/ChristopheBerthod/Piecewise.jl) adds to the [Formulas](@ref) defined in the module [Piecewise](@ref) the primitive functions corresponding to the order-``m`` Lorentz kernel ``\left[\frac{-\mathrm{Im}\,z/\pi}{(y-\mathrm{Re}\,z-x)^2+(\mathrm{Im}\,z)^2}\right]^m``, enabling the fast Lorentz transform of [`PiecewiseFunction`](@ref) objects using these formulas.

!!! warning
	The present version only implements ``m=1,2,3`` for the formulas [`POLY`](@ref), [`LOG`](@ref), [`PLS`](@ref), and [`XLOG`](@ref).
	Note that the ``m=1`` Lorentz transform is equivalent to a Hilbert transform:
	```math
	\left(L^1\circ f\right)(y, z) = -\frac{1}{\pi}\mathrm{Im}\,
	\left(H\circ f\right)(y-z).
	```

The method [`lorentz_transform`](@ref) return the Lorentz transform of a [`PiecewiseFunction`](@ref) object, as calculated using the definition above. The transform has algebraic tails at large ``y``, that may suffer from numerical errors. A numerically more stable [moment expansion](@ref momentexpansion-L) is therefore used for large values of ``|y-\mathrm{Re}\,z|``. Both the definition and the moment expansion are implemented in a [`LorentzTransform`](@ref) object created as:
```
L = LorentzTransform(f::PiecewiseFunction, m::Integer[, ground::Real])
```
The moment expansion is used if ``y-\mathrm{Re}\,z`` is outside the [support](@ref) of the piecewise function `f` and if the definition yields a result smaller than `ground` (by default `1e-10`). Once initialized, the Lorentz transform can be evaluated as `L(y, z)`.


## [Moment expansion](@id momentexpansion-L)

The expansion of the order-``m`` Lorentz transform for large ``|y-\mathrm{Re}\,z|`` is
```math
\begin{align*}
\left(L^m\circ f\right)(y, z)&=\sum_{n=2m}^{\infty}
\frac{1}{(y-\mathrm{Re}\,z)^n}
\int_{-\infty}^{\infty}dx\,f(x)P_n^{(m)}(x,\mathrm{Im}\,z)\\
P_n^{(m)}(x,a)&=\frac{1}{n!}\frac{d^n}{du^n}
\left[\frac{-a/\pi}{(1/u-x)^2+a^2}\right]^m_{u=0}\\
&=\left(\frac{a}{2\pi}\right)^m\frac{2}{(m-1)!}
\sum_{k=0}^{n-2m}\cos\left((n-k)\frac{\pi}{2}\right)\\
&\quad\times\frac{(n-1)!}{k!(n-k-2m)!!(n-k-1)!!}a^{n-k-2m}x^k.
\end{align*}
```
Dropping the terms that vanish because of the cosine, this can be recast as
```math
\begin{align*}
\left(L^m\circ f\right)(y, z)&=\left(\frac{-\mathrm{Im}\,z}{2\pi}\right)^m
\frac{2}{(m-1)!}\frac{1}{(y-\mathrm{Re}\,z)^{2m}}
\sum_{n=0}^{\infty}\frac{1}{(y-\mathrm{Re}\,z)^n}
\sum_{k=0}^{n\doteqdot 2}C_{nk}^{(m)}(\mathrm{Im}\,z)^{2k}\\
C_{nk}^{(m)}&=(-1)^k\frac{(2m+n-1)!}{(n-2k)!(2k)!!(2m+2k-1)!!}(M\circ f)(n-2k),
\end{align*}
```
where ``n\doteqdot 2`` means the integer division of ``n`` by ``2``. The structure [`LorentzTransform`](@ref) stores the piecewise function and the expansion coefficients ``C_{nk}^{(m)}``. If ``y-\mathrm{Re}\,z`` is outside the [support](@ref) of the function ``f(x)`` and if the definition of the Lorentz transform yields a result smaller than `ground`, the moment expansion is used.



## Primitives

Here, we describe the primitives added by [PiecewiseLorentz](https://github.com/ChristopheBerthod/Piecewise.jl) to the [Formulas](@ref) provided by [Piecewise](@ref). In order to handle even and odd piecewise functions, we define the primitive of ``F(x,\mathbf{a})\left[\frac{-\mathrm{Im}\,z/\pi}{(y-\mathrm{Re}\,z-sx)^2+(\mathrm{Im}\,z)^2}\right]^m`` as ``\mathcal{F}_m(s,x,\mathbf{a},y,z)`` with ``s=\pm1``. 

[`POLY`](@ref POLY-L) | 
[`LOG`](@ref LOG-L) | 
[`PLS`](@ref PLS-L)
[`XLOG`](@ref XLOG-L)

---

##### [`POLY`](@id POLY-L)

The function `POLY.value(x, a)` is
```math
F(x,\mathbf{a}) = \sum_{i=1}^k a_i x^{i-1}.
```
The primitive function for the Lorentz kernel of order ``m`` is
```math
\begin{equation*}
\mathcal{F}_m(s,x,\mathbf{a},y,z) = \left[\frac{-\mathrm{Im}\,z/\pi}{(y-\mathrm{Re}\,z)^2+(\mathrm{Im}\,z)^2}
\right]^m\sum_{i=1}^k\frac{a_ix^i}{i}
F_1\left(i,m,m,i+1,\frac{sx}{y-z},\frac{sx}{y-z^*}\right),
\end{equation*}
```
where ``F_1`` is the [Appell function](https://en.wikipedia.org/wiki/Appell_series).

Since the Appell function isn't available yet in pure `Julia` [^1], we use alternate forms for ``m=1``, ``2``, and ``3``, that involve the hypergeometric function ``{_2F_1}(a,b,c,z)``:
```math
\begin{align*}
\mathcal{F}_1(s,x,\mathbf{a},y,z) &=
\frac{-\mathrm{Im}\,z/\pi}{(y-\mathrm{Re}\,z)^2+(\mathrm{Im}\,z)^2}
\sum_{i=1}^k\frac{a_ix^i}{i}\mathrm{Re}\left\{\frac{y-z^*}{z-z^*}
\left[2\,{_2F_1}\left(1,i,i+1,\frac{sx}{y-z}\right)\right]\right\}\\
\mathcal{F}_2(s,x,\mathbf{a},y,z) &=
\frac{1}{2\pi^2[(y-\mathrm{Re}\,z)^2+(\mathrm{Im}\,z)^2]}
\sum_{i=1}^k\frac{a_ix^i}{i}\mathrm{Re}\left\{\frac{y-z^*}{z-z^*}
\left[\phantom{\frac{1}{1}}\right.\right.\\
&\quad\left.\left.\left(2+(i-1)\frac{z-z^*}{y-z}\right)
{_2F_1}\left(1,i,i+1,\frac{sx}{y-z}\right)
-i\frac{z-z^*}{y-z-sx}\right]\right\}\\
\mathcal{F}_3(s,x,\mathbf{a},y,z) &=
-\frac{1}{16\pi^3\,\mathrm{Im}\,z[(y-\mathrm{Re}\,z)^2+(\mathrm{Im}\,z)^2]}
\sum_{i=1}^k\frac{a_ix^i}{i}\mathrm{Re}\left\{\frac{y-z^*}{y-z}
\left[\phantom{\frac{1}{1}}\right.\right.\\
&\quad\left(12+6(i-1)\frac{z-z^*}{y-z}
+(i-1)(i-2)\left(\frac{z-z^*}{y-z}\right)^2\right)
{_2F_1}\left(1,i,i+1,\frac{sx}{y-z}\right)\\
&\quad\left.\left.-6i\frac{z-z^*}{y-z-sx}+i\left(1-(i-2)
\frac{y-z-sx}{y-z}\right)\left(\frac{z-z^*}{y-z-sx}\right)^2\right]\right\}.
\end{align*}
```
``{_2F_1}(a,b,c,z)`` is continuous for ``z\in\mathbb{C}``, except on the real axis ``z=r\in\mathbb{R}`` for ``r>1``. Since ``sx/(y-z)`` crosses the real axis at ``x=0``, the primitives are continuous for ``x\in\mathbb{R}``.

[^1]: It would be available through `sympy.functions.special.hyper.appellf1`.

---

##### [`LOG`](@id LOG-L)

The function `LOG.value(x, a)` is
```math
F(x,\mathbf{a}) = a_2\ln|x-a_1|.
```
The primitive function is unknown for the Lorentz kernel of general order. For the orders ``m=1``, ``2``, and ``3``, they are
```math
\begin{align*}
\mathcal{F}_1(s,x,\mathbf{a},y,z) &= \frac{sa_2}{\pi}\,\mathrm{Im}\left\{
\ln|x-a_1|\ln\left(\frac{y-z-sx}{y-z-sa_1}\right)
+\mathrm{Li}_2\left(\frac{s(x-a_1)}{y-z-sa_1}\right)\right\}\\
\mathcal{F}_2(s,x,\mathbf{a},y,z) &= -\frac{sa_2}{2\pi^2\,\mathrm{Im}\,z}
\,\mathrm{Im}\left\{\ln|x-a_1|\ln\left(\frac{y-z-sx}{y-z-sa_1}\right)
+\mathrm{Li}_2\left(\frac{s(x-a_1)}{y-z-sa_1}\right)\right.\\
&\quad\left.+i\,\mathrm{Im}\,z\left[\frac{\ln|x-a_1|}{y-z-sx}
+\frac{\ln(y-z-sx)-\ln|x-a_1|}{y-z-sa_1}\right]\right\}\\
\mathcal{F}_3(s,x,\mathbf{a},y,z) &= \frac{3sa_2}{8\pi^3(\mathrm{Im}\,z)^2}
\,\mathrm{Im}\left\{\ln|x-a_1|\ln\left(\frac{y-z-sx}{y-z-sa_1}\right)
+\mathrm{Li}_2\left(\frac{s(x-a_1)}{y-z-sa_1}\right)\right.\\
&\quad+i\,\mathrm{Im}\,z\left[\frac{\ln|x-a_1|}{y-z-sx}
+\frac{\ln(y-z-sx)-\ln|x-a_1|}{y-z-sa_1}\right]
+\frac{(\mathrm{Im}\,z)^2}{3}\left[\frac{\ln|x-a_1|}{(y-z-sx)^2}\right.\\
&\quad\left.\left.+\frac{\ln(y-z-sx)-\ln|x-a_1|}{(y-z-sa_1)^2}
-\frac{1}{(y-z-sa_1)(y-z-sx)}\right]\right\},
\end{align*}
```
where ``\mathrm{Li}_2(z)`` is the [polylogarithm](https://en.wikipedia.org/wiki/Polylogarithm). Both functions ``\ln(z)`` and ``\mathrm{Li}_2(z)`` are continuous for ``z\in\mathbb{C}``, except on the real axis, where the imaginary part has a branch cut. Since both ``(y-z-sx)/(y-z-sa_1)`` and ``s(x-a_1)/(y-z-a_1)`` cross the real axis at ``x=a_1``, the primitives are continuous in the domains where they can be used.

---

##### [`PLS`](@id PLS-L)

The function `PLS.value(x, a)` is
```math
F(x,\mathbf{a}) = a_3|x-a_1|^{a_2}.
```
The primitive function for the Lorentz kernel of order ``m`` is
```math
\begin{align*}
\mathcal{F}_m(s,x,\mathbf{a},y,z) &= \left[\frac{-\mathrm{Im}\,z/\pi}{(y-\mathrm{Re}\,z-sa_1)^2+(\mathrm{Im}\,z)^2}
\right]^ma_3|x-a_1|^{a_2}\frac{x-a_1}{a_2+1}\\
&\quad\times F_1\left(a_2+1,m,m,a_2+2,\frac{s(x-a_1)}{y-z-sa_1},
\frac{s(x-a_1)}{y-z^*-sa_1}\right),
\end{align*}
```
where ``F_1`` is the [Appell function](https://en.wikipedia.org/wiki/Appell_series).

Since the Appell function isn't available yet in pure `Julia` [^1], we use alternate forms for ``m=1``, ``2``, and ``3``, that involve the hypergeometric function ``{_2F_1}(a,b,c,z)``:
```math
\begin{align*}
\mathcal{F}_1(s,x,\mathbf{a},y,z) &= 
\frac{-\mathrm{Im}\,z/\pi}{(y-\mathrm{Re}\,z-sa_1)^2+(\mathrm{Im}\,z)^2}
a_3|x-a_1|^{a_2}\frac{x-a_1}{a_2+1}
\,\mathrm{Re}\left\{\phantom{\frac{1}{1}}\right.\\
&\quad\left.\frac{y-z^*-sa_1}{z-z^*}\left[
2\,{_2F_1}\left(1,a_2+1,a_2+2,\frac{s(x-a_1)}{y-z-sa_1}\right)\right]\right\}\\
\mathcal{F}_2(s,x,\mathbf{a},y,z) &=
\frac{1}{2\pi^2[(y-\mathrm{Re}\,z-sa_1)^2+(\mathrm{Im}\,z)^2]}
a_3|x-a_1|^{a_2}\frac{x-a_1}{a_2+1}
\,\mathrm{Re}\left\{\phantom{\frac{1}{1}}\right.\\
&\quad\frac{y-z^*-sa_1}{z-z^*}\left[\left(2+a_2\frac{z-z^*}{y-z-sa_1}\right)
{_2F_1}\left(1,a_2+1,a_2+2,\frac{s(x-a_1)}{y-z-sa_1}\right)\right.\\
&\quad\left.\left.-(a_2+1)\frac{z-z^*}{y-z-sx}\right]\right\}\\
\mathcal{F}_3(s,x,\mathbf{a},y,z) &=
-\frac{1}{16\pi^3\,\mathrm{Im}\,z[(y-\mathrm{Re}\,z-sa_1)^2+(\mathrm{Im}\,z)^2]}
a_3|x-a_1|^{a_2}\frac{x-a_1}{a_2+1}
\,\mathrm{Re}\left\{\phantom{\frac{1}{1}}\right.\\
&\quad\frac{y-z^*-sa_1}{z-z^*}\left[\left(12+6a_2\frac{z-z^*}{y-z-sa_1}
+a_2(a_2-1)\left(\frac{z-z^*}{y-z-sa_1}\right)^2\right)\right.\\
&\quad\times{_2F_1}\left(1,a_2+1,a_2+2,\frac{s(x-a_1)}{y-z-sa_1}\right)
-6(a_2+1)\frac{z-z^*}{y-z-sx}\\
&\quad\left.\left.+(a_2+1)\left(1-(a_2-1)
\frac{y-z-sx}{y-z-sa_1}\right)
\left(\frac{z-z^*}{y-z-sx}\right)^2\right]\right\}.
\end{align*}
```
``{_2F_1}(a,b,c,z)`` is continuous for ``z\in\mathbb{C}``, except on the real axis, where the imaginary part has a branch cut. Since ``s(x-a_1)/(y-z-sa_1)`` crosses the real axis at ``x=a_1``, the primitives are continuous in the domain, which must exclude ``a_1``.

---

##### [`XLOG`](@id XLOG-L)

The function `XLOG.value(x, a)` is
```math
F(x,\mathbf{a}) = a_2x\ln|x-a_1|.
```
The primitive function is unknown for the Lorentz kernel of general order. For the orders ``m=1``, ``2``, and ``3``, they are
```math
\begin{align*}
\mathcal{F}_1(s,x,\mathbf{a},y,z) &= \frac{a_2}{\pi}\,\mathrm{Im}\left\{
(y-z)\left[\ln|x-a_1|\ln\left(\frac{y-z-sx}{y-z-sa_1}\right)
+\mathrm{Li}_2\left(\frac{s(x-a_1)}{y-z-sa_1}\right)\right]\right\}\\
\mathcal{F}_2(s,x,\mathbf{a},y,z) &= -\frac{a_2}{2\pi^2\,\mathrm{Im}\,z}
\,\mathrm{Im}\left\{(y-\mathrm{Re}\,z)\left[
\ln|x-a_1|\ln\left(\frac{y-z-sx}{y-z-sa_1}\right)
+\mathrm{Li}_2\left(\frac{s(x-a_1)}{y-z-sa_1}\right)\right]\right.\\
&\quad\left.+i\,(y-z)\mathrm{Im}\,z\left[\frac{\ln|x-a_1|}{y-z-sx}
+\frac{\ln(y-z-sx)-\ln|x-a_1|}{y-z-sa_1}\right]\right\}\\
\mathcal{F}_3(s,x,\mathbf{a},y,z) &= \frac{3a_2}{8\pi^3(\mathrm{Im}\,z)^2}
\,\mathrm{Im}\left\{(y-\mathrm{Re}\,z)\left[
\ln|x-a_1|\ln\left(\frac{y-z-sx}{y-z-sa_1}\right)
+\mathrm{Li}_2\left(\frac{s(x-a_1)}{y-z-sa_1}\right)\right]\right.\\
&\quad+i\,\left(y-z+\frac{2i}{3}\mathrm{Im}\,z\right)
\mathrm{Im}\,z\left[\frac{\ln|x-a_1|}{y-z-sx}
+\frac{\ln(y-z-sx)-\ln|x-a_1|}{y-z-sa_1}\right]\\
&\quad+\frac{(\mathrm{Im}\,z)^2}{3}(y-z)\left[\frac{\ln|x-a_1|}{(y-z-sx)^2}
+\frac{\ln(y-z-sx)-\ln|x-a_1|}{(y-z-sa_1)^2}\right.\\
&\quad\left.\left.-\frac{1}{(y-z-sa_1)(y-z-sx)}\right]\right\},
\end{align*}
```
where ``\mathrm{Li}_2(z)`` is the [polylogarithm](https://en.wikipedia.org/wiki/Polylogarithm). Both functions ``\ln(z)`` and ``\mathrm{Li}_2(z)`` are continuous for ``z\in\mathbb{C}``, except on the real axis, where the imaginary part has a branch cut. Since both ``(y-z-sx)/(y-z-sa_1)`` and ``s(x-a_1)/(y-z-a_1)`` cross the real axis at ``x=a_1``, the primitives are continuous in the domains where they can be used.



## Public interface

### Index

```@index
Pages   = ["lorentz.md"]
Modules = [PiecewiseLorentz]
Order   = [:type, :function]
```


### Type

```@docs
LorentzTransform
```


### Methods

```@docs
lorentz_transform
```
