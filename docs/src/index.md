# Piecewise

A piecewise function of a real variable ``x`` returns a value computed with a formula that depends on the interval in which ``x`` lies. The [Julia](https://julialang.org) module [Piecewise](@ref) represents such a function as a collection of *pieces*. Each piece is an object of type [`Piece`](@ref) that contains an interval and a rule to compute a value, given ``x``. The rule is expressed by *formulas*, that are contained in objects of type [`Formula`](@ref).

In mathematical terms, a piecewise function can be written as
```math
	f(x) = \sum_i
	\underbrace{\theta\left(x-x_i^{\min}\right)\theta\left(x_i^{\max}-x\right)
	\underbrace{\sum_j
	\underbrace{F_{ij}(x, \mathbf{a}_{ij})}
	_{\text{formula}}}
	_{\text{rule}}}
	_{\text{piece}},
```
where ``\theta(x)`` is the [Heaviside step function](https://en.wikipedia.org/wiki/Heaviside_step_function). The functions ``F(x, \mathbf{a})`` are the *formulas*, that take a set of parameters ``\mathbf{a}``. The [`Formula`](@ref) object holds the function ``F`` together with constraints regarding its applicability with the parameters ``\mathbf{a}`` in a given interval. The quantity ``\sum_jF_{ij}(x, \mathbf{a}_{ij})`` is the rule that computes the value, given ``x`` [^1]. The parameters ``\mathbf{a}`` are in general different in each piece, but the function ``F`` may be the same. For instance, the cubic-spline interpolation of a dataset could be represented as a [`PiecewiseFunction`](@ref) object in which all formulas are 3rd-order polynomials ``F(x, \mathbf{a}_i) = a_{i1}+a_{i2}x+a_{i3}x^2+a_{i4}x^3``. The *piece* contains the rule and the interval defined by ``x_i^{\min}`` and ``x_i^{\max}``.

[^1]: The use of several formulas in a piece is therefore limited to function types for which a `+` method exists.



## Typical use

##### Encoding computer-intensive functions

Consider a function ``f(x)`` that is computationally intensive. If formulas ``F(x, \mathbf{a})`` exist, that can approximate ``f(x)`` in restricted domains to a high accuracy, then ``f(x)`` can be encoded in a fast [`PiecewiseFunction`](@ref) object. This is useful if the function ``f(x)`` is just one part of a bigger calculation that requires computing ``f(x)`` many times. The use of arbitrary functions as formulas allows one to represent behaviors, such as power laws or neighborhood of singularities, that are not well described by polynomial interpolations.

!!! tip
	The module [Piecewise](@ref) provides the method [`piecewisefit`](@ref piecewisefitsection) for constructing a piecewise approximation of a real-valued function ``f(x)`` to a user-specified accuracy, using user-specified formulas.


##### Fast integral transforms

A piecewise representation of ``f(x)`` with well-chosen formulas enables one to perform fast [**integral transforms**](https://en.wikipedia.org/wiki/Integral_transform) of ``f``. An integral transform with kernel ``K`` is defined as
```math
	(K\circ f)(\mathbf{X}) = \int_{-\infty}^{\infty}dx\,f(x) K(x, \mathbf{X}),
```
where ``\mathbf{X}`` represents a variable or a set of variables and ``K(x, \mathbf{X})`` is the kernel. If ``f(x)`` can be accurately piecewise-approximated by formulas ``F(x, \mathbf{a})`` for which continuous [primitive functions](https://en.wikipedia.org/wiki/Antiderivative) ``\mathcal{F}_K(x, \mathbf{a}, \mathbf{X})`` over the kernel ``K`` are known, i.e,
```math
	\frac{d}{dx}\mathcal{F}_K(x, \mathbf{a}, \mathbf{X})
	= F(x, \mathbf{a}) K(x, \mathbf{X}),
```
then the integral transform can be immediately evaluated as
```math
	(K\circ f)(\mathbf{X}) = \sum_i
	\theta\left(x-x_i^{\min}\right)\theta\left(x_i^{\max}-x\right)
	\sum_j\left[{\mathcal{F}_K}_{ij}(x_i^{\max}, \mathbf{a}_{ij}, \mathbf{X}) - {\mathcal{F}_K}_{ij}(x_i^{\min},  \mathbf{a}_{ij}, \mathbf{X})\right].
```
A simple example of kernel is ``K(x,n)=x^n``, which provides the ``n``-th moment of the function ``f(x)`` defined as ``(M\circ f)(n) = \int_{-\infty}^{\infty}dx\,f(x)x^n``. Other important examples include the Fourier transform with ``K(x, y)=e^{ixy}`` for ``y\in\mathbb{R}``, the Laplace transform with ``K(x, y)=\theta(x)e^{-xy}`` for ``y\in\mathbb{R}``, or the Hilbert transform with kernel ``K(x, y)=1/(y-x)`` for ``y\in\mathbb{C}\setminus\mathbb{R}``, which also yields the [Kramers-Kronig transform](https://en.wikipedia.org/wiki/Kramers–Kronig_relations).

!!! tip
	The module [Piecewise](@ref) provides a generic integral transform that can work with user-defined primitives ``\mathcal{F}_K(x, \mathbf{a}, \mathbf{X})``. It also provides several [`Formula`](@ref) objects with primitives for the kernel ``x^n`` (see [Formulas](@ref)). The module [PiecewiseHilbert](@ref) adds to these formulas the primitives needed for the Hilbert transform. The module [PiecewiseLorentz](@ref) adds the primitives needed for what we call the Lorentz transform.



## Types

##### `PiecewiseFunction`

A piecewise function can be initialized as:
```
f = PiecewiseFunction([parity,] pieces)
```
Once initialized, the piecewise function can be evaluated as `f(x)`. The optional argument `parity` can be either `:none` (default), `:even`, or `:odd`. If `:even` or `:odd` parities are specified, the piecewise function evaluates according to ``f(x<0) = \pm f(-x)``, respectively [^2]. The argument `pieces` is an object of type [`Piece`](@ref) or an array of such objects.

[^2]: This limits the use of `:odd` to function types for which a `-` method exists.

Arithmetic operations with scalars are possible, e.g., `s + f` or `s * f` yield new piecewise functions appropriately transformed. Two piecewise functions can be merged by adding them with `+` (or `sum()` for an array).


##### [`Piece`](@id piecesection)

A piece is characterized by a domain and a rule to return a value. This rule is an object of type [`Formula`](@ref) or an array of such objects, accompanied by parameters:
```
p = Piece(domain, [included,] rule, parameters)
```
The argument `domain::Tuple{Real, Real}` with `domain[1] <= domain[2]` specifies the domain of the piece. The optional argument `included::Tuple{Bool, Bool}`, by default `(true, true)`,  indicates whether the domain boundaries are included in the domain (see [Domains and boundaries](@ref)). If the argument `rule` is of type [`Formula`](@ref), the argument `parameters` must be of type `Vector{Any}`, while if `rule` is of type `Vector{Formula}`, `parameters` must be of type `Vector{Vector{Any}}`. When several formulas are provided, the rule is the sum of all formulas [^1].

A piece can also be initialized with a single function as
```
p = Piece(domain, [included,] function)
```
The argument `function` is a function with interface `function(::Real)`. It can also be passed as a string representing an anonymous function. This method of initialization is equivalent to `p = Piece(domain, [included,] Formula(2, (x, a) -> a[2] * function(a[1] * x)), [1.0, 1.0])`. The parameters `a[1]` and `a[2]` are added for the function to behave correctly under the default scaling and mirroring (see [`Formula`](@ref formulasection)).


##### [`Formula`](@id formulasection)

An object of type [`Formula`](@ref) is created as
```
F = Formula([name,] params, value [, check] [, scale] [, mirror])
```
The optional argument `name::String` is used for printout purposes. In particular, it allows the method [`format`](@ref) to replace the function by its name. The argument `params::Integer` specifies the number of parameters in the formula. If it is positive, the formula takes exactly `params` parameters; if it is negative, the formula takes at most `-param` parameters. This allows one to define formulas using functions with an unspecified number of parameters, as e.g. polynomials of various orders (see [`POLY`](@ref POLY) and [`piecewisefit`](@ref piecewisefitsection)). The argument `value::Function` is a function with interface `value(x::Real, a::Vector{Any})` that returns a value given `x` and the array of parameters `a`. If `value` is an anonymous function, the code defining that function is not stored in the `Formula` object and thus cannot be printed by the method [`format`](@ref). For that reason, it is also possible to pass `value` as a string representing an anonymous function.

The optional argument `check::Function` is a function defined as `check(a, domain, included, danger, fatal)`. This function must return `true` if the parameters `a` can be used in the function `value(x, a)` for `x` in the domain specified by the arguments `domain` and `included` (see [`Piece`](@ref piecesection)), and `false` otherwise. The arguments `danger::Bool` and `fatal::Bool` should be used to switch on and off warning and error messages, respectively (see **Example** in [`Formula`](@ref)). By default, `check = (a, domain, included, danger, fatal) -> true`, which means that no check is performed.

The optional argument `scale::Function` is a function with interface `scale(a::Vector{Any}, s::Number)` that returns the parameters of the formula after multiplication by the scalar `s`. By default, `scale = (a, s) -> [a[1:end-1]..., a[end] * s]`, which means that only the last parameter is scaled. The optional argument `mirror::Function` is a function with interface `mirror(a::Vector{Any})` that returns the parameters of the formula after even reflection through ``x=0``. By default, `mirror = a -> [-a[1], a[2:end]...]`, which means that the first parameter is negated and the others are unchanged.


!!! warning
	If the arguments `value`, `check`, `scale`, and `mirror` are function names (rather than anonymous functions), the function name is stored in the structure [`Formula`](@ref), not the function definition. If the function is redefined after the [`Formula`](@ref) object has been created, that object will use the new function.



## Domains and boundaries

The following rules determine which value is returned by the method `PiecewiseFunction(x)`.

* The value `0` is returned if `x` does not belong to any of the domains defined in `PiecewiseFunction.pieces`, with two exceptions:
    1. The value `Inf` is returned if `x` coincides with the common boundary of two domains (this boundary being excluded from both domains) and if the rules for `x-ϵ` and `x+ϵ` both yield positive values, where `ϵ = 10 *  eps(Float64)` [^3].
    2. The value `-Inf` is returned in the same situation if the rules for `x-ϵ` and `x+ϵ` both yield negative values [^3].
* If `x` coincides with the common boundary of two domains and if this boundary is included in both domains, the rule of the leftmost domain is used.
* If `x` falls inside the domain of a piece, the rule of that piece is used.

[^3]: This limits the treatment of singularities to function types for which a `sign` method exists.



## Main methods

##### `integraltransform`

A function returning the integral transform of the piecewise function `f` may be created as:
```
Kof(X) = integraltransform(f::PiecewiseFunction, X::Any)
```
This assumes that for each [`Formula`](@ref) object `F` used in the pieces of `f`, there is a method `F.value(x::Real, a::Vector{Any}, X::Any)` returning the primitive of `F.value(x, a) * K(x, X)` for the kernel `K` of interest (see [Fast integral transforms](@ref)). If the function `f` is `:even` or `:odd`, the method [`integraltransform`](@ref) requires instead a method `F.value(s, x::Real, a::Vector{Any}, X::Any)` with `s = 1` or `s = -1`, respectively, which returns the primitive of `F.value(x, a) * (K(x, X) + s * K(-x, X))`.


##### `moment`

The moment of order ``n`` of the piecewise function `f` may be obtained as
```
m = moment(f::PiecewiseFunction, n::Integer)
```
This assumes that for each [`Formula`](@ref) object `F` used in the pieces of `f`, there is a method `F.value(x::Real, a::Vector{Any}, n::Integer)` returning the primitive of `F.value(x, a) * x^n`. The argument `n` must be a non-negative integer.


##### [`piecewisefit`](@id piecewisefitsection)

Given a function ``g(x)``, a [`PiecewiseFunction`](@ref) object approximating this function can be constructed as:
```
f = piecewisefit(g, domain, formulas; kwargs...)
```
The argument `g(::Real)` is the function to approximate. The algorithm tries to fit any linear superposition of the formulas given in the argument `formulas::Vector{Formula}` to ``g(x)`` in the domain specified by the argument `domain::Tuple{Real, Real}`. The fit is successful if it matches the function within a specified tolerance and if the fitted parameters pass the tests implemented in the method `F.check` of each [`Formula`](@ref) object used. If the fit is unsuccessful, the domain is divided in two and the algorithm continues recursively in each sub-domain. The recursion uses the available threads, although the parallelism isn't optimal yet.

The argument `formulas` can hold several [`Formula`](@ref) objects with fixed number of parameters and at most one [`Formula`](@ref) object with a variable number of parameters. The formulas with fixed number of parameters are tried first, then the one with variable number of parameters, if any, progressively increasing the number of parameters up to the maximum number allowed. Then pairwise linear combinations are tried, and so on until one combination succeeds.

If the function ``g(x)`` is noiseless, the algorithm should converge with relatively large sub-domains. If ``g(x)`` has numerical noise, the algorithm will likely be trapped trying to piecewise-fit that noise. To avoid this, a minimal sub-domain size may be specified with the optional argument `grain`. If `grain > eps(Float64)`, the algorithm returns the best possible fit in sub-domains of typical size `grain`. If no fit is successful, or if `grain = eps(Float64)`, the algorithm returns a [`Piece`](@ref) object that uses either the function ``g(x)`` if the optional argument `loop` is `true`, or a linear interpolation of ``g(x)`` across the domain if it is `false` (default). On the other hand, relevant rapid variations of ``g(x)`` may be missed by the algorithm, because it tries to minimize the number of calls to `g(x)`. The optional argument `resolution` allows one to control how finely the function ``g(x)`` is sampled.

The optional keyword arguments are:

*  `parity::Symbol`: imposes a given parity (`:even` or `:odd`) to the [`PiecewiseFunction`](@ref) object (by default `parity = :none`)

* `singularities::Vector{Real}`: values of ``x`` that are treated as singularities, i.e., excluded from any sub-domain generated during the recursion (by default `singularities = []`)

* `cuts::Vector{Real}`: values of ``x`` that are forced to be piece boundaries (by default `cuts = []`)

* `grain::Real`: sub-domains of size smaller than `grain` are not split further by the algorithm (by default `grain = eps(Float64)`)

* `resolution::Real`: no fit is accepted before sampling the function ``g(x)`` in each sub-domain with at most a distance `resolution` between successive points (by default `resolution = Inf`)

* `rtol::Real`, `atol::Real`: a fit is successful if `abs(dg) < rtol * abs(g) + atol`, where `dg` are the residuals and `g` are the function values (by default `rtol = 0`, `atol=eps(Float64)`)

* `loop::Bool`: if `true`, the returned [`PiecewiseFunction`](@ref) object uses the function ``g(x)`` in sub-domains where the fit fails, otherwise it uses a linear interpolation of ``g(x)`` across the sub-domain (by default `loop = false`).



## Example

The following example adds [`PiecewiseFunction`](@ref) objects to display the [Cantor set](https://en.wikipedia.org/wiki/Cantor_set).

```@example
ENV["GKSwstype"] = "100" # hide
using ..Piecewise # hide
# Rule to construct the set
cut(d) = [(d[1], d[1] + 1//3 * (d[2] - d[1])),
    (d[1] + 2//3 * (d[2] - d[1]), d[2])]
set(n) = n == 0 ? [(0, 1)] : vcat(cut.(set(n - 1))...)

# Piecewise function at order n
cantor(n) = PiecewiseFunction([Piece(d, x -> exp(-n / 3)) for d in set(n)])

# Plot sum of piecewise functions
using Plots
# Define array x with all breakpoints
x = vcat(map(b -> [b - eps(Float64), b + eps(Float64)],
    vcat(map(d -> [d[1], d[2]], set(6))...))...)
plot(x, sum([cantor(n) for n = 0:6]).(x),
    f=(0, 0, :black), linewidth=0, axis=false, grid=false, legend=:none)
savefig("cantor.svg"); nothing # hide
```

![](cantor.svg)



## Formulas

The module [Piecewise](@ref) provides several [`Formula`](@ref) objects with appropriate primitives for computing the moments of a piecewise function. In the descriptions below, ``x_{\min}`` and ``x_{\max}`` refer to the boundaries of the domain in which the formula is used and ``{_2F_1}(a,b,c,z)`` refers to the [hypergeometric function](https://en.wikipedia.org/wiki/Hypergeometric_function). This function has a branch cut on the real axis ``z=x\in\mathbb{R}`` for ``x>1``. In all cases considered here, ``{_2F_1}(a,b,c,x)`` refers to the value below the cut, i.e., ``{_2F_1}(a,b,c,x-i0)``.

[`POLY`](@ref POLY) | 
[`TAIL`](@ref TAIL) | 
[`LOG`](@ref LOG) |
[`ISRS`](@ref ISRS) |
[`PLS`](@ref PLS) |
[`XLOG`](@ref XLOG) |
[`XISRS`](@ref XISRS)

---

##### [`POLY`](@id POLY)

* *Polynomial of varying order*
* *1 to 13 parameters*

The function `POLY.value(x, a)` is
```math
F(x,\mathbf{a}) = \sum_{i=1}^k a_i x^{i-1},
```
where the number of parameters can vary in the range ``1\leqslant k\leqslant 13``. The order of the polynomial is limited to 12 in order to reduce numerical instabilities. The parameters ``a_i`` are unrestricted.

The primitive function for the moment of order ``n`` is
```math
\mathcal{F}(x,\mathbf{a},n) = \sum_{i=1}^k \frac{a_i}{i+n} x^{i+n}.
```
This primitive is continuous for ``x\in\mathbb{R}``.

---

##### [`TAIL`](@id TAIL)

* *Rational function approaching zero as ``1/x`` or ``1/x^2`` at infinity*
* *5 parameters*

The function `TAIL.value(x, a)` is
```math
F(x,\mathbf{a}) = \frac{a_1+a_2x}{a_3+a_4x+a_5x^2}.
```
The parameters must satisfy either ``a_4^2-4a_3a_5 < 0``, such that the zeros of the denominator are not on the real axis, or they must ensure that the two zeros ``x_{\pm}=\frac{1}{2a_5}\left(-a_4\pm\sqrt{a_4^2-4a_3a_5}\right)`` lie outside the domain.

The primitive function for the moment of order ``n`` is
```math
\begin{align*}
\mathcal{F}(x,\mathbf{a},n) &= \frac{x^{n+1}}{2(n+1)a_3\Delta}\left\{
\left[2a_2a_3-a_1\left(a_4-\Delta\right)\right]
{_2F_1}\left(1,n+1,n+2,\frac{-2a_5x}{a_4+\Delta}\right)\right.\\
&\quad\left.-\left[2a_2a_3-a_1\left(a_4+\Delta\right)\right]
{_2F_1}\left(1,n+1,n+2,\frac{-2a_5x}{a_4-\Delta}\right)\right\}.
\end{align*}
```
where ``\Delta=\sqrt{a_4^2-4a_3a_5}``. Simpler expressions are used if, e.g., ``a_3=0``. This primitive is continuous for ``x\in\mathbb{R}``, except at the zeros ``x=x_{\pm}``, if they are on the real axis.

---

##### [`LOG`](@id LOG)

* *Logarithmic singularity*
* *2 parameters*

The function `LOG.value(x, a)` is
```math
F(x,\mathbf{a}) = a_2\ln|x-a_1|.
```
The domain should not include the point ``x=a_1``.

The primitive function for the moment of order ``n`` is
```math
\mathcal{F}(x,\mathbf{a},n) = \frac{a_2x^{n+1}}{(n+1)(n+2)a_1}\left[
(n+2)a_1\ln|x-a_1|+{_2F_1}\left(1,n+2,n+3,\frac{x}{a_1}\right)x\right].
```
A simpler expression is used if ``a_1=0``. This primitive is continuous in any domain excluding ``x=a_1``.

---

##### [`ISRS`](@id ISRS)

* *Inverse square-root singularity*
* *2 parameters*

The function `ISRS.value(x, a)` is
```math
F(x,\mathbf{a}) = \frac{a_2}{\sqrt{|x^2-a_1^2|}}.
```
The domain should not include any of the points ``x=\pm a_1``.

The primitive function for the moment of order ``n`` is
```math
\mathcal{F}(x,\mathbf{a},n) = \frac{a_2x^{n+1}}{(n+1)\sqrt{|x^2-a_1^2|}}
\sqrt{1-\left(\frac{x}{a_1}\right)^2}{_2F_1}\left(\frac{1}{2},\frac{n+1}{2},
\frac{n+3}{2},\left(\frac{x}{a_1}\right)^2\right).
```
A simpler expression is used if ``a_1=0``. This primitive is continuous in any domain excluding ``x=\pm a_1``.

---

##### [`PLS`](@id PLS)

* *Power-law singularity*
* *3 parameters*

The function `PLS.value(x, a)` is
```math
F(x,\mathbf{a}) = a_3|x-a_1|^{a_2}.
```
The parameter ``a_1`` must be outside the domain, i.e. ``a_1\leqslant x_{\min}`` or ``a_1\geqslant x_{\max}``, with the equal sign allowed only if ``a_2\geqslant0``. The exponent ``a_2`` must be in the range ``[-12, 12]`` (to reduce numerical instabilities like in [`POLY`](@ref POLY)).

The primitive function for the moment of order ``n`` is
```math
\mathcal{F}(x,\mathbf{a},n) = \frac{a_3|x-a_1|^{a_2}x^{n+1}}{n+1}
\left(1-\frac{x}{a_1}\right)^{-a_2}
{_2F_1}\left(n+1,-a_2,n+2,\frac{x}{a_1}\right).
```
A simpler expression is used if ``a_1=0``. This primitive is continuous in any domain excluding ``x=a_1``.

---

##### [`XLOG`](@id XLOG)

* *Logarithmic singularity times ``x``*
* *2 parameters*

The function `XLOG.value(x, a)` is
```math
F(x,\mathbf{a}) = a_2x\ln|x-a_1|.
```
The domain should not include the point ``x=a_1``.

The primitive function for the moment of order ``n`` is
```math
\mathcal{F}(x,\mathbf{a},n) = \frac{a_2x^{n+2}}{(n+2)(n+3)a_1}\left[
(n+3)a_1\ln|x-a_1|+{_2F_1}\left(1,n+3,n+4,\frac{x}{a_1}\right)x\right].
```
A simpler expression is used if ``a_1=0``. This primitive is continuous in any domain excluding ``x=a_1``.

---
 
##### [`XISRS`](@id XISRS)

* *Inverse square-root singularity times ``x``*
* *2 parameters*

The function `XISRS.value(x, a)` is
```math
F(x,\mathbf{a}) = \frac{a_2x}{\sqrt{|x^2-a_1^2|}}.
```
The domain should not include any of the points ``x=\pm a_1``.

The primitive function for the moment of order ``n`` is
```math
\mathcal{F}(x,\mathbf{a},n) = \frac{a_2x^{n+2}}{(n+2)\sqrt{|x^2-a_1^2|}}
\sqrt{1-\left(\frac{x}{a_1}\right)^2}{_2F_1}\left(\frac{1}{2},\frac{n}{2}+1,
\frac{n}{2}+2,\left(\frac{x}{a_1}\right)^2\right).
```
A simpler expression is used if ``a_1=0``. This primitive is continuous in any domain excluding ``x=\pm a_1``.


## Public interface

### Index

```@index
Pages   = ["index.md"]
Modules = [Piecewise]
Order   = [:type, :constant, :function]
```


### Types

```@docs
Formula
Piece
PiecewiseFunction
```


### Constants

```@docs
POLY
TAIL
LOG
ISRS
PLS
XLOG
XISRS
```


### Methods

```@docs
domains
intervals
support
singularities
formulas
integraltransform
moment
piecewisefit
format
```

