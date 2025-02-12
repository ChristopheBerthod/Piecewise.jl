

"""
    Formula([name::String,] params::Integer, value::Function
        [, check::Function] [, scale::Function] [, mirror::Function])

Constructor of a `Formula` object with an optional name and a number of
parameters that is exactly `params` if `params >= 0` and at most `-params` if
`params < 0`. The value of the formula is set by the function `value(::Real,
::Vector{Any})`. The optional function `check(a::Vector{Any},
domain::Tuple{Real, Real}, included::Tuple{Bool, Bool}, danger::Bool,
fatal::Bool)` must return `true` or `false` depending on whether the function
`value(x, a)` is valid in the domain specified by `domain` and `included` (see
[`Piece`](@ref)) with the parameters `a`. Warnings must be issued by `check`
only if `danger` is `true` and errors must be thrown only if `fatal` is `true`.
The optional function `scale(a::Vector{Any}, s::Number)` must return the
parameters after multiplication of the formula by `s`. The optional function
`mirror(a::Vector{Any})` must return the parameters after even reflection of the
formula through ``x=0``.

## Fields

* `name::String`
* `params::Integer`
* `value::Function`
* `check::Function`
* `scale::Function`
* `mirror::Function`

## Example

This creates a square-root singularity:
```julia-repl
julia> srs(x, a) = a[2] * sqrt(x - a[1])
julia> function srs(a, domain, included, danger, fatal)
           t = domain[1] > a[1] || (domain[1] == a[1] && ! included[1])
           !t && fatal && throw(ArgumentError("Singularity must be at left of domain."))
           t || return false
           t = a[2] >= 0
           !t && danger && @warn "Negative singularity in domain \$(domain)."
           return true
       end
julia> F = Formula("SRS", 2, srs, srs)
Formula("SRS", 2, srs, srs)

julia> F.check([1, 1], (0, 2), (true, true), true, false)
false

julia> F.check([-1, -1], (0, 2), (true, true), true, false)
┌ Warning: Negative singularity in domain (0, 2).
└ @ Main REPL[3]:6
true
```
"""
struct Formula

    name::String
    params::Int
    value::Function
    check::Function
    scale::Function
    mirror::Function
    constructor::String

    # Type constructor
    function Formula(name::String, params::Int, value::Function,
        check::Function, scale::Function, mirror::Function, constructor=missing)
        ismissing(constructor) && (constructor = "Formula(\"$(name)\", "*
            "$(params), $(value), $(check), $(scale), $(mirror))")
        hasmethod(value, (Real, Vector{T} where T<:Any)) || throw(ArgumentError(
            "The function $(value) must have a method $(value)(::Real, "*
            "::Vector{Any})."))
        hasmethod(check, (Vector{T} where T<:Any, Tuple{Real, Real},
            Tuple{Bool, Bool}, Bool, Bool)) || throw(ArgumentError(
            "The function $(check) must have a method $(check)(::Vector{Any}, "*
            "::Tuple{Real, Real}, ::Tuple{Bool, Bool}, ::Bool, ::Bool)."))
        hasmethod(scale, (Vector{T} where T<:Any, Real)) || throw(ArgumentError(
            "The function $(scale) must have a method $(scale)(::Vector{Any}, "*
            "::Number)."))
        hasmethod(mirror, (Vector{T} where T<:Any,)) || throw(ArgumentError(
            "The function $(mirror) must have a method $(mirror)(::Vector{Any})."))
        return new(name, params, value, check, scale, mirror, constructor)
    end

    # Initializations with missing functions
    function Formula(name::String, params::Int, value::Function, f1::Function,
        f2::Function, constructor=missing)
        ismissing(constructor) && (constructor = "Formula(\"$(name)\", "*
            "$(params), $(value), $(f1), $(f2))")
        applicable(f1, [0], (0.0, 0.0), (true, true), true, true) &&
        applicable(f2, [0], 0.0) && return Formula(name, params, value,
            f1, f2, a -> [-a[1], a[2:end]...], constructor)
        applicable(f1, [0], (0.0, 0.0), (true, true), true, true) &&
        applicable(f2, [0]) && return Formula(name, params, value,
            f1, (a, s) -> [a[1:end-1]..., a[end] * s], f2, constructor)
        applicable(f1, [0], 0.0) && applicable(f2, [0]) &&
            return Formula(name, params, value,
            (a, domain, included, danger, fatal) -> true, f1, f2, constructor)
        throw(ArgumentError("Illegal function argument $(f1) and/or $(f2)"))
    end

    function Formula(name::String, params::Int, value::Function, f1::Function,
        constructor=missing)
        ismissing(constructor) && (constructor = "Formula(\"$(name)\", "*
            "$(params), $(value), $(f1))")
        applicable(f1, [0], (0.0, 0.0), (true, true), true, true) &&
            return Formula(name, params, value, f1,
            (a, s) -> [a[1:end-1]..., a[end] * s], a -> [-a[1], a[2:end]...],
            constructor)
        applicable(f1, [0], 0.0) && return Formula(name, params, value,
            (a, domain, included, danger, fatal) -> true,
            f1, a -> [-a[1], a[2:end]...], constructor)
        applicable(f1, [0]) && return Formula(name, params, value,
            (a, domain, included, danger, fatal) -> true,
            (a, s) -> [a[1:end-1]..., a[end] * s], f1, constructor)
        throw(ArgumentError("Illegal function argument $(f1)"))
    end

    function Formula(name::String, params::Int, value::Function,
        constructor=missing)
        ismissing(constructor) && (constructor = "Formula(\"$(name)\", "*
            "$(params), $(value))")
        return Formula(name, params, value,
            (a, domain, included, danger, fatal) -> true,
            (a, s) -> [a[1:end-1]..., a[end] * s], a -> [-a[1], a[2:end]...],
            constructor)
    end

    # Initializations with missing name
    Formula(params::Int, value::Function, f1::Function, f2::Function,
        f3::Function) = Formula("", params, value, f1, f2, f3,
        "Formula($(params), $(value), $(f1), $(f2), $(f3))")
    Formula(params::Int, value::Function, f1::Function, f2::Function) =
        Formula("", params, value, f1, f2,
        "Formula($(params), $(value), $(f1), $(f2))")
    Formula(params::Int, value::Function, f1::Function) =
        Formula("", params, value, f1, "Formula($(params), $(value), $(f1))")
    Formula(params::Int, value::Function) = Formula("", params, value,
        "Formula($(params), $(value))")

    # Initialization with value::String instead of value::Function
    function Formula(name::String, params::Int, value::String,
        check::Function, scale::Function, mirror::Function, constructor=missing)
        invokelatest(hasmethod, eval(Meta.parse(value)), (Real, Vector{T} where T<:Any)) ||
            throw(ArgumentError("The function $(value) must have a method "*
            "with interface (::Real, ::Vector{Any})."))
        ismissing(constructor) && (constructor = "Formula(\"$(name)\", "*
            "$(params), \"$(value)\", $(check), $(scale), $(mirror))")
        f = (x, a) -> invokelatest(eval(Meta.parse(value)), x, a)
        return Formula(name, params, f, check, scale, mirror, constructor)
    end
    function Formula(name::String, params::Int, value::String, f1::Function,
        f2::Function, constructor=missing)
        invokelatest(hasmethod, eval(Meta.parse(value)), (Real, Vector{T} where T<:Any)) ||
            throw(ArgumentError("The function $(value) must have a method "*
            "with interface (::Real, ::Vector{Any})."))
        ismissing(constructor) && (constructor = "Formula(\"$(name)\", "*
            "$(params), \"$(value)\", $(f1), $(f2))")
        f = (x, a) -> invokelatest(eval(Meta.parse(value)), x, a)
        return Formula(name, params, f, f1, f2, constructor)
    end
    function Formula(name::String, params::Int, value::String, f1::Function,
        constructor=missing)
        invokelatest(hasmethod, eval(Meta.parse(value)), (Real, Vector{T} where T<:Any)) ||
            throw(ArgumentError("The function $(value) must have a method "*
            "with interface (::Real, ::Vector{Any})."))
        ismissing(constructor) && (constructor = "Formula(\"$(name)\", "*
            "$(params), \"$(value)\", $(f1))")
        f = (x, a) -> invokelatest(eval(Meta.parse(value)), x, a)
        return Formula(name, params, f, f1, constructor)
    end
    function Formula(name::String, params::Int, value::String,
        constructor=missing)
        invokelatest(hasmethod, eval(Meta.parse(value)), (Real, Vector{T} where T<:Any)) ||
            throw(ArgumentError("The function $(value) must have a method "*
            "with interface (::Real, ::Vector{Any})."))
        ismissing(constructor) && (constructor = "Formula(\"$(name)\", "*
            "$(params), \"$(value)\")")
        f = (x, a) -> invokelatest(eval(Meta.parse(value)), x, a)
        return Formula(name, params, f, constructor)
    end
    Formula(params::Int, value::String, f1::Function, f2::Function,
        f3::Function) = Formula("", params, value, f1, f2, f3,
        "Formula($(params), \"$(value)\", $(f1), $(f2), $(f3))")
    Formula(params::Int, value::String, f1::Function, f2::Function) =
        Formula("", params, value, f1, f2,
        "Formula($(params), \"$(value)\", $(f1), $(f2))")
    Formula(params::Int, value::String, f1::Function) =
        Formula("", params, value, f1, "Formula($(params), \"$(value)\", $(f1))")
    Formula(params::Int, value::String) = Formula("", params, value,
        "Formula($(params), \"$(value)\")")

    # Initialization with a string
    function Formula(value::String, constructor=missing)
        invokelatest(hasmethod, eval(Meta.parse(value)), (Real,)) || throw(ArgumentError(
            "The function $(value) must have a method with interface (::Real)."))
        ismissing(constructor) && (constructor = "Formula(\"$(value)\")")
        # Insert two parameters for default scaling and mirroring
        f = (x, a) -> a[2] * invokelatest(eval(Meta.parse(value)), a[1] * x)
        return Formula("", 2, f, constructor)
    end

end


"""
    Piece(domain, [included,] rule, parameters)

Constructor of a `Piece` object in the domain defined by the argument
`domain::Tuple{Real, Real}`. The argument `included::Tuple{Bool, Bool}`, by
default `(true, true)`, tells whether the domain boundaries belong to the
domain. The arguments `rule::Vector{Formula}` and
`parameters::Vector{Vector{Any}}` specify the rule used to evaluate the value of
the piece (see [`Formula`](@ref)). It is also possible to pass
`rule::Formula` and `parameters::Vector{Any}` for a single formula.

    Piece(domain, [included,] function)

Initialization with a single function. This is equivalent to
`Piece(domain, [included,] Formula(0, (x, a)->function(x)), Any[])`.

## Fields

* `domain::Tuple{Real, Real}`
* `included::Tuple{Bool, Bool}`
* `rule::Vector{Formula}`
* `parameters::Vector{Vector{Any}}`

## Examples

This represents ``1/x`` for strictly positive numbers:
```julia-repl
julia> Piece((0, Inf), (false, true), x -> 1 / x)
< Piece of type unnamed over the domain ]0.0, Inf] >

```
This represents ``-\\log|x|`` for ``x\\in ]0, 1]``:
```julia-repl
julia> Piece((0, 1), (false, true), Formula("LOG", 1, (x, a) -> a[1] * log(abs(x))), [-1])
< Piece of type LOG over the domain ]0.0, 1.0] >

```
"""
struct Piece

    domain::Tuple{Real, Real}
    included::Tuple{Bool, Bool}
    rule::Vector{Formula}
    parameters::Vector{Vector{T}} where T<:Any

    # Type constructor
    function Piece(domain::Tuple{Real, Real}, included::Tuple{Bool, Bool},
        rule::Vector{Formula}, parameters::Vector{Vector{T}} where T<:Any)
        # Domain boundaries must be ordered
        domain[2] >= domain[1] || throw(ArgumentError(
            "Illegal domain: $(domain)"))
        # Is there a formula?
        length(rule) > 0 || throw(ArgumentError(
            "A Piece needs at least one formula."))
        # Are formulas and parameters arrays of the same length?
        length(rule) == length(parameters) || throw(ArgumentError(
            "Formula and parameter arrays must have the same length."))
        for (i, F) in enumerate(rule)
            label = F.name == "" ? F : F.name
            # Is the number of parameters consistent with the formula?
            F.params >= 0 && length(parameters[i]) != F.params &&
                throw(ArgumentError(
                "The function $(label) requires $(F.params) parameters."))
            F.params < 0 && length(parameters[i]) == 0 &&
                throw(ArgumentError(
                "The function $(label) requires at least 1 parameter."))
            F.params < 0 && length(parameters[i]) > -F.params &&
                throw(ArgumentError(
                "The function $(label) accepts at most $(-F.params) parameters."))
            # Are the parameters valid in the domain
            F.check(parameters[i], domain, included, true, true)
        end
        return new(domain, included, rule, parameters)
    end

    # Initialization with missing argument
    Piece(domain::Tuple{Real, Real}, rule::Vector{Formula},
        parameters::Vector{Vector{T}} where T<:Any) =
        Piece(domain, (true, true), rule, parameters)

    # Initialization with a single formula
    Piece(domain::Tuple{Real, Real}, included::Tuple{Bool, Bool},
        F::Formula, parameters::Vector{T} where T<:Any) =
        Piece(domain, included, [F], [parameters])

    Piece(domain::Tuple{Real, Real}, F::Formula, parameters::Vector{T}
        where T<:Any) = Piece(domain, (true, true), [F], [parameters])

    # Initialization with a function
    function Piece(domain::Tuple{Real, Real}, included::Tuple{Bool, Bool},
        f::Function)
        applicable(f, 0.0) || throw(ArgumentError(
            "The function $(f) must have a method $(f)(::Real)."))
        # Insert two parameters for default scaling and mirroring
        return Piece(domain, included,
            Formula(2, (x, a) -> a[2] * f(a[1] * x)), [1.0, 1.0])
    end

    Piece(domain::Tuple{Real, Real}, f::Function) = Piece(domain, (true, true), f)

    # Initialization with f::String instead of f::Function
    Piece(domain::Tuple{Real, Real}, included::Tuple{Bool, Bool}, f::String) =
        # Insert two parameters for default scaling and mirroring
        # The actual modification of f is done by the Formula constructor
        Piece(domain, included, Formula(f), [1.0, 1.0])
    Piece(domain::Tuple{Real, Real}, f::String) = Piece(domain, (true, true), f)

end


"""
    PiecewiseFunction([parity::Symbol,] pieces::Vector{Piece})

Constructor of a `PiecewiseFunction`. The optional argument `parity` is `:none`
(default), `:even`, or `:odd`. The argument `pieces` contains the various pieces
(see [`Piece`](@ref)).

## Fields

* `parity::Symbol`
* `pieces::Vector{Piece}`

## Example

This represents the Cauchy principal value of ``1/x``:
```julia-repl
julia> f = PiecewiseFunction(:odd, Piece((0, Inf), (false, true), x -> 1 / x))
< Piecewise odd function with 1 piece and support [-Inf, Inf] >

julia> f.([-1, 0, 1])
3-element Vector{Real}:
 -1.0
  0
  1.0

```
"""
struct PiecewiseFunction

    parity::Symbol
    pieces::Vector{Piece}

    # Type constructor
    function PiecewiseFunction(parity::Symbol, pieces::Vector{Piece})
        # Check arguments
        parity in (:none, :even, :odd) || throw(ArgumentError(
            "Illegal parity: :$(parity)"))
        length(pieces) > 0 || throw(ArgumentError(
            "A PiecewiseFunction must have at least one Piece."))
        # Sort pieces
        ps = sort(pieces, by = p -> p.domain[1])
        if length(ps) > 1
            prod([ps[i].domain[2] <= ps[i+1].domain[1] for i in 1:(length(ps)-1)]
                ) || throw(ArgumentError("Pieces have overlapping domains."))
        end
        if parity in (:even, :odd)
            ps[1].domain[1] >= 0 || throw(ArgumentError("All domains of an "*
                "$(parity) PiecewiseFunction must be on the positive real axis."))
        end
        return new(parity, ps)
    end

    # Initialization with missing argument
    PiecewiseFunction(pieces::Vector{Piece}) = PiecewiseFunction(:none, pieces)

    # Initialization with a single piece
    PiecewiseFunction(pieces::Piece) = PiecewiseFunction(:none, [pieces])

    PiecewiseFunction(parity::Symbol, pieces::Piece) =
        PiecewiseFunction(parity, [pieces])

end
