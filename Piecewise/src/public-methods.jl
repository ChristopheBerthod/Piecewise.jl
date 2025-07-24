
# Base.show method for objects of type Formula
Base.show(io::IO, ::MIME"text/plain", F::Formula) = write(io, F.constructor)


# Base.print method for objects of type Formula
Base.print(io::IO, F::Formula) = (write(io, F.constructor); nothing)
Base.print(F::Formula) = print(stdout, F)


# Base.show method for objects of type Piece
function Base.show(io::IO, ::MIME"text/plain", p::Piece)
    p.included[1] ? (d = "[") : (d = "]")
    d *= "$(round(p.domain[1], digits=4)), $(round(p.domain[2], digits=4))"
    p.included[2] ? (d *= "]") : (d *= "[")
    label(F) = F.name == "" ? "unnamed" : F.name
    t = "$(label(p.rule[1]))"
    [t *= " + $(label(p.rule[i]))" for i=2:length(p.rule)]
    write(io, "< Piece of type $(t) in the domain $(d) >")
end


# Base.print method for objects of type Piece
Base.print(io::IO, p::Piece) = (write(io, _format(p)); nothing)
Base.print(p::Piece) = print(stdout, p)


# Base.show method for objects of type PiecewiseFunction
function Base.show(io::IO, ::MIME"text/plain", f::PiecewiseFunction)
    f.parity === :none ? (p = "") : (p = "$(f.parity) ")
    d = "[$(round(support(f)[1], digits=4)), $(round(support(f)[2], digits=4))]"
    n = length(f.pieces); n > 1 ? (s = "s") : (s = "")
    write(io, "< Piecewise $(p)function with $(n) piece$(s) and support $(d) >")
end


# Base.print method for objects of type PiecewiseFunction
Base.print(io::IO, f::PiecewiseFunction) = (write(io, _format(f)); nothing)
Base.print(f::PiecewiseFunction) = print(stdout, f)


# Objects of type Formula, Piece, and PiecewiseFunction are not iterable
# and behave like scalars
Base.broadcastable(F::Formula) = Ref(F)
Base.broadcastable(p::Piece) = Ref(p)
Base.broadcastable(f::PiecewiseFunction) = Ref(f)


# Addition operator for a piecewise function and a scalar
function +(f::PiecewiseFunction, s::Number)
    if f.parity === :even
        return f + PiecewiseFunction(:even, Piece((0, Inf), POLY, [s]))
    else
        return _unfold(f) + PiecewiseFunction(Piece((-Inf, Inf), POLY, [s]))
    end
end    
+(s::Number, f::PiecewiseFunction) = f + s
-(f::PiecewiseFunction, s::Number) = f + (-s)
-(s::Number, f::PiecewiseFunction) = s + (-1) * f


# Multiplication operator for a piece and a scalar
*(p::Piece, s::Number) = Piece(p.domain, p.included, p.rule, Vector{Vector{Any}}(
    [F.scale(a, s) for (F, a) in zip(p.rule, p.parameters)]))
*(s::Number, p::Piece) = p * s
/(p::Piece, s::Number) = p * (1 / s)


# Multiplication operator for a piecewise function and a scalar
*(f::PiecewiseFunction, s::Number) = PiecewiseFunction(f.parity,
    [p * s for p in f.pieces])
*(s::Number, f::PiecewiseFunction) = f * s
/(f::PiecewiseFunction, s::Number) = f * (1 / s)


# Addition operator for two piecewise functions
function +(f1::PiecewiseFunction, f2::PiecewiseFunction)
    # If the functions have different parities, they need to be "unfolded"
    # before being summed
    f1.parity !== f2.parity && return _unfold(f1) + _unfold(f2)
    # Sorted list of the new boundaries
    b = sort(unique(vcat([[p.domain[1], p.domain[2]]
        for p in union(f1.pieces, f2.pieces)]...)))
    # isequal(-0.0, 0.0) is false in Julia: we must remove -0.0 is 0.0 is there
    -0.0 in b && 0.0 in b && filter!(x -> !isequal(x, -0.0), b)
    # Build the new pieces
    pieces = Array{Piece}(undef, 0)
    for i = 1:length(b)-1
        domain = (b[i], b[i + 1])
        included = (_in(f1, f2, b[i], :left), _in(f1, f2, b[i + 1], :right))
        rule, parameters = Array{Formula}(undef, 0), Array{Vector{Any}}(undef, 0)
        j = _in(f1, domain)
        j > 0 && (push!(rule, f1.pieces[j].rule...);
            push!(parameters, f1.pieces[j].parameters...))
        j = _in(f2, domain)
        j > 0 && (push!(rule, f2.pieces[j].rule...);
            push!(parameters, f2.pieces[j].parameters...))
        length(rule) > 0 &&
            push!(pieces, Piece(domain, included, rule, parameters))
    end
    return PiecewiseFunction(f1.parity, pieces)
end
-(f1::PiecewiseFunction, f2::PiecewiseFunction) = f1 + (-1) * f2


# Base.sum for an array of piecewise functions
function Base.sum(f::Array{PiecewiseFunction})
    length(f) > 0 || throw(ArgumentError(
        "sum(f::Array{PiecewiseFunction}) called with empty array"))
    s = f[1]
    for i = 2:length(f)
        s += f[i]
    end
    return s
end


# Method to evaluate a PiecewiseFunction object
(f::PiecewiseFunction)(x::Real) = _evaluatepiecewise(f, x)


"""
    domains(f::PiecewiseFunction)

Return the domains of the piecewise function `f` as a `Vector{Tuple{Real, Real}}`.
"""
domains(f::PiecewiseFunction) = [p.domain for p in f.pieces]


"""
    intervals(f::PiecewiseFunction)

Return a `Vector{String}` representing the domains and boundaries of the piecewise function `f`.
"""
function intervals(f::PiecewiseFunction)
    i = String[]
    for p in f.pieces
        p.included[1] ? (d = "[") : (d = "]")
        d *= "$(p.domain[1]), $(p.domain[2])"
        p.included[2] ? (d *= "]") : (d *= "[")
        push!(i, d)
    end
    return i
end


"""
    support(f::PiecewiseFunction)

Return the support of the piecewise function `f` as a `Tuple{Real, Real}`.
"""
function support(f::PiecewiseFunction)
    if f.parity in (:even, :odd)
        return (-f.pieces[end].domain[2], f.pieces[end].domain[2])
    else
        return (f.pieces[1].domain[1], f.pieces[end].domain[2])
    end
end


"""
    singularities(f::PiecewiseFunction)

Return the singularities of the piecewise function `f` as a `Vector{Real}`.
"""
function singularities(f::PiecewiseFunction)
    s = Real[]
    ((f.parity in (:even, :odd)) && ! f.pieces[1].included[1]) && append!(s, 0)
    for i=1:length(f.pieces)-1
        if(! f.pieces[i].included[2] && ! f.pieces[i + 1].included[1])
            append!(s, f.pieces[i].domain[2])
        end
    end
    return s
end


"""
    formulas(f::PiecewiseFunction)

Return the formula names used in the piecewise function `f` as a `Vector{String}`. 
"""
formulas(f::PiecewiseFunction) =
    unique(vcat([[F.name for F in p.rule] for p in f.pieces]...))


"""
Polynomial of varying order (1 to 13 parameters)
```math
F(x,\\mathbf{a}) = \\sum_{i=1}^n a_i x^{i-1}
```
"""
POLY = Formula("POLY", -13, POLY_F, POLY_F, POLY_F)


"""
Rational function approaching zero as ``1/x`` or ``1/x^2`` at infinity (5 parameters)
```math
F(x,\\mathbf{a}) = \\frac{a_1+a_2x}{a_3+a_4x+a_5x^2}
```
"""
TAIL = Formula("TAIL", 5, TAIL_F, TAIL_F, TAIL_F, TAIL_F)


"""
Logarithmic singularity (2 parameters)
```math
F(x,\\mathbf{a}) = a_2\\ln|x-a_1|
```
"""
LOG = Formula("LOG", 2, LOG_F, LOG_F)


"""
Inverse square-root singularity (2 parameters)
```math
F(x,\\mathbf{a}) = \\frac{a_2}{\\sqrt{|x^2-a_1^2|}}
```
"""
ISRS = Formula("ISRS", 2, ISRS_F, ISRS_F)


"""
Power-law singularity (3 parameters)
```math
F(x,\\mathbf{a}) = a_3|x-a_1|^{a_2}
```
"""
PLS = Formula("PLS", 3, PLS_F, PLS_F)


"""
Logarithmic singularity times ``x`` (2 parameters)
```math
F(x,\\mathbf{a}) = a_2x\\ln|x-a_1|
```
"""
XLOG = Formula("XLOG", 2, XLOG_F, XLOG_F, XLOG_F)


"""
Inverse square-root singularity times ``x`` (2 parameters)
```math
F(x,\\mathbf{a}) = \\frac{a_2x}{\\sqrt{|x^2-a_1^2|}}
```
"""
XISRS = Formula("XISRS", 2, XISRS_F, XISRS_F, XISRS_F)


"""
    integraltransform(f::PiecewiseFunction, X::Any)
    
Return the integral transform of the piecewise function `f`, defined as
```math
(K\\circ f)(\\mathbf{X}) = \\int_{-\\infty}^{\\infty}dx\\,f(x)K(x,\\mathbf{X}).
```

This method assumes that each function `F(::Real, ::Vector{Any})` used in the piecewise
function `f` has a method `F(::Real, ::Vector{Any}, ::Any)` that returns the primitive of
``F(x, \\mathbf{a})`` multiplied by the kernel ``K(x, \\mathbf{X})``, i.e., `d/dx F(x, a, X)
= F(x, a) * K(x, X)`. If `f.parity` is `:even` or `:odd`, `integraltransform` requires a
method `F(::Integer, ::Real, ::Vector{Any}, ::Any)` such that
`d/dx F(1, x, a, X) = F(x, a) * (K(x, X) + K(-x, X))` and
`d/dx F(-1, x, a, X) = F(x, a) * (K(x, X) - K(-x, X))`.

## Example
This define a `Formula` object `LIN` that can be used to represent any piecewise linear
function, and such that `integraltransform` provides its Fourier transform:
```jldocs
julia> F(x, a) = a[1] + a[2] * x;

julia> LIN = Formula("LIN", 2, F, (a, s) -> a * s, a -> [a[1], -a[2]]);

julia> F(x, a, k) = (a[1] * k + a[2] * (k * x + im)) * exp(im * k * x)/(im * k^2);

julia> F(s, x, a, k) = 2 * (s == 1 ? real(F(x, a, k)) : im * imag(F(x, a, k)));

julia> f = PiecewiseFunction(:even, [Piece((0, π), LIN, [π, -1])])
< Piecewise even function with 1 piece and support [-3.1416, 3.1416] >

julia> integraltransform(f, 1)
4.0
```
"""
integraltransform(f::PiecewiseFunction, X::Any) =
    sum([_evaluateprimitive(p, f.parity, X) for p in f.pieces])


"""
    moment(f::PiecewiseFunction, n::Integer)

Return the moment of order `n` of the piecewise function `f`, defined as
```math
(M\\circ f)(n) = \\int_{-\\infty}^{\\infty}dx\\,f(x)x^n.
```
"""
function moment(f::PiecewiseFunction, n::Int)
    n >= 0 || throw(DomainError(n,
        "Argument n must be a non-negative integer."))
    f.parity === :even && isodd(n) && return 0.0
    f.parity === :odd && iseven(n) && return 0.0
    return real(sum([_evaluatemoment(p, f.parity, n) for p in f.pieces]))
end


"""
    piecewisefit(f::Function, domain::Tuple{Real, Real}, formulas::Vector{Formula}; kwargs...)

Return a piecewise approximation of the real-valued function `f(::Real)` in the domain
`domain`, using the formulas given in the array `formulas` (see [`Formula`](@ref)).

## Optional keyword arguments
* `parity` : Impose a parity (`:even` or `:odd`, default `:none`) to
   the piecewise function
* `singularities` : Points excluded from the domains (default `Real[]`)
* `cuts` : Points forced to be piece boundaries (default `Real[]`)
* `grain` : Minimal domain size, default `eps(Float64)`
* `resolution` : Maximal distance between sampled points, default `Inf`
* `rtol` : Relative tolerance, default `0.0`
* `atol` : Absolute tolerance, default `eps(Float64)`
* `loop` : Whether to return `f` in case of failure in a domain, default `false`

## Example
Here is a one-piece approximation to the function ``\\sin^{-1}(x)``. The power law at
``x=1`` is represented exactly and a polynomial is fit to the rest:
```jldocs
julia> f = PiecewiseFunction(:odd, Piece((0, 1), PLS, [1, 1 / 2, -sqrt(2)])) +
           piecewisefit(x -> asin(x) + sqrt(2 - 2 * x), (0, 1), [POLY],
           parity=:odd, atol=1e-5)
< Piecewise odd function with 1 piece and support [-1.0, 1.0] >

julia> max([abs(f(x) - asin(x)) for x in -1:0.1:1]...) < 1e-5
true
```
"""
function piecewisefit(f::Function, domain::Tuple{Real, Real},
    formulas::Vector{Formula}; parity::Symbol = :none,
    singularities::Vector{T} where T<:Real = Real[],
    cuts::Vector{T} where T<:Real = Real[], grain::Real = eps(Float64),
    resolution::Real = Inf, rtol::Real = 0.0, atol::Real = eps(Float64),
    loop::Bool=false)

    # Check arguments
    applicable(f, 0.0) || throw(ArgumentError(
        "The function $(f) must have a method $(f)(::Real)."))
    typeof(f(domain[1] + eps(Float64))) <: Real || throw(ArgumentError(
        "The function $(f) must return a real number."))
    domain[2] >= domain[1] || throw(ArgumentError(
        "Illegal domain: $(domain)"))
    (domain[1] == -Inf || domain[2] == Inf) && throw(ArgumentError(
        "Fit in unbounded domains not yet implemented"))
    length(unique(formulas)) == length(formulas) || throw(ArgumentError(
        "The list of formulas must contain different formulas."))
    for F in formulas
        typeof(F.value(domain[1] + eps(Float64), zeros(abs(F.params)))) <: Real ||
        throw(ArgumentError("The function $(F.value) must return a real number."))
    end
    parity in (:none, :even, :odd) || throw(ArgumentError(
        "Illegal parity: :$(parity)"))

    # Build domains and boundaries, considering parity, excluding singularities
    # and taking into account cuts
    b = Vector{Real}([domain[1], domain[2]])
    parity in (:even, :odd) && (b[1] = max(0.0, domain[1]))
    for s in singularities
        if b[1] <= s <= b[2]
            push!(b, s)
        else
            @info "Ignoring singularity at $(s)"
        end
    end
    for c in cuts
        if b[1] <= c <= b[2]
            push!(b, c)
        else
            @info "Ignoring cut at $(c)"
        end
    end
    b = unique(sort(b))
    domains = map(d -> [(d[1], d[2]),
        (!(d[1] in singularities), !(d[2] in singularities))],
        [(b[i], b[i + 1]) for i = 1:length(b)-1])

    # Build array with all combinations of formulas
    C = Array{Vector{Formula}}(undef, 0)
    for F in formulas
        push!(C, [F])
        for i = 1:length(C)-1
            push!(C, [C[i]; F])
        end
    end
    # Drop combinations with more than one formula with variable number of parameters
    filter!(C -> count(<(0), [F.params for F in C]) <= 1, C)
    # Sort combinations by increasing length
    sort!(C, by = C -> length(C))
    # Sort formulas in each combination by increasing number of parameters,
    # moving formula with variable number of parameters to the end
    [sort!(C, by = F -> F.params >=0 ? F.params : Inf) for C in C]

    # Build the list of pieces
    pieces = Array{Piece}(undef, 0)
    # not sure about a possible race condition on the array pieces (?)
    @threads for d in domains
        [push!(pieces, p) for p in _fitpieceorsplitdomain(f, d[1], d[2], C,
        abs(grain), abs(resolution), abs(rtol), abs(atol), loop,
        Vector{Tuple{Real, Real}}(undef, 0))]
    end
    # Return the PiecewiseFunction object
    if parity in (:even, :odd)
        return PiecewiseFunction(parity, pieces)
    else
        return PiecewiseFunction(pieces)
    end
end


"""
    format(f::PiecewiseFunction)

Return a string holding the constructor for the `PiecewiseFunction` object `f`. For
[`Formula`](@ref) objects than have a name, the name is used instead of the constructor of
the `Formula` object. Note that `print(f)` is equivalent to `print(format(f))`.
"""
Printf.format(f::PiecewiseFunction) = _format(f)
