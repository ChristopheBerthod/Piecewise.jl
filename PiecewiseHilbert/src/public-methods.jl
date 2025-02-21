
# Base.show method for objects of type HilbertTransform
function Base.show(io::IO, ::MIME"text/plain", H::HilbertTransform)
    H.f.parity === :none ? (p = "") : (p = "$(H.f.parity) ")
    d = "[$(round(support(H.f)[1], digits=4)), $(round(support(H.f)[2], digits=4))]"
	write(io, "< Hilbert transform of $(p)piecewise function with support $(d) >")
end


# Base.print method for objects of type HilbertTransform
function Base.print(io::IO, H::HilbertTransform)
    write(io, "HilbertTransform(")
    write(io, format(H.f))
    write(io, ", $(H.radius))")
    nothing
end
Base.print(H::HilbertTransform) = print(stdout, H)


# Objects of type HilbertTransform are not iterable and behave like scalars
Base.broadcastable(H::HilbertTransform) = Ref(H)


# Method to evaluate a HilbertTransform object
(H::HilbertTransform)(z::Complex) = _evaluatehilberttransform(H, z)


"""
    hilbert_transform(f::PiecewiseFunction, z::Complex)

Return the Hilbert transform of the piecewise function `f` at the complex number `z`:
```math
(H\\circ f)(z) = \\int_{-\\infty}^{\\infty}dx\\,\\frac{f(x)}{z-x}
```

## Example
```jldocs
julia> hilbert_transform(PiecewiseFunction([Piece((-1, 1), POLY, [1])]), im)
0.0 - 1.5707963267948966im
```

See also [`HilbertTransform`](@ref).
"""
function hilbert_transform(f::PiecewiseFunction, z::Complex)
    # Check the formulas used by f
	[[hasmethod(F.value, (Real, Vector{T} where T<:Any, Complex)) ||
	    throw(ArgumentError(
	    "The function $(F.value) used by the formula $(F.name) must have a "*
	    "method $(F.value)(::Real, ::Vector{Any}, ::Complex)"))
	    for F in p.rule] for p in f.pieces]
    imag(z) != 0 || throw(DomainError(z,
        "Argument z must have finite imaginary part."))
    return integraltransform(f, z)
end
