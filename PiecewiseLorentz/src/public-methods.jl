
# Base.show method for objects of type LorentzTransform
function Base.show(io::IO, ::MIME"text/plain", L::LorentzTransform)
    L.f.parity === :none ? (p = "") : (p = "$(L.f.parity) ")
    d = "[$(round(support(L.f)[1], digits=4)), $(round(support(L.f)[2], digits=4))]"
	write(io, "< L^$(L.m) transform of $(p)piecewise function with support $(d) >")
end


# Base.print method for objects of type L::LorentzTransform
function Base.print(io::IO, L::LorentzTransform)
    write(io, "LorentzTransform(")
    write(io, format(L.f))
    write(io, ", $(L.m), $(L.ground))")
    nothing
end
Base.print(L::LorentzTransform) = print(stdout, L)


# Objects of type LorentzTransform are not iterable and behave like scalars
Base.broadcastable(L::LorentzTransform) = Ref(L)


# Method to evaluate a LorentzTransform object
(L::LorentzTransform)(y::Real, z::Complex) = _evaluatelorentztransform(L, y, z)


"""
    lorentz_transform(f::PiecewiseFunction, m::Integer, y::Real, z::Complex)
    
``Lorentz transform of the piecewise function `f`:
```math
(L^m\\circ f)(y, z) = \\int_{-\\infty}^{\\infty}dx\\,f(x)\\left[\\frac{\\mathrm{Im}
\\,z/\\pi}{(y-\\mathrm{Re}\\,z-x)^2+(\\mathrm{Im}\\,z)^2}\\right]^m
```

## Example
```julia-repl
julia> lorentz_transform(PiecewiseFunction([Piece((-1, 1), POLY, [1])]), 2, 0, -im)
0.13023806336711655
```
"""
function lorentz_transform(f::PiecewiseFunction, m::Int, y::Real, z::Complex)
    # Check the formulas used by f
    [[hasmethod(F.value, (Real, Vector{T} where T<:Any,
        Tuple{Int, Real, Complex})) || throw(ArgumentError(
	    "The function $(F.value) used by the formula $(F.name) must have a "*
	    "method $(F.value)(::Real, ::Vector{Any}, ::Tuple{Int, Real, Complex})"))
	    for F in p.rule] for p in f.pieces]
    imag(z) < 0 || throw(DomainError(z,
        "Argument z must have negative imaginary part."))
    return integraltransform(f, (m, y, z))
end
