

"""
    HilbertTransform(f::PiecewiseFunction[, radius::Real])

Return a `HilbertTransform` object for the piecewise function `f`.

The `HilbertTransform` object behaves as a function with argument `(z::Complex)`. A moment
expansion is used beyond ``|z|`` = radius in the complex plane. If not specified, the radius
is set automatically. The field `error` gives an estimate of the absolute error at ``|z|`` =
radius.

## Fields

* `f::PiecewiseFunction`
* `radius::Real`
* `moments::Vector{Real}`
* `error::Real`

## Example

Hilbert transform of a box:
```jldocs
julia> H = HilbertTransform(PiecewiseFunction(Piece((-1, 1), POLY, [1])))
< Hilbert transform of piecewise function with support [-1.0, 1.0] >

julia> H.([1im, 10im, 100im])
3-element Vector{ComplexF64}:
 0.0 - 1.5707963267948966im
 0.0 - 0.19933730476190478im
 0.0 - 0.019999333333333334im
```
"""
struct HilbertTransform

    f::PiecewiseFunction
    radius::Real
    moments::Vector{Real}
    error::Real

    # Type constructor
    function HilbertTransform(f::PiecewiseFunction, radius::Real)
        # Check the formulas used by f
	    [[hasmethod(F.value, (Real, Vector{T} where T<:Any, Complex)) ||
	        throw(ArgumentError(
	        "The function $(F.value) used by the formula $(F.name) must have a "*
	        "method $(F.value)(::Real, ::Vector{Any}, ::Complex)"))
	        for F in p.rule] for p in f.pieces]
        # Calculate the moments and the typical error at radius
        moments = _moments(f)
        error = _error(f, abs(radius), moments)
        return new(f, abs(radius), moments, error)
    end

    # Initialization with missing radius
    function HilbertTransform(f::PiecewiseFunction)
        # Check the formulas used by f
	    [[hasmethod(F.value, (Real, Vector{T} where T<:Any, Complex)) ||
	        throw(ArgumentError(
	        "The function $(F.value) used by the formula $(F.name) must have a "*
	        "method $(F.value)(::Real, ::Vector{Any}, ::Complex)"))
	        for F in p.rule] for p in f.pieces]
        # Calculate the moments, the radius, and the error
        moments = _moments(f)
        radius = _automatic_radius(f, moments)
        error = _error(f, radius, moments)
        return new(f, radius, moments, error)
    end

end
