

"""
    LorentzTransform(f::PiecewiseFunction, m::Integer [, ground::Real=1e-10])

Constructor of a `LorentzTransform` object for the piecewise function `f`.
`m` is the order of the transform. The moment expansion is used if the result
of the transform is smaller than `ground`.

## Fields

* `f::PiecewiseFunction`
* `m::Integer`
* `ground::Real`
* `moments::Vector{Real}`

## Example

``L^2`` transform of a box:
```julia-repl
julia> L = LorentzTransform(PiecewiseFunction(Piece((-1, 1), POLY, [1])), 2)
< L^2 transform of piecewise function with support [-1.0, 1.0] >

julia> L(0, -im)
0.13023806336711655
```
"""
struct LorentzTransform

    f::PiecewiseFunction
    m::Int
    ground::Real
    C::Vector{Vector{Real}}

    # Type constructor
    function LorentzTransform(f::PiecewiseFunction, m::Int, ground::Real)
        # Check the formulas used by f
        [[hasmethod(F.value, (Real, Vector{T} where T<:Any,
            Tuple{Int, Real, Complex})) || throw(ArgumentError(
	        "The function $(F.value) used by the formula $(F.name) must have a "*
	        "method $(F.value)(::Real, ::Vector{Any}, ::Tuple{Int, Real, Complex})"))
	        for F in p.rule] for p in f.pieces]
        # Compute the coefficients of the moment expansion
        C = _coefficients(f, m)
        return new(f, m, abs(ground), C)
    end

    # Initialization with missing ground
    LorentzTransform(f::PiecewiseFunction, m::Int) =
        LorentzTransform(f, m, 1e-10)

end
