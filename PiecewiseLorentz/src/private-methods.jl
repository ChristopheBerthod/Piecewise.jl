
# Hypergeometric function 2F1, taking the value below the cut
_H2F1(a::Real, b::Real, c::Real, z::Number) =
    imag(z) ≈ 0 ? _₂F₁(a, b, c, z - im * eps(Float64)) : _₂F₁(a, b, c, z)


## Piecewise.POLY

# Primitive
function Piecewise.POLY_F(x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f01(1, x, a, X)
    X[1] == 1 && return _f011(1, x, a, X[2:3])
    X[1] == 2 && return _f012(1, x, a, X[2:3])
    X[1] == 3 && return _f013(1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for POLY with m = $(X[1])."))
end

# Even and odd primitives
function Piecewise.POLY_F(s::Int, x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f01(1, x, a, X) + s * _f01(-1, x, a, X)
    X[1] == 1 && return _f011(1, x, a, X[2:3]) + s * _f011(-1, x, a, X[2:3])
    X[1] == 2 && return _f012(1, x, a, X[2:3]) + s * _f012(-1, x, a, X[2:3])
    X[1] == 3 && return _f013(1, x, a, X[2:3]) + s * _f013(-1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for POLY with m = $(X[1])."))
end

#function _f01(s::Int, x::Real, a::Vector{T} where T<:Any,
#    X::Tuple{Int, Real, Complex})
# requires the AppellF1 function
# sympy.functions.special.hyper.appellf1
#end

function _f011(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X; u, p = y - z, 0.0
    c0 = s * x / u
    c1 = 2 * conj(u) / (z - conj(z))
    for i = 1:length(a)
        p += a[i] * x^i / i * real(c1 * _H2F1(1, i, i + 1, c0))
    end
    return  -imag(1 / u) / π * p
end

function _f012(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X; u, v, p = y - z, z - conj(z), 0.0
    c0 = s * x / u
    c1 = 1 + 2 * u / v - v / u
    c2 = 1 + v / u
    c3 = -conj(u) / (u - s * x)
    for i = 1:length(a)
        p += a[i] * x^i / i * real((c1 + c2 * i) * _H2F1(1, i, i + 1, c0) + c3 * i)
    end
    return 1 / (2 * π^2 * (real(u)^2 + imag(z)^2)) * p
end

function _f013(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X; u, ux, v, p = y - z, y - z -s * x, z - conj(z), 0.0
    c0 = s * x / u
    c1 = conj(u) / v
    c2 = 2 * (6 - 3 * v / u + (v / u)^2)
    c3 = 3 * v / u * (2 - v / u)
    c4 = (v / u)^2
    c5 = (-6 + 2 * v / u + v / ux) * v / ux
    c6 = -v^2 / u / ux
    for i = 1:length(a)
        p += a[i] * x^i / i * real(c1 * ((c2 + c3 * i + c4 * i^2) *
            _H2F1(1, i, i + 1, c0) + c5 * i + c6 * i^2))
    end
    return -1 / (16 * π^3 * imag(z) * (real(u)^2 + imag(z)^2)) * p
end


## Piecewise.TAIL

# Primitive
function Piecewise.TAIL_F(x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f02(1, x, a, X)
#    X[1] == 1 && return _f021(1, x, a, X[2:3])
#    X[1] == 2 && return _f022(1, x, a, X[2:3])
#    X[1] == 3 && return _f023(1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for TAIL with m = $(X[1])."))
end

# Even and odd primitives
function Piecewise.TAIL_F(s::Int, x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f02(1, x, a, X) + s * _f02(-1, x, a, X)
#    X[1] == 1 && return _f021(1, x, a, X[2:3]) + s * _f021(-1, x, a, X[2:3])
#    X[1] == 2 && return _f022(1, x, a, X[2:3]) + s * _f022(-1, x, a, X[2:3])
#    X[1] == 3 && return _f023(1, x, a, X[2:3]) + s * _f023(-1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for TAIL with m = $(X[1])."))
end


## Piecewise.LOG

# Primitive
function Piecewise.LOG_F(x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f03(1, x, a, X)
    X[1] == 1 && return _f031(1, x, a, X[2:3])
    X[1] == 2 && return _f032(1, x, a, X[2:3])
    X[1] == 3 && return _f033(1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for LOG with m = $(X[1])."))
end

# Even and odd primitives
function Piecewise.LOG_F(s::Int, x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f03(1, x, a, X) + s * _f03(-1, x, a, X)
    X[1] == 1 && return _f031(1, x, a, X[2:3]) + s * _f031(-1, x, a, X[2:3])
    X[1] == 2 && return _f032(1, x, a, X[2:3]) + s * _f032(-1, x, a, X[2:3])
    X[1] == 3 && return _f033(1, x, a, X[2:3]) + s * _f033(-1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for LOG with m = $(X[1])."))
end

function _f031(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X
    u, ux, w, l = y - z - s * a[1], y - z - s * x, x - a[1], log(abs(x - a[1]))
    return s * a[2] / π * imag(l * log(ux / u) + li2(s * w / u))
end

function _f032(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X
    u, ux, w, l = y - z - s * a[1], y - z - s * x, x - a[1], log(abs(x - a[1]))
    return -s * a[2] / (2 * π^2 * imag(z)) * imag(l * log(ux / u) + li2(s * w / u)
        + im * imag(z) * (l / ux + (log(ux) - l) / u))
end

function _f033(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X
    u, ux, w, l = y - z - s * a[1], y - z - s * x, x - a[1], log(abs(x - a[1]))
    return 3 * s * a[2] / (8 * π^3 * imag(z)^2) *  imag(l * log(ux / u)
        + li2(s * w / u) + im * imag(z) * (l / ux + (log(ux) - l) / u)
        + imag(z)^2 / 3 * (l / ux^2 + (log(ux) - l) / u^2 - 1 / (u * ux)))
end


## Piecewise.ISRS

# Primitive
function Piecewise.ISRS_F(x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f04(1, x, a, X)
#    X[1] == 1 && return _f041(1, x, a, X[2:3])
#    X[1] == 2 && return _f042(1, x, a, X[2:3])
#    X[1] == 3 && return _f043(1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for ISRS with m = $(X[1])."))
end

# Even and odd primitives
function Piecewise.ISRS_F(s::Int, x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f04(1, x, a, X) + s * _f04(-1, x, a, X)
#    X[1] == 1 && return _f041(1, x, a, X[2:3]) + s * _f041(-1, x, a, X[2:3])
#    X[1] == 2 && return _f042(1, x, a, X[2:3]) + s * _f042(-1, x, a, X[2:3])
#    X[1] == 3 && return _f043(1, x, a, X[2:3]) + s * _f043(-1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for ISRS with m = $(X[1])."))
end


## Piecewise.PLS

# Primitive
function Piecewise.PLS_F(x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f05(1, x, a, X)
    X[1] == 1 && return _f051(1, x, a, X[2:3])
    X[1] == 2 && return _f052(1, x, a, X[2:3])
    X[1] == 3 && return _f053(1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for PLS with m = $(X[1])."))
end

# Even and odd primitives
function Piecewise.PLS_F(s::Int, x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f05(1, x, a, X) + s * _f05(-1, x, a, X)
    X[1] == 1 && return _f051(1, x, a, X[2:3]) + s * _f051(-1, x, a, X[2:3])
    X[1] == 2 && return _f052(1, x, a, X[2:3]) + s * _f052(-1, x, a, X[2:3])
    X[1] == 3 && return _f053(1, x, a, X[2:3]) + s * _f053(-1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for PLS with m = $(X[1])."))
end

#function _f05(s::Int, x::Real, a::Vector{T} where T<:Any,
#    X::Tuple{Int, Real, Complex})
# requires the AppellF1 function
# sympy.functions.special.hyper.appellf1
#end

function _f051(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X; u, v, w = y - z - s * a[1], z - conj(z), x - a[1]
    return -imag(1 / u) / π * a[3] * abs(w)^a[2] * w / (a[2] + 1) *
        real(conj(u) / v * 2 * _H2F1(1, a[2] + 1, a[2] + 2, s * w / u))
end

function _f052(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X
    u, ux, v, w = y - z - s * a[1], y - z - s * x, z - conj(z), x - a[1]
    return 1 / (2 * π^2 * (real(u)^2 + imag(z)^2)) *
        a[3] * abs(w)^a[2] * w / (a[2] + 1) * real(conj(u) / v * (
        (2 + a[2] * v / u) * _H2F1(1, a[2] + 1, a[2] + 2, s * w / u)
        - (a[2] + 1) * v / ux))
end

function _f053(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X
    u, ux, v, w = y - z - s * a[1], y - z - s * x, z - conj(z), x - a[1]
    return -1 / (16 * π^3 * imag(z) * (real(u)^2 + imag(z)^2)) *
        a[3] * abs(w)^a[2] * w / (a[2] + 1) * real(conj(u) / v * (
        (12 + 6 * a[2] * v / u + a[2] * (a[2] - 1) * (v / u)^2) *
        _H2F1(1, a[2] + 1, a[2] + 2, s * w / u) - 6 * (a[2] + 1) * v / ux
        + (a[2] + 1) * (1 - (a[2] - 1) * ux / u) * (v / ux)^2))
end


## Piecewise.XLOG

# Primitive
function Piecewise.XLOG_F(x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f06(1, x, a, X)
    X[1] == 1 && return _f061(1, x, a, X[2:3])
    X[1] == 2 && return _f062(1, x, a, X[2:3])
    X[1] == 3 && return _f063(1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for XLOG with m = $(X[1])."))
end

# Even and odd primitives
function Piecewise.XLOG_F(s::Int, x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f06(1, x, a, X) + s * _f06(-1, x, a, X)
    X[1] == 1 && return _f061(1, x, a, X[2:3]) + s * _f061(-1, x, a, X[2:3])
    X[1] == 2 && return _f062(1, x, a, X[2:3]) + s * _f062(-1, x, a, X[2:3])
    X[1] == 3 && return _f063(1, x, a, X[2:3]) + s * _f063(-1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for XLOG with m = $(X[1])."))
end

function _f061(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X
    u, ux, w, l = y - z - s * a[1], y - z - s * x, x - a[1], log(abs(x - a[1]))
    return a[2] / π * imag((y - z) * (l * log(ux / u) + li2(s * w / u)))
end

function _f062(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X
    u, ux, w, l = y - z - s * a[1], y - z - s * x, x - a[1], log(abs(x - a[1]))
    return -a[2] / (2 * π^2 * imag(z)) * imag((y - real(z)) * (l * log(ux / u) +
        li2(s * w / u)) + im * (y - z) * imag(z) * (l / ux + (log(ux) - l) / u))
end

function _f063(s::Int, x::Real, a::Vector{T} where T<:Any, X::Tuple{Real, Complex})
    y, z = X
    u, ux, w, l = y - z - s * a[1], y - z - s * x, x - a[1], log(abs(x - a[1]))
    return 3 * a[2] / (8 * π^3 * imag(z)^2) * imag((y - real(z)) * (l * log(ux / u)
        + li2(s * w / u)) + im * (y - z + im * 2 * imag(z) / 3) * imag(z) *
        (l / ux + (log(ux) - l) / u) + imag(z)^2 / 3 * (y - z) * (l / ux^2 +
        (log(ux) - l) / u^2 - 1 / (u * ux)))
end


## Piecewise.XISRS

# Primitive
function Piecewise.XISRS_F(x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f07(1, x, a, X)
#    X[1] == 1 && return _f071(1, x, a, X[2:3])
#    X[1] == 2 && return _f072(1, x, a, X[2:3])
#    X[1] == 3 && return _f073(1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for XISRS with m = $(X[1])."))
end

# Even and odd primitives
function Piecewise.XISRS_F(s::Int, x::Real, a::Vector{T} where T<:Any,
    X::Tuple{Int, Real, Complex})
#    return _f07(1, x, a, X) + s * _f07(-1, x, a, X)
#    X[1] == 1 && return _f071(1, x, a, X[2:3]) + s * _f071(-1, x, a, X[2:3])
#    X[1] == 2 && return _f072(1, x, a, X[2:3]) + s * _f072(-1, x, a, X[2:3])
#    X[1] == 3 && return _f073(1, x, a, X[2:3]) + s * _f073(-1, x, a, X[2:3])
    throw(ArgumentError(
        "Lorentz transform not implemented for XISRS with m = $(X[1])."))
end


# Coefficients of the moment expansion
function _coefficients(f::PiecewiseFunction, m::Int)
    # Compute the moments (at maximum to order 64)
    n, moments = 0, Vector{Real}(undef, 0)
    while n <= 64 && ! isnan(moment(f, n))
        push!(moments, moment(f, n))
        n += 1
    end
    # Compute the coefficients
    C = Vector{Vector{Real}}(undef, 0)
    for n in 0:(length(moments) - 1)
        push!(C, [(-1)^k * prod(2*m+n-1:-1:1) * moments[n-2*k+1] / (prod(n-2*k:-1:1) *
            prod(2*k:-2:1) * prod(2*m+2*k-1:-2:1))  for k in 0:n÷2])
    end
    return C
end


# Lorentz transform within function support or if larger than ground,
# moment expansion otherwise
function _evaluatelorentztransform(L::LorentzTransform, y::Real, z::Complex)
    lt = lorentz_transform(L.f, L.m, y, z)
    abs(y - real(z)) <= max(abs.(support(L.f))...) && return lt
    ! isnan(lt) && abs(lt) > L.ground && return lt
    a, series, n, term, Y = imag(z), 0.0, 0, L.C[1][1], float(y - real(z))
    n += 1
    term += Y^-n * sum([L.C[n+1][k+1] * a^(2*k) for k in 0:n÷2])
    while ! isnan(term) && ! (series + term ≈ series) && n <= length(L.C) - 3
        series += term
        n += 1
        term = Y^-n * sum([L.C[n+1][k+1] * a^(2*k) for k in 0:n÷2])
        n += 1
        term += Y^-n * sum([L.C[n+1][k+1] * a^(2*k) for k in 0:n÷2])
    end
#    series + term ≈ series || @warn "Moment expansion not converged at (y, z) = "*
#        "($(y), $(z))."
    return (-a / (2 * π))^L.m * 2 / prod(L.m-1:-1:1) / Y^(2 * L.m) * series
end
