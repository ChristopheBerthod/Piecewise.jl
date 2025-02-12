
# Hypergeometric function 2F1, taking the value below the cut
_H2F1(a::Real, b::Real, c::Real, z::Number) =
    imag(z) ≈ 0 ? _₂F₁(a, b, c, z - im * eps(Float64)) : _₂F₁(a, b, c, z)

## Piecewise.POLY

# Primitive
function Piecewise.POLY_F(x::Real, a::Vector{T} where T<:Any, z::Complex)
    if length(a) == 1
        return -a[1] * log(x - z)
    else
#         return -log(x - z) * sum(a .* z .^ (0:length(a)-1)) -
#             sum([z^i * sum([a[j + 2] * x^(j-i+1) / (j-i+1) for j=i:length(a)-2])
#             for i=0:length(a)-2])

# The primitive commented above leads to numerical noise when z is large compared
# with xₘₐₓ - xₘᵢₙ but still beyond the radius: we face the difference of two
# nearly identical large numbers (see function _evaluateprimitive() in the module
# Piecewise). The following minimal code illustrates the problem:

#   include.(["src/Piecewise.jl", "src/PiecewiseHilbert.jl"])
#   using .Piecewise, .PiecewiseHilbert, DelimitedFiles
#   f = HilbertTransform(
#       PiecewiseFunction([Piece((0, 1), POLY, [0, 0, 0, 0, 0, 0, 0, 1])]), 100)
#   writedlm("noise.txt", [hcat(x, imag(f(x + im))...) for x in 0:0.1:100])

# The following primitive is stable at large z:
        return sum([a[i] * x^i * _H2F1(1, i, i + 1, x / z) / i
            for i=1:length(a)]) / z
    end
end

# Even and odd primitives
Piecewise.POLY_F(s::Int, x::Real, a::Vector{T} where T<:Any, z::Complex) =
    Piecewise.POLY_F(x, a, z) - s * Piecewise.POLY_F(x, a, -z)


## Piecewise.TAIL

# Primitive
function Piecewise.TAIL_F(x::Real, a::Vector{T} where T<:Any, z::Complex)
    Δ = sqrt(complex(a[4]^2 - 4 * a[3] * a[5]))
    return -((a[1] + a[2] * z) * (log(z - x) - log(sqrt(complex(a[3] + a[4] * x
        + a[5] * x^2)))) - (2 * a[2] * a[3] - a[1] * a[4] + (a[2] * a[4] - 2 *
        a[1] * a[5]) * z) * atanh((a[4] + 2 * a[5] * x) / Δ) / Δ) /
        (a[3] + a[4] * z + a[5] * z^2)
end

# Even and odd primitives
Piecewise.TAIL_F(s::Int, x::Real, a::Vector{T} where T<:Any, z::Complex) =
    Piecewise.TAIL_F(x, a, z) - s * Piecewise.TAIL_F(x, a, -z)


## Piecewise.LOG

# Primitive
Piecewise.LOG_F(x::Real, a::Vector{T} where T<:Any, z::Complex) =
    -a[2] * (log(abs(x - a[1])) * log((z - x) / (z - a[1])) +
    li2((x - a[1]) / (z - a[1])))

# Even and odd primitives
Piecewise.LOG_F(s::Int, x::Real, a::Vector{T} where T<:Any, z::Complex) =
    Piecewise.LOG_F(x, a, z) - s * Piecewise.LOG_F(x, a, -z)


## Piecewise.ISRS

# Primitive
function Piecewise.ISRS_F(x::Real, a::Vector{T} where T<:Any, z::Complex)
    s1, s2 = sqrt(complex(x^2 - a[1]^2)), sqrt(z^2 - a[1]^2)
    return -a[2] / sqrt(abs(x^2 - a[1]^2)) * s1 / s2 *
        atanh((a[1]^2 - x * z)/(s1 * s2))
end

# Even and odd primitives
Piecewise.ISRS_F(s::Int, x::Real, a::Vector{T} where T<:Any, z::Complex) =
    Piecewise.ISRS_F(x, a, z) - s * Piecewise.ISRS_F(x, a, -z)


## Piecewise.PLS

# Primitive
function Piecewise.PLS_F(x::Real, a::Vector{T} where T<:Any, z::Complex)
    s = (x - a[1]) / (z - a[1])
    return a[3] * abs(x - a[1])^a[2] / (1 + a[2]) * s *
        _H2F1(1, a[2] + 1, a[2] + 2, s)
end

# Even and odd primitives
Piecewise.PLS_F(s::Int, x::Real, a::Vector{T} where T<:Any, z::Complex) =
    Piecewise.PLS_F(x, a, z) - s * Piecewise.PLS_F(x, a, -z)


## Piecewise.XLOG

# Primitive
Piecewise.XLOG_F(x::Real, a::Vector{T} where T<:Any, z::Complex) =
    a[2] * (x - log(abs(x - a[1])) * (x - a[1] +
    z * log((z - x) / (z - a[1]))) - z * li2((x - a[1]) / (z - a[1])))

# Even and odd primitives
Piecewise.XLOG_F(s::Int, x::Real, a::Vector{T} where T<:Any, z::Complex) =
    Piecewise.XLOG_F(x, a, z) - s * Piecewise.XLOG_F(x, a, -z)


## Piecewise.XISRS

# Primitive
function Piecewise.XISRS_F(x::Real, a::Vector{T} where T<:Any, z::Complex)
    s1, s2 = sqrt(complex(x^2 - a[1]^2)), sqrt(z^2 - a[1]^2)
    return -a[2] * s1 / sqrt(abs(x^2 - a[1]^2)) *
        (atanh(x / s1) + z * atanh((a[1]^2 - x * z) / (s1 * s2)) / s2)
end

# Even and odd primitives
Piecewise.XISRS_F(s::Int, x::Real, a::Vector{T} where T<:Any, z::Complex) =
    Piecewise.XISRS_F(x, a, z) - s * Piecewise.XISRS_F(x, a, -z)


# Return the moments of the function f up to order 64, if they are not NaN
function _moments(f::PiecewiseFunction)
    n, moments = 0, Vector{Real}(undef, 0)
    while n <= 64 && ! isnan(moment(f, n))
        push!(moments, moment(f, n))
        n += 1
    end
    return moments
end


# Return the typical absolute error of the moment expansion at radius
function _error(f::PiecewiseFunction, radius::Real, moments::Vector{Real})
    z = radius .* exp.(im * [i * π / 11 for i=1:10])
    return max(abs.([_momentexpansion(z, moments, true)[1] -
        hilbert_transform(f, z) for z in z])...)
end


# Return the moment expansion and a flag indicating convergence
function _momentexpansion(z::Complex, moments::Vector{Real}, silent::Bool)
    n, series, term = 0, complex(0), moments[1]
    # Grouping terms by pairs, as half the moments vanish for even and odd functions
    n += 1
    term += moments[n + 1] * float(z)^-n
    while ! isnan(term) && ! (series + term ≈ series) && n <= length(moments) - 3
        series += term
        n += 1
        term = moments[n + 1] * float(z)^-n
        n += 1
        term += moments[n + 1] * float(z)^-n
    end
    silent || series + term ≈ series || @warn "Moment expansion not converged "*
        "at z = $(z). Try a larger radius."
    return series / z, series + term ≈ series
end


# The ratio of term n to term n-1 in the moment expansion is typically equal
# to n/(n+1)*s/radius, where s is the support of f. Since this is smaller than
# unity for radius > s, the expansion typically converges if radius > s.

# Return the smallest radius at which the moment expansion converges
function _automatic_radius(f::PiecewiseFunction, moments::Vector{Real})
    # We start at 10 * support if support < Inf, at 10 otherwise
    s = max(abs.(support(f))...); R = s < Inf ? 10 * s : 1.0
    # If the expansion converges at the starting radius, we divide it by 1.1,
    # otherwise we multiply by 1.1
    z = R .* exp.(im * [i * π / 11 for i=1:10])
    S = all([_momentexpansion(z, moments, true)[2] for z in z]) ? 1 / 1.1 : 1.1
    # We don't go below twice the support and above 1e3
    while (s < Inf ? 2 * s : 1) < R < 1e3
        all([_momentexpansion(z, moments, true)[2] for z in z]) || break
        R *= S
    end
    return R /= S
end


# Hilbert transform if |z| < radius and moment expansion otherwise
function _evaluatehilberttransform(H::HilbertTransform, z::Complex)
    abs(z) < H.radius && return hilbert_transform(H.f, z)
    return _momentexpansion(z, H.moments, false)[1]
end
