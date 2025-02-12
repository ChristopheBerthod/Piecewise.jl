
# Hypergeometric function 2F1, taking the value below the cut
_H2F1(a::Real, b::Real, c::Real, z::Number) =
    imag(z) ≈ 0 ? _₂F₁(a, b, c, z - im * eps(Float64)) : _₂F₁(a, b, c, z)


## POLY: Polynomial

# Function
POLY_F(x::Real, a::Vector{T} where T<:Any) =
    sum(a .* x .^ (0:length(a)-1))

# Scaling
POLY_F(a::Vector{T} where T<:Any, s::Number) = a * s

# Mirroring
POLY_F(a::Vector{T} where T<:Any) = a .* (-1) .^ (0:length(a)-1)

# Moments
POLY_F(x::Real, a::Vector{T} where T<:Any, n::Int) =
    sum(a ./ (n+1:n+length(a)) .* x .^ (n+1:n+length(a)))


## TAIL: Rational function with 1/x or 1/x^2 decay

# Function
TAIL_F(x::Real, a::Vector{T} where T<:Any) =
    (a[1] + a[2] * x) / (a[3] + a[4] * x + a[5] * x^2)

# Checks
function TAIL_F(a, domain, included, danger, fatal)
    Δ = a[4]^2 - 4 * a[3] * a[5]
    t = Δ < 0 || (!(domain[1] <= (-a[4] + sqrt(Δ)) / (2 * a[5]) <= domain[2]) && 
        !(domain[1] <= (-a[4] - sqrt(Δ)) / (2 * a[5]) <= domain[2]))
    !t && fatal && throw(ArgumentError("Singularities must be outside domain."))
    t || return false
    return true
end

# Scaling
TAIL_F(a::Vector{T} where T<:Any, s::Number) = [a[1] * s, a[2] * s, a[3:5]...]

# Mirroring
TAIL_F(a::Vector{T} where T<:Any) = [a[1], -a[2], a[3], -a[4], a[5]]

# Moments
function TAIL_F(x::Real, a::Vector{T} where T<:Any, n::Int)
    Δ = sqrt(complex(a[4]^2 - 4 * a[3] * a[5]))
    if a[3] == a[4] == 0 && n == 0
        return (-a[1] / x + a[2] * log(abs(x))) / a[5]
    elseif a[3] == a[4] == 0 && n == 1
        return (x * a[2] + a[1] * log(abs(x))) / a[5]
    elseif a[3] == a[4] == 0
        return x^n * (a[1]/((n - 1) * x) + a[2] / n) / a[5]
    elseif a[3] == 0 && n == 0
        return (a[1] * log(abs(x)) + (a[2] * a[4] / a[5] - a[1]) *
            log(abs(a[4] + a[5] * x))) / a[4]
    elseif a[3] == 0
        return x^n * real((n + 1) * a[1] * a[4] + n * (a[2] * a[4] - a[1] * a[5])
            * x * _H2F1(1, n + 1, n + 2, -a[5] * x / a[4])) / (n * (n + 1) * a[4]^2)
    elseif a[4] == a[5] == 0
        return (x^(n + 1) * a[1]  / (n + 1) + x^(n + 2) * a[2]  / (n + 2)) / a[3]
    elseif a[4] + Δ ≈ 0 || a[4] - Δ ≈ 0
        return x^(n + 1) * real(a[2] * a[3] + (a[1] * a[4] - a[2] * a[3]) *
            _H2F1(1, n + 1, n + 2, -a[4] * x / a[3])) / ((n + 1) * a[3] * a[4])
    elseif Δ == 0
        return x^(n + 1) * real((a[1] * a[4] - 2 * a[2] * a[3]) /
            (1 + a[4] * x / (2 * a[3])) - (n * a[1] * a[4] / (n + 1)
            - 2 * a[2] * a[3]) * _H2F1(1, n + 1, n + 2, -a[4] * x / (2 * a[3]))) /
            (a[3] * a[4])
    else
        return x^(n+1) * real(((2 * a[2] * a[3] - a[1] * (a[4] - Δ))
            * _H2F1(1, n + 1, n + 2, -2 * a[5] * x / (a[4] + Δ))
            - (2 * a[2] * a[3] - a[1] * (a[4] + Δ))
            * _H2F1(1, n + 1, n + 2, -2 * a[5] * x / (a[4] - Δ))) /
            (2 * (n + 1) * a[3] * Δ))
    end
end


## LOG: Logarithmic singularity

# Function
LOG_F(x::Real, a::Vector{T} where T<:Any) = a[2] * log(abs(x - a[1]))

# Checks
function LOG_F(a, domain, included, danger, fatal)
    t = !(domain[1] < a[1] < domain[2])
    !t && fatal && throw(ArgumentError("Singularity must be outside domain."))
    t || return false
    t = a[1] != domain[1] || ! included[1]
    !t && fatal && throw(ArgumentError("Use (false, $(included[2])) to exclude "*
        "the singularity."))
    t || return false
    t = a[1] != domain[2] || ! included[2]
    !t && fatal && throw(ArgumentError("Use ($(included[1]), false) to exclude "*
        "the singularity."))
    t || return false
    return true
end

# Moments
function LOG_F(x::Real, a::Vector{T} where T<:Any, n::Int)
    if a[1] == 0
        return -a[2] * x^(n + 1) * (1 / (n + 1) - log(abs(x))) / (n + 1)
    else
        return a[2] * x^(n + 1) * real((n + 2) * a[1] * log(abs(x - a[1])) +
            x * _H2F1(1, n + 2, n + 3, x / a[1])) / ((n + 1) * (n + 2) * a[1])
    end
end


## ISRS: Inverse square-root singularity

# Function
ISRS_F(x::Real, a::Vector{T} where T<:Any) = a[2] / sqrt(abs(x^2 - a[1]^2))

# Checks
function ISRS_F(a, domain, included, danger, fatal)
    t = !(domain[1] < a[1] < domain[2] || domain[1] < -a[1] < domain[2])
    !t && fatal && throw(ArgumentError("Singularity must be outside domain."))
    t || return false
    t = abs(a[1]) != abs(domain[1]) || ! included[1]
    !t && fatal && throw(ArgumentError("Use (false, $(included[2])) to exclude "*
        "the singularity."))
    t || return false
    t = abs(a[1]) != abs(domain[2]) || ! included[2]
    !t && fatal && throw(ArgumentError("Use ($(included[1]), false) to exclude "*
        "the singularity."))
    t || return false
    return true
end

# Moments
function ISRS_F(x::Real, a::Vector{T} where T<:Any, n::Int)
    if a[1] == 0 && n == 0
        return a[2] * sign(x) * log(abs(x))
    elseif a[1] == 0
        return a[2] * x^(n + 1) / (n * abs(x))
    else
        return a[2] * x^(n + 1) * real(sqrt(complex(1 - (x / a[1])^2)) *
            _H2F1(1/2, (n + 1)/2, (n + 3)/2, (x / a[1])^2)) /
            ((n + 1) * sqrt(abs(x^2 - a[1]^2)))
    end
end


## PLS: Power-law singularity

# Function
PLS_F(x::Real, a::Vector{T} where T<:Any) = a[3] * abs(x - a[1])^a[2]

# Checks
function PLS_F(a, domain, included, danger, fatal)
    t = !(domain[1] < a[1] < domain[2])
    !t && fatal && throw(ArgumentError("Singularity must be outside domain."))
    t || return false
    t = a[2] >= 0 || (a[1] != domain[1] || ! included[1])
    !t && fatal && throw(ArgumentError("Use (false, $(included[2])) to exclude "*
        "the singularity."))
    t || return false
    t = a[2] >= 0 || (a[1] != domain[2] || ! included[2])
    !t && fatal && throw(ArgumentError("Use ($(included[1]), false) to exclude "*
        "the singularity."))
    t || return false
    t = -12 <= a[2] <= 12
    !t && fatal && throw(ArgumentError("Exponent must be in the range [-12, 12]."))
    t || return false
    return true
end

# Moments
function PLS_F(x::Real, a::Vector{T} where T<:Any, n::Int)
    if a[1] == 0 && a[2] == -(n + 1)
        return a[3] * log(abs(x)) * sign(x)^(n + 1)
    elseif a[1] == 0
        return a[3] * abs(x)^a[2] * x^(n + 1) / (n + 1 + a[2])
    else
        return a[3] * abs(x - a[1])^a[2] * x^(n + 1) * real(
            complex(1 - x / a[1])^(-a[2]) *
            _H2F1(n + 1, -a[2], n + 2, x / a[1])) / (n + 1)
    end
end


## XLOG: Logarithmic singularity times x

# Function
XLOG_F(x::Real, a::Vector{T} where T<:Any) = a[2] * x * log(abs(x - a[1]))

# Checks
XLOG_F(a, domain, included, danger, fatal) =
    LOG_F(a, domain, included, danger, fatal)

# Mirroring
XLOG_F(a::Vector{T} where T<:Any) = [-a[1], -a[2]]

# Moments
function XLOG_F(x::Real, a::Vector{T} where T<:Any, n::Int)
    if a[1] == 0
        return -a[2] * x^(n + 2) * (1 / (n + 2) - log(abs(x))) / (n + 2)
    else
        return a[2] * x^(n + 2) * real((n + 3) * a[1] * log(abs(x - a[1])) +
            x * _H2F1(1, n + 3, n + 4, x / a[1])) / ((n + 2) * (n + 3) * a[1])
    end
end


## XISRS: Inverse square-root singularity times x

# Function
XISRS_F(x::Real, a::Vector{T} where T<:Any) = a[2] * x / sqrt(abs(x^2 - a[1]^2))

# Checks
XISRS_F(a, domain, included, danger, fatal) =
    ISRS_F(a, domain, included, danger, fatal)

# Mirroring
XISRS_F(a::Vector{T} where T<:Any) = [-a[1], -a[2]]

# Moments
function XISRS_F(x::Real, a::Vector{T} where T<:Any, n::Int)
    if a[1] == 0
        return a[2] * x^(n + 2) / ((n + 1) * abs(x))
    else
        return a[2] * x^(n + 2) * real(sqrt(complex(1 - (x / a[1])^2)) *
            _H2F1(1/2, n/2 + 1, n/2 + 2, (x / a[1])^2)) /
            ((n + 2) * sqrt(abs(x^2 - a[1]^2)))
    end
end

