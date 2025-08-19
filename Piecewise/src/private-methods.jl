

# Evaluate the rule of a piece
_evaluatepiece(x::Real, piece::Piece) = sum([F.value(x, a)
    for (F, a) in zip(piece.rule, piece.parameters)])


# Evaluate the primitive of a piece
function _evaluateprimitive(piece::Piece, parity::Symbol, X::Any)
    # Avoid the domain boundaries if they are excluded
    xₘᵢₙ = piece.domain[1]; ! piece.included[1] && (xₘᵢₙ += eps(Float64))
    xₘₐₓ = piece.domain[2]; ! piece.included[2] && (xₘₐₓ -= eps(Float64))
    if parity === :none
        return sum([F.value(xₘₐₓ, a, X) - F.value(xₘᵢₙ, a, X)
            for (F, a) in zip(piece.rule, piece.parameters)])
    else
        s = parity === :even ? 1 : -1
        return sum([F.value(s, xₘₐₓ, a, X) - F.value(s, xₘᵢₙ, a, X)
            for (F, a) in zip(piece.rule, piece.parameters)])
    end
end


# Evaluate the moment of order n of a piece
function _evaluatemoment(piece::Piece, parity::Symbol, n::Int)
    # Avoid the domain boundaries if they are excluded
    # We use 100*eps(Float64) here, otherwise LOG_F and XLOG_F are unstable
    xₘᵢₙ = piece.domain[1]; ! piece.included[1] && (xₘᵢₙ += 100*eps(Float64))
    xₘₐₓ = piece.domain[2]; ! piece.included[2] && (xₘₐₓ -= 100*eps(Float64))
    if parity === :none
        return sum([F.value(xₘₐₓ, a, n) - F.value(xₘᵢₙ, a, n)
            for (F, a) in zip(piece.rule, piece.parameters)])
    else
        # This is only correct for nonvanishing moments; moments that vanish by
        # symmetry are not evaluated here (see moment method)
        return 2 * sum([F.value(xₘₐₓ, a, n) - F.value(xₘᵢₙ, a, n)
            for (F, a) in zip(piece.rule, piece.parameters)])
    end
end


# Evaluate a piecewise function
function _evaluatepiecewise(f::PiecewiseFunction, x::Real)
    f.parity === :even && x < 0 &&
        return  _evaluatepiecewise(f, -x)
    f.parity === :odd  && x < 0 &&
        return -_evaluatepiecewise(f, -x)
    # Find the relevant domain
    i = 1
    while true
        # x before left edge of first domain or in a gap
        x < f.pieces[i].domain[1] && return 0
        # x at included left edge of domain
        x == f.pieces[i].domain[1] && f.pieces[i].included[1] && break
        # x at excluded left edge of domain
        if x == f.pieces[i].domain[1] && ! f.pieces[i].included[1]
            # case x = 0 for even or odd function
            x == 0.0 && f.parity === :even && return Inf *
                sign(_evaluatepiece(x + 10 * eps(Float64), f.pieces[i]))
            x == 0.0 && f.parity === :odd && return 0
            # case of first domain
            i == 1 && return 0
            # case with a gap at the left of piece[i]
            i > 1 && f.pieces[i - 1].domain[2] != x && return 0
            # If we get here, x is on a singularity
            return Inf * (sign(_evaluatepiece(x - 10 * eps(Float64), f.pieces[i-1]))
                + sign(_evaluatepiece(x + 10 * eps(Float64), f.pieces[i])))
        end
        # x in domain of piece[i]
        f.pieces[i].domain[1] < x < f.pieces[i].domain[2] && break
        # x at the included right edge of domain
        x == f.pieces[i].domain[2] && f.pieces[i].included[2] && break
        # x at the excluded right edge of domain or beyond
        i == length(f.pieces) && return 0
        i += 1
    end
    return _evaluatepiece(x, f.pieces[i])
end


# Return a piecewise function of :none parity that is equivalent to the
# piecewise function f
function _unfold(f::PiecewiseFunction)
    f.parity === :none && return f
    if f.parity === :even
        return PiecewiseFunction(f.pieces) + PiecewiseFunction([
            Piece((-p.domain[2], -p.domain[1]), reverse(p.included), p.rule,
            Vector{Vector{Any}}([F.mirror(a)
            for (F, a) in zip(p.rule, p.parameters)])) for p in f.pieces])
    else
        return PiecewiseFunction(f.pieces) + PiecewiseFunction([
            Piece((-p.domain[2], -p.domain[1]), reverse(p.included), p.rule,
            Vector{Vector{Any}}([F.scale(F.mirror(a), -1)
            for (F, a) in zip(p.rule, p.parameters)])) for p in f.pieces])
    end
end


# Return 0 if a domain is not entirely contained in a piece of f,
# otherwise the piece number
function _in(f::PiecewiseFunction, domain::Tuple{Real, Real})
    for (i, p) in enumerate(f.pieces)
        p.domain[1] <= domain[1] && domain[2] <= p.domain[2] && return i
    end
    return 0
end


# Determine whether a boundary is included in the merge of f1 and f2
function _in(f1::PiecewiseFunction, f2::PiecewiseFunction, boundary::Real,
    side::Symbol)
    # Work in progress: true by default
    # We just enumerate the cases where it should be false
    # Boundary excluded from any of the two functions
    length(findall(x -> (x[1] == boundary && ! x[2]),
        vcat([[(p.domain[1], p.included[1]), (p.domain[2], p.included[2])]
        for p in union(f1.pieces, f2.pieces)]...))) > 0 && return false
    # Boundary on a singularity of any of the two functions
    (abs(f1(boundary)) == Inf || abs(f2(boundary)) == Inf) && return false
    return true
end


# Return a method for evaluating a combination of formula for a list of x values
function _evaluatecombination(C::Vector{Formula})
    if length(C) == 1
        return (x, a) -> [C[1].value(x, a) for x in x]
    else
        k = vcat(0, [sum([abs(C[j].params) for j=1:i]) for i=1:length(C)-1])
        return (x, a) -> [(
            sum([F.value(x, a[k[i]+1:k[i + 1]]) for (i, F) in enumerate(C[1:end-1])])
            + C[end].value(x, a[k[end]+1:length(a)])
            ) for x in x]
    end
end


# Return a [Piece] object in case of successful fit or split the domain
# and call itself recursively in each sub-domain
function _fitpieceorsplitdomain(f::Function, domain::Tuple{Real, Real},
    included::Tuple{Bool, Bool}, combinations::Vector{Vector{Formula}},
    grain::Real, resolution::Real, rtol::Real, atol::Real, loop::Bool,
    pairs::Vector{Tuple{Real, Real}})

    # Empty arrays for the winning combination
    combination, parameters = Vector{Formula}[], Array{Vector{Any}}(undef, 0)

    # Array for the sorted list of already computed (x, y) pairs
    memory = pairs[:]

    # If it is the first pass in this domain and if the domain size is larger
    # than the resolution, we precompute values to take advantage of the threads
    if isempty(memory) && domain[2] - domain[1] > resolution
        n = round(Int, (domain[2] - domain[1]) / resolution)
        resize!(memory, n)
        @threads for i = 1:n
            x = domain[1] + i * (domain[2] - domain[1]) / (n + 1)
            memory[i] = (x, f(x))
        end
    end

    # Try all combinations of formula in turn
    Success = false
    for C in combinations

        if C[end].params >= 0
            # Case of fixed number of parameters
            nₘᵢₙ = nₘₐₓ = sum([F.params for F in C])
        else
            # Case of variable number of parameters
            nₘᵢₙ = sum([F.params for F in C[1:end-1]]) + 1
            nₘₐₓ = nₘᵢₙ + abs(C[end].params) - 1
        end

        for n in nₘᵢₙ:nₘₐₓ
            # Fit _evaluatecombination(C)(x, a[1:n]) to f(x) in domain
            ok, a = _fitcombination!(f, domain, _evaluatecombination(C), n,
                resolution, rtol, atol, memory)
            if ok
                # Build array of parameters
                for F in C
                    if F.params >= 0
                        append!(parameters, [a[1:F.params]])
                        deleteat!(a, 1:F.params)
                    else
                        append!(parameters, [a[1:end]])
                    end
                end
                # Exit if fitted parameters are acceptable
                if all([F.check(parameters[i], domain, included, false, false)
                    for (i, F) in enumerate(C)])
                    combination = C
                    Success = true
                    break
                else
                    empty!(parameters)
                end
            end
        end

        Success && break
    end

    if Success
        return [Piece(domain, included, combination, parameters)]
    elseif domain[2] - domain[1] > grain
        # Launch new tasks in two different threads
        p1 = @spawn vcat(
            _fitpieceorsplitdomain(f, (domain[1], (domain[1] + domain[2]) / 2),
            (included[1], true), combinations, grain, resolution, rtol, atol,
            loop, memory[findall(
            x -> domain[1] <= x[1] <= (domain[1] + domain[2]) / 2, memory)])...)
        p2 = @spawn vcat(
            _fitpieceorsplitdomain(f, ((domain[1] + domain[2]) / 2, domain[2]),
            (true, included[2]), combinations, grain, resolution, rtol, atol,
            loop, memory[findall(
            x -> (domain[1] + domain[2]) / 2 <= x[1] <= domain[2], memory)])...)
        return vcat([fetch(p1), fetch(p2)]...)
    else
        @warn "Reached discretisation limit of $(grain) in domain $(domain)"
        if grain > eps(Float64)
            try
                return _bestfit(f, domain, included, combinations, memory)
            catch
                if loop
                    return [Piece(domain, included, x -> f(x))]
                else
                    # Return a linear interpolation across the domain
                    return [Piece(domain, included,
                        POLY, (domain[1], domain[2], f(domain[1] + eps(Float64)),
                        f(domain[2] - eps(Float64))) |> x -> [
                        x[3] + (x[4] - x[3]) / (x[2] - x[1]) * (-x[1]),
                        (x[4] - x[3]) / (x[2] - x[1])])]
                end
            end
        else
            if loop
                return [Piece(domain, included, x -> f(x))]
            else
                    return [Piece(domain, included,
                        POLY, (domain[1], domain[2], f(domain[1] + eps(Float64)),
                        f(domain[2] - eps(Float64))) |> x -> [
                        x[3] + (x[4] - x[3]) / (x[2] - x[1]) * (-x[1]),
                        (x[4] - x[3]) / (x[2] - x[1])])]
            end
        end
    end
end


# Return lists x[] and y[] for at least n values of f(x) in domain
# The new values are stored in the sorted array memory
function _values!(f::Function, domain::Tuple{Real, Real}, n::Int,
    memory::Vector{Tuple{Real, Real}})

    # Case when there are enough values in memory
    length(memory) >= n && return getfield.(memory, 1), getfield.(memory, 2)

    # Case of empty memory
    if isempty(memory)
        resize!(memory, n)
        @threads for i = 1:n
            x = domain[1] + i * (domain[2] - domain[1]) / (n + 1)
            memory[i] = (x, f(x))
        end
        return getfield.(memory, 1), getfield.(memory, 2)
    end

    # Case of missing values
    while length(memory) < n
        # Add a new (x, y) pair in the largest interval
        dx = vcat(getfield.(memory, 1), domain[2]...) .-
            vcat(domain[1], getfield.(memory, 1)...)
        i = argmax(dx)
        x = domain[1] + sum(dx[1:i-1]) + dx[i] / 2
        insert!(memory, i, (x, f(x)))
    end
    return getfield.(memory, 1), getfield.(memory, 2)

end


# Success criterion for the fit
_success(y, dy, rtol, atol) = all(@. abs(dy) < rtol * abs(y) + atol)


# Return true and parameters from fitting model to the function f in domain
# within tolerance. Return (false, Nothing) if unsuccessful
function _fitcombination!(f::Function, domain::Tuple{Real, Real},
    model::Function, n::Int, resolution::Real, rtol::Real, atol::Real,
    memory::Vector{Tuple{Real, Real}})

    # Todo: the case +-Inf needs a special treatment here
    (domain[1] == -Inf || domain[2] == Inf) && throw(ArgumentError(
        "Fit in unbounded domains not yet implemented"))
    # Prepare at least as many data points as there are parameters
    x, y = _values!(f, domain, n, memory)

    # Name for the LsqFit.LsqFitResult object
    fit = 0

    # Try to fit at most n times with random initial parameters
    Success = false
    for i = 1:n
        try
            fit = curve_fit(model, x, y, 2 * rand(n) .- 1)
        catch
            continue
        end
        if _success(y, fit.resid, rtol, atol)
            Success = true
            break
        end
    end
    # If unsuccessful here, there is no need to try further
    Success || return false, Nothing

    # Add points near the boundaries to (try and) track singularities and
    # enforce continuity (unless they already exist)
    if memory[1][1] - domain[1] > 2 * eps(Float32)
        x = domain[1] + eps(Float32)
        insert!(memory, 1, (x, f(x)))
    end
    if domain[2] - memory[end][1] > 2 * eps(Float32)
        x = domain[2] - eps(Float32)
        insert!(memory, length(memory) + 1, (x, f(x)))
    end
    # Add further points at intermediate positions and fit again
    x, y = _values!(f, domain, max(2 * n + 3,
        round(Int, (domain[2] - domain[1]) / resolution)), memory)
    try
        fit = curve_fit(model, x, y, fit.param)
    catch
        return false, Nothing
    end
    if _success(y, fit.resid, rtol, atol)
        return true, fit.param
    else
        return false, Nothing
    end

end


# Return a [Piece] object with the best fit or throw an error if unsuccessful
function _bestfit(f::Function, domain::Tuple{Real, Real},
    included::Tuple{Bool, Bool}, combinations::Vector{Vector{Formula}},
    memory::Vector{Tuple{Real, Real}})

    # Empty arrays for the winning combination
    combination, parameters = Vector{Formula}[], Array{Vector{Any}}(undef, 0)

    # Try all combinations of formula in turn
    Err = Inf
    for C in combinations

        if C[end].params >= 0
            # Case of fixed number of parameters
            nₘᵢₙ = nₘₐₓ = sum([F.params for F in C])
        else
            # Case of variable number of parameters
            nₘᵢₙ = sum([F.params for F in C[1:end-1]]) + 1
            nₘₐₓ = nₘᵢₙ + abs(C[end].params) - 1
        end

        for n in nₘᵢₙ:nₘₐₓ
            # Fit _evaluatecombination(C)(x, a[1:n]) to f(x) in domain
            err, a = Err, Array{Vector{Any}}(undef, 0)
            try
                err, a = _bestfitcombination!(f, domain,
                    _evaluatecombination(C), n, nₘₐₓ, memory)
            catch
                continue
            end
            if err < Err
                # Build array of parameters
                parameters = Array{Vector{Any}}(undef, 0)
                for F in C
                    if F.params >= 0
                        append!(parameters, [a[1:F.params]])
                        deleteat!(a, 1:F.params)
                    else
                        append!(parameters, [a[1:end]])
                    end
                end
                # Save if fitted parameters are acceptable
                if all([F.check(parameters[i], domain, included, false, false)
                    for (i, F) in enumerate(C)])
                    combination = C
                    Err = err
                end
            end
        end

    end

    Err < Inf || throw(ErrorException("All fit attempts failed."))

    return [Piece(domain, included, combination, parameters)]
end


# Return error and parameters from fitting model to the function f in domain
# or throw an error if unsuccessful
function _bestfitcombination!(f::Function, domain::Tuple{Real, Real},
    model::Function, n::Int, np::Int, memory)

    # Todo: the case +-Inf needs a special treatment here
    (domain[1] == -Inf || domain[2] == Inf) && throw(ArgumentError(
        "Fit in unbounded domains not yet implemented"))
    # Prepare np data points (same number of points for all concurrent fits)
    x, y = _values!(f, domain, np, memory)

    # Names for the LsqFit.LsqFitResult and LsqFit.LsqFitResult.param objects
    fit, a = 0, 0

    # Try to fit at most n times with random initial parameters
    Err = Inf
    for i = 1:n
        try
            fit = curve_fit(model, x, y, 2 * rand(n) .- 1)
        catch
            continue
        end
        err = max(abs.(fit.resid)...)
        if err < Err
            a = fit.param
            Err = err
        end
    end
    Err < Inf || throw(ErrorException("All fit attempts failed."))

    # Add points near the boundaries and at intermediate positions and fit again
    if memory[1][1] - domain[1] > 2 * eps(Float32)
        x = domain[1] + eps(Float32)
        insert!(memory, 1, (x, f(x)))
    end
    if domain[2] - memory[end][1] > 2 * eps(Float32)
        x = domain[2] - eps(Float32)
        insert!(memory, length(memory) + 1, (x, f(x)))
    end
    # Add further points at intermediate positions and fit again
    x, y = _values!(f, domain, 2 * np + 3, memory)
    try
        fit = curve_fit(model, x, y, a)
    catch
        throw(ErrorException("All fit attempts failed."))
    end
    return max(abs.(fit.resid)...), fit.param

end


# Return a string with the constructor of p
function _format(p::Piece; indent::Int=0)
    # First line
    s = " "^indent * "Piece($(p.domain), $(p.included), ["
    # Formulas are replaced by their name, if any
    for (j, F) in enumerate(p.rule)
        j > 1 && (s *= ", ")
        s *= (F.name == "" ? F.constructor : F.name)
    end
    s *= "],\n"
    # Generate a string with all formatted parameters (indented)
    params = " "^(indent + 4) * "["
    for (j, a) in enumerate(p.parameters)
        j > 1 && (params *= ", ")
        params *= "["
        for (k, n) in enumerate(a)
            k > 1 && (params *= ", ")
            params *= @sprintf("%.15e", a[k])
        end
        params *= "]"
    end
    params *= "])"
    # Break into 81-characters lines (typically three numbers per line)
    while length(params) > 81
        cut = findlast(' ', params[1:81])
        s *= params[1:cut-1] * "\n"
        params = " "^(indent + 4) * params[cut+1:end]
    end
    s *= params
end


# Return a string with the constructor of f
function _format(f::PiecewiseFunction)
    # First line
    s = "PiecewiseFunction("
    f.parity == :none ? (s *= "[\n") : (s *= ":$(f.parity), [\n")
    # Loop over pieces
    for (i, p) in enumerate(f.pieces)
        s *= _format(p, indent=4)
        i < length(f.pieces) && (s *= ",\n")
    end
    # Last line
    s *= "\n])"
end
