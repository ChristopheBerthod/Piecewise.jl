
@testset "Print and show piecewise functions" begin

    @test typeof(show(devnull, "text/plain", POLY)) == Int64
    @test typeof(print(POLY)) == Nothing
    p = Piece((0, 1), POLY, [1])
    @test typeof(show(devnull, "text/plain", p)) == Int64
    @test typeof(print(p)) == Nothing
    f = PiecewiseFunction(:even, p)
    @test typeof(show(devnull, "text/plain", f)) == Int64
    @test typeof(print(f)) == Nothing

    # Formula, Piece, and PiecewiseFunction objects are not iterable
    @test typeof(show.([devnull, devnull], "text/plain", POLY)) == Vector{Int64}
    @test typeof(.*([1, 2], p)) == Vector{Piece}
    @test typeof(.*([1, 2], f)) == Vector{PiecewiseFunction}
    
end


@testset "Arithmetic operations on piecewise functions" begin

    # Addition of a scalar to piecewise function
    f = PiecewiseFunction(:even, Piece((0, 1), POLY, [1]))

    @test (f + 1)(1) == 2
    @test (1 + f)(1) == 2
    @test (f - 2)(1) == -1
    @test (2 - f)(1) == 1

    f = PiecewiseFunction(:odd, Piece((0, 1), POLY, [0, 1]))

    @test (f + 1)(-1) == 1
    @test (1 + f)(-1) == 1
    @test (f - 1)(-1) == -1
    @test (1 - f)(-1) == 1

    # Multiplication of piece by a scalar
    @test typeof(Piece((0, 1), POLY, [1]) * 2) == Piece
    @test typeof(2 * Piece((0, 1), POLY, [1])) == Piece
    @test typeof(Piece((0, 1), POLY, [1]) / 2) == Piece

    # Multiplication of piecewise function by a scalar
    @test (f * 2)(-1) == -2
    @test (2 * f)(1) == 2
    @test (f / 2)(-1) == -0.5

    # Addition of two piecewise funtctions
    @test typeof(f + f) == PiecewiseFunction
    @test typeof(f - f) == PiecewiseFunction
    @test (f + PiecewiseFunction(:even, Piece((0, 1), POLY, [1])))(1) == 2
    @test_throws ArgumentError sum(Vector{PiecewiseFunction}(undef, 0))
    @test typeof(sum([f, f])) == PiecewiseFunction

end


@testset "Public methods for piecewise functions" begin

    f = PiecewiseFunction(:even, Piece((0, 1), POLY, [1]))

    @test domains(f) == [(0, 1)]
    @test intervals(f) == ["[0, 1]"]
    @test support(f) == (-1, 1)
    @test support(PiecewiseFunction(Piece((0, 1), POLY, [1]))) == (0, 1)
    @test singularities(f) == []
    @test singularities(PiecewiseFunction(:even,
        Piece((0, 1), (false, true), LOG, [0, 1]))) == [0]
    @test singularities(PiecewiseFunction([Piece((0, 1), (true, false), POLY, [1]),
        Piece((1, 2), (false, true), POLY, [1])])) == [1]
    @test formulas(f) == ["POLY"]

    @test_throws MethodError integraltransform(f, 1)
    F(x, a) = a[1] + a[2] * x
    LIN = Formula(2, F, (a, s) -> a * s, a -> [a[1], -a[2]])
    F(x, a, k) = (a[1] * k + a[2] * (k * x + im)) * exp(im * k * x)/(im * k^2)
    F(s, x, a, k) = 2 * (s == 1 ? real(F(x, a, k)) : im * imag(F(x, a, k)))
    g = PiecewiseFunction(:even, [Piece((0, π), LIN, [π, -1])])
    @test integraltransform(g, 1) == 4.0
    g = PiecewiseFunction([Piece((0, π), LIN, [π, -1])])
    @test integraltransform(g, 1) == 2.0 + 3.141592653589793im

    @test_throws DomainError moment(f, -1)
    @test moment(f, 1) == 0.0
    @test moment(PiecewiseFunction(:odd, Piece((0, 1), POLY, [0, 1])), 2) == 0.0
    @test moment(f, 0) == 2.0
    @test moment(PiecewiseFunction(Piece((0, 1), POLY, [0, 1])), 0) == 0.5
    
    @test_throws ArgumentError piecewisefit((x, y) -> x, (0, 1), [POLY])
    @test_throws ArgumentError piecewisefit(x -> "x", (0, 1), [POLY])
    @test_throws ArgumentError piecewisefit(x -> x, (0, -1), [POLY])
    @test_throws ArgumentError piecewisefit(x -> x, (-Inf, 0), [POLY])
    @test_throws ArgumentError piecewisefit(x -> x, (0, Inf), [POLY])
    @test_throws ArgumentError piecewisefit(x -> x, (0, 1), [POLY, POLY])
    @test_throws ArgumentError piecewisefit(x -> x, (0, 1), [Formula("x -> im")])
    @test_throws ArgumentError piecewisefit(x -> x, (0, 1), [POLY], parity=:any)
    @test typeof(piecewisefit(x -> x, (0, 3), [POLY],
        singularities=[1, 4], cuts=[2, 4], resolution=0.1, rtol=1e-1)) == PiecewiseFunction
    @test typeof(piecewisefit(x -> x, (0, 3), [ISRS, LOG], parity=:even,
        singularities=[1], cuts=[2], rtol=1e-1)) == PiecewiseFunction
    @test typeof(piecewisefit(x -> x, (0, 1), [LOG], grain=0.5)) == PiecewiseFunction

    @test typeof(format(f)) == String
    @test typeof(format(PiecewiseFunction(Piece((0, 1), POLY, ones(13))))) == String

end

nothing
