
@testset "Private methods for piecewise functions" begin

    f = PiecewiseFunction(:even, Piece((0, 1), POLY, [1]))
    @test f(-1) == f(1)
    f = PiecewiseFunction(:odd, Piece((0, 1), POLY, [1]))
    @test f(-1) == -f(1)

    @test f(-2) == 0
    @test f(-1) == -1

    f = PiecewiseFunction(:even, Piece((0, 1), (false, true), POLY, [1]))
    @test f(0) == Inf
    f = PiecewiseFunction(:odd, Piece((0, 1), (false, true), POLY, [1]))
    @test f(0) == 0
    f = PiecewiseFunction(Piece((0, 1), (false, true), POLY, [1]))
    @test f(0) == 0
    f = PiecewiseFunction([Piece((-2, -1), POLY, [1]), Piece((0, 1), POLY, [1])])
    @test f(-0.5) == 0
    f = PiecewiseFunction([Piece((0, 1), (true, false), POLY, [1]),
        Piece((1, 2), (false, true), POLY, [1])])
    @test f(1) == Inf

    @test f(0.5) == 1.0
    @test f(2) == 1.0
    @test f(3) == 0

end

nothing
