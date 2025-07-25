using Piecewise, PiecewiseHilbert
using Test

@testset "HilbertTransform" begin

    # Initialization with or without radius
    f = PiecewiseFunction(Piece((0, 1), POLY, [1]))
    @test HilbertTransform == typeof(HilbertTransform(f, 10))
    @test HilbertTransform == typeof(HilbertTransform(f))

    # Evaluation on the real axis and moment expansion
    @test_throws DomainError HilbertTransform(f)(0im)
    @test HilbertTransform(f)(5 * im) ≈ -0.019610356906666668 - 0.19739556165079367im

    # Initialization with function lacking primitive
    f = PiecewiseFunction(Piece((0, 1), x -> x))
    @test_throws ArgumentError HilbertTransform(f, 10)
    @test_throws ArgumentError HilbertTransform(f)

end

@testset "Print and show HilbertTransform" begin

    H = HilbertTransform(PiecewiseFunction(Piece((0, 1), POLY, [1])))
    @test typeof(show(devnull, "text/plain", H)) == Int64
    @test typeof(print(H)) == Nothing

    # HilbertTransform objects are not iterable
    @test typeof(show.([devnull, devnull], "text/plain", H)) == Vector{Int64}
    
end

@testset "Formula" begin

    Hp(p) = HilbertTransform(PiecewiseFunction(p))
    He(p) = HilbertTransform(PiecewiseFunction(:even, p))

    p = Piece((0, 1), POLY, [1, 2])
    @test Hp(p)(im) ≈ -0.7757772634850759 - 1.4785453439573946im
    @test He(p)(im) ≈ 0.0 - 2.957090687914789im

    p = Piece((0, 1), TAIL, [1, 2, 3, 4, 5])
    @test Hp(p)(im) ≈ -0.10423268685047321 - 0.25063530013595703im
    @test He(p)(im) ≈ 0.0 - 0.5012706002719141im

    p = Piece((0, 1), LOG, [2, 1])
    @test Hp(p)(im) ≈ -0.10648681180992137 - 0.3350927975710153im
    @test He(p)(im) ≈ 0.0 - 0.6701855951420306im

    p = Piece((0, 1), ISRS, [2, 1])
    @test Hp(p)(im) ≈ -0.184214856111667 - 0.4077417592590041im
    @test He(p)(im) ≈ 0.0 - 0.8154835185180082im

    p = Piece((0, 1), PLS, [2, -1/2, 1])
    @test Hp(p)(im) ≈ -0.29835871733798514 - 0.6373923908392799im
    @test He(p)(im) ≈ 0.0 - 1.2747847816785598im

    p = Piece((0, 1), XLOG, [2, 1])
    @test Hp(p)(im) ≈ -0.05120156354887517 - 0.10648681180992137im
    @test He(p)(im) ≈ 0.0 - 0.21297362361984273im

    p = Piece((0, 1), XISRS, [2, 1])
    @test Hp(p)(im) ≈ -0.11585701633929485 - 0.184214856111667im
    @test He(p)(im) ≈ 0.0 - 0.368429712223334im

end

