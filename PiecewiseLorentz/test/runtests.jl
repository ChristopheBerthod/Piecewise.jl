using Piecewise, PiecewiseLorentz
using Test

@testset "LorentzTransform" begin

    # Initialization with or without ground
    f = PiecewiseFunction(Piece((0, 1), POLY, [1]))
    @test LorentzTransform == typeof(LorentzTransform(f, 2, 1e-10))
    @test LorentzTransform == typeof(LorentzTransform(f, 2))

    # Evaluation with negative imaginary part and moment expansion
    @test_throws DomainError LorentzTransform(f, 2)(0, im)
    @test LorentzTransform(f, 2, 1)(1, 5 - im) ≈ 0.00023287042221692184

    # Initialization with function lacking primitive
    f = PiecewiseFunction(Piece((0, 1), x -> x))
    @test_throws ArgumentError LorentzTransform(f, 2, 1e-10)
    @test_throws ArgumentError LorentzTransform(f, 2)

end

@testset "Print and show LorentzTransform" begin

    L = LorentzTransform(PiecewiseFunction(Piece((0, 1), POLY, [1])), 2)
    @test typeof(show(devnull, "text/plain", L)) == Int64
    @test typeof(print(L)) == Nothing

    # LorentzTransform objects are not iterable
    @test typeof(show.([devnull, devnull], "text/plain", L)) == Vector{Int64}
    
end

@testset "Formula" begin

    Lp(p, m) = LorentzTransform(PiecewiseFunction(p), m)
    Le(p, m) = LorentzTransform(PiecewiseFunction(:even, p), m)

    p = Piece((0, 1), POLY, [1, 2])
    @test Lp(p, 1)(0, -im) ≈ 0.47063560015265193
    @test Le(p, 1)(0, -im) ≈ 0.9412712003053039
    @test Lp(p, 2)(0, -im) ≈ 0.11577962350472716
    @test Le(p, 2)(0, -im) ≈ 0.23155924700945432
    @test Lp(p, 3)(0, -im) ≈ 0.029656069987218852
    @test Le(p, 3)(0, -im) ≈ 0.059312139974437704
    @test_throws ArgumentError Lp(p, 4)(0, -im)
    @test_throws ArgumentError Le(p, 4)(0, -im)

    p = Piece((0, 1), TAIL, [1, 2, 3, 4, 5])
    @test_throws ArgumentError Lp(p, 1)(0, -im)
    @test_throws ArgumentError Le(p, 1)(0, -im)

    p = Piece((0, 1), LOG, [2, 1])
    @test Lp(p, 1)(0, -im) ≈ 0.1066633502558379
    @test Le(p, 1)(0, -im) ≈ 0.2133267005116758
    @test Lp(p, 2)(0, -im) ≈ 0.030087400117168876
    @test Le(p, 2)(0, -im) ≈ 0.06017480023433775
    @test Lp(p, 3)(0, -im) ≈ 0.008636588104053535
    @test Le(p, 3)(0, -im) ≈ 0.01727317620810707
    @test_throws ArgumentError Lp(p, 4)(0, -im)
    @test_throws ArgumentError Le(p, 4)(0, -im)

    p = Piece((0, 1), ISRS, [2, 1])
    @test_throws ArgumentError Lp(p, 1)(0, -im)
    @test_throws ArgumentError Le(p, 1)(0, -im)

    p = Piece((0, 1), PLS, [2, -1/2, 1])
    @test Lp(p, 1)(0, -im) ≈ 0.2028882993824653
    @test Le(p, 1)(0, -im) ≈ 0.4057765987649306
    @test Lp(p, 2)(0, -im) ≈ 0.05189177586576684
    @test Le(p, 2)(0, -im) ≈ 0.10378355173153368
    @test Lp(p, 3)(0, -im) ≈ 0.013779563659096715
    @test Le(p, 3)(0, -im) ≈ 0.02755912731819343
    @test_throws ArgumentError Lp(p, 4)(0, -im)
    @test_throws ArgumentError Le(p, 4)(0, -im)

    p = Piece((0, 1), XLOG, [2, 1])
    @test Lp(p, 1)(0, -im) ≈ 0.0338958049472908
    @test Le(p, 1)(0, -im) ≈ 0.0677916098945816
    @test Lp(p, 2)(0, -im) ≈ 0.008665178161249475
    @test Le(p, 2)(0, -im) ≈ 0.01733035632249895
    @test Lp(p, 3)(0, -im) ≈ 0.002270885882873601
    @test Le(p, 3)(0, -im) ≈ 0.004541771765747201
    @test_throws ArgumentError Lp(p, 4)(0, -im)
    @test_throws ArgumentError Le(p, 4)(0, -im)

    p = Piece((0, 1), XISRS, [2, 1])
    @test_throws ArgumentError Lp(p, 1)(0, -im)
    @test_throws ArgumentError Le(p, 1)(0, -im)

end

