
@testset "POLY" begin

    @test Piecewise.POLY_F(1, [1, 2, 3, 4, 5, 6]) ≈ 21
    @test Piecewise.POLY_F([1, 2, 3], 2) == [2, 4, 6]
    @test Piecewise.POLY_F([1, 2, 3, 4]) == [1, -2, 3, -4]
    @test Piecewise.POLY_F(1, [1, 2, 3, 4, 5, 6], 1) ≈ 4.4071428571428575
    @test typeof(POLY) == Formula

end


@testset "TAIL" begin

    @test Piecewise.TAIL_F(1, [1, 2, 3, 4, 5]) ≈ 0.25
    @test_throws ArgumentError Piecewise.TAIL_F([1, 2, -3, 4, 5], (0, 1), (true, true), true, true)
    @test Piecewise.TAIL_F([1, 2, -3, 4, 5], (0, 1), (true, true), true, false) == false
    @test Piecewise.TAIL_F([1, 2, 3, 4, 5], (0, 1), (true, true), true, true) == true
    @test Piecewise.TAIL_F([1, 2, -3, 4, 5], (0, 1), (true, true), true, false) == false
    @test Piecewise.TAIL_F([1, 2, 3, 4, 5], 2) == [2, 4, 3, 4, 5]
    @test Piecewise.TAIL_F([1, 2, 3, 4, 5]) == [1, -2, 3, -4, 5]
    @test Piecewise.TAIL_F(1, [1, 2, 0, 0, 5], 0) ≈ -0.2
    @test Piecewise.TAIL_F(1, [1, 2, 0, 0, 5], 1) ≈ 0.4
    @test Piecewise.TAIL_F(1, [1, 2, 0, 0, 5], 5) ≈ 0.13
    @test Piecewise.TAIL_F(1, [1, 2, 0, 4, 5], 0) ≈ 0.32958368660043297
    @test Piecewise.TAIL_F(1, [1, 2, 0, 4, 5], 1) ≈ 0.30268837405404053
    @test Piecewise.TAIL_F(1, [1, 2, 3, 4, 0], 1) ≈ 0.2044341744113004
    @test Piecewise.TAIL_F(1, [1, 2, 3, 0, 0], 1) ≈ 0.38888888888888884
    @test Piecewise.TAIL_F(1, [1, 2, 3, 4, 0], 1) ≈ 0.2044341744113004
    @test Piecewise.TAIL_F(1, [1, 2, 1, 2, 1], 1) ≈ 0.4205584583201638
    @test Piecewise.TAIL_F(1, [1, 2, 3, 4, 5], 1) ≈ 0.14729513605358416
    @test typeof(TAIL) == Formula

end


@testset "LOG" begin

    @test Piecewise.LOG_F(0.5, [1, 1]) ≈ -0.6931471805599453
    @test_throws ArgumentError Piecewise.LOG_F([1, 1], (0, 2), (true, true), true, true)
    @test_throws ArgumentError Piecewise.LOG_F([1, 1], (1, 2), (true, true), true, true)
    @test_throws ArgumentError Piecewise.LOG_F([1, 1], (0, 1), (true, true), true, true)
    @test Piecewise.LOG_F([1, 1], (0, 2), (true, true), true, false) == false
    @test Piecewise.LOG_F([1, 1], (1, 2), (true, true), true, false) == false
    @test Piecewise.LOG_F([1, 1], (0, 1), (true, true), true, false) == false
    @test Piecewise.LOG_F([1, 1], (0, 1), (true, false), true, true) == true
    @test Piecewise.LOG_F([1, 1], (1, 2), (false, true), true, true) == true
    @test Piecewise.LOG_F(1, [0, 1], 0) ≈ -1.0
    @test Piecewise.LOG_F(1, [2, 1], 0) ≈ 0.38629436111988963
    @test typeof(LOG) == Formula

end


@testset "ISRS" begin

    @test Piecewise.ISRS_F(0.5, [1, 1]) ≈ 1.1547005383792517
    @test_throws ArgumentError Piecewise.ISRS_F([1, 1], (0, 2), (true, true), true, true)
    @test_throws ArgumentError Piecewise.ISRS_F([1, 1], (1, 2), (true, true), true, true)
    @test_throws ArgumentError Piecewise.ISRS_F([1, 1], (0, 1), (true, true), true, true)
    @test Piecewise.ISRS_F([1, 1], (0, 2), (true, true), true, false) == false
    @test Piecewise.ISRS_F([1, 1], (1, 2), (true, true), true, false) == false
    @test Piecewise.ISRS_F([1, 1], (0, 1), (true, true), true, false) == false
    @test Piecewise.ISRS_F([1, 1], (0, 1), (true, false), true, true) == true
    @test Piecewise.ISRS_F([1, 1], (1, 2), (false, true), true, true) == true
    @test Piecewise.ISRS_F(1, [0, 1], 0) ≈ 0.0
    @test Piecewise.ISRS_F(1, [0, 1], 1) ≈ 1.0
    @test Piecewise.ISRS_F(1, [2, 1], 0) ≈ 0.5235987755982989
    @test typeof(ISRS) == Formula

end


@testset "PLS" begin

    @test Piecewise.PLS_F(0.5, [1, -1/2, 1]) ≈ 1.4142135623730951
    @test_throws ArgumentError Piecewise.PLS_F([1, -1/2, 1], (0, 2), (true, true), true, true)
    @test_throws ArgumentError Piecewise.PLS_F([1, -1/2, 1], (1, 2), (true, true), true, true)
    @test_throws ArgumentError Piecewise.PLS_F([1, -1/2, 1], (0, 1), (true, true), true, true)
    @test_throws ArgumentError Piecewise.PLS_F([2, -13, 1], (0, 1), (true, true), true, true)
    @test Piecewise.PLS_F([1, -1/2, 1], (0, 2), (true, true), true, false) == false
    @test Piecewise.PLS_F([1, -1/2, 1], (1, 2), (true, true), true, false) == false
    @test Piecewise.PLS_F([1, -1/2, 1], (0, 1), (true, true), true, false) == false
    @test Piecewise.PLS_F([2, -13, 1], (0, 1), (true, true), true, false) == false
    @test Piecewise.PLS_F([1, -1/2, 1], (0, 1), (true, false), true, true) == true
    @test Piecewise.PLS_F([1, -1/2, 1], (1, 2), (false, true), true, true) == true
    @test Piecewise.PLS_F(0.5, [0, -2, 1], 1) ≈ -0.6931471805599453
    @test Piecewise.PLS_F(0.5, [0, -2, 1], 2) ≈ 0.5
    @test Piecewise.PLS_F(0.5, [1, -1/2, 1], 2) ≈ 0.05314694696594851
    @test typeof(PLS) == Formula

end


@testset "XLOG" begin

    @test Piecewise.XLOG_F(0.5, [1, 1]) ≈ -0.34657359027997264
    @test_throws ArgumentError Piecewise.XLOG_F([1, 1], (0, 2), (true, true), true, true)
    @test_throws ArgumentError Piecewise.XLOG_F([1, 1], (1, 2), (true, true), true, true)
    @test_throws ArgumentError Piecewise.XLOG_F([1, 1], (0, 1), (true, true), true, true)
    @test Piecewise.XLOG_F(1, [0, 1], 0) ≈ -0.25
    @test Piecewise.XLOG_F(1, [2, 1], 0) ≈ 0.1362943611198904
    @test Piecewise.XLOG_F([1, 1]) == [-1, -1]
    @test Piecewise.XLOG_F(1, [0, 1], 0) ≈ -0.25
    @test Piecewise.XLOG_F(1, [2, 1], 0) ≈ 0.1362943611198904
    @test typeof(XLOG) == Formula

end


@testset "XISRS" begin

    @test Piecewise.XISRS_F(0.5, [1, 1]) ≈ 0.5773502691896258
    @test_throws ArgumentError Piecewise.XISRS_F([1, 1], (0, 2), (true, true), true, true)
    @test_throws ArgumentError Piecewise.XISRS_F([1, 1], (1, 2), (true, true), true, true)
    @test_throws ArgumentError Piecewise.XISRS_F([1, 1], (0, 1), (true, true), true, true)
    @test Piecewise.XISRS_F([1, 1]) == [-1, -1]
    @test Piecewise.XISRS_F(1, [0, 1], 0) ≈ 1.0
    @test Piecewise.XISRS_F(1, [2, 1], 0) ≈ 0.2679491924311225
    @test typeof(XISRS) == Formula

end

nothing
