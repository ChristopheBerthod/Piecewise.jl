
@testset "Formula" begin

    value  = (x, a) -> x
    check  = (a, b, c, d, e) -> true
    scale  = (a, s) -> x
    mirror = a -> a
    fail   = (a, b, c) -> true 

    # Initialization with name and four functions
    @test Formula == typeof(   Formula("f", 1, value, check, scale, mirror))
    @test_throws ArgumentError Formula("f", 1,  fail, check, scale, mirror)
    @test_throws ArgumentError Formula("f", 1, value,  fail, scale, mirror)
    @test_throws ArgumentError Formula("f", 1, value, check,  fail, mirror)
    @test_throws ArgumentError Formula("f", 1, value, check, scale,   fail)

    # Initialization with name and three functions
    @test Formula == typeof(   Formula("f", 1, value, check, scale        ))
    @test Formula == typeof(   Formula("f", 1, value, check,        mirror))
    @test Formula == typeof(   Formula("f", 1, value,        scale, mirror))
    @test_throws ArgumentError Formula("f", 1,  fail, check, scale        )
    @test_throws ArgumentError Formula("f", 1,  fail, check,        mirror)
    @test_throws ArgumentError Formula("f", 1,  fail,        scale, mirror)
    @test_throws ArgumentError Formula("f", 1, value, check,  fail        )
    @test_throws ArgumentError Formula("f", 1, value,         fail, mirror)
    @test_throws ArgumentError Formula("f", 1, value, check,          fail)
    @test_throws ArgumentError Formula("f", 1, value,        scale,   fail)

    # Initialization with name and two functions
    @test Formula == typeof(   Formula("f", 1, value, check               ))
    @test Formula == typeof(   Formula("f", 1, value,        scale        ))
    @test Formula == typeof(   Formula("f", 1, value,               mirror))
    @test_throws ArgumentError Formula("f", 1,  fail, check               )
    @test_throws ArgumentError Formula("f", 1,  fail,        scale        )
    @test_throws ArgumentError Formula("f", 1,  fail,               mirror)
    @test_throws ArgumentError Formula("f", 1, value,  fail               )

    # Initialization with name and one function
    @test Formula == typeof(   Formula("f", 1, value                      ))
    @test_throws ArgumentError Formula("f", 1,  fail                      )

    # Initialization without name and four functions
    @test Formula == typeof(   Formula(1, value, check, scale, mirror))
    @test_throws ArgumentError Formula(1,  fail, check, scale, mirror)
    @test_throws ArgumentError Formula(1, value,  fail, scale, mirror)
    @test_throws ArgumentError Formula(1, value, check,  fail, mirror)
    @test_throws ArgumentError Formula(1, value, check, scale,   fail)

    # Initialization without name and three functions
    @test Formula == typeof(   Formula(1, value, check, scale        ))
    @test Formula == typeof(   Formula(1, value, check,        mirror))
    @test Formula == typeof(   Formula(1, value,        scale, mirror))
    @test_throws ArgumentError Formula(1,  fail, check, scale        )
    @test_throws ArgumentError Formula(1,  fail, check,        mirror)
    @test_throws ArgumentError Formula(1,  fail,        scale, mirror)
    @test_throws ArgumentError Formula(1, value, check,  fail        )
    @test_throws ArgumentError Formula(1, value,         fail, mirror)
    @test_throws ArgumentError Formula(1, value, check,          fail)
    @test_throws ArgumentError Formula(1, value,        scale,   fail)

    # Initialization without name and two functions
    @test Formula == typeof(   Formula(1, value, check               ))
    @test Formula == typeof(   Formula(1, value,        scale        ))
    @test Formula == typeof(   Formula(1, value,               mirror))
    @test_throws ArgumentError Formula(1,  fail, check               )
    @test_throws ArgumentError Formula(1,  fail,        scale        )
    @test_throws ArgumentError Formula(1,  fail,               mirror)
    @test_throws ArgumentError Formula(1, value,  fail               )

    # Initialization without name and one function
    @test Formula == typeof(   Formula(1, value                      ))
    @test_throws ArgumentError Formula(1,  fail                      )

    value  = "(x, a) -> x"
    fails  = "x -> x"

    # Initialization with name, string, and four functions
    @test Formula == typeof(   Formula("f", 1, value, check, scale, mirror))
    @test_throws ArgumentError Formula("f", 1, fails, check, scale, mirror)
    @test_throws ArgumentError Formula("f", 1, value,  fail, scale, mirror)
    @test_throws ArgumentError Formula("f", 1, value, check,  fail, mirror)
    @test_throws ArgumentError Formula("f", 1, value, check, scale,   fail)

    # Initialization with name, string, and three functions
    @test Formula == typeof(   Formula("f", 1, value, check, scale        ))
    @test Formula == typeof(   Formula("f", 1, value, check,        mirror))
    @test Formula == typeof(   Formula("f", 1, value,        scale, mirror))
    @test_throws ArgumentError Formula("f", 1, fails, check, scale        )
    @test_throws ArgumentError Formula("f", 1, fails, check,        mirror)
    @test_throws ArgumentError Formula("f", 1, fails,        scale, mirror)
    @test_throws ArgumentError Formula("f", 1, value, check,  fail        )
    @test_throws ArgumentError Formula("f", 1, value,         fail, mirror)
    @test_throws ArgumentError Formula("f", 1, value, check,          fail)
    @test_throws ArgumentError Formula("f", 1, value,        scale,   fail)

    # Initialization with name, string, and two functions
    @test Formula == typeof(   Formula("f", 1, value, check               ))
    @test Formula == typeof(   Formula("f", 1, value,        scale        ))
    @test Formula == typeof(   Formula("f", 1, value,               mirror))
    @test_throws ArgumentError Formula("f", 1, fails, check               )
    @test_throws ArgumentError Formula("f", 1, fails,        scale        )
    @test_throws ArgumentError Formula("f", 1, fails,               mirror)
    @test_throws ArgumentError Formula("f", 1, value,  fail               )

    # Initialization with name and string
    @test Formula == typeof(   Formula("f", 1, value                      ))
    @test_throws ArgumentError Formula("f", 1,  fail                      )

    # Initialization without name, with string, and four functions
    @test Formula == typeof(   Formula(1, value, check, scale, mirror))
    @test_throws ArgumentError Formula(1, fails, check, scale, mirror)
    @test_throws ArgumentError Formula(1, value,  fail, scale, mirror)
    @test_throws ArgumentError Formula(1, value, check,  fail, mirror)
    @test_throws ArgumentError Formula(1, value, check, scale,   fail)

    # Initialization without name, with string, and three functions
    @test Formula == typeof(   Formula(1, value, check, scale        ))
    @test Formula == typeof(   Formula(1, value, check,        mirror))
    @test Formula == typeof(   Formula(1, value,        scale, mirror))
    @test_throws ArgumentError Formula(1, fails, check, scale        )
    @test_throws ArgumentError Formula(1, fails, check,        mirror)
    @test_throws ArgumentError Formula(1, fails,        scale, mirror)
    @test_throws ArgumentError Formula(1, value, check,  fail        )
    @test_throws ArgumentError Formula(1, value,         fail, mirror)
    @test_throws ArgumentError Formula(1, value, check,          fail)
    @test_throws ArgumentError Formula(1, value,        scale,   fail)

    # Initialization without name, with string, and two functions
    @test Formula == typeof(   Formula(1, value, check               ))
    @test Formula == typeof(   Formula(1, value,        scale        ))
    @test Formula == typeof(   Formula(1, value,               mirror))
    @test_throws ArgumentError Formula(1, fails, check               )
    @test_throws ArgumentError Formula(1, fails,        scale        )
    @test_throws ArgumentError Formula(1, fails,               mirror)
    @test_throws ArgumentError Formula(1, value,  fail               )

    # Initialization without name and string
    @test Formula == typeof(   Formula(1, value                      ))
    @test_throws ArgumentError Formula(1,  fail                      )

    value  = "x -> x"
    fails  = "(x, y) -> x"

    # Initialization with string
    @test Formula == typeof(   Formula(value))
    @test_throws ArgumentError Formula(fails)

end


@testset "Piece" begin

    domain = (0, 1)
    included = (true, true)
    rule = [POLY]
    parameters = [[1]]

    # Initialization with a list of rules
    @test Piece == typeof(     Piece(domain, included, rule, parameters))
    @test_throws ArgumentError Piece((1, 0), included, rule, parameters)
    @test_throws ArgumentError Piece(domain, included, Vector{Formula}(undef, 0), parameters)
    @test_throws ArgumentError Piece(domain, included, rule, [[1], [1]])
    @test_throws ArgumentError Piece(domain, included, [LOG], [[1]])
    @test_throws ArgumentError Piece(domain, included, rule, [[]])
    @test_throws ArgumentError Piece(domain, included, rule, [ones(14)])

    # Initialization with a list of rules and without included
    @test Piece == typeof(Piece(domain, rule, parameters))

    rule = POLY
    parameters = [1]

    # Initialization with a single rule
    @test Piece == typeof(Piece(domain, included, rule, parameters))
    @test Piece == typeof(Piece(domain,           rule, parameters))

    # Initialization with a function
    @test Piece == typeof(     Piece(domain, included, x->x))
    @test_throws ArgumentError Piece(domain, included, (x, a)->x)
    @test Piece == typeof(     Piece(domain,           x->x))

	# Initialization with a string
    @test Piece == typeof(     Piece(domain, included, "x->x"))
    @test_throws ArgumentError Piece(domain, included, "(x, a)->x")
    @test Piece == typeof(     Piece(domain,           "x->x"))

end


@testset "PiecewiseFunction" begin

    parity = :even
    pieces = [Piece((0, 1), POLY, [1])]

	# Initialization with parity
    @test PiecewiseFunction == typeof(PiecewiseFunction(parity, pieces))
    @test_throws ArgumentError        PiecewiseFunction(:any,   pieces)
    @test_throws ArgumentError        PiecewiseFunction(parity, Vector{Piece}(undef, 0))
    @test_throws ArgumentError        PiecewiseFunction(parity, [pieces[1], pieces[1]])
    @test_throws ArgumentError        PiecewiseFunction(parity, [Piece((-1, 1), POLY, [1])])

	# Initialization with parity with a single piece
    @test PiecewiseFunction == typeof(PiecewiseFunction(parity, pieces[1]))

	# Initialization without parity
    @test PiecewiseFunction == typeof(PiecewiseFunction(pieces))
    @test PiecewiseFunction == typeof(PiecewiseFunction(pieces[1]))

end

nothing
