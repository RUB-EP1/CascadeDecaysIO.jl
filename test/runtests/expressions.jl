@testset "expression contract" begin
    f = ExpressionFunction("1/(m^2-s-complex(0,1)*m*w)", ["s"]; parameters = Dict("m" => 1.2, "w" => 0.15))
    @test f(0.8) ≈ 1 / (1.2^2 - 0.8 - 0.18im)
    @test ExpressionFunction("sqrt(s)", ["s"])(-1.0) ≈ 1im
    @test ExpressionFunction("2", ["a", "b", "c"])(1, 2, 3) == 2
    @test_throws ArgumentError f(1, 2)
    for expr in ["run(`ls`)", "Base.exit(0)", "(x=1; x)", "unknown + s", "sqrt(s,s)", "[s]"]
        @test_throws ArgumentError ExpressionFunction(expr, ["s"])
    end
    @test_throws ArgumentError ExpressionFunction("s", ["s"]; parameters = Dict("s" => 1))
end
