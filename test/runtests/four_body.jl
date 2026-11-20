@testset "four-body adapter and regeneration" begin
    ms = SystemMasses(0.2, 0.3, 0.4, 0.5; m0 = 3.0)
    top = DecayTopology(((1, 2), (3, 4)))
    spins = SystemSpins(0, 0, 0, 0; two_h0 = 0)
    chain = CascadeDecays.DecayChain(
        top, spins;
        propagators = ((1, 2) => Propagator(0, BreitWigner(1.0, 0.1)), (3, 4) => Propagator(0, BreitWigner(1.2, 0.2))),
        vertices = (
            ((1, 2), (3, 4)) => Vertex(TBD.RecouplingLS((0, 0))),
            (1, 2) => Vertex(TBD.RecouplingLS((0, 0))), (3, 4) => Vertex(TBD.RecouplingLS((0, 0))),
        )
    )
    model = CascadeDecay((chain,), top; couplings = (1.0,), names = ("four_body",))
    doc = document(model; system = ms, validation = true)
    @test all(c -> c.passed, validate(doc))
    loaded = build_model(doc)
    @test all(c -> c.passed, validate(document(loaded; validation = true)))
    @test isempty(document(loaded)["misc"])
    @test_throws ArgumentError build_model(doc, ThreeBodyBackend())
    @test document(loaded) == document(build_model(document(loaded)))
end
