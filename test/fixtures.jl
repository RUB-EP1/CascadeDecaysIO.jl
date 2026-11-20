function fixture()
    masses = SystemMasses(0.2, 0.3, 0.4; m0 = 2.5)
    spins = SystemSpins(0, 0, 0; two_h0 = 0)
    bw = BreitWigner(1.2, 0.15, 0.3, 0.4, 0, 1.5)
    ff = ExpressionFunction("1 + (p-a-b)/scale", ["p", "a", "b"]; parameters = Dict("scale" => 10.0))
    chains = map((1, 2)) do k
        i, j = TBD.ij_from_k(k)
        top = DecayTopology(((i, j), k))
        CascadeDecays.DecayChain(
            top, spins;
            propagators = ((i, j) => Propagator(0, bw),),
            vertices = ((i, j) => Vertex(TBD.RecouplingLS((0, 0)), ff), ((i, j), k) => Vertex(TBD.RecouplingLS((0, 0))))
        )
    end
    model = CascadeDecay(Tuple(chains), DecayTopology(((2, 3), 1)); couplings = (1.0, 0.3 + 0.2im), names = ("same", "same"))
    context = ModelIOContext()
    context.names[Site(1, (2, 3), :propagator)] = "shared_bw"
    context.names[Site(2, (3, 1), :propagator)] = "shared_bw"
    return model, masses, context
end
