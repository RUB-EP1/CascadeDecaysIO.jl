@testset "reference-frame transport" begin
    doc = read_document(joinpath(@__DIR__, "..", "reference_models", "lc2ppik-lhcb-2683025.json"))
    dd = only(doc["distributions"])["decay_description"]
    dd["chains"] = dd["chains"][[1, 13, 19]] # all three channels, half-integer parent
    doc["misc"] = Dict{String, Any}()
    for reference in 1:3
        i, j = TBD.ij_from_k(reference)
        dd["reference_topology"] = [[i, j], reference]
        c = build_model(doc);t = build_model(doc, ThreeBodyBackend())
        @test [ch["weight"] for ch in only(document(t)["distributions"])["decay_description"]["chains"]] ==
            [CDIO._complex_string(CDIO._parse_complex(ch["weight"])) for ch in dd["chains"]]
        for seed in 21:23
            event = phaseSpaceEvent(Xoshiro(seed), c.system.masses, c.model.reference_topology)
            inv = TBD.Invariants(t.system.ms; σ1 = CDIO._mass_squared(event[2] + event[3]), σ2 = CDIO._mass_squared(event[3] + event[1]))
            @test TBD.unpolarized_intensity(t.model, inv) ≈ CascadeDecays.unpolarized_intensity(c.model, KinematicPoint(kinematicTask(c.model), event)) rtol = 1.0e-9
        end
    end
end
