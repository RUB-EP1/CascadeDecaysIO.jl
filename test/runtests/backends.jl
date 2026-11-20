@testset "both backends and canonical documents" begin
    model, masses, context = fixture()
    doc = document(model; system = masses, context)
    @test length(doc["functions"]) == 3
    @test all(f -> !haskey(f, "x"), doc["functions"])
    loaded = build_model(doc)
    @test loaded.model.chains[1].propagators[1] isa BreitWigner
    @test loaded.model.chains[1].vertices[1].ff isa Union{ExpressionFunction, TBD.NoFormFactor}
    @test document(loaded) == doc
    task = kinematicTask(model)
    for seed in 1:5
        event = phaseSpaceEvent(Xoshiro(seed), masses, model.reference_topology)
        point = KinematicPoint(task, event)
        @test CascadeDecays.unpolarized_intensity(loaded.model, point) ≈ CascadeDecays.unpolarized_intensity(model, point)
    end
    three = build_model(doc, ThreeBodyBackend())
    @test three.model.chains[1].Xlineshape isa BreitWigner
    td = document(three)
    @test document(build_model(td, ThreeBodyBackend())) == td
    # Scalar helicities permit direct comparison without frame-convention differences.
    for seed in 6:10
        event = phaseSpaceEvent(Xoshiro(seed), masses, model.reference_topology)
        pairmass(a, b) = begin
            p = event[a] + event[b]; p.E^2 - p.px^2 - p.py^2 - p.pz^2
        end
        inv = TBD.Invariants(three.system.ms; σ1 = pairmass(2, 3), σ2 = pairmass(3, 1))
        @test TBD.unpolarized_intensity(three.model, inv) ≈ CascadeDecays.unpolarized_intensity(model, KinematicPoint(task, event)) rtol = 1.0e-10
    end
    mktempdir() do dir
        path = joinpath(dir, "model.json")
        write_document(path, loaded)
        first = read(path, String)
        write_document(path, build_model(read_document(path)))
        @test read(path, String) == first
    end
    @test_throws ArgumentError document(model; system = masses, context = ModelIOContext(names = Dict(Site(9, (1, 2), :propagator) => "bad")))
end
