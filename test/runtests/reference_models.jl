@testset "reference physics and both readers" begin
    for file in sort(readdir(joinpath(@__DIR__, "..", "reference_models"); join = true))
        endswith(file, ".json") || continue
        @testset "$(basename(file))" begin
            doc = read_document(file)
            checks = validate(doc)
            @test !isempty(checks)
            @test all(c -> c.passed, checks)
            for distribution in 1:length(doc["distributions"])
                cascade = build_model(doc; distribution)
                three = build_model(doc, ThreeBodyBackend(); distribution)
                @test document(build_model(document(cascade))) == document(cascade)
                @test document(build_model(document(three), ThreeBodyBackend())) == document(three)
                for point in doc["parameter_points"]
                    event = CDIO._four_vectors(point)
                    pairmass(a, b) = CDIO._mass_squared(event[a] + event[b])
                    inv = TBD.Invariants(three.system.ms; σ1 = pairmass(2, 3), σ2 = pairmass(3, 1))
                    ci = CascadeDecays.unpolarized_intensity(cascade.model, KinematicPoint(kinematicTask(cascade.model), event))
                    @test TBD.unpolarized_intensity(three.model, inv) ≈ ci rtol = 1.0e-8
                end
            end
        end
    end
end
