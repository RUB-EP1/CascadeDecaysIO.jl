@testset "references and edits" begin
    model, masses, context = fixture()
    doc = document(model; system = masses, context)
    bw = only(f for f in doc["functions"] if f["name"] == "shared_bw")
    bw["mass"] = "resonance_mass"
    doc["parameter_points"] = [Dict("name" => "nominal", "parameters" => [Dict("name" => "resonance_mass", "value" => 1.2)])]
    loaded = build_model(doc)
    updated = set_parameters(loaded, Dict("resonance_mass" => 1.3))
    @test all(ch -> ch.propagators[1].m == 1.3, updated.model.chains)
    result = document(updated)
    @test only(f for f in result["functions"] if f["name"] == "shared_bw")["mass"] == "resonance_mass"
    @test only(only(result["parameter_points"])["parameters"])["value"] == 1.3
    @test_throws ArgumentError set_parameters(loaded, Dict("unknown" => 1))
    bad = deepcopy(doc);bad["functions"][1]["x"] = "m12"
    @test_throws ArgumentError validate_document(bad)
    bad = deepcopy(doc);push!(bad["functions"], deepcopy(bad["functions"][1]))
    @test_throws ArgumentError validate_document(bad)
    bad = deepcopy(doc);bad["format_version"] = "99"
    @test_throws ArgumentError validate_document(bad)
    bad = deepcopy(doc);bad["parameter_points"] = []
    @test_throws ArgumentError validate_document(bad)
    bad = deepcopy(doc);only(bad["distributions"])["decay_description"]["chains"][1]["propagators"][1]["parametrization"] = "absent"
    @test_throws ArgumentError validate_document(bad)
    reversed = CascadeDecay(reverse(loaded.model.chains), loaded.model.reference_topology; couplings = reverse(loaded.model.couplings), names = reverse(loaded.model.names))
    @test_throws ArgumentError document(reversed; system = masses, context = loaded.context)
    remapped = remap_context(loaded.context, Dict(1 => 2, 2 => 1))
    @test validate_document(document(reversed; system = masses, context = remapped)) isa AbstractDict
end
