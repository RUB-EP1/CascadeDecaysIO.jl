@testset "divergence, arity and validation errors" begin
    m, ms, ctx = fixture()
    doc = document(m; system = ms, context = ctx)
    bad = deepcopy(doc)
    ff = first(f for f in bad["functions"] if f["type"] == "expression_function")
    ff["arguments"] = ["p", "a"];ff["parameters"]["b"] = 0.1
    @test_throws ArgumentError validate_document(bad)
    bad = deepcopy(doc);only(bad["domains"])["m0"] = 99
    @test_throws ArgumentError validate_document(bad)
    bad = deepcopy(doc);only(f for f in bad["functions"] if f["type"] == "BreitWigner")["width"] = -1
    @test_throws ArgumentError validate_document(bad)
    # One site changes while another retains a shared name: never silently choose one.
    changed = deepcopy(doc)
    fn = deepcopy(only(f for f in changed["functions"] if f["name"] == "shared_bw"))
    fn["name"] = "other_bw";fn["mass"] = 1.4;push!(changed["functions"], fn)
    only(changed["distributions"])["decay_description"]["chains"][2]["propagators"][1]["parametrization"] = "other_bw"
    altered = build_model(changed)
    @test_throws ArgumentError document(altered.model; system = ms, context = build_model(doc).context)
    # Independent equal values retain different generated definitions.
    independent = document(m; system = ms)
    @test count(f -> f["type"] == "BreitWigner", independent["functions"]) == 2
end
