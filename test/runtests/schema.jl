import JSONSchema
@testset "independent JSON schema" begin
    schema = JSONSchema.Schema(JSON.parsefile(joinpath(@__DIR__, "..", "..", "schema.json")))
    for path in readdir(joinpath(@__DIR__, "..", "reference_models"); join = true)
        endswith(path, ".json") || continue
        doc = JSON.parsefile(path)
        @test JSONSchema.validate(schema, doc) === nothing
        bad = deepcopy(doc);bad["functions"][1]["x"] = "m12"
        @test JSONSchema.validate(schema, bad) !== nothing
    end
    model, ms, ctx = fixture()
    doc = document(model; system = ms, context = ctx)
    @test JSONSchema.validate(schema, doc) === nothing
end
