using CascadeDecaysIO

# Run from any working directory: julia --project examples/read_write.jl
input = joinpath(@__DIR__, "..", "test", "reference_models", "lb2pkg-lhcb-2765817.json")
doc = read_document(input)
cascade = build_model(doc, CascadeBackend())
threebody = build_model(doc, ThreeBodyBackend())
@assert all(result.passed for result in validate(doc))
mktempdir() do dir
    path = joinpath(dir, "model.json")
    write_document(path, cascade)
    reread = build_model(read_document(path))
    @assert document(reread) == document(cascade)
    println("Read both backends, checked physics, and verified a canonical round trip.")
end
