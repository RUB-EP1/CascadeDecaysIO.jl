# Read, edit and write a model

Both backends read the same validated document. Exporting a loaded model writes
its selected distribution. Writing a document dictionary preserves all distributions.

```@example workflow
using CascadeDecaysIO
path = joinpath(pkgdir(CascadeDecaysIO), "test", "reference_models",
    "lb2pkg-lhcb-2765817.json")
doc = read_document(path)
loaded = build_model(doc, CascadeBackend())
threebody = build_model(doc, ThreeBodyBackend())
@assert all(result.passed for result in validate(doc))
nothing # hide
```

Canonical writing validates before opening the output file. Old checksums are
omitted; use `document(loaded; validation=true)` for fresh Cascade checks.

```@example workflow
mktempdir() do dir
    output = joinpath(dir, "model.json")
    write_document(output, loaded)
    reread = build_model(read_document(output))
    @assert document(reread) == document(loaded)
end
```

Generic expressions declare local arguments independently of parameter values:

```@example workflow
shape = ExpressionFunction("1/(mass^2-sigma-complex(0,1)*mass*width)",
    ["sigma"]; parameters=Dict("mass" => 1.52, "width" => 0.015))
shape(2.0)
```

For hand-built models use `document(model; system)` and optionally a
`ModelIOContext`. Assign names using `context.names[Site(chain, node, role)]`,
where role is `:propagator` or `:formfactor`. Shared names must have identical
function definitions. For a chain permutation explicitly use `remap_context`;
for topology changes rebuild the context.

Use `set_parameters(loaded, Dict("parameter_name" => value))` to rebuild all
functions bound to a named parameter. This requires a named parameter point in
the input document; literal numeric fields are independent values.
