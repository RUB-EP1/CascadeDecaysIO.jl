# CascadeDecaysIO.jl

Unified prototype amplitude-model IO for CascadeDecays and ThreeBodyDecays.
One JSON contract, function codec layer, expression evaluator, parameter-binding
implementation and canonical writer serve both physics backends.

This is a direct prototype redesign. There are no legacy readers or named physics
wrappers. The package retains the CascadeDecaysIO name, UUID and repository history.
Version 0.2 introduces the new format with both physics backends.

## Read, build, evaluate, write

```julia
using CascadeDecaysIO

doc = read_document("model.json")
loaded = build_model(doc, CascadeBackend(); distribution=1)
three = build_model(doc, ThreeBodyBackend(); distribution=1)

# loaded.model contains plain physics callables, loaded.system physical masses/spins.
# loaded.context holds names and parameter-reference paths.
write_document("model-out.json", loaded)
updated = set_parameters(loaded, Dict("resonance_mass" => 1.53))
write_document("fit.json", updated)
```

For hand-built models use `document(model; system)` or
`write_document(path, model; system)`. Optionally supply `context=ModelIOContext()`
and assign `context.names[Site(chain_index, node_tuple, :propagator)] = "name"`.
Use `:formfactor` for vertex functions. Names are never attached to physics objects.
`serializeToDict` is the extension point for analysis-owned callable types; emit a
supported built-in definition or an `expression_function` definition.

`read_document` validates the whole document. `build_model` builds one selected
distribution (index or name). Export of a loaded model exports that selected model,
not the other distributions. `write_document(path, doc)` writes a validated full document.

## Function contract

Built-ins define their own calling convention: lineshapes receive one invariant
mass squared; form factors receive parent, first-child, second-child masses squared.
The topology supplies those values. Function definitions never contain `x`.

```julia
f = ExpressionFunction("1/(mass^2-sigma-complex(0,1)*mass*width)",
    ["sigma"]; parameters=Dict("mass"=>1.52, "width"=>0.015))
f(2.0)
```

Generic expressions declare ordered local `arguments` and separate `parameters`.
Parameters are numbers in the callable; JSON parameter fields can instead reference
named values in a selected parameter point. `set_parameters` rebuilds the numerical
model consistently from those bindings. It never mutates a global workspace.

Only whitelisted expression syntax is interpreted: no `eval`, property access,
assignments, imports, or arbitrary function calls. See [the contract](docs/format.md).

## Editing and sharing

Explicit site names survive parameter-only edits. Assigning one name to several
sites requests one shared definition; divergent definitions or shared parameter
values are errors. Equal independent objects are not automatically deduplicated.
Generated names use chain index, topology address, role and function type, never hashes.

For a chain permutation use `remap_context(context, Dict(old_index=>new_index, ...))`.
For insertions, deletions, or topology changes supply a fresh context or rebuild the
mapping explicitly. Identical-label, identical-topology swaps cannot be inferred:
remapping is a caller responsibility even when structural checks cannot detect it.

Exports omit old checksums and validation points. For Cascade models,
`document(loaded; validation=true)` generates deterministic fresh four-vector checks.
`validate(doc)` independently evaluates stored checksums. It returns per-point
results; it does not silently turn a failing check into success.

## Development

Julia 1.11+. The prototype pins the CascadeDecays source revision; both backends
are direct dependencies because CascadeDecays itself depends on ThreeBodyDecays.

```sh
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

For a separate analysis environment, add the unregistered backend before this package:

```julia
using Pkg
Pkg.add(url="https://github.com/RUB-EP1/CascadeDecays.jl.git",
        rev="9c9c814f9fc300bc12756c4e6b86ac9126e71ab8")
Pkg.develop(path="../CascadeDecaysIO.jl")
```

Tests cover expression scope, canonical round trips, shared parameters, structural
edits, three- and four-body models, and reference physics fixtures. The JSON schema
is mirrored from the coordinated amplitude-serialization change; keep both copies
identical until the package is published and a single distribution mechanism is chosen.

## Provenance

Bare-function codecs, cascade construction and phase-space validation utilities
were adapted from RUB-EP1/CascadeDecaysIO.jl commit
`69720159976b91f633aa003565ae37922ebbeedf` (source retained with attribution).
Reference four-vector fixtures originate from amplitude-serialization revision
`16fdcf25fa0a4ffea1f61f8ad4ddcbabef080642`, via that IO repository.
The new schema removes argument bindings and adds explicit local callable interfaces.
