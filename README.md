# CascadeDecaysIO.jl

`CascadeDecaysIO.jl` writes `CascadeDecays.jl` cascade models to
amplitude-serialization-style dictionaries and JSON files. The current
implementation is writer-only: it serializes kinematics, topologies,
propagators, vertices, referenced functions, and the standard top-level JSON
sections used by amplitude-serialization documents.

> **Status:** this is a writer-first prototype. There is no JSON reader yet,
> and full write-read round-trip support is not implemented.

## What It Does

- Serializes a `CascadeSystem` into a `kinematics` dictionary.
- Serializes a concrete `DecayChain` into `topology`, `propagators`, and
  `vertices` entries.
- Serializes weighted chain lists into a `decay_description` dictionary.
- Collects referenced lineshapes and form factors into an `appendix`
  dictionary.
- Builds full amplitude-serialization documents with `distributions`,
  `functions`, and optional `domains`, `misc`, and `parameter_points` sections.
- Adds validation bundles and variable information through small mutating setter
  functions. Validation follows the ThreeBodyDecaysIO pattern:
  `misc.amplitude_model_checksums` plus `parameter_points`.
- Writes formatted JSON files with `writeJson`.
- Uses nested array nodes such as `[[[1, 2], 3], 4]` for four-body cascade
  topology addresses.

The leaves of a topology node are final-state particle indices. Nested arrays
represent successive two-body clusterings, so `[[[1, 2], 3], 4]` means particles
`1` and `2` form a subsystem, that subsystem combines with `3`, and the result
combines with `4`.

## Development Setup

This package is currently intended for local development, not registered package
installation. With `CascadeDecays.jl` and `CascadeDecaysIO.jl` cloned next to
each other, activate and test the package from the `CascadeDecaysIO.jl` folder:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.test()
```

The project pins the currently needed unregistered dependency revisions in
`Project.toml` under `[sources]`.

## Minimal Four-Body Example

The following example is executable and mirrors the test setup in
`test/runtests.jl`.

```julia
using CascadeDecays
using CascadeDecaysIO
using HadronicLineshapes
using ThreeBodyDecays: RecouplingLS

system = CascadeSystem(
    SystemSpins(0, 0, 0, 0; two_h0 = 0),
    SystemMasses(1.0, 1.1, 1.2, 1.3; m0 = 5.0),
)

topology = DecayTopology((((1, 2), 3), 4))

chain = DecayChain(
    topology;
    propagators = (
        ((1, 2), 3) => Propagator(2, BreitWigner(2.4, 0.1)),
        (1, 2) => Propagator(0, ConstantLineshape(1.0 + 0.0im)),
    ),
    vertices = (
        (((1, 2), 3), 4) => Vertex(RecouplingLS((0, 2))),
        ((1, 2), 3) => Vertex(RecouplingLS((0, 2))),
        (1, 2) => Vertex(RecouplingLS((0, 0))),
    ),
)

decay_description, appendix = serializeToDict(
    system,
    ["test_chain" => (1.0 + 0.0im, chain)];
    particle_labels = ("D0", "pi+", "D-", "K+", "B+"),
)

document = amplitudeSerializationDict(
    system,
    ["test_chain" => (1.0 + 0.0im, chain)];
    particle_labels = ("D0", "pi+", "D-", "K+", "B+"),
    name = "minimal_four_body_model",
    variables = ["m12sq", "m123sq"],
    domains = ["phase_space"],
    misc = Dict("source" => "minimal example"),
)

writeJson("minimal_four_body_model.json", document)
```

`decay_description` contains the model structure. `appendix` contains named
function definitions referenced by propagators and form factors.

`amplitudeSerializationDict` wraps the decay description and appendix into a
complete JSON-ready document. `writeJson` writes either an already constructed
document or a `CascadeSystem` plus weighted chains.

## Output Shape

The core serializer returns a pair:

```julia
decay_description, appendix = serializeToDict(system, weighted_chains)
```

For a complete amplitude-serialization document, use:

```julia
document = amplitudeSerializationDict(system, weighted_chains)
```

The document contains:

```json
{
  "distributions": [
    {
      "type": "HadronicUnpolarizedIntensity",
      "name": "cascade_model",
      "decay_description": {
        "kinematics": {},
        "reference_topology": [[[1, 2], 3], 4],
        "chains": []
      }
    }
  ],
  "functions": []
}
```

Optional document fields and distribution metadata can be inserted during
construction or later:

```julia
setVariables!(document, ["m12sq", "m123sq"])
setDomains!(document, ["phase_space"])
setParameterPoints!(document, Dict("nominal" => Dict()))
setMisc!(document, Dict("analysis" => "B2DxDK"))
setValidation!(
    document,
    Dict(
        "misc" => Dict(
            "amplitude_model_checksums" => [
                Dict("distribution" => "cascade_model", "point" => "validation_point", "value" => 1.0),
            ],
        ),
        "parameter_points" => [
            Dict("name" => "validation_point", "parameters" => []),
        ],
    ),
)
appendFunction!(document, "constant_one", Dict("type" => "ConstantLineshape", "value" => 1.0))
setSection!(document, "custom_section", Dict("kept" => true))
```

The specific setters keep the usual amplitude-serialization sections easy to
discover. `setSection!` is intentionally generic so experimental or future
schema sections can be added without changing the package.

The `decay_description` dictionary follows this structure:

```json
{
  "kinematics": {
    "initial_state": {"name": "B+", "mass": 5.0, "index": 0, "spin": 0},
    "final_state": [
      {"name": "D0", "mass": 1.0, "index": 1, "spin": 0}
    ]
  },
  "reference_topology": [[[1, 2], 3], 4],
  "chains": [
    {
      "name": "test_chain",
      "weight": "1.0 + 0.0i",
      "topology": [[[1, 2], 3], 4],
      "propagators": [
        {"node": [[1, 2], 3], "spin": 1, "parametrization": "BW_1_2_3"}
      ],

      "vertices": [
        {"node": [[[1, 2], 3], 4], "type": "ls", "l": 0, "s": 1, "formfactor": ""}
      ]
    }
  ]
}
```

The `appendix` dictionary maps generated function names to serializable
function dictionaries, for example `BreitWigner`, `MultichannelBreitWigner`,
`BlattWeisskopf`, or `ConstantLineshape`.

The `spin`, `l`, and `s` values in the JSON use ordinary angular-momentum units.
Internally, `CascadeDecays.jl` and `ThreeBodyDecays.jl` often store doubled
integer values, so `Propagator(2, ...)` is emitted as `"spin": 1`.

## Supported Serializers

| Julia object | Emitted type or section |
| --- | --- |
| `CascadeSystem` | `kinematics` |
| `DecayChain` | `topology`, `propagators`, `vertices` |
| weighted chain list | `decay_description` |
| `amplitudeSerializationDict` | full JSON-ready document |
| `BreitWigner` | `BreitWigner` function |
| `MultichannelBreitWigner` | `MultichannelBreitWigner` function |
| `ConstantLineshape` | `ConstantLineshape` function |
| `BlattWeisskopf` | `BlattWeisskopf` form factor |
| `MomentumPower` | `MomentumPower` form factor |
| `RecouplingLS` | vertex with `"type": "ls"` |
| `NoRecoupling` | vertex with `"type": "helicity"` |
| `ParityRecoupling` | vertex with `"type": "parity"` |

Unsupported lineshapes raise an `ArgumentError` that names the missing type.

## Running Tests

From the package root:

```powershell
julia --project=. -e "using Pkg; Pkg.test()"
```

If Julia/libgit2 inherits an incompatible `SSL_CERT_FILE`, clear it before
running tests:

```powershell
$env:JULIA_SSL_CA_ROOTS_PATH=''
Remove-Item Env:SSL_CERT_FILE -ErrorAction SilentlyContinue
julia --project=. -e "using Pkg; Pkg.test()"
```

## Building Documentation

The documentation site uses Quarto for the workflow tutorial and Documenter.jl
for the HTML site. From the package root:

```powershell
julia --project=docs -e "using Pkg; Pkg.instantiate()"
julia --project=docs docs/make.jl
```

The GitHub Actions workflow in `.github/workflows/Docs.yml` installs Quarto and
builds the site automatically for pushes, pull requests, and manual dispatches.

## Current Limitations

- Reader support is not implemented yet.
- The schema follows the `ThreeBodyDecaysIO.jl` field names, generalized to
  multi-step cascade topologies.
- Domains, variables, parameter points, metadata, and validation values are
  user-supplied. The writer does not infer physical phase-space boundaries or
  numerical reference amplitudes automatically.
- The API should be treated as experimental while the reader and round-trip
  tests are developed.

## Roadmap

- Add a reader for the emitted schema.
- Add writer-readback round-trip tests.
- Add a Documenter.jl site with usage, schema, and API reference pages.
- Expand supported lineshape and form-factor serializers as needed by analyses.
