# CascadeDecaysIO.jl

`CascadeDecaysIO.jl` writes [`CascadeDecays.jl`](https://github.com/RUB-EP1/CascadeDecays.jl) models to the framework-agnostic [**amplitude-serialization**](https://github.com/RUB-EP1/amplitude-serialization) JSON format and reads them back. The output has the same structure as the reference models of the standard; events for validation are given as final-state four-vectors ([amplitude-serialization#111](https://github.com/RUB-EP1/amplitude-serialization/pull/111)).

---

## Core Functionality

- **Writer** (`amplitudeSerializationDict`, `serializeToDict`, `writeJson`): converts a `CascadeDecay` plus the external masses into a complete document with the five required sections `distributions`, `functions`, `domains`, `misc`, `parameter_points`.
- **Reader** (`readJson`, `dict2instance`): rebuilds a live `CascadeDecay` and `CascadeSystem` from a file or dictionary.
- **Validation** (`validationSection`, `validationEvents`, `phaseSpaceEvent`, `validate`, `kinematicTask`): generates on-shell events at which the model has a substantial intensity, stores them as four-vector points with intensity checksums, and recomputes them from a file.
- **No lineshapes of its own**: lineshapes come from `CascadeDecays.jl` and `HadronicLineshapes.jl`. Anything else is serialized as a `generic_function` (an analytic expression) — either with [`GenericFunction`](#analysis-specific-lineshapes) or by the analysis package that defines the lineshape.
- **Strict reading**: unknown function types, missing form-factor references, missing masses, malformed spins and unsupported distribution types are errors, never silent defaults.
- **Round trip**: functions keep their names on reading, so a model read from a file can be written again unchanged. The reference models of the standard (`test/reference_models/`) are read, validated and re-written in the tests.

---

## Installation & Setup

With `CascadeDecays.jl` and `CascadeDecaysIO.jl` cloned side by side:

```julia
using Pkg
Pkg.activate("CascadeDecaysIO.jl")
Pkg.instantiate()
Pkg.test()
```

---

## Quick Start

```julia
using CascadeDecays, CascadeDecaysIO, HadronicLineshapes
using CascadeDecays.ThreeBodyDecays: RecouplingLS

spins = SystemSpins(0, 0, 0; two_h0 = 0)
masses = SystemMasses(0.14, 0.14, 0.14; m0 = 1.5)      # final states..., m0 = initial state
topology = DecayTopology(((1, 2), 3))

chain = DecayChain(topology, spins;
    propagators = ((1, 2) => Propagator(2, NamedLineshape("rho", BreitWigner(0.775, 0.15, 0.14, 0.14, 1, 1.5))),),
    vertices = (
        ((1, 2), 3) => Vertex(RecouplingLS((2, 2)), BlattWeisskopf{1}(5.0)),
        (1, 2) => Vertex(RecouplingLS((2, 0)), BlattWeisskopf{1}(1.5)),
    ))
model = CascadeDecay((chain,), topology; couplings = (1.0 + 0.0im,), names = ("rho",))

# write; three validation points (four-vectors and intensity checksums) are generated
writeJson("model.json", model, masses;
    name = "model_intensity", particle_labels = ("pi1", "pi2", "pi3", "X"))

# read back and recompute the checksums
model2, system, document = readJson("model.json")
all(r -> r.passed, validate(document, model2))
```

`examples/dummy_analysis.jl` is a complete four-body example with two topologies, generated on-shell events and the full round trip.

---

## Writing

- `amplitudeSerializationDict(cascade, masses; name, particle_labels, validation, domains, propagator_names, reference_topology)` builds the document; `writeJson(path, document)` writes it. `writeJson(path, cascade, masses; kwargs...)` does both.
- **Masses are required**: they go only into the default `phase_space` domain (`phaseSpaceDomain(masses)`), each with its particle index.
- **Spins, helicities, `l` and `s`** are written as strings (`"0"`, `"1"`, `"1/2"`), weights as complex strings (`"0.5 - 0.1i"`).
- **Names**: wrap lineshapes or form factors in `NamedLineshape("name", f)`, or pass `propagator_names` (a vector per chain or a `Dict` from node to name). Default form-factor names are `BlattWeisskopf_l<L>` and `MomentumPower_l<L>`. Two different functions with the same name are an error.
- **Argument names**: `BreitWigner`, `MultichannelBreitWigner` and `Polynomial` carry `"x"`, as in the reference models; `generic_function` has none. The name follows the `m_ij_sq` pattern of the reference models, extended to nested nodes: `m_12_sq`, `m_12_3_sq`, `m_34_sq`. A propagator is always evaluated at the invariant mass squared of its node; `x` is a label (a function used on several nodes keeps the `x` of its first use).
- **ASCII only** by default (`writeJson(...; ascii_only = false)` to switch off), so that readers in any language need no Unicode handling.

## Analysis-Specific Lineshapes

A lineshape that is not a `CascadeDecays` / `HadronicLineshapes` type has no standard representation. The analysis package that defines it adds a serializer, usually to a `generic_function`:

```julia
CascadeDecaysIO.serializeToDict(f::MyLineshape) = (
    CascadeDecaysIO.LittleDict{String,Any}("type" => "generic_function", "expression" => my_expression(f)),
    Dict{String,Any}(),
)
```

or the model uses `GenericFunction` directly:

```julia
nr = NamedLineshape("nonresonant", GenericFunction("exp(-(0.15 + 0.25 * i) * (sigma - 14.44))"))
```

A `generic_function` expression has exactly one free variable (by convention `sigma`, the invariant mass squared of the propagator node), numbers, `+ - * / ^`, elementary functions and the imaginary unit `i`. `GenericFunction` checks the expression before compiling it, so a JSON file cannot run arbitrary code. Writing an object without a serializer raises an `ArgumentError` naming the type.

## Reading

```julia
model, system, document = readJson("model.json")
model, system, document = readJson("model.json"; distribution = "name")   # by name or index; default 1
model = dict2instance(CascadeDecay, document)          # from a parsed dictionary
```

Every function comes back wrapped in a `NamedLineshape` with its name from the file. Only `HadronicUnpolarizedIntensity` distributions are supported.

Standard types are rebuilt as `HadronicLineshapes` / `CascadeDecays` objects and `generic_function` as `GenericFunction`. Other function types must be bound by name:

```julia
model, system, document = readJson("model.json"; custom_functions = Dict("lineshape_X" => σ -> my_shape(σ)))
```

## Validation

```julia
validation = validationSection(model, masses; distribution = "model_intensity")   # generate 3 points
validation = validationSection(model, masses; distribution = "model_intensity",
    node_masses = Dict((1, 2) => m12))   # node (1, 2) at a fixed mass m12
validation = validationSection(model, events; distribution = "model_intensity")   # given events
results = validate(document)                 # every checksum of every distribution
results = validate(document, model)          # use `model` for the first distribution
results = validate("model.json")
```

`amplitudeSerializationDict(model, masses; ...)` and `writeJson(path, model, masses; ...)` generate three points with `validationSection(model, masses; distribution = name)` unless `validation` is given (`validation = nothing`: no points). The generator (`validationEvents`, `phaseSpaceEvent`) builds events by successive isotropic two-body decays along the reference topology; internal node masses are uniform between their kinematic limits or fixed (`node_masses = Dict(node => m)`, e.g. a narrow resonance in the narrow-width approximation). An event is accepted if the intensity there is at least `min_relative_intensity` (default `0.1`) times the mean intensity of `pilot` (default `1000`) events drawn the same way; a model that vanishes on all of them, typically because a lineshape is narrower than the sampled mass range, is an error that asks for the node mass. A given event with vanishing intensity is also an error.

Each event becomes a point `{"name", "domain", "final_state": [{"index", "p": {"E", "px", "py", "pz"}}, ...]}` with a checksum `{"point", "distribution", "value"}`, where `value` is the unpolarized intensity `Σ|A|²`. `validate` recomputes every checksum with the distribution it names (any other name is an error), matching particles by index, and checks that the four-vectors lie on their `phase_space` domain: `|p² − m²|` in GeV² and the total four-momentum `(m0, 0, 0, 0)` in GeV, both within `atol_phase_space` (default `1e-8`). `kinematicTask(model)` returns the `KinematicTask` for evaluating the model on four-vectors.

---

## Supported Types & JSON Mapping

| Julia object | JSON | Read back as |
| --- | --- | --- |
| `CascadeDecay` + `SystemMasses` | `distributions`, `functions`, `domains`, `misc`, `parameter_points` | `readJson` / `dict2instance(CascadeDecay, dict)` |
| `CascadeSystem` | `kinematics` (`index`, `name`, `spin`) and `phase_space` domain (`index`, `mass`) | `dict2instance(CascadeSystem, dict)` |
| `DecayChain` | `name`, `topology`, `propagators`, `vertices`, `weight` | part of `CascadeDecay` |
| `BreitWigner` | `"BreitWigner"` | `HadronicLineshapes.BreitWigner` |
| `MultichannelBreitWigner`, `Flatte` | `"MultichannelBreitWigner"` | `HadronicLineshapes.MultichannelBreitWigner` |
| `ConstantLineshape` (real) | `"Polynomial"`, `"coefficients": [c]` | `CascadeDecays.ConstantLineshape` |
| `GenericFunction` / analysis serializer | `"generic_function"` | `GenericFunction` |
| `BlattWeisskopf{L}(R)` | `"BlattWeisskopf"`, `"l"`, `"radius"` | `HadronicLineshapes.BlattWeisskopf{L}(R)` |
| `MomentumPower{L}()` | `"MomentumPower"` | `HadronicLineshapes.MomentumPower{L}()` |
| `NoFormFactor` | `"formfactor": ""` | `ThreeBodyDecays.NoFormFactor` |
| `RecouplingLS((2L, 2S))` | `"type": "ls"`, `"l": "L"`, `"s": "S"` | `ThreeBodyDecays.RecouplingLS` |
| `NoRecoupling(2λa, 2λb)` | `"type": "helicity"`, `"helicities"` | `ThreeBodyDecays.NoRecoupling` |
| `ParityRecoupling(2λa, 2λb, ±)` | `"type": "parity"`, `"helicities"`, `"parity_factor"` | `ThreeBodyDecays.ParityRecoupling` |
