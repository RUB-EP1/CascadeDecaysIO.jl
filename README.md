# CascadeDecaysIO.jl

`CascadeDecaysIO.jl` is the serialization and deserialization package for [`CascadeDecays.jl`](https://github.com/RUB-EP1/CascadeDecays.jl). It converts multi-body cascade amplitude models between live Julia objects (`CascadeDecay`, `DecayChain`, `SystemSpins`) and the framework-agnostic [**Amplitude Serialization Format**](https://github.com/RUB-EP1/amplitude-serialization) dictionaries and JSON files.

---

## Core Functionality

- **Writer API (`serializeToDict`, `amplitudeSerializationDict`, `writeJson`)**:
  Exports any `CascadeDecay` model (or a spin system plus weighted chains) into an amplitude-serialization dictionary or formatted JSON file.
- **Reader API (`readJson`, `dict2instance`)**:
  Reconstructs a live `CascadeDecay` model, `CascadeSystem`, and `SystemSpins` directly from a JSON file or parsed dictionary so amplitude evaluations can be restarted from the JSON alone.
- **Clean Function Naming (`NamedLineshape`, `propagator_names`)**:
  Lets you assign human-readable names to lineshapes in `"parametrization"` and `"functions"` without structural hashes.
- **Document Modifiers (`setVariables!`, `setDomains!`, `setValidation!`, `setParameterPoints!`, `setMisc!`, `appendFunction!`, `setSection!`)**:
  Mutating helpers for attaching metadata, phase-space domains, and verification points (`four_vectors` and `amplitude_model_checksums`).
- **Custom Lineshape Support (`custom_functions`)**:
  Serializes custom lineshapes under `"type": "custom"` and allows binding arbitrary Julia callables by name during deserialization.

---

## Installation & Setup

With `CascadeDecays.jl` and `CascadeDecaysIO.jl` cloned in the same workspace:

```julia
using Pkg
Pkg.activate("CascadeDecaysIO.jl")
Pkg.instantiate()
Pkg.test()
```

---

## Tutorial: How to Set Up the JSON Writer and Reader

### Part 1: Setting Up the JSON Writer

#### 1. Naming Lineshapes (`NamedLineshape` or `propagator_names`)

When `CascadeDecaysIO` serializes a `DecayChain`, every propagator's lineshape and vertex form factor is placed in the top-level `"functions"` list and referenced by name in the chain.

You can control propagator function names in two ways:
- **Option A — Wrap lineshapes in `NamedLineshape` when building your chains**:
  ```julia
  using CascadeDecays, CascadeDecaysIO, HadronicLineshapes

  my_bw = NamedLineshape("lineshape_resonance_1", BreitWigner(m0, width, ma, mb, l, d))
  ```
- **Option B — Pass `propagator_names` during serialization** (without modifying existing chains):
  ```julia
  # Pass a vector of names per propagator in the chain, or a Dict mapping topology nodes to names:
  propagator_names = Dict((1, 2) => "lineshape_sub", ((1, 2), 3) => "lineshape_main")
  ```
*(Vertex form factors such as `BlattWeisskopf{L}(R)` and `MomentumPower{L}()` are named automatically as `BlattWeisskopf_L<L>_d<R>` and `MomentumPower_l<L>`.)*

#### 2. Building the Serialization Dictionary (`amplitudeSerializationDict` & `serializeToDict`)

If you have a `CascadeDecay` object `model` (or a `SystemSpins` / `CascadeSystem` and a list of `"chain_name" => (weight, chain)` pairs):

- **Full JSON document (`amplitudeSerializationDict`)**:
  Builds the complete top-level dictionary (`"distributions"`, `"functions"`, and any optional sections):
  ```julia
  document = amplitudeSerializationDict(
      model;
      name               = "my_distribution_name",
      particle_labels    = ("p1", "p2", "p3", "p4", "P0"), # finals..., initial
      variables          = ["m1_2sq", "m1_2_3sq"],
      propagator_names   = nothing,                        # optional override
      reference_topology = nothing,                        # defaults to model.reference_topology
  )
  ```
  *(You can also call `amplitudeSerializationDict(system_or_spins, weighted_chains; ...)`.)*

- **Low-level decay description & appendix (`serializeToDict`)**:
  If you only need the inner `decay_description` dictionary and the `appendix` of collected functions:
  ```julia
  decay_description, appendix = serializeToDict(
      model;
      particle_labels = ("p1", "p2", "p3", "p4", "P0"),
  )
  ```

#### 3. Adding Domains, Metadata, and Verification Points (`set...` Helpers)

You can pass optional sections directly as keyword arguments to `amplitudeSerializationDict`, or attach/update them on `document` using mutating setter functions:

```julia
# Set kinematic variable names on distributions[1]
setVariables!(document, ["m1_2sq", "m1_2_3sq"])

# Define phase-space domains
setDomains!(document, [
    Dict(
        "name" => "phase_space",
        "type" => "product_domain",
        "axes" => [
            Dict("name" => "m1_2sq",   "min" => 0.5, "max" => 9.0),
            Dict("name" => "m1_2_3sq", "min" => 2.0, "max" => 20.0),
        ],
    ),
])

# Attach free-form metadata
setMisc!(document, Dict("description" => "Example cascade model", "generator" => "CascadeDecaysIO.jl"))

# Attach verification parameter points (e.g. sampled four-vectors) and reference checksums
setValidation!(document, Dict(
    "parameter_points" => [
        Dict(
            "name" => "point_1",
            "four_vectors" => Dict(
                "p1" => Dict("E" => E1, "px" => px1, "py" => py1, "pz" => pz1),
                "p2" => Dict("E" => E2, "px" => px2, "py" => py2, "pz" => pz2),
                "p3" => Dict("E" => E3, "px" => px3, "py" => py3, "pz" => pz3),
                "p4" => Dict("E" => E4, "px" => px4, "py" => py4, "pz" => pz4),
            ),
        ),
    ],
    "misc" => Dict(
        "amplitude_model_checksums" => [
            Dict(
                "distribution"   => "my_distribution_name",
                "point"          => "point_1",
                "value"          => abs2(amp), # unpolarized intensity |A|^2
                "amplitude_real" => real(amp),
                "amplitude_imag" => imag(amp),
            ),
        ],
    ),
))

# Append an extra function definition or custom top-level section if needed
appendFunction!(document, "extra_fn", Dict("type" => "ConstantLineshape", "value" => "1.0 + 0.0i"))
setSection!(document, "custom_section", Dict("version" => 1))
```

#### 4. Writing to Disk (`writeJson`)

Write an already built `document` dictionary, or serialize and write a `CascadeDecay` model in a single step:

```julia
# From an existing document dictionary:
writeJson("model.json", document; indent = 4)

# Or directly from a CascadeDecay model (accepts all amplitudeSerializationDict keyword arguments):
writeJson(
    "model.json",
    model;
    name            = "my_distribution_name",
    particle_labels = ("p1", "p2", "p3", "p4", "P0"),
)
```

---

### Part 2: Setting Up the JSON Reader

#### 1. Reading a JSON File (`readJson`)

`readJson` parses a JSON file and reconstructs the runnable `CascadeDecay` model, the `CascadeSystem`, and the raw parsed dictionary:

```julia
model, system, document = readJson("model.json")
```

- `model::CascadeDecay`: Ready to evaluate on any `KinematicPoint` via `amplitude(model, point)` or `unpolarized_intensity(model, point)`.
- `system::CascadeSystem`: Contains the external `SystemSpins` (`system.quantum`) and `SystemMasses` (`system.masses`; defaults to `0.0` when masses are omitted from `kinematics`).
- `document::LittleDict{String,Any}`: The parsed JSON dictionary (allowing access to `document["parameter_points"]`, `document["misc"]`, etc.).

#### 2. Deserializing from an Existing Dictionary (`dict2instance`)

If you already have a parsed JSON dictionary in memory, use `dict2instance` to reconstruct specific types:

```julia
model  = dict2instance(CascadeDecay, document)
system = dict2instance(CascadeSystem, document)
spins  = dict2instance(SystemSpins, document)
```

#### 3. Providing Custom Lineshape Functions (`custom_functions`)

When a JSON file contains `"type": "custom"` functions that are not built-in types, pass a `custom_functions` dictionary mapping the JSON function `"name"` to any callable Julia object `f(σ)`:

```julia
my_custom_lineshapes = Dict{String,Any}(
    "my_custom_propagator" => (σ -> 1.0 / (3.5^2 - σ - 0.2im)),
)

model, system, document = readJson("model.json"; custom_functions = my_custom_lineshapes)
# or:
model = dict2instance(CascadeDecay, document; custom_functions = my_custom_lineshapes)
```

---

## Minimal Copy-Paste Template (Writer & Reader)

For a fuller multi-topology example with four-vector verification points, run [`examples/dummy_analysis.jl`](examples/dummy_analysis.jl):

```powershell
julia --project=CascadeDecaysIO.jl CascadeDecaysIO.jl/examples/dummy_analysis.jl
```

Below is a minimal self-contained template showing both writing and reading:

```julia
using CascadeDecays
using CascadeDecays.ThreeBodyDecays: RecouplingLS
using CascadeDecaysIO
using HadronicLineshapes

# 1. Build a minimal CascadeDecay model
spins    = SystemSpins(0, 0, 0, 0; two_h0 = 0)
topology = DecayTopology((((1, 2), 3), 4))

chain = DecayChain(
    topology,
    spins;
    propagators = (
        (1, 2)      => Propagator(2, NamedLineshape("const_12", ConstantLineshape(1.0 + 0.0im))),
        ((1, 2), 3) => Propagator(2, NamedLineshape("bw_123",   BreitWigner(3.5, 0.1, 1.8, 1.5, 1, 3.0))),
    ),
    vertices = (
        (((1, 2), 3), 4) => Vertex(RecouplingLS((2, 2)), BlattWeisskopf{1}(3.0)),
        ((1, 2), 3)      => Vertex(RecouplingLS((0, 2)), BlattWeisskopf{0}(3.0)),
        (1, 2)           => Vertex(RecouplingLS((2, 0))),
    ),
)

model = CascadeDecay((chain,), topology; couplings = (1.0 + 0.0im,), names = ("chain_1",))

# 2. Write to JSON
writeJson(
    "example_model.json",
    model;
    name            = "example_intensity",
    particle_labels = ("p1", "p2", "p3", "p4", "P0"),
    variables       = ["m1_2sq", "m1_2_3sq"],
)

# 3. Read back from JSON
reloaded_model, reloaded_system, reloaded_doc = readJson("example_model.json")
```

---

## Supported Types & JSON Mapping

| Julia Object / Function | Emitted JSON Representation | Reader / Round-Trip Support |
| --- | --- | --- |
| `CascadeDecay` | Top-level `distributions` + `functions` | `readJson(path)` / `dict2instance(CascadeDecay, dict)` |
| `CascadeSystem` / `SystemSpins` | `kinematics` (`initial_state`, `final_state`) | `dict2instance(CascadeSystem, dict)` / `dict2instance(SystemSpins, dict)` |
| `DecayChain` | `topology`, `propagators`, `vertices`, `weight` | Reconstructed inside `CascadeDecay` |
| `NamedLineshape(name, ls)` | Exact `"parametrization": name` & `"functions"` entry | Reconstructed into underlying lineshape |
| `BreitWigner` | `"type": "BreitWigner"` | `HadronicLineshapes.BreitWigner` |
| `MultichannelBreitWigner` | `"type": "custom"`, `"subtype": "MultichannelBreitWigner"` | `HadronicLineshapes.MultichannelBreitWigner` |
| `ConstantLineshape` | `"type": "ConstantLineshape"` | `CascadeDecays.ConstantLineshape` |
| `TFPWAMultichannelBreitWigner` | `"type": "custom"`, `"subtype": "TFPWAMultichannelBreitWigner"` | `CascadeDecaysIO.TFPWAMultichannelBreitWigner` |
| `NRExpLineshape` | `"type": "custom"`, `"subtype": "NRExpLineshape"` | `CascadeDecaysIO.NRExpLineshape` |
| `BlattWeisskopf{L}(R)` | `"type": "BlattWeisskopf"` (`BlattWeisskopf_L<L>_d<R>`) | `HadronicLineshapes.BlattWeisskopf{L}(R)` |
| `MomentumPower{L}()` | `"type": "MomentumPower"` | `HadronicLineshapes.MomentumPower{L}()` |
| `RecouplingLS((2L, 2S))` | Vertex with `"type": "ls"`, `"l": L`, `"s": S` | `ThreeBodyDecays.RecouplingLS` |
| `NoRecoupling(2λa, 2λb)` | Vertex with `"type": "helicity"` | `ThreeBodyDecays.NoRecoupling` |
| `ParityRecoupling(2λa, 2λb, ±)` | Vertex with `"type": "parity"` | `ThreeBodyDecays.ParityRecoupling` |
