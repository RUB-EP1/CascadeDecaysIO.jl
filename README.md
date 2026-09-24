# CascadeDecaysIO.jl

`CascadeDecaysIO.jl` is the input/output and serialization package for [`CascadeDecays.jl`](https://github.com/RUB-EP1/CascadeDecays.jl). It converts multi-body hadronic cascade decay models between live Julia `CascadeDecay` objects and the framework-agnostic **Amplitude Serialization Format** (`HSF-aje-2025-02`) dictionaries and JSON files.

Its goal is to allow amplitude analyses built in `CascadeDecays.jl` to be exported to a human-readable, self-contained JSON document—complete with topologies, spin couplings, lineshapes, form factors, and four-vector verification points—and to reconstruct and evaluate those models bit-exactly from the JSON file alone.

---

## Core Functionality

1. **Model Serialization (Writer)**
   - Serializes a `CascadeDecay` model directly (`amplitudeSerializationDict(model; ...)`, `serializeToDict(model; ...)`, `writeJson(path, model; ...)`), or a `CascadeSystem` / `SystemSpins` paired with weighted `DecayChain`s.
   - Emits `kinematics` (`initial_state` and `final_state` with `name`, `index`, and `spin`), omitting hardcoded `mass` fields by default so event-dependent masses are determined by four-vectors.
   - Collects all referenced propagators and vertex form factors into the top-level `"functions"` array.

2. **Model Deserialization (Reader & Round-Trip Execution)**
   - Reads an amplitude-serialization JSON file via `readJson(path)` (or a parsed dictionary via `dict2instance(CascadeDecay, dict)`), returning `(model::CascadeDecay, system::CascadeSystem, document::LittleDict)`.
   - Reconstructs all topologies, LS/helicity/parity vertices, form factors, and standard or custom lineshapes so the computation can be restarted strictly from the JSON file.
   - Supports passing `custom_functions = Dict("name" => fn)` to bind custom callable functions during deserialization.

3. **Clean, Hash-Free Function Naming (`NamedLineshape`)**
   - Wrapping any lineshape in `NamedLineshape("lineshape_ResA", ls)` (or passing `propagator_names`) ensures that `"parametrization"` references in `"chains"` and definitions in `"functions"` use exact, human-readable names without structural hashes.
   - Form factors automatically receive clean deterministic names such as `BlattWeisskopf_L1_d3.0`.

4. **Standard & Custom Physics Lineshapes**
   - **Standard types**: `BreitWigner`, `MultichannelBreitWigner` (including multi-wave partial widths in a single decay channel), `ConstantLineshape`, `BlattWeisskopf`, `MomentumPower`, and `NoFormFactor`.
   - **Custom types (`"type": "custom"`)**: `NRExpLineshape` (non-resonant exponential), `TFPWAMultichannelBreitWigner` (TF-PWA running-width convention), and general custom expression containers.

5. **Validation & Four-Vector Verification Points**
   - Provides mutating helpers (`setValidation!`, `setParameterPoints!`, `setDomains!`, `setVariables!`, `setMisc!`, `appendFunction!`, `setSection!`) to attach phase-space domains, sampled four-vectors (`parameter_points`), and expected intensity/amplitude checksums (`misc.amplitude_model_checksums`).

---

## Setup & Running Tests

With `CascadeDecays.jl` and `CascadeDecaysIO.jl` cloned in the same workspace, activate and test `CascadeDecaysIO.jl`:

```julia
using Pkg
Pkg.activate("CascadeDecaysIO.jl")
Pkg.instantiate()
Pkg.test()
```

To run the complete standalone dummy analysis and JSON round-trip example from the terminal:

```powershell
julia --project=CascadeDecaysIO.jl CascadeDecaysIO.jl/examples/dummy_analysis.jl
```

---

## Step-by-Step Tutorial: Building, Writing, and Reading a Cascade Model

For a complete multi-topology script (`(((1, 2), 3), 4)` and `((1, 2), (3, 4))`) with multiple partial waves and 3 four-vector verification events, see [`examples/dummy_analysis.jl`](examples/dummy_analysis.jl). Below is a complete walkthrough of the package workflow.

### Step 1: Define External Spins, Topology, and Named Lineshapes

Spin arguments in `SystemSpins`, `Propagator`, and `RecouplingLS` use doubled integer units (`2J`, `2L`, `2S`) internally, which `CascadeDecaysIO.jl` automatically converts to standard angular-momentum units (`J`, `L`, `S`) in the JSON.

```julia
using CascadeDecays
using CascadeDecays.ThreeBodyDecays: RecouplingLS
using CascadeDecaysIO
using FourVectors
using HadronicLineshapes

# External spins for P0(J=0) -> p1(0) p2(0) p3(0) p4(0)
spins    = SystemSpins(0, 0, 0, 0; two_h0 = 0)
topology = DecayTopology((((1, 2), 3), 4))

# NamedLineshape assigns clean, hash-free names in the JSON "parametrization" and "functions" sections
ls_v12  = NamedLineshape("constant_V12", ConstantLineshape(1.0 + 0.0im))
ls_resA = NamedLineshape(
    "lineshape_ResA",
    MultichannelBreitWigner(
        3.50,
        [
            (; gsq = 0.45, ma = 1.80, mb = 1.50, l = 0, d = 3.0),
            (; gsq = 0.20, ma = 1.80, mb = 1.50, l = 2, d = 3.0),
        ],
    ),
)
```

### Step 2: Build a `DecayChain` and `CascadeDecay` Model

```julia
chain = DecayChain(
    topology,
    spins;
    propagators = (
        (1, 2)      => Propagator(2, ls_v12),  # J = 1 intermediate state on (1, 2)
        ((1, 2), 3) => Propagator(2, ls_resA), # J = 1 resonance on ((1, 2), 3)
    ),
    vertices = (
        (((1, 2), 3), 4) => Vertex(RecouplingLS((2, 2)), BlattWeisskopf{1}(3.0)), # L=1, S=1
        ((1, 2), 3)      => Vertex(RecouplingLS((0, 2)), BlattWeisskopf{0}(3.0)), # L=0, S=1
        (1, 2)           => Vertex(RecouplingLS((2, 0))),                         # L=1, S=0
    ),
)

model = CascadeDecay(
    (chain,),
    topology;
    couplings = (0.85 - 0.15im,),
    names     = ("ResA_L1_d0",),
)
```

### Step 3: Evaluate on a Four-Vector Event & Build the Serialization Document

```julia
# Evaluate complex amplitude A and unpolarized intensity |A|^2 on a 4-vector point
task = KinematicTask((topology,))
ev   = (
    p1 = FourVector( 0.42, -0.31,  0.18; E = 1.52),
    p2 = FourVector(-0.08,  0.14, -0.05; E = 0.35),
    p3 = FourVector(-0.55,  0.29, -0.41; E = 1.65),
    p4 = FourVector( 0.21, -0.12,  0.28; E = 0.68),
)
pt  = KinematicPoint(task, (ev.p1, ev.p2, ev.p3, ev.p4))
amp = only(amplitude(model, pt))
val = abs2(amp) # unpolarized intensity |A|^2

# Build the JSON-ready dictionary
document = amplitudeSerializationDict(
    model;
    particle_labels  = ("p1", "p2", "p3", "p4", "P0"),
    name             = "minimal_cascade_model",
    variables        = ["m1_2sq", "m1_2_3sq"],
    parameter_points = [
        Dict(
            "name" => "verification_point_1",
            "four_vectors" => Dict(
                "p1" => Dict("E" => ev.p1.E, "px" => ev.p1.px, "py" => ev.p1.py, "pz" => ev.p1.pz),
                "p2" => Dict("E" => ev.p2.E, "px" => ev.p2.px, "py" => ev.p2.py, "pz" => ev.p2.pz),
                "p3" => Dict("E" => ev.p3.E, "px" => ev.p3.px, "py" => ev.p3.py, "pz" => ev.p3.pz),
                "p4" => Dict("E" => ev.p4.E, "px" => ev.p4.px, "py" => ev.p4.py, "pz" => ev.p4.pz),
            ),
        ),
    ],
    misc = Dict(
        "generator" => "CascadeDecaysIO.jl",
        "amplitude_model_checksums" => [
            Dict(
                "distribution"   => "minimal_cascade_model",
                "point"          => "verification_point_1",
                "value"          => val,
                "amplitude_real" => real(amp),
                "amplitude_imag" => imag(amp),
            ),
        ],
    ),
)

# Write formatted JSON to disk
writeJson("minimal_cascade_model.json", document)
```

### Step 4: Read the JSON File and Restart the Calculation

Using `readJson`, the entire `CascadeDecay` model is reconstructed strictly from the JSON file:

```julia
reloaded_model, reloaded_system, reloaded_doc = readJson("minimal_cascade_model.json")

# Re-evaluate at the verification point stored in the JSON
fv = reloaded_doc["parameter_points"][1]["four_vectors"]
q1 = FourVector(fv["p1"]["px"], fv["p1"]["py"], fv["p1"]["pz"]; E = fv["p1"]["E"])
q2 = FourVector(fv["p2"]["px"], fv["p2"]["py"], fv["p2"]["pz"]; E = fv["p2"]["E"])
q3 = FourVector(fv["p3"]["px"], fv["p3"]["py"], fv["p3"]["pz"]; E = fv["p3"]["E"])
q4 = FourVector(fv["p4"]["px"], fv["p4"]["py"], fv["p4"]["pz"]; E = fv["p4"]["E"])

reloaded_topologies = Tuple(unique(ch.topology for ch in reloaded_model.chains))
reloaded_pt         = KinematicPoint(KinematicTask(reloaded_topologies), (q1, q2, q3, q4))
reloaded_amp        = only(amplitude(reloaded_model, reloaded_pt))

@assert reloaded_amp == amp
```

---

## JSON Document Structure

A serialized document built with `amplitudeSerializationDict` contains:

```json
{
  "distributions": [
    {
      "type": "HadronicUnpolarizedIntensity",
      "name": "minimal_cascade_model",
      "variables": ["m1_2sq", "m1_2_3sq"],
      "decay_description": {
        "kinematics": {
          "initial_state": {"name": "P0", "index": 0, "spin": 0},
          "final_state": [
            {"name": "p1", "index": 1, "spin": 0},
            {"name": "p2", "index": 2, "spin": 0},
            {"name": "p3", "index": 3, "spin": 0},
            {"name": "p4", "index": 4, "spin": 0}
          ]
        },
        "reference_topology": [[[1, 2], 3], 4],
        "chains": [
          {
            "name": "ResA_L1_d0",
            "weight": "0.85 - 0.15i",
            "topology": [[[1, 2], 3], 4],
            "propagators": [
              {"node": [1, 2], "spin": 1, "parametrization": "constant_V12"},
              {"node": [[1, 2], 3], "spin": 1, "parametrization": "lineshape_ResA"}
            ],
            "vertices": [
              {"node": [[[1, 2], 3], 4], "type": "ls", "l": 1, "s": 1, "formfactor": "BlattWeisskopf_L1_d3.0"},
              {"node": [[1, 2], 3], "type": "ls", "l": 0, "s": 1, "formfactor": "BlattWeisskopf_L0_d3.0"},
              {"node": [1, 2], "type": "ls", "l": 1, "s": 0, "formfactor": ""}
            ]
          }
        ]
      }
    }
  ],
  "functions": [
    {
      "type": "BlattWeisskopf",
      "l": 1,
      "radius": 3.0,
      "name": "BlattWeisskopf_L1_d3.0"
    },
    {
      "type": "ConstantLineshape",
      "value": "1.0 + 0.0i",
      "name": "constant_V12"
    }
  ],
  "parameter_points": [],
  "misc": {
    "amplitude_model_checksums": []
  }
}
```

### Mutating Document Setters

Optional sections can be added or updated on an existing document dictionary at any time:

```julia
setVariables!(document, ["m1_2sq", "m1_2_3sq"])
setDomains!(document, ["phase_space"])
setParameterPoints!(document, parameter_points)
setMisc!(document, Dict("analysis" => "B2DxDK"))
setValidation!(
    document,
    Dict(
        "misc" => Dict("amplitude_model_checksums" => checksums),
        "parameter_points" => parameter_points,
    ),
)
appendFunction!(document, "constant_one", Dict("type" => "ConstantLineshape", "value" => 1.0))
setSection!(document, "custom_section", Dict("kept" => true))
```

---

## Supported Types & Functions

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

---

## Building Documentation

The documentation site uses Quarto for the workflow tutorial (`docs/writer_workflow.qmd`) and `Documenter.jl` for the HTML site:

```powershell
julia --project=docs -e "using Pkg; Pkg.instantiate()"
julia --project=docs docs/make.jl
```
