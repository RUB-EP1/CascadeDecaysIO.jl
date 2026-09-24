# CascadeDecaysIO.jl

`CascadeDecaysIO.jl` serializes and deserializes `CascadeDecays.jl` models to and from amplitude-serialization dictionaries and JSON files (`HSF-aje-2025-02`).

It supports full write-and-read round-trip workflows: exporting mass-free kinematics, multi-topology decay chains, propagators, vertices, standard and custom lineshapes, and four-vector verification points, and reconstructing live `CascadeDecay` objects via `readJson`.

## What It Provides

- `serializeToDict` and `amplitudeSerializationDict` for serializing `CascadeDecay` models or `(system, weighted_chains)` collections.
- `NamedLineshape` and `propagator_names` for clean, hash-free function naming in `"parametrization"` and `"functions"`.
- Section setters (`setDomains!`, `setVariables!`, `setValidation!`, `setParameterPoints!`, `setMisc!`, `appendFunction!`, `setSection!`).
- `writeJson` for formatted JSON export.
- `readJson` and `dict2instance` for reconstructing `CascadeDecay`, `CascadeSystem`, and `SystemSpins` from JSON files or dictionaries.

## Documentation Map

- [Serialization & Round-Trip Workflow](@ref) shows the recommended path from a `CascadeDecays.jl` model to a JSON file and back.
- [Schema Notes](@ref) explains how `CascadeDecaysIO.jl` maps cascade objects into the amplitude-serialization document layout.
- [API Reference](@ref) lists public types, functions, and docstrings.
