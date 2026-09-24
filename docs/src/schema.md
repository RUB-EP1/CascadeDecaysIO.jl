# Schema Notes

<<<<<<< HEAD
`CascadeDecaysIO.jl` follows the [Amplitude Serialization Format](https://github.com/RUB-EP1/amplitude-serialization) and `ThreeBodyDecaysIO.jl` document conventions, extending the topology representation to binary cascade trees with arbitrary numbers of final-state particles.
=======
`CascadeDecaysIO.jl` follows the ComPWA Amplitude Serialization Format (`HSF-aje-2025-02`) and `ThreeBodyDecaysIO.jl` document conventions, extending the topology representation to binary cascade trees with arbitrary numbers of final-state particles.
>>>>>>> 7e2544c628140697fa09b46cc89ff46ffd07f705

## Document Sections

A complete document built by `amplitudeSerializationDict` contains:

- `distributions`: model entries such as `HadronicUnpolarizedIntensity`.
- `functions`: named lineshapes and form factors referenced by decay chains.

Additional sections can be supplied by keyword arguments or setter helpers:

- `domains`: allowed ranges for kinematic variables.
- `misc`: free-form metadata and `amplitude_model_checksums`.
- `parameter_points`: named numerical points (including `four_vectors`) used for verification.

## Decay Description

Each distribution contains a `decay_description` dictionary with:

- `kinematics`: `initial_state` and `final_state` particle names, indices, and spins (masses are omitted by default so event-dependent masses are provided via four-vectors).
- `reference_topology`: nested array representation of the reference cascade tree.
- `chains`: weighted chain descriptions (`name`, `weight`, `topology`, `propagators`, `vertices`).

## Topology Nodes

Topology nodes use final-state indices as leaves. Nested arrays describe intermediate two-body subsystems:

```json
[[[1, 2], 3], 4]
```

means that particles `1` and `2` form a subsystem, that subsystem combines with particle `3`, and the resulting subsystem combines with particle `4`.

## Spins, LS Values, and Function Names

`CascadeDecays.jl` and `ThreeBodyDecays.jl` store angular momenta internally as doubled integers (`2J`, `2L`, `2S`). `CascadeDecaysIO.jl` converts automatically between doubled integers and standard angular-momentum units (`J`, `L`, `S`) when writing and reading:

- `Propagator(2, ...)` $\leftrightarrow$ `"spin": 1`.
- `RecouplingLS((0, 2))` $\leftrightarrow$ `"l": 0`, `"s": 1`.
- Wrapping a lineshape in `NamedLineshape("lineshape_ResA", ls)` emits `"parametrization": "lineshape_ResA"` in the chain and `"name": "lineshape_ResA"` in `"functions"`.

## Validation & Verification Points

Validation entries are stored in top-level sections:

- `misc.amplitude_model_checksums` stores expected unpolarized intensities (`"value": |A|^2`) as well as `"amplitude_real"` and `"amplitude_imag"`.
- `parameter_points` stores the corresponding sampled `four_vectors` for each final-state particle.
