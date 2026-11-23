# Schema Notes

`CascadeDecaysIO.jl` follows the [amplitude-serialization](https://github.com/RUB-EP1/amplitude-serialization) format and produces the same structure as its reference models. Events for validation are final-state four-vectors as introduced in [amplitude-serialization#111](https://github.com/RUB-EP1/amplitude-serialization/pull/111). Binary cascade topologies with any number of final-state particles are supported.

## Document Sections

`amplitudeSerializationDict` always writes the five sections the schema requires:

- `distributions`: one `HadronicUnpolarizedIntensity` with its `decay_description`.
- `functions`: named lineshapes and form factors referenced by the chains.
- `domains`: by default one `phase_space` domain, the only place of the masses: `initial_state` `{"index": 0, "mass"}` and `final_state` `[{"index", "mass"}, ...]`.
- `misc`: `amplitude_model_checksums` only.
- `parameter_points`: validation points with the final-state four-momenta.

## Decay Description

- `kinematics`: `initial_state` (index 0) and `final_state` entries with `index`, `name` and `spin`.
- `reference_topology`: nested-array cascade tree that fixes the frames of the external helicities.
- `chains`: `name`, `weight`, `topology`, `propagators`, `vertices`.

## Topology Nodes

Leaves are final-state indices; nested arrays are two-body subsystems:

```json
[[[1, 2], 3], 4]
```

## Propagators, Vertices and Functions

- A **propagator** sits on an internal node and references a lineshape of the node's invariant mass squared (`"parametrization"`). Its numerator is 1.
- A **vertex** sits on every node and holds the recoupling (`ls`, `helicity`, `parity`) and a `"formfactor"` reference (`""` for none). The form factor is evaluated at `(m_X², m_A², m_B²)` with the invariant masses of the event and is not normalized.
- `x` names the argument of `BreitWigner`, `MultichannelBreitWigner` and `Polynomial`: `m_12_sq` for node `[1, 2]`, `m_12_3_sq` for `[[1, 2], 3]`, `m_34_sq` for `[3, 4]`.
- `generic_function` expressions use one free variable, conventionally `sigma`.

## Spins and Angular Momenta

`CascadeDecays.jl` and `ThreeBodyDecays.jl` store `2J`, `2L`, `2S` and `2λ` internally. The JSON uses strings in ordinary units:

- `Propagator(2, ...)` ↔ `"spin": "1"`
- `RecouplingLS((4, 2))` ↔ `"l": "2"`, `"s": "1"`
- `NoRecoupling(1, -1)` ↔ `"helicities": ["1/2", "-1/2"]`

## Validation

- `parameter_points`: `{"name", "domain", "final_state": [{"index", "p": {"E", "px", "py", "pz"}}, ...]}`; particles are matched by index, the masses and the initial state `(m0, 0, 0, 0)` are those of the domain.
- `misc.amplitude_model_checksums`: `{"point", "distribution", "value"}` with the unpolarized intensity `Σ|A|²` of a distribution.
