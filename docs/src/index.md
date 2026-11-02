# CascadeDecaysIO.jl

`CascadeDecaysIO.jl` serializes `CascadeDecays.jl` models to [amplitude-serialization](https://github.com/RUB-EP1/amplitude-serialization) JSON documents and reads them back, with the same structure as the reference models of the standard and four-vector validation points.

## What It Provides

- `amplitudeSerializationDict`, `serializeToDict` and `writeJson` for writing a `CascadeDecay` together with its external masses.
- `readJson` and `dict2instance` for rebuilding `CascadeDecay`, `CascadeSystem` and `SystemSpins`.
- `validationSection`, `validationEvents`, `phaseSpaceEvent`, `validate` and `kinematicTask` for generated four-vector validation points and intensity checksums.
- `NamedLineshape` for function names and `GenericFunction` for `generic_function` expressions; analysis packages serialize their own lineshapes by extending `serializeToDict`.
- Section setters (`setDomains!`, `setValidation!`, `setParameterPoints!`, `setMisc!`, `setFunctions!`, `appendFunction!`, `setSection!`).

## Documentation Map

- [Writer & Reader Tutorial](@ref) shows the path from a `CascadeDecays.jl` model to a JSON file and back.
- [Schema Notes](@ref) explains how cascade objects map onto the document layout.
- [API Reference](@ref) lists public types, functions, and docstrings.
