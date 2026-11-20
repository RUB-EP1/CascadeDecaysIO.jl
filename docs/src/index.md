# CascadeDecaysIO.jl

One amplitude-model JSON contract for CascadeDecays and ThreeBodyDecays.
The IO context holds names and parameter bindings separately from physics callables.

Version 0.2 is a prototype format redesign: generic functions declare ordered
`arguments` and separate `parameters`. There is no legacy reader or `x` binding.

See the [writer workflow](writer-workflow.md), [format contract](schema.md),
and [API reference](api-reference.md).
