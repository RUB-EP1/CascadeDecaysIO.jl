# Prototype format 0.1

This version is an amplitude-serialization extension, not a claim of HS3 conformance.
It follows FlatPPL's distinction between a callable definition and its application.

A document has `format_version: "0.1"`, `functions`, `distributions`, `domains`,
`parameter_points`, and `misc`. Unknown semantic fields are rejected by the supported
schema. Function names and distribution names must each be unique in their category.
Function names are ASCII identifiers. Chain labels need not be unique.

A built-in definition contains `name`, `type` and its specified parameters.
Supported types: BreitWigner, MultichannelBreitWigner, Polynomial,
BlattWeisskopf, MomentumPower. A generic definition contains:

```json
{
  "name": "simple_BW",
  "type": "expression_function",
  "arguments": ["sigma"],
  "parameters": {"mass": 1.52, "width": 0.015},
  "expression": "1/(mass^2-sigma-complex(0,1)*mass*width)"
}
```

Arguments are ordered local names, nonempty and unique. Parameter names are local
and disjoint from arguments and reserved operators/constants. Built-in argument
interfaces are specified by type; they are not repeated in JSON. A propagator
applies its callable to the invariant of `node`; a form factor receives the squared
masses of the parent and the two ordered children. No `x` field is accepted.

Numeric function fields can contain finite real literals or string references to
named numerical parameters. References resolve in an explicitly selected parameter
point when more than one is present. Reference cycles and computed parameters are
not supported: the referenced entry must contain a numerical value. Polynomial
coefficients and multichannel fields support the same numeric-reference rule.
External particle masses currently remain numerical literals.

A shared parameter has one value across all its sites. A shared function name has
one definition across all applications. Neither sharing rule follows from accidental
numerical equality of unrelated parameters or objects.

## Expressions

Supported expressions use finite real numeric literals, local identifiers, PI,
EULER, parentheses, `+ - * / ^`, and calls to `exp log sqrt sin cos tan abs abs2
conj real imag complex`. All functions except `complex(re, im)` are unary.
Arithmetic follows scalar Julia semantics; division is floating-point. `log` and
`sqrt` promote to complex and use principal branches. Non-integral powers use the
principal complex power. Integral powers retain the usual integer-power meaning.
There are no implicit globals, inferred arguments, statement blocks, attributes,
indexing, user-code execution, or arbitrary Julia calls. Length and nesting limits
bound parsing. Expressions may be constant in any of their declared arguments.

Definitions are serializable data. The evaluator is currently an interpreter;
compilation/performance optimization can preserve this same interface later.

## Output and validation

Canonical exports sort object keys and function entries. Numeric formatting follows
JSON.jl; byte stability is tested within the supported runtime. Array ordering that
carries physics meaning is retained; adapter-owned vertex ordering is canonicalized.
Read/export/read/export is idempotent. Reproducing input formatting is not a goal.

One phase-space domain defines the initial and ordered final masses. Numerical
validation uses four-vector points and intensity checksums. Parameters and validation
points occupy distinct point records. Export does not retain stale checksums:
validation is either regenerated or omitted explicitly by the default `validation=false`.

## Backend reference-frame convention

The ThreeBodyDecays 0.16 adapter translates coupling coordinates for a half-integer
parent when a channel and the reference are 2 and 3 in either order. Its signed
principal-angle parent rotations differ from composition via channel 1 by a 2π
rotation, whose spin representation is `(-1)^(2j_parent)`. The adapter multiplies
those coupling coordinates by that sign on import and export. It does not alter
function definitions or the stored document weights. This is a reversible basis
convention translation, not a fitted correction. Regression tests compare the
same physical models in both backends at multiple points and reference choices.
