# Implementation checkpoint

## Implemented

- One package with Cascade and ThreeBody model adapters.
- Explicit expression arguments and parameter bindings; no named physics wrappers or x.
- Literal and named numerical parameters, shared-reference consistency checks.
- Site names, deterministic generated names, explicit chain-permutation remapping.
- Multi-distribution input with explicit backend/distribution selection.
- Canonical exports; stale validation data omitted or regenerated.
- New schema and three rewritten reference documents (six distributions total).
- Updated Julia documentation examples in the coordinated format branch.

## Physics findings

Cross-backend tests found a parent-spin reference-frame phase difference between
ThreeBodyDecays 0.16.0 and the pinned Cascade revision. The ThreeBody adapter now
translates coupling coordinates reversibly, with tests across all reference choices.
The physics packages themselves are unchanged. See format.md for the convention.

Comparing the revised documents against the format repository's original content
confirms that physics parameters, mathematical expressions, chains and couplings are
unchanged after accounting for the explicit syntax/name transformations. The fixtures
use independently supplied four-vector validation points and intensity checksums.

## Remaining integration work

- Review and publish the CascadeDecaysIO upgrade, then coordinate retirement of
  ThreeBodyDecaysIO. The original CascadeDecaysIO identity and history are retained.
- Update downstream analysis environments and their custom serializers. Some active
  analysis checkouts have uncommitted edits; this checkpoint does not modify them.
- Upgrade the external ampform-dpd Python reader. Its documentation execution is
  explicitly disabled on the prototype format branch; Python compatibility is not claimed.
- Render and execute the full plotting documentation environment. The core examples,
  schema, and physics regressions are tested separately.
- Evaluate interpreter performance on fitting workloads before promising throughput.

No compatibility adapters or historical readers were added. This checkpoint provides
the unified implementation and coordinated schema/reference changes; it is not a
claim that every downstream environment has already switched packages.

## Verified at this checkpoint

Julia 1.11.5, clean dependency installation, followed by
`Pkg.test(; julia_args=["--compile=min", "-O0"])`: 108 assertions passed.
The Documenter build passes, including the executable reader/writer tutorial.
CI is split between `Test.yml` and `ReusableTest.yml`, covering Julia 1.11 and
current stable on Linux, macOS and Windows; hosted CI has not run yet.
Independent JSON Schema validation accepts all three reference documents and rejects
function-level x. Both repository diffs pass whitespace checks.
