# Changelog

## 0.2.0

Breaking prototype upgrade: one IO package supports CascadeDecays and ThreeBodyDecays.
The CascadeDecaysIO package name and UUID are unchanged.

- Replace named physics wrappers with plain callables and a separate model IO context.
- Remove function-level `x`; generic expressions declare ordered `arguments` and separate `parameters`.
- Add checked expression evaluation, named parameter updates, explicit sharing and context remapping.
- Add validated document reading, backend selection and canonical writing with fresh optional physics checks.
- Update the schema and reference JSON models; legacy formats are intentionally unsupported.
- Split source and tests into focused files; add reusable CI and executable documentation.

Consumers must regenerate or update their model JSON and migrate to `read_document`,
`build_model`, `document` and `write_document`. See `docs/format.md` for the contract
and `docs/implementation-status.md` for remaining downstream integration work.
