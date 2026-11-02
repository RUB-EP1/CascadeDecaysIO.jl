"""
    NamedLineshape(name, f)

Attach an explicit name to a propagator lineshape or a vertex form factor. The
serializer uses `name` as the entry in the `functions` list and as the reference in
`parametrization` / `formfactor`; the wrapped object is evaluated unchanged.

CascadeDecaysIO defines no lineshapes of its own: `f` is a CascadeDecays or
HadronicLineshapes object, a [`GenericFunction`](@ref), or an analysis type for
which the analysis package defines `CascadeDecaysIO.serializeToDict`.
"""
struct NamedLineshape{F} <: HadronicLineshapes.AbstractFlexFunc
    name::String
    lineshape::F
end

@inline (nl::NamedLineshape)(σ) = nl.lineshape(σ)
@inline (nl::NamedLineshape)(args...) = nl.lineshape(args...)
@inline Base.getproperty(nl::NamedLineshape, s::Symbol) =
    (s === :name || s === :lineshape) ? getfield(nl, s) : getproperty(getfield(nl, :lineshape), s)
