const _FUNCTION_FIELDS = Dict(
    "BreitWigner" => ["mass", "width", "ma", "mb", "l", "d"],
    "MultichannelBreitWigner" => ["mass", "channels"],
    "Polynomial" => ["coefficients"],
    "BlattWeisskopf" => ["radius", "l"], "MomentumPower" => ["l"],
    "expression_function" => ["arguments", "parameters", "expression"]
)

_channels(fn) = [
    (; gsq = Float64(c["gsq"]), ma = Float64(c["ma"]), mb = Float64(c["mb"]), l = _parse_int(c["l"]), d = Float64(c["d"]))
        for c in collect(fn["channels"])
]

# constant -> CascadeDecays.ConstantLineshape; otherwise Polynomials.Polynomial (ascending powers of sigma)
function _polynomial(fn)
    coefficients = Float64.(fn["coefficients"])
    length(coefficients) == 1 && return ConstantLineshape(ComplexF64(only(coefficients)))
    return Polynomials.Polynomial(coefficients)
end

"""
    _build_function_workspace(functions; custom_functions)

Map every entry of the `functions` section to a callable: standard types to
CascadeDecays / HadronicLineshapes objects, `expression_function` to an
[`ExpressionFunction`](@ref). Anything else must be bound by name via
`custom_functions`; otherwise reading fails instead of guessing.
"""
function _build_function_workspace(functions_list; custom_functions = Dict{String, Any}())
    workspace = Dict{String, Any}()
    for fn in functions_list
        name = fn["name"]
        haskey(workspace, name) && throw(ArgumentError("function \"$name\" is defined twice"))
        # Names stay in this lookup table, never on the physics callable.
        workspace[name] = _function_instance(fn; custom_functions)
    end
    return workspace
end

function _function_instance(fn; custom_functions = Dict{String, Any}())
    name = fn["name"]
    haskey(custom_functions, name) && return custom_functions[name]
    type = get(fn, "type", "")
    type == "BlattWeisskopf" && return HadronicLineshapes.BlattWeisskopf{_parse_int(fn["l"])}(Float64(fn["radius"]))
    type == "MomentumPower" && return HadronicLineshapes.MomentumPower{_parse_int(fn["l"])}()
    type == "BreitWigner" && return HadronicLineshapes.BreitWigner(
        Float64(fn["mass"]), Float64(fn["width"]), Float64(fn["ma"]), Float64(fn["mb"]),
        _parse_int(fn["l"]), Float64(fn["d"])
    )
    type == "MultichannelBreitWigner" && return HadronicLineshapes.MultichannelBreitWigner(Float64(fn["mass"]), _channels(fn))
    type == "Polynomial" && return _polynomial(fn)
    type == "expression_function" && return ExpressionFunction(fn["expression"], fn["arguments"]; parameters = fn["parameters"])
    throw(
        ArgumentError(
            "function \"$name\" has type \"$type\", which CascadeDecaysIO cannot evaluate; " *
                "bind it with custom_functions = Dict(\"$name\" => f)"
        )
    )
end

"""
    serializeToDict(f)

Fallback for objects without an amplitude-serialization mapping. Lineshapes that
are not part of CascadeDecays or HadronicLineshapes must be serialized by the package
that defines them, by adding a method

```julia
CascadeDecaysIO.serializeToDict(f::MyLineshape) =
    serializeToDict(ExpressionFunction("1/(mass^2-sigma)", ["sigma"]; parameters=Dict("mass"=>f.mass)))
```

or be represented by an [`ExpressionFunction`](@ref).
"""
serializeToDict(f) = throw(
    ArgumentError(
        "no amplitude-serialization mapping for $(typeof(f)); define CascadeDecaysIO.serializeToDict " *
            "for it (built-in or expression_function) or use ExpressionFunction"
    )
)

function serializeToDict(lineshape::HadronicLineshapes.BreitWigner)
    return LittleDict{String, Any}(
            "type" => "BreitWigner",
            "mass" => lineshape.m,
            "width" => lineshape.Γ,
            "ma" => lineshape.ma,
            "mb" => lineshape.mb,
            "l" => lineshape.l,
            "d" => lineshape.d,
        ),
        _Appendix()
end

_channel_dict(c) = LittleDict{String, Any}(
    "gsq" => Float64(c.gsq), "ma" => Float64(c.ma), "mb" => Float64(c.mb), "l" => Int(c.l), "d" => Float64(c.d)
)

function serializeToDict(lineshape::HadronicLineshapes.MultichannelBreitWigner)
    return LittleDict{String, Any}(
            "type" => "MultichannelBreitWigner",
            "mass" => lineshape.m,
            "channels" => Any[_channel_dict(c) for c in lineshape.channels],
        ),
        _Appendix()
end

# HadronicLineshapes evaluates Flatte as a two-channel MultichannelBreitWigner with l = 0, d = 1
function serializeToDict(lineshape::HadronicLineshapes.Flatte)
    channels = (
        (; gsq = lineshape.gsq1, ma = lineshape.ma1, mb = lineshape.mb1, l = 0, d = 1.0),
        (; gsq = lineshape.gsq2, ma = lineshape.ma2, mb = lineshape.mb2, l = 0, d = 1.0),
    )
    return LittleDict{String, Any}(
            "type" => "MultichannelBreitWigner",
            "mass" => lineshape.m,
            "channels" => [_channel_dict(c) for c in channels],
        ),
        _Appendix()
end

function serializeToDict(lineshape::ConstantLineshape)
    value = lineshape.value
    iszero(imag(value)) || throw(
        ArgumentError(
            "ConstantLineshape($value) has a phase; move the constant into the chain weight"
        )
    )
    return LittleDict{String, Any}("type" => "Polynomial", "coefficients" => [Float64(real(value))]), _Appendix()
end

serializeToDict(p::Polynomials.Polynomial) =
    (LittleDict{String, Any}("type" => "Polynomial", "coefficients" => Float64.(Polynomials.coeffs(p))), _Appendix())

serializeToDict(f::ExpressionFunction) = (
    Dict{String, Any}(
        "type" => "expression_function", "expression" => f.expression,
        "arguments" => collect(f.arguments), "parameters" => copy(f.parameters)
    ), _Appendix(),
)

# --- vertex form factors ---

serializeToDict(ff::ThreeBodyDecays.NoFormFactor) = ("", _Appendix())

serializeToDict(ff::HadronicLineshapes.BlattWeisskopf) = (
    LittleDict{String, Any}(
        "type" => "BlattWeisskopf", "l" => HadronicLineshapes.orbital_momentum(ff), "radius" => ff.d
    ),
    _Appendix(),
)

serializeToDict(ff::HadronicLineshapes.MomentumPower) = (
    LittleDict{String, Any}("type" => "MomentumPower", "l" => HadronicLineshapes.orbital_momentum(ff)),
    _Appendix(),
)
