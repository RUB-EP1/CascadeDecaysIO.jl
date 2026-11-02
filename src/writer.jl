const _Appendix = Dict{String,Any}

# Spins, helicities, l and s are strings in the standard ("0", "1", "1/2", "-3/2")
_spin_string(two_j::Integer) = iseven(two_j) ? string(div(two_j, 2)) : string(two_j, "/2")

# Complex numbers as "a + bi", matching the ComplexNumber pattern of the schema
function _complex_string(z::Number)
    c = ComplexF64(z)
    return string(real(c), signbit(imag(c)) ? " - " : " + ", abs(imag(c)), "i")
end

_lambda_a(obj) = getproperty(obj, Symbol("two_", Char(0x03bb), "a"))
_lambda_b(obj) = getproperty(obj, Symbol("two_", Char(0x03bb), "b"))
_parity_phase_is_plus(obj) =
    getproperty(obj, Symbol(Char(0x03b7), Char(0x03b7), Char(0x03b7), "phaseisplus"))

function _json_node(address)
    address isa Tuple || return address
    return [_json_node(address[1]), _json_node(address[2])]
end

function _normalize_json_value(value)
    if value isa AbstractDict
        return LittleDict{String,Any}(string(k) => _normalize_json_value(v) for (k, v) in value)
    elseif value isa Union{AbstractVector,Tuple}
        return [_normalize_json_value(v) for v in value]
    elseif value isa Complex
        return _complex_string(value)
    else
        return value
    end
end

# Node label in the m_ij_sq convention: (1,2) -> "12", ((1,2),3) -> "12_3", ((1,2),(3,4)) -> "12_34"
function _node_label(address)
    address isa Integer && return string(address)
    a, b = address
    (a isa Integer && b isa Integer) && return string(a, b)
    return string(_node_label(a), "_", _node_label(b))
end

_mass_variable(address) = "m_$(_node_label(address))_sq"

function _line_address(topology::DecayTopology, line_ind::Integer)
    if isfinal_line_ind(topology, line_ind)
        return Int(line_ind)
    end
    vertex_ind = consumed_by(topology, line_ind)
    vertex_ind === nothing && error("could not recover bracket address for line $line_ind")
    children = child_line_inds(topology, vertex_ind)
    return (_line_address(topology, children[1]), _line_address(topology, children[2]))
end

# Function types that carry the argument name `x` in the reference models
const _TYPES_WITH_X = ("BreitWigner", "MultichannelBreitWigner", "Polynomial")

# `x` only labels the argument; a function used on several nodes keeps the `x` of its first use
_without_x(fn) = filter(kv -> first(kv) != "x", fn)

function _add_function!(appendix::AbstractDict, name::AbstractString, fn::AbstractDict)
    if haskey(appendix, name)
        _without_x(appendix[name]) == _without_x(fn) || throw(ArgumentError(
            "two different functions are named \"$name\"; give them distinct names with NamedLineshape"))
    else
        appendix[name] = fn
    end
    return appendix
end

function _merge_appendix!(target, source)
    for (name, fn) in source
        _add_function!(target, name, fn)
    end
    return target
end

# ---------------------------------------------------------------------------
# Functions: one method per standard type. Analysis packages extend
# `serializeToDict` for their own lineshape types (usually as `generic_function`).
# ---------------------------------------------------------------------------

"""
    serializeToDict(f)

Fallback for objects without an amplitude-serialization mapping. Lineshapes that
are not part of CascadeDecays or HadronicLineshapes must be serialized by the package
that defines them, by adding a method

```julia
CascadeDecaysIO.serializeToDict(f::MyLineshape) =
    (LittleDict{String,Any}("type" => "generic_function", "expression" => "..."), Dict{String,Any}())
```

or be wrapped in a [`GenericFunction`](@ref).
"""
serializeToDict(f) = throw(ArgumentError(
    "no amplitude-serialization mapping for $(typeof(f)); define CascadeDecaysIO.serializeToDict " *
    "for it (standard type or generic_function) or wrap it in GenericFunction"))

function serializeToDict(lineshape::HadronicLineshapes.BreitWigner)
    return LittleDict{String,Any}(
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

_channel_dict(c) = LittleDict{String,Any}(
    "gsq" => Float64(c.gsq), "ma" => Float64(c.ma), "mb" => Float64(c.mb), "l" => Int(c.l), "d" => Float64(c.d))

function serializeToDict(lineshape::HadronicLineshapes.MultichannelBreitWigner)
    return LittleDict{String,Any}(
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
    return LittleDict{String,Any}(
        "type" => "MultichannelBreitWigner",
        "mass" => lineshape.m,
        "channels" => [_channel_dict(c) for c in channels],
    ),
    _Appendix()
end

function serializeToDict(lineshape::ConstantLineshape)
    value = lineshape.value
    iszero(imag(value)) || throw(ArgumentError(
        "ConstantLineshape($value) has a phase; move the constant into the chain weight"))
    return LittleDict{String,Any}("type" => "Polynomial", "coefficients" => [Float64(real(value))]), _Appendix()
end

serializeToDict(p::Polynomials.Polynomial) =
    (LittleDict{String,Any}("type" => "Polynomial", "coefficients" => Float64.(Polynomials.coeffs(p))), _Appendix())

# expressions read with the Unicode names σ and π are written in ASCII
serializeToDict(f::GenericFunction) = (
    LittleDict{String,Any}("type" => "generic_function", "expression" => replace(f.expression, "σ" => "sigma", "π" => "pi")),
    _Appendix(),
)

# --- vertex form factors ---

serializeToDict(ff::ThreeBodyDecays.NoFormFactor) = ("", _Appendix())

serializeToDict(ff::HadronicLineshapes.BlattWeisskopf) = (
    LittleDict{String,Any}(
        "type" => "BlattWeisskopf", "l" => HadronicLineshapes.orbital_momentum(ff), "radius" => ff.d),
    _Appendix(),
)

serializeToDict(ff::HadronicLineshapes.MomentumPower) = (
    LittleDict{String,Any}("type" => "MomentumPower", "l" => HadronicLineshapes.orbital_momentum(ff)),
    _Appendix(),
)

_formfactor_name(ff::HadronicLineshapes.BlattWeisskopf) = "BlattWeisskopf_l$(HadronicLineshapes.orbital_momentum(ff))"
_formfactor_name(ff::HadronicLineshapes.MomentumPower) = "MomentumPower_l$(HadronicLineshapes.orbital_momentum(ff))"
_formfactor_name(ff::NamedLineshape) = ff.name
_formfactor_name(ff) = "formfactor_$(nameof(typeof(ff)))"

_unwrap(f::NamedLineshape) = f.lineshape
_unwrap(f) = f

function _serialize_formfactor(ff)
    ff isa ThreeBodyDecays.NoFormFactor && return "", _Appendix()
    fn, appendix = serializeToDict(_unwrap(ff))
    name = _formfactor_name(ff)
    fn["name"] = name
    return name, _add_function!(appendix, name, fn)
end

# --- vertex recouplings ---

function serializeToDict(h::ThreeBodyDecays.RecouplingLS)
    l, s = _spin_string.(h.two_ls)
    return LittleDict{String,Any}("type" => "ls", "l" => l, "s" => s), _Appendix()
end

function serializeToDict(h::ThreeBodyDecays.NoRecoupling)
    helicities = [_spin_string(_lambda_a(h)), _spin_string(_lambda_b(h))]
    return LittleDict{String,Any}("type" => "helicity", "helicities" => helicities), _Appendix()
end

function serializeToDict(h::ThreeBodyDecays.ParityRecoupling)
    helicities = [_spin_string(_lambda_a(h)), _spin_string(_lambda_b(h))]
    return LittleDict{String,Any}(
        "type" => "parity",
        "helicities" => helicities,
        "parity_factor" => _parity_phase_is_plus(h) ? "+" : "-",
    ),
    _Appendix()
end

"""
    serializeToDict(vertex::ThreeBodyDecays.VertexFunction)

Serialize a vertex: recoupling fields plus a `"formfactor"` reference (`""` when the
vertex has no form factor). Non-trivial form factors are returned in the appendix.
"""
function serializeToDict(vertex::ThreeBodyDecays.VertexFunction)
    formfactor_name, appendix = _serialize_formfactor(vertex.ff)
    recoupling_dict, _ = serializeToDict(vertex.h)
    vertex_dict = LittleDict{String,Any}("type" => recoupling_dict["type"])
    merge!(vertex_dict, recoupling_dict)
    vertex_dict["formfactor"] = formfactor_name
    return vertex_dict, appendix
end

# --- propagator lineshapes ---

_lineshape_name(lineshape::NamedLineshape, address) = lineshape.name
_lineshape_name(lineshape, address) = "propagator_$(_node_label(address))_$(nameof(typeof(lineshape)))"

function _serialize_lineshape(lineshape, address; name = nothing)
    fn, appendix = serializeToDict(_unwrap(lineshape))
    fn["type"] in _TYPES_WITH_X && (fn["x"] = _mass_variable(address))
    prop_name = isnothing(name) ? _lineshape_name(lineshape, address) : string(name)
    fn["name"] = prop_name
    return prop_name, _add_function!(appendix, prop_name, fn)
end

function _json_topology(topology::DecayTopology)
    return _json_node(_line_address(topology, root_line_ind(topology)))
end

_json_topology(topology) = _json_node(topology)

function _propagator_name(propagator_names, i, topology, line_ind)
    isnothing(propagator_names) && return nothing
    if propagator_names isa AbstractVector || propagator_names isa Tuple
        return (1 <= i <= length(propagator_names)) ? propagator_names[i] : nothing
    elseif propagator_names isa AbstractDict
        address = _line_address(topology, line_ind)
        return get(propagator_names, address, get(propagator_names, line_ind, nothing))
    end
    return string(propagator_names)
end

"""
    serializeToDict(chain::DecayChain; name="cascade_chain", propagator_names=nothing)

Serialize one chain: `"name"`, `"topology"`, one `"propagators"` entry per internal
line and one `"vertices"` entry per binary vertex. Returns `(chain_dict, appendix)`
where `appendix` maps function names to their definitions.
"""
function serializeToDict(
    chain::DecayChain;
    name::AbstractString = "cascade_chain",
    propagator_names = nothing,
)
    appendix = _Appendix()
    topology = chain.topology
    two_js = propagator_two_js(chain)
    line_inds = propagating_line_inds(chain)

    propagators = map(eachindex(chain.propagators)) do i
        line_ind = line_inds[i]
        address = _line_address(topology, line_ind)
        prop_name, fn_appendix = _serialize_lineshape(chain.propagators[i], address;
            name = _propagator_name(propagator_names, i, topology, line_ind))
        _merge_appendix!(appendix, fn_appendix)
        LittleDict{String,Any}(
            "node" => _json_node(address),
            "spin" => _spin_string(two_js[i]),
            "parametrization" => prop_name,
        )
    end

    vertices = map(1:nvertices(chain)) do vertex_ind
        vertex_dict, vertex_appendix = serializeToDict(chain.vertices[vertex_ind])
        _merge_appendix!(appendix, vertex_appendix)
        vertex_dict["node"] = _json_node(_line_address(topology, incoming_line_ind(topology, vertex_ind)))
        vertex_dict
    end

    chain_dict = LittleDict{String,Any}(
        "name" => name,
        "topology" => _json_topology(topology),
        "propagators" => collect(propagators),
        "vertices" => collect(vertices),
    )
    return chain_dict, appendix
end

# ---------------------------------------------------------------------------
# Kinematics and decay description
# ---------------------------------------------------------------------------

"""
    serializeToDict(system::CascadeSystem; particle_labels=nothing)

Serialize the external particles: name, index (0 for the decaying particle), spin
and mass. Returns `(kinematics, appendix)`.
"""
function serializeToDict(system::CascadeSystem; particle_labels = nothing)
    spins = system.quantum isa SystemSpinParities ? system.quantum.spins : system.quantum
    n_final = length(spins.finals)
    labels = isnothing(particle_labels) ? vcat(["p$i" for i in 1:n_final], ["X"]) : collect(particle_labels)
    length(labels) == n_final + 1 ||
        throw(ArgumentError("particle_labels must contain final-state labels followed by the initial-state label"))

    kinematics = LittleDict{String,Any}(
        "initial_state" => LittleDict{String,Any}(
            "index" => 0,
            "name" => labels[end],
            "spin" => _spin_string(spins.two_h0),
            "mass" => root_mass(system),
        ),
        "final_state" => [
            LittleDict{String,Any}(
                "index" => i,
                "name" => labels[i],
                "spin" => _spin_string(spins.finals[i]),
                "mass" => final_masses(system)[i],
            ) for i in 1:n_final
        ],
    )
    return kinematics, _Appendix()
end

serializeToDict(::Union{SystemSpins, SystemSpinParities}; kwargs...) = throw(ArgumentError(
    "the amplitude-serialization standard requires masses: use CascadeSystem(spins, SystemMasses(...))"))

"""
    serializeToDict(system::CascadeSystem, weighted_chains; particle_labels, reference_topology, propagator_names)

Serialize `name => (weight, chain)` pairs into a decay description with
`"kinematics"`, `"reference_topology"` and `"chains"`. Returns
`(decay_description, appendix)`.
"""
function serializeToDict(
    system::CascadeSystem,
    weighted_chains;
    particle_labels = nothing,
    reference_topology = nothing,
    propagator_names = nothing,
)
    kinematics, appendix = serializeToDict(system; particle_labels)
    chains = map(collect(weighted_chains)) do (name, (weight, chain))
        chain_names = propagator_names isa AbstractDict ? get(propagator_names, string(name), propagator_names) : propagator_names
        chain_dict, chain_appendix = serializeToDict(chain; name = string(name), propagator_names = chain_names)
        _merge_appendix!(appendix, chain_appendix)
        chain_dict["weight"] = _complex_string(weight)
        chain_dict
    end
    isempty(chains) && throw(ArgumentError("weighted_chains must not be empty"))
    topology = isnothing(reference_topology) ? chains[1]["topology"] : _json_topology(reference_topology)
    decay_description = LittleDict{String,Any}(
        "kinematics" => kinematics,
        "reference_topology" => topology,
        "chains" => chains,
    )
    return decay_description, appendix
end

function _cascade_spins_and_chains(cascade::CascadeDecay)
    first_chain = first(cascade.chains)
    finals_idx = final_line_inds(first_chain)
    spins = SystemSpins(
        ntuple(i -> first_chain.line_two_js[finals_idx[i]], length(finals_idx))...;
        two_h0 = first_chain.line_two_js[root_line_ind(first_chain)],
    )
    weighted_chains = [
        cascade.names[i] => (cascade.couplings[i], cascade.chains[i])
        for i in eachindex(cascade.chains)
    ]
    return spins, weighted_chains
end

"""
    serializeToDict(cascade::CascadeDecay, masses::SystemMasses; particle_labels, reference_topology, propagator_names)

Serialize a `CascadeDecay` to a decay description.
"""
function serializeToDict(
    cascade::CascadeDecay,
    masses::SystemMasses;
    reference_topology = nothing,
    kwargs...,
)
    spins, weighted_chains = _cascade_spins_and_chains(cascade)
    reference = isnothing(reference_topology) ? cascade.reference_topology : reference_topology
    return serializeToDict(CascadeSystem(spins, masses), weighted_chains; reference_topology = reference, kwargs...)
end

# ---------------------------------------------------------------------------
# Document
# ---------------------------------------------------------------------------

_appendix_functions(appendix) = [appendix[name] for name in sort!(collect(keys(appendix)))]

"""
    phaseSpaceDomain(masses::SystemMasses; name="default")

`phase_space` domain of the decay: initial-state mass `m0` and the final-state
masses in index order.
"""
phaseSpaceDomain(masses::SystemMasses; name::AbstractString = "default") = LittleDict{String,Any}(
    "name" => name,
    "type" => "phase_space",
    "m0" => masses.m0,
    "final_state_masses" => collect(masses.finals),
)

"""
    amplitudeSerializationDict(system::CascadeSystem, weighted_chains; kwargs...)
    amplitudeSerializationDict(cascade::CascadeDecay, masses::SystemMasses; kwargs...)

Build a complete amplitude-serialization document with the five sections the
standard requires: `distributions`, `functions`, `domains`, `misc` and
`parameter_points`.

# Keyword arguments
- `name`, `distribution_type`: name and type of the distribution.
- `particle_labels`: final-state labels followed by the initial-state label.
- `reference_topology`, `propagator_names`: see [`serializeToDict`](@ref).
- `domains`: defaults to the `phase_space` domain of the masses.
- `validation`: output of [`validationSection`](@ref); adds the four-vector points
  to `parameter_points` and the checksums to `misc.amplitude_model_checksums`. For a
  `cascade` the default `:generate` writes three validation points from
  `validationSection(cascade, masses)`; a model that is non-zero only at a fixed node
  mass needs `validation = validationSection(cascade, masses; distribution = name,
  node_masses)`. `nothing` writes no validation points.
"""
function amplitudeSerializationDict(
    system::CascadeSystem,
    weighted_chains;
    particle_labels = nothing,
    reference_topology = nothing,
    propagator_names = nothing,
    name::AbstractString = "cascade_model",
    distribution_type::AbstractString = "HadronicUnpolarizedIntensity",
    domains = [phaseSpaceDomain(system.masses)],
    validation = nothing,
)
    decay_description, appendix = serializeToDict(
        system, weighted_chains; particle_labels, reference_topology, propagator_names)
    document = LittleDict{String,Any}(
        "distributions" => [LittleDict{String,Any}(
            "name" => name,
            "type" => distribution_type,
            "decay_description" => decay_description,
        )],
        "functions" => _appendix_functions(appendix),
        "domains" => _normalize_json_value(collect(domains)),
        "misc" => LittleDict{String,Any}("amplitude_model_checksums" => Any[]),
        "parameter_points" => Any[],
    )
    isnothing(validation) || setValidation!(document, validation)
    return document
end

function amplitudeSerializationDict(
    cascade::CascadeDecay,
    masses::SystemMasses;
    reference_topology = nothing,
    name::AbstractString = "cascade_model",
    validation = :generate,
    kwargs...,
)
    spins, weighted_chains = _cascade_spins_and_chains(cascade)
    reference = isnothing(reference_topology) ? cascade.reference_topology : reference_topology
    generate = validation === :generate
    document = amplitudeSerializationDict(CascadeSystem(spins, masses), weighted_chains;
        reference_topology = reference, name, validation = generate ? nothing : validation, kwargs...)
    generate && setValidation!(document, validationSection(cascade, masses; distribution = name))
    return document
end

amplitudeSerializationDict(::CascadeDecay; kwargs...) = throw(ArgumentError(
    "the amplitude-serialization standard requires masses: amplitudeSerializationDict(cascade, SystemMasses(...))"))

# ---------------------------------------------------------------------------
# Section setters
# ---------------------------------------------------------------------------

"""
    setSection!(document, section, value)

Set or replace a top-level section. The standard allows only `distributions`,
`functions`, `domains`, `misc` and `parameter_points`.
"""
function setSection!(document::AbstractDict, section::AbstractString, value)
    document[section] = _normalize_json_value(value)
    return document
end

"""
    setDistributions!(document, distributions)

Set the `distributions` section; a single distribution is wrapped in a vector.
"""
function setDistributions!(document::AbstractDict, distributions)
    document["distributions"] = _normalize_json_value(distributions isa AbstractVector ? distributions : [distributions])
    return document
end

"""
    setFunctions!(document, functions_or_appendix)

Set the `functions` section from a list of function dictionaries or from an appendix
keyed by function name.
"""
function setFunctions!(document::AbstractDict, functions_or_appendix)
    document["functions"] = functions_or_appendix isa AbstractVector ?
        _normalize_json_value(functions_or_appendix) :
        _appendix_functions(LittleDict{String,Any}(string(k) => merge(_normalize_json_value(v), Dict("name" => string(k)))
            for (k, v) in functions_or_appendix))
    return document
end

"""
    appendFunction!(document, name, function_dict)

Append one function definition; it is an error if a different function of the same
name already exists.
"""
function appendFunction!(document::AbstractDict, name::AbstractString, function_dict)
    functions = get!(document, "functions", Any[])
    fn = _normalize_json_value(function_dict)
    fn["name"] = name
    existing = findfirst(f -> f["name"] == name, functions)
    if isnothing(existing)
        push!(functions, fn)
    elseif functions[existing] != fn
        throw(ArgumentError("two different functions are named \"$name\""))
    end
    return document
end

"""
    setDomains!(document, domains)

Set the `domains` section, e.g. `[phaseSpaceDomain(masses)]`.
"""
function setDomains!(document::AbstractDict, domains)
    document["domains"] = _normalize_json_value(collect(domains))
    return document
end

"""
    setMisc!(document, misc)

Replace the `misc` section. The reference models keep only
`amplitude_model_checksums` there.
"""
function setMisc!(document::AbstractDict, misc)
    document["misc"] = _normalize_json_value(misc)
    return document
end

"""
    setParameterPoints!(document, parameter_points)

Replace the `parameter_points` section.
"""
function setParameterPoints!(document::AbstractDict, parameter_points)
    document["parameter_points"] = _normalize_json_value(collect(parameter_points))
    return document
end

"""
    setValidation!(document, validation)

Add the output of [`validationSection`](@ref): its `parameter_points` are appended to
the `parameter_points` section and its `checksums` to
`misc.amplitude_model_checksums`.
"""
function setValidation!(document::AbstractDict, validation)
    points = get!(document, "parameter_points", Any[])
    append!(points, _normalize_json_value(validation.parameter_points))
    misc = get!(document, "misc", LittleDict{String,Any}())
    checksums = get!(misc, "amplitude_model_checksums", Any[])
    append!(checksums, _normalize_json_value(validation.checksums))
    return document
end

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

"""
    writeJson(path, document; indent=4, ascii_only=true)

Write an amplitude-serialization document to `path` and return `path`. With
`ascii_only` (default) the file must not contain non-ASCII characters, so that
readers in any language can parse names and expressions without Unicode handling.
"""
function writeJson(path::AbstractString, document::AbstractDict; indent::Integer = 4, ascii_only::Bool = true)
    text = sprint(io -> JSON.print(io, document, indent))
    if ascii_only && !isascii(text)
        offending = unique(filter(!isascii, collect(text)))
        throw(ArgumentError("document contains non-ASCII characters $(join(offending, " ")); " *
            "rename them or pass ascii_only = false"))
    end
    open(path, "w") do io
        print(io, text)
        println(io)
    end
    return path
end

"""
    writeJson(path, cascade::CascadeDecay, masses::SystemMasses; kwargs...)
    writeJson(path, system::CascadeSystem, weighted_chains; kwargs...)

Build the document with [`amplitudeSerializationDict`](@ref) and write it.
"""
writeJson(path::AbstractString, cascade::CascadeDecay, masses::SystemMasses; indent::Integer = 4, ascii_only::Bool = true, kwargs...) =
    writeJson(path, amplitudeSerializationDict(cascade, masses; kwargs...); indent, ascii_only)

writeJson(path::AbstractString, system::CascadeSystem, weighted_chains; indent::Integer = 4, ascii_only::Bool = true, kwargs...) =
    writeJson(path, amplitudeSerializationDict(system, weighted_chains; kwargs...); indent, ascii_only)
