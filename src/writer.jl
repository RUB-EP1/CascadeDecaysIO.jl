const _Appendix = Dict{String,Any}

_half_label(two_j::Integer) = iseven(two_j) ? div(two_j, 2) : two_j / 2
_complex_string(z::Complex) = replace(string(z), "im" => "i")
_json_scalar(x) = x isa Complex ? _complex_string(x) : x
_weight_string(x) = x isa Complex ? _complex_string(x) : string(x)
_gamma(obj) = getproperty(obj, Symbol(Char(0x0393)))
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
    else
        return _json_scalar(value)
    end
end

_node_label(address) = address isa Tuple ? join(_node_label.(address), "_") : string(address)
_mass_variable(address) = "m$(_node_label(address))sq"

function _line_address(topology::DecayTopology, line_ind::Integer)
    if isfinal_line_ind(topology, line_ind)
        return Int(line_ind)
    end
    vertex_ind = consumed_by(topology, line_ind)
    vertex_ind === nothing && error("could not recover bracket address for line $line_ind")
    children = child_line_inds(topology, vertex_ind)
    return (_line_address(topology, children[1]), _line_address(topology, children[2]))
end

function _function_name(prefix, address, payload)
    return "$(prefix)_$(_node_label(address))_$(hash(payload))"
end

_payload_name(prefix, payload) = "$(prefix)_$(hash(payload))"

function _merge_appendix!(target, source)
    merge!(target, source)
    return target
end

"""
    serializeToDict(system::CascadeSystem; particle_labels=nothing)

Serialize the external particles of a `CascadeSystem` in the same kinematics
shape used by ThreeBodyDecaysIO, generalized from three final states to `N`.

Returns `(kinematics, appendix)`, where `kinematics` contains
`"initial_state"` and `"final_state"` entries and `appendix` is empty.

# Keyword Arguments

- `particle_labels`: optional labels for all final-state particles followed by
  the initial-state label. If omitted, labels default to `p1`, `p2`, ..., `X`.

# Notes

Spin values are emitted in ordinary angular-momentum units. Internally,
`CascadeDecays.jl` stores doubled integer spins, so `two_h0 = 2` is emitted as
`"spin" => 1`.
"""
function serializeToDict(
    system::CascadeSystem;
    particle_labels = nothing,
)
    masses = system.masses
    spins = system.quantum isa SystemSpinParities ? system.quantum.spins : system.quantum
    n_final = length(masses.finals)
    labels = isnothing(particle_labels) ? vcat(["p$i" for i in 1:n_final], ["X"]) : collect(particle_labels)
    length(labels) == n_final + 1 ||
        throw(ArgumentError("particle_labels must contain final-state labels followed by the initial-state label"))

    system_dict = LittleDict{String,Any}(
        "initial_state" => LittleDict{String,Any}(
            "name" => labels[end],
            "mass" => masses.m0,
            "index" => 0,
            "spin" => _half_label(spins.two_h0),
        ),
        "final_state" => [
            LittleDict{String,Any}(
                "name" => labels[i],
                "mass" => masses.finals[i],
                "index" => i,
                "spin" => _half_label(spins.finals[i]),
            ) for i in 1:n_final
        ],
    )
    return system_dict, _Appendix()
end

"""
    serializeToDict(lineshape::ConstantLineshape)

Serialize the constant lineshape used by CascadeDecays examples and by
non-propagating placeholder lines in simple cascade tests.
"""
function serializeToDict(lineshape::ConstantLineshape)
    return LittleDict{String,Any}(
        "type" => "ConstantLineshape",
        "value" => _json_scalar(lineshape.value),
    ),
    _Appendix()
end

"""
    serializeToDict(ff::ThreeBodyDecays.NoFormFactor)

Serialize the absence of a vertex form factor.

Returns the empty form-factor reference `""` and an empty appendix.
"""
function serializeToDict(ff::ThreeBodyDecays.NoFormFactor)
    return "", _Appendix()
end

"""
    serializeToDict(ff::HadronicLineshapes.BlattWeisskopf)

Serialize a HadronicLineshapes Blatt-Weisskopf form factor.

The emitted dictionary contains `"type" => "BlattWeisskopf"`, orbital momentum
`"l"`, and radius `"radius"`.
"""
function serializeToDict(ff::HadronicLineshapes.BlattWeisskopf)
    return LittleDict{String,Any}(
        "type" => "BlattWeisskopf",
        "l" => HadronicLineshapes.orbital_momentum(ff),
        "radius" => ff.d,
    ),
    _Appendix()
end

"""
    serializeToDict(ff::HadronicLineshapes.MomentumPower)

Serialize a HadronicLineshapes threshold momentum-power form factor.
"""
function serializeToDict(ff::HadronicLineshapes.MomentumPower)
    return LittleDict{String,Any}(
        "type" => "MomentumPower",
        "l" => HadronicLineshapes.orbital_momentum(ff),
    ),
    _Appendix()
end

"""
    serializeToDict(h::ThreeBodyDecays.RecouplingLS)

Serialize an LS recoupling as a vertex dictionary fragment.

The emitted `"l"` and `"s"` values use ordinary angular-momentum units, not the
doubled integer convention used internally by `ThreeBodyDecays.jl`.
"""
function serializeToDict(h::ThreeBodyDecays.RecouplingLS)
    l, s = _half_label.(h.two_ls)
    return LittleDict{String,Any}("type" => "ls", "l" => l, "s" => s), _Appendix()
end

"""
    serializeToDict(h::ThreeBodyDecays.NoRecoupling)

Serialize a fixed-helicity vertex recoupling.
"""
function serializeToDict(h::ThreeBodyDecays.NoRecoupling)
    helicities = _half_label.((_lambda_a(h), _lambda_b(h)))
    return LittleDict{String,Any}("type" => "helicity", "helicities" => collect(helicities)), _Appendix()
end

"""
    serializeToDict(h::ThreeBodyDecays.ParityRecoupling)

Serialize a parity-constrained helicity recoupling.
"""
function serializeToDict(h::ThreeBodyDecays.ParityRecoupling)
    helicities = _half_label.((_lambda_a(h), _lambda_b(h)))
    parity_factor = _parity_phase_is_plus(h) ? '+' : '-'
    return LittleDict{String,Any}(
        "type" => "parity",
        "helicities" => collect(helicities),
        "parity_factor" => parity_factor,
    ),
    _Appendix()
end

function _serialize_formfactor(ff)
    ff_name, ff_appendix = serializeToDict(ff)
    ff_name == "" && return ff_name, ff_appendix

    name = _payload_name("formfactor", ff)
    ff_name["name"] = name
    return name, _Appendix(name => ff_name)
end

"""
    serializeToDict(vertex::ThreeBodyDecays.VertexFunction)

Serialize a `ThreeBodyDecays.VertexFunction` into a vertex dictionary fragment.

The returned dictionary contains the recoupling fields and a `"formfactor"`
reference. Non-trivial form factors are added to the returned appendix and
referenced by generated name.
"""
function serializeToDict(vertex::ThreeBodyDecays.VertexFunction)
    formfactor_name, formfactor_appendix = _serialize_formfactor(vertex.ff)
    recoupling_dict, recoupling_appendix = serializeToDict(vertex.h)
    vertex_dict = LittleDict{String,Any}("type" => recoupling_dict["type"], "formfactor" => formfactor_name)
    merge!(vertex_dict, recoupling_dict)
    appendix = _Appendix()
    _merge_appendix!(appendix, formfactor_appendix)
    _merge_appendix!(appendix, recoupling_appendix)
    return vertex_dict, appendix
end

function _serialize_lineshape(lineshape::HadronicLineshapes.BreitWigner, variable_name)
    return LittleDict{String,Any}(
        "type" => "BreitWigner",
        "mass" => lineshape.m,
        "width" => _gamma(lineshape),
        "ma" => lineshape.ma,
        "mb" => lineshape.mb,
        "l" => lineshape.l,
        "d" => lineshape.d,
        "x" => variable_name,
    ),
    _Appendix()
end

function _serialize_multichannel_bw(lineshape, variable_name; type_name = "MultichannelBreitWigner")
    channels = [
        LittleDict{String,Any}(
            "gsq" => _json_scalar(channel.gsq),
            "ma" => channel.ma,
            "mb" => channel.mb,
            "l" => channel.l,
            "d" => channel.d,
        ) for channel in lineshape.channels
    ]
    return LittleDict{String,Any}(
        "type" => type_name,
        "mass" => lineshape.m,
        "channels" => channels,
        "x" => variable_name,
    ),
    _Appendix()
end

_serialize_lineshape(lineshape::HadronicLineshapes.MultichannelBreitWigner, variable_name) =
    _serialize_multichannel_bw(lineshape, variable_name; type_name = "MultichannelBreitWigner")

_serialize_lineshape(lineshape::TFPWAMultichannelBreitWigner, variable_name) =
    _serialize_multichannel_bw(lineshape, variable_name; type_name = "TFPWAMultichannelBreitWigner")

function _serialize_lineshape(lineshape, variable_name)
    T = typeof(lineshape)
    if hasfield(T, :channels) && hasfield(T, :m)
        is_tfpwa = contains(string(nameof(T)), "TFPWA")
        return _serialize_multichannel_bw(lineshape, variable_name; type_name = is_tfpwa ? "TFPWAMultichannelBreitWigner" : "MultichannelBreitWigner")
    end
    if hasfield(T, :αβ) && hasfield(T, :m0)
        αβ = lineshape.αβ
        m0 = lineshape.m0
        α = real(αβ)
        β = imag(αβ)
        expr = "-exp(-(($α) + i*($β)) * ($variable_name - $(m0^2)))"
        return LittleDict{String,Any}(
            "type" => "NRExpLineshape",
            "expression" => expr,
            "alpha" => α,
            "beta" => β,
            "m0" => m0,
            "x" => variable_name,
        ),
        _Appendix()
    end
    if hasfield(T, :expression)
        return LittleDict{String,Any}(
            "type" => "custom",
            "expression" => string(lineshape.expression),
            "x" => variable_name,
        ),
        _Appendix()
    end
    throw(ArgumentError("unsupported lineshape type $(typeof(lineshape)) for variable $variable_name"))
end

function _serialize_named_lineshape(lineshape, address)
    name = _function_name("propagator", address, lineshape)
    if lineshape isa ConstantLineshape
        fn_dict, appendix = serializeToDict(lineshape)
    else
        fn_dict, appendix = _serialize_lineshape(lineshape, _mass_variable(address))
    end
    fn_dict["name"] = name
    return name, fn_dict, appendix
end

function _json_topology(topology::DecayTopology)
    return _json_node(_line_address(topology, root_line_ind(topology)))
end

function _json_topology(topology)
    return _json_node(topology)
end

function _serialize_propagator(topology::DecayTopology, lineshape, two_j::Integer, line_ind::Integer)
    address = _line_address(topology, line_ind)
    name, fn_dict, appendix = _serialize_named_lineshape(lineshape, address)
    prop_dict = LittleDict{String,Any}(
        "node" => _json_node(address),
        "spin" => _half_label(two_j),
        "parametrization" => name,
    )
    return prop_dict, name => fn_dict, appendix
end

function _serialize_vertex(topology::DecayTopology, vertex, vertex_ind::Integer)
    address = _line_address(topology, incoming_line_ind(topology, vertex_ind))
    vertex_dict, appendix = serializeToDict(vertex)
    vertex_dict["node"] = _json_node(address)
    return vertex_dict, appendix
end

"""
    serializeToDict(chain::DecayChain; name="cascade_chain")

Serialize one concrete CascadeDecays chain using the same chain fields as
ThreeBodyDecaysIO, generalized to all internal propagator lines and all binary
vertices of the cascade topology.

Returns `(chain_dict, appendix)`.

The `chain_dict` contains:

- `"name"`: the chain name.
- `"topology"`: nested array representation of the cascade topology.
- `"propagators"`: one entry for each internal line.
- `"vertices"`: one entry for each binary decay vertex.

The `appendix` contains named function dictionaries referenced by propagators
and non-trivial form factors.

This method writes dictionaries only. It does not write JSON files and it does
not provide read-back support.
"""
function serializeToDict(chain::DecayChain; name::AbstractString = "cascade_chain")
    appendix = _Appendix()
    functions = LittleDict{String,Any}()
    topology = chain.topology

    propagator_entries = map(zip(chain.propagators, propagator_two_js(chain), propagating_line_inds(chain))) do (lineshape, two_j, line_ind)
        prop_dict, fn_pair, fn_appendix = _serialize_propagator(topology, lineshape, two_j, line_ind)
        functions[fn_pair.first] = fn_pair.second
        _merge_appendix!(appendix, fn_appendix)
        prop_dict
    end

    vertex_entries = map(1:nvertices(chain)) do vertex_ind
        vertex_dict, vertex_appendix = _serialize_vertex(topology, chain.vertices[vertex_ind], vertex_ind)
        _merge_appendix!(appendix, vertex_appendix)
        vertex_dict
    end

    merge!(appendix, functions)
    chain_dict = LittleDict{String,Any}(
        "vertices" => collect(vertex_entries),
        "propagators" => collect(propagator_entries),
        "topology" => _json_topology(topology),
        "name" => name,
    )
    return chain_dict, appendix
end

function _chain_name_weight_model(entry)
    name, payload = entry
    weight, chain = payload
    return string(name), weight, chain
end

"""
    serializeToDict(system::CascadeSystem, weighted_chains; particle_labels=nothing, reference_topology=nothing)

Serialize a weighted list of CascadeDecays chains into a decay-description
dictionary. `weighted_chains` should be an iterable of
`name => (weight, chain)` pairs.

Returns `(decay_description, appendix)`.

The `decay_description` contains `"kinematics"`, `"reference_topology"`, and
`"chains"`. Each chain weight is emitted as a string to match the
`ThreeBodyDecaysIO.jl` style for complex coefficients.

# Keyword Arguments

- `particle_labels`: optional final-state labels followed by the initial-state
  label.
- `reference_topology`: optional topology override. This may be a
  `DecayTopology` or a nested tuple/array node such as `(((1, 2), 3), 4)`.

# Limitations

This method writes dictionaries only. It does not write JSON files and it does
not provide read-back support.
"""
function serializeToDict(
    system::CascadeSystem,
    weighted_chains;
    particle_labels = nothing,
    reference_topology = nothing,
)
    appendix = _Appendix()
    kinematics, kin_appendix = serializeToDict(system; particle_labels)
    _merge_appendix!(appendix, kin_appendix)

    chain_entries = map(weighted_chains) do entry
        name, weight, chain = _chain_name_weight_model(entry)
        chain_dict, chain_appendix = serializeToDict(chain; name)
        _merge_appendix!(appendix, chain_appendix)
        chain_dict["weight"] = _weight_string(weight)
        chain_dict
    end

    chains = collect(chain_entries)
    isempty(chains) && throw(ArgumentError("weighted_chains must not be empty"))
    topology =
        isnothing(reference_topology) ? chains[1]["topology"] : _json_topology(reference_topology)

    decay_description = LittleDict{String,Any}(
        "kinematics" => kinematics,
        "reference_topology" => topology,
        "chains" => chains,
    )
    return decay_description, appendix
end

"""
    serializeToDict(cascade::CascadeDecay, masses::SystemMasses; particle_labels=nothing, reference_topology=nothing)

Serialize a concrete `CascadeDecay` and external `masses` to a decay-description dictionary.
"""
function serializeToDict(
    cascade::CascadeDecay,
    masses::SystemMasses;
    particle_labels = nothing,
    reference_topology = nothing,
)
    first_chain = first(cascade.chains)
    spins = SystemSpins(
        first_chain.line_two_js[final_line_inds(first_chain)]...;
        two_h0 = first_chain.line_two_js[root_line_ind(first_chain)],
    )
    system = CascadeSystem(spins, masses)
    weighted_chains = [
        cascade.names[i] => (cascade.couplings[i], cascade.chains[i])
        for i in eachindex(cascade.chains)
    ]
    ref_top = isnothing(reference_topology) ? cascade.reference_topology : reference_topology
    return serializeToDict(system, weighted_chains; particle_labels, reference_topology = ref_top)
end

function serializeToDict(
    cascade::CascadeDecay;
    masses::Union{SystemMasses, Nothing} = nothing,
    particle_labels = nothing,
    reference_topology = nothing,
)
    isnothing(masses) && throw(ArgumentError("`masses` must be provided when serializing a CascadeDecay"))
    return serializeToDict(cascade, masses; particle_labels, reference_topology)
end

function _as_vector(value)
    value isa AbstractVector && return value
    return [value]
end

function _named_function_dict(name, payload)
    fn = _normalize_json_value(payload)
    fn["name"] = get(fn, "name", name)
    return fn
end

function _appendix_functions(appendix)
    names = sort!(collect(keys(appendix)); by = string)
    return [_named_function_dict(string(name), appendix[name]) for name in names]
end

function _default_distribution(decay_description; name, type, variables)
    distribution = LittleDict{String,Any}(
        "type" => type,
        "name" => name,
        "decay_description" => decay_description,
    )
    isnothing(variables) || (distribution["variables"] = variables)
    return distribution
end

function _selected_distribution(document, distribution_index)
    haskey(document, "distributions") ||
        throw(ArgumentError("document has no distributions section"))
    distributions = document["distributions"]
    1 <= distribution_index <= length(distributions) ||
        throw(ArgumentError("distribution_index $distribution_index is outside the distributions section"))
    return distributions[distribution_index]
end

function _merge_parameter_points(existing, incoming)
    isnothing(existing) && return incoming
    if existing isa AbstractVector && incoming isa AbstractVector
        return vcat(existing, incoming)
    elseif existing isa AbstractDict && incoming isa AbstractDict
        merged = LittleDict{String,Any}(existing)
        merge!(merged, incoming)
        return merged
    end
    throw(ArgumentError("cannot merge validation parameter_points with an existing incompatible parameter_points section"))
end

"""
    setSection!(document, section, value)

Set or replace a top-level amplitude-serialization section in `document`.

This generic setter is intentionally permissive so that new schema sections can
be added without changing the writer. Prefer the more specific setters for
standard sections when available.
"""
function setSection!(document::AbstractDict, section::AbstractString, value)
    document[section] = _normalize_json_value(value)
    return document
end

"""
    setDistributions!(document, distributions)

Set the top-level `"distributions"` section. A single distribution dictionary is
accepted and wrapped in a one-element vector.
"""
function setDistributions!(document::AbstractDict, distributions)
    document["distributions"] = _normalize_json_value(_as_vector(distributions))
    return document
end

"""
    setFunctions!(document, functions_or_appendix)

Set the top-level `"functions"` section.

If `functions_or_appendix` is a dictionary keyed by function name, it is
converted to the amplitude-serialization list form. Entries without a `"name"`
field receive their dictionary key as the name. If a vector is supplied, it is
used directly.
"""
function setFunctions!(document::AbstractDict, functions_or_appendix)
    functions =
        if functions_or_appendix isa AbstractVector
            _normalize_json_value(functions_or_appendix)
        else
            _appendix_functions(functions_or_appendix)
        end
    document["functions"] = functions
    return document
end

"""
    appendFunction!(document, name, function_dict)

Append one named function to the top-level `"functions"` section.

The function dictionary is copied before insertion. If it does not already
contain a `"name"` field, `name` is inserted.
"""
function appendFunction!(document::AbstractDict, name::AbstractString, function_dict)
    functions = get!(document, "functions", Any[])
    push!(functions, _named_function_dict(name, function_dict))
    return document
end

"""
    setDomains!(document, domains)

Set the top-level `"domains"` section. The writer does not infer physical
boundaries automatically; callers should pass the domain definitions used by
their analysis.
"""
function setDomains!(document::AbstractDict, domains)
    document["domains"] = _normalize_json_value(domains)
    return document
end

"""
    setMisc!(document, misc)

Set the top-level `"misc"` section for free-form metadata such as software
versions, analysis labels, or checksums.
"""
function setMisc!(document::AbstractDict, misc)
    document["misc"] = _normalize_json_value(misc)
    return document
end

"""
    setParameterPoints!(document, parameter_points)

Set the top-level `"parameter_points"` section. This is where named numerical
parameter sets can be stored when they are available independently of the model
structure.
"""
function setParameterPoints!(document::AbstractDict, parameter_points)
    document["parameter_points"] = _normalize_json_value(parameter_points)
    return document
end

"""
    setVariables!(document, variables; distribution_index=1)

Set `"variables"` on one distribution entry.

Variables belong to a concrete distribution rather than to the top-level
document, so this setter modifies `document["distributions"][distribution_index]`.
"""
function setVariables!(
    document::AbstractDict,
    variables;
    distribution_index::Integer = 1,
)
    distribution = _selected_distribution(document, distribution_index)
    distribution["variables"] = _normalize_json_value(variables)
    return document
end

"""
    setValidation!(document, validation)

Insert a ThreeBodyDecaysIO-style validation bundle into `document`.

The expected shape is a dictionary that may contain `"misc"` and
`"parameter_points"`, as returned by ThreeBodyDecaysIO's `validation_section`.
`"misc"` is merged with the document's top-level `"misc"` section, and
`"parameter_points"` replaces the top-level `"parameter_points"` section.

For convenience, a dictionary containing `"amplitude_model_checksums"` directly
is treated as the `"misc"` payload.
"""
function setValidation!(document::AbstractDict, validation)
    validation_dict = _normalize_json_value(validation)

    if haskey(validation_dict, "misc")
        misc = get(document, "misc", LittleDict{String,Any}())
        merge!(misc, validation_dict["misc"])
        document["misc"] = misc
    elseif haskey(validation_dict, "amplitude_model_checksums")
        misc = get(document, "misc", LittleDict{String,Any}())
        misc["amplitude_model_checksums"] = validation_dict["amplitude_model_checksums"]
        document["misc"] = misc
    end

    if haskey(validation_dict, "parameter_points")
        existing_points = get(document, "parameter_points", nothing)
        document["parameter_points"] = _merge_parameter_points(existing_points, validation_dict["parameter_points"])
    end

    (haskey(validation_dict, "misc") ||
     haskey(validation_dict, "amplitude_model_checksums") ||
     haskey(validation_dict, "parameter_points")) ||
        throw(ArgumentError("validation must contain misc, amplitude_model_checksums, or parameter_points"))

    return document
end

"""
    amplitudeSerializationDict(
        system::CascadeSystem,
        weighted_chains;
        particle_labels=nothing,
        reference_topology=nothing,
        name="cascade_model",
        distribution_type="HadronicUnpolarizedIntensity",
        variables=nothing,
        domains=nothing,
        misc=nothing,
        parameter_points=nothing,
        validation=nothing,
    )

Build a complete amplitude-serialization-style JSON dictionary from a
`CascadeSystem` and weighted `DecayChain` collection.

This is the document-level wrapper around `serializeToDict`. It creates the
standard `"distributions"` section, converts the collected appendix into the
top-level `"functions"` list, and inserts optional sections when keyword values
are provided.

# Keyword Placement

- `variables`: attached to the generated distribution.
- `validation`: merged into top-level `"misc"` and `"parameter_points"`.
- `domains`: written as the top-level `"domains"` section.
- `misc`: written as the top-level `"misc"` section.
- `parameter_points`: written as the top-level `"parameter_points"` section.
"""
function amplitudeSerializationDict(
    system::CascadeSystem,
    weighted_chains;
    particle_labels = nothing,
    reference_topology = nothing,
    name::AbstractString = "cascade_model",
    distribution_type::AbstractString = "HadronicUnpolarizedIntensity",
    variables = nothing,
    domains = nothing,
    misc = nothing,
    parameter_points = nothing,
    validation = nothing,
)
    decay_description, appendix = serializeToDict(
        system,
        weighted_chains;
        particle_labels,
        reference_topology,
    )
    document = LittleDict{String,Any}(
        "distributions" => [
            _default_distribution(
                decay_description;
                name,
                type = distribution_type,
                variables,
            ),
        ],
        "functions" => _appendix_functions(appendix),
    )
    isnothing(domains) || setDomains!(document, domains)
    isnothing(misc) || setMisc!(document, misc)
    isnothing(parameter_points) || setParameterPoints!(document, parameter_points)
    isnothing(validation) || setValidation!(document, validation)
    return document
end

"""
    writeJson(path, document; indent=4)

Write an amplitude-serialization dictionary to `path` as formatted JSON and
return `path`.
"""
function writeJson(path::AbstractString, document::AbstractDict; indent::Integer = 4)
    open(path, "w") do io
        JSON.print(io, document, indent)
        println(io)
    end
    return path
end

"""
    writeJson(path, system::CascadeSystem, weighted_chains; kwargs...)

Build a complete amplitude-serialization dictionary with
`amplitudeSerializationDict(system, weighted_chains; kwargs...)`, write it to
`path`, and return `path`.
"""
function writeJson(
    path::AbstractString,
    system::CascadeSystem,
    weighted_chains;
    kwargs...,
)
    return writeJson(path, amplitudeSerializationDict(system, weighted_chains; kwargs...))
end

"""
    amplitudeSerializationDict(cascade::CascadeDecay, masses::SystemMasses; kwargs...)

Build a complete amplitude-serialization dictionary from a `CascadeDecay` container and its external `masses`.
"""
function amplitudeSerializationDict(
    cascade::CascadeDecay,
    masses::SystemMasses;
    particle_labels = nothing,
    reference_topology = nothing,
    kwargs...,
)
    first_chain = first(cascade.chains)
    spins = SystemSpins(
        first_chain.line_two_js[final_line_inds(first_chain)]...;
        two_h0 = first_chain.line_two_js[root_line_ind(first_chain)],
    )
    system = CascadeSystem(spins, masses)
    weighted_chains = [
        cascade.names[i] => (cascade.couplings[i], cascade.chains[i])
        for i in eachindex(cascade.chains)
    ]
    ref_top = isnothing(reference_topology) ? cascade.reference_topology : reference_topology
    return amplitudeSerializationDict(
        system,
        weighted_chains;
        particle_labels,
        reference_topology = ref_top,
        kwargs...,
    )
end

"""
    writeJson(path, cascade::CascadeDecay, masses::SystemMasses; kwargs...)

Build a complete amplitude-serialization dictionary for `cascade` and write it to `path`.
"""
function writeJson(
    path::AbstractString,
    cascade::CascadeDecay,
    masses::SystemMasses;
    kwargs...,
)
    return writeJson(path, amplitudeSerializationDict(cascade, masses; kwargs...))
end
