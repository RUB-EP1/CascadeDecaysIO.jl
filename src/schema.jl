"Require the supported fields, rejecting unknown semantic content."
function _check_fields(object, required, optional = String[])
    object isa AbstractDict || throw(ArgumentError("expected an object"))
    missing = setdiff(required, collect(keys(object)))
    isempty(missing) || throw(ArgumentError("missing fields: $missing"))
    extra = setdiff(collect(keys(object)), vcat(required, optional))
    isempty(extra) || throw(ArgumentError("unsupported fields: $extra"))
    return nothing
end

function _nodes(node)
    if node isa Integer && !(node isa Bool) && node > 0
        return Int[node], Any[]
    end
    node isa AbstractVector && length(node) == 2 ||
        throw(ArgumentError("topology must be a binary tree"))
    left, left_nodes = _nodes(node[1])
    right, right_nodes = _nodes(node[2])
    isempty(intersect(left, right)) ||
        throw(ArgumentError("duplicate final-state index in topology"))
    return vcat(left, right), vcat(Any[node], left_nodes, right_nodes)
end

function _signature(definition)
    type = definition["type"]
    type == "expression_function" && return length(definition["arguments"])
    return type in ("BlattWeisskopf", "MomentumPower") ? 3 : 1
end

function _validate_functions(definitions, parameters)
    functions = Dict{String, Any}()
    for definition in definitions
        name = definition["name"]
        _identifier(name) || throw(ArgumentError("invalid function name $name"))
        haskey(functions, name) && throw(ArgumentError("duplicate function $name"))
        resolved = _resolved_function(definition, parameters)
        _function_instance(resolved) # Check the numerical fields and expression grammar.
        functions[name] = resolved
    end
    return functions
end

function _validate_kinematics(kinematics, domain)
    _check_fields(kinematics, ["initial_state", "final_state"])
    initial = kinematics["initial_state"]
    finals = _sorted_final_states(kinematics)
    length(finals) >= 2 || throw(ArgumentError("at least two final particles required"))
    for state in vcat([initial], finals)
        _check_fields(state, ["index", "name", "spin", "mass"])
        _finite(state["mass"]) && state["mass"] >= 0 || throw(ArgumentError("invalid external mass"))
        _spin_to_two_j(state["spin"]) >= 0 || throw(ArgumentError("negative spin"))
        state["index"] isa Integer && !(state["index"] isa Bool) ||
            throw(ArgumentError("particle indices must be integers"))
    end
    initial["index"] == 0 || throw(ArgumentError("initial index must be zero"))
    [s["index"] for s in finals] == collect(1:length(finals)) || throw(ArgumentError("invalid final indices"))
    same_masses = domain["m0"] == initial["mass"] &&
        domain["final_state_masses"] == [s["mass"] for s in finals]
    same_masses || throw(ArgumentError("phase-space domain disagrees with kinematics"))
    return length(finals)
end

function _validate_vertex(vertex)
    type = vertex["type"]
    fields = if type == "ls"
        ["l", "s"]
    elseif type == "helicity"
        ["helicities"]
    elseif type == "parity"
        ["helicities", "parity_factor"]
    else
        throw(ArgumentError("unknown vertex type $type"))
    end
    _check_fields(vertex, vcat(["node", "type", "formfactor"], fields))
    if haskey(vertex, "helicities")
        length(vertex["helicities"]) == 2 || throw(ArgumentError("two helicities required"))
    end
    if type == "parity"
        vertex["parity_factor"] in ("+", "-") || throw(ArgumentError("invalid parity factor"))
    end
    _recoupling(vertex)
    return nothing
end

function _validate_sites(entries, expected_nodes, functions, role)
    propagator = role == :propagator
    reference_key = propagator ? "parametrization" : "formfactor"
    arity = propagator ? 1 : 3
    expected = Set(_json_to_tuple.(expected_nodes))
    seen = Set{Tuple}()
    for entry in entries
        if propagator
            _check_fields(entry, ["node", "spin", reference_key])
            _spin_to_two_j(entry["spin"]) >= 0 || throw(ArgumentError("negative propagator spin"))
        else
            _validate_vertex(entry)
        end
        node = _json_to_tuple(entry["node"])
        node in expected && !(node in seen) || throw(ArgumentError("invalid or duplicate $role node $node"))
        push!(seen, node)
        name = entry[reference_key]
        !propagator && name == "" && continue
        haskey(functions, name) || throw(ArgumentError("unresolved function $name"))
        _signature(functions[name]) == arity || throw(ArgumentError("$name has wrong arity at $node"))
    end
    seen == expected || throw(ArgumentError("missing $role nodes"))
    return nothing
end

function _validate_distribution(distribution, functions, domain)
    _check_fields(distribution, ["name", "type", "decay_description"])
    distribution["type"] == "HadronicUnpolarizedIntensity" || throw(ArgumentError("unsupported distribution type"))
    decay = distribution["decay_description"]
    _check_fields(decay, ["kinematics", "reference_topology", "chains"])
    n = _validate_kinematics(decay["kinematics"], domain)
    sort(first(_nodes(decay["reference_topology"]))) == collect(1:n) || throw(ArgumentError("invalid reference topology"))
    isempty(decay["chains"]) && throw(ArgumentError("empty decay model"))
    for chain in decay["chains"]
        _check_fields(chain, ["name", "weight", "topology", "propagators", "vertices"])
        isfinite(_parse_complex(chain["weight"])) || throw(ArgumentError("nonfinite coupling"))
        leaves, nodes = _nodes(chain["topology"])
        sort(leaves) == collect(1:n) || throw(ArgumentError("topology does not match final states"))
        _validate_sites(chain["propagators"], nodes[2:end], functions, :propagator)
        _validate_sites(chain["vertices"], nodes, functions, :formfactor)
    end
    return nothing
end

"""
    validate_document(document; parameter_point=nothing)

Validate supported fields, function definitions, parameter references, topology and
call-site arity. Shared definitions are checked once, then each distribution is checked
against that workspace. Return the input document without modifying it.
"""
function validate_document(doc; parameter_point = nothing)
    _check_fields(doc, ["format_version", "distributions", "functions", "domains", "parameter_points", "misc"])
    doc["format_version"] == FORMAT_VERSION || throw(ArgumentError("unsupported amplitude format version"))
    isempty(doc["distributions"]) && throw(ArgumentError("empty distributions"))
    names = [d["name"] for d in doc["distributions"]]
    allunique(names) || throw(ArgumentError("duplicate distribution names"))
    _validate_auxiliary(doc)
    parameters = _parameter_values(doc, parameter_point)
    functions = _validate_functions(doc["functions"], parameters)
    for distribution in doc["distributions"]
        _validate_distribution(distribution, functions, only(doc["domains"]))
    end
    return doc
end
