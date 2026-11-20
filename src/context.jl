_structure(decay) = Tuple{String, Tuple}[
    (chain["name"], _json_to_tuple(chain["topology"])) for chain in decay["chains"]
]

function _read_context(doc, distribution; parameter_point)
    decay = distribution["decay_description"]
    kinematics = decay["kinematics"]
    labels = [state["name"] for state in _sorted_final_states(kinematics)]
    push!(labels, kinematics["initial_state"]["name"])
    points = [point for point in doc["parameter_points"] if haskey(point, "parameters")]
    point_name = isempty(points) ? "nominal" : parameter_point === nothing ? only(points)["name"] : parameter_point
    context = ModelIOContext(
        structure = _structure(decay),
        labels = labels,
        distribution = distribution["name"],
        domain = only(doc["domains"])["name"],
        parameter_point = point_name,
        definitions = Dict(fn["name"] => deepcopy(fn) for fn in doc["functions"]),
    )
    for (index, chain) in enumerate(decay["chains"])
        for propagator in chain["propagators"]
            context.names[Site(index, propagator["node"], :propagator)] = propagator["parametrization"]
        end
        for vertex in chain["vertices"]
            isempty(vertex["formfactor"]) && continue
            context.names[Site(index, vertex["node"], :formfactor)] = vertex["formfactor"]
        end
    end
    return context
end

function _check_context(context, model, decay)
    if context.structure !== nothing && context.structure != _structure(decay)
        throw(ArgumentError("model structure changed; remap context explicitly or use a fresh context"))
    end
    valid_sites = Set(first(site) for site in _model_sites(model))
    issubset(keys(context.names), valid_sites) || throw(ArgumentError("context contains invalid sites"))
    return nothing
end

function _default_name(site, type, functions, reserved)
    base = "chain_$(site.chain)_$(site.role)_$(_node_label(site.node))_$type"
    name = base
    suffix = 1
    while name in reserved || haskey(functions, name)
        suffix += 1
        name = "$(base)_$suffix"
    end
    return name
end

function _emit!(functions, parameters, context, site, object, reserved)
    object isa ThreeBodyDecays.NoFormFactor && return ""
    definition, _ = serializeToDict(object)
    name = get(context.names, site) do
        _default_name(site, definition["type"], functions, reserved)
    end
    if haskey(context.definitions, name)
        template = context.definitions[name]
        template["type"] == definition["type"] ||
            throw(ArgumentError("function type changed at $site; assign a new name or context"))
        for field in _FUNCTION_FIELDS[definition["type"]]
            field in ("arguments", "expression") && continue
            definition[field] = _restore_refs(definition[field], template[field], parameters, "$name.$field")
        end
    end
    definition["name"] = name
    if haskey(functions, name) && functions[name] != definition
        throw(ArgumentError("shared function $name diverged at $site"))
    end
    functions[name] = definition
    return name
end

"Remap a context using a complete old-index => new-index chain permutation."
function remap_context(context::ModelIOContext, permutation)
    context.structure === nothing && throw(ArgumentError("context has no structure"))
    indices = collect(eachindex(context.structure))
    valid = sort(collect(keys(permutation))) == indices && sort(collect(values(permutation))) == indices
    valid || throw(ArgumentError("provide a complete chain permutation"))
    result = deepcopy(context)
    result.names = Dict(Site(permutation[s.chain], s.node, s.role) => name for (s, name) in context.names)
    for (old, new) in permutation
        result.structure[new] = context.structure[old]
    end
    return result
end
