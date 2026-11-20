"""
    _select_distribution(dict, distribution)

The distribution `distribution` (index or name) of a document. Only
`HadronicUnpolarizedIntensity` is supported.
"""
function _select_distribution(dict::AbstractDict, distribution)
    haskey(dict, "distributions") || return nothing
    distributions = dict["distributions"]
    selected = if distribution isa Integer
        1 <= distribution <= length(distributions) ||
            throw(ArgumentError("document has $(length(distributions)) distributions, not $distribution"))
        distributions[distribution]
    else
        index = findfirst(d -> d["name"] == distribution, distributions)
        index === nothing && throw(ArgumentError("document has no distribution named \"$distribution\""))
        distributions[index]
    end
    type = get(selected, "type", "")
    type == "HadronicUnpolarizedIntensity" ||
        throw(ArgumentError("distribution \"$(selected["name"])\" has type \"$type\"; only HadronicUnpolarizedIntensity is supported"))
    return selected
end

function _extract_kinematics_dict(dict::AbstractDict; distribution = 1)
    haskey(dict, "kinematics") && return dict["kinematics"]
    haskey(dict, "decay_description") && return dict["decay_description"]["kinematics"]
    haskey(dict, "distributions") && return _select_distribution(dict, distribution)["decay_description"]["kinematics"]
    return dict
end

@inline function _sorted_final_states(kin::AbstractDict)
    raw_finals = kin["final_state"]
    issorted(raw_finals, by = x -> x["index"]) && return raw_finals
    return sort(collect(raw_finals), by = x -> x["index"])
end

function _state_mass(state)
    haskey(state, "mass") || throw(ArgumentError("kinematics entry \"$(state["name"])\" has no mass"))
    state["mass"] isa Real || throw(
        ArgumentError(
            "mass of \"$(state["name"])\" refers to a parameter ($(repr(state["mass"]))); only numbers are supported"
        )
    )
    return Float64(state["mass"])
end

"""
    dict2instance(::Type{SystemSpins}, dict::AbstractDict)

Deserialize the external spins from an amplitude-serialization document or its
`kinematics` section.
"""
function dict2instance(::Type{SystemSpins}, dict::AbstractDict; distribution = 1)
    kin = _extract_kinematics_dict(dict; distribution)
    finals = _sorted_final_states(kin)
    two_h0 = _spin_to_two_j(kin["initial_state"]["spin"])
    final_two_j = ntuple(i -> _spin_to_two_j(finals[i]["spin"]), length(finals))
    return SystemSpins(final_two_j...; two_h0 = two_h0)
end

"""
    dict2instance(::Type{CascadeSystem}, dict::AbstractDict; distribution=1)

Deserialize spins and masses of the external particles. Masses are required.
"""
function dict2instance(::Type{CascadeSystem}, dict::AbstractDict; distribution = 1)
    kin = _extract_kinematics_dict(dict; distribution)
    finals = _sorted_final_states(kin)
    masses = SystemMasses(map(_state_mass, finals)...; m0 = _state_mass(kin["initial_state"]))
    return CascadeSystem(dict2instance(SystemSpins, kin), masses)
end

function _recoupling(v::AbstractDict)
    type = v["type"]
    if type == "ls"
        return ThreeBodyDecays.RecouplingLS((2 * _parse_int(v["l"]), _spin_to_two_j(v["s"])))
    elseif type == "helicity"
        h = v["helicities"]
        return ThreeBodyDecays.NoRecoupling(_spin_to_two_j(h[1]), _spin_to_two_j(h[2]))
    elseif type == "parity"
        h = v["helicities"]
        return ThreeBodyDecays.ParityRecoupling(_spin_to_two_j(h[1]), _spin_to_two_j(h[2]), v["parity_factor"] == "+")
    end
    throw(ArgumentError("unknown vertex type \"$type\""))
end

function _formfactor(v::AbstractDict, workspace)
    name = get(v, "formfactor", "")
    isempty(name) && return ThreeBodyDecays.NoFormFactor()
    haskey(workspace, name) || throw(ArgumentError("vertex $(v["node"]) refers to unknown form factor \"$name\""))
    return workspace[name]
end

function _lineshape(p::AbstractDict, workspace)
    name = p["parametrization"]
    haskey(workspace, name) || throw(ArgumentError("propagator $(p["node"]) refers to unknown function \"$name\""))
    return workspace[name]
end

"""
    dict2instance(::Type{CascadeDecay}, dict::AbstractDict; distribution=1, custom_functions=Dict(), workspace=nothing)

Deserialize a `CascadeDecay` from an amplitude-serialization document. `distribution`
selects the distribution by index or name. Functions of a type CascadeDecaysIO cannot
evaluate must be bound by name in `custom_functions`. Lineshapes and form factors are plain physics callables.
"""
function dict2instance(
        ::Type{CascadeDecay}, dict::AbstractDict;
        distribution = 1, custom_functions = Dict{String, Any}(), workspace = nothing
    )
    workspace = isnothing(workspace) ? _build_function_workspace(get(dict, "functions", Any[]); custom_functions) : workspace

    decay_desc = if haskey(dict, "decay_description")
        dict["decay_description"]
    elseif haskey(dict, "distributions")
        _select_distribution(dict, distribution)["decay_description"]
    else
        dict
    end

    spins = dict2instance(SystemSpins, decay_desc["kinematics"])
    ref_topology = DecayTopology(_json_to_tuple(decay_desc["reference_topology"]))
    chains_data = decay_desc["chains"]

    chains = map(chains_data) do ch
        props = Tuple(
            _json_to_tuple(p["node"]) => Propagator(_spin_to_two_j(p["spin"]), _lineshape(p, workspace))
                for p in ch["propagators"]
        )
        verts = Tuple(
            _json_to_tuple(v["node"]) => Vertex(_recoupling(v), _formfactor(v, workspace))
                for v in ch["vertices"]
        )
        DecayChain(DecayTopology(_json_to_tuple(ch["topology"])), spins; propagators = props, vertices = verts)
    end

    return CascadeDecay(
        Tuple(chains),
        ref_topology;
        couplings = Tuple(_parse_complex(ch["weight"]) for ch in chains_data),
        names = Tuple(string(ch["name"]) for ch in chains_data),
    )
end
