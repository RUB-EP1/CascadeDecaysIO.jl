function _json_to_tuple(node)
    if node isa AbstractVector
        return Tuple(_json_to_tuple(x) for x in node)
    elseif node isa Integer
        return Int(node)
    else
        return node
    end
end

# Spins, helicities and s: strings "n" or "n/2" as in the standard (pattern of the
# schema's Spin); numbers are accepted when they are integer or half-integer
const _SPIN_PATTERN = r"^(?![-+]0)([-+]?\d+)(/2)?$"

function _spin_to_two_j(s)
    s isa Integer && return 2 * Int(s)
    if s isa Real
        isinteger(2s) || throw(ArgumentError("spin $s is not a multiple of 1/2"))
        return Int(2s)
    end
    s isa AbstractString || throw(ArgumentError("could not parse spin from $(repr(s))"))
    m = match(_SPIN_PATTERN, strip(s))
    m === nothing && throw(ArgumentError("spin $(repr(s)) is not of the form \"n\" or \"n/2\""))
    n = parse(Int, m.captures[1])
    return m.captures[2] === nothing ? 2n : n
end

# orbital angular momenta: non-negative integers, as number or string
function _parse_int(x)
    value = x isa Integer ? Int(x) : x isa AbstractString ? tryparse(Int, strip(x)) : nothing
    (value === nothing || value < 0) && throw(ArgumentError("expected a non-negative integer, got $(repr(x))"))
    return value
end

function _parse_complex(s)
    s isa Number && return ComplexF64(s)
    str = strip(string(s))
    if endswith(str, 'i')
        body = @view str[1:prevind(str, lastindex(str))]
        idx = lastindex(body)
        while idx > firstindex(body)
            c = body[idx]
            if (c == '+' || c == '-') && !(body[prevind(body, idx)] in ('e', 'E'))
                re_val = tryparse(Float64, strip(@view body[1:prevind(body, idx)]))
                im_val = tryparse(Float64, strip(@view body[nextind(body, idx):end]))
                if re_val !== nothing && im_val !== nothing
                    return ComplexF64(re_val, c == '-' ? -im_val : im_val)
                end
                break
            end
            idx = prevind(body, idx)
        end
        im_only = tryparse(Float64, strip(body))
        im_only !== nothing && return ComplexF64(0.0, im_only)
    else
        re_only = tryparse(Float64, str)
        re_only !== nothing && return ComplexF64(re_only, 0.0)
    end
    throw(ArgumentError("could not parse complex number $(repr(str))"))
end

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
CascadeDecays / HadronicLineshapes objects, `generic_function` to a
[`GenericFunction`](@ref). Anything else must be bound by name via
`custom_functions`; otherwise reading fails instead of guessing.
"""
function _build_function_workspace(functions_list; custom_functions = Dict{String,Any}())
    workspace = Dict{String,Any}()
    for fn in functions_list
        name = fn["name"]
        haskey(workspace, name) && throw(ArgumentError("function \"$name\" is defined twice"))
        # every function keeps its name, so that a model read from a file can be written again
        workspace[name] = NamedLineshape(name, _function_instance(fn; custom_functions))
    end
    return workspace
end

function _function_instance(fn; custom_functions = Dict{String,Any}())
    name = fn["name"]
    haskey(custom_functions, name) && return custom_functions[name]
    type = get(fn, "type", "")
    type == "BlattWeisskopf" && return HadronicLineshapes.BlattWeisskopf{_parse_int(fn["l"])}(Float64(fn["radius"]))
    type == "MomentumPower" && return HadronicLineshapes.MomentumPower{_parse_int(fn["l"])}()
    type == "BreitWigner" && return HadronicLineshapes.BreitWigner(
        Float64(fn["mass"]), Float64(fn["width"]), Float64(fn["ma"]), Float64(fn["mb"]),
        _parse_int(fn["l"]), Float64(fn["d"]))
    type == "MultichannelBreitWigner" && return HadronicLineshapes.MultichannelBreitWigner(Float64(fn["mass"]), _channels(fn))
    type == "Polynomial" && return _polynomial(fn)
    type == "generic_function" && return GenericFunction(fn["expression"])
    throw(ArgumentError("function \"$name\" has type \"$type\", which CascadeDecaysIO cannot evaluate; " *
        "bind it with custom_functions = Dict(\"$name\" => f)"))
end

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
    state["mass"] isa Real || throw(ArgumentError(
        "mass of \"$(state["name"])\" refers to a parameter ($(repr(state["mass"]))); only numbers are supported"))
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
evaluate must be bound by name in `custom_functions`. Every lineshape and form factor
is wrapped in a [`NamedLineshape`](@ref) carrying its name from the file.
"""
function dict2instance(::Type{CascadeDecay}, dict::AbstractDict;
        distribution = 1, custom_functions = Dict{String,Any}(), workspace = nothing)
    workspace = something(workspace, _build_function_workspace(get(dict, "functions", Any[]); custom_functions))

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
        props = Tuple(_json_to_tuple(p["node"]) => Propagator(_spin_to_two_j(p["spin"]), _lineshape(p, workspace))
                      for p in ch["propagators"])
        verts = Tuple(_json_to_tuple(v["node"]) => Vertex(_recoupling(v), _formfactor(v, workspace))
                      for v in ch["vertices"])
        DecayChain(DecayTopology(_json_to_tuple(ch["topology"])), spins; propagators = props, vertices = verts)
    end

    return CascadeDecay(
        Tuple(chains),
        ref_topology;
        couplings = Tuple(_parse_complex(ch["weight"]) for ch in chains_data),
        names = Tuple(string(ch["name"]) for ch in chains_data),
    )
end

"""
    readJson(path; distribution=1, custom_functions=Dict())

Read an amplitude-serialization JSON file and return `(cascade, system, document)` for
the distribution selected by index or name.
"""
function readJson(path::AbstractString; distribution = 1, custom_functions = Dict{String,Any}())
    dict = JSON.parsefile(path; dicttype = LittleDict{String,Any})
    cascade = dict2instance(CascadeDecay, dict; distribution, custom_functions)
    system = dict2instance(CascadeSystem, dict; distribution)
    return cascade, system, dict
end
