const _COMPLEX_RE = r"^([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\s*([+-])\s*([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)i$"
const _IMAG_RE = r"^([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)i$"

function _json_to_tuple(node)
    if node isa AbstractVector
        return Tuple(_json_to_tuple(x) for x in node)
    elseif node isa Integer
        return Int(node)
    else
        return node
    end
end

function _spin_to_two_j(s)
    if s isa Integer
        return 2 * Int(s)
    elseif s isa Real
        return round(Int, 2 * s)
    elseif s isa AbstractString
        str = strip(s)
        slash_idx = findfirst('/', str)
        if slash_idx !== nothing
            num = parse(Float64, @view str[1:prevind(str, slash_idx)])
            den = parse(Float64, @view str[nextind(str, slash_idx):end])
            return round(Int, 2 * (num / den))
        else
            return round(Int, 2 * parse(Float64, str))
        end
    end
    error("Could not parse spin from: $s")
end

function _parse_complex(s)
    s isa Number && return ComplexF64(s)
    str = strip(string(s))
    # Fast path for "re + im*i" or "re - im*i" without regex allocations
    if endswith(str, 'i')
        body = @view str[1:prevind(str, lastindex(str))]
        split_idx = nothing
        idx = lastindex(body)
        first_idx = firstindex(body)
        while idx > first_idx
            c = body[idx]
            if c == '+' || c == '-'
                prev_c = body[prevind(body, idx)]
                if prev_c != 'e' && prev_c != 'E'
                    split_idx = idx
                    break
                end
            end
            idx = prevind(body, idx)
        end
        if split_idx !== nothing
            re_str = strip(@view body[1:prevind(body, split_idx)])
            im_str = strip(@view body[nextind(body, split_idx):end])
            re_val = tryparse(Float64, re_str)
            im_val = tryparse(Float64, im_str)
            if re_val !== nothing && im_val !== nothing
                sign = body[split_idx] == '-' ? -1.0 : 1.0
                return ComplexF64(re_val, sign * im_val)
            end
        else
            im_only = tryparse(Float64, strip(body))
            im_only !== nothing && return ComplexF64(0.0, im_only)
        end
    else
        re_only = tryparse(Float64, str)
        re_only !== nothing && return ComplexF64(re_only, 0.0)
    end

    # Fallback regex path
    m = match(_COMPLEX_RE, str)
    if m !== nothing
        re = parse(Float64, m.captures[1])
        sign = m.captures[2] == "-" ? -1.0 : 1.0
        im = parse(Float64, m.captures[3])
        return ComplexF64(re, sign * im)
    end
    m_im = match(_IMAG_RE, replace(str, " " => ""))
    if m_im !== nothing
        return ComplexF64(0.0, parse(Float64, m_im.captures[1]))
    end
    error("Could not parse complex number: $str")
end

"""
    CustomExpressionLineshape

Fallback container for custom analytic lineshape expressions read from amplitude-serialization files.
"""
struct CustomExpressionLineshape <: HadronicLineshapes.AbstractFlexFunc
    expression::String
    parameters::Dict{String,Any}
end
@inline (ls::CustomExpressionLineshape)(σ) = 1.0 + 0.0im

function _parse_channels(raw_channels)
    return [
        (;
            gsq = Float64(ch["gsq"]),
            ma = Float64(ch["ma"]),
            mb = Float64(ch["mb"]),
            l = Int(ch["l"]),
            d = Float64(ch["d"]),
        ) for ch in raw_channels
    ]
end

function _build_function_workspace(functions_list; custom_functions = Dict{String,Any}())
    workspace = Dict{String,Any}()
    sizehint!(workspace, length(functions_list))
    for fn in functions_list
        name = fn["name"]
        if haskey(custom_functions, name)
            workspace[name] = custom_functions[name]
            continue
        end
        fn_type = get(fn, "type", "")
        subtype = get(fn, "subtype", "")
        if fn_type == "BlattWeisskopf"
            l = Int(fn["l"])
            r = Float64(fn["radius"])
            workspace[name] = HadronicLineshapes.BlattWeisskopf{l}(r)
        elseif fn_type == "MomentumPower"
            l = Int(fn["l"])
            workspace[name] = HadronicLineshapes.MomentumPower{l}()
        elseif fn_type == "ConstantLineshape"
            val = _parse_complex(fn["value"])
            workspace[name] = ConstantLineshape(val)
        elseif fn_type == "BreitWigner"
            m = Float64(fn["mass"])
            w = Float64(fn["width"])
            ma = Float64(fn["ma"])
            mb = Float64(fn["mb"])
            l = Int(fn["l"])
            d = Float64(fn["d"])
            workspace[name] = HadronicLineshapes.BreitWigner(m, w, ma, mb, l, d)
        elseif fn_type == "MultichannelBreitWigner" || (fn_type == "custom" && subtype == "MultichannelBreitWigner")
            m = Float64(fn["mass"])
            workspace[name] = HadronicLineshapes.MultichannelBreitWigner(m, _parse_channels(fn["channels"]))
        elseif fn_type == "TFPWAMultichannelBreitWigner" || (fn_type == "custom" && (subtype == "TFPWAMultichannelBreitWigner" || haskey(fn, "channels")))
            m = Float64(fn["mass"])
            workspace[name] = TFPWAMultichannelBreitWigner(m, _parse_channels(fn["channels"]))
        elseif fn_type == "NRExpLineshape" || (fn_type == "custom" && (subtype == "NRExpLineshape" || (haskey(fn, "alpha") && haskey(fn, "beta"))))
            alpha = Float64(fn["alpha"])
            beta = Float64(fn["beta"])
            m0 = Float64(fn["m0"])
            workspace[name] = NRExpLineshape(alpha + 1im * beta, m0)
        elseif fn_type == "custom"
            expr = get(fn, "expression", "")
            workspace[name] = CustomExpressionLineshape(expr, Dict{String,Any}())
        else
            workspace[name] = fn
        end
    end
    return workspace
end

@inline function _extract_kinematics_dict(dict::AbstractDict)
    haskey(dict, "kinematics") && return dict["kinematics"]
    haskey(dict, "decay_description") && return dict["decay_description"]["kinematics"]
    haskey(dict, "distributions") && return dict["distributions"][1]["decay_description"]["kinematics"]
    return dict
end

@inline function _sorted_final_states(kin::AbstractDict)
    raw_finals = kin["final_state"]
    issorted(raw_finals, by = x -> x["index"]) && return raw_finals
    return sort(collect(raw_finals), by = x -> x["index"])
end

"""
    dict2instance(::Type{CascadeSystem}, dict::AbstractDict)

Deserialize a `CascadeSystem` from an amplitude-serialization dictionary.
If masses are absent from kinematics, they default to 0.0.
"""
function dict2instance(::Type{CascadeSystem}, dict::AbstractDict)
    kin = _extract_kinematics_dict(dict)
    ini = kin["initial_state"]
    finals = _sorted_final_states(kin)

    m0 = haskey(ini, "mass") ? Float64(ini["mass"]) : 0.0
    two_h0 = _spin_to_two_j(ini["spin"])

    final_m = ntuple(i -> haskey(finals[i], "mass") ? Float64(finals[i]["mass"]) : 0.0, length(finals))
    final_two_j = ntuple(i -> _spin_to_two_j(finals[i]["spin"]), length(finals))

    masses = SystemMasses(final_m...; m0 = m0)
    spins = SystemSpins(final_two_j...; two_h0 = two_h0)
    return CascadeSystem(spins, masses)
end

"""
    dict2instance(::Type{SystemSpins}, dict::AbstractDict)

Deserialize `SystemSpins` from an amplitude-serialization kinematics dictionary.
"""
function dict2instance(::Type{SystemSpins}, dict::AbstractDict)
    kin = _extract_kinematics_dict(dict)
    ini = kin["initial_state"]
    finals = _sorted_final_states(kin)
    two_h0 = _spin_to_two_j(ini["spin"])
    final_two_j = ntuple(i -> _spin_to_two_j(finals[i]["spin"]), length(finals))
    return SystemSpins(final_two_j...; two_h0 = two_h0)
end

"""
    dict2instance(::Type{CascadeDecay}, dict::AbstractDict; custom_functions=Dict(), workspace=nothing)

Deserialize a `CascadeDecay` model container from an amplitude-serialization document dictionary.
Custom functions can be passed via `custom_functions` to be set by reasonable functions.
"""
function dict2instance(::Type{CascadeDecay}, dict::AbstractDict; custom_functions = Dict{String,Any}(), workspace = nothing)
    functions_list = get(dict, "functions", Any[])
    if workspace === nothing
        workspace = _build_function_workspace(functions_list; custom_functions)
    end

    decay_desc = if haskey(dict, "decay_description")
        dict["decay_description"]
    elseif haskey(dict, "distributions")
        dict["distributions"][1]["decay_description"]
    else
        dict
    end

    kin = decay_desc["kinematics"]
    spins = dict2instance(SystemSpins, kin)

    ref_topology = DecayTopology(_json_to_tuple(decay_desc["reference_topology"]))
    chains_data = decay_desc["chains"]
    n_chains = length(chains_data)

    chains = Vector{DecayChain}(undef, n_chains)
    couplings = Vector{ComplexF64}(undef, n_chains)
    names = Vector{String}(undef, n_chains)

    @inbounds for (idx, ch_data) in enumerate(chains_data)
        names[idx] = string(ch_data["name"])
        couplings[idx] = _parse_complex(ch_data["weight"])
        top = DecayTopology(_json_to_tuple(ch_data["topology"]))

        raw_props = ch_data["propagators"]
        props = ntuple(length(raw_props)) do j
            p = raw_props[j]
            node_tuple = _json_to_tuple(p["node"])
            spin_two_j = _spin_to_two_j(p["spin"])
            lineshape = workspace[p["parametrization"]]
            node_tuple => Propagator(spin_two_j, lineshape)
        end

        raw_verts = ch_data["vertices"]
        verts = ntuple(length(raw_verts)) do j
            v = raw_verts[j]
            node_tuple = _json_to_tuple(v["node"])
            ff_name = get(v, "formfactor", "")
            ff = (isempty(ff_name) || !haskey(workspace, ff_name)) ? ThreeBodyDecays.NoFormFactor() : workspace[ff_name]
            v_type = v["type"]
            recoupling = if v_type == "ls"
                ThreeBodyDecays.RecouplingLS((_spin_to_two_j(v["l"]), _spin_to_two_j(v["s"])))
            elseif v_type == "helicity"
                hels = v["helicities"]
                ThreeBodyDecays.NoRecoupling(_spin_to_two_j(hels[1]), _spin_to_two_j(hels[2]))
            elseif v_type == "parity"
                hels = v["helicities"]
                ThreeBodyDecays.ParityRecoupling(_spin_to_two_j(hels[1]), _spin_to_two_j(hels[2]), v["parity_factor"] == "+")
            else
                error("Unknown vertex recoupling type: $v_type")
            end
            node_tuple => Vertex(recoupling, ff)
        end

        chains[idx] = DecayChain(top, spins; propagators = props, vertices = verts)
    end

    return CascadeDecay(
        Tuple(chains),
        ref_topology;
        couplings = Tuple(couplings),
        names = Tuple(names),
    )
end

"""
    readJson(path::AbstractString; custom_functions=Dict{String,Any}())

Read an amplitude-serialization JSON file and return `(cascade::CascadeDecay, system::CascadeSystem, document::LittleDict)`.
`custom_functions` can be provided to set custom lineshapes to reasonable functions.
"""
function readJson(path::AbstractString; custom_functions = Dict{String,Any}())
    dict = JSON.parsefile(path; dicttype = LittleDict{String,Any})
    cascade = dict2instance(CascadeDecay, dict; custom_functions)
    system = dict2instance(CascadeSystem, dict)
    return cascade, system, dict
end
