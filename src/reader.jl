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
        if contains(str, "/")
            parts = split(str, "/")
            num = parse(Float64, parts[1])
            den = parse(Float64, parts[2])
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
    m = match(
        r"^([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\s*([+-])\s*([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)i$",
        str,
    )
    if m !== nothing
        re = parse(Float64, m.captures[1])
        sign = m.captures[2] == "-" ? -1.0 : 1.0
        im = parse(Float64, m.captures[3])
        return ComplexF64(re, sign * im)
    end
    m_im = match(r"^([+-]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)i$", replace(str, " " => ""))
    if m_im !== nothing
        return ComplexF64(0.0, parse(Float64, m_im.captures[1]))
    end
    m_re = tryparse(Float64, str)
    if m_re !== nothing
        return ComplexF64(m_re, 0.0)
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
(ls::CustomExpressionLineshape)(σ) = 1.0 + 0.0im

function _build_function_workspace(functions_list)
    workspace = Dict{String,Any}()
    for fn in functions_list
        name = fn["name"]
        fn_type = get(fn, "type", "")
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
        elseif fn_type == "MultichannelBreitWigner"
            m = Float64(fn["mass"])
            channels = [
                (;
                    gsq = Float64(ch["gsq"]),
                    ma = Float64(ch["ma"]),
                    mb = Float64(ch["mb"]),
                    l = Int(ch["l"]),
                    d = Float64(ch["d"]),
                ) for ch in fn["channels"]
            ]
            workspace[name] = HadronicLineshapes.MultichannelBreitWigner(m, channels)
        elseif fn_type == "TFPWAMultichannelBreitWigner"
            m = Float64(fn["mass"])
            channels = [
                (;
                    gsq = Float64(ch["gsq"]),
                    ma = Float64(ch["ma"]),
                    mb = Float64(ch["mb"]),
                    l = Int(ch["l"]),
                    d = Float64(ch["d"]),
                ) for ch in fn["channels"]
            ]
            workspace[name] = TFPWAMultichannelBreitWigner(m, channels)
        elseif fn_type == "NRExpLineshape"
            alpha = Float64(fn["alpha"])
            beta = Float64(fn["beta"])
            m0 = Float64(fn["m0"])
            workspace[name] = NRExpLineshape(alpha + 1im * beta, m0)
        elseif fn_type == "custom"
            expr = fn["expression"]
            workspace[name] = CustomExpressionLineshape(expr, Dict{String,Any}())
        else
            workspace[name] = fn
        end
    end
    return workspace
end

"""
    dict2instance(::Type{CascadeSystem}, dict::AbstractDict)

Deserialize a `CascadeSystem` (external kinematics: spins and masses) from an
amplitude-serialization dictionary.
"""
function dict2instance(::Type{CascadeSystem}, dict::AbstractDict)
    kin = if haskey(dict, "kinematics")
        dict["kinematics"]
    elseif haskey(dict, "decay_description")
        dict["decay_description"]["kinematics"]
    elseif haskey(dict, "distributions")
        dict["distributions"][1]["decay_description"]["kinematics"]
    else
        dict
    end

    ini = kin["initial_state"]
    finals = sort(collect(kin["final_state"]), by = x -> x["index"])

    m0 = Float64(ini["mass"])
    two_h0 = _spin_to_two_j(ini["spin"])

    final_m = [Float64(f["mass"]) for f in finals]
    final_two_j = [_spin_to_two_j(f["spin"]) for f in finals]

    masses = SystemMasses(final_m...; m0 = m0)
    spins = SystemSpins(final_two_j...; two_h0 = two_h0)
    return CascadeSystem(spins, masses)
end

"""
    dict2instance(::Type{CascadeDecay}, dict::AbstractDict; workspace=nothing)

Deserialize a `CascadeDecay` model container from an amplitude-serialization document dictionary.
"""
function dict2instance(::Type{CascadeDecay}, dict::AbstractDict; workspace = nothing)
    functions_list = get(dict, "functions", Any[])
    if workspace === nothing
        workspace = _build_function_workspace(functions_list)
    end

    decay_desc = if haskey(dict, "decay_description")
        dict["decay_description"]
    elseif haskey(dict, "distributions")
        dict["distributions"][1]["decay_description"]
    else
        dict
    end

    kin = decay_desc["kinematics"]
    system = dict2instance(CascadeSystem, kin)
    spins = system.quantum

    ref_topology = DecayTopology(_json_to_tuple(decay_desc["reference_topology"]))
    chains_data = decay_desc["chains"]

    chains = DecayChain[]
    couplings = ComplexF64[]
    names = String[]

    for ch_data in chains_data
        push!(names, string(ch_data["name"]))
        push!(couplings, _parse_complex(ch_data["weight"]))
        top = DecayTopology(_json_to_tuple(ch_data["topology"]))

        # Propagators
        props = Pair[]
        for p in ch_data["propagators"]
            node_tuple = _json_to_tuple(p["node"])
            spin_two_j = _spin_to_two_j(p["spin"])
            p_name = p["parametrization"]
            lineshape = workspace[p_name]
            push!(props, node_tuple => Propagator(spin_two_j, lineshape))
        end

        # Vertices
        verts = Pair[]
        for v in ch_data["vertices"]
            node_tuple = _json_to_tuple(v["node"])
            ff_name = get(v, "formfactor", "")
            ff = (ff_name == "" || !haskey(workspace, ff_name)) ? ThreeBodyDecays.NoFormFactor() : workspace[ff_name]
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
            push!(verts, node_tuple => Vertex(recoupling, ff))
        end

        chain = DecayChain(top, spins; propagators = Tuple(props), vertices = Tuple(verts))
        push!(chains, chain)
    end

    return CascadeDecay(
        Tuple(chains),
        ref_topology;
        couplings = Tuple(couplings),
        names = Tuple(names),
    )
end

"""
    readJson(path::AbstractString)

Read an amplitude-serialization JSON file and return `(cascade::CascadeDecay, system::CascadeSystem, document::LittleDict)`.
"""
function readJson(path::AbstractString)
    dict = JSON.parsefile(path; dicttype = LittleDict{String,Any})
    cascade = dict2instance(CascadeDecay, dict)
    system = dict2instance(CascadeSystem, dict)
    return cascade, system, dict
end
