const _Appendix = Dict{String, Any}

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
        return LittleDict{String, Any}(string(k) => _normalize_json_value(v) for (k, v) in value)
    elseif value isa Union{AbstractVector, Tuple}
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

function serializeToDict(h::ThreeBodyDecays.RecouplingLS)
    l, s = _spin_string.(h.two_ls)
    return LittleDict{String, Any}("type" => "ls", "l" => l, "s" => s), _Appendix()
end

function serializeToDict(h::ThreeBodyDecays.NoRecoupling)
    helicities = [_spin_string(_lambda_a(h)), _spin_string(_lambda_b(h))]
    return LittleDict{String, Any}("type" => "helicity", "helicities" => helicities), _Appendix()
end

function serializeToDict(h::ThreeBodyDecays.ParityRecoupling)
    helicities = [_spin_string(_lambda_a(h)), _spin_string(_lambda_b(h))]
    return LittleDict{String, Any}(
            "type" => "parity",
            "helicities" => helicities,
            "parity_factor" => _parity_phase_is_plus(h) ? "+" : "-",
        ),
        _Appendix()
end

function serializeToDict(system::CascadeSystem; particle_labels = nothing)
    spins = system.quantum isa SystemSpinParities ? system.quantum.spins : system.quantum
    n_final = length(spins.finals)
    labels = isnothing(particle_labels) ? vcat(["p$i" for i in 1:n_final], ["X"]) : collect(particle_labels)
    length(labels) == n_final + 1 ||
        throw(ArgumentError("particle_labels must contain final-state labels followed by the initial-state label"))

    kinematics = LittleDict{String, Any}(
        "initial_state" => LittleDict{String, Any}(
            "index" => 0,
            "name" => labels[end],
            "spin" => _spin_string(spins.two_h0),
            "mass" => root_mass(system),
        ),
        "final_state" => [
            LittleDict{String, Any}(
                "index" => i,
                "name" => labels[i],
                "spin" => _spin_string(spins.finals[i]),
                "mass" => final_masses(system)[i],
            ) for i in 1:n_final
        ],
    )
    return kinematics, _Appendix()
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

phaseSpaceDomain(masses::SystemMasses; name::AbstractString = "default") = LittleDict{String, Any}(
    "name" => name,
    "type" => "phase_space",
    "m0" => masses.m0,
    "final_state_masses" => collect(masses.finals),
)

function setValidation!(document::AbstractDict, validation)
    points = get!(document, "parameter_points", Any[])
    append!(points, _normalize_json_value(validation.parameter_points))
    misc = get!(document, "misc", LittleDict{String, Any}())
    checksums = get!(misc, "amplitude_model_checksums", Any[])
    append!(checksums, _normalize_json_value(validation.checksums))
    return document
end


_json_topology(topology::DecayTopology) = _json_node(_line_address(topology, root_line_ind(topology)))
