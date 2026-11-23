"""
    fourVectorPoint(name, four_vectors; domain="default")

Validation point in the four-vector format: `four_vectors` are the final-state
four-momenta (`FourVector`s) of the particles 1, 2, ..., written per particle index as
`{"index": i, "p": {"E", "px", "py", "pz"}}`. The masses and the initial state
`(m0, 0, 0, 0)` are those of the referenced `phase_space` domain.
"""
function fourVectorPoint(name::AbstractString, four_vectors; domain::AbstractString = "default")
    return LittleDict{String,Any}(
        "name" => name,
        "domain" => domain,
        "final_state" => [
            LittleDict{String,Any}("index" => i,
                "p" => LittleDict{String,Any}("E" => p.E, "px" => p.px, "py" => p.py, "pz" => p.pz))
            for (i, p) in enumerate(four_vectors)
        ],
    )
end

"""
    kinematicTask(cascade::CascadeDecay)

`KinematicTask` covering every chain topology of `cascade`, with its reference
topology and Wigner alignment for all final-state particles with spin.
"""
function kinematicTask(cascade::CascadeDecay)
    topologies = Tuple(unique(chain.topology for chain in cascade.chains))
    spins, _ = _cascade_spins_and_chains(cascade)
    spinful = Tuple(i for (i, two_j) in enumerate(spins.finals) if two_j != 0)
    return KinematicTask(topologies; reference_topology = cascade.reference_topology, wigner_finals = spinful)
end

# ---------------------------------------------------------------------------
# Phase-space events
# ---------------------------------------------------------------------------

_node_tree(topology::DecayTopology) = _line_address(topology, root_line_ind(topology))
_node_tree(tree::Tuple) = tree

_internal_nodes(::Integer) = ()
_internal_nodes(node::Tuple) = (node, _internal_nodes(node[1])..., _internal_nodes(node[2])...)

_threshold(node, masses, node_masses) =
    _minimum_mass(node[1], masses, node_masses) + _minimum_mass(node[2], masses, node_masses)

# lowest mass of a node: its fixed mass or the threshold of its decay products
function _minimum_mass(node, masses, node_masses)
    node isa Integer && return Float64(masses.finals[node])
    haskey(node_masses, node) && return Float64(node_masses[node])
    return _threshold(node, masses, node_masses)
end

function _check_node_masses(tree, masses, node_masses)
    nodes = (_internal_nodes(tree[1])..., _internal_nodes(tree[2])...)
    for (node, m) in node_masses
        node in nodes || throw(ArgumentError("node $node is not an internal node of the topology $tree"))
        m isa Real && m >= _threshold(node, masses, node_masses) || throw(ArgumentError(
            "mass of node $node must be a number above the threshold of its decay products, got $m"))
    end
    _threshold(tree, masses, node_masses) <= masses.m0 ||
        throw(ArgumentError("the node masses exceed the initial-state mass $(masses.m0) GeV"))
    return nothing
end

# fixed, or uniform between the minimum mass and `upper` (the mass left by the parent)
function _node_mass(rng, node, upper, masses, node_masses)
    lower = _minimum_mass(node, masses, node_masses)
    (node isa Integer || haskey(node_masses, node)) && return lower
    return lower + (upper - lower) * rand(rng)
end

# `p` from the rest frame of a system of mass `M` into the frame where it has momentum `P`
function _boost_from_rest(p, P, M)
    Pp = P.px * p.px + P.py * p.py + P.pz * p.pz
    k = Pp / (M * (P.E + M)) + p.E / M
    return FourVector(p.px + k * P.px, p.py + k * P.py, p.pz + k * P.pz; E = (P.E * p.E + Pp) / M)
end

# successive isotropic two-body decays of `node` (mass `M`, four-momentum `P`)
function _decay_node!(final_vectors, rng, node, P, M, masses, node_masses)
    node isa Integer && return (final_vectors[node] = P)
    a, b = node
    m_a = _node_mass(rng, a, M - _minimum_mass(b, masses, node_masses), masses, node_masses)
    m_b = _node_mass(rng, b, M - m_a, masses, node_masses)
    q = HadronicLineshapes.breakup(M, m_a, m_b)
    cosθ, ϕ = 2 * rand(rng) - 1, 2π * rand(rng)
    n = q .* (sqrt(1 - cosθ^2) * cos(ϕ), sqrt(1 - cosθ^2) * sin(ϕ), cosθ)
    _decay_node!(final_vectors, rng, a, _boost_from_rest(FourVector(n...; E = sqrt(q^2 + m_a^2)), P, M), m_a, masses, node_masses)
    _decay_node!(final_vectors, rng, b, _boost_from_rest(FourVector((.-n)...; E = sqrt(q^2 + m_b^2)), P, M), m_b, masses, node_masses)
    return final_vectors
end

"""
    phaseSpaceEvent(rng, masses, topology; node_masses=Dict())

One point of the `phase_space` domain of `masses`: the final-state four-vectors, in
index order, in the rest frame of the decaying particle. The event is built by
successive two-body decays along `topology` (a `DecayTopology` or a nested tuple such
as `(((1, 2), 3), 4)`); at every vertex the first child is emitted isotropically in the
rest frame of the parent. The mass of an internal node is `node_masses[node]` if given
(e.g. a particle with fixed mass), otherwise uniform between its kinematic limits.

Nodes are addressed as in the topology, e.g. `(1, 2)` or `((1, 2), 3)`. The events
serve as validation points; they are not distributed according to phase space.
"""
function phaseSpaceEvent(rng::AbstractRNG, masses::SystemMasses, topology; node_masses = Dict())
    tree = _node_tree(topology)
    _check_node_masses(tree, masses, node_masses)
    m0 = Float64(masses.m0)
    final_vectors = Vector{Any}(undef, length(masses.finals))
    _decay_node!(final_vectors, rng, tree, FourVector(0.0, 0.0, 0.0; E = m0), m0, masses, node_masses)
    return Tuple(final_vectors)
end

"""
    validationEvents(cascade, masses; n=3, node_masses=Dict(), topology=cascade.reference_topology,
                     rng=Xoshiro(1), min_relative_intensity=0.1, pilot=1000, max_tries=100_000)

`n` events of the phase space of `masses` (see [`phaseSpaceEvent`](@ref)) at which
`cascade` has a substantial intensity. An event is accepted if the unpolarized
intensity of `cascade` there is at least `min_relative_intensity` times the mean
intensity of `pilot` events drawn the same way.

A model whose intensity vanishes on all pilot events is an error. Its intensity is then
non-zero only on a part of the phase space the generator does not reach, typically a
lineshape narrower than the sampled mass range (e.g. a particle with fixed mass),
whose node mass has to be fixed in `node_masses`.
"""
function validationEvents(
    cascade::CascadeDecay,
    masses::SystemMasses;
    n::Integer = 3,
    node_masses = Dict(),
    topology = cascade.reference_topology,
    rng::AbstractRNG = Xoshiro(1),
    min_relative_intensity::Real = 0.1,
    pilot::Integer = 1000,
    max_tries::Integer = 100_000,
    task = kinematicTask(cascade),
)
    draw() = phaseSpaceEvent(rng, masses, topology; node_masses)
    intensity(event) = unpolarized_intensity(cascade, KinematicPoint(task, event))
    scale = sum(intensity(draw()) for _ in 1:pilot) / pilot
    scale > 0 && isfinite(scale) || throw(ArgumentError(
        "mean intensity $scale over $pilot events with node masses $node_masses: the model is non-zero only " *
        "where the generator does not sample, e.g. on a lineshape narrower than the sampled mass range; " *
        "fix that node mass with node_masses = Dict(node => mass)"))
    return map(1:n) do _
        for _ in 1:max_tries
            event = draw()
            intensity(event) >= min_relative_intensity * scale && return event
        end
        throw(ArgumentError("no event with intensity >= $(min_relative_intensity) x mean in $max_tries tries"))
    end
end

"""
    validationSection(cascade, events; distribution, domain="default", prefix="validation_point", task=kinematicTask(cascade))
    validationSection(cascade, masses::SystemMasses; distribution, n=3, kwargs...)

Validation entries for a list of events, each event a tuple of final-state
`FourVector`s. Returns `(; parameter_points, checksums)`: one four-vector point per
event and one checksum with the unpolarized intensity `Σ|A|²` of `cascade` at that
point. A point with vanishing intensity tests nothing and is an error. Pass the result
to [`setValidation!`](@ref) or to `amplitudeSerializationDict(...; validation)`.

With `masses` instead of events, `n` events are generated by [`validationEvents`](@ref),
which takes the remaining keywords (`node_masses`, `topology`, `rng`,
`min_relative_intensity`, ...).
"""
function validationSection(
    cascade::CascadeDecay,
    events;
    distribution::AbstractString,
    domain::AbstractString = "default",
    prefix::AbstractString = "validation_point",
    task = kinematicTask(cascade),
)
    parameter_points = Any[]
    checksums = Any[]
    for (n, event) in enumerate(events)
        name = "$prefix$n"
        intensity = unpolarized_intensity(cascade, KinematicPoint(task, Tuple(event)))
        intensity > 0 && isfinite(intensity) || throw(ArgumentError(
            "intensity $intensity at validation point $n: a validation point needs a non-zero amplitude"))
        push!(parameter_points, fourVectorPoint(name, event; domain))
        push!(checksums, LittleDict{String,Any}("point" => name, "distribution" => distribution, "value" => intensity))
    end
    return (; parameter_points, checksums)
end

function validationSection(
    cascade::CascadeDecay,
    masses::SystemMasses;
    distribution::AbstractString,
    domain::AbstractString = "default",
    prefix::AbstractString = "validation_point",
    task = kinematicTask(cascade),
    kwargs...,
)
    events = validationEvents(cascade, masses; task, kwargs...)
    return validationSection(cascade, events; distribution, domain, prefix, task)
end

_mass_squared(p) = p.E^2 - p.px^2 - p.py^2 - p.pz^2

_by_index(states) = sort(collect(states); by = state -> Int(state["index"]))

# final-state four-vectors of a point, in index order
function _four_vectors(point::AbstractDict)
    momenta = (state["p"] for state in _by_index(point["final_state"]))
    return Tuple(FourVector(Float64(p["px"]), Float64(p["py"]), Float64(p["pz"]); E = Float64(p["E"])) for p in momenta)
end

function _phase_space_domain(point, domains)
    name = get(point, "domain", nothing)
    haskey(domains, name) && domains[name]["type"] == "phase_space" ||
        throw(ArgumentError("point \"$(point["name"])\" must refer to a phase_space domain, got $(repr(name))"))
    return domains[name]
end

# deviation of a point from its phase_space domain
function _phase_space_offsets(point, domain::AbstractDict)
    states, finals = _by_index(point["final_state"]), _by_index(domain["final_state"])
    [state["index"] for state in states] == [state["index"] for state in finals] || throw(ArgumentError(
        "point \"$(point["name"])\" and domain \"$(domain["name"])\" list different final-state particles"))
    four_vectors = _four_vectors(point)
    # squared masses: well conditioned also for massless particles
    mass2_offset = maximum(abs(_mass_squared(p) - Float64(state["mass"])^2) for (p, state) in zip(four_vectors, finals))
    total(component) = sum(p -> getproperty(p, component), four_vectors)
    m0 = Float64(domain["initial_state"]["mass"])
    momentum_offset = max(abs(total(:E) - m0), abs(total(:px)), abs(total(:py)), abs(total(:pz)))
    return (; mass2_offset, momentum_offset)
end

"""
    validate(document; custom_functions=Dict(), rtol=1e-10, atol_phase_space=1e-8)
    validate(document, cascade; distribution=1, kwargs...)
    validate(path; kwargs...)

Recompute every entry of `misc.amplitude_model_checksums`: the unpolarized intensity
of the named distribution at the four-vector point. The point is also checked against
its `phase_space` domain: on-shell final states (|p² − m²| in GeV²) and total
four-momentum `(m0, 0, 0, 0)` (in GeV), both within `atol_phase_space`.

With `cascade` given, it is used for the distribution selected by `distribution`; all
other distributions are built from the document.
Returns one named tuple per checksum with the expected and computed value and `passed`.
"""
function validate(
    document::AbstractDict,
    cascade::Union{CascadeDecay, Nothing} = nothing;
    distribution = 1,
    custom_functions = Dict{String,Any}(),
    rtol::Real = 1e-10,
    atol_phase_space::Real = 1e-8,
)
    points = Dict(p["name"] => p for p in document["parameter_points"])
    domains = Dict(d["name"] => d for d in document["domains"])
    distribution_names = [d["name"] for d in document["distributions"]]
    workspace = _build_function_workspace(document["functions"]; custom_functions)
    models = Dict{String,Any}()
    if cascade !== nothing
        given = _select_distribution(document, distribution)["name"]
        models[given] = (cascade, kinematicTask(cascade))
    end
    model(name) = get!(models, name) do
        c = dict2instance(CascadeDecay, document; distribution = name, workspace)
        (c, kinematicTask(c))
    end

    results = []
    for checksum in document["misc"]["amplitude_model_checksums"]
        haskey(points, checksum["point"]) ||
            throw(ArgumentError("checksum refers to unknown point \"$(checksum["point"])\""))
        point = points[checksum["point"]]
        target = checksum["distribution"]
        expected = _parse_complex(checksum["value"])
        target in distribution_names ||
            throw(ArgumentError("checksum refers to \"$target\", which is not a distribution"))
        c, task = model(target)
        computed = unpolarized_intensity(c, KinematicPoint(task, _four_vectors(point)))
        offsets = _phase_space_offsets(point, _phase_space_domain(point, domains))
        on_domain = offsets.mass2_offset <= atol_phase_space && offsets.momentum_offset <= atol_phase_space
        relative_difference = abs(computed - expected) / max(abs(expected), floatmin(Float64))
        push!(results, (;
            point = checksum["point"],
            distribution = target,
            expected,
            computed,
            relative_difference,
            offsets.mass2_offset,
            offsets.momentum_offset,
            passed = relative_difference <= rtol && on_domain,
        ))
    end
    return results
end

validate(path::AbstractString; kwargs...) = validate(JSON.parsefile(path; dicttype = LittleDict{String,Any}); kwargs...)
