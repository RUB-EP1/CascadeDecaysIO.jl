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
        m isa Real && m >= _threshold(node, masses, node_masses) || throw(
            ArgumentError(
                "mass of node $node must be a number above the threshold of its decay products, got $m"
            )
        )
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
    scale > 0 && isfinite(scale) || throw(
        ArgumentError(
            "mean intensity $scale over $pilot events with node masses $node_masses: the model is non-zero only " *
                "where the generator does not sample, e.g. on a lineshape narrower than the sampled mass range; " *
                "fix that node mass with node_masses = Dict(node => mass)"
        )
    )
    return map(1:n) do _
        for _ in 1:max_tries
            event = draw()
            intensity(event) >= min_relative_intensity * scale && return event
        end
        throw(ArgumentError("no event with intensity >= $(min_relative_intensity) x mean in $max_tries tries"))
    end
end
