_backend(::ThreeBodyDecays.ThreeBodyDecay) = ThreeBodyBackend()

# ThreeBodyDecays uses signed principal-angle parent rotations. Transport between
# channels 2 and 3 differs from transport via canonical channel 1 by a 2π turn.
# A spin-j parent acquires (-1)^(2j). Cascade's full rotations use the consistent
# lift. Translate coupling coordinates at the adapter boundary, in both directions.
_threebody_phase(k, reference, two_h0) =
    ((k == 2 && reference == 3) || (k == 3 && reference == 2)) && isodd(two_h0) ? -1 : 1

function _threebody_k(topology)
    t = _json_to_tuple(topology)
    for k in 1:3
        i, j = ThreeBodyDecays.ij_from_k(k)
        t == ((i, j), k) && return k
    end
    throw(ArgumentError("ThreeBodyBackend requires cyclic ordered topology ((i,j),k)"))
end
function _build_model(doc, ::ThreeBodyBackend, workspace, context)
    dd = only(doc["distributions"])["decay_description"]
    kin = dd["kinematics"]; finals = _sorted_final_states(kin)
    length(finals) == 3 || throw(ArgumentError("ThreeBodyBackend requires three final particles"))
    ms = ThreeBodyDecays.ThreeBodyMasses((s["mass"] for s in finals)...; m0 = kin["initial_state"]["mass"])
    spins = ThreeBodyDecays.ThreeBodySpins((_spin_to_two_j(s["spin"]) for s in finals)...; two_h0 = _spin_to_two_j(kin["initial_state"]["spin"]))
    sys = ThreeBodyDecays.ThreeBodySystem(ms, spins)
    chains = map(dd["chains"]) do ch
        k = _threebody_k(ch["topology"]);i, j = ThreeBodyDecays.ij_from_k(k)
        p = only(ch["propagators"])
        vertex(node) = begin
            v = only(v for v in ch["vertices"] if _json_to_tuple(v["node"]) == node)
            ThreeBodyDecays.Vertex(_recoupling(v), _formfactor(v, workspace))
        end
        ThreeBodyDecays.DecayChain(;
            k, two_j = _spin_to_two_j(p["spin"]), Xlineshape = workspace[p["parametrization"]],
            Hij = vertex((i, j)), HRk = vertex(((i, j), k)), tbs = sys
        )
    end
    model = ThreeBodyDecays.ThreeBodyDecay(
        Tuple(chains), Tuple(_threebody_phase(_threebody_k(ch["topology"]), _threebody_k(dd["reference_topology"]), spins[4]) * _parse_complex(ch["weight"]) for ch in dd["chains"]),
        Tuple(ch["name"] for ch in dd["chains"]); reference_k = _threebody_k(dd["reference_topology"])
    )
    return LoadedModel(model, sys, context)
end
_system_masses(sys::ThreeBodyDecays.ThreeBodySystem) = SystemMasses(sys.ms)
function _model_sites(model::ThreeBodyDecays.ThreeBodyDecay)
    out = Pair{Site, Any}[]
    for (c, ch) in enumerate(model.chains)
        i, j = ThreeBodyDecays.ij_from_k(ch.k)
        push!(
            out, Site(c, (i, j), :propagator) => ch.Xlineshape,
            Site(c, (i, j), :formfactor) => ch.Hij.ff,
            Site(c, ((i, j), ch.k), :formfactor) => ch.HRk.ff
        )
    end
    return out
end
function _decay_document(model::ThreeBodyDecays.ThreeBodyDecay, system, context, emit)
    all(ch -> ch.tbs == system, model.chains) || throw(ArgumentError("system disagrees with chain kinematics"))
    spins = SystemSpins(Tuple(system.two_js)[1:3]...; two_h0 = system.two_js[4])
    kin, _ = serializeToDict(CascadeSystem(spins, _system_masses(system)); particle_labels = isempty(context.labels) ? nothing : context.labels)
    chains = Any[]
    for (c, ch) in enumerate(model.chains)
        i, j = ThreeBodyDecays.ij_from_k(ch.k)
        verts = Any[]
        for (node, v) in [((i, j), ch.Hij), (((i, j), ch.k), ch.HRk)]
            row, _ = serializeToDict(v.h)
            row["node"] = _json_node(node);row["formfactor"] = emit(Site(c, node, :formfactor), v.ff)
            push!(verts, row)
        end
        props = [Dict("spin" => _spin_string(ch.two_j), "node" => [i, j], "parametrization" => emit(Site(c, (i, j), :propagator), ch.Xlineshape))]
        push!(chains, Dict("name" => model.names[c], "weight" => _complex_string(_threebody_phase(ch.k, ThreeBodyDecays.reference_k(model), system.two_js[4]) * model.couplings[c]), "topology" => [[i, j], ch.k], "propagators" => props, "vertices" => verts))
    end
    k = ThreeBodyDecays.reference_k(model);i, j = ThreeBodyDecays.ij_from_k(k)
    return Dict("kinematics" => kin, "chains" => chains, "reference_topology" => [[i, j], k])
end
