_backend(::CascadeDecay) = CascadeBackend()

function _build_model(doc, ::CascadeBackend, workspace, context)
    return LoadedModel(dict2instance(CascadeDecay, doc; workspace), dict2instance(CascadeSystem, doc), context)
end

function _model_sites(model::CascadeDecay)
    out = Pair{Site, Any}[]
    for (c, ch) in enumerate(model.chains)
        for (i, line) in enumerate(propagating_line_inds(ch))
            push!(out, Site(c, _line_address(ch.topology, line), :propagator) => ch.propagators[i])
        end
        for i in 1:nvertices(ch)
            push!(out, Site(c, _line_address(ch.topology, incoming_line_ind(ch.topology, i)), :formfactor) => ch.vertices[i].ff)
        end
    end
    return out
end
function _decay_document(model::CascadeDecay, system, context, emit)
    spins, _ = _cascade_spins_and_chains(model)
    sys = system isa SystemMasses ? CascadeSystem(spins, system) : system
    final_two_js(sys) == final_two_js(spins) && root_two_j(sys) == root_two_j(spins) || throw(ArgumentError("system spins disagree with model"))
    kin, _ = serializeToDict(sys; particle_labels = isempty(context.labels) ? nothing : context.labels)
    chains = Any[]
    for (c, ch) in enumerate(model.chains)
        props = Any[]; verts = Any[]
        for (i, line) in enumerate(propagating_line_inds(ch))
            node = _line_address(ch.topology, line)
            push!(props, Dict("node" => _json_node(node), "spin" => _spin_string(propagator_two_js(ch)[i]), "parametrization" => emit(Site(c, node, :propagator), ch.propagators[i])))
        end
        for i in 1:nvertices(ch)
            node = _line_address(ch.topology, incoming_line_ind(ch.topology, i))
            v, _ = serializeToDict(ch.vertices[i].h)
            v["node"] = _json_node(node);v["formfactor"] = emit(Site(c, node, :formfactor), ch.vertices[i].ff)
            push!(verts, v)
        end
        push!(chains, Dict("name" => model.names[c], "weight" => _complex_string(model.couplings[c]), "topology" => _json_topology(ch.topology), "propagators" => props, "vertices" => verts))
    end
    return Dict("kinematics" => kin, "chains" => chains, "reference_topology" => _json_topology(model.reference_topology))
end

_system_masses(system::CascadeSystem) = system.masses
_system_masses(system::SystemMasses) = system
