"Read and validate a callable-format document; no backend is constructed yet."
function read_document(path::AbstractString; parameter_point = nothing)
    doc = JSON.parsefile(path; dicttype = Dict{String, Any})
    return validate_document(doc; parameter_point)
end

"""
    build_model(document, backend=CascadeBackend(); distribution=1, parameter_point=nothing)

Build one distribution with plain physics callables. The result separates the model,
physical system and export context. Select a parameter point explicitly when the
input contains several numerical parameter snapshots.
"""
function build_model(doc::AbstractDict, backend = CascadeBackend(); distribution = 1, parameter_point = nothing)
    validate_document(doc; parameter_point)
    selected = _select_distribution(doc, distribution)
    parameters = _parameter_values(doc, parameter_point)
    workspace = Dict(
        fn["name"] => _function_instance(_resolved_function(fn, parameters))
            for fn in doc["functions"]
    )
    context = _read_context(doc, selected; parameter_point)
    selected_document = merge(doc, Dict("distributions" => [selected]))
    return _build_model(selected_document, backend, workspace, context)
end

function _parameter_points(parameters, name)
    isempty(parameters) && return Any[]
    entries = [Dict("name" => key, "value" => parameters[key]) for key in sort(collect(keys(parameters)))]
    return [Dict("name" => name, "parameters" => entries)]
end

"""
    document(model; system, context=ModelIOContext(), validation=false)
    document(loaded; validation=false)

Export the current numerical model, restoring explicit names and parameter bindings
from its context. Structural changes require remapping or a fresh context. Old
checksums are omitted; `validation=true` generates fresh checks for Cascade models.
"""
function document(model; system, context = ModelIOContext(), validation = false)
    functions = Dict{String, Any}()
    parameters = Dict{String, Any}()
    reserved = Set(values(context.names))
    emit = (site, object) -> _emit!(functions, parameters, context, site, object, reserved)
    decay = _decay_document(model, system, context, emit)
    _check_context(context, model, decay)
    distribution = Dict(
        "name" => context.distribution,
        "type" => "HadronicUnpolarizedIntensity",
        "decay_description" => decay,
    )
    output = Dict{String, Any}(
        "format_version" => FORMAT_VERSION,
        "distributions" => [distribution],
        "functions" => [functions[key] for key in sort(collect(keys(functions)))],
        "domains" => [phaseSpaceDomain(_system_masses(system); name = context.domain)],
        "parameter_points" => _parameter_points(parameters, context.parameter_point),
        "misc" => Dict{String, Any}(),
    )
    if validation
        model isa CascadeDecay || throw(ArgumentError("automatic validation generation currently requires CascadeBackend"))
        checks = validationSection(model, _system_masses(system); distribution = context.distribution, domain = context.domain)
        setValidation!(output, checks)
    end
    return validate_document(output)
end

function document(loaded::LoadedModel; validation = false)
    return document(loaded.model; system = loaded.system, context = loaded.context, validation)
end

function _canonical(value)
    value isa AbstractDict && return OrderedDict(key => _canonical(value[key]) for key in sort(collect(keys(value))))
    value isa AbstractVector && return _canonical.(value)
    return value
end

"Validate and serialize completely before opening the destination for writing."
function write_document(path, value; kwargs...)
    doc = value isa AbstractDict ? validate_document(value; kwargs...) : document(value; kwargs...)
    text = sprint(io -> JSON.print(io, _canonical(doc), 4)) * "\n"
    open(path, "w") do io
        write(io, text)
    end
    return path
end
