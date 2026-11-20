function _select_point(points, selection)
    matches = [p for p in points if p["name"] == selection]
    length(matches) == 1 || throw(ArgumentError("unknown or duplicate parameter point $selection"))
    return only(matches)
end

function _parameter_values(doc, selection)
    points = [p for p in doc["parameter_points"] if haskey(p, "parameters")]
    isempty(points) && return Dict{String, Any}()
    selected = if selection === nothing
        length(points) == 1 || throw(ArgumentError("select a parameter point explicitly"))
        only(points)
    else
        _select_point(points, selection)
    end
    result = Dict{String, Any}()
    for p in selected["parameters"]
        _check_fields(p, ["name", "value"])
        haskey(result, p["name"]) && throw(ArgumentError("duplicate parameter $(p["name"])"))
        _identifier(p["name"]) && _finite(p["value"]) || throw(ArgumentError("invalid parameter"))
        result[p["name"]] = p["value"]
    end
    return result
end

function _resolve_value(v, values)
    if v isa AbstractString
        haskey(values, v) || throw(ArgumentError("unresolved parameter $v"))
        return values[v]
    elseif v isa AbstractDict
        return Dict{String, Any}(k => _resolve_value(x, values) for (k, x) in v)
    elseif v isa AbstractVector
        return [_resolve_value(x, values) for x in v]
    end
    _finite(v) || throw(ArgumentError("expected finite numeric parameter or reference"))
    return v
end
function _resolved_function(fn, values)
    type = fn["type"]
    haskey(_FUNCTION_FIELDS, type) || throw(ArgumentError("unsupported function type $type"))
    _check_fields(fn, vcat(["name", "type"], _FUNCTION_FIELDS[type]))
    out = Dict{String, Any}(fn)
    for k in _FUNCTION_FIELDS[type]
        k in ("expression", "arguments") && continue
        out[k] = _resolve_value(fn[k], values)
    end
    for k in ("mass", "width", "ma", "mb", "d", "radius")
        haskey(out, k) && out[k] < 0 && throw(ArgumentError("negative $k in $(fn["name"])"))
    end
    if haskey(out, "channels")
        isempty(out["channels"]) && throw(ArgumentError("empty channels"))
        for channel in out["channels"]
            _check_fields(channel, ["gsq", "ma", "mb", "l", "d"])
            all(v -> v >= 0, Base.values(channel)) || throw(ArgumentError("negative channel parameter"))
        end
    end
    haskey(out, "coefficients") && isempty(out["coefficients"]) && throw(ArgumentError("empty coefficients"))
    return out
end

function set_parameters(loaded::LoadedModel, updates::AbstractDict)
    doc = document(loaded)
    isempty(doc["parameter_points"]) && throw(ArgumentError("no named parameters"))
    entries = only(doc["parameter_points"])["parameters"]
    issubset(keys(updates), Set(p["name"] for p in entries)) || throw(ArgumentError("unknown parameter update"))
    for p in entries
        haskey(updates, p["name"]) && (p["value"] = updates[p["name"]])
    end
    return build_model(doc, _backend(loaded.model))
end

# Parameter references are restored only at their original, explicitly bound paths.
function _restore_refs(current, template, values, path = "")
    if template isa AbstractString
        _finite(current) || throw(ArgumentError("parameter $template no longer numeric at $path"))
        haskey(values, template) && values[template] != current && throw(ArgumentError("shared parameter $template diverged at $path"))
        values[template] = current
        return template
    elseif template isa AbstractDict
        current isa AbstractDict && Set(keys(current)) == Set(keys(template)) || throw(ArgumentError("parameter structure changed at $path; use a fresh context"))
        return Dict(k => _restore_refs(current[k], v, values, "$path.$k") for (k, v) in template)
    elseif template isa AbstractVector
        length(current) == length(template) || throw(ArgumentError("parameter array changed at $path"))
        return [_restore_refs(c, t, values, "$path[$i]") for (i, (c, t)) in enumerate(zip(current, template))]
    end
    return current
end
