function _validate_auxiliary(doc)
    length(doc["domains"]) == 1 || throw(ArgumentError("one phase-space domain required"))
    dom = only(doc["domains"])
    _check_fields(dom, ["name", "type", "m0", "final_state_masses"])
    dom["type"] == "phase_space" || throw(ArgumentError("unsupported domain type"))
    all(x -> _finite(x)&&x >= 0, [dom["m0"];dom["final_state_masses"]]) || throw(ArgumentError("invalid domain masses"))
    point_names = Set{String}(); vector_points = Set{String}()
    for p in doc["parameter_points"]
        p["name"] in point_names && throw(ArgumentError("duplicate point $(p["name"])"))
        push!(point_names, p["name"])
        if haskey(p, "parameters")
            _check_fields(p, ["name", "parameters"])
            names = Set{String}()
            for v in p["parameters"]
                _check_fields(v, ["name", "value"])
                _identifier(v["name"]) && _finite(v["value"]) || throw(ArgumentError("invalid parameter point"))
                v["name"] in names && throw(ArgumentError("duplicate parameter"))
                push!(names, v["name"])
            end
        else
            _check_fields(p, ["name", "domain", "four_vectors"])
            p["domain"] == dom["name"] || throw(ArgumentError("unknown point domain"))
            push!(vector_points, p["name"])
            ids = Int[]
            for v in p["four_vectors"]
                _check_fields(v, ["index", "E", "px", "py", "pz"])
                all(k -> _finite(v[k]), ["E", "px", "py", "pz"]) || throw(ArgumentError("nonfinite four-vector"))
                push!(ids, v["index"])
            end
            sort(ids) == collect(1:length(dom["final_state_masses"])) || throw(ArgumentError("invalid four-vector indices"))
        end
    end
    _check_fields(doc["misc"], String[], ["amplitude_model_checksums"])
    for c in get(doc["misc"], "amplitude_model_checksums", [])
        _check_fields(c, ["point", "distribution", "value"])
        c["point"] in vector_points || throw(ArgumentError("checksum needs a four-vector point"))
        c["distribution"] in [d["name"] for d in doc["distributions"]] || throw(ArgumentError("unknown checksum distribution"))
        _finite(c["value"]) && c["value"] >= 0 || throw(ArgumentError("invalid checksum"))
    end
    return
end
