function _json_to_tuple(node)
    if node isa AbstractVector
        return Tuple(_json_to_tuple(x) for x in node)
    elseif node isa Integer
        return Int(node)
    else
        return node
    end
end

# Spins, helicities and s: strings "n" or "n/2" as in the standard (pattern of the
# schema's Spin); numbers are accepted when they are integer or half-integer
const _SPIN_PATTERN = r"^(?![-+]0)([-+]?\d+)(/2)?$"

function _spin_to_two_j(s)
    s isa Integer && return 2 * Int(s)
    if s isa Real
        isinteger(2s) || throw(ArgumentError("spin $s is not a multiple of 1/2"))
        return Int(2s)
    end
    s isa AbstractString || throw(ArgumentError("could not parse spin from $(repr(s))"))
    m = match(_SPIN_PATTERN, strip(s))
    m === nothing && throw(ArgumentError("spin $(repr(s)) is not of the form \"n\" or \"n/2\""))
    n = parse(Int, m.captures[1])
    return m.captures[2] === nothing ? 2n : n
end

# orbital angular momenta: non-negative integers, as number or string
function _parse_int(x)
    value = x isa Integer ? Int(x) : x isa AbstractString ? tryparse(Int, strip(x)) : nothing
    (value === nothing || value < 0) && throw(ArgumentError("expected a non-negative integer, got $(repr(x))"))
    return value
end

function _parse_complex(s)
    s isa Number && return ComplexF64(s)
    str = strip(string(s))
    if endswith(str, 'i')
        body = @view str[1:prevind(str, lastindex(str))]
        idx = lastindex(body)
        while idx > firstindex(body)
            c = body[idx]
            if (c == '+' || c == '-') && !(body[prevind(body, idx)] in ('e', 'E'))
                re_val = tryparse(Float64, strip(@view body[1:prevind(body, idx)]))
                im_val = tryparse(Float64, strip(@view body[nextind(body, idx):end]))
                if re_val !== nothing && im_val !== nothing
                    return ComplexF64(re_val, c == '-' ? -im_val : im_val)
                end
                break
            end
            idx = prevind(body, idx)
        end
        im_only = tryparse(Float64, strip(body))
        im_only !== nothing && return ComplexF64(0.0, im_only)
    else
        re_only = tryparse(Float64, str)
        re_only !== nothing && return ComplexF64(re_only, 0.0)
    end
    throw(ArgumentError("could not parse complex number $(repr(str))"))
end
