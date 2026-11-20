# A deliberately small expression language. Parsing never evaluates Julia code.
const _OPS = Dict{Symbol, Any}(
    :+ => +, :- => -, :* => *, :/ => /, :^ => ^,
    :exp => exp, :log => log, :sqrt => sqrt, :sin => sin, :cos => cos,
    :tan => tan, :abs => abs, :abs2 => abs2, :conj => conj, :real => real,
    :imag => imag, :complex => complex
)
const _CONSTANTS = Dict(:PI => Float64(pi), :EULER => exp(1.0))
_identifier(s) = s isa AbstractString && occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", s)
_finite(x) = x isa Real && !(x isa Bool) && isfinite(x)

"""
    ExpressionFunction(expression, arguments; parameters=Dict())

Pure callable with ordered local arguments and explicitly bound numerical parameters.
Supports arithmetic, complex(re, im), PI/EULER and selected mathematical functions.
Complex log/sqrt use principal branches. No Julia eval, assignment, lookup or I/O.
"""
struct ExpressionFunction <: HadronicLineshapes.AbstractFlexFunc
    expression::String
    arguments::Tuple{Vararg{String}}
    parameters::Dict{String, Any}
    ast::Any
end

function ExpressionFunction(expression::AbstractString, arguments; parameters = Dict())
    0 < length(expression) <= 10_000 || throw(ArgumentError("invalid expression length"))
    occursin(r"^[A-Za-z0-9_ .+*/^(),\s-]+$", expression) || throw(ArgumentError("unsupported expression characters"))
    depth = 0
    for c in expression
        c == '(' && (depth += 1)
        c == ')' && (depth -= 1)
        0 <= depth <= 100 || throw(ArgumentError("expression nesting exceeds limit"))
    end
    count(==('^'), expression) <= 100 || throw(ArgumentError("too many powers"))
    args = Tuple(String.(arguments))
    !isempty(args) && length(unique(args)) == length(args) || throw(ArgumentError("arguments must be nonempty and unique"))
    pars = Dict{String, Any}(String(k) => v for (k, v) in parameters)
    names = [args...; collect(keys(pars))]
    all(_identifier, names) || throw(ArgumentError("invalid local identifier"))
    length(unique(names)) == length(names) || throw(ArgumentError("arguments and parameters overlap"))
    any(n -> haskey(_OPS, Symbol(n)) || haskey(_CONSTANTS, Symbol(n)), names) && throw(ArgumentError("reserved local identifier"))
    all(_finite, values(pars)) || throw(ArgumentError("expression parameters must be finite real numbers"))
    ast = try
        Meta.parse(expression)
    catch
        throw(ArgumentError("invalid expression syntax"))
    end
    _check_ast(ast, Set(Symbol.(names)), 0)
    return ExpressionFunction(String(expression), args, pars, ast)
end

function _check_ast(ex, names, depth)
    depth <= 100 || throw(ArgumentError("expression tree exceeds depth limit"))
    _finite(ex) && return
    ex isa Symbol && (ex in names || haskey(_CONSTANTS, ex)) && return
    if ex isa Expr && ex.head == :call && first(ex.args) isa Symbol && haskey(_OPS, first(ex.args))
        op = first(ex.args); n = length(ex.args) - 1
        valid = op == :+ ? n >= 1 : op == :* ? n >= 2 : op == :- ? n in (1, 2) : op in (:/, :^, :complex) ? n == 2 : n == 1
        valid || throw(ArgumentError("wrong arity for $op"))
        foreach(a -> _check_ast(a, names, depth + 1), ex.args[2:end])
        return
    end
    throw(ArgumentError("undeclared symbol or unsupported expression syntax: $ex"))
end

_eval_ast(ex::Number, env) = float(ex)
_eval_ast(ex::Symbol, env) = haskey(_CONSTANTS, ex) ? _CONSTANTS[ex] : env[ex]
function _eval_ast(ex::Expr, env)
    op = first(ex.args)
    args = Any[_eval_ast(a, env) for a in ex.args[2:end]]
    if op in (:sqrt, :log)
        args[1] = complex(args[1])
    elseif op == :^ && !(args[2] isa Real && isinteger(args[2]))
        args[1] = complex(args[1])
    end
    return _OPS[op](args...)
end
function (f::ExpressionFunction)(xs::Number...)
    length(xs) == length(f.arguments) || throw(ArgumentError("expected $(length(f.arguments)) arguments, got $(length(xs))"))
    env = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in f.parameters)
    for (k, v) in zip(f.arguments, xs)
        env[Symbol(k)] = v
    end
    return _eval_ast(f.ast, env)
end
