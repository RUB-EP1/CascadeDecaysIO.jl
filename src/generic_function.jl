const _GENERIC_OPERATORS = (:+, :-, :*, :/, :^)
const _GENERIC_FUNCTIONS = (
    :sqrt, :exp, :log, :log10, :sin, :cos, :tan, :asin, :acos, :atan,
    :sinh, :cosh, :tanh, :abs, :abs2, :conj, :real, :imag,
)
# functions whose real domain is restricted: evaluated with a complex argument
const _GENERIC_COMPLEX_DOMAIN = (:sqrt, :log, :log10, :asin, :acos)
# `i` (amplitude-serialization) and `im` (Julia) are the imaginary unit; `π` is
# accepted on reading, the writer only emits ASCII
const _GENERIC_CONSTANTS = (:i, :im, :pi, :π)
# limits checked before parsing: Meta.parse recurses on nesting and on chains of `^`
const _GENERIC_MAX_LENGTH = 10_000
const _GENERIC_MAX_DEPTH = 200
const _GENERIC_MAX_POWERS = 200

"""
    GenericFunction(expression)

Evaluator for the amplitude-serialization `generic_function` type: an analytic
expression in at most one free variable (by convention `sigma`, the invariant mass
squared of the propagator node), e.g.

```julia
GenericFunction("1 / (4.0 - sigma - i * 2.0 * 0.1)")
```

Grammar: numbers, `+ - * / ^`, parentheses, the functions
$(join(_GENERIC_FUNCTIONS, ", ")), the imaginary unit `i` (or `im`) and `pi`.
The expression is checked before it is compiled, so a JSON file cannot run
arbitrary code. Integer literals are read as floating-point numbers; `sqrt`, `log`,
`log10`, `asin`, `acos` and powers with non-integer exponents are evaluated in
complex arithmetic, so the function continues analytically below thresholds.
"""
struct GenericFunction{F} <: HadronicLineshapes.AbstractFlexFunc
    expression::String
    variable::Symbol
    f::F
end

function GenericFunction(expression::AbstractString)
    length(expression) <= _GENERIC_MAX_LENGTH ||
        throw(ArgumentError("generic_function longer than $(_GENERIC_MAX_LENGTH) characters"))
    _paren_depth(expression) <= _GENERIC_MAX_DEPTH ||
        throw(ArgumentError("generic_function nested deeper than $(_GENERIC_MAX_DEPTH) parentheses"))
    count(==('^'), expression) <= _GENERIC_MAX_POWERS ||
        throw(ArgumentError("generic_function with more than $(_GENERIC_MAX_POWERS) powers"))
    parsed = Meta.parse(expression; raise = true)
    free = Set{Symbol}()
    _check_generic_expression(parsed, free, expression)
    length(free) <= 1 || throw(ArgumentError(
        "generic_function \"$expression\" has more than one free variable: $(join(sort!(collect(free)), ", "))"))
    variable = isempty(free) ? :sigma : only(free)
    body = _to_julia(parsed)
    f = Core.eval(@__MODULE__, Expr(:->, variable, body))
    return GenericFunction(String(expression), variable, f)
end

(g::GenericFunction)(σ::Number) = ComplexF64(Base.invokelatest(g.f, complex(float(σ))))

function _paren_depth(expression)
    depth = max_depth = 0
    for c in expression
        c == '(' && (depth += 1; max_depth = max(max_depth, depth))
        c == ')' && (depth -= 1)
    end
    return max_depth
end

function _check_generic_expression(ex, free::Set{Symbol}, expression)
    ex isa Number && return nothing
    if ex isa Symbol
        ex in _GENERIC_CONSTANTS && return nothing
        (ex in _GENERIC_OPERATORS || ex in _GENERIC_FUNCTIONS) && throw(ArgumentError(
            "generic_function \"$expression\": `$ex` is a function, not a variable"))
        push!(free, ex)
        return nothing
    end
    if ex isa Expr && ex.head === :call
        fname = first(ex.args)
        fname isa Symbol && (fname in _GENERIC_OPERATORS || fname in _GENERIC_FUNCTIONS) ||
            throw(ArgumentError("generic_function \"$expression\": function `$fname` is not allowed"))
        foreach(arg -> _check_generic_expression(arg, free, expression), ex.args[2:end])
        return nothing
    end
    throw(ArgumentError("generic_function \"$expression\": unsupported syntax `$ex`"))
end

# translate the checked expression into Julia: `i` -> `im`, integer literals -> Float64,
# complex arguments where the real domain is restricted
_to_julia(ex) = ex === :i ? :im : ex
_to_julia(ex::Integer) = Float64(ex)
function _to_julia(ex::Expr)
    fname, args = ex.args[1], map(_to_julia, ex.args[2:end])
    if fname in _GENERIC_COMPLEX_DOMAIN
        args = [:(complex($a)) for a in args]
    elseif fname === :^ && !(ex.args[3] isa Integer)
        args = [:(complex($(args[1]))), args[2]]
    elseif fname === :^
        args = [args[1], ex.args[3]]   # keep integer exponents exact (x^2 == x*x)
    end
    return Expr(:call, fname, args...)
end
