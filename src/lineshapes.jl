import HadronicLineshapes: breakup, BlattWeisskopf, BW

@inline function _blatt_weisskopf_sq(l::Int, d::Float64, p)
    l == 0 && return true
    l == 1 && return BlattWeisskopf{1}(d)(p)^2
    l == 2 && return BlattWeisskopf{2}(d)(p)^2
    l == 3 && return BlattWeisskopf{3}(d)(p)^2
    l == 4 && return BlattWeisskopf{4}(d)(p)^2
    return BlattWeisskopf{l}(d)(p)^2
end

"""
    TFPWAMultichannelBreitWigner

TFPWA-aligned multichannel Breit–Wigner lineshape where the running width carries an extra factor of `σ`.
"""
struct TFPWAMultichannelBreitWigner{N} <: HadronicLineshapes.AbstractFlexFunc
    m::Float64
    channels::SVector{N,<:NamedTuple{(:gsq, :ma, :mb, :l, :d)}}
end

function TFPWAMultichannelBreitWigner(
    m::Real,
    channels::Vector{<:NamedTuple{(:gsq, :ma, :mb, :l, :d)}},
)
    N = length(channels)
    return TFPWAMultichannelBreitWigner(Float64(m), SVector{N}(channels...))
end

@inline function (bw::TFPWAMultichannelBreitWigner)(σ::Number)
    m0 = bw.m
    sqrt_σ = sqrt(σ)
    mΓ = zero(typeof(σ * 1.0))
    @inbounds for channel in bw.channels
        p = breakup(sqrt_σ, channel.ma, channel.mb)
        ff_sq = _blatt_weisskopf_sq(channel.l, channel.d, p)
        mΓ += channel.gsq * σ * 2p / sqrt_σ * ff_sq
    end
    return BW(σ, m0, mΓ / m0)
end
@inline (bw::TFPWAMultichannelBreitWigner)(σ::Real) = bw(σ + 1im * eps())

"""
    NRExpLineshape

Nonresonant exponential lineshape: `-exp(-αβ * (σ - m0²))`.
"""
struct NRExpLineshape <: HadronicLineshapes.AbstractFlexFunc
    αβ::ComplexF64
    m0::Float64
end

@inline (ls::NRExpLineshape)(σ::Number) = -exp(-ls.αβ * (σ - ls.m0 * ls.m0))
@inline (ls::NRExpLineshape)(σ::Real) = ls(σ + 1im * eps())

"""
    NamedLineshape

Wrapper associating an explicit name with an underlying lineshape.
When serialized by `CascadeDecaysIO`, the lineshape will use this name in
the `parametrization` field of propagators and in the `functions` list.
"""
struct NamedLineshape{F} <: HadronicLineshapes.AbstractFlexFunc
    name::String
    lineshape::F
end

@inline (nl::NamedLineshape)(σ) = nl.lineshape(σ)
@inline (nl::NamedLineshape)(args...) = nl.lineshape(args...)
@inline Base.getproperty(nl::NamedLineshape, s::Symbol) =
    (s === :name || s === :lineshape) ? getfield(nl, s) : getproperty(getfield(nl, :lineshape), s)
