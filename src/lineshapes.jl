import HadronicLineshapes: breakup, BlattWeisskopf, BW

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

function (bw::TFPWAMultichannelBreitWigner)(σ::Number)
    m0 = bw.m
    mΓ = sum(bw.channels) do channel
        gsq, ma, mb, l, d = channel.gsq, channel.ma, channel.mb, channel.l, channel.d
        FF = BlattWeisskopf{l}(d)
        p = breakup(sqrt(σ), ma, mb)
        gsq * σ * 2p / sqrt(σ) * FF(p)^2
    end
    BW(σ, m0, mΓ / m0)
end
(bw::TFPWAMultichannelBreitWigner)(σ::Real) = bw(σ + 1im * eps())

"""
    NRExpLineshape

Nonresonant exponential lineshape: `-exp(-αβ * (σ - m0²))`.
"""
struct NRExpLineshape <: HadronicLineshapes.AbstractFlexFunc
    αβ::ComplexF64
    m0::Float64
end

(ls::NRExpLineshape)(σ::Number) = -exp(-ls.αβ * (σ - ls.m0^2))
(ls::NRExpLineshape)(σ::Real) = ls(σ + 1im * eps())

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

(nl::NamedLineshape)(args...) = nl.lineshape(args...)
Base.getproperty(nl::NamedLineshape, s::Symbol) =
    (s === :name || s === :lineshape) ? getfield(nl, s) : getproperty(getfield(nl, :lineshape), s)

