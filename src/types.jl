"""
    SystemMasses(m1, m2, ...; m0)

Masses of the external particles in GeV. Final-state masses are positional; the
mass of the decaying particle is given as keyword `m0`. The amplitude-serialization
standard requires them in `kinematics` and in the `phase_space` domain.
"""
struct SystemMasses{Nf, T <: Real}
    finals::SVector{Nf, T}
    m0::T
end

function SystemMasses(ms...; m0)
    mass_tuple = Tuple(ms)
    length(mass_tuple) >= 1 ||
        throw(ArgumentError("provide at least one final-state mass before `m0`"))
    T = promote_type(typeof(m0), map(typeof, mass_tuple)...)
    return SystemMasses{length(mass_tuple), T}(
        SVector{length(mass_tuple), T}(mass_tuple),
        convert(T, m0),
    )
end

function SystemMasses(ms::ThreeBodyDecays.MassTuple)
    return SystemMasses(ms.m1, ms.m2, ms.m3; m0 = ms.m0)
end

Base.:(==)(a::SystemMasses, b::SystemMasses) =
    a.m0 == b.m0 && a.finals == b.finals

"""
    CascadeSystem(spins, masses)

External information (spins and masses) for cascade decay model kinematics serialization.
`spins` can be a `SystemSpins` or `SystemSpinParities` from `CascadeDecays`.
`masses` is a `SystemMasses`.
"""
struct CascadeSystem{Nf, Tm, Q <: Union{SystemSpins{Nf}, SystemSpinParities{Nf}}}
    quantum::Q
    masses::SystemMasses{Nf, Tm}
end

function CascadeSystem(spins::SystemSpins{Nf}, masses::SystemMasses{Nf, Tm}) where {Nf, Tm}
    return CascadeSystem{Nf, Tm, typeof(spins)}(spins, masses)
end

function CascadeSystem(quantum::SystemSpinParities{Nf}, masses::SystemMasses{Nf, Tm}) where {Nf, Tm}
    return CascadeSystem{Nf, Tm, typeof(quantum)}(quantum, masses)
end

final_two_js(system::CascadeSystem) = final_two_js(system.quantum)
root_two_j(system::CascadeSystem) = root_two_j(system.quantum)
final_masses(system::CascadeSystem) = system.masses.finals
root_mass(system::CascadeSystem) = system.masses.m0

const FORMAT_VERSION = "0.1"
struct CascadeBackend end
struct ThreeBodyBackend end

struct Site
    chain::Int
    node::Tuple
    role::Symbol
    function Site(chain, node, role)
        chain > 0 || throw(ArgumentError("chain index must be positive"))
        role in (:propagator, :formfactor) || throw(ArgumentError("invalid site role"))
        return new(chain, _json_to_tuple(node), role)
    end
end
Base.:(==)(a::Site, b::Site) = a.chain == b.chain && a.node == b.node && a.role == b.role
Base.isequal(a::Site, b::Site) = a == b
Base.hash(s::Site, h::UInt) = hash((s.chain, s.node, s.role), h)

Base.@kwdef mutable struct ModelIOContext
    names::Dict{Site, String} = Dict{Site, String}()
    # Original definitions preserve parameter-reference paths, never object identities.
    definitions::Dict{String, Any} = Dict{String, Any}()
    structure::Union{Nothing, Vector{Tuple{String, Tuple}}} = nothing
    labels::Vector{String} = String[]
    distribution::String = "amplitude_model"
    parameter_point::String = "nominal"
    domain::String = "default"
end
struct LoadedModel{M, S}
    model::M
    system::S
    context::ModelIOContext
end
