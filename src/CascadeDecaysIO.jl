module CascadeDecaysIO

using CascadeDecays
using HadronicLineshapes
using JSON
using OrderedCollections
using StaticArrays
using ThreeBodyDecays

import CascadeDecays:
    CascadeDecay,
    ConstantLineshape,
    DecayChain,
    DecayTopology,
    SystemSpinParities,
    SystemSpins,
    SystemSpinsOrSpinParities,
    child_line_inds,
    consumed_by,
    final_line_inds,
    final_two_js,
    incoming_line_ind,
    isfinal_line_ind,
    line_two_js,
    nfinal,
    nvertices,
    propagating_line_inds,
    propagator_two_js,
    root_line_ind,
    root_two_j

export SystemMasses,
    CascadeSystem

"""
    SystemMasses(m1, m2, ...; m0)

External system mass descriptor for amplitude serialization kinematics.
Final-state masses are positional; the root mass is always supplied as keyword `m0`.
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

export appendFunction!,
    amplitudeSerializationDict,
    serializeToDict,
    setDomains!,
    setDistributions!,
    setFunctions!,
    setMisc!,
    setParameterPoints!,
    setSection!,
    setValidation!,
    setVariables!,
    writeJson,
    dict2instance,
    readJson,
    TFPWAMultichannelBreitWigner,
    NRExpLineshape

include("lineshapes.jl")
include("writer.jl")
include("reader.jl")

end
