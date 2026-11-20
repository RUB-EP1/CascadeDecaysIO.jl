module CascadeDecaysIO

using CascadeDecays
using FourVectors
using HadronicLineshapes
using JSON
using OrderedCollections
import Polynomials
using Random: AbstractRNG, Xoshiro
using StaticArrays
using ThreeBodyDecays

import CascadeDecays:
    CascadeDecay,
    ConstantLineshape,
    DecayChain,
    DecayTopology,
    KinematicPoint,
    KinematicTask,
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
    root_two_j,
    unpolarized_intensity

export SystemMasses, CascadeSystem
export serializeToDict, dict2instance, phaseSpaceDomain, fourVectorPoint,
    kinematicTask, phaseSpaceEvent, validationEvents, validationSection, validate,
    ExpressionFunction, Site, ModelIOContext, LoadedModel, CascadeBackend,
    ThreeBodyBackend, read_document, build_model, write_document, document,
    validate_document, remap_context, set_parameters

include("types.jl")
include("parsing.jl")
include("expression.jl")
include("codecs.jl")
include("writer.jl")
include("reader.jl")
include("parameters.jl")
include("context.jl")
include("checks.jl")
include("schema.jl")
include("phase_space.jl")
include("validation.jl")
include("cascade.jl")
include("threebody.jl")
include("document.jl")

end
