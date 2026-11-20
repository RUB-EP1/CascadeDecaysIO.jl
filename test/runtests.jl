using Test, CascadeDecaysIO, CascadeDecays, HadronicLineshapes, Random, JSON
import ThreeBodyDecays as TBD
const CDIO = CascadeDecaysIO

include("fixtures.jl")

@testset "CascadeDecaysIO" begin
    include("runtests/expressions.jl")
    include("runtests/backends.jl")
    include("runtests/parameters.jl")
    include("runtests/reference_models.jl")
    include("runtests/reference_frames.jl")
    include("runtests/four_body.jl")
    include("runtests/errors.jl")
    include("runtests/schema.jl")
end
