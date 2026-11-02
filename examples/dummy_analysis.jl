# ==============================================================================
# Minimal cascade-decay analysis: write to amplitude-serialization JSON and
# restart the amplitude calculation from the file alone
# ==============================================================================
#   1. Build a 1 -> 4 cascade model with the B2DxDK topology and spin structure,
#      using generic parameters.
#   2. Generate three on-shell validation events (final-state four-vectors).
#   3. Write the model with four-vector validation points and intensity checksums.
#   4. Read the JSON back, recompute the checksums and compare amplitudes.
#
# Run from the package root:  julia --project=. examples/dummy_analysis.jl
# ==============================================================================

using CascadeDecays
using CascadeDecays.ThreeBodyDecays: RecouplingLS
using CascadeDecaysIO
using HadronicLineshapes
using Test

# ------------------------------------------------------------------------------
# 1. External system and topologies
# ------------------------------------------------------------------------------
# P0(J=0) -> p1 p2 p3 p4, all spin 0 (SystemSpins takes 2J)
external_spins = SystemSpins(0, 0, 0, 0; two_h0 = 0)
masses = SystemMasses(1.5, 0.3, 1.5, 0.5; m0 = 5.5)
particle_labels = ("p1", "p2", "p3", "p4", "P0")

#   seq:    P0 -> R123 p4,  R123 -> V12 p3,  V12 -> p1 p2
#   paired: P0 -> V12 R34,  V12 -> p1 p2,    R34 -> p3 p4
seq_topology = DecayTopology((((1, 2), 3), 4))
paired_topology = DecayTopology(((1, 2), (3, 4)))

# ------------------------------------------------------------------------------
# 2. Lineshapes and form factors: HadronicLineshapes types or GenericFunction
# ------------------------------------------------------------------------------
const R_BW = 3.0  # Blatt-Weisskopf radius [GeV^-1]
const M_V12 = 1.85

v12 = NamedLineshape("constant_V12", ConstantLineshape(1.0 + 0.0im))
res_a = NamedLineshape("lineshape_ResA", MultichannelBreitWigner(3.50, [
    (; gsq = 0.45, ma = M_V12, mb = 1.50, l = 0, d = R_BW),
    (; gsq = 0.20, ma = M_V12, mb = 1.50, l = 2, d = R_BW),
]))
res_b = NamedLineshape("lineshape_ResB", BreitWigner(3.65, 0.08, M_V12, 1.50, 1, R_BW))
# non-standard shape: written as "generic_function" (one free variable, ASCII)
nr_exp = NamedLineshape("lineshape_NR_Exp", GenericFunction("exp(-(0.15 + 0.25 * i) * (sigma - 14.44))"))
res_c = NamedLineshape("lineshape_ResC", BreitWigner(2.40, 0.06, 1.50, 0.50, 1, R_BW))

# ------------------------------------------------------------------------------
# 3. Chains and model
# ------------------------------------------------------------------------------
function seq_chain(lineshape, two_j, prod_two_ls, dec_two_ls)
    return DecayChain(seq_topology, external_spins;
        propagators = ((1, 2) => Propagator(2, v12), ((1, 2), 3) => Propagator(two_j, lineshape)),
        vertices = (
            (((1, 2), 3), 4) => Vertex(RecouplingLS(prod_two_ls), BlattWeisskopf{div(prod_two_ls[1], 2)}(R_BW)),
            ((1, 2), 3) => Vertex(RecouplingLS(dec_two_ls), BlattWeisskopf{div(dec_two_ls[1], 2)}(R_BW)),
            (1, 2) => Vertex(RecouplingLS((2, 0))),
        ))
end

function paired_chain(lineshape, two_j, root_two_ls, sub_two_ls)
    return DecayChain(paired_topology, external_spins;
        propagators = ((1, 2) => Propagator(2, v12), (3, 4) => Propagator(two_j, lineshape)),
        vertices = (
            ((1, 2), (3, 4)) => Vertex(RecouplingLS(root_two_ls), BlattWeisskopf{div(root_two_ls[1], 2)}(R_BW)),
            (3, 4) => Vertex(RecouplingLS(sub_two_ls), BlattWeisskopf{div(sub_two_ls[1], 2)}(R_BW)),
            (1, 2) => Vertex(RecouplingLS((2, 0))),
        ))
end

model = CascadeDecay(
    (
        seq_chain(res_a, 2, (2, 2), (0, 2)),
        seq_chain(res_a, 2, (2, 2), (4, 2)),
        seq_chain(res_b, 0, (0, 0), (2, 2)),
        seq_chain(nr_exp, 0, (0, 0), (2, 2)),
        paired_chain(res_c, 2, (2, 2), (2, 0)),
    ),
    seq_topology;
    couplings = (0.85 - 0.15im, -0.30 + 0.40im, 0.50 + 0.25im, -0.20 + 0.10im, 0.35 - 0.45im),
    names = ("ResA_l1_l0", "ResA_l1_l2", "ResB_l0_l1", "NR_Exp_l0_l1", "ResC_l1_l1"),
)

# ------------------------------------------------------------------------------
# 4. Three validation events: on-shell final-state four-vectors in the rest frame
#    of P0, V12 at its mass (like a D* in the narrow-width approximation), each with
#    an intensity of at least 0.1 x the mean intensity
# ------------------------------------------------------------------------------
events = validationEvents(model, masses; node_masses = Dict((1, 2) => M_V12))

# ------------------------------------------------------------------------------
# 5. Write: model + validation (four-vector points and intensity checksums)
# ------------------------------------------------------------------------------
json_path = joinpath(@__DIR__, "dummy_analysis_model.json")
validation = validationSection(model, events; distribution = "dummy_cascade_intensity")
document = amplitudeSerializationDict(model, masses;
    particle_labels, name = "dummy_cascade_intensity", validation)
writeJson(json_path, document)
println("Wrote $json_path ($(filesize(json_path)) bytes)")

# ------------------------------------------------------------------------------
# 6. Restart from the JSON file alone
# ------------------------------------------------------------------------------
reloaded, system, reloaded_document = readJson(json_path)
@test system.masses == masses

for result in validate(reloaded_document, reloaded)
    println("  $(result.point): |A|^2 = $(result.computed) (stored $(real(result.expected)))  passed = $(result.passed)")
    @test result.passed
end

task = kinematicTask(model)
for ev in events
    point = KinematicPoint(task, ev)
    @test amplitude(reloaded, point) ≈ amplitude(model, point) rtol = 1e-13
end
println("Done.")
