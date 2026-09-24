# ==============================================================================
# Minimal Dummy Cascade Decay Analysis & JSON Round-Trip Example
# ==============================================================================
# This standalone script demonstrates the complete workflow of:
#   1. Building a minimal 1 -> 4 cascade amplitude model mirroring the B2DxDK
#      topology and spin/lineshape structure using generic, decay-independent
#      physics parameters.
#   2. Evaluating complex amplitudes and intensities on sampled 4-vector events.
#   3. Serializing the full analysis (including custom lineshapes and 4-vector
#      verification points) to a ComPWA Amplitude Serialization JSON file.
#   4. Restarting the amplitude calculation purely from the JSON file via
#      `readJson` and verifying bit-exact agreement.
# ==============================================================================

using CascadeDecays
using CascadeDecays.ThreeBodyDecays: RecouplingLS
using CascadeDecaysIO
using FourVectors
using HadronicLineshapes
using JSON
using Test

# ------------------------------------------------------------------------------
# 1. External System & Topologies
# ------------------------------------------------------------------------------
# Generic 1 -> 4 decay: P0(J=0) -> p1(J=0) p2(J=0) p3(J=0) p4(J=0)
# (Spin values in SystemSpins are given as 2*J)
external_spins = SystemSpins(0, 0, 0, 0; two_h0 = 0)
particle_labels = ("p1", "p2", "p3", "p4", "P0")

# Two competing topologies (mirroring B2DxDK's sequential DxD and paired DK trees):
#   - seq_topology:    P0 -> R_123 + p4,  R_123 -> V_12 + p3,  V_12 -> p1 + p2
#   - paired_topology: P0 -> V_12 + R_34, V_12  -> p1 + p2,    R_34 -> p3 + p4
seq_topology    = DecayTopology((((1, 2), 3), 4))
paired_topology = DecayTopology(((1, 2), (3, 4)))
kinematic_task  = KinematicTask((seq_topology, paired_topology))

# ------------------------------------------------------------------------------
# 2. Lineshapes & Form Factors (Generic / Random Physics Parameters)
# ------------------------------------------------------------------------------
const R_BW = 3.0 # Blatt-Weisskopf barrier radius [GeV^-1]

# Shared intermediate vector state V_12 on node (1, 2) with constant propagator
v12_lineshape = NamedLineshape("constant_V12", ConstantLineshape(1.0 + 0.0im))

# (a) Physical multi-wave Breit-Wigner (S- and D-wave partial widths) on node ((1, 2), 3)
ls_resA = NamedLineshape(
    "lineshape_ResA",
    MultichannelBreitWigner(
        3.50,
        [
            (; gsq = 0.45, ma = 1.80, mb = 1.50, l = 0, d = R_BW),
            (; gsq = 0.20, ma = 1.80, mb = 1.50, l = 2, d = R_BW),
        ],
    ),
)

# (b) Standard single-channel relativistic Breit-Wigner on node ((1, 2), 3)
ls_resB = NamedLineshape(
    "lineshape_ResB",
    BreitWigner(3.65, 0.08, 1.80, 1.50, 1, R_BW),
)

# (c) Non-resonant exponential background on node ((1, 2), 3) -> serialized as "type": "custom"
ls_nr_exp = NamedLineshape(
    "lineshape_NR_Exp",
    NRExpLineshape(0.15 + 0.25im, 3.80),
)

# (d) Paired-channel Breit-Wigner on node (3, 4)
ls_resC = NamedLineshape(
    "lineshape_ResC",
    BreitWigner(2.40, 0.06, 1.50, 0.50, 1, R_BW),
)

# ------------------------------------------------------------------------------
# 3. Build Decay Chains & Combine into a CascadeDecay Model
# ------------------------------------------------------------------------------
# Helper for sequential chains: (((1, 2), 3), 4)
function make_seq_chain(res_lineshape, res_two_j, prod_two_ls, dec_two_ls)
    l_prod = div(prod_two_ls[1], 2)
    l_dec  = div(dec_two_ls[1], 2)
    return DecayChain(
        seq_topology,
        external_spins;
        propagators = (
            (1, 2)       => Propagator(2, v12_lineshape), # J = 1 (two_j = 2)
            ((1, 2), 3)  => Propagator(res_two_j, res_lineshape),
        ),
        vertices = (
            (((1, 2), 3), 4) => Vertex(RecouplingLS(prod_two_ls), BlattWeisskopf{l_prod}(R_BW)),
            ((1, 2), 3)      => Vertex(RecouplingLS(dec_two_ls),  BlattWeisskopf{l_dec}(R_BW)),
            (1, 2)           => Vertex(RecouplingLS((2, 0))),     # P-wave V_12 -> p1 p2
        ),
    )
end

# Helper for paired chains: ((1, 2), (3, 4))
function make_paired_chain(res_lineshape, res_two_j, root_two_ls, sub_two_ls)
    l_root = div(root_two_ls[1], 2)
    l_sub  = div(sub_two_ls[1], 2)
    return DecayChain(
        paired_topology,
        external_spins;
        propagators = (
            (1, 2) => Propagator(2, v12_lineshape),
            (3, 4) => Propagator(res_two_j, res_lineshape),
        ),
        vertices = (
            ((1, 2), (3, 4)) => Vertex(RecouplingLS(root_two_ls), BlattWeisskopf{l_root}(R_BW)),
            (3, 4)           => Vertex(RecouplingLS(sub_two_ls),  BlattWeisskopf{l_sub}(R_BW)),
            (1, 2)           => Vertex(RecouplingLS((2, 0))),
        ),
    )
end

# Construct 4 representative chains (S/D-wave multichannel, P-wave BW, NR exp, and paired sub-channel)
chains = (
    make_seq_chain(ls_resA,   2, (2, 2), (0, 2)), # J=1 ResA (S-wave decay)
    make_seq_chain(ls_resA,   2, (2, 2), (4, 2)), # J=1 ResA (D-wave decay)
    make_seq_chain(ls_resB,   0, (0, 0), (2, 2)), # J=0 ResB (P-wave decay)
    make_seq_chain(ls_nr_exp, 0, (0, 0), (2, 2)), # J=0 Non-resonant exponential
    make_paired_chain(ls_resC, 2, (2, 2), (2, 0)), # J=1 ResC in (3, 4) sub-channel
)

chain_names = (
    "ResA_L1_d0",
    "ResA_L1_d2",
    "ResB_L0_d1",
    "NR_Exp_L0_d1",
    "ResC_L1_d1",
)

# Random / generic complex couplings for each decay chain
chain_couplings = (
     0.85 - 0.15im,
    -0.30 + 0.40im,
     0.50 + 0.25im,
    -0.20 + 0.10im,
     0.35 - 0.45im,
)

model = CascadeDecay(
    chains,
    seq_topology;
    couplings = chain_couplings,
    names     = chain_names,
)

# ------------------------------------------------------------------------------
# 4. Generate 3 Sampled Four-Vector Events & Evaluate Amplitudes
# ------------------------------------------------------------------------------
# Synthetic 4-momenta (px, py, pz; E) for 3 phase-space events in the P0 rest frame
raw_events = [
    (
        p1 = FourVector( 0.42, -0.31,  0.18; E = 1.52),
        p2 = FourVector(-0.08,  0.14, -0.05; E = 0.35),
        p3 = FourVector(-0.55,  0.29, -0.41; E = 1.65),
        p4 = FourVector( 0.21, -0.12,  0.28; E = 0.68),
    ),
    (
        p1 = FourVector(-0.38,  0.45, -0.22; E = 1.56),
        p2 = FourVector( 0.12, -0.09,  0.11; E = 0.34),
        p3 = FourVector( 0.48, -0.52,  0.33; E = 1.67),
        p4 = FourVector(-0.22,  0.16, -0.22; E = 0.63),
    ),
    (
        p1 = FourVector( 0.19,  0.51, -0.36; E = 1.58),
        p2 = FourVector(-0.15, -0.11,  0.09; E = 0.36),
        p3 = FourVector(-0.28, -0.60,  0.44; E = 1.69),
        p4 = FourVector( 0.24,  0.20, -0.17; E = 0.62),
    ),
]

println("=== Step 1: Evaluating In-Memory Dummy Model ===")
original_amplitudes  = ComplexF64[]
original_intensities = Float64[]
parameter_points     = Dict{String,Any}[]
checksums            = Dict{String,Any}[]

for (i, ev) in enumerate(raw_events)
    pt  = KinematicPoint(kinematic_task, (ev.p1, ev.p2, ev.p3, ev.p4))
    amp = only(amplitude(model, pt))
    val = abs2(amp)

    push!(original_amplitudes, amp)
    push!(original_intensities, val)

    pt_name = "verification_point_$i"
    println("  $pt_name: A = $amp, |A|^2 = $val")

    # Record 4-vectors for the JSON parameter_points section
    push!(parameter_points, Dict(
        "name" => pt_name,
        "four_vectors" => Dict(
            "p1" => Dict("E" => ev.p1.E, "px" => ev.p1.px, "py" => ev.p1.py, "pz" => ev.p1.pz),
            "p2" => Dict("E" => ev.p2.E, "px" => ev.p2.px, "py" => ev.p2.py, "pz" => ev.p2.pz),
            "p3" => Dict("E" => ev.p3.E, "px" => ev.p3.px, "py" => ev.p3.py, "pz" => ev.p3.pz),
            "p4" => Dict("E" => ev.p4.E, "px" => ev.p4.px, "py" => ev.p4.py, "pz" => ev.p4.pz),
        ),
    ))

    # Record reference amplitude & intensity checksums
    push!(checksums, Dict(
        "distribution"   => "dummy_cascade_intensity",
        "point"          => pt_name,
        "value"          => val,
        "amplitude_real" => real(amp),
        "amplitude_imag" => imag(amp),
    ))
end

# ------------------------------------------------------------------------------
# 5. Write Model & Verification Points to JSON
# ------------------------------------------------------------------------------
json_path = joinpath(@__DIR__, "dummy_analysis_model.json")

doc = amplitudeSerializationDict(
    model;
    particle_labels,
    name             = "dummy_cascade_intensity",
    variables        = ["m1_2sq", "m1_2_3sq", "m3_4sq"],
    parameter_points = parameter_points,
    misc             = Dict(
        "description"               => "Minimal dummy analysis testing B2DxDK structure",
        "generator"                 => "CascadeDecaysIO.jl",
        "amplitude_model_checksums" => checksums,
    ),
)

writeJson(json_path, doc)
println("\n=== Step 2: Wrote Serialized JSON Model ===")
println("  File: $json_path ($(filesize(json_path)) bytes)")

# ------------------------------------------------------------------------------
# 6. Restart Computation Using ONLY the JSON File via Reader
# ------------------------------------------------------------------------------
println("\n=== Step 3: Restarting Calculation Strictly from JSON ===")

# Read model and JSON document from disk (independent of in-memory objects above)
reloaded_model, reloaded_system, reloaded_doc = readJson(json_path)

# Reconstruct kinematic task from topologies present in the reloaded model
reloaded_topologies = unique([ch.topology for ch in reloaded_model.chains])
reloaded_task = KinematicTask(Tuple(reloaded_topologies))

# Extract verification 4-vectors and expected checksums from the JSON document
json_points    = reloaded_doc["parameter_points"]
json_checksums = reloaded_doc["misc"]["amplitude_model_checksums"]

for (i, pt_dict) in enumerate(json_points)
    fv = pt_dict["four_vectors"]
    q1 = FourVector(fv["p1"]["px"], fv["p1"]["py"], fv["p1"]["pz"]; E = fv["p1"]["E"])
    q2 = FourVector(fv["p2"]["px"], fv["p2"]["py"], fv["p2"]["pz"]; E = fv["p2"]["E"])
    q3 = FourVector(fv["p3"]["px"], fv["p3"]["py"], fv["p3"]["pz"]; E = fv["p3"]["E"])
    q4 = FourVector(fv["p4"]["px"], fv["p4"]["py"], fv["p4"]["pz"]; E = fv["p4"]["E"])

    pt_reloaded  = KinematicPoint(reloaded_task, (q1, q2, q3, q4))
    amp_reloaded = only(amplitude(reloaded_model, pt_reloaded))
    val_reloaded = abs2(amp_reloaded)

    # Compare against expected JSON checksums and original evaluation
    expected_amp = ComplexF64(json_checksums[i]["amplitude_real"], json_checksums[i]["amplitude_imag"])
    expected_val = Float64(json_checksums[i]["value"])

    delta_amp = abs(amp_reloaded - expected_amp)
    delta_val = abs(val_reloaded - expected_val)

    println("  $(pt_dict["name"]): Reloaded A = $amp_reloaded, |A|^2 = $val_reloaded  (Δ|A|^2 = $delta_val, ΔA = $delta_amp)")

    @test amp_reloaded ≈ expected_amp atol = 1e-14
    @test val_reloaded ≈ expected_val atol = 1e-14
    @test amp_reloaded == original_amplitudes[i]
end

println("Done.")
