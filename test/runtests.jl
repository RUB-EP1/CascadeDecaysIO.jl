using Test
using CascadeDecays
using CascadeDecaysIO
using HadronicLineshapes
using JSON
using OrderedCollections
using Random: Xoshiro
using CascadeDecays.ThreeBodyDecays: RecouplingLS

struct UnknownShape end

const MASSES = SystemMasses(1.0, 0.2, 1.2, 0.5; m0 = 5.0)
const LABELS = ("a", "b", "c", "d", "X")
const SPINS = SystemSpins(0, 0, 0, 0; two_h0 = 0)
const SEQ = DecayTopology((((1, 2), 3), 4))
const PAIRED = DecayTopology(((1, 2), (3, 4)))

# on-shell events with m12 fixed, like a D* in the narrow-width approximation
const EVENTS = [phaseSpaceEvent(Xoshiro(k), MASSES, SEQ; node_masses = Dict((1, 2) => 1.25)) for k in 1:3]

seq_chain(lineshape, two_j, prod_two_ls, dec_two_ls; d12 = ConstantLineshape(1.0 + 0.0im)) = DecayChain(
    SEQ, SPINS;
    propagators = ((1, 2) => Propagator(2, d12), ((1, 2), 3) => Propagator(two_j, lineshape)),
    vertices = (
        (((1, 2), 3), 4) => Vertex(RecouplingLS(prod_two_ls), BlattWeisskopf{div(prod_two_ls[1], 2)}(3.0)),
        ((1, 2), 3) => Vertex(RecouplingLS(dec_two_ls), BlattWeisskopf{div(dec_two_ls[1], 2)}(3.0)),
        (1, 2) => Vertex(RecouplingLS((2, 0))),
    ),
)

paired_chain(lineshape) = DecayChain(
    PAIRED, SPINS;
    propagators = ((1, 2) => Propagator(2, ConstantLineshape(1.0 + 0.0im)), (3, 4) => Propagator(2, lineshape)),
    vertices = (
        ((1, 2), (3, 4)) => Vertex(RecouplingLS((2, 2)), BlattWeisskopf{1}(3.0)),
        (3, 4) => Vertex(RecouplingLS((2, 0)), MomentumPower{1}()),
        (1, 2) => Vertex(RecouplingLS((2, 0))),
    ),
)

function test_model()
    mcbw = NamedLineshape("lineshape_R1", MultichannelBreitWigner(3.1, [
        (; gsq = 0.4, ma = 1.25, mb = 1.2, l = 0, d = 3.0),
        (; gsq = 0.2, ma = 1.25, mb = 1.2, l = 2, d = 3.0)]))
    bw = NamedLineshape("lineshape_R2", BreitWigner(3.3, 0.1, 1.25, 1.2, 1, 3.0))
    expo = NamedLineshape("lineshape_NR", GenericFunction("exp(-(0.1 + 0.2 * i) * (sigma - 9.0))"))
    flatte = NamedLineshape("lineshape_R3", Flatte(1.9, 0.3, 1.2, 0.5, 0.1, 0.9, 0.5))
    chains = (
        seq_chain(mcbw, 2, (2, 2), (0, 2)),
        seq_chain(mcbw, 2, (2, 2), (4, 2)),
        seq_chain(bw, 0, (0, 0), (2, 2)),
        seq_chain(expo, 0, (0, 0), (2, 2)),
        paired_chain(flatte),
    )
    return CascadeDecay(chains, SEQ;
        couplings = (0.8 - 0.1im, -0.3 + 0.4im, 0.5, -0.2 + 0.1im, 0.35 - 0.45im),
        names = ("R1_S", "R1_D", "R2", "NR", "R3"))
end

function write_and_read(document; custom_functions = Dict{String,Any}(), ascii_only = true)
    mktempdir() do dir
        path = joinpath(dir, "model.json")
        writeJson(path, document; ascii_only)
        return readJson(path; custom_functions)..., read(path, String)
    end
end

@testset "CascadeDecaysIO" begin
    model = test_model()
    validation = validationSection(model, EVENTS; distribution = "test_model")
    document = amplitudeSerializationDict(model, MASSES; particle_labels = LABELS, name = "test_model", validation)
    dd = document["distributions"][1]["decay_description"]
    functions = Dict(f["name"] => f for f in document["functions"])

    @testset "document structure" begin
        @test collect(keys(document)) == ["distributions", "functions", "domains", "misc", "parameter_points"]
        @test !haskey(document["distributions"][1], "variables")
        @test collect(keys(document["misc"])) == ["amplitude_model_checksums"]
        @test document["domains"] == [phaseSpaceDomain(MASSES)]
        @test document["domains"][1]["type"] == "phase_space"
        domain = document["domains"][1]
        @test domain["initial_state"] == LittleDict("index" => 0, "mass" => 5.0)
        @test [(f["index"], f["mass"]) for f in domain["final_state"]] == [(1, 1.0), (2, 0.2), (3, 1.2), (4, 0.5)]
    end

    @testset "kinematics: string spins, masses only in the domain" begin
        kin = dd["kinematics"]
        @test kin["initial_state"] == LittleDict("index" => 0, "name" => "X", "spin" => "0")
        @test !any(haskey(f, "mass") for f in kin["final_state"])
        @test all(f["spin"] isa String for f in kin["final_state"])
        @test_throws ArgumentError serializeToDict(SPINS)
        @test_throws ArgumentError amplitudeSerializationDict(model)
    end

    @testset "chains: string spins, l and s; weights" begin
        chain = dd["chains"][2]
        @test chain["name"] == "R1_D"
        @test chain["topology"] == [[[1, 2], 3], 4]
        @test all(p["spin"] isa String for p in chain["propagators"])
        decay_vertex = only(v for v in chain["vertices"] if v["node"] == [[1, 2], 3])
        @test decay_vertex["type"] == "ls"
        @test decay_vertex["l"] == "2" && decay_vertex["s"] == "1"
        @test decay_vertex["formfactor"] == "BlattWeisskopf_l2"
        @test only(v for v in chain["vertices"] if v["node"] == [1, 2])["formfactor"] == ""
        pattern = r"^-?\d+(\.\d+)?(e[-+]?\d+)? *[+-] *\d+(\.\d+)?(e[-+]?\d+)?i$"
        @test all(occursin(pattern, c["weight"]) for c in dd["chains"])
        @test dd["chains"][3]["weight"] == "0.5 + 0.0i"
    end

    @testset "functions: standard types, names, x" begin
        @test functions["BlattWeisskopf_l1"] == LittleDict("type" => "BlattWeisskopf", "l" => 1, "radius" => 3.0, "name" => "BlattWeisskopf_l1")
        @test functions["MomentumPower_l1"]["type"] == "MomentumPower"
        @test functions["lineshape_R1"]["type"] == "MultichannelBreitWigner"
        @test functions["lineshape_R1"]["x"] == "m_12_3_sq"
        @test functions["lineshape_R2"]["type"] == "BreitWigner"
        @test functions["lineshape_R3"]["type"] == "MultichannelBreitWigner"   # Flatte
        @test functions["lineshape_R3"]["x"] == "m_34_sq"
        @test functions["lineshape_NR"] == LittleDict("type" => "generic_function",
            "expression" => "exp(-(0.1 + 0.2 * i) * (sigma - 9.0))", "name" => "lineshape_NR")
        constant = only(f for f in document["functions"] if f["type"] == "Polynomial")
        @test constant["coefficients"] == [1.0] && constant["x"] == "m_12_sq"
        @test !any(haskey(f, "subtype") for f in document["functions"])
    end

    @testset "writer rejects what it cannot express" begin
        bad = CascadeDecay((seq_chain(UnknownShape(), 0, (0, 0), (2, 2)),), SEQ; couplings = (1.0,), names = ("bad",))
        @test_throws ArgumentError amplitudeSerializationDict(bad, MASSES)
        phased = CascadeDecay((seq_chain(ConstantLineshape(1.0im), 0, (0, 0), (2, 2)),), SEQ; couplings = (1.0,), names = ("phase",))
        @test_throws ArgumentError amplitudeSerializationDict(phased, MASSES)
        clash = CascadeDecay(
            (seq_chain(NamedLineshape("same", BreitWigner(3.3, 0.1, 1.25, 1.2, 1, 3.0)), 0, (0, 0), (2, 2)),
             seq_chain(NamedLineshape("same", BreitWigner(3.5, 0.1, 1.25, 1.2, 1, 3.0)), 0, (0, 0), (2, 2))),
            SEQ; couplings = (1.0, 1.0), names = ("a", "b"))
        @test_throws ArgumentError amplitudeSerializationDict(clash, MASSES)
        mktempdir() do dir
            unicode = amplitudeSerializationDict(model, MASSES; particle_labels = ("a", "b", "c", "d", "Λb"))
            @test_throws ArgumentError writeJson(joinpath(dir, "u.json"), unicode)
        end
    end

    @testset "GenericFunction" begin
        g = GenericFunction("1 / (4.0 - sigma - i * 2.0 * 0.1)")
        @test g(3.0) ≈ 1 / (4.0 - 3.0 - 0.2im)
        @test GenericFunction("1 / (4.0 - σ - 1im * 0.2)")(3.0) ≈ g(3.0)          # Lc reference model style
        @test GenericFunction("sqrt(sigma - 2.0)")(1.0) ≈ 1.0im                    # continued below threshold
        @test GenericFunction("2.5")(1.0) == 2.5 + 0.0im
        @test_throws ArgumentError GenericFunction("run(`ls`)")
        @test_throws ArgumentError GenericFunction("sigma + x")
        @test_throws ArgumentError GenericFunction("Base.exit(0)")
        @test_throws ArgumentError GenericFunction("(a = 1; a)")
        @test GenericFunction("10^19")(1.0) == 1.0e19 + 0.0im                       # no integer overflow
        @test GenericFunction("sqrt(-4.0) + 0 * sigma")(1.0) ≈ 2.0im                # constant, continued
        @test GenericFunction("sigma^0.5")(-4.0) ≈ 2.0im
        @test_throws ArgumentError GenericFunction("("^500 * "sigma" * ")"^500)    # depth limit before parsing
        @test_throws ArgumentError GenericFunction(join(fill("sigma", 3000), "^"))  # power chain limit before parsing
    end

    read_model, system, read_document, text = write_and_read(document)
    task = kinematicTask(model)

    @testset "round trip" begin
        @test isascii(text)
        @test system.masses == MASSES
        @test read_model.names == model.names
        @test all(read_model.couplings .≈ model.couplings)
        @test any(p -> p isa NamedLineshape{<:GenericFunction}, read_model.chains[4].propagators)
        @test all(p -> p isa NamedLineshape, read_model.chains[1].propagators)    # names kept
        for event in EVENTS
            point = KinematicPoint(task, event)
            @test amplitude(read_model, point) ≈ amplitude(model, point) rtol = 1e-13
        end
    end

    @testset "validation" begin
        @test length(read_document["parameter_points"]) == 3
        point = read_document["parameter_points"][1]
        @test collect(keys(point)) == ["name", "domain", "final_state"]
        @test [v["index"] for v in point["final_state"]] == [1, 2, 3, 4]
        @test collect(keys(point["final_state"][1]["p"])) == ["E", "px", "py", "pz"]
        @test [c["point"] for c in read_document["misc"]["amplitude_model_checksums"]] ==
              ["validation_point1", "validation_point2", "validation_point3"]
        results = validate(read_document, read_model)
        @test all(r -> r.passed, results)
        @test maximum(r -> r.mass2_offset, results) < 1e-12

        tampered = deepcopy(read_document)
        tampered["misc"]["amplitude_model_checksums"][2]["value"] *= 1.001
        @test !validate(tampered, read_model)[2].passed

        off_shell = deepcopy(read_document)
        off_shell["parameter_points"][1]["final_state"][1]["p"]["E"] += 1e-3
        @test !validate(off_shell, read_model)[1].passed

        reversed = deepcopy(read_document)                 # particles are matched by index
        reverse!(reversed["parameter_points"][1]["final_state"])
        @test validate(reversed, read_model)[1].passed

        missing_particle = deepcopy(read_document)         # the point must list the particles of its domain
        pop!(missing_particle["parameter_points"][1]["final_state"])
        @test_throws ArgumentError validate(missing_particle, read_model)

        lineshape_target = deepcopy(read_document)         # every checksum is a distribution value
        lineshape_target["misc"]["amplitude_model_checksums"][1]["distribution"] = "lineshape_R2"
        @test_throws ArgumentError validate(lineshape_target, read_model)
    end

    @testset "reader is strict" begin
        unknown = deepcopy(read_document)
        for f in unknown["functions"]
            f["name"] == "lineshape_NR" && (f["type"] = "SomeAnalysisShape")
        end
        @test_throws ArgumentError CascadeDecaysIO.dict2instance(CascadeDecay, unknown)
        bound = CascadeDecaysIO.dict2instance(CascadeDecay, unknown;
            custom_functions = Dict("lineshape_NR" => σ -> exp(-(0.1 + 0.2im) * (σ - 9.0))))
        point = KinematicPoint(task, EVENTS[1])
        @test amplitude(bound, point) ≈ amplitude(model, point) rtol = 1e-13

        missing_ff = deepcopy(read_document)
        missing_ff["distributions"][1]["decay_description"]["chains"][1]["vertices"][1]["formfactor"] = "nope"
        @test_throws ArgumentError CascadeDecaysIO.dict2instance(CascadeDecay, missing_ff)

        massless = deepcopy(read_document)
        delete!(massless["domains"][1]["initial_state"], "mass")
        @test_throws ArgumentError CascadeDecaysIO.dict2instance(CascadeSystem, massless)
        @test_throws ArgumentError CascadeDecaysIO.dict2instance(CascadeSystem, LittleDict("domains" => [],
            "distributions" => read_document["distributions"]))
    end

    @testset "read -> write -> read" begin
        rewritten = amplitudeSerializationDict(read_model, system.masses; particle_labels = LABELS, name = "test_model")
        @test sort([f["name"] for f in rewritten["functions"]]) == sort([f["name"] for f in read_document["functions"]])
        again, _, _, _ = write_and_read(rewritten)
        for event in EVENTS
            point = KinematicPoint(task, event)
            @test amplitude(again, point) ≈ amplitude(model, point) rtol = 1e-13
        end
    end

    @testset "strict parsing" begin
        @test CascadeDecaysIO._spin_to_two_j("-3/2") == -3
        @test CascadeDecaysIO._spin_to_two_j(0.5) == 1
        for bad in ("1/3", "0.5", "2/4", 0.3, "x")
            @test_throws ArgumentError CascadeDecaysIO._spin_to_two_j(bad)
        end
        polarized = deepcopy(read_document)
        polarized["distributions"][1]["type"] = "HadronicPolarizedIntensity"
        @test_throws ArgumentError CascadeDecaysIO.dict2instance(CascadeDecay, polarized)
        unknown_target = deepcopy(read_document)
        unknown_target["misc"]["amplitude_model_checksums"][1]["distribution"] = "no_such_distribution"
        @test_throws ArgumentError validate(unknown_target, read_model)
    end

    @testset "polynomial lineshape round trip" begin
        poly = NamedLineshape("lineshape_poly", CascadeDecaysIO.Polynomials.Polynomial([1.0, -0.2, 0.03]))
        m = CascadeDecay((seq_chain(poly, 0, (0, 0), (2, 2)),), SEQ; couplings = (1.0,), names = ("poly",))
        doc = amplitudeSerializationDict(m, MASSES)
        fn = only(f for f in doc["functions"] if f["name"] == "lineshape_poly")
        @test fn["type"] == "Polynomial" && fn["coefficients"] == [1.0, -0.2, 0.03] && fn["x"] == "m_12_3_sq"
        m2, _, _, _ = write_and_read(doc)
        point = KinematicPoint(kinematicTask(m), EVENTS[2])
        @test amplitude(m2, point) ≈ amplitude(m, point) rtol = 1e-14
    end

    @testset "generated validation points" begin
        on_domain(event) = all(values(CascadeDecaysIO._phase_space_offsets(fourVectorPoint("p", event), phaseSpaceDomain(MASSES))) .< 1e-12)
        mass(p) = sqrt(CascadeDecaysIO._mass_squared(p))
        rng = Xoshiro(7)
        for topology in (SEQ, PAIRED, (((1, 2), 3), 4), ((3, 1), (4, 2)))
            @test all(on_domain(phaseSpaceEvent(rng, MASSES, topology)) for _ in 1:20)
        end
        fixed = [phaseSpaceEvent(rng, MASSES, SEQ; node_masses = Dict((1, 2) => 1.25, ((1, 2), 3) => 3.05))
                 for _ in 1:20]
        @test all(on_domain, fixed)
        @test all(e -> abs(mass(e[1] + e[2]) - 1.25) < 1e-12, fixed)
        @test all(e -> abs(mass(e[1] + e[2] + e[3]) - 3.05) < 1e-12, fixed)
        @test_throws ArgumentError phaseSpaceEvent(rng, MASSES, SEQ; node_masses = Dict((3, 4) => 2.0))
        @test_throws ArgumentError phaseSpaceEvent(rng, MASSES, SEQ; node_masses = Dict((1, 2) => 1.1))
        @test_throws ArgumentError phaseSpaceEvent(rng, MASSES, SEQ; node_masses = Dict(((1, 2), 3) => (3.0, 3.1)))

        events = validationEvents(model, MASSES)
        @test length(events) == 3 && all(on_domain, events)
        @test events == validationEvents(model, MASSES)                        # reproducible default seed
        mean = sum(unpolarized_intensity(model, KinematicPoint(task, phaseSpaceEvent(rng, MASSES, SEQ))) for _ in 1:1000) / 1000
        @test all(unpolarized_intensity(model, KinematicPoint(task, e)) > 0.05 * mean for e in events)

        # a particle with fixed mass: the model vanishes off m12 = 1.25
        delta = NamedLineshape("lineshape_delta", GenericFunction("exp(-0.5 * ((sqrt(sigma) - 1.25) / 1.0e-6)^16)"))
        narrow = CascadeDecay((seq_chain(NamedLineshape("lineshape_R2", BreitWigner(3.3, 0.1, 1.25, 1.2, 1, 3.0)),
            0, (0, 0), (2, 2); d12 = delta),), SEQ; couplings = (1.0,), names = ("R2",))
        @test_throws ArgumentError validationEvents(narrow, MASSES)
        @test_throws ArgumentError amplitudeSerializationDict(narrow, MASSES)
        on_shell = validationSection(narrow, MASSES; distribution = "narrow", node_masses = Dict((1, 2) => 1.25))
        @test all(c -> c["value"] > 0, on_shell.checksums)
        off_shell = phaseSpaceEvent(rng, MASSES, SEQ; node_masses = Dict((1, 2) => 1.26))
        @test_throws ArgumentError validationSection(narrow, [off_shell]; distribution = "narrow")

        # written by default when the document is created, and they validate
        doc = amplitudeSerializationDict(narrow, MASSES; name = "narrow",
            validation = validationSection(narrow, MASSES; distribution = "narrow", node_masses = Dict((1, 2) => 1.25)))
        @test length(doc["parameter_points"]) == 3 && all(r -> r.passed, validate(doc))
        generated = amplitudeSerializationDict(model, MASSES; name = "generated")
        @test length(generated["parameter_points"]) == 3
        @test all(r -> r.passed, validate(generated))
        @test isempty(amplitudeSerializationDict(model, MASSES; validation = nothing)["parameter_points"])
    end

    @testset "propagator names" begin
        chain_dict, appendix = serializeToDict(paired_chain(BreitWigner(2.4, 0.1));
            name = "override", propagator_names = ["constant_D", "lineshape_my"])
        @test Set(p["parametrization"] for p in chain_dict["propagators"]) == Set(["constant_D", "lineshape_my"])
        for p in chain_dict["propagators"]
            @test appendix[p["parametrization"]]["x"] == (p["node"] == [3, 4] ? "m_34_sq" : "m_12_sq")
        end
    end
end

# The reference models of the standard: half-integer spins, helicity and parity
# vertices, Flatte, polynomials, generic functions, several distributions.
@testset "reference models of amplitude-serialization" begin
    for file in ("lc2ppik-lhcb-2683025.json", "lb2pkg-lhcb-2765817.json", "x2pipipi-compass-1391643.json")
        path = joinpath(@__DIR__, "reference_models", file)
        @testset "$file" begin
            results = validate(path)
            @test !isempty(results)
            # CascadeDecays misses this point (all momenta in the x-z plane, pi exactly along -x)
            # by 17%; rotating the event by 1e-9 rad reproduces the checksum
            known = file == "lc2ppik-lhcb-2683025.json" ? ("validation_point",) : ()
            for r in results
                r.point in known ? (@test_broken r.passed) : (@test r.passed)
            end
            # write the first distribution again, read it back, same intensity at its points
            model, system, document = readJson(path)
            labels = (map(f -> f["name"], sort(document["distributions"][1]["decay_description"]["kinematics"]["final_state"]; by = f -> f["index"]))...,
                      document["distributions"][1]["decay_description"]["kinematics"]["initial_state"]["name"])
            rewritten = amplitudeSerializationDict(model, system.masses; particle_labels = labels)
            again, _, _, text = write_and_read(rewritten)        # σ (lc2ppik) is written as sigma
            @test isascii(text)
            @test length(rewritten["parameter_points"]) == 3     # generated validation points
            @test all(r -> r.passed, validate(rewritten))
            task = kinematicTask(model)
            name = document["distributions"][1]["name"]
            for checksum in document["misc"]["amplitude_model_checksums"]
                checksum["distribution"] == name || continue
                point = only(p for p in document["parameter_points"] if p["name"] == checksum["point"])
                kp = KinematicPoint(task, CascadeDecaysIO._four_vectors(point))
                @test unpolarized_intensity(again, kp) ≈ unpolarized_intensity(model, kp) rtol = 1e-13
            end
        end
    end
end
