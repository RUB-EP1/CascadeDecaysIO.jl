using Test
using CascadeDecays
using CascadeDecaysIO
using HadronicLineshapes
using JSON
using ThreeBodyDecays: RecouplingLS

@testset "CascadeDecaysIO writer" begin
    spins = SystemSpins(0, 0, 0, 0; two_h0 = 0)
    masses = SystemMasses(1.0, 1.1, 1.2, 1.3; m0 = 5.0)
    system = CascadeSystem(spins, masses)
    topology = DecayTopology((((1, 2), 3), 4))
    chain = DecayChain(
        topology,
        spins;
        propagators = (
            ((1, 2), 3) => Propagator(2, BreitWigner(2.4, 0.1)),
            (1, 2) => Propagator(0, ConstantLineshape(1.0 + 0.0im)),
        ),
        vertices = (
            (((1, 2), 3), 4) => Vertex(RecouplingLS((0, 2))),
            ((1, 2), 3) => Vertex(RecouplingLS((0, 2))),
            (1, 2) => Vertex(RecouplingLS((0, 0))),
        ),
    )

    kinematics, _ = serializeToDict(system; particle_labels = ("D0", "pi+", "D-", "K+", "B+"))
    @test kinematics["initial_state"]["index"] == 0
    @test length(kinematics["final_state"]) == 4
    @test kinematics["final_state"][1]["name"] == "D0"

    chain_dict, appendix = serializeToDict(chain; name = "test_chain")
    @test chain_dict["topology"] == [[[1, 2], 3], 4]
    @test length(chain_dict["propagators"]) == 2
    @test length(chain_dict["vertices"]) == 3
    @test all(haskey(v, "formfactor") for v in chain_dict["vertices"])
    @test all(haskey(p, "parametrization") for p in chain_dict["propagators"])
    @test length(appendix) >= 2

    decay_description, full_appendix = serializeToDict(
        system,
        ["test_chain" => (1.0 + 0.0im, chain)];
        particle_labels = ("D0", "pi+", "D-", "K+", "B+"),
    )
    @test decay_description["reference_topology"] == [[[1, 2], 3], 4]
    @test decay_description["chains"][1]["weight"] == "1.0 + 0.0i"
    @test length(full_appendix) >= length(appendix)

    # Test CascadeDecay container serialization
    cascade = CascadeDecay((chain,), topology; couplings = (1.0 + 0.0im,), names = ("test_chain",))
    cd_desc, cd_app = serializeToDict(cascade, masses; particle_labels = ("D0", "pi+", "D-", "K+", "B+"))
    @test cd_desc["reference_topology"] == [[[1, 2], 3], 4]
    @test cd_desc["chains"][1]["name"] == "test_chain"

    document = amplitudeSerializationDict(
        system,
        ["test_chain" => (1.0 + 0.0im, chain)];
        particle_labels = ("D0", "pi+", "D-", "K+", "B+"),
        name = "test_model",
        variables = ["m12sq", "m123sq"],
        domains = ["phase_space"],
        misc = Dict("generator" => "CascadeDecaysIO.jl test"),
        parameter_points = [Dict("name" => "nominal", "parameters" => [])],
        validation = Dict(
            :misc => Dict(
                :amplitude_model_checksums => [
                    Dict(:distribution => "test_model", :point => "validation_point", :value => 1.0),
                ],
            ),
            :parameter_points => [
                Dict(:name => "validation_point", :parameters => [Dict(:name => "m1_2", :value => 2.0)]),
            ],
        ),
    )
    @test haskey(document, "distributions")
    @test haskey(document, "functions")
    @test document["distributions"][1]["name"] == "test_model"
    @test document["distributions"][1]["variables"] == ["m12sq", "m123sq"]
    @test document["misc"]["generator"] == "CascadeDecaysIO.jl test"
    @test document["misc"]["amplitude_model_checksums"][1]["distribution"] == "test_model"
    @test length(document["parameter_points"]) == 2
    @test document["parameter_points"][1]["name"] == "nominal"
    @test document["parameter_points"][2]["name"] == "validation_point"
    @test all(haskey(fn, "name") for fn in document["functions"])
    @test document["domains"] == ["phase_space"]
    @test haskey(document, "parameter_points")

    setVariables!(document, ["m12sq"])
    setDomains!(document, ["updated_domain"])
    setMisc!(document, Dict("reviewed" => true))
    setParameterPoints!(document, Dict("alternate" => Dict("test_chain" => "0.5 + 0.0i")))
    setValidation!(document, Dict("amplitude_model_checksums" => [Dict("distribution" => "test_model", "value" => 2.0)]))
    setSection!(document, "tuple_section", (Dict(:name => "tuple_value"),))
    appendFunction!(document, "constant_one", Dict("type" => "ConstantLineshape", "value" => 1.0))
    setSection!(document, "custom_section", Dict("kept" => true))
    @test document["distributions"][1]["variables"] == ["m12sq"]
    @test document["domains"] == ["updated_domain"]
    @test document["misc"]["reviewed"] == true
    @test document["misc"]["amplitude_model_checksums"][1]["value"] == 2.0
    @test haskey(document["parameter_points"], "alternate")
    @test document["tuple_section"][1]["name"] == "tuple_value"
    @test document["functions"][end]["name"] == "constant_one"
    @test document["custom_section"]["kept"] == true

    # Paired topology test ((1, 2), (3, 4))
    dk_top = DecayTopology(((1, 2), (3, 4)))
    dk_chain = DecayChain(
        dk_top,
        spins;
        propagators = (
            (1, 2) => Propagator(0, ConstantLineshape(1.0 + 0.0im)),
            (3, 4) => Propagator(2, BreitWigner(2.4, 0.1)),
        ),
        vertices = (
            ((1, 2), (3, 4)) => Vertex(RecouplingLS((0, 2))),
            (3, 4) => Vertex(RecouplingLS((0, 2))),
            (1, 2) => Vertex(RecouplingLS((0, 0))),
        ),
    )
    dk_dict, _ = serializeToDict(dk_chain; name = "dk_test")
    @test dk_dict["topology"] == [[1, 2], [3, 4]]

    mktempdir() do dir
        output_path = joinpath(dir, "model.json")
        returned_path = writeJson(output_path, document)
        @test returned_path == output_path
        parsed = JSON.parsefile(output_path)
        @test parsed["distributions"][1]["name"] == "test_model"
        @test parsed["functions"][end]["name"] == "constant_one"

        # Test reader and round-trip
        read_cascade, read_system, read_doc = readJson(output_path)
        @test length(read_cascade.chains) == 1
        @test read_cascade.names[1] == "test_chain"
        @test read_cascade.couplings[1] ≈ 1.0 + 0.0im
        @test read_system.masses.m0 == 5.0
        @test read_system.masses.finals == [1.0, 1.1, 1.2, 1.3]
    end
end
