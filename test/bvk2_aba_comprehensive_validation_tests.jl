using Test
using CSV
using DataFrames
using Statistics

@testset "Comprehensive ABA publication packet contract" begin
    source = read(joinpath(@__DIR__, "..", "scripts",
        "write_bvk2_aba_comprehensive_validation.jl"), String)
    @test occursin("const ABA_CV_C2 = 0.16", source)
    @test occursin("bvk2_c2_calibration_nx128", source)
    @test occursin("selected_model.toml", source)
    @test occursin("aba_unified_benchmark.svg", source)
    @test occursin("aba_cell_stress_mechanism.svg", source)
    @test occursin("aba_ldis_crossing.svg", source)
    @test occursin("scft_reference_summary.csv", source)
    @test occursin("f0.50_chiN20", source)
    @test occursin("f0.50_chiN45", source)
    @test occursin("aba_cv_scft_fixed_cell_scan", source)
    @test occursin("Polyorder.gradient_wrt_cell", source)
    @test occursin("absolute D/Rg", source)
    @test occursin("SCFT stress-free D", source)
    @test !occursin("D / D0(model)", source)
    @test occursin("continuous_RPA_bifurcation", source)
    @test occursin("18.25, 18.5, 19.0, 20.0, 21.0, 22.0", source)
    @test occursin("DFM_SKIP_BVK2_ABA_COMPREHENSIVE_MAIN", source)
    @test occursin("model == \"bvk2_fixed\"", source)
    @test occursin("adaptive=model == \"bvk2\"", source)
    @test occursin("label=\"AQCE\"", source)
    boundary_source = split(
        split(source, "function aba_cv_write_boundary")[2],
        "function aba_cv_polyline")[1]
    @test !occursin("S/C", boundary_source)
end

@testset "Fixed-stiffness ABA transfer control" begin
    root = joinpath(@__DIR__, "..", "results",
        "bvk2_aba_fixed_stiffness_validation")
    summary = CSV.read(joinpath(root, "summary.csv"), DataFrame)
    profiles = CSV.read(joinpath(root, "profiles.csv"), DataFrame)
    @test nrow(summary) == 7
    @test nrow(profiles) == 7 * 256
    @test all(String.(summary.schema) .==
        "bvk2-aba-fixed-stiffness-validation-v1")
    @test all(String.(summary.model) .== "bvk2_fixed")
    @test all(.!Bool.(summary.adaptive))
    @test all(Float64.(summary.c2) .== 0.16)
    @test all(Int.(summary.nx) .== 128)
    @test all(Bool.(summary.field_gate_pass))
    @test all(Bool.(summary.composition_gate_pass))
    @test all(Bool.(summary.morphology_gate_pass))
    @test all(Bool.(summary.cell_gate_pass))
    @test all(Bool.(summary.accepted))
    @test maximum(Float64.(summary.projected_force_norm)) <= 1.0e-4
    @test maximum(Float64.(summary.projected_force_maxabs)) <= 2.0e-4
    @test isapprox(mean(abs.(Float64.(summary.period_relative_error))) * 100,
        3.83118255562417; atol=1.0e-12, rtol=0.0)
    @test isapprox(mean(Float64.(summary.profile_rms)),
        0.022672651655082292; atol=1.0e-14, rtol=0.0)
end

@testset "Expanded ABA publication artifacts" begin
    root = joinpath(@__DIR__, "..", "results",
        "bvk2_aba_comprehensive_validation")
    benchmark = CSV.read(
        joinpath(root, "benchmark_summary.csv"), DataFrame)
    scan = CSV.read(joinpath(root, "cell_stress_scan.csv"), DataFrame)
    mechanism = CSV.read(
        joinpath(root, "cell_stress_summary.csv"), DataFrame)
    crossing = CSV.read(joinpath(root, "ldis_crossing.csv"), DataFrame)
    @test nrow(benchmark) == 7 * 5
    @test all(Bool.(benchmark.accepted))
    @test sort(unique(Float64.(benchmark.chiN))) ==
        [20.0, 22.0, 25.0, 30.0, 35.0, 40.0, 45.0]
    @test nrow(scan) == 2 * 5 * 14
    @test count(String.(scan.model) .== "scft") == 2 * 14
    @test nrow(mechanism) == 2 * 5
    @test all(Bool.(mechanism.accepted))
    @test nrow(crossing) == 12
    @test sort(unique(Float64.(crossing.chiN))) ==
        [18.25, 18.5, 19.0, 20.0, 21.0, 22.0]
    chi21 = crossing[Float64.(crossing.chiN) .== 21.0, :]
    @test Set(String.(chi21.model)) == Set(("scft", "bvk2"))
    @test all(Float64.(chi21.delta_free_energy) .< 0.0)
    @test all(Float64.(chi21.amplitude) .> 0.0)
    benchmark_svg = read(
        joinpath(root, "macromolecules_figure7.svg"), String)
    mechanism_svg = read(
        joinpath(root, "aba_cell_stress_mechanism.svg"), String)
    @test occursin("period error (\\%)", benchmark_svg)
    @test occursin("profile RMS", benchmark_svg)
    @test occursin("QCE", benchmark_svg)
    @test occursin("AQCE", benchmark_svg)
    @test occursin("phi_A", benchmark_svg)
    @test occursin(raw"$D/R_g$", mechanism_svg)
    @test occursin("cell stress", mechanism_svg)
    @test !occursin("D / D0(model)", mechanism_svg)
end
