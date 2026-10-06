using CSV
using DataFrames
using Statistics
using Test

const FIXEDCELL_PERF_ROOT = normpath(joinpath(@__DIR__, ".."))
const FIXEDCELL_PERF_SCRIPT = joinpath(
    FIXEDCELL_PERF_ROOT, "scripts",
    "write_bvk2_scft_fixedcell_performance_benchmark.jl")
const FIXEDCELL_PERF_RESULTS = joinpath(
    FIXEDCELL_PERF_ROOT, "results",
    "bvk2_scft_fixedcell_performance")

@testset "BVK2/SCFT fixed-cell performance source contract" begin
    source = read(FIXEDCELL_PERF_SCRIPT, String)
    @test occursin("Polyorder.q!(scft)", source)
    @test occursin("BVK2_WARM_MAX_ITERATIONS = 50", source)
    @test occursin("Polyorder.reset(scft, next_system)", source)
    @test occursin("fixed_cell=true", source)
    @test occursin("warm_seed=true", source)
    @test occursin("continuation_accuracy_audit.csv", source)
    @test !occursin("cell_solve!", source)
end

@testset "BVK2/SCFT kernel benchmark artifacts" begin
    required = (
        "kernel_raw.csv",
        "kernel_summary.csv",
        "kernel_ratios.csv",
        "kernel_benchmark.svg",
        "environment.toml",
        "README.md",
    )
    @test all(isfile(joinpath(FIXEDCELL_PERF_RESULTS, name))
        for name in required)
    ratios = CSV.read(
        joinpath(FIXEDCELL_PERF_RESULTS, "kernel_ratios.csv"), DataFrame)
    @test nrow(ratios) == 12
    @test Set(ratios.dimension) == Set([1, 3])
    @test all(ratios.scft_over_bvk2 .> 1.0)
    @test maximum(ratios.grid_points[ratios.dimension .== 3]) == 48^3
    @test minimum(ratios.scft_over_bvk2) > 3.0
end

@testset "BVK2/SCFT warm-continuation artifacts" begin
    required = (
        "continuation_raw.csv",
        "continuation_summary.csv",
        "continuation_ratios.csv",
        "continuation_accuracy_audit.csv",
        "warm_continuation_benchmark.svg",
    )
    @test all(isfile(joinpath(FIXEDCELL_PERF_RESULTS, name))
        for name in required)
    raw = CSV.read(
        joinpath(FIXEDCELL_PERF_RESULTS, "continuation_raw.csv"),
        DataFrame)
    @test nrow(raw) == 120
    @test all(raw.accepted)
    @test all(raw.fixed_cell)
    @test all(raw.warm_seed)
    @test Set(raw.method) == Set(["bvk2", "scft"])
    @test Set(raw.architecture) == Set(["AB", "ABA"])
    @test Set(raw.trajectory) == Set(1:3)

    ratios = CSV.read(
        joinpath(FIXEDCELL_PERF_RESULTS, "continuation_ratios.csv"),
        DataFrame)
    @test nrow(ratios) == 20
    @test all(ratios.scft_over_bvk2 .> 1.0)
    @test median(ratios.scft_over_bvk2[
        ratios.architecture .== "AB"]) > 4.0
    @test median(ratios.scft_over_bvk2[
        ratios.architecture .== "ABA"]) > 3.0

    audit = CSV.read(
        joinpath(FIXEDCELL_PERF_RESULTS,
            "continuation_accuracy_audit.csv"),
        DataFrame)
    @test nrow(audit) == 2
    @test all(audit.passed)
    @test maximum(audit.absolute_energy_density_difference) <= 1.0e-10
    @test maximum(audit.profile_rms_difference) <= 1.0e-5
    @test maximum(audit.capped_residual_l2) <= 1.0e-4
    @test maximum(audit.capped_residual_linf) <= 2.0e-4
end
