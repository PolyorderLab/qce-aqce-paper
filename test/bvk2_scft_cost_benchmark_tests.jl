using Test
using CSV
using DataFrames
using TOML

@testset "Controlled BVK2/SCFT cost benchmark contract" begin
    runner = read(joinpath(@__DIR__, "..", "scripts",
        "run_bvk2_scft_cost_case.jl"), String)
    assembler = read(joinpath(@__DIR__, "..", "scripts",
        "write_bvk2_scft_cost_benchmark.jl"), String)
    @test occursin("JULIA_NUM_THREADS", assembler)
    @test occursin("taskset", assembler)
    @test occursin("/usr/bin/time", assembler)
    @test occursin("tol_stress=1.0e-5", runner)
    @test occursin("tol=1.0e-6", runner)
    @test occursin("maxΔx=0.04", runner)
    @test occursin("Random.Xoshiro(7301)", runner)
    @test occursin("single_period", runner)
    @test occursin("compilation warm-up", assembler)
end

@testset "Controlled BVK2/SCFT cost artifacts" begin
    root = joinpath(@__DIR__, "..", "results",
        "bvk2_scft_cost_benchmark")
    raw = CSV.read(joinpath(root, "raw_runs.csv"), DataFrame)
    summary = CSV.read(joinpath(root, "summary.csv"), DataFrame)
    comparison = CSV.read(joinpath(root, "comparison.csv"), DataFrame)
    environment = TOML.parsefile(
        joinpath(root, "environment.toml"))
    @test nrow(raw) == 4 * 2 * 3
    @test all(Bool.(raw.accepted))
    @test all(Bool.(raw.single_period))
    @test nrow(summary) == 4 * 2
    @test all(Bool.(summary.accepted))
    @test all(Bool.(summary.single_period))
    @test all(Float64.(summary.period_span_rg) .== 0.0)
    @test maximum(Float64.(summary.relative_wall_span)) < 0.30
    @test nrow(comparison) == 4
    @test all(Bool.(comparison.accepted))
    @test all(Float64.(
        comparison.bvk2_over_scft_solver_time) .> 1.0)
    @test environment["julia_threads"] == 1
    @test environment["blas_threads"] == 1
    @test environment["fftw_threads"] == 1
    svg = read(joinpath(
        root, "cost_memory_benchmark.svg"), String)
    @test occursin("solver time (s)", svg)
    @test occursin("peak RSS", svg)
end
