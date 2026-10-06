#!/usr/bin/env julia

# Isolated end-to-end timing unit for the controlled BVK2/Polyorder comparison.
# One process imports only the selected solver stack, performs one compilation
# warm-up, then measures three fresh seeded stress-free lamellar solves.

using FFTW
using LinearAlgebra
using Logging
using Random
using Statistics
using TOML

function arg(name::AbstractString, default=nothing)
    prefix = String(name) * "="
    for value in ARGS
        startswith(value, prefix) &&
            return split(value, "="; limit=2)[2]
    end
    default === nothing &&
        error("missing required argument $(name)")
    return default
end

const METHOD = lowercase(arg("--method"))
const ARCHITECTURE = uppercase(arg("--architecture"))
const CHIN = parse(Float64, arg("--chiN"))
const INITIAL_PERIOD_RG = parse(Float64, arg("--initial-period-rg"))
const OUTPUT = abspath(arg("--output"))
const REPEATS = parse(Int, arg("--repeats", "3"))
const C2 = 0.16

METHOD in ("bvk2", "scft") ||
    error("--method must be bvk2 or scft")
ARCHITECTURE in ("AB", "ABA") ||
    error("--architecture must be AB or ABA")
REPEATS >= 1 || error("--repeats must be positive")

BLAS.set_num_threads(1)
FFTW.set_num_threads(1)

if METHOD == "bvk2"
    @eval using DFMMonteCarlo
elseif METHOD == "scft"
    @eval using Polyorder
end

function aba_system(f::Real, chiN::Real)
    sA, sB = KuhnSegment(:A), KuhnSegment(:B)
    eb1, eb2 = BranchPoint(:EB1), BranchPoint(:EB2)
    a1 = PolymerBlock(:A1, sA, Float64(f) / 2, FreeEnd(:A1), eb1)
    b = PolymerBlock(:B, sB, 1.0 - Float64(f), eb1, eb2)
    a2 = PolymerBlock(:A2, sA, Float64(f) / 2, eb2, FreeEnd(:A2))
    chain = BlockCopolymer(:ABA, [a1, b, a2])
    return PolymerSystem(
        [Component(chain, 1.0, 1.0)],
        Dict([:A, :B] => Float64(chiN)),
    )
end

function polyorder_io_config(outdir)
    save_density = Symbol("save_", Char(0x03d5))
    kwargs = Dict{Symbol,Any}(
        :verbosity => -1,
        :progress_solve => false,
        :progress_cell => false,
        :base_dir => outdir,
        :save_config => false,
        :save_w => false,
        :save_summary => false,
        :save_trace => false,
    )
    kwargs[save_density] = false
    return Polyorder.IOConfig(; kwargs...)
end

function polyorder_density(scft, specie::Symbol)
    density_function = getproperty(Polyorder, Symbol(Char(0x03d5)))
    return vec(getfield(density_function(scft, 1, specie), :data))
end

function run_scft(repetition::Integer)
    return mktempdir() do workdir
        cd(workdir) do
            io = polyorder_io_config(workdir)
            config = Polyorder.Config(;
                io=io,
                scft=Polyorder.SCFTConfig(;
                    min_iter=10,
                    max_iter=4000,
                    tolmode=:Residual,
                    tol=1.0e-6,
                    maxΔx=0.04,
                ),
                cellopt=Polyorder.CellOptConfig(;
                    changeNx=true,
                    tol_stress=1.0e-5,
                ),
            )
            system = ARCHITECTURE == "AB" ?
                AB_system(χN=CHIN, fA=0.5) :
                aba_system(0.5, CHIN)
            scft = NoncyclicChainSCFT(
                system,
                BravaisLattice(UnitCell(INITIAL_PERIOD_RG)),
                0.005;
                mde=Polyorder.OSF,
                spacing=0.04,
                updater=Polyorder.Anderson(
                    Polyorder.SD(1.0); warmup=150, αw=0.2),
                init=:randn,
                rng=Random.Xoshiro(7301),
            )
            field_status = with_logger(NullLogger()) do
                Polyorder.solve!(scft, config)
            end
            cell_updater = Polyorder.VariableCell(
                Polyorder.BB(1.0),
                Polyorder.Anderson(
                    Polyorder.SD(1.0); warmup=80, αw=0.2);
                block=30,
            )
            cell_status, _ = with_logger(NullLogger()) do
                Polyorder.cell_solve!(
                    scft, cell_updater, config)
            end
            phi_a = polyorder_density(scft, :A)
            phi_b = polyorder_density(scft, :B)
            free_energy = Float64(Polyorder.F(scft))
            delta_free_energy =
                free_energy - Float64(Polyorder.F_DIS(system))
            stress = norm(Polyorder.gradient_wrt_cell(scft))
            period_rg =
                Float64(Polyorder.unitcell(scft).edges[1])
            single_period =
                abs(period_rg / INITIAL_PERIOD_RG - 1.0) <= 0.10
            accepted =
                occursin("Successful", string(field_status)) &&
                occursin("Successful", string(cell_status)) &&
                stress <= 1.0e-4 &&
                delta_free_energy < -1.0e-6 &&
                maximum(phi_a)-minimum(phi_a) > 1.0e-3 &&
                sqrt(mean(abs2, phi_a .+ phi_b .- 1.0)) <= 1.0e-6 &&
                single_period
            return Dict{String,Any}(
                "period_rg" => period_rg,
                "energy_density" => free_energy,
                "delta_free_energy" => delta_free_energy,
                "cell_stress" => stress,
                "grid_points" => length(phi_a),
                "field_status" => string(field_status),
                "cell_status" => string(cell_status),
                "field_converged" =>
                    occursin("Successful", string(field_status)),
                "cell_converged" =>
                    occursin("Successful", string(cell_status)),
                "single_period" => single_period,
                "accepted" => accepted,
            )
        end
    end
end

function bvk2_reference_period()
    params = ARCHITECTURE == "AB" ?
        diblock_bvk2_parameters(; f=0.5, c2=C2) :
        triblock_aba_bvk2_parameters(; f=0.5, c2=C2)
    return 2.0 * pi / params.k_star
end

function run_bvk2(repetition::Integer)
    reference_period = bvk2_reference_period()
    target_factor =
        INITIAL_PERIOD_RG / (sqrt(6.0) * reference_period)
    common = (
        f=0.5,
        chiN=CHIN,
        nx=128,
        mode_count=16,
        initial_amplitudes=(0.35,),
        c2=C2,
        max_iterations=4000,
        max_period_iterations=16,
        bootstrap_period_factors=(target_factor,),
        bootstrap_window=0.12,
        local_check_fraction=0.01,
        force_tolerance=1.0e-4,
        force_maxabs_tolerance=2.0e-4,
    )
    result = if ARCHITECTURE == "AB"
        minimize_diblock_bvk2_lamella_stress_free(; common...)
    else
        minimize_triblock_aba_bvk2_lamella_stress_free(;
            common..., period_strategy=:bootstrap_local)
    end
    field = result.result
    single_period =
        abs(result.period_rg / INITIAL_PERIOD_RG - 1.0) <= 0.10
    accepted =
        result.local_minimum_check_pass &&
        field.converged &&
        field.energy < field.homogeneous_energy - 1.0e-8 &&
        field.maximum_phi-field.minimum_phi > 1.0e-3 &&
        single_period
    return Dict{String,Any}(
        "period_rg" => Float64(result.period_rg),
        "energy_density" => Float64(result.objective),
        "delta_free_energy" =>
            Float64(field.energy-field.homogeneous_energy),
        "cell_stress" => NaN,
        "grid_points" => length(field.phi_a),
        "field_status" => string(field.converged),
        "cell_status" => string(result.local_minimum_check_pass),
        "field_converged" => Bool(field.converged),
        "cell_converged" =>
            Bool(result.local_minimum_check_pass),
        "single_period" => single_period,
        "accepted" => accepted,
    )
end

run_once(repetition) = METHOD == "bvk2" ?
    run_bvk2(repetition) : run_scft(repetition)

mkpath(dirname(OUTPUT))
warmup = @timed run_once(0)
warmup.value["accepted"] ||
    error("warm-up solve failed its acceptance contract")
GC.gc()

runs = Dict{String,Any}[]
for repetition in 1:REPEATS
    GC.gc()
    live_before = Base.gc_live_bytes()
    measurement = @timed run_once(repetition)
    result = measurement.value
    result["accepted"] ||
        error("measured repetition $(repetition) failed acceptance")
    result["repetition"] = repetition
    result["wall_s"] = measurement.time
    result["gc_s"] = measurement.gctime
    result["allocated_bytes"] = measurement.bytes
    result["gc_live_bytes_before"] = live_before
    result["gc_live_bytes_after"] = Base.gc_live_bytes()
    push!(runs, result)
end

wall = Float64[run["wall_s"] for run in runs]
allocations = Float64[run["allocated_bytes"] for run in runs]
periods = Float64[run["period_rg"] for run in runs]
payload = Dict{String,Any}(
    "schema" => "bvk2-scft-controlled-cost-case-v1",
    "method" => METHOD,
    "architecture" => ARCHITECTURE,
    "chiN" => CHIN,
    "fA" => 0.5,
    "c2" => METHOD == "bvk2" ? C2 : "not_applicable",
    "solver_package_version" => METHOD == "bvk2" ?
        string(pkgversion(DFMMonteCarlo)) :
        string(pkgversion(Polyorder)),
    "julia_threads" => Threads.nthreads(),
    "blas_threads" => BLAS.get_num_threads(),
    "fftw_threads" => 1,
    "repeats" => REPEATS,
    "rng_seed" => 7301,
    "warmup_wall_s" => warmup.time,
    "warmup_allocated_bytes" => warmup.bytes,
    "median_solver_wall_s" => median(wall),
    "minimum_solver_wall_s" => minimum(wall),
    "maximum_solver_wall_s" => maximum(wall),
    "median_allocated_bytes" => median(allocations),
    "median_period_rg" => median(periods),
    "period_span_rg" => maximum(periods)-minimum(periods),
    "runs" => runs,
)
open(OUTPUT, "w") do io
    TOML.print(io, payload; sorted=true)
end
println("wrote ", OUTPUT)
