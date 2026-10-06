const MONTE_CARLO_PROJECT = normpath(joinpath(@__DIR__, ".."))
MONTE_CARLO_PROJECT in LOAD_PATH || pushfirst!(LOAD_PATH, MONTE_CARLO_PROJECT)

# Reuse the validated Polyorder cell solve and periodic-profile alignment used by
# the existing own-cell lamellar benchmark.  Its entry point is guarded, so
# including it does not start the legacy calculation.
include(joinpath(@__DIR__, "write_diblock_stress_free_lamella_comparison.jl"))
isdefined(@__MODULE__, :BurpTiLamellarPeriodSettings) ||
    include(joinpath(@__DIR__, "burp_ti_lamellar_period_workflow.jl"))

const LPC_SCHEMA = "liu2019-stress-free-period-map-v2"
const LPC_SCFT_RESIDUAL_TOL = 1.0e-6
const LPC_UD_PROTOCOL = "mean-constrained fused-spectral L-BFGS with profile continuation"
const LPC_UD_RESIDUAL_TOL = 1.0e-6
const LPC_UD_MAX_RESIDUAL_TOL = 1.0e-5
const LPC_BVK2_PROTOCOL =
    "oversampled filtered UD-theta L-BFGS with cosine-Newton polish"
const LPC_BVK2_PERIOD_OPTIMIZER =
    "FilteredUDThetaAnalyticCellStressRoot"
const LPC_BVK2_RESIDUAL_TOL = 1.0e-7
const LPC_BVK2_MAX_RESIDUAL_TOL = 5.0e-7
const LPC_BVK2_STRESS_TOL = 1.0e-8
const LPC_BVK2_OVERSAMPLE = 2
const LPC_BVK2_SENSOR_FILTER_RATIO = 12.0
const LPC_SCFT_ANDERSON_DEPTH = 10
const LPC_SCFT_WARMUP = 100
const LPC_SCFT_SD_STEP = 0.1
const LPC_SCFT_ANDERSON_WARMUP_STEP = 0.2
const LPC_SPINODALS = [
    (f=0.20, chiN_spinodal=24.613),
    (f=0.25, chiN_spinodal=18.172),
    (f=0.30, chiN_spinodal=14.635),
    (f=0.35, chiN_spinodal=12.562),
    (f=0.40, chiN_spinodal=11.344),
    (f=0.45, chiN_spinodal=10.698),
    (f=0.50, chiN_spinodal=10.495),
]
const LPC_MODEL_ORDER = ["scft", "liu2019_opf", "ohta_kawasaki",
    "uneyama_doi",
    "burp_ti", "bvk2"]
const LPC_MODEL_LABELS = Dict(
    "scft" => "SCFT",
    "liu2019_opf" => "OPF (force/stress mapped)",
    "ohta_kawasaki" => "Ohta-Kawasaki (recomputed)",
    "uneyama_doi" => "Uneyama-Doi",
    "burp_ti" => "BURP-TI",
    "bvk2" => "BVK2",
)

const LPC_USAGE = """
Recompute the Liu et al. Figure 2 stress-free lamellar-period map.

Usage:
  julia --project=. scripts/write_liu2019_stress_free_period_map.jl [options]

State selection:
  --compositions=LIST             Comma-separated fA values (default: 0.20,...,0.50)
  --cases=LIST                    Explicit fA:chiN pairs; overrides --compositions
  --max-chiN=35                  Upper bound for default cases and explicit-case gate
  --chi-order=ascending|descending
  --models=LIST                   liu2019_opf,ohta_kawasaki,uneyama_doi,burp_ti,bvk2
  --allow-ok-extrapolation=false  Explicitly permit the Liu-2019 OK coefficient
                                  mapping above its published chiN <= 35 domain

Numerics:
  --nx=128                        OPF, OK, and BURP-TI grid size
  --bvk2-nx=128                   BVK2 grid size
  --nquad=64                      BURP-TI production quadrature order
  --burp-max-newton-iterations=256
  --ok-gradient-tol=8e-6         OK nonlinear-CG stopping tolerance
  --bvk2-c2=0.16                  Frozen BVK2 nonlinear coefficient
  --bvk2-oversample=2              UD-theta quadrature refinement factor
  --bvk2-sensor-filter-ratio=12    Grid-validated sensor cutoff / k_star
  --bvk2-energy-tol=0              Fixed-cell relative-energy stop tolerance
  --bvk2-force-tol=1e-7            BVK2 projected-force RMS tolerance
  --bvk2-force-maxabs-tol=5e-7     BVK2 projected-force maximum tolerance
  --bvk2-cell-stress-tol=1e-8      BVK2 analytic cell-stress tolerance
  --bvk2-log-period-tol=1e-7       BVK2 log-period bracket tolerance
  --model-max-iterations=10000
  --model-max-period-iterations=24
  --lower-factor=0.75             Lower period/reference-period ratio
  --upper-factor=1.8              Upper period/reference-period ratio
  --period-strategy=bootstrap_local  Strategy for non-OPF reduced models
  --model-seed-policy=history        history or reference; reference resets the
                                     period seed to the model's weak-response
                                     period while retaining field continuation
  --bootstrap-window=0.18
  --local-check-fraction=0.01
  --spacing=0.1                   Requested Polyorder spacing in Rg
  --ds=0.01                       Polyorder contour step
  --scft-presolve-max-iterations=3000
  --scft-cell-max-iterations=8000
  --scft-cell-block=5
  --scft-solver-residual-tol=1e-6
  --scft-cell-bb-type=2
  --scft-stress-tol=1e-5
  --sample-count=256              Common aligned-profile grid
  --seed=3083

Checkpoint and output control:
  --outdir=PATH                   Campaign directory
  --resume=true                   Reuse fully accepted requested checkpoints
  --reuse-existing-scft=false     Load accepted SCFT data when rerunning models
  --assemble-only=false           Merge case tables without solving
  --merge-outputs=true            Rebuild root tables after case writes
  --write-contract=true           Write the campaign README

OPF always uses per-state Eq. 18 force/stress mapping (lambda=100, c2=-1),
unbounded mean-constrained Fourier NCG field relaxation, and direct Eq. S2
stress-root bisection. BURP-TI uses bounded-KKT Newton-GMRES field relaxation
and roots its analytic cell stress; --period-strategy applies to BVK2 only.

Examples:
  # Resumable full campaign
  julia --project=. scripts/write_liu2019_stress_free_period_map.jl

  # Rerun OPF only from accepted SCFT checkpoints, then assemble
  julia --project=. scripts/write_liu2019_stress_free_period_map.jl \\
    --models=liu2019_opf --reuse-existing-scft=true --resume=false \\
    --merge-outputs=false --write-contract=false
  julia --project=. scripts/write_liu2019_stress_free_period_map.jl \\
    --assemble-only=true
"""

function lpc_parse_strings(text::AbstractString)
    return [strip(item) for item in split(text, ",") if !isempty(strip(item))]
end

function lpc_default_cases(compositions::AbstractVector{<:Real};
        chi_order::Symbol=:ascending, max_chiN::Real=35.0)
    chi_order in (:ascending, :descending) || throw(ArgumentError(
        "chi_order must be ascending or descending"))
    requested = Set(round(Float64(f); digits=8) for f in compositions)
    cases = NamedTuple[]
    for entry in reverse(LPC_SPINODALS)
        round(entry.f; digits=8) in requested || continue
        upper = floor(Int, Float64(max_chiN))
        chi_values = collect((floor(Int, entry.chiN_spinodal) + 1):upper)
        chi_order == :descending && reverse!(chi_values)
        for chiN in chi_values
            push!(cases, (f=entry.f, chiN=Float64(chiN),
                chiN_spinodal=entry.chiN_spinodal))
        end
    end
    isempty(cases) && throw(ArgumentError(
        "no requested compositions occur in the Liu et al. domain"))
    return cases
end

function lpc_explicit_cases(text::AbstractString; chi_order::Symbol=:ascending,
        max_chiN::Real=35.0)
    parsed = _parse_cases(text)
    out = NamedTuple[]
    for case in parsed
        entry = only(filter(item -> isapprox(item.f, case.f; atol=1.0e-10),
            LPC_SPINODALS))
        case.chiN > entry.chiN_spinodal || throw(DomainError(case,
            "chiN must lie above the RPA spinodal"))
        case.chiN <= Float64(max_chiN) || throw(DomainError(case,
            "chiN exceeds the requested --max-chiN=$(max_chiN)"))
        push!(out, (f=case.f, chiN=case.chiN,
            chiN_spinodal=entry.chiN_spinodal))
    end
    sort!(out; by=case -> (-case.f,
        chi_order == :ascending ? case.chiN : -case.chiN))
    return out
end

function lpc_case_dir(outdir::AbstractString, case)
    return joinpath(outdir, "cases", _case_id(case.f, case.chiN))
end

function lpc_bool(value)
    value isa Bool && return value
    return lowercase(strip(String(value))) in ("true", "t", "yes", "y", "1")
end

function lpc_bvk2_row_gate(row;
        oversample::Integer=LPC_BVK2_OVERSAMPLE,
        sensor_filter_ratio::Real=LPC_BVK2_SENSOR_FILTER_RATIO,
        residual_tol::Real=LPC_BVK2_RESIDUAL_TOL,
        maximum_residual_tol::Real=LPC_BVK2_MAX_RESIDUAL_TOL,
        stress_tol::Real=LPC_BVK2_STRESS_TOL)
    hasproperty(row, :reduced_field_protocol) &&
        !ismissing(row.reduced_field_protocol) &&
        String(row.reduced_field_protocol) == LPC_BVK2_PROTOCOL ||
        return false
    hasproperty(row, :period_optimizer) &&
        !ismissing(row.period_optimizer) &&
        String(row.period_optimizer) == LPC_BVK2_PERIOD_OPTIMIZER || return false
    hasproperty(row, :oversample_factor) &&
        !ismissing(row.oversample_factor) &&
        Int(row.oversample_factor) == Int(oversample) || return false
    hasproperty(row, :sensor_filter_ratio) &&
        !ismissing(row.sensor_filter_ratio) &&
        Float64(row.sensor_filter_ratio) == Float64(sensor_filter_ratio) ||
        return false
    hasproperty(row, :discretization_schema) &&
        !ismissing(row.discretization_schema) &&
        String(row.discretization_schema) ==
            DFMMonteCarlo._bvk2_ud_theta_discretization_schema(
                sensor_filter_ratio) || return false
    hasproperty(row, :reduced_field_residual_norm) &&
        !ismissing(row.reduced_field_residual_norm) &&
        isfinite(Float64(row.reduced_field_residual_norm)) &&
        Float64(row.reduced_field_residual_norm) <= Float64(residual_tol) ||
        return false
    hasproperty(row, :reduced_field_residual_maxabs) &&
        !ismissing(row.reduced_field_residual_maxabs) &&
        isfinite(Float64(row.reduced_field_residual_maxabs)) &&
        Float64(row.reduced_field_residual_maxabs) <=
            Float64(maximum_residual_tol) || return false
    hasproperty(row, :stress_norm) && !ismissing(row.stress_norm) &&
        isfinite(Float64(row.stress_norm)) &&
        abs(Float64(row.stress_norm)) <= Float64(stress_tol) || return false
    return true
end

function lpc_existing_case_accepted(case_dir::AbstractString,
        required_models::AbstractVector{<:AbstractString};
        burp_nquad::Integer=64, burp_max_newton_iterations::Integer=256,
        bvk2_oversample::Integer=LPC_BVK2_OVERSAMPLE,
        bvk2_sensor_filter_ratio::Real=LPC_BVK2_SENSOR_FILTER_RATIO,
        bvk2_residual_tol::Real=LPC_BVK2_RESIDUAL_TOL,
        bvk2_maximum_residual_tol::Real=LPC_BVK2_MAX_RESIDUAL_TOL,
        bvk2_stress_tol::Real=LPC_BVK2_STRESS_TOL)
    summary_path = joinpath(case_dir, "summary.csv")
    profiles_path = joinpath(case_dir, "profiles.csv")
    native_path = joinpath(case_dir, "native_profiles.csv")
    all(isfile, (summary_path, profiles_path, native_path)) || return false
    table = try
        CSV.read(summary_path, DataFrame)
    catch
        return false
    end
    all(name -> name in propertynames(table), (:schema, :model, :status)) ||
        return false
    all(String(row.schema) == LPC_SCHEMA for row in eachrow(table)) ||
        return false
    present = Set(String.(table.model))
    all(model -> model in present, ["scft"; String.(required_models)]) ||
        return false
    requested = ["scft"; String.(required_models)]
    for model in requested
        rows = [row for row in eachrow(table) if String(row.model) == model]
        length(rows) == 1 || return false
        row = only(rows)
        String(row.status) == "accepted" || return false
        if model == "uneyama_doi"
            hasproperty(row, :reduced_field_protocol) &&
                String(row.reduced_field_protocol) == LPC_UD_PROTOCOL || return false
            hasproperty(row, :reduced_field_residual_norm) &&
                !ismissing(row.reduced_field_residual_norm) &&
                isfinite(Float64(row.reduced_field_residual_norm)) &&
                Float64(row.reduced_field_residual_norm) <= LPC_UD_RESIDUAL_TOL ||
                return false
            hasproperty(row, :reduced_field_residual_maxabs) &&
                !ismissing(row.reduced_field_residual_maxabs) &&
                isfinite(Float64(row.reduced_field_residual_maxabs)) &&
                Float64(row.reduced_field_residual_maxabs) <=
                    LPC_UD_MAX_RESIDUAL_TOL || return false
            hasproperty(row, :morphology_gate_pass) &&
                lpc_bool(row.morphology_gate_pass) || return false
        end
        if model == "burp_ti"
            hasproperty(row, :reduced_field_protocol) &&
                String(row.reduced_field_protocol) ==
                "bounded-KKT Newton-GMRES" || return false
            hasproperty(row, :period_optimizer) &&
                String(row.period_optimizer) ==
                "BoundedKKTNewtonGMRES+AnalyticStressRoot" || return false
            hasproperty(row, :quadrature_order) &&
                !ismissing(row.quadrature_order) &&
                Int(row.quadrature_order) == Int(burp_nquad) || return false
            hasproperty(row, :field_iteration_cap) &&
                !ismissing(row.field_iteration_cap) &&
                Int(row.field_iteration_cap) ==
                Int(burp_max_newton_iterations) || return false
            hasproperty(row, :branch_dominant_mode) &&
                !ismissing(row.branch_dominant_mode) &&
                Int(row.branch_dominant_mode) == 1 || return false
            hasproperty(row, :stress_norm) && !ismissing(row.stress_norm) &&
                isfinite(Float64(row.stress_norm)) &&
                abs(Float64(row.stress_norm)) <= 1.0e-4 || return false
        end
        model == "bvk2" && !lpc_bvk2_row_gate(row;
            oversample=bvk2_oversample,
            sensor_filter_ratio=bvk2_sensor_filter_ratio,
            residual_tol=bvk2_residual_tol,
            maximum_residual_tol=bvk2_maximum_residual_tol,
            stress_tol=bvk2_stress_tol) && return false
    end
    return true
end

function lpc_load_existing_scft(case_dir::AbstractString)
    summary_path = joinpath(case_dir, "summary.csv")
    native_path = joinpath(case_dir, "native_profiles.csv")
    all(isfile, (summary_path, native_path)) || return nothing
    summary = CSV.read(summary_path, DataFrame)
    candidates = [row for row in eachrow(summary)
        if String(row.schema) == LPC_SCHEMA && String(row.model) == "scft" &&
            String(row.status) == "accepted"]
    length(candidates) == 1 || return nothing
    row = only(candidates)
    native = CSV.read(native_path, DataFrame)
    native = native[[String(item.schema) == LPC_SCHEMA &&
        String(item.model) == "scft" for item in eachrow(native)], :]
    isempty(native) && return nothing
    sort!(native, :native_index)
    phi_a = Float64.(native.phi_a)
    phi_b = Float64.(native.phi_b)
    return (
        convergence=String(row.scft_convergence),
        presolve_convergence="reused_validated_checkpoint",
        cell_convergence=String(row.scft_convergence),
        phi_a=phi_a, phi_b=phi_b, free_energy=Float64(row.energy),
        initial_period_rg=Float64(row.period_rg),
        period_rg=Float64(row.period_rg), grid_count=length(phi_a),
        cell_iterations=ismissing(row.period_iterations) ? 0 :
            Int(row.period_iterations), field_evaluations=missing,
        stress_norm=Float64(row.stress_norm),
        mean_phi_a=Float64(row.mean_phi_a),
        minimum_phi_a=Float64(row.minimum_phi_a),
        maximum_phi_a=Float64(row.maximum_phi_a),
        incompressibility_rms=Float64(row.scft_incompressibility_rms),
        scft_model=nothing,
        residual_norm=Float64(row.scft_residual_norm),
        requested_spacing_rg=Float64(row.scft_requested_spacing_rg),
        actual_spacing_rg=Float64(row.scft_actual_spacing_rg),
        contour_step=Float64(row.scft_contour_step),
    )
end

function lpc_load_existing_model_profile(case_dir::AbstractString,
        model::AbstractString, expected_count::Integer;
        resample_if_needed::Bool=false)
    native_path = joinpath(case_dir, "native_profiles.csv")
    isfile(native_path) || return nothing
    native = try
        CSV.read(native_path, DataFrame)
    catch
        return nothing
    end
    required = (:schema, :model, :native_index, :phi_a)
    all(name -> name in propertynames(native), required) || return nothing
    rows = native[[String(row.schema) == LPC_SCHEMA &&
        String(row.model) == String(model) for row in eachrow(native)], :]
    isempty(rows) && return nothing
    sort!(rows, :native_index)
    values = Float64.(rows.phi_a)
    all(isfinite, values) || return nothing
    length(values) == Int(expected_count) && return values
    resample_if_needed || return nothing
    return _resample(values, Int(expected_count))
end

function lpc_load_existing_model_period_factor(case_dir::AbstractString,
        model::AbstractString)
    summary_path = joinpath(case_dir, "summary.csv")
    isfile(summary_path) || return nothing
    summary = try
        CSV.read(summary_path, DataFrame)
    catch
        return nothing
    end
    rows = [row for row in eachrow(summary)
        if hasproperty(row, :schema) && String(row.schema) == LPC_SCHEMA &&
            String(row.model) == String(model)]
    length(rows) == 1 || return nothing
    factor = only(rows).period_factor
    (ismissing(factor) || !isfinite(Float64(factor))) && return nothing
    return Float64(factor)
end

function lpc_existing_model_accepted(case_dir::AbstractString,
        model::AbstractString)
    summary_path = joinpath(case_dir, "summary.csv")
    isfile(summary_path) || return false
    summary = try
        CSV.read(summary_path, DataFrame)
    catch
        return false
    end
    rows = [row for row in eachrow(summary)
        if hasproperty(row, :schema) && String(row.schema) == LPC_SCHEMA &&
            String(row.model) == String(model)]
    length(rows) == 1 || return false
    row = only(rows)
    return hasproperty(row, :status) && String(row.status) == "accepted"
end

function lpc_model_history(outdir::AbstractString,
        requested_models::AbstractVector{<:AbstractString};
        burp_nquad=nothing, burp_max_newton_iterations=nothing,
        bvk2_oversample::Integer=LPC_BVK2_OVERSAMPLE,
        bvk2_sensor_filter_ratio::Real=LPC_BVK2_SENSOR_FILTER_RATIO,
        bvk2_residual_tol::Real=LPC_BVK2_RESIDUAL_TOL,
        bvk2_maximum_residual_tol::Real=LPC_BVK2_MAX_RESIDUAL_TOL,
        bvk2_stress_tol::Real=LPC_BVK2_STRESS_TOL)
    history = Dict{Tuple{String,Float64},Vector{NamedTuple}}()
    cases_dir = joinpath(outdir, "cases")
    isdir(cases_dir) || return history
    for directory in readdir(cases_dir; join=true)
        path = joinpath(directory, "summary.csv")
        isfile(path) || continue
        table = try
            CSV.read(path, DataFrame)
        catch
            continue
        end
        for row in eachrow(table)
            model = String(row.model)
            model in requested_models || continue
            hasproperty(row, :schema) && String(row.schema) == LPC_SCHEMA ||
                continue
            hasproperty(row, :status) && String(row.status) == "accepted" ||
                continue
            if model == "uneyama_doi"
                hasproperty(row, :reduced_field_protocol) &&
                    String(row.reduced_field_protocol) == LPC_UD_PROTOCOL || continue
                hasproperty(row, :reduced_field_residual_norm) &&
                    !ismissing(row.reduced_field_residual_norm) &&
                    isfinite(Float64(row.reduced_field_residual_norm)) &&
                    Float64(row.reduced_field_residual_norm) <= LPC_UD_RESIDUAL_TOL ||
                    continue
                hasproperty(row, :reduced_field_residual_maxabs) &&
                    !ismissing(row.reduced_field_residual_maxabs) &&
                    isfinite(Float64(row.reduced_field_residual_maxabs)) &&
                    Float64(row.reduced_field_residual_maxabs) <=
                        LPC_UD_MAX_RESIDUAL_TOL || continue
                hasproperty(row, :morphology_gate_pass) &&
                    lpc_bool(row.morphology_gate_pass) || continue
            end
            if model == "burp_ti"
                hasproperty(row, :reduced_field_protocol) &&
                    String(row.reduced_field_protocol) ==
                    "bounded-KKT Newton-GMRES" || continue
                hasproperty(row, :period_optimizer) &&
                    String(row.period_optimizer) ==
                    "BoundedKKTNewtonGMRES+AnalyticStressRoot" || continue
                hasproperty(row, :branch_dominant_mode) &&
                    !ismissing(row.branch_dominant_mode) &&
                    Int(row.branch_dominant_mode) == 1 || continue
                burp_nquad === nothing ||
                    (hasproperty(row, :quadrature_order) &&
                     !ismissing(row.quadrature_order) &&
                     Int(row.quadrature_order) == Int(burp_nquad)) || continue
                burp_max_newton_iterations === nothing ||
                    (hasproperty(row, :field_iteration_cap) &&
                     !ismissing(row.field_iteration_cap) &&
                     Int(row.field_iteration_cap) ==
                        Int(burp_max_newton_iterations)) || continue
            end
            model == "bvk2" && !lpc_bvk2_row_gate(row;
                oversample=bvk2_oversample,
                sensor_filter_ratio=bvk2_sensor_filter_ratio,
                residual_tol=bvk2_residual_tol,
                maximum_residual_tol=bvk2_maximum_residual_tol,
                stress_tol=bvk2_stress_tol) && continue
            hasproperty(row, :period_factor) || continue
            factor = row.period_factor
            (ismissing(factor) || !isfinite(Float64(factor))) && continue
            key = (model, Float64(row.f))
            push!(get!(history, key, NamedTuple[]),
                (chiN=Float64(row.chiN), factor=Float64(factor)))
        end
    end
    for rows in values(history)
        unique_rows = Dict(row.chiN => row for row in rows)
        empty!(rows)
        append!(rows, values(unique_rows))
        sort!(rows; by=row -> row.chiN)
    end
    return history
end


function lpc_model_profile_history(outdir::AbstractString,
        requested_models::AbstractVector{<:AbstractString}, nx::Integer;
        burp_nquad::Integer=64, burp_max_newton_iterations::Integer=256)
    profiles = Dict{Tuple{String,Float64},Vector{NamedTuple}}()
    cases_dir = joinpath(outdir, "cases")
    isdir(cases_dir) || return profiles
    for directory in readdir(cases_dir; join=true)
        lpc_existing_case_accepted(directory, requested_models;
            burp_nquad=burp_nquad,
            burp_max_newton_iterations=burp_max_newton_iterations) || continue
        summary = try
            CSV.read(joinpath(directory, "summary.csv"), DataFrame)
        catch
            continue
        end
        all(name -> name in propertynames(summary), (:model, :f, :chiN)) ||
            continue
        for model in requested_models
            model in ("burp_ti", "uneyama_doi") || continue
            row = only(row for row in eachrow(summary)
                if String(row.model) == model)
            profile = lpc_load_existing_model_profile(directory, model, nx)
            profile === nothing && continue
            key = (String(model), Float64(row.f))
            push!(get!(profiles, key, NamedTuple[]),
                (chiN=Float64(row.chiN), profile=profile))
        end
    end
    for rows in values(profiles)
        sort!(rows; by=row -> row.chiN)
    end
    return profiles
end

function lpc_nearest_lower_model_profile(profile_history, key, chiN::Real)
    rows = get(profile_history, key, NamedTuple[])
    lower = [row for row in rows if row.chiN < Float64(chiN)]
    isempty(lower) && return nothing
    return last(sort(lower; by=row -> row.chiN)).profile
end


function lpc_nearest_model_profile(profile_history, key, chiN::Real)
    rows = get(profile_history, key, NamedTuple[])
    isempty(rows) && return nothing
    target = Float64(chiN)
    return rows[argmin(abs(row.chiN - target) for row in rows)].profile
end

function lpc_model_bootstrap_period_factors(history, key, chiN::Real;
        lower_factor::Real, upper_factor::Real,
        prefer_same_state::Bool=false)
    rows = get(history, key, NamedTuple[])
    if prefer_same_state
        same_state = [row for row in rows
            if isapprox(row.chiN, Float64(chiN); atol=1.0e-10, rtol=0.0)]
        if !isempty(same_state)
            factor = clamp(last(same_state).factor,
                Float64(lower_factor), Float64(upper_factor))
            return Float64[factor]
        end
    end
    return _bootstrap_period_factors(rows, chiN;
        lower_factor=lower_factor, upper_factor=upper_factor)
end

function lpc_scft_period_history(outdir::AbstractString)
    history = Dict{Float64,Vector{NamedTuple}}()
    cases_dir = joinpath(outdir, "cases")
    isdir(cases_dir) || return history
    for directory in readdir(cases_dir; join=true)
        path = joinpath(directory, "summary.csv")
        isfile(path) || continue
        table = try
            CSV.read(path, DataFrame)
        catch
            continue
        end
        for row in eachrow(table)
            hasproperty(row, :schema) && String(row.schema) == LPC_SCHEMA ||
                continue
            String(row.model) == "scft" || continue
            period = Float64(row.period_rg)
            amplitude = Float64(row.maximum_phi_a) - Float64(row.minimum_phi_a)
            isfinite(period) && amplitude > 1.0e-4 || continue
            push!(get!(history, Float64(row.f), NamedTuple[]),
                (chiN=Float64(row.chiN), period_rg=period))
        end
    end
    for values in values(history)
        sort!(values; by=value -> value.chiN)
    end
    return history
end

function lpc_nearest_scft_period(history, f::Real, chiN::Real, fallback::Real)
    values = get(history, Float64(f), NamedTuple[])
    isempty(values) && return Float64(fallback)
    return values[argmin(abs(value.chiN - Float64(chiN)) for value in values)].period_rg
end

function lpc_scft_initial_period(kernel_period_rg::Real, chiN::Real,
        chiN_spinodal::Real)
    reduced_distance = clamp(
        (Float64(chiN) - Float64(chiN_spinodal)) /
        (35.0 - Float64(chiN_spinodal)), 0.0, 1.0)
    return Float64(kernel_period_rg) *
        (1.0 + 0.45 * sqrt(reduced_distance))
end

function lpc_model_calibration_role(model::AbstractString, f::Real, chiN::Real)
    model == "liu2019_opf" && return "per_state_force_stress_mapping"
    model == "ohta_kawasaki" && return Float64(chiN) > 35.0 ?
        "liu2019_ok_asymptotic_quadratic_extrapolated_c3_c4" :
        "liu2019_ok_asymptotic_quadratic_mapped_c3_c4"
    model == "bvk2" || return "not_applicable"
    training = ((0.5, 15.0), (0.5, 25.0), (0.35, 30.0))
    holdout = ((0.5, 12.0), (0.5, 20.0), (0.5, 30.0))
    any(isapprox(f, point[1]; atol=1.0e-12) &&
        isapprox(chiN, point[2]; atol=1.0e-12) for point in training) &&
        return "training"
    any(isapprox(f, point[1]; atol=1.0e-12) &&
        isapprox(chiN, point[2]; atol=1.0e-12) for point in holdout) &&
        return "holdout"
    return "frozen_c2_evaluation"
end

function lpc_seed_lamellar_fields!(scft, chiN::Real)
    length(scft.wfields) >= 2 || throw(ArgumentError(
        "AB lamellar initialization requires at least two auxiliary fields"))
    grid_count = length(scft.wfields[1])
    phase = 2.0 * pi .* (0:(grid_count - 1)) ./ grid_count
    amplitude = clamp(0.25 * Float64(chiN), 2.0, 10.0)
    seed = amplitude .* cos.(phase)
    scft.wfields[1].data .= reshape(-seed, size(scft.wfields[1].data))
    scft.wfields[2].data .= reshape(seed, size(scft.wfields[2].data))
    for field in scft.wfields[3:end]
        fill!(field.data, 0.0)
    end
    return scft
end

function lpc_polyorder_cell_solved_profile(; f::Real, chiN::Real,
        initial_Lrg::Real, spacing::Real, ds::Real, seed::Integer,
        presolve_max_iter::Integer, cell_max_iter::Integer,
        cell_block::Integer, stress_tol::Real, previous_scft=nothing,
        cell_bb_step_max::Real=1.0, cell_bb_initial_step::Real=0.2,
        cell_bb_type::Integer=2,
        solver_residual_tol::Real=LPC_SCFT_RESIDUAL_TOL,
        anderson_droptol::Real=0.0,
        cell_algorithm::Symbol=:bb)
    cell_algorithm == :bb || throw(ArgumentError(
        "the production SCFT protocol requires BB cell relaxation"))
    cell_bb_type in (1, 2, 3) || throw(ArgumentError(
        "cell_bb_type must be 1, 2, or 3"))
    0.0 < solver_residual_tol <= LPC_SCFT_RESIDUAL_TOL ||
        throw(ArgumentError("solver_residual_tol must lie in (0, " *
            "$(LPC_SCFT_RESIDUAL_TOL)]"))
    io = _polyorder_io_config(mktempdir())
    presolve_config = Polyorder.Config(; io=io,
        scft=Polyorder.SCFTConfig(; min_iter=10,
            max_iter=Int(presolve_max_iter), tolmode=:Residual,
            tol=Float64(solver_residual_tol),
            maxΔx=Float64(spacing)),
        cellopt=Polyorder.CellOptConfig(; changeNx=false,
            tol_stress=Float64(stress_tol)))
    cell_config = Polyorder.Config(; io=io,
        scft=Polyorder.SCFTConfig(; min_iter=10,
            max_iter=Int(cell_max_iter), tolmode=:Residual,
            tol=Float64(solver_residual_tol),
            maxΔx=Float64(spacing)),
        cellopt=Polyorder.CellOptConfig(; changeNx=false,
            tol_stress=Float64(stress_tol)))
    system = _polyorder_ab_system(; chiN=chiN, fA=f)
    lattice = Polyorder.BravaisLattice(Polyorder.UnitCell(Float64(initial_Lrg)))
    field_updater = Polyorder.Anderson(Polyorder.SD(LPC_SCFT_SD_STEP);
        m=LPC_SCFT_ANDERSON_DEPTH, warmup=LPC_SCFT_WARMUP,
        αw=LPC_SCFT_ANDERSON_WARMUP_STEP,
        droptol=Float64(anderson_droptol))
    scft = if previous_scft === nothing
        fresh = Polyorder.NoncyclicChainSCFT(system, lattice, Float64(ds);
            mde=Polyorder.OSF, spacing=Float64(spacing),
            updater=field_updater, init=:zeros,
            rng=Random.Xoshiro(seed))
        lpc_seed_lamellar_fields!(fresh, chiN)
    else
        continued = Polyorder.NoncyclicChainSCFT(system, lattice, Float64(ds);
            mde=Polyorder.OSF, spacing=Float64(spacing),
            updater=field_updater, init=:none,
            rng=Random.Xoshiro(seed))
        Polyorder.initialize!(continued, previous_scft)
        continued
    end
    presolve_convergence = with_logger(NullLogger()) do
        Polyorder.solve!(scft, presolve_config)
    end
    cell_field_updater = Polyorder.Anderson(Polyorder.SD(LPC_SCFT_SD_STEP);
        m=LPC_SCFT_ANDERSON_DEPTH, warmup=LPC_SCFT_WARMUP,
        αw=LPC_SCFT_ANDERSON_WARMUP_STEP,
        droptol=Float64(anderson_droptol))
    cell_updater = Polyorder.VariableCell(
        Polyorder.BB(Float64(cell_bb_step_max);
            α=Float64(cell_bb_initial_step), type=Int(cell_bb_type)),
        cell_field_updater; block=Int(cell_block))
    cell_convergence, _scft_opt = with_logger(NullLogger()) do
        Polyorder.cell_solve!(scft, cell_updater, cell_config)
    end
    phi_a = _polyorder_density(scft, :A)
    phi_b = _polyorder_density(scft, :B)
    incompressibility_rms = sqrt(mean((phi_a .+ phi_b .- 1.0).^2))
    stress_norm = try
        norm(Polyorder.gradient_wrt_cell(scft))
    catch
        NaN
    end
    residual_norm = try
        Polyorder.residual(scft, cell_config)
    catch
        NaN
    end
    return (
        convergence=string(cell_convergence),
        presolve_convergence=string(presolve_convergence),
        cell_convergence=string(cell_convergence),
        phi_a=phi_a, phi_b=phi_b, free_energy=Polyorder.F(scft),
        initial_period_rg=Float64(initial_Lrg),
        period_rg=Polyorder.unitcell(scft).edges[1],
        grid_count=length(phi_a), cell_iterations=cell_updater.n,
        field_evaluations=isempty(cell_updater.evals) ? missing :
            cell_updater.evals[end],
        stress_norm=stress_norm, mean_phi_a=mean(phi_a),
        minimum_phi_a=minimum(phi_a), maximum_phi_a=maximum(phi_a),
        incompressibility_rms=incompressibility_rms, scft_model=scft,
        residual_norm=residual_norm,
        requested_spacing_rg=Float64(spacing),
        actual_spacing_rg=Polyorder.unitcell(scft).edges[1] / length(phi_a),
        contour_step=Float64(ds),
    )
end

function lpc_scft_cell_solve(; f::Real, chiN::Real,
        reference_period_rg::Real, spacing::Real, ds::Real, seed::Integer,
        presolve_max_iter::Integer, cell_max_iter::Integer,
        cell_block::Integer, stress_tol::Real, previous_scft=nothing,
        cell_bb_type::Integer=2,
        solver_residual_tol::Real=LPC_SCFT_RESIDUAL_TOL)
    errors = String[]
    if previous_scft !== nothing
        initial_period = Polyorder.unitcell(previous_scft).edges[1]
        continued = try
            lpc_polyorder_cell_solved_profile(; f=f, chiN=chiN,
                initial_Lrg=initial_period, spacing=spacing, ds=ds, seed=seed,
                presolve_max_iter=presolve_max_iter,
                cell_max_iter=cell_max_iter, cell_block=cell_block,
                stress_tol=stress_tol, previous_scft=previous_scft,
                cell_bb_type=cell_bb_type,
                solver_residual_tol=solver_residual_tol)
        catch err
            push!(errors, "continuation: " * sprint(showerror, err))
            nothing
        end
        if continued !== nothing &&
                _polyorder_success(continued.cell_convergence) &&
                continued.stress_norm <= Float64(stress_tol) &&
                continued.residual_norm < LPC_SCFT_RESIDUAL_TOL
            return continued
        end
    end

    attempts = NamedTuple[]
    for factor in (1.03, 1.0, 1.08, 1.15)
        result = try
            lpc_polyorder_cell_solved_profile(; f=f, chiN=chiN,
                initial_Lrg=reference_period_rg * factor,
                spacing=spacing, ds=ds, seed=seed,
                presolve_max_iter=presolve_max_iter,
                cell_max_iter=cell_max_iter, cell_block=cell_block,
                stress_tol=stress_tol, cell_bb_type=cell_bb_type,
                solver_residual_tol=solver_residual_tol)
        catch err
            push!(errors, "factor=$(factor): " * sprint(showerror, err))
            nothing
        end
        result === nothing || push!(attempts,
            (result=result, score=_scft_attempt_score(result)))
        result !== nothing && _polyorder_success(result.cell_convergence) &&
            result.stress_norm <= Float64(stress_tol) &&
            result.residual_norm < LPC_SCFT_RESIDUAL_TOL && return result
    end
    finite = [attempt for attempt in attempts if isfinite(attempt.score)]
    if isempty(finite)
        throw(ErrorException(
            "all SCFT cell solves failed; " * join(errors, " | ")))
    end
    return finite[argmin([attempt.score for attempt in finite])].result
end

function lpc_solve_model(model::AbstractString; f::Real, chiN::Real,
        nx::Integer, bvk2_nx::Integer, nquad::Integer, c2::Real,
        bvk2_oversample::Integer, bvk2_sensor_filter_ratio::Real,
        bvk2_energy_tolerance::Real,
        bvk2_force_tolerance::Real,
        bvk2_force_maxabs_tolerance::Real,
        bvk2_cell_stress_tolerance::Real,
        bvk2_log_period_tolerance::Real,
        max_iterations::Integer, max_period_iterations::Integer,
        burp_max_newton_iterations::Integer,
        ok_gradient_tolerance::Real,
        bootstrap_period_factors, lower_factor::Real, upper_factor::Real,
        bootstrap_window::Real, local_check_fraction::Real,
        progress_label::AbstractString, period_strategy::Symbol,
        opf_coefficients=nothing, ok_coefficients=nothing,
        initial_profile=nothing)
    initial_amplitudes = model == "liu2019_opf" ? (0.1,) :
        model == "ohta_kawasaki" ? (0.1, 0.2, 0.5) : (0.1, 0.2, 0.5, 1.0)
    common = (f=f, chiN=chiN, initial_amplitudes=initial_amplitudes,
        max_iterations=max_iterations,
        max_period_iterations=max_period_iterations,
        lower_factor=lower_factor, upper_factor=upper_factor,
        bootstrap_period_factors=bootstrap_period_factors,
        bootstrap_window=bootstrap_window,
        local_check_fraction=local_check_fraction,
        progress_label=progress_label)
    if model == "liu2019_opf"
        return minimize_diblock_opf_lamella_stress_free(; common...,
            nx=nx, mode_count=5, period_strategy=:stress_root,
            gradient_tolerance=1.0e-6, opf_coefficients=opf_coefficients)
    elseif model == "ohta_kawasaki"
        return minimize_diblock_opf_lamella_stress_free(; common...,
            nx=nx, mode_count=5, period_strategy=:stress_root,
            gradient_tolerance=Float64(ok_gradient_tolerance), check_domain=false,
            opf_coefficients=ok_coefficients)
    elseif model == "burp_ti"
        return minimize_diblock_burp_ti_lamella_stress_free_production(
            ; f=f, chiN=chiN, nx=nx, nquad=nquad,
            initial_profile=initial_profile,
            max_newton_iterations=burp_max_newton_iterations,
            bootstrap_period_factors=bootstrap_period_factors,
            lower_factor=lower_factor, upper_factor=upper_factor,
            local_check_fraction=min(Float64(local_check_fraction), 0.005),
            progress_label=progress_label)
    elseif model == "uneyama_doi"
        return minimize_diblock_uneyama_doi_lamella_stress_free(; common...,
            nx=nx, mode_count=5, relaxation_step=2.0e-2,
            tolerance=1.0e-8, gradient_tolerance=LPC_UD_RESIDUAL_TOL,
            initial_profile=initial_profile, warm_start_profiles=true,
            period_strategy=period_strategy)
    elseif model == "bvk2"
        return minimize_diblock_bvk2_ud_theta_lamella_stress_free(; common...,
            nx=bvk2_nx, mode_count=3, c2=c2,
            initial_profile=initial_profile,
            oversample=bvk2_oversample,
            sensor_filter_ratio=bvk2_sensor_filter_ratio,
            energy_tolerance=bvk2_energy_tolerance,
            force_tolerance=bvk2_force_tolerance,
            force_maxabs_tolerance=bvk2_force_maxabs_tolerance,
            cell_stress_tolerance=bvk2_cell_stress_tolerance,
            log_period_tolerance=bvk2_log_period_tolerance)
    end
    throw(ArgumentError("unsupported model $(model)"))
end

function lpc_profile_statistics(reference::AbstractVector{<:Real},
        aligned::AbstractVector{<:Real})
    residual_sum = sum(abs2, aligned .- reference)
    centered_sum = sum(abs2, reference .- mean(reference))
    return (
        rms=sqrt(residual_sum / length(reference)),
        r2=centered_sum > eps(Float64) ?
            1.0 - residual_sum / centered_sum : NaN,
        correlation=cor(reference, aligned),
    )
end

function lpc_scft_summary_row(case, scft, case_id::AbstractString,
        profile_path::AbstractString)
    field_pass = _polyorder_success(scft.cell_convergence) &&
        all(isfinite, scft.phi_a) && all(isfinite, scft.phi_b) &&
        isfinite(scft.residual_norm) &&
        scft.residual_norm < LPC_SCFT_RESIDUAL_TOL
    cell_pass = isfinite(scft.stress_norm) && scft.stress_norm <= 1.0e-5
    resolution_pass = isfinite(scft.requested_spacing_rg) &&
        scft.requested_spacing_rg > 0.0 &&
        0.8 * scft.requested_spacing_rg <= scft.actual_spacing_rg <=
            1.2 * scft.requested_spacing_rg
    composition_pass = abs(scft.mean_phi_a - case.f) <= 1.0e-6 &&
        scft.incompressibility_rms <= 1.0e-6
    morphology_pass = scft.maximum_phi_a - scft.minimum_phi_a > 1.0e-4
    accepted = field_pass && cell_pass && resolution_pass && composition_pass &&
        morphology_pass
    return (
        schema=LPC_SCHEMA, claim_kind="stress_free_lamellar_observable",
        case_id=case_id, f=case.f, chiN=case.chiN,
        chiN_spinodal=case.chiN_spinodal, model="scft",
        model_label=LPC_MODEL_LABELS["scft"], calibration_role="reference",
        period_rg=scft.period_rg, scft_period_rg=scft.period_rg,
        signed_period_error=0.0, abs_period_error=0.0,
        profile_rms=0.0, profile_r2=1.0, profile_correlation=1.0,
        alignment_shift_fraction=0.0, alignment_swapped=false,
        period_factor=missing, period_optimizer="Polyorder.cell_solve!",
        period_iterations=scft.cell_iterations,
        period_local_minimum_check_pass=cell_pass,
        period_boundary_limited=false,
        field_gate_pass=field_pass, cell_gate_pass=cell_pass,
        resolution_gate_pass=resolution_pass,
        composition_gate_pass=composition_pass,
        morphology_gate_pass=morphology_pass,
        minimum_phi_a=scft.minimum_phi_a, maximum_phi_a=scft.maximum_phi_a,
        mean_phi_a=scft.mean_phi_a, energy=scft.free_energy,
        reduced_field_residual_norm=missing,
        reduced_field_residual_maxabs=missing,
        reduced_field_protocol="not_applicable",
        stress_norm=scft.stress_norm,
        scft_convergence=scft.cell_convergence,
        scft_residual_norm=scft.residual_norm,
        scft_incompressibility_rms=scft.incompressibility_rms,
        scft_requested_spacing_rg=scft.requested_spacing_rg,
        scft_actual_spacing_rg=scft.actual_spacing_rg,
        scft_contour_step=scft.contour_step,
        scft_field_updater="Anderson(SD(0.1); m=10, warmup=100, αw=0.2)",
        scft_cell_updater="VariableCell(BB(1.0), field_updater)",
        grid_count=scft.grid_count, quadrature_order=missing,
        oversample_factor=missing, sensor_filter_ratio=missing,
        discretization_schema="Polyorder pseudospectral SCFT",
        field_iteration_cap=missing, branch_dominant_mode=missing, c2=missing,
        opf_c2=missing, opf_c3=missing, opf_c4=missing, opf_c5=missing,
        opf_c6=missing, opf_mapping_residual=missing,
        native_profile=profile_path,
        status=accepted ? "accepted" : "provisional",
        status_reason=accepted ? "reference_gates_pass" :
            "one_or_more_reference_gates_failed",
        error_message="",
    )
end

function lpc_model_summary_row(case, case_id::AbstractString,
        model::AbstractString, optimized, scft, alignment, stats,
        profile_path::AbstractString, c2::Real; opf_coefficients=nothing,
        ok_coefficients=nothing, burp_nquad::Integer=64,
        burp_max_newton_iterations::Integer=256,
        bvk2_oversample::Integer=LPC_BVK2_OVERSAMPLE,
        bvk2_sensor_filter_ratio::Real=LPC_BVK2_SENSOR_FILTER_RATIO,
        bvk2_residual_tolerance::Real=LPC_BVK2_RESIDUAL_TOL,
        bvk2_maximum_residual_tolerance::Real=LPC_BVK2_MAX_RESIDUAL_TOL,
        bvk2_stress_tolerance::Real=LPC_BVK2_STRESS_TOL)
    result = optimized.result
    polynomial_phase_field = model in ("liu2019_opf", "ohta_kawasaki")
    residual_tolerance = model == "uneyama_doi" ? LPC_UD_RESIDUAL_TOL :
        model == "bvk2" ? Float64(bvk2_residual_tolerance) :
        polynomial_phase_field ? 1.0e-5 : 5.0e-4
    maximum_residual_tolerance = model == "uneyama_doi" ?
        LPC_UD_MAX_RESIDUAL_TOL :
        model == "bvk2" ? Float64(bvk2_maximum_residual_tolerance) :
        polynomial_phase_field ? 1.0e-4 : 4.0e-3
    residual_gate = model in ("liu2019_opf", "ohta_kawasaki", "burp_ti",
        "uneyama_doi", "bvk2") ?
        isfinite(result.projected_force_norm) &&
        result.projected_force_norm <= residual_tolerance &&
        isfinite(result.projected_force_maxabs) &&
        result.projected_force_maxabs <= maximum_residual_tolerance : true
    # For polynomial phase-field and UD models, the campaign contract is stated
    # in the physical-space projected force. An optimizer's internal flag must
    # not veto a state that already satisfies the declared physical residual.
    optimizer_gate = polynomial_phase_field || model == "uneyama_doi" ?
        true : result.converged
    field_pass = optimizer_gate && residual_gate &&
        all(isfinite, result.phi_a) && isfinite(result.energy)
    composition_pass = abs(result.mean_phi - case.f) <= 1.0e-8
    branch_dominant_mode = model == "burp_ti" ?
        Int(optimized.dominant_mode) :
        model == "uneyama_doi" ?
        burp_ti_lamellar_dominant_mode(result.phi_a, case.f) : missing
    morphology_pass = result.maximum_phi - result.minimum_phi > 1.0e-4 &&
        (model == "burp_ti" ? branch_dominant_mode == 1 :
         model == "uneyama_doi" ?
            DFMMonteCarlo._diblock_one_period_lamellar_profile(result.phi_a) :
            true)
    signed_error = optimized.period_rg / scft.period_rg - 1.0
    active_coefficients = model == "liu2019_opf" ? opf_coefficients :
        model == "ohta_kawasaki" ? ok_coefficients : nothing
    stress_norm = polynomial_phase_field ? abs(diblock_opf_isotropic_stress_1d(
        result.phi_a; f=case.f, chiN=case.chiN, L=result.L,
        check_domain=false, opf_coefficients=active_coefficients)) :
        model == "burp_ti" ? abs(Float64(optimized.cell_stress)) :
        model == "bvk2" ? abs(
            diblock_bvk2_ud_theta_cell_scale_gradient_density_nd(
                asin.(sqrt.(clamp.(result.phi_a, 0.0, 1.0)));
                f=case.f, chiN=case.chiN, lengths=(result.L,), c2=c2,
                oversample=bvk2_oversample,
                sensor_filter_ratio=bvk2_sensor_filter_ratio)) : missing
    cell_pass = optimized.local_minimum_check_pass &&
        !optimized.boundary_limited &&
        (!(polynomial_phase_field || model in ("burp_ti", "bvk2")) ||
            stress_norm <= (model == "bvk2" ?
                Float64(bvk2_stress_tolerance) : 1.0e-4))
    failure_reason = if model == "uneyama_doi" && field_pass &&
            composition_pass && morphology_pass && !cell_pass &&
            !optimized.local_minimum_check_pass
        optimized.boundary_limited ?
            "no_bracketed_primitive_period_minimum" :
            "no_two_sided_primitive_period_minimum"
    else
        "one_or_more_model_gates_failed"
    end
    accepted = field_pass && cell_pass && composition_pass && morphology_pass
    return (
        schema=LPC_SCHEMA, claim_kind="stress_free_lamellar_observable",
        case_id=case_id, f=case.f, chiN=case.chiN,
        chiN_spinodal=case.chiN_spinodal, model=model,
        model_label=LPC_MODEL_LABELS[model],
        calibration_role=lpc_model_calibration_role(model, case.f, case.chiN),
        period_rg=optimized.period_rg, scft_period_rg=scft.period_rg,
        signed_period_error=signed_error, abs_period_error=abs(signed_error),
        profile_rms=stats.rms, profile_r2=stats.r2,
        profile_correlation=stats.correlation,
        alignment_shift_fraction=alignment.shift / length(alignment.values),
        alignment_swapped=alignment.swapped,
        period_factor=optimized.period_factor,
        period_optimizer=optimized.optimizer,
        period_iterations=optimized.iterations,
        period_local_minimum_check_pass=optimized.local_minimum_check_pass,
        period_boundary_limited=optimized.boundary_limited,
        field_gate_pass=field_pass, cell_gate_pass=cell_pass,
        resolution_gate_pass=true,
        composition_gate_pass=composition_pass,
        morphology_gate_pass=morphology_pass,
        minimum_phi_a=result.minimum_phi, maximum_phi_a=result.maximum_phi,
        mean_phi_a=result.mean_phi, energy=result.energy,
        reduced_field_residual_norm=result.projected_force_norm,
        reduced_field_residual_maxabs=result.projected_force_maxabs,
        reduced_field_protocol=model == "bvk2" ? LPC_BVK2_PROTOCOL :
            polynomial_phase_field ?
            "mean-constrained Fourier nonlinear conjugate gradient" :
            model == "burp_ti" ?
            "bounded-KKT Newton-GMRES" :
            model == "uneyama_doi" ?
            LPC_UD_PROTOCOL :
            "fixed-composition tangent-projected BB descent",
        stress_norm=stress_norm, scft_convergence=scft.cell_convergence,
        scft_residual_norm=scft.residual_norm,
        scft_incompressibility_rms=scft.incompressibility_rms,
        scft_requested_spacing_rg=scft.requested_spacing_rg,
        scft_actual_spacing_rg=scft.actual_spacing_rg,
        scft_contour_step=scft.contour_step,
        scft_field_updater="Anderson(SD(0.1); m=10, warmup=100, αw=0.2)",
        scft_cell_updater="VariableCell(BB(1.0), field_updater)",
        grid_count=length(result.phi_a),
        quadrature_order=model == "burp_ti" ? Int(burp_nquad) : missing,
        oversample_factor=model == "bvk2" ? Int(bvk2_oversample) : missing,
        sensor_filter_ratio=model == "bvk2" ?
            Float64(bvk2_sensor_filter_ratio) : missing,
        discretization_schema=model == "bvk2" ?
            DFMMonteCarlo._bvk2_ud_theta_discretization_schema(
                bvk2_sensor_filter_ratio) : missing,
        field_iteration_cap=model == "burp_ti" ?
            Int(burp_max_newton_iterations) : missing,
        branch_dominant_mode=branch_dominant_mode,
        c2=model == "bvk2" ? c2 : missing,
        opf_c2=polynomial_phase_field ? active_coefficients.c2 : missing,
        opf_c3=polynomial_phase_field ? active_coefficients.c3 : missing,
        opf_c4=polynomial_phase_field ? active_coefficients.c4 : missing,
        opf_c5=polynomial_phase_field ? active_coefficients.c5 : missing,
        opf_c6=polynomial_phase_field ? active_coefficients.c6 : missing,
        opf_mapping_residual=model == "liu2019_opf" ?
            opf_coefficients.force_stress_residual_norm :
            missing,
        native_profile=profile_path,
        status=accepted ? "accepted" : "provisional",
        status_reason=accepted ? "model_gates_pass" :
            failure_reason,
        error_message="",
    )
end

function lpc_failed_model_row(case, case_id::AbstractString,
        model::AbstractString, scft, message::AbstractString, c2::Real;
        bvk2_oversample::Integer=LPC_BVK2_OVERSAMPLE,
        bvk2_sensor_filter_ratio::Real=LPC_BVK2_SENSOR_FILTER_RATIO)
    return (
        schema=LPC_SCHEMA, claim_kind="stress_free_lamellar_observable",
        case_id=case_id, f=case.f, chiN=case.chiN,
        chiN_spinodal=case.chiN_spinodal, model=model,
        model_label=LPC_MODEL_LABELS[model],
        calibration_role=lpc_model_calibration_role(model, case.f, case.chiN),
        period_rg=NaN, scft_period_rg=scft.period_rg,
        signed_period_error=NaN, abs_period_error=NaN,
        profile_rms=NaN, profile_r2=NaN, profile_correlation=NaN,
        alignment_shift_fraction=NaN, alignment_swapped=false,
        period_factor=NaN, period_optimizer="failed",
        period_iterations=missing,
        period_local_minimum_check_pass=false,
        period_boundary_limited=false,
        field_gate_pass=false, cell_gate_pass=false,
        resolution_gate_pass=true,
        composition_gate_pass=false, morphology_gate_pass=false,
        minimum_phi_a=NaN, maximum_phi_a=NaN, mean_phi_a=NaN,
        energy=NaN, reduced_field_residual_norm=NaN,
        reduced_field_residual_maxabs=NaN,
        reduced_field_protocol=model == "bvk2" ?
            LPC_BVK2_PROTOCOL :
            model == "burp_ti" ? "bounded-KKT Newton-GMRES" :
            model == "uneyama_doi" ? LPC_UD_PROTOCOL :
            "fixed-composition tangent-projected BB descent",
        stress_norm=missing,
        scft_convergence=scft.cell_convergence,
        scft_residual_norm=scft.residual_norm,
        scft_incompressibility_rms=scft.incompressibility_rms,
        scft_requested_spacing_rg=scft.requested_spacing_rg,
        scft_actual_spacing_rg=scft.actual_spacing_rg,
        scft_contour_step=scft.contour_step,
        scft_field_updater="Anderson(SD(0.1); m=10, warmup=100, αw=0.2)",
        scft_cell_updater="VariableCell(BB(1.0), field_updater)",
        grid_count=missing, quadrature_order=missing,
        oversample_factor=model == "bvk2" ? Int(bvk2_oversample) : missing,
        sensor_filter_ratio=model == "bvk2" ?
            Float64(bvk2_sensor_filter_ratio) : missing,
        discretization_schema=model == "bvk2" ?
            DFMMonteCarlo._bvk2_ud_theta_discretization_schema(
                bvk2_sensor_filter_ratio) : missing,
        field_iteration_cap=missing, branch_dominant_mode=missing,
        c2=model == "bvk2" ? c2 : missing,
        opf_c2=missing, opf_c3=missing, opf_c4=missing, opf_c5=missing,
        opf_c6=missing, opf_mapping_residual=missing,
        native_profile="", status="rejected",
        status_reason="solver_failure", error_message=String(message),
    )
end

function lpc_push_common_profile!(rows, native_rows, case, case_id, model,
        period_rg, raw_values, aligned_values; phi_b_raw=nothing,
        sample_count::Integer=256)
    raw = _resample(raw_values, sample_count)
    aligned = Float64.(aligned_values)
    raw_b = phi_b_raw === nothing ? 1.0 .- raw :
        _resample(phi_b_raw, sample_count)
    aligned_b = 1.0 .- aligned
    for index in 1:sample_count
        s = (index - 1) / sample_count
        push!(rows, (
            schema=LPC_SCHEMA, case_id=case_id, f=case.f, chiN=case.chiN,
            model=model, model_label=LPC_MODEL_LABELS[model],
            period_rg=period_rg, s=s, x_rg=s * period_rg,
            phi_a=aligned[index], phi_b=aligned_b[index],
            phi_a_raw=raw[index], phi_b_raw=raw_b[index],
        ))
    end
    native_b = phi_b_raw === nothing ? 1.0 .- Float64.(raw_values) :
        Float64.(phi_b_raw)
    for index in eachindex(raw_values)
        push!(native_rows, (
            schema=LPC_SCHEMA, case_id=case_id, f=case.f, chiN=case.chiN,
            model=model, model_label=LPC_MODEL_LABELS[model],
            period_rg=period_rg, native_index=index,
            native_grid_count=length(raw_values),
            s=(index - 1) / length(raw_values),
            phi_a=Float64(raw_values[index]), phi_b=native_b[index],
        ))
    end
end

function lpc_write_case(case_dir::AbstractString, summary_rows,
        profile_rows, native_rows)
    mkpath(case_dir)
    function atomic_csv_write(path::AbstractString, table)
        temporary_path, stream = mktemp(dirname(path))
        close(stream)
        try
            CSV.write(temporary_path, table)
            mv(temporary_path, path; force=true)
        finally
            isfile(temporary_path) && rm(temporary_path; force=true)
        end
        return path
    end
    function merged_table(filename::AbstractString, rows)
        fresh = DataFrame(rows)
        path = joinpath(case_dir, filename)
        isfile(path) || return fresh
        existing = try
            CSV.read(path, DataFrame)
        catch
            return fresh
        end
        (:schema in propertynames(existing) && :model in propertynames(existing)) ||
            return fresh
        replaced = Set(String.(fresh.model))
        keep = [String(row.schema) == LPC_SCHEMA &&
            !(String(row.model) in replaced) for row in eachrow(existing)]
        preserved = existing[keep, :]
        isempty(preserved) && return fresh
        return vcat(fresh, preserved; cols=:union)
    end
    summary_path = atomic_csv_write(joinpath(case_dir, "summary.csv"),
        merged_table("summary.csv", summary_rows))
    profiles_path = atomic_csv_write(joinpath(case_dir, "profiles.csv"),
        merged_table("profiles.csv", profile_rows))
    native_path = atomic_csv_write(joinpath(case_dir, "native_profiles.csv"),
        merged_table("native_profiles.csv", native_rows))
    return summary_path, profiles_path, native_path
end

function lpc_merge_cases(outdir::AbstractString)
    cases_dir = joinpath(outdir, "cases")
    isdir(cases_dir) || return nothing
    directories = sort([dir for dir in readdir(cases_dir; join=true)
        if isfile(joinpath(dir, "summary.csv"))])
    isempty(directories) && return nothing
    summaries = vcat([CSV.read(joinpath(dir, "summary.csv"), DataFrame)
        for dir in directories]...; cols=:union)
    profiles = vcat([CSV.read(joinpath(dir, "profiles.csv"), DataFrame)
        for dir in directories]...; cols=:union)
    native = vcat([CSV.read(joinpath(dir, "native_profiles.csv"), DataFrame)
        for dir in directories]...; cols=:union)
    order = Dict(model => index for (index, model) in enumerate(LPC_MODEL_ORDER))
    sort!(summaries, [:f, :chiN, :model]; rev=[true, false, false],
        by=[identity, identity, value -> get(order, String(value), typemax(Int))])
    sort!(profiles, [:f, :chiN, :model, :s]; rev=[true, false, false, false])
    sort!(native, [:f, :chiN, :model, :native_index];
        rev=[true, false, false, false])
    CSV.write(joinpath(outdir, "summary.csv"), summaries)
    CSV.write(joinpath(outdir, "profiles.csv"), profiles)
    CSV.write(joinpath(outdir, "native_profiles.csv"), native)
    return (summaries=summaries, profiles=profiles, native=native)
end

function lpc_write_contract(outdir::AbstractString, settings, cases)
    extension = maximum(case.chiN for case in cases) > 35.0
    path = joinpath(outdir, "README.md")
    open(path, "w") do io
        println(io, extension ?
            "# Strong-segregation lamellar benchmark extension" :
            "# Liu et al. stress-free lamellar period comparison")
        println(io)
        println(io, "Review and reuse runbook: " *
            "`docs/liu2019_stress_free_period_map.md`.")
        println(io)
        println(io, "This dataset independently recomputes the domain-spacing comparison " *
            "of Liu et al. (Macromolecules 2019) using their per-state Eq. 18 " *
            "force-and-stress OPF mapping and " *
            "the Ohta-Kawasaki model, and adds Uneyama-Doi, BURP-TI, and BVK2 " *
            "(frozen c2 = $(settings.c2)). Published Figure 2 traces, when shown, are " *
            "stored separately as digitized literature data and are not conflated with " *
            "these calculations.")
        println(io)
        println(io, "## Campaign contract")
        println(io)
        println(io, "- Objective: stress-free lamellar periods and independently relaxed " *
            (extension ?
                "density profiles for a separately reported strong-segregation extension." :
                "density profiles over the published OPF domain."))
        println(io, "- Claim kind: stress-free lamellar observable at one (fA, chiN) state.")
        println(io, "- Accepted outputs: `summary.csv`, `profiles.csv`, " *
            "`native_profiles.csv`, and the case directories under `cases/`.")
        println(io, "- SCFT acceptance: successful Polyorder convergence, residual norm " *
            "below 1e-6, cell-stress norm at most 1e-5, incompressibility and " *
            "mean-composition errors at most 1e-6, " *
            "and a nonuniform lamellar profile.")
        println(io, "- OPF acceptance: mean-constrained Fourier nonlinear conjugate " *
            "gradient with field-residual RMS at most 1e-5 and maximum absolute " *
            "component at most 1e-4, direct zero-stress cell selection, and a " *
            "two-sided local-energy check. BURP-TI acceptance: bounded-KKT " *
            "Newton-GMRES residual norm at most 5e-4 and maximum component at most " *
            "4e-3, analytic cell-stress magnitude at most 1e-4, and a two-sided " *
            "relaxed-energy check. UD acceptance requires projected-force RMS " *
            "at most $(LPC_UD_RESIDUAL_TOL), maximum component at most " *
            "$(LPC_UD_MAX_RESIDUAL_TOL), fused-spectral constrained L-BFGS with " *
            "profile continuation, one contiguous A-rich domain per periodic " *
            "cell, " *
            "and a two-sided relaxed-energy check. BVK2 acceptance requires the " *
            "oversampled filtered UD-theta functional, projected-force RMS at most " *
            "$(settings.bvk2_force_tolerance), maximum component at most " *
            "$(settings.bvk2_force_maxabs_tolerance), analytic cell-stress magnitude " *
            "at most $(settings.bvk2_cell_stress_tolerance), and a two-sided " *
            "relaxed-energy check. Other reduced models use their " *
            "documented field solvers and tolerances. All models require mean-composition " *
            "error at most 1e-8, finite energy, and a nonuniform profile.")
        println(io, "- Status vocabulary: accepted, provisional, or rejected. A zero exit " *
            "status alone does not promote a row.")
        requested_chi = sort(unique(case.chiN for case in cases))
        println(io, "- Parameter domain: requested fA values and chiN = " *
            join(requested_chi, ", ") * ". " *
            (extension ?
                "These rows are an extension and do not alter the chiN <= 35 Liu benchmark cohort." :
                "All rows lie within the published comparison domain."))
        println(io, "- fA = 0.15 is absent because its reported spinodal, chiN = 38.038, " *
            "lies outside the plotted upper limit.")
        println(io, "- Retry budget: two numerically distinct retries per rejected or " *
            "provisional state; accepted checkpoints are reused.")
        println(io, "- Stop condition: every requested state has one accepted SCFT row and " *
            "one accepted row for each requested reduced model, followed by a successful " *
            "merge and figure render.")
        println(io)
        println(io, "## Numerical settings")
        println(io)
        println(io, "- Requested states: $(length(cases))")
        println(io, "- Reduced-model grids: nx = $(settings.nx); BVK2 nx = $(settings.bvk2_nx)")
        println(io, "- Requested aligned-profile grid for newly computed rows: " *
            "$(settings.sample_count) points")
        println(io, "- Polyorder spacing = $(settings.spacing) and ds = $(settings.ds) " *
            "by default; accepted reused references retain their state-specific " *
            "audited settings. The refined fA = 0.20, chiN = 25--27 references " *
            "use spacing = 0.05 and ds = 0.005. The stress tolerance is " *
            "$(settings.scft_stress_tol) throughout.")
        println(io, "- Polyorder field update = Anderson(10) accelerated SD(0.1), " *
            "with 100 warmup iterations and alpha_w = 0.2; cell optimization uses " *
            "cell_solve! with VariableCell and BB(1.0) cell relaxation.")
        println(io, "- BURP-TI quadrature order = $(settings.nquad); bounded-KKT " *
            "Newton cap = $(settings.burp_max_newton_iterations)")
        println(io, "- BVK2 discretization = oversampled UD-theta functional with " *
            "oversample factor $(settings.bvk2_oversample), one global physical " *
            "sensor cutoff kc/kstar = $(settings.bvk2_sensor_filter_ratio), " *
            "cosine-Newton field polish, and filter-consistent analytic cell stress")
        println(io, "- OPF cell strategy = direct Eq. S2 stress-root bisection")
        println(io, "- UD cell strategy = fused-spectral mean-constrained L-BFGS with " *
            "accepted-profile continuation in chiN and period, competing seed " *
            "branches, and a two-sided relaxed-energy audit")
        println(io, "- Default remaining reduced-model period strategy = " *
            "$(settings.period_strategy)")
        println(io, "- Models: $(join(settings.models, ", "))")
        println(io)
        println(io, "Profiles are stored twice. `profiles.csv` contains phase-aligned " *
            "profiles on a normalized periodic coordinate for RMS comparison; " *
            "`native_profiles.csv` preserves the solver grids before resampling or alignment.")
    end
    return path
end

function lpc_main()
    if any(argument -> argument in ("--help", "-h"), ARGS)
        print(LPC_USAGE)
        return nothing
    end
    outdir = normpath(_arg_value("--outdir", joinpath(MONTE_CARLO_PROJECT,
        "results", "liu2019_stress_free_period_map")))
    compositions = [parse(Float64, value) for value in
        lpc_parse_strings(_arg_value("--compositions", "0.20,0.25,0.30,0.35,0.40,0.45,0.50"))]
    chi_order = Symbol(lowercase(_arg_value("--chi-order", "ascending")))
    chi_order in (:ascending, :descending) || throw(ArgumentError(
        "--chi-order must be ascending or descending"))
    explicit_cases = strip(_arg_value("--cases", ""))
    max_chiN = parse(Float64, _arg_value("--max-chiN", "35"))
    isfinite(max_chiN) && max_chiN > 0.0 || throw(ArgumentError(
        "--max-chiN must be finite and positive"))
    cases = isempty(explicit_cases) ?
        lpc_default_cases(compositions; chi_order=chi_order,
            max_chiN=max_chiN) :
        lpc_explicit_cases(explicit_cases; chi_order=chi_order,
            max_chiN=max_chiN)
    models = lpc_parse_strings(_arg_value("--models",
        "liu2019_opf,ohta_kawasaki,uneyama_doi,burp_ti,bvk2"))
    all(model -> model in LPC_MODEL_ORDER[2:end], models) ||
        throw(ArgumentError("unsupported model selection"))
    allow_ok_extrapolation = _parse_bool(
        _arg_value("--allow-ok-extrapolation", "false"))
    if "ohta_kawasaki" in models &&
            any(case -> case.chiN > 35.0, cases) &&
            !allow_ok_extrapolation
        throw(ArgumentError(
            "OK states above chiN=35 require --allow-ok-extrapolation=true"))
    end
    nx = parse(Int, _arg_value("--nx", "128"))
    bvk2_nx = parse(Int, _arg_value("--bvk2-nx", "128"))
    nquad = parse(Int, _arg_value("--nquad", "64"))
    burp_max_newton_iterations = parse(Int,
        _arg_value("--burp-max-newton-iterations", "256"))
    ok_gradient_tolerance = parse(Float64,
        _arg_value("--ok-gradient-tol", "8e-6"))
    0.0 < ok_gradient_tolerance <= 8.0e-6 || throw(ArgumentError(
        "--ok-gradient-tol must lie in (0, 8e-6]"))
    c2 = parse(Float64, _arg_value("--bvk2-c2", "0.16"))
    bvk2_oversample = parse(Int, _arg_value(
        "--bvk2-oversample", string(LPC_BVK2_OVERSAMPLE)))
    bvk2_sensor_filter_ratio = parse(Float64,
        _arg_value("--bvk2-sensor-filter-ratio",
            string(LPC_BVK2_SENSOR_FILTER_RATIO)))
    bvk2_energy_tolerance = parse(Float64,
        _arg_value("--bvk2-energy-tol", "0.0"))
    bvk2_force_tolerance = parse(Float64,
        _arg_value("--bvk2-force-tol", string(LPC_BVK2_RESIDUAL_TOL)))
    bvk2_force_maxabs_tolerance = parse(Float64,
        _arg_value("--bvk2-force-maxabs-tol",
            string(LPC_BVK2_MAX_RESIDUAL_TOL)))
    bvk2_cell_stress_tolerance = parse(Float64,
        _arg_value("--bvk2-cell-stress-tol", string(LPC_BVK2_STRESS_TOL)))
    bvk2_log_period_tolerance = parse(Float64,
        _arg_value("--bvk2-log-period-tol", "1e-7"))
    bvk2_oversample >= 1 || throw(ArgumentError(
        "--bvk2-oversample must be at least 1"))
    isfinite(bvk2_sensor_filter_ratio) && bvk2_sensor_filter_ratio > 0.0 ||
        throw(ArgumentError(
            "--bvk2-sensor-filter-ratio must be finite and positive"))
    bvk2_energy_tolerance >= 0.0 || throw(ArgumentError(
        "--bvk2-energy-tol must be nonnegative"))
    bvk2_force_tolerance > 0.0 || throw(ArgumentError(
        "--bvk2-force-tol must be positive"))
    bvk2_force_maxabs_tolerance > 0.0 || throw(ArgumentError(
        "--bvk2-force-maxabs-tol must be positive"))
    bvk2_cell_stress_tolerance > 0.0 || throw(ArgumentError(
        "--bvk2-cell-stress-tol must be positive"))
    bvk2_log_period_tolerance > 0.0 || throw(ArgumentError(
        "--bvk2-log-period-tol must be positive"))
    max_iterations = parse(Int, _arg_value("--model-max-iterations", "10000"))
    max_period_iterations = parse(Int,
        _arg_value("--model-max-period-iterations", "24"))
    period_strategy = Symbol(replace(lowercase(
        _arg_value("--period-strategy", "bootstrap_local")), "-" => "_"))
    period_strategy in (:bootstrap_local, :brent) || throw(ArgumentError(
        "--period-strategy must be bootstrap_local or brent"))
    model_seed_policy = Symbol(replace(lowercase(
        _arg_value("--model-seed-policy", "history")), "-" => "_"))
    model_seed_policy in (:history, :reference) || throw(ArgumentError(
        "--model-seed-policy must be history or reference"))
    lower_factor = parse(Float64, _arg_value("--lower-factor", "0.75"))
    upper_factor = parse(Float64, _arg_value("--upper-factor", "1.8"))
    bootstrap_window = parse(Float64, _arg_value("--bootstrap-window", "0.18"))
    local_check_fraction = parse(Float64,
        _arg_value("--local-check-fraction", "0.01"))
    spacing = parse(Float64, _arg_value("--spacing", "0.1"))
    ds = parse(Float64, _arg_value("--ds", "0.01"))
    scft_presolve_max_iter = parse(Int,
        _arg_value("--scft-presolve-max-iterations", "3000"))
    scft_cell_max_iter = parse(Int,
        _arg_value("--scft-cell-max-iterations", "8000"))
    scft_cell_block = parse(Int, _arg_value("--scft-cell-block", "5"))
    scft_solver_residual_tol = parse(Float64,
        _arg_value("--scft-solver-residual-tol", "1e-6"))
    0.0 < scft_solver_residual_tol <= LPC_SCFT_RESIDUAL_TOL ||
        throw(ArgumentError("--scft-solver-residual-tol must lie in " *
            "(0, $(LPC_SCFT_RESIDUAL_TOL)]"))
    scft_cell_bb_type = parse(Int,
        _arg_value("--scft-cell-bb-type", "2"))
    scft_cell_bb_type in (1, 2, 3) || throw(ArgumentError(
        "--scft-cell-bb-type must be 1, 2, or 3"))
    scft_stress_tol = parse(Float64, _arg_value("--scft-stress-tol", "1e-5"))
    sample_count = parse(Int, _arg_value("--sample-count", "256"))
    seed = parse(Int, _arg_value("--seed", "3083"))
    resume = _parse_bool(_arg_value("--resume", "true"))
    reuse_existing_scft = _parse_bool(
        _arg_value("--reuse-existing-scft", "false"))
    assemble_only = _parse_bool(_arg_value("--assemble-only", "false"))
    merge_outputs = _parse_bool(_arg_value("--merge-outputs", "true"))
    write_contract = _parse_bool(_arg_value("--write-contract", "true"))
    mkpath(joinpath(outdir, "cases"))
    settings = (; nx, bvk2_nx, nquad, burp_max_newton_iterations,
        ok_gradient_tolerance, c2,
        bvk2_oversample, bvk2_sensor_filter_ratio,
        bvk2_energy_tolerance,
        bvk2_force_tolerance, bvk2_force_maxabs_tolerance,
        bvk2_cell_stress_tolerance, bvk2_log_period_tolerance,
        spacing, ds, scft_stress_tol,
        sample_count, models, period_strategy, allow_ok_extrapolation)
    write_contract && lpc_write_contract(outdir, settings, cases)
    assemble_only && return lpc_merge_cases(outdir)

    history = model_seed_policy == :history ?
        lpc_model_history(outdir, models; burp_nquad=nquad,
            burp_max_newton_iterations=burp_max_newton_iterations,
            bvk2_oversample=bvk2_oversample,
            bvk2_sensor_filter_ratio=bvk2_sensor_filter_ratio,
            bvk2_residual_tol=bvk2_force_tolerance,
            bvk2_maximum_residual_tol=bvk2_force_maxabs_tolerance,
            bvk2_stress_tol=bvk2_cell_stress_tolerance) :
        Dict{Tuple{String,Float64},Vector{NamedTuple}}()
    profile_history = lpc_model_profile_history(outdir, models, nx;
        burp_nquad=nquad,
        burp_max_newton_iterations=burp_max_newton_iterations)
    scft_period_history = lpc_scft_period_history(outdir)
    previous_scft = Dict{Float64,Any}()
    for (case_index, case) in enumerate(cases)
        case_id = _case_id(case.f, case.chiN)
        case_dir = lpc_case_dir(outdir, case)
        if resume && lpc_existing_case_accepted(case_dir, models;
                burp_nquad=nquad,
                burp_max_newton_iterations=burp_max_newton_iterations,
                bvk2_oversample=bvk2_oversample,
                bvk2_sensor_filter_ratio=bvk2_sensor_filter_ratio,
                bvk2_residual_tol=bvk2_force_tolerance,
                bvk2_maximum_residual_tol=bvk2_force_maxabs_tolerance,
                bvk2_stress_tol=bvk2_cell_stress_tolerance)
            println("[", case_id, "] accepted checkpoint reused")
            continue
        end
        scft = reuse_existing_scft ? lpc_load_existing_scft(case_dir) : nothing
        if scft === nothing
            println("[", case_id, "] SCFT start")
            flush(stdout)
            kernel_reference_period_rg = 2.0 * pi /
                sqrt(DFMMonteCarlo._diblock_kernel_minimum_k2(; f=case.f)) * sqrt(6.0)
            estimated_period_rg = lpc_scft_initial_period(
                kernel_reference_period_rg, case.chiN, case.chiN_spinodal)
            reference_period_rg = lpc_nearest_scft_period(scft_period_history,
                case.f, case.chiN, estimated_period_rg)
            scft = try
                lpc_scft_cell_solve(; f=case.f, chiN=case.chiN,
                    reference_period_rg=reference_period_rg,
                    spacing=spacing, ds=ds,
                    seed=seed + round(Int, 100 * case.f),
                    presolve_max_iter=scft_presolve_max_iter,
                    cell_max_iter=scft_cell_max_iter, cell_block=scft_cell_block,
                    stress_tol=scft_stress_tol,
                    cell_bb_type=scft_cell_bb_type,
                    solver_residual_tol=scft_solver_residual_tol,
                    previous_scft=get(previous_scft, Float64(case.f), nothing))
            catch err
                println("[", case_id, "] SCFT failed: ", sprint(showerror, err))
                flush(stdout)
                continue
            end
        else
            println("[", case_id, "] SCFT accepted checkpoint loaded")
        end
        println(@sprintf(
            "[%s] SCFT done L/Rg=%.8g dx/Rg=%.4g residual=%.3g stress=%.3g status=%s",
            case_id, scft.period_rg, scft.actual_spacing_rg,
            scft.residual_norm, scft.stress_norm, scft.cell_convergence))
        flush(stdout)
        if _polyorder_success(scft.cell_convergence) &&
                scft.stress_norm <= scft_stress_tol &&
                scft.residual_norm < LPC_SCFT_RESIDUAL_TOL
            if scft.scft_model !== nothing
                previous_scft[Float64(case.f)] = scft.scft_model
            end
            push!(get!(scft_period_history, Float64(case.f), NamedTuple[]),
                (chiN=Float64(case.chiN), period_rg=Float64(scft.period_rg)))
        end

        summary_rows = NamedTuple[]
        profile_rows = NamedTuple[]
        native_rows = NamedTuple[]
        native_scft_path = "cases/$(case_id)/native_profiles.csv#model=scft"
        push!(summary_rows, lpc_scft_summary_row(case, scft, case_id,
            native_scft_path))
        scft_common = _resample(scft.phi_a, sample_count)
        lpc_push_common_profile!(profile_rows, native_rows, case, case_id,
            "scft", scft.period_rg, scft.phi_a, scft_common;
            phi_b_raw=scft.phi_b, sample_count=sample_count)

        allow_swap = isapprox(case.f, 0.5; atol=1.0e-12)
        for model in models
            opf_coefficients = model == "liu2019_opf" ?
                fit_diblock_opf_coefficients_1d(scft.phi_a; f=case.f,
                    L=scft.period_rg / sqrt(6.0), lambda=100.0) : nothing
            ok_coefficients = model == "ohta_kawasaki" ?
                diblock_liu2019_ok_coefficients(
                    ; f=case.f, chiN=case.chiN,
                    check_domain=!allow_ok_extrapolation) : nothing
            key = (model, Float64(case.f))
            bootstrap = model_seed_policy == :reference ? Float64[1.0] :
                lpc_model_bootstrap_period_factors(history, key, case.chiN;
                    lower_factor=lower_factor, upper_factor=upper_factor,
                    prefer_same_state=model == "bvk2")
            model == "liu2019_opf" && period_strategy == :brent &&
                (bootstrap = Float64[])
            label = LPC_MODEL_LABELS[model]
            same_state_profile = if model == "bvk2" &&
                    lpc_existing_case_accepted(case_dir, [model];
                        bvk2_oversample=bvk2_oversample,
                        bvk2_sensor_filter_ratio=bvk2_sensor_filter_ratio,
                        bvk2_residual_tol=bvk2_force_tolerance,
                        bvk2_maximum_residual_tol=
                            bvk2_force_maxabs_tolerance,
                        bvk2_stress_tol=bvk2_cell_stress_tolerance)
                lpc_load_existing_model_profile(case_dir, model, bvk2_nx;
                    resample_if_needed=true)
            elseif model in ("burp_ti", "uneyama_doi") &&
                    lpc_existing_model_accepted(case_dir, model)
                lpc_load_existing_model_profile(case_dir, model, nx)
            else
                nothing
            end
            if model_seed_policy == :history &&
                    model in ("burp_ti", "uneyama_doi") &&
                    same_state_profile !== nothing
                same_state_factor = lpc_load_existing_model_period_factor(
                    case_dir, model)
                same_state_factor === nothing ||
                    (bootstrap = Float64[same_state_factor])
            end
            initial_profile = same_state_profile === nothing ?
                (model_seed_policy == :reference ?
                    lpc_nearest_model_profile(profile_history, key, case.chiN) :
                    lpc_nearest_lower_model_profile(profile_history, key,
                        case.chiN)) : same_state_profile
            if model == "bvk2"
                seed_description = initial_profile === nothing ? "sinusoidal" :
                    "profile[n=$(length(initial_profile))]"
                println("[", case_id, "] ", label, " seed factors=",
                    join(bootstrap, ";"), " field=", seed_description)
            end
            println("[", case_id, "] ", label, " start")
            flush(stdout)
            optimized = try
                lpc_solve_model(model; f=case.f, chiN=case.chiN,
                    nx=nx, bvk2_nx=bvk2_nx, nquad=nquad, c2=c2,
                    bvk2_oversample=bvk2_oversample,
                    bvk2_sensor_filter_ratio=bvk2_sensor_filter_ratio,
                    bvk2_energy_tolerance=bvk2_energy_tolerance,
                    bvk2_force_tolerance=bvk2_force_tolerance,
                    bvk2_force_maxabs_tolerance=
                        bvk2_force_maxabs_tolerance,
                    bvk2_cell_stress_tolerance=
                        bvk2_cell_stress_tolerance,
                    bvk2_log_period_tolerance=
                        bvk2_log_period_tolerance,
                    max_iterations=max_iterations,
                    max_period_iterations=max_period_iterations,
                    burp_max_newton_iterations=burp_max_newton_iterations,
                    ok_gradient_tolerance=ok_gradient_tolerance,
                    bootstrap_period_factors=bootstrap,
                    lower_factor=lower_factor, upper_factor=upper_factor,
                    bootstrap_window=bootstrap_window,
                    local_check_fraction=local_check_fraction,
                    progress_label="$(case_id) $(label)",
                    period_strategy=period_strategy,
                    opf_coefficients=opf_coefficients,
                    ok_coefficients=ok_coefficients,
                    initial_profile=initial_profile)
            catch err
                message = sprint(showerror, err)
                println("[", case_id, "] ", label, " failed: ", message)
                push!(summary_rows, lpc_failed_model_row(case, case_id,
                    model, scft, message, c2;
                    bvk2_oversample=bvk2_oversample,
                    bvk2_sensor_filter_ratio=bvk2_sensor_filter_ratio))
                flush(stdout)
                continue
            end
            result = optimized.result
            target, alignment = _align_periodic(scft.phi_a, result.phi_a;
                sample_count=sample_count, allow_swap=allow_swap)
            stats = lpc_profile_statistics(target, alignment.values)
            native_path = "cases/$(case_id)/native_profiles.csv#model=$(model)"
            push!(summary_rows, lpc_model_summary_row(case, case_id, model,
                optimized, scft, alignment, stats, native_path, c2;
                opf_coefficients=opf_coefficients,
                ok_coefficients=ok_coefficients, burp_nquad=nquad,
                burp_max_newton_iterations=burp_max_newton_iterations,
                bvk2_oversample=bvk2_oversample,
                bvk2_sensor_filter_ratio=bvk2_sensor_filter_ratio,
                bvk2_residual_tolerance=bvk2_force_tolerance,
                bvk2_maximum_residual_tolerance=
                    bvk2_force_maxabs_tolerance,
                bvk2_stress_tolerance=bvk2_cell_stress_tolerance))
            lpc_push_common_profile!(profile_rows, native_rows, case, case_id,
                model, optimized.period_rg, result.phi_a, alignment.values;
                sample_count=sample_count)
            if last(summary_rows).status == "accepted"
                push!(get!(history, key, NamedTuple[]),
                    (chiN=Float64(case.chiN), factor=optimized.period_factor))
                push!(get!(profile_history, key, NamedTuple[]),
                    (chiN=Float64(case.chiN),
                        profile=Float64.(result.phi_a)))
            end
            println(@sprintf("[%s] %s done L/Rg=%.8g eta=%+.5g RMS=%.5g status=%s",
                case_id, label, optimized.period_rg,
                optimized.period_rg / scft.period_rg - 1.0,
                stats.rms, last(summary_rows).status))
            flush(stdout)
        end
        lpc_write_case(case_dir, summary_rows, profile_rows, native_rows)
        merge_outputs && lpc_merge_cases(outdir)
        println("[", case_id, "] checkpoint written")
        flush(stdout)
    end
    merged = merge_outputs ? lpc_merge_cases(outdir) : nothing
    merged === nothing || println("Wrote ", joinpath(outdir, "summary.csv"))
    return merged
end

if abspath(PROGRAM_FILE) == @__FILE__
    lpc_main()
end
