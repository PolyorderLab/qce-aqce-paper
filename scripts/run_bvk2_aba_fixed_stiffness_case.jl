#!/usr/bin/env julia

const ABA_FIXED_ROOT = normpath(joinpath(@__DIR__, ".."))
ABA_FIXED_ROOT in LOAD_PATH || pushfirst!(LOAD_PATH, ABA_FIXED_ROOT)

using CSV
using DataFrames
using DFMMonteCarlo
using FFTW
using Printf
using SHA

const ABA_FIXED_PREVIOUS_SKIP = get(
    ENV, "DFM_SKIP_BVK2_ABA_COMPREHENSIVE_MAIN", nothing)
ENV["DFM_SKIP_BVK2_ABA_COMPREHENSIVE_MAIN"] = "1"
include(joinpath(@__DIR__, "write_bvk2_aba_comprehensive_validation.jl"))
if ABA_FIXED_PREVIOUS_SKIP === nothing
    delete!(ENV, "DFM_SKIP_BVK2_ABA_COMPREHENSIVE_MAIN")
else
    ENV["DFM_SKIP_BVK2_ABA_COMPREHENSIVE_MAIN"] = ABA_FIXED_PREVIOUS_SKIP
end

const ABA_FIXED_SCHEMA = "bvk2-aba-fixed-stiffness-validation-v1"
const ABA_FIXED_C2 = 0.16
const ABA_FIXED_NX = 128
const ABA_FIXED_PROFILE_SAMPLES = 256
const ABA_FIXED_CASES = Dict(
    "f0.50_chiN20" => (f=0.5, chiN=20.0),
    "f0.50_chiN22" => (f=0.5, chiN=22.0),
    "f0.50_chiN25" => (f=0.5, chiN=25.0),
    "f0.50_chiN30" => (f=0.5, chiN=30.0),
    "f0.50_chiN35" => (f=0.5, chiN=35.0),
    "f0.50_chiN40" => (f=0.5, chiN=40.0),
    "f0.50_chiN45" => (f=0.5, chiN=45.0),
)

function aba_fixed_arg(name::AbstractString, default=nothing)
    prefix = String(name) * "="
    for argument in ARGS
        startswith(argument, prefix) &&
            return split(argument, "=", limit=2)[2]
    end
    default === nothing && throw(ArgumentError("missing $(name)=..."))
    return default
end

aba_fixed_sha256(path::AbstractString) = bytes2hex(sha256(read(path)))

function aba_fixed_domain_count(profile::AbstractVector{<:Real}, threshold::Real)
    mask = Float64.(profile) .> Float64(threshold)
    (all(mask) || all(.!mask)) && return 0
    return count(index -> mask[index] &&
        !mask[mod1(index - 1, length(mask))], eachindex(mask))
end

function aba_fixed_dominant_mode(profile::AbstractVector{<:Real}, mean_value::Real)
    spectrum = abs.(rfft(Float64.(profile) .- Float64(mean_value)))
    length(spectrum) >= 2 || return 0
    return argmax(@view spectrum[2:end])
end

function aba_fixed_source(case_id::AbstractString)
    directory = joinpath(ABA_FIXED_ROOT, "results",
        "bvk2_aba_comprehensive_validation")
    summary_path = joinpath(directory, "benchmark_summary.csv")
    profile_path = joinpath(directory, "benchmark_profiles.csv")
    summary = CSV.read(summary_path, DataFrame)
    profiles = CSV.read(profile_path, DataFrame)
    scft = summary[
        (String.(summary.case_id) .== case_id) .&
        (String.(summary.model) .== "scft"), :]
    adaptive = summary[
        (String.(summary.case_id) .== case_id) .&
        (String.(summary.model) .== "bvk2"), :]
    nrow(scft) == 1 && nrow(adaptive) == 1 || throw(ArgumentError(
        "expected unique accepted SCFT and BVK2 sources for $(case_id)"))
    Bool(scft.accepted[1]) && Bool(adaptive.accepted[1]) ||
        throw(ArgumentError("source rows are not accepted for $(case_id)"))
    scft_profile = profiles[
        (String.(profiles.case_id) .== case_id) .&
        (String.(profiles.model) .== "scft"), :]
    adaptive_profile = profiles[
        (String.(profiles.case_id) .== case_id) .&
        (String.(profiles.model) .== "bvk2"), :]
    sort!(scft_profile, :s)
    sort!(adaptive_profile, :s)
    nrow(scft_profile) == ABA_FIXED_PROFILE_SAMPLES &&
        nrow(adaptive_profile) == ABA_FIXED_PROFILE_SAMPLES ||
        throw(ArgumentError("source profiles are incomplete for $(case_id)"))
    return (
        scft_row=scft[1, :],
        adaptive_row=adaptive[1, :],
        scft_profile=Float64.(scft_profile.phi_a),
        adaptive_profile=Float64.(adaptive_profile.phi_a),
        summary_path=summary_path,
        profile_path=profile_path,
        summary_sha256=aba_fixed_sha256(summary_path),
        profile_sha256=aba_fixed_sha256(profile_path),
    )
end

function aba_fixed_validate(root_row, profile_rows)
    root_row.schema == ABA_FIXED_SCHEMA || error("wrong root schema")
    root_row.model == "bvk2_fixed" && !root_row.adaptive ||
        error("wrong fixed-stiffness model metadata")
    root_row.c2 == ABA_FIXED_C2 && root_row.nx == ABA_FIXED_NX ||
        error("wrong fixed-stiffness numerical settings")
    root_row.accepted && root_row.status == "accepted" ||
        error("fixed-stiffness state did not pass every gate")
    length(profile_rows) == ABA_FIXED_PROFILE_SAMPLES ||
        error("wrong fixed-stiffness profile sample count")
    all(isfinite(row.phi_a) for row in profile_rows) ||
        error("fixed-stiffness profile contains nonfinite values")
    return true
end

function aba_fixed_publish(outdir::AbstractString, root_row, profile_rows)
    mkpath(outdir)
    stage = mktempdir(dirname(outdir); prefix=".aba-fixed-stage-")
    staged_root = joinpath(stage, "root.csv")
    staged_profile = joinpath(stage, "profile.csv")
    try
        CSV.write(staged_root, DataFrame([root_row]); newline='\n')
        CSV.write(staged_profile, DataFrame(profile_rows); newline='\n')
        aba_fixed_validate(root_row, profile_rows)
        mv(staged_profile, joinpath(outdir, "profile.csv"); force=true)
        mv(staged_root, joinpath(outdir, "root.csv"); force=true)
    finally
        rm(stage; recursive=true, force=true)
    end
end

function main()
    case_id = String(aba_fixed_arg("--case"))
    haskey(ABA_FIXED_CASES, case_id) ||
        throw(ArgumentError("unsupported ABA case $(case_id)"))
    case = ABA_FIXED_CASES[case_id]
    outdir = normpath(String(aba_fixed_arg("--outdir", joinpath(
        ABA_FIXED_ROOT, "results", "bvk2_aba_fixed_stiffness_validation",
        "cases", case_id))))
    source = aba_fixed_source(case_id)

    reference = aba_cv_model_reference("bvk2_fixed", case.f)
    target_factor = Float64(source.scft_row.period_rg) /
        (sqrt(6.0) * reference)
    seed = _resample(source.adaptive_profile, ABA_FIXED_NX)
    root = minimize_triblock_aba_bvk2_lamella_stress_free(;
        f=case.f, chiN=case.chiN, nx=ABA_FIXED_NX,
        mode_count=ABA_CV_MODES, initial_amplitudes=(0.35,),
        initial_profile=seed, warm_start_profiles=true,
        c2=ABA_FIXED_C2, adaptive=false, max_iterations=4000,
        max_period_iterations=20,
        bootstrap_period_factors=(target_factor,),
        period_strategy=:bootstrap_local, bootstrap_window=0.12,
        local_check_fraction=0.01, local_check_tolerance=1.0e-8,
        physical_residual_tolerance=1.0e-4,
        physical_max_residual_tolerance=2.0e-4,
        force_tolerance=1.0e-4, force_maxabs_tolerance=2.0e-4,
        progress_label=case_id * " fixed stiffness")

    field = Float64.(root.result.phi_a)
    reference_profile, alignment = _align(source.scft_profile, field)
    sampled = Float64.(alignment.values)
    force_pass = root.result.projected_force_norm <= 1.0e-4 &&
        root.result.projected_force_maxabs <= 2.0e-4
    composition_pass = abs(mean(field) - case.f) <= 2.0e-8
    morphology_pass = maximum(field) - minimum(field) > 1.0e-4 &&
        root.result.energy < root.result.homogeneous_energy - 1.0e-10 &&
        aba_fixed_domain_count(field, case.f) == 1 &&
        aba_fixed_dominant_mode(field, case.f) == 1
    accepted = root.result.converged && force_pass && composition_pass &&
        morphology_pass && root.local_minimum_check_pass &&
        !root.boundary_limited && isfinite(root.objective)
    root_row = (
        schema=ABA_FIXED_SCHEMA,
        claim_kind="stress_free_lamellar_observable",
        case_id=case_id, fA=case.f, chiN=case.chiN,
        model="bvk2_fixed", label="BVK, fixed stiffness",
        adaptive=false, c2=ABA_FIXED_C2, nx=ABA_FIXED_NX,
        period_rg=root.period_rg,
        scft_period_rg=Float64(source.scft_row.period_rg),
        period_relative_error=root.period_rg /
            Float64(source.scft_row.period_rg) - 1.0,
        profile_rms=alignment.rms,
        profile_correlation=alignment.correlation,
        energy_density=root.objective,
        field_iterations=root.result.iterations,
        projected_force_norm=root.result.projected_force_norm,
        projected_force_maxabs=root.result.projected_force_maxabs,
        mean_phi=mean(field), minimum_phi=minimum(field),
        maximum_phi=maximum(field),
        local_left_objective=root.local_left_objective,
        local_right_objective=root.local_right_objective,
        field_gate_pass=root.result.converged && force_pass,
        composition_gate_pass=composition_pass,
        morphology_gate_pass=morphology_pass,
        cell_gate_pass=root.local_minimum_check_pass &&
            !root.boundary_limited,
        local_minimum_check_pass=root.local_minimum_check_pass,
        accepted=accepted, status=accepted ? "accepted" : "rejected",
        source_summary=relpath(source.summary_path, ABA_FIXED_ROOT),
        source_profiles=relpath(source.profile_path, ABA_FIXED_ROOT),
        source_summary_sha256=source.summary_sha256,
        source_profiles_sha256=source.profile_sha256,
    )
    profile_rows = [(
        schema=ABA_FIXED_SCHEMA, case_id=case_id, fA=case.f,
        chiN=case.chiN, model="bvk2_fixed", adaptive=false,
        c2=ABA_FIXED_C2, nx=ABA_FIXED_NX, period_rg=root.period_rg,
        index=index, s=(index - 1) / ABA_FIXED_PROFILE_SAMPLES,
        phi_a=sampled[index])
        for index in eachindex(sampled)]
    aba_fixed_publish(outdir, root_row, profile_rows)
    println(@sprintf(
        "status=%s period/Rg=%.10g error=%.5g RMS=%.5g force=%.3g maxforce=%.3g",
        root_row.status, root_row.period_rg,
        100.0 * root_row.period_relative_error, root_row.profile_rms,
        root_row.projected_force_norm, root_row.projected_force_maxabs))
end

if get(ENV, "DFM_SKIP_BVK2_ABA_FIXED_MAIN", "0") != "1"
    main()
end
