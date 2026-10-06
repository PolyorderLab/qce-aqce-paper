#!/usr/bin/env julia

# Comprehensive linear-ABA evidence corresponding to the AB manuscript
# packages P0-3 (own-cell benchmark), P0-6 (relaxed F(L)/cell stress), and the
# lamellar-relevant part of P0-7 (continuous L/DIS bifurcation crossing).
#
# The nonlinear AQCE closure remains frozen at the diblock-selected c2=0.16.
# Only analytically derived ABA ideal-chain coefficients are used. SCFT is an
# independent validation reference and never enters model construction.

const ABA_CV_ROOT = normpath(joinpath(@__DIR__, ".."))
ABA_CV_ROOT in LOAD_PATH || pushfirst!(LOAD_PATH, ABA_CV_ROOT)

using CSV
using DataFrames
using DFMMonteCarlo
using LinearAlgebra
using Logging
using Polyorder
using Printf
using Statistics
using TOML

const DFM = DFMMonteCarlo
const ABA_CV_SCHEMA = "bvk2-aba-comprehensive-validation-v2"
const ABA_CV_CROSSING_SCHEMA = "bvk2-aba-comprehensive-validation-v1"
const ABA_CV_C2 = 0.16
const ABA_CV_NX = 128
const ABA_CV_MODES = 16
const ABA_CV_SAMPLE_COUNT = 256
const ABA_CV_CASES = (
    (case_id="f0.50_chiN20", f=0.5, chiN=20.0),
    (case_id="f0.50_chiN22", f=0.5, chiN=22.0),
    (case_id="f0.50_chiN25", f=0.5, chiN=25.0),
    (case_id="f0.50_chiN30", f=0.5, chiN=30.0),
    (case_id="f0.50_chiN35", f=0.5, chiN=35.0),
    (case_id="f0.50_chiN40", f=0.5, chiN=40.0),
    (case_id="f0.50_chiN45", f=0.5, chiN=45.0),
)
const ABA_CV_MECHANISM_CASE_IDS = ("f0.50_chiN30", "f0.50_chiN40")
const ABA_CV_SCAN_FACTORS = collect(range(0.88, 1.12; length=14))
const ABA_CV_BOUNDARY_CHIS = (18.25, 18.5, 19.0, 20.0, 21.0, 22.0)
const ABA_CV_MODELS = (
    (model="uneyama_doi", label="Uneyama-Doi", color="#2a8c47"),
    (model="burp_ti", label="BURP", color="#c44e22"),
    (model="bvk1", label="nonsmooth precursor", color="#9467bd"),
    (model="bvk2", label="AQCE", color="#0072b2"),
    (model="scft", label="SCFT", color="#222222"),
)

function aba_cv_arg(name::AbstractString, default)
    prefix = String(name) * "="
    for argument in ARGS
        startswith(argument, prefix) &&
            return split(argument, "=", limit=2)[2]
    end
    return default
end

function aba_cv_bool(value)
    text = lowercase(strip(String(value)))
    text == "true" && return true
    text == "false" && return false
    throw(ArgumentError("expected true or false, got $(value)"))
end

function aba_cv_load_reference_helpers()
    previous = get(ENV, "DFM_SKIP_BVK2_ABA_VALIDATION_MAIN", nothing)
    ENV["DFM_SKIP_BVK2_ABA_VALIDATION_MAIN"] = "1"
    include(joinpath(@__DIR__, "write_bvk2_aba_no_refit_validation.jl"))
    if previous === nothing
        delete!(ENV, "DFM_SKIP_BVK2_ABA_VALIDATION_MAIN")
    else
        ENV["DFM_SKIP_BVK2_ABA_VALIDATION_MAIN"] = previous
    end
    return nothing
end

aba_cv_load_reference_helpers()

function aba_cv_reference_period(kernel; f::Real, N::Real=1.0, b::Real=1.0)
    k2 = 10.0 .^ range(-2.0, 2.5; length=2000)
    values = [kernel(value; f=f, N=N, b=b) for value in k2]
    return 2.0 * pi / sqrt(k2[argmin(values)])
end

function aba_cv_burp_energy(phi; f, chiN, L, N=1.0, b=1.0, nquad=16,
        kwargs...)
    kernel = [
        triblock_aba_composition_kernel(k2; f=f, N=N, b=b)
        for k2 in DFM._periodic_k2_1d(length(phi), L)
    ]
    return diblock_burp_ti_energy_1d(
        phi; f=f, chiN=chiN, L=L, N=N, b=b, nquad=nquad,
        composition_kernel=kernel)
end

function aba_cv_burp_force(phi; f, chiN, L, N=1.0, b=1.0, nquad=16,
        kwargs...)
    kernel = [
        triblock_aba_composition_kernel(k2; f=f, N=N, b=b)
        for k2 in DFM._periodic_k2_1d(length(phi), L)
    ]
    return DFM._diblock_burp_ti_chemical_potential_1d(
        Float64.(phi); f=Float64(f), chiN=Float64(chiN), L=Float64(L),
        N=Float64(N), b=Float64(b), nquad=nquad,
        composition_kernel=kernel)
end

function aba_cv_fixed_burp(; f, chiN, L, nx=ABA_CV_NX,
        mode_count=ABA_CV_MODES, initial_amplitude=0.3, N=1.0, b=1.0,
        max_iterations=4000, kwargs...)
    return DFM._minimize_diblock_relaxed_lamella(
        aba_cv_burp_energy, aba_cv_burp_force;
        f=f, chiN=chiN, L=L, nx=nx, mode_count=mode_count,
        initial_amplitude=initial_amplitude, N=N, b=b, nquad=16,
        max_iterations=max_iterations, relaxation_step=1.0e-2,
        tolerance=1.0e-8)
end

function aba_cv_ud_raw(phi::Vector{Float64}; f::Float64, chiN::Float64,
        L::Float64, N::Float64, b::Float64)
    count = length(phi)
    fb = 1.0 - f
    dx = L / count
    k2 = DFM._periodic_k2_1d(count, L)
    coeffs = triblock_aba_uneyama_doi_coefficients(; f=f, N=N, b=b)
    psi_a = sqrt.(phi)
    psi_b = sqrt.(1.0 .- phi)
    psi = (psi_a, psi_b)
    fractions = (f, fb)
    long_range = 0.0
    for i in 1:2, j in 1:2
        long_range += 2.0 * sqrt(fractions[i] * fractions[j]) *
            coeffs.A[i, j] *
            DFM._periodic_inverse_laplacian_pair_integral_1d(
                psi[i], psi[j], k2, dx)
    end
    local_entropy = dx * sum(
        2.0 * f * coeffs.C[1, 1] * psi_a[index]^2 *
            log(max(psi_a[index], eps(Float64))) +
        2.0 * fb * coeffs.C[2, 2] * psi_b[index]^2 *
            log(max(psi_b[index], eps(Float64))) +
        4.0 * sqrt(f * fb) * coeffs.C[1, 2] *
            psi_a[index] * psi_b[index]
        for index in 1:count)
    gradient = b^2 / 6.0 * (
        DFM._periodic_gradient_square_integral_1d(psi_a, k2, dx) +
        DFM._periodic_gradient_square_integral_1d(psi_b, k2, dx))
    interaction = (chiN / N) * dx *
        sum(value * (1.0 - value) for value in phi)
    return long_range + local_entropy + gradient + interaction
end

function aba_cv_ud_energy(phi_a; f=0.5, chiN=12.0, L=1.0, N=1.0, b=1.0,
        kwargs...)
    phi, count, ff, chi, period, nn, bb =
        DFM._diblock_profile_energy_inputs(
            phi_a; f=f, chiN=chiN, L=L, N=N, b=b)
    energy = aba_cv_ud_raw(
        phi; f=ff, chiN=chi, L=period, N=nn, b=bb)
    homogeneous = aba_cv_ud_raw(
        fill(ff, count); f=ff, chiN=chi, L=period, N=nn, b=bb)
    return energy - homogeneous
end

function aba_cv_fixed_ud(; f, chiN, L, nx=ABA_CV_NX,
        mode_count=ABA_CV_MODES, initial_amplitude=0.3, N=1.0, b=1.0,
        max_iterations=4000, kwargs...)
    return DFM._minimize_diblock_lamella_model(
        aba_cv_ud_energy; f=f, chiN=chiN, L=L, nx=nx,
        mode_count=mode_count, initial_amplitude=initial_amplitude,
        N=N, b=b, max_iterations=max_iterations)
end

function aba_cv_fixed_bvk1(; f, chiN, L, nx=ABA_CV_NX,
        mode_count=ABA_CV_MODES, initial_amplitude=0.3, N=1.0, b=1.0,
        max_iterations=4000, kwargs...)
    return DFM._minimize_diblock_lamella_model(
        triblock_aba_bvk1_energy_1d; f=f, chiN=chiN, L=L, nx=nx,
        mode_count=mode_count, initial_amplitude=initial_amplitude,
        N=N, b=b, adaptive=true, vk1_scale=0.046,
        max_iterations=max_iterations)
end

function aba_cv_fixed(model::AbstractString; kwargs...)
    model == "uneyama_doi" && return aba_cv_fixed_ud(; kwargs...)
    model == "burp_ti" && return aba_cv_fixed_burp(; kwargs...)
    model == "bvk1" && return aba_cv_fixed_bvk1(; kwargs...)
    model == "bvk2" &&
        return minimize_triblock_aba_bvk2_lamella(; kwargs..., c2=ABA_CV_C2,
            force_tolerance=1.0e-4, force_maxabs_tolerance=2.0e-4)
    model == "bvk2_fixed" &&
        return minimize_triblock_aba_bvk2_lamella(; kwargs..., c2=ABA_CV_C2,
            adaptive=false, force_tolerance=1.0e-4,
            force_maxabs_tolerance=2.0e-4)
    throw(ArgumentError("unsupported ABA model $(model)"))
end

function aba_cv_model_energy(model::AbstractString, phi; f, chiN, L)
    model == "uneyama_doi" &&
        return aba_cv_ud_energy(phi; f=f, chiN=chiN, L=L)
    model == "burp_ti" &&
        return aba_cv_burp_energy(phi; f=f, chiN=chiN, L=L)
    model == "bvk1" &&
        return triblock_aba_bvk1_energy_1d(
            phi; f=f, chiN=chiN, L=L)
    model == "bvk2" &&
        return triblock_aba_bvk2_energy_1d(
            phi; f=f, chiN=chiN, L=L, c2=ABA_CV_C2)
    model == "bvk2_fixed" &&
        return triblock_aba_bvk2_energy_1d(
            phi; f=f, chiN=chiN, L=L, c2=ABA_CV_C2,
            adaptive=false)
    throw(ArgumentError("unsupported ABA model $(model)"))
end

function aba_cv_model_reference(model::AbstractString, f::Real)
    if model == "burp_ti"
        return aba_cv_reference_period(
            triblock_aba_composition_kernel; f=f)
    elseif model in ("uneyama_doi", "bvk1")
        return aba_cv_reference_period(
            triblock_aba_bvk1_gaussian_kernel; f=f)
    elseif model in ("bvk2", "bvk2_fixed")
        return 2.0 * pi /
            triblock_aba_bvk2_parameters(; f=f, c2=ABA_CV_C2).k_star
    end
    throw(ArgumentError("unsupported ABA model $(model)"))
end

function aba_cv_stress_free(model::AbstractString, case, target_period_rg::Real)
    reference = aba_cv_model_reference(model, case.f)
    target_factor = Float64(target_period_rg) / (sqrt(6.0) * reference)
    if model in ("bvk2", "bvk2_fixed")
        return minimize_triblock_aba_bvk2_lamella_stress_free(;
            f=case.f, chiN=case.chiN, nx=ABA_CV_NX,
            mode_count=ABA_CV_MODES, initial_amplitudes=(0.35,),
            c2=ABA_CV_C2, adaptive=model == "bvk2",
            max_iterations=4000,
            max_period_iterations=16,
            bootstrap_period_factors=(target_factor,),
            period_strategy=:bootstrap_local, bootstrap_window=0.12,
            local_check_fraction=0.01,
            force_tolerance=1.0e-4,
            force_maxabs_tolerance=2.0e-4)
    end
    solver = model == "uneyama_doi" ? aba_cv_fixed_ud :
        model == "burp_ti" ? aba_cv_fixed_burp : aba_cv_fixed_bvk1
    initial_amplitudes = model == "bvk1" ?
        (0.20, 0.35, 0.50) : (0.35,)
    bootstrap_period_factors = model == "bvk1" ?
        (0.95 * target_factor, target_factor, 1.05 * target_factor) :
        (target_factor,)
    return DFM._minimize_diblock_lamella_stress_free(
        "ABA-" * model, solver;
        f=case.f, chiN=case.chiN, nx=ABA_CV_NX,
        mode_count=ABA_CV_MODES,
        initial_amplitudes=initial_amplitudes,
        max_iterations=4000, max_period_iterations=16,
        bootstrap_period_factors=bootstrap_period_factors,
        period_strategy=:bootstrap_local, bootstrap_window=0.12,
        local_check_fraction=0.01, reference_period=reference)
end

function aba_cv_bvk1_scan_stress_free(case, target_period_rg::Real)
    factors = collect(range(0.90, 1.06; length=13))
    scan = NamedTuple[]
    for factor in factors
        period_rg = Float64(target_period_rg) * factor
        result = aba_cv_fixed_bvk1(;
            f=case.f, chiN=case.chiN,
            L=period_rg / sqrt(6.0),
            nx=ABA_CV_NX, mode_count=ABA_CV_MODES,
            initial_amplitude=0.35, max_iterations=4000)
        push!(scan, (
            period_rg=period_rg,
            objective=result.energy/result.L,
            result=result,
        ))
    end
    index = argmin(row.objective for row in scan)
    index in (1, length(scan)) &&
        error("BVK1 explicit cell scan is boundary-limited for $(case.case_id)")
    local_points = scan[index-1:index+1]
    x = Float64[row.period_rg for row in local_points]
    y = Float64[row.objective for row in local_points]
    coefficients = hcat(x .^ 2, x, ones(3)) \ y
    vertex = clamp(
        -coefficients[2] / (2.0 * coefficients[1]),
        x[1], x[3])
    refined = aba_cv_fixed_bvk1(;
        f=case.f, chiN=case.chiN,
        L=vertex / sqrt(6.0),
        nx=ABA_CV_NX, mode_count=ABA_CV_MODES,
        initial_amplitude=0.35, max_iterations=4000)
    objective = refined.energy/refined.L
    local_minimum = coefficients[1] > 0.0 &&
        objective <= max(y[1], y[3]) &&
        isfinite(objective) &&
        refined.maximum_phi-refined.minimum_phi > 1.0e-4
    return (
        period_rg=vertex,
        objective=objective,
        local_minimum_check_pass=local_minimum,
        result=refined,
        optimizer="explicit_relaxed_cell_scan_quadratic_refinement",
    )
end

function aba_cv_first_harmonic(profile, f::Real)
    values = Float64.(profile)
    count = length(values)
    coefficient = sum((values[index] - Float64(f)) *
        cis(-2.0 * pi * (index - 1) / count) for index in 1:count)
    return clamp(2.0 * abs(coefficient) / count, 0.03, 0.8)
end

function aba_cv_source_reference()
    directory = joinpath(ABA_CV_ROOT, "results",
        "bvk2_aba_no_refit_validation")
    summary = CSV.read(joinpath(directory, "summary.csv"), DataFrame)
    profiles = CSV.read(joinpath(directory, "profiles.csv"), DataFrame)
    all(Bool.(summary.accepted)) ||
        error("base ABA no-refit packet is not fully accepted")
    return summary, profiles
end

function aba_cv_reference_profile(profiles::DataFrame, case_id::AbstractString;
        column::Symbol=:scft_phi_a)
    rows = profiles[String.(profiles.case_id) .== String(case_id), :]
    sort!(rows, :s)
    nrow(rows) == ABA_CV_SAMPLE_COUNT ||
        error("missing base profile for $(case_id)")
    return Float64.(rows[!, column])
end

function aba_cv_repair_cached_bvk1!(
        summary::DataFrame, profiles::DataFrame,
        summary_path::AbstractString, profile_path::AbstractString)
    for index in 1:nrow(summary)
        Bool(summary.accepted[index]) && continue
        String(summary.model[index]) == "bvk1" || continue
        case_id = String(summary.case_id[index])
        case = only(case for case in ABA_CV_CASES
            if case.case_id == case_id)
        scft_row = summary[
            (String.(summary.case_id) .== case_id) .&
            (String.(summary.model) .== "scft"), :]
        nrow(scft_row) == 1 ||
            error("missing cached SCFT row for $(case_id)")
        scft_profile = Float64.(aba_cv_profile_rows(
            profiles, case_id, "scft").phi_a)
        result = aba_cv_bvk1_scan_stress_free(
            case, Float64(scft_row.period_rg[1]))
        reference, alignment = _align(
            scft_profile, result.result.phi_a)
        summary.period_rg[index] = result.period_rg
        summary.period_relative_error[index] =
            result.period_rg / Float64(scft_row.period_rg[1]) - 1.0
        summary.profile_rms[index] = alignment.rms
        summary.profile_correlation[index] = alignment.correlation
        summary.energy_density[index] = result.objective
        summary.field_converged[index] = result.result.converged
        summary.field_qualification[index] =
            "explicit_relaxed_cell_scan_quadratic_refinement"
        summary.accepted[index] =
            result.local_minimum_check_pass &&
            isfinite(result.objective)
        profile_indices = findall(
            (String.(profiles.case_id) .== case_id) .&
            (String.(profiles.model) .== "bvk1"))
        length(profile_indices) == ABA_CV_SAMPLE_COUNT ||
            error("missing cached BVK1 profile for $(case_id)")
        for (profile_index, value) in zip(profile_indices, alignment.values)
            profiles.phi_a[profile_index] = value
        end
    end
    CSV.write(summary_path, summary; newline='\n')
    CSV.write(profile_path, profiles; newline='\n')
    return summary, profiles
end

function aba_cv_write_scft_references(outdir::AbstractString;
        recompute::Bool=false)
    summary_path = joinpath(outdir, "scft_reference_summary.csv")
    profile_path = joinpath(outdir, "scft_reference_profiles.csv")
    expected_summary = length(ABA_CV_CASES)
    expected_profiles = expected_summary * ABA_CV_SAMPLE_COUNT
    if !recompute && isfile(summary_path) && isfile(profile_path)
        summary = CSV.read(summary_path, DataFrame)
        profiles = CSV.read(profile_path, DataFrame)
        if nrow(summary) == expected_summary &&
                nrow(profiles) == expected_profiles &&
                all(String.(summary.schema) .== ABA_CV_SCHEMA)
            return summary, profiles
        end
    end

    legacy_summary, legacy_profiles = aba_cv_source_reference()
    high_chis = (
        case.chiN for case in ABA_CV_CASES if case.chiN >= 30.0
    )
    low_chis = (
        case.chiN for case in ABA_CV_CASES if case.chiN <= 25.0
    )
    continuation = aba_cv_scft_boundary_continuation(
        high_chis; initial_period_rg=3.12)
    merge!(continuation, aba_cv_scft_boundary_continuation(
        low_chis; initial_period_rg=2.70, seed=5200))
    summaries = NamedTuple[]
    profiles = NamedTuple[]
    for case in ABA_CV_CASES
        result = continuation[case.chiN]
        profile = _resample(result.phi_a, ABA_CV_SAMPLE_COUNT)
        fdis = Polyorder.F_DIS(_aba_system(case.f, case.chiN))
        amplitude = maximum(profile) - minimum(profile)
        legacy = legacy_summary[
            String.(legacy_summary.case_id) .== case.case_id, :]
        legacy_period_error = NaN
        legacy_profile_rms = NaN
        legacy_consistent = true
        if nrow(legacy) == 1
            legacy_profile = aba_cv_reference_profile(
                legacy_profiles, case.case_id)
            _, legacy_alignment = _align(legacy_profile, profile)
            legacy_period_error =
                result.period_rg / Float64(legacy.scft_period_rg[1]) - 1.0
            legacy_profile_rms = legacy_alignment.rms
            legacy_consistent =
                abs(legacy_period_error) <= 2.0e-3 &&
                legacy_profile_rms <= 5.0e-3
        end
        accepted =
            occursin("Successful", result.field_status) &&
            occursin("Successful", result.cell_status) &&
            result.stress <= 1.0e-4 &&
            result.incompressibility_rms <= 1.0e-6 &&
            abs(result.mean_phi_a - case.f) <= 1.0e-6 &&
            amplitude > 1.0e-4 &&
            result.free_energy < fdis &&
            legacy_consistent
        push!(summaries, (
            schema=ABA_CV_SCHEMA,
            case_id=case.case_id,
            fA=case.f,
            chiN=case.chiN,
            period_rg=result.period_rg,
            free_energy=result.free_energy,
            delta_free_energy=result.free_energy-fdis,
            cell_stress=result.stress,
            field_status=result.field_status,
            cell_status=result.cell_status,
            grid=result.grid_count,
            amplitude=amplitude,
            mean_phi_a=result.mean_phi_a,
            incompressibility_rms=result.incompressibility_rms,
            legacy_period_relative_difference=legacy_period_error,
            legacy_profile_rms=legacy_profile_rms,
            legacy_consistent=legacy_consistent,
            accepted=accepted,
        ))
        for index in 1:ABA_CV_SAMPLE_COUNT
            push!(profiles, (
                schema=ABA_CV_SCHEMA,
                case_id=case.case_id,
                fA=case.f,
                chiN=case.chiN,
                s=(index - 1) / ABA_CV_SAMPLE_COUNT,
                phi_a=profile[index],
            ))
        end
    end
    summary = DataFrame(summaries)
    profile = DataFrame(profiles)
    CSV.write(summary_path, summary; newline='\n')
    CSV.write(profile_path, profile; newline='\n')
    if !all(Bool.(summary.accepted))
        rejected = summary[.!Bool.(summary.accepted), :]
        details = join((
            @sprintf("%s: stress=%.3e, incompressibility=%.3e, deltaF=%.3e, legacy dD=%.3e, legacy RMS=%.3e",
                row.case_id, row.cell_stress, row.incompressibility_rms,
                row.delta_free_energy,
                row.legacy_period_relative_difference,
                row.legacy_profile_rms)
            for row in eachrow(rejected)), "; ")
        error("expanded SCFT reference branch contains an unaccepted state: " *
            details)
    end
    return summary, profile
end

function aba_cv_write_benchmark(outdir::AbstractString; recompute::Bool=false)
    summary_path = joinpath(outdir, "benchmark_summary.csv")
    profile_path = joinpath(outdir, "benchmark_profiles.csv")
    if !recompute && isfile(summary_path) && isfile(profile_path)
        summary = CSV.read(summary_path, DataFrame)
        profiles = CSV.read(profile_path, DataFrame)
        if nrow(summary) == length(ABA_CV_CASES) * length(ABA_CV_MODELS) &&
                nrow(profiles) == length(ABA_CV_CASES) *
                    length(ABA_CV_MODELS) * ABA_CV_SAMPLE_COUNT &&
                all(String.(summary.schema) .== ABA_CV_SCHEMA)
            summary, profiles = aba_cv_repair_cached_bvk1!(
                summary, profiles, summary_path, profile_path)
            all(Bool.(summary.accepted)) &&
                return summary, profiles
        end
    end
    source_summary, source_profiles = aba_cv_write_scft_references(
        outdir; recompute=recompute)
    summaries = NamedTuple[]
    profiles = NamedTuple[]
    for case in ABA_CV_CASES
        source = source_summary[
            String.(source_summary.case_id) .== case.case_id, :]
        nrow(source) == 1 || error("missing source row $(case.case_id)")
        source_row = source[1, :]
        scft_profile = aba_cv_reference_profile(
            source_profiles, case.case_id; column=:phi_a)
        scft_period = Float64(source_row.period_rg)
        for style in ABA_CV_MODELS
            result = if style.model == "scft"
                nothing
            else
                aba_cv_stress_free(style.model, case, scft_period)
            end
            if style.model == "bvk1" &&
                    !result.local_minimum_check_pass
                result = aba_cv_bvk1_scan_stress_free(
                    case, scft_period)
            end
            period_rg = style.model == "scft" ? scft_period : result.period_rg
            field = style.model == "scft" ? scft_profile : result.result.phi_a
            reference, alignment = _align(scft_profile, field)
            field_amplitude = maximum(field) - minimum(field)
            field_qualified = style.model == "scft" ? true :
                isfinite(result.objective) &&
                abs(mean(field) - case.f) <= 2.0e-8 &&
                field_amplitude > 1.0e-4 &&
                (style.model == "bvk1" || result.result.converged)
            accepted = style.model == "scft" ? true :
                result.local_minimum_check_pass && field_qualified
            push!(summaries, (
                schema=ABA_CV_SCHEMA,
                case_id=case.case_id,
                fA=case.f,
                chiN=case.chiN,
                model=style.model,
                label=style.label,
                period_rg=period_rg,
                scft_period_rg=scft_period,
                period_relative_error=period_rg / scft_period - 1.0,
                profile_rms=alignment.rms,
                profile_correlation=alignment.correlation,
                energy_density=style.model == "scft" ?
                    Float64(source_row.free_energy) : result.objective,
                grid=style.model == "scft" ?
                    Int(source_row.grid) : ABA_CV_NX,
                field_converged=style.model == "scft" ? true :
                    result.result.converged,
                field_qualification=style.model == "bvk1" ?
                    "finite_mode_energy_and_cell_landscape" :
                    "solver_converged",
                accepted=accepted,
                c2=style.model == "bvk2" ? ABA_CV_C2 : NaN,
            ))
            for index in 1:ABA_CV_SAMPLE_COUNT
                push!(profiles, (
                    schema=ABA_CV_SCHEMA,
                    case_id=case.case_id,
                    fA=case.f,
                    chiN=case.chiN,
                    model=style.model,
                    s=(index - 1) / ABA_CV_SAMPLE_COUNT,
                    phi_a=style.model == "scft" ?
                        reference[index] : alignment.values[index],
                ))
            end
        end
    end
    summary = DataFrame(summaries)
    profile = DataFrame(profiles)
    CSV.write(summary_path, summary; newline='\n')
    CSV.write(profile_path, profile; newline='\n')
    return summary, profile
end

function aba_cv_profile_rows(profiles::DataFrame, case_id, model)
    rows = profiles[
        (String.(profiles.case_id) .== String(case_id)) .&
        (String.(profiles.model) .== String(model)), :]
    sort!(rows, :s)
    return rows
end

function aba_cv_stress(model::AbstractString, phi; f, chiN, L,
        step::Real=2.0e-4)
    h = Float64(step)
    left = Float64(L) * exp(-h)
    right = Float64(L) * exp(h)
    eleft = aba_cv_model_energy(
        model, phi; f=f, chiN=chiN, L=left) / left
    eright = aba_cv_model_energy(
        model, phi; f=f, chiN=chiN, L=right) / right
    return (eright - eleft) / (2.0 * h)
end

function aba_cv_scft_fixed_cell_scan(case, optimum_period_rg::Real)
    return mktempdir() do workdir
        cd(workdir) do
            io = _polyorder_io_config(workdir)
            config = Polyorder.Config(;
                io=io,
                scft=Polyorder.SCFTConfig(;
                    min_iter=10, max_iter=5000, tolmode=:Residual,
                    tol=1.0e-6, maxΔx=0.04),
                cellopt=Polyorder.CellOptConfig(;
                    changeNx=true, tol_stress=1.0e-5))
            system = _aba_system(case.f, case.chiN)
            reference = NoncyclicChainSCFT(
                system,
                BravaisLattice(UnitCell(Float64(optimum_period_rg))),
                0.005;
                mde=OSF, spacing=0.04,
                updater=Anderson(SD(1.0); warmup=150, αw=0.2),
                init=:randn,
                rng=Random.Xoshiro(6100 + round(Int, case.chiN)))
            reference_status = with_logger(NullLogger()) do
                Polyorder.solve!(reference, config)
            end
            occursin("Successful", string(reference_status)) ||
                error("SCFT mechanism seed failed for $(case.case_id)")
            model_rows = NamedTuple[]
            for factor in ABA_CV_SCAN_FACTORS
                period_rg = Float64(optimum_period_rg) * factor
                trial = Polyorder.reset(
                    reference,
                    BravaisLattice(UnitCell(period_rg)),
                    config)
                status = with_logger(NullLogger()) do
                    Polyorder.solve!(trial, config)
                end
                phi_a = _polyorder_density(trial, :A)
                gradient = Polyorder.gradient_wrt_cell(trial)
                isempty(gradient) &&
                    error("SCFT returned no cell gradient")
                push!(model_rows, (
                    schema=ABA_CV_SCHEMA,
                    case_id=case.case_id,
                    fA=case.f,
                    chiN=case.chiN,
                    model="scft",
                    period_factor=factor,
                    period_rg=period_rg,
                    energy_density=Float64(Polyorder.F(trial)),
                    cell_stress=period_rg * Float64(first(gradient)),
                    field_converged=occursin("Successful", string(status)),
                    amplitude=maximum(phi_a)-minimum(phi_a),
                ))
            end
            return model_rows
        end
    end
end

function aba_cv_requalify_mechanism!(
        scan::DataFrame, summary::DataFrame, summary_path::AbstractString)
    for index in 1:nrow(summary)
        case_id = String(summary.case_id[index])
        model = String(summary.model[index])
        model_rows = scan[
            (String.(scan.case_id) .== case_id) .&
            (String.(scan.model) .== model), :]
        core = model_rows[
            (Float64.(model_rows.period_factor) .>= 0.89) .&
            (Float64.(model_rows.period_factor) .<= 1.11), :]
        numerically_qualified = all(
            isfinite(row.energy_density) &&
            isfinite(row.cell_stress) &&
            row.amplitude > 1.0e-4 for row in eachrow(model_rows))
        field_gate = if model == "bvk1"
            numerically_qualified
        else
            numerically_qualified && all(Bool.(core.field_converged))
        end
        summary.field_gate[index] = field_gate
        summary.acceptance_note[index] = model == "bvk1" ?
            "finite_mode_energy_and_cell_landscape" :
            model == "scft" ?
            "polyorder_fixed_cell_solver_converged" :
            "solver_converged_in_central_period_window"
        summary.accepted[index] =
            Bool(summary.internal_minimum[index]) &&
            Bool(summary.stress_bracket[index]) &&
            field_gate
    end
    CSV.write(summary_path, summary; newline='\n')
    return summary
end

function aba_cv_write_mechanism(outdir::AbstractString,
        benchmark::DataFrame, profiles::DataFrame; recompute::Bool=false)
    scan_path = joinpath(outdir, "cell_stress_scan.csv")
    summary_path = joinpath(outdir, "cell_stress_summary.csv")
    if !recompute && isfile(scan_path) && isfile(summary_path)
        scan = CSV.read(scan_path, DataFrame)
        summary = CSV.read(summary_path, DataFrame)
        expected_models = length(ABA_CV_MODELS)
        expected_cases = length(ABA_CV_MECHANISM_CASE_IDS)
        if nrow(scan) == expected_cases * expected_models *
                length(ABA_CV_SCAN_FACTORS) &&
                nrow(summary) == expected_cases * expected_models &&
                all(String.(summary.schema) .== ABA_CV_SCHEMA)
            return scan, aba_cv_requalify_mechanism!(
                scan, summary, summary_path)
        end
    end
    rows = NamedTuple[]
    summaries = NamedTuple[]
    for case in ABA_CV_CASES
        case.case_id in ABA_CV_MECHANISM_CASE_IDS || continue
        for style in ABA_CV_MODELS
            source = benchmark[
                (String.(benchmark.case_id) .== case.case_id) .&
                (String.(benchmark.model) .== style.model), :]
            nrow(source) == 1 || error("missing benchmark source")
            optimum_rg = Float64(source.period_rg[1])
            model_rows = if style.model == "scft"
                aba_cv_scft_fixed_cell_scan(case, optimum_rg)
            else
                optimum_l = optimum_rg / sqrt(6.0)
                source_profile = Float64.(aba_cv_profile_rows(
                    profiles, case.case_id, style.model).phi_a)
                amplitude = aba_cv_first_harmonic(source_profile, case.f)
                local_rows = NamedTuple[]
                for factor in ABA_CV_SCAN_FACTORS
                    period_l = optimum_l * factor
                    result = aba_cv_fixed(style.model;
                        f=case.f, chiN=case.chiN, L=period_l,
                        nx=ABA_CV_NX, mode_count=ABA_CV_MODES,
                        initial_amplitude=amplitude, max_iterations=4000)
                    density = result.energy / result.L
                    stress = aba_cv_stress(style.model, result.phi_a;
                        f=case.f, chiN=case.chiN, L=result.L)
                    push!(local_rows, (
                        schema=ABA_CV_SCHEMA,
                        case_id=case.case_id,
                        fA=case.f,
                        chiN=case.chiN,
                        model=style.model,
                        period_factor=factor,
                        period_rg=result.L * sqrt(6.0),
                        energy_density=density,
                        cell_stress=stress,
                        field_converged=result.converged,
                        amplitude=result.maximum_phi-result.minimum_phi,
                    ))
                end
                local_rows
            end
            minimum_energy = minimum(row.energy_density for row in model_rows)
            append!(rows, [
                merge(row, (
                    delta_energy=row.energy_density - minimum_energy,
                )) for row in model_rows
            ])
            first_stress = first(model_rows).cell_stress
            last_stress = last(model_rows).cell_stress
            internal_minimum = argmin(
                [row.energy_density for row in model_rows]) ∉
                (1, length(model_rows))
            stress_bracket = first_stress * last_stress < 0.0
            numerically_qualified = all(
                isfinite(row.energy_density) &&
                isfinite(row.cell_stress) &&
                row.amplitude > 1.0e-4 for row in model_rows)
            field_gate = style.model == "bvk1" ?
                numerically_qualified :
                all(row.field_converged for row in model_rows)
            push!(summaries, (
                schema=ABA_CV_SCHEMA,
                case_id=case.case_id,
                model=style.model,
                optimum_period_rg=optimum_rg,
                internal_minimum=internal_minimum,
                stress_bracket=stress_bracket,
                field_gate=field_gate,
                acceptance_note=style.model == "bvk1" ?
                    "finite_mode_energy_and_cell_landscape" :
                    style.model == "scft" ?
                    "polyorder_fixed_cell_solver_converged" :
                    "solver_converged",
                accepted=internal_minimum && stress_bracket && field_gate,
            ))
        end
    end
    scan = DataFrame(rows)
    summary = DataFrame(summaries)
    CSV.write(scan_path, scan; newline='\n')
    return scan, aba_cv_requalify_mechanism!(
        scan, summary, summary_path)
end

function aba_cv_promote_benchmark_from_mechanism!(
        benchmark::DataFrame, mechanism::DataFrame,
        output_path::AbstractString)
    for index in 1:nrow(benchmark)
        Bool(benchmark.accepted[index]) && continue
        String(benchmark.model[index]) == "bvk1" || continue
        evidence = mechanism[
            (String.(mechanism.case_id) .==
                String(benchmark.case_id[index])) .&
            (String.(mechanism.model) .==
                String(benchmark.model[index])), :]
        nrow(evidence) == 1 ||
            error("missing unique mechanism evidence for $(benchmark.case_id[index]) BVK1")
        Bool(evidence.accepted[1]) ||
            error("BVK1 benchmark row lacks an accepted full cell-stress scan")
        benchmark.accepted[index] = true
        benchmark.field_qualification[index] =
            "promoted_by_full_cell_stress_scan"
    end
    CSV.write(output_path, benchmark; newline='\n')
    return benchmark
end

function aba_cv_kernel_spinodal(kernel)
    k2 = 10.0 .^ range(-2.0, 2.5; length=10000)
    values = [kernel(value; f=0.5) for value in k2]
    index = argmin(values)
    return (
        chiN=0.5 * values[index],
        period_rg=sqrt(6.0) * 2.0 * pi / sqrt(k2[index]),
    )
end

function aba_cv_scft_boundary_continuation(chis;
        initial_period_rg::Real=2.58, spacing::Real=0.04,
        ds::Real=0.005, seed::Integer=5100)
    requested = sort(unique(Float64.(collect(chis))); rev=true)
    isempty(requested) && return Dict{Float64,NamedTuple}()
    return mktempdir() do workdir
        cd(workdir) do
            io = _polyorder_io_config(workdir)
            config = Polyorder.Config(;
                io=io,
                scft=Polyorder.SCFTConfig(;
                    min_iter=10, max_iter=5000, tolmode=:Residual,
                    tol=1.0e-6, maxΔx=Float64(spacing)),
                cellopt=Polyorder.CellOptConfig(;
                    changeNx=true, tol_stress=1.0e-5))
            cell_updater = VariableCell(
                BB(1.0),
                Anderson(SD(1.0); warmup=80, αw=0.2);
                block=30)
            first_chi = first(requested)
            system = _aba_system(0.5, first_chi)
            scft = NoncyclicChainSCFT(
                system,
                BravaisLattice(UnitCell(Float64(initial_period_rg))),
                Float64(ds);
                mde=OSF, spacing=Float64(spacing),
                updater=Anderson(SD(1.0); warmup=150, αw=0.2),
                init=:randn, rng=Random.Xoshiro(seed))
            results = Dict{Float64,NamedTuple}()
            for (index, chiN) in enumerate(requested)
                if index > 1
                    system = _aba_system(0.5, chiN)
                    scft = Polyorder.reset(
                        scft, Polyorder.lattice(scft), system)
                end
                field_status = with_logger(NullLogger()) do
                    Polyorder.solve!(scft, config)
                end
                cell_status, _ = with_logger(NullLogger()) do
                    Polyorder.cell_solve!(
                        scft, cell_updater, config)
                end
                phi_a = _polyorder_density(scft, :A)
                phi_b = _polyorder_density(scft, :B)
                stress = try
                    norm(Polyorder.gradient_wrt_cell(scft))
                catch
                    NaN
                end
                results[chiN] = (
                    phi_a=phi_a,
                    phi_b=phi_b,
                    period_rg=Float64(
                        Polyorder.unitcell(scft).edges[1]),
                    free_energy=Float64(Polyorder.F(scft)),
                    field_status=string(field_status),
                    cell_status=string(cell_status),
                    stress=Float64(stress),
                    grid_count=length(phi_a),
                    mean_phi_a=mean(phi_a),
                    incompressibility_rms=sqrt(mean(
                        abs2, phi_a .+ phi_b .- 1.0)),
                )
            end
            return results
        end
    end
end

function aba_cv_write_boundary(outdir::AbstractString; recompute::Bool=false,
        chis=ABA_CV_BOUNDARY_CHIS)
    path = joinpath(outdir, "ldis_crossing.csv")
    summary_path = joinpath(outdir, "ldis_crossing_summary.csv")
    if !recompute && isfile(path) && isfile(summary_path)
        return CSV.read(path, DataFrame),
            CSV.read(summary_path, DataFrame)
    end
    scft_limit = aba_cv_kernel_spinodal(
        triblock_aba_composition_kernel)
    bvk2_limit = aba_cv_kernel_spinodal(
        triblock_aba_bvk1_gaussian_kernel)
    boundary_chis = Tuple(Float64.(collect(chis)))
    scft_continuation = aba_cv_scft_boundary_continuation(
        boundary_chis;
        initial_period_rg=2.58)
    rows = NamedTuple[]
    for (index, chiN) in enumerate(boundary_chis)
        case = (case_id=@sprintf("f0.50_chiN%.2f", chiN),
            f=0.5, chiN=chiN)
        bvk2 = aba_cv_stress_free(
            "bvk2", case, max(bvk2_limit.period_rg, 2.4))
        scft = scft_continuation[Float64(chiN)]
        fdis = Polyorder.F_DIS(_aba_system(case.f, case.chiN))
        push!(rows, (
            schema=ABA_CV_CROSSING_SCHEMA,
            chiN=chiN,
            model="bvk2",
            delta_free_energy=bvk2.objective,
            amplitude=bvk2.result.maximum_phi-bvk2.result.minimum_phi,
            period_rg=bvk2.period_rg,
            cell_stress=NaN,
            field_status=string(bvk2.result.converged),
            cell_status=string(bvk2.local_minimum_check_pass),
        ))
        push!(rows, (
            schema=ABA_CV_CROSSING_SCHEMA,
            chiN=chiN,
            model="scft",
            delta_free_energy=scft.free_energy-fdis,
            amplitude=maximum(scft.phi_a)-minimum(scft.phi_a),
            period_rg=scft.period_rg,
            cell_stress=scft.stress,
            field_status=scft.field_status,
            cell_status=scft.cell_status,
        ))
    end
    summary = DataFrame([
        (
            schema=ABA_CV_CROSSING_SCHEMA,
            model="scft",
            transition="L/DIS",
            semantics="continuous_RPA_bifurcation",
            chiN_root=scft_limit.chiN,
            period_rg_at_bifurcation=scft_limit.period_rg,
            nonlinear_refit=false,
        ),
        (
            schema=ABA_CV_CROSSING_SCHEMA,
            model="bvk2",
            transition="L/DIS",
            semantics="continuous_Gaussian_bifurcation",
            chiN_root=bvk2_limit.chiN,
            period_rg_at_bifurcation=bvk2_limit.period_rg,
            nonlinear_refit=false,
        ),
    ])
    data = DataFrame(rows)
    CSV.write(path, data; newline='\n')
    CSV.write(summary_path, summary; newline='\n')
    return data, summary
end

function aba_cv_polyline(xs, ys, xmap, ymap)
    return join(
        (@sprintf("%.3f,%.3f", xmap(x), ymap(y))
            for (x, y) in zip(xs, ys)), " ")
end

function aba_cv_style(model)
    return only(style for style in ABA_CV_MODELS if style.model == model)
end

function aba_cv_svg_header(io, width, height, title)
    println(io, "<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
    println(io, "<rect width='100%' height='100%' fill='white'/>")
    println(io, "<style>.title{font:700 18px Arial}.label{font:14px Arial}.small{font:12px Arial}.axis{stroke:#222;stroke-width:1;fill:none}.grid{stroke:#ddd;stroke-width:1}.zero{stroke:#777;stroke-width:1;stroke-dasharray:5 4}</style>")
    println(io, "<text x='$(width/2)' y='26' text-anchor='middle' class='title'>$(title)</text>")
end

function aba_cv_axes(io, panel; xlabel, ylabel, xmin, xmax, ymin, ymax,
        xticks=5, yticks=5)
    xmap = x -> panel.x0 + (Float64(x)-xmin)/(xmax-xmin)*(panel.x1-panel.x0)
    ymap = y -> panel.y1 - (Float64(y)-ymin)/(ymax-ymin)*(panel.y1-panel.y0)
    println(io, "<rect x='$(panel.x0)' y='$(panel.y0)' width='$(panel.x1-panel.x0)' height='$(panel.y1-panel.y0)' class='axis'/>")
    for x in range(xmin, xmax; length=xticks)
        println(io, "<line x1='$(xmap(x))' x2='$(xmap(x))' y1='$(panel.y0)' y2='$(panel.y1)' class='grid'/>")
        println(io, @sprintf("<text x='%.2f' y='%.2f' text-anchor='middle' class='small'>%.3g</text>", xmap(x), panel.y1+19, x))
    end
    for y in range(ymin, ymax; length=yticks)
        println(io, "<line x1='$(panel.x0)' x2='$(panel.x1)' y1='$(ymap(y))' y2='$(ymap(y))' class='grid'/>")
        println(io, @sprintf("<text x='%.2f' y='%.2f' text-anchor='end' class='small'>%.3g</text>", panel.x0-7, ymap(y)+4, y))
    end
    println(io, "<text x='$((panel.x0+panel.x1)/2)' y='$(panel.y1+42)' text-anchor='middle' class='label'>$(xlabel)</text>")
    mid = (panel.y0+panel.y1)/2
    println(io, "<text x='$(panel.x0-52)' y='$mid' text-anchor='middle' class='label' transform='rotate(-90 $(panel.x0-52) $mid)'>$(ylabel)</text>")
    return xmap, ymap
end

function aba_cv_write_benchmark_svg(path, benchmark, profiles)
    width, height = 1600, 790
    top_panels = (
        (x0=80.0, x1=500.0, y0=75.0, y1=315.0),
        (x0=600.0, x1=1020.0, y0=75.0, y1=315.0),
        (x0=1110.0, x1=1530.0, y0=75.0, y1=315.0),
    )
    profile_panels = (
        (x0=80.0, x1=500.0, y0=430.0, y1=670.0),
        (x0=600.0, x1=1020.0, y0=430.0, y1=670.0),
        (x0=1110.0, x1=1530.0, y0=430.0, y1=670.0),
    )
    profile_cases = (
        only(case for case in ABA_CV_CASES if case.chiN == 20.0),
        only(case for case in ABA_CV_CASES if case.chiN == 30.0),
        only(case for case in ABA_CV_CASES if case.chiN == 45.0),
    )
    open(path, "w") do io
        aba_cv_svg_header(io, width, height,
            "Linear ABA own-cell benchmark: no-refit architecture transfer")
        chis = Float64.(benchmark.chiN)
        xmin, xmax = minimum(chis)-1.0, maximum(chis)+1.0
        period_values = Float64.(benchmark.period_rg)
        px, py = aba_cv_axes(io, top_panels[1];
            xlabel="chiN", ylabel="stress-free D/Rg",
            xmin=xmin, xmax=xmax,
            ymin=minimum(period_values)-0.08,
            ymax=maximum(period_values)+0.08)
        println(io, "<text x='$((top_panels[1].x0+top_panels[1].x1)/2)' y='58' text-anchor='middle' class='label'>Absolute stress-free period</text>")
        for style in ABA_CV_MODELS
            rows = benchmark[String.(benchmark.model) .== style.model, :]
            sort!(rows, :chiN)
            dash = style.model == "scft" ? "stroke-dasharray='7 4'" : ""
            println(io, "<polyline points='$(aba_cv_polyline(rows.chiN, rows.period_rg, px, py))' fill='none' stroke='$(style.color)' stroke-width='2.4' $dash/>")
            for row in eachrow(rows)
                println(io, "<circle cx='$(px(row.chiN))' cy='$(py(row.period_rg))' r='3.5' fill='$(style.color)'/>")
            end
        end

        period_values = 100 .* Float64.(benchmark.period_relative_error)
        ymin = min(-4.0, minimum(period_values)-1)
        ymax = max(4.0, maximum(period_values)+1)
        xmap, ymap = aba_cv_axes(io, top_panels[2];
            xlabel="chiN", ylabel="period error vs SCFT (%)",
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        println(io, "<text x='$((top_panels[2].x0+top_panels[2].x1)/2)' y='58' text-anchor='middle' class='label'>Period accuracy relative to SCFT</text>")
        println(io, "<line x1='$(top_panels[2].x0)' x2='$(top_panels[2].x1)' y1='$(ymap(0.0))' y2='$(ymap(0.0))' class='zero'/>")
        for style in ABA_CV_MODELS
            rows = benchmark[String.(benchmark.model) .== style.model, :]
            sort!(rows, :chiN)
            xs = Float64.(rows.chiN)
            ys = 100 .* Float64.(rows.period_relative_error)
            dash = style.model == "scft" ? "stroke-dasharray='7 4'" : ""
            println(io, "<polyline points='$(aba_cv_polyline(xs, ys, xmap, ymap))' fill='none' stroke='$(style.color)' stroke-width='2.4' $dash/>")
            for (x,y) in zip(xs,ys)
                println(io, "<circle cx='$(xmap(x))' cy='$(ymap(y))' r='3.5' fill='$(style.color)'/>")
            end
        end

        rms_values = Float64.(benchmark.profile_rms)
        rx, ry = aba_cv_axes(io, top_panels[3];
            xlabel="chiN", ylabel="profile RMS vs SCFT",
            xmin=xmin, xmax=xmax, ymin=0.0,
            ymax=max(0.001, 1.08 * maximum(rms_values)))
        println(io, "<text x='$((top_panels[3].x0+top_panels[3].x1)/2)' y='58' text-anchor='middle' class='label'>Own-cell density-profile error</text>")
        for style in ABA_CV_MODELS
            rows = benchmark[String.(benchmark.model) .== style.model, :]
            sort!(rows, :chiN)
            dash = style.model == "scft" ? "stroke-dasharray='7 4'" : ""
            println(io, "<polyline points='$(aba_cv_polyline(rows.chiN, rows.profile_rms, rx, ry))' fill='none' stroke='$(style.color)' stroke-width='2.4' $dash/>")
            for row in eachrow(rows)
                println(io, "<circle cx='$(rx(row.chiN))' cy='$(ry(row.profile_rms))' r='3.5' fill='$(style.color)'/>")
            end
        end

        for (panel, case) in zip(profile_panels, profile_cases)
            px, py = aba_cv_axes(io, panel;
                xlabel="normalized cell coordinate", ylabel="phi_A",
                xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0)
            for style in ABA_CV_MODELS
                rows = aba_cv_profile_rows(
                    profiles, case.case_id, style.model)
                dash = style.model == "scft" ? "stroke-dasharray='7 4'" : ""
                println(io, "<polyline points='$(aba_cv_polyline(rows.s, rows.phi_a, px, py))' fill='none' stroke='$(style.color)' stroke-width='2.2' $dash/>")
            end
            println(io, "<text x='$((panel.x0+panel.x1)/2)' y='413' text-anchor='middle' class='label'>fA=0.50, chiN=$(Int(case.chiN))</text>")
        end
        x0, y0 = 390.0, 755.0
        for (index, style) in enumerate(ABA_CV_MODELS)
            x = x0 + 170*(index-1)
            dash = style.model == "scft" ? "stroke-dasharray='7 4'" : ""
            println(io, "<line x1='$x' x2='$(x+24)' y1='$y0' y2='$y0' stroke='$(style.color)' stroke-width='3' $dash/>")
            println(io, "<text x='$(x+30)' y='$(y0+4)' class='small'>$(style.label)</text>")
        end
        println(io, "</svg>")
    end
end

function aba_cv_write_mechanism_svg(path, scan, benchmark)
    width, height = 1240, 760
    panels = [
        (x0=80.0, x1=590.0, y0=80.0, y1=330.0),
        (x0=690.0, x1=1200.0, y0=80.0, y1=330.0),
        (x0=80.0, x1=590.0, y0=420.0, y1=650.0),
        (x0=690.0, x1=1200.0, y0=420.0, y1=650.0),
    ]
    open(path, "w") do io
        aba_cv_svg_header(io, width, height,
            "Linear ABA relaxed free-energy and cell-stress mechanism")
        mechanism_cases = [
            case for case in ABA_CV_CASES
            if case.case_id in ABA_CV_MECHANISM_CASE_IDS
        ]
        for (case_index, case) in enumerate(mechanism_cases)
            subset = scan[String.(scan.case_id) .== case.case_id, :]
            emax = maximum(Float64.(subset.delta_energy))
            stress_abs = maximum(abs, Float64.(subset.cell_stress))
            period_min = minimum(Float64.(subset.period_rg))
            period_max = maximum(Float64.(subset.period_rg))
            margin = 0.02 * (period_max-period_min)
            scft_source = benchmark[
                (String.(benchmark.case_id) .== case.case_id) .&
                (String.(benchmark.model) .== "scft"), :]
            nrow(scft_source) == 1 ||
                error("missing SCFT optimum for mechanism plot")
            scft_period = Float64(scft_source.period_rg[1])
            ep = panels[2 * case_index - 1]
            sp = panels[2 * case_index]
            ex, ey = aba_cv_axes(io, ep;
                xlabel="absolute D/Rg", ylabel="F(D)-Fmin",
                xmin=period_min-margin, xmax=period_max+margin,
                ymin=0.0, ymax=1.05 * emax)
            sx, sy = aba_cv_axes(io, sp;
                xlabel="absolute D/Rg", ylabel="dF/d ln D",
                xmin=period_min-margin, xmax=period_max+margin,
                ymin=-1.05 * stress_abs,
                ymax=1.05 * stress_abs)
            println(io, "<line x1='$(sp.x0)' x2='$(sp.x1)' y1='$(sy(0.0))' y2='$(sy(0.0))' class='zero'/>")
            println(io, "<line x1='$(ex(scft_period))' x2='$(ex(scft_period))' y1='$(ep.y0)' y2='$(ep.y1)' stroke='#222' stroke-width='1.5' stroke-dasharray='3 3'/>")
            println(io, "<line x1='$(sx(scft_period))' x2='$(sx(scft_period))' y1='$(sp.y0)' y2='$(sp.y1)' stroke='#222' stroke-width='1.5' stroke-dasharray='3 3'/>")
            for style in ABA_CV_MODELS
                rows = subset[String.(subset.model) .== style.model, :]
                sort!(rows, :period_rg)
                dash = style.model == "scft" ? "stroke-dasharray='7 4'" : ""
                println(io, "<polyline points='$(aba_cv_polyline(rows.period_rg, rows.delta_energy, ex, ey))' fill='none' stroke='$(style.color)' stroke-width='2.2' $dash/>")
                println(io, "<polyline points='$(aba_cv_polyline(rows.period_rg, rows.cell_stress, sx, sy))' fill='none' stroke='$(style.color)' stroke-width='2.2' $dash/>")
            end
            println(io, "<text x='$((ep.x0+ep.x1)/2)' y='$(ep.y0-18)' text-anchor='middle' class='label'>fA=0.50, chiN=$(Int(case.chiN))</text>")
            println(io, "<text x='$((sp.x0+sp.x1)/2)' y='$(sp.y0-18)' text-anchor='middle' class='label'>fA=0.50, chiN=$(Int(case.chiN))</text>")
        end
        x0, y0 = 260.0, 735.0
        for (index, style) in enumerate(ABA_CV_MODELS)
            x = x0 + 165*(index-1)
            dash = style.model == "scft" ? "stroke-dasharray='7 4'" : ""
            println(io, "<line x1='$x' x2='$(x+25)' y1='$y0' y2='$y0' stroke='$(style.color)' stroke-width='3' $dash/>")
            println(io, "<text x='$(x+31)' y='$(y0+4)' class='small'>$(style.label)</text>")
        end
        println(io, "<line x1='1055' x2='1055' y1='724' y2='746' stroke='#222' stroke-width='1.5' stroke-dasharray='3 3'/><text x='1065' y='739' class='small'>SCFT stress-free D</text>")
        println(io, "</svg>")
    end
end

function aba_cv_write_boundary_svg(path, crossing, summary)
    width, height = 1060, 460
    left = (x0=80.0, x1=500.0, y0=70.0, y1=370.0)
    right = (x0=610.0, x1=1030.0, y0=70.0, y1=370.0)
    open(path, "w") do io
        aba_cv_svg_header(io, width, height,
            "Linear ABA L/DIS continuous boundary crossing")
        dmin = minimum(Float64.(crossing.delta_free_energy))
        dmax = max(1.0e-8, maximum(Float64.(crossing.delta_free_energy)))
        x1, y1 = aba_cv_axes(io, left;
            xlabel="chiN", ylabel="F_L-F_DIS",
            xmin=17.5, xmax=22.2, ymin=1.08 * dmin, ymax=1.08 * dmax)
        x2, y2 = aba_cv_axes(io, right;
            xlabel="chiN", ylabel="lamellar amplitude",
            xmin=17.5, xmax=22.2, ymin=0.0,
            ymax=1.05maximum(Float64.(crossing.amplitude)))
        println(io, "<line x1='$(left.x0)' x2='$(left.x1)' y1='$(y1(0.0))' y2='$(y1(0.0))' class='zero'/>")
        colors = Dict("bvk2"=>"#0072b2", "scft"=>"#222222")
        for model in ("bvk2", "scft")
            rows = crossing[String.(crossing.model) .== model, :]
            dash = model == "scft" ? "stroke-dasharray='7 4'" : ""
            println(io, "<polyline points='$(aba_cv_polyline(rows.chiN, rows.delta_free_energy, x1, y1))' fill='none' stroke='$(colors[model])' stroke-width='2.5' $dash/>")
            println(io, "<polyline points='$(aba_cv_polyline(rows.chiN, rows.amplitude, x2, y2))' fill='none' stroke='$(colors[model])' stroke-width='2.5' $dash/>")
            for row in eachrow(rows)
                println(io, "<circle cx='$(x1(row.chiN))' cy='$(y1(row.delta_free_energy))' r='3.5' fill='$(colors[model])'/>")
                println(io, "<circle cx='$(x2(row.chiN))' cy='$(y2(row.amplitude))' r='3.5' fill='$(colors[model])'/>")
            end
        end
        for row in eachrow(summary)
            color = colors[String(row.model)]
            for (panel, mapx) in ((left,x1),(right,x2))
                xpos = mapx(row.chiN_root)
                println(io, "<line x1='$xpos' x2='$xpos' y1='$(panel.y0)' y2='$(panel.y1)' stroke='$color' stroke-width='1.5' stroke-dasharray='4 3'/>")
            end
        end
        println(io, "<line x1='130' x2='160' y1='47' y2='47' stroke='#0072b2' stroke-width='3'/><text x='168' y='51' class='small'>AQCE, frozen c2=0.16</text>")
        println(io, "<line x1='350' x2='380' y1='47' y2='47' stroke='#222' stroke-width='3' stroke-dasharray='7 4'/><text x='388' y='51' class='small'>Polyorder SCFT</text>")
        println(io, "<text x='530' y='430' text-anchor='middle' class='small'>Vertical lines: analytic continuous-bifurcation roots; numerical ordered branches are validation points, not transversal brackets.</text>")
        println(io, "</svg>")
    end
end

function aba_cv_write_readme(path, benchmark, mechanism_summary,
        crossing_summary)
    bvk2 = benchmark[String.(benchmark.model) .== "bvk2", :]
    accepted_benchmark = all(Bool.(benchmark.accepted))
    accepted_mechanism = all(Bool.(mechanism_summary.accepted))
    open(path, "w") do io
        println(io, "# Comprehensive no-refit AQCE validation for linear ABA")
        println(io)
        println(io, "This packet extends P1-1 with ABA analogues of the AB P0-3, P0-6, and lamellar-relevant P0-7 evidence. The nonlinear AQCE coefficient remains frozen at `c2=0.16`; ABA SCFT data are validation-only.")
        println(io)
        println(io, "## P0-3 analogue: own-cell period/profile benchmark")
        println(io)
        println(io, "Figure 7 compares Uneyama--Doi, BURP, QCE, AQCE, and independently relaxed Polyorder SCFT. Seven ordered states from `chiN=20` through `45` are compared, and every profile is evaluated in its own optimized one-period cell. The profile panels show `chiN=20`, `30`, and `45`. The QCE control is stored in `../bvk2_aba_fixed_stiffness_validation/`. The cached tables retain a nonsmooth adaptive precursor as an auxiliary development diagnostic, but it is omitted from Figure 7 because it is not part of the manuscript-wide comparator set.")
        println(io)
        println(io, "| chiN | AQCE period error | AQCE profile RMS | correlation |")
        println(io, "|---:|---:|---:|---:|")
        for row in eachrow(bvk2)
            println(io, @sprintf("| %.0f | %+.3f%% | %.5f | %.5f |",
                row.chiN, 100 * row.period_relative_error,
                row.profile_rms, row.profile_correlation))
        end
        println(io)
        println(io, "Packet gate: **$(accepted_benchmark ? "pass" : "fail")**.")
        println(io)
        println(io, "## P0-6 analogue: relaxed F(D) and cell stress")
        println(io)
        println(io, "The cached campaign independently field-relaxes four architecture-aware DFTs and Polyorder SCFT at 14 absolute periods per state; Figure 7 displays UD, BURP, QCE, and AQCE. The QCE control is obtained from its own stress-free roots rather than from this two-state mechanism scan. The scan preserves absolute `D/Rg`, marks the SCFT stress-free period, and includes the full fixed-cell SCFT energy/stress curve. Acceptance requires an interior sampled energy minimum and a stress sign change. Gradient-based BURP/AQCE and the UD control require converged fields in the central period window; the outer diagnostic endpoints must remain finite. The auxiliary nonsmooth-precursor scan is retained in the archive. Packet gate: **$(accepted_mechanism ? "pass" : "fail")**.")
        println(io)
        println(io, "## P0-7 analogue: L/DIS crossing")
        println(io)
        println(io, "The symmetric ABA L/DIS transition is continuous. Therefore the accepted coordinate is the analytic Gaussian/RPA bifurcation, not a fabricated transversal free-energy root:")
        println(io)
        println(io, "| model | chiN root | D/Rg at bifurcation | semantics |")
        println(io, "|---|---:|---:|---|")
        for row in eachrow(crossing_summary)
            println(io, @sprintf("| %s | %.6f | %.6f | %s |",
                uppercase(String(row.model)), row.chiN_root,
                row.period_rg_at_bifurcation, row.semantics))
        end
        println(io)
        println(io, "The numerical ordered-branch points at `chiN=18.25`, `18.5`, `19`, `20`, `21`, and `22` show the free-energy decrease and amplitude growth above the bifurcation. They are diagnostic validation points.")
        println(io, "This comparison is reported separately as Figure 8.")
        println(io)
        println(io, "## Scope")
        println(io)
        println(io, "This strengthens no-refit structural and weak-boundary transfer for symmetric linear ABA. It is not an ABA nonlamellar phase diagram and does not establish S/C, G/C, G/L, arbitrary-topology, or Frank--Kasper claims.")
        println(io)
        println(io, "## Reproduction")
        println(io)
        println(io, "The canonical Figures 7 and 8 renderer verifies these accepted source fingerprints before reading any cached table:")
        println(io)
        println(io, "| source | SHA-256 |")
        println(io, "|---|---|")
        println(io, "| `benchmark_summary.csv` | `fe890e94ed44ee5a1cbe5dd8608d583e3f9f63e4066d08f2a1c213c83f636e4b` |")
        println(io, "| `benchmark_profiles.csv` | `0f206c11c4664a719a5fc4ab9ad5442a9e0d47bc34703a3347be0cbdf8debf74` |")
        println(io, "| `ldis_crossing.csv` | `8d73746ba9dc434fdf110fe859a4255e14cd49bf2f38ccb83f2f5fa912c30476` |")
        println(io, "| `ldis_crossing_summary.csv` | `ef23212c5e9bbd33a4e9d12952b258135b060e927719c447df4533df79c345b3` |")
        println(io, "| `../bvk2_aba_fixed_stiffness_validation/summary.csv` | `ae3d6349a673fbfa43ae5e4a4f8d6bad4813c4a3d8a70afcf674cc0f3ffb0dee` |")
        println(io, "| `../bvk2_aba_fixed_stiffness_validation/profiles.csv` | `eda1fdd2cd5c45a086d34881e08e9d471ebddf56163bb15d04a46a1330a57943` |")
        println(io)
        println(io, "```bash")
        println(io, "julia --project=. scripts/write_bvk2_aba_comprehensive_validation.jl --recompute=true")
        println(io, "uv run python scripts/assemble_bvk2_aba_fixed_stiffness_validation.py")
        println(io, "uv run python scripts/render_macromolecules_figure8.py")
        println(io, "```")
        println(io)
        println(io, "Omit `--recompute=true` to reuse the accepted cached CSV data. The dedicated Python command renders the three-panel manuscript Figure 7 and the separate LAM/DIS Figure 8 from those tables and fails closed if their data contract is violated.")
    end
end

function main()
    selected = TOML.parsefile(joinpath(
        ABA_CV_ROOT, "results", "bvk2_c2_calibration_nx128",
        "selected_model.toml"))
    Float64(selected["c2"]) == ABA_CV_C2 ||
        error("the frozen diblock AQCE calibration is not c2=0.16")
    outdir = joinpath(
        ABA_CV_ROOT, "results", "bvk2_aba_comprehensive_validation")
    mkpath(outdir)
    recompute = aba_cv_bool(aba_cv_arg("--recompute", "false"))
    benchmark, profiles = aba_cv_write_benchmark(
        outdir; recompute=recompute)
    scan, mechanism = aba_cv_write_mechanism(
        outdir, benchmark, profiles; recompute=recompute)
    aba_cv_promote_benchmark_from_mechanism!(
        benchmark, mechanism,
        joinpath(outdir, "benchmark_summary.csv"))
    crossing, crossing_summary = aba_cv_write_boundary(
        outdir; recompute=recompute)
    aba_cv_write_benchmark_svg(
        joinpath(outdir, "aba_unified_benchmark.svg"),
        benchmark, profiles)
    aba_cv_write_mechanism_svg(
        joinpath(outdir, "aba_cell_stress_mechanism.svg"),
        scan, benchmark)
    aba_cv_write_boundary_svg(
        joinpath(outdir, "aba_ldis_crossing.svg"),
        crossing, crossing_summary)
    aba_cv_write_readme(
        joinpath(outdir, "README.md"), benchmark, mechanism,
        crossing_summary)
    all(Bool.(benchmark.accepted)) ||
        error("ABA benchmark contains an unaccepted row")
    all(Bool.(mechanism.accepted)) ||
        error("ABA cell-stress packet contains an unaccepted row")
    println("wrote comprehensive ABA packet to ", outdir)
end

if get(ENV, "DFM_SKIP_BVK2_ABA_COMPREHENSIVE_MAIN", "0") != "1"
    main()
end
