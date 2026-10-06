const SCFT_AUDIT_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "write_liu2019_stress_free_period_map.jl"))

function audit_arg(name::AbstractString, default=nothing)
    prefix = String(name) * "="
    for argument in ARGS
        startswith(argument, prefix) && return split(argument, "="; limit=2)[2]
    end
    default === nothing && throw(ArgumentError("missing required $(name)=..."))
    return default
end

function audit_domain_count(profile::AbstractVector{<:Real}, threshold::Real)
    mask = Float64.(profile) .> Float64(threshold)
    (all(mask) || !any(mask)) && return 0
    domain_count = 0
    for index in eachindex(mask)
        previous = index == firstindex(mask) ? lastindex(mask) : index - 1
        domain_count += mask[index] && !mask[previous]
    end
    return domain_count
end

function audit_dominant_mode(profile::AbstractVector{<:Real})
    values = Float64.(profile) .- mean(profile)
    count_values = length(values)
    amplitudes = [abs(sum(values[index] *
        cis(-2pi * mode * (index - 1) / count_values)
        for index in eachindex(values))) for mode in 1:(count_values ÷ 2)]
    return argmax(amplitudes)
end

function audit_main()
    f = parse(Float64, audit_arg("--f"))
    chiN = parse(Float64, audit_arg("--chiN"))
    outdir = normpath(audit_arg("--outdir"))
    spacing = parse(Float64, audit_arg("--spacing", "0.05"))
    ds = parse(Float64, audit_arg("--ds", "0.005"))
    seed = parse(Int, audit_arg("--seed", "3083"))
    cell_block = parse(Int, audit_arg("--cell-block", "5"))
    cell_max_iterations = parse(Int,
        audit_arg("--cell-max-iterations", "8000"))
    mkpath(outdir)

    canonical_rows = DataFrame(CSV.File(joinpath(SCFT_AUDIT_ROOT, "results",
        "liu2019_stress_free_period_map", "summary.csv")))
    candidates = filter(row -> row.model == "scft" &&
        isapprox(Float64(row.f), f; atol=1.0e-12) &&
        isapprox(Float64(row.chiN), chiN; atol=1.0e-12),
        eachrow(canonical_rows))
    length(candidates) == 1 || error("expected one canonical SCFT row")
    canonical = only(candidates)
    canonical.status == "accepted" || error("canonical SCFT row is not accepted")
    reference_period = Float64(canonical.period_rg)

    attempt_rows = NamedTuple[]
    results = NamedTuple[]
    for factor in (0.98, 1.02)
        result = lpc_polyorder_cell_solved_profile(; f=f, chiN=chiN,
            initial_Lrg=reference_period * factor, spacing=spacing, ds=ds,
            seed=seed, presolve_max_iter=3000,
            cell_max_iter=cell_max_iterations,
            cell_block=cell_block, stress_tol=1.0e-5, cell_algorithm=:bb)
        domain_count = audit_domain_count(result.phi_a, f)
        dominant_mode = audit_dominant_mode(result.phi_a)
        accepted = _polyorder_success(result.cell_convergence) &&
            isfinite(result.free_energy) &&
            result.residual_norm < 1.0e-6 &&
            result.stress_norm <= 1.0e-5 &&
            result.incompressibility_rms <= 1.0e-6 &&
            abs(result.mean_phi_a - f) <= 1.0e-6 &&
            maximum(result.phi_a) - minimum(result.phi_a) > 1.0e-4 &&
            domain_count == 1 && dominant_mode == 1
        push!(results, result)
        push!(attempt_rows, (
            schema="scft-low-fa-resolution-audit-v1",
            claim_kind="stress_free_lamellar_scft_observable",
            f=f, chiN=chiN, initial_factor=factor,
            requested_spacing_rg=spacing,
            actual_spacing_rg=result.actual_spacing_rg,
            contour_step=ds, grid_count=result.grid_count,
            cell_block=cell_block,
            cell_max_iterations=cell_max_iterations,
            period_rg=result.period_rg, free_energy=result.free_energy,
            residual_norm=result.residual_norm,
            stress_norm=result.stress_norm,
            incompressibility_rms=result.incompressibility_rms,
            mean_phi_a=result.mean_phi_a,
            minimum_phi_a=minimum(result.phi_a),
            maximum_phi_a=maximum(result.phi_a),
            domain_count=domain_count, dominant_mode=dominant_mode,
            convergence=result.cell_convergence,
            accepted=accepted,
        ))
    end
    CSV.write(joinpath(outdir, "attempts.csv"), DataFrame(attempt_rows))

    accepted_indices = findall(row -> row.accepted, attempt_rows)
    isempty(accepted_indices) && error("no refined start passed the SCFT gates")
    selected_index = accepted_indices[argmin([
        results[index].free_energy for index in accepted_indices
    ])]
    selected = results[selected_index]
    period_spread = maximum(result.period_rg for result in results) /
        minimum(result.period_rg for result in results) - 1.0
    all_starts_accepted = all(row.accepted for row in attempt_rows)

    root = DataFrame([(
        schema="scft-low-fa-resolution-audit-v1",
        claim_kind="stress_free_lamellar_scft_observable",
        case_id=@sprintf("f%.4g_chiN%.4g", f, chiN),
        f=f, chiN=chiN, canonical_period_rg=reference_period,
        period_rg=selected.period_rg,
        requested_spacing_rg=spacing,
        actual_spacing_rg=selected.actual_spacing_rg,
        contour_step=ds, grid_count=selected.grid_count,
        cell_block=cell_block,
        cell_max_iterations=cell_max_iterations,
        free_energy=selected.free_energy,
        residual_norm=selected.residual_norm,
        stress_norm=selected.stress_norm,
        incompressibility_rms=selected.incompressibility_rms,
        mean_phi_a=selected.mean_phi_a,
        minimum_phi_a=minimum(selected.phi_a),
        maximum_phi_a=maximum(selected.phi_a),
        domain_count=audit_domain_count(selected.phi_a, f),
        dominant_mode=audit_dominant_mode(selected.phi_a),
        convergence=selected.cell_convergence,
        selected_initial_factor=attempt_rows[selected_index].initial_factor,
        all_starts_accepted=all_starts_accepted,
        period_spread_fraction=period_spread,
        accepted=all_starts_accepted && period_spread <= 5.0e-4,
    )])
    CSV.write(joinpath(outdir, "root.csv"), root)
    profile = DataFrame(index=collect(1:length(selected.phi_a)),
        s=(0:(length(selected.phi_a) - 1)) ./ length(selected.phi_a),
        phi_a=Float64.(selected.phi_a), phi_b=Float64.(selected.phi_b))
    CSV.write(joinpath(outdir, "profile.csv"), profile)
    root.accepted[1] || error("refined SCFT starts did not agree")
    println(@sprintf("accepted f=%.2f chiN=%.0f D/Rg=%.10f nx=%d residual=%.3g stress=%.3g",
        f, chiN, selected.period_rg, selected.grid_count,
        selected.residual_norm, selected.stress_norm))
end

audit_main()
