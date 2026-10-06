const DIBLOCK_BURP_BVK2_HYBRID_FORMULA_SCHEMA =
    "burp-ti-plus-bvk2-excess-stiffness-v1"
const DIBLOCK_BURP_BVK2_HYBRID_DISCRETIZATION_SCHEMA =
    "burp-ti-full-rpa-plus-central-eta-forward-edge-excess-v1"

"""
    diblock_burp_bvk2_hybrid_fingerprint(; f=0.5, N=1, b=1, c2,
                                          adaptive=true, nquad=16)

Return a stable fingerprint for the conservative BURP-TI+BVK2 hybrid.  The
hybrid retains the complete BURP-TI free energy and adds only the BVK2 excess
stiffness `Delta K = K_psi,2 - K_psi,0`; it never adds the full BVK2 gradient
term.  Consequently the BURP-TI/full-RPA quadratic response is not counted
twice.  `c2` remains required so a no-refit calculation cannot silently adopt
another coefficient.
"""
function diblock_burp_bvk2_hybrid_fingerprint(; f::Real=0.5, N::Real=1.0,
        b::Real=1.0, c2::Real, adaptive::Bool=true, nquad::Integer=16)
    params = diblock_bvk2_parameters(; f=f, N=N, b=b, c2=c2)
    quadrature = Int(nquad)
    quadrature >= 2 || throw(ArgumentError("nquad must be at least 2"))
    return @sprintf(
        "model=BURP-TI+BVK2|schema=%s|law_mode=%s|discretization=%s|c2=%.17g|f=%.17g|N=%.17g|b=%.17g|K_psi0=%.17g|eta_star2=%.17g|k_star=%.17g|nquad=%d|k_star_definition=%s|correction=Kpsi-Kpsi0",
        DIBLOCK_BURP_BVK2_HYBRID_FORMULA_SCHEMA,
        adaptive ? "adaptive_excess" : "zero_excess_control",
        DIBLOCK_BURP_BVK2_HYBRID_DISCRETIZATION_SCHEMA,
        params.c2, params.f, params.N, params.b, params.K_psi0,
        params.eta_star2, params.k_star, quadrature,
        params.k_star_definition)
end

function _hybrid_require_interior(phi)
    minimum(phi) > 0.0 && maximum(phi) < 1.0 || throw(DomainError(
        (minimum(phi), maximum(phi)),
        "BURP-TI+BVK2 hybrid densities must remain strictly inside (0, 1)"))
    return phi
end

function _diblock_burp_bvk2_excess_gradient_1d(
        phi::AbstractVector{Float64}, f::Float64, period::Float64,
        N::Float64, b::Float64; adaptive::Bool=true, c2::Real)
    adaptive || return 0.0
    params = diblock_bvk2_parameters(; f=f, N=N, b=b, c2=c2)
    kpsi = _diblock_bvk2_kpsi_1d(phi, f, period, N, b;
        adaptive=true, c2=params.c2, parameters=params)
    delta_k = kpsi .- params.K_psi0
    return _diblock_bvk2_gradient_integral_1d(
        phi, delta_k, period / length(phi))
end

function _diblock_burp_bvk2_excess_gradient_nd(
        phi::AbstractArray{Float64}, f::Float64,
        lengths::AbstractVector{Float64}, N::Float64, b::Float64,
        dV::Float64; adaptive::Bool=true, c2::Real)
    adaptive || return 0.0
    params = diblock_bvk2_parameters(; f=f, N=N, b=b, c2=c2)
    kpsi = _diblock_bvk2_kpsi_nd(phi, f, lengths, N, b;
        adaptive=true, c2=params.c2)
    delta_k = kpsi .- params.K_psi0
    return _diblock_bvk2_gradient_integral_nd(phi, delta_k, lengths, dV)
end

"""
    diblock_burp_bvk2_hybrid_energy_1d(phi; f, chiN, L, c2, ...)
    diblock_burp_bvk2_hybrid_energy_nd(phi; f, chiN, lengths, c2, ...)

Evaluate the conservative no-double-counting hybrid

`F_hyb = F_BURP-TI + integral (K_psi,2-K_psi,0) |grad(phi)|^2/[phi(1-phi)]`.

The correction vanishes identically when `c2=0` or `adaptive=false`.
"""
function diblock_burp_bvk2_hybrid_energy_1d(
        phi_a::AbstractVector{<:Real}; f::Real=0.5, chiN::Real=12.0,
        L::Real=1.0, N::Real=1.0, b::Real=1.0, c2::Real,
        adaptive::Bool=true, nquad::Integer=16,
        mean_tolerance::Real=1.0e-8, composition_kernel=nothing,
        workspace=nothing)
    phi, _count, ff, chi, period, nn, bb =
        _diblock_profile_energy_inputs(phi_a; f=f, chiN=chiN, L=L,
            N=N, b=b, mean_tolerance=mean_tolerance)
    _hybrid_require_interior(phi)
    base = diblock_burp_ti_energy_1d(phi; f=ff, chiN=chi, L=period,
        N=nn, b=bb, nquad=nquad, mean_tolerance=mean_tolerance,
        composition_kernel=composition_kernel, workspace=workspace)
    return base + _diblock_burp_bvk2_excess_gradient_1d(
        phi, ff, period, nn, bb; adaptive=adaptive, c2=c2)
end

function diblock_burp_bvk2_hybrid_energy_nd(
        phi_a::AbstractArray{<:Real}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, N::Real=1.0, b::Real=1.0, c2::Real,
        adaptive::Bool=true, nquad::Integer=16,
        mean_tolerance::Real=1.0e-8, composition_kernel=nothing,
        workspace=nothing)
    phi, _dims, ff, chi, cell_lengths, nn, bb, dV =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    _hybrid_require_interior(phi)
    base = diblock_burp_ti_energy_nd(phi; f=ff, chiN=chi,
        lengths=cell_lengths, N=nn, b=bb, nquad=nquad,
        mean_tolerance=mean_tolerance, composition_kernel=composition_kernel,
        workspace=workspace)
    return base + _diblock_burp_bvk2_excess_gradient_nd(
        phi, ff, cell_lengths, nn, bb, dV; adaptive=adaptive, c2=c2)
end

function _diblock_burp_bvk2_excess_chemical_potential_nd(
        phi::AbstractArray{Float64}, lengths::AbstractVector{Float64},
        params; adaptive::Bool=true)
    adaptive || return zeros(Float64, size(phi))
    adaptive_gradient = _diblock_bvk2_gradient_chemical_potential_nd(
        phi, lengths, params; adaptive=true)
    fixed_gradient = _diblock_bvk2_gradient_chemical_potential_nd(
        phi, lengths, params; adaptive=false)
    adaptive_gradient .-= fixed_gradient
    return adaptive_gradient
end

"""
    diblock_burp_bvk2_hybrid_chemical_potential_1d(phi; f, chiN, L, c2, ...)
    diblock_burp_bvk2_hybrid_chemical_potential_nd(phi; f, chiN, lengths, c2, ...)

Return the exact discrete chemical potential of the implemented hybrid:
the BURP-TI derivative plus the adaptive-minus-fixed BVK2 gradient reverse
pass.  The mean constraint is not projected.
"""
function diblock_burp_bvk2_hybrid_chemical_potential_1d(
        phi_a::AbstractVector{<:Real}; f::Real=0.5, chiN::Real=12.0,
        L::Real=1.0, N::Real=1.0, b::Real=1.0, c2::Real,
        adaptive::Bool=true, nquad::Integer=16,
        mean_tolerance::Real=1.0e-8, composition_kernel=nothing,
        workspace=nothing)
    phi, _count, ff, chi, period, nn, bb =
        _diblock_profile_energy_inputs(phi_a; f=f, chiN=chiN, L=L,
            N=N, b=b, mean_tolerance=mean_tolerance)
    _hybrid_require_interior(phi)
    base = _diblock_burp_ti_chemical_potential_1d(phi; f=ff, chiN=chi,
        L=period, N=nn, b=bb, nquad=nquad,
        composition_kernel=composition_kernel, workspace=workspace)
    params = diblock_bvk2_parameters(; f=ff, N=nn, b=bb, c2=c2)
    base .+= _diblock_burp_bvk2_excess_chemical_potential_nd(
        phi, [period], params; adaptive=adaptive)
    return base
end

function diblock_burp_bvk2_hybrid_chemical_potential_nd(
        phi_a::AbstractArray{<:Real}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, N::Real=1.0, b::Real=1.0, c2::Real,
        adaptive::Bool=true, nquad::Integer=16,
        mean_tolerance::Real=1.0e-8, composition_kernel=nothing,
        workspace=nothing)
    phi, _dims, ff, chi, cell_lengths, nn, bb, _dV =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    _hybrid_require_interior(phi)
    base = diblock_burp_ti_chemical_potential_nd(phi; f=ff, chiN=chi,
        lengths=cell_lengths, N=nn, b=bb, nquad=nquad,
        mean_tolerance=mean_tolerance, composition_kernel=composition_kernel,
        workspace=workspace)
    params = diblock_bvk2_parameters(; f=ff, N=nn, b=bb, c2=c2)
    base .+= _diblock_burp_bvk2_excess_chemical_potential_nd(
        phi, cell_lengths, params; adaptive=adaptive)
    return base
end

"""
    diblock_burp_bvk2_hybrid_cell_scale_gradient_density_nd(phi; ...)

Return `d(F/V)/ds` at `s=1` when every cell length is changed to
`s*lengths[d]` at fixed nodal density values.  BURP-TI's derivative with
respect to its scalar base length is converted to this dimensionless scale
derivative before the adaptive-minus-fixed BVK2 gradient contribution is
added. When the energy is evaluated with a custom `composition_kernel`, pass
the matching scalar `kernel_derivative(k2)` here so the BURP-TI stress remains
the derivative of that same kernel.
"""
function diblock_burp_bvk2_hybrid_cell_scale_gradient_density_nd(
        phi_a::AbstractArray{<:Real}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, N::Real=1.0, b::Real=1.0, c2::Real,
        adaptive::Bool=true, nquad::Integer=16,
        mean_tolerance::Real=1.0e-8, kernel_derivative=nothing)
    phi, _dims, ff, chi, cell_lengths, nn, bb, _dV =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    _hybrid_require_interior(phi)
    burp_da = diblock_burp_ti_cell_scale_gradient_density_nd(phi; f=ff,
        chiN=chi, lengths=cell_lengths, N=nn, b=bb, nquad=nquad,
        mean_tolerance=mean_tolerance,
        kernel_derivative=kernel_derivative)
    burp_ds = cell_lengths[1] * burp_da
    adaptive || return burp_ds
    adaptive_stress = diblock_bvk2_cell_scale_gradient_density_nd(phi;
        f=ff, chiN=chi, lengths=cell_lengths, N=nn, b=bb, adaptive=true,
        c2=c2, mean_tolerance=mean_tolerance)
    fixed_stress = diblock_bvk2_cell_scale_gradient_density_nd(phi;
        f=ff, chiN=chi, lengths=cell_lengths, N=nn, b=bb, adaptive=false,
        c2=c2, mean_tolerance=mean_tolerance)
    return burp_ds + adaptive_stress - fixed_stress
end

"""
    minimize_diblock_burp_bvk2_hybrid_lamella(; f, chiN, L, c2, ...)

Relax a one-dimensional hybrid lamella in the positivity-preserving theta
chart used by the existing conservative density-functional solvers.
"""
function minimize_diblock_burp_bvk2_hybrid_lamella(; f::Real=0.5,
        chiN::Real=12.0, L=nothing, nx::Integer=128,
        mode_count::Integer=3, initial_amplitude::Real=0.2,
        initial_profile=nothing,
        N::Real=1.0, b::Real=1.0, c2::Real, adaptive::Bool=true,
        nquad::Integer=16, max_iterations::Integer=10_000,
        lbfgs_memory::Integer=10, gradient_tolerance::Real=1.0e-10,
        force_tolerance::Real=1.0e-5,
        force_maxabs_tolerance::Real=2.0e-5,
        newton_polish_iterations::Integer=12,
        newton_difference_step::Real=1.0e-5,
        newton_trust_radius::Real=1.0)
    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    period = L === nothing ?
        2.0 * pi / sqrt(_diblock_kernel_minimum_k2(; f=ff, N=nn, b=bb)) :
        Float64(L)
    count = Int(nx)
    kernel = [diblock_composition_kernel(k2; f=ff, N=nn, b=bb)
        for k2 in _periodic_k2_1d(count, period)]
    energy_evaluator = (field, _params, f_value, chi_value, period_value,
            n_value, b_value, adaptive_value) ->
        diblock_burp_bvk2_hybrid_energy_1d(field; f=f_value,
            chiN=chi_value, L=period_value, N=n_value, b=b_value, c2=c2,
            adaptive=adaptive_value, nquad=nquad, composition_kernel=kernel)
    chemical_potential_evaluator = (field, _params, f_value, chi_value,
            period_value, n_value, b_value, adaptive_value) ->
        diblock_burp_bvk2_hybrid_chemical_potential_1d(field; f=f_value,
            chiN=chi_value, L=period_value, N=n_value, b=b_value, c2=c2,
            adaptive=adaptive_value, nquad=nquad, composition_kernel=kernel)
    return _minimize_diblock_bvk2_logit_lamella(; f=ff, chiN=chiN,
        L=period, nx=count, mode_count=mode_count,
        initial_amplitude=initial_amplitude, initial_profile=initial_profile,
        N=nn, b=bb, adaptive=adaptive, c2=c2,
        energy_evaluator=energy_evaluator,
        chemical_potential_evaluator=chemical_potential_evaluator,
        max_iterations=max_iterations, lbfgs_memory=lbfgs_memory,
        gradient_tolerance=gradient_tolerance,
        force_tolerance=force_tolerance,
        force_maxabs_tolerance=force_maxabs_tolerance,
        newton_polish_iterations=newton_polish_iterations,
        newton_difference_step=newton_difference_step,
        newton_trust_radius=newton_trust_radius)
end

"""
    minimize_diblock_burp_bvk2_hybrid_lamella_stress_free(; f, chiN, c2, ...)

Minimize the hybrid lamellar energy density over the period without changing
the supplied BVK2 coefficient.
"""
function minimize_diblock_burp_bvk2_hybrid_lamella_stress_free(;
        f::Real=0.5, chiN::Real=12.0, nx::Integer=128,
        mode_count::Integer=3, initial_amplitudes=(0.2, 0.5, 1.0),
        initial_profile=nothing,
        N::Real=1.0, b::Real=1.0, c2::Real, adaptive::Bool=true,
        nquad::Integer=16, max_iterations::Integer=10_000,
        max_period_iterations::Integer=24, lower_factor::Real=0.75,
        upper_factor::Real=1.5, lbfgs_memory::Integer=10,
        gradient_tolerance::Real=1.0e-10,
        force_tolerance::Real=1.0e-5,
        force_maxabs_tolerance::Real=2.0e-5,
        newton_polish_iterations::Integer=12,
        newton_difference_step::Real=1.0e-5,
        newton_trust_radius::Real=1.0, bootstrap_period_factors=(),
        bootstrap_window::Real=0.12, local_check_fraction::Real=0.01,
        local_check_tolerance::Real=1.0e-8,
        cell_stress_tolerance::Real=1.0e-6,
        log_period_tolerance::Real=1.0e-5,
        progress_label::AbstractString="")
    amplitudes = _checked_initial_amplitudes(initial_amplitudes)
    lo = _checked_positive_factor(lower_factor, "lower_factor")
    hi = _checked_positive_factor(upper_factor, "upper_factor")
    lo < hi || throw(ArgumentError(
        "lower_factor must be smaller than upper_factor"))
    period_iterations = Int(max_period_iterations)
    period_iterations > 0 || throw(ArgumentError(
        "max_period_iterations must be positive"))
    window = Float64(bootstrap_window)
    0.0 < window < 1.0 && isfinite(window) || throw(ArgumentError(
        "bootstrap_window must be finite and between zero and one"))
    check_fraction = Float64(local_check_fraction)
    0.0 < check_fraction < 1.0 && isfinite(check_fraction) ||
        throw(ArgumentError(
            "local_check_fraction must be finite and between zero and one"))
    check_tolerance = Float64(local_check_tolerance)
    check_tolerance >= 0.0 && isfinite(check_tolerance) ||
        throw(ArgumentError(
            "local_check_tolerance must be nonnegative and finite"))
    stress_tolerance = Float64(cell_stress_tolerance)
    stress_tolerance > 0.0 && isfinite(stress_tolerance) ||
        throw(ArgumentError(
            "cell_stress_tolerance must be positive and finite"))
    log_tolerance = Float64(log_period_tolerance)
    log_tolerance > 0.0 && isfinite(log_tolerance) ||
        throw(ArgumentError(
            "log_period_tolerance must be positive and finite"))

    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    chi = Float64(chiN)
    isfinite(chi) || throw(ArgumentError("chiN must be finite"))
    reference_period =
        2.0 * pi / sqrt(_diblock_kernel_minimum_k2(; f=ff, N=nn, b=bb))
    rg_scale = sqrt(6.0 / (nn * bb^2))
    bootstrap_factors = _checked_period_factor_list(
        bootstrap_period_factors, "bootstrap_period_factors")
    center_factor = isempty(bootstrap_factors) ? 1.0 : bootstrap_factors[end]
    center_factor = clamp(center_factor, lo, hi)
    evaluations = Dict{Float64,NamedTuple}()

    function evaluate_factor(factor::Real)
        value = clamp(Float64(factor), lo, hi)
        haskey(evaluations, value) && return evaluations[value]
        period = reference_period * value
        accepted_rows = [row for row in values(evaluations) if row.valid]
        warm_profile = isempty(accepted_rows) ? nothing :
            accepted_rows[argmin(abs(log(row.factor / value))
                for row in accepted_rows)].result.phi_a
        attempts = Tuple{Union{Nothing,Vector{Float64}},Float64}[]
        warm_profile !== nothing && push!(attempts,
            (Vector{Float64}(warm_profile), amplitudes[1]))
        initial_profile !== nothing && push!(attempts,
            (Float64.(collect(initial_profile)), amplitudes[1]))
        append!(attempts, [(nothing, amplitude) for amplitude in amplitudes])
        best = nothing
        attempt_errors = String[]
        designated_center_profile = initial_profile !== nothing &&
            isempty(evaluations)
        for (attempt_index, (profile, amplitude)) in enumerate(attempts)
            result = try
                minimize_diblock_burp_bvk2_hybrid_lamella(; f=ff,
                    chiN=chi, L=period, nx=nx, mode_count=mode_count,
                    initial_amplitude=amplitude, initial_profile=profile,
                    N=nn, b=bb, c2=c2, adaptive=adaptive, nquad=nquad,
                    max_iterations=max_iterations,
                    lbfgs_memory=lbfgs_memory,
                    gradient_tolerance=gradient_tolerance,
                    force_tolerance=force_tolerance,
                    force_maxabs_tolerance=force_maxabs_tolerance,
                    newton_polish_iterations=newton_polish_iterations,
                    newton_difference_step=newton_difference_step,
                    newton_trust_radius=newton_trust_radius)
            catch error
                message = sprint(showerror, error)
                push!(attempt_errors, message)
                _period_progress(progress_label, @sprintf(
                    "hybrid fixed-period attempt %d failed at factor %.8g: %s",
                    attempt_index, value, message))
                designated_center_profile && attempt_index == 1 &&
                    throw(ErrorException(
                        "designated continuation seed failed: " * message))
                nothing
            end
            result === nothing && continue
            identity_pass =
                result.energy < result.homogeneous_energy - 1.0e-10 &&
                result.maximum_phi - result.minimum_phi > 1.0e-4
            valid = result.converged && identity_pass && isfinite(result.energy)
            stress = valid ?
                diblock_burp_bvk2_hybrid_cell_scale_gradient_density_nd(
                    result.phi_a; f=ff, chiN=chi, lengths=(period,),
                    N=nn, b=bb, c2=c2, adaptive=adaptive,
                    nquad=nquad) : NaN
            valid &= isfinite(stress)
            row = (factor=value, log_factor=log(value), period=period,
                objective=result.energy / period, stress=stress,
                result=result, amplitude=amplitude, valid=valid,
                identity_pass=identity_pass)
            designated_center_profile && attempt_index == 1 && !row.valid &&
                throw(ErrorException(
                    "designated continuation seed did not reach a stationary lamella"))
            if best === nothing || (row.valid && !best.valid) ||
                    (row.valid == best.valid && row.objective < best.objective)
                best = row
            end
            row.valid && break
        end
        best === nothing && throw(ErrorException(@sprintf(
            "all hybrid fixed-period attempts failed at factor %.8g: %s",
            value, join(unique(attempt_errors), " | "))))
        evaluations[value] = best
        _period_progress(progress_label, @sprintf(
            "hybrid cell root factor=%.8g L/Rg=%.8g stress=%.8g field=%s",
            value, period * rg_scale, best.stress,
            string(best.result.converged)))
        return best
    end

    center = evaluate_factor(center_factor)
    center.valid || throw(ErrorException(
        "hybrid bootstrap period did not produce a stationary lamella"))
    direction = center.stress > 0.0 ? -1.0 : 1.0
    search_step = min(0.02, window / 4.0)
    current = center
    bracket = nothing
    for _ in 1:60
        candidate_factor = current.factor * exp(direction * search_step)
        lo < candidate_factor < hi || break
        trial = evaluate_factor(candidate_factor)
        if !trial.valid
            search_step *= 0.5
            search_step >= 1.0e-4 || break
            continue
        end
        if current.stress < 0.0 < trial.stress
            bracket = (current, trial)
            break
        elseif trial.stress < 0.0 < current.stress
            bracket = (trial, current)
            break
        end
        current = trial
    end
    bracket === nothing && throw(ErrorException(
        "hybrid analytic cell stress has no valid oriented continuation bracket"))
    lower, upper = bracket
    selected = abs(lower.stress) <= abs(upper.stress) ? lower : upper
    for _ in 1:period_iterations
        width = upper.log_factor - lower.log_factor
        abs(selected.stress) <= stress_tolerance && width <= log_tolerance &&
            break
        secant = (lower.log_factor * upper.stress -
            upper.log_factor * lower.stress) /
            (upper.stress - lower.stress)
        guard = 0.1 * width
        trial_log = clamp(secant, lower.log_factor + guard,
            upper.log_factor - guard)
        trial = evaluate_factor(exp(trial_log))
        trial.valid || break
        abs(trial.stress) < abs(selected.stress) && (selected = trial)
        if trial.stress < 0.0
            lower = trial
        else
            upper = trial
        end
    end
    width_pass = upper.log_factor - lower.log_factor <= log_tolerance
    stress_pass = abs(selected.stress) <= stress_tolerance

    left_factor = max(lo, selected.factor * (1.0 - check_fraction))
    right_factor = min(hi, selected.factor * (1.0 + check_fraction))
    left = evaluate_factor(left_factor)
    right = evaluate_factor(right_factor)
    boundary_limited = left_factor == lo || right_factor == hi
    orientation_pass = left.valid && right.valid &&
        left.stress < 0.0 < right.stress
    energy_pass = left.valid && right.valid &&
        left.objective >= selected.objective - check_tolerance &&
        right.objective >= selected.objective - check_tolerance
    local_pass = !boundary_limited && selected.valid && stress_pass &&
        width_pass && orientation_pass && energy_pass

    return DiblockLamellaPeriodOptimizationResult("BURP-TI+BVK2",
        selected.result, reference_period, selected.factor, selected.period,
        selected.period * rg_scale, selected.objective,
        "AnalyticCellStressRoot", local_pass, length(evaluations),
        selected.amplitude, lo, hi, !isempty(bootstrap_factors),
        center_factor, window, selected.factor, selected.objective,
        left.factor, left.objective, right.factor, right.objective,
        check_fraction, local_pass, boundary_limited)
end
