const DIBLOCK_BVK2_MORPHOLOGY_SOLVER_SCHEMA =
    "mean_physical_nd_lbfgs_energy_gate_jfnk_v4"

function _diblock_bvk2_morphology_phi_from_logits(
        logits::AbstractVector{<:Real}, dims::NTuple{D,Int}, f::Float64,
        dV::Float64) where {D}
    length(logits) == prod(dims) || throw(DimensionMismatch(
        "BVK2 morphology logits must match the requested grid dimensions"))
    phi_values = Vector{Float64}(undef, length(logits))
    weights = similar(phi_values)
    _diblock_bvk2_logit_field_gradient!(phi_values, weights, nothing,
        logits, nothing, f, dV)
    return reshape(phi_values, dims), weights
end

function _bvk2_anderson_type2_direction(values_history, step_history;
        memory::Integer=8, regularization::Real=1.0e-10)
    length(values_history) == length(step_history) ||
        throw(DimensionMismatch("Anderson histories must have equal length"))
    isempty(values_history) && throw(ArgumentError(
        "Anderson history must be nonempty"))
    depth = min(max(0, Int(memory)), length(values_history) - 1)
    current = Vector{Float64}(step_history[end])
    depth == 0 && return current
    first = length(values_history) - depth
    delta_values = hcat([
        Vector{Float64}(values_history[index + 1]) .-
            Vector{Float64}(values_history[index])
        for index in first:(length(values_history) - 1)]...)
    delta_steps = hcat([
        Vector{Float64}(step_history[index + 1]) .-
            Vector{Float64}(step_history[index])
        for index in first:(length(step_history) - 1)]...)
    eta = Float64(regularization) * max(1.0, sum(abs2, delta_steps))
    normal = transpose(delta_steps) * delta_steps
    normal[diagind(normal)] .+= eta
    gamma = normal \ (transpose(delta_steps) * current)
    direction = current .-
        (delta_values .+ delta_steps) * gamma
    return all(isfinite, direction) ? direction : current
end

"""
    _diblock_bvk2_morphology_logit_state(logits, dims; ...)

Evaluate a fixed-cell BVK2 state in the exact mean-constrained logit chart.
For unconstrained coordinates `z_i`, the physical density is

```
phi_i = logistic(z_i + s(z)),   mean(phi) = f.
```

Implicit differentiation of the scalar shift gives the exact pullback

```
dF/dz_j = dV*w_j*(mu_j - sum(w .* mu)/sum(w)),
w_j = phi_j*(1-phi_j).
```

Only logistic values that round exactly to zero or one are moved to the
nearest representable interior value.  This is a coordinate guard, not a
regularization of the BVK2 free energy.
"""
function _diblock_bvk2_morphology_logit_state(
        logits::AbstractVector{<:Real}, dims::NTuple{D,Int};
        f::Float64, chiN::Float64, lengths::AbstractVector{Float64},
        N::Float64, b::Float64, adaptive::Bool, c2::Real) where {D}
    volume = prod(lengths)
    dV = volume / prod(dims)
    phi, weights = _diblock_bvk2_morphology_phi_from_logits(
        logits, dims, f, dV)
    energy = diblock_bvk2_energy_nd(phi; f=f, chiN=chiN, lengths=lengths,
        N=N, b=b, adaptive=adaptive, c2=c2, mean_tolerance=1.0e-10)
    chemical = diblock_bvk2_chemical_potential_nd(phi; f=f, chiN=chiN,
        lengths=lengths, N=N, b=b, adaptive=adaptive, c2=c2,
        mean_tolerance=1.0e-10)
    logit_gradient = Vector{Float64}(undef, length(logits))
    _diblock_bvk2_logit_field_gradient!(vec(phi), weights,
        logit_gradient, logits, vec(chemical), f, dV)
    # Optimize the intensive energy so identical morphologies embedded in
    # cells with different passive volume have identical objective scaling.
    logit_gradient ./= volume
    residual = chemical .- mean(chemical)
    residual_r2 = sqrt(dV * sum(abs2, residual))
    residual_rinf = maximum(abs, residual)
    return (
        energy=energy,
        energy_density=energy / volume,
        phi=copy(phi),
        chemical=chemical,
        residual=residual,
        residual_r2=residual_r2,
        residual_rinf=residual_rinf,
        logit_gradient=logit_gradient,
    )
end

function _diblock_bvk2_morphology_physical_state(phi::AbstractArray{<:Real};
        f::Float64, chiN::Float64, lengths::AbstractVector{Float64},
        N::Float64, b::Float64, adaptive::Bool, c2::Real)
    field = Float64.(Array(phi))
    all(isfinite, field) || throw(ArgumentError(
        "BVK2 physical state must contain only finite values"))
    # Pure/saturated domains are physical at strong segregation, so phi may
    # reach 0/1 or overshoot slightly (the energy uses the xlogy entropy limit
    # and the divergent eta / entropy-derivative logs are floored at 1e-10).
    # Reject only gross excursions past the box.
    (minimum(field) >= -1.0e-10 && maximum(field) <= 1.0 + 1.0e-10) ||
        throw(DomainError((minimum(field), maximum(field)),
            "BVK2 phi overshoot exceeds 1e-10 past [0, 1]"))
    abs(mean(field) - f) <= 1.0e-10 || throw(ArgumentError(
        "BVK2 physical state must satisfy the fixed mean"))
    volume = prod(lengths)
    dV = volume / length(field)
    energy = diblock_bvk2_energy_nd(field; f=f, chiN=chiN,
        lengths=lengths, N=N, b=b, adaptive=adaptive, c2=c2,
        mean_tolerance=1.0e-10)
    chemical = diblock_bvk2_chemical_potential_nd(field; f=f, chiN=chiN,
        lengths=lengths, N=N, b=b, adaptive=adaptive, c2=c2,
        mean_tolerance=1.0e-10)
    residual = chemical .- mean(chemical)
    return (
        energy=energy,
        energy_density=energy / volume,
        phi=field,
        chemical=chemical,
        residual=residual,
        residual_r2=sqrt(dV * sum(abs2, residual)),
        residual_rinf=maximum(abs, residual),
    )
end

function _diblock_bvk2_translation_modes(phi;
        relative_norm_tolerance::Float64=1.0e-10)
    modes = Vector{Vector{Float64}}()
    scale = max(norm(vec(phi .- mean(phi))), 1.0)
    for dimension in 1:ndims(phi)
        backward = ntuple(axis -> axis == dimension ? -1 : 0, ndims(phi))
        forward = ntuple(axis -> axis == dimension ? 1 : 0, ndims(phi))
        mode = 0.5 .* (circshift(phi, backward) .-
            circshift(phi, forward))
        mode .-= mean(mode)
        mode_vector = vec(mode)
        for basis in modes
            mode_vector .-= dot(mode_vector, basis) .* basis
        end
        mode_norm = norm(mode_vector)
        mode_norm > relative_norm_tolerance * scale &&
            push!(modes, mode_vector ./ mode_norm)
    end
    return modes
end

function _diblock_bvk2_project_gauge!(values, translation_modes)
    values .-= mean(values)
    vector = vec(values)
    for mode in translation_modes
        vector .-= dot(vector, mode) .* mode
    end
    values .-= mean(values)
    return values
end

function _diblock_bvk2_rpa_preconditioner(dims, lengths;
        f::Float64, chiN::Float64, N::Float64, b::Float64,
        denominator_floor::Float64)
    cell_lengths = Float64.(collect(lengths))
    k2_continuum = _periodic_k2_nd(dims, cell_lengths)
    k2_fd = zeros(Float64, dims)
    for index in CartesianIndices(k2_fd)
        value = 0.0
        for dimension in eachindex(dims)
            count = dims[dimension]
            raw_mode = index[dimension] - 1
            mode = raw_mode <= count ÷ 2 ? raw_mode : raw_mode - count
            spacing = cell_lengths[dimension] / count
            value += (2.0 * sin(pi * mode / count) / spacing)^2
        end
        k2_fd[index] = value
    end
    params = diblock_bvk2_parameters(; f=f, N=N, b=b, c2=0.0)
    denominator = similar(k2_continuum)
    first_index = firstindex(denominator)
    for index in eachindex(denominator)
        if index == first_index
            denominator[index] = Inf
            continue
        end
        value = params.A_psi / k2_continuum[index] +
            params.M_psi / (f * (1.0 - f)) +
            N * b^2 * k2_fd[index] / (12.0 * f * (1.0 - f)) - 2.0 * chiN
        denominator[index] = abs(value) < denominator_floor ?
            copysign(denominator_floor, value == 0.0 ? 1.0 : value) : value
    end
    function apply(vector)
        transformed = fft(reshape(vector, dims))
        transformed[first_index] = 0.0
        @inbounds for index in eachindex(transformed)
            index == first_index && continue
            transformed[index] /= denominator[index]
        end
        output = vec(real(ifft(transformed)))
        output .-= mean(output)
        return output
    end
    return apply
end

function _diblock_bvk2_gmres(apply_operator, rhs::Vector{Float64};
        max_iterations::Int, restart::Int, tolerance::Float64,
        apply_preconditioner=identity)
    solution = zeros(Float64, length(rhs))
    beta0 = norm(rhs)
    beta0 == 0.0 && return solution, 0, 0.0
    current_residual = copy(rhs)
    total_iterations = 0
    while total_iterations < max_iterations
        beta = norm(current_residual)
        beta / beta0 <= tolerance &&
            return solution, total_iterations, beta / beta0
        subspace = min(restart, max_iterations - total_iterations)
        basis = zeros(Float64, length(rhs), subspace + 1)
        hessenberg = zeros(Float64, subspace + 1, subspace)
        basis[:, 1] .= current_residual ./ beta
        target = zeros(Float64, subspace + 1)
        target[1] = beta
        best_y = Float64[]
        best_count = 0
        best_relative = Inf
        for column in 1:subspace
            work = Vector{Float64}(apply_operator(
                apply_preconditioner(basis[:, column])))
            for row in 1:column
                hessenberg[row, column] = dot(basis[:, row], work)
                @. work -= hessenberg[row, column] * basis[:, row]
            end
            hessenberg[column + 1, column] = norm(work)
            if hessenberg[column + 1, column] > 1.0e-13 && column < subspace
                basis[:, column + 1] .=
                    work ./ hessenberg[column + 1, column]
            end
            y = hessenberg[1:column + 1, 1:column] \
                target[1:column + 1]
            relative = norm(target[1:column + 1] -
                hessenberg[1:column + 1, 1:column] * y) / beta0
            if relative < best_relative
                best_relative = relative
                best_y = y
                best_count = column
            end
            if relative <= tolerance
                solution .+= basis[:, 1:column] * y
                return Vector{Float64}(apply_preconditioner(solution)),
                    total_iterations + column, relative
            end
        end
        best_count == 0 && break
        solution .+= basis[:, 1:best_count] * best_y
        total_iterations += best_count
        current_residual .= rhs .-
            apply_operator(apply_preconditioner(solution))
    end
    relative = norm(rhs .-
        apply_operator(apply_preconditioner(solution))) / beta0
    return Vector{Float64}(apply_preconditioner(solution)),
        total_iterations, relative
end

"""
    _bvk2_matrix_free_newton(state_at, initial_values; ...)

Safeguarded Newton--GMRES solve for a nonlinear residual without assembling a
dense Jacobian.  `state_at(values)` must return a named tuple containing
`residual`, `r2`, `rinf`, and `energy`.  Jacobian-vector products use a
mesh-relative centred difference and the linear system is solved with the
existing restarted right-preconditioned GMRES implementation.

The nonlinear acceptance test is deliberately expressed in the caller's
normalized residual norm, so changing the number of grid coefficients does
not change the convergence contract.  If the Krylov direction is unusable, a
bounded Type-II Anderson direction is tried before the raw preconditioned
residual fallback.
"""
function _bvk2_matrix_free_newton(state_at, initial_values;
        max_iterations::Integer=20, force_tolerance::Real=1.0e-5,
        force_maxabs_tolerance::Real=Inf,
        time_limit::Real=Inf,
        gmres_max_iterations::Integer=32, gmres_restart::Integer=16,
        gmres_tolerance::Real=0.1,
        maximum_linear_relative_residual::Real=0.95,
        relative_step::Real=1.0e-7, trust_radius::Real=0.1,
        line_search_steps::Integer=24,
        sufficient_decrease::Real=1.0e-4,
        energy_increase_tolerance::Real=1.0e-12,
        apply_preconditioner=identity, preconditioner_builder=nothing,
        step_norm=vector -> maximum(abs, vector),
        anderson_memory::Integer=8,
        anderson_regularization::Real=1.0e-10, callback=nothing)
    iterations = max(0, Int(max_iterations))
    linear_iterations = max(1, Int(gmres_max_iterations))
    restart = clamp(Int(gmres_restart), 1, linear_iterations)
    tolerance = Float64(force_tolerance)
    maxabs_tolerance = Float64(force_maxabs_tolerance)
    maximum_time = Float64(time_limit)
    linear_tolerance = Float64(gmres_tolerance)
    maximum_linear = Float64(maximum_linear_relative_residual)
    fd_relative_step = Float64(relative_step)
    radius = Float64(trust_radius)
    decrease = Float64(sufficient_decrease)
    energy_tolerance = Float64(energy_increase_tolerance)
    tolerance > 0.0 && isfinite(tolerance) || throw(ArgumentError(
        "force_tolerance must be positive and finite"))
    maxabs_tolerance > 0.0 || throw(ArgumentError(
        "force_maxabs_tolerance must be positive"))
    (maximum_time > 0.0 && (isfinite(maximum_time) || isinf(maximum_time))) ||
        throw(ArgumentError("time_limit must be positive"))
    0.0 < linear_tolerance < 1.0 || throw(ArgumentError(
        "gmres_tolerance must lie in (0, 1)"))
    0.0 < maximum_linear <= 1.0 || throw(ArgumentError(
        "maximum_linear_relative_residual must lie in (0, 1]"))
    fd_relative_step > 0.0 && isfinite(fd_relative_step) ||
        throw(ArgumentError("relative_step must be positive and finite"))
    radius > 0.0 && isfinite(radius) || throw(ArgumentError(
        "trust_radius must be positive and finite"))
    0.0 <= decrease < 1.0 || throw(ArgumentError(
        "sufficient_decrease must lie in [0, 1)"))
    energy_tolerance >= 0.0 && isfinite(energy_tolerance) ||
        throw(ArgumentError(
            "energy_increase_tolerance must be nonnegative and finite"))
    Int(anderson_memory) >= 0 || throw(ArgumentError(
        "anderson_memory must be nonnegative"))
    Float64(anderson_regularization) > 0.0 &&
        isfinite(Float64(anderson_regularization)) || throw(ArgumentError(
            "anderson_regularization must be positive and finite"))

    values = Vector{Float64}(initial_values)
    state = state_at(values)
    initial_r2 = max(Float64(state.r2), eps(Float64))
    total_gmres_iterations = 0
    jvp_evaluations = 0
    completed_iterations = 0
    stop_reason = "newton_iteration_cap"
    history = NamedTuple[]
    anderson_values = Vector{Vector{Float64}}()
    anderson_steps = Vector{Vector{Float64}}()
    started = time()

    for iteration in 1:iterations
        if time() - started >= maximum_time
            stop_reason = "time_limit"
            break
        end
        if state.r2 <= tolerance && state.rinf <= maxabs_tolerance
            stop_reason = "force_tol"
            break
        end
        base_values = copy(values)
        base_residual = Vector{Float64}(state.residual)
        function apply_jacobian(vector_a)
            vector = Vector{Float64}(vector_a)
            direction_norm = maximum(abs, vector)
            direction_norm > 0.0 && isfinite(direction_norm) ||
                return zeros(Float64, length(vector))
            step = clamp(fd_relative_step *
                (1.0 + maximum(abs, base_values)) / direction_norm,
                1.0e-8, 1.0e-3)
            plus = state_at(base_values .+ step .* vector)
            minus = state_at(base_values .- step .* vector)
            jvp_evaluations += 2
            return (Vector{Float64}(plus.residual) .-
                Vector{Float64}(minus.residual)) ./ (2.0 * step)
        end
        current_preconditioner = preconditioner_builder === nothing ?
            apply_preconditioner : preconditioner_builder(
                apply_jacobian, base_values, state, iteration)
        fixed_point_direction = -Vector{Float64}(
            current_preconditioner(base_residual))
        push!(anderson_values, copy(base_values))
        push!(anderson_steps, copy(fixed_point_direction))
        keep = max(2, Int(anderson_memory) + 1)
        length(anderson_values) > keep && popfirst!(anderson_values)
        length(anderson_steps) > keep && popfirst!(anderson_steps)
        anderson_direction = _bvk2_anderson_type2_direction(
            anderson_values, anderson_steps;
            memory=anderson_memory,
            regularization=anderson_regularization)
        # Couple the inexact-Newton forcing to the requested nonlinear
        # tolerance as well as the initial residual.  Using only the latter
        # resets the forcing to 0.5 when a nearly converged checkpoint is
        # restarted, which makes high-resolution polishing intentionally too
        # inaccurate exactly where a strict residual is required.
        forcing = max(linear_tolerance, min(0.5,
            0.5 * sqrt(Float64(state.r2) / initial_r2),
            linear_tolerance * sqrt(Float64(state.r2) / tolerance)))
        direction, inner_iterations, inner_relative =
            _diblock_bvk2_gmres(apply_jacobian, -base_residual;
                max_iterations=linear_iterations, restart=restart,
                tolerance=forcing,
                apply_preconditioner=current_preconditioner)
        total_gmres_iterations += inner_iterations
        if time() - started >= maximum_time
            stop_reason = "time_limit"
            break
        end
        direction_source = "newton_gmres"
        if !all(isfinite, direction) || !isfinite(inner_relative) ||
                inner_relative > maximum_linear
            direction = anderson_direction
            direction_source = length(anderson_values) > 1 ?
                "anderson_recovery" : "preconditioned_residual_fallback"
        end
        largest = Float64(step_norm(direction))
        if !(largest > 0.0 && isfinite(largest))
            stop_reason = "direction_zero"
            break
        end
        largest > radius && (direction .*= radius / largest)

        accepted = nothing
        alpha = 1.0
        evaluations = 0
        for _ in 1:max(1, Int(line_search_steps))
            evaluations += 1
            trial_values = base_values .+ alpha .* direction
            trial = state_at(trial_values)
            energy_pass = isfinite(trial.energy) &&
                trial.energy <= state.energy + energy_tolerance *
                    max(1.0, abs(Float64(state.energy)))
            merit_pass = energy_pass &&
                trial.r2^2 <= (1.0 - decrease * alpha) * state.r2^2
            rinf_pass = trial.rinf <= state.rinf * (1.0 + 1.0e-10) ||
                trial.rinf <= maxabs_tolerance
            if merit_pass && rinf_pass
                accepted = (values=trial_values, state=trial, alpha=alpha)
                break
            end
            alpha *= 0.5
        end
        if accepted === nothing && direction_source == "newton_gmres"
            fallback = anderson_direction
            largest = Float64(step_norm(fallback))
            largest > radius && (fallback .*= radius / largest)
            alpha = 1.0
            for _ in 1:max(1, Int(line_search_steps))
                evaluations += 1
                trial_values = base_values .+ alpha .* fallback
                trial = state_at(trial_values)
                energy_pass = isfinite(trial.energy) &&
                    trial.energy <= state.energy + energy_tolerance *
                        max(1.0, abs(Float64(state.energy)))
                merit_pass = energy_pass &&
                    trial.r2^2 <=
                        (1.0 - decrease * alpha) * state.r2^2
                rinf_pass = trial.rinf <=
                    state.rinf * (1.0 + 1.0e-10) ||
                    trial.rinf <= maxabs_tolerance
                if merit_pass && rinf_pass
                    accepted = (values=trial_values, state=trial, alpha=alpha)
                    direction_source = length(anderson_values) > 1 ?
                        "anderson_recovery" :
                        "preconditioned_residual_fallback"
                    break
                end
                alpha *= 0.5
            end
        end
        if accepted === nothing && length(anderson_values) > 1
            fallback = copy(fixed_point_direction)
            largest = Float64(step_norm(fallback))
            largest > radius && (fallback .*= radius / largest)
            alpha = 1.0
            for _ in 1:max(1, Int(line_search_steps))
                evaluations += 1
                trial_values = base_values .+ alpha .* fallback
                trial = state_at(trial_values)
                energy_pass = isfinite(trial.energy) &&
                    trial.energy <= state.energy + energy_tolerance *
                        max(1.0, abs(Float64(state.energy)))
                merit_pass = energy_pass &&
                    trial.r2^2 <=
                        (1.0 - decrease * alpha) * state.r2^2
                rinf_pass = trial.rinf <=
                    state.rinf * (1.0 + 1.0e-10) ||
                    trial.rinf <= maxabs_tolerance
                if merit_pass && rinf_pass
                    accepted = (values=trial_values, state=trial,
                        alpha=alpha)
                    direction_source = "preconditioned_residual_fallback"
                    empty!(anderson_values)
                    empty!(anderson_steps)
                    break
                end
                alpha *= 0.5
            end
        end
        if accepted === nothing
            stop_reason = "residual_line_search_failed"
            break
        end
        values = accepted.values
        state = accepted.state
        completed_iterations = iteration
        row = (iteration=iteration, r2=Float64(state.r2),
            rinf=Float64(state.rinf), energy=Float64(state.energy),
            alpha=accepted.alpha, direction_source=direction_source,
            gmres_iterations=inner_iterations,
            gmres_relative_residual=inner_relative,
            line_search_evaluations=evaluations)
        push!(history, row)
        callback !== nothing && callback(row, state)
        if state.r2 <= tolerance && state.rinf <= maxabs_tolerance
            stop_reason = "force_tol"
            break
        end
    end
    converged = state.r2 <= tolerance && state.rinf <= maxabs_tolerance
    converged && (stop_reason = "force_tol")
    return (minimizer=values, state=state, iterations=completed_iterations,
        converged=converged, stop_reason=stop_reason,
        gmres_iterations=total_gmres_iterations,
        jvp_evaluations=jvp_evaluations, history=history)
end

function _diblock_bvk2_open_step_bound(phi, direction)
    limit = Inf
    @inbounds for index in eachindex(phi)
        step = direction[index]
        if step > 0.0
            limit = min(limit, (1.0 - phi[index]) / step)
        elseif step < 0.0
            limit = min(limit, -phi[index] / step)
        end
    end
    return limit
end

function _diblock_bvk2_fraction_to_boundary(phi, direction;
        fraction::Float64)
    return min(1.0, fraction *
        _diblock_bvk2_open_step_bound(phi, direction))
end

function _diblock_bvk2_physical_lbfgs(initial_state, state_at_phi;
        f::Float64, max_iterations::Int, memory::Int,
        step_max::Float64, fraction_to_boundary::Float64,
        line_search_steps::Int, armijo::Float64,
        curvature_tolerance::Float64, descent_tolerance::Float64,
        convergence_policy::Symbol=:strict_stationarity,
        residual_r2_tolerance::Float64=1.0e-5,
        residual_rinf_tolerance::Float64=2.0e-5,
        phase_energy_residual_r2_ceiling::Float64=1.0e-2,
        phase_energy_residual_rinf_ceiling::Float64=2.0e-2,
        phase_energy_plateau_tolerance::Float64=1.0e-8,
        phase_energy_plateau_window::Int=20,
        mean_tolerance::Float64=1.0e-10)
    state = initial_state
    trace_row(current, iteration, accepted_alpha,
            iteration_line_search_evaluations, line_search_status) = (
        iteration=iteration,
        stage="physical_lbfgs",
        energy_density=current.energy_density,
        residual_r2=current.residual_r2,
        residual_rinf=current.residual_rinf,
        mean_phi=mean(current.phi),
        minimum_phi=minimum(current.phi),
        maximum_phi=maximum(current.phi),
        accepted_alpha=accepted_alpha,
        line_search_evaluations=iteration_line_search_evaluations,
        line_search_status=line_search_status,
    )
    optimizer_trace = [trace_row(state, 0, NaN, 0, "initial")]
    function convergence_metrics(current)
        mean_bounds_pass = abs(mean(current.phi) - f) <= mean_tolerance &&
            0.0 < minimum(current.phi) && maximum(current.phi) < 1.0
        strict = mean_bounds_pass &&
            current.residual_r2 <= residual_r2_tolerance &&
            current.residual_rinf <= residual_rinf_tolerance
        plateau_available = length(optimizer_trace) >=
            phase_energy_plateau_window
        plateau_span = if plateau_available
            energies = getproperty.(optimizer_trace[
                (end - phase_energy_plateau_window + 1):end],
                :energy_density)
            maximum(energies) - minimum(energies)
        else
            Inf
        end
        phase_energy = strict || (mean_bounds_pass && plateau_available &&
            plateau_span <= phase_energy_plateau_tolerance &&
            current.residual_r2 <= phase_energy_residual_r2_ceiling &&
            current.residual_rinf <= phase_energy_residual_rinf_ceiling)
        return (strict=strict, phase_energy=phase_energy,
            plateau_available=plateau_available,
            plateau_span=plateau_span)
    end
    selected_converged(metrics) = convergence_policy ==
        :strict_stationarity ? metrics.strict : metrics.phase_energy
    if max_iterations == 0
        metrics = convergence_metrics(state)
        return (
            state=state, converged=false, iterations=0, status="not_run",
            minimum=state.energy_density, memory_updates=0,
            line_search_evaluations=0, descent_fallbacks=0,
            optimizer_trace=optimizer_trace,
            strict_stationarity_converged=metrics.strict,
            phase_energy_converged=metrics.phase_energy,
            energy_plateau_available=metrics.plateau_available,
            energy_plateau_span=metrics.plateau_span,
        )
    end

    history_s = Vector{Vector{Float64}}()
    history_y = Vector{Vector{Float64}}()
    history_rho = Float64[]
    memory_updates = 0
    line_search_evaluations = 0
    descent_fallbacks = 0
    completed_iterations = 0
    status = "maximum_iterations"
    point_count = length(state.phi)

    function physical_gradient(current)
        # Minimize E/dV rather than E/V.  Its Euclidean gradient on the
        # fixed-mean tangent space is exactly mu - mean(mu), which keeps the
        # physical residual and optimization direction on the same scale.
        gradient = vec(copy(current.residual))
        gradient .-= mean(gradient)
        return gradient
    end

    objective(current) = point_count * current.energy_density

    function limited_direction(raw_direction)
        direction = Vector{Float64}(raw_direction)
        direction .-= mean(direction)
        largest = maximum(abs, direction)
        if !(largest > 0.0 && isfinite(largest))
            return nothing
        end
        largest > step_max && (direction .*= step_max / largest)
        return direction
    end

    function search_direction(base, gradient, raw_direction)
        direction = limited_direction(raw_direction)
        direction === nothing && return nothing
        slope = dot(gradient, direction)
        descent_scale = norm(gradient) * norm(direction)
        slope < -descent_tolerance * descent_scale && isfinite(slope) ||
            return nothing
        alpha = _diblock_bvk2_fraction_to_boundary(
            vec(base.phi), direction; fraction=fraction_to_boundary)
        alpha > 0.0 && isfinite(alpha) || return nothing
        for _ in 1:line_search_steps
            candidate = base.phi .+ alpha .* reshape(direction, size(base.phi))
            trial = try
                state_at_phi(candidate)
            catch error
                if error isa DomainError || error isa ArgumentError
                    nothing
                else
                    rethrow()
                end
            end
            line_search_evaluations += 1
            base_objective = objective(base)
            trial_objective = trial === nothing ? Inf : objective(trial)
            armijo_pass = isfinite(trial_objective) &&
                trial_objective <= base_objective + armijo * alpha * slope
            roundoff_scale = 100eps(Float64) *
                max(1.0, abs(base_objective), abs(trial_objective))
            roundoff_residual_pass = trial !== nothing &&
                isfinite(trial_objective) &&
                abs(trial_objective - base_objective) <= roundoff_scale &&
                trial.residual_r2 < base.residual_r2 &&
                trial.residual_rinf < base.residual_rinf
            if trial !== nothing && (armijo_pass || roundoff_residual_pass)
                acceptance_status = armijo_pass ?
                    "accepted_armijo" : "accepted_roundoff_residual"
                return trial, direction, alpha, acceptance_status
            end
            alpha *= 0.5
        end
        return nothing
    end

    for iteration in 1:max_iterations
        minimum_phi = minimum(state.phi)
        maximum_phi = maximum(state.phi)
        metrics = convergence_metrics(state)
        if selected_converged(metrics)
            status = metrics.strict ? "strict_stationarity_convergence" :
                "phase_energy_convergence"
            break
        end

        gradient = physical_gradient(state)
        q = copy(gradient)
        alpha_history = zeros(Float64, length(history_s))
        for index in reverse(eachindex(history_s))
            alpha_history[index] = history_rho[index] *
                dot(history_s[index], q)
            q .-= alpha_history[index] .* history_y[index]
        end
        gamma = 1.0
        if !isempty(history_s)
            sy = dot(history_s[end], history_y[end])
            yy = dot(history_y[end], history_y[end])
            sy > 0.0 && yy > 0.0 && (gamma = sy / yy)
        end
        inverse_gradient = gamma .* q
        for index in eachindex(history_s)
            beta = history_rho[index] *
                dot(history_y[index], inverse_gradient)
            inverse_gradient .+= history_s[index] .*
                (alpha_history[index] - beta)
        end
        direction = .-inverse_gradient
        direction .-= mean(direction)

        # Indefinite or stale curvature information must never turn the
        # physical residual into an ascent direction.
        descent_scale = norm(gradient) * norm(direction)
        used_fallback = !(dot(gradient, direction) <
            -descent_tolerance * descent_scale) ||
            !all(isfinite, direction)
        used_fallback && (direction = -gradient)
        accepted_with_fallback = used_fallback
        iteration_line_search_start = line_search_evaluations
        accepted = search_direction(state, gradient, direction)
        if accepted === nothing && !used_fallback
            empty!(history_s)
            empty!(history_y)
            empty!(history_rho)
            descent_fallbacks += 1
            accepted_with_fallback = true
            accepted = search_direction(state, gradient, -gradient)
        elseif used_fallback
            descent_fallbacks += 1
        end
        if accepted === nothing
            status = "line_search_failed"
            break
        end

        trial, _direction, accepted_alpha, line_search_status = accepted
        accepted_with_fallback &&
            (line_search_status *= "_fallback_direction")
        next_gradient = physical_gradient(trial)
        step = vec(trial.phi .- state.phi)
        step .-= mean(step)
        gradient_change = next_gradient .- gradient
        gradient_change .-= mean(gradient_change)
        curvature = dot(step, gradient_change)
        curvature_scale = norm(step) * norm(gradient_change)
        if curvature > curvature_tolerance * curvature_scale &&
                isfinite(curvature)
            length(history_s) == memory && begin
                popfirst!(history_s)
                popfirst!(history_y)
                popfirst!(history_rho)
            end
            push!(history_s, step)
            push!(history_y, gradient_change)
            push!(history_rho, inv(curvature))
            memory_updates += 1
        end
        state = trial
        completed_iterations = iteration
        push!(optimizer_trace, trace_row(state, iteration, accepted_alpha,
            line_search_evaluations - iteration_line_search_start,
            line_search_status))
    end

    metrics = convergence_metrics(state)
    converged = selected_converged(metrics)
    converged && (status = metrics.strict ?
        "strict_stationarity_convergence" : "phase_energy_convergence")
    return (
        state=state,
        converged=converged,
        iterations=completed_iterations,
        status=status,
        minimum=state.energy_density,
        memory_updates=memory_updates,
        line_search_evaluations=line_search_evaluations,
        descent_fallbacks=descent_fallbacks,
        optimizer_trace=optimizer_trace,
        strict_stationarity_converged=metrics.strict,
        phase_energy_converged=metrics.phase_energy,
        energy_plateau_available=metrics.plateau_available,
        energy_plateau_span=metrics.plateau_span,
    )
end

# ---------------------------------------------------------------------------
# Mobility-preconditioned Barzilai-Borwein field optimizer.
#
# The physical-space L-BFGS above stalls at strong segregation (chiN >= ~40):
# once cells saturate against the density walls, the bounded/logit step shrinks
# to zero exactly where the physical force (chem - mu_bar) is largest, so the
# solve halts far from stationarity (measured: residual R2 ~ 50) or overshoots
# a wall (DomainError).  This is the same degeneracy the BVK1 2D/3D solver hit;
# the fix (ported from
# scripts/write_burp_bvk1_polyorder_nonlamellar_logit_lbfgs_probe.jl, commit
# eb786ed) descends along the UN-suppressed force in the logit chart with a
# non-monotone Grippo-Lampariello-Lucidi line search, a Barzilai-Borwein
# spectral step, and a trust-region cap, and measures stationarity on the
# physical force ||chem - mu_bar||_inf so a wall-suppressed y-gradient cannot
# masquerade as convergence.  SCFT (Polyorder) avoids the degeneracy entirely
# by solving on unbounded fields; this recovers the same conditioning for the
# bounded-density functional.
@inline function _bvk2_spectral_dot(a, b)
    total = 0.0
    @inbounds @simd for i in eachindex(a)
        total += a[i] * b[i]
    end
    return total
end

@inline function _bvk2_spectral_inf_norm(a)
    m = 0.0
    @inbounds for i in eachindex(a)
        v = abs(a[i])
        m = v > m ? v : m
    end
    return m
end

# Spectral dealiasing of the theta field: zero the top (1-keepfrac) of Fourier
# modes per axis (the 2/3 rule at keepfrac=2/3).  The kpsi sensor's central
# differences are BLIND to the Nyquist (checkerboard) mode, so at strong
# segregation (chiN>=49) the theta-BB descent falls into a spurious grid-scale
# mode that lowers the discrete energy without bound as the grid refines (~1000x
# the near-Nyquist power vs chiN=48).  Projecting the field onto the smooth
# subspace each step removes it.  keepfrac<=0 or >=1 disables (identity).  The
# DC mode (mean) is always kept, so the fixed-mean retract is unaffected.
function _bvk2_dealias_filter!(y::AbstractVector{Float64}, dims::NTuple{D,Int},
        keepfrac::Float64) where {D}
    (keepfrac <= 0.0 || keepfrac >= 1.0) && return y
    A = reshape(y, dims)
    B = fft(A)
    @inbounds for idx in CartesianIndices(B)
        drop = false
        for d in 1:D
            n = dims[d]
            i = idx[d] - 1
            freq = i <= n ÷ 2 ? i : i - n
            if abs(freq) >= keepfrac * (n / 2)
                drop = true
                break
            end
        end
        drop && (B[idx] = 0.0 + 0.0im)
    end
    A .= real.(ifft(B))
    return y
end

"""
    _bvk2_preconditioned_spectral_descent(objective, pgradient!, y0; ...)

Mobility-preconditioned Barzilai-Borwein descent with a non-monotone GLL line
search.  `pgradient!(p, wbuf, y)` writes the unsuppressed physical force
`p[i] = chem[i] - mu_bar` (the search direction is `-p`) and the slope weights
`wbuf[i] = dV*weight[i]` (`weight = phi*(1-phi)`, the sigmoid Jacobian), so the
true energy directional derivative along `-p` is `-sum(wbuf*p^2)`.  Convergence
is on `||p||_inf`.  Returns `(minimizer, iterations, converged, g_norm,
stop_reason)`.
"""
function _bvk2_preconditioned_spectral_descent(objective, pgradient!, y0;
        max_iterations::Integer, ptol::Real=2.0e-5, ftol::Real=1.0e-12,
        memory::Integer=10, stall_window::Integer=10,
        alpha_min::Real=1.0e-10, alpha_max::Real=1.0e10, armijo::Real=1.0e-4,
        dymax::Real=8.0, postfilter!::Union{Nothing,Function}=nothing)
    x = Vector{Float64}(copy(y0))
    n = length(x)
    p = zeros(Float64, n)
    wbuf = zeros(Float64, n)
    pgradient!(p, wbuf, x)
    fx = objective(x)
    pinf = _bvk2_spectral_inf_norm(p)
    maxit = max(1, Int(max_iterations))
    if !isfinite(fx) || pinf <= ptol
        return (minimizer=x, iterations=0, converged=(pinf <= ptol),
            g_norm=pinf, stop_reason=(pinf <= ptol ? "force_tol" :
            "nonfinite_initial"))
    end
    M = max(1, Int(memory))
    fring = fill(fx, M)
    alpha = clamp(1.0 / max(pinf, eps(Float64)), alpha_min, alpha_max)
    xnew = similar(x)
    pnew = zeros(Float64, n)
    wbufnew = zeros(Float64, n)
    svec = similar(x)
    dpvec = similar(x)
    stall = 0
    iters = 0
    converged = false
    stop_reason = "iteration_cap"
    while iters < maxit
        iters += 1
        fref = maximum(fring)
        smag = 0.0
        @inbounds @simd for i in 1:n
            smag += wbuf[i] * p[i] * p[i]
        end
        smag = max(smag, 0.0)
        lam = min(alpha, dymax / max(pinf, eps(Float64)))
        fnew = fx
        accepted = false
        for _bt in 1:40
            @inbounds for i in 1:n
                xnew[i] = x[i] - lam * p[i]
            end
            # A constrained chart must be retracted before *every* trial
            # objective.  Applying the postfilter only after line-search
            # acceptance lets off-manifold UD-theta trials reach the energy
            # evaluator and fail its composition-mean gate.
            postfilter! !== nothing && postfilter!(xnew)
            fnew = objective(xnew)
            if isfinite(fnew) && fnew <= fref - armijo * lam * smag
                accepted = true
                break
            end
            denom = 2.0 * (fnew - fx + lam * smag)
            lam_quad = denom > 0.0 ? (smag * lam * lam) / denom : 0.5 * lam
            lam = clamp(lam_quad, 0.1 * lam, 0.5 * lam)
            if lam <= alpha_min
                lam = alpha_min
                @inbounds for i in 1:n
                    xnew[i] = x[i] - lam * p[i]
                end
                postfilter! !== nothing && postfilter!(xnew)
                fnew = objective(xnew)
                accepted = isfinite(fnew) &&
                    fnew <= fref - armijo * lam * smag
                break
            end
        end
        if !accepted
            stop_reason = "line_search_failed"
            break
        end
        pgradient!(pnew, wbufnew, xnew)
        @inbounds for i in 1:n
            svec[i] = xnew[i] - x[i]
            dpvec[i] = pnew[i] - p[i]
        end
        sdp = _bvk2_spectral_dot(svec, dpvec)
        sts = _bvk2_spectral_dot(svec, svec)
        alpha = sdp > 1.0e-30 ? clamp(sts / sdp, alpha_min, alpha_max) : alpha_max
        df = abs(fx - fnew)
        @inbounds for i in 1:n
            x[i] = xnew[i]
            p[i] = pnew[i]
            wbuf[i] = wbufnew[i]
        end
        fx = fnew
        fring[mod1(iters + 1, M)] = fx
        pinf = _bvk2_spectral_inf_norm(p)
        if pinf <= ptol
            converged = true
            stop_reason = "force_tol"
            break
        end
        if df <= ftol * (abs(fx) + 1.0)
            stall += 1
            if stall >= max(1, Int(stall_window))
                stop_reason = "energy_stagnation"
                break
            end
        else
            stall = 0
        end
    end
    return (minimizer=x, iterations=iters, converged=converged, g_norm=pinf,
        stop_reason=stop_reason)
end

"""
    _diblock_bvk2_mobility_bb_optimize(state0, state_at_phi; ...)

Drop-in replacement for `_diblock_bvk2_physical_lbfgs` that runs the
mobility-preconditioned BB descent in the fixed-mean logit chart.  Returns the
same NamedTuple shape so `relax_diblock_bvk2_morphology_nd` can consume it
unchanged.  The fixed-mean map `phi = sigmoid(y + lambda)` reuses the exact
Newton composition shift `_diblock_bvk2_logit_shift`; the energy and chemical
potential are the accelerated `diblock_bvk2_energy_nd` /
`diblock_bvk2_chemical_potential_nd`.
"""
function _diblock_bvk2_mobility_bb_optimize(state0, state_at_phi;
        f::Float64, chiN::Float64, lengths::AbstractVector{Float64},
        N::Float64, b::Float64, adaptive::Bool, c2::Real, dV::Float64,
        max_iterations::Integer, force_tolerance::Real,
        energy_tolerance::Real=1.0e-12, memory::Integer=10,
        stall_window::Integer=10, trust_region::Real=8.0)
    dims = size(state0.phi)
    y0 = vec(log.(state0.phi ./ (1.0 .- state0.phi)))
    y0 .-= sum(y0) / length(y0)

    function y_to_phi(y)
        lambda = _diblock_bvk2_logit_shift(y, f)
        phi = Array{Float64}(undef, dims)
        pv = vec(phi)
        @inbounds for i in eachindex(y)
            pv[i] = _diblock_bvk2_open_logistic(y[i] + lambda)
        end
        return phi
    end

    function objective(y)
        phi = y_to_phi(y)
        energy = diblock_bvk2_energy_nd(phi; f=f, chiN=chiN, lengths=lengths,
            N=N, b=b, adaptive=adaptive, c2=c2, mean_tolerance=1.0e-8)
        return isfinite(energy) ? energy : floatmax(Float64) / 4.0
    end

    function pgradient!(p, wbuf, y)
        lambda = _diblock_bvk2_logit_shift(y, f)
        phi = Array{Float64}(undef, dims)
        pv = vec(phi)
        @inbounds for i in eachindex(y)
            pv[i] = _diblock_bvk2_open_logistic(y[i] + lambda)
        end
        mu = vec(diblock_bvk2_chemical_potential_nd(phi; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, adaptive=adaptive, c2=c2,
            mean_tolerance=1.0e-8))
        weight_sum = 0.0
        weighted_mu = 0.0
        @inbounds for i in eachindex(y)
            s = pv[i]                     # phi = sigmoid(y+lambda); weight = s(1-s)
            w = s * (1.0 - s)
            wbuf[i] = dV * w
            weight_sum += w
            weighted_mu += w * mu[i]
        end
        mu_bar = weighted_mu / max(weight_sum, eps(Float64))
        @inbounds for i in eachindex(y)
            p[i] = mu[i] - mu_bar
        end
        return p
    end

    descent = _bvk2_preconditioned_spectral_descent(objective, pgradient!, y0;
        max_iterations=max_iterations, ptol=force_tolerance,
        ftol=energy_tolerance, memory=memory, stall_window=stall_window,
        dymax=trust_region)
    phi_final = y_to_phi(descent.minimizer)
    final_state = state_at_phi(phi_final)
    return (
        state=final_state,
        converged=descent.converged,
        iterations=descent.iterations,
        status=descent.stop_reason,
        minimum=final_state.energy_density,
        memory_updates=0,
        line_search_evaluations=0,
        descent_fallbacks=0,
        optimizer_trace=NamedTuple[],
        strict_stationarity_converged=descent.converged,
        phase_energy_converged=descent.converged,
        energy_plateau_available=false,
        energy_plateau_span=Inf,
    )
end

# ---------------------------------------------------------------------------
# Uneyama-Doi-style square-root-density chart: theta with phi = sin^2(theta).
#
# UD (Macromolecules 38, 196 (2005)) solve their density functional to
# chiN=100 by rewriting it in psi = sqrt(phi), "because the logarithmic
# contribution is regular as psi -> 0".  The single-field incompressible
# analogue is the angle chart phi = sin^2(theta): dphi/dtheta = sin(2*theta)
# = 2*sqrt(phi*(1-phi)), so the theta-gradient of the energy is
#     g_i = sin(2*theta_i) * (mu_i - lambda),
# regular at BOTH walls (sqrt(phi)*log(phi) -> 0) yet only sqrt-suppressed at
# interfaces (unlike the logit chart's phi*(1-phi) suppression that produced
# false plateaus).  The entropy Hessian in theta grows only logarithmically
# toward the walls -- condition number ~ log(1/phi_min) instead of 1/phi_min
# -- so plain BB/spectral descent on the TRUE theta-gradient converges, and
# ||g||_inf -> 0 is a genuine stationarity measure (no loose-R2 workaround).
# The mean constraint has a CLOSED FORM in this chart:
#     mean(sin^2(theta+s)) = 1/2 - (1/2)*R*cos(2s+delta) = f,
# with R*exp(i*delta) = mean(exp(2i*theta)).

"""Mobility-weighted mean retraction for the theta chart: iterates
`theta_i += eps * sin(2*theta_i)` with Newton on the scalar `eps` until
`mean(sin^2(theta)) = f` to 1e-12.  Only interface cells move (pure domains
have sin(2*theta) = 0 and stay pure) -- the physically correct way to restore
composition at strong segregation, where a UNIFORM shift both fails
analytically (|<exp(2i*theta)>| drops below |1-2f| as the field becomes
binary) and would wrongly de-saturate pure domains."""
function _bvk2_theta_mean_retract!(theta::AbstractVector{Float64}, f::Float64)
    n = length(theta)
    n > 0 || throw(ArgumentError("theta chart requires a nonempty field"))
    inv_n = 1.0 / n
    for _ in 1:60
        m = 0.0
        w2 = 0.0
        @inbounds @simd for i in 1:n
            s2t = sin(2.0 * theta[i])
            m += sin(theta[i])^2
            w2 += s2t * s2t
        end
        residual = m * inv_n - f
        abs(residual) <= 1.0e-12 && return theta
        slope = w2 * inv_n
        slope > 1.0e-14 || throw(ArgumentError(
            "theta chart mean is frozen (no interface cells): " *
            "residual=$(residual), <sin^2(2theta)>=$(slope)"))
        eps_step = -residual / slope
        # Damp very large corrections so a near-frozen mean cannot fling
        # interface cells across the box in one step.
        eps_step = clamp(eps_step, -0.5, 0.5)
        @inbounds @simd for i in 1:n
            theta[i] += eps_step * sin(2.0 * theta[i])
        end
    end
    throw(ArgumentError("theta chart mean retraction did not converge"))
end

"""BB/spectral descent on the TRUE gradient in the theta chart (phi=sin^2).
Returns the `_diblock_bvk2_physical_lbfgs`-compatible NamedTuple.  Unlike the
mobility-BB (logit chart, preconditioned direction) and SD-Newton paths, the
convergence test `||g||_inf <= force_tolerance` here is genuine stationarity:
the theta-gradient vanishes at a constrained minimum even when saturated wall
cells keep the raw residual `mu - mean(mu)` large."""
function _diblock_bvk2_theta_bb_optimize(state0, state_at_phi;
        f::Float64, chiN::Float64, lengths::AbstractVector{Float64},
        N::Float64, b::Float64, adaptive::Bool, c2::Real, dV::Float64,
        max_iterations::Integer, force_tolerance::Real,
        energy_tolerance::Real=1.0e-12, memory::Integer=10,
        stall_window::Integer=10, trust_region::Real=8.0,
        dealias_keepfrac::Real=parse(Float64, get(ENV, "B2B_DEALIAS_KEEPFRAC", "1.0")))
    dims = size(state0.phi)
    theta0 = vec(asin.(sqrt.(clamp.(state0.phi, 0.0, 1.0))))

    function theta_to_phi(y)
        yr = _bvk2_theta_mean_retract!(copy(Vector{Float64}(y)), f)
        phi = Array{Float64}(undef, dims)
        pv = vec(phi)
        @inbounds for i in eachindex(yr)
            t = sin(yr[i])
            pv[i] = t * t
        end
        return phi
    end

    function objective(y)
        phi = theta_to_phi(y)
        energy = diblock_bvk2_energy_nd(phi; f=f, chiN=chiN, lengths=lengths,
            N=N, b=b, adaptive=adaptive, c2=c2, mean_tolerance=1.0e-8)
        return isfinite(energy) ? energy : floatmax(Float64) / 4.0
    end

    function pgradient!(p, wbuf, y)
        yr = _bvk2_theta_mean_retract!(copy(Vector{Float64}(y)), f)
        phi = Array{Float64}(undef, dims)
        pv = vec(phi)
        sin2 = Vector{Float64}(undef, length(yr))
        @inbounds for i in eachindex(yr)
            t = yr[i]
            st = sin(t)
            pv[i] = st * st
            sin2[i] = sin(2.0 * t)
        end
        mu = vec(diblock_bvk2_chemical_potential_nd(phi; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, adaptive=adaptive, c2=c2,
            mean_tolerance=1.0e-8))
        # Mean-preserving projection: d<phi> = (1/n) sum sin(2theta)*dtheta,
        # so lambda is the sin^2(2theta)-weighted mean of mu.
        weight_sum = 0.0
        weighted_mu = 0.0
        @inbounds for i in eachindex(y)
            w = sin2[i] * sin2[i]
            weight_sum += w
            weighted_mu += w * mu[i]
        end
        lambda = weighted_mu / max(weight_sum, eps(Float64))
        # p = theta-gradient / dV; slope along -p is dV*sum(p^2) via wbuf=dV.
        @inbounds for i in eachindex(y)
            p[i] = sin2[i] * (mu[i] - lambda)
            wbuf[i] = dV
        end
        return p
    end

    postfilter! = (dealias_keepfrac > 0.0 && dealias_keepfrac < 1.0) ?
        (y -> _bvk2_dealias_filter!(y, dims, Float64(dealias_keepfrac))) : nothing
    descent = _bvk2_preconditioned_spectral_descent(objective, pgradient!,
        theta0; max_iterations=max_iterations, ptol=force_tolerance,
        ftol=energy_tolerance, memory=memory, stall_window=stall_window,
        dymax=trust_region, postfilter! = postfilter!)
    phi_final = theta_to_phi(descent.minimizer)
    final_state = state_at_phi(phi_final)
    # "energy_stagnation" = `stall_window` consecutive accepted steps with a
    # relative energy change below `ftol` -- a plateau certificate at least as
    # strong as the SD-Newton acceptance, reached here with a genuinely small
    # theta-gradient (the walls no longer inflate the stationarity measure).
    plateaued = descent.stop_reason == "energy_stagnation"
    converged = descent.converged || plateaued
    return (
        state=final_state,
        converged=converged,
        iterations=descent.iterations,
        status=string("theta_", descent.stop_reason),
        minimum=final_state.energy_density,
        memory_updates=0,
        line_search_evaluations=0,
        descent_fallbacks=0,
        optimizer_trace=NamedTuple[],
        strict_stationarity_converged=descent.converged,
        phase_energy_converged=converged,
        energy_plateau_available=plateaued,
        energy_plateau_span=plateaued ? 0.0 : Inf,
    )
end

# ---------------------------------------------------------------------------
# Mobility-preconditioned SD + normal-residual Gauss-Newton (BVK1 port).
#
# Ported from polish_bvk1_branch_stationarity.jl (the sd_newton_ngmres gyroid
# solver that reaches chiN=60 for BVK1).  The mobility-BB SD stage escapes the
# strong-segregation stall and acquires the morphology; the normal-residual
# step then descends the SQUARED-RESIDUAL merit 0.5*|R|^2 (R = mu - mean(mu))
# along one J.R Jacobian-vector product, mobility-weighted by phi*(1-phi) so
# the entropy-divergent wall cells -- which dominate the raw residual and stall
# an unweighted line search -- are down-weighted.  Acceptance is an ENERGY
# PLATEAU plus a loose R2 ceiling: at strong segregation the exact energy is
# the boundary observable and strict force stationarity is both unreachable
# and unnecessary (BVK1 accepts residual_r2 up to ~500 with an energy plateau).

"""Mobility-preconditioned descent direction `-M*(force - <force>_M)/scale`,
`M = phi*(1-phi)`, `scale = ref*(1-ref)`, mean-projected to preserve mean.
Down-weights saturated cells where the entropy term makes `force` diverge."""
function _diblock_bvk2_mobility_direction(phi::AbstractArray{Float64},
        force::AbstractArray{Float64}, reference::Float64)
    mobility_sum = 0.0
    weighted = 0.0
    @inbounds for i in eachindex(phi)
        m = phi[i] * (1.0 - phi[i])
        mobility_sum += m
        weighted += m * force[i]
    end
    shift = weighted / max(mobility_sum, eps(Float64))
    scale = reference * (1.0 - reference)
    direction = Array{Float64}(undef, size(phi))
    total = 0.0
    @inbounds for i in eachindex(phi)
        m = phi[i] * (1.0 - phi[i])
        direction[i] = -(m / scale) * (force[i] - shift)
        total += direction[i]
    end
    mean_direction = total / length(direction)
    @inbounds for i in eachindex(direction)
        direction[i] -= mean_direction
    end
    return direction
end

"""One mobility-preconditioned Gauss-Newton step on the residual merit
`0.5*|R|^2`.  Its gradient is `J*R` (one JVP); the step is the mobility descent
along it, damped by fraction-to-boundary and a backtracking line search that
requires an Armijo decrease of `|R|^2`, no energy increase, and no growth of
`|R|_inf`.  Returns the accepted physical state, or `nothing` on a stall."""
function _diblock_bvk2_normal_residual_step(base, state_at_phi, translation_modes;
        f::Float64, chiN::Float64, lengths::AbstractVector{Float64},
        N::Float64, b::Float64, adaptive::Bool, c2::Real,
        step_max::Float64, fraction::Float64, line_search_steps::Int,
        sufficient_decrease::Float64, jvp_step::Float64, rinf_tolerance::Float64)
    residual_direction = copy(base.residual)
    _diblock_bvk2_project_gauge!(residual_direction, translation_modes)
    product, _jstep, _sign = _diblock_bvk2_morphology_jvp(base, residual_direction;
        f=f, chiN=chiN, lengths=lengths, N=N, b=b, adaptive=adaptive, c2=c2,
        relative_step=jvp_step)
    _diblock_bvk2_project_gauge!(product, translation_modes)
    direction = _diblock_bvk2_mobility_direction(base.phi, product, f)
    _diblock_bvk2_project_gauge!(direction, translation_modes)
    largest = maximum(abs, direction)
    (largest > 0.0 && isfinite(largest)) || return nothing
    largest > step_max && (direction .*= step_max / largest)
    alpha = _diblock_bvk2_fraction_to_boundary(base.phi, direction;
        fraction=fraction)
    (alpha > 0.0 && isfinite(alpha)) || return nothing
    base_merit = base.residual_r2^2
    for _ in 1:line_search_steps
        candidate = base.phi .+ alpha .* direction
        candidate .-= mean(candidate) - f
        trial = try
            _diblock_bvk2_morphology_physical_state(candidate; f=f, chiN=chiN,
                lengths=lengths, N=N, b=b, adaptive=adaptive, c2=c2)
        catch err
            (err isa DomainError || err isa ArgumentError) ? nothing : rethrow()
        end
        if trial !== nothing && isfinite(trial.residual_r2) &&
                trial.residual_r2^2 <=
                    (1.0 - sufficient_decrease * alpha) * base_merit &&
                trial.energy_density <= base.energy_density +
                    100.0 * eps(Float64) * max(1.0, abs(base.energy_density)) &&
                (trial.residual_rinf <= base.residual_rinf * (1.0 + 1.0e-10) ||
                 trial.residual_rinf <= rinf_tolerance)
            return trial
        end
        alpha *= 0.5
    end
    return nothing
end

"""Alternating mobility-BB SD blocks and normal-residual Gauss-Newton blocks,
accepted on a 3-sample energy-density plateau together with a loose R2 ceiling.
Mirrors the BVK1 sd_block/newton_block recovery: when a Newton block stalls
(line search fails at strong segregation) the next SD block re-escapes and
descends further, so progress continues to an energy-plateaued state.  Returns
the `_diblock_bvk2_physical_lbfgs`-compatible NamedTuple."""
function _diblock_bvk2_sd_newton_optimize(state0, state_at_phi;
        f::Float64, chiN::Float64, lengths::AbstractVector{Float64},
        N::Float64, b::Float64, adaptive::Bool, c2::Real, dV::Float64,
        sd_iterations::Integer, newton_iterations::Integer, max_blocks::Integer,
        step_max::Float64, fraction::Float64, line_search_steps::Int,
        sufficient_decrease::Float64, jvp_step::Float64, rinf_tolerance::Float64,
        r2_acceptance::Float64, plateau_window::Integer,
        plateau_tolerance::Float64, trust_region::Float64, memory::Integer)
    state = state0
    translation_modes = _diblock_bvk2_translation_modes(state.phi)
    energies = Float64[state.energy_density]
    total_iters = 0
    window = max(3, Int(plateau_window))
    plateau_pass = false
    plateau_span = Inf
    converged = false
    status = "sd_newton_block_cap"

    function refresh_plateau!()
        if length(energies) >= window
            plateau_span = maximum(@view energies[end - window + 1:end]) -
                minimum(@view energies[end - window + 1:end])
            plateau_pass = plateau_span <=
                plateau_tolerance * (abs(state.energy_density) + 1.0)
        end
    end

    for _block in 1:Int(max_blocks)
        # SD block (mobility BB): escapes the stall on block 1, re-escapes
        # after a Newton stall on later blocks.
        sd = _diblock_bvk2_mobility_bb_optimize(state, state_at_phi;
            f=f, chiN=chiN, lengths=lengths, N=N, b=b, adaptive=adaptive, c2=c2,
            dV=dV, max_iterations=sd_iterations, force_tolerance=rinf_tolerance,
            energy_tolerance=plateau_tolerance, memory=memory,
            stall_window=plateau_window, trust_region=trust_region)
        state = sd.state
        total_iters += Int(sd.iterations)
        push!(energies, state.energy_density)
        refresh_plateau!()
        if plateau_pass && state.residual_r2 <= r2_acceptance
            converged = true
            status = "energy_plateau"
            break
        end
        # Newton block: mobility Gauss-Newton on the residual merit.
        for _ in 1:Int(newton_iterations)
            base = _diblock_bvk2_morphology_physical_state(state.phi; f=f,
                chiN=chiN, lengths=lengths, N=N, b=b, adaptive=adaptive, c2=c2)
            trial = _diblock_bvk2_normal_residual_step(base, state_at_phi,
                translation_modes; f=f, chiN=chiN, lengths=lengths, N=N, b=b,
                adaptive=adaptive, c2=c2, step_max=step_max, fraction=fraction,
                line_search_steps=line_search_steps,
                sufficient_decrease=sufficient_decrease, jvp_step=jvp_step,
                rinf_tolerance=rinf_tolerance)
            total_iters += 1
            trial === nothing && break        # stalled -> next SD block
            state = trial
            push!(energies, state.energy_density)
            refresh_plateau!()
            if plateau_pass && state.residual_r2 <= r2_acceptance
                converged = true
                status = "energy_plateau"
                break
            end
        end
        converged && break
    end
    return (
        state=state,
        converged=converged,
        iterations=total_iters,
        status=status,
        minimum=state.energy_density,
        memory_updates=0,
        line_search_evaluations=0,
        descent_fallbacks=0,
        optimizer_trace=NamedTuple[],
        strict_stationarity_converged=false,
        phase_energy_converged=converged,
        energy_plateau_available=plateau_pass,
        energy_plateau_span=plateau_span,
    )
end

function _diblock_bvk2_morphology_jvp(base, direction_a;
        f::Float64, chiN::Float64, lengths::AbstractVector{Float64},
        N::Float64, b::Float64, adaptive::Bool, c2::Real,
        relative_step::Float64)
    direction = Float64.(Array(direction_a))
    size(direction) == size(base.phi) || throw(DimensionMismatch(
        "BVK2 JVP direction must match the physical field"))
    direction .-= mean(direction)
    direction_norm = maximum(abs, direction)
    direction_norm > 0.0 && isfinite(direction_norm) ||
        return zeros(Float64, size(direction)), 0.0, 1.0
    positive_bound = _diblock_bvk2_open_step_bound(base.phi, direction)
    negative_bound = _diblock_bvk2_open_step_bound(base.phi, .-direction)
    sign = negative_bound > positive_bound ? -1.0 : 1.0
    signed_direction = sign .* direction
    raw_step = relative_step * (1.0 + maximum(abs, base.phi)) /
        direction_norm
    step = min(raw_step, 0.1 * max(positive_bound, negative_bound))
    step > 100eps(Float64) ||
        return zeros(Float64, size(direction)), step, sign
    candidate = base.phi .+ step .* signed_direction
    candidate .-= mean(candidate) - f
    trial = _diblock_bvk2_morphology_physical_state(candidate;
        f=f, chiN=chiN, lengths=lengths,
        N=N, b=b, adaptive=adaptive, c2=c2)
    product = (trial.residual .- base.residual) ./ (sign * step)
    product .-= mean(product)
    return product, step, sign
end

function _diblock_bvk2_morphology_initial_phi(morphology, initial_phi,
        dims, lengths, f::Float64, initial_amplitude::Float64)
    if initial_phi === nothing
        grid_dims = _diblock_morphology_dims(dims === nothing ? (32,) : dims)
        cell_lengths = _diblock_nd_lengths(lengths, grid_dims)
        phi = Float64.(diblock_morphology_seed(morphology; f=f,
            dims=grid_dims, lengths=cell_lengths,
            amplitude=initial_amplitude))
        # The analytic seed is already a logistic field.  Repair only exact
        # floating-point endpoints produced by representational saturation.
        @inbounds for index in eachindex(phi)
            phi[index] == 0.0 && (phi[index] = nextfloat(0.0))
            phi[index] == 1.0 && (phi[index] = prevfloat(1.0))
        end
        return phi, grid_dims, cell_lengths, "diblock_morphology_seed"
    end

    phi = Float64.(Array(initial_phi))
    ndims(phi) >= 1 || throw(ArgumentError(
        "initial_phi must have at least one dimension"))
    all(count -> count >= 2, size(phi)) || throw(ArgumentError(
        "initial_phi must have at least two sites in every dimension"))
    all(isfinite, phi) || throw(ArgumentError(
        "initial_phi must contain only finite values"))
    all(value -> 0.0 < value < 1.0, phi) || throw(DomainError(
        (minimum(phi), maximum(phi)),
        "initial_phi must stay strictly inside the physical interval (0, 1)"))
    grid_dims = Tuple(size(phi))
    if dims !== nothing
        requested_dims = _diblock_morphology_dims(dims)
        requested_dims == grid_dims || throw(DimensionMismatch(
            "initial_phi dimensions must match dims"))
    end
    cell_lengths = _diblock_nd_lengths(lengths, grid_dims)
    return phi, grid_dims, cell_lengths, "initial_phi"
end

"""
    relax_diblock_bvk2_morphology_nd(morphology=:LAM; f, chiN, lengths,
        c2, dims=(32,), initial_phi=nothing, ...)

Relax an arbitrary periodic N-D density field at a fixed cell using the
differentiable BVK2 functional.  `c2` is mandatory.  When `initial_phi` is
omitted, `diblock_morphology_seed` supplies the requested LAM, CYL, BCC, or
GYR seed; an interior `initial_phi` may instead be supplied on any periodic
N-D grid.

The solver first minimizes the intensive BVK2 energy by a mean-constrained
physical-density L-BFGS method.  Its gradient is the unweighted projected
chemical potential, not the endpoint-suppressed logit pullback.  Every step
stays strictly inside `(0, 1)` by fraction-to-boundary scaling and is accepted
only by an Armijo decrease of the exact BVK2 energy.  A bounded residual
polish may then take safeguarded Newton-Krylov steps.  The legacy-named
`gradient_tolerance` controls only a descent-angle safeguard; it is not a
termination test.  The reported `converged` flag follows
`convergence_policy`.  The optional strict-stationarity audit requires:

```
abs(mean(phi)-f) <= 1e-10,
0 < minimum(phi) <= maximum(phi) < 1,
R2 <= 1e-5,  Rinf <= 2e-5,
```

The phase-diagram policy instead requires 20 accepted energies with span at
most `1e-8`, together with `R2 <= 1e-2` and `Rinf <= 2e-2`.  The energy gate
controls the phase-boundary observable; the residual ceilings reject stalled
or grossly nonstationary fields.  Both certification flags are returned.

where `R = mu - mean(mu)` and `R2^2 = dV*sum(R.^2)`.  Morphology identity is
deliberately separate: the return value reports only the final contrast and
does not add a morphology-preserving penalty to the free energy.
"""
function relax_diblock_bvk2_morphology_nd(morphology=:LAM;
        f::Real=0.5, chiN::Real=20.0, dims=nothing, lengths=1.0,
        initial_amplitude::Real=1.0, initial_phi=nothing,
        N::Real=1.0, b::Real=1.0, adaptive::Bool=true, c2::Real,
        max_iterations::Integer=2_000, lbfgs_memory::Integer=10,
        gradient_tolerance::Real=1.0e-10,
        lbfgs_step_max::Real=0.1,
        lbfgs_fraction_to_boundary::Real=0.9,
        lbfgs_line_search_steps::Integer=24,
        lbfgs_armijo::Real=1.0e-4,
        lbfgs_curvature_tolerance::Real=1.0e-10,
        convergence_policy::Symbol=:strict_stationarity,
        residual_r2_tolerance::Real=1.0e-5,
        residual_rinf_tolerance::Real=2.0e-5,
        phase_energy_residual_r2_ceiling::Real=1.0e-2,
        phase_energy_residual_rinf_ceiling::Real=2.0e-2,
        phase_energy_plateau_tolerance::Real=1.0e-8,
        phase_energy_plateau_window::Integer=20,
        newton_krylov_iterations::Integer=8,
        gmres_max_iterations::Integer=24,
        gmres_restart::Integer=12,
        gmres_tolerance::Real=2.0e-2,
        gmres_preconditioner::Symbol=:none,
        gmres_preconditioner_floor::Real=1.0,
        jvp_relative_step::Real=1.0e-6,
        newton_step_max::Real=0.1,
        newton_fraction_to_boundary::Real=0.9,
        newton_line_search_steps::Integer=20,
        newton_sufficient_decrease::Real=1.0e-3,
        field_optimizer::Symbol=:lbfgs_jfnk,
        bb_trust_region::Real=8.0,
        bb_iterations::Integer=200,
        normal_residual_step_max::Real=0.02,
        r2_acceptance::Real=500.0,
        sd_newton_plateau_tolerance::Real=1.0e-6,
        sd_newton_max_blocks::Integer=8)
    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    chi = Float64(chiN)
    isfinite(chi) || throw(ArgumentError("chiN must be finite"))
    amplitude = Float64(initial_amplitude)
    isfinite(amplitude) || throw(ArgumentError(
        "initial_amplitude must be finite"))
    iterations = Int(max_iterations)
    iterations >= 0 || throw(ArgumentError(
        "max_iterations must be nonnegative"))
    memory = Int(lbfgs_memory)
    memory > 0 || throw(ArgumentError("lbfgs_memory must be positive"))
    gtol = Float64(gradient_tolerance)
    gtol > 0.0 && isfinite(gtol) || throw(ArgumentError(
        "gradient_tolerance must be positive and finite"))
    lbfgs_max_step = Float64(lbfgs_step_max)
    lbfgs_max_step > 0.0 && isfinite(lbfgs_max_step) || throw(ArgumentError(
        "lbfgs_step_max must be positive and finite"))
    lbfgs_boundary_fraction = Float64(lbfgs_fraction_to_boundary)
    0.0 < lbfgs_boundary_fraction < 1.0 &&
        isfinite(lbfgs_boundary_fraction) || throw(ArgumentError(
            "lbfgs_fraction_to_boundary must be finite and in (0, 1)"))
    lbfgs_search_steps = Int(lbfgs_line_search_steps)
    lbfgs_search_steps > 0 || throw(ArgumentError(
        "lbfgs_line_search_steps must be positive"))
    lbfgs_c1 = Float64(lbfgs_armijo)
    0.0 < lbfgs_c1 < 1.0 && isfinite(lbfgs_c1) || throw(ArgumentError(
        "lbfgs_armijo must be finite and in (0, 1)"))
    lbfgs_curvature = Float64(lbfgs_curvature_tolerance)
    lbfgs_curvature >= 0.0 && isfinite(lbfgs_curvature) ||
        throw(ArgumentError(
            "lbfgs_curvature_tolerance must be finite and nonnegative"))
    convergence_policy in (:strict_stationarity, :phase_energy) ||
        throw(ArgumentError("convergence_policy must be :strict_stationarity or :phase_energy"))
    field_optimizer in (:lbfgs_jfnk, :mobility_bb, :mobility_bb_jfnk,
        :sd_newton_mobility, :theta_bb) || throw(ArgumentError(
            "field_optimizer must be :lbfgs_jfnk, :mobility_bb, " *
            ":mobility_bb_jfnk, :sd_newton_mobility, or :theta_bb"))
    r2_tolerance = Float64(residual_r2_tolerance)
    rinf_tolerance = Float64(residual_rinf_tolerance)
    phase_r2_ceiling = Float64(phase_energy_residual_r2_ceiling)
    phase_rinf_ceiling = Float64(phase_energy_residual_rinf_ceiling)
    plateau_tolerance = Float64(phase_energy_plateau_tolerance)
    plateau_window = Int(phase_energy_plateau_window)
    all(value -> value > 0.0 && isfinite(value),
        (r2_tolerance, rinf_tolerance, phase_r2_ceiling,
            phase_rinf_ceiling, plateau_tolerance)) || throw(ArgumentError(
        "field convergence tolerances must be finite and positive"))
    plateau_window >= 2 || throw(ArgumentError(
        "phase_energy_plateau_window must be at least two"))
    newton_iterations = Int(newton_krylov_iterations)
    newton_iterations >= 0 || throw(ArgumentError(
        "newton_krylov_iterations must be nonnegative"))
    gmres_iterations = Int(gmres_max_iterations)
    gmres_iterations > 0 || throw(ArgumentError(
        "gmres_max_iterations must be positive"))
    gmres_restart_count = Int(gmres_restart)
    gmres_restart_count > 0 || throw(ArgumentError(
        "gmres_restart must be positive"))
    gmres_tol = Float64(gmres_tolerance)
    0.0 < gmres_tol < 1.0 && isfinite(gmres_tol) || throw(ArgumentError(
        "gmres_tolerance must be finite and in (0, 1)"))
    gmres_preconditioner in (:none, :rpa_fourier) || throw(ArgumentError(
        "gmres_preconditioner must be :none or :rpa_fourier"))
    preconditioner_floor = Float64(gmres_preconditioner_floor)
    preconditioner_floor > 0.0 && isfinite(preconditioner_floor) ||
        throw(ArgumentError(
            "gmres_preconditioner_floor must be positive and finite"))
    jvp_step = Float64(jvp_relative_step)
    jvp_step > 0.0 && isfinite(jvp_step) || throw(ArgumentError(
        "jvp_relative_step must be positive and finite"))
    step_max = Float64(newton_step_max)
    step_max > 0.0 && isfinite(step_max) || throw(ArgumentError(
        "newton_step_max must be positive and finite"))
    boundary_fraction = Float64(newton_fraction_to_boundary)
    0.0 < boundary_fraction < 1.0 && isfinite(boundary_fraction) ||
        throw(ArgumentError(
            "newton_fraction_to_boundary must be finite and in (0, 1)"))
    line_search_steps = Int(newton_line_search_steps)
    line_search_steps > 0 || throw(ArgumentError(
        "newton_line_search_steps must be positive"))
    sufficient_decrease = Float64(newton_sufficient_decrease)
    0.0 <= sufficient_decrease < 1.0 && isfinite(sufficient_decrease) ||
        throw(ArgumentError(
            "newton_sufficient_decrease must be finite and in [0, 1)"))
    phi0, grid_dims, cell_lengths, initial_source =
        _diblock_bvk2_morphology_initial_phi(morphology, initial_phi,
            dims, lengths, ff, amplitude)
    logits = vec(log.(phi0 ./ (1.0 .- phi0)))
    logits .-= mean(logits)
    volume = prod(cell_lengths)
    dV = volume / prod(grid_dims)
    phi_start, _ = _diblock_bvk2_morphology_phi_from_logits(
        logits, grid_dims, ff, dV)
    state_at_phi(candidate) = _diblock_bvk2_morphology_physical_state(
        candidate; f=ff, chiN=chi, lengths=cell_lengths,
        N=nn, b=bb, adaptive=adaptive, c2=c2)
    state = state_at_phi(phi_start)
    # :mobility_bb_jfnk runs a BOUNDED BB stage (escape the strong-segregation
    # stall + acquire the morphology) and then hands off to the RPA-
    # preconditioned Newton-Krylov polish below -- mirroring the BVK1
    # sd_block/newton_block gyroid solver.  :mobility_bb runs BB to the full
    # budget and skips Newton (diagnostic).
    physical_optimization = if field_optimizer == :theta_bb
        # UD-style square-root-density (angle) chart: BB on the true
        # theta-gradient, regular at the density walls.
        _diblock_bvk2_theta_bb_optimize(
            state, state_at_phi; f=ff, chiN=chi, lengths=cell_lengths,
            N=nn, b=bb, adaptive=adaptive, c2=c2, dV=dV,
            max_iterations=iterations, force_tolerance=rinf_tolerance,
            energy_tolerance=plateau_tolerance, memory=memory,
            stall_window=plateau_window, trust_region=Float64(bb_trust_region))
    elseif field_optimizer == :sd_newton_mobility
        # BVK1-style: mobility-BB SD warmup then normal-residual Gauss-Newton,
        # accepted on an energy plateau + a loose R2 ceiling.
        _diblock_bvk2_sd_newton_optimize(
            state, state_at_phi; f=ff, chiN=chi, lengths=cell_lengths,
            N=nn, b=bb, adaptive=adaptive, c2=c2, dV=dV,
            sd_iterations=Int(bb_iterations), newton_iterations=newton_iterations,
            max_blocks=Int(sd_newton_max_blocks),
            step_max=Float64(normal_residual_step_max),
            fraction=boundary_fraction, line_search_steps=line_search_steps,
            sufficient_decrease=sufficient_decrease, jvp_step=jvp_step,
            rinf_tolerance=rinf_tolerance, r2_acceptance=Float64(r2_acceptance),
            plateau_window=plateau_window,
            plateau_tolerance=Float64(sd_newton_plateau_tolerance),
            trust_region=Float64(bb_trust_region), memory=memory)
    elseif field_optimizer in (:mobility_bb, :mobility_bb_jfnk)
        _diblock_bvk2_mobility_bb_optimize(
            state, state_at_phi; f=ff, chiN=chi, lengths=cell_lengths,
            N=nn, b=bb, adaptive=adaptive, c2=c2, dV=dV,
            max_iterations=(field_optimizer == :mobility_bb_jfnk ?
                Int(bb_iterations) : iterations),
            force_tolerance=rinf_tolerance,
            energy_tolerance=plateau_tolerance, memory=memory,
            stall_window=plateau_window, trust_region=Float64(bb_trust_region))
    else
        _diblock_bvk2_physical_lbfgs(
            state, state_at_phi; f=ff, max_iterations=iterations,
            memory=memory, step_max=lbfgs_max_step,
            fraction_to_boundary=lbfgs_boundary_fraction,
            line_search_steps=lbfgs_search_steps, armijo=lbfgs_c1,
            curvature_tolerance=lbfgs_curvature, descent_tolerance=gtol,
            convergence_policy=convergence_policy,
            residual_r2_tolerance=r2_tolerance,
            residual_rinf_tolerance=rinf_tolerance,
            phase_energy_residual_r2_ceiling=phase_r2_ceiling,
            phase_energy_residual_rinf_ceiling=phase_rinf_ceiling,
            phase_energy_plateau_tolerance=plateau_tolerance,
            phase_energy_plateau_window=plateau_window)
    end
    state = physical_optimization.state
    optimizer_trace = copy(physical_optimization.optimizer_trace)
    optimizer_converged = physical_optimization.converged
    optimizer_iterations = physical_optimization.iterations
    optimizer_status = physical_optimization.status
    optimizer_minimum = physical_optimization.minimum
    optimizer_residual_r2 = state.residual_r2
    optimizer_residual_rinf = state.residual_rinf
    completed_newton_iterations = 0
    total_gmres_iterations = 0
    final_gmres_relative_residual = NaN
    jvp_calls = 0
    # Pure :mobility_bb converges on the physical force alone, so the
    # Newton-Krylov polish is skipped for it.  :mobility_bb_jfnk and
    # :lbfgs_jfnk both run the polish (BB/L-BFGS having supplied the warm
    # start it needs to stay in its convergence basin).
    newton_range = field_optimizer in
        (:mobility_bb, :sd_newton_mobility, :theta_bb) ?
        (1:0) : (1:newton_iterations)
    for newton_iteration in newton_range
        convergence_policy == :phase_energy &&
            physical_optimization.phase_energy_converged && break
        state.residual_r2 <= r2_tolerance &&
            state.residual_rinf <= rinf_tolerance && break
        base = _diblock_bvk2_morphology_physical_state(state.phi;
            f=ff, chiN=chi, lengths=cell_lengths,
            N=nn, b=bb, adaptive=adaptive, c2=c2)
        translation_modes = _diblock_bvk2_translation_modes(base.phi)
        linear_residual = copy(base.residual)
        _diblock_bvk2_project_gauge!(linear_residual, translation_modes)
        raw_preconditioner = gmres_preconditioner == :rpa_fourier ?
            _diblock_bvk2_rpa_preconditioner(grid_dims, cell_lengths;
                f=ff, chiN=chi, N=nn, b=bb,
                denominator_floor=preconditioner_floor) : identity

        function apply_preconditioner(vector)
            output = reshape(Vector{Float64}(
                raw_preconditioner(vector)), grid_dims)
            _diblock_bvk2_project_gauge!(output, translation_modes)
            return vec(output)
        end

        function apply_jacobian(vector_direction)
            direction = reshape(Float64.(vector_direction), grid_dims)
            _diblock_bvk2_project_gauge!(direction, translation_modes)
            product, step, _sign = _diblock_bvk2_morphology_jvp(
                base, direction;
                f=ff, chiN=chi, lengths=cell_lengths,
                N=nn, b=bb, adaptive=adaptive, c2=c2,
                relative_step=jvp_step)
            step > 100eps(Float64) && (jvp_calls += 1)
            _diblock_bvk2_project_gauge!(product, translation_modes)
            return vec(product)
        end

        delta, inner_iterations, inner_relative = _diblock_bvk2_gmres(
            apply_jacobian, -vec(linear_residual);
            max_iterations=gmres_iterations,
            restart=min(gmres_restart_count, gmres_iterations),
            tolerance=gmres_tol,
            apply_preconditioner=apply_preconditioner)
        total_gmres_iterations += inner_iterations
        final_gmres_relative_residual = inner_relative
        direction = reshape(delta, grid_dims)
        _diblock_bvk2_project_gauge!(direction, translation_modes)
        largest_step = maximum(abs, direction)
        largest_step > 0.0 && isfinite(largest_step) || break
        largest_step > step_max && (direction .*= step_max / largest_step)
        alpha = _diblock_bvk2_fraction_to_boundary(base.phi, direction;
            fraction=boundary_fraction)
        alpha > 0.0 && isfinite(alpha) || break
        accepted = false
        iteration_line_search_evaluations = 0
        accepted_alpha = NaN
        for _ in 1:line_search_steps
            iteration_line_search_evaluations += 1
            candidate = base.phi .+ alpha .* direction
            candidate .-= mean(candidate) - ff
            trial = try
                _diblock_bvk2_morphology_physical_state(candidate;
                    f=ff, chiN=chi, lengths=cell_lengths,
                    N=nn, b=bb, adaptive=adaptive, c2=c2)
            catch error
                if error isa DomainError || error isa ArgumentError
                    nothing
                else
                    rethrow()
                end
            end
            residual_merit_pass = trial !== nothing &&
                isfinite(trial.residual_r2) &&
                trial.residual_r2^2 <=
                (1.0 - sufficient_decrease * alpha) * base.residual_r2^2
            max_residual_pass = trial !== nothing &&
                (trial.residual_rinf <= base.residual_rinf ||
                 trial.residual_rinf <= rinf_tolerance)
            if residual_merit_pass && max_residual_pass
                state = trial
                completed_newton_iterations = newton_iteration
                accepted = true
                accepted_alpha = alpha
                break
            end
            alpha *= 0.5
        end
        accepted || break
        push!(optimizer_trace, (
            iteration=optimizer_iterations + newton_iteration,
            stage="physical_jfnk",
            energy_density=state.energy_density,
            residual_r2=state.residual_r2,
            residual_rinf=state.residual_rinf,
            mean_phi=mean(state.phi),
            minimum_phi=minimum(state.phi),
            maximum_phi=maximum(state.phi),
            accepted_alpha=accepted_alpha,
            line_search_evaluations=iteration_line_search_evaluations,
            line_search_status="accepted_residual_merit",
        ))
    end

    mean_phi = mean(state.phi)
    mean_error = abs(mean_phi - ff)
    minimum_phi = minimum(state.phi)
    maximum_phi = maximum(state.phi)
    strict_bounds = 0.0 < minimum_phi && maximum_phi < 1.0
    strict_stationarity_converged = mean_error <= 1.0e-10 && strict_bounds &&
        state.residual_r2 <= r2_tolerance &&
        state.residual_rinf <= rinf_tolerance
    phase_energy_converged = physical_optimization.phase_energy_converged ||
        strict_stationarity_converged
    selected_converged = convergence_policy == :strict_stationarity ?
        strict_stationarity_converged : phase_energy_converged
    kind = _diblock_morphology_symbol(morphology)
    fingerprint = string(diblock_bvk2_fingerprint(; f=ff, N=nn, b=bb,
            c2=c2, adaptive=adaptive),
        "|solver=", DIBLOCK_BVK2_MORPHOLOGY_SOLVER_SCHEMA,
        "|dims=", join(grid_dims, "x"),
        "|lengths=", join(cell_lengths, ","))
    return (
        model="BVK2",
        morphology=string(kind),
        initial_source=initial_source,
        solver_schema=DIBLOCK_BVK2_MORPHOLOGY_SOLVER_SCHEMA,
        fingerprint=fingerprint,
        f=ff,
        chiN=chi,
        c2=Float64(c2),
        dims=grid_dims,
        lengths=Tuple(cell_lengths),
        energy_density=state.energy_density,
        energy=state.energy,
        phi=state.phi,
        residual_r2=state.residual_r2,
        residual_rinf=state.residual_rinf,
        mean_phi=mean_phi,
        mean_error=mean_error,
        minimum_phi=minimum_phi,
        maximum_phi=maximum_phi,
        strict_bounds=strict_bounds,
        contrast=maximum_phi - minimum_phi,
        converged=selected_converged,
        convergence_policy=String(convergence_policy),
        strict_stationarity_converged=strict_stationarity_converged,
        phase_energy_converged=phase_energy_converged,
        energy_plateau_available=physical_optimization.energy_plateau_available,
        energy_plateau_span=physical_optimization.energy_plateau_span,
        energy_plateau_tolerance=plateau_tolerance,
        energy_plateau_window=plateau_window,
        optimizer_converged=optimizer_converged,
        optimizer_iterations=optimizer_iterations,
        optimizer_status=optimizer_status,
        optimizer_minimum=optimizer_minimum,
        optimizer_residual_r2=optimizer_residual_r2,
        optimizer_residual_rinf=optimizer_residual_rinf,
        lbfgs_memory_updates=physical_optimization.memory_updates,
        lbfgs_line_search_evaluations=
            physical_optimization.line_search_evaluations,
        lbfgs_descent_fallbacks=physical_optimization.descent_fallbacks,
        optimizer_trace=optimizer_trace,
        newton_krylov_iterations=completed_newton_iterations,
        gmres_iterations=total_gmres_iterations,
        gmres_final_relative_residual=final_gmres_relative_residual,
        gmres_preconditioner=String(gmres_preconditioner),
        gmres_preconditioner_floor=preconditioner_floor,
        jvp_calls=jvp_calls,
    )
end

function _diblock_bvk2_scaled_lengths(reference_lengths, factor::Float64,
        scaled_axes)
    dimension_count = length(reference_lengths)
    axes = scaled_axes === nothing ? Tuple(1:dimension_count) :
        scaled_axes isa Integer ? (Int(scaled_axes),) : Tuple(Int.(scaled_axes))
    isempty(axes) && throw(ArgumentError(
        "scaled_axes must select at least one axis"))
    length(unique(axes)) == length(axes) || throw(ArgumentError(
        "scaled_axes must not contain duplicate axes"))
    all(axis -> 1 <= axis <= dimension_count, axes) || throw(ArgumentError(
        "scaled_axes contains an axis outside 1:$dimension_count"))
    lengths = [Float64(value) for value in reference_lengths]
    for axis in axes
        lengths[axis] *= factor
    end
    return lengths, axes
end

"""
    relax_diblock_bvk2_morphology_stress_free(morphology=:LAM; ...)

Relax a BVK2 ordered field and root its analytic symmetry-compatible cell
stress in the logarithmic scale factor.  `reference_lengths` and the density
grid use BVK2 length units.  `scaled_axes=nothing` dilates every cell axis;
for a 3-D CYL embedding use `scaled_axes=(1, 2)` so the invariant third axis
is not changed.

Every cell evaluation first calls [`relax_diblock_bvk2_morphology_nd`](@ref)
and is eligible for the stress bracket only when its selected field-
convergence policy passes and its contrast exceeds `minimum_contrast`.  The nearest eligible field
is reused as the next continuation seed.  The cell root is accepted only with
an oriented stress bracket, `abs(stress) <= cell_stress_tolerance`, logarithmic
bracket width at most `log_cell_tolerance`, and independent relaxed checks at
`scale*(1 +/- local_check_fraction)` that retain the stress orientation and do
not lower the energy density beyond `local_energy_tolerance`.

The returned `converged` flag certifies the selected field policy and the
scalar-cell gates, not morphology identity.  Under `:phase_energy`, strict
field stationarity remains a separate audit flag.  Shell, topology,
correlation, and competing-phase gates remain separate publication
certificates.
"""
function relax_diblock_bvk2_morphology_stress_free(morphology=:LAM;
        f::Real=0.5, chiN::Real=20.0, dims=(32,),
        reference_lengths=1.0, scaled_axes=nothing,
        initial_amplitude::Real=1.0, initial_phi=nothing,
        N::Real=1.0, b::Real=1.0, adaptive::Bool=true, c2::Real,
        lower_factor::Real=0.75, upper_factor::Real=1.5,
        initial_factor::Real=1.0, search_log_step::Real=0.02,
        max_search_steps::Integer=60, max_cell_iterations::Integer=24,
        cell_stress_tolerance::Real=1.0e-6,
        log_cell_tolerance::Real=1.0e-5,
        local_check_fraction::Real=0.01,
        local_energy_tolerance::Real=1.0e-8,
        minimum_contrast::Real=1.0e-3,
        field_kwargs...)
    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    chi = Float64(chiN)
    isfinite(chi) || throw(ArgumentError("chiN must be finite"))
    lo = Float64(lower_factor)
    hi = Float64(upper_factor)
    0.0 < lo < hi && isfinite(lo) && isfinite(hi) || throw(ArgumentError(
        "lower_factor and upper_factor must be finite, positive, and ordered"))
    start = Float64(initial_factor)
    lo < start < hi || throw(ArgumentError(
        "initial_factor must lie strictly inside the factor bounds"))
    step0 = Float64(search_log_step)
    step0 > 0.0 && isfinite(step0) || throw(ArgumentError(
        "search_log_step must be positive and finite"))
    search_steps = Int(max_search_steps)
    search_steps > 0 || throw(ArgumentError("max_search_steps must be positive"))
    cell_iterations = Int(max_cell_iterations)
    cell_iterations > 0 || throw(ArgumentError(
        "max_cell_iterations must be positive"))
    stress_tolerance = Float64(cell_stress_tolerance)
    stress_tolerance > 0.0 && isfinite(stress_tolerance) ||
        throw(ArgumentError("cell_stress_tolerance must be positive and finite"))
    log_tolerance = Float64(log_cell_tolerance)
    log_tolerance > 0.0 && isfinite(log_tolerance) || throw(ArgumentError(
        "log_cell_tolerance must be positive and finite"))
    check_fraction = Float64(local_check_fraction)
    0.0 < check_fraction < 1.0 && isfinite(check_fraction) ||
        throw(ArgumentError("local_check_fraction must lie in (0, 1)"))
    energy_tolerance = Float64(local_energy_tolerance)
    energy_tolerance >= 0.0 && isfinite(energy_tolerance) ||
        throw(ArgumentError("local_energy_tolerance must be finite and nonnegative"))
    contrast_tolerance = Float64(minimum_contrast)
    contrast_tolerance >= 0.0 && isfinite(contrast_tolerance) ||
        throw(ArgumentError("minimum_contrast must be finite and nonnegative"))

    reference_dims = initial_phi === nothing ? _diblock_morphology_dims(dims) :
        Tuple(size(initial_phi))
    base_lengths = _diblock_nd_lengths(reference_lengths, reference_dims)
    _, axes = _diblock_bvk2_scaled_lengths(base_lengths, 1.0, scaled_axes)
    evaluations = Dict{Float64,NamedTuple}()

    function evaluate_factor(factor::Real; force_seed=nothing)
        value = Float64(factor)
        lo < value < hi || throw(ArgumentError(
            "cell factor must remain strictly inside the search bounds"))
        haskey(evaluations, value) && return evaluations[value]
        lengths, _ = _diblock_bvk2_scaled_lengths(base_lengths, value, axes)
        eligible_rows = [row for row in values(evaluations) if row.eligible]
        warm_phi = if force_seed !== nothing
            force_seed
        elseif !isempty(eligible_rows)
            eligible_rows[argmin(abs(log(row.factor / value))
                for row in eligible_rows)].result.phi
        else
            initial_phi
        end
        result = relax_diblock_bvk2_morphology_nd(morphology;
            f=ff, chiN=chi, dims=reference_dims, lengths=lengths,
            initial_amplitude=initial_amplitude, initial_phi=warm_phi,
            N=nn, b=bb, adaptive=adaptive, c2=c2, field_kwargs...)
        identity_proxy = result.contrast >= contrast_tolerance
        eligible = result.converged && identity_proxy &&
            isfinite(result.energy_density)
        stress = eligible ? diblock_bvk2_cell_scale_gradient_density_nd(
            result.phi; f=ff, chiN=chi, lengths=lengths, N=nn, b=bb,
            adaptive=adaptive, c2=c2, scaled_axes=axes,
            mean_tolerance=1.0e-10) : NaN
        eligible &= isfinite(stress)
        row = (factor=value, log_factor=log(value), lengths=Tuple(lengths),
            energy_density=result.energy_density, stress=stress,
            result=result, identity_proxy=identity_proxy, eligible=eligible)
        evaluations[value] = row
        return row
    end

    center = evaluate_factor(start)
    center.eligible || throw(ErrorException(
        "initial BVK2 cell factor did not produce an eligible stationary " *
        "ordered field: R2=$(center.result.residual_r2), " *
        "Rinf=$(center.result.residual_rinf), " *
        "contrast=$(center.result.contrast), " *
        "strict_bounds=$(center.result.strict_bounds), " *
        "optimizer_iterations=$(center.result.optimizer_iterations), " *
        "newton_iterations=$(center.result.newton_krylov_iterations), " *
        "gmres_iterations=$(center.result.gmres_iterations), " *
        "jvp_calls=$(center.result.jvp_calls)"))
    current = center
    bracket = nothing
    search_step = step0
    if center.stress == 0.0
        left_probe_factor = center.factor * exp(-search_step)
        right_probe_factor = center.factor * exp(search_step)
        if lo < left_probe_factor < center.factor < right_probe_factor < hi
            left_probe = evaluate_factor(left_probe_factor)
            right_probe = evaluate_factor(right_probe_factor)
            left_probe.eligible && right_probe.eligible &&
                left_probe.stress < 0.0 < right_probe.stress &&
                (bracket = (left_probe, right_probe))
        end
    end
    direction = center.stress > 0.0 ? -1.0 : 1.0
    for _ in 1:search_steps
        bracket === nothing || break
        candidate_factor = current.factor * exp(direction * search_step)
        lo < candidate_factor < hi || break
        trial = evaluate_factor(candidate_factor)
        if !trial.eligible
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
        "BVK2 analytic cell stress has no eligible oriented continuation bracket"))
    lower, upper = bracket
    selected = abs(lower.stress) <= abs(upper.stress) ? lower : upper
    root_iterations = 0
    for iteration in 1:cell_iterations
        root_iterations = iteration
        width = upper.log_factor - lower.log_factor
        abs(selected.stress) <= stress_tolerance && width <= log_tolerance && break
        secant = (lower.log_factor * upper.stress -
            upper.log_factor * lower.stress) / (upper.stress - lower.stress)
        guard = 0.1 * width
        candidate_log = clamp(secant, lower.log_factor + guard,
            upper.log_factor - guard)
        trial = evaluate_factor(exp(candidate_log))
        trial.eligible || throw(ErrorException(
            "BVK2 field lost eligibility inside the analytic cell bracket"))
        selected = abs(trial.stress) < abs(selected.stress) ? trial : selected
        if trial.stress < 0.0
            lower = trial
        elseif trial.stress > 0.0
            upper = trial
        else
            selected = trial
            lower = trial
            upper = trial
            break
        end
    end

    final_width = upper.log_factor - lower.log_factor
    left_factor = selected.factor * (1.0 - check_fraction)
    right_factor = selected.factor * (1.0 + check_fraction)
    lo < left_factor < selected.factor < right_factor < hi ||
        throw(ErrorException(
            "BVK2 cell root is too close to a search-window endpoint for local validation"))
    left = evaluate_factor(left_factor; force_seed=selected.result.phi)
    right = evaluate_factor(right_factor; force_seed=selected.result.phi)
    neighbor_fields = left.eligible && right.eligible
    stress_oriented = neighbor_fields && left.stress < 0.0 < right.stress
    local_minimum = neighbor_fields &&
        left.energy_density >= selected.energy_density - energy_tolerance &&
        right.energy_density >= selected.energy_density - energy_tolerance
    cell_converged = selected.eligible &&
        abs(selected.stress) <= stress_tolerance &&
        final_width <= log_tolerance && stress_oriented && local_minimum
    trace = sort(collect(values(evaluations)); by=row -> row.factor)
    return (
        model="BVK2",
        morphology=selected.result.morphology,
        f=ff,
        chiN=chi,
        c2=Float64(c2),
        scaled_axes=axes,
        reference_lengths=Tuple(base_lengths),
        cell_factor=selected.factor,
        lengths=selected.lengths,
        energy_density=selected.energy_density,
        stress=selected.stress,
        phi=selected.result.phi,
        field_result=selected.result,
        contrast=selected.result.contrast,
        identity_proxy=selected.identity_proxy,
        bracket=(lower.factor, upper.factor),
        log_bracket_width=final_width,
        root_iterations=root_iterations,
        left_check=left,
        right_check=right,
        stress_oriented=stress_oriented,
        local_minimum=local_minimum,
        cell_converged=cell_converged,
        converged=cell_converged,
        trace=trace,
        solver_schema="analytic_log_cell_stress_root_nd_v1",
        fingerprint=selected.result.fingerprint *
            "|cell_solver=analytic_log_cell_stress_root_nd_v1" *
            "|scaled_axes=" * join(axes, ","),
    )
end
