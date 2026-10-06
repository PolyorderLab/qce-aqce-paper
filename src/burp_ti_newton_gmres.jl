"""
    DiblockBurpTiForceWorkspaceND(dims)

Preallocated scratch buffers and FFTW plans for repeated BURP-TI
chemical-potential evaluations on a fixed N-D periodic grid.
"""
struct DiblockBurpTiForceWorkspaceND{D,PF,PI}
    dims::NTuple{D,Int}
    psi::Array{Float64,D}
    phi_lambda::Array{Float64,D}
    eta::Array{Float64,D}
    kernel_psi::Array{Float64,D}
    kernel_eta::Array{Float64,D}
    field_c::Array{ComplexF64,D}
    field_hat::Array{ComplexF64,D}
    kernel_hat::Array{ComplexF64,D}
    kernel_time::Array{ComplexF64,D}
    fplan::PF
    iplan::PI
end

function DiblockBurpTiForceWorkspaceND(dims::NTuple{D,Int}) where {D}
    all(count -> count >= 2, dims) ||
        throw(ArgumentError("BURP-TI N-D force workspace requires at least two grid points per dimension"))
    dummy = zeros(ComplexF64, dims...)
    fplan = plan_fft(dummy)
    iplan = plan_ifft(dummy)
    return DiblockBurpTiForceWorkspaceND{D,typeof(fplan),typeof(iplan)}(dims,
        zeros(Float64, dims...), zeros(Float64, dims...), zeros(Float64, dims...),
        zeros(Float64, dims...), zeros(Float64, dims...),
        zeros(ComplexF64, dims...), zeros(ComplexF64, dims...),
        zeros(ComplexF64, dims...), zeros(ComplexF64, dims...), fplan, iplan)
end

function DiblockBurpTiForceWorkspaceND(dims::AbstractVector{<:Integer})
    return DiblockBurpTiForceWorkspaceND(Tuple(Int.(dims)))
end

function _burp_ti_kernel_apply_nd!(out::AbstractArray{Float64,D},
        values::AbstractArray{Float64,D}, kern::AbstractArray{<:Real,D},
        ws::DiblockBurpTiForceWorkspaceND{D}) where {D}
    size(values) == ws.dims ||
        throw(DimensionMismatch("values dims $(size(values)) do not match workspace dims $(ws.dims)"))
    size(out) == ws.dims ||
        throw(DimensionMismatch("output dims $(size(out)) do not match workspace dims $(ws.dims)"))
    size(kern) == ws.dims ||
        throw(DimensionMismatch("kernel dims $(size(kern)) do not match workspace dims $(ws.dims)"))
    field_c = ws.field_c
    @inbounds for idx in eachindex(values)
        field_c[idx] = values[idx]
    end
    mul!(ws.field_hat, ws.fplan, field_c)
    first = firstindex(ws.kernel_hat)
    @inbounds for idx in eachindex(ws.kernel_hat)
        ws.kernel_hat[idx] = idx == first ? 0.0 : kern[idx] * ws.field_hat[idx]
    end
    mul!(ws.kernel_time, ws.iplan, ws.kernel_hat)
    @inbounds for idx in eachindex(out)
        out[idx] = real(ws.kernel_time[idx])
    end
    return out
end

@inline _burp_ti_eta_prime(value::Real, f::Real, fb::Real) =
    0.5 * (f / value + fb / (1.0 - value))

@inline _burp_ti_eta_double_prime(value::Real, f::Real, fb::Real) =
    0.5 * (-f / (value * value) + fb / ((1.0 - value) * (1.0 - value)))

"""
    diblock_burp_ti_chemical_potential_nd(phi_a; f, chiN, lengths, ...)

Return the unconstrained BURP-TI chemical potential on an N-D periodic density
grid. Morphology-specific choices are represented only by `phi_a` and
`lengths`; the implementation is shared by all phases.
"""
function diblock_burp_ti_chemical_potential_nd(phi_a::AbstractArray{<:Real};
        f::Real=0.5, chiN::Real=12.0, lengths=1.0, N::Real=1.0, b::Real=1.0,
        nquad::Integer=16, mean_tolerance::Real=1.0e-8,
        composition_kernel=nothing, workspace=nothing)
    phi, dims, ff, chi, cell_lengths, nn, bb, _ =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    qcount = Int(nquad)
    qcount >= 2 || throw(ArgumentError("nquad must be at least 2"))
    fb = 1.0 - ff
    kern = if composition_kernel === nothing
        k2_values = _periodic_k2_nd(dims, cell_lengths)
        map(k2 -> diblock_composition_kernel(k2; f=ff, N=nn, b=bb), k2_values)
    else
        size(composition_kernel) == dims ||
            throw(DimensionMismatch("precomputed composition kernel dims must match phi dims"))
        composition_kernel
    end

    ws = workspace === nothing ? DiblockBurpTiForceWorkspaceND(dims) : workspace
    ws.dims == dims ||
        throw(DimensionMismatch("workspace dims $(ws.dims) do not match phi dims $(dims)"))

    psi = ws.psi
    @inbounds for idx in eachindex(phi)
        psi[idx] = phi[idx] - ff
    end
    _burp_ti_kernel_apply_nd!(ws.kernel_psi, psi, kern, ws)

    chemical = zeros(Float64, dims...)
    for q in 1:qcount
        lambda = (q - 0.5) / qcount
        @inbounds for idx in eachindex(phi)
            ws.phi_lambda[idx] = ff + lambda * psi[idx]
            ws.eta[idx] = diblock_burp_eta_psi(ws.phi_lambda[idx], ff)
        end
        _burp_ti_kernel_apply_nd!(ws.kernel_eta, ws.eta, kern, ws)
        @inbounds for idx in eachindex(phi)
            eta_prime = _burp_ti_eta_prime(ws.phi_lambda[idx], ff, fb)
            chemical[idx] += ws.kernel_eta[idx] +
                lambda * eta_prime * ws.kernel_psi[idx]
        end
    end
    chemical ./= qcount
    @. chemical += chi * (1.0 - 2.0 * phi)
    return chemical
end

struct DiblockBurpTiProjectedForce{D}
    residual::Array{Float64,D}
    chemical_potential::Array{Float64,D}
    norm::Float64
    maxabs::Float64
end

struct DiblockBurpTiBoundedProjectedForce{D}
    residual::Array{Float64,D}
    chemical_potential::Array{Float64,D}
    norm::Float64
    maxabs::Float64
    multiplier::Float64
    lower_active_count::Int
    upper_active_count::Int
    free_count::Int
end

function _diblock_projected_force_result(chemical::Array{Float64,D},
        dV::Float64) where {D}
    residual = chemical .- mean(chemical)
    return DiblockBurpTiProjectedForce{D}(residual, chemical,
        sqrt(dV * sum(abs2, residual)), maximum(abs, residual))
end

function _diblock_bounded_projected_force_result(chemical::Array{Float64,D},
        phi::Array{Float64,D}, dV::Float64; lower::Real, upper::Real,
        active_tolerance::Real) where {D}
    lo = Float64(lower)
    hi = Float64(upper)
    active_tol = Float64(active_tolerance)
    lower_active = phi .<= lo + active_tol
    upper_active = phi .>= hi - active_tol
    lower_active_count = count(lower_active)
    upper_active_count = count(upper_active)

    if lower_active_count == 0 && upper_active_count == 0
        multiplier = mean(chemical)
        residual = chemical .- multiplier
        return DiblockBurpTiBoundedProjectedForce{D}(residual, chemical,
            sqrt(dV * sum(abs2, residual)), maximum(abs, residual), multiplier,
            0, 0, length(chemical))
    end

    free = .!(lower_active .| upper_active)

    function violation_sum(multiplier::Float64)
        total = 0.0
        @inbounds for idx in eachindex(chemical)
            raw = chemical[idx] - multiplier
            if lower_active[idx]
                total += min(raw, 0.0)
            elseif upper_active[idx]
                total += max(raw, 0.0)
            else
                total += raw
            end
        end
        return total
    end

    spread = maximum(chemical) - minimum(chemical)
    pad = max(spread, 1.0)
    low = minimum(chemical) - pad
    high = maximum(chemical) + pad
    for _ in 1:100
        mid = 0.5 * (low + high)
        if violation_sum(mid) > 0.0
            low = mid
        else
            high = mid
        end
    end
    multiplier = 0.5 * (low + high)

    residual = similar(chemical)
    @inbounds for idx in eachindex(chemical)
        raw = chemical[idx] - multiplier
        if lower_active[idx]
            residual[idx] = min(raw, 0.0)
        elseif upper_active[idx]
            residual[idx] = max(raw, 0.0)
        else
            residual[idx] = raw
        end
    end
    return DiblockBurpTiBoundedProjectedForce{D}(residual, chemical,
        sqrt(dV * sum(abs2, residual)), maximum(abs, residual), multiplier,
        lower_active_count, upper_active_count, count(free))
end

function _diblock_burp_ti_bounded_newton_movable_mask(
        bounded::DiblockBurpTiBoundedProjectedForce{D},
        phi::AbstractArray{<:Real,D}; lower::Real, upper::Real,
        active_tolerance::Real) where {D}
    size(phi) == size(bounded.residual) ||
        throw(DimensionMismatch("phi and bounded residual dimensions must match"))
    lo = Float64(lower)
    hi = Float64(upper)
    active_tol = Float64(active_tolerance)
    movable = trues(size(phi))
    @inbounds for idx in eachindex(phi)
        at_bound = phi[idx] <= lo + active_tol ||
            phi[idx] >= hi - active_tol
        movable[idx] = !at_bound || bounded.residual[idx] != 0.0
    end
    return movable
end

function _diblock_burp_ti_project_bounded_tangent!(values::AbstractArray,
        movable::AbstractArray{Bool})
    size(values) == size(movable) ||
        throw(DimensionMismatch("values and movable mask dimensions must match"))
    movable_count = count(movable)
    movable_count == 0 && return fill!(values, zero(eltype(values)))
    movable_sum = zero(eltype(values))
    @inbounds for idx in eachindex(values)
        if movable[idx]
            movable_sum += values[idx]
        else
            values[idx] = zero(eltype(values))
        end
    end
    movable_mean = movable_sum / movable_count
    @inbounds for idx in eachindex(values)
        movable[idx] && (values[idx] -= movable_mean)
    end
    return values
end

"""
    diblock_burp_ti_projected_force_nd(phi_a; f, chiN, lengths, ...)

Return the fixed-composition stationarity residual `mu - mean(mu)`, its
cell-volume weighted L2 norm, and its max absolute component.
"""
function diblock_burp_ti_projected_force_nd(phi_a::AbstractArray{<:Real};
        f::Real=0.5, chiN::Real=12.0, lengths=1.0, N::Real=1.0, b::Real=1.0,
        nquad::Integer=16, mean_tolerance::Real=1.0e-8,
        composition_kernel=nothing, workspace=nothing)
    chemical = diblock_burp_ti_chemical_potential_nd(phi_a; f=f, chiN=chiN,
        lengths=lengths, N=N, b=b, nquad=nquad,
        mean_tolerance=mean_tolerance, composition_kernel=composition_kernel,
        workspace=workspace)
    cell_lengths = _diblock_nd_lengths(lengths, size(chemical))
    dV = prod(cell_lengths) / length(chemical)
    return _diblock_projected_force_result(chemical, dV)
end

"""
    diblock_burp_ti_bounded_projected_force_nd(phi_a; f, chiN, lengths, ...)

Return the KKT residual for fixed-composition stationarity with box bounds.
Interior cells use `mu - lambda`; lower-bound cells only count negative
violations, and upper-bound cells only count positive violations. The scalar
multiplier `lambda` is chosen to minimize this bounded residual.
"""
function diblock_burp_ti_bounded_projected_force_nd(
        phi_a::AbstractArray{<:Real}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, N::Real=1.0, b::Real=1.0, nquad::Integer=16,
        mean_tolerance::Real=1.0e-8, density_floor::Real=1.0e-6,
        active_tolerance::Real=10.0 * Float64(density_floor),
        composition_kernel=nothing, workspace=nothing)
    phi = Float64.(Array(phi_a))
    chemical = diblock_burp_ti_chemical_potential_nd(phi; f=f, chiN=chiN,
        lengths=lengths, N=N, b=b, nquad=nquad,
        mean_tolerance=mean_tolerance, composition_kernel=composition_kernel,
        workspace=workspace)
    cell_lengths = _diblock_nd_lengths(lengths, size(chemical))
    dV = prod(cell_lengths) / length(chemical)
    floor = Float64(density_floor)
    return _diblock_bounded_projected_force_result(chemical, phi, dV;
        lower=floor, upper=1.0 - floor,
        active_tolerance=Float64(active_tolerance))
end

function _diblock_burp_ti_context(dims::NTuple{D,Int},
        cell_lengths::AbstractVector{Float64}; f::Float64, chiN::Float64,
        N::Float64, b::Float64, nquad::Integer) where {D}
    qcount = Int(nquad)
    qcount >= 2 || throw(ArgumentError("nquad must be at least 2"))
    k2_values = _periodic_k2_nd(dims, cell_lengths)
    return (
        f=f,
        fb=1.0 - f,
        chiN=chiN,
        nquad=qcount,
        lengths=Float64.(collect(cell_lengths)),
        dims=dims,
        dV=prod(cell_lengths) / prod(dims),
        composition_kernel=map(k2 -> diblock_composition_kernel(k2; f=f,
            N=N, b=b), k2_values),
        energy_workspace=DiblockBurpTiWorkspaceND(dims),
        force_workspace=DiblockBurpTiForceWorkspaceND(dims),
        jv_workspace=DiblockBurpTiForceWorkspaceND(dims),
        preconditioner_workspace=DiblockBurpTiForceWorkspaceND(dims),
        kpsi=zeros(Float64, dims...),
        kv=zeros(Float64, dims...),
        tmp=zeros(Float64, dims...),
        jout=zeros(Float64, dims...),
    )
end

function _diblock_burp_ti_projected_jv!(out::AbstractArray{Float64,D},
        phi::AbstractArray{Float64,D}, direction::AbstractArray{Float64,D},
        ctx) where {D}
    ws = ctx.jv_workspace
    ff = ctx.f
    fb = ctx.fb
    @inbounds for idx in eachindex(phi)
        ws.psi[idx] = phi[idx] - ff
    end
    _burp_ti_kernel_apply_nd!(ctx.kpsi, ws.psi, ctx.composition_kernel, ws)
    _burp_ti_kernel_apply_nd!(ctx.kv, direction, ctx.composition_kernel, ws)
    fill!(out, 0.0)
    for q in 1:ctx.nquad
        lambda = (q - 0.5) / ctx.nquad
        @inbounds for idx in eachindex(phi)
            phi_lambda = ff + lambda * ws.psi[idx]
            ctx.tmp[idx] = _burp_ti_eta_prime(phi_lambda, ff, fb) * direction[idx]
        end
        _burp_ti_kernel_apply_nd!(ws.kernel_eta, ctx.tmp,
            ctx.composition_kernel, ws)
        @inbounds for idx in eachindex(phi)
            phi_lambda = ff + lambda * ws.psi[idx]
            eta_prime = _burp_ti_eta_prime(phi_lambda, ff, fb)
            eta_double = _burp_ti_eta_double_prime(phi_lambda, ff, fb)
            out[idx] += lambda * ws.kernel_eta[idx] +
                lambda * eta_prime * ctx.kv[idx] +
                lambda * lambda * eta_double * direction[idx] * ctx.kpsi[idx]
        end
    end
    @inbounds for idx in eachindex(out)
        out[idx] = out[idx] / ctx.nquad - 2.0 * ctx.chiN * direction[idx]
    end
    out .-= mean(out)
    return out
end

function _diblock_burp_ti_rpa_preconditioner!(
        out::AbstractArray{Float64,D}, values::AbstractArray{Float64,D},
        ctx; floor::Real=1.0) where {D}
    ws = ctx.preconditioner_workspace
    size(values) == ws.dims ||
        throw(DimensionMismatch("values dims $(size(values)) do not match workspace dims $(ws.dims)"))
    cutoff = Float64(floor)
    cutoff > 0.0 ||
        throw(ArgumentError("preconditioner floor must be positive"))
    @inbounds for idx in eachindex(values)
        ws.field_c[idx] = values[idx]
    end
    mul!(ws.field_hat, ws.fplan, ws.field_c)
    first = firstindex(ws.kernel_hat)
    @inbounds for idx in eachindex(ws.kernel_hat)
        if idx == first
            ws.kernel_hat[idx] = 0.0
        else
            denom = ctx.composition_kernel[idx] - 2.0 * ctx.chiN
            if abs(denom) < cutoff
                denom = denom < 0.0 ? -cutoff : cutoff
            end
            ws.kernel_hat[idx] = ws.field_hat[idx] / denom
        end
    end
    mul!(ws.kernel_time, ws.iplan, ws.kernel_hat)
    @inbounds for idx in eachindex(out)
        out[idx] = real(ws.kernel_time[idx])
    end
    out .-= mean(out)
    return out
end

function _diblock_burp_ti_apply_preconditioner!(
        out::AbstractArray{Float64,D}, values::AbstractArray{Float64,D},
        ctx; preconditioner::Symbol=:none, floor::Real=1.0) where {D}
    if preconditioner == :none
        out .= values
        out .-= mean(out)
        return out
    elseif preconditioner == :rpa_fourier
        return _diblock_burp_ti_rpa_preconditioner!(out, values, ctx;
            floor=floor)
    end
    throw(ArgumentError("unsupported BURP-TI GMRES preconditioner $(preconditioner)"))
end

"""
    diblock_burp_ti_projected_jacobian_vector_product_nd(phi, direction; ...)

Apply the exact Jacobian of the fixed-composition BURP-TI residual to a
zero-mean perturbation. This is the matrix-free operator used by Newton-GMRES.
"""
function diblock_burp_ti_projected_jacobian_vector_product_nd(
        phi_a::AbstractArray{<:Real}, direction_a::AbstractArray{<:Real};
        f::Real=0.5, chiN::Real=12.0, lengths=1.0, N::Real=1.0, b::Real=1.0,
        nquad::Integer=16, mean_tolerance::Real=1.0e-8,
        composition_kernel=nothing, workspace=nothing)
    phi, dims, ff, chi, cell_lengths, nn, bb, _ =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    size(direction_a) == dims ||
        throw(DimensionMismatch("direction dims $(size(direction_a)) do not match phi dims $(dims)"))
    direction = Float64.(Array(direction_a))
    direction .-= mean(direction)
    ctx = _diblock_burp_ti_context(dims, cell_lengths; f=ff, chiN=chi,
        N=nn, b=bb, nquad=nquad)
    if composition_kernel !== nothing
        size(composition_kernel) == dims ||
            throw(DimensionMismatch("precomputed composition kernel dims must match phi dims"))
        ctx = merge(ctx, (composition_kernel=composition_kernel,))
    end
    if workspace !== nothing
        workspace.dims == dims ||
            throw(DimensionMismatch("workspace dims $(workspace.dims) do not match phi dims $(dims)"))
        ctx = merge(ctx, (jv_workspace=workspace,))
    end
    out = zeros(Float64, dims...)
    return _diblock_burp_ti_projected_jv!(out, phi, direction, ctx)
end

Base.@kwdef struct DiblockBurpTiNewtonOptions
    force_tolerance::Float64 = 1.0e-4
    max_newton_iterations::Int = 14
    gmres_max_iterations::Int = 20
    gmres_restart::Int = 20
    gmres_tolerance::Float64 = 1.0e-2
    step_max::Float64 = 0.1
    density_floor::Float64 = 1.0e-6
    line_search_max_steps::Int = 20
    line_search_sufficient_decrease::Float64 = 1.0e-3
    nquad::Int = 64
    mean_tolerance::Float64 = 1.0e-8
    energy_density_tolerance::Float64 = 0.0
    energy_delta_force_tolerance::Float64 = Inf
    force_stagnation_iterations::Int = 0
    force_stagnation_relative_tolerance::Float64 = 0.0
    force_stagnation_energy_density_tolerance::Float64 = 0.0
    gmres_preconditioner::Symbol = :none
    gmres_preconditioner_floor::Float64 = 1.0
    progress_interval::Int = 0
    progress_callback::Any = nothing
    force_metric_callback::Any = nothing
    stationarity_mode::Symbol = :projected
    active_tolerance_factor::Float64 = 10.0
end

struct DiblockBurpTiNewtonPhasePolicy
    morphology::Symbol
    dimensionality::Int
    seed_strategy::String
    seed_amplitude::Float64
    max_newton_iterations::Int
    gmres_max_iterations::Int
    gmres_restart::Int
    gmres_tolerance::Float64
    step_max::Float64
    line_search_max_steps::Int
    line_search_sufficient_decrease::Float64
    gmres_preconditioner::Symbol
    gmres_preconditioner_floor::Float64
end

function diblock_burp_ti_newton_phase_policy(morphology)
    kind = _diblock_morphology_symbol(morphology)
    if kind == :LAM
        return DiblockBurpTiNewtonPhasePolicy(:LAM, 1,
            "warm_seed_or_lamellar_primitive", 0.75, 18, 12, 12, 1.0e-2,
            0.08, 24, 1.0e-3, :none, 1.0)
    elseif kind == :CYL
        return DiblockBurpTiNewtonPhasePolicy(:CYL, 2,
            "warm_seed_or_cylindrical_primitive", 0.65, 14, 20, 20, 1.0e-2,
            0.10, 20, 1.0e-3, :none, 1.0)
    elseif kind == :BCC
        return DiblockBurpTiNewtonPhasePolicy(:BCC, 3,
            "warm_seed_or_conservative_bcc_primitive", 0.50, 36, 24, 24,
            3.0e-3, 0.060, 32, 1.0e-3, :none, 1.0)
    elseif kind == :GYR
        return DiblockBurpTiNewtonPhasePolicy(:GYR, 3,
            "warm_seed_or_conservative_gyroid_primitive", 0.35, 48, 96, 64,
            1.0e-4, 0.080, 36, 1.0e-3, :none, 1.0)
    end
    throw(ArgumentError("unsupported BURP-TI Newton-GMRES phase $(morphology)"))
end

_burp_ti_policy_override(value, default) = value === nothing ? default : value

"""
    diblock_burp_ti_newton_options_for_phase(morphology; ...)

Return production-oriented Newton-GMRES options for a morphology.  CYL keeps
the benchmark settings.  BCC/GYR use larger Krylov subspaces because production
3D grids otherwise make each Newton step too inexact, then control the basin
with step caps and line search.
"""
function diblock_burp_ti_newton_options_for_phase(morphology;
        force_tolerance::Real=1.0e-4, nquad::Integer=64,
        density_floor::Real=1.0e-6, mean_tolerance::Real=1.0e-8,
        max_newton_iterations=nothing, gmres_max_iterations=nothing,
        gmres_restart=nothing, gmres_tolerance=nothing, step_max=nothing,
        gmres_preconditioner=nothing, gmres_preconditioner_floor=nothing,
        line_search_max_steps=nothing,
        line_search_sufficient_decrease=nothing, progress_interval::Integer=0,
        progress_callback=nothing, energy_density_tolerance::Real=0.0,
        energy_delta_force_tolerance::Real=Inf,
        force_stagnation_iterations::Integer=0,
        force_stagnation_relative_tolerance::Real=0.0,
        force_stagnation_energy_density_tolerance::Real=0.0,
        force_metric_callback=nothing,
        stationarity_mode::Symbol=:projected,
        active_tolerance_factor::Real=10.0)
    policy = diblock_burp_ti_newton_phase_policy(morphology)
    gmres_max = Int(_burp_ti_policy_override(gmres_max_iterations,
        policy.gmres_max_iterations))
    restart_default = gmres_max_iterations === nothing ?
        policy.gmres_restart : gmres_max
    return DiblockBurpTiNewtonOptions(;
        force_tolerance=Float64(force_tolerance),
        max_newton_iterations=Int(_burp_ti_policy_override(
            max_newton_iterations, policy.max_newton_iterations)),
        gmres_max_iterations=gmres_max,
        gmres_restart=Int(_burp_ti_policy_override(gmres_restart,
            restart_default)),
        gmres_tolerance=Float64(_burp_ti_policy_override(gmres_tolerance,
            policy.gmres_tolerance)),
        step_max=Float64(_burp_ti_policy_override(step_max,
            policy.step_max)),
        density_floor=Float64(density_floor),
        line_search_max_steps=Int(_burp_ti_policy_override(
            line_search_max_steps, policy.line_search_max_steps)),
        line_search_sufficient_decrease=Float64(_burp_ti_policy_override(
            line_search_sufficient_decrease,
            policy.line_search_sufficient_decrease)),
        nquad=Int(nquad),
        mean_tolerance=Float64(mean_tolerance),
        energy_density_tolerance=Float64(energy_density_tolerance),
        energy_delta_force_tolerance=Float64(energy_delta_force_tolerance),
        force_stagnation_iterations=Int(force_stagnation_iterations),
        force_stagnation_relative_tolerance=Float64(
            force_stagnation_relative_tolerance),
        force_stagnation_energy_density_tolerance=Float64(
            force_stagnation_energy_density_tolerance),
        gmres_preconditioner=Symbol(_burp_ti_policy_override(
            gmres_preconditioner, policy.gmres_preconditioner)),
        gmres_preconditioner_floor=Float64(_burp_ti_policy_override(
            gmres_preconditioner_floor, policy.gmres_preconditioner_floor)),
        progress_interval=Int(progress_interval),
        progress_callback=progress_callback,
        force_metric_callback=force_metric_callback,
        stationarity_mode=Symbol(stationarity_mode),
        active_tolerance_factor=Float64(active_tolerance_factor))
end

"""
    diblock_burp_ti_prepare_newton_seed(morphology; f, dims, lengths, seed, ...)

Prepare a fixed-mean, box-bounded density for Newton-GMRES.  If `seed` is
provided, it is treated as the preferred warm/in-model seed and only projected
to the requested composition and density box.  Otherwise a conservative
phase-specific primitive seed is generated.
"""
function diblock_burp_ti_prepare_newton_seed(morphology; f::Real,
        dims=nothing, lengths=1.0, seed=nothing, amplitude=nothing,
        density_floor::Real=1.0e-6)
    ff = Float64(f)
    0.0 < ff < 1.0 ||
        throw(ArgumentError("diblock fraction f must be in (0, 1)"))
    floor = Float64(density_floor)
    0.0 <= floor < min(ff, 1.0 - ff) ||
        throw(ArgumentError("density_floor must be nonnegative and below both f and 1-f"))
    if seed !== nothing
        raw = Float64.(Array(seed))
        if dims !== nothing
            target_dims = _diblock_morphology_dims(dims)
            size(raw) == target_dims ||
                throw(DimensionMismatch("seed dims $(size(raw)) do not match requested dims $(target_dims)"))
        end
        return _project_box_mean_density(raw, ff; lo=floor, hi=1.0 - floor)
    end

    dims === nothing &&
        throw(ArgumentError("dims is required when no warm seed is supplied"))
    policy = diblock_burp_ti_newton_phase_policy(morphology)
    amp = Float64(_burp_ti_policy_override(amplitude, policy.seed_amplitude))
    raw = diblock_morphology_seed(morphology; f=ff, dims=dims,
        lengths=lengths, amplitude=amp)
    return _project_box_mean_density(raw, ff; lo=floor, hi=1.0 - floor)
end

struct DiblockBurpTiNewtonResult{D}
    phi::Array{Float64,D}
    initial_energy::Float64
    final_energy::Float64
    initial_projected_force_norm::Float64
    initial_projected_force_maxabs::Float64
    final_projected_force_norm::Float64
    final_projected_force_maxabs::Float64
    converged::Bool
    stop_reason::String
    newton_iterations::Int
    gmres_iterations::Int
    jv_calls::Int
    residual_calls::Int
    energy_calls::Int
    gmres_final_relative_residual::Float64
    energy_history::Vector{Float64}
    final_neighbor_energy_delta::Float64
    final_neighbor_energy_density_delta::Float64
    density_floor::Float64
    nquad::Int
    lengths::Vector{Float64}
end

function _validate_burp_ti_newton_options(options::DiblockBurpTiNewtonOptions,
        f::Float64)
    isfinite(options.force_tolerance) && options.force_tolerance > 0.0 ||
        throw(ArgumentError("force_tolerance must be positive and finite"))
    options.max_newton_iterations >= 0 ||
        throw(ArgumentError("max_newton_iterations must be nonnegative"))
    options.gmres_max_iterations >= 1 ||
        throw(ArgumentError("gmres_max_iterations must be positive"))
    options.gmres_restart >= 1 ||
        throw(ArgumentError("gmres_restart must be positive"))
    isfinite(options.gmres_tolerance) && options.gmres_tolerance > 0.0 ||
        throw(ArgumentError("gmres_tolerance must be positive and finite"))
    options.gmres_preconditioner in (:none, :rpa_fourier) ||
        throw(ArgumentError("unsupported gmres_preconditioner $(options.gmres_preconditioner)"))
    isfinite(options.gmres_preconditioner_floor) &&
        options.gmres_preconditioner_floor > 0.0 ||
        throw(ArgumentError("gmres_preconditioner_floor must be positive and finite"))
    isfinite(options.step_max) && options.step_max > 0.0 ||
        throw(ArgumentError("step_max must be positive and finite"))
    0.0 <= options.density_floor < min(f, 1.0 - f) ||
        throw(ArgumentError("density_floor must be nonnegative and below both f and 1-f"))
    options.line_search_max_steps >= 1 ||
        throw(ArgumentError("line_search_max_steps must be positive"))
    isfinite(options.line_search_sufficient_decrease) &&
        options.line_search_sufficient_decrease >= 0.0 ||
        throw(ArgumentError("line_search_sufficient_decrease must be nonnegative and finite"))
    options.nquad >= 2 || throw(ArgumentError("nquad must be at least 2"))
    isfinite(options.mean_tolerance) && options.mean_tolerance > 0.0 ||
        throw(ArgumentError("mean_tolerance must be positive and finite"))
    isfinite(options.energy_density_tolerance) &&
        options.energy_density_tolerance >= 0.0 ||
        throw(ArgumentError("energy_density_tolerance must be nonnegative and finite"))
    options.energy_delta_force_tolerance > 0.0 ||
        throw(ArgumentError("energy_delta_force_tolerance must be positive"))
    options.force_stagnation_iterations >= 0 ||
        throw(ArgumentError("force_stagnation_iterations must be nonnegative"))
    isfinite(options.force_stagnation_relative_tolerance) &&
        options.force_stagnation_relative_tolerance >= 0.0 ||
        throw(ArgumentError("force_stagnation_relative_tolerance must be nonnegative and finite"))
    isfinite(options.force_stagnation_energy_density_tolerance) &&
        options.force_stagnation_energy_density_tolerance >= 0.0 ||
        throw(ArgumentError("force_stagnation_energy_density_tolerance must be nonnegative and finite"))
    options.progress_interval >= 0 ||
        throw(ArgumentError("progress_interval must be nonnegative"))
    options.stationarity_mode in (:projected, :bounded_kkt) ||
        throw(ArgumentError("unsupported stationarity_mode $(options.stationarity_mode)"))
    isfinite(options.active_tolerance_factor) &&
        options.active_tolerance_factor >= 1.0 ||
        throw(ArgumentError("active_tolerance_factor must be finite and at least one"))
    return nothing
end

function _diblock_burp_ti_newton_progress(options::DiblockBurpTiNewtonOptions,
        event)
    options.progress_callback === nothing && return nothing
    options.progress_callback(event)
    return nothing
end

function _diblock_burp_ti_newton_force_metric(
        options::DiblockBurpTiNewtonOptions, phi, force; f::Float64,
        chiN::Float64, lengths, stationarity_force=force)
    options.force_metric_callback === nothing &&
        return (norm=Float64(stationarity_force.norm),
            maxabs=Float64(stationarity_force.maxabs),
            name=options.stationarity_mode == :bounded_kkt ?
                :bounded_kkt_force : :projected_force,
            raw=stationarity_force)
    metric = options.force_metric_callback((phi=phi, force=force, f=f,
        stationarity_force=stationarity_force, chiN=chiN, lengths=lengths,
        density_floor=options.density_floor, nquad=options.nquad))
    if metric isa Real
        return (norm=Float64(metric), maxabs=NaN, name=:custom, raw=metric)
    end
    norm_value = Float64(getproperty(metric, :norm))
    maxabs_value = hasproperty(metric, :maxabs) ?
        Float64(getproperty(metric, :maxabs)) : NaN
    metric_name = hasproperty(metric, :name) ? getproperty(metric, :name) :
        :custom
    return merge((norm=norm_value, maxabs=maxabs_value, name=metric_name),
        (raw=metric,))
end

function _diblock_burp_ti_newton_stationarity(force, phi,
        options::DiblockBurpTiNewtonOptions, dV::Float64)
    options.stationarity_mode == :projected &&
        return (force=force, movable=nothing)
    active_tolerance = options.active_tolerance_factor * options.density_floor
    bounded = _diblock_bounded_projected_force_result(
        force.chemical_potential, phi, dV;
        lower=options.density_floor, upper=1.0 - options.density_floor,
        active_tolerance=active_tolerance)
    movable = _diblock_burp_ti_bounded_newton_movable_mask(bounded, phi;
        lower=options.density_floor, upper=1.0 - options.density_floor,
        active_tolerance=active_tolerance)
    return (force=bounded, movable=movable)
end

function _project_box_mean_density(values::AbstractArray{<:Real}, f::Float64;
        lo::Float64, hi::Float64)
    lo <= f <= hi ||
        throw(ArgumentError("target mean f must lie inside the projection box"))
    input = Float64.(Array(values))
    all(isfinite, input) ||
        throw(ArgumentError("density field must contain only finite values"))
    target = f * length(input)
    low = lo - maximum(input)
    high = hi - minimum(input)
    for _ in 1:80
        mid = 0.5 * (low + high)
        total = 0.0
        @inbounds for value in input
            total += clamp(value + mid, lo, hi)
        end
        if total > target
            high = mid
        else
            low = mid
        end
    end
    shift = 0.5 * (low + high)
    out = similar(input)
    @inbounds for idx in eachindex(input)
        out[idx] = clamp(input[idx] + shift, lo, hi)
    end
    return out
end

function _project_box_mean_density_bounded_step(
        values::AbstractArray{<:Real}, current::AbstractArray{<:Real},
        movable::AbstractArray{Bool}, f::Float64; lo::Float64, hi::Float64)
    size(values) == size(current) == size(movable) ||
        throw(DimensionMismatch("values, current field, and movable mask dimensions must match"))
    movable_count = count(movable)
    movable_count > 0 || return Float64.(Array(current))
    fixed_sum = 0.0
    movable_min = Inf
    movable_max = -Inf
    @inbounds for idx in eachindex(values)
        if movable[idx]
            value = Float64(values[idx])
            movable_min = min(movable_min, value)
            movable_max = max(movable_max, value)
        else
            fixed_sum += Float64(current[idx])
        end
    end
    movable_target = f * length(values) - fixed_sum
    lo * movable_count - 100eps(Float64) <= movable_target <=
        hi * movable_count + 100eps(Float64) ||
        throw(ArgumentError("fixed active bounds leave an infeasible composition target"))
    low = lo - movable_max
    high = hi - movable_min
    for _ in 1:80
        mid = 0.5 * (low + high)
        total = 0.0
        @inbounds for idx in eachindex(values)
            movable[idx] || continue
            total += clamp(Float64(values[idx]) + mid, lo, hi)
        end
        if total > movable_target
            high = mid
        else
            low = mid
        end
    end
    shift = 0.5 * (low + high)
    out = Float64.(Array(current))
    @inbounds for idx in eachindex(out)
        movable[idx] || continue
        out[idx] = clamp(Float64(values[idx]) + shift, lo, hi)
    end
    return out
end

function _diblock_burp_ti_gmres(phi::Array{Float64,D},
        residual::Array{Float64,D}, ctx; maxiter::Integer, restart::Integer,
        tolerance::Real, preconditioner::Symbol=:none,
        preconditioner_floor::Real=1.0, counters,
        movable_mask=nothing) where {D}
    rhs = -vec(residual)
    n = length(rhs)
    solution = zeros(Float64, n)
    beta0 = norm(rhs)
    beta0 == 0.0 && return solution, 0, 0.0
    preconditioned = zeros(Float64, ctx.dims...)

    function apply_preconditioner(vec_direction)
        direction = reshape(vec_direction, ctx.dims)
        _diblock_burp_ti_apply_preconditioner!(preconditioned, direction, ctx;
            preconditioner=preconditioner, floor=preconditioner_floor)
        movable_mask === nothing ||
            _diblock_burp_ti_project_bounded_tangent!(preconditioned,
                movable_mask)
        return preconditioned
    end

    function apply_jacobian(vec_direction)
        direction = apply_preconditioner(vec_direction)
        _diblock_burp_ti_projected_jv!(ctx.jout, phi, direction, ctx)
        movable_mask === nothing ||
            _diblock_burp_ti_project_bounded_tangent!(ctx.jout,
                movable_mask)
        counters.jv_calls[] += 1
        return vec(copy(ctx.jout))
    end

    current_residual = copy(rhs)
    total_iterations = 0
    final_relres = Inf
    while total_iterations < maxiter
        beta = norm(current_residual)
        final_relres = beta / beta0
        final_relres < tolerance &&
            return solution, total_iterations, final_relres
        m = min(Int(restart), Int(maxiter) - total_iterations)
        basis = zeros(Float64, n, m + 1)
        hessenberg = zeros(Float64, m + 1, m)
        basis[:, 1] .= current_residual ./ beta
        target = zeros(Float64, m + 1)
        target[1] = beta
        best_y = zeros(Float64, 0)
        best_j = 0
        best_residual = Inf
        for j in 1:m
            work = apply_jacobian(basis[:, j])
            for i in 1:j
                hessenberg[i, j] = dot(basis[:, i], work)
                @. work -= hessenberg[i, j] * basis[:, i]
            end
            hessenberg[j + 1, j] = norm(work)
            if hessenberg[j + 1, j] > 1.0e-13 && j < m
                basis[:, j + 1] .= work ./ hessenberg[j + 1, j]
            end
            y = hessenberg[1:j + 1, 1:j] \ target[1:j + 1]
            relres = norm(target[1:j + 1] -
                hessenberg[1:j + 1, 1:j] * y) / beta0
            if relres < best_residual
                best_residual = relres
                best_y = y
                best_j = j
            end
            if relres < tolerance
                solution .+= basis[:, 1:j] * y
                delta = vec(copy(apply_preconditioner(solution)))
                return delta, total_iterations + j, relres
            end
        end
        if best_j == 0
            delta = vec(copy(apply_preconditioner(solution)))
            return delta, total_iterations, final_relres
        end
        solution .+= basis[:, 1:best_j] * best_y
        total_iterations += best_j
        current_residual = rhs - apply_jacobian(solution)
    end
    final_relres = norm(rhs - apply_jacobian(solution)) / beta0
    delta = vec(copy(apply_preconditioner(solution)))
    return delta, total_iterations, final_relres
end

"""
    solve_diblock_burp_ti_newton_gmres(initial_phi; f, chiN, lengths, options, ...)

Solve the fixed-composition BURP-TI stationarity equations with Newton steps
and a matrix-free exact-Jacobian GMRES linear solve. The solver is
phase-agnostic: it accepts any periodic density array and matching cell
lengths. Set `options.stationarity_mode=:bounded_kkt` to solve in the movable
tangent space of the current box-constrained KKT active set.
"""
function solve_diblock_burp_ti_newton_gmres(initial_phi::AbstractArray{<:Real};
        f::Real=0.5, chiN::Real=12.0, lengths=1.0, N::Real=1.0, b::Real=1.0,
        options::DiblockBurpTiNewtonOptions=DiblockBurpTiNewtonOptions(),
        composition_kernel=nothing)
    raw = Float64.(Array(initial_phi))
    dims = size(raw)
    length(raw) >= 8 ||
        throw(ArgumentError("BURP-TI Newton-GMRES requires at least eight grid points"))
    all(count -> count >= 2, dims) ||
        throw(ArgumentError("BURP-TI Newton-GMRES requires at least two points per dimension"))
    all(isfinite, raw) ||
        throw(ArgumentError("initial density field must contain only finite values"))
    ff = Float64(f)
    0.0 < ff < 1.0 || throw(ArgumentError("diblock fraction f must be in (0, 1)"))
    chi = Float64(chiN)
    isfinite(chi) || throw(ArgumentError("chiN must be finite"))
    _, nn, bb, _ = _diblock_validate(ff, N, b, 1.0)
    _validate_burp_ti_newton_options(options, ff)
    cell_lengths = _diblock_nd_lengths(lengths, dims)

    phi = _project_box_mean_density(raw, ff; lo=options.density_floor,
        hi=1.0 - options.density_floor)
    ctx = _diblock_burp_ti_context(dims, cell_lengths; f=ff, chiN=chi,
        N=nn, b=bb, nquad=options.nquad)
    if composition_kernel !== nothing
        size(composition_kernel) == dims ||
            throw(DimensionMismatch("precomputed composition kernel dims must match phi dims"))
        ctx = merge(ctx, (composition_kernel=composition_kernel,))
    end

    counters = (jv_calls=Ref(0), residual_calls=Ref(0), energy_calls=Ref(0))
    function residual_with_count(field)
        counters.residual_calls[] += 1
        chemical = diblock_burp_ti_chemical_potential_nd(field; f=ff,
            chiN=chi, lengths=cell_lengths, N=nn, b=bb, nquad=options.nquad,
            mean_tolerance=options.mean_tolerance,
            composition_kernel=ctx.composition_kernel,
            workspace=ctx.force_workspace)
        return _diblock_projected_force_result(chemical, ctx.dV)
    end
    function energy_with_count(field)
        counters.energy_calls[] += 1
        return diblock_burp_ti_energy_nd(field; f=ff, chiN=chi,
            lengths=cell_lengths, N=nn, b=bb, nquad=options.nquad,
            mean_tolerance=options.mean_tolerance,
            composition_kernel=ctx.composition_kernel,
            workspace=ctx.energy_workspace)
    end

    force = residual_with_count(phi)
    stationarity = _diblock_burp_ti_newton_stationarity(force, phi, options,
        ctx.dV)
    force_metric = _diblock_burp_ti_newton_force_metric(options, phi, force;
        f=ff, chiN=chi, lengths=cell_lengths,
        stationarity_force=stationarity.force)
    initial_force_norm = force.norm
    initial_force_max = force.maxabs
    initial_energy = energy_with_count(phi)
    cell_volume = prod(cell_lengths)
    energy_history = Float64[Float64(initial_energy)]
    last_energy = Float64(initial_energy)
    final_neighbor_energy_delta = Inf
    final_neighbor_energy_density_delta = Inf
    stop_reason = force_metric.norm <= options.force_tolerance ?
        "force_tol" : "iteration_cap"
    newton_iterations = 0
    gmres_iterations = 0
    gmres_final_relative_residual = NaN
    stagnation_reference_iter = 0
    stagnation_reference_force = Float64(force_metric.norm)
    solve_started = time()
    _diblock_burp_ti_newton_progress(options, (stage=:initial,
        iteration=0, force_norm=Float64(force_metric.norm),
        force_maxabs=Float64(force_metric.maxabs),
        raw_force_norm=Float64(force.norm),
        raw_force_maxabs=Float64(force.maxabs),
        force_metric=String(force_metric.name), gmres_iterations=0,
        gmres_total_iterations=0, gmres_final_relative_residual=NaN,
        alpha=0.0, step_max=Float64(options.step_max),
        energy=Float64(initial_energy), energy_delta=NaN,
        energy_density_delta=NaN,
        stop_reason=String(stop_reason), elapsed_s=0.0,
        jv_calls=0, residual_calls=Int(counters.residual_calls[]),
        energy_calls=Int(counters.energy_calls[])))

    for newton_iter in 1:options.max_newton_iterations
        force_metric.norm <= options.force_tolerance && break
        delta, inner_iterations, inner_residual = _diblock_burp_ti_gmres(phi,
            stationarity.force.residual, ctx;
            maxiter=options.gmres_max_iterations,
            restart=options.gmres_restart, tolerance=options.gmres_tolerance,
            preconditioner=options.gmres_preconditioner,
            preconditioner_floor=options.gmres_preconditioner_floor,
            counters=counters, movable_mask=stationarity.movable)
        newton_iterations = newton_iter
        gmres_iterations += inner_iterations
        gmres_final_relative_residual = inner_residual

        direction = reshape(delta, dims)
        if stationarity.movable === nothing
            direction .-= mean(direction)
        else
            _diblock_burp_ti_project_bounded_tangent!(direction,
                stationarity.movable)
        end
        max_delta = maximum(abs, direction)
        scaled_max_delta = max_delta
        if !isfinite(max_delta) || max_delta == 0.0
            stop_reason = "newton_direction_zero"
            break
        elseif max_delta > options.step_max
            direction .*= options.step_max / max_delta
            scaled_max_delta = options.step_max
        end

        best_phi = phi
        best_force = force
        best_stationarity = stationarity
        best_force_metric = force_metric
        best_alpha = 0.0
        alpha = 1.0
        for _ in 1:options.line_search_max_steps
            proposed = phi .+ alpha .* direction
            candidate = if stationarity.movable === nothing
                _project_box_mean_density(proposed, ff;
                    lo=options.density_floor,
                    hi=1.0 - options.density_floor)
            else
                _project_box_mean_density_bounded_step(proposed, phi,
                    stationarity.movable, ff; lo=options.density_floor,
                    hi=1.0 - options.density_floor)
            end
            trial_force = residual_with_count(candidate)
            trial_stationarity = _diblock_burp_ti_newton_stationarity(
                trial_force, candidate, options, ctx.dV)
            trial_force_metric = _diblock_burp_ti_newton_force_metric(options,
                candidate, trial_force; f=ff, chiN=chi,
                lengths=cell_lengths,
                stationarity_force=trial_stationarity.force)
            if isfinite(trial_force_metric.norm) &&
                    trial_force_metric.norm < best_force_metric.norm
                best_phi = candidate
                best_force = trial_force
                best_stationarity = trial_stationarity
                best_force_metric = trial_force_metric
                best_alpha = alpha
            end
            if trial_force_metric.norm <=
                    (1.0 - options.line_search_sufficient_decrease * alpha) *
                    force_metric.norm
                break
            end
            alpha *= 0.5
        end

        phi = best_phi
        force = best_force
        stationarity = best_stationarity
        force_metric = best_force_metric
        current_energy = energy_with_count(phi)
        push!(energy_history, Float64(current_energy))
        final_neighbor_energy_delta = abs(Float64(current_energy) - last_energy)
        final_neighbor_energy_density_delta =
            final_neighbor_energy_delta / cell_volume
        last_energy = Float64(current_energy)
        if force_metric.norm <= options.force_tolerance
            stop_reason = "force_tol"
            break
        elseif best_alpha == 0.0
            stop_reason = "line_search_stalled"
            break
        elseif options.energy_density_tolerance > 0.0 &&
                force_metric.norm <= options.energy_delta_force_tolerance &&
                final_neighbor_energy_density_delta <=
                options.energy_density_tolerance
            stop_reason = "energy_delta_tol"
            break
        elseif options.force_stagnation_iterations > 0 &&
                newton_iter - stagnation_reference_iter >=
                options.force_stagnation_iterations
            relative_improvement =
                (stagnation_reference_force - force_metric.norm) /
                max(stagnation_reference_force, eps(Float64))
            if relative_improvement <=
                    options.force_stagnation_relative_tolerance &&
                    final_neighbor_energy_density_delta <=
                    options.force_stagnation_energy_density_tolerance
                stop_reason = "force_stagnation"
                break
            end
            stagnation_reference_iter = newton_iter
            stagnation_reference_force = Float64(force_metric.norm)
        else
            stop_reason = "iteration_cap"
        end
        if options.progress_interval > 0 &&
                (newton_iter % options.progress_interval == 0 ||
                 force_metric.norm <= options.force_tolerance ||
                 best_alpha == 0.0 ||
                 newton_iter == options.max_newton_iterations)
            _diblock_burp_ti_newton_progress(options, (stage=:iteration,
                iteration=Int(newton_iter),
                force_norm=Float64(force_metric.norm),
                force_maxabs=Float64(force_metric.maxabs),
                raw_force_norm=Float64(force.norm),
                raw_force_maxabs=Float64(force.maxabs),
                force_metric=String(force_metric.name),
                gmres_iterations=Int(inner_iterations),
                gmres_total_iterations=Int(gmres_iterations),
                gmres_final_relative_residual=Float64(inner_residual),
                alpha=Float64(best_alpha),
                step_max=Float64(scaled_max_delta),
                energy=Float64(current_energy),
                energy_delta=Float64(final_neighbor_energy_delta),
                energy_density_delta=Float64(final_neighbor_energy_density_delta),
                stop_reason=String(stop_reason),
                elapsed_s=Float64(time() - solve_started),
                jv_calls=Int(counters.jv_calls[]),
                residual_calls=Int(counters.residual_calls[]),
                energy_calls=Int(counters.energy_calls[])))
        end
    end

    final_energy = last_energy
    _diblock_burp_ti_newton_progress(options, (stage=:final,
        iteration=Int(newton_iterations), force_norm=Float64(force_metric.norm),
        force_maxabs=Float64(force_metric.maxabs),
        raw_force_norm=Float64(force.norm),
        raw_force_maxabs=Float64(force.maxabs),
        force_metric=String(force_metric.name),
        gmres_iterations=0, gmres_total_iterations=Int(gmres_iterations),
        gmres_final_relative_residual=Float64(gmres_final_relative_residual),
        alpha=0.0, step_max=Float64(options.step_max),
        energy=Float64(final_energy),
        energy_delta=Float64(final_neighbor_energy_delta),
        energy_density_delta=Float64(final_neighbor_energy_density_delta),
        stop_reason=String(stop_reason),
        elapsed_s=Float64(time() - solve_started),
        jv_calls=Int(counters.jv_calls[]),
        residual_calls=Int(counters.residual_calls[]),
        energy_calls=Int(counters.energy_calls[])))
    energy_delta_converged = options.energy_density_tolerance > 0.0 &&
        force_metric.norm <= options.energy_delta_force_tolerance &&
        final_neighbor_energy_density_delta <= options.energy_density_tolerance
    converged = force_metric.norm <= options.force_tolerance ||
        energy_delta_converged
    return DiblockBurpTiNewtonResult{length(dims)}(phi, Float64(initial_energy),
        Float64(final_energy), Float64(initial_force_norm),
        Float64(initial_force_max), Float64(force.norm), Float64(force.maxabs),
        Bool(converged), String(stop_reason),
        Int(newton_iterations), Int(gmres_iterations), Int(counters.jv_calls[]),
        Int(counters.residual_calls[]), Int(counters.energy_calls[]),
        Float64(gmres_final_relative_residual), Float64.(energy_history),
        Float64(final_neighbor_energy_delta),
        Float64(final_neighbor_energy_density_delta),
        Float64(options.density_floor),
        Int(options.nquad), Float64.(collect(cell_lengths)))
end
