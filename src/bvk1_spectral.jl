const DIBLOCK_BVK1_SPECTRAL_DISCRETIZATION_SCHEMA =
    "bvk1-spectral-continuum-v1"

const DIBLOCK_BVK1_DEALIASED_OBJECTIVE_SCHEMA =
    "bvk1-spectral-theta-phase-averaged-v3"

const _BVK1_GOLDEN_TRANSLATION = (sqrt(5.0)-1.0)/2.0

@inline function _bvk1_fft_mode(index::Int, count::Int)
    zero_based = index - 1
    return zero_based <= count ÷ 2 ? zero_based : zero_based - count
end

@inline _bvk1_fft_index(mode::Int, count::Int) =
    mode >= 0 ? mode + 1 : count + mode + 1

@inline function _bvk1_fft_mode_representable(mode::Int, count::Int)
    half = count ÷ 2
    return iseven(count) ? abs(mode) < half : abs(mode) <= half
end

function _bvk1_spectral_shift(values::AbstractArray{<:Real,D},
        shifts::NTuple{D,Float64}) where {D}
    field = Float64.(values)
    dims = size(field)
    transformed = fft(field)
    for index in CartesianIndices(field)
        modes = ntuple(axis -> _bvk1_fft_mode(index[axis], dims[axis]), D)
        if any(axis -> iseven(dims[axis]) &&
                abs(modes[axis]) == dims[axis] ÷ 2, 1:D)
            transformed[index] = 0.0
            continue
        end
        phase = cis(-2pi * sum(modes[axis] * shifts[axis] /
            dims[axis] for axis in 1:D))
        transformed[index] *= phase
    end
    return real.(ifft(transformed))
end

function _bvk1_translation_choices(rank::Integer, count::Integer)
    phases = Int(count)
    phases >= 1 || throw(ArgumentError(
        "translation_phases must be positive"))
    offsets = phases == 1 ? (0.0,) : Tuple(mod(
        (index + _BVK1_GOLDEN_TRANSLATION) / phases, 1.0)
        for index in 0:phases-1)
    return collect(Iterators.product(ntuple(_ -> offsets, Int(rank))...))
end

"""Cell-centred Fourier interpolation used by the dealiased objective.

Even-grid Nyquist hyperplanes are removed.  Their sign is ambiguous under a
fractional translation and a well-resolved ordered field must not carry
physical information there.  Removing them makes the interpolation map and
its adjoint exact, while also preventing the optimizer from exploiting a
grid-scale checkerboard mode.
"""
function _bvk1_cell_centered_resample(values::AbstractArray{<:Real,D},
        target_dims::NTuple{D,Int}) where {D}
    source = Float64.(values)
    all(>(0), target_dims) || throw(ArgumentError(
        "target grid dimensions must be positive"))
    size(source) == target_dims && return copy(source)
    source_dims = size(source)
    source_count = prod(source_dims)
    target_count = prod(target_dims)
    transformed = fft(source)
    padded = zeros(ComplexF64, target_dims)
    ratio = target_count / source_count
    for source_index in CartesianIndices(source)
        modes = ntuple(axis -> _bvk1_fft_mode(
            source_index[axis], source_dims[axis]), D)
        any(axis -> iseven(source_dims[axis]) &&
            abs(modes[axis]) == source_dims[axis] ÷ 2, 1:D) && continue
        any(axis -> !_bvk1_fft_mode_representable(
            modes[axis], target_dims[axis]), 1:D) && continue
        target_index = CartesianIndex(ntuple(axis ->
            _bvk1_fft_index(modes[axis], target_dims[axis]), D))
        phase = cis(pi * sum(modes[axis] *
            (inv(target_dims[axis]) - inv(source_dims[axis])) for axis in 1:D))
        padded[target_index] = ratio * phase * transformed[source_index]
    end
    return real.(ifft(padded))
end

"""Euclidean adjoint of `_bvk1_cell_centered_resample`.

The reverse pass uses the same retained-mode set and conjugate cell-centre
phase as the forward interpolation.  This is deliberately not implemented as
a reverse resampling call because even-grid Nyquist splitting would then make
the two maps only approximately adjoint.
"""
function _bvk1_cell_centered_resample_adjoint(
        values::AbstractArray{<:Real,D}, source_dims::NTuple{D,Int}) where {D}
    target_dims = size(values)
    all(>(0), source_dims) || throw(ArgumentError(
        "source grid dimensions must be positive"))
    target_dims == source_dims && return Float64.(values)
    source_count = prod(source_dims)
    target_count = prod(target_dims)
    target_adjoint = fft(Float64.(values)) ./ target_count
    source_adjoint = zeros(ComplexF64, source_dims)
    ratio = target_count / source_count
    for source_index in CartesianIndices(source_adjoint)
        modes = ntuple(axis -> _bvk1_fft_mode(
            source_index[axis], source_dims[axis]), D)
        any(axis -> iseven(source_dims[axis]) &&
            abs(modes[axis]) == source_dims[axis] ÷ 2, 1:D) && continue
        any(axis -> !_bvk1_fft_mode_representable(
            modes[axis], target_dims[axis]), 1:D) && continue
        target_index = CartesianIndex(ntuple(axis ->
            _bvk1_fft_index(modes[axis], target_dims[axis]), D))
        phase = cis(pi * sum(modes[axis] *
            (inv(target_dims[axis]) - inv(source_dims[axis])) for axis in 1:D))
        source_adjoint[source_index] = ratio * conj(phase) *
            target_adjoint[target_index]
    end
    return real.(source_count .* ifft(source_adjoint))
end

function _bvk1_theta_uniform_mean_retract!(
        theta::AbstractArray{Float64}, f::Real;
        tolerance::Real=1.0e-13, max_iterations::Integer=80)
    target = Float64(f)
    0.0 < target < 1.0 || throw(ArgumentError("f must be in (0, 1)"))
    inv_count = inv(Float64(length(theta)))
    for _ in 1:Int(max_iterations)
        density_sum = 0.0
        slope_sum = 0.0
        @inbounds @simd for index in eachindex(theta)
            sine = sin(theta[index])
            sine2 = sin(2theta[index])
            density_sum += sine * sine
            slope_sum += sine2
        end
        residual = density_sum * inv_count - target
        abs(residual) <= Float64(tolerance) && return theta
        slope = slope_sum * inv_count
        abs(slope) > 1.0e-14 || throw(ArgumentError(
            "theta mean is frozen because the field has no interface cells"))
        correction = clamp(-residual / slope, -0.25, 0.25)
        @inbounds @simd for index in eachindex(theta)
            theta[index] += correction
        end
    end
    throw(ArgumentError("theta mean retraction did not converge"))
end

"""Signed Fourier wavenumbers for a real periodic collocation grid.

The even-grid Nyquist coefficient is set to zero for first derivatives.  Its
sign is ambiguous for a real collocation field, and zero makes the discrete
first-derivative operator exactly skew-adjoint.  The Laplacian continues to
use the full squared wavenumber supplied by `_periodic_k2_nd`.
"""
function _bvk1_spectral_wavenumbers_1d(n::Integer, L::Real)
    count = Int(n)
    count >= 2 || throw(ArgumentError("periodic grid requires at least two points"))
    length_value = Float64(L)
    length_value > 0.0 || throw(ArgumentError("periodic length must be positive"))
    half = count ÷ 2
    values = Vector{Float64}(undef, count)
    for index in 1:count
        mode = index <= half + 1 ? index - 1 : index - 1 - count
        if iseven(count) && mode == half
            mode = 0
        end
        values[index] = 2.0pi * mode / length_value
    end
    return values
end

function _bvk1_spectral_derivative(values::AbstractArray{Float64,D},
        lengths::AbstractVector{Float64}, axis::Integer) where {D}
    1 <= Int(axis) <= D || throw(ArgumentError("spectral derivative axis out of range"))
    length(lengths) == D || throw(ArgumentError(
        "lengths must have one entry per array dimension"))
    wave = _bvk1_spectral_wavenumbers_1d(size(values, axis), lengths[axis])
    transformed = fft(values)
    differentiated = similar(transformed)
    for index in CartesianIndices(values)
        differentiated[index] = im * wave[index[axis]] * transformed[index]
    end
    return real.(ifft(differentiated))
end

# The nonlinear solver calls the spectral objective sequentially within one
# Julia process; campaign parallelism uses separate processes.  Cache the two
# complex FFT buffers and their in-place plans at each grid size, matching the
# established BVK2 workspace contract.  Physical-space outputs are always
# copied into fresh real arrays before returning, so no caller can alias this
# scratch between objective evaluations.
mutable struct _BVK1SpectralFFTWorkspace{D}
    lengths::NTuple{D,Float64}
    waves::NTuple{D,Vector{Float64}}
    k2::Array{Float64,D}
    source_hat::Array{ComplexF64,D}
    work_hat::Array{ComplexF64,D}
    source_forward::Any
    source_inverse::Any
    work_forward::Any
    work_inverse::Any
end

function _bvk1_spectral_fft_workspace(dims::NTuple{D,Int}, lengths) where {D}
    len = NTuple{D,Float64}(Float64(value) for value in lengths)
    waves = ntuple(axis -> _bvk1_spectral_wavenumbers_1d(
        dims[axis], len[axis]), D)
    k2 = _periodic_k2_nd(dims, collect(len))
    source_hat = zeros(ComplexF64, dims)
    work_hat = zeros(ComplexF64, dims)
    return _BVK1SpectralFFTWorkspace{D}(len, waves, k2,
        source_hat, work_hat, plan_fft!(source_hat),
        plan_ifft!(source_hat), plan_fft!(work_hat), plan_ifft!(work_hat))
end

const _BVK1_SPECTRAL_FFT_WS_CACHE = Dict{Tuple{Vararg{Int}},Any}()
const _BVK1_SPECTRAL_FFT_WS_CACHE_LOCK = ReentrantLock()

function _bvk1_spectral_fft_workspace_cached(
        dims::NTuple{D,Int}, lengths) where {D}
    key = (Threads.threadid(), dims...)
    workspace = lock(_BVK1_SPECTRAL_FFT_WS_CACHE_LOCK) do
        get!(() -> _bvk1_spectral_fft_workspace(dims, lengths),
            _BVK1_SPECTRAL_FFT_WS_CACHE, key)
    end::_BVK1SpectralFFTWorkspace{D}
    len = NTuple{D,Float64}(Float64(value) for value in lengths)
    if len != workspace.lengths
        workspace.lengths = len
        workspace.waves = ntuple(axis -> _bvk1_spectral_wavenumbers_1d(
            dims[axis], len[axis]), D)
        workspace.k2 .= _periodic_k2_nd(dims, collect(len))
    end
    return workspace
end

function _bvk1_map_translation_choices(operation::Function, choices)
    count = length(choices)
    if count <= 1 || Threads.nthreads() == 1
        return map(operation, choices)
    end
    results = Vector{Any}(undef, count)
    Threads.@threads :static for index in eachindex(choices)
        results[index] = operation(choices[index])
    end
    return results
end

@inline function _bvk1_copy_real!(destination::AbstractArray{Float64},
        source::AbstractArray{ComplexF64})
    @inbounds @simd for index in eachindex(destination, source)
        destination[index] = real(source[index])
    end
    return destination
end

"""Return spectral gradients and, optionally, the spectral Laplacian.

All components share the same forward transform of `values`.  This matters in
the three-dimensional relaxation hot path: computing each component through
`_bvk1_spectral_derivative` would repeat an identical forward FFT `D` times.
The inverse transforms cannot be shared because the component fields differ.
"""
function _bvk1_spectral_gradient_laplacian_impl(
        values::AbstractArray{Float64,D},
        lengths::AbstractVector{Float64}, ::Val{L}) where {D,L}
    length(lengths) == D || throw(ArgumentError(
        "lengths must have one entry per array dimension"))
    workspace = _bvk1_spectral_fft_workspace_cached(size(values), lengths)
    @inbounds @simd for index in eachindex(workspace.source_hat, values)
        workspace.source_hat[index] = complex(values[index])
    end
    workspace.source_forward * workspace.source_hat
    gradients = Vector{Array{Float64,D}}(undef, D)
    for axis in 1:D
        wave = workspace.waves[axis]
        @inbounds for index in CartesianIndices(values)
            workspace.work_hat[index] = im * wave[index[axis]] *
                workspace.source_hat[index]
        end
        workspace.work_inverse * workspace.work_hat
        gradients[axis] = _bvk1_copy_real!(
            Array{Float64,D}(undef, size(values)), workspace.work_hat)
    end
    if !L
        return gradients, nothing
    end
    @inbounds @simd for index in eachindex(workspace.work_hat)
        workspace.work_hat[index] = -workspace.k2[index] *
            workspace.source_hat[index]
    end
    workspace.work_inverse * workspace.work_hat
    laplacian_values = _bvk1_copy_real!(
        Array{Float64,D}(undef, size(values)), workspace.work_hat)
    return gradients, laplacian_values
end

function _bvk1_spectral_gradient_laplacian(
        values::AbstractArray{Float64,D},
        lengths::AbstractVector{Float64};
        laplacian::Bool=true) where {D}
    return _bvk1_spectral_gradient_laplacian_impl(
        values, lengths, Val(laplacian))
end

"""Apply the divergence of a vector field with one fused inverse FFT."""
function _bvk1_spectral_divergence(
        components::AbstractVector{<:AbstractArray{Float64,D}},
        lengths::AbstractVector{Float64}) where {D}
    length(components) == D || throw(DimensionMismatch(
        "divergence requires one component per array dimension"))
    length(lengths) == D || throw(ArgumentError(
        "lengths must have one entry per array dimension"))
    dims = size(first(components))
    all(component -> size(component) == dims, components) ||
        throw(DimensionMismatch("divergence component sizes must agree"))
    workspace = _bvk1_spectral_fft_workspace_cached(dims, lengths)
    fill!(workspace.source_hat, zero(ComplexF64))
    for axis in 1:D
        component = components[axis]
        @inbounds @simd for index in eachindex(workspace.work_hat, component)
            workspace.work_hat[index] = complex(component[index])
        end
        workspace.work_forward * workspace.work_hat
        wave = workspace.waves[axis]
        @inbounds for index in CartesianIndices(workspace.source_hat)
            workspace.source_hat[index] += im * wave[index[axis]] *
                workspace.work_hat[index]
        end
    end
    workspace.source_inverse * workspace.source_hat
    return _bvk1_copy_real!(Array{Float64,D}(undef, dims),
        workspace.source_hat)
end

"""Apply `laplacian(laplacian_adjoint) - div(grad_adjoint)`.

The transpose contributions share one accumulated Fourier spectrum and one
inverse transform.  Each input still needs its own forward transform; linearity
then makes combining all adjoints before the inverse exact to roundoff.
"""
function _bvk1_spectral_laplacian_minus_divergence(
        components::AbstractVector{<:AbstractArray{Float64,D}},
        laplacian_adjoint::AbstractArray{Float64,D},
        lengths::AbstractVector{Float64}) where {D}
    length(components) == D || throw(DimensionMismatch(
        "gradient adjoint requires one component per array dimension"))
    length(lengths) == D || throw(ArgumentError(
        "lengths must have one entry per array dimension"))
    dims = size(laplacian_adjoint)
    all(component -> size(component) == dims, components) ||
        throw(DimensionMismatch("gradient-adjoint component sizes must agree"))
    workspace = _bvk1_spectral_fft_workspace_cached(dims, lengths)
    @inbounds @simd for index in eachindex(
            workspace.source_hat, laplacian_adjoint)
        workspace.source_hat[index] = complex(laplacian_adjoint[index])
    end
    workspace.source_forward * workspace.source_hat
    @inbounds @simd for index in eachindex(workspace.source_hat)
        workspace.source_hat[index] *= -workspace.k2[index]
    end
    for axis in 1:D
        component = components[axis]
        @inbounds @simd for index in eachindex(workspace.work_hat, component)
            workspace.work_hat[index] = complex(component[index])
        end
        workspace.work_forward * workspace.work_hat
        wave = workspace.waves[axis]
        @inbounds for index in CartesianIndices(workspace.source_hat)
            workspace.source_hat[index] -= im * wave[index[axis]] *
                workspace.work_hat[index]
        end
    end
    workspace.source_inverse * workspace.source_hat
    return _bvk1_copy_real!(Array{Float64,D}(undef, dims),
        workspace.source_hat)
end

function _bvk1_spectral_derivatives(values::AbstractArray{Float64,D},
        lengths::AbstractVector{Float64}) where {D}
    return _bvk1_spectral_gradient_laplacian_impl(
        values, lengths, Val(true))
end

_bvk1_spectral_gradient_components(values::AbstractArray{Float64,D},
    lengths::AbstractVector{Float64}) where {D} = first(
        _bvk1_spectral_gradient_laplacian_impl(
            values, lengths, Val(false)))

@inline function _bvk1_spectral_eta_prime(phi::Real, f::Real)
    raw = Float64(phi)
    (raw <= 1.0e-10 || raw >= 1.0 - 1.0e-10) && return 0.0
    value = raw
    ff = Float64(f)
    return 0.5 * (ff / value + (1.0 - ff) / (1.0 - value))
end

function _bvk1_spectral_state(phi::AbstractArray{Float64,D}, f::Float64,
        lengths::AbstractVector{Float64}, N::Float64, b::Float64;
        adaptive::Bool=true, vk1_scale::Real=0.046,
        eta_floor::Real=1.0e-12, denominator_floor::Real=1.0e-12,
        laplacian_smoothing::Real=0.0,
        composition_gradient::Bool=true) where {D}
    scale = N * Float64(vk1_scale)
    eta_floor_value = Float64(eta_floor)
    denominator_floor_value = Float64(denominator_floor)
    smoothing = Float64(laplacian_smoothing)
    scale >= 0.0 && isfinite(scale) || throw(ArgumentError(
        "vk1_scale must be nonnegative and finite"))
    eta_floor_value > 0.0 && isfinite(eta_floor_value) || throw(ArgumentError(
        "eta_floor must be positive and finite"))
    denominator_floor_value > 0.0 && isfinite(denominator_floor_value) ||
        throw(ArgumentError("denominator_floor must be positive and finite"))
    smoothing >= 0.0 && isfinite(smoothing) || throw(ArgumentError(
        "laplacian_smoothing must be nonnegative and finite"))

    eta = map(value -> diblock_bvk1_eta_psi(value; f=f), phi)
    eta_gradients, eta_laplacian = _bvk1_spectral_derivatives(
        eta, lengths)
    grad_eta2 = zeros(Float64, size(phi))
    for axis in 1:D
        grad_eta2 .+= eta_gradients[axis].^2
    end
    eta_denominator = max.(abs.(eta), eta_floor_value)
    abs_laplacian = similar(eta_laplacian)
    abs_laplacian_derivative = similar(eta_laplacian)
    for index in eachindex(eta_laplacian)
        pair = _diblock_bvk1_regularized_abs_pair(
            eta_laplacian[index], smoothing)
        abs_laplacian[index] = pair[1]
        abs_laplacian_derivative[index] = pair[2]
    end
    sensor = scale .* (2.0 .* abs_laplacian ./ eta_denominator .+
        grad_eta2 ./ eta_denominator.^2)
    activation = tanh.(eta.^2 ./ (f * (1.0 - f)))
    k0 = N * b^2 / 24.0
    kpsi = adaptive ?
        k0 .* (1.0 .+ 0.5 .* activation .* tanh.(sensor)) :
        fill(k0, size(phi))
    phi_gradients = composition_gradient ?
        _bvk1_spectral_gradient_components(phi, lengths) : nothing
    grad_phi2 = composition_gradient ? zeros(Float64, size(phi)) : nothing
    if composition_gradient
        for component in phi_gradients
            grad_phi2 .+= component.^2
        end
    end
    raw_denominator = composition_gradient ? phi .* (1.0 .- phi) : nothing
    phi_denominator = composition_gradient ?
        max.(raw_denominator, denominator_floor_value) : nothing
    quotient = composition_gradient ? grad_phi2 ./ phi_denominator : nothing
    return (; eta, eta_gradients, eta_laplacian, grad_eta2,
        phi_gradients, grad_phi2, eta_denominator, abs_laplacian,
        abs_laplacian_derivative, sensor, activation, kpsi,
        raw_denominator, phi_denominator, quotient, k0, scale,
        eta_floor=eta_floor_value, denominator_floor=denominator_floor_value)
end

"""Continuum-collocation BVK1 energy with spectral spatial derivatives.

This evaluates the same continuum adaptive-stiffness law as BVK1, while
removing the registry-dependent mixture of central node sensors and forward
edge gradients used by the legacy finite-difference discretization.
"""
function diblock_bvk1_spectral_energy_nd(phi_a::AbstractArray{<:Real};
        f::Real=0.5, chiN::Real=12.0, lengths=1.0, N::Real=1.0,
        b::Real=1.0, adaptive::Bool=true, vk1_scale::Real=0.046,
        eta_floor::Real=1.0e-12, denominator_floor::Real=1.0e-12,
        laplacian_smoothing::Real=0.0, mean_tolerance::Real=1.0e-8)
    phi, dims, ff, chi, cell_lengths, nn, bb, dV =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    k2_values = _periodic_k2_nd(dims, cell_lengths)
    coefficients = diblock_bvk1_coefficients(; f=ff, N=nn, b=bb)
    psi = phi .- ff
    connectivity = 0.5 * coefficients.A_psi *
        _periodic_inverse_laplacian_pair_integral_nd(
            psi, psi, k2_values, dV)
    local_entropy = coefficients.M_psi * dV * sum(
        _diblock_relative_entropy(value, ff) for value in phi)
    state = _bvk1_spectral_state(phi, ff, cell_lengths, nn, bb;
        adaptive=adaptive, vk1_scale=vk1_scale, eta_floor=eta_floor,
        denominator_floor=denominator_floor,
        laplacian_smoothing=laplacian_smoothing)
    gradient = dV * sum(state.kpsi .* state.quotient)
    interaction = chi * dV * sum(
        value * (1.0 - value) - ff * (1.0 - ff) for value in phi)
    return connectivity + local_entropy + gradient + interaction
end

function diblock_bvk1_spectral_energy_1d(phi_a::AbstractVector{<:Real};
        f::Real=0.5, chiN::Real=12.0, L::Real=1.0, kwargs...)
    return diblock_bvk1_spectral_energy_nd(phi_a; f=f, chiN=chiN,
        lengths=(Float64(L),), kwargs...)
end

"""Exact nodal chemical potential of `diblock_bvk1_spectral_energy_nd`.

The reverse pass differentiates both the direct composition-gradient quotient
and the complete adaptive dependence of `K_psi` through eta, its spectral
gradient, and its spectral Laplacian.  It returns `(1/dV) dF/dphi`.
"""
function diblock_bvk1_spectral_chemical_potential_nd(
        phi_a::AbstractArray{<:Real}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, N::Real=1.0, b::Real=1.0,
        adaptive::Bool=true, vk1_scale::Real=0.046,
        eta_floor::Real=1.0e-12, denominator_floor::Real=1.0e-12,
        laplacian_smoothing::Real=0.0, mean_tolerance::Real=1.0e-8)
    phi, dims, ff, chi, cell_lengths, nn, bb, _dV =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    coefficients = diblock_bvk1_coefficients(; f=ff, N=nn, b=bb)
    k2_values = _periodic_k2_nd(dims, cell_lengths)
    state = _bvk1_spectral_state(phi, ff, cell_lengths, nn, bb;
        adaptive=adaptive, vk1_scale=vk1_scale, eta_floor=eta_floor,
        denominator_floor=denominator_floor,
        laplacian_smoothing=laplacian_smoothing)

    chemical = coefficients.A_psi .* _periodic_inverse_laplacian_apply_nd(
        phi .- ff, k2_values)
    for index in eachindex(phi)
        value = phi[index]
        entropy_derivative = value <= 0.0 ? -Inf : value >= 1.0 ? Inf :
            log(value * (1.0 - ff) / ((1.0 - value) * ff))
        chemical[index] += coefficients.M_psi * entropy_derivative +
            chi * (1.0 - 2.0 * value)
    end

    # Direct derivative of K |grad(phi)|^2 / [phi(1-phi)].
    active_denominator = state.raw_denominator .> state.denominator_floor
    chemical .+= ifelse.(active_denominator,
        -state.kpsi .* state.grad_phi2 .* (1.0 .- 2.0 .* phi) ./
            state.phi_denominator.^2,
        0.0)
    gradient_adjoints = [2.0 .* state.kpsi .* component ./
        state.phi_denominator for component in state.phi_gradients]
    chemical .-= _bvk1_spectral_divergence(
        gradient_adjoints, cell_lengths)

    if adaptive
        tanh_sensor = tanh.(state.sensor)
        sensor_adjoint = state.quotient .* (0.5 * state.k0) .*
            state.activation .* (1.0 .- tanh_sensor.^2)
        activation_adjoint = state.quotient .* (0.5 * state.k0) .*
            tanh_sensor
        activation_scale = ff * (1.0 - ff)
        eta_adjoint = activation_adjoint .* (1.0 .- state.activation.^2) .*
            (2.0 .* state.eta ./ activation_scale)

        laplacian_adjoint = sensor_adjoint .* state.scale .* 2.0 ./
            state.eta_denominator .* state.abs_laplacian_derivative
        denominator_adjoint = sensor_adjoint .* state.scale .* (
            -2.0 .* state.abs_laplacian ./ state.eta_denominator.^2 .-
            2.0 .* state.grad_eta2 ./ state.eta_denominator.^3)
        eta_adjoint .+= ifelse.(abs.(state.eta) .> state.eta_floor,
            denominator_adjoint .* sign.(state.eta), 0.0)
        eta_gradient_adjoints = [sensor_adjoint .* state.scale .* 2.0 .*
            component ./ state.eta_denominator.^2
            for component in state.eta_gradients]
        eta_adjoint .+= _bvk1_spectral_laplacian_minus_divergence(
            eta_gradient_adjoints, laplacian_adjoint, cell_lengths)
        for index in eachindex(phi)
            chemical[index] += eta_adjoint[index] *
                _bvk1_spectral_eta_prime(phi[index], ff)
        end
    end
    return chemical
end

function _bvk1_spectral_theta_energy_gradient_fine(
        theta::AbstractArray{Float64,D}; f::Real, chiN::Real, lengths,
        N::Real, b::Real, adaptive::Bool, vk1_scale::Real,
        eta_floor::Real, denominator_floor::Real,
        laplacian_smoothing::Real, gradient::Bool) where {D}
    phi = sin.(theta).^2
    physical, dims, ff, chi, cell_lengths, nn, bb, dV =
        _diblock_profile_energy_inputs_nd(phi; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=2.0e-12)
    coefficients = diblock_bvk1_coefficients(; f=ff, N=nn, b=bb)
    k2_values = _periodic_k2_nd(dims, cell_lengths)
    psi = physical .- ff
    connectivity = 0.5 * coefficients.A_psi *
        _periodic_inverse_laplacian_pair_integral_nd(
            psi, psi, k2_values, dV)
    local_entropy = coefficients.M_psi * dV * sum(
        _diblock_relative_entropy(value, ff) for value in physical)
    interaction = chi * dV * sum(
        value * (1.0-value) - ff * (1.0-ff) for value in physical)
    state = _bvk1_spectral_state(physical, ff, cell_lengths, nn, bb;
        adaptive=adaptive, vk1_scale=vk1_scale, eta_floor=eta_floor,
        denominator_floor=denominator_floor,
        laplacian_smoothing=laplacian_smoothing,
        composition_gradient=false)
    theta_gradients = _bvk1_spectral_gradient_components(theta, cell_lengths)
    grad_theta2 = zeros(Float64, dims)
    for component in theta_gradients
        grad_theta2 .+= component.^2
    end
    # In the UD angle chart phi=sin(theta)^2, so the BVK1 outer gradient
    # quotient is exactly 4|grad(theta)|^2 in the continuum.  Using this form
    # avoids an artificial 0/0 wall and a discrete chain-rule alias.
    gradient_energy = 4.0 * dV * sum(state.kpsi .* grad_theta2)
    total = connectivity + local_entropy + interaction + gradient_energy
    gradient || return (energy=total, theta_gradient=nothing,
        phi=physical, cell_lengths=cell_lengths)

    sine2 = sin.(2.0 .* theta)
    theta_gradient = (coefficients.A_psi .*
        _periodic_inverse_laplacian_apply_nd(psi, k2_values) .+
        chi .* (1.0 .- 2.0 .* physical)) .* sine2
    for index in eachindex(physical)
        value = physical[index]
        # The composite wall limit is theta*log(theta^2) -> 0.  Evaluating
        # this product directly preserves the derivative of the xlogy entropy
        # without an inconsistent density clamp.
        if 0.0 < value < 1.0
            theta_gradient[index] += coefficients.M_psi * sine2[index] *
                log(value * (1.0-ff) / ((1.0-value)*ff))
        end
    end
    theta_gradient .-= 8.0 .* _bvk1_spectral_divergence(
        [state.kpsi .* component for component in theta_gradients],
        cell_lengths)

    if adaptive
        k_adjoint = 4.0 .* grad_theta2
        tanh_sensor = tanh.(state.sensor)
        sensor_adjoint = k_adjoint .* (0.5 * state.k0) .*
            state.activation .* (1.0 .- tanh_sensor.^2)
        activation_adjoint = k_adjoint .* (0.5 * state.k0) .* tanh_sensor
        activation_scale = ff * (1.0-ff)
        eta_adjoint = activation_adjoint .* (1.0 .- state.activation.^2) .*
            (2.0 .* state.eta ./ activation_scale)
        laplacian_adjoint = sensor_adjoint .* state.scale .* 2.0 ./
            state.eta_denominator .* state.abs_laplacian_derivative
        denominator_adjoint = sensor_adjoint .* state.scale .* (
            -2.0 .* state.abs_laplacian ./ state.eta_denominator.^2 .-
            2.0 .* state.grad_eta2 ./ state.eta_denominator.^3)
        eta_adjoint .+= ifelse.(abs.(state.eta) .> state.eta_floor,
            denominator_adjoint .* sign.(state.eta), 0.0)
        eta_gradient_adjoints = [sensor_adjoint .* state.scale .* 2.0 .*
            component ./ state.eta_denominator.^2
            for component in state.eta_gradients]
        eta_adjoint .+= _bvk1_spectral_laplacian_minus_divergence(
            eta_gradient_adjoints, laplacian_adjoint, cell_lengths)
        for index in eachindex(physical)
            theta_gradient[index] += eta_adjoint[index] *
                _bvk1_spectral_eta_prime(physical[index], ff) * sine2[index]
        end
    end
    return (energy=total, theta_gradient=theta_gradient,
        phi=physical, cell_lengths=cell_lengths)
end

function _bvk1_dealiased_theta_state(theta::AbstractArray{<:Real,D};
        f::Real, lengths, target_dims, oversampling::Real) where {D}
    coarse = Float64.(theta)
    dims = if target_dims === nothing
        factor = Float64(oversampling)
        factor >= 1.0 && isfinite(factor) || throw(ArgumentError(
            "oversampling must be finite and at least one"))
        ntuple(axis -> max(size(coarse, axis),
            ceil(Int, factor * size(coarse, axis))), D)
    else
        Tuple(Int.(target_dims))
    end
    length(dims) == D || throw(DimensionMismatch(
        "target_dims must have one entry per theta dimension"))
    all(axis -> dims[axis] >= size(coarse, axis), 1:D) ||
        throw(ArgumentError("target_dims may not undersample theta"))
    cell_lengths = _diblock_nd_lengths(lengths, size(coarse))
    fine_theta = _bvk1_cell_centered_resample(coarse, dims)
    _bvk1_theta_uniform_mean_retract!(fine_theta, f)
    return (coarse=coarse, dims=dims, cell_lengths=cell_lengths,
        fine_theta=fine_theta, fine_phi=sin.(fine_theta).^2)
end

"""Free-energy density of the padded, cell-centred BVK1 theta objective."""
function diblock_bvk1_dealiased_theta_energy(
        theta::AbstractArray{<:Real}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, target_dims=nothing, oversampling::Real=2.0,
        translation_phases::Integer=1,
        N::Real=1.0, b::Real=1.0, adaptive::Bool=true,
        vk1_scale::Real=0.046, eta_floor::Real=1.0e-12,
        denominator_floor::Real=1.0e-12,
        laplacian_smoothing::Real=0.0)
    choices = _bvk1_translation_choices(ndims(theta), translation_phases)
    if length(choices) > 1
        energies = _bvk1_map_translation_choices(choices) do choice
            shifted = _bvk1_spectral_shift(theta, Tuple(choice))
            diblock_bvk1_dealiased_theta_energy(shifted; f=f, chiN=chiN,
                lengths=lengths, target_dims=target_dims,
                oversampling=oversampling, translation_phases=1, N=N, b=b,
                adaptive=adaptive, vk1_scale=vk1_scale,
                eta_floor=eta_floor, denominator_floor=denominator_floor,
                laplacian_smoothing=laplacian_smoothing)
        end
        return mean(energies)
    end
    state = _bvk1_dealiased_theta_state(theta; f=f, lengths=lengths,
        target_dims=target_dims, oversampling=oversampling)
    fine = _bvk1_spectral_theta_energy_gradient_fine(state.fine_theta;
        f=f, chiN=chiN, lengths=state.cell_lengths, N=N, b=b,
        adaptive=adaptive, vk1_scale=vk1_scale, eta_floor=eta_floor,
        denominator_floor=denominator_floor,
        laplacian_smoothing=laplacian_smoothing, gradient=false)
    return fine.energy / prod(state.cell_lengths)
end

"""Evaluate the padded BVK1 objective and its exact coarse-theta gradient.

The optimization variable is a coarse-grid Uneyama--Doi angle field,
`phi=sin(theta)^2`.  It is Fourier-interpolated at cell centres to
`target_dims`, retracted to the requested composition, and only then passed
to the spectral BVK1 functional.  The returned objective is a free-energy
density.  Its gradient includes the fixed-composition tangent projection on
the fine grid and the exact adjoint of the Fourier interpolation.

This padded objective removes the nonlinear collocation alias that otherwise
acts as a grid-pinning potential for translated ordered phases.  Richardson
correction is still required across independently relaxed coarse grids.
"""
function diblock_bvk1_dealiased_theta_energy_gradient(
        theta::AbstractArray{<:Real,D}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, target_dims=nothing, oversampling::Real=2.0,
        translation_phases::Integer=1,
        N::Real=1.0, b::Real=1.0, adaptive::Bool=true,
        vk1_scale::Real=0.046, eta_floor::Real=1.0e-12,
        denominator_floor::Real=1.0e-12,
        laplacian_smoothing::Real=0.0) where {D}
    choices = _bvk1_translation_choices(D, translation_phases)
    if length(choices) > 1
        results = _bvk1_map_translation_choices(choices) do choice
            shifted = _bvk1_spectral_shift(theta, Tuple(choice))
            diblock_bvk1_dealiased_theta_energy_gradient(shifted;
                f=f, chiN=chiN, lengths=lengths,
                target_dims=target_dims, oversampling=oversampling,
                translation_phases=1, N=N, b=b, adaptive=adaptive,
                vk1_scale=vk1_scale, eta_floor=eta_floor,
                denominator_floor=denominator_floor,
                laplacian_smoothing=laplacian_smoothing)
        end
        gradient = zeros(Float64, size(theta))
        for (choice, result) in zip(choices, results)
            gradient .+= _bvk1_spectral_shift(result.gradient,
                ntuple(axis -> -Float64(choice[axis]), D))
        end
        gradient ./= length(results)
        energies = getproperty.(results, :energy_density)
        return merge(first(results), (
            energy_density=mean(energies), gradient=gradient,
            translation_phase_count=length(results),
            translation_energy_spread=maximum(energies)-minimum(energies),
            schema=DIBLOCK_BVK1_DEALIASED_OBJECTIVE_SCHEMA))
    end
    state = _bvk1_dealiased_theta_state(theta; f=f, lengths=lengths,
        target_dims=target_dims, oversampling=oversampling)
    fine = _bvk1_spectral_theta_energy_gradient_fine(state.fine_theta;
        f=f, chiN=chiN, lengths=state.cell_lengths, N=N, b=b,
        adaptive=adaptive, vk1_scale=vk1_scale, eta_floor=eta_floor,
        denominator_floor=denominator_floor,
        laplacian_smoothing=laplacian_smoothing, gradient=true)

    sine2 = sin.(2.0 .* state.fine_theta)
    normal_sum = sum(sine2)
    abs(normal_sum) > eps(Float64) || throw(ArgumentError(
        "theta gradient is frozen because the field has no interface cells"))
    multiplier = sum(fine.theta_gradient) / normal_sum
    # Pull back the exact uniform-shift composition retraction, then convert
    # total energy to energy density.
    fine_gradient = (fine.theta_gradient .- multiplier .* sine2) ./
        prod(state.dims)
    coarse_gradient = _bvk1_cell_centered_resample_adjoint(
        fine_gradient, size(state.coarse))
    volume = prod(state.cell_lengths)
    return (energy_density=fine.energy / volume,
        gradient=coarse_gradient, fine_theta=state.fine_theta,
        fine_phi=state.fine_phi, target_dims=state.dims,
        composition_multiplier=multiplier,
        mean_error=mean(state.fine_phi)-Float64(f),
        translation_phase_count=1,
        translation_energy_spread=0.0,
        schema=DIBLOCK_BVK1_DEALIASED_OBJECTIVE_SCHEMA)
end
