# ===========================================================================
# Oversampled ("honest") quadrature for the BVK2 energy.
#
# WHY.  The θ-chart field solver represents the density as φ = sin²θ on a
# periodic grid.  Even when θ is band-limited to the grid Nyquist, φ = sin²θ =
# (1−cos2θ)/2 carries Fourier content up to *twice* the Nyquist, so the native
# grid quadrature of the strongly nonlinear energy terms (the singular-weight
# gradient integral ∫K(∇φ)²/[φ(1−φ)], the RPA connectivity, the entropy, the
# interaction) ALIASES that content.  Above χN≈49 the aliased near-Nyquist band
# is mis-priced (a measured sign flip in the gradient term), and the θ-BB
# minimiser slides down a spurious grid-scale mode: the energy diverges downward
# with grid refinement instead of converging (see
# docs/manuscript/macromolecules/working_note.md §7.5).
#
# FIX.  Evaluate the energy on an m×-refined grid.  Upsample θ by exact Fourier
# interpolation to the fine grid, form φ = sin²(θ_fine) there, and integrate the
# *existing* energy on the fine grid.  The fine grid resolves φ's content up to
# m× the coarse Nyquist, so the quadrature no longer aliases and the energy
# grid-converges.  The optimisation variable stays the coarse θ; the θ-gradient
# is obtained by the EXACT adjoint of the upsampler, so it is finite-difference
# consistent with the oversampled energy (unlike a truncation-based downsample).
#
# COST.  Inherently ~mᵈ× the native FFT work (the fine grid has mᵈ× the points);
# the fine-grid energy/chemical-potential FFTs dominate and are already cached
# (`_bvk2_spectral_workspace_cached`, `_bvk2_chempot_workspace_cached`).  The
# resampler here caches its own plans and reuses preallocated buffers so it adds
# ZERO allocation per call (a FourierTools.resample round trip allocated ~57 MB
# per 48³ call — ~TBs of GC over a full solve).
# ===========================================================================

const DIBLOCK_BVK2_UD_THETA_DISCRETIZATION_SCHEMA =
    "ud-theta-spectral-adaptive-k-oversampled-v1"
const DIBLOCK_BVK2_UD_THETA_FILTERED_DISCRETIZATION_SCHEMA =
    "ud-theta-spectral-adaptive-k-physical-sensor-filter-oversampled-v2"

# The adaptive stiffness sensor contains two derivatives of eta.  Without a
# physical ultraviolet scale, a spectrally tiny high-k tail is amplified by
# k^2 and changes the discrete constitutive law as the grid is refined.  Use a
# fixed continuum cutoff k_c = ratio*k_star and the screened-biharmonic
# transfer H(k)=1/(1+(k/k_c)^4), equivalently
# (I+k_c^-4*Delta^2)*eta_sensor=eta. With ratio=12 the ordering band is
# essentially unchanged (H(k_star)=0.99995) while high-k curvature decays as
# k^-2. The global ratio 12 is the largest cutoff that passes the complete
# lamellar period/profile nx=128/256 audit; it is not promoted for competing-
# morphology thermodynamics. `Inf` selects the legacy sensor.
const DIBLOCK_BVK2_UD_SENSOR_FILTER_RATIO = 12.0

function _bvk2_ud_sensor_filter_ratio(value::Real)
    ratio = Float64(value)
    (isfinite(ratio) && ratio > 0.0) || ratio == Inf ||
        throw(ArgumentError(
            "sensor_filter_ratio must be positive or Inf, got $value"))
    return ratio
end

_bvk2_ud_theta_discretization_schema(sensor_filter_ratio::Real) =
    _bvk2_ud_sensor_filter_ratio(sensor_filter_ratio) == Inf ?
    DIBLOCK_BVK2_UD_THETA_DISCRETIZATION_SCHEMA :
    DIBLOCK_BVK2_UD_THETA_FILTERED_DISCRETIZATION_SCHEMA

"""Apply the real-even physical sensor filter; this operator is self-adjoint."""
function _bvk2_ud_sensor_filter(values::AbstractArray{<:Real,D}, lengths,
        k_star2::Real, sensor_filter_ratio::Real) where {D}
    ratio = _bvk2_ud_sensor_filter_ratio(sensor_filter_ratio)
    ratio == Inf && return values isa Array{Float64,D} ? values :
        Array{Float64,D}(values)
    workspace = _bvk1_spectral_fft_workspace_cached(size(values), lengths)
    @inbounds @. workspace.source_hat = complex(values)
    workspace.source_forward * workspace.source_hat
    cutoff2 = ratio^2 * Float64(k_star2)
    @inbounds for index in eachindex(workspace.source_hat)
        reduced_k2 = workspace.k2[index] / cutoff2
        workspace.source_hat[index] /= 1.0 + reduced_k2^2
    end
    workspace.source_inverse * workspace.source_hat
    return _bvk1_copy_real!(Array{Float64,D}(undef, size(values)),
        workspace.source_hat)
end

# Per-axis coarse→fine Fourier-index map for an even coarse length `n` upsampled
# to fine length `M = m*n`.  A coarse mode at frequency `q` lands at the same
# integer frequency on the fine grid; the coarse Nyquist mode (q = n/2) is real
# and split half to +Nyquist and half to −Nyquist so the interpolant stays real
# and the operator is its own transpose's partner.
struct _BVK2OSAxis
    fmap::Vector{Int}   # coarse index k -> primary fine index (Nyq -> +Nyq)
    nyq_c::Int          # coarse Nyquist index (n÷2 + 1)
    pos::Int            # +Nyquist fine index (n÷2 + 1)
    neg::Int            # −Nyquist fine index (M − n÷2 + 1)
end

function _bvk2_os_axis(n::Int, M::Int)
    iseven(n) || throw(ArgumentError("oversampled quadrature needs even grid lengths, got $n"))
    fmap = Vector{Int}(undef, n)
    @inbounds for k in 1:n
        q = (k - 1) <= n ÷ 2 ? (k - 1) : (k - 1 - n)
        fmap[k] = q >= 0 ? q + 1 : M + q + 1
    end
    return _BVK2OSAxis(fmap, n ÷ 2 + 1, n ÷ 2 + 1, M - n ÷ 2 + 1)
end

# Cached workspace for a fixed (coarse dims, factor).  Holds the four FFT plans
# (coarse/fine forward+inverse) and the coarse/fine complex + real buffers.
mutable struct _BVK2OversampleWS{D}
    dims::NTuple{D,Int}
    fdims::NTuple{D,Int}
    factor::Int
    ax::NTuple{D,_BVK2OSAxis}
    coarse_spec::Array{ComplexF64,D}   # coarse spectrum buffer
    fine_spec::Array{ComplexF64,D}     # fine spectrum buffer
    fine_real::Array{Float64,D}        # real fine output (upsample)
    coarse_real::Array{Float64,D}      # real coarse output (adjoint)
    plan_fwd_coarse::Any
    plan_inv_fine::Any
    plan_fwd_fine::Any
    plan_inv_coarse::Any
    fine_lin::LinearIndices{D,NTuple{D,Base.OneTo{Int}}}
end

function _bvk2_oversample_ws(dims::NTuple{D,Int}, factor::Int) where {D}
    factor >= 2 || throw(ArgumentError("oversample factor must be ≥ 2, got $factor"))
    fdims = ntuple(d -> factor * dims[d], D)
    ax = ntuple(d -> _bvk2_os_axis(dims[d], fdims[d]), D)
    coarse_spec = zeros(ComplexF64, dims)
    fine_spec = zeros(ComplexF64, fdims)
    return _BVK2OversampleWS{D}(dims, fdims, factor, ax,
        coarse_spec, fine_spec,
        Array{Float64,D}(undef, fdims), Array{Float64,D}(undef, dims),
        plan_fft!(coarse_spec), plan_ifft!(fine_spec),
        plan_fft!(copy(fine_spec)), plan_ifft!(copy(coarse_spec)),
        LinearIndices(fdims))
end

# Per-process cache keyed by (dims, factor).  Safe without a lock for the same
# reason as `_bvk2_spectral_workspace_cached`: the morphology solver drives
# resampling sequentially on its own thread; campaign parallelism is
# process-level.
const _BVK2_OVERSAMPLE_WS_CACHE = Dict{Tuple{Any,Int},Any}()

function _bvk2_oversample_ws_cached(dims::NTuple{D,Int}, factor::Int) where {D}
    return get!(() -> _bvk2_oversample_ws(dims, factor),
        _BVK2_OVERSAMPLE_WS_CACHE, (dims, factor))::_BVK2OversampleWS{D}
end

# Per-axis (index, weight) targets for a coarse index `k` on one axis: a single
# target for a regular mode, or the two split Nyquist targets.
@inline function _bvk2_os_targets(a::_BVK2OSAxis, k::Int)
    k == a.nyq_c && return (a.pos, 0.5, a.neg, 0.5, 2)
    return (a.fmap[k], 1.0, 0, 0.0, 1)
end

# Scatter the coarse spectrum into the (zeroed) fine spectrum with Nyquist
# splitting.  D ≤ 3; the Nyquist hyperplanes are the only indices with >1 target.
function _bvk2_os_scatter!(ws::_BVK2OversampleWS{D}) where {D}
    fill!(ws.fine_spec, zero(ComplexF64))
    ax = ws.ax; dst = ws.fine_spec; src = ws.coarse_spec; flin = ws.fine_lin
    @inbounds for ci in CartesianIndices(ws.dims)
        val = src[ci]
        if D == 1
            t = _bvk2_os_targets(ax[1], ci[1])
            dst[t[1]] += t[2] * val
            t[5] == 2 && (dst[t[3]] += t[4] * val)
        elseif D == 2
            a1 = _bvk2_os_targets(ax[1], ci[1]); a2 = _bvk2_os_targets(ax[2], ci[2])
            for u in 1:a1[5], v in 1:a2[5]
                iu = u == 1 ? a1[1] : a1[3]; wu = u == 1 ? a1[2] : a1[4]
                iv = v == 1 ? a2[1] : a2[3]; wv = v == 1 ? a2[2] : a2[4]
                dst[flin[iu, iv]] += (wu * wv) * val
            end
        else
            a1 = _bvk2_os_targets(ax[1], ci[1]); a2 = _bvk2_os_targets(ax[2], ci[2])
            a3 = _bvk2_os_targets(ax[3], ci[3])
            for u in 1:a1[5], v in 1:a2[5], w in 1:a3[5]
                iu = u == 1 ? a1[1] : a1[3]; wu = u == 1 ? a1[2] : a1[4]
                iv = v == 1 ? a2[1] : a2[3]; wv = v == 1 ? a2[2] : a2[4]
                iw = w == 1 ? a3[1] : a3[3]; ww = w == 1 ? a3[2] : a3[4]
                dst[flin[iu, iv, iw]] += (wu * wv * ww) * val
            end
        end
    end
    return ws.fine_spec
end

# Gather (exact adjoint of scatter): coarse[ci] = Σ weight * fine[target].
function _bvk2_os_gather!(ws::_BVK2OversampleWS{D}) where {D}
    ax = ws.ax; fine = ws.fine_spec; dst = ws.coarse_spec; flin = ws.fine_lin
    @inbounds for ci in CartesianIndices(ws.dims)
        acc = zero(ComplexF64)
        if D == 1
            t = _bvk2_os_targets(ax[1], ci[1])
            acc += t[2] * fine[t[1]]
            t[5] == 2 && (acc += t[4] * fine[t[3]])
        elseif D == 2
            a1 = _bvk2_os_targets(ax[1], ci[1]); a2 = _bvk2_os_targets(ax[2], ci[2])
            for u in 1:a1[5], v in 1:a2[5]
                iu = u == 1 ? a1[1] : a1[3]; wu = u == 1 ? a1[2] : a1[4]
                iv = v == 1 ? a2[1] : a2[3]; wv = v == 1 ? a2[2] : a2[4]
                acc += (wu * wv) * fine[flin[iu, iv]]
            end
        else
            a1 = _bvk2_os_targets(ax[1], ci[1]); a2 = _bvk2_os_targets(ax[2], ci[2])
            a3 = _bvk2_os_targets(ax[3], ci[3])
            for u in 1:a1[5], v in 1:a2[5], w in 1:a3[5]
                iu = u == 1 ? a1[1] : a1[3]; wu = u == 1 ? a1[2] : a1[4]
                iv = v == 1 ? a2[1] : a2[3]; wv = v == 1 ? a2[2] : a2[4]
                iw = w == 1 ? a3[1] : a3[3]; ww = w == 1 ? a3[2] : a3[4]
                acc += (wu * wv * ww) * fine[flin[iu, iv, iw]]
            end
        end
        dst[ci] = acc
    end
    return ws.coarse_spec
end

"""
    _bvk2_fourier_upsample!(ws, x) -> ws.fine_real

Exact factor-`m` Fourier interpolation of the real coarse array `x` onto the
fine grid, reproducing `FourierTools.resample(x, m.*size(x))` to roundoff.
Reuses `ws` buffers; allocates nothing.
"""
function _bvk2_fourier_upsample!(ws::_BVK2OversampleWS{D},
        x::AbstractArray{Float64,D}) where {D}
    size(x) == ws.dims || throw(DimensionMismatch("upsample input size mismatch"))
    @inbounds @. ws.coarse_spec = complex(x)
    ws.plan_fwd_coarse * ws.coarse_spec
    _bvk2_os_scatter!(ws)
    ws.plan_inv_fine * ws.fine_spec
    scale = prod(ws.fdims) / prod(ws.dims)      # amplitude-preserving normalization
    @inbounds @. ws.fine_real = scale * real(ws.fine_spec)
    return ws.fine_real
end

"""
    _bvk2_fourier_upsample_adjoint!(ws, g) -> ws.coarse_real

Exact adjoint of [`_bvk2_fourier_upsample!`](@ref): given a real fine array `g`,
returns the coarse array `h` with `⟨upsample(x), g⟩ == ⟨x, h⟩` for all `x`.
Reuses `ws` buffers; allocates nothing.
"""
function _bvk2_fourier_upsample_adjoint!(ws::_BVK2OversampleWS{D},
        g::AbstractArray{Float64,D}) where {D}
    size(g) == ws.fdims || throw(DimensionMismatch("adjoint input size mismatch"))
    @inbounds @. ws.fine_spec = complex(g)
    ws.plan_fwd_fine * ws.fine_spec
    _bvk2_os_gather!(ws)
    ws.plan_inv_coarse * ws.coarse_spec
    @inbounds @. ws.coarse_real = real(ws.coarse_spec)
    return ws.coarse_real
end

"""
    diblock_bvk2_energy_oversampled_nd(theta_a; f, chiN, lengths, c2,
                                       oversample=2, adaptive=true,
                                       mean_tolerance=1e-4) -> Float64

Honest (`oversample`×-quadrature) BVK2 energy density of the coarse θ-chart
field `theta_a` (so `φ = sin²θ`).  θ is Fourier-upsampled by `oversample`, then
the standard [`diblock_bvk2_energy_nd`](@ref) is evaluated on `φ = sin²(θ_fine)`
at the SAME physical `lengths`.  `oversample=1` bypasses upsampling and returns
the native energy of `sin²(theta_a)` exactly.  Composition is NOT re-projected;
pass a θ whose `sin²` already integrates to `f` (the solver mean-retracts).
"""
function diblock_bvk2_energy_oversampled_nd(theta_a::AbstractArray{<:Real};
        f::Real=0.5, chiN::Real=12.0, lengths=1.0, N::Real=1.0, b::Real=1.0,
        adaptive::Bool=true, c2::Real, oversample::Integer=2,
        mean_tolerance::Real=1.0e-4)
    dims = size(theta_a)
    len = _bvk2_oversample_lengths(lengths, length(dims))
    if oversample <= 1
        phi = sin.(Float64.(theta_a)).^2
        return diblock_bvk2_energy_nd(phi; f=f, chiN=chiN, lengths=len,
            N=N, b=b, adaptive=adaptive, c2=c2, mean_tolerance=mean_tolerance)
    end
    ws = _bvk2_oversample_ws_cached(dims, Int(oversample))
    theta_fine = _bvk2_fourier_upsample!(ws, Array{Float64}(theta_a))
    phi_fine = sin.(theta_fine).^2
    return diblock_bvk2_energy_nd(phi_fine; f=f, chiN=chiN, lengths=len,
        N=N, b=b, adaptive=adaptive, c2=c2, mean_tolerance=mean_tolerance)
end

"""
    diblock_bvk2_theta_gradient_oversampled_nd(theta_a; f, chiN, lengths, c2,
        oversample=2, adaptive=true, mean_tolerance=1e-4,
        project_mean=false) -> Array

Exact θ-gradient of [`diblock_bvk2_energy_oversampled_nd`](@ref): the array
`dE/dθ` of the same shape as `theta_a`.  Equals
`Uᵀ(dV_fine · μ_fine ⊙ sin(2 θ_fine))`, with `μ_fine` the fine-grid chemical
potential and `Uᵀ` the exact adjoint of the Fourier upsampler — so it is
finite-difference consistent with the oversampled energy.  With
`project_mean=true` the mean-preserving projection (orthogonal to `sin(2θ)`) is
applied, giving the constrained gradient the θ-BB solver descends.
"""
function diblock_bvk2_theta_gradient_oversampled_nd(theta_a::AbstractArray{<:Real};
        f::Real=0.5, chiN::Real=12.0, lengths=1.0, N::Real=1.0, b::Real=1.0,
        adaptive::Bool=true, c2::Real, oversample::Integer=2,
        mean_tolerance::Real=1.0e-4, project_mean::Bool=false)
    dims = size(theta_a)
    len = _bvk2_oversample_lengths(lengths, length(dims))
    theta = Array{Float64}(theta_a)
    if oversample <= 1
        phi = sin.(theta).^2
        mu = diblock_bvk2_chemical_potential_nd(phi; f=f, chiN=chiN, lengths=len,
            N=N, b=b, adaptive=adaptive, c2=c2, mean_tolerance=mean_tolerance)
        dV = prod(Float64.(len)) / prod(dims)
        grad = dV .* mu .* sin.(2.0 .* theta)
        project_mean && _bvk2_project_mean_direction!(grad, theta)
        return grad
    end
    ws = _bvk2_oversample_ws_cached(dims, Int(oversample))
    theta_fine = copy(_bvk2_fourier_upsample!(ws, theta))
    phi_fine = sin.(theta_fine).^2
    mu_fine = diblock_bvk2_chemical_potential_nd(phi_fine; f=f, chiN=chiN,
        lengths=len, N=N, b=b, adaptive=adaptive, c2=c2, mean_tolerance=mean_tolerance)
    dV_fine = prod(Float64.(len)) / prod(ws.fdims)
    g_fine = (dV_fine .* mu_fine) .* sin.(2.0 .* theta_fine)
    grad = copy(_bvk2_fourier_upsample_adjoint!(ws, g_fine))
    project_mean && _bvk2_project_mean_direction!(grad, theta)
    return grad
end

"""
    _diblock_bvk2_ud_theta_fine(theta; f, chiN, lengths, N, b, adaptive, c2,
                                gradient)

Evaluate BVK2 after applying the Uneyama--Doi square-root/angle
transformation *before* spatial discretization.  The continuum identity

```
K(phi) * |grad(phi)|^2 / (phi * (1 - phi))
    = 4 * K(sin(theta)^2) * |grad(theta)|^2
```

is represented directly with spectral derivatives.  This differs from merely
using `phi=sin(theta)^2` as an optimizer chart for the legacy finite-difference
edge quotient.  The returned `theta_gradient_density` is the exact derivative
of the discrete total energy divided by the uniform cell volume `dV`.
"""
function _diblock_bvk2_ud_theta_fine(
        theta::AbstractArray{Float64,D}; f::Real, chiN::Real, lengths,
        N::Real, b::Real, adaptive::Bool, c2::Real,
        gradient::Bool, sensor_filter_ratio::Real=
            Inf) where {D}
    filter_ratio = _bvk2_ud_sensor_filter_ratio(sensor_filter_ratio)
    phi = sin.(theta).^2
    physical, dims, ff, chi, cell_lengths, nn, bb, dV =
        _diblock_profile_energy_inputs_nd(phi; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=1.0e-4)
    params = diblock_bvk2_parameters(; f=ff, N=nn, b=bb, c2=c2)
    k2_values = _periodic_k2_nd(dims, cell_lengths)
    psi = physical .- ff

    connectivity = 0.5 * params.A_psi *
        _periodic_inverse_laplacian_pair_integral_nd(
            psi, psi, k2_values, dV)
    local_entropy = params.M_psi * dV * sum(
        _diblock_relative_entropy(value, ff) for value in physical)
    interaction = chi * dV * sum(
        value * (1.0 - value) - ff * (1.0 - ff) for value in physical)

    eta = map(value -> diblock_bvk1_eta_psi(value; f=ff), physical)
    sensor_eta = adaptive && filter_ratio != Inf ?
        _bvk2_ud_sensor_filter(
            eta, cell_lengths, params.k_star2, filter_ratio) : eta
    eta_gradients, eta_laplacian =
        _bvk1_spectral_derivatives(sensor_eta, cell_lengths)
    eta_grad2 = zeros(Float64, dims)
    for component in eta_gradients
        eta_grad2 .+= component.^2
    end
    eta_radius2 = sensor_eta.^2 .+ params.eta_star2
    reduced_q = eta_laplacian ./
        (params.k_star2 .* sqrt.(eta_radius2))
    reduced_p = eta_grad2 ./ (params.k_star2 .* eta_radius2)
    curvature = reduced_q.^2 ./ sqrt.(1.0 .+ reduced_q.^2)
    sensor = params.c2 .* (2.0 .* curvature .+ reduced_p)
    kpsi = adaptive ?
        params.K_psi0 .* (1.0 .+ 0.5 .* tanh.(sensor)) :
        fill(params.K_psi0, dims)

    theta_gradients =
        _bvk1_spectral_gradient_components(theta, cell_lengths)
    theta_grad2 = zeros(Float64, dims)
    for component in theta_gradients
        theta_grad2 .+= component.^2
    end
    gradient_energy = 4.0 * dV * sum(kpsi .* theta_grad2)
    total = connectivity + local_entropy + interaction + gradient_energy
    gradient || return (
        energy=total,
        theta_gradient_density=nothing,
        phi=physical,
        cell_lengths=cell_lengths,
        kpsi=kpsi,
        sensor_eta=sensor_eta,
        reduced_q=reduced_q,
        reduced_p=reduced_p,
        sensor=sensor,
        sensor_filter_ratio=filter_ratio,
        schema=_bvk2_ud_theta_discretization_schema(filter_ratio),
    )

    sine2 = sin.(2.0 .* theta)
    theta_gradient_density = (
        params.A_psi .* _periodic_inverse_laplacian_apply_nd(
            psi, k2_values) .+
        chi .* (1.0 .- 2.0 .* physical)
    ) .* sine2
    for index in eachindex(physical)
        value = physical[index]
        if 0.0 < value < 1.0
            theta_gradient_density[index] += params.M_psi * sine2[index] *
                log(value * (1.0 - ff) / ((1.0 - value) * ff))
        end
    end

    # Direct derivative of 4*K*|grad(theta)|^2.  The spectral derivative is
    # skew-adjoint, hence d/dtheta = -8*div(K*grad(theta)).
    theta_gradient_density .-= 8.0 .* _bvk1_spectral_divergence(
        [kpsi .* component for component in theta_gradients],
        cell_lengths)

    if adaptive
        # Complete reverse pass through
        # K -> sensor -> (Q,P) -> (laplacian(eta), grad(eta), eta) -> phi.
        k_adjoint = 4.0 .* theta_grad2
        tanh_sensor = tanh.(sensor)
        sensor_adjoint = k_adjoint .* (0.5 * params.K_psi0) .*
            (1.0 .- tanh_sensor.^2)
        curvature_prime = reduced_q .* (2.0 .+ reduced_q.^2) ./
            (1.0 .+ reduced_q.^2).^(1.5)
        q_adjoint = sensor_adjoint .* params.c2 .* 2.0 .*
            curvature_prime
        p_adjoint = sensor_adjoint .* params.c2
        laplacian_adjoint = q_adjoint ./
            (params.k_star2 .* sqrt.(eta_radius2))
        radius2_adjoint = -0.5 .* q_adjoint .* reduced_q ./ eta_radius2 .-
            p_adjoint .* reduced_p ./ eta_radius2
        sensor_eta_adjoint = 2.0 .* sensor_eta .* radius2_adjoint
        grad2_adjoint = p_adjoint ./ (params.k_star2 .* eta_radius2)
        eta_gradient_adjoints = [
            2.0 .* component .* grad2_adjoint
            for component in eta_gradients
        ]
        sensor_eta_adjoint .+= _bvk1_spectral_laplacian_minus_divergence(
            eta_gradient_adjoints, laplacian_adjoint, cell_lengths)
        eta_adjoint = filter_ratio == Inf ? sensor_eta_adjoint :
            _bvk2_ud_sensor_filter(sensor_eta_adjoint,
                cell_lengths, params.k_star2, filter_ratio)
        for index in eachindex(physical)
            value = physical[index]
            eta_prime = _bvk1_spectral_eta_prime(value, ff)
            if eta_prime != 0.0
                theta_gradient_density[index] +=
                    eta_adjoint[index] * eta_prime * sine2[index]
            end
        end
    end

    return (
        energy=total,
        theta_gradient_density=theta_gradient_density,
        phi=physical,
        cell_lengths=cell_lengths,
        kpsi=kpsi,
        sensor_eta=sensor_eta,
        reduced_q=reduced_q,
        reduced_p=reduced_p,
        sensor=sensor,
        sensor_filter_ratio=filter_ratio,
        schema=_bvk2_ud_theta_discretization_schema(filter_ratio),
    )
end

"""
    diblock_bvk2_ud_theta_energy_nd(theta; ..., oversample=2,
                                    sensor_filter_ratio=Inf)

Total BVK2 free energy with the Uneyama--Doi angle transformation applied
before discretization. `theta` is Fourier-interpolated first, then the regular
spectral form `4*K(phi)*|grad(theta)|^2` is evaluated. The adaptive curvature
sensor can use a fixed-physical-scale screened-biharmonic filter with cutoff
`sensor_filter_ratio*k_star`. The default `Inf` retains the legacy raw sensor;
passing a finite ratio opts into the v2 filtered schema. The filter choice and
ratio must be held fixed across competing phases.
"""
function diblock_bvk2_ud_theta_energy_nd(
        theta_a::AbstractArray{<:Real,D}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, N::Real=1.0, b::Real=1.0, adaptive::Bool=true,
        c2::Real, oversample::Integer=2, sensor_filter_ratio::Real=
            Inf) where {D}
    theta = theta_a isa Array{Float64,D} ?
        theta_a : Array{Float64,D}(theta_a)
    cell_lengths = _bvk2_oversample_lengths(lengths, D)
    fine_theta = if oversample <= 1
        theta
    else
        ws = _bvk2_oversample_ws_cached(size(theta), Int(oversample))
        _bvk2_fourier_upsample!(ws, theta)
    end
    return _diblock_bvk2_ud_theta_fine(fine_theta; f=f, chiN=chiN,
        lengths=cell_lengths, N=N, b=b, adaptive=adaptive, c2=c2,
        gradient=false, sensor_filter_ratio=sensor_filter_ratio).energy
end

"""
    diblock_bvk2_ud_theta_value_gradient_nd(theta; ..., oversample=2,
                                            project_mean=false)

Evaluate the filtered UD-theta energy and its exact coarse-grid derivative in
one fine pass. The fused path shares Fourier interpolation and every energy
intermediate already produced by the analytic reverse pass, including the
self-adjoint sensor-filter pullback.
"""
function diblock_bvk2_ud_theta_value_gradient_nd(
        theta_a::AbstractArray{<:Real,D}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, N::Real=1.0, b::Real=1.0, adaptive::Bool=true,
        c2::Real, oversample::Integer=2,
        project_mean::Bool=false, sensor_filter_ratio::Real=
            Inf) where {D}
    theta = theta_a isa Array{Float64,D} ?
        theta_a : Array{Float64,D}(theta_a)
    cell_lengths = _bvk2_oversample_lengths(lengths, D)
    if oversample <= 1
        fine = _diblock_bvk2_ud_theta_fine(theta; f=f, chiN=chiN,
            lengths=cell_lengths, N=N, b=b, adaptive=adaptive, c2=c2,
            gradient=true, sensor_filter_ratio=sensor_filter_ratio)
        dV = prod(cell_lengths) / length(theta)
        grad = dV .* fine.theta_gradient_density
        project_mean && _bvk2_project_mean_direction!(grad, theta)
        return (energy=fine.energy, gradient=grad)
    end

    ws = _bvk2_oversample_ws_cached(size(theta), Int(oversample))
    fine_theta = _bvk2_fourier_upsample!(ws, theta)
    fine = _diblock_bvk2_ud_theta_fine(fine_theta; f=f, chiN=chiN,
        lengths=cell_lengths, N=N, b=b, adaptive=adaptive, c2=c2,
        gradient=true, sensor_filter_ratio=sensor_filter_ratio)
    dV_fine = prod(cell_lengths) / length(fine_theta)
    fine_gradient = fine.theta_gradient_density
    fine_gradient .*= dV_fine
    grad = copy(_bvk2_fourier_upsample_adjoint!(ws, fine_gradient))
    if project_mean
        fine_normal = sin.(2.0 .* fine_theta)
        coarse_normal = copy(
            _bvk2_fourier_upsample_adjoint!(ws, fine_normal))
        normal2 = sum(abs2, coarse_normal)
        multiplier = dot(grad, coarse_normal) / max(normal2, eps(Float64))
        grad .-= multiplier .* coarse_normal
    end
    return (energy=fine.energy, gradient=grad)
end

"""
    diblock_bvk2_ud_theta_gradient_nd(theta; ..., oversample=2,
                                      project_mean=false)

Exact coarse-grid derivative of [`diblock_bvk2_ud_theta_energy_nd`](@ref).
The fine-grid functional is differentiated analytically and pulled back with
the exact adjoint of the Fourier interpolator.
"""
function diblock_bvk2_ud_theta_gradient_nd(
        theta_a::AbstractArray{<:Real,D}; kwargs...) where {D}
    return diblock_bvk2_ud_theta_value_gradient_nd(
        theta_a; kwargs...).gradient
end

function _bvk2_ud_selected_axes(rank::Integer, scaled_axes)
    axis_values = scaled_axes === nothing ? collect(1:rank) :
        scaled_axes isa Integer ? [scaled_axes] : collect(scaled_axes)
    isempty(axis_values) &&
        throw(ArgumentError("scaled_axes must select at least one axis"))
    all(axis -> axis isa Integer && !(axis isa Bool), axis_values) ||
        throw(ArgumentError("scaled_axes must contain integer axis indices"))
    selected = falses(rank)
    for axis_value in axis_values
        axis = Int(axis_value)
        1 <= axis <= rank ||
            throw(ArgumentError("scaled axis $axis is outside 1:$rank"))
        selected[axis] &&
            throw(ArgumentError("scaled_axes must not contain duplicate axes"))
        selected[axis] = true
    end
    return selected
end

"""
    diblock_bvk2_ud_theta_cell_scale_gradient_density_nd(theta; ...,
        oversample=2, scaled_axes=nothing)

Exact frozen-field derivative of the direct UD-theta energy density under a
cell dilation of the selected axes. The same Fourier-interpolated theta
samples, physical sensor filter, spectral derivatives, and adaptive
`K_psi(Q,P)` law as
[`diblock_bvk2_ud_theta_energy_nd`](@ref) are differentiated; no legacy
density-edge quotient enters this stress.
"""
function diblock_bvk2_ud_theta_cell_scale_gradient_density_nd(
        theta_a::AbstractArray{<:Real,D}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, N::Real=1.0, b::Real=1.0, adaptive::Bool=true,
        c2::Real, oversample::Integer=2, scaled_axes=nothing,
        sensor_filter_ratio::Real=
            Inf) where {D}
    theta = Array{Float64,D}(theta_a)
    cell_lengths = _bvk2_oversample_lengths(lengths, D)
    selected = _bvk2_ud_selected_axes(D, scaled_axes)
    fine_theta = if oversample <= 1
        theta
    else
        ws = _bvk2_oversample_ws_cached(size(theta), Int(oversample))
        copy(_bvk2_fourier_upsample!(ws, theta))
    end
    phi = sin.(fine_theta).^2
    physical, dims, ff, _chi, fine_lengths, nn, bb, _dV =
        _diblock_profile_energy_inputs_nd(phi; f=f, chiN=chiN,
            lengths=cell_lengths, N=N, b=b, mean_tolerance=1.0e-4)
    params = diblock_bvk2_parameters(; f=ff, N=nn, b=bb, c2=c2)
    point_count = length(physical)

    values_by_dim = [
        _periodic_k2_1d(dims[axis], fine_lengths[axis])
        for axis in 1:D
    ]
    k2_grid = zeros(Float64, dims)
    selected_k2_grid = zeros(Float64, dims)
    psi_hat = fft(physical .- ff)
    connectivity_scale_sum = 0.0
    for index in CartesianIndices(psi_hat)
        coordinates = Tuple(index)
        k2 = 0.0
        selected_k2 = 0.0
        for axis in 1:D
            directional_k2 = values_by_dim[axis][coordinates[axis]]
            k2 += directional_k2
            selected[axis] && (selected_k2 += directional_k2)
        end
        k2_grid[index] = k2
        selected_k2_grid[index] = selected_k2
        k2 == 0.0 && continue
        connectivity_scale_sum +=
            abs2(psi_hat[index]) * selected_k2 / k2^2
    end
    derivative_density =
        params.A_psi * connectivity_scale_sum / point_count^2

    eta = map(value -> diblock_bvk1_eta_psi(value; f=ff), physical)
    ratio = _bvk2_ud_sensor_filter_ratio(sensor_filter_ratio)
    eta_hat = fft(eta)
    filter_multiplier = ratio == Inf ? ones(Float64, dims) :
        @. 1.0 /
            (1.0 + (k2_grid / (ratio^2 * params.k_star2))^2)
    sensor_eta_hat = eta_hat .* filter_multiplier
    sensor_eta = real.(ifft(sensor_eta_hat))

    # Under a logarithmic dilation of the selected axes, k^2 changes by
    # -2*k_selected^2.  Differentiate the physical filter itself; omitting
    # this term would make the cell stress inconsistent with the energy.
    sensor_eta_scale_hat = if adaptive && ratio != Inf
        cutoff2 = ratio^2 * params.k_star2
        @. eta_hat * 4.0 * filter_multiplier^2 * k2_grid *
            selected_k2_grid / cutoff2^2
    else
        zero.(sensor_eta_hat)
    end
    sensor_eta_scale = real.(ifft(sensor_eta_scale_hat))
    eta_gradients, eta_laplacian =
        _bvk1_spectral_derivatives(sensor_eta, fine_lengths)
    eta_scale_gradients, eta_scale_laplacian =
        _bvk1_spectral_derivatives(sensor_eta_scale, fine_lengths)
    theta_gradients =
        _bvk1_spectral_gradient_components(fine_theta, fine_lengths)
    eta_grad2 = zeros(Float64, dims)
    eta_grad2_scale = zeros(Float64, dims)
    theta_grad2 = zeros(Float64, dims)
    selected_theta_grad2 = zeros(Float64, dims)
    for axis in 1:D
        eta_component2 = eta_gradients[axis].^2
        theta_component2 = theta_gradients[axis].^2
        eta_grad2 .+= eta_component2
        theta_grad2 .+= theta_component2
        eta_component_scale = eta_scale_gradients[axis]
        if selected[axis]
            eta_component_scale = eta_component_scale .- eta_gradients[axis]
            selected_theta_grad2 .+= theta_component2
        end
        eta_grad2_scale .+=
            2.0 .* eta_gradients[axis] .* eta_component_scale
    end

    selected_laplacian_hat =
        @. -selected_k2_grid * sensor_eta_hat
    selected_eta_laplacian = real.(ifft(selected_laplacian_hat))
    eta_laplacian_scale =
        eta_scale_laplacian .- 2.0 .* selected_eta_laplacian

    eta_radius2 = sensor_eta.^2 .+ params.eta_star2
    eta_radius2_scale = 2.0 .* sensor_eta .* sensor_eta_scale
    reduced_q = eta_laplacian ./
        (params.k_star2 .* sqrt.(eta_radius2))
    reduced_p = eta_grad2 ./ (params.k_star2 .* eta_radius2)
    reduced_q_scale = eta_laplacian_scale ./
        (params.k_star2 .* sqrt.(eta_radius2)) .-
        0.5 .* reduced_q .* eta_radius2_scale ./ eta_radius2
    reduced_p_scale = eta_grad2_scale ./
        (params.k_star2 .* eta_radius2) .-
        reduced_p .* eta_radius2_scale ./ eta_radius2
    curvature = reduced_q.^2 ./ sqrt.(1.0 .+ reduced_q.^2)
    sensor = params.c2 .* (2.0 .* curvature .+ reduced_p)
    kpsi = adaptive ?
        params.K_psi0 .* (1.0 .+ 0.5 .* tanh.(sensor)) :
        fill(params.K_psi0, dims)
    kpsi_scale_derivative = zeros(Float64, dims)
    if adaptive
        curvature_prime = reduced_q .* (2.0 .+ reduced_q.^2) ./
            (1.0 .+ reduced_q.^2).^(1.5)
        sensor_scale_derivative = params.c2 .*
            (2.0 .* curvature_prime .* reduced_q_scale .+
             reduced_p_scale)
        kpsi_scale_derivative .=
            0.5 .* params.K_psi0 .* (1.0 .- tanh.(sensor).^2) .*
            sensor_scale_derivative
    end

    derivative_density += mean(
        4.0 .* kpsi_scale_derivative .* theta_grad2 .-
        8.0 .* kpsi .* selected_theta_grad2)
    return derivative_density
end

"""Fine-quadrature composition and its exact coarse-theta normal."""
function _bvk2_ud_theta_composition_normal(
        theta::AbstractArray{Float64,D}; oversample::Integer=2) where {D}
    if oversample <= 1
        fine_theta = theta
        normal = sin.(2.0 .* theta) ./ length(theta)
        return mean(sin.(theta).^2), normal
    end
    ws = _bvk2_oversample_ws_cached(size(theta), Int(oversample))
    fine_theta = _bvk2_fourier_upsample!(ws, theta)
    fine_normal = sin.(2.0 .* fine_theta)
    normal = copy(_bvk2_fourier_upsample_adjoint!(ws, fine_normal))
    normal ./= length(fine_theta)
    return mean(sin.(fine_theta).^2), normal
end

"""Fill a reusable coarse-grid composition normal and return the composition."""
function _bvk2_ud_theta_composition_normal!(
        normal::AbstractArray{Float64,D},
        fine_normal::AbstractArray{Float64,D},
        theta::AbstractArray{Float64,D}; oversample::Integer=2) where {D}
    size(normal) == size(theta) ||
        throw(DimensionMismatch("composition-normal buffer size mismatch"))
    factor = Int(oversample)
    if factor <= 1
        inverse_count = inv(Float64(length(theta)))
        composition_sum = 0.0
        @inbounds @simd for index in eachindex(theta, normal)
            angle = theta[index]
            sine = sin(angle)
            composition_sum += sine * sine
            normal[index] = sin(2.0 * angle) * inverse_count
        end
        return composition_sum * inverse_count
    end
    ws = _bvk2_oversample_ws_cached(size(theta), factor)
    size(fine_normal) == ws.fdims ||
        throw(DimensionMismatch("fine composition-normal buffer size mismatch"))
    fine_theta = _bvk2_fourier_upsample!(ws, theta)
    composition_sum = 0.0
    @inbounds @simd for index in eachindex(fine_theta, fine_normal)
        angle = fine_theta[index]
        sine = sin(angle)
        composition_sum += sine * sine
        fine_normal[index] = sin(2.0 * angle)
    end
    coarse_normal = _bvk2_fourier_upsample_adjoint!(ws, fine_normal)
    inverse_count = inv(Float64(length(fine_theta)))
    @inbounds @simd for index in eachindex(normal, coarse_normal)
        normal[index] = coarse_normal[index] * inverse_count
    end
    return composition_sum * inverse_count
end

function _bvk2_ud_theta_mean_retract_buffered!(
        theta::AbstractArray{Float64,D}, f::Real,
        normal::AbstractArray{Float64,D},
        fine_normal::AbstractArray{Float64,D};
        oversample::Integer=2, tolerance::Real=1.0e-12,
        max_iterations::Integer=60) where {D}
    target = Float64(f)
    0.0 < target < 1.0 ||
        throw(ArgumentError("target composition must lie in (0, 1)"))
    for _ in 1:max(1, Int(max_iterations))
        composition = _bvk2_ud_theta_composition_normal!(
            normal, fine_normal, theta; oversample=oversample)
        residual = composition - target
        abs(residual) <= Float64(tolerance) && return theta
        slope = sum(abs2, normal)
        slope > 1.0e-20 || throw(ArgumentError(
            "UD-theta composition is frozen because the field has no interface"))
        coefficient = -residual / slope
        largest_normal = maximum(abs, normal)
        max_coefficient = 0.5 / max(largest_normal, eps(Float64))
        coefficient = clamp(coefficient, -max_coefficient, max_coefficient)
        @inbounds @simd for index in eachindex(theta, normal)
            theta[index] += coefficient * normal[index]
        end
    end
    composition = _bvk2_ud_theta_composition_normal!(
        normal, fine_normal, theta; oversample=oversample)
    throw(ArgumentError(
        "UD-theta fine-quadrature mean retraction did not converge: " *
        "error=$(composition-target)"))
end

"""
    _bvk2_ud_theta_mean_retract!(theta, f; oversample=2)

Restore the target composition measured by the same fine quadrature used by
the UD-theta objective.  Newton steps follow the exact coarse-theta constraint
normal, so saturated bulk nodes remain effectively fixed while interface
degrees of freedom carry the correction.
"""
function _bvk2_ud_theta_mean_retract!(
        theta::AbstractArray{Float64,D}, f::Real;
        oversample::Integer=2, tolerance::Real=1.0e-12,
        max_iterations::Integer=60) where {D}
    target = Float64(f)
    0.0 < target < 1.0 ||
        throw(ArgumentError("target composition must lie in (0, 1)"))
    for _ in 1:max(1, Int(max_iterations))
        composition, normal = _bvk2_ud_theta_composition_normal(
            theta; oversample=oversample)
        residual = composition - target
        abs(residual) <= Float64(tolerance) && return theta
        slope = sum(abs2, normal)
        slope > 1.0e-20 || throw(ArgumentError(
            "UD-theta composition is frozen because the field has no interface"))
        coefficient = -residual / slope
        largest_normal = maximum(abs, normal)
        max_coefficient = 0.5 / max(largest_normal, eps(Float64))
        coefficient = clamp(coefficient, -max_coefficient, max_coefficient)
        theta .+= coefficient .* normal
    end
    composition, _ = _bvk2_ud_theta_composition_normal(
        theta; oversample=oversample)
    throw(ArgumentError(
        "UD-theta fine-quadrature mean retraction did not converge: " *
        "error=$(composition-target)"))
end

"""Composition manifold shared by every fixed-cell UD-theta morphology.

The manifold depends only on the target composition and grid shape.  It has no
morphology-specific logic, so the same implementation is used for LAM, CYL,
BCC, GYR, FCC, O70, and any future periodic seed represented in the UD-theta
chart.
"""
struct DiblockBVK2UDThetaMean{D} <: Optim.Manifold
    f::Float64
    constraint_oversample::Int
    dims::NTuple{D,Int}
    normal::Array{Float64,D}
    fine_normal::Array{Float64,D}
end

function DiblockBVK2UDThetaMean(f::Real, constraint_oversample::Integer,
        dims::NTuple{D,Int}) where {D}
    factor = Int(constraint_oversample)
    fine_dims = factor <= 1 ? dims : ntuple(axis -> factor * dims[axis], D)
    return DiblockBVK2UDThetaMean{D}(Float64(f), factor, dims,
        Array{Float64,D}(undef, dims),
        Array{Float64,D}(undef, fine_dims))
end

function Optim.retract!(manifold::DiblockBVK2UDThetaMean, theta)
    _bvk2_ud_theta_mean_retract_buffered!(
        reshape(theta, manifold.dims), manifold.f,
        manifold.normal, manifold.fine_normal;
        oversample=manifold.constraint_oversample)
    return theta
end

function Optim.project_tangent!(manifold::DiblockBVK2UDThetaMean,
    gradient, theta)
    field = reshape(theta, manifold.dims)
    _bvk2_ud_theta_composition_normal!(manifold.normal,
        manifold.fine_normal, field;
        oversample=manifold.constraint_oversample)
    normal2 = sum(abs2, manifold.normal)
    normal2 > eps(Float64) ||
        throw(ArgumentError("UD-theta composition manifold is frozen"))
    normal_dot_gradient = 0.0
    @inbounds @simd for index in eachindex(manifold.normal, gradient)
        normal_dot_gradient += manifold.normal[index] * gradient[index]
    end
    coefficient = normal_dot_gradient / normal2
    @inbounds @simd for index in eachindex(manifold.normal, gradient)
        gradient[index] -= coefficient * manifold.normal[index]
    end
    return gradient
end

function _bvk2_ud_theta_energy_or_infinity(evaluate::Function)
    try
        return evaluate()
    catch exception
        message = sprint(showerror, exception)
        if exception isa ArgumentError && occursin(
                "mean phi_A must match f within mean_tolerance", message)
            return Inf
        end
        rethrow()
    end
end

"""Lean phase-neutral value/gradient callback for fixed-cell UD-theta solves.

Energy-only line-search calls evaluate only the energy.  Gradient calls use the
fused analytic value/gradient evaluator, while mean retraction and tangent
projection are handled by [`DiblockBVK2UDThetaMean`](@ref).  Expensive stress,
topology, and identity diagnostics therefore remain outside the optimization
hot loop.
"""
mutable struct DiblockBVK2UDThetaHotLoop{D,L,C}
    dims::NTuple{D,Int}
    f::Float64
    chiN::Float64
    lengths::L
    N::Float64
    b::Float64
    adaptive::Bool
    c2::Float64
    oversample::Int
    sensor_filter_ratio::Float64
    on_gradient::C
    calls::Int
    energy_calls::Int
    gradient_calls::Int
    last_energy::Float64
end

function DiblockBVK2UDThetaHotLoop(theta::AbstractArray{<:Real,D};
        f::Real=0.5, chiN::Real=12.0, lengths=1.0, N::Real=1.0,
        b::Real=1.0, adaptive::Bool=true, c2::Real,
        oversample::Integer=2, sensor_filter_ratio::Real=Inf,
        on_gradient=nothing) where {D}
    cell_lengths = _bvk2_oversample_lengths(lengths, D)
    return DiblockBVK2UDThetaHotLoop{D,typeof(cell_lengths),typeof(on_gradient)}(
        size(theta), Float64(f), Float64(chiN), cell_lengths, Float64(N),
        Float64(b), adaptive, Float64(c2), Int(oversample),
        _bvk2_ud_sensor_filter_ratio(sensor_filter_ratio), on_gradient,
        0, 0, 0, NaN)
end

function (hotloop::DiblockBVK2UDThetaHotLoop)(F, G, values)
    hotloop.calls += 1
    theta = reshape(values, hotloop.dims)
    energy = NaN
    gradient = nothing
    if G === nothing
        if F !== nothing
            energy = _bvk2_ud_theta_energy_or_infinity() do
                diblock_bvk2_ud_theta_energy_nd(theta;
                    f=hotloop.f, chiN=hotloop.chiN,
                    lengths=hotloop.lengths, N=hotloop.N, b=hotloop.b,
                    adaptive=hotloop.adaptive, c2=hotloop.c2,
                    oversample=hotloop.oversample,
                    sensor_filter_ratio=hotloop.sensor_filter_ratio)
            end
        end
    else
        fused = diblock_bvk2_ud_theta_value_gradient_nd(theta;
            f=hotloop.f, chiN=hotloop.chiN, lengths=hotloop.lengths,
            N=hotloop.N, b=hotloop.b, adaptive=hotloop.adaptive,
            c2=hotloop.c2, oversample=hotloop.oversample,
            project_mean=false,
            sensor_filter_ratio=hotloop.sensor_filter_ratio)
        energy = fused.energy
        gradient = fused.gradient
        G .= vec(gradient)
        hotloop.gradient_calls += 1
        hotloop.on_gradient === nothing || hotloop.on_gradient((
            energy=energy,
            gradient=gradient,
            theta=theta,
            gradient_call=hotloop.gradient_calls,
        ))
    end
    if F !== nothing
        hotloop.energy_calls += 1
        isfinite(energy) && (hotloop.last_energy = energy)
    end
    return F === nothing ? nothing : energy
end

"""
    optimize_diblock_bvk2_ud_theta_hotloop(theta0; kwargs...)

Run the shared fixed-cell L-BFGS hot loop for an arbitrary-dimensional
UD-theta morphology.  The result includes the retracted final field, the raw
`Optim` result, callback counters, elapsed time, and a stable stop label.
"""
function optimize_diblock_bvk2_ud_theta_hotloop(
        theta0::AbstractArray{<:Real,D}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, N::Real=1.0, b::Real=1.0, adaptive::Bool=true,
        c2::Real, oversample::Integer=2, constraint_oversample::Integer=1,
        sensor_filter_ratio::Real=Inf, max_iterations::Integer=2000,
        time_limit::Real=3600.0, lbfgs_memory::Integer=20,
        energy_tolerance::Real=1.0e-8, gradient_tolerance::Real=1.0e-12,
        on_gradient=nothing) where {D}
    theta = Array{Float64,D}(theta0)
    _bvk2_ud_theta_mean_retract!(theta, f;
        oversample=constraint_oversample)
    hotloop = DiblockBVK2UDThetaHotLoop(theta; f=f, chiN=chiN,
        lengths=lengths, N=N, b=b, adaptive=adaptive, c2=c2,
        oversample=oversample, sensor_filter_ratio=sensor_filter_ratio,
        on_gradient=on_gradient)
    manifold = DiblockBVK2UDThetaMean(
        Float64(f), Int(constraint_oversample), size(theta))
    started = time()
    optimization = Optim.optimize(
        Optim.NLSolversBase.only_fg!(hotloop), vec(theta),
        Optim.LBFGS(m=Int(lbfgs_memory), manifold=manifold,
            linesearch=Optim.LineSearches.BackTracking(order=3)),
        Optim.Options(iterations=Int(max_iterations),
            g_tol=Float64(gradient_tolerance),
            f_reltol=Float64(energy_tolerance),
            time_limit=Float64(time_limit), show_trace=false,
            store_trace=false, extended_trace=false))
    elapsed = time() - started
    final_theta = reshape(copy(Vector{Float64}(Optim.minimizer(optimization))),
        size(theta))
    _bvk2_ud_theta_mean_retract!(final_theta, f;
        oversample=constraint_oversample)
    stop_reason = Optim.converged(optimization) ? "hotloop_converged" :
        Optim.iteration_limit_reached(optimization) ?
        "hotloop_iteration_cap" :
        elapsed >= 0.95 * Float64(time_limit) ? "hotloop_time_limit" :
        "hotloop_stopped"
    return (
        theta=final_theta,
        optimization=optimization,
        iterations=Optim.iterations(optimization),
        converged=Optim.converged(optimization),
        calls=hotloop.calls,
        energy_calls=hotloop.energy_calls,
        gradient_calls=hotloop.gradient_calls,
        elapsed_seconds=elapsed,
        stop_reason=stop_reason,
    )
end

"""
    minimize_diblock_bvk2_ud_theta_lamella(; f, chiN, L, c2, ...)

Relax a one-period lamellar field at fixed period with the oversampled
UD-`theta` BVK2 functional.  A finite `sensor_filter_ratio` applies the same
physical sensor filter in the energy, exact field gradient, and subsequent
cell-stress calculation.  Convergence is certified from the projected
coarse-`theta` force after optimization rather than from the optimizer's
internal stopping flag.
"""
function minimize_diblock_bvk2_ud_theta_lamella(; f::Real=0.5,
        chiN::Real=12.0, L=nothing, nx::Integer=128,
        mode_count::Integer=3, initial_amplitude::Real=0.2,
        initial_profile=nothing, N::Real=1.0, b::Real=1.0,
        adaptive::Bool=true, c2::Real, oversample::Integer=2,
        sensor_filter_ratio::Real=DIBLOCK_BVK2_UD_SENSOR_FILTER_RATIO,
        max_iterations::Integer=4000, lbfgs_memory::Integer=20,
        energy_tolerance::Real=1.0e-10,
        force_tolerance::Real=1.0e-5,
        force_maxabs_tolerance::Real=2.0e-5,
        polish_iterations::Integer=8,
        polish_difference_step::Real=1.0e-6,
        polish_trust_radius::Real=0.1,
        time_limit::Real=1800.0)
    count = Int(nx)
    count >= 16 && iseven(count) || throw(ArgumentError(
        "UD-theta lamellar grid size must be even and at least 16"))
    modes = Int(mode_count)
    modes >= 1 || throw(ArgumentError("mode_count must be positive"))
    factor = Int(oversample)
    factor >= 1 || throw(ArgumentError("oversample must be positive"))
    iterations = Int(max_iterations)
    iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    force_tol = Float64(force_tolerance)
    force_max_tol = Float64(force_maxabs_tolerance)
    force_tol > 0.0 && isfinite(force_tol) || throw(ArgumentError(
        "force_tolerance must be positive and finite"))
    force_max_tol > 0.0 && isfinite(force_max_tol) || throw(ArgumentError(
        "force_maxabs_tolerance must be positive and finite"))

    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    chi = Float64(chiN)
    isfinite(chi) || throw(ArgumentError("chiN must be finite"))
    reference_period = 2.0 * pi /
        sqrt(_diblock_kernel_minimum_k2(; f=ff, N=nn, b=bb))
    period = L === nothing ? reference_period : Float64(L)
    period > 0.0 && isfinite(period) || throw(ArgumentError(
        "lamella period L must be positive and finite"))
    ratio = _bvk2_ud_sensor_filter_ratio(sensor_filter_ratio)
    polish_count = Int(polish_iterations)
    polish_count >= 0 || throw(ArgumentError(
        "polish_iterations must be nonnegative"))
    difference_step = Float64(polish_difference_step)
    difference_step > 0.0 && isfinite(difference_step) ||
        throw(ArgumentError(
            "polish_difference_step must be positive and finite"))
    trust_radius = Float64(polish_trust_radius)
    trust_radius > 0.0 && isfinite(trust_radius) ||
        throw(ArgumentError(
            "polish_trust_radius must be positive and finite"))

    phi0 = if initial_profile === nothing
        last(diblock_lamella_profile_from_modes(
            [Float64(initial_amplitude); zeros(Float64, modes - 1)];
            f=ff, nx=count, L=period))
    else
        candidate = Float64.(collect(initial_profile))
        length(candidate) == count || throw(DimensionMismatch(
            "initial_profile length must equal nx"))
        all(isfinite, candidate) || throw(ArgumentError(
            "initial_profile must contain only finite values"))
        clamp.(candidate, 1.0e-12, 1.0 - 1.0e-12)
    end
    theta0 = asin.(sqrt.(phi0))
    _bvk2_ud_theta_mean_retract!(theta0, ff; oversample=factor)
    hotloop = optimize_diblock_bvk2_ud_theta_hotloop(theta0;
        f=ff, chiN=chi, lengths=(period,), N=nn, b=bb,
        adaptive=adaptive, c2=c2, oversample=factor,
        constraint_oversample=factor, sensor_filter_ratio=ratio,
        max_iterations=iterations, time_limit=time_limit,
        lbfgs_memory=lbfgs_memory, energy_tolerance=energy_tolerance,
        gradient_tolerance=1.0e-12)
    theta = vec(hotloop.theta)
    dx = period / count
    completed_polish_iterations = 0
    if polish_count > 0
        seed = _diblock_bvk2_reflection_symmetrize(theta)
        _bvk2_ud_theta_mean_retract!(seed, ff; oversample=factor)
        basis = [cos(2.0 * pi * mode * (index - 0.5) / count)
            for index in 1:count, mode in 1:(count ÷ 2 - 1)]
        anchor = mean(seed)
        coefficients = (2.0 / count) .* (transpose(basis) * seed)
        function cosine_state(values)
            trial_theta = anchor .+ basis * values
            _bvk2_ud_theta_mean_retract!(
                trial_theta, ff; oversample=factor)
            fused = diblock_bvk2_ud_theta_value_gradient_nd(trial_theta;
                f=ff, chiN=chi, lengths=(period,), N=nn, b=bb,
                adaptive=adaptive, c2=c2, oversample=factor,
                project_mean=true, sensor_filter_ratio=ratio)
            force = fused.gradient ./ dx
            residual = (2.0 / count) .* (transpose(basis) * force)
            return (theta=trial_theta, energy=fused.energy, force=force,
                residual=residual, r2=sqrt(mean(abs2, force)),
                rinf=maximum(abs, force))
        end
        state = cosine_state(coefficients)
        for iteration in 1:polish_count
            state.r2 <= force_tol && state.rinf <= force_max_tol && break
            coefficient_count = length(coefficients)
            jacobian = Matrix{Float64}(undef,
                coefficient_count, coefficient_count)
            for column in 1:coefficient_count
                step = difference_step *
                    (1.0 + abs(coefficients[column]))
                plus = copy(coefficients)
                minus = copy(coefficients)
                plus[column] += step
                minus[column] -= step
                jacobian[:, column] .=
                    (cosine_state(plus).residual .-
                     cosine_state(minus).residual) ./ (2.0 * step)
            end
            all(isfinite, jacobian) || break
            direction = try
                -(pinv(jacobian; rtol=1.0e-11) * state.residual)
            catch
                break
            end
            all(isfinite, direction) || break
            physical_step = basis * direction
            largest_step = maximum(abs, physical_step)
            largest_step > trust_radius &&
                (direction .*= trust_radius / largest_step)
            accepted = nothing
            alpha = 1.0
            for _ in 1:24
                trial_coefficients = coefficients .+ alpha .* direction
                trial = cosine_state(trial_coefficients)
                if isfinite(trial.energy) && trial.r2 < state.r2
                    accepted = (coefficients=trial_coefficients,
                        state=trial)
                    break
                end
                alpha *= 0.5
            end
            accepted === nothing && break
            coefficients = accepted.coefficients
            state = accepted.state
            completed_polish_iterations = iteration
        end
        theta = state.theta
    end
    evaluation = diblock_bvk2_ud_theta_value_gradient_nd(theta;
        f=ff, chiN=chi, lengths=(period,), N=nn, b=bb,
        adaptive=adaptive, c2=c2, oversample=factor,
        project_mean=true, sensor_filter_ratio=ratio)
    projected_force = evaluation.gradient ./ dx
    projected_force_norm = sqrt(mean(abs2, projected_force))
    projected_force_maxabs = maximum(abs, projected_force)
    converged = projected_force_norm <= force_tol &&
        projected_force_maxabs <= force_max_tol
    phi = sin.(theta).^2
    composition, _ = _bvk2_ud_theta_composition_normal(
        theta; oversample=factor)
    homogeneous_theta = fill(asin(sqrt(ff)), count)
    homogeneous_energy = diblock_bvk2_ud_theta_energy_nd(
        homogeneous_theta; f=ff, chiN=chi, lengths=(period,), N=nn,
        b=bb, adaptive=adaptive, c2=c2, oversample=factor,
        sensor_filter_ratio=ratio)
    x = [(index - 0.5) * dx for index in 1:count]
    coefficients = _diblock_cosine_coefficients_from_profile(phi; f=ff,
        L=period, mode_count=modes)
    return DiblockBurpLamellaResult(ff, chi, period, count, modes, x,
        phi, coefficients, evaluation.energy, homogeneous_energy,
        converged, hotloop.iterations + completed_polish_iterations,
        minimum(phi), maximum(phi),
        composition, projected_force_norm, projected_force_maxabs)
end

"""
    minimize_diblock_bvk2_ud_theta_lamella_stress_free(; f, chiN, c2, ...)

Root the analytic cell stress of the oversampled, optionally sensor-filtered
UD-`theta` BVK2 functional.  Every field, objective, and cell-stress sample
uses the same `oversample` and `sensor_filter_ratio`; accepted roots also pass
two-sided stress orientation and energy checks.
"""
function minimize_diblock_bvk2_ud_theta_lamella_stress_free(;
        f::Real=0.5, chiN::Real=12.0, nx::Integer=128,
        mode_count::Integer=3, initial_amplitudes=(0.2, 0.5, 1.0),
        initial_profile=nothing, N::Real=1.0, b::Real=1.0,
        adaptive::Bool=true, c2::Real, oversample::Integer=2,
        sensor_filter_ratio::Real=DIBLOCK_BVK2_UD_SENSOR_FILTER_RATIO,
        max_iterations::Integer=4000, max_period_iterations::Integer=24,
        lower_factor::Real=0.75, upper_factor::Real=2.0,
        lbfgs_memory::Integer=20, energy_tolerance::Real=1.0e-10,
        force_tolerance::Real=1.0e-5,
        force_maxabs_tolerance::Real=2.0e-5,
        time_limit::Real=1800.0, bootstrap_period_factors=(),
        bootstrap_window::Real=0.12,
        local_check_fraction::Real=0.01,
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
    reference_period = 2.0 * pi /
        sqrt(_diblock_kernel_minimum_k2(; f=ff, N=nn, b=bb))
    rg_scale = sqrt(6.0 / (nn * bb^2))
    ratio = _bvk2_ud_sensor_filter_ratio(sensor_filter_ratio)
    factor_os = Int(oversample)
    factor_os >= 1 || throw(ArgumentError("oversample must be positive"))
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
        for (profile, amplitude) in attempts
            result = try
                minimize_diblock_bvk2_ud_theta_lamella(; f=ff, chiN=chi,
                    L=period, nx=nx, mode_count=mode_count,
                    initial_amplitude=amplitude, initial_profile=profile,
                    N=nn, b=bb, adaptive=adaptive, c2=c2,
                    oversample=factor_os, sensor_filter_ratio=ratio,
                    max_iterations=max_iterations,
                    lbfgs_memory=lbfgs_memory,
                    energy_tolerance=energy_tolerance,
                    force_tolerance=force_tolerance,
                    force_maxabs_tolerance=force_maxabs_tolerance,
                    time_limit=time_limit)
            catch
                nothing
            end
            result === nothing && continue
            identity_pass = result.energy <
                result.homogeneous_energy - 1.0e-10 &&
                result.maximum_phi - result.minimum_phi > 1.0e-4 &&
                _diblock_one_period_lamellar_profile(result.phi_a)
            valid = result.converged && identity_pass &&
                isfinite(result.energy)
            stress = valid ?
                diblock_bvk2_ud_theta_cell_scale_gradient_density_nd(
                    asin.(sqrt.(clamp.(result.phi_a, 0.0, 1.0)));
                    f=ff, chiN=chi, lengths=(period,), N=nn, b=bb,
                    adaptive=adaptive, c2=c2, oversample=factor_os,
                    sensor_filter_ratio=ratio) : NaN
            valid &= isfinite(stress)
            row = (factor=value, log_factor=log(value), period=period,
                objective=result.energy / period, stress=stress,
                result=result, amplitude=amplitude, valid=valid,
                identity_pass=identity_pass)
            if best === nothing || (row.valid && !best.valid) ||
                    (row.valid == best.valid &&
                     row.objective < best.objective)
                best = row
            end
            row.valid && break
        end
        best === nothing && throw(ErrorException(@sprintf(
            "all filtered BVK2 fixed-period attempts failed at factor %.8g",
            value)))
        evaluations[value] = best
        _period_progress(progress_label, @sprintf(
            "filtered BVK2 cell root factor=%.8g L/Rg=%.8g stress=%.8g field=%s force=%.4g maxforce=%.4g identity=%s",
            value, period * rg_scale, best.stress,
            string(best.result.converged),
            best.result.projected_force_norm,
            best.result.projected_force_maxabs,
            string(best.identity_pass)))
        return best
    end

    center = evaluate_factor(center_factor)
    center.valid || throw(ErrorException(
        "filtered BVK2 bootstrap period did not produce a stationary lamella"))
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
        "filtered BVK2 analytic cell stress has no valid oriented bracket"))
    lower, upper = bracket
    selected = abs(lower.stress) <= abs(upper.stress) ? lower : upper
    root_steps = 0
    for iteration in 1:period_iterations
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
        root_steps = iteration
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

    return DiblockLamellaPeriodOptimizationResult("BVK2", selected.result,
        reference_period, selected.factor, selected.period,
        selected.period * rg_scale, selected.objective,
        "FilteredUDThetaAnalyticCellStressRoot", local_pass,
        length(evaluations), selected.amplitude, lo, hi,
        !isempty(bootstrap_factors), center_factor, window,
        selected.factor, selected.objective,
        left.factor, left.objective, right.factor, right.objective,
        check_fraction, local_pass, boundary_limited)
end

# Project out the mean-changing direction w = sin(2θ): with dφ = sin(2θ) dθ,
# the mean constraint d⟨φ⟩ = 0 is ⟨w, dθ⟩ = 0, so remove the w-component of the
# gradient.  Matches the θ-BB solver's λ-projection.
function _bvk2_project_mean_direction!(grad::AbstractArray{Float64},
        theta::AbstractArray{Float64})
    ws_ = 0.0; wg = 0.0
    @inbounds for i in eachindex(grad)
        w = sin(2.0 * theta[i])
        ws_ += w * w
        wg += w * grad[i]
    end
    lambda = wg / max(ws_, eps(Float64))
    @inbounds for i in eachindex(grad)
        grad[i] -= lambda * sin(2.0 * theta[i])
    end
    return grad
end

_bvk2_oversample_lengths(lengths, rank::Int) =
    lengths isa Number ? fill(Float64(lengths), rank) : collect(Float64.(lengths))
