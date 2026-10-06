const DIBLOCK_BVK2_FORMULA_SCHEMA = "bvk2-rpa-eta-v1"
const DIBLOCK_BVK2_DISCRETIZATION_SCHEMA =
    "central_eta_forward_edge_exact_v1"
const DIBLOCK_BVK2_SOLVER_SCHEMA = "mean_logit_cosine_lbfgs_residual_newton_v2"

"""
    diblock_bvk2_parameters(; f=0.5, N=1.0, b=1.0, c2)

Return the continuum parameters of the differentiable BVK2 constitutive law.
The intrinsic ordering scale is fixed by the Uneyama--Doi Gaussian response,
`k_star^4 = 12f(1-f)A_psi/(N*b^2)`; `c2` is the sole free constitutive number.
"""
function _bvk2_parameters(base, model::Symbol, architecture::Symbol;
        f::Real, N::Real, b::Real, c2::Real)
    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    coefficient = Float64(c2)
    coefficient >= 0.0 && isfinite(coefficient) ||
        throw(ArgumentError("c2 must be nonnegative and finite"))
    eta_star2 = ff * (1.0 - ff)
    k_star2 = sqrt(12.0 * eta_star2 * base.A_psi / (nn * bb^2))
    return (
        model=model,
        architecture=architecture,
        formula_schema=DIBLOCK_BVK2_FORMULA_SCHEMA,
        coordinate=:eta_psi,
        coordinate_units=:lengths_and_b_same_units,
        f=ff,
        N=nn,
        b=bb,
        c2=coefficient,
        A_psi=base.A_psi,
        M_psi=base.M_psi,
        K_psi0=base.K_psi0,
        eta_star2=eta_star2,
        k_star2=k_star2,
        k_star=sqrt(k_star2),
        k_star_definition="k_star^4=12*f*(1-f)*A_psi/(N*b^2)",
    )
end

function diblock_bvk2_parameters(; f::Real=0.5, N::Real=1.0,
        b::Real=1.0, c2::Real)
    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    base = diblock_bvk1_coefficients(; f=ff, N=nn, b=bb)
    return _bvk2_parameters(base, :BVK2, :AB;
        f=ff, N=nn, b=bb, c2=c2)
end

"""
    triblock_aba_bvk2_parameters(; f=0.5, N=1.0, b=1.0, c2)

Return BVK2 parameters for a symmetric linear `A_(f/2)-B_(1-f)-A_(f/2)`
triblock. Only the analytically derived ideal-chain coefficients change from
the diblock model; the nonlinear closure coefficient `c2` is passed through
unchanged.
"""
function triblock_aba_bvk2_parameters(; f::Real=0.5, N::Real=1.0,
        b::Real=1.0, c2::Real)
    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    base = triblock_aba_bvk1_coefficients(; f=ff, N=nn, b=bb)
    return _bvk2_parameters(base, :BVK2, :ABA;
        f=ff, N=nn, b=bb, c2=c2)
end

"""
    diblock_bvk2_fingerprint(; f=0.5, N=1.0, b=1.0, c2)

Return a stable, human-readable model fingerprint.  It records every choice
that affects the BVK2 constitutive law and deliberately identifies BVK2
separately from the frozen BVK1 model.
"""
function _bvk2_fingerprint(params; adaptive::Bool=true)
    law_mode = adaptive ? "adaptive" : "nonadaptive_control"
    model_identity = params.architecture == :AB ?
        "model=BVK2" : "model=BVK2|architecture=$(params.architecture)"
    return @sprintf(
        "%s|schema=%s|law_mode=%s|discretization=%s|coordinate=eta_psi|units=lengths_and_b_same_units|c2=%.17g|f=%.17g|N=%.17g|b=%.17g|k_star_definition=%s|k_star=%.17g",
        model_identity, params.formula_schema, law_mode,
        DIBLOCK_BVK2_DISCRETIZATION_SCHEMA,
        params.c2, params.f, params.N, params.b, params.k_star_definition,
        params.k_star)
end

function diblock_bvk2_fingerprint(; adaptive::Bool=true, kwargs...)
    return _bvk2_fingerprint(diblock_bvk2_parameters(; kwargs...);
        adaptive=adaptive)
end

"""
    triblock_aba_bvk2_fingerprint(; adaptive=true, f=0.5, N=1, b=1, c2)

Fingerprint the no-refit ABA BVK2 model, including the architecture label and
the unchanged nonlinear closure coefficient.
"""
function triblock_aba_bvk2_fingerprint(; adaptive::Bool=true, kwargs...)
    return _bvk2_fingerprint(triblock_aba_bvk2_parameters(; kwargs...);
        adaptive=adaptive)
end

@inline function _diblock_bvk2_reduced_sensor(eta::Float64, grad2::Float64,
        laplacian::Float64, params)
    eta_radius2 = eta^2 + params.eta_star2
    q = laplacian / (params.k_star2 * sqrt(eta_radius2))
    p = grad2 / (params.k_star2 * eta_radius2)
    analytic_curvature = q^2 / sqrt(1.0 + q^2)
    return params.c2 * (2.0 * analytic_curvature + p)
end

function _diblock_bvk2_kpsi_1d(phi::AbstractVector{Float64}, f::Float64,
        period::Float64, N::Float64, b::Float64; adaptive::Bool=true,
        c2::Real, parameters=nothing)
    count = length(phi)
    params = parameters === nothing ?
        diblock_bvk2_parameters(; f=f, N=N, b=b, c2=c2) : parameters
    !adaptive && return fill(params.K_psi0, count)
    dx = period / count
    eta = [diblock_bvk1_eta_psi(value; f=f) for value in phi]
    out = Vector{Float64}(undef, count)
    for idx in eachindex(phi)
        left = idx == 1 ? count : idx - 1
        right = idx == count ? 1 : idx + 1
        derivative = (eta[right] - eta[left]) / (2.0 * dx)
        laplacian = (eta[right] - 2.0 * eta[idx] + eta[left]) / dx^2
        sensor = _diblock_bvk2_reduced_sensor(
            eta[idx], derivative^2, laplacian, params)
        out[idx] = params.K_psi0 * (1.0 + 0.5 * tanh(sensor))
    end
    return out
end

function _diblock_bvk2_kpsi_nd(phi::AbstractArray{Float64}, f::Float64,
        lengths::AbstractVector{Float64}, N::Float64, b::Float64;
        adaptive::Bool=true, c2::Real)
    dims = size(phi)
    length(lengths) == length(dims) ||
        throw(ArgumentError("lengths must have one entry per array dimension"))
    params = diblock_bvk2_parameters(; f=f, N=N, b=b, c2=c2)
    !adaptive && return fill(params.K_psi0, dims)
    spacings = [lengths[dim] / dims[dim] for dim in eachindex(dims)]
    eta = similar(phi)
    tforeach(CartesianIndices(phi)) do index
        @inbounds eta[index] = diblock_bvk1_eta_psi(phi[index]; f=f)
    end
    out = Array{Float64}(undef, dims)
    # Per-node write, no cross-node dependence (reads neighbour eta only).
    tforeach(CartesianIndices(phi)) do index
        grad2 = 0.0
        laplacian = 0.0
        @inbounds for dim in eachindex(dims)
            left = _periodic_shift_index(index, dims, dim, -1)
            right = _periodic_shift_index(index, dims, dim, 1)
            spacing = spacings[dim]
            derivative = (eta[right] - eta[left]) / (2.0 * spacing)
            grad2 += derivative^2
            laplacian += (eta[right] - 2.0 * eta[index] + eta[left]) / spacing^2
        end
        @inbounds sensor = _diblock_bvk2_reduced_sensor(
            eta[index], grad2, laplacian, params)
        @inbounds out[index] = params.K_psi0 * (1.0 + 0.5 * tanh(sensor))
    end
    return out
end

"""
    diblock_bvk2_kpsi_profile_1d(phi_a; f, L, N=1, b=1, adaptive=true, c2)

Evaluate the analytic RPA-scaled BVK2 stiffness on a periodic one-dimensional
interior density field.
"""
function diblock_bvk2_kpsi_profile_1d(phi_a::AbstractVector{<:Real};
        f::Real=0.5, L::Real=1.0, N::Real=1.0, b::Real=1.0,
        adaptive::Bool=true, c2::Real, mean_tolerance::Real=1.0e-8)
    phi, _count, ff, _chi, period, nn, bb = _diblock_profile_energy_inputs(phi_a;
        f=f, chiN=0.0, L=L, N=N, b=b, mean_tolerance=mean_tolerance)
    return _diblock_bvk2_kpsi_1d(phi, ff, period, nn, bb;
        adaptive=adaptive, c2=c2)
end

"""
    diblock_bvk2_kpsi_profile_nd(phi_a; f, lengths, N=1, b=1, adaptive=true, c2)

Evaluate the analytic RPA-scaled BVK2 stiffness on a periodic N-dimensional
interior density field.
"""
function diblock_bvk2_kpsi_profile_nd(phi_a::AbstractArray{<:Real};
        f::Real=0.5, lengths=1.0, N::Real=1.0, b::Real=1.0,
        adaptive::Bool=true, c2::Real, mean_tolerance::Real=1.0e-8)
    phi, _dims, ff, _chi, cell_lengths, nn, bb, _dV =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=0.0,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    return _diblock_bvk2_kpsi_nd(phi, ff, cell_lengths, nn, bb;
        adaptive=adaptive, c2=c2)
end

function _diblock_bvk2_gradient_integral_1d(phi::AbstractVector{Float64},
        kpsi::AbstractVector{Float64}, dx::Float64)
    count = length(phi)
    length(kpsi) == count || throw(DimensionMismatch("K_psi length must match phi"))
    total = 0.0
    for idx in eachindex(phi)
        right = idx == count ? 1 : idx + 1
        edge_phi = 0.5 * (phi[idx] + phi[right])
        denominator = edge_phi * (1.0 - edge_phi)
        gradient = (phi[right] - phi[idx]) / dx
        edge_k = 0.5 * (kpsi[idx] + kpsi[right])
        total += dx * edge_k * gradient^2 / denominator
    end
    return total
end

function _diblock_bvk2_gradient_integral_nd(phi::AbstractArray{Float64},
        kpsi::AbstractArray{Float64}, lengths::AbstractVector{Float64},
        dV::Float64)
    dims = size(phi)
    size(kpsi) == dims || throw(DimensionMismatch("K_psi size must match phi"))
    length(lengths) == length(dims) ||
        throw(ArgumentError("lengths must have one entry per array dimension"))
    spacings = [lengths[dim] / dims[dim] for dim in eachindex(dims)]
    # Threaded reduction over the low endpoint of every edge (each edge is
    # counted once, exactly as the serial sweep).  Reproduces the serial
    # total to floating-point roundoff (summation order differs).
    total = tmapreduce(+, CartesianIndices(phi)) do index
        acc = 0.0
        @inbounds for dim in eachindex(dims)
            right = _periodic_shift_index(index, dims, dim, 1)
            edge_phi = 0.5 * (phi[index] + phi[right])
            denominator = edge_phi * (1.0 - edge_phi)
            gradient = (phi[right] - phi[index]) / spacings[dim]
            edge_k = 0.5 * (kpsi[index] + kpsi[right])
            acc += dV * edge_k * gradient^2 / denominator
        end
        acc
    end
    return total
end

@inline function _diblock_bvk2_eta_prime(phi::Float64, f::Float64)
    # Floor into the open box so d(eta)/d(phi) stays finite at saturated cells.
    p = clamp(phi, 1.0e-10, 1.0 - 1.0e-10)
    return 0.5 * (f / p + (1.0 - f) / (1.0 - p))
end

# Reusable scratch for the chemical-potential reverse pass.  The JFNK inner
# loop re-enters `_diblock_bvk2_gradient_chemical_potential_nd` thousands of
# times at fixed dims; allocating ~16 grid arrays per call dominates memory
# traffic (~280 MB/call at 128^3) and stalls the threaded loops on GC.  A
# per-dims workspace is cached and reused so the hot path allocates nothing
# but the returned `chemical` array.  Every scratch field is fully overwritten
# each call (`=`, never stale `+=`), so no zeroing is required.
struct _BVK2ChemPotWorkspace{D}
    eta::Array{Float64,D}
    eta_prime::Array{Float64,D}
    gradient_components::Vector{Array{Float64,D}}
    eta_radius2::Array{Float64,D}
    reduced_q::Array{Float64,D}
    reduced_p::Array{Float64,D}
    sensor::Array{Float64,D}
    kpsi::Array{Float64,D}
    k_adjoint::Array{Float64,D}
    eta_adjoint::Array{Float64,D}
    laplacian_adjoint::Array{Float64,D}
    gradient_component_adjoint::Vector{Array{Float64,D}}
end

function _bvk2_chempot_workspace(dims::NTuple{D,Int}) where {D}
    z() = Array{Float64,D}(undef, dims)
    return _BVK2ChemPotWorkspace{D}(
        z(), z(), [z() for _ in 1:D], z(), z(), z(), z(), z(),
        z(), z(), z(), [z() for _ in 1:D])
end

# Per-process cache keyed by grid dims.  Safe without a lock: the morphology
# solver drives chemical-potential calls sequentially on its own thread, and
# the in-call OhMyThreads regions never re-enter this function.  Campaign
# parallelism is process-level (separate julia workers), not shared-thread.
const _BVK2_CHEMPOT_WS_CACHE = Dict{Tuple{Vararg{Int}},Any}()

function _bvk2_chempot_workspace_cached(dims::NTuple{D,Int}) where {D}
    return get!(() -> _bvk2_chempot_workspace(dims),
        _BVK2_CHEMPOT_WS_CACHE, dims)::_BVK2ChemPotWorkspace{D}
end

# Reusable spectral scratch for the RPA inverse-Laplacian term.  Each
# chemical-potential call previously did an out-of-place `fft`/`ifft` round
# trip -- ~4 fresh complex arrays (~130 MB/call at 128^3) plus a re-planned
# transform every time.  Cache one in-place plan pair, one complex buffer,
# and the k^2 grid (rebuilt only when the cell changes), so the term reuses
# them and the FFT runs on the FFTW thread pool (`B2B_FFTW_THREADS`).
mutable struct _BVK2SpectralWorkspace{D}
    lengths::NTuple{D,Float64}
    k2::Array{Float64,D}
    buffer::Array{ComplexF64,D}
    forward::Any
    inverse::Any
end

function _bvk2_spectral_workspace(dims::NTuple{D,Int}, lengths) where {D}
    buffer = zeros(ComplexF64, dims)
    forward = plan_fft!(buffer)
    inverse = plan_ifft!(buffer)
    len = NTuple{D,Float64}(Float64(l) for l in lengths)
    k2 = _periodic_k2_nd(dims, collect(Float64, lengths))
    return _BVK2SpectralWorkspace{D}(len, k2, buffer, forward, inverse)
end

const _BVK2_SPECTRAL_WS_CACHE = Dict{Tuple{Vararg{Int}},Any}()

function _bvk2_spectral_workspace_cached(dims::NTuple{D,Int}, lengths) where {D}
    ws = get!(() -> _bvk2_spectral_workspace(dims, lengths),
        _BVK2_SPECTRAL_WS_CACHE, dims)::_BVK2SpectralWorkspace{D}
    len = NTuple{D,Float64}(Float64(l) for l in lengths)
    if len != ws.lengths          # cell changed -> rebuild the k^2 grid only
        ws.k2 .= _periodic_k2_nd(dims, collect(Float64, lengths))
        ws.lengths = len
    end
    return ws
end

# Function barrier: `forward`/`inverse` enter concretely typed so the FFT
# `mul!` specialises.  Applies A_psi * (-Laplacian)^{-1} psi into `out`,
# reproducing the out-of-place path to floating-point roundoff.
function _bvk2_apply_inverse_laplacian!(out::AbstractArray{Float64},
        psi::AbstractArray{Float64}, buffer::AbstractArray{ComplexF64},
        k2::AbstractArray{Float64}, forward, inverse, a_psi::Float64)
    @inbounds @. buffer = complex(psi)
    forward * buffer                      # in-place forward FFT
    tforeach(eachindex(buffer)) do i
        @inbounds buffer[i] =
            iszero(k2[i]) ? zero(ComplexF64) : buffer[i] / k2[i]
    end
    inverse * buffer                      # in-place inverse FFT (normalised)
    @inbounds @. out = a_psi * real(buffer)
    return out
end

"""
    _diblock_bvk2_gradient_chemical_potential_nd(phi, lengths, params; adaptive)

Differentiate the *discrete* BVK2 edge-gradient energy.  The return value is
the nodal chemical potential, i.e. the derivative after division by the
uniform cell volume `dV`.  This is deliberately a reverse accumulation of the
same central-difference sensor and forward-edge quotient used by
`_diblock_bvk2_kpsi_nd` and `_diblock_bvk2_gradient_integral_nd`:

```
E_grad/dV = sum_edges (K_i + K_j)/2 * ((phi_j-phi_i)/h)^2 /
                                      (phi_e*(1-phi_e)).
```

Consequently the force includes both the direct edge-quotient derivative and
the full constitutive derivative through `K(M(Q,P))`, `Q`, `P`, the central
gradient/Laplacian stencils, and `eta(phi)`.
"""
function _diblock_bvk2_gradient_chemical_potential_nd(
        phi::AbstractArray{Float64}, lengths::AbstractVector{Float64},
        params; adaptive::Bool=true)
    dims = size(phi)
    length(lengths) == length(dims) ||
        throw(ArgumentError("lengths must have one entry per array dimension"))
    spacings = [lengths[dim] / dims[dim] for dim in eachindex(dims)]

    ws = _bvk2_chempot_workspace_cached(dims)
    eta = ws.eta
    eta_prime = ws.eta_prime
    tforeach(CartesianIndices(phi)) do index
        @inbounds eta[index] = diblock_bvk1_eta_psi(phi[index]; f=params.f)
        @inbounds eta_prime[index] = _diblock_bvk2_eta_prime(phi[index], params.f)
    end

    # Store the node-local sensor intermediates needed by the reverse pass.
    gradient_components = ws.gradient_components
    eta_radius2 = ws.eta_radius2
    reduced_q = ws.reduced_q
    reduced_p = ws.reduced_p
    sensor = ws.sensor
    kpsi = ws.kpsi
    # Per-node write (each node owns its slot in every output array); reads
    # only neighbour eta, so no cross-node dependence.
    tforeach(CartesianIndices(phi)) do index
        grad2 = 0.0
        lap = 0.0
        @inbounds for dim in eachindex(dims)
            left = _periodic_shift_index(index, dims, dim, -1)
            right = _periodic_shift_index(index, dims, dim, 1)
            spacing = spacings[dim]
            derivative = (eta[right] - eta[left]) / (2.0 * spacing)
            gradient_components[dim][index] = derivative
            grad2 += derivative^2
            lap += (eta[right] - 2.0 * eta[index] + eta[left]) / spacing^2
        end
        @inbounds begin
            radius2 = eta[index]^2 + params.eta_star2
            q = lap / (params.k_star2 * sqrt(radius2))
            p = grad2 / (params.k_star2 * radius2)
            curvature = q^2 / sqrt(1.0 + q^2)
            node_sensor = params.c2 * (2.0 * curvature + p)
            eta_radius2[index] = radius2
            reduced_q[index] = q
            reduced_p[index] = p
            sensor[index] = node_sensor
            kpsi[index] = adaptive ?
                params.K_psi0 * (1.0 + 0.5 * tanh(node_sensor)) : params.K_psi0
        end
    end

    # `chemical` is returned to the caller, so it stays a fresh allocation
    # (never the reused workspace) to prevent aliasing across calls.
    chemical = zeros(Float64, dims)
    k_adjoint = ws.k_adjoint

    # Direct derivative of each exact edge quotient and its K endpoint
    # adjoints.  Working with E_grad/dV makes `chemical` the continuum-style
    # nodal chemical potential without carrying dV through the reverse pass.
    # Gather form of the edge scatter: each node sums the contributions of
    # the two edges it caps along each axis -- the right edge (node is the
    # low endpoint) and the left edge (node is the high endpoint).  Every
    # edge is evaluated from both of its endpoints, doubling the edge
    # arithmetic, so each node writes only its own slot and the loop is
    # data-parallel.  Reproduces the scatter result to floating-point
    # roundoff (per-node accumulation order differs).
    tforeach(CartesianIndices(phi)) do index
        chem = 0.0
        kadj = 0.0
        @inbounds for dim in eachindex(dims)
            right = _periodic_shift_index(index, dims, dim, 1)
            left = _periodic_shift_index(index, dims, dim, -1)
            spacing = spacings[dim]
            # Right edge (index, right): `index` is the low endpoint.
            edge_phi = 0.5 * (phi[index] + phi[right])
            denominator = edge_phi * (1.0 - edge_phi)
            difference = (phi[right] - phi[index]) / spacing
            edge_k = 0.5 * (kpsi[index] + kpsi[right])
            common = 0.5 * difference^2 * (1.0 - 2.0 * edge_phi) / denominator^2
            chem += edge_k *
                (-2.0 * difference / (spacing * denominator) - common)
            kadj += 0.5 * difference^2 / denominator
            # Left edge (left, index): `index` is the high endpoint.
            edge_phi_l = 0.5 * (phi[left] + phi[index])
            denominator_l = edge_phi_l * (1.0 - edge_phi_l)
            difference_l = (phi[index] - phi[left]) / spacing
            edge_k_l = 0.5 * (kpsi[left] + kpsi[index])
            common_l =
                0.5 * difference_l^2 * (1.0 - 2.0 * edge_phi_l) / denominator_l^2
            chem += edge_k_l *
                (2.0 * difference_l / (spacing * denominator_l) - common_l)
            kadj += 0.5 * difference_l^2 / denominator_l
        end
        @inbounds chemical[index] = chem
        @inbounds k_adjoint[index] = kadj
    end

    adaptive || return chemical

    eta_adjoint = ws.eta_adjoint
    laplacian_adjoint = ws.laplacian_adjoint
    gradient_component_adjoint = ws.gradient_component_adjoint

    # K -> sensor -> (Q,P) -> (laplacian, |grad eta|^2, eta-radius).
    # Per-node write (eta_adjoint[index] += lands on this node's own zeroed
    # slot); the neighbour scatter into eta_adjoint happens in the next loop.
    tforeach(CartesianIndices(phi)) do index
        @inbounds begin
            tanh_sensor = tanh(sensor[index])
            sensor_adjoint = k_adjoint[index] * 0.5 * params.K_psi0 *
                (1.0 - tanh_sensor^2)
            q = reduced_q[index]
            p = reduced_p[index]
            radius2 = eta_radius2[index]
            curvature_prime = q * (2.0 + q^2) / (1.0 + q^2)^(1.5)
            q_adjoint = sensor_adjoint * params.c2 * 2.0 * curvature_prime
            p_adjoint = sensor_adjoint * params.c2
            laplacian_adjoint[index] =
                q_adjoint / (params.k_star2 * sqrt(radius2))
            radius2_adjoint = -0.5 * q_adjoint * q / radius2 -
                p_adjoint * p / radius2
            # `=` not `+=`: eta_adjoint is a reused workspace array (each node
            # written exactly once here); the neighbour scatter accumulates
            # onto this in the following gather loop.
            eta_adjoint[index] = 2.0 * eta[index] * radius2_adjoint
            grad2_adjoint = p_adjoint / (params.k_star2 * radius2)
            for dim in eachindex(dims)
                gradient_component_adjoint[dim][index] =
                    2.0 * gradient_components[dim][index] * grad2_adjoint
            end
        end
    end

    # Gather transpose of the periodic central-gradient and Laplacian
    # stencils: each node collects the adjoint mass its neighbours would
    # scatter onto it (reading their component/Laplacian adjoints) instead
    # of writing into them, and adds it to the eta_adjoint already seeded by
    # the sensor loop.  Central-gradient transpose is antisymmetric
    # (+left/-right); Laplacian transpose is the symmetric [1,-2,1] stencil.
    tforeach(CartesianIndices(phi)) do index
        acc = 0.0
        @inbounds for dim in eachindex(dims)
            left = _periodic_shift_index(index, dims, dim, -1)
            right = _periodic_shift_index(index, dims, dim, 1)
            spacing = spacings[dim]
            component_adjoint = gradient_component_adjoint[dim]
            acc += (component_adjoint[left] - component_adjoint[right]) /
                (2.0 * spacing)
            acc += (laplacian_adjoint[left] - 2.0 * laplacian_adjoint[index] +
                    laplacian_adjoint[right]) / spacing^2
        end
        @inbounds eta_adjoint[index] += acc
    end

    @. chemical += eta_adjoint * eta_prime
    return chemical
end

function _diblock_bvk2_chemical_potential_validated(
        phi::AbstractArray{Float64}, cell_lengths::AbstractVector{Float64},
        ff::Float64, chi::Float64, nn::Float64, bb::Float64;
        adaptive::Bool, c2::Real, parameters=nothing)
    dims = size(phi)
    params = parameters === nothing ?
        diblock_bvk2_parameters(; f=ff, N=nn, b=bb, c2=c2) : parameters
    sw = _bvk2_spectral_workspace_cached(dims, cell_lengths)
    psi = phi .- ff
    chemical = similar(phi)
    _bvk2_apply_inverse_laplacian!(chemical, psi, sw.buffer, sw.k2,
        sw.forward, sw.inverse, params.A_psi)
    entropy_prefactor = params.M_psi
    fb = 1.0 - ff
    tforeach(eachindex(chemical)) do index
        @inbounds begin
            # Floor the entropy-derivative log at phi=1e-10 so the force stays
            # bounded at saturated cells (those cells are mobility-frozen in
            # the SD/Newton step); the linear interaction term keeps the true
            # phi.
            pc = clamp(phi[index], 1.0e-10, 1.0 - 1.0e-10)
            chemical[index] += entropy_prefactor *
                log(pc * fb / (ff * (1.0 - pc))) +
                chi * (1.0 - 2.0 * phi[index])
        end
    end
    chemical .+= _diblock_bvk2_gradient_chemical_potential_nd(
        phi, cell_lengths, params; adaptive=adaptive)
    return chemical
end

function _bvk2_energy_1d_validated(phi::Vector{Float64}, count::Int,
        ff::Float64, chi::Float64, period::Float64, nn::Float64,
        bb::Float64, params; adaptive::Bool)
    dx = period / count
    k2_values = _periodic_k2_1d(count, period)
    psi = phi .- ff
    connectivity = 0.5 * params.A_psi *
        _periodic_inverse_laplacian_pair_integral_1d(psi, psi, k2_values, dx)
    local_entropy = params.M_psi * dx * sum(
        _diblock_relative_entropy(value, ff) for value in phi)
    kpsi = _diblock_bvk2_kpsi_1d(phi, ff, period, nn, bb;
        adaptive=adaptive, c2=params.c2, parameters=params)
    gradient = _diblock_bvk2_gradient_integral_1d(phi, kpsi, dx)
    interaction = chi * dx * sum(
        value * (1.0 - value) - ff * (1.0 - ff) for value in phi)
    return connectivity + local_entropy + gradient + interaction
end

"""
    diblock_bvk2_energy_1d(phi_a; f, chiN, L, N=1, b=1, adaptive=true, c2)

Evaluate the BVK2 excess free energy on a periodic one-dimensional density
field.  The edge quotient uses its exact open-domain denominator
`phi_e*(1-phi_e)`; endpoint densities are rejected by the shared input gate.
"""
function diblock_bvk2_energy_1d(phi_a::AbstractVector{<:Real};
        f::Real=0.5, chiN::Real=12.0, L::Real=1.0, N::Real=1.0,
        b::Real=1.0, adaptive::Bool=true, c2::Real,
        mean_tolerance::Real=1.0e-8)
    phi, count, ff, chi, period, nn, bb = _diblock_profile_energy_inputs(phi_a;
        f=f, chiN=chiN, L=L, N=N, b=b, mean_tolerance=mean_tolerance)
    params = diblock_bvk2_parameters(; f=ff, N=nn, b=bb, c2=c2)
    return _bvk2_energy_1d_validated(
        phi, count, ff, chi, period, nn, bb, params; adaptive=adaptive)
end

"""
    triblock_aba_bvk2_energy_1d(phi_a; f, chiN, L, N=1, b=1,
                                adaptive=true, c2)

Evaluate the frozen-closure BVK2 energy for a symmetric linear ABA triblock.
The architecture-specific `A_psi`, `M_psi`, and ordering scale are derived
from the ABA ideal-chain correlation; the nonlinear `c2` law is unchanged.
"""
function triblock_aba_bvk2_energy_1d(phi_a::AbstractVector{<:Real};
        f::Real=0.5, chiN::Real=12.0, L::Real=1.0, N::Real=1.0,
        b::Real=1.0, adaptive::Bool=true, c2::Real,
        mean_tolerance::Real=1.0e-8)
    phi, count, ff, chi, period, nn, bb = _diblock_profile_energy_inputs(phi_a;
        f=f, chiN=chiN, L=L, N=N, b=b, mean_tolerance=mean_tolerance)
    params = triblock_aba_bvk2_parameters(; f=ff, N=nn, b=bb, c2=c2)
    return _bvk2_energy_1d_validated(
        phi, count, ff, chi, period, nn, bb, params; adaptive=adaptive)
end

"""
    diblock_bvk2_energy_nd(phi_a; f, chiN, lengths, N=1, b=1, adaptive=true, c2)

Evaluate the BVK2 excess free energy on a periodic N-dimensional density
field using the exact open-domain edge denominator.
"""
function diblock_bvk2_energy_nd(phi_a::AbstractArray{<:Real};
        f::Real=0.5, chiN::Real=12.0, lengths=1.0, N::Real=1.0,
        b::Real=1.0, adaptive::Bool=true, c2::Real,
        mean_tolerance::Real=1.0e-8)
    phi, dims, ff, chi, cell_lengths, nn, bb, dV =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    k2_values = _periodic_k2_nd(dims, cell_lengths)
    params = diblock_bvk2_parameters(; f=ff, N=nn, b=bb, c2=c2)
    psi = phi .- ff
    connectivity = 0.5 * params.A_psi *
        _periodic_inverse_laplacian_pair_integral_nd(psi, psi, k2_values, dV)
    local_entropy = params.M_psi * dV *
        tmapreduce(value -> _diblock_relative_entropy(value, ff), +, phi)
    kpsi = _diblock_bvk2_kpsi_nd(phi, ff, cell_lengths, nn, bb;
        adaptive=adaptive, c2=params.c2)
    gradient = _diblock_bvk2_gradient_integral_nd(phi, kpsi, cell_lengths, dV)
    interaction = chi * dV *
        tmapreduce(value -> value * (1.0 - value) - ff * (1.0 - ff), +, phi)
    return connectivity + local_entropy + gradient + interaction
end

"""
    diblock_bvk2_cell_scale_gradient_density_nd(phi_a; f, chiN, lengths,
        N=1, b=1, adaptive=true, c2, scaled_axes=nothing)

Return the analytic frozen-field derivative of the *discrete* BVK2 excess
free-energy density under a cell dilation of the selected axes,

```
d/ds [F(phi; lengths(s)) / prod(lengths(s))] at s = 1,
```

where `lengths(s)[d] = s * lengths[d]` for `d in scaled_axes` and is unchanged
otherwise.  `scaled_axes=nothing` selects every axis and therefore preserves
the original isotropic derivative.  A single integer or a collection of
distinct one-based axis indices may be supplied.

The nodal density samples `phi` are held fixed.  Local entropy and contact
interaction densities are invariant under this dilation.  For each Fourier
mode, the spectral inverse-Laplacian derivative is weighted by the selected
part of its squared wave number.  Only selected forward-edge quotients have
the explicit `s^-2` derivative.  In adaptive BVK2, the stiffness also changes
through the selected-axis parts of the reduced curvature and gradient
sensors.  This includes the complete analytic chain `K(M(Q,P))`, including
the smooth curvature response at `Q=0`.

At a converged fixed-composition field the envelope theorem identifies this
frozen-field derivative with the scalar cell stationarity condition.  Away
from field stationarity it remains an exact diagnostic derivative of the
implemented discrete energy, but is not by itself a stress-free-cell proof.
"""
function diblock_bvk2_cell_scale_gradient_density_nd(
        phi_a::AbstractArray{<:Real}; f::Real=0.5, chiN::Real=12.0,
        lengths=1.0, N::Real=1.0, b::Real=1.0,
        adaptive::Bool=true, c2::Real,
        scaled_axes=nothing,
        mean_tolerance::Real=1.0e-8)
    phi, dims, ff, _chi, cell_lengths, nn, bb, _dV =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    params = diblock_bvk2_parameters(; f=ff, N=nn, b=bb, c2=c2)
    point_count = length(phi)
    dimension_count = length(dims)
    axis_values = scaled_axes === nothing ? collect(1:dimension_count) :
        scaled_axes isa Integer ? [scaled_axes] : collect(scaled_axes)
    isempty(axis_values) &&
        throw(ArgumentError("scaled_axes must select at least one axis"))
    all(axis -> axis isa Integer && !(axis isa Bool), axis_values) ||
        throw(ArgumentError("scaled_axes must contain integer axis indices"))
    selected = falses(dimension_count)
    for axis_value in axis_values
        axis = Int(axis_value)
        1 <= axis <= dimension_count ||
            throw(ArgumentError("scaled axis $axis is outside 1:$dimension_count"))
        selected[axis] &&
            throw(ArgumentError("scaled_axes must not contain duplicate axes"))
        selected[axis] = true
    end

    # With fixed Fourier coefficients, k_d^2 -> k_d^2/s^2 on each selected
    # axis.  Thus d(1/k^2)/ds = 2*k_selected^2/k^4 at s=1.  Division by the
    # dilated volume cancels the Fourier quadrature volume factor, so this is
    # directly the connectivity-density derivative.
    values_by_dim = [
        _periodic_k2_1d(dims[dim], cell_lengths[dim])
        for dim in 1:dimension_count]
    psi = phi .- ff
    psi_hat = fft(psi)
    connectivity_scale_sum = 0.0
    for index in CartesianIndices(psi_hat)
        coordinates = Tuple(index)
        k2 = 0.0
        selected_k2 = 0.0
        for dim in 1:dimension_count
            directional_k2 = values_by_dim[dim][coordinates[dim]]
            k2 += directional_k2
            selected[dim] && (selected_k2 += directional_k2)
        end
        k2 == 0.0 && continue
        connectivity_scale_sum += abs2(psi_hat[index]) * selected_k2 / k2^2
    end
    # Julia's unnormalized FFT contributes one factor of `point_count` in
    # Parseval's identity, while the inverse-Laplacian integral and the
    # volume normalization contribute two inverse factors in total.
    derivative_density =
        params.A_psi * connectivity_scale_sum / point_count^2

    spacings = [cell_lengths[dim] / dims[dim] for dim in eachindex(dims)]
    eta = similar(phi)
    for index in CartesianIndices(phi)
        eta[index] = diblock_bvk1_eta_psi(phi[index]; f=ff)
    end

    kpsi = fill(params.K_psi0, dims)
    kpsi_scale_derivative = zeros(Float64, dims)
    if adaptive
        for index in CartesianIndices(phi)
            grad2 = 0.0
            selected_grad2 = 0.0
            laplacian = 0.0
            selected_laplacian = 0.0
            for dim in eachindex(dims)
                left = _periodic_shift_index(index, dims, dim, -1)
                right = _periodic_shift_index(index, dims, dim, 1)
                spacing = spacings[dim]
                eta_gradient = (eta[right] - eta[left]) / (2.0 * spacing)
                directional_grad2 = eta_gradient^2
                directional_laplacian =
                    (eta[right] - 2.0 * eta[index] + eta[left]) / spacing^2
                grad2 += directional_grad2
                laplacian += directional_laplacian
                if selected[dim]
                    selected_grad2 += directional_grad2
                    selected_laplacian += directional_laplacian
                end
            end

            radius2 = eta[index]^2 + params.eta_star2
            q = laplacian / (params.k_star2 * sqrt(radius2))
            p = grad2 / (params.k_star2 * radius2)
            selected_q =
                selected_laplacian / (params.k_star2 * sqrt(radius2))
            selected_p = selected_grad2 / (params.k_star2 * radius2)
            curvature = q^2 / sqrt(1.0 + q^2)
            sensor = params.c2 * (2.0 * curvature + p)
            tanh_sensor = tanh(sensor)
            kpsi[index] =
                params.K_psi0 * (1.0 + 0.5 * tanh_sensor)

            # The selected-axis stencil pieces obey Q'=-2*Q_selected and
            # P'=-2*P_selected.  For C(Q)=Q^2/sqrt(1+Q^2),
            # C'(Q)=Q*(2+Q^2)/(1+Q^2)^(3/2), which vanishes continuously
            # at Q=0.  Therefore M'=-2*c2*(2*C'(Q)*Q_selected+P_selected).
            curvature_prime = q * (2.0 + q^2) / (1.0 + q^2)^(1.5)
            sensor_scale_derivative = -2.0 * params.c2 *
                (2.0 * selected_q * curvature_prime + selected_p)
            kpsi_scale_derivative[index] = 0.5 * params.K_psi0 *
                (1.0 - tanh_sensor^2) * sensor_scale_derivative
        end
    end

    # F_grad/V is the arithmetic mean of the directional edge densities.
    # A selected-axis quotient contributes -2*K*q_edge from its explicit
    # s^-2 scaling; every edge also receives K'*q_edge from the adaptive
    # constitutive response.
    for index in CartesianIndices(phi)
        for dim in eachindex(dims)
            right = _periodic_shift_index(index, dims, dim, 1)
            spacing = spacings[dim]
            edge_phi = 0.5 * (phi[index] + phi[right])
            denominator = edge_phi * (1.0 - edge_phi)
            edge_gradient = (phi[right] - phi[index]) / spacing
            quotient = edge_gradient^2 / denominator
            edge_k = 0.5 * (kpsi[index] + kpsi[right])
            edge_k_scale_derivative = 0.5 *
                (kpsi_scale_derivative[index] +
                 kpsi_scale_derivative[right])
            derivative_density +=
                (edge_k_scale_derivative -
                 (selected[dim] ? 2.0 * edge_k : 0.0)) *
                quotient / point_count
        end
    end
    return derivative_density
end

"""
    diblock_bvk2_chemical_potential_1d(phi_a; f, chiN, L, N=1, b=1,
                                       adaptive=true, c2)

Return the analytic derivative of the discrete BVK2 energy with respect to a
one-dimensional periodic density field.  With `dx=L/length(phi_a)`, the
normalization is `dF = dx*sum(mu .* dphi)`.  The mean constraint is not
projected here; subtract `mean(mu)` when a fixed-composition residual is
required.
"""
function diblock_bvk2_chemical_potential_1d(phi_a::AbstractVector{<:Real};
        f::Real=0.5, chiN::Real=12.0, L::Real=1.0, N::Real=1.0,
        b::Real=1.0, adaptive::Bool=true, c2::Real,
        mean_tolerance::Real=1.0e-8)
    phi, _count, ff, chi, period, nn, bb = _diblock_profile_energy_inputs(phi_a;
        f=f, chiN=chiN, L=L, N=N, b=b, mean_tolerance=mean_tolerance)
    return _diblock_bvk2_chemical_potential_validated(
        phi, [period], ff, chi, nn, bb; adaptive=adaptive, c2=c2)
end

"""
    triblock_aba_bvk2_chemical_potential_1d(phi_a; f, chiN, L, N=1, b=1,
                                            adaptive=true, c2)

Return the analytic discrete chemical potential of the symmetric ABA BVK2
energy. Its normalization and fixed-mean projection convention match
[`diblock_bvk2_chemical_potential_1d`].
"""
function triblock_aba_bvk2_chemical_potential_1d(
        phi_a::AbstractVector{<:Real};
        f::Real=0.5, chiN::Real=12.0, L::Real=1.0, N::Real=1.0,
        b::Real=1.0, adaptive::Bool=true, c2::Real,
        mean_tolerance::Real=1.0e-8)
    phi, _count, ff, chi, period, nn, bb = _diblock_profile_energy_inputs(phi_a;
        f=f, chiN=chiN, L=L, N=N, b=b, mean_tolerance=mean_tolerance)
    params = triblock_aba_bvk2_parameters(; f=ff, N=nn, b=bb, c2=c2)
    return _diblock_bvk2_chemical_potential_validated(
        phi, [period], ff, chi, nn, bb; adaptive=adaptive, c2=c2,
        parameters=params)
end

"""
    diblock_bvk2_chemical_potential_nd(phi_a; f, chiN, lengths, N=1, b=1,
                                       adaptive=true, c2)

Return the analytic derivative of the discrete BVK2 energy on an N-D
periodic grid.  For uniform cell volume `dV`, the normalization is
`dF = dV*sum(mu .* dphi)`.  It includes the full derivative of the adaptive
stiffness and the exact open-domain edge denominator.  The result is
unconstrained; fixed-mean solvers should project out its spatial mean.
"""
function diblock_bvk2_chemical_potential_nd(phi_a::AbstractArray{<:Real};
        f::Real=0.5, chiN::Real=12.0, lengths=1.0, N::Real=1.0,
        b::Real=1.0, adaptive::Bool=true, c2::Real,
        mean_tolerance::Real=1.0e-8)
    phi, _dims, ff, chi, cell_lengths, nn, bb, _dV =
        _diblock_profile_energy_inputs_nd(phi_a; f=f, chiN=chiN,
            lengths=lengths, N=N, b=b, mean_tolerance=mean_tolerance)
    return _diblock_bvk2_chemical_potential_validated(
        phi, cell_lengths, ff, chi, nn, bb; adaptive=adaptive, c2=c2)
end

@inline function _diblock_bvk2_open_logistic(value::Real)
    phi = _logistic_stable(value)
    phi == 0.0 && return nextfloat(0.0)
    phi == 1.0 && return prevfloat(1.0)
    return phi
end

# Composition sum for the logit-shift Newton step: returns the shifted
# logistic mean and its s-derivative mean in ONE serial pairwise pass, packed
# as a complex accumulator whose real part reproduces `Statistics.mean`
# bit-for-bit (the mean constraint is verified at ~2e-15, so this reduction
# must not be threaded).  `s` is an argument, so the closure is type-stable and
# allocation-free -- calling `sum(do ...)` on the reassigned loop variable
# directly would box it.
function _diblock_bvk2_logit_mean_slope(logits::AbstractVector{<:Real},
        s::Float64, inv_n::Float64)
    acc = sum(logits) do value
        p = _diblock_bvk2_open_logistic(Float64(value) + s)
        complex(p, p * (1.0 - p))
    end
    return real(acc) * inv_n, imag(acc) * inv_n
end

function _diblock_bvk2_logit_shift(logits::AbstractVector{<:Real}, f::Float64)
    isempty(logits) && throw(ArgumentError("BVK2 logits must not be empty"))
    lo = -maximum(Float64, logits) - 80.0
    hi = -minimum(Float64, logits) + 80.0
    inv_n = 1.0 / length(logits)
    # Solve mean(sigmoid(logit + s)) = f for the composition-matching shift.
    # SAFEGUARDED NEWTON: the mean is monotone increasing in s, so each step
    # gathers the mean AND its derivative in a SINGLE serial pairwise pass --
    # packed as a complex accumulator whose real part reproduces the pairwise
    # sum that `Statistics.mean` uses downstream bit-for-bit (the composition
    # constraint is checked at ~2e-15, so the reduction cannot be threaded, as
    # a parallel summation order shifts the mean by more than that tolerance).
    # Newton reaches machine precision in ~6-8 evaluations versus ~55 bisection
    # sweeps; the maintained [lo, hi] bracket falls back to bisection whenever
    # a Newton step would leave it, keeping bisection-grade robustness.
    s = 0.5 * (lo + hi)
    for _ in 1:80
        # Function barrier: `s` enters as an argument (never boxed), so the
        # pairwise sum closure stays type-stable and allocation-free.  Passing
        # the reassigned loop variable directly into a `sum(do ...)` here would
        # box it -- the classic captured-and-mutated-variable trap.
        mean_value, slope = _diblock_bvk2_logit_mean_slope(logits, s, inv_n)
        residual = mean_value - f
        residual < 0.0 ? (lo = s) : (hi = s)
        s_previous = s
        candidate = slope > 0.0 ? s - residual / slope : 0.5 * (lo + hi)
        s = (lo < candidate < hi) ? candidate : 0.5 * (lo + hi)
        # Newton lands on the root to full precision in a handful of steps;
        # stop once the step is at the ulp scale (further steps only limit-cycle
        # in the last bit).  `s` is already at the root, so the mean constraint
        # is satisfied exactly regardless of this tolerance.
        abs(s - s_previous) <= 1.0e-13 * (1.0 + abs(s)) && break
    end
    return s
end

"""
    _diblock_bvk2_logit_field_gradient!(phi, weights, gradient, logits,
                                        chemical, f, dx)

Map unconstrained logits to the open density interval while enforcing
`mean(phi) == f`, then optionally pull the physical chemical potential back
to the logit chart.  If `s(logits)` is the scalar mean-enforcing shift,

```
phi_i = logistic(logits_i + s),
sum_i phi_i / n = f,
```

implicit differentiation gives

```
ds/dlogits_j = -w_j / sum_i w_i,  w_i = phi_i(1-phi_i),
dF/dlogits_j = dx*w_j*(mu_j - sum_i(w_i*mu_i)/sum_i(w_i)).
```

The final subtraction of the gradient mean only removes roundoff in the
analytically null uniform-logit gauge.  `chemical` and `gradient` may be
`nothing` when only the feasible field is needed.
"""
function _diblock_bvk2_logit_field_gradient!(
        phi::Vector{Float64}, weights::Vector{Float64}, gradient,
        logits::AbstractVector{<:Real}, chemical, f::Float64, dx::Float64)
    length(phi) == length(weights) == length(logits) ||
        throw(DimensionMismatch("BVK2 logit work vectors must have equal length"))
    shift = _diblock_bvk2_logit_shift(logits, f)
    tforeach(eachindex(phi)) do index
        @inbounds value = _diblock_bvk2_open_logistic(Float64(logits[index]) + shift)
        @inbounds phi[index] = value
        @inbounds weights[index] = value * (1.0 - value)
    end
    chemical === nothing && return shift
    gradient === nothing && return shift
    length(chemical) == length(phi) == length(gradient) ||
        throw(DimensionMismatch("BVK2 chemical potential and gradient must match logits"))
    weight_sum = sum(weights)
    weight_sum > 0.0 || throw(ArgumentError("BVK2 logit chart lost all movable density weight"))
    weighted_chemical = dot(weights, chemical) / weight_sum
    tforeach(eachindex(gradient)) do index
        @inbounds gradient[index] = dx * weights[index] *
            (Float64(chemical[index]) - weighted_chemical)
    end
    # The map is invariant under logits -> logits + constant.  Enforce the
    # corresponding exact null gradient against floating-point accumulation.
    gradient .-= mean(gradient)
    return shift
end

function _diblock_bvk2_lamellar_logit_basis(count::Int)
    iseven(count) || throw(ArgumentError(
        "BVK2 lamellar cosine solver requires an even nx"))
    mode_count = count ÷ 2
    return [cos(2.0 * pi * mode * (index - 0.5) / count)
        for index in 1:count, mode in 1:mode_count]
end

function _diblock_bvk2_reflection_symmetrize(values::Vector{Float64})
    best = copy(values)
    best_defect = Inf
    for shift in 0:(length(values) - 1)
        candidate = circshift(values, shift)
        defect = sum(abs2, candidate .- reverse(candidate))
        if defect < best_defect
            best = candidate
            best_defect = defect
        end
    end
    return 0.5 .* (best .+ reverse(best))
end

function _minimize_diblock_bvk2_logit_lamella(; f::Real=0.5,
        chiN::Real=12.0, L=nothing, nx::Integer=128,
        mode_count::Integer=3, initial_amplitude::Real=0.2,
        initial_profile=nothing,
        N::Real=1.0, b::Real=1.0, adaptive::Bool=true, c2::Real,
        parameter_builder=diblock_bvk2_parameters,
        energy_evaluator=nothing, chemical_potential_evaluator=nothing,
        reference_period=nothing,
        max_iterations::Integer=10_000, lbfgs_memory::Integer=10,
        gradient_tolerance::Real=1.0e-10,
        force_tolerance::Real=1.0e-5,
        force_maxabs_tolerance::Real=2.0e-5,
        newton_polish_iterations::Integer=12,
        newton_difference_step::Real=1.0e-5,
        newton_trust_radius::Real=1.0)
    modes = Int(mode_count)
    modes >= 1 || throw(ArgumentError("mode_count must be positive"))
    count = Int(nx)
    count >= 8 || throw(ArgumentError("nx must be at least 8"))
    iterations = Int(max_iterations)
    iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    memory = Int(lbfgs_memory)
    memory > 0 || throw(ArgumentError("lbfgs_memory must be positive"))
    gtol = Float64(gradient_tolerance)
    gtol > 0.0 && isfinite(gtol) ||
        throw(ArgumentError("gradient_tolerance must be positive and finite"))
    force_tol = Float64(force_tolerance)
    force_max_tol = Float64(force_maxabs_tolerance)
    force_tol > 0.0 && isfinite(force_tol) ||
        throw(ArgumentError("force_tolerance must be positive and finite"))
    force_max_tol > 0.0 && isfinite(force_max_tol) ||
        throw(ArgumentError("force_maxabs_tolerance must be positive and finite"))
    polish_iterations = Int(newton_polish_iterations)
    polish_iterations >= 0 || throw(ArgumentError(
        "newton_polish_iterations must be nonnegative"))
    difference_step = Float64(newton_difference_step)
    difference_step > 0.0 && isfinite(difference_step) || throw(ArgumentError(
        "newton_difference_step must be positive and finite"))
    trust_radius = Float64(newton_trust_radius)
    trust_radius > 0.0 && isfinite(trust_radius) || throw(ArgumentError(
        "newton_trust_radius must be positive and finite"))

    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    chi = Float64(chiN)
    isfinite(chi) || throw(ArgumentError("chiN must be finite"))
    params = parameter_builder(; f=ff, N=nn, b=bb, c2=c2)
    default_period = reference_period === nothing ?
        2.0 * pi / sqrt(_diblock_kernel_minimum_k2(; f=ff, N=nn, b=bb)) :
        Float64(reference_period)
    period = L === nothing ? default_period : Float64(L)
    period > 0.0 && isfinite(period) ||
        throw(ArgumentError("lamella period L must be positive and finite"))
    dx = period / count

    initial_phi = if initial_profile === nothing
        last(diblock_lamella_profile_from_modes(
            [Float64(initial_amplitude); zeros(Float64, modes - 1)];
            f=ff, nx=count, L=period))
    else
        candidate = Float64.(collect(initial_profile))
        length(candidate) == count || throw(DimensionMismatch(
            "initial_profile length must equal nx"))
        all(isfinite, candidate) || throw(ArgumentError(
            "initial_profile must contain only finite values"))
        all(value -> 0.0 < value < 1.0, candidate) || throw(DomainError(
            (minimum(candidate), maximum(candidate)),
            "initial_profile must stay strictly inside (0, 1)"))
        candidate
    end
    logits = log.(initial_phi ./ (1.0 .- initial_phi))
    logits .-= mean(logits)
    logits = _diblock_bvk2_reflection_symmetrize(logits)
    basis = _diblock_bvk2_lamellar_logit_basis(count)
    logit_coefficients = (2.0 / count) .* (transpose(basis) * logits)
    logits = similar(initial_phi)
    phi = similar(logits)
    weights = similar(logits)
    chemical = similar(logits)
    logit_gradient = similar(logits)

    model_energy = energy_evaluator === nothing ?
        (field -> _bvk2_energy_1d_validated(
            field, count, ff, chi, period, nn, bb, params;
            adaptive=adaptive)) :
        (field -> energy_evaluator(
            field, params, ff, chi, period, nn, bb, adaptive))
    model_force = chemical_potential_evaluator === nothing ?
        (field -> _diblock_bvk2_chemical_potential_validated(
            field, [period], ff, chi, nn, bb; adaptive=adaptive, c2=c2,
            parameters=params)) :
        (field -> chemical_potential_evaluator(
            field, params, ff, chi, period, nn, bb, adaptive))

    function objective_gradient!(F, G, coefficients)
        mul!(logits, basis, coefficients)
        _diblock_bvk2_logit_field_gradient!(
            phi, weights, nothing, logits, nothing, ff, dx)
        energy = model_energy(phi)
        if G !== nothing
            chemical .= model_force(phi)
            _diblock_bvk2_logit_field_gradient!(
                phi, weights, logit_gradient, logits, chemical, ff, dx)
            mul!(G, transpose(basis), logit_gradient)
        end
        return F === nothing ? nothing : energy
    end

    optimization = Optim.optimize(
        Optim.NLSolversBase.only_fg!(objective_gradient!), logit_coefficients,
        Optim.LBFGS(m=memory,
            linesearch=Optim.LineSearches.BackTracking(order=3)),
        Optim.Options(iterations=iterations,
            g_tol=gtol, f_reltol=0.0, show_trace=false, store_trace=false,
            extended_trace=false))

    function physical_state(coefficients)
        state_logits = basis * coefficients
        state_phi = zeros(Float64, count)
        state_weights = similar(state_phi)
        _diblock_bvk2_logit_field_gradient!(state_phi, state_weights,
            nothing, state_logits, nothing, ff, dx)
        state_energy = model_energy(state_phi)
        state_chemical = model_force(state_phi)
        state_force = state_chemical .- mean(state_chemical)
        state_norm = sqrt(dx * sum(abs2, state_force))
        state_maxabs = maximum(abs, state_force)
        state_residual_coefficients =
            (2.0 / count) .* (transpose(basis) * state_force)
        return (coefficients=copy(coefficients), phi=state_phi,
            energy=state_energy, force=state_force,
            force_norm=state_norm, force_maxabs=state_maxabs,
            residual_coefficients=state_residual_coefficients)
    end

    state = physical_state(Vector{Float64}(Optim.minimizer(optimization)))
    completed_polish_iterations = 0
    for polish_iteration in 1:polish_iterations
        state.force_norm <= force_tol && state.force_maxabs <= force_max_tol &&
            break
        coefficient_count = length(state.coefficients)
        jacobian = Matrix{Float64}(undef, coefficient_count, coefficient_count)
        for column in 1:coefficient_count
            step = difference_step * (1.0 + abs(state.coefficients[column]))
            plus = copy(state.coefficients)
            minus = copy(state.coefficients)
            plus[column] += step
            minus[column] -= step
            plus_residual = physical_state(plus).residual_coefficients
            minus_residual = physical_state(minus).residual_coefficients
            jacobian[:, column] .= (plus_residual .- minus_residual) ./ (2.0 * step)
        end
        all(isfinite, jacobian) || break
        step_direction = try
            -(pinv(jacobian; rtol=1.0e-10) * state.residual_coefficients)
        catch
            break
        end
        all(isfinite, step_direction) || break
        largest_step = maximum(abs, step_direction)
        largest_step > trust_radius &&
            (step_direction .*= trust_radius / largest_step)
        accepted = false
        alpha = 1.0
        for _ in 1:24
            trial = physical_state(
                state.coefficients .+ alpha .* step_direction)
            if isfinite(trial.energy) && trial.force_norm < state.force_norm
                state = trial
                accepted = true
                completed_polish_iterations = polish_iteration
                break
            end
            alpha *= 0.5
        end
        accepted || break
    end

    phi .= state.phi
    energy = state.energy
    projected_force_norm = state.force_norm
    projected_force_maxabs = state.force_maxabs
    accepted = projected_force_norm <= force_tol &&
        projected_force_maxabs <= force_max_tol
    x = [(index - 0.5) * dx for index in 1:count]
    coefficients = _diblock_cosine_coefficients_from_profile(phi; f=ff,
        L=period, mode_count=modes)
    homogeneous = model_energy(fill(ff, count))
    return DiblockBurpLamellaResult(ff, chi, period, count, modes, x,
        copy(phi), coefficients, energy, homogeneous, accepted,
        Optim.iterations(optimization) + completed_polish_iterations,
        minimum(phi), maximum(phi), mean(phi),
        projected_force_norm, projected_force_maxabs)
end

"""
    minimize_diblock_bvk2_lamella(; f, chiN, L, c2, ...)

Relax a one-dimensional BVK2 lamellar profile at fixed period using a smooth,
mean-constrained logit chart and analytic L-BFGS gradient.  `c2` is required so
calibration and production calculations cannot silently use a provisional
coefficient.  The returned `converged` flag is fail-closed: it is true only
when both projected physical chemical-potential tolerances pass, independent
of the optimizer's internal stopping status.
"""
function minimize_diblock_bvk2_lamella(; f::Real=0.5, chiN::Real=12.0,
        L=nothing, nx::Integer=128, mode_count::Integer=3,
        initial_amplitude::Real=0.2, initial_profile=nothing,
        N::Real=1.0, b::Real=1.0,
        adaptive::Bool=true, c2::Real, max_iterations::Integer=10_000,
        lbfgs_memory::Integer=10, gradient_tolerance::Real=1.0e-10,
        force_tolerance::Real=1.0e-5,
        force_maxabs_tolerance::Real=2.0e-5,
        newton_polish_iterations::Integer=12,
        newton_difference_step::Real=1.0e-5,
        newton_trust_radius::Real=1.0)
    return _minimize_diblock_bvk2_logit_lamella(; f=f, chiN=chiN, L=L,
        nx=nx, mode_count=mode_count, initial_amplitude=initial_amplitude,
        initial_profile=initial_profile,
        N=N, b=b, adaptive=adaptive, c2=c2,
        max_iterations=max_iterations, lbfgs_memory=lbfgs_memory,
        gradient_tolerance=gradient_tolerance,
        force_tolerance=force_tolerance,
        force_maxabs_tolerance=force_maxabs_tolerance,
        newton_polish_iterations=newton_polish_iterations,
        newton_difference_step=newton_difference_step,
        newton_trust_radius=newton_trust_radius)
end

"""
    minimize_triblock_aba_bvk2_lamella(; f, chiN, L, c2, ...)

Relax a symmetric ABA BVK2 lamella at fixed period. The solver and nonlinear
closure are identical to the diblock path; only the analytically derived
ideal-chain parameter builder is replaced.
"""
function minimize_triblock_aba_bvk2_lamella(; f::Real=0.5,
        chiN::Real=12.0, L=nothing, N::Real=1.0, b::Real=1.0,
        c2::Real, kwargs...)
    params = triblock_aba_bvk2_parameters(; f=f, N=N, b=b, c2=c2)
    return _minimize_diblock_bvk2_logit_lamella(; f=f, chiN=chiN, L=L,
        N=N, b=b, c2=c2, parameter_builder=triblock_aba_bvk2_parameters,
        reference_period=2.0 * pi / params.k_star, kwargs...)
end

"""
    minimize_triblock_aba_bvk2_lamella_stress_free(; f, chiN, c2, ...)

Find the stress-free single-period ABA BVK2 lamella with the shared
field-relaxed scalar-period workflow. The search is centered on the ABA BVK2
ordering scale and never changes the supplied `c2`.
"""
function minimize_triblock_aba_bvk2_lamella_stress_free(; f::Real=0.5,
        chiN::Real=12.0, N::Real=1.0, b::Real=1.0, c2::Real,
        kwargs...)
    params = triblock_aba_bvk2_parameters(; f=f, N=N, b=b, c2=c2)
    return _minimize_diblock_lamella_stress_free("ABA-BVK2",
        minimize_triblock_aba_bvk2_lamella; f=f, chiN=chiN, N=N, b=b,
        c2=c2, reference_period=2.0 * pi / params.k_star, kwargs...)
end

"""
    minimize_diblock_bvk2_lamella_stress_free(; f, chiN, c2, ...)

Minimize the BVK2 lamellar energy density over the period using the shared
stress-free-period workflow and the analytic fixed-period BVK2 relaxation.
"""
function minimize_diblock_bvk2_lamella_stress_free(; f::Real=0.5,
        chiN::Real=12.0, nx::Integer=128, mode_count::Integer=3,
        initial_amplitudes=(0.2, 0.5, 1.0), initial_profile=nothing,
        N::Real=1.0, b::Real=1.0,
        adaptive::Bool=true, c2::Real, max_iterations::Integer=10_000,
        max_period_iterations::Integer=24, lower_factor::Real=0.75,
        upper_factor::Real=2.0, lbfgs_memory::Integer=10,
        gradient_tolerance::Real=1.0e-10,
        force_tolerance::Real=1.0e-5,
        force_maxabs_tolerance::Real=2.0e-5,
        newton_polish_iterations::Integer=12,
        newton_difference_step::Real=1.0e-5,
        newton_trust_radius::Real=1.0,
        bootstrap_period_factors=(),
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
        for (profile, amplitude) in attempts
            result = try
                minimize_diblock_bvk2_lamella(; f=ff, chiN=chi, L=period,
                    nx=nx, mode_count=mode_count,
                    initial_amplitude=amplitude, initial_profile=profile,
                    N=nn, b=bb, adaptive=adaptive, c2=c2,
                    max_iterations=max_iterations, lbfgs_memory=lbfgs_memory,
                    gradient_tolerance=gradient_tolerance,
                    force_tolerance=force_tolerance,
                    force_maxabs_tolerance=force_maxabs_tolerance,
                    newton_polish_iterations=newton_polish_iterations,
                    newton_difference_step=newton_difference_step,
                    newton_trust_radius=newton_trust_radius)
            catch
                nothing
            end
            result === nothing && continue
            identity_pass = result.energy < result.homogeneous_energy - 1.0e-10 &&
                result.maximum_phi - result.minimum_phi > 1.0e-4
            valid = result.converged && identity_pass && isfinite(result.energy)
            stress = valid ? diblock_bvk2_cell_scale_gradient_density_nd(
                result.phi_a; f=ff, chiN=chi, lengths=(period,), N=nn,
                b=bb, adaptive=adaptive, c2=c2) : NaN
            valid &= isfinite(stress)
            row = (factor=value, log_factor=log(value), period=period,
                objective=result.energy / period, stress=stress,
                result=result, amplitude=amplitude, valid=valid,
                identity_pass=identity_pass)
            if best === nothing || (row.valid && !best.valid) ||
                    (row.valid == best.valid && row.objective < best.objective)
                best = row
            end
            row.valid && break
        end
        best === nothing && throw(ErrorException(@sprintf(
            "all BVK2 fixed-period attempts failed at factor %.8g", value)))
        evaluations[value] = best
        _period_progress(progress_label, @sprintf(
            "BVK2 cell root factor=%.8g L/Rg=%.8g stress=%.8g field=%s force=%.4g maxforce=%.4g identity=%s",
            value, period * rg_scale, best.stress,
            string(best.result.converged),
            best.result.projected_force_norm,
            best.result.projected_force_maxabs,
            string(best.identity_pass)))
        return best
    end

    center = evaluate_factor(center_factor)
    center.valid || throw(ErrorException(
        "BVK2 bootstrap period did not produce a stationary lamella"))
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
        "BVK2 analytic cell stress has no valid oriented continuation bracket"))
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
        "AnalyticCellStressRoot", local_pass,
        length(evaluations), selected.amplitude, lo, hi,
        !isempty(bootstrap_factors), center_factor, window,
        selected.factor, selected.objective,
        left.factor, left.objective, right.factor, right.objective,
        check_fraction, local_pass, boundary_limited)
end
