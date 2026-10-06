"""
Liu et al. optimized phase-field (OPF) model for AB diblock melts.

The implementation follows Macromolecules 52, 2878--2888 (2019),
Eqs. 13 and 20--24 and Tables 1--2. Lengths in the published functional
are measured in `R_g`; the public API retains this package's convention that
`L` is measured in the same units as `b`, and converts internally.
"""

# Table 2 includes all four replacements in the Liu et al. correction
# (DOI 10.1021/acs.macromol.9b01332): B2[1][3], B4[2][1], B5[2][1], B6[2][1].
const _DIBLOCK_OPF_B2 = (
    (5.920, -14.31, -398.5),
    (2.025, -4.285, 39.47),
    (0.005522, -0.1915, -1.003),
)

const _DIBLOCK_OPF_B3 = (
    (9.741, -39.46, -999.0),
    (9.224, 14.69, -510.3),
    (0.06433, -1.281, 15.70),
)

const _DIBLOCK_OPF_B4 = (
    (9.686, 53.00, -1775.0),
    (3.6, 0.0, 0.0),
    (0.02068, -0.2385, -0.4559),
)

const _DIBLOCK_OPF_B5 = (
    (0.7853, -5.654, -16.22),
    (0.2126, -1.170, 3.659),
    (0.1185, -0.7423, 5.481),
)

const _DIBLOCK_OPF_B6 = (
    (0.5, 0.0, 0.0),
    (0.2227, -1.956, 7.147),
    (0.0006666, -0.02858, 0.05316),
)

@inline function _diblock_opf_polynomial(table, x::Float64, g::Float64;
        odd_g::Bool=false)
    total = 0.0
    for j in 0:2, k in 0:2
        exponent = 2 * k + (odd_g ? 1 : 0)
        total += table[j + 1][k + 1] * x^j * g^exponent
    end
    return total
end

"""
    diblock_opf_spinodal(; f=0.5, N=1, b=1)

Numerically evaluate the RPA spinodal `χN_s = min_k Γ₂(k)/2`. The result
reproduces Liu et al. Table 1 (for example `10.495` at `f=0.5`).
"""
function diblock_opf_spinodal(; f::Real=0.5, N::Real=1.0, b::Real=1.0)
    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    k2_star = _diblock_kernel_minimum_k2(; f=ff, N=nn, b=bb)
    return 0.5 * diblock_composition_kernel(k2_star; f=ff, N=nn, b=bb)
end

"""
    diblock_opf_coefficients(; f=0.5, chiN=20, N=1, b=1,
                              check_domain=true)

Return the published OPF coefficients `c2,...,c6`. The regression is certified
only for `0.2 ≤ f ≤ 0.8`, `χN_s < χN ≤ 35`.
"""
function diblock_opf_coefficients(; f::Real=0.5, chiN::Real=20.0,
        N::Real=1.0, b::Real=1.0, check_domain::Bool=true)
    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    chi = Float64(chiN)
    isfinite(chi) || throw(ArgumentError("chiN must be finite"))
    spinodal = diblock_opf_spinodal(; f=ff, N=nn, b=bb)
    if check_domain
        0.2 <= ff <= 0.8 || throw(DomainError(ff,
            "Liu2019 OPF regression is valid only for 0.2 <= f <= 0.8"))
        spinodal < chi <= 35.0 || throw(DomainError(chi,
            "Liu2019 OPF regression requires chiN_s < chiN <= 35"))
    end

    g = 0.5 - ff
    x = chi - spinodal
    c2 = -_diblock_opf_polynomial(_DIBLOCK_OPF_B2, x, g)
    c3 = -_diblock_opf_polynomial(_DIBLOCK_OPF_B3, x, g; odd_g=true)
    c4 = _diblock_opf_polynomial(_DIBLOCK_OPF_B4, x, g)
    numerator5 = sum(
        (_DIBLOCK_OPF_B5[1][k + 1] + _DIBLOCK_OPF_B5[2][k + 1] * x) *
        g^(2 * k) for k in 0:2)
    denominator5 = 1.0 + sum(
        _DIBLOCK_OPF_B5[3][k + 1] * x * g^(2 * k) for k in 0:2)
    denominator6 = _diblock_opf_polynomial(_DIBLOCK_OPF_B6, x, g)
    c5 = numerator5 / denominator5
    c6 = -c2 / denominator6
    all(isfinite, (c2, c3, c4, c5, c6)) ||
        throw(ErrorException("non-finite Liu2019 OPF regression coefficient"))
    c4 > 0.0 || throw(DomainError(c4, "OPF quartic coefficient must be positive"))
    c5 > 0.0 || throw(DomainError(c5, "OPF square-gradient coefficient must be positive"))
    c6 > 0.0 || throw(DomainError(c6, "OPF nonlocal coefficient must be positive"))
    return (c2=c2, c3=c3, c4=c4, c5=c5, c6=c6,
        f=ff, chiN=chi, chiN_spinodal=spinodal, x=x, g=g,
        N=nn, b=bb, regression="Liu2019_Eqs20-24_Table2")
end

function _diblock_opf_coefficients_for_evaluation(opf_coefficients;
        f::Real, chiN::Real, N::Real, b::Real, check_domain::Bool)
    opf_coefficients === nothing && return diblock_opf_coefficients(
        ; f=f, chiN=chiN, N=N, b=b, check_domain=check_domain)
    values = (c2=Float64(opf_coefficients.c2),
        c3=Float64(opf_coefficients.c3), c4=Float64(opf_coefficients.c4),
        c5=Float64(opf_coefficients.c5), c6=Float64(opf_coefficients.c6))
    all(isfinite, values) || throw(ArgumentError("OPF coefficients must be finite"))
    values.c4 > 0.0 || throw(DomainError(values.c4,
        "OPF quartic coefficient must be positive"))
    values.c5 > 0.0 || throw(DomainError(values.c5,
        "OPF square-gradient coefficient must be positive"))
    values.c6 > 0.0 || throw(DomainError(values.c6,
        "OPF nonlocal coefficient must be positive"))
    return values
end

"""
    fit_diblock_opf_coefficients_1d(phi_a; f, L, lambda=100, N=1, b=1)

Fit the OPF basis coefficients to a stress-free SCFT lamellar profile using
the force- and stress-matching metric of Liu et al. Eq. 18. The scale
degeneracy is removed by setting `c2=-1`; this common positive scale does not
affect the equilibrium profile or stress-free period.
"""
function fit_diblock_opf_coefficients_1d(phi_a::AbstractVector{<:Real};
        f::Real=0.5, L::Real, lambda::Real=100.0, N::Real=1.0,
        b::Real=1.0, mean_tolerance::Real=1.0e-8)
    phi, count, ff, _, period, nn, bb = _diblock_profile_energy_inputs(phi_a;
        f=f, chiN=0.0, L=L, N=N, b=b, mean_tolerance=mean_tolerance,
        enforce_density_bounds=true)
    weight = Float64(lambda)
    weight >= 0.0 && isfinite(weight) || throw(ArgumentError(
        "lambda must be nonnegative and finite"))
    rg = bb * sqrt(nn / 6.0)
    period_rg = period / rg
    dx_rg = period_rg / count
    k2_rg = _periodic_k2_1d(count, period_rg)
    delta = phi .- ff
    delta .-= mean(delta)
    laplacian = _periodic_laplacian_apply_1d(delta, k2_rg)
    inverse_laplacian = _periodic_inverse_laplacian_apply_1d(delta, k2_rg)
    force_basis = hcat(2.0 .* delta, 3.0 .* delta.^2,
        4.0 .* delta.^3, -2.0 .* laplacian,
        2.0 .* inverse_laplacian)
    gradient = _periodic_gradient_square_integral_1d(delta, k2_rg, dx_rg)
    nonlocal = _periodic_inverse_laplacian_pair_integral_1d(
        delta, delta, k2_rg, dx_rg)
    stress_basis = [0.0, 0.0, 0.0, -2.0 * gradient / period_rg,
        2.0 * nonlocal / period_rg]
    design = vcat(force_basis ./ sqrt(count),
        sqrt(weight) .* transpose(stress_basis))
    fitted = design[:, 2:5] \ design[:, 1]
    coefficients = (c2=-1.0, c3=fitted[1], c4=fitted[2],
        c5=fitted[3], c6=fitted[4])
    _diblock_opf_coefficients_for_evaluation(coefficients; f=ff, chiN=0.0,
        N=nn, b=bb, check_domain=false)
    residual = design * collect(coefficients)
    return merge(coefficients, (f=ff, lambda=weight,
        force_stress_residual_norm=norm(residual),
        mapping="Liu2019_Eq18_c2_fixed_-1"))
end

"""
    diblock_liu2019_ok_coefficients(; f=0.5, chiN=20, N=1, b=1)

Construct the Ohta--Kawasaki comparator used in Liu et al. Figure 2.  The
quadratic coefficients `c2`, `c5`, and `c6` follow their Eqs. 15--17.  Liu et
al. did not use phenomenological double-well coefficients: `c3` and `c4` are
inherited from their phase-field mapping, represented here by the published
Eqs. 21--22 regression (including the 2019 correction to Table 2).
"""
function diblock_liu2019_ok_coefficients(; f::Real=0.5, chiN::Real=20.0,
        N::Real=1.0, b::Real=1.0, check_domain::Bool=true)
    mapped = diblock_opf_coefficients(; f=f, chiN=chiN, N=N, b=b,
        check_domain=check_domain)
    ff = mapped.f
    chi = mapped.chiN
    nn = mapped.N
    bb = mapped.b
    spinodal = mapped.chiN_spinodal
    c5 = 1.0 / (4.0 * ff * (1.0 - ff))
    c6 = 3.0 / (4.0 * ff^2 * (1.0 - ff)^2)
    k_star = (c6 / c5)^0.25
    c2 = -chi + spinodal - c5 * k_star^2 - c6 / k_star^2
    coefficients = (c2=c2, c3=mapped.c3, c4=mapped.c4, c5=c5, c6=c6)
    _diblock_opf_coefficients_for_evaluation(coefficients; f=ff, chiN=chi,
        N=nn, b=bb, check_domain=false)
    return merge(coefficients, (f=ff, chiN=chi, chiN_spinodal=spinodal,
        N=nn, b=bb,
        mapping="Liu2019_OK_Eqs15-17_with_mapped_c3_c4_Eqs21-22"))
end

"""
    diblock_liu2019_ok_gaussian_kernel(k2; f=0.5, N=1, b=1)

Return the ideal-chain Gaussian vertex implied by the Ohta--Kawasaki
coefficients in Liu et al., Eqs. 15--17. The local quadratic coefficient is
chosen so that the interacting Hessian becomes unstable at the RPA spinodal,
whereas the preferred wavevector is set by the published `c5` and `c6`
coefficients.

The input `k2` is the squared physical wavevector. The returned per-chain
vertex uses the same normalization as `diblock_burp_ti_gaussian_kernel`, so
subtracting `2chiN` gives the homogeneous-state quadratic Hessian.
"""
function diblock_liu2019_ok_gaussian_kernel(k2::Real; f::Real=0.5,
        N::Real=1.0, b::Real=1.0)
    ff, nn, bb, _ = _diblock_validate(f, N, b, 0.0)
    kk = Float64(k2)
    isfinite(kk) || throw(ArgumentError("k2 must be finite"))
    kk >= 0.0 || throw(ArgumentError("k2 must be nonnegative"))
    kk == 0.0 && return Inf

    q2 = kk * nn * bb^2 / 6.0
    c5 = 1.0 / (4.0 * ff * (1.0 - ff))
    c6 = 3.0 / (4.0 * ff^2 * (1.0 - ff)^2)
    q2_star = sqrt(c6 / c5)
    spinodal = diblock_opf_spinodal(; f=ff, N=nn, b=bb)
    return 2.0 * (spinodal + c5 * (q2 - q2_star) +
        c6 * (inv(q2) - inv(q2_star)))
end

"""
    diblock_opf_reference_period(; f=0.5, chiN=20, N=1, b=1)

Return the single-mode OPF period in the package's physical length units.
The corresponding value in `R_g` units is `2π/(c6/c5)^(1/4)`.
"""
function diblock_opf_reference_period(; f::Real=0.5, chiN::Real=20.0,
        N::Real=1.0, b::Real=1.0, check_domain::Bool=true,
        opf_coefficients=nothing)
    coefficients = _diblock_opf_coefficients_for_evaluation(opf_coefficients;
        f=f, chiN=chiN, N=N, b=b, check_domain=check_domain)
    period_rg = 2.0 * pi / (coefficients.c6 / coefficients.c5)^0.25
    rg = Float64(b) * sqrt(Float64(N) / 6.0)
    return period_rg * rg
end

function diblock_opf_energy_1d(phi_a::AbstractVector{<:Real};
        f::Real=0.5, chiN::Real=20.0, L::Real=1.0, N::Real=1.0,
        b::Real=1.0, mean_tolerance::Real=1.0e-8,
        check_domain::Bool=true, opf_coefficients=nothing)
    phi, count, ff, chi, period, nn, bb = _diblock_profile_energy_inputs(phi_a;
        f=f, chiN=chiN, L=L, N=N, b=b, mean_tolerance=mean_tolerance,
        enforce_density_bounds=false)
    coefficients = _diblock_opf_coefficients_for_evaluation(opf_coefficients;
        f=ff, chiN=chi, N=nn, b=bb, check_domain=check_domain)
    rg = bb * sqrt(nn / 6.0)
    period_rg = period / rg
    dx_rg = period_rg / count
    k2_rg = _periodic_k2_1d(count, period_rg)
    delta = phi .- ff
    local_energy = dx_rg * sum(
        coefficients.c2 * value^2 +
        coefficients.c3 * value^3 +
        coefficients.c4 * value^4 for value in delta)
    gradient = coefficients.c5 *
        _periodic_gradient_square_integral_1d(delta, k2_rg, dx_rg)
    nonlocal = coefficients.c6 *
        _periodic_inverse_laplacian_pair_integral_1d(delta, delta, k2_rg, dx_rg)
    return local_energy + gradient + nonlocal
end

function _diblock_opf_chemical_potential_1d(phi_a::AbstractVector{Float64};
        f::Float64, chiN::Float64, L::Float64, N::Float64, b::Float64,
        check_domain::Bool=true, opf_coefficients=nothing)
    coefficients = _diblock_opf_coefficients_for_evaluation(opf_coefficients;
        f=f, chiN=chiN, N=N, b=b, check_domain=check_domain)
    rg = b * sqrt(N / 6.0)
    period_rg = L / rg
    k2_rg = _periodic_k2_1d(length(phi_a), period_rg)
    delta = phi_a .- f
    laplacian = _periodic_laplacian_apply_1d(delta, k2_rg)
    inverse_laplacian = _periodic_inverse_laplacian_apply_1d(delta, k2_rg)
    return @. 2.0 * coefficients.c2 * delta +
        3.0 * coefficients.c3 * delta^2 +
        4.0 * coefficients.c4 * delta^3 -
        2.0 * coefficients.c5 * laplacian +
        2.0 * coefficients.c6 * inverse_laplacian
end

"""
    diblock_opf_isotropic_stress_1d(phi_a; f, chiN, L, N=1, b=1)

Evaluate the derivative of the intensive OPF energy with respect to isotropic
cell strain, as given by Liu et al. Supporting Information Eq. S2. The profile
is held fixed in fractional cell coordinates. A stress-free cell has zero
stress.
"""
function diblock_opf_isotropic_stress_1d(phi_a::AbstractVector{<:Real};
        f::Real=0.5, chiN::Real=20.0, L::Real=1.0, N::Real=1.0,
        b::Real=1.0, mean_tolerance::Real=1.0e-8,
        check_domain::Bool=true, opf_coefficients=nothing)
    phi, count, ff, chi, period, nn, bb = _diblock_profile_energy_inputs(phi_a;
        f=f, chiN=chiN, L=L, N=N, b=b, mean_tolerance=mean_tolerance,
        enforce_density_bounds=false)
    coefficients = _diblock_opf_coefficients_for_evaluation(opf_coefficients;
        f=ff, chiN=chi, N=nn, b=bb, check_domain=check_domain)
    rg = bb * sqrt(nn / 6.0)
    period_rg = period / rg
    dx_rg = period_rg / count
    k2_rg = _periodic_k2_1d(count, period_rg)
    delta = phi .- ff
    gradient = _periodic_gradient_square_integral_1d(delta, k2_rg, dx_rg)
    nonlocal = _periodic_inverse_laplacian_pair_integral_1d(
        delta, delta, k2_rg, dx_rg)
    return (-2.0 * coefficients.c5 * gradient +
        2.0 * coefficients.c6 * nonlocal) / period_rg
end

function _minimize_diblock_opf_lamella_unconstrained(; f::Real=0.5,
        chiN::Real=20.0, L::Real, nx::Integer=128, mode_count::Integer=5,
        initial_amplitude::Real=0.2, N::Real=1.0, b::Real=1.0,
        max_iterations::Integer=10_000, tolerance::Real=1.0e-8,
        gradient_tolerance::Real=1.0e-8, check_domain::Bool=true,
        opf_coefficients=nothing)
    count = Int(nx)
    count >= 8 || throw(ArgumentError("nx must be at least 8"))
    iseven(count) || throw(ArgumentError("the OPF pseudospectral grid must be even"))
    reported_modes = Int(mode_count)
    reported_modes >= 1 || throw(ArgumentError("mode_count must be positive"))
    ff, nn, bb, _ = _diblock_validate(f, N, b, 1.0)
    chi = Float64(chiN)
    period = Float64(L)
    period > 0.0 && isfinite(period) || throw(ArgumentError(
        "lamella period L must be positive and finite"))
    amplitude = Float64(initial_amplitude)
    isfinite(amplitude) || throw(ArgumentError("initial_amplitude must be finite"))

    # Eq. 13 of Liu et al. defines a conserved polynomial phase field.  It is
    # not a volume-fraction parametrization and is therefore not constrained
    # pointwise to [0,1].  Expand the one-dimensional lamella in every resolved
    # cosine mode while omitting the constant mode, which enforces <phi_A>=f
    # exactly and fixes the otherwise arbitrary translational phase.
    resolved_modes = count ÷ 2
    s = [(index - 1) / count for index in 1:count]
    basis = Matrix{Float64}(undef, count, resolved_modes)
    for mode in 1:resolved_modes, index in 1:count
        basis[index, mode] = cospi(2.0 * mode * s[index])
    end
    coefficients0 = zeros(Float64, resolved_modes)
    coefficients0[1] = amplitude
    phi = zeros(Float64, count)
    chemical = similar(phi)
    projected = similar(phi)
    rg = bb * sqrt(nn / 6.0)
    dx_rg = (period / rg) / count

    function objective_gradient!(F, G, coefficients)
        mul!(phi, basis, coefficients)
        @. phi += ff
        energy = diblock_opf_energy_1d(phi; f=ff, chiN=chi, L=period,
            N=nn, b=bb, check_domain=check_domain,
            opf_coefficients=opf_coefficients)
        if G !== nothing
            chemical .= _diblock_opf_chemical_potential_1d(phi; f=ff,
                chiN=chi, L=period, N=nn, b=bb,
                check_domain=check_domain, opf_coefficients=opf_coefficients)
            chemical .-= mean(chemical)
            mul!(G, transpose(basis), chemical)
            G .*= dx_rg
        end
        return F === nothing ? nothing : energy
    end

    ncg_optimization = Optim.optimize(
        Optim.NLSolversBase.only_fg!(objective_gradient!), coefficients0,
        Optim.ConjugateGradient(
            linesearch=Optim.LineSearches.StrongWolfe(c_1=1.0e-4, c_2=0.1)),
        Optim.Options(iterations=Int(max_iterations),
            g_tol=Float64(gradient_tolerance),
            f_reltol=0.0, f_abstol=0.0, x_reltol=0.0, x_abstol=0.0,
            show_trace=false,
            store_trace=false, extended_trace=false))

    function physical_state(optimization)
        coefficients = Vector{Float64}(Optim.minimizer(optimization))
        mul!(phi, basis, coefficients)
        @. phi += ff
        energy = diblock_opf_energy_1d(phi; f=ff, chiN=chi, L=period,
            N=nn, b=bb, check_domain=check_domain,
            opf_coefficients=opf_coefficients)
        chemical .= _diblock_opf_chemical_potential_1d(phi; f=ff,
            chiN=chi, L=period, N=nn, b=bb, check_domain=check_domain,
            opf_coefficients=opf_coefficients)
        projected .= chemical .- mean(chemical)
        return (coefficients=coefficients, phi=copy(phi), energy=energy,
            residual_norm=sqrt(mean(abs2, projected)),
            residual_maxabs=maximum(abs, projected))
    end

    optimization = ncg_optimization
    state = physical_state(ncg_optimization)
    iteration_count = Optim.iterations(ncg_optimization)
    if state.residual_norm > Float64(gradient_tolerance)
        polish = Optim.optimize(
            Optim.NLSolversBase.only_fg!(objective_gradient!), state.coefficients,
            Optim.LBFGS(m=10,
                linesearch=Optim.LineSearches.StrongWolfe(c_1=1.0e-4, c_2=0.1)),
            Optim.Options(iterations=min(Int(max_iterations), 2_000),
                g_tol=0.1 * Float64(gradient_tolerance),
                f_reltol=0.0, f_abstol=0.0, x_reltol=0.0, x_abstol=0.0,
                show_trace=false, store_trace=false, extended_trace=false))
        polished_state = physical_state(polish)
        iteration_count += Optim.iterations(polish)
        energy_tolerance = 1.0e-12 * max(1.0, abs(state.energy))
        if polished_state.energy <= state.energy + energy_tolerance &&
                polished_state.residual_norm < state.residual_norm
            optimization = polish
            state = polished_state
        end
    end
    if state.residual_norm > Float64(gradient_tolerance)
        fallback = Optim.optimize(
            Optim.NLSolversBase.only_fg!(objective_gradient!), state.coefficients,
            Optim.LBFGS(m=20,
                linesearch=Optim.LineSearches.BackTracking(order=3)),
            Optim.Options(iterations=min(Int(max_iterations), 4_000),
                g_tol=0.05 * Float64(gradient_tolerance),
                f_reltol=0.0, f_abstol=0.0, x_reltol=0.0, x_abstol=0.0,
                show_trace=false, store_trace=false, extended_trace=false))
        fallback_state = physical_state(fallback)
        iteration_count += Optim.iterations(fallback)
        energy_tolerance = 1.0e-12 * max(1.0, abs(state.energy))
        if fallback_state.energy <= state.energy + energy_tolerance &&
                fallback_state.residual_norm < state.residual_norm
            optimization = fallback
            state = fallback_state
        end
    end
    phi .= state.phi
    energy = state.energy
    projected_force_norm = state.residual_norm
    projected_force_maxabs = state.residual_maxabs
    reported_coefficients = [
        2.0 / count * sum((phi[index] - ff) *
            cospi(2.0 * mode * s[index]) for index in 1:count)
        for mode in 1:reported_modes
    ]
    x = [index * period / count for index in 0:(count - 1)]
    homogeneous = diblock_opf_energy_1d(fill(ff, count); f=ff, chiN=chi,
        L=period, N=nn, b=bb, check_domain=check_domain,
        opf_coefficients=opf_coefficients)
    # `Optim` applies `gradient_tolerance` to the Fourier-coefficient gradient.
    # A line-search termination can nevertheless leave the stricter physical-
    # space residual below the same tolerance, which is also a valid stationary
    # solution. The physical residual remains reported and campaign-gated.
    converged = Optim.converged(optimization) ||
        projected_force_norm <= Float64(gradient_tolerance)
    return DiblockBurpLamellaResult(ff, chi, period, count, reported_modes,
        x, copy(phi), reported_coefficients, energy, homogeneous, converged,
        iteration_count, minimum(phi), maximum(phi), mean(phi),
        projected_force_norm, projected_force_maxabs)
end

function minimize_diblock_opf_lamella(; f::Real=0.5, chiN::Real=20.0,
        L=nothing, nx::Integer=128, mode_count::Integer=5,
        initial_amplitude::Real=0.2, N::Real=1.0, b::Real=1.0,
        max_iterations::Integer=10_000, relaxation_step::Real=1.0e-2,
        tolerance::Real=1.0e-8, gradient_tolerance::Real=1.0e-5,
        check_domain::Bool=true, opf_coefficients=nothing)
    period = L === nothing ?
        diblock_opf_reference_period(; f=f, chiN=chiN, N=N, b=b,
            check_domain=check_domain,
            opf_coefficients=opf_coefficients) : Float64(L)
    return _minimize_diblock_opf_lamella_unconstrained(; f=f, chiN=chiN,
        L=period,
        nx=nx, mode_count=mode_count, initial_amplitude=initial_amplitude,
        N=N, b=b, max_iterations=max_iterations,
        tolerance=tolerance,
        gradient_tolerance=gradient_tolerance,
        check_domain=check_domain, opf_coefficients=opf_coefficients)
end

function _select_diblock_opf_stress_root_candidate(candidates)
    isempty(candidates) && throw(ArgumentError(
        "at least one audited stress-root candidate is required"))
    local_minima = [candidate for candidate in candidates
        if candidate.local_minimum_check_pass]
    pool = isempty(local_minima) ? candidates : local_minima
    converged = [candidate for candidate in pool if candidate.result.converged]
    isempty(converged) || (pool = converged)
    return pool[argmin([candidate.objective for candidate in pool])]
end

function _minimize_diblock_opf_lamella_stress_free_root(; f::Real,
        chiN::Real, nx::Integer, mode_count::Integer, initial_amplitudes,
        N::Real, b::Real, max_iterations::Integer,
        max_period_iterations::Integer, lower_factor::Real,
        upper_factor::Real, tolerance::Real, gradient_tolerance::Real,
        local_check_fraction::Real, local_check_tolerance::Real,
        check_domain::Bool, opf_coefficients)
    reference = diblock_opf_reference_period(; f=f, chiN=chiN, N=N, b=b,
        check_domain=check_domain, opf_coefficients=opf_coefficients)
    lo = Float64(lower_factor)
    hi = Float64(upper_factor)
    0.0 < lo < hi || throw(ArgumentError(
        "lower_factor and upper_factor must define a positive interval"))
    amplitudes = _checked_initial_amplitudes(initial_amplitudes)
    evaluations = NamedTuple[]

    function evaluate(factor::Float64, amplitude::Float64)
        result = minimize_diblock_opf_lamella(; f=f, chiN=chiN,
            L=reference * factor, nx=nx, mode_count=mode_count,
            initial_amplitude=amplitude, N=N, b=b,
            max_iterations=max_iterations, tolerance=tolerance,
            gradient_tolerance=gradient_tolerance, check_domain=check_domain,
            opf_coefficients=opf_coefficients)
        stress = diblock_opf_isotropic_stress_1d(result.phi_a; f=f,
            chiN=chiN, L=result.L, N=N, b=b, check_domain=check_domain,
            opf_coefficients=opf_coefficients)
        item = (factor=factor, amplitude=amplitude, result=result,
            stress=stress, objective=result.energy / result.L)
        push!(evaluations, item)
        return item
    end

    root_candidates = NamedTuple[]
    scan_count = max(17, min(65, Int(max_period_iterations) + 9))
    for amplitude in amplitudes
        scan = [evaluate(factor, amplitude)
            for factor in range(lo, hi; length=scan_count)]
        for index in 1:(length(scan) - 1)
            left = scan[index]
            right = scan[index + 1]
            all(isfinite, (left.stress, right.stress)) || continue
            left.stress * right.stress <= 0.0 || continue
            for _ in 1:Int(max_period_iterations)
                center = evaluate((left.factor + right.factor) / 2.0, amplitude)
                if abs(center.stress) <= 1.0e-7 ||
                        right.factor - left.factor <= 1.0e-8
                    left = center
                    right = center
                    break
                elseif left.stress * center.stress <= 0.0
                    right = center
                else
                    left = center
                end
            end
            candidate = abs(left.stress) <= abs(right.stress) ? left : right
            abs(candidate.stress) <= 1.0e-4 && push!(root_candidates, candidate)
        end
    end
    isempty(root_candidates) && throw(ErrorException(
        "no zero-stress OPF lamellar root was found in the period interval"))
    fraction = Float64(local_check_fraction)
    audited_candidates = map(root_candidates) do candidate
        left_factor = candidate.factor * (1.0 - fraction)
        right_factor = candidate.factor * (1.0 + fraction)
        boundary_limited = !(lo < left_factor && right_factor < hi)
        left = boundary_limited ? nothing : evaluate(left_factor, candidate.amplitude)
        right = boundary_limited ? nothing : evaluate(right_factor, candidate.amplitude)
        objective_tolerance = Float64(local_check_tolerance) *
            max(1.0, abs(candidate.objective))
        local_pass = !boundary_limited &&
            left.objective > candidate.objective + objective_tolerance &&
            right.objective > candidate.objective + objective_tolerance
        merge(candidate, (left_factor=left_factor, right_factor=right_factor,
            left_objective=boundary_limited ? Inf : left.objective,
            right_objective=boundary_limited ? Inf : right.objective,
            local_minimum_check_pass=local_pass,
            boundary_limited=boundary_limited))
    end
    selected = _select_diblock_opf_stress_root_candidate(audited_candidates)
    rg_scale = sqrt(6.0 / (Float64(N) * Float64(b)^2))
    return DiblockLamellaPeriodOptimizationResult("Liu2019 OPF",
        selected.result, reference, selected.factor, selected.result.L,
        selected.result.L * rg_scale, selected.objective,
        "StressRootBisection", selected.result.converged &&
            abs(selected.stress) <= 1.0e-4,
        length(evaluations), selected.amplitude, lo, hi, false, NaN, NaN,
        selected.factor, selected.objective,
        selected.boundary_limited ? NaN : selected.left_factor,
        selected.left_objective,
        selected.boundary_limited ? NaN : selected.right_factor,
        selected.right_objective,
        fraction, selected.local_minimum_check_pass, selected.boundary_limited)
end

function minimize_diblock_opf_lamella_stress_free(; f::Real=0.5,
        chiN::Real=20.0, nx::Integer=128, mode_count::Integer=5,
        initial_amplitudes=(0.1, 0.2, 0.5), N::Real=1.0, b::Real=1.0,
        max_iterations::Integer=10_000, max_period_iterations::Integer=24,
        lower_factor::Real=0.75, upper_factor::Real=1.5,
        relaxation_step::Real=1.0e-2, tolerance::Real=1.0e-8,
        gradient_tolerance::Real=1.0e-5,
        bootstrap_period_factors=(), bootstrap_window::Real=0.12,
        local_check_fraction::Real=0.01, local_check_tolerance::Real=1.0e-8,
        period_strategy=:stress_root, progress_label::AbstractString="",
        check_domain::Bool=true, opf_coefficients=nothing)
    if period_strategy == :stress_root
        return _minimize_diblock_opf_lamella_stress_free_root(; f=f,
            chiN=chiN, nx=nx, mode_count=mode_count,
            initial_amplitudes=initial_amplitudes, N=N, b=b,
            max_iterations=max_iterations,
            max_period_iterations=max_period_iterations,
            lower_factor=lower_factor, upper_factor=upper_factor,
            tolerance=tolerance, gradient_tolerance=gradient_tolerance,
            local_check_fraction=local_check_fraction,
            local_check_tolerance=local_check_tolerance,
            check_domain=check_domain, opf_coefficients=opf_coefficients)
    end
    reference = diblock_opf_reference_period(; f=f, chiN=chiN, N=N, b=b,
        check_domain=check_domain, opf_coefficients=opf_coefficients)
    return _minimize_diblock_lamella_stress_free("Liu2019 OPF",
        minimize_diblock_opf_lamella; f=f, chiN=chiN, nx=nx,
        mode_count=mode_count, initial_amplitudes=initial_amplitudes, N=N, b=b,
        lower_factor=lower_factor, upper_factor=upper_factor,
        max_iterations=max_iterations, max_period_iterations=max_period_iterations,
        bootstrap_period_factors=bootstrap_period_factors,
        bootstrap_window=bootstrap_window, local_check_fraction=local_check_fraction,
        local_check_tolerance=local_check_tolerance, period_strategy=period_strategy,
        progress_label=progress_label, reference_period=reference,
        relaxation_step=relaxation_step, tolerance=tolerance,
        gradient_tolerance=gradient_tolerance,
        check_domain=check_domain, opf_coefficients=opf_coefficients)
end
