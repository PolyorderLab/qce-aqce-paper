using Test
using Statistics
using DFMMonteCarlo

@testset "Liu2019 OPF coefficient regressions" begin
    @test DFMMonteCarlo._DIBLOCK_OPF_B4[2][1] == 3.6
    @test diblock_opf_spinodal(; f=0.5) ≈ 10.494868245363874 atol=1.0e-10
    @test diblock_opf_spinodal(; f=0.4) ≈ 11.343968614645185 atol=1.0e-10

    coefficients = diblock_opf_coefficients(; f=0.5, chiN=20.0)
    @test coefficients.c2 ≈ -25.666790861992418 rtol=1.0e-12
    @test coefficients.c3 == 0.0
    @test coefficients.c4 ≈ 45.77286123032753 rtol=1.0e-12
    @test coefficients.c5 ≈ 1.3196700000718236 rtol=1.0e-12
    @test coefficients.c6 ≈ 9.587827209148502 rtol=1.0e-12

    left = diblock_opf_coefficients(; f=0.4, chiN=30.0)
    right = diblock_opf_coefficients(; f=0.6, chiN=30.0)
    @test left.c2 ≈ right.c2 rtol=1.0e-12
    @test left.c3 ≈ -right.c3 rtol=1.0e-12
    @test left.c4 ≈ right.c4 rtol=1.0e-12
    @test left.c5 ≈ right.c5 rtol=1.0e-12
    @test left.c6 ≈ right.c6 rtol=1.0e-12

    @test_throws DomainError diblock_opf_coefficients(; f=0.19, chiN=30.0)
    @test_throws DomainError diblock_opf_coefficients(; f=0.5, chiN=40.0)
    @test_throws DomainError diblock_opf_coefficients(; f=0.5, chiN=10.0)
end

@testset "Liu2019 Ohta-Kawasaki comparison coefficients" begin
    f = 0.30
    chiN = 20.0
    coefficients = diblock_liu2019_ok_coefficients(; f=f, chiN=chiN)
    mapped = diblock_opf_coefficients(; f=f, chiN=chiN)
    @test coefficients.c5 ≈ 1.0 / (4.0 * f * (1.0 - f))
    @test coefficients.c6 ≈ 3.0 / (4.0 * f^2 * (1.0 - f)^2)
    k_star = (coefficients.c6 / coefficients.c5)^0.25
    @test coefficients.c2 ≈ -chiN + diblock_opf_spinodal(; f=f) -
        coefficients.c5 * k_star^2 - coefficients.c6 / k_star^2
    @test coefficients.c3 == mapped.c3
    @test coefficients.c4 == mapped.c4
    @test coefficients.mapping ==
        "Liu2019_OK_Eqs15-17_with_mapped_c3_c4_Eqs21-22"

    q2_star = sqrt(coefficients.c6 / coefficients.c5)
    k2_star = 6.0 * q2_star
    @test diblock_liu2019_ok_gaussian_kernel(k2_star; f=f) ≈
        2.0 * diblock_opf_spinodal(; f=f)

    q2 = 5.0
    k2 = 6.0 * q2
    full_hessian = 2.0 * (coefficients.c2 + coefficients.c5 * q2 +
        coefficients.c6 / q2)
    @test diblock_liu2019_ok_gaussian_kernel(k2; f=f) - 2.0 * chiN ≈
        full_hessian
    @test isinf(diblock_liu2019_ok_gaussian_kernel(0.0; f=f))
    @test_throws ArgumentError diblock_liu2019_ok_gaussian_kernel(-1.0; f=f)

    N = 8.0
    b = 1.7
    k2_star_scaled = 6.0 * q2_star / (N * b^2)
    @test diblock_liu2019_ok_gaussian_kernel(
        k2_star_scaled; f=f, N=N, b=b) ≈
        2.0 * diblock_opf_spinodal(; f=f, N=N, b=b)

    for chain_length in (2.0, 4.0)
        @test diblock_opf_spinodal(; f=f, N=chain_length, b=b) ≈
            diblock_opf_spinodal(; f=f, N=1.0, b=b) rtol=3.0e-10
        physical_k2 = 6.0 * q2 / (chain_length * b^2)
        @test diblock_liu2019_ok_gaussian_kernel(physical_k2;
            f=f, N=chain_length, b=b) ≈
            diblock_liu2019_ok_gaussian_kernel(6.0 * q2 / b^2;
                f=f, N=1.0, b=b) rtol=2.0e-14

        count = 128
        period_rg = 2.0pi / sqrt(q2)
        period = period_rg * b * sqrt(chain_length / 6.0)
        phase = 2.0pi .* (0:count-1) ./ count
        amplitude = 2.0e-5
        plus = f .+ amplitude .* cos.(phase)
        minus = f .- amplitude .* cos.(phase)
        ok_coefficients = diblock_liu2019_ok_coefficients(;
            f=f, chiN=chiN, N=chain_length, b=b)
        measured = 2.0 * (
            diblock_opf_energy_1d(plus; f=f, chiN=chiN, L=period,
                N=chain_length, b=b, opf_coefficients=ok_coefficients) +
            diblock_opf_energy_1d(minus; f=f, chiN=chiN, L=period,
                N=chain_length, b=b, opf_coefficients=ok_coefficients)
        ) / (period_rg * amplitude^2)
        expected = diblock_liu2019_ok_gaussian_kernel(physical_k2;
            f=f, N=chain_length, b=b) - 2.0 * chiN
        @test measured ≈ expected rtol=3.0e-6 atol=3.0e-6
    end
end

@testset "Liu2019 force/stress coefficient mapping" begin
    count = 96
    f = 0.35
    period_rg = 4.1
    period = period_rg / sqrt(6.0)
    s = (0:(count - 1)) ./ count
    phi = @. f + 0.28 * cos(2.0 * pi * s) +
        0.04 * cos(4.0 * pi * s) - 0.015 * cos(6.0 * pi * s)
    mapped = fit_diblock_opf_coefficients_1d(phi;
        f=f, L=period, lambda=100.0)
    @test mapped.c2 == -1.0
    @test all(isfinite, (mapped.c3, mapped.c4, mapped.c5, mapped.c6))
    @test mapped.c4 > 0.0
    @test mapped.c5 > 0.0
    @test mapped.c6 > 0.0
    expected_reference = 2.0 * pi / (mapped.c6 / mapped.c5)^0.25
    @test diblock_opf_reference_period(; f=f, chiN=30.0,
        opf_coefficients=mapped) * sqrt(6.0) ≈ expected_reference
end

@testset "Liu2019 OPF variational implementation" begin
    count = 64
    f = 0.4
    chiN = 30.0
    period_rg = 4.2
    period = period_rg / sqrt(6.0)
    x = [(index - 0.5) / count for index in 1:count]
    phi = @. f + 0.18 * cos(2.0 * pi * x) + 0.03 * cos(4.0 * pi * x)
    direction = @. sin(2.0 * pi * x) - 0.3 * sin(6.0 * pi * x)
    direction .-= mean(direction)
    direction ./= maximum(abs, direction)
    epsilon = 1.0e-6

    energy_plus = diblock_opf_energy_1d(phi .+ epsilon .* direction;
        f=f, chiN=chiN, L=period)
    energy_minus = diblock_opf_energy_1d(phi .- epsilon .* direction;
        f=f, chiN=chiN, L=period)
    finite_difference = (energy_plus - energy_minus) / (2.0 * epsilon)
    chemical = DFMMonteCarlo._diblock_opf_chemical_potential_1d(phi;
        f=f, chiN=chiN, L=period, N=1.0, b=1.0)
    analytic = (period_rg / count) * sum(chemical .* direction)
    @test analytic ≈ finite_difference rtol=2.0e-7 atol=1.0e-8

    strain = 1.0e-5
    density_plus = diblock_opf_energy_1d(phi; f=f, chiN=chiN,
        L=period * (1.0 + strain)) / (period_rg * (1.0 + strain))
    density_minus = diblock_opf_energy_1d(phi; f=f, chiN=chiN,
        L=period * (1.0 - strain)) / (period_rg * (1.0 - strain))
    stress_difference = (density_plus - density_minus) / (2.0 * strain)
    stress = diblock_opf_isotropic_stress_1d(phi;
        f=f, chiN=chiN, L=period)
    @test stress ≈ stress_difference rtol=2.0e-8 atol=1.0e-9
    @test diblock_opf_energy_1d(fill(f, count);
        f=f, chiN=chiN, L=period) == 0.0
end

@testset "Liu2019 OPF stress-free lamella smoke" begin
    result = minimize_diblock_opf_lamella_stress_free(; f=0.5, chiN=20.0,
        nx=48, mode_count=4, initial_amplitudes=(0.1,),
        max_iterations=3000, max_period_iterations=12,
        local_check_fraction=0.02, period_strategy=:stress_root,
        tolerance=1.0e-7, gradient_tolerance=1.0e-5)
    @test result.result.converged
    @test result.local_minimum_check_pass
    @test !result.boundary_limited
    @test 3.7 < result.period_rg < 4.4
    @test abs(result.result.mean_phi - 0.5) < 1.0e-10
    @test abs(diblock_opf_isotropic_stress_1d(result.result.phi_a;
        f=0.5, chiN=20.0, L=result.result.L)) < 1.0e-4
end

@testset "Liu2019 OPF order parameter is not pointwise bounded" begin
    result = minimize_diblock_opf_lamella_stress_free(; f=0.2, chiN=35.0,
        nx=64, mode_count=5, initial_amplitudes=(0.1, 0.2, 0.5),
        max_iterations=5000, max_period_iterations=16,
        period_strategy=:stress_root, tolerance=1.0e-7,
        gradient_tolerance=1.0e-5)
    @test result.result.converged
    @test result.local_minimum_check_pass
    @test !result.boundary_limited
    @test 3.30 < result.period_rg < 3.52
    @test result.result.minimum_phi < -0.15
    @test result.result.maximum_phi < 1.0
    @test abs(result.result.mean_phi - 0.2) < 1.0e-10
end

@testset "fixed-composition lamellar relaxation uses a tangent descent" begin
    result = minimize_diblock_ohta_kawasaki_lamella(; f=0.25, chiN=20.0,
        L=3.4 / sqrt(6.0), nx=64, mode_count=3, initial_amplitude=1.0,
        max_iterations=10_000, relaxation_step=1.0e-2,
        tolerance=1.0e-8)
    @test result.converged
    @test result.iterations > 1
    @test result.projected_force_norm < 5.0e-4
    @test result.projected_force_maxabs < 4.0e-3
end

@testset "stress-root selection prioritizes verified cell minima" begin
    candidates = [
        (objective=-2.0, local_minimum_check_pass=false,
            result=(converged=true,), label=:lower_saddle),
        (objective=-1.0, local_minimum_check_pass=true,
            result=(converged=true,), label=:verified_minimum),
        (objective=-0.5, local_minimum_check_pass=true,
            result=(converged=false,), label=:unconverged_minimum),
    ]
    selected = DFMMonteCarlo._select_diblock_opf_stress_root_candidate(candidates)
    @test selected.label == :verified_minimum

    no_minimum = [
        (objective=-2.0, local_minimum_check_pass=false,
            result=(converged=true,), label=:lower),
        (objective=-1.0, local_minimum_check_pass=false,
            result=(converged=true,), label=:higher),
    ]
    @test DFMMonteCarlo._select_diblock_opf_stress_root_candidate(no_minimum).label ==
        :lower
end
