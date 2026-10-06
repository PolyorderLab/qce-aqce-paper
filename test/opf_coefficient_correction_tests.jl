using Test
using DFMMonteCarlo

function corrected_liu2019_opf_coefficients(; f, chiN)
    spinodal = diblock_opf_spinodal(; f=f)
    x = chiN - spinodal
    g = 0.5 - f
    g2 = g^2
    g4 = g^4

    c2 = -(5.920 - 14.31g2 - 398.5g4 +
        x * (2.025 - 4.285g2 + 39.47g4) +
        x^2 * (0.005522 - 0.1915g2 - 1.003g4))
    c3 = -(g * (9.741 - 39.46g2 - 999.0g4) +
        x * g * (9.224 + 14.69g2 - 510.3g4) +
        x^2 * g * (0.06433 - 1.281g2 + 15.70g4))
    c4 = 9.686 + 53.00g2 - 1775.0g4 + 3.6x +
        x^2 * (0.02068 - 0.2385g2 - 0.4559g4)
    c5 = (0.7853 - 5.654g2 - 16.22g4 +
        x * (0.2126 - 1.170g2 + 3.659g4)) /
        (1.0 + x * (0.1185 - 0.7423g2 + 5.481g4))
    denominator6 = 0.5 + x * (0.2227 - 1.956g2 + 7.147g4) +
        x^2 * (0.0006666 - 0.02858g2 + 0.05316g4)
    return (c2=c2, c3=c3, c4=c4, c5=c5, c6=-c2 / denominator6)
end

@testset "Liu2019 corrected OPF Table 2 coefficients" begin
    @test DFMMonteCarlo._DIBLOCK_OPF_B2[1][3] == -398.5
    @test DFMMonteCarlo._DIBLOCK_OPF_B4[2][1] == 3.6
    @test DFMMonteCarlo._DIBLOCK_OPF_B5[2][1] == 0.2126
    @test DFMMonteCarlo._DIBLOCK_OPF_B6[2][1] == 0.2227

    for parameters in ((f=0.5, chiN=20.0), (f=0.4, chiN=30.0))
        actual = diblock_opf_coefficients(; parameters...)
        expected = corrected_liu2019_opf_coefficients(; parameters...)
        for coefficient in (:c2, :c3, :c4, :c5, :c6)
            @test getproperty(actual, coefficient) ≈ getproperty(expected, coefficient) rtol=1.0e-14
        end
    end
end

@testset "explicit OPF coefficients bypass the fitted regression" begin
    explicit = (c2=-1.0, c3=0.25, c4=2.0, c5=3.0, c6=4.0)
    evaluated = DFMMonteCarlo._diblock_opf_coefficients_for_evaluation(explicit;
        f=0.1, chiN=100.0, N=1.0, b=1.0, check_domain=true)
    @test evaluated == explicit
    @test diblock_opf_reference_period(; f=0.1, chiN=100.0,
        opf_coefficients=explicit) * sqrt(6.0) ≈ 2.0pi / (4.0 / 3.0)^0.25

    count = 96
    f = 0.35
    s = (0:(count - 1)) ./ count
    phi = @. f + 0.28 * cos(2.0pi * s) + 0.04 * cos(4.0pi * s) -
        0.015 * cos(6.0pi * s)
    fitted = fit_diblock_opf_coefficients_1d(phi;
        f=f, L=4.1 / sqrt(6.0), lambda=100.0)
    @test diblock_opf_reference_period(; f=f, chiN=30.0,
        opf_coefficients=fitted) * sqrt(6.0) ≈
        2.0pi / (fitted.c6 / fitted.c5)^0.25
end

@testset "OPF correction leaves Liu2019 OK coefficients unchanged" begin
    for parameters in ((f=0.3, chiN=20.0), (f=0.5, chiN=20.0),
            (f=0.7, chiN=20.0))
        f = parameters.f
        chiN = parameters.chiN
        ok = diblock_liu2019_ok_coefficients(; parameters...)
        mapped = corrected_liu2019_opf_coefficients(; parameters...)
        expected_c5 = 1.0 / (4.0 * f * (1.0 - f))
        expected_c6 = 3.0 / (4.0 * f^2 * (1.0 - f)^2)
        q2_star = sqrt(expected_c6 / expected_c5)
        expected_c2 = -chiN + diblock_opf_spinodal(; f=f) -
            expected_c5 * q2_star - expected_c6 / q2_star
        @test ok.c2 ≈ expected_c2 rtol=1.0e-14
        @test ok.c3 ≈ mapped.c3 rtol=1.0e-14
        @test ok.c4 ≈ mapped.c4 rtol=1.0e-14
        @test ok.c5 ≈ expected_c5 rtol=1.0e-14
        @test ok.c6 ≈ expected_c6 rtol=1.0e-14
    end
end

@testset "corrected OPF validation gates pass representative domain points" begin
    for f in (0.2, 0.35, 0.5, 0.65, 0.8)
        spinodal = diblock_opf_spinodal(; f=f)
        for chiN in (spinodal + 0.5, 35.0)
            coefficients = diblock_opf_coefficients(; f=f, chiN=chiN)
            @test all(isfinite,
                (coefficients.c2, coefficients.c3, coefficients.c4,
                    coefficients.c5, coefficients.c6))
            @test coefficients.c4 > 0.0
            @test coefficients.c5 > 0.0
            @test coefficients.c6 > 0.0
        end
    end
end
