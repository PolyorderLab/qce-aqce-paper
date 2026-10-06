const MONTE_CARLO_PROJECT = normpath(joinpath(@__DIR__, ".."))
if !(MONTE_CARLO_PROJECT in LOAD_PATH)
    pushfirst!(LOAD_PATH, MONTE_CARLO_PROJECT)
end

using CSV
using DataFrames
using DFMMonteCarlo
using Optim
using Printf
using Statistics

function _arg_value(name, default)
    prefix = name * "="
    for arg in ARGS
        startswith(arg, prefix) && return split(arg, "=", limit=2)[2]
    end
    return default
end

function _parse_float_list(text::AbstractString)
    values = Float64[]
    for item in split(text, ",")
        stripped = strip(item)
        isempty(stripped) && continue
        push!(values, parse(Float64, stripped))
    end
    isempty(values) && throw(ArgumentError("list argument must contain at least one value"))
    return values
end

_rg(; N::Real=1.0, b::Real=1.0) = sqrt(Float64(N)) * Float64(b) / sqrt(6.0)
_k2_from_qrg(qrg::Real; N::Real=1.0, b::Real=1.0) = (Float64(qrg) / _rg(; N=N, b=b))^2
_qrg_from_k2(k2::Real; N::Real=1.0, b::Real=1.0) = sqrt(Float64(k2)) * _rg(; N=N, b=b)

function _ud_kernel(k2; f, N, b)
    return DFMMonteCarlo.diblock_uneyama_doi_composition_kernel(k2; f=f, N=N, b=b)
end

function _model_specs()
    return [
        (id="exact_rpa", label="Full RPA kernel (BURP-TI Gaussian limit)",
            color="#1f5fbf", dash="", kernel=(k2, f, N, b) ->
                DFMMonteCarlo.diblock_burp_ti_gaussian_kernel(k2; f=f, N=N, b=b)),
        (id="uneyama_doi", label="Literal nonlinear UD Hessian",
            color="#d97706", dash="2 3", kernel=(k2, f, N, b) ->
                DFMMonteCarlo.diblock_uneyama_doi_composition_kernel(
                    k2; f=f, N=N, b=b)),
        (id="bvk1_bvk2", label="Nominal reduced kernel shared by BVK1 and BVK2",
            color="#2a8c47", dash="6 4", kernel=(k2, f, N, b) ->
                DFMMonteCarlo.diblock_bvk1_gaussian_kernel(k2; f=f, N=N, b=b)),
        (id="ok", label="OK",
            color="#7a4fb0", dash="2 4", kernel=(k2, f, N, b) ->
                DFMMonteCarlo.diblock_liu2019_ok_gaussian_kernel(
                    k2; f=f, N=N, b=b)),
    ]
end

function _kernel_minimum(kernel, f::Real, N::Real, b::Real, qrg_min::Real, qrg_max::Real)
    k2_min = _k2_from_qrg(qrg_min; N=N, b=b)
    k2_max = _k2_from_qrg(qrg_max; N=N, b=b)
    objective(log_k2) = kernel(exp(log_k2), f, N, b)
    result = Optim.optimize(objective, log(k2_min), log(k2_max), Optim.GoldenSection())
    k2_star = exp(Optim.minimizer(result))
    kernel_star = kernel(k2_star, f, N, b)
    qrg_star = _qrg_from_k2(k2_star; N=N, b=b)
    return (
        k2_star=k2_star,
        qrg_star=qrg_star,
        period_rg=2.0 * pi / qrg_star,
        kernel_star=kernel_star,
        chiN_spinodal=0.5 * Float64(N) * kernel_star,
        converged=Optim.converged(result),
        iterations=Optim.iterations(result),
    )
end

function _fmt(value)
    v = Float64(value)
    if abs(v) >= 1.0e4 || (abs(v) < 1.0e-3 && v != 0.0)
        return @sprintf("%.4e", v)
    end
    return @sprintf("%.6g", v)
end

function _svg_x(x0, width, xmin, xmax, value)
    return x0 + width * (Float64(value) - xmin) / (xmax - xmin)
end

function _svg_y(y0, height, ymin, ymax, value)
    return y0 + height * (1.0 - (Float64(value) - ymin) / (ymax - ymin))
end

function _polyline(points)
    return join([@sprintf("%.3f,%.3f", point[1], point[2]) for point in points], " ")
end

function _nice_step(raw_step::Real)
    raw = Float64(raw_step)
    raw > 0.0 && isfinite(raw) || return 1.0
    exponent = floor(log10(raw))
    base = 10.0^exponent
    fraction = raw / base
    nice_fraction = if fraction <= 1.0
        1.0
    elseif fraction <= 2.0
        2.0
    elseif fraction <= 2.5
        2.5
    elseif fraction <= 5.0
        5.0
    else
        10.0
    end
    return nice_fraction * base
end

function _ticks(vmin::Real, vmax::Real; count::Integer=6)
    lo = Float64(vmin)
    hi = Float64(vmax)
    span = max(hi - lo, 1.0e-12)
    step = _nice_step(span / max(Int(count), 2))
    start = ceil(lo / step) * step
    stop = floor(hi / step) * step
    major = collect(start:step:stop)
    isempty(major) && (major = [lo, hi])
    minor_step = step / 2.0
    minor = [tick for tick in collect((start - minor_step):minor_step:(stop + minor_step))
        if lo <= tick <= hi && all(!isapprox(tick, m; atol=step * 1.0e-8) for m in major)]
    return (major=major, minor=minor)
end

function _write_kernel_svg(path, rows, f_values, qrg_min, qrg_max)
    x_min = qrg_min^2
    x_max = qrg_max^2
    panel_width = 360.0
    panel_height = 220.0
    margin_left = 78.0
    margin_top = 72.0
    gap_x = 78.0
    gap_y = 84.0
    width = margin_left + 2.0 * panel_width + gap_x + 250.0
    height = margin_top + length(f_values) * panel_height +
        (length(f_values) - 1) * gap_y + 92.0
    specs = _model_specs()
    color_by_model = Dict(spec.id => spec.color for spec in specs)
    dash_by_model = Dict(spec.id => spec.dash for spec in specs)
    label_by_model = Dict(spec.id => spec.label for spec in specs)

    open(path, "w") do io
        println(io, @sprintf("<svg xmlns='http://www.w3.org/2000/svg' width='%.0f' height='%.0f' viewBox='0 0 %.0f %.0f'>",
            width, height, width, height))
        println(io, "<rect width='100%' height='100%' fill='white'/>")
        println(io, "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1c1c1c}.axis{stroke:#202020;stroke-width:1.2}.grid{stroke:#dedede;stroke-width:0.75}.tick{stroke:#202020;stroke-width:1}.minor{stroke:#aaa;stroke-width:0.75}.line{fill:none;stroke-width:2.2;stroke-linejoin:round;stroke-linecap:round}</style>")
        println(io, "<text x='$(width / 2)' y='30' text-anchor='middle' font-size='18' font-weight='700'>Diblock Gaussian kernel and susceptibility comparison</text>")
        println(io, "<text x='$(width / 2)' y='51' text-anchor='middle' font-size='12'>BURP-TI preserves exact RPA at Gaussian order; BVK1 shares the Uneyama-Doi Gaussian skeleton.</text>")

        x_ticks = _ticks(x_min, x_max; count=7)
        for (row_index, f) in enumerate(f_values)
            row_y = margin_top + (row_index - 1) * (panel_height + gap_y)
            case_rows = [row for row in rows if isapprox(row.f, f; atol=1.0e-12)]
            exact_rows = [row for row in case_rows if row.model == "exact_rpa"]
            exact_smax = maximum(row.susceptibility for row in exact_rows)

            for col in 1:2
                x0 = margin_left + (col - 1) * (panel_width + gap_x)
                y0 = row_y
                y_label = col == 1 ? "S(k) / max S_RPA" : "(K_model - K_RPA) / K_RPA"
                panel_title = col == 1 ? "normalized susceptibility" : "relative inverse-kernel error"
                y_min = col == 1 ? 0.0 : -0.55
                y_max = col == 1 ? 1.08 : 1.05
                y_ticks = _ticks(y_min, y_max; count=6)

                println(io, @sprintf("<text x='%.3f' y='%.3f' font-size='13' font-weight='700'>fA = %.3g, %s</text>",
                    x0, y0 - 18.0, f, panel_title))

                for tick in x_ticks.minor
                    x = _svg_x(x0, panel_width, x_min, x_max, tick)
                    println(io, @sprintf("<line class='grid' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                        x, y0, x, y0 + panel_height))
                end
                for tick in y_ticks.minor
                    y = _svg_y(y0, panel_height, y_min, y_max, tick)
                    println(io, @sprintf("<line class='grid' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                        x0, y, x0 + panel_width, y))
                end
                for tick in x_ticks.major
                    x = _svg_x(x0, panel_width, x_min, x_max, tick)
                    println(io, @sprintf("<line class='tick' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                        x, y0 + panel_height, x, y0 + panel_height + 6.0))
                    println(io, @sprintf("<text x='%.3f' y='%.3f' text-anchor='middle' font-size='10'>%s</text>",
                        x, y0 + panel_height + 20.0, _fmt(tick)))
                end
                for tick in y_ticks.major
                    y = _svg_y(y0, panel_height, y_min, y_max, tick)
                    println(io, @sprintf("<line class='tick' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                        x0 - 6.0, y, x0, y))
                    println(io, @sprintf("<text x='%.3f' y='%.3f' text-anchor='end' dominant-baseline='middle' font-size='10'>%s</text>",
                        x0 - 10.0, y, _fmt(tick)))
                end
                zero_y = _svg_y(y0, panel_height, y_min, y_max, 0.0)
                if y0 <= zero_y <= y0 + panel_height
                    println(io, @sprintf("<line class='minor' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                        x0, zero_y, x0 + panel_width, zero_y))
                end
                println(io, @sprintf("<line class='axis' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x0, y0 + panel_height, x0 + panel_width, y0 + panel_height))
                println(io, @sprintf("<line class='axis' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x0, y0, x0, y0 + panel_height))
                println(io, @sprintf("<text x='%.3f' y='%.3f' text-anchor='middle' font-size='12'>k^2 R_g^2</text>",
                    x0 + panel_width / 2.0, y0 + panel_height + 42.0))
                println(io, @sprintf("<text transform='translate(%.3f %.3f) rotate(-90)' text-anchor='middle' font-size='12'>%s</text>",
                    x0 - 50.0, y0 + panel_height / 2.0, y_label))

                for spec in specs
                    model_rows = [row for row in case_rows if row.model == spec.id]
                    sort!(model_rows; by=row -> row.qRg)
                    points = Tuple{Float64,Float64}[]
                    for row in model_rows
                        raw_y = col == 1 ? row.susceptibility / exact_smax :
                            row.relative_kernel_error
                        clipped = clamp(raw_y, y_min, y_max)
                        push!(points, (_svg_x(x0, panel_width, x_min, x_max, row.qRg^2),
                            _svg_y(y0, panel_height, y_min, y_max, clipped)))
                    end
                    dash = isempty(dash_by_model[spec.id]) ? "" :
                        " stroke-dasharray='$(dash_by_model[spec.id])'"
                    println(io, "<polyline class='line' stroke='$(color_by_model[spec.id])'$(dash) points='$(_polyline(points))'/>")
                end
            end
        end

        legend_x = width - 225.0
        legend_y = margin_top + 8.0
        println(io, @sprintf("<text x='%.3f' y='%.3f' font-size='13' font-weight='700'>Gaussian models</text>",
            legend_x, legend_y - 18.0))
        for (idx, spec) in enumerate(specs)
            y = legend_y + 23.0 * (idx - 1)
            dash = isempty(spec.dash) ? "" : " stroke-dasharray='$(spec.dash)'"
            println(io, @sprintf("<line class='line' stroke='%s'%s x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                spec.color, dash, legend_x, y, legend_x + 36.0, y))
            println(io, @sprintf("<text x='%.3f' y='%.3f' dominant-baseline='middle' font-size='11'>%s</text>",
                legend_x + 44.0, y, spec.label))
        end
        println(io, "</svg>")
    end
    return path
end

function _write_note(path, summary_rows, settings)
    open(path, "w") do io
        println(io, "# Diblock Gaussian kernel comparison")
        println(io)
        println(io, "This artifact is the Uneyama-Doi-style weak-response companion to the stress-free lamella profile comparison. It compares the chi-independent parts of the second-order composition vertices for an incompressible AB diblock.")
        println(io)
        println(io, "The plotted models are:")
        println(io)
        println(io, "- `Exact RPA / BURP-TI Gaussian`: `K_RPA(k) = v^T h(k)^{-1} v`, with `v=(1,-1)`. BURP-TI reduces to this kernel because `eta_psi(phi)=phi-f+O((phi-f)^2)`.")
        println(io, "- `Literal nonlinear Uneyama-Doi Hessian`: `A_psi/k^2 + C_UD^(2) + b^2 k^2/[12 f(1-f)]`, where `C_UD^(2) = C_AA + C_BB - C_AB/[f(1-f)]`.")
        println(io, "- `BVK1/BVK2 reduced Gaussian`: `A_psi/k^2 + C_psi + b^2 k^2/[12 f(1-f)]`. The adaptive stiffness activates beyond second order.")
        println(io, "- `OK`: the Ohta--Kawasaki comparator with the published local, gradient, and nonlocal coefficients from Liu et al., Eqs. 15--17. The local coefficient reproduces the RPA spinodal, while the gradient and nonlocal coefficients determine an approximate preferred wavevector.")
        println(io)
        println(io, "## Minimum-kernel summary")
        println(io)
        println(io, "| fA | model | kRg at min | L/Rg | chiN spinodal | relative spinodal error |")
        println(io, "| ---: | --- | ---: | ---: | ---: | ---: |")
        for row in summary_rows
            println(io, @sprintf("| %.3g | %s | %s | %s | %s | %s |",
                row.f, row.model_label, _fmt(row.qRg_star), _fmt(row.period_rg),
                _fmt(row.chiN_spinodal), _fmt(row.relative_spinodal_error)))
        end
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, "1. BURP-TI is not a new weak-amplitude approximation: at Gaussian order it is exactly the full diblock RPA kernel.  Its approximation enters through the nonlinear log-density completion and thermodynamic integration.")
        println(io, "2. The literal nonlinear Uneyama-Doi continuation and the nominal reduced kernel have the same nonlocal and gradient coefficients but different local curvatures.  Their preferred weak-segregation wavelength is therefore the same, whereas their spinodals differ.")
        println(io, "3. BVK1 and BVK2 retain the nominal reduced kernel by construction. Their finite-amplitude differences arise from the adaptive gradient term, not from a hidden RPA fit.")
        println(io, "4. The OK model reproduces the RPA spinodal through its local quadratic coefficient, but its published gradient and nonlocal coefficients do not exactly reproduce the RPA minimizing wavevector. Spinodal agreement alone is therefore not sufficient evidence for accurate finite-amplitude profiles or periods.")
        println(io, "5. This motivates a manuscript structure with two evidence layers: first the Gaussian response shown here, then finite-amplitude stress-free lamella profiles and periods against SCFT.")
        println(io)
        println(io, "## Outputs")
        println(io)
        println(io, "- `kernel_samples.csv`: sampled inverse kernels and normalized susceptibilities.")
        println(io, "- `kernel_minima.csv`: minimum kernel, predicted weak-segregation period, and spinodal values.")
        println(io, "- `diblock_kernel_comparison.svg`: publication-style comparison figure.")
        println(io)
        println(io, "Settings: `N=$(settings.N)`, `b=$(settings.b)`, `kRg_min=$(settings.qrg_min)`, `kRg_max=$(settings.qrg_max)`, `sample_count=$(settings.sample_count)`.")
    end
    return path
end

function main()
    outdir = normpath(_arg_value("--outdir",
        joinpath(MONTE_CARLO_PROJECT, "results", "diblock_kernel_comparison")))
    mkpath(outdir)
    f_values = _parse_float_list(_arg_value("--f-values", "0.1,0.25,0.3,0.35,0.5"))
    sample_count = parse(Int, _arg_value("--sample-count", "420"))
    qrg_min = parse(Float64, _arg_value("--qrg-min", "0.05"))
    qrg_max = parse(Float64, _arg_value("--qrg-max", "5.0"))
    N = parse(Float64, _arg_value("--N", "1.0"))
    b = parse(Float64, _arg_value("--b", "1.0"))
    sample_count >= 64 || throw(ArgumentError("sample-count must be at least 64"))
    0.0 < qrg_min < qrg_max || throw(ArgumentError("expected 0 < qrg-min < qrg-max"))

    specs = _model_specs()
    qrg_values = collect(range(qrg_min, qrg_max; length=sample_count))
    sample_rows = NamedTuple[]
    summary_rows = NamedTuple[]

    for f in f_values
        exact_values = [
            DFMMonteCarlo.diblock_composition_kernel(_k2_from_qrg(qrg; N=N, b=b);
                f=f, N=N, b=b)
            for qrg in qrg_values
        ]
        exact_susceptibility = 1.0 ./ exact_values
        exact_smax = maximum(exact_susceptibility)
        exact_min = _kernel_minimum(specs[1].kernel, f, N, b, qrg_min, qrg_max)
        for spec in specs
            model_min = _kernel_minimum(spec.kernel, f, N, b, qrg_min, qrg_max)
            rel_spinodal = (model_min.chiN_spinodal - exact_min.chiN_spinodal) /
                exact_min.chiN_spinodal
            push!(summary_rows, (
                f=f,
                N=N,
                b=b,
                model=spec.id,
                model_label=spec.label,
                qRg_star=model_min.qrg_star,
                period_rg=model_min.period_rg,
                kernel_star=model_min.kernel_star,
                chiN_spinodal=model_min.chiN_spinodal,
                exact_qRg_star=exact_min.qrg_star,
                exact_period_rg=exact_min.period_rg,
                exact_chiN_spinodal=exact_min.chiN_spinodal,
                relative_spinodal_error=rel_spinodal,
                converged=model_min.converged,
                iterations=model_min.iterations,
            ))

            for (idx, qrg) in enumerate(qrg_values)
                k2 = _k2_from_qrg(qrg; N=N, b=b)
                exact_kernel = exact_values[idx]
                kernel = spec.kernel(k2, f, N, b)
                susceptibility = 1.0 / kernel
                push!(sample_rows, (
                    f=f,
                    N=N,
                    b=b,
                    qRg=qrg,
                    k2=k2,
                    model=spec.id,
                    model_label=spec.label,
                    kernel=kernel,
                    exact_kernel=exact_kernel,
                    relative_kernel_error=(kernel - exact_kernel) / exact_kernel,
                    susceptibility=susceptibility,
                    exact_susceptibility=exact_susceptibility[idx],
                    susceptibility_over_exact_smax=susceptibility / exact_smax,
                ))
            end
        end
    end

    samples_csv = CSV.write(joinpath(outdir, "kernel_samples.csv"),
        DataFrame(sample_rows))
    minima_csv = CSV.write(joinpath(outdir, "kernel_minima.csv"),
        DataFrame(summary_rows))
    svg_path = _write_kernel_svg(joinpath(outdir, "diblock_kernel_comparison.svg"),
        sample_rows, f_values, qrg_min, qrg_max)
    note_path = _write_note(joinpath(outdir, "kernel_comparison.md"), summary_rows,
        (N=N, b=b, qrg_min=qrg_min, qrg_max=qrg_max, sample_count=sample_count))

    println("Wrote $(samples_csv)")
    println("Wrote $(minima_csv)")
    println("Wrote $(svg_path)")
    println("Wrote $(note_path)")
end

main()
