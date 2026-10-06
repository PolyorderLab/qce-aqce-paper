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
    isempty(values) && throw(ArgumentError("float list cannot be empty"))
    return values
end

function _f_values()
    raw = _arg_value("--f-values", "")
    if !isempty(strip(raw))
        return _parse_float_list(raw)
    end
    fmin = parse(Float64, _arg_value("--f-min", "0.10"))
    fmax = parse(Float64, _arg_value("--f-max", "0.50"))
    count = parse(Int, _arg_value("--f-count", "81"))
    count >= 5 || throw(ArgumentError("f-count must be at least 5"))
    0.0 < fmin < fmax <= 0.5 || throw(ArgumentError("expected 0 < f-min < f-max <= 0.5"))
    return collect(range(fmin, fmax; length=count))
end

_rg(; N::Real=1.0, b::Real=1.0) = sqrt(Float64(N)) * Float64(b) / sqrt(6.0)
_k2_from_qrg(qrg::Real; N::Real=1.0, b::Real=1.0) =
    (Float64(qrg) / _rg(; N=N, b=b))^2
_qrg_from_k2(k2::Real; N::Real=1.0, b::Real=1.0) =
    sqrt(Float64(k2)) * _rg(; N=N, b=b)

function _model_specs()
    return [
        (id="exact_rpa", label="Full RPA kernel (BURP-TI Gaussian limit)",
            short_label="Full RPA (BURP-TI)", color="#1f5fbf", dash="",
            kernel=(k2, f, N, b) ->
                DFMMonteCarlo.diblock_burp_ti_gaussian_kernel(k2; f=f, N=N, b=b)),
        (id="uneyama_doi", label="Literal nonlinear UD Hessian",
            short_label="Nonlinear UD Hessian", color="#d97706", dash="2 3",
            kernel=(k2, f, N, b) ->
                DFMMonteCarlo.diblock_uneyama_doi_composition_kernel(
                    k2; f=f, N=N, b=b)),
        (id="bvk1_bvk2", label="Nominal reduced kernel shared by BVK1 and BVK2",
            short_label="Reduced kernel (BVK1; BVK2)", color="#2a8c47", dash="6 4",
            kernel=(k2, f, N, b) ->
                DFMMonteCarlo.diblock_bvk1_gaussian_kernel(k2; f=f, N=N, b=b)),
        (id="ok", label="OK",
            short_label="OK", color="#7a4fb0", dash="2 4",
            kernel=(k2, f, N, b) ->
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
    qrg_star = _qrg_from_k2(k2_star; N=N, b=b)
    kernel_star = kernel(k2_star, f, N, b)
    margin = (log(k2_star) - log(k2_min)) / (log(k2_max) - log(k2_min))
    return (
        k2_star=k2_star,
        qRg_star=qrg_star,
        period_rg=2.0 * pi / qrg_star,
        kernel_star=kernel_star,
        chiN_spinodal=0.5 * Float64(N) * kernel_star,
        converged=Optim.converged(result),
        iterations=Optim.iterations(result),
        boundary_limited=margin < 1.0e-3 || margin > 1.0 - 1.0e-3,
    )
end

function _fmt(value)
    v = Float64(value)
    if abs(v) >= 1.0e4 || (abs(v) < 1.0e-3 && v != 0.0)
        return @sprintf("%.3e", v)
    end
    return @sprintf("%.5g", v)
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
    lo == hi && (hi = lo + 1.0)
    span = hi - lo
    step = _nice_step(span / max(Int(count), 2))
    start = floor(lo / step) * step
    stop = ceil(hi / step) * step
    return collect(start:step:stop)
end

_svg_x(x0, width, xmin, xmax, value) =
    x0 + width * (Float64(value) - xmin) / (xmax - xmin)
_svg_y(y0, height, ymin, ymax, value) =
    y0 + height * (1.0 - (Float64(value) - ymin) / (ymax - ymin))

function _polyline(points)
    return join([@sprintf("%.3f,%.3f", point[1], point[2]) for point in points], " ")
end

function _field_value(row, field)
    return Float64(getproperty(row, field))
end

function _write_svg(path, rows, f_values)
    specs = _model_specs()
    color_by_model = Dict(spec.id => spec.color for spec in specs)
    dash_by_model = Dict(spec.id => spec.dash for spec in specs)
    label_by_model = Dict(spec.id => spec.short_label for spec in specs)
    panels = [
        (field=:chiN_spinodal, title="spinodal", ylabel="chi N at spinodal",
            error=false),
        (field=:period_rg, title="weak-segregation period", ylabel="L / R_g",
            error=false),
        (field=:relative_spinodal_error, title="spinodal error", ylabel="relative error",
            error=true),
        (field=:relative_period_error, title="period error", ylabel="relative error",
            error=true),
    ]

    panel_width = 315.0
    panel_height = 220.0
    margin_left = 82.0
    margin_top = 76.0
    gap_x = 76.0
    gap_y = 82.0
    width = margin_left + 2.0 * panel_width + gap_x + 230.0
    height = margin_top + 2.0 * panel_height + gap_y + 96.0
    xmin = minimum(f_values)
    xmax = maximum(f_values)
    x_ticks = _ticks(xmin, xmax; count=5)

    open(path, "w") do io
        println(io, @sprintf("<svg xmlns='http://www.w3.org/2000/svg' width='%.0f' height='%.0f' viewBox='0 0 %.0f %.0f'>",
            width, height, width, height))
        println(io, "<rect width='100%' height='100%' fill='white'/>")
        println(io, "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1c1c1c}.axis{stroke:#202020;stroke-width:1.2}.grid{stroke:#dedede;stroke-width:0.75}.tick{stroke:#202020;stroke-width:1}.zero{stroke:#777;stroke-width:0.9;stroke-dasharray:3 3}.line{fill:none;stroke-width:2.2;stroke-linejoin:round;stroke-linecap:round}</style>")
        println(io, "<text x='$(width / 2)' y='30' text-anchor='middle' font-size='18' font-weight='700'>Diblock weak-response map</text>")
        println(io, "<text x='$(width / 2)' y='51' text-anchor='middle' font-size='12'>Spinodal and weak-period predictions across composition; fA shown only to 0.5 by A/B symmetry.</text>")

        for (panel_index, panel) in enumerate(panels)
            col = (panel_index - 1) % 2
            row_index = (panel_index - 1) ÷ 2
            x0 = margin_left + col * (panel_width + gap_x)
            y0 = margin_top + row_index * (panel_height + gap_y)
            panel_rows = panel.error ? [row for row in rows if row.model != "exact_rpa"] : rows
            values = [_field_value(row, panel.field) for row in panel_rows]
            ymin = minimum(values)
            ymax = maximum(values)
            if panel.error
                limit = max(abs(ymin), abs(ymax), 1.0e-3)
                ymin = -1.15 * limit
                ymax = 1.15 * limit
            else
                padding = 0.08 * (ymax - ymin)
                ymin -= padding
                ymax += padding
            end
            y_ticks = _ticks(ymin, ymax; count=5)
            ymin = minimum(y_ticks)
            ymax = maximum(y_ticks)

            println(io, @sprintf("<text x='%.3f' y='%.3f' text-anchor='middle' font-size='13' font-weight='700'>%s</text>",
                x0 + panel_width / 2.0, y0 - 18.0, panel.title))
            for tick in x_ticks
                x = _svg_x(x0, panel_width, xmin, xmax, tick)
                println(io, @sprintf("<line class='grid' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x, y0, x, y0 + panel_height))
                println(io, @sprintf("<line class='tick' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x, y0 + panel_height, x, y0 + panel_height + 6.0))
                println(io, @sprintf("<text x='%.3f' y='%.3f' text-anchor='middle' font-size='10'>%s</text>",
                    x, y0 + panel_height + 20.0, _fmt(tick)))
            end
            for tick in y_ticks
                y = _svg_y(y0, panel_height, ymin, ymax, tick)
                println(io, @sprintf("<line class='grid' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x0, y, x0 + panel_width, y))
                println(io, @sprintf("<line class='tick' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x0 - 6.0, y, x0, y))
                println(io, @sprintf("<text x='%.3f' y='%.3f' text-anchor='end' dominant-baseline='middle' font-size='10'>%s</text>",
                    x0 - 10.0, y, _fmt(tick)))
            end
            if ymin < 0.0 < ymax
                y = _svg_y(y0, panel_height, ymin, ymax, 0.0)
                println(io, @sprintf("<line class='zero' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x0, y, x0 + panel_width, y))
            end
            println(io, @sprintf("<line class='axis' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                x0, y0 + panel_height, x0 + panel_width, y0 + panel_height))
            println(io, @sprintf("<line class='axis' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                x0, y0, x0, y0 + panel_height))
            println(io, @sprintf("<text x='%.3f' y='%.3f' text-anchor='middle' font-size='12'>fA</text>",
                x0 + panel_width / 2.0, y0 + panel_height + 42.0))
            println(io, @sprintf("<text transform='translate(%.3f %.3f) rotate(-90)' text-anchor='middle' font-size='12'>%s</text>",
                x0 - 52.0, y0 + panel_height / 2.0, panel.ylabel))

            for spec in specs
                panel.error && spec.id == "exact_rpa" && continue
                model_rows = [row for row in rows if row.model == spec.id]
                sort!(model_rows; by=row -> row.f)
                points = Tuple{Float64,Float64}[]
                for model_row in model_rows
                    value = _field_value(model_row, panel.field)
                    push!(points, (_svg_x(x0, panel_width, xmin, xmax, model_row.f),
                        _svg_y(y0, panel_height, ymin, ymax, value)))
                end
                dash = isempty(dash_by_model[spec.id]) ? "" :
                    " stroke-dasharray='$(dash_by_model[spec.id])'"
                println(io, "<polyline class='line' stroke='$(color_by_model[spec.id])'$(dash) points='$(_polyline(points))'/>")
            end
        end

        legend_x = width - 198.0
        legend_y = margin_top + 8.0
        println(io, @sprintf("<text x='%.3f' y='%.3f' font-size='13' font-weight='700'>Models</text>",
            legend_x, legend_y - 18.0))
        for (idx, spec) in enumerate(specs)
            y = legend_y + 23.0 * (idx - 1)
            dash = isempty(spec.dash) ? "" : " stroke-dasharray='$(spec.dash)'"
            println(io, @sprintf("<line class='line' stroke='%s'%s x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                spec.color, dash, legend_x, y, legend_x + 34.0, y))
            println(io, @sprintf("<text x='%.3f' y='%.3f' dominant-baseline='middle' font-size='11'>%s</text>",
                legend_x + 42.0, y, label_by_model[spec.id]))
        end
        println(io, @sprintf("<text x='%.3f' y='%.3f' font-size='11'>Error panels omit the exact</text>",
            legend_x, legend_y + 90.0))
        println(io, @sprintf("<text x='%.3f' y='%.3f' font-size='11'>RPA/BURP-TI curve.</text>",
            legend_x, legend_y + 106.0))
        println(io, "</svg>")
    end
    return path
end

function _write_note(path, rows, aggregate_rows, settings)
    open(path, "w") do io
        println(io, "# Diblock weak-response map")
        println(io)
        println(io, "This artifact extends the representative Gaussian-kernel panels into a continuous composition scan.  It is the lamellar/saddle analogue of the Uneyama-Doi weak-response and phase-boundary checks: the question is whether the real-density constructions reproduce the full diblock RPA spinodal and weak-segregation wavelength before any finite-amplitude saddle claim is made.")
        println(io)
        println(io, "The scan uses `fA <= 0.5` because AB exchange symmetry covers the other half of composition space.")
        println(io)
        println(io, "## Aggregate weak-response errors")
        println(io)
        println(io, "| model | max abs spinodal error | mean abs spinodal error | max abs period error | mean abs period error | boundary-limited minima | all converged |")
        println(io, "| --- | ---: | ---: | ---: | ---: | ---: | --- |")
        for row in aggregate_rows
            println(io, @sprintf("| %s | %s | %s | %s | %s | %d | `%s` |",
                row.model_label, _fmt(row.max_abs_spinodal_error),
                _fmt(row.mean_abs_spinodal_error), _fmt(row.max_abs_period_error),
                _fmt(row.mean_abs_period_error), row.boundary_limited_count,
                string(row.all_converged)))
        end
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, "1. BURP-TI is exactly the full diblock RPA kernel at Gaussian order for every scanned composition, so its spinodal and weak-period errors are numerically zero.")
        println(io, "2. The literal nonlinear Uneyama-Doi functional has a larger local Hessian than the nominal reduced vertex. The shift changes its spinodal but not its preferred weak-segregation wavevector.")
        println(io, "3. BVK1 and BVK2 share the nominal reduced kernel. Their finite-amplitude differences therefore arise beyond Gaussian order.")
        println(io, "4. Here OK denotes the Ohta--Kawasaki comparator with the published coefficients of Liu et al. Its local quadratic coefficient reproduces the RPA spinodal, but its gradient and nonlocal coefficients leave a composition-dependent error in the preferred weak-segregation period. This distinction helps separate instability-threshold matching from finite-amplitude accuracy.")
        println(io)
        println(io, "Outputs:")
        println(io)
        println(io, "- `weak_response_map.csv`: minima and relative errors for every model/composition pair.")
        println(io, "- `weak_response_summary.csv`: aggregate error metrics by model.")
        println(io, "- `diblock_weak_response_map.svg`: spinodal, weak-period, and error curves.")
        println(io)
        println(io, "Settings: `N=$(settings.N)`, `b=$(settings.b)`, `kRg_min=$(settings.qrg_min)`, `kRg_max=$(settings.qrg_max)`, `f_count=$(settings.f_count)`.")
    end
    return path
end

function main()
    outdir = normpath(_arg_value("--outdir",
        joinpath(MONTE_CARLO_PROJECT, "results", "diblock_weak_response_map")))
    mkpath(outdir)
    f_values = _f_values()
    qrg_min = parse(Float64, _arg_value("--qrg-min", "0.05"))
    qrg_max = parse(Float64, _arg_value("--qrg-max", "8.0"))
    N = parse(Float64, _arg_value("--N", "1.0"))
    b = parse(Float64, _arg_value("--b", "1.0"))
    0.0 < qrg_min < qrg_max || throw(ArgumentError("expected 0 < qRg-min < qRg-max"))

    specs = _model_specs()
    rows = NamedTuple[]
    for f in f_values
        exact = _kernel_minimum(specs[1].kernel, f, N, b, qrg_min, qrg_max)
        for spec in specs
            model_min = spec.id == "exact_rpa" ? exact :
                _kernel_minimum(spec.kernel, f, N, b, qrg_min, qrg_max)
            push!(rows, (
                f=f,
                N=N,
                b=b,
                model=spec.id,
                model_label=spec.label,
                qRg_star=model_min.qRg_star,
                period_rg=model_min.period_rg,
                kernel_star=model_min.kernel_star,
                chiN_spinodal=model_min.chiN_spinodal,
                exact_qRg_star=exact.qRg_star,
                exact_period_rg=exact.period_rg,
                exact_chiN_spinodal=exact.chiN_spinodal,
                relative_qRg_error=(model_min.qRg_star - exact.qRg_star) / exact.qRg_star,
                relative_period_error=(model_min.period_rg - exact.period_rg) /
                    exact.period_rg,
                relative_spinodal_error=(model_min.chiN_spinodal - exact.chiN_spinodal) /
                    exact.chiN_spinodal,
                converged=model_min.converged,
                iterations=model_min.iterations,
                boundary_limited=model_min.boundary_limited,
            ))
        end
    end

    aggregate_rows = NamedTuple[]
    for spec in specs
        model_rows = [row for row in rows if row.model == spec.id]
        push!(aggregate_rows, (
            model=spec.id,
            model_label=spec.label,
            f_count=length(model_rows),
            max_abs_spinodal_error=maximum(abs(row.relative_spinodal_error) for row in model_rows),
            mean_abs_spinodal_error=mean(abs(row.relative_spinodal_error) for row in model_rows),
            max_abs_period_error=maximum(abs(row.relative_period_error) for row in model_rows),
            mean_abs_period_error=mean(abs(row.relative_period_error) for row in model_rows),
            boundary_limited_count=count(row.boundary_limited for row in model_rows),
            all_converged=all(row.converged for row in model_rows),
        ))
    end

    map_csv = CSV.write(joinpath(outdir, "weak_response_map.csv"), DataFrame(rows))
    summary_csv = CSV.write(joinpath(outdir, "weak_response_summary.csv"),
        DataFrame(aggregate_rows))
    svg_path = _write_svg(joinpath(outdir, "diblock_weak_response_map.svg"),
        rows, f_values)
    note_path = _write_note(joinpath(outdir, "weak_response_map.md"), rows,
        aggregate_rows, (N=N, b=b, qrg_min=qrg_min, qrg_max=qrg_max,
            f_count=length(f_values)))

    println("Wrote $(map_csv)")
    println("Wrote $(summary_csv)")
    println("Wrote $(svg_path)")
    println("Wrote $(note_path)")
end

main()
