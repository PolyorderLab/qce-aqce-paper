const MONTE_CARLO_PROJECT = normpath(joinpath(@__DIR__, ".."))
MONTE_CARLO_PROJECT in LOAD_PATH || pushfirst!(LOAD_PATH, MONTE_CARLO_PROJECT)

using CSV
using DataFrames
using DFMMonteCarlo
using Printf
using Statistics

const ULB_SCHEMA = "unified-own-cell-lamella-v1"
const ULB_CASES = [
    (case_id="f0.5_chiN12", f=0.50, chiN=12.0, bvk2_split="holdout"),
    (case_id="f0.5_chiN15", f=0.50, chiN=15.0, bvk2_split="train"),
    (case_id="f0.5_chiN20", f=0.50, chiN=20.0, bvk2_split="holdout"),
    (case_id="f0.5_chiN25", f=0.50, chiN=25.0, bvk2_split="train"),
    (case_id="f0.5_chiN30", f=0.50, chiN=30.0, bvk2_split="holdout"),
    (case_id="f0.35_chiN30", f=0.35, chiN=30.0, bvk2_split="train"),
]
const ULB_MODELS = [
    (model="scft", label="SCFT", color="#222222", dash="7 4"),
    (model="ohta_kawasaki", label="Ohta-Kawasaki", color="#7a4fb0", dash=""),
    (model="uneyama_doi", label="Uneyama-Doi", color="#2a8c47", dash=""),
    (model="liu2019_opf", label="Liu-OPF", color="#2878b5", dash=""),
    (model="burp_ti", label="BURP-TI", color="#c44e22", dash=""),
    (model="bvk2", label="BVK2", color="#e69f00", dash=""),
]
const ULB_APPROXIMATION_MODELS = [style for style in ULB_MODELS
    if style.model != "scft"]

function ulb_arg_value(name::AbstractString, default)
    prefix = String(name) * "="
    for argument in ARGS
        startswith(argument, prefix) &&
            return split(argument, "=", limit=2)[2]
    end
    return default
end

function ulb_periodic_linear(values::AbstractVector{<:Real}, coordinate::Real)
    count = length(values)
    count > 0 || throw(ArgumentError("cannot interpolate an empty profile"))
    scaled = mod(Float64(coordinate), 1.0) * count
    base = floor(Int, scaled)
    left = mod1(base + 1, count)
    right = mod1(left + 1, count)
    fraction = scaled - base
    return (1.0 - fraction) * Float64(values[left]) +
        fraction * Float64(values[right])
end

function ulb_resample(values::AbstractVector{<:Real}, sample_count::Integer=256)
    count = Int(sample_count)
    count > 1 || throw(ArgumentError("sample_count must exceed one"))
    return [ulb_periodic_linear(values, (index - 1) / count)
        for index in 1:count]
end

function ulb_align_periodic(reference::AbstractVector{<:Real},
        model::AbstractVector{<:Real}; sample_count::Integer=256,
        allow_swap::Bool=false)
    target = ulb_resample(reference, sample_count)
    model_values = ulb_resample(model, sample_count)
    candidates = allow_swap ? (model_values, 1.0 .- model_values) :
        (model_values,)
    best = (rms=Inf, shift=0, swapped=false, values=copy(model_values))
    for (candidate_index, candidate) in enumerate(candidates)
        for shift in 0:(sample_count - 1)
            aligned = [candidate[mod1(index + shift, sample_count)]
                for index in 1:sample_count]
            rms = sqrt(mean(abs2, aligned .- target))
            if rms < best.rms
                best = (rms=rms, shift=shift,
                    swapped=candidate_index == 2, values=aligned)
            end
        end
    end
    return target, best
end

function ulb_profile_statistics(reference, model)
    target = Float64.(reference)
    values = Float64.(model)
    length(target) == length(values) ||
        throw(DimensionMismatch("profile vectors must have the same length"))
    residual_sum = sum(abs2, values .- target)
    centered_sum = sum(abs2, target .- mean(target))
    return (
        rms=sqrt(residual_sum / length(target)),
        r2=centered_sum > 0.0 ? 1.0 - residual_sum / centered_sum : 1.0,
        correlation=cor(target, values),
        minimum=minimum(values),
        maximum=maximum(values),
        mean=mean(values),
    )
end

function ulb_case_profiles(profiles::DataFrame, case_id::AbstractString,
        model::AbstractString)
    rows = profiles[(String.(profiles.case_id) .== String(case_id)) .&
        (String.(profiles.model) .== String(model)), :]
    nrow(rows) > 1 || throw(ArgumentError(
        "missing $(model) profile for $(case_id)"))
    sort!(rows, :s)
    return Float64.(rows.phi_a)
end

function ulb_selected_bvk2_rows(training::DataFrame, holdout::DataFrame)
    combined = vcat(training, holdout; cols=:union)
    selected = combined[
        isapprox.(Float64.(combined.c2), 0.16; atol=1.0e-12, rtol=0.0) .&
        Bool.(combined.selection_eligible), :]
    case_ids = Set(case.case_id for case in ULB_CASES)
    selected = selected[in.(String.(selected.case_id), Ref(case_ids)), :]
    for case in ULB_CASES
        count(String.(selected.case_id) .== case.case_id) == 1 ||
            throw(ArgumentError("expected one accepted BVK2 row for $(case.case_id)"))
    end
    sort!(selected, :chiN)
    return selected
end

function ulb_validate_bvk2_cache(cache::DataFrame, selected::DataFrame;
        sample_count::Integer=256)
    all(name -> name in propertynames(cache),
        (:case_id, :model, :period_rg, :s, :phi_a, :nx, :c2)) || return false
    for row in eachrow(selected)
        rows = cache[(String.(cache.case_id) .== String(row.case_id)) .&
            (String.(cache.model) .== "bvk2"), :]
        nrow(rows) == sample_count || return false
        all(Int.(rows.nx) .== 128) || return false
        all(isapprox.(Float64.(rows.c2), 0.16; atol=1.0e-12, rtol=0.0)) ||
            return false
        all(isapprox.(Float64.(rows.period_rg), Float64(row.model_period_rg);
            atol=1.0e-10, rtol=0.0)) || return false
    end
    return true
end

function ulb_materialize_bvk2_profiles(base_summary::DataFrame,
        base_profiles::DataFrame, selected::DataFrame, cache_path::AbstractString;
        sample_count::Integer=256, reuse_cache::Bool=true)
    if reuse_cache && isfile(cache_path)
        cached = CSV.read(cache_path, DataFrame)
        ulb_validate_bvk2_cache(cached, selected; sample_count=sample_count) &&
            return cached
    end

    rows = NamedTuple[]
    for selected_row in eachrow(selected)
        case_id = String(selected_row.case_id)
        case = only(case for case in ULB_CASES if case.case_id == case_id)
        reference_rows = base_summary[String.(base_summary.case_id) .== case_id, :]
        nrow(reference_rows) > 0 ||
            throw(ArgumentError("missing reference summary for $(case_id)"))
        reference_row = reference_rows[1, :]
        scft_profile = ulb_case_profiles(base_profiles, case_id, "scft")
        n_value = Float64(reference_row.N)
        b_value = Float64(reference_row.b)
        period_rg = Float64(selected_row.model_period_rg)
        rg_scale = sqrt(6.0 / (n_value * b_value^2))
        period_b = period_rg / rg_scale

        coarse = minimize_diblock_bvk2_lamella(; f=case.f, chiN=case.chiN,
            L=period_b, nx=64, mode_count=3,
            initial_profile=ulb_resample(scft_profile, 64),
            N=n_value, b=b_value, c2=0.16, max_iterations=4_000,
            gradient_tolerance=1.0e-10, force_tolerance=1.0e-5,
            force_maxabs_tolerance=2.0e-5)
        coarse.converged || throw(ErrorException(
            "BVK2 coarse profile replay failed for $(case_id)"))
        fine = minimize_diblock_bvk2_lamella(; f=case.f, chiN=case.chiN,
            L=period_b, nx=128, mode_count=3,
            initial_profile=ulb_resample(coarse.phi_a, 128),
            N=n_value, b=b_value, c2=0.16, max_iterations=4_000,
            gradient_tolerance=1.0e-10, force_tolerance=1.0e-5,
            force_maxabs_tolerance=2.0e-5)
        fine.converged || throw(ErrorException(
            "BVK2 fine profile replay failed for $(case_id)"))
        target, alignment = ulb_align_periodic(scft_profile, fine.phi_a;
            sample_count=sample_count, allow_swap=isapprox(case.f, 0.5;
                atol=1.0e-12))
        abs(alignment.rms - Float64(selected_row.profile_rms)) <= 2.0e-4 ||
            throw(ErrorException(@sprintf(
                "BVK2 replay RMS mismatch for %s: %.8g vs frozen %.8g",
                case_id, alignment.rms, Float64(selected_row.profile_rms))))
        for index in 1:sample_count
            push!(rows, (
                schema=ULB_SCHEMA,
                case_id=case_id,
                f=case.f,
                chiN=case.chiN,
                bvk2_split=case.bvk2_split,
                model="bvk2",
                period_rg=period_rg,
                s=(index - 1) / sample_count,
                phi_a=alignment.values[index],
                nx=128,
                c2=0.16,
                alignment_shift_fraction=alignment.shift / sample_count,
                alignment_swapped=alignment.swapped,
            ))
        end
    end
    cache = DataFrame(rows)
    CSV.write(cache_path, cache)
    return cache
end

function ulb_base_summary_row(row, case, model_label)
    model_amplitude = Float64(row.model_maximum_phi_a) -
        Float64(row.model_minimum_phi_a)
    scft_amplitude = Float64(row.scft_maximum_phi_a) -
        Float64(row.scft_minimum_phi_a)
    signed_period_error = Float64(row.model_period_rg) /
        Float64(row.scft_period_rg) - 1.0
    field_gate = Bool(row.model_converged)
    cell_gate = Bool(row.period_local_minimum_check_pass) &&
        !Bool(row.period_boundary_limited)
    return (
        schema=ULB_SCHEMA, case_id=case.case_id, f=case.f, chiN=case.chiN,
        bvk2_split=case.bvk2_split, model=String(row.model),
        model_label=model_label, model_calibration_role="not_applicable",
        period_rg=Float64(row.model_period_rg),
        scft_period_rg=Float64(row.scft_period_rg),
        signed_period_error=signed_period_error,
        abs_period_error=abs(signed_period_error),
        profile_rms=Float64(row.profile_rms),
        profile_r2=Float64(row.profile_r2),
        profile_correlation=Float64(row.profile_correlation),
        alignment_shift_fraction=Float64(row.alignment_shift_fraction),
        alignment_swapped=Bool(row.alignment_swapped),
        min_phi_a=Float64(row.model_minimum_phi_a),
        max_phi_a=Float64(row.model_maximum_phi_a),
        mean_phi_a=Float64(row.model_mean_phi_a),
        amplitude_error=abs(model_amplitude - scft_amplitude),
        field_gate_pass=field_gate, cell_gate_pass=cell_gate,
        accepted=field_gate && cell_gate, status_reason=field_gate && cell_gate ?
            "accepted" : "source_gate_failed",
        source_artifact="results/diblock_stress_free_lamella_comparison/summary.csv",
        source_fingerprint="canonical-own-cell-$(String(row.model))",
    )
end

function ulb_summary_rows(base_summary::DataFrame, base_profiles::DataFrame,
        opf_summary::DataFrame, selected_bvk2::DataFrame,
        bvk2_profiles::DataFrame)
    rows = NamedTuple[]
    for case in ULB_CASES
        source_rows = base_summary[String.(base_summary.case_id) .== case.case_id, :]
        nrow(source_rows) > 0 ||
            throw(ArgumentError("missing base case $(case.case_id)"))
        reference = source_rows[1, :]
        scft_field = occursin("Successful", String(reference.scft_convergence))
        scft_cell = Float64(reference.scft_stress_norm) <= 1.0e-5
        push!(rows, (
            schema=ULB_SCHEMA, case_id=case.case_id, f=case.f, chiN=case.chiN,
            bvk2_split=case.bvk2_split, model="scft", model_label="SCFT",
            model_calibration_role="reference",
            period_rg=Float64(reference.scft_period_rg),
            scft_period_rg=Float64(reference.scft_period_rg),
            signed_period_error=0.0, abs_period_error=0.0, profile_rms=0.0,
            profile_r2=1.0, profile_correlation=1.0,
            alignment_shift_fraction=0.0, alignment_swapped=false,
            min_phi_a=Float64(reference.scft_minimum_phi_a),
            max_phi_a=Float64(reference.scft_maximum_phi_a),
            mean_phi_a=Float64(reference.scft_mean_phi_a),
            amplitude_error=0.0, field_gate_pass=scft_field,
            cell_gate_pass=scft_cell, accepted=scft_field && scft_cell,
            status_reason=scft_field && scft_cell ? "reference" :
                "scft_gate_failed",
            source_artifact="results/diblock_stress_free_lamella_comparison/summary.csv",
            source_fingerprint="Polyorder-cell-solved-SCFT",
        ))

        for model in ("ohta_kawasaki", "uneyama_doi", "burp_ti")
            matches = source_rows[String.(source_rows.model) .== model, :]
            nrow(matches) == 1 ||
                throw(ArgumentError("missing $(model) row for $(case.case_id)"))
            label = only(style.label for style in ULB_MODELS
                if style.model == model)
            push!(rows, ulb_base_summary_row(matches[1, :], case, label))
        end

        opf_matches = opf_summary[String.(opf_summary.case_id) .== case.case_id, :]
        nrow(opf_matches) == 1 ||
            throw(ArgumentError("missing in-domain OPF row for $(case.case_id)"))
        opf = opf_matches[1, :]
        scft_amplitude = Float64(opf.scft_maximum_phi) -
            Float64(opf.scft_minimum_phi)
        opf_amplitude = Float64(opf.opf_maximum_phi) -
            Float64(opf.opf_minimum_phi)
        opf_field = Bool(opf.solver_converged)
        opf_cell = Bool(opf.period_local_minimum_check_pass) &&
            !Bool(opf.period_boundary_limited)
        push!(rows, (
            schema=ULB_SCHEMA, case_id=case.case_id, f=case.f, chiN=case.chiN,
            bvk2_split=case.bvk2_split, model="liu2019_opf",
            model_label="Liu-OPF", model_calibration_role="not_applicable",
            period_rg=Float64(opf.opf_period_rg),
            scft_period_rg=Float64(opf.scft_period_rg),
            signed_period_error=Float64(opf.signed_period_error),
            abs_period_error=abs(Float64(opf.signed_period_error)),
            profile_rms=Float64(opf.profile_rms),
            profile_r2=Float64(opf.profile_r2),
            profile_correlation=Float64(opf.profile_correlation),
            alignment_shift_fraction=Float64(opf.alignment_shift_fraction),
            alignment_swapped=Bool(opf.alignment_swapped),
            min_phi_a=Float64(opf.opf_minimum_phi),
            max_phi_a=Float64(opf.opf_maximum_phi),
            mean_phi_a=Float64(opf.opf_mean_phi),
            amplitude_error=abs(opf_amplitude - scft_amplitude),
            field_gate_pass=opf_field, cell_gate_pass=opf_cell,
            accepted=Bool(opf.accepted),
            status_reason=Bool(opf.accepted) ? "accepted" : "source_gate_failed",
            source_artifact="results/liu2019_opf_baseline/summary.csv",
            source_fingerprint="Liu2019-Eqs20-24-Table2",
        ))

        bvk2_matches = selected_bvk2[
            String.(selected_bvk2.case_id) .== case.case_id, :]
        nrow(bvk2_matches) == 1 ||
            throw(ArgumentError("missing selected BVK2 row for $(case.case_id)"))
        bvk2 = bvk2_matches[1, :]
        bvk2_phi = ulb_case_profiles(bvk2_profiles, case.case_id, "bvk2")
        scft_phi = ulb_case_profiles(base_profiles, case.case_id, "scft")
        profile_stats = ulb_profile_statistics(scft_phi, bvk2_phi)
        abs(profile_stats.rms - Float64(bvk2.profile_rms)) <= 2.0e-4 ||
            throw(ErrorException("BVK2 cached profile does not reproduce frozen RMS"))
        scft_amplitude = Float64(reference.scft_maximum_phi_a) -
            Float64(reference.scft_minimum_phi_a)
        bvk2_amplitude = profile_stats.maximum - profile_stats.minimum
        push!(rows, (
            schema=ULB_SCHEMA, case_id=case.case_id, f=case.f, chiN=case.chiN,
            bvk2_split=case.bvk2_split, model="bvk2", model_label="BVK2",
            model_calibration_role=case.bvk2_split,
            period_rg=Float64(bvk2.model_period_rg),
            scft_period_rg=Float64(bvk2.scft_period_rg),
            signed_period_error=Float64(bvk2.signed_period_error),
            abs_period_error=Float64(bvk2.period_abs_error),
            profile_rms=Float64(bvk2.profile_rms),
            profile_r2=profile_stats.r2,
            profile_correlation=profile_stats.correlation,
            alignment_shift_fraction=Float64(bvk2.alignment_shift_fraction),
            alignment_swapped=Bool(bvk2.alignment_swapped),
            min_phi_a=profile_stats.minimum, max_phi_a=profile_stats.maximum,
            mean_phi_a=profile_stats.mean,
            amplitude_error=abs(bvk2_amplitude - scft_amplitude),
            field_gate_pass=Bool(bvk2.stationarity_pass),
            cell_gate_pass=Bool(bvk2.cell_stationarity_pass),
            accepted=Bool(bvk2.selection_eligible),
            status_reason=Bool(bvk2.selection_eligible) ? "accepted_frozen_c2" :
                "source_gate_failed",
            source_artifact=case.bvk2_split == "train" ?
                "results/bvk2_c2_calibration_nx128/stage2_rows.csv" :
                "results/bvk2_c2_calibration_holdout_nx128/stage2_rows.csv",
            source_fingerprint=String(bvk2.fingerprint),
        ))
    end
    model_order = Dict(style.model => index for (index, style) in
        enumerate(ULB_MODELS))
    case_order = Dict(case.case_id => index for (index, case) in
        enumerate(ULB_CASES))
    sort!(rows; by=row -> (case_order[row.case_id], model_order[row.model]))
    return rows
end

function ulb_aggregate_metrics(summary_rows)
    rows = [row for row in summary_rows if row.model != "scft" && row.accepted]
    aggregates = NamedTuple[]
    for style in ULB_APPROXIMATION_MODELS
        model_rows = [row for row in rows if row.model == style.model]
        length(model_rows) == length(ULB_CASES) || throw(ArgumentError(
            "$(style.model) does not cover the six-case common panel"))
        period_wins = 0
        profile_wins = 0
        for case in ULB_CASES
            case_rows = [row for row in rows if row.case_id == case.case_id]
            model_row = only(row for row in case_rows if row.model == style.model)
            isapprox(model_row.abs_period_error,
                minimum(row.abs_period_error for row in case_rows);
                atol=1.0e-14) && (period_wins += 1)
            isapprox(model_row.profile_rms,
                minimum(row.profile_rms for row in case_rows);
                atol=1.0e-14) && (profile_wins += 1)
        end
        push!(aggregates, (
            model=style.model, model_label=style.label,
            case_count=length(model_rows),
            mean_abs_period_error=mean(row.abs_period_error for row in model_rows),
            max_abs_period_error=maximum(row.abs_period_error for row in model_rows),
            mean_profile_rms=mean(row.profile_rms for row in model_rows),
            max_profile_rms=maximum(row.profile_rms for row in model_rows),
            period_wins=period_wins, profile_wins=profile_wins,
        ))
    end
    return aggregates
end

function ulb_unified_profiles(base_profiles::DataFrame,
        opf_profiles::DataFrame, bvk2_profiles::DataFrame,
        summary_rows; sample_count::Integer=256)
    period_by_key = Dict((row.case_id, row.model) => row.period_rg
        for row in summary_rows)
    rows = NamedTuple[]
    for case in ULB_CASES
        for model in ("scft", "ohta_kawasaki", "uneyama_doi", "burp_ti")
            values = ulb_case_profiles(base_profiles, case.case_id, model)
            values = ulb_resample(values, sample_count)
            for index in 1:sample_count
                push!(rows, (
                    schema=ULB_SCHEMA, case_id=case.case_id, f=case.f,
                    chiN=case.chiN, bvk2_split=case.bvk2_split, model=model,
                    period_rg=period_by_key[(case.case_id, model)],
                    s=(index - 1) / sample_count, phi_a=values[index],
                ))
            end
        end
        opf_values = ulb_case_profiles(opf_profiles, case.case_id,
            "liu2019_opf")
        for index in 1:sample_count
            push!(rows, (
                schema=ULB_SCHEMA, case_id=case.case_id, f=case.f,
                chiN=case.chiN, bvk2_split=case.bvk2_split,
                model="liu2019_opf",
                period_rg=period_by_key[(case.case_id, "liu2019_opf")],
                s=(index - 1) / sample_count, phi_a=opf_values[index],
            ))
        end
        cached = bvk2_profiles[
            (String.(bvk2_profiles.case_id) .== case.case_id) .&
            (String.(bvk2_profiles.model) .== "bvk2"), :]
        sort!(cached, :s)
        nrow(cached) == sample_count ||
            throw(ArgumentError("invalid BVK2 cache for $(case.case_id)"))
        for row in eachrow(cached)
            push!(rows, (
                schema=ULB_SCHEMA, case_id=case.case_id, f=case.f,
                chiN=case.chiN, bvk2_split=case.bvk2_split, model="bvk2",
                period_rg=period_by_key[(case.case_id, "bvk2")],
                s=Float64(row.s), phi_a=Float64(row.phi_a),
            ))
        end
    end
    return rows
end

ulb_svg_x(x0, width, value) = x0 + width * Float64(value)
ulb_svg_y(y0, height, value, ymin, ymax) =
    y0 + height * (1.0 - (Float64(value) - ymin) / (ymax - ymin))

function ulb_polyline(xs, ys)
    return join((@sprintf("%.3f,%.3f", x, y) for (x, y) in zip(xs, ys)), " ")
end

function ulb_write_profiles_svg(path::AbstractString, profiles::DataFrame)
    panel_width, panel_height = 370.0, 220.0
    left, top, gap_x, gap_y = 74.0, 94.0, 54.0, 64.0
    width = left + 2panel_width + gap_x + 34.0
    height = top + 3panel_height + 2gap_y + 74.0
    open(path, "w") do io
        println(io, @sprintf(
            "<svg xmlns='http://www.w3.org/2000/svg' width='%.0f' height='%.0f' viewBox='0 0 %.0f %.0f'>",
            width, height, width, height))
        println(io, "<rect width='100%' height='100%' fill='white'/>")
        println(io, "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202020}.axis{stroke:#222;stroke-width:1.1}.grid{stroke:#e2e2e2;stroke-width:.8}.curve{fill:none;stroke-width:2;stroke-linejoin:round;stroke-linecap:round}</style>")
        println(io, "<text x='$(width/2)' y='28' text-anchor='middle' font-size='18' font-weight='700'>Stress-free lamellar profiles in each model's own cell</text>")
        legend_y = 57.0
        legend_x = 45.0
        for (index, style) in enumerate(ULB_MODELS)
            x = legend_x + (index - 1) * 136.0
            dash = isempty(style.dash) ? "" :
                " stroke-dasharray='$(style.dash)'"
            println(io, @sprintf(
                "<line class='curve' stroke='%s'%s x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f'/>",
                style.color, dash, x, legend_y, x + 25, legend_y))
            println(io, @sprintf(
                "<text x='%.1f' y='%.1f' dominant-baseline='middle' font-size='11'>%s</text>",
                x + 31, legend_y, style.label))
        end
        for (case_index, case) in enumerate(ULB_CASES)
            column = mod(case_index - 1, 2)
            row_index = div(case_index - 1, 2)
            x0 = left + column * (panel_width + gap_x)
            y0 = top + row_index * (panel_height + gap_y)
            for tick in 0.0:0.25:1.0
                x = ulb_svg_x(x0, panel_width, tick)
                y = ulb_svg_y(y0, panel_height, tick, 0.0, 1.0)
                println(io, @sprintf(
                    "<line class='grid' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x, y0, x, y0 + panel_height))
                println(io, @sprintf(
                    "<line class='grid' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x0, y, x0 + panel_width, y))
            end
            println(io, @sprintf(
                "<line class='axis' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                x0, y0 + panel_height, x0 + panel_width, y0 + panel_height))
            println(io, @sprintf(
                "<line class='axis' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                x0, y0, x0, y0 + panel_height))
            println(io, @sprintf(
                "<text x='%.3f' y='%.3f' font-size='13' font-weight='700'>fA=%.2f, chiN=%.0f</text>",
                x0 + 8, y0 + 18, case.f, case.chiN))
            badge = case.bvk2_split == "train" ? "BVK2 train" :
                "BVK2 held-out"
            println(io, @sprintf(
                "<text x='%.3f' y='%.3f' text-anchor='end' font-size='10' fill='#555'>%s</text>",
                x0 + panel_width - 8, y0 + 18, badge))
            for style in ULB_MODELS
                rows = profiles[
                    (String.(profiles.case_id) .== case.case_id) .&
                    (String.(profiles.model) .== style.model), :]
                sort!(rows, :s)
                xs = [ulb_svg_x(x0, panel_width, s) for s in Float64.(rows.s)]
                ys = [ulb_svg_y(y0, panel_height, phi, 0.0, 1.0)
                    for phi in Float64.(rows.phi_a)]
                dash = isempty(style.dash) ? "" :
                    " stroke-dasharray='$(style.dash)'"
                println(io, "<polyline class='curve' stroke='$(style.color)'$(dash) points='$(ulb_polyline(xs,ys))'/>")
            end
            println(io, @sprintf(
                "<text x='%.3f' y='%.3f' text-anchor='middle' font-size='11'>x / L_model</text>",
                x0 + panel_width / 2, y0 + panel_height + 30))
            println(io, @sprintf(
                "<text transform='translate(%.3f %.3f) rotate(-90)' text-anchor='middle' font-size='11'>phi_A</text>",
                x0 - 42, y0 + panel_height / 2))
        end
        println(io, "</svg>")
    end
    return path
end

function ulb_marker(io, x, y, color, model, split)
    if model == "bvk2"
        half = 4.8
        points = @sprintf("%.3f,%.3f %.3f,%.3f %.3f,%.3f %.3f,%.3f",
            x, y-half, x-half, y, x, y+half, x+half, y)
        fill = split == "train" ? color : "white"
        println(io, "<polygon points='$(points)' fill='$(fill)' stroke='$(color)' stroke-width='1.8'/>")
    else
        println(io, @sprintf(
            "<circle cx='%.3f' cy='%.3f' r='4.1' fill='%s' stroke='white' stroke-width='1'/>",
            x, y, color))
    end
end

function ulb_write_errors_svg(path::AbstractString, summary_rows)
    rows = [row for row in summary_rows if row.model != "scft"]
    panel_width, panel_height = 500.0, 330.0
    left, top, gap = 80.0, 78.0, 92.0
    width = left + 2panel_width + gap + 42.0
    height = top + panel_height + 132.0
    period_values = 100.0 .* [row.signed_period_error for row in rows]
    period_min = floor((minimum(period_values) - 1.0) / 5.0) * 5.0
    period_max = ceil((maximum(period_values) + 1.0) / 5.0) * 5.0
    profile_max = ceil(maximum(row.profile_rms for row in rows) * 100.0) / 100.0
    metrics = [
        (field=:signed_period_error, scale=100.0, ymin=period_min,
            ymax=period_max, title="Signed period error", ylabel="100(L/L_SCFT - 1)"),
        (field=:profile_rms, scale=1.0, ymin=0.0, ymax=profile_max,
            title="Normalized-profile RMS", ylabel="RMS vs SCFT"),
    ]
    model_index = Dict(style.model => index for (index, style) in
        enumerate(ULB_APPROXIMATION_MODELS))
    open(path, "w") do io
        println(io, @sprintf(
            "<svg xmlns='http://www.w3.org/2000/svg' width='%.0f' height='%.0f' viewBox='0 0 %.0f %.0f'>",
            width, height, width, height))
        println(io, "<rect width='100%' height='100%' fill='white'/>")
        println(io, "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202020}.axis{stroke:#222;stroke-width:1.1}.grid{stroke:#e0e0e0;stroke-width:.8}.zero{stroke:#555;stroke-width:1.1;stroke-dasharray:5 4}</style>")
        println(io, "<text x='$(width/2)' y='28' text-anchor='middle' font-size='18' font-weight='700'>Six-model own-cell lamellar errors</text>")
        for (metric_index, metric) in enumerate(metrics)
            x0 = left + (metric_index - 1) * (panel_width + gap)
            y0 = top
            for fraction in 0.0:0.25:1.0
                value = metric.ymin + fraction * (metric.ymax - metric.ymin)
                y = ulb_svg_y(y0, panel_height, value, metric.ymin, metric.ymax)
                println(io, @sprintf(
                    "<line class='grid' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x0, y, x0 + panel_width, y))
                println(io, @sprintf(
                    "<text x='%.3f' y='%.3f' text-anchor='end' dominant-baseline='middle' font-size='10'>%.3g</text>",
                    x0 - 8, y, value))
            end
            if metric.ymin <= 0.0 <= metric.ymax
                zero_y = ulb_svg_y(y0, panel_height, 0.0, metric.ymin,
                    metric.ymax)
                println(io, @sprintf(
                    "<line class='zero' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                    x0, zero_y, x0 + panel_width, zero_y))
            end
            println(io, @sprintf(
                "<line class='axis' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                x0, y0 + panel_height, x0 + panel_width, y0 + panel_height))
            println(io, @sprintf(
                "<line class='axis' x1='%.3f' y1='%.3f' x2='%.3f' y2='%.3f'/>",
                x0, y0, x0, y0 + panel_height))
            println(io, @sprintf(
                "<text x='%.3f' y='%.3f' text-anchor='middle' font-size='14' font-weight='700'>%s</text>",
                x0 + panel_width / 2, y0 - 20, metric.title))
            for (case_index, case) in enumerate(ULB_CASES)
                center = x0 + panel_width * (case_index - 0.5) /
                    length(ULB_CASES)
                println(io, @sprintf(
                    "<text x='%.3f' y='%.3f' text-anchor='middle' font-size='9'>%.2f/%g</text>",
                    center, y0 + panel_height + 18, case.f, case.chiN))
                for style in ULB_APPROXIMATION_MODELS
                    row = only(row for row in rows
                        if row.case_id == case.case_id &&
                           row.model == style.model)
                    offset = (model_index[style.model] -
                        (length(ULB_APPROXIMATION_MODELS) + 1) / 2) * 7.5
                    value = Float64(getproperty(row, metric.field)) *
                        metric.scale
                    y = ulb_svg_y(y0, panel_height, value, metric.ymin,
                        metric.ymax)
                    ulb_marker(io, center + offset, y, style.color,
                        style.model, case.bvk2_split)
                end
            end
            println(io, @sprintf(
                "<text transform='translate(%.3f %.3f) rotate(-90)' text-anchor='middle' font-size='11'>%s</text>",
                x0 - 52, y0 + panel_height / 2, metric.ylabel))
            println(io, @sprintf(
                "<text x='%.3f' y='%.3f' text-anchor='middle' font-size='10'>fA / chiN</text>",
                x0 + panel_width / 2, y0 + panel_height + 39))
        end
        legend_y = height - 42.0
        legend_x = 95.0
        for (index, style) in enumerate(ULB_APPROXIMATION_MODELS)
            x = legend_x + (index - 1) * 175.0
            ulb_marker(io, x, legend_y, style.color, style.model, "train")
            println(io, @sprintf(
                "<text x='%.3f' y='%.3f' dominant-baseline='middle' font-size='11'>%s</text>",
                x + 11, legend_y, style.label))
        end
        println(io, @sprintf(
            "<text x='%.3f' y='%.3f' font-size='10'>BVK2: filled diamond=train; open diamond=held-out. Split applies only to BVK2.</text>",
            width - 430, height - 15))
        println(io, "</svg>")
    end
    return path
end

function ulb_write_readme(path, aggregates, summary_path, profiles_path,
        profile_svg, error_svg)
    open(path, "w") do io
        println(io, "# Unified stress-free lamellar benchmark")
        println(io)
        println(io, "This is the main-text six-model comparison for SCFT, ")
        println(io, "Ohta--Kawasaki, Uneyama--Doi, Liu2019 OPF, BURP-TI, and ")
        println(io, "the frozen BVK2 model (`c2=0.16`).")
        println(io)
        println(io, "Every profile is solved in its model's own stress-free cell. ")
        println(io, "Profiles are mapped to `s=x/L_model`, periodically aligned to ")
        println(io, "SCFT, and sampled at 256 points. Period and profile errors are ")
        println(io, "therefore separate observables; the overlay does not hide period ")
        println(io, "errors by imposing the SCFT cell.")
        println(io)
        println(io, "The common panel contains six state points inside Liu et al.'s ")
        println(io, "published regression domain (`0.2<=f_A<=0.8`, ")
        println(io, "`chiN_s(f_A)<chiN<=35`). The two existing `chiN=50` SCFT/BVK2 ")
        println(io, "cases are excluded rather than extrapolating OPF.")
        println(io)
        println(io, "| model | cases | mean absolute period error | maximum absolute period error | mean profile RMS | maximum profile RMS | period wins | profile wins |")
        println(io, "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for row in aggregates
            println(io, @sprintf(
                "| %s | %d | %.5f | %.5f | %.5f | %.5f | %d | %d |",
                row.model_label, row.case_count, row.mean_abs_period_error,
                row.max_abs_period_error, row.mean_profile_rms,
                row.max_profile_rms, row.period_wins, row.profile_wins))
        end
        println(io)
        println(io, "BVK2 training cases are `(0.50,15)`, `(0.50,25)`, and ")
        println(io, "`(0.35,30)`; the other three common cases are held out. ")
        println(io, "The split labels apply only to BVK2. BVK1 is excluded from this ")
        println(io, "main comparison because it is the mechanism-discovery precursor; ")
        println(io, "its eight-case ablation remains in the earlier benchmark.")
        println(io)
        println(io, "Outputs:")
        println(io)
        println(io, "- `$(basename(summary_path))`: 36 accepted case/model rows.")
        println(io, "- `$(basename(profiles_path))`: aligned own-cell profiles.")
        println(io, "- `$(basename(profile_svg))`: six-panel profile figure.")
        println(io, "- `$(basename(error_svg))`: period/profile error figure.")
        println(io)
        println(io, "Reproduce with:")
        println(io)
        println(io, "```bash")
        println(io, "julia --project=. scripts/write_liu2019_opf_baseline.jl")
        println(io, "julia --project=. scripts/write_unified_lamellar_benchmark.jl")
        println(io, "```")
    end
    return path
end

function main()
    base_dir = normpath(ulb_arg_value("--base-dir", joinpath(
        MONTE_CARLO_PROJECT, "results",
        "diblock_stress_free_lamella_comparison")))
    opf_dir = normpath(ulb_arg_value("--opf-dir", joinpath(
        MONTE_CARLO_PROJECT, "results", "liu2019_opf_baseline")))
    training_dir = normpath(ulb_arg_value("--bvk2-training-dir", joinpath(
        MONTE_CARLO_PROJECT, "results", "bvk2_c2_calibration_nx128")))
    holdout_dir = normpath(ulb_arg_value("--bvk2-holdout-dir", joinpath(
        MONTE_CARLO_PROJECT, "results",
        "bvk2_c2_calibration_holdout_nx128")))
    outdir = normpath(ulb_arg_value("--outdir", joinpath(
        MONTE_CARLO_PROJECT, "results", "unified_lamellar_benchmark")))
    sample_count = parse(Int, ulb_arg_value("--sample-count", "256"))
    reuse_cache = lowercase(ulb_arg_value("--reuse-bvk2-profiles", "true")) in
        ("true", "1", "yes")
    mkpath(outdir)

    base_summary = CSV.read(joinpath(base_dir, "summary.csv"), DataFrame)
    base_profiles = CSV.read(joinpath(base_dir, "profiles.csv"), DataFrame)
    opf_summary = CSV.read(joinpath(opf_dir, "summary.csv"), DataFrame)
    opf_profiles = CSV.read(joinpath(opf_dir, "profiles.csv"), DataFrame)
    selected_bvk2 = ulb_selected_bvk2_rows(
        CSV.read(joinpath(training_dir, "stage2_rows.csv"), DataFrame),
        CSV.read(joinpath(holdout_dir, "stage2_rows.csv"), DataFrame))
    bvk2_cache_path = joinpath(outdir, "bvk2_profiles_nx128.csv")
    bvk2_profiles = ulb_materialize_bvk2_profiles(base_summary, base_profiles,
        selected_bvk2, bvk2_cache_path; sample_count=sample_count,
        reuse_cache=reuse_cache)

    summary_rows = ulb_summary_rows(base_summary, base_profiles, opf_summary,
        selected_bvk2, bvk2_profiles)
    all(row.accepted for row in summary_rows) ||
        throw(ErrorException("one or more unified benchmark rows failed gates"))
    profiles = ulb_unified_profiles(base_profiles, opf_profiles,
        bvk2_profiles, summary_rows; sample_count=sample_count)
    aggregates = ulb_aggregate_metrics(summary_rows)

    summary_path = CSV.write(joinpath(outdir, "summary.csv"),
        DataFrame(summary_rows))
    profiles_path = CSV.write(joinpath(outdir, "profiles.csv"),
        DataFrame(profiles))
    profile_svg = ulb_write_profiles_svg(joinpath(outdir,
        "stress_free_lamella_profiles.svg"), DataFrame(profiles))
    error_svg = ulb_write_errors_svg(joinpath(outdir,
        "stress_free_lamella_errors.svg"), summary_rows)
    readme = ulb_write_readme(joinpath(outdir, "README.md"), aggregates,
        summary_path, profiles_path, profile_svg, error_svg)

    println("Wrote ", summary_path)
    println("Wrote ", profiles_path)
    println("Wrote ", profile_svg)
    println("Wrote ", error_svg)
    println("Wrote ", readme)
end

if get(ENV, "DFM_SKIP_UNIFIED_LAMELLAR_MAIN", "0") != "1"
    main()
end
