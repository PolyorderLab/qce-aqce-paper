#!/usr/bin/env julia

# Frozen-closure architecture-transfer packet for the Macromolecules P1-1
# claim. BVK2 keeps the diblock-calibrated c2=0.16; only the ideal-chain ABA
# coefficients are changed. Polyorder SCFT is an independently relaxed
# reference and never enters the BVK2 parameter construction.

using CSV
using DataFrames
using DFMMonteCarlo
using LinearAlgebra
using Logging
using Polyorder
using Polyorder: Anderson, BB, SD, VariableCell
using Printf
using Random
using Statistics
using TOML

const C2 = 0.16
const CASES = ((f=0.5, chiN=30.0), (f=0.5, chiN=40.0))
const NX = 128
const MODES = 16
const SAMPLE_COUNT = 256
const PERIOD_RELATIVE_TOLERANCE = 0.05
const PROFILE_RMS_TOLERANCE = 0.05
const PROFILE_CORRELATION_TOLERANCE = 0.98
const SCFT_STRESS_TOLERANCE = 1.0e-4

function _density_symbol()
    return Symbol(Char(0x03d5))
end

function _save_density_symbol()
    return Symbol("save_", Char(0x03d5))
end

function _polyorder_io_config(outdir)
    kwargs = Dict{Symbol,Any}(
        :verbosity => -1,
        :progress_solve => false,
        :progress_cell => false,
        :base_dir => outdir,
        :save_config => false,
        :save_w => false,
        :save_summary => false,
        :save_trace => false,
    )
    kwargs[_save_density_symbol()] = false
    return Polyorder.IOConfig(; kwargs...)
end

function _polyorder_density(scft, specie::Symbol)
    density_function = getproperty(Polyorder, _density_symbol())
    return vec(getfield(density_function(scft, 1, specie), :data))
end

function _aba_system(f::Real, chiN::Real)
    sA, sB = KuhnSegment(:A), KuhnSegment(:B)
    eb1, eb2 = BranchPoint(:EB1), BranchPoint(:EB2)
    a1 = PolymerBlock(:A1, sA, Float64(f) / 2, FreeEnd(:A1), eb1)
    b = PolymerBlock(:B, sB, 1.0 - Float64(f), eb1, eb2)
    a2 = PolymerBlock(:A2, sA, Float64(f) / 2, eb2, FreeEnd(:A2))
    chain = BlockCopolymer(:ABA, [a1, b, a2])
    return PolymerSystem(
        [Component(chain, 1.0, 1.0)],
        Dict([:A, :B] => Float64(chiN)),
    )
end

function _scft_profile(f::Real, chiN::Real, initial_period_rg::Real;
        spacing::Real=0.04, ds::Real=0.005, seed::Integer=3083)
    return mktempdir() do workdir
        cd(workdir) do
            io = _polyorder_io_config(workdir)
            field_config = Polyorder.Config(;
                io=io,
                scft=Polyorder.SCFTConfig(;
                    min_iter=10,
                    max_iter=4000,
                    tolmode=:Residual,
                    tol=1.0e-6,
                    maxΔx=Float64(spacing),
                ),
                cellopt=Polyorder.CellOptConfig(;
                    changeNx=true,
                    tol_stress=1.0e-5,
                ),
            )
            system = _aba_system(f, chiN)
            lattice = BravaisLattice(UnitCell(Float64(initial_period_rg)))
            scft = NoncyclicChainSCFT(
                system,
                lattice,
                Float64(ds);
                mde=OSF,
                spacing=Float64(spacing),
                updater=Anderson(SD(1.0); warmup=150, αw=0.2),
                init=:randn,
                rng=Random.Xoshiro(seed),
            )
            field_status = with_logger(NullLogger()) do
                Polyorder.solve!(scft, field_config)
            end
            cell_updater = VariableCell(
                BB(1.0),
                Anderson(SD(1.0); warmup=80, αw=0.2);
                block=30,
            )
            cell_status, _ = with_logger(NullLogger()) do
                Polyorder.cell_solve!(scft, cell_updater, field_config)
            end
            phi_a = _polyorder_density(scft, :A)
            phi_b = _polyorder_density(scft, :B)
            stress = try
                norm(Polyorder.gradient_wrt_cell(scft))
            catch
                NaN
            end
            period = Polyorder.unitcell(scft).edges[1]
            return (
                phi_a=phi_a,
                phi_b=phi_b,
                period_rg=Float64(period),
                free_energy=Float64(Polyorder.F(scft)),
                field_status=string(field_status),
                cell_status=string(cell_status),
                stress=Float64(stress),
                grid_count=length(phi_a),
                mean_phi_a=mean(phi_a),
                incompressibility_rms=
                    sqrt(mean(abs2, phi_a .+ phi_b .- 1.0)),
            )
        end
    end
end

function _periodic_linear(values::AbstractVector{<:Real}, s::Real)
    n = length(values)
    x = mod(Float64(s), 1.0) * n
    left0 = floor(Int, x)
    fraction = x - left0
    left = mod1(left0 + 1, n)
    right = mod1(left + 1, n)
    return (1.0 - fraction) * Float64(values[left]) +
        fraction * Float64(values[right])
end

function _resample(values, count::Integer=SAMPLE_COUNT)
    return [_periodic_linear(values, (index - 1) / count)
        for index in 1:count]
end

function _align(reference, model)
    ref = _resample(reference)
    candidate = _resample(model)
    best = (rms=Inf, correlation=-Inf, shift=0, values=candidate)
    for shift in 0:(SAMPLE_COUNT - 1)
        aligned = circshift(candidate, shift)
        rms = sqrt(mean(abs2, aligned .- ref))
        correlation = cor(ref, aligned)
        if rms < best.rms
            best = (
                rms=rms,
                correlation=correlation,
                shift=shift,
                values=aligned,
            )
        end
    end
    return ref, best
end

function _svg_polyline(xs, ys, xmap, ymap)
    return join(
        (@sprintf("%.3f,%.3f", xmap(x), ymap(y)) for (x, y) in zip(xs, ys)),
        " ",
    )
end

function _write_svg(path, summaries, profiles)
    width, height = 1240, 430
    left = (x0=75.0, x1=390.0, y0=70.0, y1=355.0)
    panels = [
        (x0=470.0, x1=800.0, y0=70.0, y1=355.0),
        (x0=870.0, x1=1200.0, y0=70.0, y1=355.0),
    ]
    chis = Float64.(summaries.chiN)
    periods = vcat(Float64.(summaries.bvk2_period_rg),
        Float64.(summaries.scft_period_rg))
    ymin, ymax = minimum(periods) - 0.08, maximum(periods) + 0.08
    xmin, xmax = minimum(chis) - 1.0, maximum(chis) + 1.0
    xmap = x -> left.x0 + (Float64(x) - xmin) / (xmax - xmin) *
        (left.x1 - left.x0)
    ymap = y -> left.y1 - (Float64(y) - ymin) / (ymax - ymin) *
        (left.y1 - left.y0)
    open(path, "w") do io
        println(io, "<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
        println(io, "<rect width='100%' height='100%' fill='white'/>")
        println(io, "<style>.t{font:16px Arial}.s{font:13px Arial}.a{stroke:#333;stroke-width:1;fill:none}.g{stroke:#ddd;stroke-width:1}.b{stroke:#1f77b4;stroke-width:2.5;fill:none}.r{stroke:#d62728;stroke-width:2.5;fill:none}</style>")
        println(io, "<text x='620' y='28' text-anchor='middle' class='t' font-weight='bold'>No-refit ABA transfer: frozen BVK2 c2=0.16 versus SCFT</text>")
        for panel in (left, panels...)
            println(io, "<rect x='$(panel.x0)' y='$(panel.y0)' width='$(panel.x1-panel.x0)' height='$(panel.y1-panel.y0)' class='a'/>")
        end
        for y in range(ymin, ymax; length=5)
            println(io, "<line x1='$(left.x0)' x2='$(left.x1)' y1='$(ymap(y))' y2='$(ymap(y))' class='g'/>")
            println(io, @sprintf("<text x='%.1f' y='%.1f' text-anchor='end' class='s'>%.2f</text>", left.x0-8, ymap(y)+5, y))
        end
        for chi in chis
            println(io, @sprintf("<text x='%.1f' y='%.1f' text-anchor='middle' class='s'>%.0f</text>", xmap(chi), left.y1+22, chi))
        end
        println(io, "<text x='$(0.5*(left.x0+left.x1))' y='395' text-anchor='middle' class='s'>χN</text>")
        println(io, "<text x='24' y='$(0.5*(left.y0+left.y1))' text-anchor='middle' class='s' transform='rotate(-90 24 $(0.5*(left.y0+left.y1)))'>stress-free D/Rg</text>")
        println(io, "<polyline points='$(_svg_polyline(chis, summaries.bvk2_period_rg, xmap, ymap))' class='b'/>")
        println(io, "<polyline points='$(_svg_polyline(chis, summaries.scft_period_rg, xmap, ymap))' class='r'/>")
        for row in eachrow(summaries)
            println(io, "<circle cx='$(xmap(row.chiN))' cy='$(ymap(row.bvk2_period_rg))' r='4' fill='#1f77b4'/>")
            println(io, "<circle cx='$(xmap(row.chiN))' cy='$(ymap(row.scft_period_rg))' r='4' fill='#d62728'/>")
        end
        println(io, "<line x1='95' x2='125' y1='88' y2='88' class='b'/><text x='132' y='93' class='s'>BVK2 ABA</text>")
        println(io, "<line x1='95' x2='125' y1='108' y2='108' class='r'/><text x='132' y='113' class='s'>SCFT ABA</text>")

        for (panel, row) in zip(panels, eachrow(summaries))
            subset = profiles[profiles.case_id .== row.case_id, :]
            ss = Float64.(subset.s)
            bvk = Float64.(subset.bvk2_phi_a)
            scft = Float64.(subset.scft_phi_a)
            pmapx = s -> panel.x0 + Float64(s) * (panel.x1 - panel.x0)
            pmapy = p -> panel.y1 - Float64(p) * (panel.y1 - panel.y0)
            for y in 0.0:0.25:1.0
                println(io, "<line x1='$(panel.x0)' x2='$(panel.x1)' y1='$(pmapy(y))' y2='$(pmapy(y))' class='g'/>")
            end
            println(io, "<polyline points='$(_svg_polyline(ss, bvk, pmapx, pmapy))' class='b'/>")
            println(io, "<polyline points='$(_svg_polyline(ss, scft, pmapx, pmapy))' class='r'/>")
            println(io, @sprintf("<text x='%.1f' y='%.1f' text-anchor='middle' class='s'>χN=%.0f; RMS=%.3f; r=%.4f</text>", 0.5*(panel.x0+panel.x1), panel.y0-12, row.chiN, row.profile_rms, row.profile_correlation))
            println(io, "<text x='$(0.5*(panel.x0+panel.x1))' y='395' text-anchor='middle' class='s'>normalized cell coordinate</text>")
        end
        println(io, "<text x='835' y='418' text-anchor='middle' class='s'>Each profile uses its independently optimized stress-free cell.</text>")
        println(io, "</svg>")
    end
    return path
end

function _write_readme(path, summaries, fingerprint)
    max_period_error = maximum(abs, summaries.period_relative_error)
    max_profile_rms = maximum(summaries.profile_rms)
    min_correlation = minimum(summaries.profile_correlation)
    accepted = all(summaries.accepted)
    open(path, "w") do io
        println(io, "# Frozen no-refit BVK2 transfer to a symmetric ABA triblock")
        println(io)
        println(io, "This packet tests P1-1 of the Macromolecules plan. The BVK2 nonlinear closure is the diblock-calibrated `c2 = 0.16`; only the ideal-chain correlation and its analytically derived ABA Gaussian coefficients change. No ABA SCFT result is used in parameter construction or fitting.")
        println(io)
        println(io, "Model fingerprint: `", fingerprint, "`")
        println(io)
        println(io, "| χN | BVK2 D/Rg | SCFT D/Rg | signed error | profile RMS | correlation | accepted |")
        println(io, "|---:|---:|---:|---:|---:|---:|:---:|")
        for row in eachrow(summaries)
            println(io, @sprintf("| %.0f | %.6f | %.6f | %+.3f%% | %.4f | %.5f | %s |",
                row.chiN, row.bvk2_period_rg, row.scft_period_rg,
                100 * row.period_relative_error, row.profile_rms,
                row.profile_correlation, row.accepted ? "yes" : "no"))
        end
        println(io)
        println(io, "Predeclared gates: `|ΔD|/D_SCFT ≤ 5%`, aligned profile RMS `≤ 0.05`, correlation `≥ 0.98`, converged/local-minimum BVK2 fields, and SCFT cell stress `≤ 1e-4`.")
        println(io)
        println(io, @sprintf("Result: **%s**. Maximum period error is %.3f%%, maximum profile RMS is %.4f, and minimum correlation is %.5f.",
            accepted ? "accepted" : "not accepted", 100 * max_period_error,
            max_profile_rms, min_correlation))
        println(io)
        println(io, "Scope: this establishes quantitative no-refit structural transfer for symmetric linear ABA lamellae at two strong-segregation state points. It does not establish an ABA phase diagram, arbitrary-topology universality, or Frank–Kasper stabilization.")
        println(io)
        println(io, "Reproduce with `julia --project=. scripts/write_bvk2_aba_no_refit_validation.jl`.")
    end
    return path
end

function main()
    outdir = normpath(joinpath(
        @__DIR__, "..", "results", "bvk2_aba_no_refit_validation"))
    mkpath(outdir)
    fingerprint = triblock_aba_bvk2_fingerprint(; f=0.5, c2=C2)
    selected_model = TOML.parsefile(joinpath(
        @__DIR__, "..", "results", "bvk2_c2_calibration_nx128",
        "selected_model.toml"))
    Float64(selected_model["c2"]) == C2 ||
        error("selected diblock BVK2 c2 is not the frozen value 0.16")

    summaries = NamedTuple[]
    profile_rows = NamedTuple[]
    for (case_index, case) in enumerate(CASES)
        case_id = @sprintf("f%.2f_chiN%.0f", case.f, case.chiN)
        println("[$case_id] relaxing ABA BVK2")
        bvk2 = minimize_triblock_aba_bvk2_lamella_stress_free(;
            f=case.f, chiN=case.chiN, nx=NX, mode_count=MODES,
            initial_amplitudes=(0.35,), c2=C2,
            max_iterations=4000, max_period_iterations=16,
            bootstrap_period_factors=(1.15, 1.25),
            period_strategy=:bootstrap_local, bootstrap_window=0.12,
            local_check_fraction=0.01,
            progress_label="ABA-BVK2 $case_id")
        println("[$case_id] relaxing independent ABA SCFT")
        scft = _scft_profile(
            case.f, case.chiN, bvk2.period_rg;
            seed=3083 + case_index,
        )
        scft_ref, aligned = _align(scft.phi_a, bvk2.result.phi_a)
        period_error = (bvk2.period_rg - scft.period_rg) / scft.period_rg
        bvk2_pass = bvk2.converged && bvk2.local_minimum_check_pass &&
            bvk2.result.converged
        scft_pass = isfinite(scft.stress) &&
            scft.stress <= SCFT_STRESS_TOLERANCE &&
            occursin("Successful", scft.cell_status)
        accepted = bvk2_pass && scft_pass &&
            abs(period_error) <= PERIOD_RELATIVE_TOLERANCE &&
            aligned.rms <= PROFILE_RMS_TOLERANCE &&
            aligned.correlation >= PROFILE_CORRELATION_TOLERANCE
        push!(summaries, (
            case_id=case_id,
            fA=case.f,
            chiN=case.chiN,
            c2=C2,
            bvk2_period_rg=bvk2.period_rg,
            scft_period_rg=scft.period_rg,
            period_relative_error=period_error,
            profile_rms=aligned.rms,
            profile_correlation=aligned.correlation,
            profile_shift=aligned.shift,
            bvk2_energy_density=bvk2.objective,
            scft_free_energy=scft.free_energy,
            bvk2_force_r2=bvk2.result.projected_force_norm,
            bvk2_force_rinf=bvk2.result.projected_force_maxabs,
            bvk2_local_minimum=bvk2.local_minimum_check_pass,
            bvk2_converged=bvk2.result.converged,
            scft_field_status=scft.field_status,
            scft_cell_status=scft.cell_status,
            scft_cell_stress=scft.stress,
            scft_incompressibility_rms=scft.incompressibility_rms,
            bvk2_grid=NX,
            scft_grid=scft.grid_count,
            accepted=accepted,
            model_fingerprint=fingerprint,
        ))
        for index in 1:SAMPLE_COUNT
            push!(profile_rows, (
                case_id=case_id,
                fA=case.f,
                chiN=case.chiN,
                s=(index - 1) / SAMPLE_COUNT,
                bvk2_phi_a=aligned.values[index],
                scft_phi_a=scft_ref[index],
            ))
        end
    end

    summary = DataFrame(summaries)
    profiles = DataFrame(profile_rows)
    CSV.write(joinpath(outdir, "summary.csv"), summary; newline='\n')
    CSV.write(joinpath(outdir, "profiles.csv"), profiles; newline='\n')
    _write_svg(joinpath(outdir, "aba_no_refit_validation.svg"),
        summary, profiles)
    _write_readme(joinpath(outdir, "README.md"), summary, fingerprint)
    all(summary.accepted) ||
        error("one or more ABA no-refit validation cases failed")
    println("accepted ", nrow(summary), " / ", nrow(summary), " cases")
end

if get(ENV, "DFM_SKIP_BVK2_ABA_VALIDATION_MAIN", "0") != "1"
    main()
end
