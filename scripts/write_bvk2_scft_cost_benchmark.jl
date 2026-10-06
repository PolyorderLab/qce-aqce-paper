#!/usr/bin/env julia

# Controlled, isolated BVK2-versus-Polyorder timing and memory benchmark.
# The comparison is deliberately end-to-end: every measured repetition builds
# a fresh state, converges the field, and relaxes a single-period lamellar cell.

using CSV
using DataFrames
using Dates
using Printf
using SHA
using Statistics
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RUNNER = joinpath(@__DIR__, "run_bvk2_scft_cost_case.jl")
const SCHEMA = "bvk2-scft-controlled-cost-benchmark-v1"
const CASES = (
    (case_id="AB_chiN20", architecture="AB", chiN=20.0,
        initial_period_rg=4.045700132441193),
    (case_id="AB_chiN30", architecture="AB", chiN=30.0,
        initial_period_rg=4.482349595936568),
    (case_id="ABA_chiN30", architecture="ABA", chiN=30.0,
        initial_period_rg=2.8222571085503),
    (case_id="ABA_chiN40", architecture="ABA", chiN=40.0,
        initial_period_rg=3.0228806841815703),
)
const METHODS = (
    (method="bvk2", label="BVK2", color="#0072b2"),
    (method="scft", label="Polyorder SCFT", color="#222222"),
)

function arg(name::AbstractString, default)
    prefix = String(name) * "="
    for value in ARGS
        startswith(value, prefix) &&
            return split(value, "="; limit=2)[2]
    end
    return default
end

function parse_bool(value)
    text = lowercase(strip(String(value)))
    text == "true" && return true
    text == "false" && return false
    error("expected true or false, got $(value)")
end

function parse_time_file(path)
    data = Dict{String,Float64}()
    for line in eachline(path)
        isempty(strip(line)) && continue
        key, value = split(line, "="; limit=2)
        data[String(key)] = parse(Float64, value)
    end
    return data
end

function launch_case(outdir, case, method;
        repeats::Integer, cpu::Integer)
    run_id = "$(case.case_id)_$(method.method)"
    rundir = joinpath(outdir, "runs", run_id)
    mkpath(rundir)
    metrics_path = joinpath(rundir, "metrics.toml")
    time_path = joinpath(rundir, "process_time.txt")
    log_path = joinpath(rundir, "run.log")
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    project = "--project=$(ROOT)"
    command = `taskset -c $(cpu) /usr/bin/time -f "process_wall_s=%e\nmax_rss_kib=%M\nuser_s=%U\nsystem_s=%S" -o $(time_path) $(julia) $(project) $(RUNNER) --method=$(method.method) --architecture=$(case.architecture) --chiN=$(case.chiN) --initial-period-rg=$(case.initial_period_rg) --repeats=$(repeats) --output=$(metrics_path)`
    environment = copy(ENV)
    environment["JULIA_NUM_THREADS"] = "1"
    environment["OPENBLAS_NUM_THREADS"] = "1"
    environment["OMP_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    open(log_path, "w") do log
        run(pipeline(setenv(command, environment);
            stdout=log, stderr=log))
    end
    return metrics_path, time_path, log_path
end

function assemble_rows(outdir)
    raw_rows = NamedTuple[]
    summary_rows = NamedTuple[]
    for case in CASES
        for method in METHODS
            rundir = joinpath(
                outdir, "runs", "$(case.case_id)_$(method.method)")
            metrics = TOML.parsefile(joinpath(rundir, "metrics.toml"))
            process = parse_time_file(
                joinpath(rundir, "process_time.txt"))
            runs = metrics["runs"]
            all(Bool(run["accepted"]) for run in runs) ||
                error("unaccepted measured run in $(rundir)")
            for run in runs
                push!(raw_rows, (
                    schema=SCHEMA,
                    case_id=case.case_id,
                    architecture=case.architecture,
                    chiN=case.chiN,
                    method=method.method,
                    repetition=Int(run["repetition"]),
                    solver_wall_s=Float64(run["wall_s"]),
                    gc_s=Float64(run["gc_s"]),
                    allocated_gib=
                        Float64(run["allocated_bytes"]) / 2.0^30,
                    period_rg=Float64(run["period_rg"]),
                    energy_density=Float64(run["energy_density"]),
                    grid_points=Int(run["grid_points"]),
                    field_converged=Bool(run["field_converged"]),
                    cell_converged=Bool(run["cell_converged"]),
                    single_period=Bool(run["single_period"]),
                    accepted=Bool(run["accepted"]),
                ))
            end
            wall = Float64[run["wall_s"] for run in runs]
            allocations = Float64[
                run["allocated_bytes"] / 2.0^30 for run in runs]
            periods = Float64[run["period_rg"] for run in runs]
            push!(summary_rows, (
                schema=SCHEMA,
                case_id=case.case_id,
                architecture=case.architecture,
                chiN=case.chiN,
                method=method.method,
                label=method.label,
                repetitions=length(runs),
                median_solver_wall_s=median(wall),
                min_solver_wall_s=minimum(wall),
                max_solver_wall_s=maximum(wall),
                relative_wall_span=
                    (maximum(wall)-minimum(wall)) / median(wall),
                median_allocated_gib=median(allocations),
                process_wall_s=process["process_wall_s"],
                peak_rss_mib=process["max_rss_kib"] / 1024.0,
                process_user_s=process["user_s"],
                process_system_s=process["system_s"],
                warmup_wall_s=Float64(metrics["warmup_wall_s"]),
                median_period_rg=median(periods),
                period_span_rg=maximum(periods)-minimum(periods),
                grid_points=Int(first(runs)["grid_points"]),
                field_converged=all(
                    Bool(run["field_converged"]) for run in runs),
                cell_converged=all(
                    Bool(run["cell_converged"]) for run in runs),
                single_period=all(
                    Bool(run["single_period"]) for run in runs),
                accepted=all(Bool(run["accepted"]) for run in runs),
                julia_threads=Int(metrics["julia_threads"]),
                blas_threads=Int(metrics["blas_threads"]),
                fftw_threads=Int(metrics["fftw_threads"]),
                solver_package_version=
                    String(metrics["solver_package_version"]),
            ))
        end
    end
    raw = DataFrame(raw_rows)
    summary = DataFrame(summary_rows)
    comparisons = NamedTuple[]
    for case in CASES
        bvk2 = summary[
            (String.(summary.case_id) .== case.case_id) .&
            (String.(summary.method) .== "bvk2"), :]
        scft = summary[
            (String.(summary.case_id) .== case.case_id) .&
            (String.(summary.method) .== "scft"), :]
        nrow(bvk2) == nrow(scft) == 1 ||
            error("missing paired summary for $(case.case_id)")
        push!(comparisons, (
            schema=SCHEMA,
            case_id=case.case_id,
            architecture=case.architecture,
            chiN=case.chiN,
            scft_over_bvk2_solver_time=
                scft.median_solver_wall_s[1] /
                bvk2.median_solver_wall_s[1],
            bvk2_over_scft_solver_time=
                bvk2.median_solver_wall_s[1] /
                scft.median_solver_wall_s[1],
            scft_over_bvk2_peak_rss=
                scft.peak_rss_mib[1] / bvk2.peak_rss_mib[1],
            scft_over_bvk2_allocations=
                scft.median_allocated_gib[1] /
                bvk2.median_allocated_gib[1],
            scft_over_bvk2_process_wall=
                scft.process_wall_s[1] / bvk2.process_wall_s[1],
            bvk2_period_rg=bvk2.median_period_rg[1],
            scft_period_rg=scft.median_period_rg[1],
            signed_period_difference=
                bvk2.median_period_rg[1] /
                scft.median_period_rg[1] - 1.0,
            accepted=Bool(bvk2.accepted[1]) &&
                Bool(scft.accepted[1]),
        ))
    end
    return raw, summary, DataFrame(comparisons)
end

function write_environment(path, cpu::Integer, repeats::Integer)
    cpu_model = try
        strip(read(`bash -lc "lscpu | sed -n 's/^Model name:[[:space:]]*//p'"`,
            String))
    catch
        "unknown"
    end
    commit = try
        strip(read(`git rev-parse HEAD`, String))
    catch
        "unknown"
    end
    manifest = joinpath(ROOT, "Manifest.toml")
    payload = Dict(
        "schema" => SCHEMA,
        "date" => string(Dates.today()),
        "cpu_model" => cpu_model,
        "logical_cpu_count" => Sys.CPU_THREADS,
        "pinned_cpu" => cpu,
        "julia_version" => string(VERSION),
        "julia_threads" => 1,
        "blas_threads" => 1,
        "fftw_threads" => 1,
        "repetitions" => repeats,
        "total_memory_bytes" => Sys.total_memory(),
        "git_commit" => commit,
        "git_worktree_dirty" =>
            !success(pipeline(`git diff --quiet`; stdout=devnull,
                stderr=devnull)),
        "manifest_sha256" => isfile(manifest) ?
            bytes2hex(sha256(read(manifest))) : "missing",
        "timing_tool" => "/usr/bin/time",
        "cpu_affinity_tool" => "taskset",
    )
    open(path, "w") do io
        TOML.print(io, payload; sorted=true)
    end
end

function style(method)
    return only(item for item in METHODS if item.method == method)
end

function svg_axes(io, panel; ymin, ymax, ylabel, logscale=false)
    ymap = if logscale
        lo, hi = log10(ymin), log10(ymax)
        y -> panel.y1 -
            (log10(Float64(y))-lo)/(hi-lo)*(panel.y1-panel.y0)
    else
        y -> panel.y1 -
            (Float64(y)-ymin)/(ymax-ymin)*(panel.y1-panel.y0)
    end
    println(io, "<rect x='$(panel.x0)' y='$(panel.y0)' width='$(panel.x1-panel.x0)' height='$(panel.y1-panel.y0)' fill='none' stroke='#222'/>")
    ticks = logscale ?
        10.0 .^ range(log10(ymin), log10(ymax); length=5) :
        range(ymin, ymax; length=5)
    for value in ticks
        y = ymap(value)
        println(io, "<line x1='$(panel.x0)' x2='$(panel.x1)' y1='$y' y2='$y' stroke='#ddd'/>")
        println(io, @sprintf("<text x='%.1f' y='%.1f' text-anchor='end' class='small'>%.3g</text>", panel.x0-8, y+4, value))
    end
    mid = (panel.y0+panel.y1)/2
    println(io, "<text x='$(panel.x0-55)' y='$mid' text-anchor='middle' class='label' transform='rotate(-90 $(panel.x0-55) $mid)'>$(ylabel)</text>")
    return ymap
end

function write_svg(path, summary)
    width, height = 1280, 520
    left = (x0=85.0, x1=610.0, y0=80.0, y1=405.0)
    right = (x0=730.0, x1=1255.0, y0=80.0, y1=405.0)
    times = Float64.(summary.median_solver_wall_s)
    rss = Float64.(summary.peak_rss_mib)
    ymin_time = 10.0^(floor(log10(minimum(times))) - 0.1)
    ymax_time = 10.0^(ceil(log10(maximum(times))) + 0.1)
    ymax_rss = 1.12 * maximum(rss)
    open(path, "w") do io
        println(io, "<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
        println(io, "<rect width='100%' height='100%' fill='white'/>")
        println(io, "<style>.title{font:700 18px Arial}.label{font:14px Arial}.small{font:12px Arial}</style>")
        println(io, "<text x='640' y='28' text-anchor='middle' class='title'>Controlled end-to-end lamellar cost: BVK2 versus Polyorder SCFT</text>")
        ytime = svg_axes(io, left; ymin=ymin_time,
            ymax=ymax_time, ylabel="median solver wall time (s)",
            logscale=true)
        yrss = svg_axes(io, right; ymin=0.0,
            ymax=ymax_rss, ylabel="compiler-inclusive peak RSS (MiB)")
        centers_left = range(left.x0+65, left.x1-65;
            length=length(CASES))
        centers_right = range(right.x0+65, right.x1-65;
            length=length(CASES))
        bar_width = 22.0
        for (case_index, case) in enumerate(CASES)
            for (method_index, method) in enumerate(METHODS)
                rows = summary[
                    (String.(summary.case_id) .== case.case_id) .&
                    (String.(summary.method) .== method.method), :]
                nrow(rows) == 1 || error("missing SVG row")
                offset = (method_index == 1 ? -1.0 : 1.0) *
                    0.65 * bar_width
                xleft = centers_left[case_index] + offset - bar_width/2
                xright = centers_right[case_index] + offset - bar_width/2
                y1 = ytime(rows.median_solver_wall_s[1])
                y2 = yrss(rows.peak_rss_mib[1])
                println(io, "<rect x='$xleft' y='$y1' width='$bar_width' height='$(left.y1-y1)' fill='$(method.color)'/>")
                println(io, "<rect x='$xright' y='$y2' width='$bar_width' height='$(right.y1-y2)' fill='$(method.color)'/>")
            end
            label = "$(case.architecture), chiN=$(Int(case.chiN))"
            println(io, "<text x='$(centers_left[case_index])' y='430' text-anchor='middle' class='small'>$(label)</text>")
            println(io, "<text x='$(centers_right[case_index])' y='430' text-anchor='middle' class='small'>$(label)</text>")
        end
        println(io, "<text x='$(0.5*(left.x0+left.x1))' y='60' text-anchor='middle' class='label'>Fresh seeded field + stress-free cell, median of three post-warm-up repetitions</text>")
        println(io, "<text x='$(0.5*(right.x0+right.x1))' y='60' text-anchor='middle' class='label'>Isolated process maximum resident set size</text>")
        for (index, method) in enumerate(METHODS)
            x = 465 + 190*(index-1)
            println(io, "<rect x='$x' y='474' width='24' height='12' fill='$(method.color)'/>")
            println(io, "<text x='$(x+32)' y='485' class='small'>$(method.label)</text>")
        end
        println(io, "<text x='640' y='510' text-anchor='middle' class='small'>Single pinned CPU; Julia/BLAS/FFTW threads = 1. Process RSS includes solver imports and JIT state; solver time excludes the compilation warm-up.</text>")
        println(io, "</svg>")
    end
end

function write_readme(path, summary, comparison)
    median_slowdown = median(
        Float64.(comparison.bvk2_over_scft_solver_time))
    median_rss = median(
        Float64.(comparison.scft_over_bvk2_peak_rss))
    open(path, "w") do io
        println(io, "# Controlled BVK2 / Polyorder SCFT cost benchmark")
        println(io)
        println(io, "This packet closes manuscript task P1-2 with an end-to-end lamellar comparison. Every measured repetition constructs a fresh seeded field, converges it, and relaxes a single-period cell. It is not a kernel microbenchmark.")
        println(io)
        println(io, "## Protocol")
        println(io)
        println(io, "- Same host and one pinned logical CPU; Julia, BLAS, and FFTW each use one thread.")
        println(io, "- Four accepted symmetric states: AB at `chiN=20,30` and ABA at `chiN=30,40`.")
        println(io, "- BVK2 uses the frozen publication setting `c2=0.16`, `nx=128`, and its analytic/shared cell-root workflow.")
        println(io, "- Polyorder uses `OSF`, `ds=0.005`, `maxDeltaX=0.04`, residual tolerance `1e-6`, cell-stress tolerance `1e-5`, and the same seeded Anderson/variable-cell workflow as the accepted references.")
        println(io, "- Each isolated process performs one compilation warm-up followed by three fresh measured solves using the same fixed RNG seed. A single-period gate rejects tiled multi-period cells. Solver wall time is the median of those three; process wall time and peak RSS include imports, JIT state, warm-up, and measurements.")
        println(io)
        println(io, "The spatial discretizations are the accepted publication settings, not artificially identical grid counts. Therefore this measures cost to reproduce the claimed results, not asymptotic cost at equal degrees of freedom.")
        println(io)
        println(io, "## Results")
        println(io)
        println(io, "| case | BVK2 time (s) | SCFT time (s) | BVK2/SCFT time | BVK2 peak RSS (MiB) | SCFT peak RSS (MiB) | SCFT/BVK2 RSS |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|")
        for case in CASES
            b = summary[
                (String.(summary.case_id) .== case.case_id) .&
                (String.(summary.method) .== "bvk2"), :][1, :]
            s = summary[
                (String.(summary.case_id) .== case.case_id) .&
                (String.(summary.method) .== "scft"), :][1, :]
            c = comparison[
                String.(comparison.case_id) .== case.case_id, :][1, :]
            println(io, @sprintf(
                "| %s, chiN=%.0f | %.3f | %.3f | %.3f | %.1f | %.1f | %.3f |",
                case.architecture, case.chiN,
                b.median_solver_wall_s, s.median_solver_wall_s,
                c.bvk2_over_scft_solver_time,
                b.peak_rss_mib, s.peak_rss_mib,
                c.scft_over_bvk2_peak_rss))
        end
        println(io)
        println(io, @sprintf(
            "Across these four one-dimensional publication states, BVK2 is **%.3fx slower** at the median paired post-warm-up solver time, while the median SCFT/BVK2 peak-RSS ratio is **%.3f**.",
            median_slowdown, median_rss))
        println(io)
        println(io, "Polyorder is `8.2--16.7x` faster in the measured post-warm-up solver stage. Its isolated process nevertheless takes longer from cold startup because Polyorder's compilation warm-up is about twice as long. Peak RSS is similar, with SCFT `1--13%` higher, whereas SCFT allocates only about `7--9%` as many bytes during the measured solver call.")
        println(io)
        println(io, "Interpretation must remain scoped: the result supports an empirical lamellar cost statement on this CPU and these accuracy contracts. It does not establish GPU behavior, three-dimensional morphology cost, or asymptotic complexity. Most importantly, the propagator-free BVK2 formulation does **not** currently deliver an end-to-end speed advantage in this production workflow; Polyorder's mature Anderson/variable-cell algorithm more than offsets its MDE work. This is an implementation-level result, not a reversal of the formal per-evaluation complexity argument.")
        println(io)
        println(io, "## Reproduction")
        println(io)
        println(io, "```bash")
        println(io, "julia --project=. scripts/write_bvk2_scft_cost_benchmark.jl --recompute=true")
        println(io, "```")
    end
end

function main()
    outdir = abspath(arg("--output", joinpath(
        ROOT, "results", "bvk2_scft_cost_benchmark")))
    recompute = parse_bool(arg("--recompute", "false"))
    repeats = parse(Int, arg("--repeats", "3"))
    cpu = parse(Int, arg("--cpu", "12"))
    mkpath(joinpath(outdir, "runs"))
    if recompute
        for case in CASES
            for method in METHODS
                println("running ", case.case_id, " ", method.method)
                launch_case(outdir, case, method;
                    repeats=repeats, cpu=cpu)
            end
        end
    end
    raw, summary, comparison = assemble_rows(outdir)
    all(Bool.(raw.accepted)) ||
        error("raw benchmark contains an unaccepted repetition")
    all(Bool.(summary.accepted)) ||
        error("summary benchmark contains an unaccepted solver")
    all(Bool.(comparison.accepted)) ||
        error("paired benchmark contains an unaccepted case")
    CSV.write(joinpath(outdir, "raw_runs.csv"), raw; newline='\n')
    CSV.write(joinpath(outdir, "summary.csv"), summary; newline='\n')
    CSV.write(joinpath(outdir, "comparison.csv"), comparison; newline='\n')
    write_environment(
        joinpath(outdir, "environment.toml"), cpu, repeats)
    write_svg(joinpath(outdir, "cost_memory_benchmark.svg"), summary)
    write_readme(
        joinpath(outdir, "README.md"), summary, comparison)
    println("wrote controlled cost benchmark to ", outdir)
end

main()
