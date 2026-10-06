#!/usr/bin/env julia

# Fixed-cell performance benchmarks for BVK2 and Polyorder SCFT.
#
# Two questions are kept separate:
# 1. Kernel cost: one BVK2 energy+chemical-potential evaluation versus one
#    Polyorder q! evaluation (MDE propagation+density integration+force) at
#    equal full-grid M.
# 2. Warm continuation: the marginal cost of an adjacent-chiN fixed-cell
#    solve seeded from the immediately preceding accepted state.

using CSV
using DataFrames
using Dates
using DFMMonteCarlo
using FFTW
using LinearAlgebra
using Logging
using Polyorder
using Printf
using Random
using Statistics
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_OUTDIR =
    joinpath(ROOT, "results", "bvk2_scft_fixedcell_performance")
const SCHEMA = "bvk2-scft-fixedcell-performance-v1"
const C2 = 0.16
const DS = 0.005
const RNG_SEED = 7301
const BVK2_WARM_MAX_ITERATIONS = 50
const KERNEL_1D_GRIDS = (32, 64, 128, 256, 512, 1024)
const KERNEL_3D_GRIDS = (12, 16, 24, 32, 40, 48)
const CONTINUATION_CASES = (
    (architecture="AB", chi_anchor=20.0, chi_values=collect(21.0:30.0),
        period_rg=4.2),
    (architecture="ABA", chi_anchor=30.0, chi_values=collect(31.0:40.0),
        period_rg=2.9),
)

function arg(name::AbstractString, default)
    prefix = String(name) * "="
    for value in ARGS
        startswith(value, prefix) &&
            return split(value, "="; limit=2)[2]
    end
    return default
end

function aba_system(f::Real, chiN::Real)
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

polymer_system(architecture::AbstractString, chiN::Real) =
    architecture == "AB" ?
    AB_system(χN=Float64(chiN), fA=0.5) :
    aba_system(0.5, Float64(chiN))

function quiet_io(outdir)
    save_density = Symbol("save_", Char(0x03d5))
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
    kwargs[save_density] = false
    return Polyorder.IOConfig(; kwargs...)
end

function scft_config(outdir; tol=1.0e-6)
    return Polyorder.Config(;
        io=quiet_io(outdir),
        scft=Polyorder.SCFTConfig(;
            min_iter=2,
            max_iter=4000,
            tolmode=:Residual,
            tol=Float64(tol),
            maxΔx=0.04,
        ),
        cellopt=Polyorder.CellOptConfig(; changeNx=false),
    )
end

function exact_template(period_rg::Real, dims::Tuple)
    lat = length(dims) == 1 ?
        BravaisLattice(UnitCell(Float64(period_rg))) :
        BravaisLattice(UnitCell(Cubic(), Float64(period_rg)))
    return AuxiliaryField(zeros(Float64, dims...), lat)
end

function make_kernel_scft(dims::Tuple; chiN=30.0, period_rg=4.2)
    template = exact_template(period_rg, dims)
    scft = NoncyclicChainSCFT(
        AB_system(χN=Float64(chiN), fA=0.5),
        template,
        DS;
        mde=Polyorder.OSF,
        updater=Polyorder.SD(0.2),
        profile=:fast,
        use_dct=false,
        symmetrize=false,
        init=:none,
    )
    nx = dims[1]
    phase = [2.0 * pi * (i - 0.5) / nx for i in 1:nx]
    if length(dims) == 1
        scft.wfields[1].data .= 2.0 .* cos.(phase)
        scft.wfields[2].data .= -2.0 .* cos.(phase)
    else
        for index in CartesianIndices(scft.wfields[1].data)
            value = 2.0 * cos(phase[index[1]])
            scft.wfields[1].data[index] = value
            scft.wfields[2].data[index] = -value
        end
    end
    Polyorder.q!(scft)
    Polyorder.q!(scft)
    return scft
end

function kernel_density(dims::Tuple)
    nx = dims[1]
    phase = [2.0 * pi * (i - 0.5) / nx for i in 1:nx]
    if length(dims) == 1
        return 0.5 .+ 0.35 .* cos.(phase)
    end
    phi = Array{Float64}(undef, dims)
    for index in CartesianIndices(phi)
        phi[index] = 0.5 + 0.35 * cos(phase[index[1]])
    end
    return phi
end

function bvk2_kernel(phi, lengths_b)
    if ndims(phi) == 1
        energy = diblock_bvk2_energy_1d(
            phi; f=0.5, chiN=30.0, L=only(lengths_b), c2=C2)
        chemical = diblock_bvk2_chemical_potential_1d(
            phi; f=0.5, chiN=30.0, L=only(lengths_b), c2=C2)
    else
        energy = diblock_bvk2_energy_nd(
            phi; f=0.5, chiN=30.0, lengths=lengths_b, c2=C2)
        chemical = diblock_bvk2_chemical_potential_nd(
            phi; f=0.5, chiN=30.0, lengths=lengths_b, c2=C2)
    end
    return Float64(energy) + sum(chemical)
end

function median_measurement(call, repeats::Integer)
    call()
    call()
    rows = NamedTuple[]
    for repetition in 1:repeats
        GC.gc()
        measurement = @timed call()
        push!(rows, (
            repetition=repetition,
            wall_s=measurement.time,
            allocated_bytes=measurement.bytes,
            checksum=Float64(measurement.value),
        ))
    end
    return rows
end

function run_kernel_benchmark()
    raw = NamedTuple[]
    for (dimension, grids, repeats) in (
            (1, KERNEL_1D_GRIDS, 9),
            (3, KERNEL_3D_GRIDS, 5))
        for n in grids
            dims = ntuple(_ -> n, dimension)
            m = prod(dims)
            phi = kernel_density(dims)
            period_rg = 4.2
            lengths_b = fill(period_rg / sqrt(6.0), dimension)
            scft = make_kernel_scft(dims; period_rg=period_rg)
            bvk2_rows = median_measurement(
                () -> bvk2_kernel(phi, lengths_b), repeats)
            scft_rows = median_measurement(() -> begin
                Polyorder.q!(scft)
                return Float64(Polyorder.F(scft)) +
                    sum(field -> sum(field.data), scft.forces)
            end, repeats)
            for (method, rows) in (("bvk2", bvk2_rows), ("scft", scft_rows))
                for row in rows
                    push!(raw, (
                        schema=SCHEMA,
                        dimension=dimension,
                        linear_grid=n,
                        grid_points=m,
                        method=method,
                        repetition=row.repetition,
                        wall_s=row.wall_s,
                        allocated_bytes=row.allocated_bytes,
                        checksum=row.checksum,
                        contour_step=method == "scft" ? DS : NaN,
                        contour_steps=method == "scft" ?
                            round(Int, 1.0 / DS) : 0,
                        full_grid=true,
                        fixed_cell=true,
                    ))
                end
            end
        end
    end
    raw_df = DataFrame(raw)
    summary = combine(
        groupby(raw_df, [:dimension, :linear_grid, :grid_points, :method]),
        :wall_s => median => :median_wall_s,
        :wall_s => minimum => :minimum_wall_s,
        :wall_s => maximum => :maximum_wall_s,
        :allocated_bytes => median => :median_allocated_bytes,
        nrow => :repetitions,
    )
    ratios = NamedTuple[]
    for group in groupby(summary, [:dimension, :linear_grid, :grid_points])
        b = only(group.median_wall_s[group.method .== "bvk2"])
        s = only(group.median_wall_s[group.method .== "scft"])
        push!(ratios, (
            schema=SCHEMA,
            dimension=group.dimension[1],
            linear_grid=group.linear_grid[1],
            grid_points=group.grid_points[1],
            bvk2_median_wall_s=b,
            scft_median_wall_s=s,
            scft_over_bvk2=s / b,
        ))
    end
    return raw_df, summary, DataFrame(ratios)
end

function make_continuation_scft(architecture, chiN, period_rg, nx, workdir)
    template = exact_template(period_rg, (nx,))
    return NoncyclicChainSCFT(
        polymer_system(architecture, chiN),
        template,
        DS;
        mde=Polyorder.OSF,
        updater=Polyorder.Anderson(
            Polyorder.SD(1.0); warmup=150, αw=0.2),
        profile=:balanced,
        init=:randn,
        rng=Random.Xoshiro(RNG_SEED),
    ), scft_config(workdir)
end

function scft_density(scft, specie::Symbol)
    density_function = getproperty(Polyorder, Symbol(Char(0x03d5)))
    return vec(getfield(density_function(scft, 1, specie), :data))
end

function scft_accepted(scft, status, system, config)
    phi_a = scft_density(scft, :A)
    phi_b = scft_density(scft, :B)
    return occursin("Successful", string(status)) &&
        Float64(Polyorder.F(scft) - Polyorder.F_DIS(system)) < -1.0e-6 &&
        maximum(phi_a) - minimum(phi_a) > 1.0e-3 &&
        sqrt(mean(abs2, phi_a .+ phi_b .- 1.0)) <= 1.0e-6 &&
        Polyorder.residual(scft, config) <= 1.1e-6
end

function bvk2_anchor(architecture, chiN, period_rg, nx)
    common = (
        f=0.5,
        chiN=Float64(chiN),
        L=Float64(period_rg) / sqrt(6.0),
        nx=nx,
        mode_count=16,
        initial_amplitude=0.35,
        c2=C2,
        max_iterations=4000,
        force_tolerance=1.0e-4,
        force_maxabs_tolerance=2.0e-4,
    )
    return architecture == "AB" ?
        minimize_diblock_bvk2_lamella(; common...) :
        minimize_triblock_aba_bvk2_lamella(; common...)
end

function bvk2_continue(architecture, chiN, period_rg, nx, profile)
    common = (
        f=0.5,
        chiN=Float64(chiN),
        L=Float64(period_rg) / sqrt(6.0),
        nx=nx,
        mode_count=16,
        initial_profile=profile,
        c2=C2,
        max_iterations=BVK2_WARM_MAX_ITERATIONS,
        force_tolerance=1.0e-4,
        force_maxabs_tolerance=2.0e-4,
    )
    return architecture == "AB" ?
        minimize_diblock_bvk2_lamella(; common...) :
        minimize_triblock_aba_bvk2_lamella(; common...)
end

function bvk2_continue_reference(
        architecture, chiN, period_rg, nx, profile)
    common = (
        f=0.5,
        chiN=Float64(chiN),
        L=Float64(period_rg) / sqrt(6.0),
        nx=nx,
        mode_count=16,
        initial_profile=profile,
        c2=C2,
        max_iterations=4000,
        force_tolerance=1.0e-4,
        force_maxabs_tolerance=2.0e-4,
    )
    return architecture == "AB" ?
        minimize_diblock_bvk2_lamella(; common...) :
        minimize_triblock_aba_bvk2_lamella(; common...)
end

function run_bvk2_warm_accuracy_audit(nx)
    rows = NamedTuple[]
    for case in CONTINUATION_CASES
        anchor = bvk2_anchor(
            case.architecture, case.chi_anchor, case.period_rg, nx)
        target_chiN = first(case.chi_values)
        capped = bvk2_continue(case.architecture, target_chiN,
            case.period_rg, nx, anchor.phi_a)
        reference = bvk2_continue_reference(case.architecture, target_chiN,
            case.period_rg, nx, anchor.phi_a)
        capped.converged && reference.converged ||
            error("BVK2 warm-cap accuracy audit did not converge")
        energy_scale = capped.L
        energy_difference =
            (capped.energy - reference.energy) / energy_scale
        profile_rms =
            sqrt(mean(abs2, capped.phi_a .- reference.phi_a))
        push!(rows, (
            schema=SCHEMA,
            architecture=case.architecture,
            chiN=target_chiN,
            period_rg=case.period_rg,
            capped_max_iterations=BVK2_WARM_MAX_ITERATIONS,
            capped_iterations=capped.iterations,
            reference_max_iterations=4000,
            reference_iterations=reference.iterations,
            capped_energy_density=capped.energy / energy_scale,
            reference_energy_density=reference.energy / energy_scale,
            absolute_energy_density_difference=abs(energy_difference),
            profile_rms_difference=profile_rms,
            capped_residual_l2=capped.projected_force_norm,
            capped_residual_linf=capped.projected_force_maxabs,
            passed=abs(energy_difference) <= 1.0e-10 &&
                profile_rms <= 1.0e-5 &&
                capped.projected_force_norm <= 1.0e-4 &&
                capped.projected_force_maxabs <= 2.0e-4,
        ))
    end
    audit = DataFrame(rows)
    all(audit.passed) ||
        error("BVK2 warm iteration cap failed the accuracy audit")
    return audit
end

function run_continuation_path(case, trajectory, nx, workdir)
    rows = NamedTuple[]
    anchor_bvk2 =
        bvk2_anchor(case.architecture, case.chi_anchor, case.period_rg, nx)
    anchor_bvk2.converged ||
        error("BVK2 anchor failed for $(case.architecture)")
    scft, config = make_continuation_scft(
        case.architecture, case.chi_anchor, case.period_rg, nx, workdir)
    anchor_status = with_logger(NullLogger()) do
        Polyorder.solve!(scft, config)
    end
    scft_accepted(scft, anchor_status, scft.system, config) ||
        error("SCFT anchor failed for $(case.architecture)")
    bvk2_profile = copy(anchor_bvk2.phi_a)
    for chiN in case.chi_values
        bvk2_measurement = @timed bvk2_continue(
            case.architecture, chiN, case.period_rg, nx, bvk2_profile)
        bvk2 = bvk2_measurement.value
        bvk2.converged ||
            error("BVK2 continuation failed at $(case.architecture) chiN=$(chiN)")
        bvk2_profile = copy(bvk2.phi_a)
        push!(rows, (
            schema=SCHEMA,
            architecture=case.architecture,
            trajectory=trajectory,
            chiN=chiN,
            delta_chiN=chiN - case.chi_anchor,
            period_rg=case.period_rg,
            grid_points=nx,
            method="bvk2",
            setup_wall_s=0.0,
            solve_wall_s=bvk2_measurement.time,
            total_wall_s=bvk2_measurement.time,
            allocated_bytes=bvk2_measurement.bytes,
            iterations=bvk2.iterations,
            function_evaluations=missing,
            residual_l2=bvk2.projected_force_norm,
            residual_linf=bvk2.projected_force_maxabs,
            energy_density=bvk2.energy / bvk2.L,
            accepted=bvk2.converged &&
                bvk2.energy < bvk2.homogeneous_energy - 1.0e-8,
            fixed_cell=true,
            warm_seed=true,
        ))

        next_system = polymer_system(case.architecture, chiN)
        reset_measurement = @timed Polyorder.reset(scft, next_system)
        scft = reset_measurement.value
        solve_measurement = @timed with_logger(NullLogger()) do
            Polyorder.solve!(scft, config)
        end
        status = solve_measurement.value
        accepted = scft_accepted(scft, status, next_system, config)
        accepted ||
            error("SCFT continuation failed at $(case.architecture) chiN=$(chiN)")
        updater = scft.updater
        evaluations = isempty(updater.evals) ? updater.n : updater.evals[end]
        push!(rows, (
            schema=SCHEMA,
            architecture=case.architecture,
            trajectory=trajectory,
            chiN=chiN,
            delta_chiN=chiN - case.chi_anchor,
            period_rg=case.period_rg,
            grid_points=length(scft_density(scft, :A)),
            method="scft",
            setup_wall_s=reset_measurement.time,
            solve_wall_s=solve_measurement.time,
            total_wall_s=reset_measurement.time + solve_measurement.time,
            allocated_bytes=reset_measurement.bytes + solve_measurement.bytes,
            iterations=updater.n,
            function_evaluations=evaluations,
            residual_l2=Float64(Polyorder.residual(scft, config)),
            residual_linf=NaN,
            energy_density=Float64(Polyorder.F(scft)),
            accepted=accepted,
            fixed_cell=true,
            warm_seed=true,
        ))
    end
    return rows
end

function run_continuation_benchmark(outdir; trajectories=3, nx=128)
    raw = NamedTuple[]
    workdir = joinpath(outdir, "scratch")
    mkpath(workdir)
    # One excluded path compiles reset/solve and the architecture-specific BVK2
    # parameter builder before any measurement is retained.
    for case in CONTINUATION_CASES
        run_continuation_path(case, 0, nx, workdir)
    end
    for trajectory in 1:trajectories
        for case in CONTINUATION_CASES
            append!(raw,
                run_continuation_path(case, trajectory, nx, workdir))
        end
    end
    raw_df = DataFrame(raw)
    summary = combine(
        groupby(raw_df, [:architecture, :chiN, :period_rg, :method]),
        :setup_wall_s => median => :median_setup_wall_s,
        :solve_wall_s => median => :median_solve_wall_s,
        :total_wall_s => median => :median_total_wall_s,
        :allocated_bytes => median => :median_allocated_bytes,
        :iterations => median => :median_iterations,
        :residual_l2 => maximum => :maximum_residual_l2,
        :accepted => all => :accepted,
        nrow => :trajectories,
    )
    ratios = NamedTuple[]
    for group in groupby(summary, [:architecture, :chiN, :period_rg])
        b = only(group.median_total_wall_s[group.method .== "bvk2"])
        s = only(group.median_total_wall_s[group.method .== "scft"])
        push!(ratios, (
            schema=SCHEMA,
            architecture=group.architecture[1],
            chiN=group.chiN[1],
            period_rg=group.period_rg[1],
            bvk2_median_total_wall_s=b,
            scft_median_total_wall_s=s,
            scft_over_bvk2=s / b,
        ))
    end
    return raw_df, summary, DataFrame(ratios)
end

function log_map(value, lo, hi, y0, y1)
    return y1 - (log10(value) - log10(lo)) /
        (log10(hi) - log10(lo)) * (y1 - y0)
end

function write_kernel_svg(path, ratios)
    width, height = 1040, 500
    panels = (
        (dimension=1, x0=85.0, x1=500.0, title="1D equal full grid"),
        (dimension=3, x0=610.0, x1=1025.0, title="3D equal full grid"),
    )
    y0, y1 = 70.0, 405.0
    ymin = max(0.5, 10.0^floor(log10(minimum(ratios.scft_over_bvk2))))
    ymax = 10.0^ceil(log10(maximum(ratios.scft_over_bvk2)))
    open(path, "w") do io
        println(io, "<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
        println(io, "<rect width='100%' height='100%' fill='white'/>")
        println(io, "<style>.title{font:700 18px Arial}.label{font:14px Arial}.small{font:12px Arial}</style>")
        println(io, "<text x='520' y='28' text-anchor='middle' class='title'>Fixed-cell kernel cost: one nonlinear evaluation</text>")
        for panel in panels
            rows = ratios[ratios.dimension .== panel.dimension, :]
            xmin, xmax = minimum(rows.grid_points), maximum(rows.grid_points)
            println(io, "<rect x='$(panel.x0)' y='$y0' width='$(panel.x1-panel.x0)' height='$(y1-y0)' fill='none' stroke='#222'/>")
            println(io, "<text x='$((panel.x0+panel.x1)/2)' y='54' text-anchor='middle' class='label'>$(panel.title)</text>")
            for ratio in (1.0, 10.0, 100.0, 1000.0)
                ymin <= ratio <= ymax || continue
                y = log_map(ratio, ymin, ymax, y0, y1)
                println(io, "<line x1='$(panel.x0)' x2='$(panel.x1)' y1='$y' y2='$y' stroke='#ddd'/>")
                println(io, @sprintf("<text x='%.1f' y='%.1f' text-anchor='end' class='small'>%.0fx</text>", panel.x0-8, y+4, ratio))
            end
            points = String[]
            for row in eachrow(sort(rows, :grid_points))
                x = panel.x0 + (log10(row.grid_points)-log10(xmin)) /
                    (log10(xmax)-log10(xmin))*(panel.x1-panel.x0)
                y = log_map(row.scft_over_bvk2, ymin, ymax, y0, y1)
                push!(points, @sprintf("%.3f,%.3f", x, y))
                println(io, "<circle cx='$x' cy='$y' r='4.5' fill='#d55e00'/>")
                println(io, "<text x='$x' y='$(y-9)' text-anchor='middle' class='small'>$(round(row.scft_over_bvk2; sigdigits=3))x</text>")
            end
            println(io, "<polyline points='$(join(points, " "))' fill='none' stroke='#d55e00' stroke-width='2'/>")
            println(io, "<text x='$((panel.x0+panel.x1)/2)' y='438' text-anchor='middle' class='label'>grid points M (log scale)</text>")
            println(io, "<text x='$(panel.x0)' y='421' text-anchor='middle' class='small'>$xmin</text>")
            println(io, "<text x='$(panel.x1)' y='421' text-anchor='middle' class='small'>$xmax</text>")
        end
        println(io, "<text x='20' y='$((y0+y1)/2)' text-anchor='middle' class='label' transform='rotate(-90 20 $((y0+y1)/2))'>SCFT q! time / BVK2 energy+force time</text>")
        println(io, "<text x='520' y='478' text-anchor='middle' class='small'>OSF ds=0.005 (about 200 contour steps); full-grid FFT paths; one CPU thread; compilation excluded.</text>")
        println(io, "</svg>")
    end
end

function write_continuation_svg(path, summary)
    width, height = 1040, 500
    panels = (
        (architecture="AB", x0=85.0, x1=500.0, title="AB, fixed D/Rg=4.2"),
        (architecture="ABA", x0=610.0, x1=1025.0, title="ABA, fixed D/Rg=2.9"),
    )
    y0, y1 = 70.0, 405.0
    times = summary.median_total_wall_s
    ymin = 10.0^floor(log10(minimum(times)))
    ymax = 10.0^ceil(log10(maximum(times)))
    styles = Dict("bvk2" => ("#0072b2", "BVK2"),
        "scft" => ("#222222", "Polyorder SCFT"))
    open(path, "w") do io
        println(io, "<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
        println(io, "<rect width='100%' height='100%' fill='white'/>")
        println(io, "<style>.title{font:700 18px Arial}.label{font:14px Arial}.small{font:12px Arial}</style>")
        println(io, "<text x='520' y='28' text-anchor='middle' class='title'>Warm continuation at fixed cell</text>")
        for panel in panels
            rows = summary[summary.architecture .== panel.architecture, :]
            xmin, xmax = minimum(rows.chiN), maximum(rows.chiN)
            println(io, "<rect x='$(panel.x0)' y='$y0' width='$(panel.x1-panel.x0)' height='$(y1-y0)' fill='none' stroke='#222'/>")
            println(io, "<text x='$((panel.x0+panel.x1)/2)' y='54' text-anchor='middle' class='label'>$(panel.title)</text>")
            for value in 10.0 .^ (floor(log10(ymin)):ceil(log10(ymax)))
                ymin <= value <= ymax || continue
                y = log_map(value, ymin, ymax, y0, y1)
                println(io, "<line x1='$(panel.x0)' x2='$(panel.x1)' y1='$y' y2='$y' stroke='#ddd'/>")
                println(io, @sprintf("<text x='%.1f' y='%.1f' text-anchor='end' class='small'>%.3g</text>", panel.x0-8, y+4, value))
            end
            for method in ("bvk2", "scft")
                method_rows = sort(rows[rows.method .== method, :], :chiN)
                color, _ = styles[method]
                points = String[]
                for row in eachrow(method_rows)
                    x = panel.x0 + (row.chiN-xmin)/(xmax-xmin)*(panel.x1-panel.x0)
                    y = log_map(row.median_total_wall_s, ymin, ymax, y0, y1)
                    push!(points, @sprintf("%.3f,%.3f", x, y))
                    println(io, "<circle cx='$x' cy='$y' r='3.5' fill='$color'/>")
                end
                println(io, "<polyline points='$(join(points, " "))' fill='none' stroke='$color' stroke-width='2'/>")
            end
            println(io, "<text x='$((panel.x0+panel.x1)/2)' y='438' text-anchor='middle' class='label'>chiN</text>")
        end
        println(io, "<text x='20' y='$((y0+y1)/2)' text-anchor='middle' class='label' transform='rotate(-90 20 $((y0+y1)/2))'>median marginal point time (s, log scale)</text>")
        for (index, method) in enumerate(("bvk2", "scft"))
            color, label = styles[method]
            x = 330 + 210 * (index - 1)
            println(io, "<line x1='$x' x2='$(x+28)' y1='470' y2='470' stroke='$color' stroke-width='3'/>")
            println(io, "<text x='$(x+36)' y='475' class='small'>$label</text>")
        end
        println(io, "</svg>")
    end
end

function write_environment(path)
    payload = Dict(
        "schema" => SCHEMA,
        "date" => string(Dates.today()),
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "fftw_threads" => 1,
        "polyorder_version" => string(pkgversion(Polyorder)),
        "dfmmontecarlo_version" => string(pkgversion(DFMMonteCarlo)),
        "contour_step" => DS,
        "rng_seed" => RNG_SEED,
        "cpu_model" => try
            strip(read(`bash -lc "lscpu | sed -n 's/^Model name:[[:space:]]*//p'"`,
                String))
        catch
            "unknown"
        end,
    )
    open(path, "w") do io
        TOML.print(io, payload; sorted=true)
    end
end

function write_readme(
        path, kernel_ratios, continuation_ratios, accuracy_audit)
    k1 = kernel_ratios[kernel_ratios.dimension .== 1, :]
    k3 = kernel_ratios[kernel_ratios.dimension .== 3, :]
    ab = continuation_ratios[
        continuation_ratios.architecture .== "AB", :]
    aba = continuation_ratios[
        continuation_ratios.architecture .== "ABA", :]
    open(path, "w") do io
        println(io, "# BVK2 / SCFT fixed-cell performance benchmark")
        println(io)
        println(io, "This packet separates intrinsic nonlinear-evaluation cost from fixed-cell continuation cost. No cell optimization, stress root, period bracket, or local-period audit is timed.")
        println(io)
        println(io, "## Kernel benchmark")
        println(io)
        println(io, "- BVK2 timing unit: one energy plus analytic chemical-potential evaluation.")
        println(io, "- SCFT timing unit: one `Polyorder.q!` evaluation, including all OSF MDE propagations, density integration, and force construction.")
        println(io, "- Equal full grids, no DCT/CFFT/symmetry compression, one CPU thread, `ds=0.005` (about 200 contour steps).")
        println(io, "- 1D grids: `$(join(KERNEL_1D_GRIDS, ", "))`; 3D grids: `$(join(KERNEL_3D_GRIDS, ", "))`.")
        println(io)
        println(io, @sprintf(
            "SCFT/BVK2 kernel-time ratios span `%.2f--%.2fx` in 1D and `%.2f--%.2fx` in 3D.",
            minimum(k1.scft_over_bvk2), maximum(k1.scft_over_bvk2),
            minimum(k3.scft_over_bvk2), maximum(k3.scft_over_bvk2)))
        println(io)
        println(io, "## Fixed-cell warm continuation")
        println(io)
        println(io, "- AB: anchor `chiN=20`, then `21:30`, fixed `D/Rg=4.2`.")
        println(io, "- ABA: anchor `chiN=30`, then `31:40`, fixed `D/Rg=2.9`.")
        println(io, "- Anchors and one compilation trajectory are excluded. Every measured point is seeded from the immediately previous accepted point.")
        println(io, "- BVK2 uses the previous density profile. Polyorder uses `reset(scft, new_system)`, which carries the converged auxiliary fields and fixed lattice forward.")
        println(io, "- BVK2 warm L-BFGS is capped at `$BVK2_WARM_MAX_ITERATIONS` iterations and then uses its existing residual-Newton polish. This is an accuracy-qualified stop, not a weakened acceptance rule.")
        for row in eachrow(accuracy_audit)
            println(io, @sprintf(
                "- %s cap audit at chiN=%.0f: `|Delta F|=%.3e`, profile RMS `%.3e`, capped residuals `(%.3e, %.3e)` versus the 4000-iteration reference.",
                row.architecture, row.chiN,
                row.absolute_energy_density_difference,
                row.profile_rms_difference,
                row.capped_residual_l2, row.capped_residual_linf))
        end
        println(io, "- Reported time is the complete marginal point cost. For SCFT it includes reset/setup plus fixed-cell `solve!`; for BVK2 it is the complete fixed-period minimizer call.")
        println(io)
        println(io, @sprintf(
            "Median SCFT/BVK2 marginal-time ratios are `%.2fx` for AB and `%.2fx` for ABA.",
            median(ab.scft_over_bvk2), median(aba.scft_over_bvk2)))
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, "The kernel ratio tests the formal removal of the `Ns` MDE factor. The continuation ratio tests whether that lower evaluation cost survives a production-like warm solve. Neither result includes cell optimization, so it must not be compared directly with the earlier stress-free end-to-end benchmark.")
        println(io)
        println(io, "## Reproduction")
        println(io)
        println(io, "```bash")
        println(io, "JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 12 \\")
        println(io, "  julia --project=. scripts/write_bvk2_scft_fixedcell_performance_benchmark.jl")
        println(io, "```")
    end
end

function main()
    outdir = abspath(arg("--outdir", DEFAULT_OUTDIR))
    trajectories = parse(Int, arg("--trajectories", "3"))
    nx = parse(Int, arg("--continuation-nx", "128"))
    mkpath(outdir)
    BLAS.set_num_threads(1)
    FFTW.set_num_threads(1)

    kernel_raw, kernel_summary, kernel_ratios = run_kernel_benchmark()
    accuracy_audit = run_bvk2_warm_accuracy_audit(nx)
    continuation_raw, continuation_summary, continuation_ratios =
        run_continuation_benchmark(outdir;
            trajectories=trajectories, nx=nx)
    all(continuation_raw.accepted) ||
        error("warm continuation contains an unaccepted point")

    CSV.write(joinpath(outdir, "kernel_raw.csv"), kernel_raw)
    CSV.write(joinpath(outdir, "kernel_summary.csv"), kernel_summary)
    CSV.write(joinpath(outdir, "kernel_ratios.csv"), kernel_ratios)
    CSV.write(joinpath(outdir, "continuation_accuracy_audit.csv"),
        accuracy_audit)
    CSV.write(joinpath(outdir, "continuation_raw.csv"), continuation_raw)
    CSV.write(joinpath(outdir, "continuation_summary.csv"),
        continuation_summary)
    CSV.write(joinpath(outdir, "continuation_ratios.csv"),
        continuation_ratios)
    write_kernel_svg(joinpath(outdir, "kernel_benchmark.svg"),
        kernel_ratios)
    write_continuation_svg(
        joinpath(outdir, "warm_continuation_benchmark.svg"),
        continuation_summary)
    write_environment(joinpath(outdir, "environment.toml"))
    write_readme(joinpath(outdir, "README.md"),
        kernel_ratios, continuation_ratios, accuracy_audit)
    rm(joinpath(outdir, "scratch"); recursive=true, force=true)
    println("wrote fixed-cell performance packet to ", outdir)
end

main()
