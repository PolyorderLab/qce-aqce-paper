#!/usr/bin/env julia

"""Compute the AB-diblock RPA stability limit used by BVK2 plots.

Run from the MonteCarlo project root with the isolated Polyorder environment:

    julia --startup-file=no --project=reproducibility/polyorder \
        scripts/write_bvk2_rpa_stability_limit.jl

The CSV stores the independently computed `f_A <= 0.5` branch. Plotters reflect
it about `f_A = 0.5`, just as they do for the numerical phase boundaries.
"""

using Optim
using Polyorder
using Roots

const SCHEMA = "polyorder-rpa-ab-diblock-stability-v1"
const DEFAULT_OUTPUT = normpath(joinpath(
    @__DIR__, "..", "results", "bvk2_publication_phase_diagram",
    "rpa_stability_limit.csv",
))

function parse_args(args::Vector{String})
    output = DEFAULT_OUTPUT
    points = 181
    f_min = 0.05
    index = 1
    while index <= length(args)
        option = args[index]
        index == length(args) && error("missing value for $option")
        value = args[index + 1]
        if option == "--output"
            output = abspath(value)
        elseif option == "--points"
            points = parse(Int, value)
        elseif option == "--f-min"
            f_min = parse(Float64, value)
        else
            error("unknown option: $option")
        end
        index += 2
    end
    0.0 < f_min < 0.5 || error("--f-min must lie strictly between 0 and 0.5")
    points >= 2 || error("--points must be at least 2")
    return output, points, f_min
end

function compute_curve(points::Int, f_min::Float64)
    rows = NamedTuple[]
    for f_a in range(f_min, 0.5; length=points)
        system = AB_system(; fA=f_a)
        chi_n_star, d_star = Polyorder.RPA.compute_stability_limit(
            system; χN0=10.0,
        )
        k_star = 2pi / d_star
        all(isfinite, (chi_n_star, d_star, k_star)) ||
            error("non-finite RPA result at f_A=$f_a")
        chi_n_star > 0.0 && d_star > 0.0 && k_star > 0.0 ||
            error("non-positive RPA result at f_A=$f_a")
        push!(rows, (
            fA=f_a,
            chiN_spinodal=chi_n_star,
            D_star_Rg=d_star,
            k_star_Rg_inv=k_star,
        ))
    end

    symmetric = rows[end].chiN_spinodal
    isapprox(symmetric, 10.4949; atol=5e-4, rtol=0.0) ||
        error("unexpected symmetric-diblock spinodal: $symmetric")
    return rows
end

function write_curve(path::String, rows)
    mkpath(dirname(path))
    version = pkgversion(Polyorder)
    open(path, "w") do io
        println(
            io,
            "schema,fA,chiN_spinodal,D_star_Rg,k_star_Rg_inv,polyorder_version",
        )
        for row in rows
            println(
                io,
                join((
                    SCHEMA,
                    repr(row.fA),
                    repr(row.chiN_spinodal),
                    repr(row.D_star_Rg),
                    repr(row.k_star_Rg_inv),
                    version,
                ), ','),
            )
        end
    end
end

output, points, f_min = parse_args(ARGS)
curve = compute_curve(points, f_min)
write_curve(output, curve)
println(
    "wrote $(length(curve)) Polyorder.RPA points to $output ",
    "(symmetric chiN*=$(curve[end].chiN_spinodal))",
)
