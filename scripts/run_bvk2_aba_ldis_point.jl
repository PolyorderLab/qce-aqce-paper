#!/usr/bin/env julia

"""Compute and validate one ABA LAM--DIS ordered-branch point.

The runner reuses the publication workflow without rebuilding the complete ABA
benchmark. Results are written to an isolated directory and require explicit
promotion into the canonical crossing table.
"""

const ROOT = normpath(joinpath(@__DIR__, ".."))
ENV["DFM_SKIP_BVK2_ABA_COMPREHENSIVE_MAIN"] = "1"
include(joinpath(@__DIR__, "write_bvk2_aba_comprehensive_validation.jl"))

using CSV
using DataFrames

function argument(name::AbstractString, default)
    prefix = String(name) * "="
    for value in ARGS
        startswith(value, prefix) && return split(value, "="; limit=2)[2]
    end
    return default
end

function main()
    chiN = parse(Float64, argument("--chiN", "21.0"))
    outdir = abspath(argument(
        "--outdir", joinpath(ROOT, "results", "bvk2_aba_ldis_chi21")))
    mkpath(outdir)
    crossing, _ = aba_cv_write_boundary(
        outdir; recompute=true, chis=(chiN,))

    nrow(crossing) == 2 || error("expected exactly two model rows")
    Set(String.(crossing.model)) == Set(("scft", "bvk2")) ||
        error("expected matched SCFT and BVK2 rows")
    all(Float64.(crossing.chiN) .== chiN) ||
        error("crossing table contains the wrong chiN")
    for row in eachrow(crossing)
        all(isfinite, (
            Float64(row.delta_free_energy), Float64(row.amplitude),
            Float64(row.period_rg))) || error("nonfinite plotted observable")
        Float64(row.delta_free_energy) < 0.0 ||
            error("ordered branch is not below DIS")
        Float64(row.amplitude) > 0.0 || error("ordered branch has zero amplitude")
        if String(row.model) == "scft"
            String(row.field_status) == "Polyorder.Successful()" ||
                error("SCFT field solve did not converge")
            String(row.cell_status) == "Polyorder.Successful()" ||
                error("SCFT cell solve did not converge")
            isfinite(Float64(row.cell_stress)) || error("SCFT stress is nonfinite")
            abs(Float64(row.cell_stress)) <= 1.0e-5 ||
                error("SCFT stress exceeds 1e-5")
        else
            lowercase(String(row.field_status)) == "true" ||
                error("BVK2 field solve did not converge")
            lowercase(String(row.cell_status)) == "true" ||
                error("BVK2 cell minimum gate failed")
        end
    end

    decision = DataFrame([(
        schema="bvk2-aba-ldis-point-decision-v1",
        claim_kind="ordered_branch_point",
        chiN=chiN,
        models="scft;bvk2",
        scientific_status="accepted",
        decision="accept",
        reason="matched finite ordered branches; field, cell, stress, and finiteness gates pass",
        source="ldis_crossing.csv",
    )])
    CSV.write(joinpath(outdir, "decision.csv"), decision; newline='\n')
    println("accepted ABA L/DIS point at chiN=", chiN, " in ", outdir)
end

main()
