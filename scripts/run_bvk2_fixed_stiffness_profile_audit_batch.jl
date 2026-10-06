const FIXED_PROFILE_BATCH_ROOT = normpath(joinpath(@__DIR__, ".."))

ENV["DFM_SKIP_BVK2_FIXED_FIGURE6_ROOT_MAIN"] = "1"
include(joinpath(@__DIR__, "run_bvk2_fixed_stiffness_figure6_root.jl"))

const FIXED_PROFILE_BATCH_USAGE = """
Compute matched fixed-stiffness BVK roots for several accepted adaptive cases.

Usage:
  julia --project=. scripts/run_bvk2_fixed_stiffness_profile_audit_batch.jl \\
    --cases=f0.5_chiN36,f0.5_chiN37 \\
    --source-results-dir=results/bvk_strong_segregation_extension \\
    --outdir=results/bvk_strong_segregation_extension/fixed_cases
"""

function fixed_profile_batch_arg(name::AbstractString, default=nothing)
    prefix = String(name) * "="
    for argument in ARGS
        startswith(argument, prefix) &&
            return split(argument, "="; limit=2)[2]
    end
    default === nothing && throw(ArgumentError("missing $(name)=..."))
    return default
end

function fixed_profile_batch_case(case_id::AbstractString)
    matched = match(r"^f([0-9]+(?:\.[0-9]+)?)_chiN([0-9]+(?:\.[0-9]+)?)$",
        String(case_id))
    matched === nothing && throw(ArgumentError(
        "case must have the form f0.25_chiN36"))
    f = parse(Float64, matched.captures[1])
    chiN = parse(Float64, matched.captures[2])
    0.0 < f <= 0.5 || throw(ArgumentError("f must lie in (0, 0.5]"))
    chiN > 0.0 || throw(ArgumentError("chiN must be positive"))
    return (f=f, chiN=chiN)
end

function fixed_profile_batch_main()
    if any(argument -> argument in ("--help", "-h"), ARGS)
        print(FIXED_PROFILE_BATCH_USAGE)
        return nothing
    end
    case_ids = [strip(value) for value in split(
        fixed_profile_batch_arg("--cases"), ",") if !isempty(strip(value))]
    isempty(case_ids) && throw(ArgumentError("--cases must not be empty"))
    length(unique(case_ids)) == length(case_ids) || throw(ArgumentError(
        "--cases contains duplicate case identifiers"))
    source_results_dir = normpath(fixed_profile_batch_arg(
        "--source-results-dir"))
    outdir = normpath(fixed_profile_batch_arg("--outdir"))
    original_args = copy(ARGS)
    try
        for case_id in case_ids
            BFA_CASES[case_id] = fixed_profile_batch_case(case_id)
            empty!(ARGS)
            append!(ARGS, [
                "--case=$(case_id)",
                "--source-results-dir=$(source_results_dir)",
                "--outdir=$(joinpath(outdir, case_id))",
            ])
            main()
        end
    finally
        empty!(ARGS)
        append!(ARGS, original_args)
    end
    return nothing
end

fixed_profile_batch_main()
