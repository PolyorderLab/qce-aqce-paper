const FIXED_PROFILE_AUDIT_ROOT = normpath(joinpath(@__DIR__, ".."))

ENV["DFM_SKIP_BVK2_FIXED_FIGURE6_ROOT_MAIN"] = "1"
include(joinpath(@__DIR__, "run_bvk2_fixed_stiffness_figure6_root.jl"))

function fixed_profile_audit_case(case_id::AbstractString)
    matched = match(r"^f([0-9]+(?:\.[0-9]+)?)_chiN([0-9]+(?:\.[0-9]+)?)$",
        String(case_id))
    matched === nothing && throw(ArgumentError(
        "case must have the form f0.25_chiN26"))
    f = parse(Float64, matched.captures[1])
    chiN = parse(Float64, matched.captures[2])
    0.0 < f <= 0.5 || throw(ArgumentError("f must lie in (0, 0.5]"))
    chiN > 0.0 || throw(ArgumentError("chiN must be positive"))
    return (f=f, chiN=chiN)
end

function fixed_profile_audit_main()
    case_id = String(bfa_arg_value("--case"))
    case = fixed_profile_audit_case(case_id)
    BFA_CASES[case_id] = case
    main()
end

fixed_profile_audit_main()
