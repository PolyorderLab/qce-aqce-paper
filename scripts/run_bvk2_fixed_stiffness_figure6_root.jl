const MONTE_CARLO_PROJECT = normpath(joinpath(@__DIR__, ".."))
if !(MONTE_CARLO_PROJECT in LOAD_PATH)
    pushfirst!(LOAD_PATH, MONTE_CARLO_PROJECT)
end

using CSV
using DataFrames
using DFMMonteCarlo
using FFTW
using Printf
using SHA

const BFA_SCHEMA = "bvk2-fixed-stiffness-figure6-root-v1"
const BFA_C2 = 0.16
const BFA_NX = 256
const BFA_OVERSAMPLE = 2
const BFA_SENSOR_FILTER_RATIO = 12.0
const BFA_CASES = Dict(
    "f0.5_chiN20" => (f=0.50, chiN=20.0),
    "f0.35_chiN30" => (f=0.35, chiN=30.0),
)

function bfa_file_fingerprint(path::AbstractString)
    return bytes2hex(sha256(read(path)))
end

function bfa_arg_value(name::AbstractString, default=nothing)
    prefix = String(name) * "="
    for arg in ARGS
        startswith(arg, prefix) && return split(arg, "=", limit=2)[2]
    end
    default === nothing && throw(ArgumentError("missing required $(name)=..."))
    return default
end

function bfa_source(case_id::AbstractString)
    source_root = normpath(bfa_arg_value("--source-results-dir", joinpath(
        MONTE_CARLO_PROJECT, "results", "bvk2_filtered_period_map_r12_nx256")))
    source_dir = joinpath(source_root, "cases", case_id)
    summary_path = joinpath(source_dir, "summary.csv")
    profile_path = joinpath(source_dir, "native_profiles.csv")
    isfile(summary_path) || throw(ArgumentError("missing $(summary_path)"))
    isfile(profile_path) || throw(ArgumentError("missing $(profile_path)"))
    summary = CSV.read(summary_path, DataFrame)
    rows = summary[String.(summary.model) .== "bvk2", :]
    nrow(rows) == 1 || throw(ArgumentError(
        "expected one production BVK2 row for $(case_id)"))
    row = rows[1, :]
    String(row.case_id) == case_id || throw(ArgumentError(
        "production BVK2 source has wrong case id"))
    String(row.model) == "bvk2" || throw(ArgumentError(
        "production BVK2 source has wrong model"))
    Float64(row.c2) == BFA_C2 || throw(ArgumentError(
        "production BVK2 source has wrong c2"))
    isfinite(Float64(row.period_rg)) || throw(ArgumentError(
        "production BVK2 source has nonfinite period"))
    String(row.status) == "accepted" || throw(ArgumentError(
        "production BVK2 source is not accepted for $(case_id)"))
    Bool(row.field_gate_pass) && Bool(row.cell_gate_pass) &&
        Bool(row.resolution_gate_pass) && Bool(row.composition_gate_pass) &&
        Bool(row.morphology_gate_pass) ||
        throw(ArgumentError("production BVK2 source gates failed for $(case_id)"))
    Int(row.grid_count) == BFA_NX || throw(ArgumentError(
        "production BVK2 source is not nx=$(BFA_NX)"))
    Int(row.oversample_factor) == BFA_OVERSAMPLE || throw(ArgumentError(
        "production BVK2 source has wrong oversample factor"))
    isapprox(Float64(row.sensor_filter_ratio), BFA_SENSOR_FILTER_RATIO;
        atol=0.0, rtol=0.0) || throw(ArgumentError(
        "production BVK2 source has wrong sensor filter"))

    profiles = CSV.read(profile_path, DataFrame)
    profiles = profiles[String.(profiles.model) .== "bvk2", :]
    sort!(profiles, :native_index)
    nrow(profiles) == BFA_NX || throw(ArgumentError(
        "expected $(BFA_NX) native source samples for $(case_id)"))
    profile = Float64.(profiles.phi_a)
    all(isfinite, profile) || throw(ArgumentError(
        "production BVK2 profile contains nonfinite values"))
    all(isapprox.(Float64.(profiles.period_rg), Float64(row.period_rg);
        atol=1.0e-12, rtol=0.0)) || throw(ArgumentError(
        "production BVK2 profile period mismatch"))
    bfa_periodic_domain_count(profile, Float64(row.f)) == 1 ||
        throw(ArgumentError("production BVK2 profile is not one-domain LAM"))
    summary_sha256 = bfa_file_fingerprint(summary_path)
    profile_sha256 = bfa_file_fingerprint(profile_path)
    length(summary_sha256) == 64 && length(profile_sha256) == 64 ||
        throw(ArgumentError("production BVK2 source provenance is incomplete"))
    return (
        row=row,
        profile=profile,
        summary_path=summary_path,
        profile_path=profile_path,
        summary_sha256=summary_sha256,
        profile_sha256=profile_sha256,
    )
end

function bfa_periodic_domain_count(profile::AbstractVector{<:Real},
        threshold::Real)
    mask = Float64.(profile) .> Float64(threshold)
    (all(mask) || all(.!mask)) && return 0
    return count(index -> mask[index] && !mask[mod1(index - 1,
        length(mask))], eachindex(mask))
end

function bfa_dominant_mode(profile::AbstractVector{<:Real}, f::Real)
    spectrum = abs.(rfft(Float64.(profile) .- Float64(f)))
    length(spectrum) >= 2 || return 0
    return argmax(@view spectrum[2:end])
end

function bfa_stress(result, case; adaptive::Bool=false)
    theta = asin.(sqrt.(clamp.(result.phi_a, 0.0, 1.0)))
    return diblock_bvk2_ud_theta_cell_scale_gradient_density_nd(theta;
        f=case.f, chiN=case.chiN, lengths=(result.L,), adaptive=adaptive,
        c2=BFA_C2, oversample=BFA_OVERSAMPLE,
        sensor_filter_ratio=BFA_SENSOR_FILTER_RATIO)
end

function bfa_relax_local(root, case, factor::Real; adaptive::Bool=false)
    period = root.period_b * Float64(factor)
    result = minimize_diblock_bvk2_ud_theta_lamella(; f=case.f,
        chiN=case.chiN, L=period, nx=BFA_NX, mode_count=3,
        initial_profile=root.result.phi_a, adaptive=adaptive, c2=BFA_C2,
        oversample=BFA_OVERSAMPLE,
        sensor_filter_ratio=BFA_SENSOR_FILTER_RATIO,
        max_iterations=4000, energy_tolerance=1.0e-10,
        force_tolerance=1.0e-5, force_maxabs_tolerance=2.0e-5,
        time_limit=1800.0)
    return (result=result, objective=result.energy / result.L,
        stress=bfa_stress(result, case; adaptive=adaptive))
end

function bfa_git_provenance()
    commit = try
        readchomp(`git -C $MONTE_CARLO_PROJECT rev-parse HEAD`)
    catch
        "unknown"
    end
    dirty = try
        !isempty(strip(read(`git -C $MONTE_CARLO_PROJECT status --porcelain`,
            String)))
    catch
        true
    end
    return (commit=commit, dirty=dirty)
end

function bfa_validate_staged_outputs(root_path::AbstractString,
        profile_path::AbstractString)
    root = CSV.read(root_path, DataFrame)
    nrow(root) == 1 || throw(ArgumentError(
        "fixed-stiffness publication requires exactly one root row"))
    row = root[1, :]
    required_root = (:schema, :case_id, :model, :adaptive, :c2, :nx,
        :oversample_factor, :sensor_filter_ratio, :source_summary,
        :source_profile, :source_summary_sha256, :source_profile_sha256,
        :field_gate_pass, :composition_gate_pass, :morphology_gate_pass,
        :cell_gate_pass, :local_minimum_check_pass,
        :stress_orientation_pass, :accepted, :status, :period_rg)
    all(name -> name in propertynames(root), required_root) ||
        throw(ArgumentError("fixed-stiffness root lacks publication columns"))
    String(row.schema) == BFA_SCHEMA || throw(ArgumentError(
        "fixed-stiffness root has wrong schema"))
    haskey(BFA_CASES, String(row.case_id)) || throw(ArgumentError(
        "fixed-stiffness root has unsupported case"))
    String(row.model) == "bvk2_fixed" && !Bool(row.adaptive) ||
        throw(ArgumentError("fixed-stiffness root has wrong model"))
    Float64(row.c2) == BFA_C2 && Int(row.nx) == BFA_NX &&
        Int(row.oversample_factor) == BFA_OVERSAMPLE &&
        Float64(row.sensor_filter_ratio) == BFA_SENSOR_FILTER_RATIO ||
        throw(ArgumentError("fixed-stiffness root has wrong production settings"))
    all(Bool(getproperty(row, name)) for name in (
        :field_gate_pass, :composition_gate_pass, :morphology_gate_pass,
        :cell_gate_pass, :local_minimum_check_pass,
        :stress_orientation_pass, :accepted)) || throw(ArgumentError(
        "fixed-stiffness root is not accepted by every gate"))
    String(row.status) == "accepted" || throw(ArgumentError(
        "fixed-stiffness root status is not accepted"))
    isfinite(Float64(row.period_rg)) || throw(ArgumentError(
        "fixed-stiffness root period is nonfinite"))
    source_paths = Dict{Symbol,String}()
    for name in (:source_summary, :source_profile)
        value = first(split(String(getproperty(row, name)), '#'))
        !isempty(value) || throw(ArgumentError(
            "fixed-stiffness source path is empty"))
        !isabspath(value) || throw(ArgumentError(
            "fixed-stiffness source paths must be repository-relative"))
        resolved = normpath(joinpath(MONTE_CARLO_PROJECT, value))
        isfile(resolved) || throw(ArgumentError(
            "fixed-stiffness source artifact is missing: $(value)"))
        source_paths[name] = resolved
    end
    for (fingerprint_name, path_name) in ((:source_summary_sha256,
            :source_summary), (:source_profile_sha256, :source_profile))
        value = String(getproperty(row, fingerprint_name))
        occursin(r"^[0-9a-f]{64}$", value) || throw(ArgumentError(
            "fixed-stiffness source fingerprint is invalid"))
        value == bfa_file_fingerprint(source_paths[path_name]) ||
            throw(ArgumentError(
                "fixed-stiffness source fingerprint does not match artifact"))
    end

    profile = CSV.read(profile_path, DataFrame)
    required_profile = (:schema, :case_id, :model, :adaptive, :c2, :nx,
        :period_rg, :index, :s, :phi_a)
    all(name -> name in propertynames(profile), required_profile) ||
        throw(ArgumentError("fixed-stiffness profile lacks publication columns"))
    nrow(profile) == BFA_NX || throw(ArgumentError(
        "fixed-stiffness profile has wrong sample count"))
    all(String.(profile.schema) .== BFA_SCHEMA) &&
        all(String.(profile.case_id) .== String(row.case_id)) &&
        all(String.(profile.model) .== "bvk2_fixed") &&
        all(.!Bool.(profile.adaptive)) &&
        all(Float64.(profile.c2) .== BFA_C2) &&
        all(Int.(profile.nx) .== BFA_NX) || throw(ArgumentError(
        "fixed-stiffness profile metadata mismatch"))
    values = Float64.(profile.phi_a)
    all(isfinite, values) && all(isfinite, Float64.(profile.s)) ||
        throw(ArgumentError("fixed-stiffness profile contains nonfinite values"))
    all(isapprox.(Float64.(profile.period_rg), Float64(row.period_rg);
        atol=1.0e-12, rtol=0.0)) || throw(ArgumentError(
        "fixed-stiffness profile period mismatch"))
    bfa_periodic_domain_count(values, BFA_CASES[String(row.case_id)].f) == 1 ||
        throw(ArgumentError("fixed-stiffness profile is not one-domain LAM"))
    return true
end

function bfa_publish_accepted(outdir::AbstractString, root_row,
        profile_rows; promotion_hook=(_stage -> nothing))
    parent = dirname(normpath(outdir))
    mkpath(parent)
    staging = mktempdir(parent; prefix=".bvk2-fixed-root-stage-")
    target_root = joinpath(outdir, "root.csv")
    target_profile = joinpath(outdir, "profile.csv")
    incoming_root = joinpath(outdir, ".root.csv.incoming")
    incoming_profile = joinpath(outdir, ".profile.csv.incoming")
    backup_root = joinpath(staging, "root.csv.backup")
    backup_profile = joinpath(staging, "profile.csv.backup")
    root_existed = isfile(target_root)
    profile_existed = isfile(target_profile)
    root_promoted = false
    profile_promoted = false
    try
        staged_root = joinpath(staging, "root.csv")
        staged_profile = joinpath(staging, "profile.csv")
        CSV.write(staged_root, DataFrame([root_row]))
        CSV.write(staged_profile, DataFrame(profile_rows))
        bfa_validate_staged_outputs(staged_root, staged_profile)
        mkpath(outdir)
        root_existed && cp(target_root, backup_root; force=true)
        profile_existed && cp(target_profile, backup_profile; force=true)
        cp(staged_root, incoming_root; force=true)
        cp(staged_profile, incoming_profile; force=true)
        promotion_hook(:before_promote)
        mv(incoming_profile, target_profile; force=true)
        profile_promoted = true
        promotion_hook(:after_profile)
        mv(incoming_root, target_root; force=true)
        root_promoted = true
        promotion_hook(:after_root)
    catch
        if profile_promoted
            profile_existed ? cp(backup_profile, target_profile; force=true) :
                rm(target_profile; force=true)
        end
        if root_promoted
            root_existed ? cp(backup_root, target_root; force=true) :
                rm(target_root; force=true)
        end
        rethrow()
    finally
        rm(incoming_root; force=true)
        rm(incoming_profile; force=true)
        rm(staging; recursive=true, force=true)
    end
    return (root=target_root, profile=target_profile)
end

function main()
    case_id = String(bfa_arg_value("--case"))
    haskey(BFA_CASES, case_id) || throw(ArgumentError(
        "unsupported case $(case_id); expected one of $(collect(keys(BFA_CASES)))"))
    case = BFA_CASES[case_id]
    outdir = normpath(bfa_arg_value("--outdir", joinpath(
        MONTE_CARLO_PROJECT, "results",
        "bvk2_fixed_stiffness_figure6_nx256", case_id)))
    mkpath(outdir)
    source = bfa_source(case_id)
    source_factor = Float64(source.row.period_factor)
    println(@sprintf(
        "[%s] fixed-stiffness root; source adaptive factor=%.12g period/Rg=%.12g",
        case_id, source_factor, Float64(source.row.period_rg)))
    flush(stdout)

    root = minimize_diblock_bvk2_ud_theta_lamella_stress_free(;
        f=case.f, chiN=case.chiN, nx=BFA_NX, mode_count=3,
        initial_profile=source.profile, adaptive=false, c2=BFA_C2,
        oversample=BFA_OVERSAMPLE,
        sensor_filter_ratio=BFA_SENSOR_FILTER_RATIO,
        max_iterations=4000, max_period_iterations=30,
        bootstrap_period_factors=(source_factor,), bootstrap_window=0.12,
        lower_factor=0.75, upper_factor=2.0,
        energy_tolerance=1.0e-10, force_tolerance=1.0e-5,
        force_maxabs_tolerance=2.0e-5,
        cell_stress_tolerance=1.0e-6, log_period_tolerance=1.0e-5,
        local_check_fraction=0.01, local_check_tolerance=1.0e-8,
        time_limit=1800.0, progress_label=case_id * " fixed Kpsi")

    center_stress = bfa_stress(root.result, case; adaptive=false)
    left = bfa_relax_local(root, case, 0.99; adaptive=false)
    right = bfa_relax_local(root, case, 1.01; adaptive=false)
    composition_pass = abs(root.result.mean_phi - case.f) <= 2.0e-8
    force_pass = root.result.projected_force_norm <= 1.0e-5 &&
        root.result.projected_force_maxabs <= 2.0e-5
    morphology_pass = root.result.maximum_phi - root.result.minimum_phi >
        1.0e-4 && root.result.energy < root.result.homogeneous_energy -
        1.0e-10 && bfa_dominant_mode(root.result.phi_a, case.f) == 1 &&
        bfa_periodic_domain_count(root.result.phi_a, case.f) == 1
    stress_pass = abs(center_stress) <= 1.0e-6
    orientation_pass = left.result.converged && right.result.converged &&
        left.stress < 0.0 < right.stress
    energy_pass = left.objective >= root.objective - 1.0e-8 &&
        right.objective >= root.objective - 1.0e-8
    accepted = root.converged && root.local_minimum_check_pass &&
        !root.boundary_limited && root.result.converged && composition_pass &&
        force_pass && morphology_pass && stress_pass && orientation_pass &&
        energy_pass
    provenance = bfa_git_provenance()
    root_row = (
        schema=BFA_SCHEMA, claim_kind="stress_free_lamellar_observable",
        case_id=case_id, f=case.f, chiN=case.chiN,
        model="bvk2_fixed", model_label="BVK2, fixed stiffness",
        adaptive=false, c2=BFA_C2, nx=BFA_NX,
        oversample_factor=BFA_OVERSAMPLE,
        sensor_filter_ratio=BFA_SENSOR_FILTER_RATIO,
        discretization_schema=
            "ud-theta-spectral-fixed-k-physical-sensor-filter-oversampled-v2",
        source_model="bvk2", source_period_rg=Float64(source.row.period_rg),
        source_period_factor=source_factor,
        source_summary=relpath(source.summary_path, MONTE_CARLO_PROJECT),
        source_profile=relpath(source.profile_path, MONTE_CARLO_PROJECT) *
            "#model=bvk2", period_rg=root.period_rg,
        source_summary_sha256=source.summary_sha256,
        source_profile_sha256=source.profile_sha256,
        period_b=root.period_b, period_factor=root.period_factor,
        objective=root.objective, optimizer=root.optimizer,
        root_iterations=root.iterations,
        field_iterations=root.result.iterations,
        mean_phi=root.result.mean_phi,
        minimum_phi=root.result.minimum_phi,
        maximum_phi=root.result.maximum_phi,
        projected_force_norm=root.result.projected_force_norm,
        projected_force_maxabs=root.result.projected_force_maxabs,
        root_cell_stress=center_stress,
        local_left_factor=0.99 * root.period_factor,
        local_left_objective=left.objective,
        local_left_stress=left.stress,
        local_right_factor=1.01 * root.period_factor,
        local_right_objective=right.objective,
        local_right_stress=right.stress,
        field_gate_pass=root.result.converged && force_pass,
        composition_gate_pass=composition_pass,
        morphology_gate_pass=morphology_pass,
        cell_gate_pass=root.converged && stress_pass && orientation_pass &&
            energy_pass && !root.boundary_limited,
        local_minimum_check_pass=root.local_minimum_check_pass && energy_pass,
        stress_orientation_pass=orientation_pass,
        accepted=accepted, status=accepted ? "accepted" : "rejected",
        status_reason=accepted ? "all_contract_gates_pass" :
            "one_or_more_contract_gates_failed",
        git_commit=provenance.commit, git_worktree_dirty=provenance.dirty,
    )
    profile_rows = [(
        schema=BFA_SCHEMA, case_id=case_id, f=case.f, chiN=case.chiN,
        model="bvk2_fixed", adaptive=false, c2=BFA_C2, nx=BFA_NX,
        oversample_factor=BFA_OVERSAMPLE,
        sensor_filter_ratio=BFA_SENSOR_FILTER_RATIO,
        period_rg=root.period_rg, index=index,
        s=(index - 0.5) / BFA_NX, phi_a=root.result.phi_a[index])
        for index in eachindex(root.result.phi_a)]
    bfa_publish_accepted(outdir, root_row, profile_rows)
    println(@sprintf(
        "status=%s period/Rg=%.12g stress=%.4g force=%.4g maxforce=%.4g left_stress=%.4g right_stress=%.4g",
        root_row.status, root.period_rg, center_stress,
        root.result.projected_force_norm,
        root.result.projected_force_maxabs, left.stress, right.stress))
end

if get(ENV, "DFM_SKIP_BVK2_FIXED_FIGURE6_ROOT_MAIN", "0") != "1"
    main()
end
