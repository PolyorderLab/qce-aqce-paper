using Test
using FFTW
using LinearAlgebra
using Statistics
using TOML
using DFMMonteCarlo
import Arianna

include("burp_ti_newton_gmres_tests.jl")
include("burp_ti_lamellar_period_workflow_tests.jl")
include("liu2019_period_map_history_tests.jl")
include("ud_low_fa_regression_tests.jl")
include("bvk1_newton_ngmres_tests.jl")
include("bvk1_gc_gl_sd_newton_adapter_tests.jl")
include("bvk1_gc_gl_sd_newton_sample_tests.jl")
include("bvk1_richardson_correction_tests.jl")
include("bvk1_spectral_tests.jl")
include("bvk1_dealiased_richardson_contract_tests.jl")
include("bvk1_spectral_theta_fixedcell_tests.jl")
include("bvk1_spectral_theta_stressfree_tests.jl")
include("bvk1_lam1024_direct_reference_tests.jl")
include("bvk1_lam1024_endpoint_campaign_tests.jl")
include("bvk1_nonsmooth_kkt_audit_tests.jl")
include("bvk2_tests.jl")
include("generic_n_scaling_tests.jl")
include("burp_bvk2_hybrid_tests.jl")
include("burp_bvk2_hybrid_ablation_tests.jl")
include("bvk2_c2_calibration_tests.jl")
include("bvk2_morphology_tests.jl")
include("bvk2_ud_theta_gl_chi45_tests.jl")
include("bvk2_lam_resolution_mechanism_tests.jl")
include("bvk2_nonlamellar_dispatch_tests.jl")
include("bvk2_morphology_pilot_runner_tests.jl")
include("bvk2_boundary_unit_tests.jl")
include("bvk2_oversample_tests.jl")
include("opf_coefficient_correction_tests.jl")
include("opf_tests.jl")
include("unified_lamellar_benchmark_tests.jl")
include("lamellar_cell_stress_mechanism_tests.jl")
include("bvk2_fixed_stiffness_figure6_tests.jl")
include("triblock_aba_tests.jl")
include("bvk2_polyorder_phase_seed_tests.jl")
include("bvk2_ud_theta_extended_phase_tests.jl")
include("bvk2_extended_cell_root_tests.jl")
include("bvk2_aba_comprehensive_validation_tests.jl")
include("bvk2_scft_cost_benchmark_tests.jl")
include("bvk2_scft_fixedcell_performance_tests.jl")
@testset "BVK2 manuscript reproducibility package" begin
    python_test = joinpath(
        @__DIR__, "test_write_bvk2_manuscript_reproducibility_package.py")
    run(`python3 $python_test`)
end
@testset "Macromolecules manuscript drafts" begin
    python_test = joinpath(
        @__DIR__, "test_macromolecules_manuscript_drafts.py")
    run(`python3 $python_test`)
end
@testset "BVK2 manuscript Figure 1 and TOC" begin
    python_test = joinpath(
        @__DIR__, "test_write_bvk2_figure1_and_toc.py")
    run(`python3 $python_test`)
end
@testset "BVK2 SI convergence and morphology packet" begin
    python_test = joinpath(
        @__DIR__, "test_write_bvk2_si_convergence_morphology.py")
    run(`python3 $python_test`)
end
include("gc_gl_boundary_campaign_tests.jl")
include("gc_gl_lowchi_continuation_tests.jl")
include("bvk1_gc_gl_fixed_f_continuation_tests.jl")
include("bvk1_sdis_fixed_chi_tests.jl")
include("certify_bvk1_sdis_endpoint_tests.jl")
include("bvk1_sc_highchi_fixed_chi_tests.jl")
include("bvk1_sc_true_fixed_chi_tests.jl")

@testset "UD phase diagram digitization source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_phase_diagram_digitization.jl")
    source = read(script_path, String)
    @test occursin("UD_AXIS_XMIN = 119.0", source)
    @test occursin("UD_AXIS_XMAX = 602.0", source)
    @test occursin("UD_AXIS_Y_CHIN0 = 606.0", source)
    @test occursin("UD_AXIS_Y_CHIN120 = 141.0", source)
    @test occursin("left_G_arrow_tip", source)
    @test occursin("current_f0.39_chiN20", source)
    @test occursin("proposed_left_g_interior_0.41_24", source)
    @test occursin("proposed_left_g_interior_0.42_24", source)
end

@testset "UD gyroid interior target source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_interior_target_scan.jl")
    source = read(script_path, String)
    target_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_target_resolution_scan.jl"), String)
    @test occursin("f0.41_chiN24_g8:0.41:24.0:8x8x8", source)
    @test occursin("f0.42_chiN24_g8:0.42:24.0:8x8x8", source)
    @test occursin("0.94,1.0,1.06,1.12,1.18", source)
    @test occursin("write_uneyama_doi_gyroid_target_resolution_scan.jl", source)
    @test occursin("Base.include(@__MODULE__", source)
    @test occursin("--lam-iterations", source)
    @test occursin("--cyl-dims", source)
    @test occursin("--group-balance-tolerance", source)
    @test occursin("stationarity_check_pass=stationarity", target_source)
    @test occursin("projected_gradient_norm <= Float64(projected_tolerance)",
        target_source)
end

@testset "UD gyroid high-grid polish source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_high_grid_polish.jl")
    source = read(script_path, String)
    @test occursin("f0.42_chiN24_g10:0.42:24.0:10x10x10", source)
    @test occursin("1.06,1.08,1.10,1.12,1.14,1.16,1.18", source)
    @test occursin("--warm-iterations\", \"1200\"", source)
    @test occursin("--final-iterations\", \"1800\"", source)
    @test occursin("--projected-tolerance\", \"0.002\"", source)
    @test occursin("write_uneyama_doi_gyroid_target_resolution_scan.jl", source)
    @test occursin("Base.include(@__MODULE__", source)
end

@testset "UD checkpointed gyroid polish source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_checkpointed_polish.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    target_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_target_resolution_scan.jl"), String)
    @test occursin("DFM_SKIP_UD_GYR_TARGET_RESOLUTION_MAIN", target_source)
    @test occursin("f0.42_chiN24_g10:0.42:24.0:10x10x10", source)
    @test occursin("0.94,1.00,1.06", source)
    @test occursin("ud_gyroid_checkpointed_polish_gyr_rows.csv", source)
    @test occursin("ud_gyroid_checkpointed_polish_competitor_rows.csv", source)
    @test occursin("ud_gyroid_checkpointed_polish_cyl_evaluations.csv", source)
    @test occursin("_completed_scale_factors", source)
    @test occursin("_completed_competitor_morphologies", source)
    @test occursin("_checkpointed_cylinder_candidate", source)
    @test occursin("_checkpoint_competitors!", source)
    @test occursin("--competitor-source", source)
    @test occursin("_load_competitor_rows", source)
    @test occursin("checkpointed UD competitor", source)
    @test occursin("checkpointed UD CYL factor", source)
    @test occursin("skip completed CYL factor", source)
    @test occursin("skip completed", source)
    @test occursin("CSV.write(gyr_csv", source)
    @test occursin("CSV.write(competitor_csv", source)
    @test occursin("CSV.write(cyl_eval_csv", source)
    @test occursin("Target._annotate_gyr_scan", source)
    @test occursin("Target._write_markdown", source)
end

@testset "UD gyroid nonlinear identity probe source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_nonlinear_identity_probe.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("f0.378_chiN28_g10:0.378:28.0:10x10x10", source)
    @test occursin("1.00,1.06,1.12,1.18,1.24,1.30,1.36,1.42",
        source)
    @test occursin("ud_gyroid_nonlinear_identity_probe_rows.csv", source)
    @test occursin("nonlinear_gyroid_correlation", source)
    @test occursin("_sigmoid_reference", source)
    @test occursin("_best_nonlinear_reference_correlation", source)
    @test occursin("accepted_by_nonlinear_identity", source)
    @test occursin("first-star", source)
    @test occursin("Target._gyr_warm_steps", source)
    packet_script = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("ud_gyroid_nonlinear_identity_probe", packet_script)
    @test occursin("SI Fig./Table 35", packet_script)
    @test occursin("SI Data 36", packet_script)
    @test occursin("polyorder_gyroid_reference_format", packet_script)
    @test occursin("hard first-star classifier", packet_script)
end

@testset "UD gyroid Polyorder reference probe source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_polyorder_reference_probe.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("Polyorder_results", source)
    @test occursin("GYR_opt.mat", source)
    @test occursin("using MAT", source)
    @test occursin("_load_polyorder_gyroid_reference", source)
    @test occursin("_downsample_periodic_mean", source)
    @test occursin("_best_oriented_shifted_correlation", source)
    @test occursin("polyorder_gyroid_reference_probe_rows.csv", source)
    @test occursin("polyorder_gyroid_reference_probe_summary.csv", source)
    @test occursin("polyorder_scf_gyr_correlation", source)
    @test occursin("analytic_sigmoid_gyr_correlation", source)
    @test occursin("analytic_sigmoid_bcc_correlation", source)
    @test occursin("field_dominant_shell", source)
    @test occursin("polyorder_reference_requires_lattice_convention_mapping",
        source)
    @test occursin("PRIMITIVE_TO_CONVENTIONAL_RECIPROCAL", source)
    @test occursin("_mapped_star_shell_label", source)
    @test occursin("_mapped_star_shell_counts", source)
    @test occursin("_mapped_shell_phase_metrics", source)
    @test occursin("maps_110_to_211", source)
    @test occursin("mapped_star_shell_counts", source)
    @test occursin("mapped_shell_phase_residual", source)
    @test occursin("mapped_shell_phase_similarity", source)
    @test occursin("_shell_invariant_gyroid_support_metrics", source)
    @test occursin("topology_label", source)
    @test occursin("shell_invariant_support_label", source)
    @test occursin("dominant_shell_effective_mode_count", source)
    @test occursin("a_rich_euler_characteristic", source)
    @test occursin("b_rich_euler_characteristic", source)
    @test occursin("negative_euler_double_network", source)
    @test occursin("fA=0.36", source)
    @test occursin("chiN=16", source)
    @test occursin("60x60x60", source)
    packet_script = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("polyorder_gyroid_reference_probe", packet_script)
    @test occursin("SI Fig./Table 37", packet_script)
    audit_script = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_extension_audit.jl"), String)
    @test occursin("has_ud_gyroid_polyorder_reference_probe", audit_script)
    @test occursin("shell-invariant gyroid classifier", audit_script)
end

@testset "UD gyroid Polyorder-initialized saddle probe source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_polyorder_initialized_probe.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("Polyorder_results", source)
    @test occursin("GYR_opt.mat", source)
    @test occursin("using MAT", source)
    @test occursin("_load_polyorder_gyroid_reference", source)
    @test occursin("_downsample_periodic_mean", source)
    @test occursin("_mean_shift_density", source)
    @test occursin("_polyorder_initial_theta", source)
    @test occursin("relax_diblock_uneyama_doi_morphology_seed_nd", source)
    @test occursin("_shell_invariant_gyroid_support_metrics", source)
    @test occursin("POLYORDER_211_TO_UD_110_LENGTH_RATIO", source)
    @test occursin("polyorder_initialized_gyr211_pass", source)
    @test occursin("polyorder_initialized_gyr_support_pass", source)
    @test occursin("ud_gyroid_polyorder_initialized_rows.csv", source)
    @test occursin("ud_gyroid_polyorder_initialized_summary.csv", source)
    @test occursin("DFM_SKIP_UD_GYR_POLYORDER_INITIALIZED_MAIN", source)
    if isfile(script_path)
        mod = Module(:UDGyroidPolyorderInitializedProbeTest)
        key = "DFM_SKIP_UD_GYR_POLYORDER_INITIALIZED_MAIN"
        previous = get(ENV, key, nothing)
        try
            ENV[key] = "1"
            Base.include(mod, script_path)
        finally
            if previous === nothing
                delete!(ENV, key)
            else
                ENV[key] = previous
            end
        end
        seed = reshape(collect(range(0.02, 0.98; length=64)), 4, 4, 4)
        mean_shift_density = Base.invokelatest(getproperty, mod, Symbol("_mean_shift_density"))
        shifted = Base.invokelatest(mean_shift_density, seed, 0.378)
        @test size(shifted) == size(seed)
        @test all(value -> 0.0 < value < 1.0, shifted)
        @test mean(shifted) ≈ 0.378 atol=1.0e-10
    end
end

@testset "BURP-BVK1 same-parameter CYL competitor source contract" begin
    reference_script = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_reference_comparison.jl")
    unrestricted_script = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_unrestricted_relaxation.jl")
    cyl_logit_script = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_cyl_logit_cell_scan.jl")
    gate_script = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_representative_phase_rank_gate.jl")

    reference_source = read(reference_script, String)
    unrestricted_source = read(unrestricted_script, String)
    cyl_logit_source = read(cyl_logit_script, String)
    gate_source = read(gate_script, String)

    @test occursin("--cyl-dir", reference_source)
    @test occursin("f=cyl_f, chiN=cyl_chiN", reference_source)
    @test occursin("--cyl-dir", unrestricted_source)
    @test occursin("f=cyl_f, chiN=cyl_chiN", unrestricted_source)
    @test occursin("haskey(selected, (MORPHOLOGY, model))",
        cyl_logit_source)
    @test occursin("start_label == \"branch_selected\" && selected_row === nothing",
        cyl_logit_source)
    @test occursin("same_parameter_cyl_competitor", gate_source)
    @test occursin("cyl_same_parameter_logit_cell_scan_summary.csv",
        gate_source)
end

@testset "BURP-BVK1 GYR phase-rank energy panel source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyr_phase_rank_energy_panel.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("gyr_phase_rank_energy_panel.svg", source)
    @test occursin("gyr_phase_rank_energy_panel_rows.csv", source)
    @test occursin("gyr_phase_rank_energy_panel_summary.csv", source)
    @test occursin("GYR_BURP_TI", source)
    @test occursin("GYR_BVK1", source)
    @test occursin("energy_margin_to_next_accepted", source)
    @test occursin("endpoint-limited BVK1 CYL", source)
    @test occursin("single representative BURP-TI GYR phase-rank pass",
        source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("gyr_phase_rank_energy_panel", packet_source)
    @test occursin("SI Fig./Table 0g", packet_source)
end

@testset "BURP-BVK1 phase-rank experiment manifest source contract" begin
    target_selector_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_nonlamellar_phase_rank_target_selector.jl")
    @test isfile(target_selector_path)
    target_selector_source = isfile(target_selector_path) ?
        read(target_selector_path, String) : ""
    @test occursin("ud_cyl_bcc_method_equity_gate.csv",
        target_selector_source)
    @test occursin("ud_method_equity_status", target_selector_source)
    @test occursin(
        "ud_cyl_bcc_representative_method_equity_pass_with_endpoint_caveat",
        target_selector_source)
    @test occursin("target_model_blocker_after_ud_control_pass",
        target_selector_source)
    @test occursin("ud_gyr_open", target_selector_source)

    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_phase_rank_experiment_manifest.jl")
    @test isfile(script_path)
    source = read(script_path, String)
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("phase_rank_experiment_manifest.csv", source)
    @test occursin("phase_rank_experiment_queue.csv", source)
    @test occursin("phase_rank_experiment_manifest.svg", source)
    @test occursin("Polyorder_results/BCC/BCC_opt.mat", source)
    @test occursin("Polyorder_results/CYL/CYL_opt.mat", source)
    @test occursin("Polyorder_results/CYL_f0.36_chiN16/CYL_opt.mat", source)
    @test occursin("source_kind=user_provided_polyorder_3d_seed", source)
    @test occursin("source_kind=project_generated_polyorder_2d_random_seed", source)
    @test occursin("ud_cyl_bcc_method_equity_gate.csv", source)
    @test occursin("ud_method_equity_status", source)
    @test occursin("UD CYL/BCC control pass", source)
    @test occursin("target-model blocker after UD control pass", source)
    @test occursin("GYR_BVK1_same_parameter_competitor_unblock", source)
    @test occursin("CYL_2D_polyorder_random_seed_route", source)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate_summary.csv",
        source)
    @test occursin("burp_bvk1_ud_neighbor_f0p41_chiN24_phase_rank_gate",
        source)
    @test occursin("burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate",
        source)
    @test occursin("phase_rank_gate_summary.csv", source)
    @test occursin("phase_rank_completed", source)
    @test occursin("GYR_ud_neighbor_f0p41_chiN24", source)
    @test occursin("BCC_user_seed_finite_contrast_route", source)
    @test occursin("stop_condition", source)
    @test occursin("write_burp_bvk1_phase_rank_experiment_manifest.jl",
        pipeline_source)
    @test occursin("phase_rank_experiment_manifest", packet_source)
    @test occursin("SI Fig./Table 0d2c", packet_source)
    manifest_csv = read(joinpath(@__DIR__, "..", "results",
        "burp_bvk1_phase_rank_experiment_manifest",
        "phase_rank_experiment_manifest.csv"), String)
    queue_csv = read(joinpath(@__DIR__, "..", "results",
        "burp_bvk1_phase_rank_experiment_manifest",
        "phase_rank_experiment_queue.csv"), String)
    manifest_lines = split(chomp(manifest_csv), '\n')
    queue_lines = split(chomp(queue_csv), '\n')
    @test occursin("GYR_ud_target_f0p39_chiN20", manifest_csv)
    @test occursin("GYR_ud_neighbor_f0p41_chiN24", manifest_csv)
    @test occursin("GYR_ud_neighbor_f0p42_chiN24", manifest_csv)
    @test occursin("CYL_2D_polyorder_random_seed_route", manifest_csv)
    @test occursin("burp_bvk1_ud_neighbor_f0p41_chiN24_phase_rank_gate",
        manifest_csv)
    @test occursin("burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate",
        manifest_csv)
    @test occursin("phase_rank_completed", manifest_csv)
    @test !occursin("GYR_ud_target_f0p39_chiN20", queue_csv)
    @test !occursin("GYR_ud_neighbor_f0p41_chiN24", queue_csv)
    @test !occursin("GYR_ud_neighbor_f0p42_chiN24", queue_csv)
    @test !occursin("CYL_2D_polyorder_random_seed_route", queue_csv)
    @test length(manifest_lines) >= 2
    @test startswith(manifest_lines[2], "0,GYR_ud_target_f0p39_chiN20")
    @test length(queue_lines) >= 2
    @test startswith(queue_lines[2], "1,BCC_user_seed_finite_contrast_route")
    @test occursin("next_bcc_finite_contrast_guardrail", queue_lines[2])
    @test occursin("historical_f0p36_gyr_bvk1_cyl_blocker", queue_csv)
    @test !occursin("next_reproducibility_target", queue_lines[2])
end

@testset "BURP-BVK1 UD-target SCFT reference packet source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_target_scft_reference_packet.jl")
    @test isfile(script_path)
    source = read(script_path, String)
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("ud_target_scft_reference_packet.csv", source)
    @test occursin("ud_target_scft_reference_packet.svg", source)
    @test occursin("Polyorder_results/CYL_f0.39_chiN20/CYL_opt.mat", source)
    @test occursin("Polyorder_results/gyroid/GYR_opt.mat", source)
    @test occursin("Polyorder_results/GYR_f0.39_chiN20/GYR_opt.mat", source)
    @test occursin("GYR_ud_target_f0p39_chiN20", source)
    @test occursin("project_generated_polyorder_2d_random_seed", source)
    @test occursin("gyroid_warm_start_from_existing_polyorder_seed", source)
    @test occursin("same_parameter_phase_rank_inputs", source)
    @test occursin("generated CYL, LAM, and DIS", source)
    @test occursin("ready_same_parameter_lam_competitors", source)
    @test occursin("ready_same_parameter_dis_baseline", source)
    @test occursin("write_burp_bvk1_ud_target_scft_reference_packet.jl",
        pipeline_source)
    @test occursin("ud_target_scft_reference_packet", packet_source)
    @test occursin("SI Fig./Table 0d2d", packet_source)
    @test occursin("Polyorder_results/CYL_f0.39_chiN20/CYL_opt.mat",
        packet_source)
    @test !occursin("../Polyorder_results/CYL_f0.39_chiN20", packet_source)
end

@testset "BURP-BVK1 GYR-neighbor reproducibility packet source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyr_neighbor_reproducibility_packet.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("gyr_neighbor_reproducibility_packet.csv", source)
    @test occursin("gyr_neighbor_reproducibility_summary.csv", source)
    @test occursin("gyr_neighbor_reproducibility_packet.svg", source)
    @test occursin("gyr_neighbor_reproducibility_packet.md", source)
    @test occursin("phase_rank_experiment_queue.csv", source)
    @test occursin("GYR_ud_neighbor_f0p41_chiN24", source)
    @test occursin("GYR_ud_neighbor_f0p42_chiN24", source)
    @test occursin("completed_f0p39_positive_control", source)
    @test occursin("Polyorder_results/GYR_f0.41_chiN24/GYR_opt.mat",
        source)
    @test occursin("Polyorder_results/CYL_f0.41_chiN24/CYL_opt.mat",
        source)
    @test occursin("polyorder_hexrect_reference_required", source)
    @test occursin("strict DIS/LAM/CYL/GYR", source)
    @test occursin("write_burp_bvk1_gyr_neighbor_reproducibility_packet.jl",
        pipeline_source)
    @test occursin("gyr_neighbor_reproducibility_packet", packet_source)
    @test occursin("SI Fig./Table 0d2f", packet_source)
    packet_csv_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_gyr_neighbor_reproducibility_packet",
        "gyr_neighbor_reproducibility_packet.csv")
    summary_csv_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_gyr_neighbor_reproducibility_packet",
        "gyr_neighbor_reproducibility_summary.csv")
    @test isfile(packet_csv_path)
    @test isfile(summary_csv_path)
    packet_csv = isfile(packet_csv_path) ? read(packet_csv_path, String) : ""
    summary_csv = isfile(summary_csv_path) ? read(summary_csv_path, String) : ""
    @test occursin("GYR_ud_neighbor_f0p41_chiN24", packet_csv)
    @test occursin("Polyorder_results/GYR_f0.41_chiN24/GYR_opt.mat",
        packet_csv)
    @test occursin("Polyorder_results/CYL_f0.41_chiN24/CYL_opt.mat",
        packet_csv)
    @test occursin("polyorder_hexrect_reference_required", packet_csv)
    @test occursin("phase_rank_completed", summary_csv)
    @test occursin("Retain as completed GYR-neighbor evidence", summary_csv)
end

@testset "BURP-BVK1 UD-target f0.39 chiN20 phase-rank gate source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_target_f0p39_chiN20_phase_rank_gate.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("GYR_ud_target_f0p39_chiN20", source)
    @test occursin("gyr_logit_cell_scan_summary.csv", source)
    @test occursin("cyl_logit_cell_scan_summary.csv", source)
    @test occursin("lam_competitors.csv", source)
    @test occursin("dis_baseline.csv", source)
    @test occursin("Polyorder_results/CYL_f0.39_chiN20/CYL_opt.mat", source)
    @test occursin("stress_free_cell_candidate", source)
    @test occursin("missing_target_model_cyl_competitor", source)
    @test occursin("cell_size_energy_minimum", source)
    @test occursin("cell_stress_proxy", source)
    @test occursin("target_cyl_status", source)
    @test occursin("selected_morphology", source)
    @test occursin("lowest_accepted_energy_morphology", source)
    @test occursin("refined_energy_density", source)
    @test occursin("refined_cell_factor", source)
    @test occursin("cyl_geometry_admissible", source)
    @test occursin("geometry_convention_invalid", source)
    @test occursin("ud_cyl_bcc_method_equity_gate.csv", source)
    @test occursin("UD-first", source)
    @test occursin("numerical-method evidence first", source)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate.csv", source)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate_summary.csv", source)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate.svg", source)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate.md", source)
end

@testset "Polyorder GYR f0.39 chiN20 warm-start source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "run_polyorder_gyr_f0p39_chiN20_warm_start.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("Polyorder_results/gyroid/GYR_opt.mat", source)
    @test occursin("Polyorder_results/GYR_f0.39_chiN20/GYR_opt.mat", source)
    @test occursin("gyroid_warm_start_from_existing_polyorder_seed", source)
    @test occursin("fA = 0.39", source)
    @test occursin("chiN = 20.0", source)
    @test occursin("Polyorder.initialize!", source)
    @test occursin("Polyorder.cell_solve!", source)
    @test occursin("Polyorder.VariableCell", source)
    @test occursin("dominant_shell", source)
    @test occursin("{211}", source)
    @test occursin("findmax(label_power)[2]", source)
    @test !occursin("first(findmax(label_power)[2])", source)
    @test occursin("\"component\"", source)
    @test occursin("\"wfields\"", source)
    @test occursin("run_polyorder_gyr_f0p39_chiN20_warm_start.jl",
        pipeline_source)
end

@testset "Polyorder f0.41 chiN24 neighbor reference runners source contract" begin
    gyr_script_path = joinpath(@__DIR__, "..", "scripts",
        "run_polyorder_gyr_f0p41_chiN24_warm_start.jl")
    cyl_script_path = joinpath(@__DIR__, "..", "scripts",
        "run_polyorder_cyl_f0p41_chiN24_hexrect.jl")
    @test isfile(gyr_script_path)
    @test isfile(cyl_script_path)
    gyr_source = isfile(gyr_script_path) ? read(gyr_script_path, String) : ""
    cyl_source = isfile(cyl_script_path) ? read(cyl_script_path, String) : ""
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyr_neighbor_reproducibility_packet.jl"), String)
    @test occursin("TARGET_FA = 0.41", gyr_source)
    @test occursin("TARGET_CHIN = 24.0", gyr_source)
    @test occursin("Polyorder_results/GYR_f0.41_chiN24/GYR_opt.mat",
        gyr_source)
    @test occursin("Polyorder_results/GYR_f0.39_chiN20/GYR_opt.mat",
        gyr_source)
    @test occursin("gyroid_neighbor_warm_start_from_f0p39_control",
        gyr_source)
    @test occursin("DEFAULT_INITIAL_A = 9.7559158", gyr_source)
    @test occursin("Polyorder.initialize!", gyr_source)
    @test occursin("Polyorder.cell_solve!", gyr_source)
    @test occursin("{211}", gyr_source)
    @test occursin("TARGET_FA = 0.41", cyl_source)
    @test occursin("TARGET_CHIN = 24.0", cyl_source)
    @test occursin("Polyorder_results/CYL_f0.41_chiN24/CYL_opt.mat",
        cyl_source)
    @test occursin("polyorder_hexrect_reference_required", cyl_source)
    @test occursin("Polyorder.UnitCell(Polyorder.HexRect()", cyl_source)
    @test occursin("DEFAULT_INITIAL_A = 4.62679", cyl_source)
    @test occursin("DEFAULT_SEED = 4101", cyl_source)
    @test occursin("Polyorder.cell_solve!", cyl_source)
    @test occursin("run_polyorder_gyr_f0p41_chiN24_warm_start.jl",
        pipeline_source)
    @test occursin("run_polyorder_cyl_f0p41_chiN24_hexrect.jl",
        pipeline_source)
    @test occursin("run_polyorder_gyr_f0p41_chiN24_warm_start.jl",
        packet_source)
    @test occursin("run_polyorder_cyl_f0p41_chiN24_hexrect.jl",
        packet_source)
end

@testset "BURP-BVK1 f0.41 chiN24 neighbor LAM/DIS competitor source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_neighbor_f0p41_chiN24_lam_dis_competitors.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyr_neighbor_reproducibility_packet.jl"), String)
    @test occursin("TARGET_ID = \"GYR_ud_neighbor_f0p41_chiN24\"", source)
    @test occursin("TARGET_FA = 0.41", source)
    @test occursin("TARGET_CHIN = 24.0", source)
    @test occursin("Polyorder_results/GYR_f0.41_chiN24/GYR_opt.mat",
        source)
    @test occursin("Polyorder_results/CYL_f0.41_chiN24/CYL_opt.mat",
        source)
    @test occursin("polyorder_hexrect_reference_required", source)
    @test occursin("minimize_diblock_burp_lamella_stress_free", source)
    @test occursin("minimize_diblock_bvk1_lamella_stress_free", source)
    @test occursin("accepted_lamellar_competitor", source)
    @test occursin("period_local_minimum_check_pass", source)
    @test occursin("energy_density_proxy", source)
    @test occursin("lam_competitors.csv", source)
    @test occursin("lam_competitors.svg", source)
    @test occursin("lam_competitors.md", source)
    @test occursin("dis_baseline.csv", source)
    @test occursin("homogeneous_disordered_baseline", source)
    @test occursin("write_burp_bvk1_ud_neighbor_f0p41_chiN24_lam_dis_competitors.jl",
        pipeline_source)
    @test occursin("results/burp_bvk1_ud_neighbor_f0p41_chiN24_lam_competitors/lam_competitors.csv",
        packet_source)
    @test occursin("results/burp_bvk1_ud_neighbor_f0p41_chiN24_dis_baseline/dis_baseline.csv",
        packet_source)

    lam_csv_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_ud_neighbor_f0p41_chiN24_lam_competitors",
        "lam_competitors.csv")
    dis_csv_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_ud_neighbor_f0p41_chiN24_dis_baseline",
        "dis_baseline.csv")
    @test isfile(lam_csv_path)
    @test isfile(dis_csv_path)
    lam_csv = isfile(lam_csv_path) ? read(lam_csv_path, String) : ""
    dis_csv = isfile(dis_csv_path) ? read(dis_csv_path, String) : ""
    @test occursin("GYR_BURP_TI", lam_csv)
    @test occursin("GYR_BVK1", lam_csv)
    @test occursin("accepted_lamellar_competitor", lam_csv)
    @test occursin("period_local_minimum_check_pass", lam_csv)
    @test occursin("GYR_ud_neighbor_f0p41_chiN24", dis_csv)
    @test occursin("homogeneous_disordered_baseline", dis_csv)
end

@testset "BURP-BVK1 f0.41 chiN24 neighbor GYR/CYL scan source contract" begin
    gyr_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_neighbor_f0p41_chiN24_gyr_logit_cell_scan.jl")
    cyl_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_neighbor_f0p41_chiN24_cyl_logit_cell_scan.jl")
    @test isfile(gyr_script_path)
    @test isfile(cyl_script_path)
    gyr_source = isfile(gyr_script_path) ? read(gyr_script_path, String) : ""
    cyl_source = isfile(cyl_script_path) ? read(cyl_script_path, String) : ""
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyr_neighbor_reproducibility_packet.jl"), String)

    @test occursin("TARGET_ID = \"GYR_ud_neighbor_f0p41_chiN24\"",
        gyr_source)
    @test occursin("TARGET_FA = 0.41", gyr_source)
    @test occursin("TARGET_CHIN = 24.0", gyr_source)
    @test occursin("Polyorder_results/GYR_f0.41_chiN24/GYR_opt.mat",
        gyr_source)
    @test occursin("--gyr-target-dims=14,14,14", gyr_source)
    @test occursin("--starts=scft_point", gyr_source)
    @test occursin("--cell-factors=0.70,0.78,0.86,0.94,1.02,1.10,1.18,1.26,1.34",
        gyr_source)
    @test occursin("write_burp_bvk1_polyorder_gyr_logit_cell_scan.jl",
        gyr_source)
    @test occursin("gyr_logit_cell_scan_summary.csv", gyr_source)

    @test occursin("TARGET_ID = \"GYR_ud_neighbor_f0p41_chiN24\"",
        cyl_source)
    @test occursin("TARGET_FA = 0.41", cyl_source)
    @test occursin("TARGET_CHIN = 24.0", cyl_source)
    @test occursin("Polyorder_results/CYL_f0.41_chiN24/CYL_opt.mat",
        cyl_source)
    @test occursin("polyorder_hexrect_reference_required", cyl_source)
    @test occursin("--cyl-target-dims=16,28", cyl_source)
    @test occursin("--starts=scft_point", cyl_source)
    @test occursin("--cell-factors=0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.00,1.02,1.04,1.06,1.08,1.10,1.12",
        cyl_source)
    @test occursin("write_burp_bvk1_polyorder_cyl_logit_cell_scan.jl",
        cyl_source)
    @test occursin("cyl_logit_cell_scan_summary.csv", cyl_source)

    @test occursin("write_burp_bvk1_ud_neighbor_f0p41_chiN24_gyr_logit_cell_scan.jl",
        pipeline_source)
    @test occursin("write_burp_bvk1_ud_neighbor_f0p41_chiN24_cyl_logit_cell_scan.jl",
        pipeline_source)
    @test occursin("results/burp_bvk1_ud_neighbor_f0p41_chiN24_gyr_logit_cell_scan/gyr_logit_cell_scan_summary.csv",
        packet_source)
    @test occursin("results/burp_bvk1_ud_neighbor_f0p41_chiN24_cyl_logit_cell_scan/cyl_logit_cell_scan_summary.csv",
        packet_source)

    gyr_summary_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_ud_neighbor_f0p41_chiN24_gyr_logit_cell_scan",
        "gyr_logit_cell_scan_summary.csv")
    cyl_summary_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_ud_neighbor_f0p41_chiN24_cyl_logit_cell_scan",
        "cyl_logit_cell_scan_summary.csv")
    @test isfile(gyr_summary_path)
    @test isfile(cyl_summary_path)
    gyr_summary = isfile(gyr_summary_path) ? read(gyr_summary_path, String) : ""
    cyl_summary = isfile(cyl_summary_path) ? read(cyl_summary_path, String) : ""
    @test occursin("Polyorder_results/GYR_f0.41_chiN24/GYR_opt.mat",
        gyr_summary)
    @test occursin("Polyorder_results/CYL_f0.41_chiN24/CYL_opt.mat",
        cyl_summary)
    @test occursin("stress_free_cell_candidate", gyr_summary)
    @test occursin("phase_rank_cyl_candidate", cyl_summary)
    @test occursin("polyorder_hexrect_sqrt3", cyl_summary)
end

@testset "BURP-BVK1 f0.41 chiN24 neighbor phase-rank gate source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_neighbor_f0p41_chiN24_phase_rank_gate.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyr_neighbor_reproducibility_packet.jl"), String)
    @test occursin("TARGET_ID = \"GYR_ud_neighbor_f0p41_chiN24\"", source)
    @test occursin("TARGET_FA = 0.41", source)
    @test occursin("TARGET_CHIN = 24.0", source)
    @test occursin("REQUIRED_COMPETITORS = (\"DIS\", \"LAM\", \"CYL\", \"GYR\")",
        source)
    @test occursin("burp_bvk1_ud_neighbor_f0p41_chiN24_dis_baseline",
        source)
    @test occursin("burp_bvk1_ud_neighbor_f0p41_chiN24_lam_competitors",
        source)
    @test occursin("burp_bvk1_ud_neighbor_f0p41_chiN24_cyl_logit_cell_scan",
        source)
    @test occursin("burp_bvk1_ud_neighbor_f0p41_chiN24_gyr_logit_cell_scan",
        source)
    @test occursin("phase_rank_cyl_candidate", source)
    @test occursin("stress_free_cell_candidate", source)
    @test occursin("polyorder_hexrect_sqrt3", source)
    @test occursin("phase_rank_gate_summary.csv", source)
    @test occursin("phase_rank_pass", source)
    @test occursin("write_burp_bvk1_ud_neighbor_f0p41_chiN24_phase_rank_gate.jl",
        pipeline_source)
    @test occursin("results/burp_bvk1_ud_neighbor_f0p41_chiN24_phase_rank_gate/phase_rank_gate_summary.csv",
        packet_source)

    summary_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_ud_neighbor_f0p41_chiN24_phase_rank_gate",
        "phase_rank_gate_summary.csv")
    gate_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_ud_neighbor_f0p41_chiN24_phase_rank_gate",
        "phase_rank_gate.csv")
    @test isfile(summary_path)
    @test isfile(gate_path)
    summary = isfile(summary_path) ? read(summary_path, String) : ""
    gate = isfile(gate_path) ? read(gate_path, String) : ""
    @test occursin("GYR_ud_neighbor_f0p41_chiN24", summary)
    @test occursin("BURP-TI", summary)
    @test occursin("BVK1", summary)
    @test occursin("phase_rank_pass", summary)
    @test occursin("selected_morphology", summary)
    @test occursin("accepted_competitor", gate)
    @test occursin("polyorder_hexrect_sqrt3", gate)
end

@testset "BURP-BVK1 f0.42 chiN24 neighbor production path contract" begin
    gyr_ref_script = joinpath(@__DIR__, "..", "scripts",
        "run_polyorder_gyr_f0p42_chiN24_warm_start.jl")
    cyl_ref_script = joinpath(@__DIR__, "..", "scripts",
        "run_polyorder_cyl_f0p42_chiN24_hexrect.jl")
    lam_dis_script = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_neighbor_f0p42_chiN24_lam_dis_competitors.jl")
    gyr_scan_script = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_neighbor_f0p42_chiN24_gyr_logit_cell_scan.jl")
    cyl_scan_script = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_neighbor_f0p42_chiN24_cyl_logit_cell_scan.jl")
    phase_gate_script = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate.jl")
    @test isfile(gyr_ref_script)
    @test isfile(cyl_ref_script)
    @test isfile(lam_dis_script)
    @test isfile(gyr_scan_script)
    @test isfile(cyl_scan_script)
    @test isfile(phase_gate_script)

    gyr_ref_source = isfile(gyr_ref_script) ? read(gyr_ref_script, String) : ""
    cyl_ref_source = isfile(cyl_ref_script) ? read(cyl_ref_script, String) : ""
    lam_dis_source = isfile(lam_dis_script) ? read(lam_dis_script, String) : ""
    gyr_scan_source = isfile(gyr_scan_script) ? read(gyr_scan_script, String) : ""
    cyl_scan_source = isfile(cyl_scan_script) ? read(cyl_scan_script, String) : ""
    phase_gate_source = isfile(phase_gate_script) ? read(phase_gate_script, String) : ""
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyr_neighbor_reproducibility_packet.jl"), String)

    @test occursin("TARGET_FA = 0.42", gyr_ref_source)
    @test occursin("TARGET_CHIN = 24.0", gyr_ref_source)
    @test occursin("Polyorder_results/GYR_f0.42_chiN24/GYR_opt.mat",
        gyr_ref_source)
    @test occursin("Polyorder_results/GYR_f0.41_chiN24/GYR_opt.mat",
        gyr_ref_source)
    @test occursin("gyroid_neighbor_warm_start_from_f0p41_control",
        gyr_ref_source)
    @test occursin("DEFAULT_INITIAL_A = 10.413026", gyr_ref_source)
    @test occursin("Polyorder.initialize!", gyr_ref_source)
    @test occursin("{211}", gyr_ref_source)

    @test occursin("TARGET_FA = 0.42", cyl_ref_source)
    @test occursin("TARGET_CHIN = 24.0", cyl_ref_source)
    @test occursin("Polyorder_results/CYL_f0.42_chiN24/CYL_opt.mat",
        cyl_ref_source)
    @test occursin("polyorder_hexrect_reference_required", cyl_ref_source)
    @test occursin("Polyorder.UnitCell(Polyorder.HexRect()", cyl_ref_source)
    @test occursin("DEFAULT_INITIAL_A = 4.8747042", cyl_ref_source)
    @test occursin("DEFAULT_SEED = 4201", cyl_ref_source)

    @test occursin("TARGET_ID = \"GYR_ud_neighbor_f0p42_chiN24\"",
        lam_dis_source)
    @test occursin("TARGET_FA = 0.42", lam_dis_source)
    @test occursin("TARGET_CHIN = 24.0", lam_dis_source)
    @test occursin("Polyorder_results/GYR_f0.42_chiN24/GYR_opt.mat",
        lam_dis_source)
    @test occursin("Polyorder_results/CYL_f0.42_chiN24/CYL_opt.mat",
        lam_dis_source)
    @test occursin("lam_competitors.csv", lam_dis_source)
    @test occursin("dis_baseline.csv", lam_dis_source)

    @test occursin("TARGET_ID = \"GYR_ud_neighbor_f0p42_chiN24\"",
        gyr_scan_source)
    @test occursin("Polyorder_results/GYR_f0.42_chiN24/GYR_opt.mat",
        gyr_scan_source)
    @test occursin("--gyr-target-dims=14,14,14", gyr_scan_source)
    @test occursin("gyr_logit_cell_scan_summary.csv", gyr_scan_source)

    @test occursin("TARGET_ID = \"GYR_ud_neighbor_f0p42_chiN24\"",
        cyl_scan_source)
    @test occursin("Polyorder_results/CYL_f0.42_chiN24/CYL_opt.mat",
        cyl_scan_source)
    @test occursin("polyorder_hexrect_reference_required", cyl_scan_source)
    @test occursin("--cyl-target-dims=15,25", cyl_scan_source)
    @test occursin("cyl_logit_cell_scan_summary.csv", cyl_scan_source)

    @test occursin("TARGET_ID = \"GYR_ud_neighbor_f0p42_chiN24\"",
        phase_gate_source)
    @test occursin("REQUIRED_COMPETITORS = (\"DIS\", \"LAM\", \"CYL\", \"GYR\")",
        phase_gate_source)
    @test occursin("burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate",
        phase_gate_source)
    @test occursin("phase_rank_gate_summary.csv", phase_gate_source)
    @test occursin("phase_rank_pass", phase_gate_source)

    @test occursin("run_polyorder_gyr_f0p42_chiN24_warm_start.jl",
        pipeline_source)
    @test occursin("run_polyorder_cyl_f0p42_chiN24_hexrect.jl",
        pipeline_source)
    @test occursin("write_burp_bvk1_ud_neighbor_f0p42_chiN24_lam_dis_competitors.jl",
        pipeline_source)
    @test occursin("write_burp_bvk1_ud_neighbor_f0p42_chiN24_gyr_logit_cell_scan.jl",
        pipeline_source)
    @test occursin("write_burp_bvk1_ud_neighbor_f0p42_chiN24_cyl_logit_cell_scan.jl",
        pipeline_source)
    @test occursin("write_burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate.jl",
        pipeline_source)
    @test occursin("Polyorder_results/GYR_f0.42_chiN24/GYR_opt.mat",
        packet_source)
    @test occursin("Polyorder_results/CYL_f0.42_chiN24/CYL_opt.mat",
        packet_source)
    @test occursin("results/burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate/phase_rank_gate_summary.csv",
        packet_source)

    summary_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate",
        "phase_rank_gate_summary.csv")
    gate_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate",
        "phase_rank_gate.csv")
    @test isfile(summary_path)
    @test isfile(gate_path)
    summary = isfile(summary_path) ? read(summary_path, String) : ""
    gate = isfile(gate_path) ? read(gate_path, String) : ""
    @test occursin("GYR_ud_neighbor_f0p42_chiN24", summary)
    @test occursin("BURP-TI", summary)
    @test occursin("BVK1", summary)
    @test occursin("phase_rank_pass", summary)
    @test occursin("selected_morphology", summary)
    @test occursin("accepted_competitor", gate)
    @test occursin("polyorder_hexrect_sqrt3", gate)
end

@testset "3D density field visual renderer" begin
    outdir = mktempdir()
    field = DFMMonteCarlo._diblock_morphology_seed_field(:GYR, (6, 6, 6);
        amplitude=1.0)
    svg_path = DFMMonteCarlo.write_density_field_visual_svg(field,
        joinpath(outdir, "gyr_density_visual.svg");
        title="GYR density visual test", threshold=0.0)
    @test isfile(svg_path)
    text = read(svg_path, String)
    @test occursin("GYR density visual test", text)
    @test occursin("isometric voxel projection", text)
    @test occursin("central XY slice", text)
    @test occursin("central XZ slice", text)
    @test occursin("central YZ slice", text)

    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_density_visuals.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("write_density_field_visual_svg", source)
    @test occursin("ud_gyr_polyorder_initialized_density.csv", source)
    @test occursin("ud_gyr_polyorder_initialized_density_visual.svg", source)
    @test occursin("bcc_reference_density_visual.svg", source)
    @test occursin("DFM_SKIP_UD_GYR_DENSITY_VISUALS_MAIN", source)

    complex_visual_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_complex_structure_density_visuals.jl")
    @test isfile(complex_visual_script_path)
    complex_visual_source = isfile(complex_visual_script_path) ?
        read(complex_visual_script_path, String) : ""
    @test occursin("_read_density_csv", complex_visual_source)
    @test occursin("_render_phase_pair!", complex_visual_source)
    @test occursin("write_density_field_visual_svg", complex_visual_source)
    @test occursin("phase=:above", complex_visual_source)
    @test occursin("phase=:below", complex_visual_source)
    @test occursin("gyr_reference", complex_visual_source)
    @test occursin("bcc_reference", complex_visual_source)
    @test occursin("burp_ti_selected_density.csv", complex_visual_source)
    @test occursin("bvk1_selected_density.csv", complex_visual_source)
    @test occursin("complex_structure_density_visuals.html", complex_visual_source)
    @test occursin("complex_structure_visual_summary.csv", complex_visual_source)
    @test occursin("DFM_SKIP_COMPLEX_STRUCTURE_DENSITY_VISUALS_MAIN",
        complex_visual_source)
    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("complex_structure_density_visuals", packet_source)
    @test occursin("SI Fig./Table 41", packet_source)

    logit_gallery_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_logit_visual_gallery.jl")
    @test isfile(logit_gallery_script_path)
    logit_gallery_source = isfile(logit_gallery_script_path) ?
        read(logit_gallery_script_path, String) : ""
    @test occursin("_load_logit_selected_rows", logit_gallery_source)
    @test occursin("_read_density_csv", logit_gallery_source)
    @test occursin("_write_slice_gallery_svg", logit_gallery_source)
    @test occursin("_render_3d_visuals!", logit_gallery_source)
    @test occursin("nonlamellar_logit_visual_gallery.svg", logit_gallery_source)
    @test occursin("nonlamellar_logit_visual_gallery.html", logit_gallery_source)
    @test occursin("nonlamellar_logit_visual_summary.csv", logit_gallery_source)
    @test occursin("CYL", logit_gallery_source)
    @test occursin("BCC", logit_gallery_source)
    @test occursin("GYR", logit_gallery_source)
    @test occursin("stress_free_cell_candidate", logit_gallery_source)
    @test occursin("endpoint_limited_cell_scan", logit_gallery_source)
    @test occursin(
        "DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_LOGIT_VISUAL_GALLERY_MAIN",
        logit_gallery_source)
    @test occursin("nonlamellar_logit_visual_gallery", packet_source)

    transfer_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_polyorder_initialized_probe.jl")
    @test isfile(transfer_script_path)
    transfer_source = isfile(transfer_script_path) ?
        read(transfer_script_path, String) : ""
    @test occursin("Polyorder_results/gyroid/GYR_opt.mat", transfer_source)
    @test occursin("BURP-TI,BVK1", transfer_source)
    @test occursin("initial_density=initial.phi", transfer_source)
    @test occursin("write_density_field_visual_svg", transfer_source)
    @test occursin("burp_bvk1_gyr_polyorder_initialized_rows.csv",
        transfer_source)
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_POLYORDER_INITIALIZED_MAIN",
        transfer_source)

    reduced_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_reduced_shape_optimizer.jl")
    @test isfile(reduced_script_path)
    reduced_source = isfile(reduced_script_path) ?
        read(reduced_script_path, String) : ""
    @test occursin("Polyorder_results/gyroid/GYR_opt.mat", reduced_source)
    @test occursin("Optim.Brent", reduced_source)
    @test occursin("reduced_branch_stationarity_pass", reduced_source)
    @test occursin("BURP-TI,BVK1", reduced_source)
    @test occursin("write_density_field_visual_svg", reduced_source)
    @test occursin("burp_bvk1_gyr_reduced_shape_rows.csv", reduced_source)
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_REDUCED_SHAPE_MAIN",
        reduced_source)

    endpoint_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_endpoint_amplitude_scan.jl")
    @test isfile(endpoint_script_path)
    endpoint_source = isfile(endpoint_script_path) ?
        read(endpoint_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_ENDPOINT_SCAN_MAIN",
        endpoint_source)
    @test occursin("amplitude_values", endpoint_source)
    @test occursin("endpoint_drive", endpoint_source)
    @test occursin("BURP-TI,BVK1", endpoint_source)
    @test occursin("Polyorder_results/gyroid/GYR_opt.mat", endpoint_source)
    @test occursin("burp_bvk1_gyr_endpoint_amplitude_rows.csv",
        endpoint_source)
    @test occursin("burp_bvk1_gyr_endpoint_amplitude_summary.csv",
        endpoint_source)
    @test occursin("burp_bvk1_gyr_endpoint_amplitude_scan.svg",
        endpoint_source)
    @test occursin("burp_bvk1_gyr_endpoint_amplitude_scan.md",
        endpoint_source)
    @test occursin("_density_from_shape_logits", endpoint_source)
    @test occursin("_energy_density", endpoint_source)

    packet_source_after_endpoint = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("endpoint_amplitude_scan", packet_source_after_endpoint)
    @test occursin("SI Fig./Table 42", packet_source_after_endpoint)

    audit_source_after_endpoint = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_extension_audit.jl"), String)
    @test occursin("has_burp_bvk1_gyroid_endpoint_amplitude_scan",
        audit_source_after_endpoint)
    @test occursin("endpoint-driven reduced branch",
        audit_source_after_endpoint)

    decomposition_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_endpoint_energy_decomposition.jl")
    @test isfile(decomposition_script_path)
    decomposition_source = isfile(decomposition_script_path) ?
        read(decomposition_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_ENDPOINT_DECOMPOSITION_MAIN",
        decomposition_source)
    @test occursin("_burp_ti_energy_terms", decomposition_source)
    @test occursin("_bvk1_energy_terms", decomposition_source)
    @test occursin("energy_closure_error", decomposition_source)
    @test occursin("interaction_drop_from_a1", decomposition_source)
    @test occursin("burp_bvk1_gyr_endpoint_energy_decomposition_rows.csv",
        decomposition_source)
    @test occursin("burp_bvk1_gyr_endpoint_energy_decomposition_summary.csv",
        decomposition_source)
    @test occursin("burp_bvk1_gyr_endpoint_energy_decomposition.svg",
        decomposition_source)
    @test occursin("burp_bvk1_gyr_endpoint_energy_decomposition.md",
        decomposition_source)

    packet_source_after_decomposition = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("endpoint_energy_decomposition",
        packet_source_after_decomposition)
    @test occursin("SI Fig./Table 43", packet_source_after_decomposition)

    regularization_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_regularized_endpoint_sensitivity.jl")
    @test isfile(regularization_script_path)
    regularization_source = isfile(regularization_script_path) ?
        read(regularization_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_REGULARIZED_ENDPOINT_MAIN",
        regularization_source)
    @test occursin("_boundary_barrier_density", regularization_source)
    @test occursin("barrier_strength_values", regularization_source)
    @test occursin("regularized_energy_density", regularization_source)
    @test occursin("minimum_interior_barrier_strength", regularization_source)
    @test occursin("burp_bvk1_gyr_regularized_endpoint_sensitivity_rows.csv",
        regularization_source)
    @test occursin("burp_bvk1_gyr_regularized_endpoint_sensitivity_summary.csv",
        regularization_source)
    @test occursin("burp_bvk1_gyr_regularized_endpoint_sensitivity.svg",
        regularization_source)
    @test occursin("burp_bvk1_gyr_regularized_endpoint_sensitivity.md",
        regularization_source)

    packet_source_after_regularization = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("regularized_endpoint_sensitivity",
        packet_source_after_regularization)
    @test occursin("SI Fig./Table 44", packet_source_after_regularization)

    barrier_calibration_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_endpoint_barrier_calibration.jl")
    @test isfile(barrier_calibration_script_path)
    barrier_calibration_source = isfile(barrier_calibration_script_path) ?
        read(barrier_calibration_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_ENDPOINT_BARRIER_CALIBRATION_MAIN",
        barrier_calibration_source)
    @test occursin("_stationarity_lambda", barrier_calibration_source)
    @test occursin("dH_da", barrier_calibration_source)
    @test occursin("dB_da", barrier_calibration_source)
    @test occursin("polyorder_target_amplitude", barrier_calibration_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_calibration_rows.csv",
        barrier_calibration_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_calibration_summary.csv",
        barrier_calibration_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_calibration.svg",
        barrier_calibration_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_calibration.md",
        barrier_calibration_source)

    packet_source_after_barrier_calibration = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("endpoint_barrier_calibration",
        packet_source_after_barrier_calibration)
    @test occursin("SI Fig./Table 45", packet_source_after_barrier_calibration)

    barrier_stability_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_endpoint_barrier_stability.jl")
    @test isfile(barrier_stability_script_path)
    barrier_stability_source = isfile(barrier_stability_script_path) ?
        read(barrier_stability_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_ENDPOINT_BARRIER_STABILITY_MAIN",
        barrier_stability_source)
    @test occursin("_nonendpoint_identity_rows", barrier_stability_source)
    @test occursin("lambda_contrast_ratio", barrier_stability_source)
    @test occursin("lambda_model_mismatch_at_a1", barrier_stability_source)
    @test occursin("contrast_sensitivity_conclusion", barrier_stability_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_stability_rows.csv",
        barrier_stability_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_stability_summary.csv",
        barrier_stability_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_stability.svg",
        barrier_stability_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_stability.md",
        barrier_stability_source)

    packet_source_after_barrier_stability = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("endpoint_barrier_stability",
        packet_source_after_barrier_stability)
    @test occursin("SI Fig./Table 46", packet_source_after_barrier_stability)

    barrier_law_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_endpoint_barrier_law_fit.jl")
    @test isfile(barrier_law_script_path)
    barrier_law_source = isfile(barrier_law_script_path) ?
        read(barrier_law_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_ENDPOINT_BARRIER_LAW_MAIN",
        barrier_law_source)
    @test occursin("_fit_shared_barrier_density_quadratic", barrier_law_source)
    @test occursin("shared_law_rmse", barrier_law_source)
    @test occursin("shared_law_r2", barrier_law_source)
    @test occursin("model_offset_rmse", barrier_law_source)
    @test occursin("law_fit_conclusion", barrier_law_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_law_fit_rows.csv",
        barrier_law_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_law_fit_summary.csv",
        barrier_law_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_law_fit.svg",
        barrier_law_source)
    @test occursin("burp_bvk1_gyr_endpoint_barrier_law_fit.md",
        barrier_law_source)

    packet_source_after_barrier_law = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("endpoint_barrier_law_fit",
        packet_source_after_barrier_law)
    @test occursin("SI Fig./Table 47", packet_source_after_barrier_law)

    force_consistency_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_nd_force_consistency.jl")
    @test isfile(force_consistency_script_path)
    force_consistency_source = isfile(force_consistency_script_path) ?
        read(force_consistency_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_ND_FORCE_CONSISTENCY_MAIN",
        force_consistency_source)
    @test occursin("_burp_ti_chemical_potential_nd",
        force_consistency_source)
    @test occursin("_bvk1_frozen_k_chemical_potential_nd",
        force_consistency_source)
    @test occursin("adaptive_k_force_gap", force_consistency_source)
    @test occursin("central_finite_difference_directional_derivative",
        force_consistency_source)
    @test occursin("burp_bvk1_nd_force_consistency_rows.csv",
        force_consistency_source)
    @test occursin("burp_bvk1_nd_force_consistency_summary.csv",
        force_consistency_source)
    @test occursin("burp_bvk1_nd_force_consistency.svg",
        force_consistency_source)
    @test occursin("burp_bvk1_nd_force_consistency.md",
        force_consistency_source)

    packet_source_after_force_consistency = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("nd_force_consistency",
        packet_source_after_force_consistency)
    @test occursin("SI Fig./Table 48", packet_source_after_force_consistency)

    analytic_relaxation_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_ti_gyroid_analytic_force_relaxation.jl")
    @test isfile(analytic_relaxation_script_path)
    analytic_relaxation_source = isfile(analytic_relaxation_script_path) ?
        read(analytic_relaxation_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_TI_GYR_ANALYTIC_FORCE_MAIN",
        analytic_relaxation_source)
    @test occursin("_burp_ti_fixed_mean_projected_descent",
        analytic_relaxation_source)
    @test occursin("_burp_ti_chemical_potential_nd",
        analytic_relaxation_source)
    @test occursin("projected_chemical_potential_norm",
        analytic_relaxation_source)
    @test occursin("endpoint_descent", analytic_relaxation_source)
    @test occursin("burp_ti_gyr_analytic_force_relaxation_rows.csv",
        analytic_relaxation_source)
    @test occursin("burp_ti_gyr_analytic_force_relaxation_summary.csv",
        analytic_relaxation_source)
    @test occursin("burp_ti_gyr_analytic_force_relaxation.svg",
        analytic_relaxation_source)
    @test occursin("burp_ti_gyr_analytic_force_relaxation.md",
        analytic_relaxation_source)

    packet_source_after_analytic_relaxation = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("burp_ti_gyr_analytic_force_relaxation",
        packet_source_after_analytic_relaxation)
    @test occursin("SI Fig./Table 49",
        packet_source_after_analytic_relaxation)

    bvk1_adaptive_force_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_bvk1_adaptive_force_consistency.jl")
    @test isfile(bvk1_adaptive_force_script_path)
    bvk1_adaptive_force_source = isfile(bvk1_adaptive_force_script_path) ?
        read(bvk1_adaptive_force_script_path, String) : ""
    @test occursin("DFM_SKIP_BVK1_ADAPTIVE_FORCE_CONSISTENCY_MAIN",
        bvk1_adaptive_force_source)
    @test occursin("_bvk1_adaptive_k_chemical_potential_nd",
        bvk1_adaptive_force_source)
    @test occursin("_accumulate_bvk1_adaptive_k_derivative!",
        bvk1_adaptive_force_source)
    @test occursin("adaptive_k_derivative_closes_gap",
        bvk1_adaptive_force_source)
    @test occursin("bvk1_adaptive_force_consistency_rows.csv",
        bvk1_adaptive_force_source)
    @test occursin("bvk1_adaptive_force_consistency_summary.csv",
        bvk1_adaptive_force_source)
    @test occursin("bvk1_adaptive_force_consistency.svg",
        bvk1_adaptive_force_source)
    @test occursin("bvk1_adaptive_force_consistency.md",
        bvk1_adaptive_force_source)

    packet_source_after_bvk1_adaptive_force = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("bvk1_adaptive_force_consistency",
        packet_source_after_bvk1_adaptive_force)
    @test occursin("SI Fig./Table 50",
        packet_source_after_bvk1_adaptive_force)

    bvk1_force_regression = Module(:BVK1AdaptiveForceRegression)
    previous_force_main = get(ENV,
        "DFM_SKIP_BVK1_ADAPTIVE_FORCE_CONSISTENCY_MAIN", nothing)
    try
        ENV["DFM_SKIP_BVK1_ADAPTIVE_FORCE_CONSISTENCY_MAIN"] = "1"
        Base.include(bvk1_force_regression, joinpath(@__DIR__, "..", "scripts",
            "write_bvk1_adaptive_force_consistency.jl"))
    finally
        if previous_force_main === nothing
            delete!(ENV, "DFM_SKIP_BVK1_ADAPTIVE_FORCE_CONSISTENCY_MAIN")
        else
            ENV["DFM_SKIP_BVK1_ADAPTIVE_FORCE_CONSISTENCY_MAIN"] =
                previous_force_main
        end
    end
    small_phi = DFMMonteCarlo._diblock_density_from_seed_field(
        DFMMonteCarlo._diblock_morphology_seed_field(:GYR, (6, 6, 6);
            amplitude=0.35), 0.37)
    small_rows = bvk1_force_regression._rows_for_force_check(small_phi,
        (3.2, 3.2, 3.2); f=0.37, chiN=16.0, direction_count=3,
        finite_difference_step=1.0e-5, seed=20260713)
    adaptive_rows = filter(row ->
        row.force_gate == "bvk1_full_adaptive_force_vs_adaptive_energy",
        small_rows)
    @test maximum(row.relative_error for row in adaptive_rows) < 1.0e-5

    # The BVK1 coefficient contains |Delta eta| and is therefore nonsmooth on
    # symmetry-enforced Laplacian-zero manifolds.  The continuation solver uses
    # a differentiable, uniformly convergent approximation only while relaxing
    # the field; the default zero smoothing must remain exactly the production
    # BVK1 functional.
    legacy_energy = DFMMonteCarlo.diblock_bvk1_energy_nd(small_phi;
        f=0.37, chiN=16.0, lengths=(3.2, 3.2, 3.2), adaptive=true)
    zero_smoothing_energy = DFMMonteCarlo.diblock_bvk1_energy_nd(small_phi;
        f=0.37, chiN=16.0, lengths=(3.2, 3.2, 3.2), adaptive=true,
        laplacian_smoothing=0.0)
    @test zero_smoothing_energy == legacy_energy

    regularized_rows = bvk1_force_regression._rows_for_force_check(small_phi,
        (3.2, 3.2, 3.2); f=0.37, chiN=16.0, direction_count=3,
        finite_difference_step=1.0e-5, seed=20260713,
        laplacian_smoothing=1.0e-3)
    regularized_adaptive_rows = filter(row ->
        row.force_gate == "bvk1_full_adaptive_force_vs_adaptive_energy",
        regularized_rows)
    @test maximum(row.relative_error for row in regularized_adaptive_rows) <
        1.0e-5

    coarse_regularized_energy = DFMMonteCarlo.diblock_bvk1_energy_nd(small_phi;
        f=0.37, chiN=16.0, lengths=(3.2, 3.2, 3.2), adaptive=true,
        laplacian_smoothing=1.0e-2)
    fine_regularized_energy = DFMMonteCarlo.diblock_bvk1_energy_nd(small_phi;
        f=0.37, chiN=16.0, lengths=(3.2, 3.2, 3.2), adaptive=true,
        laplacian_smoothing=1.0e-4)
    @test abs(fine_regularized_energy - legacy_energy) <
        abs(coarse_regularized_energy - legacy_energy)

    allocation_phi = DFMMonteCarlo._diblock_density_from_seed_field(
        DFMMonteCarlo._diblock_morphology_seed_field(:GYR, (24, 24, 24);
            amplitude=0.35), 0.37)
    bvk1_force_regression._bvk1_adaptive_k_chemical_potential_nd(
        allocation_phi; f=0.37, chiN=16.0, lengths=(3.2, 3.2, 3.2))
    GC.gc()
    force_allocations = @allocated begin
        bvk1_force_regression._bvk1_adaptive_k_chemical_potential_nd(
            allocation_phi; f=0.37, chiN=16.0,
            lengths=(3.2, 3.2, 3.2))
    end
    @test force_allocations < 60_000_000

    gyroid_period_audit_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_gyroid_period_convention_audit.jl")
    @test isfile(gyroid_period_audit_script_path)
    gyroid_period_audit_source = isfile(gyroid_period_audit_script_path) ?
        read(gyroid_period_audit_script_path, String) : ""
    @test occursin("DFM_SKIP_GYROID_PERIOD_CONVENTION_AUDIT_MAIN",
        gyroid_period_audit_source)
    @test occursin("sqrt(6.0)", gyroid_period_audit_source)
    @test occursin("legacy_110", gyroid_period_audit_source)
    @test occursin("polyorder_211", gyroid_period_audit_source)
    @test occursin("same_unit_ratio", gyroid_period_audit_source)
    @test occursin("do_not_compare_b_to_rg", gyroid_period_audit_source)
    @test occursin("gyroid_period_convention_audit.csv",
        gyroid_period_audit_source)
    @test occursin("gyroid_period_convention_audit.svg",
        gyroid_period_audit_source)
    @test occursin("gyroid_period_convention_audit.md",
        gyroid_period_audit_source)

    packet_source_after_gyroid_period_audit = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("gyroid_period_convention_audit",
        packet_source_after_gyroid_period_audit)
    @test occursin("SI Fig./Table 51",
        packet_source_after_gyroid_period_audit)

    bvk1_adaptive_relaxation_script_path = joinpath(@__DIR__, "..",
        "scripts", "write_bvk1_gyroid_adaptive_force_relaxation.jl")
    @test isfile(bvk1_adaptive_relaxation_script_path)
    bvk1_adaptive_relaxation_source =
        isfile(bvk1_adaptive_relaxation_script_path) ?
        read(bvk1_adaptive_relaxation_script_path, String) : ""
    @test occursin("DFM_SKIP_BVK1_GYR_ADAPTIVE_FORCE_RELAXATION_MAIN",
        bvk1_adaptive_relaxation_source)
    @test occursin("_bvk1_fixed_mean_projected_descent",
        bvk1_adaptive_relaxation_source)
    @test occursin("_bvk1_projected_chemical_potential",
        bvk1_adaptive_relaxation_source)
    @test occursin("_bvk1_adaptive_k_chemical_potential_nd",
        bvk1_adaptive_relaxation_source)
    @test occursin("verified_adaptive_bvk1_nd_force",
        bvk1_adaptive_relaxation_source)
    @test occursin("endpoint_descent", bvk1_adaptive_relaxation_source)
    @test occursin("bvk1_gyr_adaptive_force_relaxation_rows.csv",
        bvk1_adaptive_relaxation_source)
    @test occursin("bvk1_gyr_adaptive_force_relaxation_summary.csv",
        bvk1_adaptive_relaxation_source)
    @test occursin("bvk1_gyr_adaptive_force_relaxation.svg",
        bvk1_adaptive_relaxation_source)
    @test occursin("bvk1_gyr_adaptive_force_relaxation.md",
        bvk1_adaptive_relaxation_source)

    packet_source_after_bvk1_adaptive_relaxation = read(joinpath(@__DIR__,
        "..", "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("bvk1_gyr_adaptive_force_relaxation",
        packet_source_after_bvk1_adaptive_relaxation)
    @test occursin("SI Fig./Table 52",
        packet_source_after_bvk1_adaptive_relaxation)

    endpoint_policy_doc_path = joinpath(@__DIR__, "..", "docs",
        "nonlamellar_endpoint_admissibility_policy.md")
    @test isfile(endpoint_policy_doc_path)
    endpoint_policy_text = isfile(endpoint_policy_doc_path) ?
        read(endpoint_policy_doc_path, String) : ""
    @test occursin("endpoint distance", endpoint_policy_text)
    @test occursin("d_endpoint", endpoint_policy_text)
    @test occursin("Polyorder-logit amplitude", endpoint_policy_text)
    @test occursin("psi_i(r)=sqrt(phi_i(r))", endpoint_policy_text)
    @test occursin("saturation is physically natural",
        endpoint_policy_text)
    @test occursin("saturation_is_not_automatic_failure",
        endpoint_policy_text)
    @test occursin("grid/bound stability",
        endpoint_policy_text)
    @test occursin("branch-local barrier law", endpoint_policy_text)
    @test occursin("not a transferable closure", endpoint_policy_text)
    @test occursin("independent SCFT/UD reference profiles",
        endpoint_policy_text)
    @test occursin("admissibility gate", endpoint_policy_text)

    packet_source_after_endpoint_policy = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("nonlamellar_endpoint_admissibility_policy",
        packet_source_after_endpoint_policy)
    @test occursin("SI Note 53", packet_source_after_endpoint_policy)

    audit_source_after_endpoint_policy = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_morphology_extension_audit.jl"), String)
    @test occursin("has_nonlamellar_endpoint_admissibility_policy",
        audit_source_after_endpoint_policy)
    @test occursin("endpoint admissibility policy",
        audit_source_after_endpoint_policy)

    endpoint_chin_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_endpoint_chin_scan.jl")
    @test isfile(endpoint_chin_script_path)
    endpoint_chin_source = isfile(endpoint_chin_script_path) ?
        read(endpoint_chin_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_ENDPOINT_CHIN_SCAN_MAIN",
        endpoint_chin_source)
    @test occursin("chiN_values", endpoint_chin_source)
    @test occursin("endpoint_drive_threshold_chiN", endpoint_chin_source)
    @test occursin("lower_chiN_interior_window", endpoint_chin_source)
    @test occursin("high_chiN_endpoint_drive", endpoint_chin_source)
    @test occursin("fixed Polyorder-GYR branch", endpoint_chin_source)
    @test occursin("burp_bvk1_gyr_endpoint_chin_rows.csv",
        endpoint_chin_source)
    @test occursin("burp_bvk1_gyr_endpoint_chin_summary.csv",
        endpoint_chin_source)
    @test occursin("burp_bvk1_gyr_endpoint_chin_scan.svg",
        endpoint_chin_source)
    @test occursin("burp_bvk1_gyr_endpoint_chin_scan.md",
        endpoint_chin_source)

    packet_source_after_endpoint_chin = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("endpoint_chin_scan",
        packet_source_after_endpoint_chin)
    @test occursin("SI Fig./Table 54",
        packet_source_after_endpoint_chin)

    force_chin_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_force_chin_scan.jl")
    @test isfile(force_chin_script_path)
    force_chin_source = isfile(force_chin_script_path) ?
        read(force_chin_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_FORCE_CHIN_SCAN_MAIN",
        force_chin_source)
    @test occursin("chiN_values", force_chin_source)
    @test occursin("full-field projected descent", force_chin_source)
    @test occursin("stationary_tolerance", force_chin_source)
    @test occursin("endpoint_floor", force_chin_source)
    @test occursin("interior_stalled", force_chin_source)
    @test occursin("force_chin_result", force_chin_source)
    @test occursin("burp_bvk1_gyr_force_chin_rows.csv",
        force_chin_source)
    @test occursin("burp_bvk1_gyr_force_chin_summary.csv",
        force_chin_source)
    @test occursin("burp_bvk1_gyr_force_chin_scan.svg",
        force_chin_source)
    @test occursin("burp_bvk1_gyr_force_chin_scan.md",
        force_chin_source)

    packet_source_after_force_chin = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("force_chin_scan", packet_source_after_force_chin)
    @test occursin("SI Fig./Table 55", packet_source_after_force_chin)
    @test occursin("force_chin_long_scan", packet_source_after_force_chin)
    @test occursin("SI Fig./Table 56", packet_source_after_force_chin)

    audit_source_after_force_chin = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_morphology_extension_audit.jl"), String)
    @test occursin("has_burp_bvk1_gyroid_force_chin_scan",
        audit_source_after_force_chin)
    @test occursin("has_burp_bvk1_gyroid_force_chin_long_scan",
        audit_source_after_force_chin)
    @test occursin("full-field force chiN scan", audit_source_after_force_chin)
    @test occursin("lower-chiN long convergence", audit_source_after_force_chin)

    endpoint_kkt_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyroid_endpoint_kkt_diagnostic.jl")
    @test isfile(endpoint_kkt_script_path)
    endpoint_kkt_source = isfile(endpoint_kkt_script_path) ?
        read(endpoint_kkt_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_ENDPOINT_KKT_MAIN",
        endpoint_kkt_source)
    @test occursin("fixed-mean bound-constrained KKT", endpoint_kkt_source)
    @test occursin("active_lower_count", endpoint_kkt_source)
    @test occursin("active_upper_count", endpoint_kkt_source)
    @test occursin("free_projected_residual_norm", endpoint_kkt_source)
    @test occursin("kkt_violation_max", endpoint_kkt_source)
    @test occursin("kkt_candidate", endpoint_kkt_source)
    @test occursin("burp_bvk1_gyr_endpoint_kkt_rows.csv",
        endpoint_kkt_source)
    @test occursin("burp_bvk1_gyr_endpoint_kkt_summary.csv",
        endpoint_kkt_source)
    @test occursin("burp_bvk1_gyr_endpoint_kkt.svg", endpoint_kkt_source)
    @test occursin("burp_bvk1_gyr_endpoint_kkt.md", endpoint_kkt_source)

    packet_source_after_endpoint_kkt = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("endpoint_kkt_diagnostic",
        packet_source_after_endpoint_kkt)
    @test occursin("SI Fig./Table 57", packet_source_after_endpoint_kkt)

    audit_source_after_endpoint_kkt = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_morphology_extension_audit.jl"), String)
    @test occursin("has_burp_bvk1_gyroid_endpoint_kkt_diagnostic",
        audit_source_after_endpoint_kkt)
    @test occursin("endpoint KKT diagnostic", audit_source_after_endpoint_kkt)

    saturation_audit_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_strong_segregation_saturation_audit.jl")
    @test isfile(saturation_audit_script_path)
    saturation_audit_source = isfile(saturation_audit_script_path) ?
        read(saturation_audit_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_SATURATION_AUDIT_MAIN",
        saturation_audit_source)
    @test occursin("physical_strong_segregation_window",
        saturation_audit_source)
    @test occursin("grid_bound_sensitivity_required",
        saturation_audit_source)
    @test occursin("stationarity_required",
        saturation_audit_source)
    @test occursin("saturation_is_not_automatic_failure",
        saturation_audit_source)
    @test occursin("polyorder_reference_chiN", saturation_audit_source)
    @test occursin("burp_bvk1_strong_segregation_saturation_rows.csv",
        saturation_audit_source)
    @test occursin("burp_bvk1_strong_segregation_saturation_summary.csv",
        saturation_audit_source)
    @test occursin("burp_bvk1_strong_segregation_saturation.svg",
        saturation_audit_source)
    @test occursin("burp_bvk1_strong_segregation_saturation.md",
        saturation_audit_source)

    packet_source_after_saturation_audit = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("strong_segregation_saturation_audit",
        packet_source_after_saturation_audit)
    @test occursin("SI Fig./Table 58", packet_source_after_saturation_audit)

    audit_source_after_saturation = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_morphology_extension_audit.jl"), String)
    @test occursin("has_burp_bvk1_strong_segregation_saturation_audit",
        audit_source_after_saturation)
    @test occursin("saturation is not automatic failure",
        audit_source_after_saturation)

    bound_sensitivity_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_gyr_bound_sensitivity.jl")
    @test isfile(bound_sensitivity_script_path)
    bound_sensitivity_source = isfile(bound_sensitivity_script_path) ?
        read(bound_sensitivity_script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_GYR_BOUND_SENSITIVITY_MAIN",
        bound_sensitivity_source)
    @test occursin("bound_stability_pass", bound_sensitivity_source)
    @test occursin("endpoint_floor_values", bound_sensitivity_source)
    @test occursin("stationarity_or_kkt_pass", bound_sensitivity_source)
    @test occursin("grid_bound_sensitivity_required",
        bound_sensitivity_source)
    @test occursin("burp_bvk1_gyr_bound_sensitivity_rows.csv",
        bound_sensitivity_source)
    @test occursin("burp_bvk1_gyr_bound_sensitivity_summary.csv",
        bound_sensitivity_source)
    @test occursin("burp_bvk1_gyr_bound_sensitivity.svg",
        bound_sensitivity_source)
    @test occursin("burp_bvk1_gyr_bound_sensitivity.md",
        bound_sensitivity_source)

    packet_source_after_bound_sensitivity = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("gyr_bound_sensitivity",
        packet_source_after_bound_sensitivity)
    @test occursin("SI Fig./Table 59",
        packet_source_after_bound_sensitivity)

    audit_source_after_bound_sensitivity = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_morphology_extension_audit.jl"), String)
    @test occursin("has_burp_bvk1_gyr_bound_sensitivity",
        audit_source_after_bound_sensitivity)
    @test occursin("bound sensitivity",
        audit_source_after_bound_sensitivity)

    variational_note_path = joinpath(@__DIR__, "..", "docs",
        "strong_segregation_variational_principle.md")
    @test isfile(variational_note_path)
    variational_note_text = isfile(variational_note_path) ?
        read(variational_note_path, String) : ""
    @test occursin("fixed-composition variational inequality",
        variational_note_text)
    @test occursin("Gamma-limit interpretation", variational_note_text)
    @test occursin("floor-stability criterion", variational_note_text)
    @test occursin("KKT sign convention", variational_note_text)
    @test occursin("Uneyama-Doi square-root variable",
        variational_note_text)
    @test occursin("not evidence against saturation",
        variational_note_text)

    packet_source_after_variational_note = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("strong_segregation_variational_principle",
        packet_source_after_variational_note)
    @test occursin("SI Note 60", packet_source_after_variational_note)

    audit_source_after_variational_note = read(joinpath(@__DIR__, "..",
        "scripts", "write_burp_bvk1_morphology_extension_audit.jl"),
        String)
    @test occursin("has_strong_segregation_variational_principle",
        audit_source_after_variational_note)
    @test occursin("variational principle",
        audit_source_after_variational_note)
end

@testset "BURP-BVK1 publication readiness matrix source contract" begin
    readiness_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl")
    @test isfile(readiness_script_path)
    readiness_source = isfile(readiness_script_path) ?
        read(readiness_script_path, String) : ""
    @test occursin("manuscript_readiness_matrix.csv", readiness_source)
    @test occursin("manuscript_readiness_matrix.md", readiness_source)
    @test occursin("publish_now", readiness_source)
    @test occursin("supplement_guardrail", readiness_source)
    @test occursin("not_yet", readiness_source)
    @test occursin("claim_boundary", readiness_source)
    @test occursin("cyl_logit_cell_scan_summary.csv", readiness_source)
    @test occursin("nonlamellar_logit_visual_gallery", readiness_source)
    @test occursin("visual morphology audit", readiness_source)
    @test occursin("DIS/LAM/CYL/BCC/GYR", readiness_source)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate",
        readiness_source)
    @test occursin("limited same-target GYR phase-rank passes",
        readiness_source)
    @test occursin("fA=0.39, chiN=20", readiness_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("manuscript_readiness_matrix", packet_source)
    @test occursin("SI Table 61", packet_source)

    audit_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_extension_audit.jl"), String)
    @test occursin("has_burp_bvk1_publication_readiness_matrix",
        audit_source)
    @test occursin("publication readiness", audit_source)
end

@testset "Uneyama-Doi reference coverage source contract" begin
    coverage_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl")
    @test isfile(coverage_script_path)
    coverage_source = isfile(coverage_script_path) ?
        read(coverage_script_path, String) : ""
    @test occursin("ud_reference_coverage.csv", coverage_source)
    @test occursin("ud_reference_coverage.md", coverage_source)
    @test occursin("mean_square_distance.eps", coverage_source)
    @test occursin("S_diblock.eps", coverage_source)
    @test occursin("S_triblock.eps", coverage_source)
    @test occursin("AB_phasediagram.eps", coverage_source)
    @test occursin("AB_period.eps", coverage_source)
    @test occursin("AB_melt.eps", coverage_source)
    @test occursin("A_B_AB_blend.eps", coverage_source)
    @test occursin("AB_C_blend1.eps", coverage_source)
    @test occursin("AB_C_blend2.eps", coverage_source)
    @test occursin("covered_core_ab_diblock", coverage_source)
    @test occursin("partial_nonlamellar", coverage_source)
    @test occursin("out_of_scope_current_model", coverage_source)
    @test occursin("cyl_logit_cell_scan_summary.csv", coverage_source)
    @test occursin("nonlamellar_logit_visual_gallery", coverage_source)
    @test occursin("visual morphology audit", coverage_source)
    @test occursin("DIS/LAM/CYL/BCC/GYR", coverage_source)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate",
        coverage_source)
    @test occursin("limited same-target GYR phase-rank passes",
        coverage_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("ud_reference_coverage", packet_source)
    @test occursin("SI Table 62", packet_source)

    audit_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_extension_audit.jl"), String)
    @test occursin("has_burp_bvk1_ud_reference_coverage", audit_source)
    @test occursin("UD reference coverage", audit_source)
end

@testset "UD CYL/BCC method-equity gate source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_cyl_bcc_method_equity_gate.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_UD_CYL_BCC_METHOD_EQUITY_MAIN", source)
    @test occursin("ud_cyl_bcc_method_equity_gate.csv", source)
    @test occursin("ud_cyl_bcc_method_equity_gate_summary.csv", source)
    @test occursin("ud_cyl_bcc_method_equity_gate.svg", source)
    @test occursin("ud_cyl_bcc_method_equity_gate.md", source)
    @test occursin("f0.15_chiN50_BCC", source)
    @test occursin("f0.30_chiN50_CYL", source)
    @test occursin("ud_representative_phase_rank_gate.csv", source)
    @test occursin("ud_endpoint_refinement_gate.csv", source)
    @test occursin("UD-first CYL/BCC method-equity gate", source)
    @test occursin("representative_cyl_bcc_pass_count", source)
    @test occursin("endpoint_caveat_count", source)

    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("ud_cyl_bcc_method_equity_gate", pipeline_source)
    @test occursin("write_uneyama_doi_cyl_bcc_method_equity_gate.jl",
        pipeline_source)

    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("ud_cyl_bcc_method_equity_gate", coverage_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("UD CYL/BCC method-equity gate", packet_source)
    @test occursin("ud_cyl_bcc_method_equity_gate", packet_source)
end

@testset "BURP-BVK1 UD parity scorecard source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_parity_scorecard.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_UD_PARITY_SCORECARD_MAIN",
        source)
    @test occursin("ud_reference_coverage.csv", source)
    @test occursin("covered_core_ab_diblock", source)
    @test occursin("partial_nonlamellar", source)
    @test occursin("out_of_scope_current_model", source)
    @test occursin("equivalent_ud_claim_ready", source)
    @test occursin("ud_parity_scorecard.csv", source)
    @test occursin("ud_parity_scorecard_summary.csv", source)
    @test occursin("ud_parity_scorecard.svg", source)
    @test occursin("ud_parity_scorecard.md", source)
    @test occursin("AB diblock", source)
    @test occursin("phase diagram", source)
    @test occursin("blend/triblock", source)

    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("ud_parity_scorecard", pipeline_source)
    @test occursin("write_burp_bvk1_ud_parity_scorecard.jl",
        pipeline_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("Uneyama-Doi parity scorecard", packet_source)
    @test occursin("ud_parity_scorecard", packet_source)

    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("ud_parity_scorecard", readiness_source)
end

@testset "Polyorder nonlamellar reference provenance source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_reference_provenance.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_REFERENCE_PROVENANCE_MAIN",
        source)
    @test occursin("Polyorder_results/CYL/CYL_opt.mat", source)
    @test occursin("Polyorder_results/BCC/BCC_opt.mat", source)
    @test occursin("Polyorder_results/gyroid/GYR_opt.mat", source)
    @test occursin("Polyorder_results/CYL_f0.36_chiN16/CYL_opt.mat",
        source)
    @test occursin("same_parameter_cyl_for_gyr_rank", source)
    @test occursin("polyorder_reference_provenance.csv", source)
    @test occursin("polyorder_reference_provenance.svg", source)
    @test occursin("polyorder_reference_provenance.md", source)

    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("polyorder_reference_provenance", pipeline_source)
    @test occursin("write_burp_bvk1_polyorder_reference_provenance.jl",
        pipeline_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("Polyorder SCFT reference provenance", packet_source)
    @test occursin("polyorder_reference_provenance", packet_source)

    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("polyorder_reference_provenance", readiness_source)
end

@testset "Polyorder nonlamellar branch landscape source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_branch_landscape.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_BRANCH_LANDSCAPE_MAIN",
        source)
    @test occursin("nonlamellar_polyorder_scale_scan.csv", source)
    @test occursin("initial_curvature_proxy", source)
    @test occursin("finite_amplitude_well_depth", source)
    @test occursin("landscape_class", source)
    @test occursin("dis_stable_branch", source)
    @test occursin("finite_amplitude_well", source)
    @test occursin("endpoint_or_scan_boundary_drive", source)
    @test occursin("nonlamellar_branch_landscape.csv", source)
    @test occursin("nonlamellar_branch_landscape.svg", source)
    @test occursin("nonlamellar_branch_landscape.md", source)

    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("nonlamellar_branch_landscape", pipeline_source)
    @test occursin(
        "write_burp_bvk1_polyorder_nonlamellar_branch_landscape.jl",
        pipeline_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("Nonlamellar SCFT-branch amplitude landscape",
        packet_source)
    @test occursin("nonlamellar_branch_landscape", packet_source)

    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("nonlamellar_branch_landscape", readiness_source)
end

@testset "Polyorder nonlamellar shape fidelity source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_shape_fidelity.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_SHAPE_FIDELITY_MAIN",
        source)
    @test occursin("nonlamellar_logit_visual_summary.csv", source)
    @test occursin("density_rms", source)
    @test occursin("centered_correlation", source)
    @test occursin("contrast_ratio", source)
    @test occursin("shape_fidelity_interpretation", source)
    @test occursin("nonlamellar_shape_fidelity.csv", source)
    @test occursin("nonlamellar_shape_fidelity.svg", source)
    @test occursin("nonlamellar_shape_fidelity.md", source)

    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("nonlamellar_shape_fidelity", pipeline_source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_shape_fidelity.jl",
        pipeline_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("Nonlamellar SCFT shape fidelity", packet_source)
    @test occursin("nonlamellar_shape_fidelity", packet_source)

    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("nonlamellar_shape_fidelity", readiness_source)
end

@testset "Polyorder nonlamellar BVK1 adaptive-stiffness source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_adaptive_stiffness.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_ADAPTIVE_STIFFNESS_MAIN",
        source)
    @test occursin("nonlamellar_logit_visual_summary.csv", source)
    @test occursin("DFMMonteCarlo._diblock_bvk1_kpsi_nd", source)
    @test occursin("adaptive_stiffness_mean", source)
    @test occursin("adaptive_stiffness_p95", source)
    @test occursin("adaptive_stiffness_max", source)
    @test occursin("interface_weighted_mean", source)
    @test occursin("nonlamellar_adaptive_stiffness.csv", source)
    @test occursin("nonlamellar_adaptive_stiffness.svg", source)
    @test occursin("nonlamellar_adaptive_stiffness.md", source)

    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("nonlamellar_adaptive_stiffness", pipeline_source)
    @test occursin(
        "write_burp_bvk1_polyorder_nonlamellar_adaptive_stiffness.jl",
        pipeline_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("Nonlamellar BVK1 adaptive stiffness", packet_source)
    @test occursin("nonlamellar_adaptive_stiffness", packet_source)

    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("nonlamellar_adaptive_stiffness", readiness_source)
end

@testset "Polyorder nonlamellar cell-scale stationarity source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_cell_scale_stationarity.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_CELL_MAIN",
        source)
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_REFERENCE_MAIN",
        source)
    @test occursin("_generate_cyl_reference", source)
    @test occursin("Polyorder_results", source)
    @test occursin("BCC_opt.mat", source)
    @test occursin("GYR_opt.mat", source)
    @test occursin("_scale_density_logit", source)
    @test occursin("amplitude_scales", source)
    @test occursin("cell_factors", source)
    @test occursin("local_2d_minimum_check_pass", source)
    @test occursin("homogeneous_branch_degenerate_cell_scale", source)
    @test occursin("nonlamellar_cell_scale_references.csv", source)
    @test occursin("nonlamellar_cell_scale_grid.csv", source)
    @test occursin("nonlamellar_cell_scale_summary.csv", source)
    @test occursin("nonlamellar_cell_scale_stationarity.svg", source)
    @test occursin("nonlamellar_cell_scale_stationarity.md", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_nonlamellar_cell_scale_stationarity.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("branch-local cell-scale", doc_text)
    @test occursin("CYL/BURP-TI", doc_text)
    @test occursin("GYR | BVK1", doc_text)
    @test occursin("BCC", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("Polyorder nonlamellar branch cell-scale stationarity",
        packet_source)
    @test occursin("nonlamellar_cell_scale_stationarity",
        packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("nonlamellar_cell_scale_summary.csv", readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("polyorder_scft_nonlamellar_cell_scale_stationarity",
        coverage_source)
end

@testset "Polyorder nonlamellar unrestricted relaxation source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_unrestricted_relaxation.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_UNRESTRICTED_MAIN",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_reference_comparison.jl",
        source)
    @test occursin("write_burp_bvk1_nd_force_consistency.jl", source)
    @test occursin("write_bvk1_adaptive_force_consistency.jl", source)
    @test occursin("scft_point", source)
    @test occursin("branch_selected", source)
    @test occursin("fixed_mean_projected_descent", source)
    # The dispatch harness carries a BVK2-specific force tolerance; the
    # contract pins the exact dispatch expression so a silent revert to the
    # single-model tolerance is caught.
    @test occursin(
        "force_tolerance::Real=(model == \"BVK2\" ? 1.0e-5 : tolerance)",
        source)
    @test occursin("force_norm <= force_tol", source)
    @test occursin("nonlamellar_unrestricted_relaxation_summary.csv", source)
    @test occursin("nonlamellar_unrestricted_relaxation_paths.csv", source)
    @test occursin("nonlamellar_unrestricted_relaxation.svg", source)
    @test occursin("nonlamellar_unrestricted_relaxation.md", source)
    @test occursin("identity_preserved", source)
    @test occursin("endpoint_limited", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_nonlamellar_unrestricted_relaxation.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("unrestricted fixed-composition", doc_text)
    @test occursin("Polyorder SCFT", doc_text)
    @test occursin("BURP-TI", doc_text)
    @test occursin("BVK1", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("nonlamellar_unrestricted_relaxation",
        packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("nonlamellar_unrestricted_relaxation_summary.csv",
        readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("polyorder_scft_nonlamellar_unrestricted_relaxation",
        coverage_source)
end

@testset "BVK1 production branch solver source contract" begin
    source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_branch_solver_3d.jl"), String)
    @test occursin("_branch3d_periodic_linear_resample", source)
    @test occursin("_branch3d_prepare_warm_field", source)
    @test occursin("_branch3d_solve_bvk1_scalar", source)
    @test occursin("bvk1_alternating_projected_field_frozen_cell_bracket",
        source)
    @test occursin("warm_cell_window_fraction", source)
    @test occursin("bounded_kkt_maxabs", source)
    @test occursin("field_stationarity", source)
    @test occursin("cell_accepted", source)
    @test occursin("\"BCC\" => (32, 32, 32)", source)
    @test occursin("\"CYL\" => (40, 72)", source)
    @test occursin("\"GYR\" => (48, 48, 48)", source)
end

@testset "BVK1 fixed-fA S/DIS campaign source contract" begin
    point_driver = joinpath(@__DIR__, "..", "scripts",
        "run_bvk1_sdis_fixed_fa_point.jl")
    launcher = joinpath(@__DIR__, "..", "scripts",
        "launch_bvk1_sdis_fixed_fa_campaign.py")
    assembler = joinpath(@__DIR__, "..", "scripts",
        "assemble_bvk1_publication_phase_diagram.py")
    @test isfile(point_driver)
    @test isfile(launcher)
    @test isfile(assembler)
    point_source = read(point_driver, String)
    @test occursin("downward_warm_continuation", point_source)
    @test occursin("production_resolution", point_source)
    @test occursin("resolution_consistent", point_source)
    @test occursin("gates_lower", point_source)
    @test occursin("gates_upper", point_source)
    # Pins refreshed to the reworked fixed-fA driver: collapse detection is
    # now the "phase_collapsed" status, endpoint acceptance chains the strict
    # bracket with the lower-endpoint gates, and production dims are 36^3.
    @test occursin("phase_collapsed", point_source)
    @test occursin(
        "gates_pass = strict_bracketed && result.low_eval.gates",
        point_source)
    @test occursin("bfs_arg(\"--dims\", \"36x36x36\")", point_source)
    launcher_source = read(launcher, String)
    @test occursin("sdis_fixed_fa_roots.csv", launcher_source)
    @test occursin("default=0.09", launcher_source)
    @test occursin("bfs_arg(\"--max-grid-spacing\", \"0.09\")", point_source)
    @test occursin("accepted_fixed_fa_sdis", read(assembler, String))
    python_test = joinpath(@__DIR__,
        "test_assemble_bvk1_publication_phase_diagram.py")
    run(`python3 $python_test`)
end

@testset "BVK1 resolution-gated pairwise locator source contract" begin
    driver = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_bb_boundary_bisection.jl")
    launcher = joinpath(@__DIR__, "..", "scripts",
        "launch_bvk1_bb_resolution_campaign.py")
    assembler = joinpath(@__DIR__, "..", "scripts",
        "assemble_bvk1_publication_phase_diagram.py")
    @test isfile(driver)
    @test isfile(launcher)
    driver_source = read(driver, String)
    @test occursin("max_grid_spacing", driver_source)
    @test occursin("coarse_confirmation", driver_source)
    @test occursin("resolution_consistent", driver_source)
    @test occursin("production_resolution", driver_source)
    @test occursin("valid_a && valid_b", driver_source)
    launcher_source = read(launcher, String)
    @test occursin("--coarse-roots=", launcher_source)
    @test occursin("default=\"40x40x40\"", launcher_source)
    @test occursin("default=\"48x48x48\"", launcher_source)
    @test occursin("default=\"28x28x28\"", launcher_source)
    @test occursin("default=\"32x32x32\"", launcher_source)
    @test occursin("f\"--gyr-dims={args.locator_gyr_dims}\"",
        launcher_source)
    @test occursin("f\"--gyr-dims={args.production_gyr_dims}\"",
        launcher_source)
    @test occursin("f\"--bcc-dims={args.locator_bcc_dims}\"",
        launcher_source)
    @test occursin("f\"--bcc-dims={args.production_bcc_dims}\"",
        launcher_source)
    @test occursin("--coarse-max-grid-spacing", launcher_source)
    @test occursin("--production-max-grid-spacing", launcher_source)
    @test occursin("PRODUCTION_MAX_GRID_SPACING = 0.09", launcher_source)
    @test occursin("default=PRODUCTION_MAX_GRID_SPACING", launcher_source)
    @test occursin("accepted_pair_roots", read(assembler, String))
end

@testset "BVK1 branch stationarity translation gauge" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "polish_bvk1_branch_stationarity.jl")
    mod = Module(:TestBVK1NewtonGMRESTranslationGauge)
    previous = get(ENV, "DFM_SKIP_BVK1_BRANCH_STATIONARITY_MAIN", nothing)
    try
        ENV["DFM_SKIP_BVK1_BRANCH_STATIONARITY_MAIN"] = "1"
        Base.include(mod, script_path)
    finally
        if previous === nothing
            delete!(ENV, "DFM_SKIP_BVK1_BRANCH_STATIONARITY_MAIN")
        else
            ENV["DFM_SKIP_BVK1_BRANCH_STATIONARITY_MAIN"] = previous
        end
    end

    stable_trace = [
        (iteration=0, kkt_norm=9.0e-3, kkt_maxabs=0.3,
            energy_density=-1.0),
        (iteration=50, kkt_norm=8.0e-3, kkt_maxabs=0.4,
            energy_density=-1.0-4.0e-9),
        (iteration=100, kkt_norm=7.0e-3, kkt_maxabs=0.5,
            energy_density=-1.0-7.0e-9),
    ]
    plateau = mod.bng_r2_energy_plateau(stable_trace;
        r2_tolerance=1.0e-2, energy_tolerance=1.0e-8, window=3)
    @test plateau.pass
    @test plateau.available
    @test plateau.energy_span <= 1.0e-8
    @test plateau.max_abs_delta <= 1.0e-8
    @test !mod.bng_r2_energy_plateau(stable_trace[1:2];
        r2_tolerance=1.0e-2, energy_tolerance=1.0e-8, window=3).pass
    @test !mod.bng_r2_energy_plateau([
            stable_trace[1], merge(stable_trace[2], (kkt_norm=1.1e-2,)),
            stable_trace[3]];
        r2_tolerance=1.0e-2, energy_tolerance=1.0e-8, window=3).pass
    @test !mod.bng_r2_energy_plateau([
            stable_trace[1], merge(stable_trace[2],
                (energy_density=-1.0-2.0e-8,)), stable_trace[3]];
        r2_tolerance=1.0e-2, energy_tolerance=1.0e-8, window=3).pass
    duplicate_trace = [stable_trace[1], stable_trace[2],
        merge(stable_trace[2], (energy_density=-1.0-5.0e-9,))]
    @test !mod.bng_r2_energy_plateau(duplicate_trace;
        r2_tolerance=1.0e-2, energy_tolerance=1.0e-8, window=3).available
    @test_throws ArgumentError mod.bng_r2_energy_plateau(stable_trace;
        r2_tolerance=1.0e-2, energy_tolerance=1.0e-8, window=1)
    @test mod.bng_convergence_policy(:strict_joint_force) ==
        :strict_joint_force
    @test_throws ArgumentError mod.bng_convergence_policy(:unknown)

    dims = (8, 10, 2)
    phi = [0.36 + 0.04 * cospi(2 * (i - 1) / dims[1]) +
        0.03 * cospi(2 * (j - 1) / dims[2])
        for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3]]
    movable = trues(dims)
    modes = mod.bng_translation_modes(phi, movable)
    @test length(modes) == 2
    @test all(isapprox(norm(mode), 1.0; atol=1.0e-12) for mode in modes)
    @test isapprox(dot(modes[1], modes[2]), 0.0; atol=1.0e-12)

    direction = reshape(collect(range(-1.0, 1.0; length=prod(dims))), dims)
    mod.bng_project_gauge!(direction, movable, modes)
    @test isapprox(mean(direction), 0.0; atol=1.0e-12)
    @test all(isapprox(dot(vec(direction), mode), 0.0; atol=1.0e-12)
        for mode in modes)

    @test mod.bng_invariant_axes(phi) == (3,)
    reduced = mod.bng_reduce_invariant(phi, (3,))
    @test size(reduced) == (8, 10, 1)
    @test mod.bng_expand_invariant(reduced, dims, (3,)) == phi

    bounds = (0.98, 1.04)
    for cell in (0.981, 1.0, 1.015, 1.039)
        coordinate = mod.bng_cell_factor_coordinate(cell, bounds)
        recovered, derivative =
            mod.bng_cell_factor_from_coordinate(coordinate, bounds)
        @test recovered ≈ cell atol=1.0e-14
        @test derivative > 0.0
    end
    @test collect(mod.bng_scaled_lengths(
        (2.0, 3.0, 4.0), 1.1, (1, 2))) ≈ [2.2, 3.3, 4.0]
    @test_throws ArgumentError mod.bng_cell_factor_coordinate(0.98, bounds)

    diagonal = collect(range(1.0, 2.0; length=80))
    function apply_diagonal(vector)
        field = reshape(vector, dims)
        plane = vec(field[:, :, 1])
        return vec(repeat(reshape(diagonal .* plane, 8, 10, 1), 1, 1, 2))
    end
    rhs_plane = collect(range(-0.5, 0.5; length=80))
    rhs = vec(repeat(reshape(rhs_plane, 8, 10, 1), 1, 1, 2))
    dense = mod.bng_dense_svd_solve(apply_diagonal, rhs, dims;
        invariant_axes=(3,))
    @test dense.variable_count == 80
    @test dense.numerical_rank == 80
    @test dense.relative_residual < 1.0e-12
    @test vec(reshape(dense.solution, dims)[:, :, 1]) ≈ rhs_plane ./ diagonal

    joint = mod.polish_bvk1_joint_cell_lbfgs(phi;
        morphology="CYL", fA=0.36, chiN=16.0,
        base_lengths=(1.7, 2.9, 1.0),
        initial_cell_factor=1.0, cell_factor_bounds=(0.98, 1.02),
        scaled_axes=(1, 2), invariant_axes=(3,), laplacian_smoothing=0.1,
        max_iterations=3, force_tolerance=0.0,
        cell_stress_tolerance=0.0, force_check_interval=1,
        verbose=false)
    @test 0.98 < joint.cell_factor < 1.02
    @test collect(joint.lengths) ≈ collect(mod.bng_scaled_lengths(
        (1.7, 2.9, 1.0), joint.cell_factor, (1, 2)))
    @test isfinite(joint.final_energy_density)
    @test isfinite(joint.cell_stress)
    @test mean(joint.phi) ≈ 0.36 atol=1.0e-12
    @test joint.phi[:, :, 1] ≈ joint.phi[:, :, 2]
    @test joint.objective_calls > 0
    @test joint.gradient_calls > 0

    gyr_phi = [0.36 + 0.02 * cospi(2 * (i - 1) / 4) *
        cospi(2 * (j - 1) / 4) * cospi(2 * (k - 1) / 4)
        for i in 1:4, j in 1:4, k in 1:4]
    gyr_joint = mod.polish_bvk1_joint_cell_lbfgs(gyr_phi;
        morphology="GYR", fA=0.36, chiN=16.0,
        base_lengths=(3.6, 3.6, 3.6), initial_cell_factor=1.0,
        cell_factor_bounds=(0.98, 1.02), laplacian_smoothing=0.1,
        max_iterations=1, force_tolerance=0.0,
        cell_stress_tolerance=0.0, force_check_interval=1,
        verbose=false)
    @test size(gyr_joint.phi) == size(gyr_phi)
    @test gyr_joint.invariant_axes == ()
    @test all(isapprox(length, 3.6 * gyr_joint.cell_factor;
        atol=1.0e-12) for length in gyr_joint.lengths)
    @test isapprox(mean(gyr_joint.phi), 0.36; atol=1.0e-12)

    function mock_block_field(phi, lengths; converged::Bool,
            iterations::Int, kkt_norm::Float64, kkt_maxabs::Float64,
            energy::Float64)
        trace = [
            (iteration=0, energy_density=energy + 1.0e-5,
                kkt_norm=10.0 * kkt_norm,
                kkt_maxabs=10.0 * kkt_maxabs,
                stop_reason="initial_state"),
            (iteration=iterations, energy_density=energy,
                kkt_norm=kkt_norm, kkt_maxabs=kkt_maxabs,
                stop_reason=converged ? "force_tol" :
                    "optimizer_iteration_cap"),
        ]
        return (phi=phi, lengths=lengths, history=trace,
            initial_kkt_norm=trace[1].kkt_norm,
            initial_kkt_maxabs=trace[1].kkt_maxabs,
            final_kkt_norm=kkt_norm, final_kkt_maxabs=kkt_maxabs,
            converged=converged,
            stop_reason=converged ? "force_tol" :
                "optimizer_iteration_cap",
            optimizer_iterations=iterations,
            optimizer_converged=converged,
            objective_calls=iterations + 1, gradient_calls=iterations)
    end
    function run_block_rollback_case(reason::Symbol)
        calls = Ref(0)
        accepted_phi = fill(0.20, 2, 2, 2)
        rejected_phi = reshape(repeat([0.10, 0.30], 4), 2, 2, 2)
        energy(field, lengths) = (lengths[1] - 1.02)^2 +
            10.0 * var(vec(field))
        function field_polish(_phi; lengths, kwargs...)
            calls[] += 1
            first_call = calls[] == 1
            candidate = first_call ? accepted_phi : rejected_phi
            converged = first_call || reason == :energy_increase
            kkt_norm = converged ? 5.0e-4 : 2.0e-2
            kkt_maxabs = converged ? 7.0e-4 : 3.0e-2
            return mock_block_field(candidate, lengths;
                converged=converged,
                iterations=first_call ? 5 : 7,
                kkt_norm=kkt_norm, kkt_maxabs=kkt_maxabs,
                energy=energy(candidate, lengths))
        end
        result = mod.polish_bvk1_block_cell_lbfgs(accepted_phi;
            morphology="GYR", fA=0.20, chiN=15.0,
            base_lengths=(1.0, 1.0, 1.0), initial_cell_factor=1.0,
            cell_factor_bounds=(0.98, 1.04), max_cell_iterations=3,
            max_relative_cell_step=5.0e-3,
            force_tolerance=1.0e-3, cell_stress_tolerance=1.0e-4,
            field_polish=field_polish, energy_evaluator=energy,
            checkpoint_residual_evaluator=(_field, _lengths) ->
                (kkt_norm=1.04e-3, kkt_maxabs=4.64e-3),
            verbose=false)
        return result, calls[]
    end

    failed_trial, failed_calls = run_block_rollback_case(:field_failure)
    @test failed_calls == 2
    @test failed_trial.stop_reason ==
        "cell_trial_field_not_converged_rolled_back"
    @test failed_trial.cell_factor == 1.0
    @test failed_trial.phi == fill(0.20, 2, 2, 2)
    @test failed_trial.final_kkt_norm == 5.0e-4
    @test failed_trial.history[1].checkpoint_committed
    @test !failed_trial.history[2].checkpoint_committed
    @test failed_trial.history[2].rollback_reason == "field_not_converged"

    higher_trial, higher_calls = run_block_rollback_case(:energy_increase)
    @test higher_calls == 2
    @test higher_trial.stop_reason ==
        "cell_trial_energy_increase_rolled_back"
    @test higher_trial.cell_factor == 1.0
    @test higher_trial.phi == fill(0.20, 2, 2, 2)
    @test higher_trial.history[2].field_converged
    @test !higher_trial.history[2].checkpoint_committed
    @test higher_trial.history[2].rollback_reason == "energy_increase"
    @test propertynames(first(higher_trial.field_history)) == (
        :cumulative_iteration, :cell_iteration, :field_iteration,
        :cell_factor, :energy_density, :kkt_norm, :kkt_maxabs,
        :checkpoint_committed, :rollback_reason, :stop_reason)
    @test getproperty.(higher_trial.field_history,
        :cumulative_iteration) == [0, 0, 5, 5, 12]
    @test getproperty.(higher_trial.field_history,
        :checkpoint_committed) == [true, true, true, false, false]

    baseline_phi = fill(0.20, 2, 2, 2)
    degraded_phi = reshape(repeat([0.10, 0.30], 4), 2, 2, 2)
    baseline_energy(field, lengths) = 4.0e-6 * log(lengths[1]) +
        10.0 * var(vec(field))
    degraded_polish(_phi; lengths, kwargs...) = mock_block_field(
        degraded_phi, lengths; converged=false, iterations=7,
        kkt_norm=2.0e-2, kkt_maxabs=3.0e-2,
        energy=baseline_energy(degraded_phi, lengths))
    input_rollback = mod.polish_bvk1_block_cell_lbfgs(baseline_phi;
        morphology="GYR", fA=0.20, chiN=15.0,
        base_lengths=(1.0, 1.0, 1.0), initial_cell_factor=1.0,
        cell_factor_bounds=(0.98, 1.04), max_cell_iterations=3,
        force_tolerance=1.0e-3, cell_stress_tolerance=1.0e-4,
        field_polish=degraded_polish, energy_evaluator=baseline_energy,
        checkpoint_residual_evaluator=(_field, _lengths) ->
            (kkt_norm=1.04e-3, kkt_maxabs=4.64e-3),
        verbose=false)
    @test input_rollback.stop_reason ==
        "cell_trial_field_not_converged_rolled_back"
    @test input_rollback.phi == baseline_phi
    @test input_rollback.cell_factor == 1.0
    @test input_rollback.initial_kkt_norm == 1.04e-3
    @test input_rollback.initial_kkt_maxabs == 4.64e-3
    @test input_rollback.final_kkt_norm == 1.04e-3
    @test input_rollback.final_kkt_maxabs == 4.64e-3
    @test isapprox(input_rollback.initial_cell_stress, 4.0e-6;
        atol=1.0e-10)
    @test input_rollback.final_energy_density ==
        input_rollback.initial_energy_density

    relaxed_calls = Ref(0)
    relaxed_polish(phi; lengths, kwargs...) = begin
        relaxed_calls[] += 1
        mock_block_field(phi, lengths; converged=false, iterations=5,
            kkt_norm=5.0e-3, kkt_maxabs=8.0e-3, energy=-0.1)
    end
    relaxed_checkpoint = mod.polish_bvk1_block_cell_lbfgs(baseline_phi;
        morphology="GYR", fA=0.20, chiN=15.0,
        base_lengths=(1.0, 1.0, 1.0), initial_cell_factor=1.0,
        cell_factor_bounds=(0.98, 1.04), max_cell_iterations=3,
        force_tolerance=1.0e-3, checkpoint_force_tolerance=1.0e-2,
        cell_stress_tolerance=1.0e-4, field_polish=relaxed_polish,
        energy_evaluator=(_field, _lengths) -> -0.1,
        checkpoint_residual_evaluator=(_field, _lengths) ->
            (kkt_norm=2.0e-2, kkt_maxabs=2.0e-2), verbose=false)
    @test relaxed_calls[] == 1
    @test !relaxed_checkpoint.converged
    @test relaxed_checkpoint.stop_reason ==
        "checkpoint_force_and_cell_tol"
    @test relaxed_checkpoint.final_kkt_norm == 5.0e-3
    @test relaxed_checkpoint.final_kkt_maxabs == 8.0e-3

    smoothed_policy_seen = Ref{Any}(nothing)
    smoothed_polish(phi; lengths, convergence_policy,
            smoothed_r2_tolerance, energy_plateau_tolerance,
            energy_plateau_window, kwargs...) = begin
        smoothed_policy_seen[] = (convergence_policy,
            smoothed_r2_tolerance, energy_plateau_tolerance,
            energy_plateau_window)
        merge(mock_block_field(phi, lengths; converged=true, iterations=5,
            kkt_norm=5.0e-3, kkt_maxabs=5.0,
            energy=-0.1), (phase_energy_converged=true,
            energy_plateau_available=true, energy_plateau_span=5.0e-9,
            energy_plateau_max_abs_delta=4.0e-9))
    end
    smoothed_checkpoint = mod.polish_bvk1_block_cell_lbfgs(baseline_phi;
        morphology="GYR", fA=0.20, chiN=15.0,
        base_lengths=(1.0, 1.0, 1.0), initial_cell_factor=1.0,
        cell_factor_bounds=(0.98, 1.04), max_cell_iterations=1,
        force_tolerance=1.0e-3, checkpoint_force_tolerance=1.0e-2,
        cell_stress_tolerance=1.0e-4, field_polish=smoothed_polish,
        convergence_policy=:smoothed_r2_energy_plateau,
        smoothed_r2_tolerance=1.0e-2,
        energy_plateau_tolerance=1.0e-8, energy_plateau_window=3,
        energy_evaluator=(_field, _lengths) -> -0.1,
        checkpoint_residual_evaluator=(_field, _lengths) ->
            (kkt_norm=2.0e-2, kkt_maxabs=20.0), verbose=false)
    @test smoothed_policy_seen[] ==
        (:smoothed_r2_energy_plateau, 1.0e-2, 1.0e-8, 3)
    @test smoothed_checkpoint.converged
    @test smoothed_checkpoint.phase_energy_converged
    @test smoothed_checkpoint.energy_plateau_available
    @test smoothed_checkpoint.energy_plateau_span == 5.0e-9
    @test smoothed_checkpoint.energy_plateau_max_abs_delta == 4.0e-9
    @test smoothed_checkpoint.final_kkt_maxabs == 5.0
    @test smoothed_checkpoint.stop_reason ==
        "smoothed_r2_energy_plateau_and_cell_tol"
end

@testset "BVK1 allocation-free periodic neighbor indexing" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_nd_force_consistency.jl")
    mod = Module(:TestBVK1PeriodicNeighborIndexing)
    previous = get(ENV,
        "DFM_SKIP_BURP_BVK1_ND_FORCE_CONSISTENCY_MAIN", nothing)
    try
        ENV["DFM_SKIP_BURP_BVK1_ND_FORCE_CONSISTENCY_MAIN"] = "1"
        Base.include(mod, script_path)
    finally
        if previous === nothing
            delete!(ENV, "DFM_SKIP_BURP_BVK1_ND_FORCE_CONSISTENCY_MAIN")
        else
            ENV["DFM_SKIP_BURP_BVK1_ND_FORCE_CONSISTENCY_MAIN"] = previous
        end
    end

    dims = (3, 4, 2)
    for index in CartesianIndices(dims), dim in 1:3
        expected = collect(Tuple(index))
        expected[dim] = mod1(expected[dim] + 1, dims[dim])
        @test mod._periodic_forward_index(index, dims, dim) ==
            CartesianIndex(Tuple(expected))
    end
end

@testset "BVK1 exact nonsmooth KKT audit" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_bvk1_nonsmooth_kkt_audit.jl")
    mod = Module(:TestBVK1NonsmoothKKTAudit)
    previous = get(ENV, "DFM_SKIP_BVK1_NONSMOOTH_KKT_AUDIT_MAIN", nothing)
    try
        ENV["DFM_SKIP_BVK1_NONSMOOTH_KKT_AUDIT_MAIN"] = "1"
        Base.include(mod, script_path)
    finally
        if previous === nothing
            delete!(ENV, "DFM_SKIP_BVK1_NONSMOOTH_KKT_AUDIT_MAIN")
        else
            ENV["DFM_SKIP_BVK1_NONSMOOTH_KKT_AUDIT_MAIN"] = previous
        end
    end

    dims = (4, 6, 2)
    phi = [0.36 + 0.02 * cospi(2 * (i - 1) / dims[1]) +
        0.015 * sinpi(2 * (j - 1) / dims[2]) +
        (k == 1 ? 0.0037 : -0.0037)
        for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3]]
    lengths = (1.3, 2.1, 0.9)
    @test isapprox(mean(phi), 0.36; atol=1.0e-14)

    no_cusp = mod.bnk_nonsmooth_kkt_audit(phi; f=0.36, chiN=16.0,
        lengths=lengths, cusp_absolute_tolerance=0.0,
        cusp_relative_tolerance=0.0, max_iterations=50)
    @test isempty(no_cusp.cusp.cusp_indices)
    @test no_cusp.adjusted.norm == no_cusp.ordinary.norm

    all_cusp = mod.bnk_nonsmooth_kkt_audit(phi; f=0.36, chiN=16.0,
        lengths=lengths, cusp_absolute_tolerance=Inf,
        cusp_relative_tolerance=0.0, max_iterations=100)
    @test length(all_cusp.cusp.cusp_indices) == length(phi)
    @test all_cusp.adjusted.norm <= all_cusp.ordinary.norm + 1.0e-12

    chemical = reshape(collect(range(-0.004, 0.006; length=length(phi))),
        size(phi))
    authoritative = mod.bnk_bounded_result(vec(chemical), phi, lengths;
        density_floor=1.0e-6, active_tolerance_factor=10.0)
    partition = vec(mod.bnk_bound_partition(phi;
        density_floor=1.0e-6, active_tolerance_factor=10.0))
    fast_residual = similar(vec(chemical))
    fast_multiplier = mod.bnk_bounded_residual!(fast_residual,
        vec(chemical), partition)
    @test authoritative.lower_active_count == 0
    @test authoritative.upper_active_count == 0
    @test fast_multiplier ≈ authoritative.multiplier atol=1.0e-14
    @test fast_residual ≈ vec(authoritative.residual) atol=1.0e-14
    @test sqrt(prod(lengths) / length(phi) * sum(abs2, fast_residual)) ≈
        authoritative.norm atol=1.0e-14
    @test maximum(abs, fast_residual) ≈ authoritative.maxabs atol=1.0e-14
    @test isnothing(mod.bnk_require_physical_density(phi))
    active_phi = fill(0.36, 2, 2, 2)
    active_phi[1] = nextfloat(0.0)
    active_phi[2] = prevfloat(1.0)
    @test isnothing(mod.bnk_require_physical_density(active_phi))
    active_partition = mod.bnk_bound_partition(active_phi;
        density_floor=1.0e-6, active_tolerance_factor=10.0)
    @test active_partition[1] == 0x01
    @test active_partition[2] == 0x02
    for invalid in (0.0, 1.0, -eps(), 1.0 + eps(), NaN)
        invalid_phi = copy(phi)
        invalid_phi[1] = invalid
        @test_throws ErrorException mod.bnk_require_physical_density(
            invalid_phi)
    end

    feasibility_phi = fill(0.36, 2, 2, 2)
    feasibility_phi[1] = 5.0e-7
    feasibility_phi[2] = 1.0 - 5.0e-7
    chemical_zero = [0.002, -0.002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    basis = reshape([-0.001, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        :, 1)
    feasibility = mod.bnk_maxabs_feasibility_audit(chemical_zero, basis,
        feasibility_phi, (1.0, 1.0, 1.0), [0.0];
        maxabs_tolerance=1.0e-3, max_iterations=100)
    @test feasibility.optimizer_converged
    @test all(-1.0 .<= feasibility.subgradients .<= 1.0)
    @test feasibility.adjusted.maxabs <= 1.0e-3 + 1.0e-12
    @test feasibility.adjusted.norm <= 1.0e-3

    source = read(script_path, String)
    @test occursin("LBFGSB()", source)
    @test !occursin("Fminbox(LBFGS())", source)

    metric(norm, maxabs) = (norm=Float64(norm), maxabs=Float64(maxabs))
    feasibility_result(norm, maxabs; converged=true) = (
        adjusted=metric(norm, maxabs), subgradients=[0.25],
        optimizer_converged=converged, optimizer_iterations=3,
        objective_calls=5)
    calls = Ref(0)
    selected = mod.bnk_select_cusp_certificate(
        metric(5.0e-4, 1.2e-3), [0.0];
        feasibility=() -> begin
            calls[] += 1
            feasibility_result(6.0e-4, 8.0e-4)
        end)
    @test calls[] == 1
    @test selected.maxabs_feasibility_attempted
    @test selected.maxabs_feasibility_converged
    @test selected.maxabs_feasibility_pass
    @test selected.certificate_mode == "maxabs_feasibility"
    @test mod.bnk_cusp_gate(selected.adjusted)

    robust_certificates = [mod.bnk_select_cusp_certificate(
        metric(5.0e-4, 1.2e-3), [0.0];
        feasibility=() -> feasibility_result(6.0e-4, maxabs))
        for maxabs in (8.0e-4, 9.0e-4)]
    @test count(certificate -> mod.bnk_cusp_gate(certificate.adjusted),
        robust_certificates) >= 2
    @test all(certificate ->
        certificate.adjusted.norm <= 1.0e-3 &&
        certificate.adjusted.maxabs <= 1.0e-3, robust_certificates)

    rejected = mod.bnk_select_cusp_certificate(
        metric(5.0e-4, 1.2e-3), [0.0];
        feasibility=() -> feasibility_result(6.0e-4, 1.1e-3))
    @test rejected.maxabs_feasibility_attempted
    @test !rejected.maxabs_feasibility_pass
    @test rejected.certificate_mode == "minimum_norm"
    @test !mod.bnk_cusp_gate(rejected.adjusted)

    descent = mod.bnk_nonsmooth_kkt_descent(phi; f=0.36, chiN=16.0,
        lengths=lengths, cusp_absolute_tolerance=Inf,
        force_tolerance=Inf, max_iterations=4, audit_max_iterations=50)
    @test descent.converged
    @test descent.iterations == 0
    @test descent.phi == phi
end

@testset "BVK1 branch stationarity anneal resampling" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "run_bvk1_branch_stationarity_anneal.jl")
    mod = Module(:TestBVK1BranchStationarityAnneal)
    previous = get(ENV,
        "DFM_SKIP_BVK1_BRANCH_STATIONARITY_ANNEAL_MAIN", nothing)
    try
        ENV["DFM_SKIP_BVK1_BRANCH_STATIONARITY_ANNEAL_MAIN"] = "1"
        Base.include(mod, script_path)
    finally
        if previous === nothing
            delete!(ENV, "DFM_SKIP_BVK1_BRANCH_STATIONARITY_ANNEAL_MAIN")
        else
            ENV["DFM_SKIP_BVK1_BRANCH_STATIONARITY_ANNEAL_MAIN"] = previous
        end
    end

    @test mod.bca_default_input_dims("BCC", (32, 32, 32)) == (32, 32, 32)
    @test mod.bca_default_input_dims("GYR", (48, 48, 48)) == (48, 48, 48)
    @test mod.bca_default_input_dims("CYL", (40, 72)) == (40, 72, 2)

    source_dims = (12, 16)
    target_dims = (18, 24)
    periodic(x, y) = 0.4 + 0.07cos(2pi*x) + 0.03sin(4pi*y) +
        0.02cos(2pi*(x+y))
    source_centered = [periodic((i-0.5)/source_dims[1],
        (j-0.5)/source_dims[2]) for i in 1:source_dims[1],
        j in 1:source_dims[2]]
    expected_centered = [periodic((i-0.5)/target_dims[1],
        (j-0.5)/target_dims[2]) for i in 1:target_dims[1],
        j in 1:target_dims[2]]
    centered = mod.bca_cell_centered_resample(source_centered, target_dims)
    @test maximum(abs, centered .- expected_centered) <= 1.0e-12
    @test isapprox(mean(centered), mean(source_centered); atol=1.0e-13)

    theta_periodic(x, y) = 0.61 + 0.08cos(2pi*x) + 0.04sin(2pi*y)
    theta_source = [theta_periodic((i-0.5)/source_dims[1],
        (j-0.5)/source_dims[2]) for i in 1:source_dims[1],
        j in 1:source_dims[2]]
    theta_expected = [theta_periodic((i-0.5)/target_dims[1],
        (j-0.5)/target_dims[2]) for i in 1:target_dims[1],
        j in 1:target_dims[2]]
    theta_phi = sin.(theta_source).^2
    theta_target_mean = mean(sin.(theta_expected).^2)
    theta_resampled = mod.bca_resample_theta(theta_phi, target_dims,
        theta_target_mean)
    @test maximum(abs, theta_resampled .- sin.(theta_expected).^2) <= 1.0e-12
    @test isapprox(mean(theta_resampled), theta_target_mean; atol=1.0e-13)
    @test 0.0 < minimum(theta_resampled) < maximum(theta_resampled) < 1.0

    phi = [0.36 + 0.03 * cospi(2 * (i - 1) / 8) +
        0.02 * sinpi(2 * (j - 1) / 10) +
        (k == 1 ? 0.001 : -0.001)
        for i in 1:8, j in 1:10, k in 1:2]
    resampled = mod.bca_resample_branch(phi, (6, 9), 0.36;
        morphology="CYL", density_floor=1.0e-6)
    @test size(resampled) == (6, 9, 2)
    @test isapprox(mean(resampled), 0.36; atol=1.0e-12)
    @test maximum(abs, resampled[:, :, 1] .- resampled[:, :, 2]) == 0.0
    @test 1.0e-6 <= minimum(resampled) < maximum(resampled) <= 1.0 - 1.0e-6

    gyr = mod.bca_resample_branch(phi, (6, 7, 5), 0.36;
        morphology="GYR", density_floor=1.0e-6)
    @test size(gyr) == (6, 7, 5)
    @test isapprox(mean(gyr), 0.36; atol=1.0e-12)
    @test 1.0e-6 <= minimum(gyr) < maximum(gyr) <= 1.0 - 1.0e-6

    shifted = mod.bca_resample_branch(phi, size(phi), 0.46;
        morphology="GYR", density_floor=1.0e-6)
    linear = mod.BCA_SOLVER.DFMMonteCarlo._project_box_mean_density(phi,
        0.46; lo=1.0e-6, hi=1.0 - 1.0e-6)
    @test size(shifted) == size(phi)
    @test isapprox(mean(shifted), 0.46; atol=1.0e-12)
    @test maximum(abs, shifted .- linear) > 1.0e-6
end

@testset "Polyorder nonlamellar BCC identity gate" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_unrestricted_relaxation.jl")
    mod = Module(:TestPolyorderNonlamellarBCCIdentityGate)
    previous = get(ENV,
        "DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_UNRESTRICTED_MAIN",
        nothing)
    try
        ENV["DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_UNRESTRICTED_MAIN"] = "1"
        Base.include(mod, script_path)
    finally
        if previous === nothing
            delete!(ENV,
                "DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_UNRESTRICTED_MAIN")
        else
            ENV["DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_UNRESTRICTED_MAIN"] =
                previous
        end
    end

    bcc_reference_like = (
        dominant_shell_label="110",
        dominant_shell_power_fraction=0.72,
        negative_euler_double_network=false,
    )
    bcc_homogeneous = (
        dominant_shell_label="0",
        dominant_shell_power_fraction=0.0,
        negative_euler_double_network=false,
    )
    gyr_without_double_network = (
        dominant_shell_label="211",
        dominant_shell_power_fraction=0.90,
        negative_euler_double_network=false,
    )

    @test mod._identity_preserved("BCC", bcc_reference_like, "110";
        contrast=0.2)
    @test !mod._identity_preserved("BCC", bcc_reference_like, "110";
        contrast=1.0e-6)
    @test !mod._identity_preserved("BCC", bcc_homogeneous, "110")
    @test !mod._identity_preserved("GYR", gyr_without_double_network, "211")
end

@testset "Polyorder nonlamellar stress-free relaxation source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_stress_free_relaxation.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_STRESS_FREE_MAIN",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_unrestricted_relaxation.jl",
        source)
    @test occursin("scalar_cell_coordinate_descent", source)
    @test occursin("cell_factor_bracket", source)
    @test occursin("stress_free_relaxation_summary.csv", source)
    @test occursin("stress_free_relaxation_paths.csv", source)
    @test occursin("stress_free_relaxation_cell_scan.csv", source)
    @test occursin("stress_free_relaxation.svg", source)
    @test occursin("stress_free_relaxation.md", source)
    @test occursin("_arg_value(\"--starts\", \"branch_selected,scft_point\")",
        source)
    @test occursin("identity_preserved", source)
    @test occursin("endpoint_limited", source)
    @test occursin("cell_minimum", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_nonlamellar_stress_free_relaxation.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("stress-free", doc_text)
    @test occursin("fixed-composition", doc_text)
    @test occursin("Polyorder SCFT", doc_text)
    @test occursin("scalar cell factor", doc_text)
    @test occursin("BURP-TI", doc_text)
    @test occursin("BVK1", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("stress_free_relaxation", packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("stress_free_relaxation_summary.csv", readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("polyorder_scft_nonlamellar_stress_free_relaxation",
        coverage_source)
end

@testset "Polyorder nonlamellar GYR convergence ladder source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_gyr_stress_free_convergence_ladder.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_GYR_STRESS_FREE_LADDER_MAIN",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_stress_free_relaxation.jl",
        source)
    @test occursin("field_iteration_ladder", source)
    @test occursin("force_reduction_ratio", source)
    @test occursin("gyr_stress_free_convergence_ladder_summary.csv", source)
    @test occursin("gyr_stress_free_convergence_ladder_paths.csv", source)
    @test occursin("gyr_stress_free_convergence_ladder.svg", source)
    @test occursin("gyr_stress_free_convergence_ladder.md", source)
    @test occursin("identity_preserved", source)
    @test occursin("cell_minimum", source)
    @test occursin("stationarity_candidate", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_gyr_stress_free_convergence_ladder.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("GYR", doc_text)
    @test occursin("convergence ladder", doc_text)
    @test occursin("projected-force", doc_text)
    @test occursin("Polyorder SCFT", doc_text)
    @test occursin("BURP-TI", doc_text)
    @test occursin("BVK1", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("gyr_stress_free_convergence_ladder",
        packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("gyr_stress_free_convergence_ladder_summary.csv",
        readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("polyorder_scft_gyr_stress_free_convergence_ladder",
        coverage_source)
end

@testset "Polyorder nonlamellar BCC convergence ladder source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_bcc_stress_free_convergence_ladder.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_BCC_STRESS_FREE_LADDER_MAIN",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_stress_free_relaxation.jl",
        source)
    @test occursin("scft_point", source)
    @test occursin("field_iteration_ladder", source)
    @test occursin("force_reduction_ratio", source)
    @test occursin("bcc_stress_free_convergence_ladder_summary.csv", source)
    @test occursin("bcc_stress_free_convergence_ladder_paths.csv", source)
    @test occursin("bcc_stress_free_convergence_ladder.svg", source)
    @test occursin("bcc_stress_free_convergence_ladder.md", source)
    @test occursin("identity_preserved", source)
    @test occursin("cell_minimum", source)
    @test occursin("stationarity_candidate", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_bcc_stress_free_convergence_ladder.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("BCC", doc_text)
    @test occursin("convergence ladder", doc_text)
    @test occursin("{110}", doc_text)
    @test occursin("Polyorder SCFT", doc_text)
    @test occursin("BURP-TI", doc_text)
    @test occursin("BVK1", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("bcc_stress_free_convergence_ladder",
        packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("bcc_stress_free_convergence_ladder_summary.csv",
        readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("polyorder_scft_bcc_stress_free_convergence_ladder",
        coverage_source)
end

@testset "Polyorder nonlamellar logit LBFGS optimizer probe source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_logit_lbfgs_probe.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_LOGIT_LBFGS_MAIN",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_stress_free_relaxation.jl",
        source)
    @test occursin("Optim.LBFGS", source)
    @test occursin("fixed_mean_logit_lbfgs", source)
    @test occursin("_phi_from_logit_mean", source)
    @test occursin("_logit_mean_gradient!", source)
    @test occursin("nonlamellar_logit_lbfgs_probe_summary.csv", source)
    @test occursin("nonlamellar_logit_lbfgs_probe_paths.csv", source)
    @test occursin("nonlamellar_logit_lbfgs_probe.svg", source)
    @test occursin("nonlamellar_logit_lbfgs_probe.md", source)
    @test occursin("identity_preserved", source)
    @test occursin("stationarity_candidate", source)
    @test occursin("endpoint_limited", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_nonlamellar_logit_lbfgs_probe.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("logit", doc_text)
    @test occursin("L-BFGS", doc_text)
    @test occursin("fixed composition", doc_text)
    @test occursin("phase ranking", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("nonlamellar_logit_lbfgs_probe", packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("polyorder_scft_nonlamellar_logit_lbfgs_probe",
        readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("nonlamellar_logit_lbfgs_probe_summary.csv",
        coverage_source)
end

@testset "Polyorder GYR logit cell scan source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_gyr_logit_cell_scan.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_GYR_LOGIT_CELL_SCAN_MAIN",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_logit_lbfgs_probe.jl",
        source)
    @test occursin("Optim.LBFGS", source)
    @test occursin("cell_factor_grid", source)
    @test occursin("gyr_logit_cell_scan_rows.csv", source)
    @test occursin("gyr_logit_cell_scan_summary.csv", source)
    @test occursin("gyr_logit_cell_scan.svg", source)
    @test occursin("gyr_logit_cell_scan.md", source)
    @test occursin("stress_free_cell_candidate", source)
    @test occursin("stationarity_candidate", source)
    @test occursin("cell_stress_proxy", source)
    @test occursin("_quadratic_cell_refinement", source)
    @test occursin("refined_cell_factor", source)
    @test occursin("refined_energy_density", source)
    @test occursin("refined_cell_curvature", source)
    @test occursin("refined_cell_stress_proxy", source)
    @test occursin("identity_preserved", source)
    @test occursin("grid_boundary_limited", source)
    @test occursin("--gyr-reference", source)
    @test occursin("--gyr-target-dims", source)
    @test occursin("--gyr-f", source)
    @test occursin("--gyr-chiN", source)
    @test occursin("source_reference", source)
    @test occursin("_needs_branch_selected_start", source)
    @test occursin("Polyorder_results/GYR_f0.39_chiN20/GYR_opt.mat", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_gyr_logit_cell_scan.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("GYR", doc_text)
    @test occursin("logit", doc_text)
    @test occursin("cubic cell", doc_text)
    @test occursin("stress-free", doc_text)
    @test occursin("phase ranking", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("gyr_logit_cell_scan", packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("polyorder_scft_gyr_logit_cell_scan", readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("gyr_logit_cell_scan_summary.csv", coverage_source)
end

@testset "Polyorder BCC logit cell scan source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_bcc_logit_cell_scan.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_BCC_LOGIT_CELL_SCAN_MAIN",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_logit_lbfgs_probe.jl",
        source)
    @test occursin("Optim.LBFGS", source)
    @test occursin("BCC", source)
    @test occursin("{110}", source)
    @test occursin("bcc_logit_cell_scan_rows.csv", source)
    @test occursin("bcc_logit_cell_scan_summary.csv", source)
    @test occursin("bcc_logit_cell_scan.svg", source)
    @test occursin("bcc_logit_cell_scan.md", source)
    @test occursin("stress_free_cell_candidate", source)
    @test occursin("stationarity_candidate", source)
    @test occursin("cell_stress_proxy", source)
    @test occursin("identity_preserved", source)
    @test occursin("grid_boundary_limited", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_bcc_logit_cell_scan.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("BCC", doc_text)
    @test occursin("logit", doc_text)
    @test occursin("cubic cell", doc_text)
    @test occursin("{110}", doc_text)
    @test occursin("phase ranking", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("bcc_logit_cell_scan", packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("polyorder_scft_bcc_logit_cell_scan", readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("bcc_logit_cell_scan_summary.csv", coverage_source)
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("bcc_logit_cell_scan", pipeline_source)
end

@testset "Polyorder CYL logit cell scan source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_cyl_logit_cell_scan.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_CYL_LOGIT_CELL_SCAN_MAIN",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_logit_lbfgs_probe.jl",
        source)
    @test occursin("Optim.LBFGS", source)
    @test occursin("CYL", source)
    @test occursin("cell_stress_proxy", source)
    @test occursin("_quadratic_cell_refinement", source)
    @test occursin("refined_cell_factor", source)
    @test occursin("refined_energy_density", source)
    @test occursin("refined_cell_curvature", source)
    @test occursin("refined_cell_stress_proxy", source)
    @test occursin("cyl_logit_cell_scan_rows.csv", source)
    @test occursin("cyl_logit_cell_scan_summary.csv", source)
    @test occursin("cyl_logit_cell_scan.svg", source)
    @test occursin("cyl_logit_cell_scan.md", source)
    @test occursin("stress_free_cell_candidate", source)
    @test occursin("stationarity_candidate", source)
    @test occursin("identity_preserved", source)
    @test occursin("grid_boundary_limited", source)
    @test occursin("--cyl-reference", source)
    @test occursin("--cyl-target-dims", source)
    @test occursin("--cyl-f", source)
    @test occursin("--cyl-chiN", source)
    @test occursin("source_reference", source)
    @test occursin("_needs_branch_selected_start", source)
    @test occursin("Polyorder_results/CYL_f0.39_chiN20/CYL_opt.mat",
        source)
    @test occursin("cyl_geometry_admissible", source)
    @test occursin("special_rectangle_2_over_sqrt3", source)
    @test occursin("rhombic_hexagonal", source)
    @test occursin("geometry_convention_invalid", source)
    @test occursin("polyorder_hexrect_sqrt3", source)
    @test !occursin("not_user_convention", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_cyl_logit_cell_scan.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("CYL", doc_text)
    @test occursin("logit", doc_text)
    @test occursin("cubic cell", doc_text)
    @test occursin("endpoint", doc_text)
    @test occursin("phase ranking", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("cyl_logit_cell_scan", packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("polyorder_scft_cyl_logit_cell_scan", readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("cyl_logit_cell_scan_summary.csv", coverage_source)
    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("cyl_logit_cell_scan", pipeline_source)

    result_dir = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_polyorder_cyl_logit_cell_scan")
    summary_path = joinpath(result_dir, "cyl_logit_cell_scan_summary.csv")
    result_md_path = joinpath(result_dir, "cyl_logit_cell_scan.md")
    @test isfile(summary_path)
    @test isfile(result_md_path)
    summary = isfile(summary_path) ? read(summary_path, String) : ""
    result_md = isfile(result_md_path) ? read(result_md_path, String) : ""
    @test occursin("refined_cell_factor", summary)
    @test occursin("refined_cell_curvature", summary)
    @test occursin("cyl_geometry_admissible", summary)
    @test occursin("cyl_geometry_convention", summary)
    @test occursin("phase_rank_cyl_candidate", summary)
    @test occursin("polyorder_hexrect_sqrt3", summary)
    @test occursin("Endpoint saturation and projected-force residuals are diagnostics here, not hard acceptance gates",
        result_md)
    @test !occursin("remain away from the density endpoint", result_md)
    @test !occursin("Endpoint-limited rows are reported as `endpoint_limited_cell_scan`",
        result_md)
end

@testset "Polyorder nonlamellar KKT gate source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_kkt_gate.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_KKT_GATE_MAIN",
        source)
    @test occursin("fixed-mean bound-constrained KKT", source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_stress_free_relaxation.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_bcc_stress_free_convergence_ladder.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_gyr_stress_free_convergence_ladder.jl",
        source)
    @test occursin("bcc_stress_free_convergence_ladder_summary.csv",
        source)
    @test occursin("_latest_bcc_ladder_rows", source)
    @test occursin("bcc_ladder_", source)
    @test occursin("active_lower_count", source)
    @test occursin("active_upper_count", source)
    @test occursin("free_projected_residual_norm", source)
    @test occursin("kkt_violation_max", source)
    @test occursin("stationarity_candidate", source)
    @test occursin("accepted_saddle_candidate", source)
    @test occursin("polyorder_nonlamellar_kkt_gate_rows.csv", source)
    @test occursin("polyorder_nonlamellar_kkt_gate_summary.csv", source)
    @test occursin("polyorder_nonlamellar_kkt_gate.svg", source)
    @test occursin("polyorder_nonlamellar_kkt_gate.md", source)
    @test occursin("_summary_best_row", source)
    @test occursin("identity_preserved", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_nonlamellar_kkt_gate.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("KKT", doc_text)
    @test occursin("endpoint", doc_text)
    @test occursin("CYL", doc_text)
    @test occursin("BCC", doc_text)
    @test occursin("GYR", doc_text)
    @test occursin("phase ranking", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("polyorder_nonlamellar_kkt_gate", packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("polyorder_scft_nonlamellar_kkt_gate",
        readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("polyorder_nonlamellar_kkt_gate_summary.csv",
        coverage_source)
end

@testset "Polyorder nonlamellar KKT tolerance sensitivity source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_kkt_tolerance_sensitivity.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_KKT_TOLERANCE_MAIN",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_kkt_gate.jl",
        source)
    @test occursin("active_tolerance_values", source)
    @test occursin("kkt_tolerance", source)
    @test occursin("min_kkt_violation", source)
    @test occursin("max_kkt_violation", source)
    @test occursin("accepted_saddle_candidate_count", source)
    @test occursin("CYL", source)
    @test occursin("polyorder_nonlamellar_kkt_tolerance_rows.csv", source)
    @test occursin("polyorder_nonlamellar_kkt_tolerance_summary.csv", source)
    @test occursin("polyorder_nonlamellar_kkt_tolerance_sensitivity.svg",
        source)
    @test occursin("polyorder_nonlamellar_kkt_tolerance_sensitivity.md",
        source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_nonlamellar_kkt_tolerance_sensitivity.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("active-bound tolerance", doc_text)
    @test occursin("KKT", doc_text)
    @test occursin("CYL", doc_text)
    @test occursin("robustness", doc_text)
    @test occursin("phase ranking", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("kkt_tolerance_sensitivity", packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("polyorder_scft_nonlamellar_kkt_tolerance_sensitivity",
        readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("polyorder_nonlamellar_kkt_tolerance_summary.csv",
        coverage_source)
end

@testset "Polyorder nonlamellar energy decomposition source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_energy_decomposition.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_ENERGY_DECOMP_MAIN",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_reference_comparison.jl",
        source)
    @test occursin("_burp_ti_energy_terms", source)
    @test occursin("_bvk1_energy_terms", source)
    @test occursin("entropic_work_density", source)
    @test occursin("connectivity_density", source)
    @test occursin("local_entropy_density", source)
    @test occursin("adaptive_gradient_density", source)
    @test occursin("interaction_density", source)
    @test occursin("positive_terms_density", source)
    @test occursin("energy_closure_error", source)
    @test occursin("mechanism_label", source)
    @test occursin("CYL", source)
    @test occursin("BCC", source)
    @test occursin("GYR", source)
    @test occursin("polyorder_nonlamellar_energy_decomposition_rows.csv",
        source)
    @test occursin("polyorder_nonlamellar_energy_decomposition_summary.csv",
        source)
    @test occursin("polyorder_nonlamellar_energy_decomposition.svg",
        source)
    @test occursin("polyorder_nonlamellar_energy_decomposition.md",
        source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_nonlamellar_energy_decomposition.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("energy decomposition", doc_text)
    @test occursin("BURP-TI", doc_text)
    @test occursin("BVK1", doc_text)
    @test occursin("CYL", doc_text)
    @test occursin("BCC", doc_text)
    @test occursin("GYR", doc_text)
    @test occursin("phase ranking", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("polyorder_nonlamellar_energy_decomposition",
        packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("polyorder_scft_nonlamellar_energy_decomposition",
        readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("polyorder_nonlamellar_energy_decomposition_summary.csv",
        coverage_source)
end

@testset "Polyorder nonlamellar manuscript synthesis panel source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_polyorder_nonlamellar_manuscript_panel.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_POLYORDER_NONLAMELLAR_PANEL_MAIN",
        source)
    @test occursin("nonlamellar_manuscript_panel.csv", source)
    @test occursin("nonlamellar_manuscript_panel.svg", source)
    @test occursin("nonlamellar_manuscript_panel.md", source)
    @test occursin("polyorder_nonlamellar_energy_decomposition_summary.csv",
        source)
    @test occursin("polyorder_nonlamellar_kkt_gate_summary.csv", source)
    @test occursin("bcc_stress_free_convergence_ladder_summary.csv", source)
    @test occursin("gyr_stress_free_convergence_ladder_summary.csv", source)
    @test occursin("gyr_logit_cell_scan_summary.csv", source)
    @test occursin("gyr_logit_cell_force_norm", source)
    @test occursin("gyr_logit_cell_factor", source)
    @test occursin("gyr_logit_cell_stress_free_candidate", source)
    @test occursin("gyr_branch_stress_free_cell_candidate", source)
    @test occursin("FINITE_AMPLITUDE_CONTRAST_FLOOR", source)
    @test occursin("DIS_DEGENERACY_ENERGY_TOLERANCE", source)
    @test occursin("finite_amplitude_ordered_candidate", source)
    @test occursin("rejected_dis_degenerate", source)
    @test occursin("claim_boundary", source)
    @test occursin("phase ranking", source)
    @test occursin("DIS/LAM/CYL/BCC/GYR", source)
    @test occursin("CYL", source)
    @test occursin("BCC", source)
    @test occursin("GYR", source)
    @test occursin("mechanism_label", source)
    @test occursin("accepted_saddle_candidate", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "polyorder_scft_nonlamellar_manuscript_panel.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("manuscript panel", doc_text)
    @test occursin("CYL", doc_text)
    @test occursin("BCC", doc_text)
    @test occursin("GYR", doc_text)
    @test occursin("phase ranking", doc_text)
    @test occursin("Uneyama-Doi", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("nonlamellar_manuscript_panel", packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("polyorder_scft_nonlamellar_manuscript_panel",
        readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("nonlamellar_manuscript_panel.csv", coverage_source)
end

@testset "BURP-BVK1 UD phase-space coverage source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_phase_space_coverage.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_UD_PHASE_SPACE_COVERAGE_MAIN",
        source)
    @test occursin("ud_phase_space_coverage.csv", source)
    @test occursin("ud_phase_space_coverage.svg", source)
    @test occursin("ud_phase_space_coverage.md", source)
    @test occursin("ud_ab_phase_diagram_boundaries.csv", source)
    @test occursin("ud_ab_phase_diagram_landmarks.csv", source)
    @test occursin("diblock_stress_free_lamella_comparison", source)
    @test occursin("nonlamellar_manuscript_panel.csv", source)
    @test occursin("convergence_ladder_force_norm", source)
    @test occursin("Uneyama-Doi", source)
    @test occursin("phase-space coverage", source)
    @test occursin("phase ranking", source)
    @test occursin("LAM", source)
    @test occursin("CYL", source)
    @test occursin("BCC", source)
    @test occursin("GYR", source)
    @test occursin("accepted_saddle_candidate", source)
    @test occursin("diagnostic_only", source)
    @test occursin("not_yet", source)

    doc_path = joinpath(@__DIR__, "..", "docs",
        "burp_bvk1_ud_phase_space_coverage.md")
    @test isfile(doc_path)
    doc_text = isfile(doc_path) ? read(doc_path, String) : ""
    @test occursin("phase-space coverage", doc_text)
    @test occursin("Uneyama-Doi", doc_text)
    @test occursin("LAM", doc_text)
    @test occursin("CYL", doc_text)
    @test occursin("BCC", doc_text)
    @test occursin("GYR", doc_text)
    @test occursin("phase ranking", doc_text)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("ud_phase_space_coverage", packet_source)
    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("burp_bvk1_ud_phase_space_coverage", readiness_source)
    coverage_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_ud_reference_coverage.jl"), String)
    @test occursin("ud_phase_space_coverage.csv", coverage_source)
end

@testset "BURP-BVK1 publication packet nonlamellar claim contract" begin
    packet_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl")
    @test isfile(packet_script_path)
    packet_source = isfile(packet_script_path) ?
        read(packet_script_path, String) : ""
    @test occursin("claim_id=\"C8\"", packet_source)
    @test occursin("claim_id=\"C9\"", packet_source)
    @test occursin("claim_id=\"C9a\"", packet_source)
    @test occursin("Polyorder SCFT CYL branch gives saturated HexRect stress-free cell candidates",
        packet_source)
    @test occursin("representative CYL phase-rank pass", packet_source)
    @test occursin("cyl_logit_cell_force_norm", packet_source)
    @test occursin("cyl_logit_cell_scan_summary.csv", packet_source)
    @test occursin("claim_id=\"C10\"", packet_source)
    @test occursin("SCFT-referenced CYL/BCC/GYR", packet_source)
    @test occursin("GYR branch stress-free cubic-cell candidates", packet_source)
    @test occursin("near-DIS BCC", packet_source)
    @test occursin("finite-amplitude", packet_source)
    @test occursin("gyr_logit_cell_force_norm", packet_source)
    @test occursin("gyr_logit_cell_factor", packet_source)
    @test occursin("gyr_logit_cell_scan_summary.csv", packet_source)
    @test occursin("phase ranking remains not_yet", packet_source)
    @test occursin("claim_id=\"C10a\"", packet_source)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate_summary.csv",
        packet_source)
    @test occursin("both BURP-TI and BVK1 select GYR at fA=0.39, chiN=20",
        packet_source)
    @test occursin("claim_id=\"C10b\"", packet_source)
    @test occursin("burp_bvk1_ud_neighbor_f0p41_chiN24_phase_rank_gate",
        packet_source)
    @test occursin("both BURP-TI and BVK1 select GYR at fA=0.41, chiN=24",
        packet_source)
    @test occursin("claim_id=\"C10c\"", packet_source)
    @test occursin("burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate",
        packet_source)
    @test occursin("both BURP-TI and BVK1 select GYR at fA=0.42, chiN=24",
        packet_source)
    @test occursin("limited same-target GYR phase-rank passes", packet_source)
    @test occursin("homogeneous/disordered", packet_source)
    @test occursin("DIS/LAM/CYL/BCC/GYR", packet_source)
    @test occursin("nonlamellar_manuscript_panel.csv", packet_source)
    @test occursin("ud_phase_space_coverage.csv", packet_source)
    @test occursin("accepted_saddle_candidate", packet_source)
    @test occursin("endpoint_rejected", packet_source)
    @test occursin("homogeneous_rejected", packet_source)
    @test occursin("publish_ready_lamellar", packet_source)
    @test occursin("scalar_cell_candidate", packet_source)
    @test occursin("UD-first CYL/BCC method-equity gate", packet_source)
    @test occursin("if UD CYL/BCC reproduction fails under the same workflow, BURP-TI/BVK1 CYL/BCC mismatches are numerical-method evidence first",
        packet_source)
    claim_matrix_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_publication_packet", "claim_matrix.csv")
    figure_inventory_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_publication_packet", "figure_inventory.csv")
    @test isfile(claim_matrix_path)
    @test isfile(figure_inventory_path)
    claim_matrix = isfile(claim_matrix_path) ? read(claim_matrix_path, String) : ""
    figure_inventory = isfile(figure_inventory_path) ?
        read(figure_inventory_path, String) : ""
    @test occursin("C10b", claim_matrix)
    @test occursin("fA=0.41, chiN=24", claim_matrix)
    @test occursin("C10c", claim_matrix)
    @test occursin("fA=0.42, chiN=24", claim_matrix)
    @test occursin("phase_rank_pass 2/2", claim_matrix)
    @test occursin("burp_bvk1_ud_neighbor_f0p41_chiN24_phase_rank_gate",
        claim_matrix)
    @test occursin("UD-neighbor fA=0.41 chiN=24 phase-rank gate",
        figure_inventory)
    @test occursin("UD-neighbor fA=0.42 chiN=24 phase-rank gate",
        figure_inventory)
    @test occursin("phase_rank_gate_summary.csv", figure_inventory)
    mechanism_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_nonlamellar_mechanism_synthesis.jl"), String)
    @test occursin("UD-first CYL/BCC method-equity gate", mechanism_source)
    @test occursin("not Hamiltonian failure evidence by themselves",
        mechanism_source)
end

@testset "BURP-BVK1 manuscript blueprint current nonlamellar evidence contract" begin
    blueprint_path = joinpath(@__DIR__, "..", "docs",
        "burp_bvk1_manuscript_blueprint.md")
    @test isfile(blueprint_path)
    blueprint = isfile(blueprint_path) ? read(blueprint_path, String) : ""
    @test occursin("gyr_logit_cell_scan", blueprint)
    @test occursin("cyl_logit_cell_scan", blueprint)
    @test occursin("bcc_logit_cell_scan", blueprint)
    @test occursin("GYR branch stress-free cubic-cell candidates", blueprint)
    @test occursin("endpoint_limited_cell_scan, force=17.187", blueprint)
    @test occursin("near-DIS logit diagnostic", blueprint)
    @test occursin("contrast=7.8765e-07", blueprint)
    @test occursin("not finite-amplitude", blueprint)
    @test occursin("3.2968e-06", blueprint)
    @test occursin("3.3735e-06", blueprint)
    @test occursin("stress_free_cell_candidate 2/2", blueprint)
    @test occursin("Polyorder_results/BCC/BCC_opt.mat", blueprint)
    @test occursin("same-parameter CYL competitor", blueprint)
    @test occursin("gyr_phase_rank_energy_panel", blueprint)
    @test occursin("0.012194", blueprint)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate", blueprint)
    @test occursin("both BURP-TI and BVK1 select GYR over DIS, LAM, and same-parameter CYL at `f_A=0.39`, `chi N=20`",
        blueprint)
    @test occursin("GYR_ud_neighbor_f0p41_chiN24", blueprint)
    @test occursin("both BURP-TI and BVK1 select GYR over DIS, LAM, and same-parameter CYL at `f_A=0.41`, `chi N=24`",
        blueprint)
    @test occursin("GYR_ud_neighbor_f0p42_chiN24", blueprint)
    @test occursin("both BURP-TI and BVK1 select GYR over DIS, LAM, and same-parameter CYL at `f_A=0.42`, `chi N=24`",
        blueprint)
    @test occursin("phase_rank_completed", blueprint)
    @test occursin("7/7", blueprint)
    @test occursin("phase_rank_gate_summary.csv", blueprint)
    @test !occursin("remaining GYR warm-start SCFT input", blueprint)
    @test !occursin("remaining blockers are same-parameter LAM/DIS generation",
        blueprint)
    @test occursin("not a phase diagram", blueprint)
    @test occursin("not unrestricted stress-free nonlamellar phase ranking",
        blueprint)
    @test occursin("UD-first CYL/BCC method-equity gate", blueprint)
    @test occursin("if UD CYL/BCC reproduction fails under the same workflow",
        blueprint)
end

@testset "BURP-BVK1 nonlamellar publication pipeline source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("PIPELINE_STEPS", source)
    @test occursin("Polyorder_results/CYL/CYL_opt.mat", source)
    @test occursin("Polyorder_results/BCC/BCC_opt.mat", source)
    @test occursin("Polyorder_results/gyroid/GYR_opt.mat", source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_reference_comparison.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_cell_scale_stationarity.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_unrestricted_relaxation.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_stress_free_relaxation.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_bcc_stress_free_convergence_ladder.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_gyr_stress_free_convergence_ladder.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_logit_lbfgs_probe.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_gyr_logit_cell_scan.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_kkt_gate.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_kkt_tolerance_sensitivity.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_energy_decomposition.jl",
        source)
    @test occursin("write_burp_bvk1_polyorder_nonlamellar_manuscript_panel.jl",
        source)
    @test occursin("write_burp_bvk1_representative_phase_rank_gate.jl",
        source)
    @test occursin("write_burp_bvk1_phase_rank_competitor_matrix.jl",
        source)
    @test occursin("write_burp_bvk1_ud_phase_space_coverage.jl", source)
    @test occursin("write_burp_bvk1_publication_readiness_matrix.jl",
        source)
    @test occursin("write_burp_bvk1_ud_reference_coverage.jl", source)
    @test occursin("write_burp_bvk1_publication_packet.jl", source)
    @test occursin("--dry-run", source)
    @test occursin("--from-step", source)
    @test occursin("--to-step", source)
    @test occursin("--skip-step", source)
    @test occursin("skipped_existing", source)
    @test occursin("--continue-on-error", source)
    @test occursin("nonlamellar_publication_pipeline.csv", source)
    @test occursin("nonlamellar_publication_pipeline.md", source)
    @test occursin("nonlamellar_publication_pipeline.log", source)
    @test occursin("DIS/LAM/CYL/BCC/GYR", source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("nonlamellar_publication_pipeline", packet_source)
end

@testset "BURP-BVK1 representative phase-rank gate source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_representative_phase_rank_gate.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("burp_bvk1_representative_phase_rank_gate.csv", source)
    @test occursin("burp_bvk1_representative_phase_rank_gate_summary.csv",
        source)
    @test occursin("burp_bvk1_representative_phase_rank_gate.svg", source)
    @test occursin("burp_bvk1_representative_phase_rank_gate.md", source)
    @test occursin("nonlamellar_manuscript_panel.csv", source)
    @test occursin("same_point_lamellar_competitors.csv", source)
    @test occursin("morphology_phase_readiness_summary.csv", source)
    @test occursin("required_morphologies", source)
    @test occursin("accepted_candidate_count", source)
    @test occursin("phase_rank_pass", source)
    @test occursin("missing_accepted_competitors", source)
    @test occursin("_disordered_baseline_candidate_rows", source)
    @test occursin("homogeneous_disordered_baseline", source)
    @test occursin("energy_proxy_scope", source)
    @test occursin("model_energy_density", source)
    @test occursin("not a phase diagram", source)
    @test occursin("LAM", source)
    @test occursin("DIS", source)
    @test occursin("CYL", source)
    @test occursin("BCC", source)
    @test occursin("GYR", source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("burp_bvk1_representative_phase_rank_gate",
        packet_source)
end

@testset "BURP-BVK1 phase-rank competitor matrix source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_phase_rank_competitor_matrix.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("DFM_SKIP_BURP_BVK1_PHASE_RANK_COMPETITOR_MATRIX_MAIN",
        source)
    @test occursin("burp_bvk1_representative_phase_rank_gate.csv", source)
    @test occursin("burp_bvk1_representative_phase_rank_gate_summary.csv",
        source)
    @test occursin("polyorder_nonlamellar_reference", source)
    @test occursin("required_morphologies", source)
    @test occursin("accepted_competitor", source)
    @test occursin("rejected_competitor", source)
    @test occursin("missing_competitor", source)
    @test occursin("not_required", source)
    @test occursin("phase_rank_pass", source)
    @test occursin("nonlamellar_phase_rank_competitor_matrix.csv", source)
    @test occursin("nonlamellar_phase_rank_competitor_summary.csv", source)
    @test occursin("nonlamellar_phase_rank_competitor_matrix.svg", source)
    @test occursin("nonlamellar_phase_rank_competitor_matrix.md", source)
    @test occursin("DIS/LAM/CYL/BCC/GYR", source)

    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("phase_rank_competitor_matrix", pipeline_source)
    @test occursin("write_burp_bvk1_phase_rank_competitor_matrix.jl",
        pipeline_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("Nonlamellar phase-rank competitor matrix",
        packet_source)
    @test occursin("nonlamellar_phase_rank_competitor_matrix",
        packet_source)

    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("nonlamellar_phase_rank_competitor_matrix",
        readiness_source)
end

@testset "BURP-BVK1 nonlamellar mechanism synthesis source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_nonlamellar_mechanism_synthesis.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("nonlamellar_mechanism_synthesis.csv", source)
    @test occursin("nonlamellar_mechanism_synthesis.svg", source)
    @test occursin("nonlamellar_mechanism_synthesis.md", source)
    @test occursin("nonlamellar_phase_rank_competitor_summary.csv", source)
    @test occursin("nonlamellar_shape_fidelity.csv", source)
    @test occursin("nonlamellar_adaptive_stiffness.csv", source)
    @test occursin("polyorder_nonlamellar_energy_decomposition_summary.csv",
        source)
    @test occursin("finite_amplitude_rank_pass", source)
    @test occursin("endpoint_limited", source)
    @test occursin("near_dis_collapse", source)
    @test occursin("nonlamellar_mechanism_class", source)

    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("nonlamellar_mechanism_synthesis", pipeline_source)
    @test occursin("write_burp_bvk1_nonlamellar_mechanism_synthesis.jl",
        pipeline_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("Nonlamellar mechanism synthesis", packet_source)
    @test occursin("nonlamellar_mechanism_synthesis", packet_source)

    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("nonlamellar_mechanism_synthesis", readiness_source)
end

@testset "BURP-BVK1 publication model comparison synthesis source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_model_comparison_synthesis.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("publication_model_comparison_synthesis.csv", source)
    @test occursin("publication_model_comparison_synthesis_summary.csv",
        source)
    @test occursin("publication_model_comparison_synthesis.svg", source)
    @test occursin("publication_model_comparison_synthesis.md", source)
    @test occursin("aggregate_metrics.csv", source)
    @test occursin("weak_response_summary.csv", source)
    @test occursin("manuscript_readiness_matrix.csv", source)
    @test occursin("nonlamellar_mechanism_synthesis_summary.csv", source)
    @test occursin("ud_parity_scorecard_summary.csv", source)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate_summary.csv",
        source)
    @test occursin("weak_response", source)
    @test occursin("finite_amplitude_lamella", source)
    @test occursin("nonlamellar_transfer", source)
    @test occursin("claim_scope", source)
    @test occursin("limited same-target GYR phase-rank passes", source)
    @test occursin("fA=0.39, chiN=20", source)
    @test occursin("fA=0.41/f0.42, chiN=24", source)
    @test occursin("fA=0.42", source)
    @test occursin("bounded fA=0.39, fA=0.41, and fA=0.42 GYR rank pass",
        source)
    @test occursin("best_current_publishable_role", source)
    @test occursin("UD-first method-equity gate", source)
    @test occursin("BURP-TI/BVK1 CYL/BCC mismatches are numerical-method evidence first",
        source)

    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("publication_model_comparison_synthesis",
        pipeline_source)
    @test occursin("write_burp_bvk1_publication_model_comparison_synthesis.jl",
        pipeline_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("Publication model comparison synthesis", packet_source)
    @test occursin("publication_model_comparison_synthesis", packet_source)
    @test occursin("ud_target_f0p39_chiN20_phase_rank_gate", packet_source)
    @test occursin("burp_bvk1_ud_neighbor_f0p41_chiN24_phase_rank_gate",
        packet_source)
    @test occursin("burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate",
        packet_source)

    readiness_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_readiness_matrix.jl"), String)
    @test occursin("publication_model_comparison_synthesis",
        readiness_source)

    blueprint_path = joinpath(@__DIR__, "..", "docs",
        "burp_bvk1_manuscript_blueprint.md")
    blueprint = isfile(blueprint_path) ? read(blueprint_path, String) : ""
    @test occursin("publication_model_comparison_synthesis", blueprint)
    @test occursin("UD-neighbor f0.42 GYR phase-rank passes", blueprint)
    @test occursin(
        "results/burp_bvk1_ud_neighbor_f0p42_chiN24_phase_rank_gate/phase_rank_gate_summary.csv",
        blueprint)
    @test !occursin("f0.42`, `chi N=24` target or a second independent",
        blueprint)
    @test !occursin("`GYR_ud_neighbor_f0p42_chiN24` remains the queued follow-up",
        blueprint)

    summary_path = joinpath(@__DIR__, "..", "results",
        "burp_bvk1_publication_model_comparison_synthesis",
        "publication_model_comparison_synthesis_summary.csv")
    @test isfile(summary_path)
    summary = isfile(summary_path) ? read(summary_path, String) : ""
    @test occursin("fA=0.41/f0.42, chiN=24", summary)
    @test occursin("f0.42", summary)
    @test occursin("phase-rank passes 6/6", summary)
end

@testset "BURP-BVK1 same-point lamellar competitors source contract" begin
    script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_same_point_lamellar_competitors.jl")
    @test isfile(script_path)
    source = isfile(script_path) ? read(script_path, String) : ""
    @test occursin("same_point_lamellar_competitors.csv", source)
    @test occursin("same_point_lamellar_competitors.svg", source)
    @test occursin("same_point_lamellar_competitors.md", source)
    @test occursin("minimize_diblock_burp_lamella_stress_free", source)
    @test occursin("minimize_diblock_bvk1_lamella_stress_free", source)
    @test occursin("accepted_lamellar_competitor", source)
    @test occursin("period_local_minimum_check_pass", source)
    @test occursin("model_converged", source)
    @test occursin("Polyorder_results/CYL/CYL_opt.mat", source)
    @test occursin("Polyorder_results/BCC/BCC_opt.mat", source)
    @test occursin("Polyorder_results/gyroid/GYR_opt.mat", source)
    @test occursin("CYL", source)
    @test occursin("BCC", source)
    @test occursin("GYR", source)
    @test occursin("0.30", source)
    @test occursin("0.15", source)
    @test occursin("0.36", source)
    @test occursin("30.0", source)
    @test occursin("16.0", source)
    @test occursin("_target_scope_phrase", source)
    @test occursin("target-specific", source)
    @test count("flush(stdout)", source) >= 2

    pipeline_source = read(joinpath(@__DIR__, "..", "scripts",
        "run_burp_bvk1_nonlamellar_publication_pipeline.jl"), String)
    @test occursin("write_burp_bvk1_same_point_lamellar_competitors.jl",
        pipeline_source)

    packet_source = read(joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl"), String)
    @test occursin("same_point_lamellar_competitors", packet_source)
end

@testset "Diblock morphology density initialization" begin
    field = DFMMonteCarlo._diblock_morphology_seed_field(:GYR, (4, 4, 4);
        amplitude=0.8)
    phi = DFMMonteCarlo._diblock_density_from_seed_field(field, 0.37)
    seed = DFMMonteCarlo._diblock_seed_field_from_density(phi, 0.37)
    roundtrip = DFMMonteCarlo._diblock_density_from_seed_field(seed, 0.37)
    @test size(seed) == size(phi)
    @test mean(roundtrip) ≈ 0.37 atol=1.0e-12
    @test maximum(abs.(roundtrip .- phi)) < 1.0e-10

    result = DFMMonteCarlo.relax_diblock_morphology_seed_nd(:GYR;
        model=:BVK1, f=0.37, chiN=16.0, dims=size(phi), lengths=(1.0, 1.0, 1.0),
        initial_density=phi, max_iterations=0)
    @test result.initial_source == "density_logit"
    @test result.accepted_steps == 0
    @test maximum(abs.(result.phi .- phi)) < 1.0e-10
end

@testset "Periodic binary topology diagnostics" begin
    empty_mask = falses(4, 4, 4)
    full_mask = trues(4, 4, 4)
    single_cube = falses(4, 4, 4)
    single_cube[2, 2, 2] = true
    two_cubes = falses(4, 4, 4)
    two_cubes[1, 1, 1] = true
    two_cubes[3, 3, 3] = true

    @test DFMMonteCarlo._periodic_binary_cubical_euler(empty_mask).euler_characteristic == 0
    @test DFMMonteCarlo._periodic_binary_cubical_euler(full_mask).euler_characteristic == 0
    @test DFMMonteCarlo._periodic_binary_cubical_euler(single_cube).euler_characteristic == 1
    @test DFMMonteCarlo._periodic_binary_cubical_euler(two_cubes).euler_characteristic == 2

    empty_components = DFMMonteCarlo._periodic_binary_component_summary(empty_mask)
    full_components = DFMMonteCarlo._periodic_binary_component_summary(full_mask)
    single_components = DFMMonteCarlo._periodic_binary_component_summary(single_cube)
    @test empty_components.component_count == 0
    @test full_components.component_count == 1
    @test full_components.max_wrap_dimensions == 3
    @test single_components.component_count == 1
    @test single_components.max_wrap_dimensions == 0

    lamellar_field = reshape([i <= 2 ? 1.0 : 0.0 for i in 1:4, j in 1:4, k in 1:4],
        4, 4, 4)
    lamellar = DFMMonteCarlo._periodic_binary_topology_metrics(lamellar_field;
        threshold=0.5)
    @test lamellar.topology_label == "lamellar_or_slab"
    @test !lamellar.negative_euler_double_network
end

@testset "Shell-invariant gyroid support diagnostics" begin
    dims = (10, 10, 10)
    gyr = DFMMonteCarlo._shell_invariant_gyroid_support_metrics(
        DFMMonteCarlo._diblock_morphology_seed_field(:GYR, dims; amplitude=1.0))
    bcc = DFMMonteCarlo._shell_invariant_gyroid_support_metrics(
        DFMMonteCarlo._diblock_morphology_seed_field(:BCC, dims; amplitude=1.0))
    lam = DFMMonteCarlo._shell_invariant_gyroid_support_metrics(
        DFMMonteCarlo._diblock_morphology_seed_field(:LAM, dims; amplitude=1.0))
    cyl = DFMMonteCarlo._shell_invariant_gyroid_support_metrics(
        DFMMonteCarlo._diblock_morphology_seed_field(:CYL, dims; amplitude=1.0))

    @test gyr.shell_invariant_support_label == "gyroid_or_bcc_double_network_support"
    @test gyr.dominant_shell_effective_mode_count == 12
    @test gyr.dominant_shell_participation_ratio ≈ 12.0
    @test gyr.dominant_shell_power_fraction ≈ 1.0
    @test gyr.negative_euler_double_network

    @test bcc.shell_invariant_support_label == "gyroid_or_bcc_double_network_support"
    @test bcc.dominant_shell_effective_mode_count == 12
    @test bcc.negative_euler_double_network
    @test bcc.phase_information_required

    @test lam.shell_invariant_support_label != "gyroid_or_bcc_double_network_support"
    @test lam.dominant_shell_effective_mode_count == 2
    @test !lam.negative_euler_double_network

    @test cyl.shell_invariant_support_label != "gyroid_or_bcc_double_network_support"
    @test cyl.dominant_shell_effective_mode_count == 4
    @test !cyl.negative_euler_double_network
end

@testset "resumable blocked-coordinate 3D VK1 chunk runner" begin
    outdir = mktempdir()
    first = DFMMonteCarlo.run_slit_3d_vk1_blocked_chunk!(outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        local_potential_coeffs=(0.10, -0.05, 0.02),
        target_xi_eff=0.24, profile_measure_power=1.0,
        lateral_fraction_amplitude=5.0e-3,
        profile_transfer_amplitude=5.0e-4,
        profile_mode_amplitude=5.0e-4,
        reference_width_amplitude=2.0e-4,
        lateral_fraction_move_weight=0.50,
        profile_transfer_move_weight=0.20,
        profile_mode_move_weight=0.10,
        reference_width_move_weight=0.20,
        reference_width_reflection_move_weight=0.05,
        profile_mode_count=3,
        chunk_steps=60, first_burnin_steps=20, sample_stride=20,
        chains=2, chain_indices=1:2, parallel_chains=true,
        reference_width_xi=0.24,
        seed=1701, initial_profile=:reference_width)
    second = DFMMonteCarlo.run_slit_3d_vk1_blocked_chunk!(outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        local_potential_coeffs=(0.10, -0.05, 0.02),
        target_xi_eff=0.24, profile_measure_power=1.0,
        lateral_fraction_amplitude=5.0e-3,
        profile_transfer_amplitude=5.0e-4,
        profile_mode_amplitude=5.0e-4,
        reference_width_amplitude=2.0e-4,
        lateral_fraction_move_weight=0.50,
        profile_transfer_move_weight=0.20,
        profile_mode_move_weight=0.10,
        reference_width_move_weight=0.20,
        reference_width_reflection_move_weight=0.05,
        profile_mode_count=3,
        chunk_steps=60, first_burnin_steps=20, sample_stride=20,
        chains=2, chain_indices=1:2, parallel_chains=true,
        reference_width_xi=0.24,
        seed=1701, initial_profile=:reference_width)

    @test first.chunk_count == 1
    @test second.chunk_count == 2
    @test second.params.local_potential_coeffs == (0.10, -0.05, 0.02)
    @test second.sample_count > first.sample_count
    @test second.chain_count == 2
    @test isfile(second.summary_csv)
    @test isfile(second.profile_csv)
    @test isfile(second.diagnostics_csv)
    @test isfile(joinpath(outdir, "slit_3d_mc_observable_convergence.csv"))
    @test isfile(second.profile_svg)
    @test isfile(second.metadata_path)
    @test isfile(second.chunk_manifest_csv)
    @test all(isfile(path) for path in second.final_state_paths)

    manifest_header, manifest_rows =
        DFMMonteCarlo._read_csv_strings(second.chunk_manifest_csv)
    @test "chain" in manifest_header
    @test "chunk" in manifest_header
    @test "trace_path" in manifest_header
    @test length(manifest_rows) == 4

    diagnostic_header, diagnostic_rows =
        DFMMonteCarlo._read_csv_strings(second.diagnostics_csv)
    @test "chunked_sampler" in diagnostic_header
    @test "state_space" in diagnostic_header
    @test "lateral_fraction_acceptance_rate" in diagnostic_header
    @test "reference_width_reflection_acceptance_rate" in diagnostic_header
    @test "xi_block_drift_pass" in diagnostic_header
    @test "observable_chain_xi_sd" in diagnostic_header
    @test "observable_leave_one_chain_xi_min" in diagnostic_header
    @test "observable_central_rhat_max" in diagnostic_header
    @test "observable_calibration_pass" in diagnostic_header
    @test length(diagnostic_rows) == 1
    state_space_idx = findfirst(==("state_space"), diagnostic_header)
    chunk_count_idx = findfirst(==("chunk_count"), diagnostic_header)
    @test diagnostic_rows[1][state_space_idx] == "plane_mean_lateral_fraction"
    @test parse(Float64, diagnostic_rows[1][chunk_count_idx]) == 2.0

    metadata = TOML.parsefile(second.metadata_path)
    @test metadata["run"]["mode"] == "slit_3d_hevk1_blocked_profile_fraction_chunked"
    @test metadata["sampling"]["state_space"] == "plane_mean_lateral_fraction"
    @test metadata["sampling"]["chunk_count"] == 2
    @test metadata["sampling"]["chain_count"] == 2
    @test metadata["sampling"]["parallel_chains"] == true
    @test metadata["sampling"]["phi_min"] == 1.0e-10
    @test metadata["sampling"]["reference_width_reflection_move_weight"] == 0.05
    @test metadata["parameters"]["local_potential_a1"] == 0.10
    @test metadata["parameters"]["local_potential_a2"] == -0.05
    @test metadata["parameters"]["local_potential_a3"] == 0.02

    replica_outdir = mktempdir()
    replica_chunk = DFMMonteCarlo.run_slit_3d_vk1_blocked_chunk!(replica_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        local_potential_coeffs=(0.10, -0.05, 0.02),
        target_xi_eff=0.24, profile_measure_power=1.0,
        lateral_fraction_amplitude=5.0e-3,
        profile_transfer_amplitude=5.0e-4,
        profile_mode_amplitude=5.0e-4,
        reference_width_amplitude=2.0e-4,
        lateral_fraction_move_weight=0.50,
        profile_transfer_move_weight=0.20,
        profile_mode_move_weight=0.10,
        reference_width_move_weight=0.20,
        chunk_steps=80, first_burnin_steps=20, sample_stride=20,
        chains=1, chain_indices=1:1,
        reference_width_xi=0.24,
        seed=1801, initial_profile=:reference_width,
        replica_temperatures=[1.0, 1.4],
        replica_swap_interval=8)
    @test replica_chunk.sample_count > 0
    replica_header, replica_rows =
        DFMMonteCarlo._read_csv_strings(replica_chunk.diagnostics_csv)
    @test "replica_count" in replica_header
    @test "replica_swap_attempts" in replica_header
    @test parse(Int, replica_rows[1][findfirst(==("replica_count"),
        replica_header)]) == 2
    @test parse(Int, replica_rows[1][findfirst(==("replica_swap_attempts"),
        replica_header)]) > 0
    replica_metadata = TOML.parsefile(replica_chunk.metadata_path)
    @test replica_metadata["run"]["sampler"] ==
          "arianna_blocked_profile_fraction_replica_exchange_chunked"
    @test replica_metadata["sampling"]["replica_count"] == 2
end

@testset "resumable 3D measure-corrected chunk runner" begin
    outdir = mktempdir()
    first = DFMMonteCarlo.run_slit_3d_measure_corrected_chunk!(outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        local_potential_coeffs=(0.10, -0.05, 0.02),
        target_xi_eff=0.24, profile_measure_power=1.0,
        profile_mode_amplitude=1.0e-4, profile_mode_count=3,
        chunk_steps=60, first_burnin_steps=20, sample_stride=20,
        chains=2, chain_indices=1:2,
        initial_pair_amplitude=1.0e-3, initial_plane_amplitude=5.0e-4,
        initial_global_amplitude=2.0e-4, profile_mode_move_weight=0.25,
        initial_reference_width_amplitude=2.0e-4,
        reference_width_move_weight=0.10, reference_width_xi=0.24,
        tune_rounds=0, seed=701)
    second = DFMMonteCarlo.run_slit_3d_measure_corrected_chunk!(outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        local_potential_coeffs=(0.10, -0.05, 0.02),
        target_xi_eff=0.24, profile_measure_power=1.0,
        profile_mode_amplitude=1.0e-4, profile_mode_count=3,
        chunk_steps=60, first_burnin_steps=20, sample_stride=20,
        chains=2, chain_indices=1:2,
        initial_pair_amplitude=1.0e-3, initial_plane_amplitude=5.0e-4,
        initial_global_amplitude=2.0e-4, profile_mode_move_weight=0.25,
        initial_reference_width_amplitude=2.0e-4,
        reference_width_move_weight=0.10, reference_width_xi=0.24,
        tune_rounds=0, seed=701)

    @test first.chunk_count == 1
    @test second.chunk_count == 2
    @test second.params.local_potential_coeffs == (0.10, -0.05, 0.02)
    @test second.sample_count > first.sample_count
    @test second.chain_count == 2
    @test isfile(second.profile_csv)
    @test isfile(second.diagnostics_csv)
    @test isfile(second.profile_svg)
    @test isfile(second.metadata_path)
    @test isfile(second.chunk_manifest_csv)
    @test all(isfile(path) for path in second.final_state_paths)

    manifest_header, manifest_rows =
        DFMMonteCarlo._read_csv_strings(second.chunk_manifest_csv)
    @test "chain" in manifest_header
    @test "chunk" in manifest_header
    @test "trace_path" in manifest_header
    @test length(manifest_rows) == 4

    diagnostic_header, diagnostic_rows =
        DFMMonteCarlo._read_csv_strings(second.diagnostics_csv)
    @test "chunked_sampler" in diagnostic_header
    @test "chunk_count" in diagnostic_header
    @test "xi_block_drift_pass" in diagnostic_header
    @test length(diagnostic_rows) == 1
    chunked_idx = findfirst(==("chunked_sampler"), diagnostic_header)
    chunk_count_idx = findfirst(==("chunk_count"), diagnostic_header)
    @test diagnostic_rows[1][chunked_idx] == "true"
    @test parse(Float64, diagnostic_rows[1][chunk_count_idx]) == 2.0

    metadata = TOML.parsefile(second.metadata_path)
    @test metadata["run"]["mode"] == "slit_3d_measure_corrected_chunked"
    @test metadata["parameters"]["local_potential_a1"] == 0.10
    @test metadata["parameters"]["local_potential_a2"] == -0.05
    @test metadata["parameters"]["local_potential_a3"] == 0.02
    @test metadata["sampling"]["chunk_count"] == 2
    @test metadata["sampling"]["chain_count"] == 2

    tail = DFMMonteCarlo.run_slit_3d_measure_corrected_chunk!(outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        local_potential_coeffs=(0.10, -0.05, 0.02),
        target_xi_eff=0.24, profile_measure_power=1.0,
        profile_mode_amplitude=1.0e-4, profile_mode_count=3,
        chunk_steps=60, first_burnin_steps=20, sample_stride=20,
        chains=2, chain_indices=1:2, analysis_tail_chunks=1,
        initial_pair_amplitude=1.0e-3, initial_plane_amplitude=5.0e-4,
        initial_global_amplitude=2.0e-4, profile_mode_move_weight=0.25,
        initial_reference_width_amplitude=2.0e-4,
        reference_width_move_weight=0.10, reference_width_xi=0.24,
        tune_rounds=0, seed=701)
    @test tail.chunk_count == 1
    @test tail.completed_chunk_count == 3
    @test tail.sample_count < second.sample_count
    tail_header, tail_rows = DFMMonteCarlo._read_csv_strings(tail.diagnostics_csv)
    @test "analysis_tail_chunks" in tail_header
    @test "completed_chunk_count" in tail_header
    tail_chunk_idx = findfirst(==("chunk_count"), tail_header)
    completed_chunk_idx = findfirst(==("completed_chunk_count"), tail_header)
    @test parse(Float64, tail_rows[1][tail_chunk_idx]) == 1.0
    @test parse(Float64, tail_rows[1][completed_chunk_idx]) == 3.0
    tail_metadata = TOML.parsefile(tail.metadata_path)
    @test tail_metadata["sampling"]["chunk_count"] == 1
    @test tail_metadata["sampling"]["completed_chunk_count"] == 3
    @test tail_metadata["sampling"]["analysis_tail_chunks"] == 1

    continuation_outdir = mktempdir()
    continuation = DFMMonteCarlo.run_slit_3d_measure_corrected_chunk!(continuation_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        target_xi_eff=0.24, profile_measure_power=1.001,
        profile_mode_amplitude=1.0e-4, profile_mode_count=3,
        chunk_steps=60, first_burnin_steps=20, sample_stride=20,
        chains=2, chain_indices=1:2, initial_state_outdir=outdir,
        initial_pair_amplitude=1.0e-3, initial_plane_amplitude=5.0e-4,
        initial_global_amplitude=2.0e-4, profile_mode_move_weight=0.25,
        initial_reference_width_amplitude=2.0e-4,
        reference_width_move_weight=0.10, reference_width_xi=0.24,
        tune_rounds=0, seed=801)
    @test continuation.chunk_count == 1
    @test continuation.completed_chunk_count == 1
    @test isfile(continuation.metadata_path)
    continuation_metadata = TOML.parsefile(continuation.metadata_path)
    @test continuation_metadata["sampling"]["initial_state_outdir"] == outdir
    @test continuation_metadata["sampling"]["initial_state_mode"] == "external_chunk_state"
end

@testset "Debye and perturbation references" begin
    @test debye(0.0) == 1.0
    @test isapprox(debye(1.0), 2 * (exp(-1.0) + 1.0 - 1.0), rtol=1e-14)
    @test isapprox(dfm_pressure(300.0, 1.0, 20.0), 212.03120808345622, rtol=1e-10)
    @test isapprox(dfm_pressure_smeared(300.0, 1.0, 20.0; alpha=0.1, n=1001, kmax=100.0),
        217.75163315442185, rtol=1e-8)
    coarse_grid_pressure = dfm_pressure_discrete(300.0, 1.0, 30.0, Grid3D(8, 6.4); alpha=0.1)
    continuum_pressure = dfm_pressure_smeared(300.0, 1.0, 30.0; alpha=0.1)
    fine_grid_pressure = dfm_pressure_discrete(300.0, 1.0, 30.0, Grid3D(64, 6.4); alpha=0.1)
    mean_field = DFMMonteCarlo.mean_field_pressure(1.0, 30.0)
    @test abs(coarse_grid_pressure - mean_field) < abs(continuum_pressure - mean_field)
    @test isapprox(fine_grid_pressure, continuum_pressure; atol=0.1)
end

@testset "block-matrix URP first deterministic gates" begin
    @test isapprox(block_urp_eta(1.2, 1.0), log(1.2), rtol=1e-14)
    @test isapprox(block_urp_eta(0.7, 0.5), 0.5 * log(0.7 / 0.5), rtol=1e-14)

    f = 0.37
    eps = 1.0e-6
    central_slope = (diblock_burp_eta_psi(f + eps, f) -
                     diblock_burp_eta_psi(f - eps, f)) / (2eps)
    @test isapprox(central_slope, 1.0; rtol=1e-9, atol=1e-9)

    k2 = 3.0
    gamma = diblock_ideal_vertex(k2; f=f, N=1.0, b=1.0)
    v = [1.0, -1.0]
    @test isapprox(diblock_composition_kernel(k2; f=f), dot(v, gamma * v); rtol=1e-13)

    small_k = diblock_composition_kernel(1.0e-5; f=0.5)
    mid_k = diblock_composition_kernel(10.0; f=0.5)
    large_k = diblock_composition_kernel(1.0e5; f=0.5)
    scan_k2 = 10 .^ range(-4, 4; length=241)
    scan_values = [diblock_composition_kernel(candidate_k2; f=0.5) for candidate_k2 in scan_k2]
    min_index = argmin(scan_values)
    @test 1 < min_index < length(scan_values)
    @test scan_values[min_index] < small_k
    @test scan_values[min_index] < large_k
    @test mid_k < small_k
    @test mid_k < large_k

    lamella = minimize_diblock_burp_lamella(; f=0.5, chiN=12.0, nx=64,
        mode_count=3, nquad=8, max_iterations=5_000,
        relaxation_step=1.0e-2, tolerance=1.0e-7)
    @test lamella.converged
    @test lamella.minimum_phi > 0.0
    @test lamella.maximum_phi < 1.0
    @test isapprox(lamella.mean_phi, 0.5; atol=1.0e-12)
    @test lamella.maximum_phi - lamella.minimum_phi > 0.1
    @test lamella.energy < lamella.homogeneous_energy - 1.0e-4
    @test isapprox(
        diblock_burp_ti_energy_1d(lamella.phi_a; f=0.5, chiN=12.0, L=lamella.L, nquad=8),
        lamella.energy; rtol=1.0e-12)

    stress_free_lamella = minimize_diblock_burp_lamella_stress_free(; f=0.5,
        chiN=12.0, nx=32, mode_count=2, nquad=4, lower_factor=0.9,
        upper_factor=1.1, initial_amplitudes=[0.2], max_iterations=1200,
        max_period_iterations=5)
    @test stress_free_lamella.optimizer == "Brent"
    @test 0.9 < stress_free_lamella.period_factor < 1.1
    @test isfinite(stress_free_lamella.objective)
    @test stress_free_lamella.result.converged
    @test stress_free_lamella.result.minimum_phi > 0.0
    @test stress_free_lamella.result.maximum_phi < 1.0

    _, weak_phi = diblock_lamella_profile_from_modes([1.0e-3]; f=0.5, nx=64, L=lamella.L)
    @test isapprox(
        diblock_rpa_quadratic_energy_1d(weak_phi; f=0.5, chiN=12.0, L=lamella.L),
        diblock_burp_ti_energy_1d(weak_phi; f=0.5, chiN=12.0, L=lamella.L, nquad=8);
        rtol=1.0e-3, atol=1.0e-10)

    ud_coeffs = diblock_uneyama_doi_coefficients(; f=0.5, N=1.0, b=1.0)
    @test ud_coeffs.A ≈ ud_coeffs.A'
    @test ud_coeffs.C ≈ ud_coeffs.C'
    @test ud_coeffs.C[1, 2] ≈ -1.0
    @test ud_coeffs.s_f > 0.0
    @test ud_coeffs.s_one_minus_f > 0.0

    bvk1_coeffs = diblock_bvk1_coefficients(; f=0.35, N=1.0, b=1.0)
    @test isapprox(bvk1_coeffs.A_psi, 9.0 / (0.35^2 * 0.65^2); rtol=1.0e-12)
    @test isapprox(bvk1_coeffs.M_psi / (0.35 * 0.65), bvk1_coeffs.C_psi; rtol=1.0e-12)
    @test bvk1_coeffs.C_psi > 0.0
    @test bvk1_coeffs.M_psi > 0.0
    @test bvk1_coeffs.K_psi0 ≈ 1.0 / 24.0

    weak_nx = 512
    weak_L = 2.0
    weak_f = 0.5
    weak_chiN = 12.0
    weak_x = [(idx - 0.5) * weak_L / weak_nx for idx in 1:weak_nx]
    weak_amp = 1.0e-4
    weak_phi = weak_f .+ weak_amp .* cos.(2.0 * pi .* weak_x ./ weak_L)
    weak_coeffs = diblock_bvk1_coefficients(; f=weak_f, N=1.0, b=1.0)
    weak_ud_coeffs = diblock_uneyama_doi_coefficients(
        ; f=weak_f, N=1.0, b=1.0)
    weak_k2 = (2.0 * pi / weak_L)^2
    weak_kernel = weak_coeffs.A_psi / weak_k2 + weak_coeffs.C_psi +
        weak_coeffs.b^2 * weak_k2 / (12.0 * weak_f * (1.0 - weak_f))
    weak_ud_kernel = weak_coeffs.A_psi / weak_k2 +
        weak_ud_coeffs.C[1, 1] + weak_ud_coeffs.C[2, 2] -
        weak_ud_coeffs.C[1, 2] / (weak_f * (1.0 - weak_f)) +
        weak_coeffs.b^2 * weak_k2 / (12.0 * weak_f * (1.0 - weak_f))
    @test isapprox(
        diblock_uneyama_doi_composition_kernel(weak_k2; f=weak_f),
        weak_ud_kernel; rtol=1.0e-13)
    @test isapprox(
        diblock_bvk1_gaussian_kernel(weak_k2; f=weak_f),
        weak_kernel; rtol=1.0e-13)
    @test diblock_bvk1_gaussian_kernel(weak_k2; f=weak_f) !=
          diblock_uneyama_doi_composition_kernel(weak_k2; f=weak_f)
    @test isapprox(
        diblock_burp_ti_gaussian_kernel(weak_k2; f=weak_f),
        diblock_composition_kernel(weak_k2; f=weak_f); rtol=1.0e-13)
    @test diblock_uneyama_doi_composition_kernel(weak_k2; f=weak_f) !=
          diblock_composition_kernel(weak_k2; f=weak_f)
    weak_hat = fft(weak_phi .- weak_f)
    weak_dx = weak_L / weak_nx
    weak_expected = 0.5 * weak_dx / weak_nx * weak_kernel *
        sum(abs2, weak_hat[2:end]) +
        weak_chiN * weak_dx * sum(value * (1.0 - value) - weak_f * (1.0 - weak_f)
            for value in weak_phi)
    @test isapprox(
        diblock_bvk1_energy_1d(weak_phi; f=weak_f, chiN=weak_chiN, L=weak_L,
            adaptive=false),
        weak_expected; rtol=2.0e-4, atol=1.0e-12)
    @test isapprox(
        diblock_bvk1_energy_1d(weak_phi; f=weak_f, chiN=weak_chiN, L=weak_L),
        weak_expected; rtol=2.0e-4, atol=1.0e-12)
    weak_kpsi = diblock_bvk1_kpsi_profile_1d(weak_phi; f=weak_f, L=weak_L)
    @test length(weak_kpsi) == weak_nx
    @test minimum(weak_kpsi) >= weak_coeffs.K_psi0
    @test maximum(weak_kpsi) <= 1.5 * weak_coeffs.K_psi0 + 1.0e-12
    uniform_kpsi = diblock_bvk1_kpsi_profile_1d(fill(weak_f, weak_nx);
        f=weak_f, L=weak_L)
    @test all(isapprox(value, weak_coeffs.K_psi0; atol=1.0e-14)
        for value in uniform_kpsi)

    bvk1_lamella = minimize_diblock_bvk1_lamella(; f=0.5, chiN=12.0, nx=64,
        mode_count=3, max_iterations=3_000)
    @test bvk1_lamella.converged
    @test bvk1_lamella.minimum_phi > 0.0
    @test bvk1_lamella.maximum_phi < 1.0
    @test isapprox(bvk1_lamella.mean_phi, 0.5; atol=1.0e-12)
    @test bvk1_lamella.maximum_phi - bvk1_lamella.minimum_phi > 0.1
    @test bvk1_lamella.energy < bvk1_lamella.homogeneous_energy - 1.0e-4

    stress_free_bvk1 = minimize_diblock_bvk1_lamella_stress_free(; f=0.5,
        chiN=12.0, nx=32, mode_count=2, lower_factor=0.9, upper_factor=1.1,
        initial_amplitudes=[0.2], max_iterations=1200, max_period_iterations=5)
    @test stress_free_bvk1.model == "BVK1"
    @test stress_free_bvk1.optimizer == "Brent"
    @test 0.9 < stress_free_bvk1.period_factor < 1.1
    @test isfinite(stress_free_bvk1.objective)
    @test stress_free_bvk1.result.converged
    @test stress_free_bvk1.result.minimum_phi > 0.0
    @test stress_free_bvk1.result.maximum_phi < 1.0
    @test stress_free_bvk1.result.maximum_phi - stress_free_bvk1.result.minimum_phi > 0.1
    @test stress_free_bvk1.result.energy <
        stress_free_bvk1.result.homogeneous_energy - 1.0e-4

    bootstrapped_bvk1 = minimize_diblock_bvk1_lamella_stress_free(; f=0.5,
        chiN=12.0, nx=32, mode_count=2, lower_factor=0.9, upper_factor=1.2,
        bootstrap_period_factors=[1.05], bootstrap_window=0.08,
        initial_amplitudes=[0.2], max_iterations=1200, max_period_iterations=5,
        local_check_fraction=0.02, period_strategy=:direct_bootstrap)
    @test bootstrapped_bvk1.bootstrap_used
    @test isapprox(bootstrapped_bvk1.bootstrap_period_factor, 1.05; atol=1.0e-12)
    @test isfinite(bootstrapped_bvk1.local_left_objective)
    @test isfinite(bootstrapped_bvk1.local_right_objective)
    @test bootstrapped_bvk1.local_left_factor < bootstrapped_bvk1.period_factor <
          bootstrapped_bvk1.local_right_factor
    @test isapprox(bootstrapped_bvk1.period_b,
        bootstrapped_bvk1.reference_period_b * bootstrapped_bvk1.period_factor;
        rtol=1.0e-12)
    @test isapprox(bootstrapped_bvk1.result.L, bootstrapped_bvk1.period_b;
        rtol=1.0e-12)
    @test !bootstrapped_bvk1.boundary_limited
    @test bootstrapped_bvk1.optimizer == "DirectBootstrap"
    @test isapprox(bootstrapped_bvk1.period_factor, 1.05; atol=1.0e-12)
    @test bootstrapped_bvk1.local_minimum_check_pass ==
          (bootstrapped_bvk1.local_left_objective > bootstrapped_bvk1.objective &&
           bootstrapped_bvk1.local_right_objective > bootstrapped_bvk1.objective)
    @test isapprox(bootstrapped_bvk1.local_center_factor,
        bootstrapped_bvk1.period_factor; atol=1.0e-12)
    @test isapprox(bootstrapped_bvk1.local_center_objective,
        bootstrapped_bvk1.objective; atol=1.0e-12)

    locally_corrected_bvk1 = minimize_diblock_bvk1_lamella_stress_free(; f=0.5,
        chiN=12.0, nx=32, mode_count=2, lower_factor=0.9, upper_factor=1.2,
        bootstrap_period_factors=[1.0], bootstrap_window=0.08,
        initial_amplitudes=[0.2], max_iterations=1200, max_period_iterations=8,
        local_check_fraction=0.02, period_strategy=:bootstrap_local)
    @test locally_corrected_bvk1.optimizer == "BootstrapLocal"
    @test locally_corrected_bvk1.period_factor > 1.0
    @test isfinite(locally_corrected_bvk1.local_left_objective)
    @test isfinite(locally_corrected_bvk1.local_right_objective)
    @test isfinite(locally_corrected_bvk1.local_center_objective)

    _, embedded_phi_1d = diblock_lamella_profile_from_modes([0.7, -0.12];
        f=0.5, nx=32, L=3.2)
    embedded_phi_2d = repeat(reshape(embedded_phi_1d, :, 1), 1, 4)
    embedded_phi_3d = repeat(reshape(embedded_phi_1d, :, 1, 1), 1, 3, 2)
    burp_1d = diblock_burp_ti_energy_1d(embedded_phi_1d; f=0.5, chiN=18.0,
        L=3.2, nquad=8)
    @test isapprox(diblock_burp_ti_energy_nd(embedded_phi_2d; f=0.5,
        chiN=18.0, lengths=(3.2, 1.0), nquad=8), burp_1d; rtol=1.0e-12,
        atol=1.0e-12)
    @test isapprox(diblock_burp_ti_energy_nd(embedded_phi_3d; f=0.5,
        chiN=18.0, lengths=(3.2, 1.0, 1.0), nquad=8), burp_1d;
        rtol=1.0e-12, atol=1.0e-12)
    bvk1_1d = diblock_bvk1_energy_1d(embedded_phi_1d; f=0.5, chiN=18.0,
        L=3.2)
    @test isapprox(diblock_bvk1_energy_nd(embedded_phi_2d; f=0.5,
        chiN=18.0, lengths=(3.2, 1.0)), bvk1_1d; rtol=1.0e-12,
        atol=1.0e-12)
    @test isapprox(diblock_bvk1_energy_nd(embedded_phi_3d; f=0.5,
        chiN=18.0, lengths=(3.2, 1.0, 1.0)), bvk1_1d; rtol=1.0e-12,
        atol=1.0e-12)
    ud_1d = diblock_uneyama_doi_energy_1d(embedded_phi_1d; f=0.5, chiN=18.0,
        L=3.2)
    @test isapprox(diblock_uneyama_doi_energy_nd(embedded_phi_2d; f=0.5,
        chiN=18.0, lengths=(3.2, 1.0)), ud_1d; rtol=1.0e-12,
        atol=1.0e-12)
    @test isapprox(diblock_uneyama_doi_energy_nd(embedded_phi_3d; f=0.5,
        chiN=18.0, lengths=(3.2, 1.0, 1.0)), ud_1d; rtol=1.0e-12,
        atol=1.0e-12)

    _, one_mode_phi = diblock_lamella_profile_from_modes([0.7]; f=0.5,
        nx=32, L=3.2)
    lam_seed = diblock_morphology_seed(:LAM; f=0.5, dims=(32,),
        lengths=(3.2,), amplitude=0.7)
    @test size(lam_seed) == (32,)
    @test maximum(abs.(lam_seed .- one_mode_phi)) < 1.0e-13

    morphology_seed_specs = [
        (:CYL, (18, 18), (3.2, 3.2)),
        (:BCC, (8, 8, 8), (3.2, 3.2, 3.2)),
        (:GYR, (8, 8, 8), (3.2, 3.2, 3.2)),
    ]
    for (kind, dims, lengths) in morphology_seed_specs
        seed_phi = diblock_morphology_seed(kind; f=0.35, dims=dims,
            lengths=lengths, amplitude=0.9)
        @test size(seed_phi) == dims
        @test isapprox(mean(seed_phi), 0.35; atol=1.0e-12)
        @test 0.0 < minimum(seed_phi) < maximum(seed_phi) < 1.0
        @test maximum(seed_phi) - minimum(seed_phi) > 0.05
        @test isfinite(diblock_burp_ti_energy_nd(seed_phi; f=0.35,
            chiN=20.0, lengths=lengths, nquad=4))
        @test isfinite(diblock_bvk1_energy_nd(seed_phi; f=0.35,
            chiN=20.0, lengths=lengths))
    end

    amplitude_envelope = optimize_diblock_morphology_seed_amplitude(:CYL;
        model=:BVK1, f=0.35, chiN=20.0, dims=(12, 12),
        lengths=(3.2, 3.2), lower_amplitude=0.0, upper_amplitude=1.4,
        iterations=12)
    @test amplitude_envelope.model == "BVK1"
    @test amplitude_envelope.morphology == "CYL"
    @test 0.0 <= amplitude_envelope.amplitude <= 1.4
    @test isfinite(amplitude_envelope.energy)
    @test isfinite(amplitude_envelope.homogeneous_energy)
    @test amplitude_envelope.energy <= amplitude_envelope.homogeneous_energy + 1.0e-10
    @test isapprox(amplitude_envelope.mean_phi, 0.35; atol=1.0e-12)
    @test 0.0 < amplitude_envelope.min_phi < amplitude_envelope.max_phi < 1.0

    pilot_relaxation = relax_diblock_morphology_seed_nd(:CYL; model=:BVK1,
        f=0.35, chiN=20.0, dims=(8, 8), lengths=(3.2, 3.2),
        initial_amplitude=0.9, max_iterations=2, gradient_step=0.08,
        finite_difference_step=2.0e-4, nquad=3)
    @test pilot_relaxation.model == "BVK1"
    @test pilot_relaxation.morphology == "CYL"
    @test pilot_relaxation.energy <= pilot_relaxation.initial_energy + 1.0e-10
    @test pilot_relaxation.accepted_steps >= 1
    @test isapprox(pilot_relaxation.mean_phi, 0.35; atol=1.0e-12)
    @test 0.0 < pilot_relaxation.min_phi < pilot_relaxation.max_phi < 1.0
    @test pilot_relaxation.gradient_norm > 0.0

    stationarity_probe = probe_diblock_morphology_relaxation_stationarity(:CYL;
        model=:BVK1, f=0.35, chiN=20.0, dims=(6, 6), lengths=(3.2, 3.2),
        initial_amplitude=0.9, max_iterations=2, gradient_step=0.08,
        finite_difference_step=2.0e-4, nquad=3, probe_count=2,
        perturbation_scale=1.0e-4, seed=17)
    @test stationarity_probe.model == "BVK1"
    @test stationarity_probe.morphology == "CYL"
    @test stationarity_probe.probe_count == 2
    @test isfinite(stationarity_probe.min_forward_delta)
    @test isfinite(stationarity_probe.min_backward_delta)
    @test isfinite(stationarity_probe.min_symmetric_delta)
    @test stationarity_probe.min_symmetric_delta <= stationarity_probe.max_symmetric_delta
    @test stationarity_probe.local_stationarity_pass ==
          (stationarity_probe.min_symmetric_delta >= -stationarity_probe.energy_tolerance)
    @test stationarity_probe.relaxation.energy <= stationarity_probe.relaxation.initial_energy + 1.0e-10

    ud_psi_relaxation = relax_diblock_uneyama_doi_morphology_seed_nd(:CYL;
        f=0.35, chiN=20.0, dims=(8, 8), lengths=(3.2, 3.2),
        initial_amplitude=0.9, max_iterations=3, relaxation_step=1.0e-2,
        tolerance=1.0e-8)
    @test ud_psi_relaxation.model == "Uneyama-Doi"
    @test ud_psi_relaxation.method == "psi_theta_relaxation"
    @test ud_psi_relaxation.morphology == "CYL"
    @test ud_psi_relaxation.energy <= ud_psi_relaxation.initial_energy + 1.0e-10
    @test ud_psi_relaxation.accepted_steps >= 1
    @test haskey(ud_psi_relaxation, :projected_gradient_norm)
    @test 0.0 <= ud_psi_relaxation.projected_gradient_norm <=
          ud_psi_relaxation.gradient_norm + 1.0e-12
    @test isapprox(ud_psi_relaxation.mean_phi, 0.35; atol=1.0e-12)
    @test 0.0 < ud_psi_relaxation.min_phi < ud_psi_relaxation.max_phi < 1.0
    ud_projected_relaxation = relax_diblock_uneyama_doi_morphology_seed_nd(:CYL;
        f=0.5, chiN=30.0, dims=(8, 8), lengths=(3.0, 3.0),
        initial_amplitude=0.5, max_iterations=4, relaxation_step=1.0e-2,
        descent_direction=:projected)
    @test ud_projected_relaxation.descent_direction == "projected"
    @test ud_projected_relaxation.energy <=
          ud_projected_relaxation.initial_energy + 1.0e-10
    @test isfinite(ud_projected_relaxation.projected_gradient_norm)
    ud_projected_polish = relax_diblock_uneyama_doi_morphology_seed_nd(:CYL;
        f=0.5, chiN=30.0, dims=(8, 8), lengths=(3.0, 3.0),
        initial_theta=ud_projected_relaxation.theta, max_iterations=2,
        relaxation_step=1.0e-2, descent_direction=:projected)
    @test ud_projected_polish.initial_source == "theta"
    @test isapprox(ud_projected_polish.initial_energy,
        ud_projected_relaxation.energy; atol=1.0e-10)
    @test ud_projected_polish.energy <= ud_projected_polish.initial_energy + 1.0e-10

    ud_cell_scan = scan_diblock_uneyama_doi_morphology_cell_scales(:CYL;
        f=0.35, chiN=20.0, dims=(6, 6), base_lengths=(3.2, 3.2),
        scale_factors=[0.95, 1.0, 1.05], initial_amplitude=0.9,
        max_iterations=2, relaxation_step=1.0e-2)
    @test ud_cell_scan.model == "Uneyama-Doi"
    @test ud_cell_scan.method == "psi_cell_scale_scan"
    @test ud_cell_scan.morphology == "CYL"
    @test length(ud_cell_scan.evaluations) == 3
    @test ud_cell_scan.selected_factor in [0.95, 1.0, 1.05]
    @test isfinite(ud_cell_scan.selected_energy_density)
    @test ud_cell_scan.selected_energy_density <=
          minimum(row.energy_density for row in ud_cell_scan.evaluations) + 1.0e-12
    @test 0.0 < ud_cell_scan.selected_result.min_phi <
          ud_cell_scan.selected_result.max_phi < 1.0
    ud_adaptive_cell = optimize_diblock_uneyama_doi_morphology_cell_scale(:CYL;
        f=0.35, chiN=20.0, dims=(8, 8), base_lengths=(3.2, 3.2),
        initial_scale_factor=1.0, initial_log_step=0.1,
        max_scale_expansions=4, max_refinement_iterations=1,
        initial_amplitude=0.9, max_iterations=4, relaxation_step=1.0e-2)
    @test ud_adaptive_cell.model == "Uneyama-Doi"
    @test ud_adaptive_cell.method == "psi_adaptive_cell_scale"
    @test ud_adaptive_cell.morphology == "CYL"
    @test length(ud_adaptive_cell.evaluations) >= 3
    @test isfinite(ud_adaptive_cell.selected_energy_density)
    @test ud_adaptive_cell.selected_energy_density <=
          minimum(row.energy_density for row in ud_adaptive_cell.evaluations) + 1.0e-12
    @test !ud_adaptive_cell.boundary_limited
    @test ud_adaptive_cell.local_minimum_check_pass
    @test ud_adaptive_cell.local_left_energy_density >
          ud_adaptive_cell.selected_energy_density
    @test ud_adaptive_cell.local_right_energy_density >
          ud_adaptive_cell.selected_energy_density
    @test 0.0 < ud_adaptive_cell.selected_result.min_phi <
          ud_adaptive_cell.selected_result.max_phi < 1.0
    ud_multistart_cell = optimize_diblock_uneyama_doi_morphology_cell_scale_multistart(:CYL;
        f=0.35, chiN=20.0, dims=(8, 8), base_lengths=(3.2, 3.2),
        initial_scale_factor=1.0, initial_log_step=0.1,
        max_scale_expansions=2, max_refinement_iterations=1,
        initial_amplitudes=[0.5, 0.9], max_iterations=3,
        relaxation_step=1.0e-2)
    @test ud_multistart_cell.model == "Uneyama-Doi"
    @test ud_multistart_cell.method == "psi_adaptive_cell_scale_multistart"
    @test ud_multistart_cell.morphology == "CYL"
    @test ud_multistart_cell.initial_amplitudes == (0.5, 0.9)
    @test length(ud_multistart_cell.amplitude_results) == 2
    @test ud_multistart_cell.selected_amplitude in (0.5, 0.9)
    @test length(ud_multistart_cell.evaluations) >= 6
    @test all(row -> haskey(row, :initial_amplitude), ud_multistart_cell.evaluations)
    @test ud_multistart_cell.selected_energy_density <=
          minimum(row.selected_energy_density
              for row in ud_multistart_cell.amplitude_results) + 1.0e-12
    ud_shape_scan = scan_diblock_uneyama_doi_morphology_cell_shapes(:CYL;
        f=0.35, chiN=20.0, dims=(6, 6), base_lengths=(3.2, 3.2),
        shape_factors=[(1.0, 1.0), (1.1, 1.0 / 1.1), (1.0 / 1.1, 1.1)],
        initial_amplitudes=[0.5], max_iterations=2, relaxation_step=1.0e-2)
    @test ud_shape_scan.model == "Uneyama-Doi"
    @test ud_shape_scan.method == "psi_cell_shape_scan"
    @test ud_shape_scan.morphology == "CYL"
    @test length(ud_shape_scan.evaluations) == 3
    @test ud_shape_scan.shape_count == 3
    @test ud_shape_scan.amplitude_count == 1
    @test any(row -> row.shape_factors != (1.0, 1.0), ud_shape_scan.evaluations)
    @test all(row -> length(row.lengths) == 2, ud_shape_scan.evaluations)
    @test isfinite(ud_shape_scan.selected_energy_density)
    @test ud_shape_scan.selected_energy_density <=
          minimum(row.energy_density for row in ud_shape_scan.evaluations) + 1.0e-12
    @test 0.0 < ud_shape_scan.selected_result.min_phi <
          ud_shape_scan.selected_result.max_phi < 1.0
    ud_shape_opt = optimize_diblock_uneyama_doi_morphology_cell_shape_coordinate(:CYL;
        f=0.35, chiN=20.0, dims=(6, 6), base_lengths=(3.2, 3.2),
        shape_direction=(-1.0, 1.0), initial_amplitudes=[0.5],
        initial_alpha_step=0.08, max_shape_expansions=2,
        max_refinement_iterations=1, max_iterations=2, relaxation_step=1.0e-2)
    @test ud_shape_opt.model == "Uneyama-Doi"
    @test ud_shape_opt.method == "psi_cell_shape_coordinate_optimization"
    @test ud_shape_opt.morphology == "CYL"
    @test ud_shape_opt.initial_amplitudes == (0.5,)
    @test length(ud_shape_opt.evaluations) >= 3
    @test any(row -> abs(row.alpha) > 0.0, ud_shape_opt.evaluations)
    @test isfinite(ud_shape_opt.selected_energy_density)
    @test ud_shape_opt.selected_energy_density <=
          minimum(row.energy_density for row in ud_shape_opt.evaluations) + 1.0e-12
    @test isapprox(prod(ud_shape_opt.selected_shape_factors), 1.0; atol=1.0e-12)
    @test length(ud_shape_opt.selected_lengths) == 2
    @test isfinite(ud_shape_opt.local_left_energy_density) ||
          isfinite(ud_shape_opt.local_right_energy_density)
    @test 0.0 < ud_shape_opt.selected_result.min_phi <
          ud_shape_opt.selected_result.max_phi < 1.0
    ud_multi_shape = optimize_diblock_uneyama_doi_morphology_cell_shape_multicoordinate(
        :BCC; f=0.35, chiN=20.0, dims=(4, 4, 4),
        base_lengths=(3.2, 3.2, 3.2), initial_amplitudes=[0.5],
        coordinate_passes=1, initial_alpha_step=0.06,
        max_shape_expansions=1, max_refinement_iterations=1,
        max_iterations=2, relaxation_step=1.0e-2)
    @test ud_multi_shape.model == "Uneyama-Doi"
    @test ud_multi_shape.method == "psi_cell_shape_multicoordinate_optimization"
    @test ud_multi_shape.morphology == "BCC"
    @test ud_multi_shape.direction_count == 2
    @test ud_multi_shape.coordinate_passes == 1
    @test length(ud_multi_shape.pass_results) == 2
    @test all(row -> length(row.shape_direction) == 3, ud_multi_shape.pass_results)
    @test all(row -> haskey(row, :pass_index) && haskey(row, :direction_index),
        ud_multi_shape.evaluations)
    @test isapprox(prod(ud_multi_shape.total_shape_factors), 1.0; atol=1.0e-12)
    @test isapprox(prod(ud_multi_shape.selected_lengths), prod((3.2, 3.2, 3.2));
        rtol=1.0e-12)
    @test isfinite(ud_multi_shape.selected_energy_density)
    @test haskey(ud_multi_shape, :initial_cell_energy_density)
    @test haskey(ud_multi_shape, :shape_improvement_over_initial_cell)
    @test isapprox(ud_multi_shape.shape_improvement_over_initial_cell,
        ud_multi_shape.initial_cell_energy_density - ud_multi_shape.selected_energy_density;
        atol=1.0e-12)
    @test ud_multi_shape.selected_energy_density <=
          minimum(row.selected_energy_density for row in ud_multi_shape.pass_results) +
          1.0e-12
    @test 0.0 < ud_multi_shape.selected_result.min_phi <
          ud_multi_shape.selected_result.max_phi < 1.0
    ud_lam_projection_regression = relax_diblock_uneyama_doi_morphology_seed_nd(:LAM;
        f=0.35, chiN=20.0, dims=(32,), lengths=(3.36,),
        initial_amplitude=0.9, max_iterations=16, relaxation_step=1.0e-2)
    @test isapprox(ud_lam_projection_regression.mean_phi, 0.35; atol=1.0e-12)

    asymmetric_local_period = minimize_diblock_uneyama_doi_lamella_stress_free(;
        f=0.25, chiN=50.0, nx=48, mode_count=4, lower_factor=0.65,
        upper_factor=2.2, bootstrap_period_factors=[1.0], bootstrap_window=0.18,
        initial_amplitudes=[1.0], max_iterations=8000, max_period_iterations=24,
        local_check_fraction=0.01, period_strategy=:bootstrap_local)
    @test asymmetric_local_period.optimizer == "BootstrapLocal"
    @test asymmetric_local_period.period_factor > 1.0
    @test asymmetric_local_period.local_minimum_check_pass

    residual_converged = DiblockBurpLamellaResult(0.5, 12.0, 1.0, 4, 1,
        collect(range(0.0, 0.75; length=4)), fill(0.5, 4), [0.0],
        0.0, 0.0, false, 10_000, 0.4, 0.6, 0.5, 5.0e-7, 2.0e-6)
    @test DFMMonteCarlo._diblock_period_field_converged(
        residual_converged, 1.0e-6, 1.0e-5)
    @test !DFMMonteCarlo._diblock_period_field_converged(
        residual_converged, 1.0e-7, 1.0e-5)

    fixed_script = read(joinpath(@__DIR__, "..", "scripts",
        "write_diblock_burp_polyorder_comparison.jl"), String)
    stress_free_script = read(joinpath(@__DIR__, "..", "scripts",
        "write_diblock_stress_free_lamella_comparison.jl"), String)
    for script_text in (fixed_script, stress_free_script)
        @test occursin("id=\"bvk1\"", script_text)
        @test occursin("label=\"BVK1\"", script_text)
        @test occursin("source=\"bvk1_aligned\"", script_text)
    end
    augment_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_bvk1_stress_free_rows.jl")
    @test isfile(augment_script_path)
    augment_script = read(augment_script_path, String)
    @test occursin("minimize_diblock_bvk1_lamella_stress_free", augment_script)
    @test occursin("--source-dir", augment_script)
    @test occursin("--outdir", augment_script)
    @test occursin("stress_free_lamella_periods.svg", stress_free_script)
    @test occursin("Polyorder SCFT", stress_free_script)
    @test occursin("_period_x_value", stress_free_script)
    @test occursin(">chiN</text>", stress_free_script)
    @test occursin("Series by fA", stress_free_script)
    @test occursin("--reuse-scft-dir", stress_free_script)
    @test occursin("period_local_minimum_check_pass", stress_free_script)
    @test occursin("period_local_center_objective", stress_free_script)
    @test occursin("progress_label", stress_free_script)
    @test occursin("--period-strategy", stress_free_script)
    @test occursin("DirectBootstrap", read(joinpath(@__DIR__, "..", "src",
        "DFMMonteCarlo.jl"), String))
    @test occursin("BootstrapLocal", read(joinpath(@__DIR__, "..", "src",
        "DFMMonteCarlo.jl"), String))
    @test occursin("Float64[clamp(1.0, lo, hi)]", stress_free_script)
    @test occursin("Float64[clamp(1.0, lo, hi)]", augment_script)
    @test occursin("flush(stdout)", read(joinpath(@__DIR__, "..", "src",
        "DFMMonteCarlo.jl"), String))
    kernel_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_diblock_kernel_comparison.jl")
    @test isfile(kernel_script_path)
    kernel_script = read(kernel_script_path, String)
    @test occursin("diblock_burp_ti_gaussian_kernel", kernel_script)
    @test occursin("diblock_bvk1_gaussian_kernel", kernel_script)
    @test occursin("diblock_kernel_comparison.svg", kernel_script)
    weak_map_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_diblock_weak_response_map.jl")
    @test isfile(weak_map_script_path)
    weak_map_script = read(weak_map_script_path, String)
    @test occursin("diblock_weak_response_map.svg", weak_map_script)
    @test occursin("chiN_spinodal", weak_map_script)
    @test occursin("diblock_burp_ti_gaussian_kernel", weak_map_script)
    @test occursin("diblock_bvk1_gaussian_kernel", weak_map_script)
    finite_amplitude_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_diblock_finite_amplitude_diagnostic.jl")
    @test isfile(finite_amplitude_script_path)
    finite_amplitude_script = read(finite_amplitude_script_path, String)
    @test occursin("finite_amplitude_diagnostic.svg", finite_amplitude_script)
    @test occursin("finite_amplitude_mechanism_summary.csv", finite_amplitude_script)
    @test occursin("production stress-free cases", finite_amplitude_script)
    @test occursin("bvk1_closer_than_ud_count", finite_amplitude_script)
    @test occursin("diblock_bvk1_kpsi_profile_1d", finite_amplitude_script)
    @test occursin("amplitude_boundary", finite_amplitude_script)
    metric_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_diblock_stress_free_metric_summary.jl")
    @test isfile(metric_script_path)
    metric_script = read(metric_script_path, String)
    @test occursin("stress_free_metric_summary.svg", metric_script)
    @test occursin("mean_period_abs_error", metric_script)
    @test occursin("amplitude_error", metric_script)
    bracket_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_diblock_period_stationarity_brackets.jl")
    @test isfile(bracket_script_path)
    bracket_script = read(bracket_script_path, String)
    @test occursin("period_stationarity_brackets.svg", bracket_script)
    @test occursin("period_local_left_objective", bracket_script)
    @test occursin("period_local_right_objective", bracket_script)
    @test occursin("period_stationarity_margin", bracket_script)
    packet_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_publication_packet.jl")
    @test isfile(packet_script_path)
    packet_script = read(packet_script_path, String)
    @test occursin("publication_evidence_packet.md", packet_script)
    @test occursin("claim_matrix.csv", packet_script)
    @test occursin("figure_inventory.csv", packet_script)
    @test occursin("finite_amplitude_mechanism_summary.csv", packet_script)
    @test occursin("bvk1_closer_than_ud_count", packet_script)
    @test occursin("morphology_phase_readiness", packet_script)
    @test occursin("Uneyama-Doi", packet_script)
    @test occursin("ud_energy_density_normalization_gate", packet_script)
    @test occursin("ud_strong_form_equivalence", packet_script)
    @test occursin("ud_morphology_identity_gate", packet_script)
    @test occursin("ud_symmetry_constrained_cell_scan", packet_script)
    @test occursin("ud_endpoint_refinement_gate", packet_script)
    @test occursin("ud_representative_phase_rank_gate", packet_script)
    @test occursin("ud_gyroid_random_basin_probe", packet_script)
    @test occursin("ud_gyroid_first_star_phase_probe", packet_script)
    @test occursin("ud_gyroid_first_star_optimizer_probe", packet_script)
    @test occursin("ud_gyroid_projected_basin_shape_probe", packet_script)
    @test occursin("ud_gyroid_lattice_classifier_design.md", packet_script)
    @test occursin("ud_gyroid_lattice_classifier_probe", packet_script)
    @test occursin("ud_gyroid_cubic_continuation_probe", packet_script)
    @test occursin("ud_gyroid_first_star_cubic_continuation_probe", packet_script)
    @test occursin("ud_gyroid_neighbor_warm_start_probe", packet_script)
    @test occursin("SI Fig./Table 24", packet_script)
    @test occursin("ud_gyroid_target_resolution_scan", packet_script)
    @test occursin("SI Fig./Table 25", packet_script)
    @test occursin("ud_gyroid_wide_scale_diagnostic", packet_script)
    @test occursin("SI Fig./Table 26", packet_script)
    @test occursin("ud_gyroid_classifier_constrained_period_scan", packet_script)
    @test occursin("SI Fig./Table 27", packet_script)
    @test occursin("ud_gyroid_scale_local_phase_optimizer_scan", packet_script)
    @test occursin("SI Fig./Table 28", packet_script)
    @test occursin("simple bicontinuity", packet_script)
    @test occursin("burp_bvk1_gyr_polyorder_initialized_probe", packet_script)
    @test occursin("SI Fig./Table 39", packet_script)
    @test occursin("burp_bvk1_gyr_reduced_shape_optimizer", packet_script)
    @test occursin("SI Fig./Table 40", packet_script)
    robustness_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_diblock_saddle_robustness.jl")
    @test isfile(robustness_script_path)
    robustness_script = read(robustness_script_path, String)
    @test occursin("saddle_robustness_summary.svg", robustness_script)
    @test occursin("mode_count", robustness_script)
    @test occursin("period_strategy=:bootstrap_local", robustness_script)
    morphology_audit_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_extension_audit.jl")
    @test isfile(morphology_audit_script_path)
    morphology_audit_script = read(morphology_audit_script_path, String)
    @test occursin("morphology_extension_audit.md", morphology_audit_script)
    @test occursin("morphology_gap_matrix.csv", morphology_audit_script)
    @test occursin("phase_diagram_gap", morphology_audit_script)
    @test occursin("morphology_phase_readiness", morphology_audit_script)
    @test occursin("morphology_saturation_chin_scan", morphology_audit_script)
    @test occursin("ud_energy_density_normalization_gate", morphology_audit_script)
    @test occursin("ud_strong_form_equivalence", morphology_audit_script)
    @test occursin("has_burp_bvk1_gyroid_polyorder_initialized_probe",
        morphology_audit_script)
    @test occursin("has_burp_bvk1_gyroid_reduced_shape_optimizer",
        morphology_audit_script)
    @test occursin("ud_morphology_identity_gate", morphology_audit_script)
    @test occursin("ud_symmetry_constrained_cell_scan", morphology_audit_script)
    @test occursin("ud_endpoint_refinement_gate", morphology_audit_script)
    @test occursin("ud_representative_phase_rank_gate", morphology_audit_script)
    @test occursin("ud_gyroid_random_basin_probe", morphology_audit_script)
    @test occursin("ud_gyroid_first_star_phase_probe", morphology_audit_script)
    @test occursin("ud_gyroid_first_star_optimizer_probe", morphology_audit_script)
    @test occursin("ud_gyroid_projected_basin_shape_probe", morphology_audit_script)
    @test occursin("ud_gyroid_lattice_classifier_design.md", morphology_audit_script)
    @test occursin("ud_gyroid_lattice_classifier_probe", morphology_audit_script)
    @test occursin("ud_gyroid_cubic_continuation_probe", morphology_audit_script)
    @test occursin("ud_gyroid_first_star_cubic_continuation_probe",
        morphology_audit_script)
    @test occursin("ud_gyroid_neighbor_warm_start_probe", morphology_audit_script)
    @test occursin("ud_gyroid_target_resolution_scan", morphology_audit_script)
    @test occursin("ud_gyroid_wide_scale_diagnostic", morphology_audit_script)
    @test occursin("ud_gyroid_classifier_constrained_period_scan",
        morphology_audit_script)
    @test occursin("ud_gyroid_scale_local_phase_optimizer_scan",
        morphology_audit_script)
    @test occursin("simple bicontinuity", morphology_audit_script)
    @test occursin("target/resolution", morphology_audit_script)
    @test occursin("LAM", morphology_audit_script)
    @test occursin("CYL", morphology_audit_script)
    @test occursin("GYR", morphology_audit_script)
    nd_equivalence_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_nd_energy_equivalence.jl")
    @test isfile(nd_equivalence_script_path)
    nd_equivalence_script = read(nd_equivalence_script_path, String)
    @test occursin("nd_energy_equivalence.csv", nd_equivalence_script)
    @test occursin("diblock_burp_ti_energy_nd", nd_equivalence_script)
    @test occursin("diblock_bvk1_energy_nd", nd_equivalence_script)
    ud_nonlamellar_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_nonlamellar_reproduction_baseline.jl")
    @test isfile(ud_nonlamellar_script_path)
    ud_nonlamellar_script = read(ud_nonlamellar_script_path, String)
    @test occursin("ud_nonlamellar_reproduction_baseline.csv",
        ud_nonlamellar_script)
    @test occursin("AB_phasediagram.eps", ud_nonlamellar_script)
    @test occursin("relaxation_method", ud_nonlamellar_script)
    @test occursin("psi = sqrt(phi)", ud_nonlamellar_script)
    @test occursin("not a BURP-TI/BVK1 failure", ud_nonlamellar_script)
    ud_psi_relaxation_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_nonlamellar_psi_relaxation.jl")
    @test isfile(ud_psi_relaxation_script_path)
    ud_psi_relaxation_script = read(ud_psi_relaxation_script_path, String)
    @test occursin("ud_nonlamellar_psi_relaxation.csv", ud_psi_relaxation_script)
    @test occursin("relax_diblock_uneyama_doi_morphology_seed_nd",
        ud_psi_relaxation_script)
    @test occursin("psi = sqrt(phi)", ud_psi_relaxation_script)
    @test occursin("fixed-cell relaxation", ud_psi_relaxation_script)
    ud_cell_scan_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_nonlamellar_cell_scan.jl")
    @test isfile(ud_cell_scan_script_path)
    ud_cell_scan_script = read(ud_cell_scan_script_path, String)
    @test occursin("ud_nonlamellar_cell_scan.csv", ud_cell_scan_script)
    @test occursin("scan_diblock_uneyama_doi_morphology_cell_scales",
        ud_cell_scan_script)
    @test occursin("optimize_diblock_uneyama_doi_morphology_cell_scale",
        ud_cell_scan_script)
    @test occursin("optimize_diblock_uneyama_doi_morphology_cell_scale_multistart",
        ud_cell_scan_script)
    @test occursin("initial_amplitudes", ud_cell_scan_script)
    @test occursin("adaptive scalar", ud_cell_scan_script)
    @test occursin("wide scalar-bound", ud_cell_scan_script)
    @test occursin("_arg_value(\"--min-scale-factor\", \"0.35\")",
        ud_cell_scan_script)
    @test occursin("_arg_value(\"--max-scale-factor\", \"2.4\")",
        ud_cell_scan_script)
    @test occursin("_arg_value(\"--max-scale-expansions\", \"12\")",
        ud_cell_scan_script)
    @test occursin("stress-free-scale", ud_cell_scan_script)
    @test occursin("not a phase diagram", ud_cell_scan_script)
    ud_saturation_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_saturation_grid_check.jl")
    @test isfile(ud_saturation_script_path)
    ud_saturation_script = read(ud_saturation_script_path, String)
    @test occursin("ud_saturation_grid_check.csv", ud_saturation_script)
    @test occursin("saturation-flagged rows", ud_saturation_script)
    @test occursin("grid/iteration follow-up", ud_saturation_script)
    @test occursin("not a phase diagram", ud_saturation_script)
    @test occursin("optimize_diblock_uneyama_doi_morphology_cell_scale_multistart",
        ud_saturation_script)
    ud_shape_scan_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_nonlamellar_cell_shape_scan.jl")
    @test isfile(ud_shape_scan_script_path)
    ud_shape_scan_script = read(ud_shape_scan_script_path, String)
    @test occursin("ud_nonlamellar_cell_shape_scan.csv", ud_shape_scan_script)
    @test occursin("scan_diblock_uneyama_doi_morphology_cell_shapes",
        ud_shape_scan_script)
    @test occursin("anisotropic", ud_shape_scan_script)
    @test occursin("shape_delta", ud_shape_scan_script)
    @test occursin("not a phase diagram", ud_shape_scan_script)
    ud_shape_opt_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_cell_shape_coordinate_optimization.jl")
    @test isfile(ud_shape_opt_script_path)
    ud_shape_opt_script = read(ud_shape_opt_script_path, String)
    @test occursin("ud_cell_shape_coordinate_optimization.csv",
        ud_shape_opt_script)
    @test occursin("optimize_diblock_uneyama_doi_morphology_cell_shape_coordinate",
        ud_shape_opt_script)
    @test occursin("bracketing", ud_shape_opt_script)
    @test occursin("not a phase diagram", ud_shape_opt_script)
    ud_multi_shape_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_multicoordinate_lattice_optimization.jl")
    @test isfile(ud_multi_shape_script_path)
    ud_multi_shape_script = read(ud_multi_shape_script_path, String)
    @test occursin("ud_multicoordinate_lattice_optimization.csv",
        ud_multi_shape_script)
    @test occursin("optimize_diblock_uneyama_doi_morphology_cell_shape_multicoordinate",
        ud_multi_shape_script)
    @test occursin("multi-coordinate", ud_multi_shape_script)
    @test occursin("not a phase diagram", ud_multi_shape_script)
    ud_production_ladder_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_production_convergence_ladder.jl")
    @test isfile(ud_production_ladder_script_path)
    ud_production_ladder_script = read(ud_production_ladder_script_path, String)
    @test occursin("ud_production_convergence_ladder.csv",
        ud_production_ladder_script)
    @test occursin("projected_gradient_norm", ud_production_ladder_script)
    @test occursin("convergence ladder", ud_production_ladder_script)
    @test occursin("not a phase diagram", ud_production_ladder_script)
    seed_audit_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_seed_audit.jl")
    @test isfile(seed_audit_script_path)
    seed_audit_script = read(seed_audit_script_path, String)
    @test occursin("morphology_seed_audit.csv", seed_audit_script)
    @test occursin("diblock_morphology_seed", seed_audit_script)
    @test occursin("LAM", seed_audit_script)
    @test occursin("CYL", seed_audit_script)
    @test occursin("GYR", seed_audit_script)
    amplitude_audit_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_amplitude_envelope.jl")
    @test isfile(amplitude_audit_script_path)
    amplitude_audit_script = read(amplitude_audit_script_path, String)
    @test occursin("morphology_amplitude_envelope.csv", amplitude_audit_script)
    @test occursin("optimize_diblock_morphology_seed_amplitude", amplitude_audit_script)
    @test occursin("not a full saddle", amplitude_audit_script)
    phase_readiness_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_phase_readiness.jl")
    @test isfile(phase_readiness_script_path)
    phase_readiness_script = read(phase_readiness_script_path, String)
    @test occursin("morphology_phase_readiness_summary.csv", phase_readiness_script)
    @test occursin("energy_density", phase_readiness_script)
    @test occursin("not a phase diagram", phase_readiness_script)
    relaxation_audit_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_pilot_relaxation.jl")
    @test isfile(relaxation_audit_script_path)
    relaxation_audit_script = read(relaxation_audit_script_path, String)
    @test occursin("morphology_pilot_relaxation.csv", relaxation_audit_script)
    @test occursin("relax_diblock_morphology_seed_nd", relaxation_audit_script)
    @test occursin("pilot relaxation", relaxation_audit_script)
    stationarity_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_stationarity_probe.jl")
    @test isfile(stationarity_script_path)
    stationarity_script = read(stationarity_script_path, String)
    @test occursin("morphology_stationarity_probe.csv", stationarity_script)
    @test occursin("probe_diblock_morphology_relaxation_stationarity", stationarity_script)
    @test occursin("stationarity probe", stationarity_script)
    stress_test_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_relaxation_stress_test.jl")
    @test isfile(stress_test_script_path)
    stress_test_script = read(stress_test_script_path, String)
    @test occursin("morphology_relaxation_stress_test.csv", stress_test_script)
    @test occursin("saturation_warning", stress_test_script)
    @test occursin("not a phase diagram", stress_test_script)
    saturation_scan_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_burp_bvk1_morphology_saturation_chin_scan.jl")
    @test isfile(saturation_scan_script_path)
    saturation_scan_script = read(saturation_scan_script_path, String)
    @test occursin("morphology_saturation_chin_scan.csv", saturation_scan_script)
    @test occursin("chiN_values", saturation_scan_script)
    @test occursin("saturation_warning", saturation_scan_script)
    @test occursin("not a phase diagram", saturation_scan_script)
    ud_symmetric_gate_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_symmetric_phase_gate.jl")
    @test isfile(ud_symmetric_gate_script_path)
    ud_symmetric_gate_script = read(ud_symmetric_gate_script_path, String)
    @test occursin("ud_symmetric_phase_gate.csv", ud_symmetric_gate_script)
    @test occursin("expected_morphology", ud_symmetric_gate_script)
    @test occursin("symmetric-composition", ud_symmetric_gate_script)
    ud_energy_norm_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_energy_density_normalization_gate.jl")
    @test isfile(ud_energy_norm_script_path)
    ud_energy_norm_script = read(ud_energy_norm_script_path, String)
    @test occursin("ud_energy_density_normalization_gate.csv", ud_energy_norm_script)
    @test occursin("energy-density normalization", ud_energy_norm_script)
    @test occursin("lamellar/nonlamellar", ud_energy_norm_script)
    ud_strong_form_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_strong_form_equivalence.jl")
    @test isfile(ud_strong_form_script_path)
    ud_strong_form_script = read(ud_strong_form_script_path, String)
    @test occursin("ud_strong_form_equivalence.csv", ud_strong_form_script)
    @test occursin("freeenergy_diblock_strong", ud_strong_form_script)
    @test occursin("tau-local", ud_strong_form_script)
    @test occursin("not a phase diagram", ud_strong_form_script)
    ud_identity_gate_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_morphology_identity_gate.jl")
    @test isfile(ud_identity_gate_script_path)
    ud_identity_gate_script = read(ud_identity_gate_script_path, String)
    @test occursin("ud_morphology_identity_gate.csv", ud_identity_gate_script)
    @test occursin("target_star_q_spread", ud_identity_gate_script)
    @test occursin("morphology identity", ud_identity_gate_script)
    @test occursin("not a phase diagram", ud_identity_gate_script)
    ud_symmetry_cell_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_symmetry_constrained_cell_scan.jl")
    @test isfile(ud_symmetry_cell_script_path)
    ud_symmetry_cell_script = read(ud_symmetry_cell_script_path, String)
    @test occursin("ud_symmetry_constrained_cell_scan.csv",
        ud_symmetry_cell_script)
    @test occursin("sqrt(3.0)", ud_symmetry_cell_script)
    @test occursin("morphology-preserving", ud_symmetry_cell_script)
    @test occursin("target_star_q_spread", ud_symmetry_cell_script)
    @test occursin("not a phase diagram", ud_symmetry_cell_script)
    ud_refinement_gate_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_endpoint_refinement_gate.jl")
    @test isfile(ud_refinement_gate_script_path)
    ud_refinement_gate_script = read(ud_refinement_gate_script_path, String)
    @test occursin("ud_endpoint_refinement_gate.csv", ud_refinement_gate_script)
    @test occursin("dims_multiplier", ud_refinement_gate_script)
    @test occursin("raw_then_projected", ud_refinement_gate_script)
    @test occursin("endpoint-adjacent", ud_refinement_gate_script)
    @test occursin("not a phase diagram", ud_refinement_gate_script)
    ud_representative_gate_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_representative_phase_rank_gate.jl")
    @test isfile(ud_representative_gate_script_path)
    ud_representative_gate_script = read(ud_representative_gate_script_path, String)
    @test occursin("ud_representative_phase_rank_gate.csv",
        ud_representative_gate_script)
    @test occursin("ud_representative_phase_rank_gate_aggregate.csv",
        ud_representative_gate_script)
    @test occursin("f0.15_chiN50_BCC", ud_representative_gate_script)
    @test occursin("f0.30_chiN50_CYL", ud_representative_gate_script)
    @test occursin("f0.39_chiN20_GYR", ud_representative_gate_script)
    @test occursin("AB_phasediagram.eps", ud_representative_gate_script)
    @test occursin("endpoint_allowed_by_ud_psi", ud_representative_gate_script)
    @test occursin("morphology_identity_pass", ud_representative_gate_script)
    @test occursin("pass_count", ud_representative_gate_script)
    @test occursin("--min-scale-factor", ud_representative_gate_script)
    @test occursin("--max-scale-factor", ud_representative_gate_script)
    @test occursin("not a full phase diagram", ud_representative_gate_script)
    ud_gyroid_random_probe_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_random_basin_probe.jl")
    @test isfile(ud_gyroid_random_probe_script_path)
    ud_gyroid_random_probe_script = read(ud_gyroid_random_probe_script_path, String)
    @test occursin("ud_gyroid_random_basin_probe.csv",
        ud_gyroid_random_probe_script)
    @test occursin("f0.39_chiN20_GYR", ud_gyroid_random_probe_script)
    @test occursin("random_theta", ud_gyroid_random_probe_script)
    @test occursin("best_shift_corr", ud_gyroid_random_probe_script)
    @test occursin("star110_power_fraction", ud_gyroid_random_probe_script)
    @test occursin("not a phase diagram", ud_gyroid_random_probe_script)
    ud_gyroid_phase_probe_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_first_star_phase_probe.jl")
    @test isfile(ud_gyroid_phase_probe_script_path)
    ud_gyroid_phase_probe_script = read(ud_gyroid_phase_probe_script_path, String)
    @test occursin("ud_gyroid_first_star_phase_modes.csv",
        ud_gyroid_phase_probe_script)
    @test occursin("ud_gyroid_first_star_phase_summary.csv",
        ud_gyroid_phase_probe_script)
    @test occursin("f0.39_chiN20_GYR", ud_gyroid_phase_probe_script)
    @test occursin("star110_modes", ud_gyroid_phase_probe_script)
    @test occursin("phase_residual", ud_gyroid_phase_probe_script)
    @test occursin("not a phase diagram", ud_gyroid_phase_probe_script)
    ud_gyroid_optimizer_probe_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_first_star_optimizer_probe.jl")
    @test isfile(ud_gyroid_optimizer_probe_script_path)
    ud_gyroid_optimizer_probe_script = read(ud_gyroid_optimizer_probe_script_path,
        String)
    @test occursin("ud_gyroid_first_star_optimizer_candidates.csv",
        ud_gyroid_optimizer_probe_script)
    @test occursin("ud_gyroid_first_star_optimizer_summary.csv",
        ud_gyroid_optimizer_probe_script)
    @test occursin("Optim.NelderMead", ud_gyroid_optimizer_probe_script)
    @test occursin("positive_star110_modes", ud_gyroid_optimizer_probe_script)
    @test occursin("relax_selected_first_star_seeds", ud_gyroid_optimizer_probe_script)
    @test occursin("not a phase diagram", ud_gyroid_optimizer_probe_script)
    ud_gyroid_shape_probe_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_projected_basin_shape_probe.jl")
    @test isfile(ud_gyroid_shape_probe_script_path)
    ud_gyroid_shape_probe_script = read(ud_gyroid_shape_probe_script_path, String)
    @test occursin("ud_gyroid_projected_basin_shape_scan.csv",
        ud_gyroid_shape_probe_script)
    @test occursin("ud_gyroid_projected_basin_shape_summary.csv",
        ud_gyroid_shape_probe_script)
    @test occursin("volume_preserving_lengths", ud_gyroid_shape_probe_script)
    @test occursin("shape_direction", ud_gyroid_shape_probe_script)
    @test occursin("q2_weighted_spread", ud_gyroid_shape_probe_script)
    @test occursin("not a phase diagram", ud_gyroid_shape_probe_script)
    ud_gyroid_classifier_design_path = joinpath(@__DIR__, "..", "docs",
        "ud_gyroid_lattice_classifier_design.md")
    @test isfile(ud_gyroid_classifier_design_path)
    ud_gyroid_classifier_design = read(ud_gyroid_classifier_design_path, String)
    @test occursin("classifier-first", ud_gyroid_classifier_design)
    @test occursin("P_{110}", ud_gyroid_classifier_design)
    @test occursin("R_GYR", ud_gyroid_classifier_design)
    @test occursin("not a valid gyroid point", ud_gyroid_classifier_design)
    ud_gyroid_classifier_probe_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_lattice_classifier_probe.jl")
    @test isfile(ud_gyroid_classifier_probe_script_path)
    ud_gyroid_classifier_probe_script = read(ud_gyroid_classifier_probe_script_path,
        String)
    @test occursin("ud_gyroid_lattice_classifier_probe.csv",
        ud_gyroid_classifier_probe_script)
    @test occursin("ud_gyroid_lattice_classifier_summary.csv",
        ud_gyroid_classifier_probe_script)
    @test occursin("cubic_rotation_matrices", ud_gyroid_classifier_probe_script)
    @test occursin("phase_residual_to_reference", ud_gyroid_classifier_probe_script)
    @test occursin("q2_weighted_spread", ud_gyroid_classifier_probe_script)
    @test occursin("group_balance", ud_gyroid_classifier_probe_script)
    @test occursin("classifier-first", ud_gyroid_classifier_probe_script)
    @test occursin("not a phase diagram", ud_gyroid_classifier_probe_script)
    ud_gyroid_continuation_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_cubic_continuation_probe.jl")
    @test isfile(ud_gyroid_continuation_script_path)
    ud_gyroid_continuation_script = read(ud_gyroid_continuation_script_path,
        String)
    @test occursin("ud_gyroid_cubic_continuation_rows.csv",
        ud_gyroid_continuation_script)
    @test occursin("ud_gyroid_cubic_continuation_summary.csv",
        ud_gyroid_continuation_script)
    @test occursin("phase_start_label", ud_gyroid_continuation_script)
    @test occursin("classifier_label", ud_gyroid_continuation_script)
    @test occursin("accepted_as_gyroid", ud_gyroid_continuation_script)
    @test occursin("classifier_accepted_gyr_count", ud_gyroid_continuation_script)
    @test occursin("phase_rank_pass", ud_gyroid_continuation_script)
    @test occursin("ud_gyr_cubic_continuation_no_local_accepted_gyr",
        ud_gyroid_continuation_script)
    @test occursin("cubic GYR period/phase continuation", ud_gyroid_continuation_script)
    @test occursin("not a phase diagram", ud_gyroid_continuation_script)
    ud_gyroid_first_star_continuation_script_path = joinpath(@__DIR__, "..",
        "scripts", "write_uneyama_doi_gyroid_first_star_cubic_continuation_probe.jl")
    @test isfile(ud_gyroid_first_star_continuation_script_path)
    ud_gyroid_first_star_continuation_script = read(
        ud_gyroid_first_star_continuation_script_path, String)
    @test occursin("ud_gyroid_first_star_cubic_continuation_rows.csv",
        ud_gyroid_first_star_continuation_script)
    @test occursin("ud_gyroid_first_star_cubic_continuation_summary.csv",
        ud_gyroid_first_star_continuation_script)
    @test occursin("first-star cubic continuation", ud_gyroid_first_star_continuation_script)
    @test occursin("best_random_projection", ud_gyroid_first_star_continuation_script)
    @test occursin("optimized_first_star", ud_gyroid_first_star_continuation_script)
    @test occursin("classifier_accepted_gyr_count",
        ud_gyroid_first_star_continuation_script)
    @test occursin("accepted_for_gyr_ranking",
        ud_gyroid_first_star_continuation_script)
    @test occursin("phase_rank_pass", ud_gyroid_first_star_continuation_script)
    @test occursin("topology_label", ud_gyroid_first_star_continuation_script)
    @test occursin("bicontinuous_periodic", ud_gyroid_first_star_continuation_script)
    @test occursin("_periodic_component_summary",
        ud_gyroid_first_star_continuation_script)
    @test occursin("not a phase diagram", ud_gyroid_first_star_continuation_script)
    ud_gyroid_warm_start_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_neighbor_warm_start_probe.jl")
    @test isfile(ud_gyroid_warm_start_script_path)
    ud_gyroid_warm_start_script = read(ud_gyroid_warm_start_script_path, String)
    @test occursin("ud_gyroid_neighbor_warm_start_rows.csv",
        ud_gyroid_warm_start_script)
    @test occursin("ud_gyroid_neighbor_warm_start_summary.csv",
        ud_gyroid_warm_start_script)
    @test occursin("neighbor warm-start", ud_gyroid_warm_start_script)
    @test occursin("warm_start_label", ud_gyroid_warm_start_script)
    @test occursin("annealing_path", ud_gyroid_warm_start_script)
    @test occursin("target_f", ud_gyroid_warm_start_script)
    @test occursin("target_chiN", ud_gyroid_warm_start_script)
    @test occursin("classifier_accepted_gyr_count", ud_gyroid_warm_start_script)
    @test occursin("phase_rank_pass", ud_gyroid_warm_start_script)
    @test occursin("not a phase diagram", ud_gyroid_warm_start_script)
    ud_gyroid_target_resolution_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_target_resolution_scan.jl")
    @test isfile(ud_gyroid_target_resolution_script_path)
    ud_gyroid_target_resolution_script =
        read(ud_gyroid_target_resolution_script_path, String)
    @test occursin("ud_gyroid_target_resolution_scan_rows.csv",
        ud_gyroid_target_resolution_script)
    @test occursin("ud_gyroid_target_resolution_scan_gyr_rows.csv",
        ud_gyroid_target_resolution_script)
    @test occursin("ud_gyroid_target_resolution_scan_paths.csv",
        ud_gyroid_target_resolution_script)
    @test occursin("ud_gyroid_target_resolution_scan_summary.csv",
        ud_gyroid_target_resolution_script)
    @test occursin("ud_gyroid_target_resolution_scan_aggregate.csv",
        ud_gyroid_target_resolution_script)
    @test occursin("target/resolution", ud_gyroid_target_resolution_script)
    @test occursin("gyr_base_scale_factor", ud_gyroid_target_resolution_script)
    @test occursin("target_f", ud_gyroid_target_resolution_script)
    @test occursin("target_chiN", ud_gyroid_target_resolution_script)
    @test occursin("accepted_lam_energy_density", ud_gyroid_target_resolution_script)
    @test occursin("accepted_cyl_energy_density", ud_gyroid_target_resolution_script)
    @test occursin("accepted_gyr_energy_density", ud_gyroid_target_resolution_script)
    @test occursin("phase_rank_pass", ud_gyroid_target_resolution_script)
    @test occursin("not a phase diagram", ud_gyroid_target_resolution_script)
    ud_gyroid_wide_scale_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_wide_scale_diagnostic.jl")
    @test isfile(ud_gyroid_wide_scale_script_path)
    ud_gyroid_wide_scale_script = read(ud_gyroid_wide_scale_script_path, String)
    @test occursin("ud_gyroid_wide_scale_diagnostic_summary.csv",
        ud_gyroid_wide_scale_script)
    @test occursin("ud_gyroid_wide_scale_diagnostic.md",
        ud_gyroid_wide_scale_script)
    @test occursin("write_uneyama_doi_gyroid_target_resolution_scan.jl",
        ud_gyroid_wide_scale_script)
    @test occursin("wide cubic scale diagnostic", ud_gyroid_wide_scale_script)
    @test occursin("not_first_star_dominant", ud_gyroid_wide_scale_script)
    @test occursin("not a phase diagram", ud_gyroid_wide_scale_script)
    ud_gyroid_classifier_period_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_classifier_constrained_period_scan.jl")
    @test isfile(ud_gyroid_classifier_period_script_path)
    ud_gyroid_classifier_period_script =
        read(ud_gyroid_classifier_period_script_path, String)
    @test occursin("ud_gyroid_classifier_constrained_period_summary.csv",
        ud_gyroid_classifier_period_script)
    @test occursin("ud_gyroid_classifier_constrained_period_scan.md",
        ud_gyroid_classifier_period_script)
    @test occursin("write_uneyama_doi_gyroid_target_resolution_scan.jl",
        ud_gyroid_classifier_period_script)
    @test occursin("classifier-constrained period scan",
        ud_gyroid_classifier_period_script)
    @test occursin("best_gyroid_like_energy_density",
        ud_gyroid_classifier_period_script)
    @test occursin("transition_bracket_width", ud_gyroid_classifier_period_script)
    @test occursin("not a phase diagram", ud_gyroid_classifier_period_script)
    ud_gyroid_scale_local_phase_script_path = joinpath(@__DIR__, "..", "scripts",
        "write_uneyama_doi_gyroid_scale_local_phase_optimizer_scan.jl")
    @test isfile(ud_gyroid_scale_local_phase_script_path)
    ud_gyroid_scale_local_phase_script =
        read(ud_gyroid_scale_local_phase_script_path, String)
    @test occursin("ud_gyroid_scale_local_phase_optimizer_summary.csv",
        ud_gyroid_scale_local_phase_script)
    @test occursin("ud_gyroid_scale_local_phase_optimizer_scan.md",
        ud_gyroid_scale_local_phase_script)
    @test occursin("write_uneyama_doi_gyroid_first_star_cubic_continuation_probe.jl",
        ud_gyroid_scale_local_phase_script)
    @test occursin("scale-local first-star phase optimizer",
        ud_gyroid_scale_local_phase_script)
    @test occursin("best_phase_optimized_gyroid_energy_density",
        ud_gyroid_scale_local_phase_script)
    @test occursin("best_phase_optimized_gyroid_topology_label",
        ud_gyroid_scale_local_phase_script)
    @test occursin("bicontinuous_row_count", ud_gyroid_scale_local_phase_script)
    @test occursin("phase_variable_rank_pass", ud_gyroid_scale_local_phase_script)
    @test occursin("not a phase diagram", ud_gyroid_scale_local_phase_script)

    kernel_min = DFMMonteCarlo._diblock_kernel_minimum_k2(; f=0.5, N=1.0, b=1.0)
    ok_at_min = diblock_ohta_kawasaki_kernel(kernel_min; f=0.5)
    @test isfinite(ok_at_min)
    @test ok_at_min < diblock_ohta_kawasaki_kernel(kernel_min / 10.0; f=0.5)
    @test ok_at_min < diblock_ohta_kawasaki_kernel(kernel_min * 10.0; f=0.5)
    @test isapprox(ok_at_min, diblock_composition_kernel(kernel_min; f=0.5); rtol=1.0e-6)

    ud_lamella = minimize_diblock_uneyama_doi_lamella(; f=0.5, chiN=12.0, nx=64,
        mode_count=3, max_iterations=4_000)
    ok_lamella = minimize_diblock_ohta_kawasaki_lamella(; f=0.5, chiN=12.0, nx=64,
        mode_count=3, max_iterations=5_000, relaxation_step=1.0e-2,
        tolerance=1.0e-7)
    @test ud_lamella.converged
    @test ok_lamella.converged
    @test 0.0 < ud_lamella.minimum_phi < ud_lamella.maximum_phi < 1.0
    @test 0.0 < ok_lamella.minimum_phi < ok_lamella.maximum_phi < 1.0
    @test ud_lamella.minimum_phi < 0.40
    @test ud_lamella.maximum_phi > 0.60
    @test ud_lamella.maximum_phi - ud_lamella.minimum_phi > 0.1
    @test ok_lamella.maximum_phi - ok_lamella.minimum_phi > 0.1
end

@testset "Gaussian perturbation Monte Carlo" begin
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_gaussian_perturbation_mc.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_gaussian_analysis_comparison_suite.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_gaussian_sample_size_study.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "write_validation_audit.jl"))
    @test occursin("gaussian_field_analysis_comparison",
        read(joinpath(@__DIR__, "..", "scripts", "write_validation_audit.jl"), String))

    grid = Grid3D(16, 6.4)
    params = HEURPParams(A=300.0, B=1.0, C=20.0, alpha=0.1)
    contributions = DFMMonteCarlo._dfm_pressure_mode_contributions(params, grid)
    @test length(contributions) > 100
    @test DFMMonteCarlo.mean_field_pressure(params.B, params.C) + sum(contributions) ≈
          dfm_pressure_discrete(params.A, params.B, params.C, grid; alpha=params.alpha)
    summary = DFMMonteCarlo._gaussian_mode_pressure_summary(
        params,
        grid;
        samples=5000,
        chains=4,
        seed=17,
        pressure_tolerance=0.2,
    )
    @test :reference_error in propertynames(summary)
    @test :reference_zscore in propertynames(summary)
    @test :reference_agreement_pass in propertynames(summary)
    @test :within_chain_stderr in propertynames(summary)
    @test :between_chain_stderr in propertynames(summary)
    @test :chain_pressure_std in propertynames(summary)
    @test :required_effective_sample_size in propertynames(summary)
    @test :required_samples_per_chain in propertynames(summary)
    @test :sampling_budget_pass in propertynames(summary)
    @test summary.reference_error ≈ abs(summary.mc_pressure - summary.dfm_pressure)
    @test summary.reference_zscore ≈ summary.reference_error / summary.mc_error
    @test summary.within_chain_stderr >= 0.0
    @test summary.between_chain_stderr >= 0.0
    @test summary.chain_pressure_std >= 0.0
    @test summary.mc_error == max(summary.within_chain_stderr, summary.between_chain_stderr, 1.0e-12)
    @test summary.required_effective_sample_size >= 1.0
    @test summary.required_samples_per_chain >= 1
    @test summary.diagnostic_pass == min(summary.pressure_error_pass,
        summary.reference_agreement_pass, summary.sampling_budget_pass)

    field_summary = DFMMonteCarlo._gaussian_field_pressure_summary(
        params,
        grid;
        samples=30_000,
        chains=4,
        seed=29,
        pressure_tolerance=0.5,
    )
    @test field_summary.sampler == "gaussian_field_mc"
    @test field_summary.dfm_pressure ≈ dfm_pressure_discrete(params.A, params.B, params.C, grid; alpha=params.alpha)
    @test field_summary.effective_sample_size == 120_000.0
    @test field_summary.mc_error > 0.0
    @test field_summary.reference_error <= max(field_summary.reference_tolerance,
        field_summary.reference_sigma_tolerance * field_summary.mc_error)
    @test field_summary.diagnostic_pass == min(field_summary.pressure_error_pass,
        field_summary.reference_agreement_pass, field_summary.sampling_budget_pass)
    no_mode_params = HEURPParams(A=300.0, B=1.0, C=0.1, alpha=0.1)
    no_mode_field_summary = DFMMonteCarlo._gaussian_field_pressure_summary(
        no_mode_params,
        grid;
        samples=1000,
        chains=2,
        seed=31,
        pressure_tolerance=0.5,
    )
    @test no_mode_field_summary.stable_mode_count == 0
    @test no_mode_field_summary.mc_pressure ≈ DFMMonteCarlo.mean_field_pressure(no_mode_params.B, no_mode_params.C)
    @test no_mode_field_summary.dfm_pressure ≈ no_mode_field_summary.mc_pressure
    @test no_mode_field_summary.reference_error == 0.0
    @test no_mode_field_summary.diagnostic_pass == 1.0

    configuration_grid = Grid3D(8, 6.4)
    configuration_summary = DFMMonteCarlo._gaussian_field_configuration_pressure_summary(
        params,
        configuration_grid;
        samples=400,
        chains=3,
        seed=41,
        pressure_tolerance=Inf,
    )
    @test configuration_summary.sampler == "gaussian_field_configuration_mc"
    @test configuration_summary.dfm_pressure ≈
          dfm_pressure_discrete(params.A, params.B, params.C, configuration_grid; alpha=params.alpha)
    @test configuration_summary.stable_mode_count > 0
    @test configuration_summary.effective_sample_size == 1200.0
    @test configuration_summary.mc_error > 0.0
    @test configuration_summary.reference_error <= max(0.5,
        configuration_summary.reference_sigma_tolerance * configuration_summary.mc_error)

    no_mode_configuration_summary = DFMMonteCarlo._gaussian_field_configuration_pressure_summary(
        no_mode_params,
        grid;
        samples=1000,
        chains=2,
        seed=43,
        pressure_tolerance=0.5,
    )
    @test no_mode_configuration_summary.stable_mode_count == 0
    @test no_mode_configuration_summary.mc_pressure ≈
          DFMMonteCarlo.mean_field_pressure(no_mode_params.B, no_mode_params.C)
    @test no_mode_configuration_summary.dfm_pressure ≈ no_mode_configuration_summary.mc_pressure
    @test no_mode_configuration_summary.reference_error == 0.0
    @test no_mode_configuration_summary.diagnostic_pass == 1.0

    field_outdir = mktempdir()
    field_result = run_gaussian_perturbation_mc(field_outdir; Nx=16, C_values=[20.0],
        samples=30_000, chains=4, seed=29, pressure_tolerance=0.5, estimator=:field)
    field_metadata = TOML.parsefile(joinpath(field_outdir, "run_metadata.toml"))
    @test field_metadata["run"]["sampler"] == "gaussian_field_mc"
    field_header, field_rows = DFMMonteCarlo._read_csv_strings(field_result.diagnostics_csv)
    field_sampler_idx = findfirst(==("sampler"), field_header)
    field_pass_idx = findfirst(==("diagnostic_pass"), field_header)
    @test all(row[field_sampler_idx] == "gaussian_field_mc" for row in field_rows)
    @test all(parse(Float64, row[field_pass_idx]) == 1.0 for row in field_rows)

    configuration_outdir = mktempdir()
    configuration_result = run_gaussian_perturbation_mc(configuration_outdir; Nx=8, C_values=[20.0],
        samples=400, chains=3, seed=41, pressure_tolerance=Inf, estimator=:field_configuration)
    configuration_metadata = TOML.parsefile(joinpath(configuration_outdir, "run_metadata.toml"))
    @test configuration_metadata["run"]["sampler"] == "gaussian_field_configuration_mc"
    @test configuration_metadata["sampling"]["estimator"] == "field_configuration"
    @test configuration_metadata["estimator"]["formula_id"] ==
          "finite_grid_gaussian_field_configuration_pressure"
    configuration_header, configuration_rows =
        DFMMonteCarlo._read_csv_strings(configuration_result.diagnostics_csv)
    configuration_sampler_idx = findfirst(==("sampler"), configuration_header)
    @test all(row[configuration_sampler_idx] == "gaussian_field_configuration_mc"
        for row in configuration_rows)

    outdir = mktempdir()
    result = run_gaussian_perturbation_mc(outdir; Nx=16, C_values=[20.0],
        samples=5000, chains=4, seed=17, pressure_tolerance=0.2)
    @test isfile(result.pressure_csv)
    @test isfile(result.diagnostics_csv)
    gaussian_metadata_path = joinpath(outdir, "run_metadata.toml")
    @test isfile(gaussian_metadata_path)
    gaussian_metadata = TOML.parsefile(gaussian_metadata_path)
    @test gaussian_metadata["run"]["sampler"] == "gaussian_mode_mc"
    @test gaussian_metadata["run"]["seed"] == 17
    @test gaussian_metadata["sampling"]["samples"] == 5000
    @test gaussian_metadata["sampling"]["chains"] == 4
    @test gaussian_metadata["grid"]["Nx"] == 16
    @test gaussian_metadata["grid"]["target_Nx"] == 64
    @test gaussian_metadata["params"]["A"] == 300.0
    @test gaussian_metadata["params"]["B"] == 1.0
    @test gaussian_metadata["params"]["alpha"] == 0.1
    @test gaussian_metadata["estimator"]["formula_id"] == "finite_grid_gaussian_mode_sum"
    @test gaussian_metadata["environment"]["package"] == "DFMMonteCarlo"
    @test gaussian_metadata["outputs"]["pressure_csv"] == result.pressure_csv
    @test gaussian_metadata["outputs"]["diagnostics_csv"] == result.diagnostics_csv
    @test result.config.Nx == 16
    @test result.config.target_Nx == 64

    header, rows = DFMMonteCarlo._read_csv_numeric(result.pressure_csv)
    mc_idx = findfirst(==("mc_pressure"), header)
    err_idx = findfirst(==("mc_error"), header)
    dfm_idx = findfirst(==("dfm_pressure"), header)
    cl_idx = findfirst(==("cl_pressure_nearest"), header)
    mc_minus_dfm_idx = findfirst(==("mc_minus_dfm"), header)
    mc_dfm_zscore_idx = findfirst(==("mc_dfm_zscore"), header)
    mc_minus_cl_idx = findfirst(==("mc_minus_cl"), header)
    dfm_minus_cl_idx = findfirst(==("dfm_minus_cl"), header)
    cl_model_residual_idx = findfirst(==("cl_model_residual"), header)
    cl_context_pass_idx = findfirst(==("cl_context_pass"), header)
    reference_pass_idx = findfirst(==("reference_agreement_pass"), header)
    sampler_idx = findfirst(==("sampler"), header)
    @test abs(rows[1][mc_idx] - rows[1][dfm_idx]) <= 5.0 * rows[1][err_idx]
    @test rows[1][err_idx] <= 0.2
    @test rows[1][mc_minus_dfm_idx] ≈ rows[1][mc_idx] - rows[1][dfm_idx]
    @test rows[1][mc_dfm_zscore_idx] ≈ abs(rows[1][mc_minus_dfm_idx]) / rows[1][err_idx]
    @test rows[1][mc_minus_cl_idx] ≈ rows[1][mc_idx] - rows[1][cl_idx]
    @test rows[1][dfm_minus_cl_idx] ≈ rows[1][dfm_idx] - rows[1][cl_idx]
    @test rows[1][cl_model_residual_idx] ≈ rows[1][mc_minus_cl_idx] - rows[1][dfm_minus_cl_idx]
    @test rows[1][cl_model_residual_idx] ≈ rows[1][mc_minus_dfm_idx]
    @test rows[1][cl_context_pass_idx] == rows[1][reference_pass_idx]
    @test header[sampler_idx] == "sampler"

    diag_header, diag_rows = DFMMonteCarlo._read_csv_strings(result.diagnostics_csv)
    @test "samples" in diag_header
    @test "chain_count" in diag_header
    @test "within_chain_stderr" in diag_header
    @test "between_chain_stderr" in diag_header
    @test "chain_pressure_std" in diag_header
    @test "effective_sample_size" in diag_header
    @test "reference_error" in diag_header
    @test "reference_zscore" in diag_header
    @test "reference_agreement_pass" in diag_header
    @test "mc_minus_dfm" in diag_header
    @test "mc_dfm_zscore" in diag_header
    @test "mc_minus_cl" in diag_header
    @test "dfm_minus_cl" in diag_header
    @test "cl_model_residual" in diag_header
    @test "cl_context_pass" in diag_header
    @test "required_effective_sample_size" in diag_header
    @test "required_samples_per_chain" in diag_header
    @test "sampling_budget_pass" in diag_header
    @test "pressure_error_pass" in diag_header
    @test "diagnostic_pass" in diag_header
    @test "sampler" in diag_header
    diag_sampler_idx = findfirst(==("sampler"), diag_header)
    diag_pass_idx = findfirst(==("diagnostic_pass"), diag_header)
    diag_error_idx = findfirst(==("mc_error"), diag_header)
    diag_within_idx = findfirst(==("within_chain_stderr"), diag_header)
    diag_between_idx = findfirst(==("between_chain_stderr"), diag_header)
    @test all(row[diag_sampler_idx] == "gaussian_mode_mc" for row in diag_rows)
    @test all(parse(Float64, row[diag_pass_idx]) == 1.0 for row in diag_rows)
    @test all(parse(Float64, row[diag_error_idx]) ==
              max(parse(Float64, row[diag_within_idx]),
                  parse(Float64, row[diag_between_idx]), 1.0e-12) for row in diag_rows)

    adaptive_outdir = mktempdir()
    adaptive = run_gaussian_perturbation_mc(adaptive_outdir; Nx=16, C_values=[20.0],
        samples=10, max_samples=200, sample_growth=2.0, chains=4, seed=17,
        pressure_tolerance=0.2)
    adaptive_header, adaptive_rows = DFMMonteCarlo._read_csv_numeric(adaptive.diagnostics_csv)
    adaptive_samples_idx = findfirst(==("samples"), adaptive_header)
    adaptive_max_samples_idx = findfirst(==("max_samples"), adaptive_header)
    adaptive_attempts_idx = findfirst(==("sampling_attempts"), adaptive_header)
    adaptive_pass_idx = findfirst(==("diagnostic_pass"), adaptive_header)
    @test adaptive_rows[1][adaptive_samples_idx] > 10
    @test adaptive_rows[1][adaptive_samples_idx] <= 200
    @test adaptive_rows[1][adaptive_max_samples_idx] == 200
    @test adaptive_rows[1][adaptive_attempts_idx] > 1
    @test adaptive_rows[1][adaptive_pass_idx] == 1.0

    script_outdir = mktempdir()
    old_args = copy(ARGS)
    try
        empty!(ARGS)
        append!(ARGS, [
            "--outdir", script_outdir,
            "--nx", "16",
            "--c-values", "20.0",
            "--samples", "10",
            "--max-samples", "200",
            "--sample-growth", "2.0",
            "--chains", "4",
            "--pressure-tolerance", "0.2",
        ])
        redirect_stdout(devnull) do
            include(joinpath(@__DIR__, "..", "scripts", "run_gaussian_perturbation_mc.jl"))
        end
    finally
        empty!(ARGS)
        append!(ARGS, old_args)
    end
    @test isfile(joinpath(script_outdir, "gaussian_sampling_diagnostics.csv"))
    script_header, script_rows = DFMMonteCarlo._read_csv_numeric(joinpath(script_outdir, "gaussian_sampling_diagnostics.csv"))
    script_samples_idx = findfirst(==("samples"), script_header)
    script_pass_idx = findfirst(==("diagnostic_pass"), script_header)
    @test script_rows[1][script_samples_idx] > 10
    @test script_rows[1][script_pass_idx] == 1.0

    field_script_outdir = mktempdir()
    old_args = copy(ARGS)
    try
        empty!(ARGS)
        append!(ARGS, [
            "--outdir", field_script_outdir,
            "--nx", "16",
            "--c-values", "20.0",
            "--samples", "1000",
            "--max-samples", "1000",
            "--chains", "2",
            "--pressure-tolerance", "Inf",
            "--estimator", "field",
        ])
        redirect_stdout(devnull) do
            include(joinpath(@__DIR__, "..", "scripts", "run_gaussian_perturbation_mc.jl"))
        end
    finally
        empty!(ARGS)
        append!(ARGS, old_args)
    end
    field_script_metadata = TOML.parsefile(joinpath(field_script_outdir, "run_metadata.toml"))
    @test field_script_metadata["run"]["sampler"] == "gaussian_field_mc"
    @test field_script_metadata["sampling"]["estimator"] == "field"
    field_script_text = read(joinpath(field_script_outdir, "gaussian_mc_validation.md"), String)
    @test occursin("Gaussian field-amplitude", field_script_text)

    configuration_script_outdir = mktempdir()
    old_args = copy(ARGS)
    try
        empty!(ARGS)
        append!(ARGS, [
            "--outdir", configuration_script_outdir,
            "--nx", "8",
            "--c-values", "20.0",
            "--samples", "400",
            "--max-samples", "400",
            "--chains", "3",
            "--pressure-tolerance", "Inf",
            "--estimator", "field_configuration",
        ])
        redirect_stdout(devnull) do
            include(joinpath(@__DIR__, "..", "scripts", "run_gaussian_perturbation_mc.jl"))
        end
    finally
        empty!(ARGS)
        append!(ARGS, old_args)
    end
    configuration_script_metadata =
        TOML.parsefile(joinpath(configuration_script_outdir, "run_metadata.toml"))
    @test configuration_script_metadata["run"]["sampler"] == "gaussian_field_configuration_mc"
    @test configuration_script_metadata["sampling"]["estimator"] == "field_configuration"
    configuration_script_text =
        read(joinpath(configuration_script_outdir, "gaussian_mc_validation.md"), String)
    @test occursin("Gaussian field configurations", configuration_script_text)

    repeated_script_outdir = mktempdir()
    old_args = copy(ARGS)
    repeated_stderr_path = tempname()
    try
        empty!(ARGS)
        append!(ARGS, [
            "--outdir", repeated_script_outdir,
            "--nx", "4",
            "--c-values", "20.0",
            "--samples", "1",
            "--max-samples", "1",
            "--chains", "1",
            "--pressure-tolerance", "Inf",
        ])
        open(repeated_stderr_path, "w+") do stderr_io
            redirect_stdout(devnull) do
                redirect_stderr(stderr_io) do
                    include(joinpath(@__DIR__, "..", "scripts", "run_gaussian_perturbation_mc.jl"))
                    include(joinpath(@__DIR__, "..", "scripts", "run_gaussian_perturbation_mc.jl"))
                end
            end
        end
    finally
        empty!(ARGS)
        append!(ARGS, old_args)
    end
    @test !occursin("WARNING: Method definition", read(repeated_stderr_path, String))

    study_outdir = mktempdir()
    study = run_gaussian_sample_size_study(study_outdir; Nx=16, C=20.0,
        sample_counts=[25, 100, 400], chains=4, seed=17,
        pressure_tolerance=0.2)
    @test isfile(study.csv)
    @test isfile(study.text)
    @test study.config.Nx == 16
    @test study.config.C_values == [20.0]
    study_metadata = TOML.parsefile(joinpath(study_outdir, "run_metadata.toml"))
    @test study_metadata["run"]["mode"] == "gaussian_sample_size_study"
    @test study_metadata["sampling"]["sample_counts"] == [25, 100, 400]
    study_header, study_rows = DFMMonteCarlo._read_csv_numeric(study.csv)
    @test length(study_rows) == 3
    study_samples_idx = findfirst(==("samples"), study_header)
    study_chain_idx = findfirst(==("chain_count"), study_header)
    study_ess_idx = findfirst(==("effective_sample_size"), study_header)
    study_error_idx = findfirst(==("mc_error"), study_header)
    study_within_idx = findfirst(==("within_chain_stderr"), study_header)
    study_between_idx = findfirst(==("between_chain_stderr"), study_header)
    study_scaled_error_idx = findfirst(==("error_times_sqrt_effective_sample_size"), study_header)
    study_mc_minus_dfm_idx = findfirst(==("mc_minus_dfm"), study_header)
    study_reference_idx = findfirst(==("reference_error"), study_header)
    study_budget_idx = findfirst(==("sampling_budget_pass"), study_header)
    study_pass_idx = findfirst(==("diagnostic_pass"), study_header)
    @test [row[study_samples_idx] for row in study_rows] == [25.0, 100.0, 400.0]
    @test all(row[study_chain_idx] == 4.0 for row in study_rows)
    @test all(row[study_ess_idx] == row[study_samples_idx] * row[study_chain_idx] for row in study_rows)
    @test all(row[study_error_idx] ==
              max(row[study_within_idx], row[study_between_idx], 1.0e-12) for row in study_rows)
    @test all(row[study_scaled_error_idx] ≈
              row[study_error_idx] * sqrt(row[study_ess_idx]) for row in study_rows)
    @test all(abs(row[study_mc_minus_dfm_idx]) ≈ row[study_reference_idx] for row in study_rows)
    @test study_rows[1][study_budget_idx] == 0.0
    @test study_rows[end][study_budget_idx] == 1.0
    @test study_rows[end][study_pass_idx] == 1.0
    study_text = read(study.text, String)
    @test occursin("Gaussian sample-size convergence study", study_text)
    @test occursin("effective independent samples", study_text)
    @test occursin("1/sqrt(N)", study_text)

    configuration_study_outdir = mktempdir()
    configuration_study = run_gaussian_sample_size_study(configuration_study_outdir; Nx=16, C=20.0,
        sample_counts=[10, 40], chains=3, seed=47, pressure_tolerance=0.5,
        estimator=:field_configuration)
    configuration_study_metadata =
        TOML.parsefile(joinpath(configuration_study_outdir, "run_metadata.toml"))
    @test configuration_study_metadata["run"]["sampler"] == "gaussian_field_configuration_mc"
    @test configuration_study_metadata["sampling"]["estimator"] == "field_configuration"
    @test configuration_study_metadata["estimator"]["formula_id"] ==
          "finite_grid_gaussian_field_configuration_pressure"
    configuration_study_header, configuration_study_rows =
        DFMMonteCarlo._read_csv_strings(configuration_study.csv)
    configuration_study_sampler_idx = findfirst(==("sampler"), configuration_study_header)
    @test all(row[configuration_study_sampler_idx] == "gaussian_field_configuration_mc"
        for row in configuration_study_rows)
    configuration_study_text = read(configuration_study.text, String)
    @test occursin("field-configuration", configuration_study_text)

    study_script_outdir = mktempdir()
    old_args = copy(ARGS)
    try
        empty!(ARGS)
        append!(ARGS, [
            "--outdir", study_script_outdir,
            "--nx", "16",
            "--c", "20.0",
            "--sample-counts", "25,100",
            "--chains", "4",
            "--pressure-tolerance", "0.2",
        ])
        redirect_stdout(devnull) do
            include(joinpath(@__DIR__, "..", "scripts", "run_gaussian_sample_size_study.jl"))
        end
    finally
        empty!(ARGS)
        append!(ARGS, old_args)
    end
    @test isfile(joinpath(study_script_outdir, "gaussian_sample_size_study.csv"))
    @test isfile(joinpath(study_script_outdir, "gaussian_sample_size_study.md"))
    @test isfile(joinpath(study_script_outdir, "run_metadata.toml"))

    configuration_study_script_outdir = mktempdir()
    old_args = copy(ARGS)
    try
        empty!(ARGS)
        append!(ARGS, [
            "--outdir", configuration_study_script_outdir,
            "--nx", "16",
            "--c", "20.0",
            "--sample-counts", "10,40",
            "--chains", "3",
            "--pressure-tolerance", "0.5",
            "--estimator", "field_configuration",
        ])
        redirect_stdout(devnull) do
            include(joinpath(@__DIR__, "..", "scripts", "run_gaussian_sample_size_study.jl"))
        end
    finally
        empty!(ARGS)
        append!(ARGS, old_args)
    end
    configuration_study_script_metadata =
        TOML.parsefile(joinpath(configuration_study_script_outdir, "run_metadata.toml"))
    @test configuration_study_script_metadata["run"]["sampler"] == "gaussian_field_configuration_mc"
    @test configuration_study_script_metadata["sampling"]["estimator"] == "field_configuration"

    text_path = write_gaussian_mc_text(result, outdir)
    @test isfile(text_path)
    text = read(text_path, String)
    @test occursin("Gaussian HE-URP perturbation Monte Carlo", text)
    @test occursin("finite-grid DFM", text)
    @test occursin("diagnostic gate", text)
    @test occursin("reference-agreement", text)
    @test occursin("sample-budget", text)
    @test occursin("within-chain", text)
    @test occursin("between-chain", text)
    @test occursin("finite-grid DFM-vs-CL model offset", text)
    @test occursin("residual after subtracting", text)

    comparison_outdir = mktempdir()
    comparison_paths = run_gaussian_analysis_comparison_suite(
        comparison_outdir;
        Nx=16,
        max_points_per_figure=1,
        figure_names=["PivsC_B1_A300_50", "dPivsC_B1_A10_50", "PivsA_B1_C20"],
        samples=5000,
        max_samples=5000,
        chains=4,
        seed=31,
        pressure_tolerance=0.25,
    )
    @test length(comparison_paths) == 3
    comparison_metadata_path = joinpath(comparison_outdir, "run_metadata.toml")
    @test isfile(comparison_metadata_path)
    comparison_metadata = TOML.parsefile(comparison_metadata_path)
    @test comparison_metadata["run"]["sampler"] == "gaussian_mode_mc"
    @test comparison_metadata["run"]["mode"] == "gaussian_analysis_comparison"
    @test comparison_metadata["run"]["seed"] == 31
    @test comparison_metadata["grid"]["Nx"] == 16
    @test comparison_metadata["sampling"]["samples"] == 5000
    @test comparison_metadata["sampling"]["chains"] == 4
    @test comparison_metadata["outputs"]["comparison_csv_count"] == 3
    @test length(comparison_metadata["outputs"]["comparison_csvs"]) == 3
    for path in comparison_paths
        @test isfile(path)
        header, rows = DFMMonteCarlo._read_csv_strings(path)
        @test "mc_pressure" in header
        @test "mc_error" in header
        @test "mean_field" in header
        @test "one_loop" in header
        @test "dfm_pressure" in header
        @test "stable_mode_count" in header
        @test "samples" in header
        @test "chain_count" in header
        @test "within_chain_stderr" in header
        @test "between_chain_stderr" in header
        @test "chain_pressure_std" in header
        @test "effective_sample_size" in header
        @test "reference_error" in header
        @test "reference_zscore" in header
        @test "reference_agreement_pass" in header
        @test "mc_minus_dfm" in header
        @test "mc_dfm_zscore" in header
        @test "mc_minus_cl" in header
        @test "dfm_minus_cl" in header
        @test "cl_model_residual" in header
        @test "cl_context_pass" in header
        @test "required_effective_sample_size" in header
        @test "required_samples_per_chain" in header
        @test "sampling_budget_pass" in header
        @test "pressure_error_pass" in header
        @test "diagnostic_pass" in header
        @test "sampler" in header
        @test length(rows) == 1
        mc_idx = findfirst(==("mc_pressure"), header)
        err_idx = findfirst(==("mc_error"), header)
        dfm_idx = findfirst(==("dfm_pressure"), header)
        pass_idx = findfirst(==("diagnostic_pass"), header)
        sampler_idx = findfirst(==("sampler"), header)
        within_idx = findfirst(==("within_chain_stderr"), header)
        between_idx = findfirst(==("between_chain_stderr"), header)
        mc_minus_dfm_idx = findfirst(==("mc_minus_dfm"), header)
        cl_residual_idx = findfirst(==("cl_model_residual"), header)
        cl_context_idx = findfirst(==("cl_context_pass"), header)
        ref_pass_idx = findfirst(==("reference_agreement_pass"), header)
        mc = parse(Float64, rows[1][mc_idx])
        err = parse(Float64, rows[1][err_idx])
        dfm = parse(Float64, rows[1][dfm_idx])
        @test abs(mc - dfm) <= max(0.25, 5.0 * err)
        @test err == max(parse(Float64, rows[1][within_idx]),
            parse(Float64, rows[1][between_idx]), 1.0e-12)
        @test parse(Float64, rows[1][mc_minus_dfm_idx]) ≈ mc - dfm
        @test parse(Float64, rows[1][cl_residual_idx]) ≈ parse(Float64, rows[1][mc_minus_dfm_idx])
        @test parse(Float64, rows[1][cl_context_idx]) == parse(Float64, rows[1][ref_pass_idx])
        @test parse(Float64, rows[1][pass_idx]) == 1.0
        @test rows[1][sampler_idx] == "gaussian_mode_mc"
    end
    gaussian_svg = joinpath(comparison_outdir, "figures", "PivsC_B1_A300_50_gaussian_mc_comparison.svg")
    @test isfile(gaussian_svg)
    gaussian_svg_text = read(gaussian_svg, String)
    @test occursin("class='plot-title'", gaussian_svg_text)
    @test occursin("class='x-axis-label'", gaussian_svg_text)
    @test occursin("class='y-axis-label'", gaussian_svg_text)

    summary_path = write_gaussian_analysis_summary(comparison_outdir)
    @test isfile(summary_path)
    summary_text = read(summary_path, String)
    @test occursin("Gaussian Analysis MC validation", summary_text)
    @test occursin("finite-grid DFM", summary_text)
    @test occursin("3 comparison CSVs", summary_text)
    @test occursin("3 of 3 MC points passed", summary_text)
    @test occursin("maximum absolute MC deviation from finite-grid DFM", summary_text)
    @test occursin("maximum absolute MC deviation from CL", summary_text)
    @test occursin("finite-grid DFM-vs-CL model offset", summary_text)
    @test occursin("maximum residual after subtracting", summary_text)
    @test occursin("reference-agreement", summary_text)
    @test occursin("sample-budget", summary_text)
    @test occursin("within-chain", summary_text)
    @test occursin("between-chain", summary_text)
    @test occursin("Current status: manuscript-ready perturbation MC overlays", summary_text)

    field_comparison_outdir = mktempdir()
    field_comparison_paths = run_gaussian_analysis_comparison_suite(
        field_comparison_outdir;
        Nx=16,
        max_points_per_figure=1,
        figure_names=["PivsC_B1_A300_50"],
        samples=1000,
        max_samples=1000,
        chains=2,
        seed=53,
        pressure_tolerance=Inf,
        estimator=:field,
    )
    @test length(field_comparison_paths) == 1
    field_comparison_metadata = TOML.parsefile(joinpath(field_comparison_outdir, "run_metadata.toml"))
    @test field_comparison_metadata["run"]["sampler"] == "gaussian_field_mc"
    @test field_comparison_metadata["sampling"]["estimator"] == "field"
    @test field_comparison_metadata["estimator"]["formula_id"] == "finite_grid_gaussian_field_amplitude_pressure"
    field_header, field_rows = DFMMonteCarlo._read_csv_strings(only(field_comparison_paths))
    field_sampler_idx = findfirst(==("sampler"), field_header)
    field_pass_idx = findfirst(==("diagnostic_pass"), field_header)
    @test only(field_rows)[field_sampler_idx] == "gaussian_field_mc"
    @test parse(Float64, only(field_rows)[field_pass_idx]) == 1.0
    field_summary_path = write_gaussian_analysis_summary(field_comparison_outdir)
    field_summary_text = read(field_summary_path, String)
    @test occursin("field-amplitude", field_summary_text)
    @test !occursin("mode samples per chain", field_summary_text)
end

@testset "Slit confinement saddle validation" begin
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_slit_saddle_validation.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_slit_3d_validation.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_slit_3d_measure_corrected_scan.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_slit_3d_profile_mode_basis_ak_scan.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_slit_3d_profile_mode_a_ak_scan.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "write_slit_validation_text.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_slit_vk1_field_saddle_validation.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_alexander_katz_figures.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_alexander_katz_fig4_digitization.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "write_ak_fig4_b25_vk1_a_mapping.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "write_ak_fig4_b25_a_prediction.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "summarize_ak_fig4_b25_multiseed_anchor.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "summarize_ak_fig4_multiseed_anchor.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "write_ak_fig4_multiseed_mapping.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "write_ak_fig4_global_vk1_summary.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "write_ak_fig4_low_bc_branch_quadrature.jl"))
    @test isfile(joinpath(@__DIR__, "..", "scripts", "write_ak_fig4_b25_corrected_multiseed_mapping.jl"))

    grid = SlitGrid1D(40, 1.0)
    @test grid.Nz == 40
    @test grid.L == 1.0
    @test length(grid.z) == 40
    @test length(grid.weights) == 40
    @test DFMMonteCarlo.slit_integral(ones(40), grid) ≈ 1.0
    @test DFMMonteCarlo.normalize_slit_density(fill(2.0, 40), grid) ≈ ones(40)

    phi_low = DFMMonteCarlo.slit_low_density_profile(grid)
    @test DFMMonteCarlo.slit_integral(phi_low, grid) ≈ grid.L atol = 2.0e-3
    @test phi_low[1] < phi_low[20]
    @test phi_low[end] < phi_low[20]

    reference = DFMMonteCarlo.slit_reference_profile(grid.z; B=25.0, C=0.5, L=1.0, xi=0.37)
    fitted_xi, fit_error = DFMMonteCarlo.fit_slit_xi_eff(grid.z, reference; L=1.0)
    @test fitted_xi ≈ 0.37 atol = 5.0e-3
    @test fit_error < 1.0e-6

    sine_mode = sin.(pi .* grid.z ./ grid.L)
    operated = DFMMonteCarlo.apply_slit_debye_inverse(sine_mode, grid; alpha=0.0)
    expected = sine_mode ./ debye((pi / grid.L)^2)
    @test maximum(abs.(operated .- expected)) < 2.0e-2

    params = SlitParams(B=1.0, C=10.0, A=0.0, alpha=0.0, model=:urp_ti)
    uniform = ones(grid.Nz)
    @test isfinite(slit_heurp_ti_energy(uniform, grid, params))
    @test isfinite(slit_hevk1_energy(uniform, grid, params))

    wall_grid = SlitGrid1D(4, 4.0)
    unit_density = ones(wall_grid.Nz)
    unit_k1 = ones(wall_grid.Nz)
    @test DFMMonteCarlo._slit_vk1_gradient_integral(unit_density, unit_k1, wall_grid) ≈ 8.0

    urp_result = minimize_slit_density(:urp_ti, 1.0, 10.0, 0.0; Nz=40, L=1.0, alpha=0.0)
    vk1_result = minimize_slit_density(:vk1, 1.0, 10.0, 0.0; Nz=40, L=1.0, alpha=0.0)
    @test urp_result.model == :urp_ti
    @test vk1_result.model == :vk1
    @test DFMMonteCarlo.slit_integral(urp_result.phi, urp_result.grid) ≈ 1.0 atol = 1.0e-8
    @test DFMMonteCarlo.slit_integral(vk1_result.phi, vk1_result.grid) ≈ 1.0 atol = 1.0e-8
    @test minimum(urp_result.phi) > 0.0
    @test minimum(vk1_result.phi) > 0.0
    @test isfinite(urp_result.energy)
    @test isfinite(vk1_result.energy)
    @test urp_result.xi_eff > 0.0
    @test vk1_result.xi_eff > 0.0

    field_result = minimize_slit_vk1_field_saddle(25.0, 0.5, 0.0; Nz=40, L=1.0, max_iterations=120)
    ansatz_result = minimize_slit_density(:vk1, 25.0, 0.5, 0.0; Nz=40, L=1.0, alpha=0.0)
    @test field_result.model == :vk1
    @test field_result.solver == "full_field_log_density_lbfgs"
    @test DFMMonteCarlo.slit_integral(field_result.phi, field_result.grid) ≈ field_result.grid.L atol = 1.0e-10
    @test minimum(field_result.phi) > 0.0
    @test isfinite(field_result.energy)
    @test field_result.energy <= ansatz_result.energy + 1.0e-5
    @test field_result.rms_reference_error < 5.0e-2
    @test isfinite(field_result.projected_gradient_rms)
    @test field_result.projected_gradient_rms < 1.0

    smooth_result = minimize_slit_vk1_smooth_basis_saddle(25.0, 0.5, 0.0;
        Nz=40, L=1.0, basis_count=8, max_iterations=1000)
    @test smooth_result.model == :vk1
    @test startswith(smooth_result.solver, "smooth_log_density_cosine_nelder_mead_")
    @test smooth_result.converged
    @test DFMMonteCarlo.slit_integral(smooth_result.phi, smooth_result.grid) ≈ smooth_result.grid.L atol = 1.0e-10
    @test minimum(smooth_result.phi) > 0.0
    @test isfinite(smooth_result.energy)
    @test smooth_result.energy <= ansatz_result.energy + 1.0e-5
    @test smooth_result.rms_reference_error < 7.5e-2
    @test isfinite(smooth_result.projected_gradient_rms)

    outdir = mktempdir()
    result = run_slit_saddle_validation(outdir; Nz=40, B_values=[1.0, 25.0],
        C_values=[0.5, 10.0], A_values=[0.0, 299.0], alpha=0.0)
    @test isfile(result.profile_csv)
    @test isfile(result.summary_csv)
    @test isfile(result.profile_svg)
    @test isfile(result.xi_svg)
    @test isfile(result.markdown_path)
    validation_text = read(result.markdown_path, String)
    @test occursin("HE-VK1 as the current real density Hamiltonian", validation_text)
    @test occursin("diagnostic bulk-kernel extension", validation_text)
    summary_header, summary_rows = DFMMonteCarlo._read_csv_strings(result.summary_csv)
    @test "model" in summary_header
    @test "A" in summary_header
    @test "xi_eff" in summary_header
    @test "status" in summary_header
    @test length(summary_rows) >= 8
    profile_svg = read(result.profile_svg, String)
    xi_svg = read(result.xi_svg, String)
    @test occursin("Slit density profiles", profile_svg)
    @test occursin("z / L", profile_svg)
    @test occursin("phi(z)", profile_svg)
    @test occursin("Effective slit correlation length", xi_svg)
    @test occursin("B C", xi_svg)
    @test occursin("xi_eff", xi_svg)

    field_outdir = mktempdir()
    field_validation = run_slit_vk1_field_saddle_validation(field_outdir; B=25.0, C=0.5,
        A=0.0, Nz=40, L=1.0, max_iterations=80, smooth_basis_count=6,
        smooth_max_iterations=1000)
    @test isfile(field_validation.profile_csv)
    @test isfile(field_validation.summary_csv)
    field_profile_header, field_profile_rows = DFMMonteCarlo._read_csv_strings(field_validation.profile_csv)
    field_profile_idx = findfirst(==("profile"), field_profile_header)
    @test any(row[field_profile_idx] == "vk1_smooth_basis" for row in field_profile_rows)
    field_summary_header, field_summary_rows = DFMMonteCarlo._read_csv_strings(field_validation.summary_csv)
    field_solver_idx = findfirst(==("solver"), field_summary_header)
    field_energy_idx = findfirst(==("energy"), field_summary_header)
    @test length(field_summary_rows) == 1
    @test any(startswith(row[field_solver_idx], "smooth_log_density_cosine_nelder_mead_") for row in field_summary_rows)

    fig2_series_script = joinpath(@__DIR__, "..", "scripts", "write_ak_fig2_vk1_saddle_series.jl")
    @test isfile(fig2_series_script)
    fig2_series_outdir = mktempdir()
    old_args = copy(ARGS)
    empty!(ARGS)
    append!(ARGS, [
        "--outdir=$fig2_series_outdir",
        "--nz=32",
        "--smooth-basis-count=4",
        "--smooth-max-iterations=300",
        "--bc-values=1,10",
    ])
    try
        include(fig2_series_script)
    finally
        empty!(ARGS)
        append!(ARGS, old_args)
    end
    @test isfile(joinpath(fig2_series_outdir, "ak_fig2_vk1_saddle_profiles.csv"))
    @test isfile(joinpath(fig2_series_outdir, "ak_fig2_vk1_saddle_summary.csv"))
    @test isfile(joinpath(fig2_series_outdir, "ak_fig2_vk1_saddle_profiles.svg"))
    fig2_series_header, fig2_series_rows =
        DFMMonteCarlo._read_csv_strings(joinpath(fig2_series_outdir, "ak_fig2_vk1_saddle_summary.csv"))
    fig2_series_bc_idx = findfirst(==("BC"), fig2_series_header)
    fig2_series_solver_idx = findfirst(==("solver"), fig2_series_header)
    @test sort(unique(parse(Float64, row[fig2_series_bc_idx]) for row in fig2_series_rows)) == [1.0, 10.0]
    @test all(startswith(row[fig2_series_solver_idx], "smooth_log_density_cosine_nelder_mead_") for row in fig2_series_rows)

    field_diagnostic_outdir = mktempdir()
    field_diagnostic = run_slit_vk1_field_saddle_validation(field_diagnostic_outdir; B=25.0, C=0.5,
        A=0.0, Nz=40, L=1.0, max_iterations=80, smooth_basis_count=6,
        smooth_max_iterations=1000, include_full_field=true)
    diagnostic_header, diagnostic_rows = DFMMonteCarlo._read_csv_strings(field_diagnostic.summary_csv)
    diagnostic_solver_idx = findfirst(==("solver"), diagnostic_header)
    diagnostic_energy_idx = findfirst(==("energy"), diagnostic_header)
    @test length(diagnostic_rows) >= 2
    @test any(row[diagnostic_solver_idx] == "full_field_log_density_lbfgs" for row in diagnostic_rows)
    field_energy = only(parse(Float64, row[diagnostic_energy_idx]) for row in diagnostic_rows if row[diagnostic_solver_idx] == "full_field_log_density_lbfgs")
    smooth_energy = only(parse(Float64, row[diagnostic_energy_idx]) for row in diagnostic_rows if startswith(row[diagnostic_solver_idx], "smooth_log_density_cosine_nelder_mead_"))
    @test field_energy <= smooth_energy + 1.0e-5

    noise = run_slit_noise_profile_mc(25.0, 0.5, 0.0; Nz=40, basis_count=4,
        steps=8_000, burnin_steps=2_000, sample_stride=20,
        proposal_width=0.03, temperature=1.0, seed=91)
    @test noise.sample_count > 0
    @test 0.0 < noise.acceptance_rate < 1.0
    @test DFMMonteCarlo.slit_integral(noise.mean_phi, noise.grid) ≈ noise.grid.L atol = 1.0e-10
    @test maximum(noise.mean_phi) > maximum(noise.saddle_phi)
    @test noise.xi_eff > noise.xi_bulk

    ak_outdir = mktempdir()
    ak_validation = run_alexander_katz_figure_validation(ak_outdir; Nz=80,
        fig3_mc_steps=8_000, fig3_mc_burnin_steps=2_000, fig3_mc_sample_stride=20,
        fig3_mc_seed=91)
    @test isfile(ak_validation.fig2_csv)
    @test isfile(ak_validation.fig2_svg)
    @test isfile(ak_validation.fig3_csv)
    @test isfile(ak_validation.fig3_svg)
    @test isfile(ak_validation.summary_csv)
    @test isfile(ak_validation.markdown_path)
    ak_fig3_header, ak_fig3_rows = DFMMonteCarlo._read_csv_strings(ak_validation.fig3_csv)
    ak_series_idx = findfirst(==("series"), ak_fig3_header)
    ak_phi_idx = findfirst(==("phi"), ak_fig3_header)
    @test any(row[ak_series_idx] == "mc_noise_profile" for row in ak_fig3_rows)
    @test any(row[ak_series_idx] == "mc_noise_points" for row in ak_fig3_rows)
    @test any(row[ak_series_idx] == "mc_eq22_fit" for row in ak_fig3_rows)
    @test !any(row[ak_series_idx] == "fluctuation_fit" for row in ak_fig3_rows)
    @test !any(row[ak_series_idx] == "fluctuation_points" for row in ak_fig3_rows)
    mc_peak = maximum(parse(Float64, row[ak_phi_idx]) for row in ak_fig3_rows if row[ak_series_idx] == "mc_noise_profile")
    fit_peak = maximum(parse(Float64, row[ak_phi_idx]) for row in ak_fig3_rows if row[ak_series_idx] == "mc_eq22_fit")
    @test fit_peak > mc_peak
    ak_summary_header, ak_summary_rows = DFMMonteCarlo._read_csv_strings(ak_validation.summary_csv)
    ak_metric_idx = findfirst(==("metric"), ak_summary_header)
    ak_value_idx = findfirst(==("value"), ak_summary_header)
    summary = Dict(row[ak_metric_idx] => parse(Float64, row[ak_value_idx]) for row in ak_summary_rows)
    @test summary["fig2_bc_count"] == 4.0
    @test summary["fig2_min_bc"] == 1.0
    @test summary["fig2_max_bc"] == 1000.0
    @test summary["fig3_B"] == 25.0
    @test summary["fig3_C"] == 0.5
    @test summary["fig3_BC"] == 12.5
    @test summary["fig3_xi_bulk"] ≈ 0.2
    @test summary["fig3_noise_source_code"] == 1.0
    @test summary["fig3_mc_sample_count"] > 0.0
    @test 0.0 < summary["fig3_mc_acceptance_rate"] < 1.0
    @test summary["fig3_mc_xi_eff"] > summary["fig3_xi_bulk"]
    ak_text = read(ak_validation.markdown_path, String)
    @test occursin("Alexander-Katz Fig. 2", ak_text)
    @test occursin("Alexander-Katz Fig. 3", ak_text)
    @test occursin("Monte Carlo", ak_text)
    @test occursin("effective correlation length", ak_text)
    @test !occursin("not a new complex-Langevin or Monte Carlo", ak_text)

    digitization_outdir = mktempdir()
    digitization = run_alexander_katz_fig3_digitization(digitization_outdir;
        normalization_Nz=400)
    @test isfile(digitization.profile_csv)
    @test isfile(digitization.summary_csv)
    @test isfile(digitization.profile_svg)
    @test isfile(digitization.markdown_path)
    @test isfile(digitization.metadata_path)
    @test digitization.point_count >= 30
    @test 0.30 <= digitization.xi_eff <= 0.36
    @test digitization.xi_fit_r2 > 0.995
    digitization_header, digitization_rows =
        DFMMonteCarlo._read_csv_strings(digitization.profile_csv)
    @test "pixel_x" in digitization_header
    @test "pixel_y" in digitization_header
    @test "z_over_L" in digitization_header
    @test "phi" in digitization_header
    @test "fit_phi" in digitization_header
    @test length(digitization_rows) == digitization.point_count
    digitization_summary_header, digitization_summary_rows =
        DFMMonteCarlo._read_csv_strings(digitization.summary_csv)
    digitization_metric_idx = findfirst(==("metric"), digitization_summary_header)
    digitization_value_idx = findfirst(==("value"), digitization_summary_header)
    digitization_summary = Dict(row[digitization_metric_idx] => parse(Float64, row[digitization_value_idx])
        for row in digitization_summary_rows)
    @test digitization_summary["plot_frame_left_px"] == 61.0
    @test digitization_summary["plot_frame_right_px"] == 797.0
    @test digitization_summary["plot_frame_top_px"] == 47.0
    @test digitization_summary["plot_frame_bottom_px"] == 670.0
    digitization_metadata = TOML.parsefile(digitization.metadata_path)
    @test digitization_metadata["run"]["mode"] == "alexander_katz_fig3_digitization"
    @test digitization_metadata["fit"]["normalization_Nz"] == 400
    digitization_text = read(digitization.markdown_path, String)
    @test occursin("Alexander-Katz Fig. 3", digitization_text)
    @test occursin("digitized filled-square", digitization_text)

    fig4_outdir = mktempdir()
    fig4_digitization = run_alexander_katz_fig4_digitization(fig4_outdir)
    @test isfile(fig4_digitization.xi_csv)
    @test isfile(fig4_digitization.summary_csv)
    @test isfile(fig4_digitization.xi_svg)
    @test isfile(fig4_digitization.markdown_path)
    @test isfile(fig4_digitization.metadata_path)
    @test fig4_digitization.point_count >= 20
    fig4_header, fig4_rows = DFMMonteCarlo._read_csv_strings(fig4_digitization.xi_csv)
    @test "B" in fig4_header
    @test "BC" in fig4_header
    @test "C" in fig4_header
    @test "xi_eff" in fig4_header
    @test "xi_bulk" in fig4_header
    fig4_b_idx = findfirst(==("B"), fig4_header)
    fig4_bc_idx = findfirst(==("BC"), fig4_header)
    fig4_xi_idx = findfirst(==("xi_eff"), fig4_header)
    fig4_bulk_idx = findfirst(==("xi_bulk"), fig4_header)
    fig4_b_values = sort(unique(parse(Float64, row[fig4_b_idx]) for row in fig4_rows))
    @test fig4_b_values == [1.0, 10.0, 25.0]
    fig4_b25_anchor = [
        row for row in fig4_rows
        if parse(Float64, row[fig4_b_idx]) == 25.0 &&
           12.0 <= parse(Float64, row[fig4_bc_idx]) <= 13.0
    ]
    @test length(fig4_b25_anchor) == 1
    @test parse(Float64, only(fig4_b25_anchor)[fig4_xi_idx]) ≈ 0.35 atol = 0.02
    fig4_b25_targets = DFMMonteCarlo._read_ak_fig4_b25_digitized_targets(fig4_digitization.xi_csv)
    @test length(fig4_b25_targets) == 9
    @test DFMMonteCarlo._ak_fig4_loglog_target_xi(fig4_b25_targets, 12.5) ≈ 0.3496 atol = 0.003
    @test DFMMonteCarlo._ak_fig4_loglog_target_xi(fig4_b25_targets, 16.482431389109692) ≈
          0.2804808303704316 atol = 1.0e-10
    fig4_b10_targets = DFMMonteCarlo._read_ak_fig4_digitized_targets(fig4_digitization.xi_csv; B=10.0)
    @test length(fig4_b10_targets) == 10
    @test DFMMonteCarlo._ak_fig4_loglog_target_xi(fig4_b10_targets, 10.0) ≈
          0.25320705975988633 atol = 1.0e-12
    @test any(parse(Float64, row[fig4_b_idx]) == 25.0 &&
              40.0 <= parse(Float64, row[fig4_bc_idx]) <= 80.0 &&
              parse(Float64, row[fig4_xi_idx]) > parse(Float64, row[fig4_bulk_idx])
        for row in fig4_rows)
    fig4_summary_header, fig4_summary_rows = DFMMonteCarlo._read_csv_strings(fig4_digitization.summary_csv)
    fig4_metric_idx = findfirst(==("metric"), fig4_summary_header)
    fig4_value_idx = findfirst(==("value"), fig4_summary_header)
    fig4_summary = Dict(row[fig4_metric_idx] => parse(Float64, row[fig4_value_idx])
        for row in fig4_summary_rows)
    @test fig4_summary["plot_frame_left_px"] == 198.5
    @test fig4_summary["x_decade_px"] == 327.5
    @test fig4_summary["y_decade_px"] == 425.8
    fig4_metadata = TOML.parsefile(fig4_digitization.metadata_path)
    @test fig4_metadata["run"]["mode"] == "alexander_katz_fig4_digitization"
    fig4_text = read(fig4_digitization.markdown_path, String)
    @test occursin("Alexander-Katz Fig. 4", fig4_text)
    @test occursin("xi_eff", fig4_text)

    quoted_csv_outdir = mktempdir()
    quoted_rows = [
        (A=0.0, temperature=1.0, target_xi_eff=0.22,
            target_source="Alexander-Katz Fig. 4, B=25, BC=57.77",
            xi_eff=0.18, abs_xi_error=0.04, xi_tolerance=0.05,
            calibration_pass=0.0, xi_bulk=0.09, xi_fit_rmse=0.01,
            xi_fit_r2=0.99, xi_fit_quality_pass=1.0,
            target_profile_rmse=0.02, target_profile_r2=0.996,
            target_profile_max_abs_error=0.03, target_profile_wall_rms=0.01,
            target_profile_r2_threshold=0.995, profile_chain_rhat_max=1.01,
            profile_split_rhat_max=1.02, profile_chain_rmse_max=0.001,
            acceptance_rate=0.2, sample_count=10.0, steps=100.0,
            burnin_steps=20.0, sample_stride=5.0, chain_count=1.0,
            z_basis_count=4, proposal_width=0.8, energy_mean=1.0,
            energy_stderr=0.1, density_variance_mean=0.1,
            lateral_variance_mean=0.0, run_outdir="run0",
            profile_csv="profile0.csv", diagnostics_csv="diag0.csv",
            profile_svg="profile0.svg", selection_rank=2),
        (A=1.0, temperature=1.0, target_xi_eff=0.22,
            target_source="Alexander-Katz Fig. 4, B=25, BC=57.77",
            xi_eff=0.221, abs_xi_error=0.001, xi_tolerance=0.05,
            calibration_pass=1.0, xi_bulk=0.09, xi_fit_rmse=0.01,
            xi_fit_r2=0.99, xi_fit_quality_pass=1.0,
            target_profile_rmse=0.02, target_profile_r2=0.996,
            target_profile_max_abs_error=0.03, target_profile_wall_rms=0.01,
            target_profile_r2_threshold=0.995, profile_chain_rhat_max=1.01,
            profile_split_rhat_max=1.02, profile_chain_rmse_max=0.001,
            acceptance_rate=0.2, sample_count=10.0, steps=100.0,
            burnin_steps=20.0, sample_stride=5.0, chain_count=1.0,
            z_basis_count=4, proposal_width=0.8, energy_mean=1.0,
            energy_stderr=0.1, density_variance_mean=0.1,
            lateral_variance_mean=0.0, run_outdir="run1",
            profile_csv="profile1.csv", diagnostics_csv="diag1.csv",
            profile_svg="profile1.svg", selection_rank=1),
    ]
    quoted_csv = DFMMonteCarlo._write_dataframe(joinpath(quoted_csv_outdir,
        "profile_mode_a_scan.csv"), quoted_rows)
    quoted_header, quoted_read_rows = DFMMonteCarlo._read_csv_strings(quoted_csv)
    quoted_source_idx = findfirst(==("target_source"), quoted_header)
    quoted_xi_idx = findfirst(==("xi_eff"), quoted_header)
    @test quoted_read_rows[1][quoted_source_idx] == "Alexander-Katz Fig. 4, B=25, BC=57.77"
    @test parse(Float64, quoted_read_rows[2][quoted_xi_idx]) ≈ 0.221
    quoted_svg = DFMMonteCarlo._write_slit_3d_profile_mode_a_scan_svg(quoted_csv,
        quoted_csv_outdir; target_label="Alexander-Katz Fig. 4")
    @test occursin("Alexander-Katz Fig. 4 target xi", read(quoted_svg, String))

    target_label_outdir = mktempdir()
    markdown_path = DFMMonteCarlo._write_slit_3d_profile_mode_a_scan_text(
        target_label_outdir, quoted_rows, quoted_rows[2],
        "Alexander-Katz Fig. 4, B=25, BC=57.77";
        target_label="Alexander-Katz Fig. 4")
    @test basename(markdown_path) == "ak_fig4_3d_profile_mode_a_scan.md"
    target_label_text = read(markdown_path, String)
    @test occursin("Alexander-Katz Fig. 4", target_label_text)
    @test !occursin("AK Fig. 3 A scan", target_label_text)

    mapping_prod1 = mktempdir()
    mapping_prod2 = mktempdir()
    mapping_z = collect(range(1 / 16, 15 / 16; length=8))
    mapping_phi1 = DFMMonteCarlo.slit_reference_profile(mapping_z; B=25.0, C=1.0,
        L=1.0, xi=0.31)
    mapping_phi2 = DFMMonteCarlo.slit_reference_profile(mapping_z; B=25.0, C=2.0,
        L=1.0, xi=0.22)
    DFMMonteCarlo._write_dataframe(joinpath(mapping_prod1, "slit_3d_mc_diagnostics.csv"), [
        (sampler="profile_mode_metropolis", A=22.5, B=25.0, C=1.0,
            Nx=20, Ny=20, Nz=8, Lxy=1.0, Lz=1.0, steps=1000,
            burnin_steps=200, sample_stride=10, z_basis_count=4,
            proposal_width=0.8, temperature=1.0, chain_count=2,
            sample_count=160, profile_estimator="profile_mode_mirror_symmetrized",
            profile_symmetrized=true, attempted_moves=2000, accepted_moves=600,
            acceptance_rate=0.3, energy_mean=1.0, energy_stderr=0.1,
            density_variance_mean=0.1, density_variance_stderr=0.01,
            lateral_variance_mean=0.0, lateral_variance_stderr=0.0,
            xy_profile_roughness_mean=0.0, xi_bulk=0.1414, xi_eff=0.31,
            xi_fit_error=0.0, xi_fit_rmse=0.0, xi_fit_r2=1.0,
            xi_fit_mae=0.0, xi_fit_max_abs_error=0.0, xi_fit_wall_rms=0.0,
            xi_fit_nrmse_range=0.0, xi_fit_nrmse_mean=0.0,
            xi_fit_quality_pass=1.0, profile_chain_rhat_max=1.01,
            profile_split_rhat_max=1.02, profile_chain_rmse_max=0.001,
            profile_chain_max_abs_deviation=0.002),
    ])
    DFMMonteCarlo._write_dataframe(joinpath(mapping_prod2, "slit_3d_mc_diagnostics.csv"), [
        (sampler="profile_mode_metropolis", A=115.0, B=25.0, C=2.0,
            Nx=20, Ny=20, Nz=8, Lxy=1.0, Lz=1.0, steps=1000,
            burnin_steps=200, sample_stride=10, z_basis_count=4,
            proposal_width=0.8, temperature=1.0, chain_count=2,
            sample_count=160, profile_estimator="profile_mode_mirror_symmetrized",
            profile_symmetrized=true, attempted_moves=2000, accepted_moves=500,
            acceptance_rate=0.25, energy_mean=1.0, energy_stderr=0.1,
            density_variance_mean=0.1, density_variance_stderr=0.01,
            lateral_variance_mean=0.0, lateral_variance_stderr=0.0,
            xy_profile_roughness_mean=0.0, xi_bulk=0.1, xi_eff=0.22,
            xi_fit_error=0.0, xi_fit_rmse=0.0, xi_fit_r2=1.0,
            xi_fit_mae=0.0, xi_fit_max_abs_error=0.0, xi_fit_wall_rms=0.0,
            xi_fit_nrmse_range=0.0, xi_fit_nrmse_mean=0.0,
            xi_fit_quality_pass=1.0, profile_chain_rhat_max=1.01,
            profile_split_rhat_max=1.02, profile_chain_rmse_max=0.001,
            profile_chain_max_abs_deviation=0.002),
    ])
    DFMMonteCarlo._write_dataframe(joinpath(mapping_prod1, "slit_3d_mc_profile.csv"), [
        (series="mc_mean_3d", B=25.0, C=1.0, A=22.5, Nx=20, Ny=20,
            Nz=8, Lxy=1.0, Lz=1.0, z_over_L=mapping_z[i],
            phi=mapping_phi1[i], stderr=0.0, xi=0.31,
            sampler="profile_mode_mirror_symmetrized") for i in eachindex(mapping_z)
    ])
    DFMMonteCarlo._write_dataframe(joinpath(mapping_prod2, "slit_3d_mc_profile.csv"), [
        (series="mc_mean_3d", B=25.0, C=2.0, A=115.0, Nx=20, Ny=20,
            Nz=8, Lxy=1.0, Lz=1.0, z_over_L=mapping_z[i],
            phi=mapping_phi2[i], stderr=0.0, xi=0.22,
            sampler="profile_mode_mirror_symmetrized") for i in eachindex(mapping_z)
    ])
    mapping_outdir = mktempdir()
    mapping_result = write_ak_fig4_b25_vk1_a_mapping(mapping_outdir; targets=[
        (tag="BC25", target_xi_eff=0.31, tolerance=0.02,
            production_outdir=mapping_prod1),
        (tag="BC50", target_xi_eff=0.22, tolerance=0.02,
            production_outdir=mapping_prod2),
    ])
    @test mapping_result.point_count == 2
    @test isfile(mapping_result.summary_csv)
    @test isfile(mapping_result.a_svg)
    @test isfile(mapping_result.xi_svg)
    @test isfile(mapping_result.markdown_path)
    mapping_header, mapping_rows = DFMMonteCarlo._read_csv_strings(mapping_result.summary_csv)
    @test "target_profile_r2" in mapping_header
    @test "calibration_status" in mapping_header
    mapping_status_idx = findfirst(==("calibration_status"), mapping_header)
    @test all(row[mapping_status_idx] == "full_profile_pass" for row in mapping_rows)
    mapping_text = read(mapping_result.markdown_path, String)
    @test occursin("AK Fig. 4 B=25 VK1 A mapping", mapping_text)
    @test occursin("Passing rows are calibration evidence", mapping_text)

    prediction_mapping_csv = DFMMonteCarlo._write_dataframe(joinpath(mktempdir(),
        "mapping.csv"), [
        (tag="BC10", B=25.0, BC=10.0, C=0.4, target_xi_eff=0.5,
            xi_bulk=0.2236, A=2.0, tolerance=0.02, production_xi_eff=0.5,
            abs_xi_error=0.0, rel_xi_error=0.0, fit_rmse=0.0, fit_r2=1.0,
            target_profile_rmse=0.0, target_profile_r2=1.0,
            target_profile_r2_threshold=0.995, profile_chain_rhat_max=1.0,
            profile_split_rhat_max=1.0, acceptance_rate=0.4,
            sample_count=100.0, calibration_status="full_profile_pass",
            production_outdir="run10"),
        (tag="BC100", B=25.0, BC=100.0, C=4.0, target_xi_eff=0.15,
            xi_bulk=0.0707, A=20.0, tolerance=0.02, production_xi_eff=0.15,
            abs_xi_error=0.0, rel_xi_error=0.0, fit_rmse=0.0, fit_r2=1.0,
            target_profile_rmse=0.0, target_profile_r2=1.0,
            target_profile_r2_threshold=0.995, profile_chain_rhat_max=1.0,
            profile_split_rhat_max=1.0, acceptance_rate=0.3,
            sample_count=100.0, calibration_status="full_profile_pass",
            production_outdir="run100"),
    ])
    prediction_targets_csv = DFMMonteCarlo._write_dataframe(joinpath(mktempdir(),
        "fig4.csv"), [
        (source="ak_fig4", series="B25_filled_square", B=25.0,
            marker="filled_square", pixel_x=0.0, pixel_y=0.0,
            log10_BC=log10(5.0), BC=5.0, C=0.2,
            log10_xi_eff=log10(0.8), xi_eff=0.8, xi_bulk=1 / sqrt(10.0),
            xi_ratio=0.8 * sqrt(10.0)),
        (source="ak_fig4", series="B25_filled_square", B=25.0,
            marker="filled_square", pixel_x=0.0, pixel_y=0.0,
            log10_BC=log10(10.0), BC=10.0, C=0.4,
            log10_xi_eff=log10(0.5), xi_eff=0.5, xi_bulk=1 / sqrt(20.0),
            xi_ratio=0.5 * sqrt(20.0)),
        (source="ak_fig4", series="B25_filled_square", B=25.0,
            marker="filled_square", pixel_x=0.0, pixel_y=0.0,
            log10_BC=log10(30.0), BC=30.0, C=1.2,
            log10_xi_eff=log10(0.3), xi_eff=0.3, xi_bulk=1 / sqrt(60.0),
            xi_ratio=0.3 * sqrt(60.0)),
        (source="ak_fig4", series="B25_filled_square", B=25.0,
            marker="filled_square", pixel_x=0.0, pixel_y=0.0,
            log10_BC=log10(300.0), BC=300.0, C=12.0,
            log10_xi_eff=log10(0.08), xi_eff=0.08, xi_bulk=1 / sqrt(600.0),
            xi_ratio=0.08 * sqrt(600.0)),
    ])
    prediction_outdir = mktempdir()
    prediction_result = write_ak_fig4_b25_a_prediction(prediction_outdir;
        mapping_csv=prediction_mapping_csv, fig4_xi_csv=prediction_targets_csv)
    @test isfile(prediction_result.prediction_csv)
    @test isfile(prediction_result.prediction_svg)
    @test isfile(prediction_result.markdown_path)
    prediction_header, prediction_rows =
        DFMMonteCarlo._read_csv_strings(prediction_result.prediction_csv)
    prediction_bc_idx = findfirst(==("BC"), prediction_header)
    prediction_A_idx = findfirst(==("predicted_A"), prediction_header)
    prediction_status_idx = findfirst(==("prediction_status"), prediction_header)
    prediction_by_bc = Dict(parse(Float64, row[prediction_bc_idx]) => row for row in prediction_rows)
    @test prediction_by_bc[10.0][prediction_status_idx] == "calibrated"
    @test prediction_by_bc[30.0][prediction_status_idx] == "interpolation"
    @test prediction_by_bc[5.0][prediction_status_idx] == "low_extrapolation"
    @test prediction_by_bc[300.0][prediction_status_idx] == "high_extrapolation"
    @test parse(Float64, prediction_by_bc[30.0][prediction_A_idx]) ≈ 6.0 atol = 1.0e-10
end

@testset "3D slit real-density HE-VK1 Metropolis sampler" begin
    grid3 = SlitGrid3D(4, 8, 1.0; Lxy=1.0)
    @test grid3.Nx == 4
    @test grid3.Ny == 4
    @test grid3.Nz == 8
    @test grid3.volume ≈ 1.0
    @test length(grid3.z) == 8
    @test DFMMonteCarlo._proposals_per_sweep(grid3) == 4 * 4 * 8
    @test DFMMonteCarlo._steps_from_sweeps(grid3, 0.5) == 64
    @test DFMMonteCarlo._sweeps_from_steps(grid3, 64) ≈ 0.5

    profile = DFMMonteCarlo.slit_reference_profile(grid3.z; B=25.0, C=0.5, L=1.0)
    field = DFMMonteCarlo.slit_3d_field_from_profile(profile, grid3)
    params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0, model=:vk1)
    @test size(field) == (4, 4, 8)
    @test mean(field) ≈ 1.0 atol = 1.0e-12
    @test slit_hevk1_energy_3d(field, grid3, params) ≈
          slit_hevk1_energy(profile, grid3.slit_grid, params) atol = 1.0e-10
    profile_log_jacobian = DFMMonteCarlo._slit_3d_profile_log_jacobian(field, grid3)
    expected_log_jacobian = (grid3.Nx * grid3.Ny - 1) * sum(log.(profile))
    @test profile_log_jacobian ≈ expected_log_jacobian atol = 1.0e-10
    @test DFMMonteCarlo._slit_3d_profile_log_jacobian(ones(size(field)), grid3) ≈ 0.0 atol = 1.0e-12
    @test DFMMonteCarlo.slit_hevk1_energy_3d_profile_measure(field, grid3, params;
        profile_measure=:flat_cell) ≈ slit_hevk1_energy_3d(field, grid3, params)
    @test DFMMonteCarlo.slit_hevk1_energy_3d_profile_measure(field, grid3, params;
        profile_measure=:flat_profile) ≈
          slit_hevk1_energy_3d(field, grid3, params) + expected_log_jacobian
    @test DFMMonteCarlo.slit_hevk1_energy_3d_blocked_profile_measure(field, grid3, params;
        profile_measure_power=1.0) ≈ slit_hevk1_energy_3d(field, grid3, params)
    @test DFMMonteCarlo.slit_hevk1_energy_3d_blocked_profile_measure(field, grid3, params;
        profile_measure_power=1.25) ≈
          slit_hevk1_energy_3d(field, grid3, params) + 0.25 * expected_log_jacobian
    lateral_field = copy(field)
    for i in 1:grid3.Nx, j in 1:grid3.Ny, k in 1:grid3.Nz
        lateral_field[i, j, k] *= 1.0 + 0.10 * (isodd(i + j) ? 1.0 : -1.0)
    end
    @test DFMMonteCarlo._slit_3d_xy_profile(lateral_field, grid3) ≈ profile atol = 1.0e-12
    @test DFMMonteCarlo._slit_3d_lateral_gradient_integral(field, grid3) ≈ 0.0 atol = 1.0e-12
    lateral_gradient_integral =
        DFMMonteCarlo._slit_3d_lateral_gradient_integral(lateral_field, grid3)
    @test lateral_gradient_integral > 0.0
    lateral_gradient_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1, lateral_gradient_coeff=3.0)
    @test slit_hevk1_energy_3d(field, grid3, lateral_gradient_params) ≈
          slit_hevk1_energy_3d(field, grid3, params) atol = 1.0e-10
    @test slit_hevk1_energy_3d(lateral_field, grid3, lateral_gradient_params) ≈
          slit_hevk1_energy_3d(lateral_field, grid3, params) +
          3.0 * lateral_gradient_integral atol = 1.0e-10
    projected_entropy_reference = DFMMonteCarlo.normalize_slit_density(
        DFMMonteCarlo.slit_reference_profile(grid3.z; B=25.0, C=0.5, L=1.0),
        grid3.slit_grid)
    projected_entropy_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1, projected_entropy_coeff=-1.5,
        projected_entropy_reference=projected_entropy_reference)
    @test DFMMonteCarlo._slit_projected_entropy_integral(
        projected_entropy_reference, grid3.slit_grid,
        projected_entropy_reference) ≈ 0.0 atol = 1.0e-12
    constrained_mode = cos.(2.0 .* pi .* grid3.z ./ grid3.Lz)
    constrained_mode .-= DFMMonteCarlo.slit_integral(constrained_mode, grid3.slit_grid) / grid3.Lz
    eps_projected = 1.0e-5
    projected_plus = projected_entropy_reference .+ eps_projected .* constrained_mode
    projected_minus = projected_entropy_reference .- eps_projected .* constrained_mode
    @test DFMMonteCarlo.slit_integral(projected_plus, grid3.slit_grid) ≈ grid3.Lz atol = 1.0e-12
    @test DFMMonteCarlo.slit_integral(projected_minus, grid3.slit_grid) ≈ grid3.Lz atol = 1.0e-12
    projected_first_variation =
        (DFMMonteCarlo._slit_projected_entropy_integral(projected_plus,
             grid3.slit_grid, projected_entropy_reference) -
         DFMMonteCarlo._slit_projected_entropy_integral(projected_minus,
             grid3.slit_grid, projected_entropy_reference)) / (2.0 * eps_projected)
    @test abs(projected_first_variation) < 1.0e-8
    projected_quadratic_weight = DFMMonteCarlo._slit_projected_quadratic_weight(:wall,
        grid3.slit_grid; wall_width=0.2)
    @test length(projected_quadratic_weight) == grid3.Nz
    @test minimum(projected_quadratic_weight) > 0.0
    @test DFMMonteCarlo.slit_integral(projected_quadratic_weight, grid3.slit_grid) / grid3.Lz ≈
          1.0 atol = 1.0e-12
    projected_quadratic_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1, projected_quadratic_coeff=-0.25,
        projected_quadratic_reference=projected_entropy_reference,
        projected_quadratic_weight=projected_quadratic_weight)
    projected_center_weight = DFMMonteCarlo._slit_projected_quadratic_weight(:center,
        grid3.slit_grid; wall_width=0.2)
    projected_multi_quadratic_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1,
        projected_quadratic_coeffs=[-0.25, 0.15],
        projected_quadratic_references=[
            projected_entropy_reference,
            projected_entropy_reference,
        ],
        projected_quadratic_weights=[
            projected_quadratic_weight,
            projected_center_weight,
        ])
    @test DFMMonteCarlo._slit_projected_quadratic_integral(
        projected_entropy_reference, grid3.slit_grid, projected_entropy_reference,
        projected_quadratic_weight) ≈ 0.0 atol = 1.0e-12
    projected_quadratic_first_variation =
        (DFMMonteCarlo._slit_projected_quadratic_integral(projected_plus,
             grid3.slit_grid, projected_entropy_reference, projected_quadratic_weight) -
         DFMMonteCarlo._slit_projected_quadratic_integral(projected_minus,
             grid3.slit_grid, projected_entropy_reference, projected_quadratic_weight)) /
        (2.0 * eps_projected)
    @test abs(projected_quadratic_first_variation) < 1.0e-8
    projected_field = DFMMonteCarlo.slit_3d_field_from_profile(projected_entropy_reference, grid3)
    projected_plus_field = DFMMonteCarlo.slit_3d_field_from_profile(projected_plus, grid3)
    @test slit_hevk1_energy(projected_entropy_reference, grid3.slit_grid,
        projected_entropy_params) ≈ slit_hevk1_energy(projected_entropy_reference,
        grid3.slit_grid, params) atol = 1.0e-10
    @test slit_hevk1_energy(projected_plus, grid3.slit_grid, projected_entropy_params) -
          slit_hevk1_energy(projected_plus, grid3.slit_grid, params) ≈
          slit_hevk1_energy_3d(projected_plus_field, grid3, projected_entropy_params) -
          slit_hevk1_energy_3d(projected_plus_field, grid3, params) atol = 1.0e-10
    @test slit_hevk1_energy_3d(projected_field, grid3, projected_entropy_params) ≈
          slit_hevk1_energy_3d(projected_field, grid3, params) atol = 1.0e-10
    @test slit_hevk1_energy(projected_entropy_reference, grid3.slit_grid,
        projected_quadratic_params) ≈ slit_hevk1_energy(projected_entropy_reference,
        grid3.slit_grid, params) atol = 1.0e-10
    @test slit_hevk1_energy(projected_plus, grid3.slit_grid, projected_quadratic_params) -
          slit_hevk1_energy(projected_plus, grid3.slit_grid, params) ≈
          slit_hevk1_energy_3d(projected_plus_field, grid3, projected_quadratic_params) -
          slit_hevk1_energy_3d(projected_plus_field, grid3, params) atol = 1.0e-10
    multi_quadratic_correction_1d =
        slit_hevk1_energy(projected_plus, grid3.slit_grid, projected_multi_quadratic_params) -
        slit_hevk1_energy(projected_plus, grid3.slit_grid, params)
    expected_multi_quadratic_correction =
        -0.25 * DFMMonteCarlo._slit_projected_quadratic_integral(
            projected_plus, grid3.slit_grid, projected_entropy_reference,
            projected_quadratic_weight) +
        0.15 * DFMMonteCarlo._slit_projected_quadratic_integral(
            projected_plus, grid3.slit_grid, projected_entropy_reference,
            projected_center_weight)
    @test multi_quadratic_correction_1d ≈ expected_multi_quadratic_correction atol = 1.0e-10
    @test slit_hevk1_energy_3d(projected_plus_field, grid3,
        projected_multi_quadratic_params) - slit_hevk1_energy_3d(projected_plus_field,
        grid3, params) ≈ expected_multi_quadratic_correction atol = 1.0e-10
    @test slit_hevk1_energy(projected_entropy_reference, grid3.slit_grid,
        projected_multi_quadratic_params) ≈ slit_hevk1_energy(projected_entropy_reference,
        grid3.slit_grid, params) atol = 1.0e-10
    projected_local_coeffs = (0.40, -0.15, 0.20)
    projected_local_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1, projected_local_potential_coeffs=projected_local_coeffs,
        projected_local_potential_reference=projected_entropy_reference)
    @test DFMMonteCarlo._slit_projected_local_potential_integral(
        projected_entropy_reference, grid3.slit_grid, projected_entropy_reference,
        projected_local_coeffs) ≈ 0.0 atol = 1.0e-12
    projected_local_first_variation =
        (DFMMonteCarlo._slit_projected_local_potential_integral(projected_plus,
             grid3.slit_grid, projected_entropy_reference, projected_local_coeffs) -
         DFMMonteCarlo._slit_projected_local_potential_integral(projected_minus,
             grid3.slit_grid, projected_entropy_reference, projected_local_coeffs)) /
        (2.0 * eps_projected)
    @test abs(projected_local_first_variation) < 1.0e-8
    projected_local_correction_1d =
        slit_hevk1_energy(projected_plus, grid3.slit_grid, projected_local_params) -
        slit_hevk1_energy(projected_plus, grid3.slit_grid, params)
    expected_projected_local_correction =
        DFMMonteCarlo._slit_projected_local_potential_integral(projected_plus,
            grid3.slit_grid, projected_entropy_reference, projected_local_coeffs)
    @test projected_local_correction_1d ≈ expected_projected_local_correction atol = 1.0e-10
    @test slit_hevk1_energy_3d(projected_plus_field, grid3,
        projected_local_params) - slit_hevk1_energy_3d(projected_plus_field,
        grid3, params) ≈ expected_projected_local_correction atol = 1.0e-10
    @test slit_hevk1_energy(projected_entropy_reference, grid3.slit_grid,
        projected_local_params) ≈ slit_hevk1_energy(projected_entropy_reference,
        grid3.slit_grid, params) atol = 1.0e-10

    profile_kernel_matrix = [2.0 0.25; 0.25 0.5]
    @test DFMMonteCarlo._slit_projected_profile_kernel_matrix(profile_kernel_matrix) ≈
          profile_kernel_matrix
    @test_throws ArgumentError DFMMonteCarlo._slit_projected_profile_kernel_matrix(
        [1.0 2.0; 0.0 1.0])
    profile_kernel_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1, projected_profile_kernel_matrix=profile_kernel_matrix,
        projected_profile_kernel_reference=projected_entropy_reference)
    @test DFMMonteCarlo._slit_projected_profile_kernel_integral(
        projected_entropy_reference, grid3.slit_grid, projected_entropy_reference,
        profile_kernel_matrix) ≈ 0.0 atol = 1.0e-12
    profile_kernel_first_variation =
        (DFMMonteCarlo._slit_projected_profile_kernel_integral(projected_plus,
             grid3.slit_grid, projected_entropy_reference, profile_kernel_matrix) -
         DFMMonteCarlo._slit_projected_profile_kernel_integral(projected_minus,
             grid3.slit_grid, projected_entropy_reference, profile_kernel_matrix)) /
        (2.0 * eps_projected)
    @test abs(profile_kernel_first_variation) < 1.0e-8
    profile_kernel_coeffs =
        DFMMonteCarlo._slit_projected_profile_kernel_coefficients(projected_plus,
            grid3.slit_grid, projected_entropy_reference, profile_kernel_matrix)
    expected_profile_kernel_correction =
        0.5 * sum(profile_kernel_coeffs .* (profile_kernel_matrix * profile_kernel_coeffs))
    @test DFMMonteCarlo._slit_projected_profile_kernel_integral(projected_plus,
        grid3.slit_grid, projected_entropy_reference, profile_kernel_matrix) ≈
          expected_profile_kernel_correction atol = 1.0e-12
    @test slit_hevk1_energy(projected_plus, grid3.slit_grid, profile_kernel_params) -
          slit_hevk1_energy(projected_plus, grid3.slit_grid, params) ≈
          expected_profile_kernel_correction atol = 1.0e-10
    @test slit_hevk1_energy_3d(projected_plus_field, grid3,
        profile_kernel_params) - slit_hevk1_energy_3d(projected_plus_field,
        grid3, params) ≈ expected_profile_kernel_correction atol = 1.0e-10
    lateral_projected_plus = copy(projected_plus_field)
    for i in 1:grid3.Nx, j in 1:grid3.Ny, k in 1:grid3.Nz
        lateral_projected_plus[i, j, k] *= 1.0 + 0.05 * (isodd(i + j) ? 1.0 : -1.0)
    end
    @test DFMMonteCarlo._slit_3d_xy_profile(lateral_projected_plus, grid3) ≈
          projected_plus atol = 1.0e-12
    @test DFMMonteCarlo._slit_3d_projected_profile_kernel_integral(
        lateral_projected_plus, grid3, projected_entropy_reference,
        profile_kernel_matrix) ≈ expected_profile_kernel_correction atol = 1.0e-10
    @test slit_hevk1_energy(projected_entropy_reference, grid3.slit_grid,
        profile_kernel_params) ≈ slit_hevk1_energy(projected_entropy_reference,
        grid3.slit_grid, params) atol = 1.0e-10

    profile_kernel_coeff_target = [1.0e-4, -5.0e-5]
    profile_kernel_from_coeffs =
        DFMMonteCarlo._slit_profile_from_projected_profile_kernel_coefficients(
            profile_kernel_coeff_target, grid3.slit_grid, projected_entropy_reference)
    recovered_profile_kernel_coeffs =
        DFMMonteCarlo._slit_projected_profile_kernel_coefficients(
            profile_kernel_from_coeffs, grid3.slit_grid,
            projected_entropy_reference, [1.0 0.0; 0.0 1.0])
    @test recovered_profile_kernel_coeffs ≈ profile_kernel_coeff_target atol = 1.0e-12
    @test DFMMonteCarlo.slit_integral(profile_kernel_from_coeffs, grid3.slit_grid) ≈
          grid3.Lz atol = 1.0e-12
    @test DFMMonteCarlo._regularized_lateral_logdet([1.0, 3.0];
        curvature_floor=1.0e-12, curvature_regularization=0.0) ≈
          0.5 * log(3.0) atol = 1.0e-12
    @test DFMMonteCarlo._regularized_lateral_logdet([1.0, 3.0];
        curvature_floor=1.0e-12, curvature_regularization=2.0) ≈
          0.5 * (log(3.0) + log(5.0)) atol = 1.0e-12
    @test isinf(DFMMonteCarlo._regularized_lateral_logdet([-3.0];
        curvature_floor=1.0e-12, curvature_regularization=2.0))
    @test_throws ArgumentError DFMMonteCarlo._regularized_lateral_logdet([1.0];
        curvature_floor=1.0e-12, curvature_regularization=-1.0)
    determinant_logdet =
        DFMMonteCarlo._slit_lateral_logdet_for_projected_profile_kernel_coefficients(
            zeros(2), grid3, params, projected_entropy_reference;
            lateral_mode_count=1, lateral_z_basis_count=0,
            hessian_step=1.0e-4, curvature_floor=1.0e-12)
    @test isfinite(determinant_logdet.logdet)
    @test determinant_logdet.lateral_coefficient_count > 0
    @test determinant_logdet.negative_eigenvalue_count == 0
    regularized_determinant_logdet =
        DFMMonteCarlo._slit_lateral_logdet_for_projected_profile_kernel_coefficients(
            zeros(2), grid3, params, projected_entropy_reference;
            lateral_mode_count=1, lateral_z_basis_count=0,
            hessian_step=1.0e-4, curvature_floor=1.0e-12,
            curvature_regularization=2.0)
    @test regularized_determinant_logdet.lateral_coefficient_count ==
          determinant_logdet.lateral_coefficient_count
    @test regularized_determinant_logdet.logdet > determinant_logdet.logdet
    @test regularized_determinant_logdet.curvature_regularization == 2.0
    determinant_kernel =
        DFMMonteCarlo._slit_projected_profile_kernel_lateral_logdet_matrix(
            grid3, params, projected_entropy_reference;
            mode_count=2, lateral_mode_count=1, lateral_z_basis_count=0,
            lateral_hessian_step=1.0e-4, profile_hessian_step=1.0e-4,
            curvature_floor=1.0e-12)
    @test size(determinant_kernel.matrix) == (2, 2)
    @test determinant_kernel.matrix ≈ transpose(determinant_kernel.matrix) atol = 1.0e-10
    @test all(isfinite, determinant_kernel.matrix)
    @test determinant_kernel.lateral_coefficient_count ==
          determinant_logdet.lateral_coefficient_count
    @test determinant_kernel.reference_negative_eigenvalue_count == 0
    determinant_kernel_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1, projected_profile_kernel_matrix=determinant_kernel.matrix,
        projected_profile_kernel_reference=projected_entropy_reference)
    @test slit_hevk1_energy(projected_entropy_reference, grid3.slit_grid,
        determinant_kernel_params) ≈ slit_hevk1_energy(projected_entropy_reference,
        grid3.slit_grid, params) atol = 1.0e-10
    profile_entropy_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1, profile_entropy_coeff=1.5)
    local_a_params = SlitParams(B=25.0, C=0.5, A=1.5, alpha=0.0, model=:vk1)
    @test DFMMonteCarlo._slit_3d_profile_entropy_integral(field, grid3) ≈
          DFMMonteCarlo._slit_3d_profile_entropy_integral(lateral_field, grid3) atol = 1.0e-12
    @test slit_hevk1_energy(profile, grid3.slit_grid, profile_entropy_params) ≈
          slit_hevk1_energy(profile, grid3.slit_grid, local_a_params) atol = 1.0e-10
    @test slit_hevk1_energy_3d(field, grid3, profile_entropy_params) ≈
          slit_hevk1_energy_3d(field, grid3, local_a_params) atol = 1.0e-10
    local_entropy_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1, local_potential_coeffs=(2.0, 0.0, 0.0))
    scalar_entropy_params = SlitParams(B=25.0, C=0.5, A=2.0, alpha=0.0,
        model=:vk1)
    @test slit_hevk1_energy(profile, grid3.slit_grid, local_entropy_params) ≈
          slit_hevk1_energy(profile, grid3.slit_grid, scalar_entropy_params) atol = 1.0e-10
    @test slit_hevk1_energy_3d(field, grid3, local_entropy_params) ≈
          slit_hevk1_energy_3d(field, grid3, scalar_entropy_params) atol = 1.0e-10
    blocked_energy_system = DFMMonteCarlo.SlitHEVK13DBlockedSystem(
        copy(lateral_field), grid3, lateral_gradient_params; profile_measure_power=1.25)
    @test hasproperty(blocked_energy_system, :energy_workspace)
    blocked_workspace = blocked_energy_system.energy_workspace
    blocked_workspace_k1 = blocked_workspace.k1
    direct_blocked_energy =
        DFMMonteCarlo.slit_hevk1_energy_3d_blocked_profile_measure(
            blocked_energy_system.phi, grid3, lateral_gradient_params;
            profile_measure_power=1.25)
    @test blocked_energy_system.energy ≈ direct_blocked_energy atol = 1.0e-10
    blocked_energy_action = DFMMonteCarlo.SlitLateralFractionPairAction(1, 1, 2, 1.0e-4)
    Arianna.perform_action!(blocked_energy_system, blocked_energy_action)
    @test blocked_energy_system.energy ≈
          DFMMonteCarlo.slit_hevk1_energy_3d_blocked_profile_measure(
              blocked_energy_system.phi, grid3, lateral_gradient_params;
              profile_measure_power=1.25) atol = 1.0e-10
    @test blocked_energy_system.energy_workspace === blocked_workspace
    @test blocked_energy_system.energy_workspace.k1 === blocked_workspace_k1
    local_fit = fit_slit_vk1_local_potential(Nz=16, B=25.0, C=0.5, A=0.0)
    @test all(isfinite, local_fit.coefficients)
    @test isfinite(local_fit.chemical_potential)
    @test local_fit.force_r2 > 0.75
    corrected_force = DFMMonteCarlo._slit_hevk1_functional_derivative(
        local_fit.target_phi, local_fit.grid, local_fit.params)
    base_force = DFMMonteCarlo._slit_hevk1_functional_derivative(
        local_fit.target_phi, local_fit.grid,
        SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0, model=:vk1))
    @test std(corrected_force) < std(base_force)

    outdir = mktempdir()
    result = run_slit_3d_vk1_mc(outdir; Nx=4, Nz=8, Lz=1.0, Lxy=1.0,
        B=25.0, C=0.5, A=0.0, steps=120, burnin_steps=30,
        sample_stride=10, initial_pair_amplitude=1.0e-3,
        initial_plane_amplitude=5.0e-4, initial_global_amplitude=2.0e-4,
        tune_rounds=1, tune_steps=12, seed=602)

    @test result.sample_count > 0
    @test result.steps == 120
    @test result.burnin_steps == 30
    @test 0.0 <= result.acceptance_rate <= 1.0
    @test DFMMonteCarlo.slit_integral(result.mean_phi, result.grid.slit_grid) ≈
          result.grid.Lz atol = 1.0e-10
    @test minimum(result.mean_phi) > 0.0
    @test isfinite(result.energy_mean)
    @test result.xi_eff > 0.0
    @test result.xi_fit_rmse ≈ sqrt(result.xi_fit_error)
    @test isfinite(result.xi_fit_r2)
    @test result.xi_fit_r2 <= 1.0
    @test result.xi_fit_max_abs_error >= result.xi_fit_rmse
    @test result.xi_fit_wall_rms >= 0.0
    @test result.xi_fit_quality_pass in (0.0, 1.0)
    @test isfinite(result.lateral_variance_mean)
    @test result.sampler == "arianna_metropolis"
    @test isfile(result.profile_csv)
    @test isfile(result.diagnostics_csv)
    @test isfile(result.profile_svg)
    @test isfile(result.markdown_path)
    @test isfile(result.metadata_path)
    profile_blocks_csv = joinpath(outdir, "slit_3d_mc_profile_blocks.csv")
    xi_blocks_csv = joinpath(outdir, "slit_3d_mc_xi_blocks.csv")
    @test isfile(profile_blocks_csv)
    @test isfile(xi_blocks_csv)

    profile_header, profile_rows = DFMMonteCarlo._read_csv_strings(result.profile_csv)
    @test "series" in profile_header
    @test "phi" in profile_header
    @test "stderr" in profile_header
    series_idx = findfirst(==("series"), profile_header)
    phi_idx = findfirst(==("phi"), profile_header)
    @test any(row[series_idx] == "mc_mean_3d" for row in profile_rows)
    @test any(row[series_idx] == "mc_mean_3d_raw" for row in profile_rows)
    @test any(row[series_idx] == "eq22_fit_to_3d_mc" for row in profile_rows)
    primary_profile = [parse(Float64, row[phi_idx]) for row in profile_rows if row[series_idx] == "mc_mean_3d"]
    raw_profile = [parse(Float64, row[phi_idx]) for row in profile_rows if row[series_idx] == "mc_mean_3d_raw"]
    @test length(primary_profile) == length(raw_profile) == result.grid.Nz
    @test primary_profile ≈ reverse(primary_profile)

    block_header, block_rows = DFMMonteCarlo._read_csv_strings(profile_blocks_csv)
    @test "chain" in block_header
    @test "block" in block_header
    @test "sample_start" in block_header
    @test "sample_stop" in block_header
    @test "z_over_L" in block_header
    @test "phi" in block_header
    @test length(block_rows) >= result.grid.Nz

    xi_block_header, xi_block_rows = DFMMonteCarlo._read_csv_strings(xi_blocks_csv)
    @test "chain" in xi_block_header
    @test "block" in xi_block_header
    @test "sample_start" in xi_block_header
    @test "sample_stop" in xi_block_header
    @test "sample_count" in xi_block_header
    @test "xi_eff" in xi_block_header
    xi_block_idx = findfirst(==("xi_eff"), xi_block_header)
    @test all(isfinite(parse(Float64, row[xi_block_idx])) for row in xi_block_rows)

    diagnostic_header, diagnostic_rows = DFMMonteCarlo._read_csv_strings(result.diagnostics_csv)
    @test "sampler" in diagnostic_header
    @test "chain_count" in diagnostic_header
    @test "profile_estimator" in diagnostic_header
    @test "profile_symmetrized" in diagnostic_header
    @test "lateral_variance_mean" in diagnostic_header
    @test "xy_profile_roughness_mean" in diagnostic_header
    @test "xi_fit_rmse" in diagnostic_header
    @test "xi_fit_r2" in diagnostic_header
    @test "xi_fit_max_abs_error" in diagnostic_header
    @test "xi_fit_wall_rms" in diagnostic_header
    @test "xi_fit_quality_pass" in diagnostic_header
    @test "profile_measure" in diagnostic_header
    @test "profile_measure_log_jacobian_mean" in diagnostic_header
    @test length(diagnostic_rows) == 1
    profile_measure_idx = findfirst(==("profile_measure"), diagnostic_header)
    @test diagnostic_rows[1][profile_measure_idx] == "flat_cell"

    metadata = TOML.parsefile(result.metadata_path)
    @test metadata["run"]["mode"] == "slit_3d_hevk1_metropolis"
    @test metadata["run"]["sampler"] == "arianna_metropolis"
    @test metadata["grid"]["Nx"] == 4
    @test metadata["grid"]["Nz"] == 8
    @test metadata["parameters"]["A"] == 0.0
    @test metadata["sampling"]["chain_count"] == 1
    @test metadata["sampling"]["profile_symmetrized"] == true
    @test metadata["sampling"]["profile_measure"] == "flat_cell"
    @test haskey(metadata["diagnostics"], "xi_fit_r2")
    @test haskey(metadata["diagnostics"], "xi_fit_quality_pass")
    @test metadata["outputs"]["profile_blocks_csv"] == profile_blocks_csv
    @test metadata["outputs"]["xi_blocks_csv"] == xi_blocks_csv

    text = read(result.markdown_path, String)
    @test occursin("3D real-density HE-VK1 Metropolis", text)
    @test occursin("Arianna", text)
    @test occursin("periodic x/y", text)
    @test occursin("fit R2", text)

    sweep_outdir = mktempdir()
    sweep_result = run_slit_3d_vk1_mc(sweep_outdir; Nx=4, Nz=8, Lz=1.0, Lxy=1.0,
        B=25.0, C=0.5, A=0.0, sweeps=1.0, burnin_sweeps=0.25,
        sample_stride_sweeps=0.25, initial_pair_amplitude=1.0e-3,
        initial_plane_amplitude=5.0e-4, initial_global_amplitude=2.0e-4,
        tune_rounds=1, tune_steps=12, seed=603)
    @test sweep_result.steps == 128
    @test sweep_result.burnin_steps == 32
    @test sweep_result.sample_stride == 32
    sweep_header, sweep_rows = DFMMonteCarlo._read_csv_strings(sweep_result.diagnostics_csv)
    @test "proposals_per_sweep" in sweep_header
    @test "total_sweeps" in sweep_header
    @test "production_sweeps" in sweep_header
    @test "burnin_sweeps" in sweep_header
    @test "sample_stride_sweeps" in sweep_header
    @test "profile_chain_rhat_max" in sweep_header
    @test "profile_split_rhat_max" in sweep_header
    @test "profile_chain_rmse_max" in sweep_header
    proposals_idx = findfirst(==("proposals_per_sweep"), sweep_header)
    total_sweeps_idx = findfirst(==("total_sweeps"), sweep_header)
    production_sweeps_idx = findfirst(==("production_sweeps"), sweep_header)
    burnin_sweeps_idx = findfirst(==("burnin_sweeps"), sweep_header)
    stride_sweeps_idx = findfirst(==("sample_stride_sweeps"), sweep_header)
    @test parse(Float64, sweep_rows[1][proposals_idx]) == 128.0
    @test parse(Float64, sweep_rows[1][total_sweeps_idx]) ≈ 1.0
    @test parse(Float64, sweep_rows[1][production_sweeps_idx]) ≈ 0.75
    @test parse(Float64, sweep_rows[1][burnin_sweeps_idx]) ≈ 0.25
    @test parse(Float64, sweep_rows[1][stride_sweeps_idx]) ≈ 0.25
    sweep_metadata = TOML.parsefile(sweep_result.metadata_path)
    @test sweep_metadata["sampling"]["proposals_per_sweep"] == 128
    @test sweep_metadata["sampling"]["total_sweeps"] ≈ 1.0
    @test sweep_metadata["sampling"]["production_sweeps"] ≈ 0.75
    @test sweep_metadata["sampling"]["burnin_sweeps"] ≈ 0.25

    corrected_outdir = mktempdir()
    corrected_result = run_slit_3d_vk1_mc(corrected_outdir; Nx=4, Nz=8, Lz=1.0, Lxy=1.0,
        B=25.0, C=0.5, A=0.0, steps=90, burnin_steps=20, sample_stride=10,
        initial_pair_amplitude=1.0e-3, initial_plane_amplitude=5.0e-4,
        initial_global_amplitude=2.0e-4, tune_rounds=1, tune_steps=12,
        profile_measure=:flat_profile, initial_profile_mode_amplitude=5.0e-4,
        profile_mode_move_weight=0.30, profile_mode_count=3, seed=606)
    @test corrected_result.sample_count > 0
    corrected_header, corrected_rows =
        DFMMonteCarlo._read_csv_strings(corrected_result.diagnostics_csv)
    corrected_measure_idx = findfirst(==("profile_measure"), corrected_header)
    corrected_logj_idx = findfirst(==("profile_measure_log_jacobian_mean"), corrected_header)
    corrected_profile_acceptance_idx = findfirst(==("profile_mode_acceptance_rate"), corrected_header)
    @test "profile_mode_count" in corrected_header
    @test corrected_rows[1][corrected_measure_idx] == "flat_profile"
    @test isfinite(parse(Float64, corrected_rows[1][corrected_logj_idx]))
    @test 0.0 <= parse(Float64, corrected_rows[1][corrected_profile_acceptance_idx]) <= 1.0
    corrected_metadata = TOML.parsefile(corrected_result.metadata_path)
    @test corrected_metadata["sampling"]["profile_measure"] == "flat_profile"
    @test corrected_metadata["sampling"]["profile_measure_power"] == 1.0
    @test corrected_metadata["sampling"]["profile_mode_move_weight"] == 0.30
    @test corrected_metadata["sampling"]["profile_mode_count"] == 3

    width_move_outdir = mktempdir()
    width_move_result = run_slit_3d_vk1_mc(width_move_outdir;
        Nx=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        steps=90, burnin_steps=20, sample_stride=10,
        initial_pair_amplitude=1.0e-3, initial_plane_amplitude=5.0e-4,
        initial_global_amplitude=2.0e-4, initial_profile_mode_amplitude=5.0e-4,
        initial_reference_width_amplitude=2.0e-4, reference_width_move_weight=0.20,
        reference_width_reflection_move_weight=0.05,
        reference_width_xi=0.24, tune_rounds=0, profile_measure=:flat_profile,
        profile_mode_move_weight=0.10, profile_mode_count=3, seed=609)
    @test width_move_result.sample_count > 0
    width_header, width_rows = DFMMonteCarlo._read_csv_strings(width_move_result.diagnostics_csv)
    @test "reference_width_acceptance_rate" in width_header
    @test "reference_width_reflection_acceptance_rate" in width_header
    @test "reference_width_amplitude" in width_header
    @test "reference_width_xi" in width_header
    width_acceptance_idx = findfirst(==("reference_width_acceptance_rate"), width_header)
    width_reflection_acceptance_idx = findfirst(==("reference_width_reflection_acceptance_rate"),
        width_header)
    width_amplitude_idx = findfirst(==("reference_width_amplitude"), width_header)
    width_xi_idx = findfirst(==("reference_width_xi"), width_header)
    @test 0.0 <= parse(Float64, width_rows[1][width_acceptance_idx]) <= 1.0
    @test 0.0 <= parse(Float64, width_rows[1][width_reflection_acceptance_idx]) <= 1.0
    @test parse(Float64, width_rows[1][width_amplitude_idx]) == 2.0e-4
    @test parse(Float64, width_rows[1][width_xi_idx]) == 0.24
    width_metadata = TOML.parsefile(width_move_result.metadata_path)
    @test width_metadata["sampling"]["reference_width_move_weight"] == 0.20
    @test width_metadata["sampling"]["reference_width_reflection_move_weight"] == 0.05
    @test width_metadata["sampling"]["reference_width_xi"] == 0.24

    projected_blocked_outdir = mktempdir()
    projected_blocked_result = run_slit_3d_vk1_blocked_mc(projected_blocked_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        projected_entropy_coeff=-0.75, projected_entropy_reference=:saddle,
        steps=64, burnin_steps=16, sample_stride=16, chains=1,
        lateral_fraction_amplitude=1.0e-3, profile_transfer_amplitude=5.0e-4,
        profile_mode_amplitude=5.0e-4, reference_width_amplitude=2.0e-4,
        profile_measure_power=1.0, initial_profile=:eq22, seed=610)
    @test projected_blocked_result.sample_count > 0
    @test projected_blocked_result.params.projected_entropy_coeff == -0.75
    @test projected_blocked_result.params.projected_entropy_reference ≈
          projected_blocked_result.saddle_phi
    bare_blocked_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0, model=:vk1)
    @test DFMMonteCarlo.slit_hevk1_energy(projected_blocked_result.saddle_phi,
        projected_blocked_result.grid.slit_grid, projected_blocked_result.params) ≈
          DFMMonteCarlo.slit_hevk1_energy(projected_blocked_result.saddle_phi,
        projected_blocked_result.grid.slit_grid, bare_blocked_params) atol = 1.0e-10
    projected_header, projected_rows =
        DFMMonteCarlo._read_csv_strings(projected_blocked_result.diagnostics_csv)
    @test "projected_entropy_coeff" in projected_header
    @test "projected_entropy_reference_source" in projected_header
    projected_coeff_idx = findfirst(==("projected_entropy_coeff"), projected_header)
    projected_ref_idx = findfirst(==("projected_entropy_reference_source"), projected_header)
    @test parse(Float64, projected_rows[1][projected_coeff_idx]) == -0.75
    @test projected_rows[1][projected_ref_idx] == "saddle"
    projected_metadata = TOML.parsefile(projected_blocked_result.metadata_path)
    @test projected_metadata["parameters"]["projected_entropy_coeff"] == -0.75
    @test projected_metadata["parameters"]["projected_entropy_reference_source"] == "saddle"
    @test projected_metadata["parameters"]["A"] == 0.0

    projected_quadratic_outdir = mktempdir()
    projected_quadratic_result = run_slit_3d_vk1_blocked_mc(projected_quadratic_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        projected_quadratic_coeff=-0.30, projected_quadratic_reference=:saddle,
        projected_quadratic_weight=:wall, projected_quadratic_wall_width=0.2,
        steps=64, burnin_steps=16, sample_stride=16, chains=1,
        lateral_fraction_amplitude=1.0e-3, profile_transfer_amplitude=5.0e-4,
        profile_mode_amplitude=5.0e-4, reference_width_amplitude=2.0e-4,
        profile_measure_power=1.0, initial_profile=:eq22, seed=611)
    @test projected_quadratic_result.sample_count > 0
    @test projected_quadratic_result.params.projected_quadratic_coeff == -0.30
    @test projected_quadratic_result.params.projected_quadratic_reference ≈
          projected_quadratic_result.saddle_phi
    @test slit_hevk1_energy(projected_quadratic_result.saddle_phi,
        projected_quadratic_result.grid.slit_grid, projected_quadratic_result.params) ≈
          slit_hevk1_energy(projected_quadratic_result.saddle_phi,
        projected_quadratic_result.grid.slit_grid, bare_blocked_params) atol = 1.0e-10
    projected_quadratic_header, projected_quadratic_rows =
        DFMMonteCarlo._read_csv_strings(projected_quadratic_result.diagnostics_csv)
    @test "projected_quadratic_coeff" in projected_quadratic_header
    @test "projected_quadratic_reference_source" in projected_quadratic_header
    @test "projected_quadratic_weight_source" in projected_quadratic_header
    @test "projected_quadratic_wall_width" in projected_quadratic_header
    projected_quadratic_coeff_idx =
        findfirst(==("projected_quadratic_coeff"), projected_quadratic_header)
    projected_quadratic_ref_idx =
        findfirst(==("projected_quadratic_reference_source"), projected_quadratic_header)
    projected_quadratic_weight_idx =
        findfirst(==("projected_quadratic_weight_source"), projected_quadratic_header)
    projected_quadratic_width_idx =
        findfirst(==("projected_quadratic_wall_width"), projected_quadratic_header)
    @test parse(Float64, projected_quadratic_rows[1][projected_quadratic_coeff_idx]) == -0.30
    @test projected_quadratic_rows[1][projected_quadratic_ref_idx] == "saddle"
    @test projected_quadratic_rows[1][projected_quadratic_weight_idx] == "wall"
    @test parse(Float64, projected_quadratic_rows[1][projected_quadratic_width_idx]) == 0.2
    projected_quadratic_metadata = TOML.parsefile(projected_quadratic_result.metadata_path)
    @test projected_quadratic_metadata["parameters"]["projected_quadratic_coeff"] == -0.30
    @test projected_quadratic_metadata["parameters"]["projected_quadratic_reference_source"] == "saddle"
    @test projected_quadratic_metadata["parameters"]["projected_quadratic_weight_source"] == "wall"
    @test projected_quadratic_metadata["parameters"]["projected_quadratic_wall_width"] == 0.2

    projected_multi_outdir = mktempdir()
    projected_multi_result = run_slit_3d_vk1_blocked_mc(projected_multi_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        projected_quadratic_coeffs=[-0.30, 0.20],
        projected_quadratic_reference=:saddle,
        projected_quadratic_weights=[:center, :wall],
        projected_quadratic_wall_width=0.2,
        steps=64, burnin_steps=16, sample_stride=16, chains=1,
        lateral_fraction_amplitude=1.0e-3, profile_transfer_amplitude=5.0e-4,
        profile_mode_amplitude=5.0e-4, reference_width_amplitude=2.0e-4,
        profile_measure_power=1.0, initial_profile=:eq22, seed=612)
    @test projected_multi_result.sample_count > 0
    @test projected_multi_result.params.projected_quadratic_coeffs == [-0.30, 0.20]
    @test length(projected_multi_result.params.projected_quadratic_weights) == 2
    @test slit_hevk1_energy(projected_multi_result.saddle_phi,
        projected_multi_result.grid.slit_grid, projected_multi_result.params) ≈
          slit_hevk1_energy(projected_multi_result.saddle_phi,
        projected_multi_result.grid.slit_grid, bare_blocked_params) atol = 1.0e-10
    projected_multi_header, projected_multi_rows =
        DFMMonteCarlo._read_csv_strings(projected_multi_result.diagnostics_csv)
    @test "projected_quadratic_coeffs" in projected_multi_header
    @test "projected_quadratic_weight_sources" in projected_multi_header
    projected_multi_coeffs_idx =
        findfirst(==("projected_quadratic_coeffs"), projected_multi_header)
    projected_multi_weights_idx =
        findfirst(==("projected_quadratic_weight_sources"), projected_multi_header)
    @test projected_multi_rows[1][projected_multi_coeffs_idx] == "-0.3;0.2"
    @test projected_multi_rows[1][projected_multi_weights_idx] == "center;wall"
    projected_multi_metadata = TOML.parsefile(projected_multi_result.metadata_path)
    @test projected_multi_metadata["parameters"]["projected_quadratic_coeffs"] == "-0.3;0.2"
    @test projected_multi_metadata["parameters"]["projected_quadratic_weight_sources"] == "center;wall"

    projected_local_outdir = mktempdir()
    projected_local_result = run_slit_3d_vk1_blocked_mc(projected_local_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        projected_local_potential_coeffs=projected_local_coeffs,
        projected_local_potential_reference=:saddle,
        steps=64, burnin_steps=16, sample_stride=16, chains=1,
        lateral_fraction_amplitude=1.0e-3, profile_transfer_amplitude=5.0e-4,
        profile_mode_amplitude=5.0e-4, reference_width_amplitude=2.0e-4,
        profile_measure_power=1.0, initial_profile=:eq22, seed=613)
    @test projected_local_result.sample_count > 0
    @test projected_local_result.params.projected_local_potential_coeffs ==
          projected_local_coeffs
    @test projected_local_result.params.projected_local_potential_reference ≈
          projected_local_result.saddle_phi
    @test slit_hevk1_energy(projected_local_result.saddle_phi,
        projected_local_result.grid.slit_grid, projected_local_result.params) ≈
          slit_hevk1_energy(projected_local_result.saddle_phi,
        projected_local_result.grid.slit_grid, bare_blocked_params) atol = 1.0e-10
    projected_local_header, projected_local_rows =
        DFMMonteCarlo._read_csv_strings(projected_local_result.diagnostics_csv)
    @test "projected_local_potential_a1" in projected_local_header
    @test "projected_local_potential_a2" in projected_local_header
    @test "projected_local_potential_a3" in projected_local_header
    @test "projected_local_potential_reference_source" in projected_local_header
    projected_local_a1_idx =
        findfirst(==("projected_local_potential_a1"), projected_local_header)
    projected_local_ref_idx =
        findfirst(==("projected_local_potential_reference_source"), projected_local_header)
    @test parse(Float64, projected_local_rows[1][projected_local_a1_idx]) ==
          projected_local_coeffs[1]
    @test projected_local_rows[1][projected_local_ref_idx] == "saddle"
    projected_local_metadata = TOML.parsefile(projected_local_result.metadata_path)
    @test projected_local_metadata["parameters"]["projected_local_potential_a1"] ==
          projected_local_coeffs[1]
    @test projected_local_metadata["parameters"]["projected_local_potential_reference_source"] ==
          "saddle"

    projected_profile_kernel_outdir = mktempdir()
    projected_profile_kernel_result = run_slit_3d_vk1_blocked_mc(
        projected_profile_kernel_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        projected_profile_kernel_matrix=profile_kernel_matrix,
        projected_profile_kernel_reference=:saddle,
        steps=64, burnin_steps=16, sample_stride=16, chains=1,
        lateral_fraction_amplitude=1.0e-3, profile_transfer_amplitude=5.0e-4,
        profile_mode_amplitude=5.0e-4, reference_width_amplitude=2.0e-4,
        profile_measure_power=1.0, initial_profile=:eq22, seed=614)
    @test projected_profile_kernel_result.sample_count > 0
    @test projected_profile_kernel_result.params.projected_profile_kernel_matrix ≈
          profile_kernel_matrix
    @test projected_profile_kernel_result.params.projected_profile_kernel_reference ≈
          projected_profile_kernel_result.saddle_phi
    @test slit_hevk1_energy(projected_profile_kernel_result.saddle_phi,
        projected_profile_kernel_result.grid.slit_grid,
        projected_profile_kernel_result.params) ≈
          slit_hevk1_energy(projected_profile_kernel_result.saddle_phi,
        projected_profile_kernel_result.grid.slit_grid, bare_blocked_params) atol = 1.0e-10
    projected_profile_kernel_header, projected_profile_kernel_rows =
        DFMMonteCarlo._read_csv_strings(projected_profile_kernel_result.diagnostics_csv)
    @test "projected_profile_kernel_matrix_count" in projected_profile_kernel_header
    @test "projected_profile_kernel_matrix_trace" in projected_profile_kernel_header
    @test "projected_profile_kernel_matrix_max_abs" in projected_profile_kernel_header
    @test "projected_profile_kernel_reference_source" in projected_profile_kernel_header
    profile_kernel_count_idx =
        findfirst(==("projected_profile_kernel_matrix_count"),
            projected_profile_kernel_header)
    profile_kernel_trace_idx =
        findfirst(==("projected_profile_kernel_matrix_trace"),
            projected_profile_kernel_header)
    profile_kernel_ref_idx =
        findfirst(==("projected_profile_kernel_reference_source"),
            projected_profile_kernel_header)
    @test parse(Int, projected_profile_kernel_rows[1][profile_kernel_count_idx]) == 4
    @test parse(Float64, projected_profile_kernel_rows[1][profile_kernel_trace_idx]) ≈
          2.5
    @test projected_profile_kernel_rows[1][profile_kernel_ref_idx] == "saddle"
    projected_profile_kernel_metadata =
        TOML.parsefile(projected_profile_kernel_result.metadata_path)
    @test projected_profile_kernel_metadata["parameters"]["projected_profile_kernel_matrix_count"] ==
          4
    @test projected_profile_kernel_metadata["parameters"]["projected_profile_kernel_matrix_trace"] ≈
          2.5
    @test projected_profile_kernel_metadata["parameters"]["projected_profile_kernel_reference_source"] ==
          "saddle"

    width_system = DFMMonteCarlo.SlitHEVK13DSystem(copy(field), grid3, params;
        profile_measure=:flat_profile)
    width_action = DFMMonteCarlo.SlitReferenceWidthAction(0.0, 0.24)
    width_parameters = [0.20, 0.24]
    lower_width, upper_width = DFMMonteCarlo._slit_reference_width_delta_bounds(
        width_system.phi, width_system.grid, width_system.params, 0.24,
        width_system.phi_min; amplitude=0.20)
    @test lower_width < 0.0 < upper_width
    @test upper_width - lower_width <= 0.40
    rng = DFMMonteCarlo.StableRNG(711)
    Arianna.sample_action!(width_action, DFMMonteCarlo.SlitReferenceWidthPolicy(),
        width_parameters, width_system, rng)
    @test lower_width <= width_action.δ <= upper_width
    forward_logq = Arianna.log_proposal_density(width_action,
        DFMMonteCarlo.SlitReferenceWidthPolicy(), width_parameters, width_system)
    @test isfinite(forward_logq)
    old_energy, new_energy = Arianna.perform_action!(width_system, width_action)
    @test isfinite(old_energy)
    @test isfinite(new_energy)
    @test minimum(width_system.phi) > width_system.phi_min
    @test DFMMonteCarlo.slit_integral(
        DFMMonteCarlo._slit_3d_xy_profile(width_system.phi, grid3),
        grid3.slit_grid) ≈ grid3.Lz atol = 1.0e-10
    Arianna.invert_action!(width_action, width_system)
    backward_logq = Arianna.log_proposal_density(width_action,
        DFMMonteCarlo.SlitReferenceWidthPolicy(), width_parameters, width_system)
    @test isfinite(backward_logq)
    Arianna.revert_action!(width_system, width_action)
    @test width_system.phi ≈ field
    @test width_system.energy ≈ old_energy

    reflection_shape = DFMMonteCarlo._slit_reference_width_shape(grid3, params, 0.24)
    reflection_field = copy(field)
    for k in 1:grid3.Nz
        @views reflection_field[:, :, k] .+= 0.005 * reflection_shape[k]
    end
    reflection_system = DFMMonteCarlo.SlitHEVK13DSystem(reflection_field, grid3, params;
        profile_measure=:flat_profile)
    reflection_action = DFMMonteCarlo.SlitReferenceWidthReflectionAction(0.0, 0.24)
    reflection_parameters = [0.24]
    reflection_target = DFMMonteCarlo.slit_reference_profile(grid3.z;
        B=params.B, C=params.C, L=grid3.Lz, xi=0.24)
    reflection_profile_before = DFMMonteCarlo._slit_3d_xy_profile(
        reflection_system.phi, grid3)
    reflection_offset_before = sum((reflection_profile_before .- reflection_target) .*
                                   reflection_shape)
    Arianna.sample_action!(reflection_action,
        DFMMonteCarlo.SlitReferenceWidthReflectionPolicy(), reflection_parameters,
        reflection_system, rng)
    @test isfinite(reflection_action.δ)
    @test abs(reflection_action.δ) > 0.0
    @test Arianna.log_proposal_density(reflection_action,
        DFMMonteCarlo.SlitReferenceWidthReflectionPolicy(), reflection_parameters,
        reflection_system) == 0.0
    reflection_old_energy, reflection_new_energy =
        Arianna.perform_action!(reflection_system, reflection_action)
    @test isfinite(reflection_old_energy)
    @test isfinite(reflection_new_energy)
    @test minimum(reflection_system.phi) > reflection_system.phi_min
    @test DFMMonteCarlo.slit_integral(
        DFMMonteCarlo._slit_3d_xy_profile(reflection_system.phi, grid3),
        grid3.slit_grid) ≈ grid3.Lz atol = 1.0e-10
    reflection_profile_after = DFMMonteCarlo._slit_3d_xy_profile(
        reflection_system.phi, grid3)
    reflection_offset_after = sum((reflection_profile_after .- reflection_target) .*
                                  reflection_shape)
    @test reflection_offset_after ≈ -reflection_offset_before atol = 1.0e-10
    Arianna.invert_action!(reflection_action, reflection_system)
    @test Arianna.log_proposal_density(reflection_action,
        DFMMonteCarlo.SlitReferenceWidthReflectionPolicy(), reflection_parameters,
        reflection_system) == 0.0
    Arianna.perform_action!(reflection_system, reflection_action)
    @test reflection_system.phi ≈ reflection_field atol = 1.0e-12
    @test reflection_system.energy ≈ reflection_old_energy atol = 1.0e-10

    blocked_reflection_field = copy(field)
    for k in 1:grid3.Nz
        @views blocked_reflection_field[:, :, k] .+= 0.005 * reflection_shape[k]
    end
    blocked_reflection_system = DFMMonteCarlo.SlitHEVK13DBlockedSystem(
        blocked_reflection_field, grid3, params; profile_measure_power=1.0)
    blocked_reflection_action =
        DFMMonteCarlo.SlitBlockedReferenceWidthReflectionAction(0.0, 0.24)
    blocked_profile_before = DFMMonteCarlo._slit_3d_xy_profile(
        blocked_reflection_system.phi, grid3)
    blocked_offset_before = sum((blocked_profile_before .- reflection_target) .*
                                reflection_shape)
    Arianna.sample_action!(blocked_reflection_action,
        DFMMonteCarlo.SlitBlockedReferenceWidthReflectionPolicy(), reflection_parameters,
        blocked_reflection_system, rng)
    @test isfinite(blocked_reflection_action.δ)
    @test abs(blocked_reflection_action.δ) > 0.0
    @test Arianna.log_proposal_density(blocked_reflection_action,
        DFMMonteCarlo.SlitBlockedReferenceWidthReflectionPolicy(), reflection_parameters,
        blocked_reflection_system) == 0.0
    blocked_old_energy, blocked_new_energy =
        Arianna.perform_action!(blocked_reflection_system, blocked_reflection_action)
    @test isfinite(blocked_old_energy)
    @test isfinite(blocked_new_energy)
    @test minimum(blocked_reflection_system.phi) > blocked_reflection_system.phi_min
    @test DFMMonteCarlo.slit_integral(
        DFMMonteCarlo._slit_3d_xy_profile(blocked_reflection_system.phi, grid3),
        grid3.slit_grid) ≈ grid3.Lz atol = 1.0e-10
    blocked_profile_after = DFMMonteCarlo._slit_3d_xy_profile(
        blocked_reflection_system.phi, grid3)
    blocked_offset_after = sum((blocked_profile_after .- reflection_target) .*
                               reflection_shape)
    @test blocked_offset_after ≈ -blocked_offset_before atol = 1.0e-10
    Arianna.invert_action!(blocked_reflection_action, blocked_reflection_system)
    Arianna.perform_action!(blocked_reflection_system, blocked_reflection_action)
    @test blocked_reflection_system.phi ≈ blocked_reflection_field atol = 1.0e-12
    @test blocked_reflection_system.energy ≈ blocked_old_energy atol = 1.0e-10

    parallel_chain_outdir = mktempdir()
    parallel_chain_result = run_slit_3d_vk1_mc(parallel_chain_outdir;
        Nx=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        steps=80, burnin_steps=20, sample_stride=20, chains=2,
        initial_pair_amplitude=1.0e-3, initial_plane_amplitude=5.0e-4,
        initial_global_amplitude=2.0e-4, profile_measure=:flat_profile,
        initial_reference_width_amplitude=2.0e-4, reference_width_move_weight=0.20,
        reference_width_xi=0.24, tune_rounds=0, parallel_chains=true, seed=610)
    @test parallel_chain_result.chain_count == 2
    parallel_header, parallel_rows =
        DFMMonteCarlo._read_csv_strings(parallel_chain_result.diagnostics_csv)
    @test "parallel_chains" in parallel_header
    parallel_idx = findfirst(==("parallel_chains"), parallel_header)
    @test parallel_rows[1][parallel_idx] == "true"
    parallel_metadata = TOML.parsefile(parallel_chain_result.metadata_path)
    @test parallel_metadata["sampling"]["parallel_chains"] == true

    drift_outdir = mktempdir()
    stable_blocks_csv = DFMMonteCarlo._write_dataframe(joinpath(drift_outdir,
        "stable_xi_blocks.csv"), [
        (chain=1, block=1, sample_start=1, sample_stop=2, sample_count=2,
            xi_eff=0.240, xi_fit_error=1.0e-5),
        (chain=1, block=2, sample_start=3, sample_stop=4, sample_count=2,
            xi_eff=0.242, xi_fit_error=1.0e-5),
        (chain=1, block=3, sample_start=5, sample_stop=6, sample_count=2,
            xi_eff=0.241, xi_fit_error=1.0e-5),
        (chain=2, block=1, sample_start=1, sample_stop=2, sample_count=2,
            xi_eff=0.239, xi_fit_error=1.0e-5),
        (chain=2, block=2, sample_start=3, sample_stop=4, sample_count=2,
            xi_eff=0.240, xi_fit_error=1.0e-5),
        (chain=2, block=3, sample_start=5, sample_stop=6, sample_count=2,
            xi_eff=0.238, xi_fit_error=1.0e-5),
    ])
    drifting_blocks_csv = DFMMonteCarlo._write_dataframe(joinpath(drift_outdir,
        "drifting_xi_blocks.csv"), [
        (chain=1, block=1, sample_start=1, sample_stop=2, sample_count=2,
            xi_eff=0.240, xi_fit_error=1.0e-5),
        (chain=1, block=2, sample_start=3, sample_stop=4, sample_count=2,
            xi_eff=0.280, xi_fit_error=1.0e-5),
        (chain=1, block=3, sample_start=5, sample_stop=6, sample_count=2,
            xi_eff=0.330, xi_fit_error=1.0e-5),
        (chain=2, block=1, sample_start=1, sample_stop=2, sample_count=2,
            xi_eff=0.235, xi_fit_error=1.0e-5),
        (chain=2, block=2, sample_start=3, sample_stop=4, sample_count=2,
            xi_eff=0.285, xi_fit_error=1.0e-5),
        (chain=2, block=3, sample_start=5, sample_stop=6, sample_count=2,
            xi_eff=0.360, xi_fit_error=1.0e-5),
    ])
    stable_drift = DFMMonteCarlo._slit_xi_block_drift_diagnostics(stable_blocks_csv;
        drift_tolerance=0.05)
    drifting_drift = DFMMonteCarlo._slit_xi_block_drift_diagnostics(drifting_blocks_csv;
        drift_tolerance=0.05)
    @test stable_drift.block_drift_pass == 1.0
    @test stable_drift.max_abs_start_end_delta < 0.05
    @test stable_drift.monotone_chain_fraction == 0.0
    @test drifting_drift.block_drift_pass == 0.0
    @test drifting_drift.max_abs_start_end_delta > 0.05
    @test drifting_drift.monotone_chain_fraction == 1.0

    measure_scan_outdir = mktempdir()
    measure_scan = run_slit_3d_measure_corrected_scan(measure_scan_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        target_xi_eff=0.24, profile_measure_powers=[1.0, 1.0002],
        profile_mode_amplitudes=[1.0e-4], profile_mode_count=3,
        steps=100, burnin_steps=20, sample_stride=20, chains=1,
        initial_pair_amplitude=1.0e-3, initial_plane_amplitude=5.0e-4,
        initial_global_amplitude=2.0e-4, profile_mode_move_weight=0.25,
        initial_reference_width_amplitude=2.0e-4, reference_width_move_weight=0.10,
        xi_tolerance=0.2, rhat_threshold=1.25, r2_threshold=0.5,
        parallel_chains=true, seed=607)
    @test measure_scan.run_count == 2
    @test measure_scan.target_xi_eff == 0.24
    @test measure_scan.best_profile_measure_power in (1.0, 1.0002)
    @test measure_scan.best_profile_mode_amplitude == 1.0e-4
    @test isfile(measure_scan.summary_csv)
    @test isfile(measure_scan.scan_svg)
    @test isfile(measure_scan.markdown_path)
    @test isfile(measure_scan.metadata_path)
    @test isfile(measure_scan.best_profile_csv)
    @test isfile(measure_scan.best_diagnostics_csv)
    measure_header, measure_rows = DFMMonteCarlo._read_csv_strings(measure_scan.summary_csv)
    @test length(measure_rows) == 2
    @test "profile_measure_power" in measure_header
    @test "profile_mode_amplitude" in measure_header
    @test "profile_mode_acceptance_rate" in measure_header
    @test "reference_width_acceptance_rate" in measure_header
    @test "reference_width_move_weight" in measure_header
    @test "xi_block_drift_pass" in measure_header
    @test "xi_block_max_abs_start_end_delta" in measure_header
    @test "xi_block_monotone_chain_fraction" in measure_header
    @test "calibration_pass" in measure_header
    @test "selection_rank" in measure_header
    measure_rank_idx = findfirst(==("selection_rank"), measure_header)
    @test count(row -> row[measure_rank_idx] == "1", measure_rows) == 1
    measure_metadata = TOML.parsefile(measure_scan.metadata_path)
    @test measure_metadata["run"]["mode"] == "slit_3d_measure_corrected_scan"
    @test measure_metadata["target"]["observable"] == "AK Fig. 3 effective xi"

    measure_script_outdir = mktempdir()
    old_args = copy(ARGS)
    try
        empty!(ARGS)
        append!(ARGS, [
            "--outdir=$(measure_script_outdir)",
            "--nx=4",
            "--ny=4",
            "--nz=8",
            "--steps=80",
            "--burnin-steps=20",
            "--sample-stride=20",
            "--chains=1",
            "--target-xi-eff=0.24",
            "--profile-measure-powers=1.0",
            "--profile-mode-amplitudes=0.0001",
            "--profile-mode-count=3",
            "--profile-mode-move-weight=0.25",
            "--initial-reference-width-amplitude=0.0002",
            "--reference-width-move-weight=0.1",
            "--parallel-chains=true",
            "--initial-pair-amplitude=0.001",
            "--initial-plane-amplitude=0.0005",
            "--initial-global-amplitude=0.0002",
            "--seed=608",
        ])
        redirect_stdout(devnull) do
            include(joinpath(@__DIR__, "..", "scripts", "run_slit_3d_measure_corrected_scan.jl"))
        end
    finally
        empty!(ARGS)
        append!(ARGS, old_args)
    end
    @test isfile(joinpath(measure_script_outdir, "measure_corrected_scan.csv"))
    measure_script_metadata = TOML.parsefile(joinpath(measure_script_outdir, "run_metadata.toml"))
    @test measure_script_metadata["run"]["mode"] == "slit_3d_measure_corrected_scan"

    profile_direct_outdir = mktempdir()
    profile_direct_result = run_slit_3d_profile_direct_vk1_mc(profile_direct_outdir;
        Nx=4, Ny=4, Nz=12, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        steps=900, burnin_steps=200, sample_stride=20, proposal_width=0.012,
        chains=2, seed=611)
    @test profile_direct_result.sample_count > 0
    @test profile_direct_result.chain_count == 2
    @test profile_direct_result.sampler == "profile_direct_metropolis"
    @test profile_direct_result.lateral_variance_mean == 0.0
    @test DFMMonteCarlo.slit_integral(profile_direct_result.mean_phi,
        profile_direct_result.grid.slit_grid) ≈ profile_direct_result.grid.Lz atol = 1.0e-10
    @test minimum(profile_direct_result.mean_phi) > 0.0
    @test isfile(profile_direct_result.profile_csv)
    @test isfile(profile_direct_result.diagnostics_csv)
    @test isfile(profile_direct_result.profile_svg)
    @test isfile(profile_direct_result.markdown_path)
    @test isfile(profile_direct_result.metadata_path)
    profile_direct_header, profile_direct_rows =
        DFMMonteCarlo._read_csv_strings(profile_direct_result.diagnostics_csv)
    @test "sampler" in profile_direct_header
    @test "proposal_width" in profile_direct_header
    @test "profile_chain_rhat_max" in profile_direct_header
    @test "profile_split_rhat_max" in profile_direct_header
    pd_sampler_idx = findfirst(==("sampler"), profile_direct_header)
    @test profile_direct_rows[1][pd_sampler_idx] == "profile_direct_metropolis"
    profile_direct_metadata = TOML.parsefile(profile_direct_result.metadata_path)
    @test profile_direct_metadata["run"]["mode"] == "slit_3d_hevk1_profile_direct_metropolis"
    @test profile_direct_metadata["run"]["sampler"] == "profile_direct_metropolis"
    @test profile_direct_metadata["sampling"]["state_space"] == "positive_normalized_plane_profile"

    blocked_outdir = mktempdir()
    blocked_result = run_slit_3d_vk1_blocked_mc(blocked_outdir;
        Nx=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        local_potential_coeffs=local_fit.coefficients,
        lateral_gradient_coeff=0.75,
        steps=100, burnin_steps=20, sample_stride=10, chains=2,
        lateral_fraction_amplitude=1.0e-3, profile_transfer_amplitude=5.0e-4,
        profile_mode_amplitude=5.0e-4, reference_width_amplitude=2.0e-4,
        lateral_fraction_move_weight=0.50, profile_transfer_move_weight=0.25,
        profile_mode_move_weight=0.10, reference_width_move_weight=0.15,
        profile_measure_power=1.0, seed=612)
    @test blocked_result.sample_count > 0
    @test blocked_result.chain_count == 2
    @test blocked_result.sampler == "arianna_blocked_profile_fraction"
    @test DFMMonteCarlo.slit_integral(blocked_result.mean_phi,
        blocked_result.grid.slit_grid) ≈ blocked_result.grid.Lz atol = 1.0e-10
    @test minimum(blocked_result.mean_phi) > 0.0
    @test isfinite(blocked_result.energy_mean)
    @test isfinite(blocked_result.lateral_variance_mean)
    @test isfile(blocked_result.profile_csv)
    @test isfile(blocked_result.diagnostics_csv)
    @test isfile(blocked_result.profile_svg)
    @test isfile(blocked_result.markdown_path)
    @test isfile(blocked_result.metadata_path)
    blocked_header, blocked_rows = DFMMonteCarlo._read_csv_strings(blocked_result.diagnostics_csv)
    @test "state_space" in blocked_header
    @test "coordinate_log_jacobian_power" in blocked_header
    @test "lateral_fraction_acceptance_rate" in blocked_header
    @test "profile_transfer_acceptance_rate" in blocked_header
    @test "reference_width_acceptance_rate" in blocked_header
    @test "local_potential_a1" in blocked_header
    @test "lateral_gradient_coeff" in blocked_header
    blocked_sampler_idx = findfirst(==("sampler"), blocked_header)
    blocked_state_idx = findfirst(==("state_space"), blocked_header)
    blocked_power_idx = findfirst(==("coordinate_log_jacobian_power"), blocked_header)
    blocked_a1_idx = findfirst(==("local_potential_a1"), blocked_header)
    blocked_lateral_gradient_idx = findfirst(==("lateral_gradient_coeff"), blocked_header)
    @test blocked_rows[1][blocked_sampler_idx] == "arianna_blocked_profile_fraction"
    @test blocked_rows[1][blocked_state_idx] == "plane_mean_lateral_fraction"
    @test parse(Float64, blocked_rows[1][blocked_power_idx]) == 0.0
    @test parse(Float64, blocked_rows[1][blocked_a1_idx]) ≈ local_fit.coefficients[1]
    @test parse(Float64, blocked_rows[1][blocked_lateral_gradient_idx]) ≈ 0.75
    blocked_metadata = TOML.parsefile(blocked_result.metadata_path)
    @test blocked_metadata["run"]["mode"] == "slit_3d_hevk1_blocked_profile_fraction_metropolis"
    @test blocked_metadata["sampling"]["state_space"] == "plane_mean_lateral_fraction"
    @test blocked_metadata["sampling"]["coordinate_log_jacobian_power"] == 0.0
    @test blocked_metadata["parameters"]["local_potential_a1"] ≈ local_fit.coefficients[1]
    @test blocked_metadata["parameters"]["lateral_gradient_coeff"] ≈ 0.75

    blocked_replica_outdir = mktempdir()
    blocked_replica_result = run_slit_3d_vk1_blocked_mc(blocked_replica_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        local_potential_coeffs=local_fit.coefficients,
        steps=160, burnin_steps=40, sample_stride=20, chains=1,
        lateral_fraction_amplitude=1.0e-3, profile_transfer_amplitude=5.0e-4,
        profile_mode_amplitude=5.0e-4, reference_width_amplitude=2.0e-4,
        lateral_fraction_move_weight=0.50, profile_transfer_move_weight=0.25,
        profile_mode_move_weight=0.10, reference_width_move_weight=0.15,
        profile_measure_power=1.0, seed=613,
        replica_temperatures=[1.0, 1.5],
        replica_swap_interval=8)
    @test blocked_replica_result.sample_count > 0
    @test DFMMonteCarlo.slit_integral(blocked_replica_result.mean_phi,
        blocked_replica_result.grid.slit_grid) ≈ blocked_replica_result.grid.Lz atol = 1.0e-10
    blocked_replica_header, blocked_replica_rows =
        DFMMonteCarlo._read_csv_strings(blocked_replica_result.diagnostics_csv)
    @test "replica_count" in blocked_replica_header
    @test "replica_temperatures" in blocked_replica_header
    @test "replica_swap_attempts" in blocked_replica_header
    @test "replica_swap_acceptance_rate" in blocked_replica_header
    @test "replica_round_trips" in blocked_replica_header
    @test parse(Int, blocked_replica_rows[1][findfirst(==("replica_count"),
        blocked_replica_header)]) == 2
    @test parse(Int, blocked_replica_rows[1][findfirst(==("replica_swap_attempts"),
        blocked_replica_header)]) > 0
    blocked_replica_metadata = TOML.parsefile(blocked_replica_result.metadata_path)
    @test blocked_replica_metadata["run"]["sampler"] ==
          "arianna_blocked_profile_fraction_replica_exchange"
    @test blocked_replica_metadata["sampling"]["replica_count"] == 2

    profile_mode_outdir = mktempdir()
    profile_mode_result = run_slit_3d_profile_mode_vk1_mc(profile_mode_outdir;
        Nx=4, Ny=4, Nz=16, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        z_basis_count=4, steps=4_000, burnin_steps=1_000, sample_stride=40,
        proposal_width=0.035, chains=2, seed=604)
    @test profile_mode_result.sample_count > 0
    @test profile_mode_result.chain_count == 2
    @test profile_mode_result.sampler == "profile_mode_metropolis"
    @test profile_mode_result.xi_eff > profile_mode_result.xi_bulk
    @test profile_mode_result.xi_fit_quality_pass == 1.0
    @test DFMMonteCarlo.slit_integral(profile_mode_result.mean_phi,
        profile_mode_result.grid.slit_grid) ≈ profile_mode_result.grid.Lz atol = 1.0e-10
    @test profile_mode_result.lateral_variance_mean == 0.0
    @test isfile(profile_mode_result.profile_csv)
    @test isfile(profile_mode_result.diagnostics_csv)
    @test isfile(profile_mode_result.profile_svg)
    @test isfile(profile_mode_result.markdown_path)
    @test isfile(profile_mode_result.metadata_path)
    profile_mode_blocks_csv = joinpath(profile_mode_outdir, "slit_3d_mc_profile_blocks.csv")
    profile_mode_xi_blocks_csv = joinpath(profile_mode_outdir, "slit_3d_mc_xi_blocks.csv")
    @test isfile(profile_mode_blocks_csv)
    @test isfile(profile_mode_xi_blocks_csv)
    profile_mode_header, profile_mode_rows =
        DFMMonteCarlo._read_csv_strings(profile_mode_result.diagnostics_csv)
    @test "sampler" in profile_mode_header
    @test "z_basis_count" in profile_mode_header
    @test "profile_chain_rhat_max" in profile_mode_header
    pm_sampler_idx = findfirst(==("sampler"), profile_mode_header)
    @test profile_mode_rows[1][pm_sampler_idx] == "profile_mode_metropolis"

    spectral_basis = DFMMonteCarlo._slit_3d_spectral_basis(grid3;
        z_basis_count=2, lateral_mode_count=1, lateral_z_basis_count=1)
    @test size(spectral_basis, 1) == length(field)
    @test size(spectral_basis, 2) > 2
    spectral_zero = DFMMonteCarlo._slit_3d_field_from_spectral_coefficients(
        zeros(size(spectral_basis, 2)), profile, spectral_basis, grid3)
    @test spectral_zero ≈ field atol = 1.0e-12
    lateral_coeffs = zeros(size(spectral_basis, 2))
    lateral_column = findfirst(3:size(spectral_basis, 2)) do column
        trial_coeffs = zeros(size(spectral_basis, 2))
        trial_coeffs[column] = 0.15
        trial_field = DFMMonteCarlo._slit_3d_field_from_spectral_coefficients(
            trial_coeffs, profile, spectral_basis, grid3)
        DFMMonteCarlo._slit_3d_lateral_variance(trial_field, grid3) > 1.0e-8
    end
    @test lateral_column !== nothing
    lateral_coeffs[lateral_column] = 0.15
    spectral_lateral = DFMMonteCarlo._slit_3d_field_from_spectral_coefficients(
        lateral_coeffs, profile, spectral_basis, grid3)
    @test minimum(spectral_lateral) > 0.0
    @test DFMMonteCarlo.slit_integral(
        DFMMonteCarlo._slit_3d_xy_profile(spectral_lateral, grid3),
        grid3.slit_grid) ≈ grid3.Lz atol = 1.0e-12
    @test DFMMonteCarlo._slit_3d_lateral_variance(spectral_lateral, grid3) > 0.0
    z_basis = DFMMonteCarlo._slit_cosine_basis_matrix(grid3.slit_grid, 2)
    lateral_basis = DFMMonteCarlo._slit_3d_lateral_spectral_basis(grid3;
        lateral_mode_count=1, lateral_z_basis_count=1)
    blocked_profile_coeffs = [0.10, -0.05]
    blocked_lateral_coeffs = zeros(size(lateral_basis, 2))
    blocked_lateral_coeffs[1] = 0.20
    blocked_spectral = DFMMonteCarlo._slit_3d_field_from_blocked_spectral_coefficients(
        blocked_profile_coeffs, blocked_lateral_coeffs, profile, z_basis, lateral_basis, grid3)
    blocked_profile = DFMMonteCarlo._slit_profile_from_smooth_coefficients(
        blocked_profile_coeffs, profile, z_basis, grid3.slit_grid)
    @test DFMMonteCarlo._slit_3d_xy_profile(blocked_spectral, grid3) ≈
          blocked_profile atol = 1.0e-12
    @test DFMMonteCarlo._slit_3d_lateral_variance(blocked_spectral, grid3) > 0.0
    @test DFMMonteCarlo._slit_3d_lateral_variance_integral(field, grid3) ≈ 0.0 atol = 1.0e-12
    blocked_lateral_variance_integral =
        DFMMonteCarlo._slit_3d_lateral_variance_integral(blocked_spectral, grid3)
    @test blocked_lateral_variance_integral > 0.0
    lateral_variance_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1, lateral_variance_coeff=2.5)
    @test slit_hevk1_energy_3d(field, grid3, lateral_variance_params) ≈
          slit_hevk1_energy_3d(field, grid3, params) atol = 1.0e-10
    @test slit_hevk1_energy_3d(blocked_spectral, grid3, lateral_variance_params) ≈
          slit_hevk1_energy_3d(blocked_spectral, grid3, params) +
          2.5 * blocked_lateral_variance_integral atol = 1.0e-10
    zcoupled_lateral_coeffs = zeros(size(lateral_basis, 2))
    zcoupled_lateral_coeffs[3] = 0.20
    zcoupled_spectral = DFMMonteCarlo._slit_3d_field_from_blocked_spectral_coefficients(
        blocked_profile_coeffs, zcoupled_lateral_coeffs, profile, z_basis, lateral_basis, grid3)
    @test DFMMonteCarlo._slit_3d_lateral_fraction_z_gradient_integral(field, grid3) ≈
          0.0 atol = 1.0e-12
    @test DFMMonteCarlo._slit_3d_lateral_fraction_z_gradient_integral(blocked_spectral, grid3) ≈
          0.0 atol = 1.0e-12
    zcoupled_lateral_z_integral =
        DFMMonteCarlo._slit_3d_lateral_fraction_z_gradient_integral(zcoupled_spectral, grid3)
    @test zcoupled_lateral_z_integral > 0.0
    lateral_z_params = SlitParams(B=25.0, C=0.5, A=0.0, alpha=0.0,
        model=:vk1, lateral_fraction_z_gradient_coeff=1.75)
    @test slit_hevk1_energy_3d(field, grid3, lateral_z_params) ≈
          slit_hevk1_energy_3d(field, grid3, params) atol = 1.0e-10
    @test slit_hevk1_energy_3d(blocked_spectral, grid3, lateral_z_params) ≈
          slit_hevk1_energy_3d(blocked_spectral, grid3, params) atol = 1.0e-10
    @test slit_hevk1_energy_3d(zcoupled_spectral, grid3, lateral_z_params) ≈
          slit_hevk1_energy_3d(zcoupled_spectral, grid3, params) +
          1.75 * zcoupled_lateral_z_integral atol = 1.0e-10
    lateral_metadata = DFMMonteCarlo._slit_3d_lateral_spectral_mode_metadata(grid3;
        lateral_mode_count=1, lateral_z_basis_count=1)
    @test length(lateral_metadata) == size(lateral_basis, 2)
    @test lateral_metadata[1] == (mx=0, my=1, zmode=0, phase_kind=:cos)
    @test lateral_metadata[2] == (mx=0, my=1, zmode=0, phase_kind=:sin)
    @test lateral_metadata[3] == (mx=0, my=1, zmode=1, phase_kind=:cos)
    @test count(mode -> mode.zmode == 0, lateral_metadata) == 6
    @test count(mode -> mode.zmode == 1, lateral_metadata) == 6
    zcoupled_weights = [mode.zmode == 0 ? 0.0 : 2.5 for mode in lateral_metadata]
    @test DFMMonteCarlo._slit_lateral_mode_quadratic_energy(
        zeros(length(zcoupled_weights)), zcoupled_weights) ≈ 0.0
    test_lateral_coeffs = zeros(length(zcoupled_weights))
    test_lateral_coeffs[3] = 0.2
    @test DFMMonteCarlo._slit_lateral_mode_quadratic_energy(
        test_lateral_coeffs, zcoupled_weights) ≈ 0.5 * 2.5 * 0.2^2
    @test_throws DimensionMismatch DFMMonteCarlo._slit_lateral_mode_quadratic_energy(
        test_lateral_coeffs[1:(end - 1)], zcoupled_weights)
    full_spectral_weights = vcat([4.0, 1.5], zcoupled_weights)
    full_spectral_coeffs = zeros(length(full_spectral_weights))
    full_spectral_coeffs[1] = 0.25
    full_spectral_coeffs[5] = 0.2
    @test DFMMonteCarlo._slit_spectral_mode_quadratic_energy(
        full_spectral_coeffs, full_spectral_weights) ≈
          0.5 * (4.0 * 0.25^2 + 2.5 * 0.2^2)
    @test_throws DimensionMismatch DFMMonteCarlo._slit_spectral_mode_quadratic_energy(
        full_spectral_coeffs[1:(end - 1)], full_spectral_weights)
    quadratic_hessian = DFMMonteCarlo._finite_difference_hessian(
        x -> 0.5 * (3.0 * x[1]^2 + 5.0 * x[2]^2), zeros(2); step=1.0e-4)
    @test quadratic_hessian ≈ [3.0 0.0; 0.0 5.0] atol = 1.0e-6
    nonlinear_scalar(x) = 4.0 + 1.5 * x[1] - 0.7 * x[2] +
                          0.5 * (2.0 * x[1]^2 + 0.6 * x[1] * x[2] +
                                 3.0 * x[2]^2) + 0.25 * x[1]^3
    scalar_gradient = DFMMonteCarlo._finite_difference_gradient_vector(
        nonlinear_scalar, zeros(2); step=1.0e-5)
    @test scalar_gradient ≈ [1.5, -0.7] atol = 1.0e-8
    projected_scalar_zero = DFMMonteCarlo._projected_scalar_functional_value(
        nonlinear_scalar, zeros(2); finite_difference_step=1.0e-5)
    @test projected_scalar_zero ≈ 0.0 atol = 1.0e-10
    projected_scalar_mode = [1.0e-4, 0.0]
    projected_scalar_first_variation =
        (DFMMonteCarlo._projected_scalar_functional_value(
             nonlinear_scalar, projected_scalar_mode; finite_difference_step=1.0e-5) -
         DFMMonteCarlo._projected_scalar_functional_value(
             nonlinear_scalar, -projected_scalar_mode; finite_difference_step=1.0e-5)) / 2.0e-4
    @test abs(projected_scalar_first_variation) < 1.0e-6
    @test DFMMonteCarlo._projected_scalar_functional_value(
        nonlinear_scalar, [0.2, -0.1]; finite_difference_step=1.0e-5) ≈
          (nonlinear_scalar([0.2, -0.1]) - nonlinear_scalar(zeros(2)) -
           sum([1.5, -0.7] .* [0.2, -0.1])) atol = 1.0e-8
    @test DFMMonteCarlo._profile_coefficient_trust_region_barrier(
        [0.1, 0.2], 0.3; strength=10.0) ≈ 0.0
    @test DFMMonteCarlo._profile_coefficient_trust_region_barrier(
        [0.4, 0.0], 0.3; strength=10.0) ≈
          10.0 * (0.4 / 0.3 - 1.0)^2
    @test_throws ArgumentError DFMMonteCarlo._profile_coefficient_trust_region_barrier(
        [0.1], 0.0; strength=10.0)
    @test_throws ArgumentError DFMMonteCarlo._profile_coefficient_trust_region_barrier(
        [0.1], 0.3; strength=-1.0)
    profile_logj_zero = DFMMonteCarlo._slit_profile_coordinate_log_jacobian(
        zeros(2), profile, z_basis, grid3.slit_grid)
    @test isfinite(profile_logj_zero)
    profile_logj_shifted = DFMMonteCarlo._slit_profile_coordinate_log_jacobian(
        [0.15, -0.05], profile, z_basis, grid3.slit_grid)
    @test isfinite(profile_logj_shifted)
    projected_profile_logj_zero =
        DFMMonteCarlo._slit_projected_profile_coordinate_log_jacobian(
            zeros(2), profile, z_basis, grid3.slit_grid)
    @test projected_profile_logj_zero ≈ 0.0 atol = 1.0e-10
    profile_logj_mode = [1.0e-4, 0.0]
    projected_logj_first_variation =
        (DFMMonteCarlo._slit_projected_profile_coordinate_log_jacobian(
             profile_logj_mode, profile, z_basis, grid3.slit_grid) -
         DFMMonteCarlo._slit_projected_profile_coordinate_log_jacobian(
             -profile_logj_mode, profile, z_basis, grid3.slit_grid)) / (2.0e-4)
    @test abs(projected_logj_first_variation) < 1.0e-6
    profile_logj_energy = DFMMonteCarlo._slit_profile_coordinate_log_jacobian_energy(
        [0.15, -0.05], 0.75, profile, z_basis, grid3.slit_grid)
    @test isfinite(profile_logj_energy)
    @test DFMMonteCarlo._slit_profile_coordinate_log_jacobian_energy(
        zeros(2), 0.75, profile, z_basis, grid3.slit_grid) ≈ 0.0 atol = 1.0e-10
    profile_xi_values = [0.12, 0.20, 0.35, 0.45]
    profile_xi_potential_values = [0.20, 0.08, 0.0, 0.03]
    @test DFMMonteCarlo._slit_profile_xi_potential_coeff(1.25) == 1.25
    @test DFMMonteCarlo._slit_profile_xi_potential_value(
        0.275, profile_xi_values, profile_xi_potential_values) ≈ 0.04
    projected_xi_potential_zero =
        DFMMonteCarlo._slit_projected_profile_xi_potential(
            zeros(2), profile, z_basis, grid3.slit_grid,
            profile_xi_values, profile_xi_potential_values)
    @test projected_xi_potential_zero ≈ 0.0 atol = 1.0e-8
    profile_xi_mode = [1.0e-4, 0.0]
    projected_xi_first_variation =
        (DFMMonteCarlo._slit_projected_profile_xi_potential(
             profile_xi_mode, profile, z_basis, grid3.slit_grid,
             profile_xi_values, profile_xi_potential_values) -
         DFMMonteCarlo._slit_projected_profile_xi_potential(
             -profile_xi_mode, profile, z_basis, grid3.slit_grid,
             profile_xi_values, profile_xi_potential_values)) / (2.0e-4)
    @test abs(projected_xi_first_variation) < 1.0e-5
    profile_xi_potential_energy =
        DFMMonteCarlo._slit_profile_xi_potential_energy(
            [0.15, -0.05], 1.25, profile, z_basis, grid3.slit_grid,
            profile_xi_values, profile_xi_potential_values)
    @test isfinite(profile_xi_potential_energy)
    @test DFMMonteCarlo._slit_profile_xi_potential_energy(
        zeros(2), 1.25, profile, z_basis, grid3.slit_grid,
        profile_xi_values, profile_xi_potential_values) ≈ 0.0 atol = 1.0e-8

    profile_quadratic_matrix = [2.0 0.25; 0.25 0.5]
    @test DFMMonteCarlo._slit_profile_mode_quadratic_matrix(
        profile_quadratic_matrix, 2) ≈ profile_quadratic_matrix
    @test_throws DimensionMismatch DFMMonteCarlo._slit_profile_mode_quadratic_matrix(
        profile_quadratic_matrix, 3)
    @test_throws ArgumentError DFMMonteCarlo._slit_profile_mode_quadratic_matrix(
        [2.0 0.4; 0.25 0.5], 2)
    @test DFMMonteCarlo._slit_profile_mode_quadratic_energy(
        [0.2, -0.1], profile_quadratic_matrix) ≈ 0.0375

    replica_xi_values = [0.05, 0.20, 0.50, 3.0]
    replica_xi_potential_values = [0.5, 0.0, 0.4, 2.0]
    zero_lateral_basis = DFMMonteCarlo._slit_3d_lateral_spectral_basis(grid3;
        lateral_mode_count=0, lateral_z_basis_count=0)
    profile_quadratic_coeffs = [0.2, -0.1]
    profile_quadratic_field =
        DFMMonteCarlo._slit_3d_field_from_blocked_spectral_coefficients(
            profile_quadratic_coeffs, Float64[], profile, z_basis,
            zero_lateral_basis, grid3)
    profile_quadratic_base_energy =
        slit_hevk1_energy_3d(profile_quadratic_field, grid3, params)
    @test DFMMonteCarlo._slit_3d_blocked_spectral_energy_from_coefficients(
        profile_quadratic_coeffs, grid3, params, profile, z_basis,
        zero_lateral_basis, zeros(2);
        profile_mode_quadratic_matrix=profile_quadratic_matrix) ≈
          profile_quadratic_base_energy + 0.0375

    replica_trace =
        DFMMonteCarlo._run_slit_3d_blocked_spectral_replica_exchange_trace(
            grid3, params, profile, z_basis, zero_lateral_basis;
            steps=80, burnin_steps=10, sample_stride=10,
            proposal_width=0.8, profile_proposal_width=0.8,
            lateral_proposal_width=0.8, target_temperature=1.0e6,
            replica_temperatures=[1.0e6], replica_swap_interval=10,
            seed=2468, profile_xi_potential_coeff=1000.0,
            profile_xi_potential_xi_values=replica_xi_values,
            profile_xi_potential_values=replica_xi_potential_values)
    @test replica_trace.accepted_moves > 0
    @test replica_trace.sample_count == length(replica_trace.energy_trace)
    for sample_idx in 1:replica_trace.sample_count
        sample_coeffs = [
            replica_trace.coefficient_traces[mode][sample_idx]
            for mode in eachindex(replica_trace.coefficient_traces)
        ]
        expected_energy =
            DFMMonteCarlo._slit_3d_blocked_spectral_energy_from_coefficients(
                sample_coeffs, grid3, params, profile, z_basis, zero_lateral_basis,
                zeros(Float64, length(sample_coeffs));
                profile_xi_potential_coeff=1000.0,
                profile_xi_potential_xi_values=replica_xi_values,
                profile_xi_potential_values=replica_xi_potential_values)
        @test replica_trace.energy_trace[sample_idx] ≈ expected_energy atol = 1.0e-8
    end

    spectral_outdir = mktempdir()
    spectral_result = run_slit_3d_spectral_mode_vk1_mc(spectral_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        profile_entropy_coeff=0.5,
        lateral_variance_coeff=0.25,
        lateral_fraction_z_gradient_coeff=0.75,
        lateral_mode_quadratic_coeffs=zcoupled_weights,
        spectral_mode_quadratic_coeffs=full_spectral_weights,
        profile_mode_quadratic_matrix=profile_quadratic_matrix,
        profile_coordinate_jacobian_coeff=0.5,
        profile_xi_potential_coeff=0.4,
        profile_xi_potential_xi_values=profile_xi_values,
        profile_xi_potential_values=profile_xi_potential_values,
        z_basis_count=2, lateral_mode_count=1, lateral_z_basis_count=1,
        steps=260, burnin_steps=60, sample_stride=20,
        proposal_width=0.04, profile_proposal_width=0.08,
        lateral_proposal_width=0.03, chains=2, seed=614)
    @test spectral_result.sample_count > 0
    @test spectral_result.chain_count == 2
    @test spectral_result.sampler == "spectral_mode_metropolis"
    @test spectral_result.lateral_variance_mean > 0.0
    @test DFMMonteCarlo.slit_integral(spectral_result.mean_phi,
        spectral_result.grid.slit_grid) ≈ spectral_result.grid.Lz atol = 1.0e-10
    @test isfile(spectral_result.profile_csv)
    @test isfile(spectral_result.diagnostics_csv)
    @test isfile(spectral_result.profile_svg)
    @test isfile(spectral_result.markdown_path)
    @test isfile(spectral_result.metadata_path)
    coefficient_csv = joinpath(spectral_outdir, "slit_3d_mc_coefficients.csv")
    hessian_csv = joinpath(spectral_outdir, "slit_3d_spectral_hessian.csv")
    branch_blocks_csv = joinpath(spectral_outdir, "slit_3d_spectral_branch_blocks.csv")
    branch_summary_csv = joinpath(spectral_outdir, "slit_3d_spectral_branch_summary.csv")
    @test isfile(coefficient_csv)
    @test isfile(hessian_csv)
    @test isfile(branch_blocks_csv)
    @test isfile(branch_summary_csv)
    spectral_header, spectral_rows =
        DFMMonteCarlo._read_csv_strings(spectral_result.diagnostics_csv)
    @test "sampler" in spectral_header
    @test "spectral_basis_count" in spectral_header
    @test "lateral_mode_count" in spectral_header
    @test "profile_entropy_coeff" in spectral_header
    @test "lateral_variance_coeff" in spectral_header
    @test "lateral_fraction_z_gradient_coeff" in spectral_header
    @test "lateral_mode_quadratic_coeff_count" in spectral_header
    @test "lateral_mode_quadratic_coeff_sum" in spectral_header
    @test "lateral_mode_quadratic_coeff_max" in spectral_header
    @test "spectral_mode_quadratic_coeff_count" in spectral_header
    @test "spectral_mode_quadratic_coeff_sum" in spectral_header
    @test "spectral_mode_quadratic_coeff_max" in spectral_header
    @test "profile_mode_quadratic_coeff_sum" in spectral_header
    @test "lateral_total_mode_quadratic_coeff_sum" in spectral_header
    @test "profile_mode_quadratic_matrix_coeff_count" in spectral_header
    @test "profile_mode_quadratic_matrix_trace" in spectral_header
    @test "profile_mode_quadratic_matrix_max_abs" in spectral_header
    @test "profile_coordinate_jacobian_coeff" in spectral_header
    @test "profile_xi_potential_coeff" in spectral_header
    @test "profile_xi_potential_point_count" in spectral_header
    @test "initial_profile" in spectral_header
    @test "profile_proposal_width" in spectral_header
    @test "lateral_proposal_width" in spectral_header
    @test "plane_normalized_lateral" in spectral_header
    spectral_sampler_idx = findfirst(==("sampler"), spectral_header)
    spectral_profile_entropy_idx = findfirst(==("profile_entropy_coeff"), spectral_header)
    spectral_lateral_variance_idx = findfirst(==("lateral_variance_coeff"), spectral_header)
    spectral_lateral_z_idx = findfirst(==("lateral_fraction_z_gradient_coeff"), spectral_header)
    spectral_mode_quad_count_idx =
        findfirst(==("lateral_mode_quadratic_coeff_count"), spectral_header)
    spectral_mode_quad_sum_idx =
        findfirst(==("lateral_mode_quadratic_coeff_sum"), spectral_header)
    spectral_mode_quad_max_idx =
        findfirst(==("lateral_mode_quadratic_coeff_max"), spectral_header)
    spectral_full_quad_count_idx =
        findfirst(==("spectral_mode_quadratic_coeff_count"), spectral_header)
    spectral_full_quad_sum_idx =
        findfirst(==("spectral_mode_quadratic_coeff_sum"), spectral_header)
    spectral_profile_quad_sum_idx =
        findfirst(==("profile_mode_quadratic_coeff_sum"), spectral_header)
    spectral_lateral_total_quad_sum_idx =
        findfirst(==("lateral_total_mode_quadratic_coeff_sum"), spectral_header)
    spectral_profile_quad_matrix_count_idx =
        findfirst(==("profile_mode_quadratic_matrix_coeff_count"), spectral_header)
    spectral_profile_quad_matrix_trace_idx =
        findfirst(==("profile_mode_quadratic_matrix_trace"), spectral_header)
    spectral_profile_quad_matrix_max_idx =
        findfirst(==("profile_mode_quadratic_matrix_max_abs"), spectral_header)
    spectral_profile_logj_idx =
        findfirst(==("profile_coordinate_jacobian_coeff"), spectral_header)
    spectral_profile_xi_potential_idx =
        findfirst(==("profile_xi_potential_coeff"), spectral_header)
    spectral_profile_xi_potential_count_idx =
        findfirst(==("profile_xi_potential_point_count"), spectral_header)
    spectral_initial_profile_idx = findfirst(==("initial_profile"), spectral_header)
    spectral_profile_width_idx = findfirst(==("profile_proposal_width"), spectral_header)
    spectral_lateral_width_idx = findfirst(==("lateral_proposal_width"), spectral_header)
    spectral_plane_normalized_idx = findfirst(==("plane_normalized_lateral"), spectral_header)
    @test spectral_rows[1][spectral_sampler_idx] == "spectral_mode_metropolis"
    @test parse(Float64, spectral_rows[1][spectral_profile_entropy_idx]) ≈ 0.5
    @test parse(Float64, spectral_rows[1][spectral_lateral_variance_idx]) ≈ 0.25
    @test parse(Float64, spectral_rows[1][spectral_lateral_z_idx]) ≈ 0.75
    @test parse(Int, spectral_rows[1][spectral_mode_quad_count_idx]) == length(zcoupled_weights)
    @test parse(Float64, spectral_rows[1][spectral_mode_quad_sum_idx]) ≈ sum(zcoupled_weights)
    @test parse(Float64, spectral_rows[1][spectral_mode_quad_max_idx]) ≈ maximum(zcoupled_weights)
    @test parse(Int, spectral_rows[1][spectral_full_quad_count_idx]) ==
          length(full_spectral_weights)
    @test parse(Float64, spectral_rows[1][spectral_full_quad_sum_idx]) ≈
          sum(full_spectral_weights)
    @test parse(Float64, spectral_rows[1][spectral_profile_quad_sum_idx]) ≈ 5.5
    @test parse(Float64, spectral_rows[1][spectral_lateral_total_quad_sum_idx]) ≈
          2.0 * sum(zcoupled_weights)
    @test parse(Int, spectral_rows[1][spectral_profile_quad_matrix_count_idx]) == 4
    @test parse(Float64, spectral_rows[1][spectral_profile_quad_matrix_trace_idx]) ≈
          sum(profile_quadratic_matrix[i, i] for i in 1:size(profile_quadratic_matrix, 1))
    @test parse(Float64, spectral_rows[1][spectral_profile_quad_matrix_max_idx]) ≈
          maximum(abs.(profile_quadratic_matrix))
    @test parse(Float64, spectral_rows[1][spectral_profile_logj_idx]) ≈ 0.5
    @test parse(Float64, spectral_rows[1][spectral_profile_xi_potential_idx]) ≈ 0.4
    @test parse(Int, spectral_rows[1][spectral_profile_xi_potential_count_idx]) ==
          length(profile_xi_values)
    @test spectral_rows[1][spectral_initial_profile_idx] == "eq22"
    @test parse(Float64, spectral_rows[1][spectral_profile_width_idx]) ≈ 0.08
    @test parse(Float64, spectral_rows[1][spectral_lateral_width_idx]) ≈ 0.03
    @test spectral_rows[1][spectral_plane_normalized_idx] == "true"
    coefficient_header, coefficient_rows = DFMMonteCarlo._read_csv_strings(coefficient_csv)
    @test "chain" in coefficient_header
    @test "sample" in coefficient_header
    @test "mode_index" in coefficient_header
    @test "mode_family" in coefficient_header
    @test "coefficient" in coefficient_header
    spectral_basis_count = parse(Int, spectral_rows[1][findfirst(==("spectral_basis_count"), spectral_header)])
    @test length(coefficient_rows) == spectral_result.sample_count * spectral_basis_count
    @test any(row -> row[findfirst(==("mode_family"), coefficient_header)] == "profile_z",
        coefficient_rows)
    @test any(row -> row[findfirst(==("mode_family"), coefficient_header)] == "lateral",
        coefficient_rows)
    hessian_header, hessian_rows = DFMMonteCarlo._read_csv_strings(hessian_csv)
    @test "eigen_index" in hessian_header
    @test "eigenvalue" in hessian_header
    @test "dominant_mode_index" in hessian_header
    @test "dominant_mode_family" in hessian_header
    @test length(hessian_rows) == spectral_basis_count
    hessian_eigenvalues = [
        parse(Float64, row[findfirst(==("eigenvalue"), hessian_header)])
        for row in hessian_rows
    ]
    @test all(isfinite, hessian_eigenvalues)
    @test issorted(hessian_eigenvalues)
    branch_header, branch_rows = DFMMonteCarlo._read_csv_strings(branch_blocks_csv)
    @test "branch" in branch_header
    @test "chain" in branch_header
    @test "block" in branch_header
    @test "mode_index" in branch_header
    @test "coefficient_mean" in branch_header
    @test length(branch_rows) > 0
    @test length(branch_rows) % spectral_basis_count == 0
    branch_summary_header, branch_summary_rows =
        DFMMonteCarlo._read_csv_strings(branch_summary_csv)
    @test "branch" in branch_summary_header
    @test "block_count" in branch_summary_header
    @test "sample_count" in branch_summary_header
    @test "mean_xi_eff" in branch_summary_header
    @test length(branch_summary_rows) >= 1

    replica_outdir = mktempdir()
    replica_result = run_slit_3d_spectral_mode_vk1_mc(replica_outdir;
        Nx=4, Ny=4, Nz=8, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        profile_entropy_coeff=0.5,
        lateral_variance_coeff=0.25,
        lateral_mode_quadratic_coeffs=zcoupled_weights,
        z_basis_count=2, lateral_mode_count=1, lateral_z_basis_count=1,
        steps=260, burnin_steps=60, sample_stride=20,
        proposal_width=0.04, chains=1, seed=715,
        replica_temperatures=[1.0, 1.6],
        replica_swap_interval=5,
        write_hessian_diagnostics=false)
    @test replica_result.sample_count > 0
    @test isfile(joinpath(replica_outdir, "slit_3d_mc_coefficients.csv"))
    @test isfile(joinpath(replica_outdir, "slit_3d_spectral_branch_summary.csv"))
    replica_header, replica_rows =
        DFMMonteCarlo._read_csv_strings(replica_result.diagnostics_csv)
    @test "replica_count" in replica_header
    @test "replica_temperatures" in replica_header
    @test "replica_swap_attempts" in replica_header
    @test "replica_swap_acceptance_rate" in replica_header
    @test "replica_round_trips" in replica_header
    @test parse(Int, replica_rows[1][findfirst(==("replica_count"), replica_header)]) == 2
    @test parse(Int, replica_rows[1][findfirst(==("replica_swap_attempts"), replica_header)]) > 0

    scan_outdir = mktempdir()
    scan_result = run_slit_3d_profile_mode_ak_scan(scan_outdir;
        Nx=4, Ny=4, Nz=16, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        target_xi_eff=0.28, temperatures=[0.75, 1.25],
        z_basis_count=4, steps=1_800, burnin_steps=600, sample_stride=60,
        proposal_width=0.035, chains=1, seed=605)
    @test scan_result.run_count == 2
    @test scan_result.target_xi_eff == 0.28
    @test scan_result.best_temperature in (0.75, 1.25)
    @test isfile(scan_result.summary_csv)
    @test isfile(scan_result.scan_svg)
    @test isfile(scan_result.markdown_path)
    @test isfile(scan_result.metadata_path)
    @test isfile(scan_result.best_profile_csv)
    @test isfile(scan_result.best_diagnostics_csv)
    scan_header, scan_rows = DFMMonteCarlo._read_csv_strings(scan_result.summary_csv)
    @test length(scan_rows) == 2
    @test "temperature" in scan_header
    @test "target_xi_eff" in scan_header
    @test "abs_xi_error" in scan_header
    @test "target_profile_rmse" in scan_header
    @test "target_profile_r2" in scan_header
    @test "selection_rank" in scan_header
    error_idx = findfirst(==("abs_xi_error"), scan_header)
    rank_idx = findfirst(==("selection_rank"), scan_header)
    scan_errors = [parse(Float64, row[error_idx]) for row in scan_rows]
    @test scan_result.best_abs_xi_error ≈ minimum(scan_errors)
    @test count(row -> row[rank_idx] == "1", scan_rows) == 1
    scan_metadata = TOML.parsefile(scan_result.metadata_path)
    @test scan_metadata["run"]["mode"] == "slit_3d_profile_mode_ak_scan"
    @test scan_metadata["target"]["observable"] == "AK Fig. 3 effective xi"
    scan_text = read(scan_result.markdown_path, String)
    @test occursin("Alexander-Katz Fig. 3", scan_text)
    @test occursin("effective correlation length", scan_text)
    @test occursin("calibration", scan_text)

    basis_scan_outdir = mktempdir()
    basis_scan_result = run_slit_3d_profile_mode_basis_ak_scan(basis_scan_outdir;
        Nx=4, Ny=4, Nz=16, Lz=1.0, Lxy=1.0, B=25.0, C=0.5, A=0.0,
        target_xi_eff=0.28, temperatures=[0.75, 1.25], z_basis_counts=[2, 3],
        steps=900, burnin_steps=300, sample_stride=60, proposal_width=0.035,
        chains=1, seed=607)
    @test basis_scan_result.run_count == 4
    @test basis_scan_result.target_xi_eff == 0.28
    @test basis_scan_result.best_z_basis_count in (2, 3)
    @test basis_scan_result.best_temperature in (0.75, 1.25)
    @test isfile(basis_scan_result.summary_csv)
    @test isfile(basis_scan_result.markdown_path)
    @test isfile(basis_scan_result.metadata_path)
    @test isfile(basis_scan_result.best_profile_csv)
    basis_header, basis_rows = DFMMonteCarlo._read_csv_strings(basis_scan_result.summary_csv)
    @test length(basis_rows) == 4
    @test "z_basis_count" in basis_header
    @test "temperature" in basis_header
    @test "selection_rank" in basis_header
    basis_rank_idx = findfirst(==("selection_rank"), basis_header)
    @test count(row -> row[basis_rank_idx] == "1", basis_rows) == 1
    basis_metadata = TOML.parsefile(basis_scan_result.metadata_path)
    @test basis_metadata["run"]["mode"] == "slit_3d_profile_mode_basis_ak_scan"
    @test basis_metadata["scan"]["z_basis_counts"] == [2, 3]
    @test basis_metadata["target"]["observable"] == "AK Fig. 3 effective xi"

    basis_script_outdir = mktempdir()
    basis_old_args = copy(ARGS)
    try
        empty!(ARGS)
        append!(ARGS, [
            "--outdir=$(basis_script_outdir)",
            "--nx=4",
            "--ny=4",
            "--nz=16",
            "--Lxy=1.0",
            "--Lz=1.0",
            "--B=25.0",
            "--C=0.5",
            "--A=0.0",
            "--target-xi-eff=0.28",
            "--temperatures=0.75",
            "--z-basis-counts=2",
            "--steps=500",
            "--burnin-steps=150",
            "--sample-stride=50",
            "--proposal-width=0.035",
            "--chains=1",
            "--seed=608",
        ])
        redirect_stdout(devnull) do
            include(joinpath(@__DIR__, "..", "scripts", "run_slit_3d_profile_mode_basis_ak_scan.jl"))
        end
    finally
        empty!(ARGS)
        append!(ARGS, basis_old_args)
    end
    @test isfile(joinpath(basis_script_outdir, "profile_mode_basis_temperature_scan.csv"))
    basis_script_metadata = TOML.parsefile(joinpath(basis_script_outdir, "run_metadata.toml"))
    @test basis_script_metadata["run"]["mode"] == "slit_3d_profile_mode_basis_ak_scan"
    @test basis_script_metadata["scan"]["z_basis_counts"] == [2]

    a_scan_outdir = mktempdir()
    a_scan_result = run_slit_3d_profile_mode_a_ak_scan(a_scan_outdir;
        Nx=4, Ny=4, Nz=16, Lz=1.0, Lxy=1.0, B=25.0, C=0.5,
        target_xi_eff=0.28, A_values=[0.0, 0.5], z_basis_count=2,
        temperature=1.0, steps=700, burnin_steps=200, sample_stride=50,
        proposal_width=0.2, chains=1, seed=609)
    @test a_scan_result.run_count == 2
    @test a_scan_result.target_xi_eff == 0.28
    @test a_scan_result.best_A in (0.0, 0.5)
    @test a_scan_result.best_temperature == 1.0
    @test isfile(a_scan_result.summary_csv)
    @test isfile(a_scan_result.markdown_path)
    @test isfile(a_scan_result.metadata_path)
    @test isfile(a_scan_result.best_profile_csv)
    a_header, a_rows = DFMMonteCarlo._read_csv_strings(a_scan_result.summary_csv)
    @test length(a_rows) == 2
    @test "A" in a_header
    @test "temperature" in a_header
    @test "selection_rank" in a_header
    a_rank_idx = findfirst(==("selection_rank"), a_header)
    @test count(row -> row[a_rank_idx] == "1", a_rows) == 1
    a_metadata = TOML.parsefile(a_scan_result.metadata_path)
    @test a_metadata["run"]["mode"] == "slit_3d_profile_mode_a_ak_scan"
    @test a_metadata["scan"]["A_values"] == [0.0, 0.5]
    @test a_metadata["scan"]["temperature"] == 1.0

    a_script_outdir = mktempdir()
    a_old_args = copy(ARGS)
    try
        empty!(ARGS)
        append!(ARGS, [
            "--outdir=$(a_script_outdir)",
            "--nx=4",
            "--ny=4",
            "--nz=16",
            "--Lxy=1.0",
            "--Lz=1.0",
            "--B=25.0",
            "--C=0.5",
            "--target-xi-eff=0.28",
            "--A-values=0.0",
            "--z-basis-count=2",
            "--temperature=1.0",
            "--steps=400",
            "--burnin-steps=100",
            "--sample-stride=50",
            "--proposal-width=0.2",
            "--chains=1",
            "--seed=610",
        ])
        redirect_stdout(devnull) do
            include(joinpath(@__DIR__, "..", "scripts", "run_slit_3d_profile_mode_a_ak_scan.jl"))
        end
    finally
        empty!(ARGS)
        append!(ARGS, a_old_args)
    end
    @test isfile(joinpath(a_script_outdir, "profile_mode_a_scan.csv"))
    a_script_metadata = TOML.parsefile(joinpath(a_script_outdir, "run_metadata.toml"))
    @test a_script_metadata["run"]["mode"] == "slit_3d_profile_mode_a_ak_scan"
    @test a_script_metadata["scan"]["A_values"] == [0.0]
end

@testset "3D grid and HE-URP energy" begin
    @test isfile(joinpath(@__DIR__, "..", "scripts", "run_direct_stability_probe.jl"))

    grid = Grid3D(4, 6.4)
    @test grid.dimension == 3
    @test grid.volume ≈ 6.4^3
    @test size(wavenumber_squared(grid)) == (4, 4, 4)

    phi = ones(4, 4, 4)
    params = HEURPParams(A=300.0, B=1.0, C=20.0, alpha=0.1)
    @test params.alpha == 0.1
    @test isfinite(heurp_energy(phi, grid, params))
    @test heurp_energy(phi, grid, params) ≈ 0.5 * params.B * params.C^2 * grid.volume
    @test DFMMonteCarlo._heurp_energy_sparse(phi, grid, params) ≈ heurp_energy(phi, grid, params)
    @test isfinite(hevk1_energy(phi, grid, params))
    @test hevk1_energy(phi, grid, params) ≈ heurp_energy(phi, grid, params)
    vk1_uniform = DFMMonteCarlo._vk1_coefficient_field(phi, grid)
    @test size(vk1_uniform) == size(phi)
    @test all(isfinite, vk1_uniform)
    @test all(value -> 1 / 6 <= value <= 1 / 4, vk1_uniform)
    @test hevk1_pressure_observable(phi, grid, params) ≈ params.C + 0.5 * params.B * params.C^2

    bad_phi = copy(phi)
    bad_phi[1, 1, 1] = 0.0
    @test_throws DomainError heurp_energy(bad_phi, grid, params)
    @test_throws DomainError hevk1_energy(bad_phi, grid, params)

    function cosine_density(grid; amplitude=0.05)
        field = ones(grid.Nx, grid.Nx, grid.Nx)
        for i in 1:grid.Nx, j in 1:grid.Nx, k in 1:grid.Nx
            field[i, j, k] = 1.0 + amplitude * cos(2π * (i - 1) / grid.Nx)
        end
        return field
    end
    function spectral_excess(field, grid, params; energy=heurp_energy)
        cell = grid.volume / length(field)
        real_term = 0.5 * params.B * params.C^2 * cell * sum(abs2, field)
        return energy(field, grid, params) - real_term
    end

    coarse_grid = Grid3D(4, 6.4)
    fine_grid = Grid3D(8, 6.4)
    coarse_field = cosine_density(coarse_grid)
    fine_field = cosine_density(fine_grid)
    coarse_spectral = spectral_excess(coarse_field, coarse_grid, params)
    fine_spectral = spectral_excess(fine_field, fine_grid, params)
    @test isapprox(fine_spectral, coarse_spectral; rtol=1e-3)
    @test DFMMonteCarlo.low_mode_structure_factor(phi) == 0.0
    coarse_structure = DFMMonteCarlo.low_mode_structure_factor(coarse_field)
    fine_structure = DFMMonteCarlo.low_mode_structure_factor(fine_field)
    @test coarse_structure > 0.0
    @test isapprox(fine_structure, coarse_structure; rtol=1e-12)
    @test DFMMonteCarlo._heurp_energy_sparse(fine_field, fine_grid, params) ≈
          heurp_energy(fine_field, fine_grid, params)
    @test DFMMonteCarlo._heurp_energy_fft(fine_field, fine_grid, params) ≈
          heurp_energy(fine_field, fine_grid, params)
    @test DFMMonteCarlo._energy_backend(Grid3D(8, 6.4)) == :fft
    @test DFMMonteCarlo._energy_backend(Grid3D(64, 6.4)) == :fft

    unsmeared = HEURPParams(A=300.0, B=1.0, C=20.0, alpha=0.0)
    smeared = HEURPParams(A=300.0, B=1.0, C=20.0, alpha=0.1)
    nonuniform = frozen_pressure_fields(grid)[end]
    @test heurp_energy(nonuniform, grid, smeared) > heurp_energy(nonuniform, grid, unsmeared)
    @test isfinite(hevk1_energy(nonuniform, grid, smeared))
    @test isfinite(hevk1_pressure_observable(nonuniform, grid, smeared))
    @test DFMMonteCarlo._proposals_per_sweep(grid) == 64
    @test DFMMonteCarlo._steps_from_sweeps(grid, 0.25) == 16
    @test DFMMonteCarlo._steps_from_sweeps(grid, 1.0) == 64

    probe_outdir = mktempdir()
    probe = run_direct_stability_probe(probe_outdir; Nx=4, C=20.0,
        modes=[(1, 0, 0), (2, 0, 0)], amplitude_count=8, boundary_margin=1.0e-3)
    @test isfile(probe.csv_path)
    @test isfile(probe.summary_path)
    probe_header, probe_rows = DFMMonteCarlo._read_csv_numeric(probe.csv_path)
    @test "mode_ax" in probe_header
    @test "amplitude" in probe_header
    @test "min_phi" in probe_header
    @test "energy" in probe_header
    @test "energy_delta" in probe_header
    @test "pressure_observable" in probe_header
    @test length(probe_rows) == 16
    min_phi_idx = findfirst(==("min_phi"), probe_header)
    energy_delta_idx = findfirst(==("energy_delta"), probe_header)
    @test minimum(row[min_phi_idx] for row in probe_rows) >= 1.0e-3
    @test minimum(row[energy_delta_idx] for row in probe_rows) < -1.0e4
    probe_text = read(probe.summary_path, String)
    @test occursin("direct-density HE-URP stability probe", probe_text)
    @test occursin("positivity boundary", probe_text)
    @test occursin("full nonlinear target", probe_text)

    basin_probe = run_direct_stability_probe(mktempdir(); Nx=4, C=20.0,
        modes=[(1, 0, 0), (2, 0, 0)], amplitude_count=8, boundary_margin=1.0e-3,
        basin_variance_max=0.02)
    basin_header, basin_rows = DFMMonteCarlo._read_csv_numeric(basin_probe.csv_path)
    @test "basin_variance_max" in basin_header
    @test "inside_basin" in basin_header
    basin_variance_idx = findfirst(==("basin_variance_max"), basin_header)
    inside_basin_idx = findfirst(==("inside_basin"), basin_header)
    density_variance_idx = findfirst(==("density_variance"), basin_header)
    @test all(row[basin_variance_idx] == 0.02 for row in basin_rows)
    @test any(row[inside_basin_idx] == 1.0 for row in basin_rows)
    @test any(row[inside_basin_idx] == 0.0 for row in basin_rows)
    @test all(row[density_variance_idx] <= 0.02 + 1e-12 for row in basin_rows if row[inside_basin_idx] == 1.0)
    basin_text = read(basin_probe.summary_path, String)
    @test occursin("Basin restriction: density_variance <= 0.02", basin_text)
    @test occursin("Basin-restricted scan kept", basin_text)

    basin_scan = run_basin_variance_scan(mktempdir(); Nx=4, C=20.0,
        modes=[(1, 0, 0), (2, 0, 0)], amplitude_count=8,
        variance_caps=[0.0, 0.02], boundary_margin=1.0e-3)
    @test isfile(basin_scan.csv_path)
    @test isfile(basin_scan.summary_path)
    scan_header, scan_rows = DFMMonteCarlo._read_csv_numeric(basin_scan.csv_path)
    @test "basin_variance_max" in scan_header
    @test "in_basin_points" in scan_header
    @test "min_in_basin_energy_delta" in scan_header
    @test "basin_stability_pass" in scan_header
    @test length(scan_rows) == 2
    scan_cap_idx = findfirst(==("basin_variance_max"), scan_header)
    scan_points_idx = findfirst(==("in_basin_points"), scan_header)
    scan_delta_idx = findfirst(==("min_in_basin_energy_delta"), scan_header)
    scan_pass_idx = findfirst(==("basin_stability_pass"), scan_header)
    @test scan_rows[1][scan_cap_idx] == 0.0
    @test scan_rows[1][scan_delta_idx] == 0.0
    @test scan_rows[1][scan_pass_idx] == 1.0
    @test scan_rows[2][scan_cap_idx] == 0.02
    @test scan_rows[2][scan_points_idx] > scan_rows[1][scan_points_idx]
    @test scan_rows[2][scan_delta_idx] < 0.0
    @test scan_rows[2][scan_pass_idx] == 0.0
    scan_text = read(basin_scan.summary_path, String)
    @test occursin("basin variance scan", scan_text)
    @test occursin("largest stable variance cap", scan_text)
end

@testset "Direct phi moves preserve constraints" begin
    rng = StableRNG(123)
    phi = ones(4, 4, 4)
    moved, accepted_geometry = pair_conserving_proposal(rng, phi, 0.05; phi_min=1e-8)
    @test accepted_geometry
    @test mean(moved) ≈ 1.0
    @test minimum(moved) > 0.0

    mode_moved, ok = fourier_mode_proposal(rng, phi, 0.01; phi_min=1e-8)
    @test ok
    @test mean(mode_moved) ≈ 1.0 atol=1e-14
    @test minimum(mode_moved) > 0.0

    grid = Grid3D(4, 6.4)
    params = HEURPParams(A=300.0, B=1.0, C=20.0, alpha=0.1)
    system = DFMMonteCarlo.DensityFieldSystem(copy(phi), grid, params)
    original_energy = system.energy
    @test DFMMonteCarlo._boundary_barrier_strength(nothing) == 0.0
    @test DFMMonteCarlo._boundary_barrier_strength(12.5) == 12.5
    @test_throws ArgumentError DFMMonteCarlo._boundary_barrier_strength(-1.0)
    barrier_strength = 12.5
    barrier_system = DFMMonteCarlo.DensityFieldSystem(copy(phi), grid, params;
        boundary_barrier_strength=barrier_strength)
    @test barrier_system.boundary_barrier_strength == barrier_strength
    @test barrier_system.energy ≈ original_energy
    boundary_phi = copy(phi)
    boundary_phi[1] = 0.1
    boundary_phi[2] = 1.9
    base_boundary_energy = DFMMonteCarlo._target_energy_with_basin(
        boundary_phi, grid, params, Inf, :heurp, 0.0)
    barrier_boundary_energy = DFMMonteCarlo._target_energy_with_basin(
        boundary_phi, grid, params, Inf, :heurp, barrier_strength)
    @test barrier_boundary_energy ≈
          base_boundary_energy + DFMMonteCarlo._boundary_barrier_energy(boundary_phi, barrier_strength)
    @test barrier_boundary_energy > base_boundary_energy
    action = DFMMonteCarlo.FourierModeDensityAction(0, 0, 0, 0.0, 0.0)
    policy = DFMMonteCarlo.FourierModeDensityPolicy()
    Arianna.sample_action!(action, policy, [0.01, 2.0], system, rng)
    @test 0 <= action.ax <= 2
    @test 0 <= action.ay <= 2
    @test 0 <= action.az <= 2
    @test action.ax != 0 || action.ay != 0 || action.az != 0
    @test 0.0 <= action.phase <= 2π
    @test action.amplitude == 0.01
    old_energy, new_energy = Arianna.perform_action!(system, action)
    @test old_energy ≈ original_energy
    @test isfinite(new_energy)
    @test mean(system.phi) ≈ 1.0 atol=1e-14
    @test minimum(system.phi) > 0.0
    @test var(vec(system.phi)) > 0.0
    Arianna.invert_action!(action, system)
    Arianna.perform_action!(system, action)
    @test maximum(abs.(system.phi .- 1.0)) < 1.0e-12
    @test system.energy ≈ original_energy
    @test isfinite(Arianna.log_proposal_density(action, policy, [0.01, 2.0], system))

    global_system = DFMMonteCarlo.DensityFieldSystem(copy(phi), grid, params)
    global_original_energy = global_system.energy
    global_action = DFMMonteCarlo.GlobalGaussianDensityAction(zeros(4, 4, 4))
    global_policy = DFMMonteCarlo.GlobalGaussianDensityPolicy()
    Arianna.sample_action!(global_action, global_policy, [0.01], global_system, rng)
    @test mean(global_action.δ) ≈ 0.0 atol=1.0e-14
    @test std(vec(global_action.δ); corrected=false) ≈ 0.01 rtol=1.0e-12
    global_old_energy, global_new_energy = Arianna.perform_action!(global_system, global_action)
    @test global_old_energy ≈ global_original_energy
    @test isfinite(global_new_energy)
    @test mean(global_system.phi) ≈ 1.0 atol=1.0e-14
    @test minimum(global_system.phi) > 0.0
    @test var(vec(global_system.phi)) > 0.0
    Arianna.invert_action!(global_action, global_system)
    Arianna.perform_action!(global_system, global_action)
    @test maximum(abs.(global_system.phi .- 1.0)) < 1.0e-12
    @test global_system.energy ≈ global_original_energy
    @test isfinite(Arianna.log_proposal_density(global_action, global_policy, [0.01], global_system))

    radial_phi = copy(phi)
    radial_phi[1] += 0.08
    radial_phi[2] -= 0.08
    radial_system = DFMMonteCarlo.DensityFieldSystem(radial_phi, grid, params)
    radial_original_phi = copy(radial_system.phi)
    radial_original_energy = radial_system.energy
    radial_original_variance = var(vec(radial_system.phi))
    radial_action = DFMMonteCarlo.RadialVarianceDensityAction(0.0)
    radial_policy = DFMMonteCarlo.RadialVarianceDensityPolicy()
    Arianna.sample_action!(radial_action, radial_policy, [0.01], radial_system, rng)
    @test abs(radial_action.δ) <= 0.01
    radial_old_energy, radial_new_energy = Arianna.perform_action!(radial_system, radial_action)
    @test radial_old_energy ≈ radial_original_energy
    @test isfinite(radial_new_energy)
    @test mean(radial_system.phi) ≈ 1.0 atol=1.0e-14
    @test minimum(radial_system.phi) > 0.0
    @test var(vec(radial_system.phi)) != radial_original_variance
    Arianna.invert_action!(radial_action, radial_system)
    Arianna.perform_action!(radial_system, radial_action)
    @test maximum(abs.(radial_system.phi .- radial_original_phi)) < 1.0e-12
    @test radial_system.energy ≈ radial_original_energy
    @test isfinite(Arianna.log_proposal_density(radial_action, radial_policy, [0.01], radial_system))

    uniform_radial_system = DFMMonteCarlo.DensityFieldSystem(copy(phi), grid, params)
    uniform_radial_action = DFMMonteCarlo.RadialVarianceDensityAction(0.0)
    Arianna.sample_action!(uniform_radial_action, radial_policy, [0.01], uniform_radial_system, rng)
    uniform_radial_old_energy, uniform_radial_new_energy =
        Arianna.perform_action!(uniform_radial_system, uniform_radial_action)
    @test uniform_radial_old_energy ≈ uniform_radial_new_energy
    @test maximum(abs.(uniform_radial_system.phi .- 1.0)) < 1.0e-12

    uniform_initial = DFMMonteCarlo._initial_density_field(
        StableRNG(124),
        grid;
        initialization=:uniform,
        initial_density_amplitude=0.05,
        phi_min=1.0e-8,
    )
    @test all(uniform_initial .== 1.0)
    overdispersed_initial = DFMMonteCarlo._initial_density_field(
        StableRNG(125),
        grid;
        initialization=:overdispersed,
        initial_density_amplitude=0.05,
        phi_min=1.0e-8,
    )
    @test mean(overdispersed_initial) ≈ 1.0 atol=1.0e-14
    @test minimum(overdispersed_initial) > 1.0e-8
    @test var(vec(overdispersed_initial)) > 0.0
    @test var(vec(overdispersed_initial)) <= 0.05^2 * 1.03
    basin_initial = DFMMonteCarlo._initial_density_field(
        StableRNG(126),
        grid;
        initialization=:overdispersed,
        initial_density_amplitude=0.5,
        basin_variance_max=1.0e-4,
        phi_min=1.0e-8,
    )
    @test var(vec(basin_initial)) <= 1.0e-4 + 1.0e-14

    stable_min_k2 = 2.0 * params.A / params.C
    no_stable_action = DFMMonteCarlo.FourierModeDensityAction(0, 0, 0, 0.0, 0.0)
    @test_throws ArgumentError Arianna.sample_action!(
        no_stable_action,
        policy,
        [0.01, 2.0, stable_min_k2],
        system,
        rng,
    )

    stable_grid = Grid3D(8, 6.4)
    stable_system = DFMMonteCarlo.DensityFieldSystem(ones(8, 8, 8), stable_grid, params)
    stable_action = DFMMonteCarlo.FourierModeDensityAction(0, 0, 0, 0.0, 0.0)
    for _ in 1:20
        Arianna.sample_action!(stable_action, policy, [0.01, 4.0, stable_min_k2], stable_system, rng)
        k2 = (2π / stable_grid.L)^2 *
             (stable_action.ax^2 + stable_action.ay^2 + stable_action.az^2)
        @test k2 >= stable_min_k2 - 1.0e-12
        @test 0 <= stable_action.ax <= 4
        @test 0 <= stable_action.ay <= 4
        @test 0 <= stable_action.az <= 4
    end
    @test isfinite(Arianna.log_proposal_density(stable_action, policy,
        [0.01, 4.0, stable_min_k2], stable_system))
end

@testset "Arianna proposal tuning" begin
    grid = Grid3D(4, 6.4)
    params = HEURPParams(A=300.0, B=1.0, C=30.0)
    @test DFMMonteCarlo._default_tuning_steps(0) == 25
    @test DFMMonteCarlo._default_tuning_steps(1) == 25
    @test DFMMonteCarlo._default_tuning_steps(512) == 512
    @test DFMMonteCarlo._default_tuning_steps(12_800) == 2_000
    @test DFMMonteCarlo._adapt_metropolis_amplitude(0.01, 0.6; target_acceptance=0.3) > 0.01
    @test DFMMonteCarlo._adapt_metropolis_amplitude(0.01, 0.1; target_acceptance=0.3) < 0.01
    tuned = DFMMonteCarlo._tune_pair_move_amplitude(
        grid,
        params;
        initial_amplitude=1.0e-4,
        target_acceptance=0.44,
        rounds=8,
        steps_per_round=20,
        seed=100,
    )
    @test tuned.amplitude > 0.01
    @test 0.2 <= tuned.tuning_acceptance_rate <= 0.8
    constrained_acceptance = DFMMonteCarlo._arianna_pair_acceptance_probe(
        grid,
        params,
        0.05;
        seed=101,
        steps=20,
        basin_variance_max=0.0,
    )
    @test constrained_acceptance == 0.0
    rough_phi = ones(grid.Nx, grid.Nx, grid.Nx)
    rough_phi[1] += 0.2
    rough_phi[2] -= 0.2
    rough_start_acceptance = DFMMonteCarlo._arianna_pair_acceptance_probe(
        grid,
        params,
        0.05;
        seed=102,
        steps=20,
        initial_phi=rough_phi,
    )
    @test 0.0 <= rough_start_acceptance <= 1.0
    @test rough_phi[1] == 1.2
    @test rough_phi[2] == 0.8

    stable_grid = Grid3D(8, 6.4)
    stable_params = HEURPParams(A=300.0, B=1.0, C=20.0, alpha=0.1)
    stable_min_k2 = 2.0 * stable_params.A / stable_params.C
    small_fourier_acceptance = DFMMonteCarlo._arianna_fourier_acceptance_probe(
        stable_grid,
        stable_params,
        1.0e-4,
        4,
        stable_min_k2;
        seed=201,
        steps=20,
    )
    large_fourier_acceptance = DFMMonteCarlo._arianna_fourier_acceptance_probe(
        stable_grid,
        stable_params,
        0.1,
        4,
        stable_min_k2;
        seed=202,
        steps=20,
    )
    @test 0.0 <= large_fourier_acceptance <= small_fourier_acceptance <= 1.0
    rough_stable_phi = ones(stable_grid.Nx, stable_grid.Nx, stable_grid.Nx)
    rough_stable_phi[1] += 0.1
    rough_stable_phi[2] -= 0.1
    rough_fourier_acceptance = DFMMonteCarlo._arianna_fourier_acceptance_probe(
        stable_grid,
        stable_params,
        0.01,
        4,
        stable_min_k2;
        seed=204,
        steps=20,
        initial_phi=rough_stable_phi,
    )
    @test 0.0 <= rough_fourier_acceptance <= 1.0
    @test rough_stable_phi[1] == 1.1
    @test rough_stable_phi[2] == 0.9
    tuned_fourier = DFMMonteCarlo._tune_fourier_move_amplitude(
        stable_grid,
        stable_params;
        initial_amplitude=1.0e-4,
        fourier_max_mode=4,
        fourier_min_k2=stable_min_k2,
        target_acceptance=0.44,
        rounds=6,
        steps_per_round=20,
        seed=203,
    )
    @test tuned_fourier.amplitude > 1.0e-4
    @test isfinite(tuned_fourier.tuning_acceptance_rate)
    @test 0.0 <= tuned_fourier.tuning_acceptance_rate <= 1.0

    warm_phi = ones(grid.Nx, grid.Nx, grid.Nx)
    warmup = DFMMonteCarlo._adaptive_warmup_density_field!(
        warm_phi,
        grid,
        params;
        pair_amplitude=0.01,
        fourier_amplitude=0.0025,
        pair_weight=0.7,
        collective_amplitude_scale=0.5,
        collective_pair_count=4,
        collective_weight=0.2,
        fourier_weight=0.1,
        fourier_max_mode=2,
        fourier_min_k2=0.0,
        stable_mode_proposals=false,
        rounds=1,
        steps_per_round=3,
        target_acceptance=0.44,
        seed=301,
    )
    @test warmup.adaptive_warmup_rounds == 1
    @test warmup.adaptive_warmup_steps == 3
    @test warmup.pair_proposal_amplitude > 0.0
    @test warmup.fourier_mode_amplitude > 0.0
    @test 0.0 <= warmup.adaptive_warmup_acceptance_rate <= 1.0
    @test all(warm_phi .> 0.0)
    @test sum(warm_phi) ≈ length(warm_phi)
end

@testset "Pressure estimator contract" begin
    contract_path = joinpath(@__DIR__, "..", "docs", "pressure_estimator_contract.toml")
    contract = load_pressure_contract(contract_path)
    @test contract.measure == "phi_direct"
    @test contract.dimension == 3
    @test contract.observable_formula_id == "heurp_logvolume_fixed_cv_alpha_central_difference_v3"
    @test "C_times_V" in contract.fixed_quantities
    @test "alpha" in contract.fixed_quantities

    grid = Grid3D(4, 6.4)
    params = HEURPParams(A=300.0, B=1.0, C=20.0, alpha=0.1)
    fields = frozen_pressure_fields(grid)
    values = [pressure_observable(phi, grid, params, contract) for phi in fields]
    summary = batch_summary(values)
    @test length(values) == 3
    @test isfinite(summary.mean)
    @test summary.stderr >= 0.0
    @test pressure_observable(ones(4, 4, 4), grid, params, contract) ≈ params.C + 0.5 * params.B * params.C^2
    @test DFMMonteCarlo._pressure_observable_sparse(ones(4, 4, 4), grid, params, contract) ≈
          pressure_observable(ones(4, 4, 4), grid, params, contract)
    nonuniform = fields[end]
    eps = 1.0e-4
    volume_scale = exp(eps)
    length_scale = exp(eps / grid.dimension)
    plus_grid = Grid3D(grid.Nx, grid.L * length_scale)
    minus_grid = Grid3D(grid.Nx, grid.L / length_scale)
    plus_params = HEURPParams(A=params.A, B=params.B, C=params.C / volume_scale, alpha=params.alpha)
    minus_params = HEURPParams(A=params.A, B=params.B, C=params.C * volume_scale, alpha=params.alpha)
    dH_dlogV_fixed_chains = (heurp_energy(nonuniform, plus_grid, plus_params) -
                             heurp_energy(nonuniform, minus_grid, minus_params)) / (2.0 * eps)
    expected_pressure = params.C - dH_dlogV_fixed_chains / grid.volume
    @test pressure_observable(nonuniform, grid, params, contract) ≈ expected_pressure
    @test DFMMonteCarlo._pressure_observable_sparse(nonuniform, grid, params, contract) ≈ expected_pressure
    @test DFMMonteCarlo._pressure_observable_fft(nonuniform, grid, params, contract) ≈ expected_pressure

    bad = PressureContract(contract.measure, 2, contract.fourier_normalization,
        contract.volume_derivative, contract.fixed_quantities,
        contract.observable_formula_id, contract.cl_column_mapping)
    @test_throws ArgumentError pressure_observable(fields[1], grid, params, bad)
    bad_norm = PressureContract(contract.measure, contract.dimension, "changed_normalization",
        contract.volume_derivative, contract.fixed_quantities,
        contract.observable_formula_id, contract.cl_column_mapping)
    @test_throws ArgumentError pressure_observable(fields[1], grid, params, bad_norm)
    bad_derivative = PressureContract(contract.measure, contract.dimension, contract.fourier_normalization,
        "changed_volume_derivative", contract.fixed_quantities,
        contract.observable_formula_id, contract.cl_column_mapping)
    @test_throws ArgumentError pressure_observable(fields[1], grid, params, bad_derivative)
    missing_fixed = PressureContract(contract.measure, contract.dimension, contract.fourier_normalization,
        contract.volume_derivative, ["A", "B", "C"],
        contract.observable_formula_id, contract.cl_column_mapping)
    @test_throws ArgumentError pressure_observable(fields[1], grid, params, missing_fixed)
    bad_formula = PressureContract(contract.measure, contract.dimension, contract.fourier_normalization,
        contract.volume_derivative, contract.fixed_quantities,
        "heurp_logvolume_fixed_cv_central_difference_v2", contract.cl_column_mapping)
    @test_throws ArgumentError pressure_observable(fields[1], grid, params, bad_formula)
end

@testset "Trace uncertainty" begin
    summary = batch_summary([1.0, 2.0, 3.0, 4.0])
    @test summary.mean == 2.5
    @test isapprox(summary.stderr, 1.0, rtol=1e-14)
    @test summary.integrated_autocorrelation_time >= 0.5
    @test 1.0 <= summary.effective_sample_size <= summary.n
    alternating = batch_summary([1.0, -1.0, 1.0, -1.0, 1.0, -1.0])
    @test alternating.integrated_autocorrelation_time == 0.5
    @test alternating.effective_sample_size == alternating.n
    correlated = batch_summary(vcat(fill(1.0, 8), fill(3.0, 8)))
    @test correlated.mean == 2.0
    @test correlated.stderr > 0.0
    @test correlated.integrated_autocorrelation_time > 0.5
    @test correlated.effective_sample_size < correlated.n
    @test correlated.stderr >= correlated.independent_stderr
    @test DFMMonteCarlo._trace_drift_zscore([1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0]) == 0.0
    @test isinf(DFMMonteCarlo._trace_drift_zscore(vcat(fill(1.0, 8), fill(10.0, 8))))
    @test isnan(DFMMonteCarlo._trace_drift_zscore([1.0, 2.0, 3.0]))
    @test DFMMonteCarlo._drift_pass((trace_drift_zscore=NaN,), 3.0) == 1.0
    @test DFMMonteCarlo._drift_pass((trace_drift_zscore=2.0,), 3.0) == 1.0
    @test DFMMonteCarlo._drift_pass((trace_drift_zscore=4.0,), 3.0) == 0.0
    @test DFMMonteCarlo._drift_pass((trace_drift_zscore=Inf,), 3.0) == 0.0
    chain_offset_summaries = [batch_summary(fill(1.0, 20)), batch_summary(fill(2.0, 20))]
    @test DFMMonteCarlo._combined_chain_trace_drift_zscore(chain_offset_summaries) == 0.0
    @test isinf(DFMMonteCarlo._trace_drift_zscore(vcat(fill(1.0, 20), fill(2.0, 20))))

    drifting_traces = [[1.0, 1.0, 1.0, 10.0, 10.0, 10.0],
                       [1.0, 1.0, 1.0, 10.0, 10.0, 10.0]]
    @test DFMMonteCarlo._chain_rhat(drifting_traces) == 1.0
    @test isinf(DFMMonteCarlo._split_chain_rhat(drifting_traces))
    stable_traces = [[1.0, 2.0, 1.0, 2.0], [1.0, 2.0, 1.0, 2.0]]
    @test DFMMonteCarlo._split_chain_rhat(stable_traces) == 1.0
    @test isnan(DFMMonteCarlo._split_chain_rhat([[1.0, 2.0, 3.0]]))
    @test DFMMonteCarlo._rhat_value_pass(NaN, 1.1) == 1.0
    @test DFMMonteCarlo._rhat_value_pass(1.0, 1.1) == 1.0
    @test DFMMonteCarlo._rhat_value_pass(1.2, 1.1) == 0.0
    @test DFMMonteCarlo._rhat_value_pass(Inf, 1.1) == 0.0

    @test DFMMonteCarlo._trace_samples_pass((n=4,), 5) == 0.0
    @test DFMMonteCarlo._trace_samples_pass((n=5,), 5) == 1.0
    @test DFMMonteCarlo._pressure_error_pass((stderr=0.2,), 0.1) == 0.0
    @test DFMMonteCarlo._pressure_error_pass((stderr=0.05,), 0.1) == 1.0
    @test DFMMonteCarlo._pressure_error_pass((stderr=0.2,), Inf) == 1.0
    steady_summary = batch_summary(fill(1.0, 20))
    steady_equilibration = DFMMonteCarlo._equilibration_sensitivity_diagnostics(
        fill(1.0, 20),
        steady_summary;
        threshold=3.0,
    )
    @test steady_equilibration.equilibration_tail_fraction == 0.5
    @test steady_equilibration.equilibration_tail_sample_count == 10
    @test steady_equilibration.equilibration_tail_mean == 1.0
    @test steady_equilibration.equilibration_sensitivity_zscore == 0.0
    @test steady_equilibration.equilibration_sensitivity_pass == 1.0
    shifted_trace = vcat(fill(1.0, 20), fill(10.0, 20))
    shifted_equilibration = DFMMonteCarlo._equilibration_sensitivity_diagnostics(
        shifted_trace,
        batch_summary(shifted_trace);
        threshold=3.0,
    )
    @test shifted_equilibration.equilibration_tail_sample_count == 20
    @test shifted_equilibration.equilibration_tail_mean == 10.0
    @test shifted_equilibration.equilibration_sensitivity_zscore > 3.0
    @test shifted_equilibration.equilibration_sensitivity_pass == 0.0
    chain_a_equilibration = DFMMonteCarlo._equilibration_sensitivity_diagnostics(
        fill(1.0, 20),
        batch_summary(fill(1.0, 20));
        threshold=3.0,
    )
    chain_b_equilibration = DFMMonteCarlo._equilibration_sensitivity_diagnostics(
        fill(10.0, 20),
        batch_summary(fill(10.0, 20));
        threshold=3.0,
    )
    combined_equilibration = DFMMonteCarlo._combine_equilibration_sensitivity_diagnostics(
        [chain_a_equilibration, chain_b_equilibration],
    )
    @test combined_equilibration.equilibration_sensitivity_zscore == 0.0
    @test combined_equilibration.equilibration_sensitivity_pass == 1.0

    budget = DFMMonteCarlo._sampling_budget_estimate(
        (stderr=0.2, integrated_autocorrelation_time=2.0, effective_sample_size=25.0, n=100),
        production_steps=200,
        proposals_per_sweep=64,
        burn_fraction=0.5,
        sample_stride=2,
        ess_target=50.0,
        min_trace_samples=20,
        pressure_tolerance=0.1,
    )
    @test budget.pressure_target_effective_sample_size == 100.0
    @test budget.required_effective_sample_size == 100.0
    @test budget.required_trace_samples == 400
    @test budget.estimated_required_steps == 1600
    @test budget.estimated_required_sweeps == 25.0
    @test budget.production_budget_shortfall == 8.0
    @test budget.sampling_budget_pass == 0.0
    replicated_budget = DFMMonteCarlo._sampling_budget_estimate(
        (stderr=0.2, integrated_autocorrelation_time=2.0, effective_sample_size=25.0, n=400),
        production_steps=200,
        proposals_per_sweep=64,
        burn_fraction=0.5,
        sample_stride=2,
        ess_target=50.0,
        min_trace_samples=20,
        pressure_tolerance=0.1,
        chain_count=4,
    )
    @test replicated_budget.required_trace_samples == 400
    @test replicated_budget.estimated_required_steps == 400
    @test replicated_budget.production_budget_shortfall == 2.0
    @test replicated_budget.sampling_budget_pass == 0.0
    decorrelation = DFMMonteCarlo._decorrelation_diagnostics(
        (integrated_autocorrelation_time=2.0, n=100);
        sample_stride=3,
        production_steps=420,
        equilibration_steps=120,
        required_effective_sample_size=20.0,
    )
    @test decorrelation.decorrelation_stride_samples == 4
    @test decorrelation.decorrelation_stride_steps == 12
    @test decorrelation.decorrelated_sample_count == 25
    @test decorrelation.equilibration_autocorrelation_windows == 10.0
    @test decorrelation.production_autocorrelation_windows == 25.0
    @test decorrelation.decorrelation_pass == 1.0
    decorrelation_short = DFMMonteCarlo._decorrelation_diagnostics(
        (integrated_autocorrelation_time=2.0, n=100);
        sample_stride=3,
        production_steps=420,
        equilibration_steps=120,
        required_effective_sample_size=30.0,
    )
    @test decorrelation_short.decorrelation_pass == 0.0

    ready_budget = DFMMonteCarlo._sampling_budget_estimate(
        (stderr=0.05, integrated_autocorrelation_time=0.5, effective_sample_size=100.0, n=100),
        production_steps=200,
        proposals_per_sweep=64,
        burn_fraction=0.0,
        sample_stride=1,
        ess_target=50.0,
        min_trace_samples=20,
        pressure_tolerance=0.1,
    )
    @test ready_budget.pressure_target_effective_sample_size == 25.0
    @test ready_budget.required_effective_sample_size == 50.0
    @test ready_budget.required_trace_samples == 50
    @test ready_budget.estimated_required_steps == 50
    @test ready_budget.sampling_budget_pass == 1.0
end

@testset "CL parser and run configuration" begin
    cl_path = joinpath(@__DIR__, "..", "..", "Analysis", "GREM_CL_a0.1_B1.0_L6.4_Nx64_dt0.1ETD.dat")
    data = read_cl_pressure(cl_path)
    @test length(data.C) >= 20
    @test data.C[1] == 0.1
    @test data.pressure[1] ≈ 0.10442656253
    @test data.pressure_error[1] ≈ 2.0578125849e-5

    analysis = analysis_reference_config()
    @test analysis.mode == :analysis
    @test analysis.alpha == 0.1
    @test analysis.A == 300.0
    @test analysis.B == 1.0
    @test analysis.L == 6.4
    @test analysis.target_Nx == 64
    @test analysis.C_values == data.C[11:end]
    @test first(analysis.C_values) == 19.0

    quick = run_config(:quick)
    production = run_config(:production)
    @test quick.dimension == 3
    @test quick.manuscript_grade == false
    @test production.dimension == 3
    @test production.target_Nx == 64
    @test production.manuscript_grade == true
    @test production.alpha == 0.1
end

@testset "Analysis figure comparison suite" begin
    specs = analysis_figure_specs()
    names = [spec.name for spec in specs]
    @test "PivsC_B1_A300_50" in names
    @test "dPivsC_B1_A10_50" in names
    @test "PivsA_B1_C20" in names

    pivsc = specs[findfirst(==("PivsC_B1_A300_50"), names)]
    @test pivsc.kind == :PivsC
    @test pivsc.A == 300.0
    @test pivsc.B == 1.0
    @test pivsc.Cmax == 50.0

    outdir = mktempdir()
    results = run_analysis_comparison_suite(outdir; steps=4, Nx=2, max_points_per_figure=1,
        figure_names=["PivsC_B1_A300_50"], burn_fraction=0.25, sample_stride=2)
    @test length(results) == 1
    @test isfile(results[1])
    header, rows = DFMMonteCarlo._read_csv_numeric(results[1])
    @test "mc_pressure" in header
    @test "mean_field" in header
    @test "one_loop" in header
    @test "dfm_pressure" in header
    @test "cl_pressure" in header
    @test "collective_proposal_amplitude" in header
    @test "collective_pair_count" in header
    @test "collective_move_weight" in header
    @test "fourier_mode_amplitude" in header
    @test "fourier_max_mode" in header
    @test "fourier_move_weight" in header
    @test "ess_target" in header
    @test "converged" in header
    @test "burn_fraction" in header
    @test "sample_stride_sweeps" in header
    @test "production_steps" in header
    @test "max_production_steps" in header
    @test "proposals_per_sweep" in header
    @test "production_sweeps" in header
    @test "max_production_sweeps" in header
    @test "convergence_attempts" in header
    @test "alpha" in header
    @test length(rows) == 1
    analysis_burn_idx = findfirst(==("burn_fraction"), header)
    analysis_stride_idx = findfirst(==("sample_stride"), header)
    analysis_stride_sweeps_idx = findfirst(==("sample_stride_sweeps"), header)
    analysis_proposals_per_sweep_idx = findfirst(==("proposals_per_sweep"), header)
    @test rows[1][analysis_burn_idx] == 0.25
    @test rows[1][analysis_stride_idx] == 2.0
    @test rows[1][analysis_stride_sweeps_idx] ≈
          rows[1][analysis_stride_idx] / rows[1][analysis_proposals_per_sweep_idx]

    dpivsc_paths = run_dpivsc_all_cl_comparisons(outdir; steps=1, Nx=2,
        figure_names=["dPivsC_B1_A10_50"], burn_fraction=0.25, sample_stride=2)
    @test length(dpivsc_paths) == 1
    dp_header, dp_rows = DFMMonteCarlo._read_csv_numeric(dpivsc_paths[1])
    @test "mc_pressure" in dp_header
    @test "collective_proposal_amplitude" in dp_header
    @test "collective_pair_count" in dp_header
    @test "collective_move_weight" in dp_header
    @test "fourier_mode_amplitude" in dp_header
    @test "fourier_max_mode" in dp_header
    @test "fourier_move_weight" in dp_header
    @test "ess_target" in dp_header
    @test "converged" in dp_header
    @test "burn_fraction" in dp_header
    @test "sample_stride_sweeps" in dp_header
    @test "production_steps" in dp_header
    @test "max_production_steps" in dp_header
    @test "proposals_per_sweep" in dp_header
    @test "production_sweeps" in dp_header
    @test "max_production_sweeps" in dp_header
    @test "convergence_attempts" in dp_header
    @test "alpha" in dp_header
    dp_cl_path = joinpath(@__DIR__, "..", "..", "Analysis", "GREM_CL_a0.1_B1.0_L6.4_Nx64_dt0.1ETD.dat")
    @test length(dp_rows) == length(read_cl_pressure(dp_cl_path).C)
    dp_svg = joinpath(outdir, "figures", "dPivsC_B1_A10_50_mc_comparison.svg")
    svg_text = read(dp_svg, String)
    @test occursin("class='plot-title'", svg_text)
    @test occursin("class='x-axis-label'", svg_text)
    @test occursin("class='y-axis-label'", svg_text)
    @test occursin("class='x-tick-label'", svg_text)
    @test occursin("class='y-tick-label'", svg_text)
    @test length(collect(eachmatch(r"class='x-major-tick'", svg_text))) >= 9
    @test length(collect(eachmatch(r"class='y-major-tick'", svg_text))) >= 9
    @test length(collect(eachmatch(r"class='x-tick-label'", svg_text))) >= 9
    @test length(collect(eachmatch(r"class='y-tick-label'", svg_text))) >= 9
    @test length(collect(eachmatch(r"class='x-minor-tick'", svg_text))) >= 32
    @test length(collect(eachmatch(r"class='y-minor-tick'", svg_text))) >= 32
end

@testset "Arianna spike and artifact generation" begin
    outdir = mktempdir()
    spike = run_arianna_spike(outdir; seed=7, steps=20)
    @test spike.success
    @test isfile(spike.output_path)
    @test spike.acceptance_rate >= 0.0
    @test spike.acceptance_rate <= 1.0
    spike_header, spike_rows = DFMMonteCarlo._read_csv_strings(spike.output_path)
    @test "low_mode_structure_factor" in spike_header
    spike_structure_idx = findfirst(==("low_mode_structure_factor"), spike_header)
    @test all(parse(Float64, row[spike_structure_idx]) >= 0.0 for row in spike_rows)

    result = run_pressure_sweep(outdir; mode=:quick, seed=9, steps=6, ess_target=1.0, max_steps=12,
        chains=2, min_acceptance=0.2, max_acceptance=0.8, rhat_threshold=1.1, min_trace_samples=1,
        burn_fraction=0.25, sample_stride=2)
    @test isfile(result.pressure_csv)
    @test isfile(result.diagnostics_csv)
    metadata_path = joinpath(outdir, "run_metadata.toml")
    @test isfile(metadata_path)
    metadata = TOML.parsefile(metadata_path)
    @test metadata["run"]["sampler"] == "arianna"
    @test metadata["run"]["mode"] == "quick"
    @test metadata["run"]["seed"] == 9
    @test metadata["grid"]["dimension"] == 3
    @test metadata["grid"]["Nx"] == result.config.Nx
    @test metadata["grid"]["target_Nx"] == result.config.target_Nx
    @test metadata["params"]["A"] == 300.0
    @test metadata["params"]["B"] == 1.0
    @test metadata["params"]["alpha"] == result.config.alpha
    @test metadata["sampling"]["chains"] == 2
    @test metadata["sampling"]["steps"] == 6
    @test metadata["sampling"]["max_steps"] == 12
    @test metadata["sampling"]["burn_fraction"] == 0.25
    @test metadata["sampling"]["sample_stride"] == 2
    @test metadata["sampling"]["phi_min"] == 1.0e-8
    @test metadata["sampling"]["global_move_weight"] == 0.0
    @test metadata["sampling"]["global_amplitude_scale"] == 0.1
    @test metadata["run"]["target_model"] == "HE-URP"
    @test metadata["pressure_contract"]["observable_formula_id"] == "heurp_logvolume_fixed_cv_alpha_central_difference_v3"
    @test metadata["pressure_contract"]["measure"] == "phi_direct"
    @test metadata["environment"]["package"] == "DFMMonteCarlo"
    @test haskey(metadata["environment"], "julia_version")
    @test metadata["outputs"]["pressure_csv"] == result.pressure_csv
    @test metadata["outputs"]["diagnostics_csv"] == result.diagnostics_csv
    @test metadata["outputs"]["observable_traces_csv"] == joinpath(outdir, "observable_traces.csv")
    @test result.config.dimension == 3
    header, rows = DFMMonteCarlo._read_csv_numeric(result.pressure_csv)
    @test "mc_pressure" in header
    @test "raw_pressure_observable" in header
    @test "integrated_autocorrelation_time" in header
    @test "effective_sample_size" in header
    @test "trace_drift_zscore" in header
    @test "ess_target" in header
    @test "converged" in header
    @test "burn_fraction" in header
    @test "sample_stride_sweeps" in header
    @test "production_steps" in header
    @test "max_production_steps" in header
    @test "proposals_per_sweep" in header
    @test "production_sweeps" in header
    @test "max_production_sweeps" in header
    @test "convergence_attempts" in header
    @test "chain_count" in header
    @test "chains_converged" in header
    @test "chain_rhat" in header
    @test "split_chain_rhat" in header
    @test "chain_between_stderr" in header
    @test "chain_within_stderr" in header
    @test "acceptance_min" in header
    @test "acceptance_max" in header
    @test "acceptance_pass" in header
    @test "drift_threshold" in header
    @test "drift_pass" in header
    @test "rhat_threshold" in header
    @test "rhat_pass" in header
    @test "ess_pass" in header
    @test "nonlinear_stability_pass" in header
    @test "nonlinear_stability_min_energy_delta" in header
    @test "nonlinear_stability_min_phi" in header
    @test "min_trace_samples" in header
    @test "trace_samples_pass" in header
    @test "pressure_tolerance" in header
    @test "pressure_error_pass" in header
    @test "pressure_target_effective_sample_size" in header
    @test "required_effective_sample_size" in header
    @test "required_trace_samples" in header
    @test "estimated_required_steps" in header
    @test "estimated_required_sweeps" in header
    @test "production_budget_shortfall" in header
    @test "sampling_budget_pass" in header
    @test "adaptive_warmup_rounds" in header
    @test "adaptive_warmup_steps" in header
    @test "adaptive_warmup_acceptance_rate" in header
    @test "adaptive_warmup_pair_acceptance_rate" in header
    @test "adaptive_warmup_collective_acceptance_rate" in header
    @test "adaptive_warmup_fourier_acceptance_rate" in header
    @test "diagnostic_pass" in header
    @test "alpha" in header
    @test "target_model_code" in header
    dfm_idx = findfirst(==("dfm_pressure"), header)
    @test all(row[dfm_idx] ≈ dfm_pressure_discrete(300.0, 1.0, row[1], Grid3D(result.config.Nx, result.config.L);
        alpha=result.config.alpha) for row in rows)
    @test any(abs(row[2] - row[4]) > 1.0 for row in rows)
    alpha_idx = findfirst(==("alpha"), header)
    @test all(row[alpha_idx] == result.config.alpha for row in rows)
    target_model_idx = findfirst(==("target_model_code"), header)
    @test all(row[target_model_idx] == 1.0 for row in rows)
    diag_header, diag_rows = DFMMonteCarlo._read_csv_strings(result.diagnostics_csv)
    @test "sampler" in diag_header
    @test "equilibration_steps" in diag_header
    @test "sample_stride" in diag_header
    @test "trace_samples" in diag_header
    @test "integrated_autocorrelation_time" in diag_header
    @test "effective_sample_size" in diag_header
    @test "trace_drift_zscore" in diag_header
    @test "proposal_amplitude" in diag_header
    @test "pair_proposal_amplitude" in diag_header
    @test "pair_move_weight" in diag_header
    @test "pair_acceptance_rate" in diag_header
    @test "pair_tuning_acceptance_rate" in diag_header
    @test "collective_acceptance_rate" in diag_header
    @test "fourier_acceptance_rate" in diag_header
    @test "global_acceptance_rate" in diag_header
    @test "target_acceptance" in diag_header
    @test "tuning_rounds" in diag_header
    @test "tuning_steps" in diag_header
    @test "tuning_acceptance_rate" in diag_header
    @test "adaptive_warmup_rounds" in diag_header
    @test "adaptive_warmup_steps" in diag_header
    @test "adaptive_warmup_acceptance_rate" in diag_header
    @test "adaptive_warmup_pair_acceptance_rate" in diag_header
    @test "adaptive_warmup_collective_acceptance_rate" in diag_header
    @test "adaptive_warmup_fourier_acceptance_rate" in diag_header
    @test "collective_proposal_amplitude" in diag_header
    @test "collective_pair_count" in diag_header
    @test "collective_move_weight" in diag_header
    @test "fourier_mode_amplitude" in diag_header
    @test "fourier_max_mode" in diag_header
    @test "fourier_move_weight" in diag_header
    @test "global_field_amplitude" in diag_header
    @test "global_move_weight" in diag_header
    @test "fourier_tuning_acceptance_rate" in diag_header
    @test "ess_target" in diag_header
    @test "converged" in diag_header
    @test "burn_fraction" in diag_header
    @test "sample_stride_sweeps" in diag_header
    @test "production_steps" in diag_header
    @test "max_production_steps" in diag_header
    @test "proposals_per_sweep" in diag_header
    @test "production_sweeps" in diag_header
    @test "max_production_sweeps" in diag_header
    @test "equilibration_sweeps" in diag_header
    @test "convergence_attempts" in diag_header
    @test "chain_count" in diag_header
    @test "chains_converged" in diag_header
    @test "chain_rhat" in diag_header
    @test "split_chain_rhat" in diag_header
    @test "chain_between_stderr" in diag_header
    @test "chain_within_stderr" in diag_header
    @test "acceptance_min" in diag_header
    @test "acceptance_max" in diag_header
    @test "acceptance_pass" in diag_header
    @test "drift_threshold" in diag_header
    @test "drift_pass" in diag_header
    @test "rhat_threshold" in diag_header
    @test "rhat_pass" in diag_header
    @test "ess_pass" in diag_header
    @test "nonlinear_stability_pass" in diag_header
    @test "nonlinear_stability_min_energy_delta" in diag_header
    @test "nonlinear_stability_min_phi" in diag_header
    @test "min_trace_samples" in diag_header
    @test "trace_samples_pass" in diag_header
    @test "pressure_tolerance" in diag_header
    @test "pressure_error_pass" in diag_header
    @test "pressure_target_effective_sample_size" in diag_header
    @test "required_effective_sample_size" in diag_header
    @test "required_trace_samples" in diag_header
    @test "estimated_required_steps" in diag_header
    @test "estimated_required_sweeps" in diag_header
    @test "production_budget_shortfall" in diag_header
    @test "sampling_budget_pass" in diag_header
    @test "decorrelation_pass" in diag_header
    @test "equilibration_sensitivity_zscore" in diag_header
    @test "equilibration_sensitivity_pass" in diag_header
    @test "diagnostic_pass" in diag_header
    @test "low_mode_structure_factor" in diag_header
    @test "alpha" in diag_header
    @test "target_model_code" in diag_header
    sampler_idx = findfirst(==("sampler"), diag_header)
    @test all(row[sampler_idx] == "arianna" for row in diag_rows)
    trace_samples_idx = findfirst(==("trace_samples"), diag_header)
    ess_idx = findfirst(==("effective_sample_size"), diag_header)
    drift_idx = findfirst(==("trace_drift_zscore"), diag_header)
    amplitude_idx = findfirst(==("proposal_amplitude"), diag_header)
    target_idx = findfirst(==("target_acceptance"), diag_header)
    tuning_idx = findfirst(==("tuning_rounds"), diag_header)
    tuning_steps_idx = findfirst(==("tuning_steps"), diag_header)
    adaptive_rounds_idx = findfirst(==("adaptive_warmup_rounds"), diag_header)
    adaptive_steps_idx = findfirst(==("adaptive_warmup_steps"), diag_header)
    adaptive_acceptance_idx = findfirst(==("adaptive_warmup_acceptance_rate"), diag_header)
    adaptive_pair_acceptance_idx = findfirst(==("adaptive_warmup_pair_acceptance_rate"), diag_header)
    adaptive_collective_acceptance_idx = findfirst(==("adaptive_warmup_collective_acceptance_rate"), diag_header)
    adaptive_fourier_acceptance_idx = findfirst(==("adaptive_warmup_fourier_acceptance_rate"), diag_header)
    collective_amplitude_idx = findfirst(==("collective_proposal_amplitude"), diag_header)
    collective_pair_count_idx = findfirst(==("collective_pair_count"), diag_header)
    collective_idx = findfirst(==("collective_move_weight"), diag_header)
    fourier_amplitude_idx = findfirst(==("fourier_mode_amplitude"), diag_header)
    fourier_max_mode_idx = findfirst(==("fourier_max_mode"), diag_header)
    fourier_idx = findfirst(==("fourier_move_weight"), diag_header)
    global_acceptance_idx = findfirst(==("global_acceptance_rate"), diag_header)
    global_amplitude_idx = findfirst(==("global_field_amplitude"), diag_header)
    global_idx = findfirst(==("global_move_weight"), diag_header)
    ess_target_idx = findfirst(==("ess_target"), diag_header)
    converged_idx = findfirst(==("converged"), diag_header)
    burn_fraction_idx = findfirst(==("burn_fraction"), diag_header)
    sample_stride_idx = findfirst(==("sample_stride"), diag_header)
    sample_stride_sweeps_idx = findfirst(==("sample_stride_sweeps"), diag_header)
    production_steps_idx = findfirst(==("production_steps"), diag_header)
    max_steps_idx = findfirst(==("max_production_steps"), diag_header)
    proposals_per_sweep_idx = findfirst(==("proposals_per_sweep"), diag_header)
    production_sweeps_idx = findfirst(==("production_sweeps"), diag_header)
    max_sweeps_idx = findfirst(==("max_production_sweeps"), diag_header)
    equilibration_sweeps_idx = findfirst(==("equilibration_sweeps"), diag_header)
    attempts_idx = findfirst(==("convergence_attempts"), diag_header)
    chain_count_idx = findfirst(==("chain_count"), diag_header)
    chains_converged_idx = findfirst(==("chains_converged"), diag_header)
    chain_rhat_idx = findfirst(==("chain_rhat"), diag_header)
    split_chain_rhat_idx = findfirst(==("split_chain_rhat"), diag_header)
    chain_between_idx = findfirst(==("chain_between_stderr"), diag_header)
    chain_within_idx = findfirst(==("chain_within_stderr"), diag_header)
    acceptance_min_idx = findfirst(==("acceptance_min"), diag_header)
    acceptance_max_idx = findfirst(==("acceptance_max"), diag_header)
    acceptance_pass_idx = findfirst(==("acceptance_pass"), diag_header)
    drift_threshold_idx = findfirst(==("drift_threshold"), diag_header)
    drift_pass_idx = findfirst(==("drift_pass"), diag_header)
    rhat_threshold_idx = findfirst(==("rhat_threshold"), diag_header)
    rhat_pass_idx = findfirst(==("rhat_pass"), diag_header)
    ess_pass_idx = findfirst(==("ess_pass"), diag_header)
    stability_pass_idx = findfirst(==("nonlinear_stability_pass"), diag_header)
    stability_delta_idx = findfirst(==("nonlinear_stability_min_energy_delta"), diag_header)
    stability_min_phi_idx = findfirst(==("nonlinear_stability_min_phi"), diag_header)
    min_trace_samples_idx = findfirst(==("min_trace_samples"), diag_header)
    trace_samples_pass_idx = findfirst(==("trace_samples_pass"), diag_header)
    pressure_tolerance_idx = findfirst(==("pressure_tolerance"), diag_header)
    pressure_error_pass_idx = findfirst(==("pressure_error_pass"), diag_header)
    pressure_target_ess_idx = findfirst(==("pressure_target_effective_sample_size"), diag_header)
    required_ess_idx = findfirst(==("required_effective_sample_size"), diag_header)
    required_trace_idx = findfirst(==("required_trace_samples"), diag_header)
    required_steps_idx = findfirst(==("estimated_required_steps"), diag_header)
    required_sweeps_idx = findfirst(==("estimated_required_sweeps"), diag_header)
    shortfall_idx = findfirst(==("production_budget_shortfall"), diag_header)
    budget_pass_idx = findfirst(==("sampling_budget_pass"), diag_header)
    decorrelation_pass_idx = findfirst(==("decorrelation_pass"), diag_header)
    equilibration_sensitivity_idx = findfirst(==("equilibration_sensitivity_zscore"), diag_header)
    equilibration_sensitivity_pass_idx = findfirst(==("equilibration_sensitivity_pass"), diag_header)
    diagnostic_pass_idx = findfirst(==("diagnostic_pass"), diag_header)
    structure_idx = findfirst(==("low_mode_structure_factor"), diag_header)
    diag_alpha_idx = findfirst(==("alpha"), diag_header)
    diag_target_model_idx = findfirst(==("target_model_code"), diag_header)
    @test all(parse(Float64, row[trace_samples_idx]) >= 1.0 for row in diag_rows)
    @test all(parse(Float64, row[ess_idx]) >= 1.0 for row in diag_rows)
    @test all(parse(Float64, row[drift_idx]) >= 0.0 || isnan(parse(Float64, row[drift_idx])) for row in diag_rows)
    @test all(parse(Float64, row[ess_target_idx]) == 1.0 for row in diag_rows)
    @test all(parse(Float64, row[converged_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(parse(Float64, row[burn_fraction_idx]) == 0.25 for row in diag_rows)
    @test all(parse(Float64, row[sample_stride_idx]) == 2.0 for row in diag_rows)
    @test all(6.0 <= parse(Float64, row[production_steps_idx]) <= 12.0 for row in diag_rows)
    @test all(parse(Float64, row[max_steps_idx]) == 12.0 for row in diag_rows)
    @test all(parse(Float64, row[proposals_per_sweep_idx]) == result.config.Nx^3 for row in diag_rows)
    @test all(parse(Float64, row[production_sweeps_idx]) ≈
              parse(Float64, row[production_steps_idx]) / parse(Float64, row[proposals_per_sweep_idx]) for row in diag_rows)
    @test all(parse(Float64, row[sample_stride_sweeps_idx]) ≈
              parse(Float64, row[sample_stride_idx]) / parse(Float64, row[proposals_per_sweep_idx]) for row in diag_rows)
    @test all(parse(Float64, row[max_sweeps_idx]) ≈
              parse(Float64, row[max_steps_idx]) / parse(Float64, row[proposals_per_sweep_idx]) for row in diag_rows)
    @test all(parse(Float64, row[equilibration_sweeps_idx]) >= 0.0 for row in diag_rows)
    @test all(parse(Float64, row[attempts_idx]) >= 1.0 for row in diag_rows)
    @test all(parse(Float64, row[converged_idx]) == 0.0 ||
              parse(Float64, row[ess_idx]) >= parse(Float64, row[ess_target_idx]) for row in diag_rows)
    @test all(parse(Float64, row[chain_count_idx]) == 2.0 for row in diag_rows)
    @test all(0.0 <= parse(Float64, row[chains_converged_idx]) <= 2.0 for row in diag_rows)
    @test all(parse(Float64, row[equilibration_sensitivity_idx]) >= 0.0 ||
              isnan(parse(Float64, row[equilibration_sensitivity_idx])) for row in diag_rows)
    @test all(parse(Float64, row[equilibration_sensitivity_pass_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(parse(Float64, row[chain_rhat_idx]) >= 1.0 for row in diag_rows)
    @test all(parse(Float64, row[split_chain_rhat_idx]) >= 1.0 ||
              isnan(parse(Float64, row[split_chain_rhat_idx])) for row in diag_rows)
    @test all(parse(Float64, row[chain_between_idx]) >= 0.0 for row in diag_rows)
    @test all(parse(Float64, row[chain_within_idx]) >= 0.0 for row in diag_rows)
    @test all(parse(Float64, row[acceptance_min_idx]) == 0.2 for row in diag_rows)
    @test all(parse(Float64, row[acceptance_max_idx]) == 0.8 for row in diag_rows)
    @test all(parse(Float64, row[drift_threshold_idx]) == 3.0 for row in diag_rows)
    @test all(parse(Float64, row[rhat_threshold_idx]) == 1.1 for row in diag_rows)
    @test all(parse(Float64, row[acceptance_pass_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(parse(Float64, row[drift_pass_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(parse(Float64, row[rhat_pass_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(parse(Float64, row[ess_pass_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(parse(Float64, row[stability_pass_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(isfinite(parse(Float64, row[stability_delta_idx])) for row in diag_rows)
    @test all(parse(Float64, row[stability_min_phi_idx]) > 0.0 for row in diag_rows)
    @test all(parse(Float64, row[min_trace_samples_idx]) == 1.0 for row in diag_rows)
    @test all(parse(Float64, row[trace_samples_pass_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(isinf(parse(Float64, row[pressure_tolerance_idx])) for row in diag_rows)
    @test all(parse(Float64, row[pressure_error_pass_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(parse(Float64, row[pressure_target_ess_idx]) == 0.0 for row in diag_rows)
    @test all(parse(Float64, row[required_ess_idx]) >= parse(Float64, row[ess_target_idx]) for row in diag_rows)
    @test all(parse(Float64, row[required_trace_idx]) >= parse(Float64, row[min_trace_samples_idx]) for row in diag_rows)
    @test all(parse(Float64, row[required_steps_idx]) >= parse(Float64, row[required_trace_idx]) for row in diag_rows)
    @test all(parse(Float64, row[required_sweeps_idx]) ≈
              parse(Float64, row[required_steps_idx]) / parse(Float64, row[proposals_per_sweep_idx]) for row in diag_rows)
    @test all(parse(Float64, row[shortfall_idx]) >= 1.0 for row in diag_rows)
    @test all(parse(Float64, row[budget_pass_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(parse(Float64, row[diagnostic_pass_idx]) in (0.0, 1.0) for row in diag_rows)
    @test all(parse(Float64, row[structure_idx]) >= 0.0 for row in diag_rows)
    @test all(parse(Float64, row[diag_alpha_idx]) == result.config.alpha for row in diag_rows)
    @test all(parse(Float64, row[diag_target_model_idx]) == 1.0 for row in diag_rows)
    @test all(isnan(parse(Float64, row[global_acceptance_idx])) for row in diag_rows)
    @test all(parse(Float64, row[global_amplitude_idx]) > 0.0 for row in diag_rows)
    @test all(parse(Float64, row[global_idx]) == 0.0 for row in diag_rows)
    @test all(parse(Float64, row[adaptive_rounds_idx]) == 0.0 for row in diag_rows)
    @test all(parse(Float64, row[adaptive_steps_idx]) == 0.0 for row in diag_rows)
    @test all(isnan(parse(Float64, row[adaptive_acceptance_idx])) for row in diag_rows)
    @test all(isnan(parse(Float64, row[adaptive_pair_acceptance_idx])) for row in diag_rows)
    @test all(isnan(parse(Float64, row[adaptive_collective_acceptance_idx])) for row in diag_rows)
    @test all(isnan(parse(Float64, row[adaptive_fourier_acceptance_idx])) for row in diag_rows)

    vk1_outdir = mktempdir()
    vk1_result = run_pressure_sweep(vk1_outdir; mode=:quick, seed=10, steps=3,
        max_steps=3, C_values=[20.0], ess_target=0.0, chains=1,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=true, tune_proposal=false, target_model=:hevk1)
    vk1_metadata = TOML.parsefile(joinpath(vk1_outdir, "run_metadata.toml"))
    @test vk1_metadata["run"]["target_model"] == "HE-VK1"
    vk1_pressure_header, vk1_pressure_rows = DFMMonteCarlo._read_csv_numeric(vk1_result.pressure_csv)
    vk1_diag_header, vk1_diag_rows = DFMMonteCarlo._read_csv_numeric(vk1_result.diagnostics_csv)
    vk1_pressure_target_idx = findfirst(==("target_model_code"), vk1_pressure_header)
    vk1_diag_target_idx = findfirst(==("target_model_code"), vk1_diag_header)
    vk1_stability_idx = findfirst(==("nonlinear_stability_pass"), vk1_diag_header)
    vk1_initial_energy_idx = findfirst(==("initial_energy"), vk1_diag_header)
    vk1_energy_drift_idx = findfirst(==("energy_drift_zscore"), vk1_diag_header)
    vk1_variance_drift_idx = findfirst(==("density_variance_drift_zscore"), vk1_diag_header)
    vk1_structure_drift_idx = findfirst(==("low_mode_structure_drift_zscore"), vk1_diag_header)
    vk1_observable_drift_idx = findfirst(==("max_observable_drift_zscore"), vk1_diag_header)
    vk1_observable_pass_idx = findfirst(==("observable_stationarity_pass"), vk1_diag_header)
    vk1_pair_acc_idx = findfirst(==("pair_acceptance_rate"), vk1_diag_header)
    vk1_collective_acc_idx = findfirst(==("collective_acceptance_rate"), vk1_diag_header)
    vk1_fourier_acc_idx = findfirst(==("fourier_acceptance_rate"), vk1_diag_header)
    vk1_pair_amp_idx = findfirst(==("pair_proposal_amplitude"), vk1_diag_header)
    vk1_fourier_tune_idx = findfirst(==("fourier_tuning_acceptance_rate"), vk1_diag_header)
    vk1_adaptive_rounds_idx = findfirst(==("adaptive_warmup_rounds"), vk1_diag_header)
    vk1_adaptive_steps_idx = findfirst(==("adaptive_warmup_steps"), vk1_diag_header)
    vk1_adaptive_acceptance_idx = findfirst(==("adaptive_warmup_acceptance_rate"), vk1_diag_header)
    vk1_diagnostic_idx = findfirst(==("diagnostic_pass"), vk1_diag_header)
    @test all(row[vk1_pressure_target_idx] == 2.0 for row in vk1_pressure_rows)
    @test all(row[vk1_diag_target_idx] == 2.0 for row in vk1_diag_rows)
    @test all(row[vk1_stability_idx] == 1.0 for row in vk1_diag_rows)
    @test all(row[vk1_initial_energy_idx] ≈ hevk1_energy(ones(4, 4, 4), Grid3D(4, 6.4),
        HEURPParams(A=300.0, B=1.0, C=row[1], alpha=0.1)) for row in vk1_diag_rows)
    @test vk1_energy_drift_idx !== nothing
    @test vk1_variance_drift_idx !== nothing
    @test vk1_structure_drift_idx !== nothing
    @test vk1_observable_drift_idx !== nothing
    @test vk1_observable_pass_idx !== nothing
    @test vk1_pair_acc_idx !== nothing
    @test vk1_collective_acc_idx !== nothing
    @test vk1_fourier_acc_idx !== nothing
    @test vk1_pair_amp_idx !== nothing
    @test vk1_fourier_tune_idx !== nothing
    @test vk1_adaptive_rounds_idx !== nothing
    @test vk1_adaptive_steps_idx !== nothing
    @test vk1_adaptive_acceptance_idx !== nothing
    @test all(row[vk1_energy_drift_idx] >= 0.0 || isnan(row[vk1_energy_drift_idx]) for row in vk1_diag_rows)
    @test all(row[vk1_variance_drift_idx] >= 0.0 || isnan(row[vk1_variance_drift_idx]) for row in vk1_diag_rows)
    @test all(row[vk1_structure_drift_idx] >= 0.0 || isnan(row[vk1_structure_drift_idx]) for row in vk1_diag_rows)
    @test all(row[vk1_observable_drift_idx] >= 0.0 || isnan(row[vk1_observable_drift_idx]) for row in vk1_diag_rows)
    @test all(0.0 <= row[vk1_pair_acc_idx] <= 1.0 || isnan(row[vk1_pair_acc_idx]) for row in vk1_diag_rows)
    @test all(0.0 <= row[vk1_collective_acc_idx] <= 1.0 || isnan(row[vk1_collective_acc_idx]) for row in vk1_diag_rows)
    @test all(0.0 <= row[vk1_fourier_acc_idx] <= 1.0 || isnan(row[vk1_fourier_acc_idx]) for row in vk1_diag_rows)
    @test all(row[vk1_pair_amp_idx] > 0.0 for row in vk1_diag_rows)
    @test all(0.0 <= row[vk1_fourier_tune_idx] <= 1.0 || isnan(row[vk1_fourier_tune_idx]) for row in vk1_diag_rows)
    @test all(row[vk1_adaptive_rounds_idx] == 0.0 for row in vk1_diag_rows)
    @test all(row[vk1_adaptive_steps_idx] == 0.0 for row in vk1_diag_rows)
    @test all(isnan(row[vk1_adaptive_acceptance_idx]) for row in vk1_diag_rows)
    @test all(row[vk1_observable_pass_idx] in (0.0, 1.0) for row in vk1_diag_rows)
    @test all(row[vk1_diagnostic_idx] <= row[vk1_observable_pass_idx] for row in vk1_diag_rows)
    vk1_text = read(write_manuscript_text(vk1_result, vk1_outdir), String)
    @test occursin("HE-VK1", vk1_text)
    @test occursin("Eq. HE-VK1", vk1_text)

    adaptive_outdir = mktempdir()
    adaptive_result = run_pressure_sweep(adaptive_outdir; mode=:quick, seed=12, steps=4,
        max_steps=4, C_values=[20.0], ess_target=0.0, chains=1,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=false, tune_proposal=false, target_model=:hevk1,
        adaptive_warmup=true, adaptive_warmup_rounds=1, adaptive_warmup_steps=2,
        adaptive_warmup_gain=0.5)
    adaptive_metadata = TOML.parsefile(joinpath(adaptive_outdir, "run_metadata.toml"))
    @test adaptive_metadata["sampling"]["adaptive_warmup"] == true
    @test adaptive_metadata["sampling"]["adaptive_warmup_rounds"] == 1
    @test adaptive_metadata["sampling"]["adaptive_warmup_steps"] == 2
    @test adaptive_metadata["sampling"]["adaptive_warmup_gain"] == 0.5
    adaptive_diag_header, adaptive_diag_rows = DFMMonteCarlo._read_csv_numeric(adaptive_result.diagnostics_csv)
    adaptive_rounds_idx = findfirst(==("adaptive_warmup_rounds"), adaptive_diag_header)
    adaptive_steps_idx = findfirst(==("adaptive_warmup_steps"), adaptive_diag_header)
    adaptive_acceptance_idx = findfirst(==("adaptive_warmup_acceptance_rate"), adaptive_diag_header)
    adaptive_fourier_acceptance_idx = findfirst(==("adaptive_warmup_fourier_acceptance_rate"), adaptive_diag_header)
    @test all(row[adaptive_rounds_idx] == 1.0 for row in adaptive_diag_rows)
    @test all(row[adaptive_steps_idx] == 2.0 for row in adaptive_diag_rows)
    @test all(0.0 <= row[adaptive_acceptance_idx] <= 1.0 for row in adaptive_diag_rows)
    @test all(0.0 <= row[adaptive_fourier_acceptance_idx] <= 1.0 ||
              isnan(row[adaptive_fourier_acceptance_idx]) for row in adaptive_diag_rows)
    adaptive_chain_header, adaptive_chain_rows = DFMMonteCarlo._read_csv_numeric(
        joinpath(adaptive_outdir, "pressure_chain_diagnostics.csv"))
    adaptive_chain_rounds_idx = findfirst(==("adaptive_warmup_rounds"), adaptive_chain_header)
    adaptive_chain_steps_idx = findfirst(==("adaptive_warmup_steps"), adaptive_chain_header)
    adaptive_chain_acceptance_idx = findfirst(==("adaptive_warmup_acceptance_rate"), adaptive_chain_header)
    @test adaptive_chain_rounds_idx !== nothing
    @test adaptive_chain_steps_idx !== nothing
    @test adaptive_chain_acceptance_idx !== nothing
    @test all(row[adaptive_chain_rounds_idx] == 1.0 for row in adaptive_chain_rows)
    @test all(row[adaptive_chain_steps_idx] == 2.0 for row in adaptive_chain_rows)
    @test all(0.0 <= row[adaptive_chain_acceptance_idx] <= 1.0 for row in adaptive_chain_rows)
    adaptive_budget_plan = read(joinpath(adaptive_outdir, "direct_mc_budget_plan.md"), String)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP=true", adaptive_budget_plan)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP_ROUNDS=1", adaptive_budget_plan)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP_STEPS=2", adaptive_budget_plan)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP_GAIN=0.5", adaptive_budget_plan)
    adaptive_text = read(write_manuscript_text(adaptive_result, adaptive_outdir), String)
    @test occursin("adaptive warmup discarded", adaptive_text)
    @test occursin("proposal amplitudes were frozen", adaptive_text)

    global_move_outdir = mktempdir()
    global_move_result = run_pressure_sweep(global_move_outdir; mode=:quick, seed=14, steps=20,
        max_steps=20, C_values=[20.0], ess_target=0.0, chains=1,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=false, tune_proposal=false, target_model=:hevk1,
        global_move_weight=1.0, global_amplitude_scale=0.1)
    global_metadata = TOML.parsefile(joinpath(global_move_outdir, "run_metadata.toml"))
    @test global_metadata["sampling"]["global_move_weight"] == 1.0
    @test global_metadata["sampling"]["global_amplitude_scale"] == 0.1
    global_diag_header, global_diag_rows = DFMMonteCarlo._read_csv_numeric(global_move_result.diagnostics_csv)
    global_diag_acceptance_idx = findfirst(==("global_acceptance_rate"), global_diag_header)
    global_diag_amplitude_idx = findfirst(==("global_field_amplitude"), global_diag_header)
    global_diag_weight_idx = findfirst(==("global_move_weight"), global_diag_header)
    @test all(0.0 <= row[global_diag_acceptance_idx] <= 1.0 for row in global_diag_rows)
    @test all(row[global_diag_amplitude_idx] > 0.0 for row in global_diag_rows)
    @test all(row[global_diag_weight_idx] > 0.0 for row in global_diag_rows)

    global_warmup_outdir = mktempdir()
    global_warmup_result = run_pressure_sweep(global_warmup_outdir; mode=:quick, seed=15, steps=20,
        max_steps=20, C_values=[20.0], ess_target=0.0, chains=1,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=false, tune_proposal=false, target_model=:hevk1,
        adaptive_warmup=true, adaptive_warmup_rounds=1, adaptive_warmup_steps=2,
        global_move_weight=0.25, global_amplitude_scale=0.1)
    global_warmup_diag_header, global_warmup_diag_rows =
        DFMMonteCarlo._read_csv_numeric(global_warmup_result.diagnostics_csv)
    global_warmup_acceptance_idx = findfirst(==("global_acceptance_rate"), global_warmup_diag_header)
    global_warmup_weight_idx = findfirst(==("global_move_weight"), global_warmup_diag_header)
    @test all(0.0 <= row[global_warmup_acceptance_idx] <= 1.0 for row in global_warmup_diag_rows)
    @test all(row[global_warmup_weight_idx] > 0.0 for row in global_warmup_diag_rows)

    radial_move_outdir = mktempdir()
    radial_move_result = run_pressure_sweep(radial_move_outdir; mode=:quick, seed=16, steps=24,
        max_steps=24, C_values=[20.0], ess_target=0.0, chains=1,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=false, tune_proposal=false, target_model=:hevk1,
        initial_density_mode=:overdispersed, initial_density_amplitude=0.03,
        radial_variance_move_weight=0.5, radial_variance_amplitude_scale=0.2)
    radial_metadata = TOML.parsefile(joinpath(radial_move_outdir, "run_metadata.toml"))
    @test radial_metadata["sampling"]["radial_variance_move_weight"] == 0.5
    @test radial_metadata["sampling"]["radial_variance_amplitude_scale"] == 0.2
    radial_diag_header, radial_diag_rows = DFMMonteCarlo._read_csv_numeric(radial_move_result.diagnostics_csv)
    radial_acceptance_idx = findfirst(==("radial_variance_acceptance_rate"), radial_diag_header)
    radial_amplitude_idx = findfirst(==("radial_variance_amplitude"), radial_diag_header)
    radial_weight_idx = findfirst(==("radial_variance_move_weight"), radial_diag_header)
    @test radial_acceptance_idx !== nothing
    @test radial_amplitude_idx !== nothing
    @test radial_weight_idx !== nothing
    @test all(0.0 <= row[radial_acceptance_idx] <= 1.0 || isnan(row[radial_acceptance_idx]) for row in radial_diag_rows)
    @test all(row[radial_amplitude_idx] > 0.0 for row in radial_diag_rows)
    @test all(row[radial_weight_idx] > 0.0 for row in radial_diag_rows)

    rx_outdir = mktempdir()
    rx_result = run_direct_replica_exchange(rx_outdir; seed=21, steps=24,
        C=20.0, Nx=4, betas=[1.0, 0.5], target_model=:hevk1,
        sample_stride=2, swap_interval=4, proposal_amplitude=0.01,
        global_move_weight=0.25)
    @test isfile(rx_result.samples_csv)
    @test isfile(rx_result.diagnostics_csv)
    @test isfile(rx_result.markdown_path)
    rx_header, rx_rows = DFMMonteCarlo._read_csv_numeric(rx_result.samples_csv)
    rx_diag_header, rx_diag_rows = DFMMonteCarlo._read_csv_numeric(rx_result.diagnostics_csv)
    @test "pressure" in rx_header
    @test "density_variance" in rx_header
    @test "replica_beta" in rx_header
    @test "swap_acceptance_rate" in rx_diag_header
    @test "target_pressure_mean" in rx_diag_header
    @test "density_variance_drift_zscore" in rx_diag_header
    @test "low_mode_structure_drift_zscore" in rx_diag_header
    @test "observable_stationarity_pass" in rx_diag_header
    @test "replica_round_trips" in rx_diag_header
    @test "max_replica_round_trips" in rx_diag_header
    @test "min_replica_round_trips" in rx_diag_header
    @test !isempty(rx_rows)
    target_sample_idx = findfirst(==("target_sample_count"), rx_diag_header)
    swap_idx = findfirst(==("swap_acceptance_rate"), rx_diag_header)
    variance_drift_idx = findfirst(==("density_variance_drift_zscore"), rx_diag_header)
    structure_drift_idx = findfirst(==("low_mode_structure_drift_zscore"), rx_diag_header)
    observable_pass_idx = findfirst(==("observable_stationarity_pass"), rx_diag_header)
    round_trips_idx = findfirst(==("replica_round_trips"), rx_diag_header)
    max_round_trips_idx = findfirst(==("max_replica_round_trips"), rx_diag_header)
    min_round_trips_idx = findfirst(==("min_replica_round_trips"), rx_diag_header)
    diag_row = only(rx_diag_rows)
    @test diag_row[target_sample_idx] >= 1.0
    @test 0.0 <= diag_row[swap_idx] <= 1.0
    @test diag_row[variance_drift_idx] >= 0.0 || isnan(diag_row[variance_drift_idx])
    @test diag_row[structure_drift_idx] >= 0.0 || isnan(diag_row[structure_drift_idx])
    @test diag_row[observable_pass_idx] in (0.0, 1.0)
    @test diag_row[round_trips_idx] >= 0.0
    @test diag_row[max_round_trips_idx] >= diag_row[min_round_trips_idx]

    mixed_outdir = mktempdir()
    mixed_result = run_pressure_sweep(mixed_outdir; mode=:quick, seed=11, steps=3,
        max_steps=3, C_values=[20.0], ess_target=0.0, chains=2,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=true, tune_proposal=false, target_model=:hevk1,
        initial_density_mode=:mixed, initial_density_amplitude=0.02)
    mixed_metadata = TOML.parsefile(joinpath(mixed_outdir, "run_metadata.toml"))
    @test mixed_metadata["sampling"]["initial_density_mode"] == "mixed"
    mixed_diag_header, mixed_diag_rows = DFMMonteCarlo._read_csv_numeric(mixed_result.diagnostics_csv)
    mixed_over_idx = findfirst(==("initial_density_overdispersed"), mixed_diag_header)
    mixed_var_idx = findfirst(==("initial_density_variance"), mixed_diag_header)
    @test all(row[mixed_over_idx] == 1.0 for row in mixed_diag_rows)
    @test all(row[mixed_var_idx] > 0.0 for row in mixed_diag_rows)
    mixed_chain_header, mixed_chain_rows = DFMMonteCarlo._read_csv_numeric(
        joinpath(mixed_outdir, "pressure_chain_diagnostics.csv"))
    mixed_chain_idx = findfirst(==("chain"), mixed_chain_header)
    mixed_chain_over_idx = findfirst(==("initial_density_overdispersed"), mixed_chain_header)
    mixed_chain_var_idx = findfirst(==("initial_density_variance"), mixed_chain_header)
    mixed_chain_adaptive_rounds_idx = findfirst(==("adaptive_warmup_rounds"), mixed_chain_header)
    mixed_chain_adaptive_steps_idx = findfirst(==("adaptive_warmup_steps"), mixed_chain_header)
    mixed_chain_adaptive_acceptance_idx = findfirst(==("adaptive_warmup_acceptance_rate"), mixed_chain_header)
    @test Set(row[mixed_chain_idx] for row in mixed_chain_rows) == Set([1.0, 2.0])
    @test Set(row[mixed_chain_over_idx] for row in mixed_chain_rows) == Set([0.0, 1.0])
    @test only(row[mixed_chain_var_idx] for row in mixed_chain_rows if row[mixed_chain_over_idx] == 0.0) == 0.0
    @test only(row[mixed_chain_var_idx] for row in mixed_chain_rows if row[mixed_chain_over_idx] == 1.0) > 0.0
    @test mixed_chain_adaptive_rounds_idx !== nothing
    @test mixed_chain_adaptive_steps_idx !== nothing
    @test mixed_chain_adaptive_acceptance_idx !== nothing
    @test all(row[mixed_chain_adaptive_rounds_idx] == 0.0 for row in mixed_chain_rows)
    @test all(row[mixed_chain_adaptive_steps_idx] == 0.0 for row in mixed_chain_rows)
    @test all(isnan(row[mixed_chain_adaptive_acceptance_idx]) for row in mixed_chain_rows)

    mixed_budget_plan_path = joinpath(mixed_outdir, "direct_mc_budget_plan.md")
    @test isfile(mixed_budget_plan_path)
    mixed_budget_plan = read(mixed_budget_plan_path, String)
    @test occursin("Recommended direct MC rerun", mixed_budget_plan)
    @test occursin("estimated_required_steps", mixed_budget_plan)
    @test occursin("DFM_MC_QUICK_TARGET_MODEL=hevk1", mixed_budget_plan)
    @test occursin("DFM_MC_QUICK_C_VALUES=20", mixed_budget_plan)
    @test occursin("DFM_MC_QUICK_INITIAL_DENSITY_MODE=mixed", mixed_budget_plan)
    @test occursin("DFM_MC_QUICK_STEPS", mixed_budget_plan)

    sweep_script = read(joinpath(@__DIR__, "..", "scripts", "run_pressure_sweep.jl"), String)
    @test occursin("DFM_MC_QUICK_TARGET_MODEL", sweep_script)
    @test occursin("DFM_MC_QUICK_C_VALUES", sweep_script)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP", sweep_script)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP_ROUNDS", sweep_script)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP_STEPS", sweep_script)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP_GAIN", sweep_script)
    @test occursin("DFM_MC_QUICK_GLOBAL_MOVE_WEIGHT", sweep_script)
    @test occursin("DFM_MC_QUICK_GLOBAL_AMPLITUDE_SCALE", sweep_script)
    @test occursin("C_values=c_values", sweep_script)
    @test occursin("target_model=target_model", sweep_script)
    @test occursin("adaptive_warmup=adaptive_warmup", sweep_script)
    @test occursin("adaptive_warmup_rounds=adaptive_warmup_rounds", sweep_script)
    @test occursin("adaptive_warmup_steps=adaptive_warmup_steps", sweep_script)
    @test occursin("adaptive_warmup_gain=adaptive_warmup_gain", sweep_script)
    @test occursin("global_move_weight=global_move_weight", sweep_script)
    @test occursin("global_amplitude_scale=global_amplitude_scale", sweep_script)
    @test occursin("write_manuscript_text(result, outdir)", sweep_script)
    @test occursin("write_pressure_svg(result.pressure_csv, outdir)", sweep_script)
    @test occursin("write_diagnostics_svg(result.diagnostics_csv, outdir)", sweep_script)

    trace_csv = joinpath(outdir, "pressure_traces.csv")
    @test isfile(trace_csv)
    trace_header, trace_rows = DFMMonteCarlo._read_csv_numeric(trace_csv)
    @test trace_header == ["C", "chain", "sample_index", "pressure"]
    trace_C_idx = findfirst(==("C"), trace_header)
    trace_chain_idx = findfirst(==("chain"), trace_header)
    trace_pressure_idx = findfirst(==("pressure"), trace_header)
    @test Set(row[trace_C_idx] for row in trace_rows) == Set(result.config.C_values)
    @test Set(row[trace_chain_idx] for row in trace_rows) == Set([1.0, 2.0])
    @test all(isfinite(row[trace_pressure_idx]) for row in trace_rows)
    @test length(trace_rows) == sum(round(Int, parse(Float64, row[trace_samples_idx])) for row in diag_rows)
    for diag_row in diag_rows
        C_value = parse(Float64, diag_row[1])
        expected_count = round(Int, parse(Float64, diag_row[trace_samples_idx]))
        @test count(row -> row[trace_C_idx] == C_value, trace_rows) == expected_count
    end
    observable_trace_csv = joinpath(outdir, "observable_traces.csv")
    @test isfile(observable_trace_csv)
    observable_trace_header, observable_trace_rows = DFMMonteCarlo._read_csv_numeric(observable_trace_csv)
    @test observable_trace_header == [
        "C",
        "chain",
        "sample_index",
        "pressure",
        "energy",
        "density_variance",
        "low_mode_structure_factor",
    ]
    observable_trace_C_idx = findfirst(==("C"), observable_trace_header)
    observable_trace_chain_idx = findfirst(==("chain"), observable_trace_header)
    observable_trace_pressure_idx = findfirst(==("pressure"), observable_trace_header)
    observable_trace_energy_idx = findfirst(==("energy"), observable_trace_header)
    observable_trace_variance_idx = findfirst(==("density_variance"), observable_trace_header)
    observable_trace_structure_idx = findfirst(==("low_mode_structure_factor"), observable_trace_header)
    @test length(observable_trace_rows) == length(trace_rows)
    @test Set(row[observable_trace_C_idx] for row in observable_trace_rows) == Set(result.config.C_values)
    @test Set(row[observable_trace_chain_idx] for row in observable_trace_rows) == Set([1.0, 2.0])
    @test all(isfinite(row[observable_trace_pressure_idx]) for row in observable_trace_rows)
    @test all(isfinite(row[observable_trace_energy_idx]) for row in observable_trace_rows)
    @test all(row[observable_trace_variance_idx] >= 0.0 for row in observable_trace_rows)
    @test all(row[observable_trace_structure_idx] >= 0.0 for row in observable_trace_rows)
    trace_summary = write_pressure_trace_summary(outdir)
    @test isfile(trace_summary.csv_path)
    @test isfile(trace_summary.markdown_path)
    trace_summary_header, trace_summary_rows =
        DFMMonteCarlo._read_csv_numeric(trace_summary.csv_path)
    @test trace_summary_header == [
        "C",
        "chain_count",
        "trace_samples",
        "min_chain_samples",
        "max_chain_samples",
        "pressure_mean",
        "pressure_stderr",
        "pressure_stderr_independent",
        "integrated_autocorrelation_time",
        "effective_sample_size",
        "trace_drift_zscore",
        "chain_rhat",
        "split_chain_rhat",
        "chain_between_stderr",
        "chain_within_stderr",
    ]
    @test length(trace_summary_rows) == length(result.config.C_values)
    summary_C_idx = findfirst(==("C"), trace_summary_header)
    summary_chain_idx = findfirst(==("chain_count"), trace_summary_header)
    summary_samples_idx = findfirst(==("trace_samples"), trace_summary_header)
    summary_mean_idx = findfirst(==("pressure_mean"), trace_summary_header)
    summary_error_idx = findfirst(==("pressure_stderr"), trace_summary_header)
    summary_ess_idx = findfirst(==("effective_sample_size"), trace_summary_header)
    summary_rhat_idx = findfirst(==("chain_rhat"), trace_summary_header)
    summary_split_idx = findfirst(==("split_chain_rhat"), trace_summary_header)
    @test Set(row[summary_C_idx] for row in trace_summary_rows) == Set(result.config.C_values)
    @test all(row[summary_chain_idx] == 2.0 for row in trace_summary_rows)
    @test all(row[summary_samples_idx] ==
              count(trace_row -> trace_row[trace_C_idx] == row[summary_C_idx], trace_rows)
              for row in trace_summary_rows)
    @test all(isfinite(row[summary_mean_idx]) for row in trace_summary_rows)
    @test all(row[summary_error_idx] >= 0.0 for row in trace_summary_rows)
    @test all(row[summary_ess_idx] >= 1.0 for row in trace_summary_rows)
    @test all(row[summary_rhat_idx] >= 1.0 || isnan(row[summary_rhat_idx]) for row in trace_summary_rows)
    @test all(row[summary_split_idx] >= 1.0 || isnan(row[summary_split_idx]) for row in trace_summary_rows)
    trace_summary_text = read(trace_summary.markdown_path, String)
    @test occursin("Trace-only Arianna pressure summary", trace_summary_text)
    @test occursin("pressure_traces.csv", trace_summary_text)
    @test occursin("interrupted", lowercase(trace_summary_text))

    resume_outdir = mktempdir()
    first_resume = run_pressure_sweep(resume_outdir; mode=:quick, seed=113, steps=3,
        max_steps=3, C_values=[20.0], ess_target=0.0, chains=1,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=false, tune_proposal=false)
    resume_trace_csv = joinpath(resume_outdir, "pressure_traces.csv")
    resume_chain_csv = joinpath(resume_outdir, "pressure_chain_diagnostics.csv")
    @test isfile(resume_trace_csv)
    @test isfile(resume_chain_csv)
    resume_header_1, resume_rows_1 = DFMMonteCarlo._read_csv_numeric(resume_trace_csv)
    resume_chain_idx = findfirst(==("chain"), resume_header_1)
    @test Set(row[resume_chain_idx] for row in resume_rows_1) == Set([1.0])
    resumed = run_pressure_sweep(resume_outdir; mode=:quick, seed=113, steps=3,
        max_steps=3, C_values=[20.0], ess_target=0.0, chains=2,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=false, tune_proposal=false, resume_traces=true)
    @test resumed.pressure_csv == first_resume.pressure_csv
    resume_header_2, resume_rows_2 = DFMMonteCarlo._read_csv_numeric(resume_trace_csv)
    @test resume_header_2 == resume_header_1
    @test Set(row[resume_chain_idx] for row in resume_rows_2) == Set([1.0, 2.0])
    @test count(row -> row[resume_chain_idx] == 1.0, resume_rows_2) ==
          count(row -> row[resume_chain_idx] == 1.0, resume_rows_1)
    chain_header, chain_rows = DFMMonteCarlo._read_csv_numeric(resume_chain_csv)
    @test "acceptance_rate" in chain_header
    @test "final_energy" in chain_header
    @test "equilibration_sensitivity_pass" in chain_header
    chain_diag_idx = findfirst(==("chain"), chain_header)
    @test Set(row[chain_diag_idx] for row in chain_rows) == Set([1.0, 2.0])
    resumed_diag_header, resumed_diag_rows = DFMMonteCarlo._read_csv_numeric(resumed.diagnostics_csv)
    resumed_chain_count_idx = findfirst(==("chain_count"), resumed_diag_header)
    resumed_rhat_idx = findfirst(==("chain_rhat"), resumed_diag_header)
    @test all(row[resumed_chain_count_idx] == 2.0 for row in resumed_diag_rows)
    @test all(row[resumed_rhat_idx] >= 1.0 || isnan(row[resumed_rhat_idx]) for row in resumed_diag_rows)
    @test all(parse(Float64, row[diagnostic_pass_idx]) ==
              min(parse(Float64, row[acceptance_pass_idx]),
                  parse(Float64, row[drift_pass_idx]),
                  parse(Float64, row[rhat_pass_idx]),
                  parse(Float64, row[ess_pass_idx]),
                  parse(Float64, row[stability_pass_idx]),
                  parse(Float64, row[trace_samples_pass_idx]),
                  parse(Float64, row[pressure_error_pass_idx]),
                  parse(Float64, row[budget_pass_idx]),
                  parse(Float64, row[decorrelation_pass_idx]),
                  parse(Float64, row[equilibration_sensitivity_pass_idx])) for row in diag_rows)
    if amplitude_idx !== nothing && target_idx !== nothing && tuning_idx !== nothing
        @test all(parse(Float64, row[amplitude_idx]) > 0.0 for row in diag_rows)
        @test all(0.0 < parse(Float64, row[target_idx]) < 1.0 for row in diag_rows)
        @test all(parse(Float64, row[tuning_idx]) >= 8.0 for row in diag_rows)
        @test all(parse(Float64, row[tuning_steps_idx]) >= 25.0 for row in diag_rows)
    end
    if collective_amplitude_idx !== nothing && collective_pair_count_idx !== nothing && collective_idx !== nothing
        @test all(parse(Float64, row[collective_amplitude_idx]) > 0.0 for row in diag_rows)
        @test all(parse(Float64, row[collective_pair_count_idx]) >= 2.0 for row in diag_rows)
        @test all(0.0 < parse(Float64, row[collective_idx]) < 1.0 for row in diag_rows)
    end
    if fourier_amplitude_idx !== nothing && fourier_max_mode_idx !== nothing && fourier_idx !== nothing
        @test all(parse(Float64, row[fourier_amplitude_idx]) > 0.0 for row in diag_rows)
        @test all(parse(Float64, row[fourier_max_mode_idx]) >= 1.0 for row in diag_rows)
        @test all(0.0 < parse(Float64, row[fourier_idx]) < 1.0 for row in diag_rows)
    end

    basin_outdir = mktempdir()
    basin_cap = 1.0e-4
    basin_result = run_pressure_sweep(basin_outdir; mode=:quick, seed=11, steps=6,
        max_steps=6, C_values=[20.0], ess_target=0.0, chains=1,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=false, basin_variance_max=basin_cap)
    basin_metadata = TOML.parsefile(joinpath(basin_outdir, "run_metadata.toml"))
    @test basin_metadata["sampling"]["basin_variance_max"] == basin_cap
    basin_header, _ = DFMMonteCarlo._read_csv_numeric(basin_result.pressure_csv)
    @test "basin_variance_max" in basin_header
    @test "basin_constraint_pass" in basin_header
    @test "phi_min" in basin_header
    basin_diag_header, basin_diag_rows = DFMMonteCarlo._read_csv_strings(basin_result.diagnostics_csv)
    @test "basin_variance_max" in basin_diag_header
    @test "basin_constraint_pass" in basin_diag_header
    @test "phi_min" in basin_diag_header
    basin_cap_idx = findfirst(==("basin_variance_max"), basin_diag_header)
    basin_pass_idx = findfirst(==("basin_constraint_pass"), basin_diag_header)
    basin_density_idx = findfirst(==("density_variance"), basin_diag_header)
    basin_phi_min_idx = findfirst(==("phi_min"), basin_diag_header)
    @test all(parse(Float64, row[basin_cap_idx]) == basin_cap for row in basin_diag_rows)
    @test all(parse(Float64, row[basin_pass_idx]) == 1.0 for row in basin_diag_rows)
    @test all(parse(Float64, row[basin_density_idx]) <= basin_cap + 1.0e-12 for row in basin_diag_rows)
    @test all(parse(Float64, row[basin_phi_min_idx]) == 1.0e-8 for row in basin_diag_rows)
    basin_text_path = write_manuscript_text(basin_result, basin_outdir)
    basin_text = read(basin_text_path, String)
    @test occursin("hard-wall Gaussian-basin restriction", basin_text)
    @test occursin("constrained-target diagnostic", basin_text)
    @test !occursin("manuscript-supporting validation", basin_text)
    basin_audit = write_validation_audit(mktempdir(); direct_dirs=[basin_outdir])
    basin_audit_metadata = TOML.parsefile(basin_audit.toml_path)
    @test basin_audit_metadata["verdict"]["direct_status"] == "diagnostic"
    @test basin_audit_metadata["direct"]["finite_basin_points"] == 1
    @test occursin("basin-restricted direct outputs remain diagnostic", read(basin_audit.markdown_path, String))

    barrier_outdir = mktempdir()
    barrier_strength = 20.0
    barrier_result = run_pressure_sweep(barrier_outdir; mode=:quick, seed=12, steps=4,
        max_steps=4, C_values=[20.0], ess_target=0.0, chains=1,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=false, boundary_barrier_strength=barrier_strength)
    barrier_metadata = TOML.parsefile(joinpath(barrier_outdir, "run_metadata.toml"))
    @test barrier_metadata["sampling"]["boundary_barrier_strength"] == barrier_strength
    barrier_header, barrier_rows = DFMMonteCarlo._read_csv_numeric(barrier_result.pressure_csv)
    barrier_diag_header, barrier_diag_rows = DFMMonteCarlo._read_csv_numeric(barrier_result.diagnostics_csv)
    @test "boundary_barrier_strength" in barrier_header
    @test "boundary_barrier_strength" in barrier_diag_header
    barrier_pressure_idx = findfirst(==("boundary_barrier_strength"), barrier_header)
    barrier_diag_idx = findfirst(==("boundary_barrier_strength"), barrier_diag_header)
    @test all(row[barrier_pressure_idx] == barrier_strength for row in barrier_rows)
    @test all(row[barrier_diag_idx] == barrier_strength for row in barrier_diag_rows)
    barrier_text = read(write_manuscript_text(barrier_result, barrier_outdir), String)
    @test occursin("logarithmic boundary-barrier regularization", barrier_text)
    @test occursin("regularized direct-target diagnostic", barrier_text)
    @test !occursin("manuscript-supporting validation", barrier_text)
    barrier_audit = write_validation_audit(mktempdir(); direct_dirs=[barrier_outdir])
    barrier_audit_metadata = TOML.parsefile(barrier_audit.toml_path)
    @test barrier_audit_metadata["verdict"]["direct_status"] == "diagnostic"
    @test barrier_audit_metadata["direct"]["boundary_barrier_points"] == 1
    @test occursin("boundary-barrier direct outputs remain diagnostic",
        read(barrier_audit.markdown_path, String))

    lower_bound_outdir = mktempdir()
    lower_bound_result = run_pressure_sweep(lower_bound_outdir; mode=:quick, seed=13, steps=4,
        max_steps=4, C_values=[20.0], ess_target=0.0, chains=1,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=4, stability_probe=true, phi_min=0.05)
    lower_bound_metadata = TOML.parsefile(joinpath(lower_bound_outdir, "run_metadata.toml"))
    @test lower_bound_metadata["sampling"]["phi_min"] == 0.05
    lower_pressure_header, lower_pressure_rows = DFMMonteCarlo._read_csv_numeric(lower_bound_result.pressure_csv)
    lower_diag_header, lower_diag_rows = DFMMonteCarlo._read_csv_numeric(lower_bound_result.diagnostics_csv)
    lower_phi_pressure_idx = findfirst(==("phi_min"), lower_pressure_header)
    lower_phi_diag_idx = findfirst(==("phi_min"), lower_diag_header)
    lower_stability_phi_idx = findfirst(==("nonlinear_stability_min_phi"), lower_diag_header)
    @test lower_phi_pressure_idx !== nothing
    @test lower_phi_diag_idx !== nothing
    @test all(row[lower_phi_pressure_idx] == 0.05 for row in lower_pressure_rows)
    @test all(row[lower_phi_diag_idx] == 0.05 for row in lower_diag_rows)
    @test all(row[lower_stability_phi_idx] >= 0.05 - 1.0e-12 for row in lower_diag_rows)
    lower_bound_audit = write_validation_audit(mktempdir(); direct_dirs=[lower_bound_outdir])
    lower_bound_audit_metadata = TOML.parsefile(lower_bound_audit.toml_path)
    @test lower_bound_audit_metadata["direct"]["max_phi_min"] == 0.05
    @test occursin("density lower bound", read(lower_bound_audit.markdown_path, String))
    legacy_direct_outdir = mktempdir()
    open(joinpath(legacy_direct_outdir, "sampling_diagnostics.csv"), "w") do io
        println(io, "diagnostic_pass,nonlinear_stability_pass,sampling_budget_pass,nonlinear_stability_min_energy_delta,nonlinear_stability_min_phi,production_steps,production_sweeps,effective_sample_size,low_mode_structure_factor")
        println(io, "0,0,0,-1.0,0.001,20,0.001,1.0,0.0")
    end
    downhill_direct_audit = write_validation_audit(mktempdir(); direct_dirs=[legacy_direct_outdir])
    downhill_direct_metadata = TOML.parsefile(downhill_direct_audit.toml_path)
    @test downhill_direct_metadata["verdict"]["direct_recommendation"] == "reformulate_or_constrain_nonlinear_target"
    @test occursin("Recommended direct-MC next step: reformulate or constrain the nonlinear target",
        read(downhill_direct_audit.markdown_path, String))
    mixed_support_audit = write_validation_audit(mktempdir(); direct_dirs=[legacy_direct_outdir, lower_bound_outdir])
    mixed_support_metadata = TOML.parsefile(mixed_support_audit.toml_path)
    @test mixed_support_metadata["direct"]["phi_min_present_csv_count"] == 1
    @test mixed_support_metadata["direct"]["phi_min_missing_csv_count"] == 1
    @test occursin("missing explicit phi_min", read(mixed_support_audit.markdown_path, String))
    lower_bound_text = read(write_manuscript_text(lower_bound_result, lower_bound_outdir), String)
    @test occursin("density lower bound", lower_bound_text)

    @test isdefined(DFMMonteCarlo, :run_direct_floor_sensitivity_study)
    floor_script_path = joinpath(@__DIR__, "..", "scripts", "run_direct_floor_sensitivity_study.jl")
    @test isfile(floor_script_path)
    if isdefined(DFMMonteCarlo, :run_direct_floor_sensitivity_study) && isfile(floor_script_path)
        floor_outdir = mktempdir()
        floor_study = DFMMonteCarlo.run_direct_floor_sensitivity_study(floor_outdir;
            Nx=3, C=20.0, phi_min_values=[0.05, 0.1], steps=2, max_steps=2,
            chains=1, ess_target=0.0, min_trace_samples=1, pressure_tolerance=Inf,
            burn_fraction=0.0, sample_stride=1, tune_proposal=false, stability_probe=false)
        @test isfile(floor_study.csv_path)
        @test isfile(floor_study.markdown_path)
        @test isfile(floor_study.budget_plan_path)
        @test isfile(floor_study.svg_path)
        floor_svg = read(floor_study.svg_path, String)
        @test occursin("<title>Direct nonlinear floor sensitivity</title>", floor_svg)
        @test occursin("phi_min", floor_svg)
        @test occursin("Pi", floor_svg)
        @test occursin("x-major-tick", floor_svg)
        @test occursin("x-minor-tick", floor_svg)
        @test occursin("y-major-tick", floor_svg)
        @test occursin("y-minor-tick", floor_svg)
        @test occursin("MC direct", floor_svg)
        @test occursin("CL nearest", floor_svg)
        @test occursin("DFM grid", floor_svg)
        @test occursin("Mean field", floor_svg)
        floor_header, floor_rows = DFMMonteCarlo._read_csv_numeric(floor_study.csv_path)
        @test "phi_min" in floor_header
        @test "mc_pressure" in floor_header
        @test "effective_sample_size" in floor_header
        @test "pressure_span" in floor_header
        @test "required_effective_sample_size" in floor_header
        @test "required_trace_samples" in floor_header
        @test "estimated_required_steps" in floor_header
        @test "estimated_required_sweeps" in floor_header
        @test "production_budget_shortfall" in floor_header
        @test "diagnostic_pass" in floor_header
        @test length(floor_rows) == 2
        floor_phi_idx = findfirst(==("phi_min"), floor_header)
        floor_span_idx = findfirst(==("pressure_span"), floor_header)
        floor_required_steps_idx = findfirst(==("estimated_required_steps"), floor_header)
        floor_shortfall_idx = findfirst(==("production_budget_shortfall"), floor_header)
        @test sort([row[floor_phi_idx] for row in floor_rows]) == [0.05, 0.1]
        @test all(row[floor_span_idx] >= 0.0 for row in floor_rows)
        @test all(row[floor_required_steps_idx] >= 1.0 for row in floor_rows)
        @test all(row[floor_shortfall_idx] >= 1.0 for row in floor_rows)
        floor_text = read(floor_study.markdown_path, String)
        @test occursin("Direct nonlinear floor-sensitivity MC study", floor_text)
        @test occursin("lower-bound-regularized target", floor_text)
        @test occursin("estimated required", floor_text)
        @test occursin("not unconstrained HE-URP evidence", floor_text)
        floor_budget_plan = read(floor_study.budget_plan_path, String)
        @test occursin("Recommended direct floor-sensitivity rerun", floor_budget_plan)
        @test occursin("DFM_DIRECT_FLOOR_STEPS", floor_budget_plan)
        @test occursin("not unconstrained HE-URP validation", floor_budget_plan)
        floor_reused = DFMMonteCarlo.run_direct_floor_sensitivity_study(floor_outdir;
            Nx=3, C=20.0, phi_min_values=[0.05, 0.1], steps=2, max_steps=2,
            chains=1, ess_target=0.0, min_trace_samples=1, pressure_tolerance=Inf,
            burn_fraction=0.0, sample_stride=1, tune_proposal=false,
            stability_probe=false, reuse_existing=true)
        @test isfile(floor_reused.csv_path)
        @test isfile(floor_reused.budget_plan_path)
        @test isfile(floor_reused.svg_path)
        floor_script_text = read(floor_script_path, String)
        @test occursin("run_direct_floor_sensitivity_study", floor_script_text)
        @test occursin("DFM_DIRECT_FLOOR_PHI_MIN_VALUES", floor_script_text)
        floor_audit = write_validation_audit(mktempdir(); floor_sensitivity_dirs=[floor_outdir])
        floor_audit_metadata = TOML.parsefile(floor_audit.toml_path)
        @test floor_audit_metadata["direct"]["floor_sensitivity_points"] == 2
        @test floor_audit_metadata["direct"]["floor_sensitivity_phi_min_count"] == 2
        @test haskey(floor_audit_metadata["direct"], "max_floor_sensitivity_budget_shortfall")
        @test floor_audit_metadata["direct"]["max_floor_sensitivity_budget_shortfall"] >= 1.0
        floor_audit_text = read(floor_audit.markdown_path, String)
        @test occursin("Direct floor-sensitivity MC", floor_audit_text)
        @test occursin("maximum budget shortfall", floor_audit_text)
    end

    stable_subspace_outdir = mktempdir()
    stable_subspace_result = run_pressure_sweep(stable_subspace_outdir; mode=:quick, seed=17, steps=4,
        max_steps=4, C_values=[20.0], ess_target=0.0, chains=1,
        min_acceptance=0.0, max_acceptance=1.0, min_trace_samples=1,
        pressure_tolerance=Inf, burn_fraction=0.0, sample_stride=1,
        Nx=8, stability_probe=true, stable_mode_proposals=true,
        fourier_max_mode=4, fourier_move_weight=1.0,
        tuning_rounds=2, tuning_steps=5,
        initial_density_mode=:overdispersed, initial_density_amplitude=0.02)
    stable_metadata = TOML.parsefile(joinpath(stable_subspace_outdir, "run_metadata.toml"))
    @test stable_metadata["sampling"]["stable_mode_proposals"] == true
    @test stable_metadata["sampling"]["fourier_min_k2"] == 30.0
    @test stable_metadata["sampling"]["initial_density_mode"] == "overdispersed"
    @test stable_metadata["sampling"]["initial_density_amplitude"] == 0.02
    stable_diag_header, stable_diag_rows =
        DFMMonteCarlo._read_csv_numeric(stable_subspace_result.diagnostics_csv)
    stable_flag_idx = findfirst(==("stable_mode_proposals"), stable_diag_header)
    stable_min_k2_idx = findfirst(==("fourier_min_k2"), stable_diag_header)
    stable_collective_idx = findfirst(==("collective_move_weight"), stable_diag_header)
    stable_fourier_idx = findfirst(==("fourier_move_weight"), stable_diag_header)
    stable_stability_idx = findfirst(==("nonlinear_stability_pass"), stable_diag_header)
    stable_amplitude_idx = findfirst(==("proposal_amplitude"), stable_diag_header)
    stable_fourier_amplitude_idx = findfirst(==("fourier_mode_amplitude"), stable_diag_header)
    stable_tuning_rounds_idx = findfirst(==("tuning_rounds"), stable_diag_header)
    stable_tuning_steps_idx = findfirst(==("tuning_steps"), stable_diag_header)
    stable_tuning_acceptance_idx = findfirst(==("tuning_acceptance_rate"), stable_diag_header)
    stable_initial_flag_idx = findfirst(==("initial_density_overdispersed"), stable_diag_header)
    stable_initial_amp_idx = findfirst(==("initial_density_amplitude"), stable_diag_header)
    stable_initial_var_idx = findfirst(==("initial_density_variance"), stable_diag_header)
    stable_initial_energy_idx = findfirst(==("initial_energy"), stable_diag_header)
    @test stable_flag_idx !== nothing
    @test stable_min_k2_idx !== nothing
    @test stable_initial_flag_idx !== nothing
    @test stable_initial_amp_idx !== nothing
    @test stable_initial_var_idx !== nothing
    @test stable_initial_energy_idx !== nothing
    @test all(row[stable_flag_idx] == 1.0 for row in stable_diag_rows)
    @test all(row[stable_min_k2_idx] == 30.0 for row in stable_diag_rows)
    @test all(row[stable_initial_flag_idx] == 1.0 for row in stable_diag_rows)
    @test all(row[stable_initial_amp_idx] == 0.02 for row in stable_diag_rows)
    @test all(row[stable_initial_var_idx] > 0.0 for row in stable_diag_rows)
    @test all(isfinite(row[stable_initial_energy_idx]) for row in stable_diag_rows)
    @test all(row[stable_collective_idx] == 0.0 for row in stable_diag_rows)
    @test all(row[stable_fourier_idx] == 1.0 for row in stable_diag_rows)
    @test all(row[stable_stability_idx] in (0.0, 1.0) for row in stable_diag_rows)
    @test all(row[stable_tuning_rounds_idx] == 2.0 for row in stable_diag_rows)
    @test all(row[stable_tuning_steps_idx] == 5.0 for row in stable_diag_rows)
    @test all(isfinite(row[stable_tuning_acceptance_idx]) for row in stable_diag_rows)
    @test all(row[stable_amplitude_idx] == row[stable_fourier_amplitude_idx] for row in stable_diag_rows)
    stable_decorr_stride_idx = findfirst(==("decorrelation_stride_samples"), stable_diag_header)
    stable_decorr_steps_idx = findfirst(==("decorrelation_stride_steps"), stable_diag_header)
    stable_decorr_count_idx = findfirst(==("decorrelated_sample_count"), stable_diag_header)
    stable_decorr_pass_idx = findfirst(==("decorrelation_pass"), stable_diag_header)
    stable_equil_tail_idx = findfirst(==("equilibration_tail_sample_count"), stable_diag_header)
    stable_equil_z_idx = findfirst(==("equilibration_sensitivity_zscore"), stable_diag_header)
    stable_equil_pass_idx = findfirst(==("equilibration_sensitivity_pass"), stable_diag_header)
    @test stable_decorr_stride_idx !== nothing
    @test stable_decorr_steps_idx !== nothing
    @test stable_decorr_count_idx !== nothing
    @test stable_decorr_pass_idx !== nothing
    @test stable_equil_tail_idx !== nothing
    @test stable_equil_z_idx !== nothing
    @test stable_equil_pass_idx !== nothing
    @test all(row[stable_decorr_stride_idx] >= 1.0 for row in stable_diag_rows)
    @test all(row[stable_decorr_steps_idx] >= row[stable_decorr_stride_idx] for row in stable_diag_rows)
    @test all(row[stable_decorr_count_idx] >= 1.0 for row in stable_diag_rows)
    @test all(row[stable_decorr_pass_idx] in (0.0, 1.0) for row in stable_diag_rows)
    @test all(row[stable_equil_tail_idx] >= 1.0 for row in stable_diag_rows)
    @test all(!isnan(row[stable_equil_z_idx]) for row in stable_diag_rows)
    @test all(row[stable_equil_pass_idx] in (0.0, 1.0) for row in stable_diag_rows)
    stable_subspace_text = read(write_manuscript_text(stable_subspace_result, stable_subspace_outdir), String)
    @test occursin("stable-Fourier-subspace target diagnostic", stable_subspace_text)
    @test occursin("decorrelation stride", stable_subspace_text)
    @test occursin("burn-in sensitivity", stable_subspace_text)
    @test !occursin("mixed local pair-density moves", stable_subspace_text)
    @test isdefined(DFMMonteCarlo, :write_direct_stable_subspace_svg)
    stable_subspace_svg = DFMMonteCarlo.write_direct_stable_subspace_svg(
        stable_subspace_result.pressure_csv,
        stable_subspace_result.diagnostics_csv,
        stable_subspace_outdir,
    )
    @test isfile(stable_subspace_svg)
    stable_subspace_svg_text = read(stable_subspace_svg, String)
    @test occursin("<title>Constrained direct nonlinear MC pressure</title>", stable_subspace_svg_text)
    @test occursin("stable Fourier subspace", stable_subspace_svg_text)
    @test occursin("not unconstrained HE-URP", stable_subspace_svg_text)
    @test occursin("MC direct stable", stable_subspace_svg_text)
    @test occursin("DFM grid", stable_subspace_svg_text)
    @test occursin("CL nearest", stable_subspace_svg_text)
    @test occursin("Mean field", stable_subspace_svg_text)
    @test occursin("class='x-axis-label'", stable_subspace_svg_text)
    @test occursin("class='y-axis-label'", stable_subspace_svg_text)
    @test length(collect(eachmatch(r"class='y-major-tick'", stable_subspace_svg_text))) >= 9
    @test length(collect(eachmatch(r"class='y-minor-tick'", stable_subspace_svg_text))) >= 32
    stable_subspace_figure_script = joinpath(@__DIR__, "..", "scripts", "write_direct_stable_subspace_figure.jl")
    @test isfile(stable_subspace_figure_script)
    stable_subspace_figure_script_text = read(stable_subspace_figure_script, String)
    @test occursin("write_direct_stable_subspace_svg", stable_subspace_figure_script_text)
    @test occursin("DFM_DIRECT_STABLE_SUBSPACE_OUTDIR", stable_subspace_figure_script_text)

    stable_support_gate_outdir = mktempdir()
    stable_support_pressure_csv = DFMMonteCarlo._write_table(
        joinpath(stable_support_gate_outdir, "pressure_comparison.csv"),
        ["C", "mc_pressure", "mc_error", "dfm_pressure", "cl_pressure_nearest", "cl_error_nearest"],
        [(20.0, 218.0, 0.1, 218.0, 218.0, 0.1)],
    )
    stable_support_diagnostics_csv = DFMMonteCarlo._write_table(
        joinpath(stable_support_gate_outdir, "sampling_diagnostics.csv"),
        [
            "C",
            "stable_mode_proposals",
            "fourier_min_k2",
            "phi_min",
            "diagnostic_pass",
            "nonlinear_stability_pass",
            "sampling_budget_pass",
            "nonlinear_stability_min_energy_delta",
            "nonlinear_stability_min_phi",
            "production_steps",
            "production_sweeps",
            "effective_sample_size",
            "low_mode_structure_factor",
        ],
        [(20.0, 1.0, 30.0, 1.0e-8, 1.0, 1.0, 1.0, 0.0, 1.0, 16000.0, 31.25, 100.0, 0.0)],
    )
    stable_support_text_path = write_manuscript_text(
        SweepResult(stable_support_pressure_csv, stable_support_diagnostics_csv, run_config(:quick)),
        stable_support_gate_outdir)
    stable_support_text = read(stable_support_text_path, String)
    @test occursin("stable-Fourier-subspace target diagnostic", stable_support_text)
    @test occursin("k^2 >= 30", stable_support_text)
    @test !occursin("manuscript-supporting validation", stable_support_text)
    stable_support_audit = write_validation_audit(mktempdir(); direct_dirs=[stable_support_gate_outdir])
    stable_support_metadata = TOML.parsefile(stable_support_audit.toml_path)
    @test stable_support_metadata["verdict"]["direct_status"] == "diagnostic"
    @test stable_support_metadata["direct"]["stable_subspace_points"] == 1
    @test stable_support_metadata["direct"]["stable_subspace_passed_points"] == 1
    @test occursin("stable-Fourier-subspace direct outputs remain diagnostic",
        read(stable_support_audit.markdown_path, String))
    @test occursin("1 passed their local MC gates", read(stable_support_audit.markdown_path, String))

    sweep_outdir = mktempdir()
    sweep_result = run_pressure_sweep(sweep_outdir; mode=:quick, seed=19, sweeps=0.25, max_sweeps=0.25,
        C_values=[30.0], ess_target=0.0, chains=1, min_trace_samples=1)
    sweep_header, sweep_rows = DFMMonteCarlo._read_csv_numeric(sweep_result.diagnostics_csv)
    sweep_steps_idx = findfirst(==("production_steps"), sweep_header)
    sweep_equiv_idx = findfirst(==("production_sweeps"), sweep_header)
    sweep_max_equiv_idx = findfirst(==("max_production_sweeps"), sweep_header)
    sweep_stride_equiv_idx = findfirst(==("sample_stride_sweeps"), sweep_header)
    @test sweep_rows[1][sweep_steps_idx] == 2^3 * 0.25
    @test sweep_rows[1][sweep_equiv_idx] ≈ 0.25
    @test sweep_rows[1][sweep_max_equiv_idx] ≈ 0.25
    @test sweep_rows[1][sweep_stride_equiv_idx] ≈ 1 / 2^3

    production_pilot_outdir = mktempdir()
    production_pilot = run_pressure_sweep(production_pilot_outdir; mode=:production, Nx=3,
        C_values=[20.0], steps=2, max_steps=2, ess_target=0.0, chains=1,
        min_trace_samples=1, pressure_tolerance=Inf, burn_fraction=0.0)
    @test production_pilot.config.Nx == 3
    production_target = read(joinpath(production_pilot_outdir, "production_target.txt"), String)
    @test occursin("Manuscript target grid: Nx=64", production_target)
    @test occursin("Effective sampling grid: Nx=3", production_target)
    pilot_header, pilot_rows = DFMMonteCarlo._read_csv_numeric(production_pilot.diagnostics_csv)
    pilot_proposals_idx = findfirst(==("proposals_per_sweep"), pilot_header)
    pilot_sweeps_idx = findfirst(==("production_sweeps"), pilot_header)
    @test pilot_rows[1][pilot_proposals_idx] == 3^3
    @test pilot_rows[1][pilot_sweeps_idx] ≈ 2 / 3^3

    contract = load_pressure_contract(joinpath(@__DIR__, "..", "docs", "pressure_estimator_contract.toml"))
    guarded_grid = Grid3D(2, 6.4)
    guarded_params = HEURPParams(A=300.0, B=1.0, C=20.0, alpha=0.1)
    guarded_sample, guarded_summary = DFMMonteCarlo._run_arianna_pressure_trace_until_ess(
        guarded_grid,
        guarded_params,
        contract;
        seed=23,
        steps=2,
        ess_target=0.0,
        max_steps=12,
        min_trace_samples=5,
        pressure_tolerance=Inf,
    )
    @test guarded_summary.n >= 5
    @test guarded_sample.trace_samples_pass == 1.0
    @test guarded_sample.production_steps > 2
    @test guarded_sample.min_trace_samples == 5
    @test isinf(guarded_sample.pressure_tolerance)

    pressure_fig = write_pressure_svg(result.pressure_csv, outdir)
    diagnostics_fig = write_diagnostics_svg(result.diagnostics_csv, outdir)
    manuscript_text = write_manuscript_text(result, outdir)
    @test isfile(pressure_fig)
    @test isfile(diagnostics_fig)
    pressure_svg = read(pressure_fig, String)
    diagnostics_svg = read(diagnostics_fig, String)
    for svg_text in (pressure_svg, diagnostics_svg)
        @test occursin("class='plot-title'", svg_text)
        @test occursin("class='x-axis-label'", svg_text)
        @test occursin("class='y-axis-label'", svg_text)
        @test length(collect(eachmatch(r"class='x-major-tick'", svg_text))) >= 9
        @test length(collect(eachmatch(r"class='y-major-tick'", svg_text))) >= 9
        @test length(collect(eachmatch(r"class='x-tick-label'", svg_text))) >= 9
        @test length(collect(eachmatch(r"class='y-tick-label'", svg_text))) >= 9
        @test length(collect(eachmatch(r"class='x-minor-tick'", svg_text))) >= 32
        @test length(collect(eachmatch(r"class='y-minor-tick'", svg_text))) >= 32
    end
    @test occursin("Monte Carlo", read(manuscript_text, String))
    @test occursin("effective sample size", read(manuscript_text, String))
    @test occursin("burn fraction", read(manuscript_text, String))
    @test occursin("ESS target", read(manuscript_text, String))
    @test occursin("state points reached all active targets", read(manuscript_text, String))
    @test occursin("estimated required production budget", read(manuscript_text, String))
    @test occursin("convergence", read(manuscript_text, String))
    @test occursin("independent chains", read(manuscript_text, String))
    @test occursin("R-hat", read(manuscript_text, String))
    @test occursin("split R-hat", read(manuscript_text, String))
    @test occursin("diagnostic gate", read(manuscript_text, String))
    @test occursin("nonlinear stability", read(manuscript_text, String))
    @test occursin("positivity boundary", read(manuscript_text, String))
    @test occursin("acceptance window", read(manuscript_text, String))
    @test occursin("proposal amplitude", read(manuscript_text, String))
    @test occursin("multi-pair", read(manuscript_text, String))
    @test occursin("Fourier", read(manuscript_text, String))
    @test occursin("structure factor", read(manuscript_text, String))
    pressure_script_text = read(joinpath(@__DIR__, "..", "scripts", "run_pressure_sweep.jl"), String)
    @test occursin("DFM_MC_QUICK_TUNING_ROUNDS", pressure_script_text)
    @test occursin("DFM_MC_QUICK_TUNING_STEPS", pressure_script_text)
    @test occursin("DFM_MC_QUICK_TARGET_ACCEPTANCE", pressure_script_text)
    @test occursin("DFM_MC_QUICK_INITIAL_PROPOSAL_AMPLITUDE", pressure_script_text)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP", pressure_script_text)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP_ROUNDS", pressure_script_text)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP_STEPS", pressure_script_text)
    @test occursin("DFM_MC_QUICK_ADAPTIVE_WARMUP_GAIN", pressure_script_text)
    @test occursin("DFM_MC_QUICK_GLOBAL_MOVE_WEIGHT", pressure_script_text)
    @test occursin("DFM_MC_QUICK_GLOBAL_AMPLITUDE_SCALE", pressure_script_text)
    @test occursin("DFM_MC_QUICK_BOUNDARY_BARRIER_STRENGTH", pressure_script_text)
    trace_summary_script_text = read(joinpath(@__DIR__, "..", "scripts", "write_pressure_trace_summary.jl"), String)
    @test occursin("write_pressure_trace_summary", trace_summary_script_text)
    @test occursin("DFM_MC_TRACE_SUMMARY_OUTDIR", trace_summary_script_text)
    stability_scan_script_text = read(joinpath(@__DIR__, "..", "scripts", "run_direct_stability_scan.jl"), String)
    @test occursin("run_direct_stability_scan", stability_scan_script_text)
    @test occursin("DFM_DIRECT_STABILITY_SCAN_C_VALUES", stability_scan_script_text)
    boundary_scan_script_text = read(joinpath(@__DIR__, "..", "scripts", "run_direct_boundary_scan.jl"), String)
    @test occursin("run_direct_boundary_scan", boundary_scan_script_text)
    @test occursin("DFM_DIRECT_BOUNDARY_SCAN_MARGINS", boundary_scan_script_text)
    asymptotic_scan_script_text = read(joinpath(@__DIR__, "..", "scripts", "run_direct_boundary_asymptotic_scan.jl"), String)
    @test occursin("run_direct_boundary_asymptotic_scan", asymptotic_scan_script_text)
    @test occursin("DFM_DIRECT_BOUNDARY_ASYMPTOTIC_SCAN_C_VALUES", asymptotic_scan_script_text)
    boundary_ray_script_text = read(joinpath(@__DIR__, "..", "scripts", "run_direct_boundary_ray_mc.jl"), String)
    @test occursin("run_direct_boundary_ray_mc", boundary_ray_script_text)
    @test occursin("DFM_DIRECT_BOUNDARY_RAY_MC_EPSILON_MIN", boundary_ray_script_text)

    @test isdefined(DFMMonteCarlo, :write_direct_target_preflight)
    preflight_script_path = joinpath(@__DIR__, "..", "scripts", "write_direct_target_preflight.jl")
    @test isfile(preflight_script_path)
    if isdefined(DFMMonteCarlo, :write_direct_target_preflight) && isfile(preflight_script_path)
        preflight_boundary_dir = mktempdir()
        DFMMonteCarlo._write_table(
            joinpath(preflight_boundary_dir, "direct_boundary_asymptotic_scan.csv"),
            [
                "C",
                "mode_ax",
                "mode_ay",
                "mode_az",
                "log_boundary_slope",
                "min_cells",
                "min_eta",
                "max_eta",
                "asymptotic_unbounded_below",
            ],
            [
                (20.0, 1, 0, 0, 100.0, 1, -1.0, 1.0, 1.0),
                (20.0, 1, 1, 0, -1.0, 2, -1.0, 1.0, 0.0),
            ],
        )
        preflight_floor_dir = mktempdir()
        DFMMonteCarlo._write_table(
            joinpath(preflight_floor_dir, "direct_floor_sensitivity.csv"),
            [
                "phi_min",
                "mc_pressure",
                "estimated_required_steps",
                "estimated_required_sweeps",
                "production_budget_shortfall",
                "pressure_span",
            ],
            [
                (0.02, 1000.0, 100000.0, 1000.0, 10.0, 300.0),
                (0.1, 700.0, 20000.0, 200.0, 2.0, 300.0),
            ],
        )
        preflight_ray_dir = mktempdir()
        DFMMonteCarlo._write_table(
            joinpath(preflight_ray_dir, "direct_boundary_ray_summary.csv"),
            [
                "cutoff_dominated",
                "lower_cutoff_hit_fraction",
                "pressure_minus_dfm",
                "pressure_minus_cl",
                "jacobian_corrected_log_slope",
                "epsilon_min",
            ],
            [
                (1.0, 0.95, 1200.0, 1201.0, 100.0, 1.0e-8),
            ],
        )
        preflight = write_direct_target_preflight(mktempdir();
            boundary_asymptotic_dirs=[preflight_boundary_dir],
            floor_sensitivity_dirs=[preflight_floor_dir],
            boundary_ray_dirs=[preflight_ray_dir],
            floor_pressure_span_tolerance=10.0)
        @test isfile(preflight.markdown_path)
        @test isfile(preflight.toml_path)
        preflight_metadata = TOML.parsefile(preflight.toml_path)
        @test preflight_metadata["decision"]["status"] == "blocked_unconstrained_direct_production"
        @test preflight_metadata["decision"]["recommended_action"] == "reformulate_or_constrain_nonlinear_target"
        @test preflight_metadata["checks"]["boundary_asymptotic_unbounded_points"] == 1
        @test preflight_metadata["checks"]["floor_sensitivity_pressure_span"] == 300.0
        @test preflight_metadata["checks"]["boundary_ray_cutoff_dominated_points"] == 1
        @test "boundary_ray_cutoff_dominated" in preflight_metadata["decision"]["reasons"]
        preflight_text = read(preflight.markdown_path, String)
        @test occursin("Direct nonlinear target preflight", preflight_text)
        @test occursin("do not spend longer unconstrained direct Arianna budgets", preflight_text)
        @test occursin("asymptotically unbounded below", preflight_text)
        @test occursin("floor-sensitive", preflight_text)
        @test occursin("Boundary-ray MC check", preflight_text)
        preflight_script_text = read(preflight_script_path, String)
        @test occursin("write_direct_target_preflight", preflight_script_text)
        @test occursin("DFM_DIRECT_TARGET_PREFLIGHT_BOUNDARY_DIRS", preflight_script_text)
        @test occursin("DFM_DIRECT_TARGET_PREFLIGHT_BOUNDARY_RAY_DIRS", preflight_script_text)
        @test occursin("direct_heurp_boundary_asymptotic_scan_lowC", preflight_script_text)
    end

    audit_gaussian_outdir = mktempdir()
    gaussian_audit_result = run_gaussian_perturbation_mc(audit_gaussian_outdir; Nx=16,
        C_values=[20.0], samples=5000, max_samples=5000, chains=4, seed=37,
        pressure_tolerance=0.25)
    write_gaussian_mc_text(gaussian_audit_result, audit_gaussian_outdir)
    run_gaussian_sample_size_study(audit_gaussian_outdir; Nx=16, C=20.0,
        sample_counts=[100, 400], chains=4, seed=37, pressure_tolerance=0.25)
    field_audit_outdir = mktempdir()
    field_audit_result = run_gaussian_perturbation_mc(field_audit_outdir; Nx=16,
        C_values=[20.0], samples=1000, max_samples=1000, chains=2, seed=41,
        pressure_tolerance=Inf, estimator=:field)
    write_gaussian_mc_text(field_audit_result, field_audit_outdir)
    audit_outdir = mktempdir()
    audit = write_validation_audit(audit_outdir; gaussian_dirs=[audit_gaussian_outdir, field_audit_outdir],
        direct_dirs=[outdir])
    @test isfile(audit.markdown_path)
    @test isfile(audit.toml_path)
    audit_text = read(audit.markdown_path, String)
    @test occursin("Accepted evidence: finite-grid Gaussian perturbation MC", audit_text)
    @test occursin("Direct nonlinear Arianna status: diagnostic", audit_text)
    @test occursin("CL-context", audit_text)
    @test occursin("do not use diagnostic direct-density outputs as manuscript validation", audit_text)
    audit_metadata = TOML.parsefile(audit.toml_path)
    @test audit_metadata["verdict"]["gaussian_status"] == "accepted"
    @test audit_metadata["verdict"]["direct_status"] == "diagnostic"
    @test audit_metadata["gaussian"]["csv_count"] == 3
    @test audit_metadata["gaussian"]["total_points"] == 4
    @test audit_metadata["gaussian"]["passed_points"] == audit_metadata["gaussian"]["total_points"]
    @test audit_metadata["gaussian"]["sampler_counts"]["gaussian_mode_mc"] == 3
    @test audit_metadata["gaussian"]["sampler_counts"]["gaussian_field_mc"] == 1
    @test haskey(audit_metadata["gaussian"], "max_cl_model_residual")
    @test haskey(audit_metadata["gaussian"], "cl_context_passed_points")
    @test audit_metadata["gaussian"]["max_cl_model_residual"] >= 0.0
    @test audit_metadata["gaussian"]["cl_context_passed_points"] == audit_metadata["gaussian"]["finite_cl_points"]
    @test audit_metadata["direct"]["passed_points"] < audit_metadata["direct"]["total_points"]
    @test audit_metadata["readiness"]["overall_status"] == "partial"
    @test audit_metadata["readiness"]["perturbative_mc"] == "accepted"
    @test audit_metadata["readiness"]["unconstrained_direct_mc"] == "diagnostic"
    @test audit_metadata["readiness"]["manuscript_evidence"] == "finite_grid_gaussian_perturbation_mc"
    @test audit_metadata["readiness"]["full_objective_complete"] == false
    @test occursin("Full objective readiness: partial", audit_text)
    @test occursin("accepted perturbative MC evidence", audit_text)
    @test occursin("Gaussian samplers:", audit_text)
    @test occursin("gaussian_field_mc=1", audit_text)

    accepted_direct_outdir = mktempdir()
    DFMMonteCarlo._write_table(
        joinpath(accepted_direct_outdir, "sampling_diagnostics.csv"),
        [
            "C",
            "diagnostic_pass",
            "nonlinear_stability_pass",
            "sampling_budget_pass",
            "production_steps",
            "production_sweeps",
            "effective_sample_size",
            "low_mode_structure_factor",
            "target_model_code",
            "basin_variance_max",
            "stable_mode_proposals",
            "fourier_min_k2",
            "phi_min",
            "acceptance_rate",
            "radial_variance_acceptance_rate",
            "pressure_stderr",
            "pressure_tolerance",
            "chain_rhat",
            "split_chain_rhat",
            "density_variance_drift_zscore",
            "low_mode_structure_drift_zscore",
        ],
        [(20.0, 1.0, 1.0, 1.0, 32768.0, 64.0, 921.0, 1.36e-6, 2.0, Inf, 0.0, 0.0,
            1.0e-8, 0.492, 0.485, 0.00267, 0.01, 1.0, 1.00095, 0.683, 2.09)],
    )
    accepted_direct_audit = write_validation_audit(mktempdir();
        gaussian_dirs=[audit_gaussian_outdir, field_audit_outdir],
        direct_dirs=[outdir, accepted_direct_outdir])
    accepted_direct_text = read(accepted_direct_audit.markdown_path, String)
    accepted_direct_metadata = TOML.parsefile(accepted_direct_audit.toml_path)
    @test accepted_direct_metadata["verdict"]["direct_status"] == "accepted"
    @test accepted_direct_metadata["direct"]["accepted_unconstrained_points"] == 1
    @test accepted_direct_metadata["direct"]["accepted_unconstrained_hevk1_points"] == 1
    @test accepted_direct_metadata["direct"]["accepted_unconstrained_heurp_points"] == 0
    @test accepted_direct_metadata["readiness"]["unconstrained_direct_mc"] == "accepted"
    @test accepted_direct_metadata["readiness"]["unconstrained_heurp_direct_mc"] == "diagnostic"
    @test accepted_direct_metadata["readiness"]["unconstrained_hevk1_direct_mc"] == "accepted"
    @test accepted_direct_metadata["readiness"]["heurp_direct_required"] == true
    @test accepted_direct_metadata["readiness"]["required_direct_target_complete"] == false
    @test accepted_direct_metadata["readiness"]["full_objective_complete"] == false
    @test accepted_direct_metadata["readiness"]["overall_status"] == "partial"
    @test occursin("Direct nonlinear Arianna status: accepted", accepted_direct_text)
    @test occursin("accepted unconstrained direct state point", accepted_direct_text)

    accepted_heurp_direct_outdir = mktempdir()
    DFMMonteCarlo._write_table(
        joinpath(accepted_heurp_direct_outdir, "sampling_diagnostics.csv"),
        [
            "C",
            "diagnostic_pass",
            "nonlinear_stability_pass",
            "sampling_budget_pass",
            "production_steps",
            "production_sweeps",
            "effective_sample_size",
            "low_mode_structure_factor",
            "target_model_code",
            "basin_variance_max",
            "stable_mode_proposals",
            "fourier_min_k2",
            "phi_min",
        ],
        [(20.0, 1.0, 1.0, 1.0, 32768.0, 64.0, 921.0, 1.36e-6, 1.0, Inf, 0.0, 0.0,
            1.0e-8)],
    )
    accepted_heurp_direct_audit = write_validation_audit(mktempdir();
        direct_dirs=[accepted_direct_outdir, accepted_heurp_direct_outdir])
    accepted_heurp_direct_text = read(accepted_heurp_direct_audit.markdown_path, String)
    accepted_heurp_direct_metadata = TOML.parsefile(accepted_heurp_direct_audit.toml_path)
    @test accepted_heurp_direct_metadata["direct"]["accepted_unconstrained_points"] == 2
    @test accepted_heurp_direct_metadata["direct"]["accepted_unconstrained_hevk1_points"] == 1
    @test accepted_heurp_direct_metadata["direct"]["accepted_unconstrained_heurp_points"] == 1
    @test accepted_heurp_direct_metadata["readiness"]["unconstrained_heurp_direct_mc"] == "accepted"
    @test accepted_heurp_direct_metadata["readiness"]["unconstrained_hevk1_direct_mc"] == "accepted"
    @test accepted_heurp_direct_metadata["readiness"]["required_direct_target_complete"] == true
    @test accepted_heurp_direct_metadata["readiness"]["full_objective_complete"] == false
    @test occursin("HE-URP=1", accepted_heurp_direct_text)
    @test occursin("HE-VK1=1", accepted_heurp_direct_text)

    manuscript_audit_outdir = mktempdir()
    manuscript_audit = write_manuscript_evidence_text(manuscript_audit_outdir;
        gaussian_dirs=[audit_gaussian_outdir, field_audit_outdir], direct_dirs=[outdir])
    @test basename(manuscript_audit.markdown_path) == "mc_validation.md"
    @test basename(manuscript_audit.toml_path) == "mc_validation.toml"
    manuscript_audit_text = read(manuscript_audit.markdown_path, String)
    @test occursin("Accepted evidence: finite-grid Gaussian perturbation MC", manuscript_audit_text)
    @test occursin("Full objective readiness: partial", manuscript_audit_text)
    @test occursin("direct nonlinear arianna sampler remains a separate target", lowercase(manuscript_audit_text))
    manuscript_script_text = read(joinpath(@__DIR__, "..", "scripts", "write_manuscript_text.jl"), String)
    validation_script_text = read(joinpath(@__DIR__, "..", "scripts", "write_validation_audit.jl"), String)
    @test occursin("write_manuscript_evidence_text", manuscript_script_text)
    @test occursin("gaussian_field_configuration_perturbation_mc_Nx64_C20", manuscript_script_text)
    @test occursin("gaussian_field_configuration_perturbation_mc_Nx64_C20", validation_script_text)
    @test occursin("gaussian_field_configuration_sample_size_study_Nx64_C20", manuscript_script_text)
    @test occursin("gaussian_field_configuration_sample_size_study_Nx64_C20", validation_script_text)
    @test occursin("gaussian_field_configuration_analysis_comparison", manuscript_script_text)
    @test occursin("gaussian_field_configuration_analysis_comparison", validation_script_text)
    @test occursin("direct_nonlinear_unconstrained_Nx8_C20_probe", manuscript_script_text)
    @test occursin("direct_nonlinear_unconstrained_Nx8_C20_probe", validation_script_text)
    @test occursin("direct_heurp_arianna_Nx8_C20_radial_variance_attack_32768", manuscript_script_text)
    @test occursin("direct_heurp_arianna_Nx8_C20_radial_variance_attack_32768", validation_script_text)
    @test occursin("direct_heurp_unconstrained_nonlinear_attack_Nx8_C20_65536", manuscript_script_text)
    @test occursin("direct_heurp_unconstrained_nonlinear_attack_Nx8_C20_65536", validation_script_text)
    @test occursin("direct_heurp_boundary_barrier_lambda140_Nx8_C20_589824_burn85", manuscript_script_text)
    @test occursin("direct_heurp_boundary_barrier_lambda140_Nx8_C20_589824_burn85", validation_script_text)
    @test occursin("direct_heurp_boundary_ray_mc_Nx8_C20_eps1e-8_20000", manuscript_script_text)
    @test occursin("direct_heurp_boundary_ray_mc_Nx8_C20_eps1e-8_20000", validation_script_text)
    @test occursin("direct_heurp_basin_constrained_Nx8_C20_v1e-6_burn75_100000", manuscript_script_text)
    @test occursin("direct_heurp_basin_constrained_Nx8_C20_v1e-6_burn75_100000", validation_script_text)
    @test occursin("direct_vk1_arianna_Nx8_C20_radial_variance_attack_32768", manuscript_script_text)
    @test occursin("direct_vk1_arianna_Nx8_C20_radial_variance_attack_32768", validation_script_text)
    @test occursin("direct_boundary_asymptotic_scan_Nx8_representative", manuscript_script_text)
    @test occursin("direct_boundary_asymptotic_scan_Nx8_representative", validation_script_text)
    @test occursin("direct_heurp_boundary_asymptotic_scan_lowC", manuscript_script_text)
    @test occursin("direct_heurp_boundary_asymptotic_scan_lowC", validation_script_text)
    @test occursin("direct_floor_sensitivity_Nx4_C20_budgeted", manuscript_script_text)
    @test occursin("direct_floor_sensitivity_Nx4_C20_budgeted", validation_script_text)
    @test occursin("stable_subspace_direct_Nx8_C20_resumed_50000", manuscript_script_text)
    @test occursin("stable_subspace_direct_Nx8_C20_resumed_50000", validation_script_text)
    @test occursin("stable_subspace_direct_Nx8_C20_tuned_25000_chain1", manuscript_script_text)
    @test occursin("stable_subspace_direct_Nx8_C20_tuned_25000_chain1", validation_script_text)
    @test !occursin("joinpath(root, \"results\", \"pressure\")", manuscript_script_text)

    probe_outdir = mktempdir()
    DFMMonteCarlo._write_table(
        joinpath(probe_outdir, "direct_stability_probe.csv"),
        [
            "mode_ax",
            "mode_ay",
            "mode_az",
            "amplitude_index",
            "amplitude",
            "min_phi",
            "max_phi",
            "density_variance",
            "energy",
            "energy_delta",
            "pressure_observable",
            "alpha",
        ],
        [
            (1, 0, 0, 1, 0.0, 1.0, 1.0, 0.0, 10.0, 0.0, 220.0, 0.1),
            (2, 0, 0, 2, 0.5, 0.001, 1.5, 0.25, 5.0, -5.0, 230.0, 0.1),
        ],
    )
    audit_with_probe = write_validation_audit(mktempdir(); gaussian_dirs=[audit_gaussian_outdir],
        direct_dirs=[outdir], stability_probe_dirs=[probe_outdir])
    probe_text = read(audit_with_probe.markdown_path, String)
    @test occursin("Standalone stability probe", probe_text)
    @test occursin("Recommended direct-MC next step: reformulate or constrain the nonlinear target", probe_text)
    probe_metadata = TOML.parsefile(audit_with_probe.toml_path)
    @test probe_metadata["direct"]["stability_probe_points"] == 2
    @test probe_metadata["direct"]["stability_probe_min_energy_delta"] == -5.0
    @test probe_metadata["direct"]["stability_probe_min_phi"] == 0.001
    @test probe_metadata["verdict"]["direct_recommendation"] == "reformulate_or_constrain_nonlinear_target"
    audit_with_basin_and_probe = write_validation_audit(mktempdir(); gaussian_dirs=[audit_gaussian_outdir],
        direct_dirs=[basin_outdir], stability_probe_dirs=[probe_outdir])
    basin_probe_metadata = TOML.parsefile(audit_with_basin_and_probe.toml_path)
    @test basin_probe_metadata["verdict"]["direct_recommendation"] == "reformulate_or_constrain_nonlinear_target"

    stability_scan_outdir = mktempdir()
    stability_scan = run_direct_stability_scan(stability_scan_outdir; Nx=4,
        C_values=[20.0, 30.0], amplitude_count=5,
        modes=[(1, 0, 0), (2, 0, 0)], boundary_margin=0.001)
    @test isfile(stability_scan.csv_path)
    @test isfile(stability_scan.markdown_path)
    scan_header, scan_rows = DFMMonteCarlo._read_csv_numeric(stability_scan.csv_path)
    @test scan_header == [
        "C",
        "stability_pass",
        "min_energy_delta",
        "min_phi",
        "mode_ax",
        "mode_ay",
        "mode_az",
        "points",
        "boundary_margin",
        "phi_min",
        "energy_tolerance",
        "stable_fourier_min_k2",
        "mean_field_pressure",
        "dfm_pressure",
        "pressure_at_min_energy",
        "pressure_minus_mean_field",
        "pressure_minus_dfm",
        "basin_variance_max",
        "alpha",
    ]
    @test length(scan_rows) == 2
    scan_C_idx = findfirst(==("C"), scan_header)
    scan_delta_idx = findfirst(==("min_energy_delta"), scan_header)
    scan_phi_idx = findfirst(==("min_phi"), scan_header)
    scan_k2_idx = findfirst(==("stable_fourier_min_k2"), scan_header)
    scan_dfm_idx = findfirst(==("dfm_pressure"), scan_header)
    @test Set(row[scan_C_idx] for row in scan_rows) == Set([20.0, 30.0])
    @test all(isfinite(row[scan_delta_idx]) for row in scan_rows)
    @test all(row[scan_phi_idx] > 0.0 for row in scan_rows)
    @test all(row[scan_k2_idx] > 0.0 for row in scan_rows)
    @test all(isfinite(row[scan_dfm_idx]) for row in scan_rows)
    scan_text = read(stability_scan.markdown_path, String)
    @test occursin("Direct nonlinear stability scan", scan_text)
    @test occursin("full nonlinear HE-URP target", scan_text)

    vk1_scan_outdir = mktempdir()
    vk1_scan = run_direct_vk1_stability_scan(vk1_scan_outdir; Nx=4,
        C_values=[20.0, 40.0], amplitude_count=5,
        modes=[(1, 0, 0), (2, 0, 0)], boundary_margin=0.001)
    @test isfile(vk1_scan.csv_path)
    @test isfile(vk1_scan.markdown_path)
    vk1_header, vk1_rows = DFMMonteCarlo._read_csv_numeric(vk1_scan.csv_path)
    @test vk1_header == [
        "C",
        "stability_pass",
        "min_energy_delta",
        "min_phi",
        "mode_ax",
        "mode_ay",
        "mode_az",
        "points",
        "boundary_margin",
        "phi_min",
        "energy_tolerance",
        "mean_field_pressure",
        "dfm_pressure",
        "pressure_at_min_energy",
        "pressure_minus_mean_field",
        "pressure_minus_dfm",
        "alpha",
        "target_model_code",
    ]
    @test length(vk1_rows) == 2
    vk1_C_idx = findfirst(==("C"), vk1_header)
    vk1_delta_idx = findfirst(==("min_energy_delta"), vk1_header)
    vk1_pressure_idx = findfirst(==("pressure_at_min_energy"), vk1_header)
    vk1_target_idx = findfirst(==("target_model_code"), vk1_header)
    @test Set(row[vk1_C_idx] for row in vk1_rows) == Set([20.0, 40.0])
    @test all(isfinite(row[vk1_delta_idx]) for row in vk1_rows)
    @test all(isfinite(row[vk1_pressure_idx]) for row in vk1_rows)
    @test all(row[vk1_target_idx] == 2.0 for row in vk1_rows)
    vk1_text = read(vk1_scan.markdown_path, String)
    @test occursin("Direct nonlinear HE-VK1 stability scan", vk1_text)
    @test occursin("Eq. HE-VK1", vk1_text)

    boundary_scan_outdir = mktempdir()
    boundary_scan = run_direct_boundary_scan(boundary_scan_outdir; Nx=4, C=20.0,
        boundary_margins=[0.1, 0.01, 0.001], amplitude_count=5,
        modes=[(1, 0, 0), (2, 0, 0)])
    @test isfile(boundary_scan.csv_path)
    @test isfile(boundary_scan.markdown_path)
    boundary_header, boundary_rows = DFMMonteCarlo._read_csv_numeric(boundary_scan.csv_path)
    @test boundary_header == [
        "boundary_margin",
        "stability_pass",
        "min_energy_delta",
        "min_phi",
        "mode_ax",
        "mode_ay",
        "mode_az",
        "points",
        "phi_min",
        "energy_tolerance",
        "C",
        "stable_fourier_min_k2",
        "mean_field_pressure",
        "dfm_pressure",
        "pressure_at_min_energy",
        "pressure_minus_mean_field",
        "pressure_minus_dfm",
        "basin_variance_max",
        "alpha",
    ]
    @test length(boundary_rows) == 3
    margin_idx = findfirst(==("boundary_margin"), boundary_header)
    delta_idx = findfirst(==("min_energy_delta"), boundary_header)
    min_phi_idx = findfirst(==("min_phi"), boundary_header)
    pressure_idx = findfirst(==("pressure_at_min_energy"), boundary_header)
    pressure_mf_idx = findfirst(==("pressure_minus_mean_field"), boundary_header)
    pressure_dfm_idx = findfirst(==("pressure_minus_dfm"), boundary_header)
    mean_field_idx = findfirst(==("mean_field_pressure"), boundary_header)
    dfm_idx = findfirst(==("dfm_pressure"), boundary_header)
    @test [row[margin_idx] for row in boundary_rows] == [0.1, 0.01, 0.001]
    @test all(row[min_phi_idx] > 0.0 for row in boundary_rows)
    @test boundary_rows[end][delta_idx] <= boundary_rows[1][delta_idx]
    @test all(isfinite(row[pressure_idx]) for row in boundary_rows)
    @test all(row[pressure_mf_idx] ≈ row[pressure_idx] - row[mean_field_idx] for row in boundary_rows)
    @test all(row[pressure_dfm_idx] ≈ row[pressure_idx] - row[dfm_idx] for row in boundary_rows)
    boundary_text = read(boundary_scan.markdown_path, String)
    @test occursin("Direct nonlinear boundary scan", boundary_text)
    @test occursin("positivity wall", boundary_text)
    @test occursin("pressure observable at the same minimum-energy points", boundary_text)

    asymptotic_scan_outdir = mktempdir()
    asymptotic_scan = run_direct_boundary_asymptotic_scan(asymptotic_scan_outdir; Nx=4,
        C_values=[20.0, 40.0], modes=[(2, 0, 0)])
    @test isfile(asymptotic_scan.csv_path)
    @test isfile(asymptotic_scan.markdown_path)
    asymptotic_header, asymptotic_rows = DFMMonteCarlo._read_csv_numeric(asymptotic_scan.csv_path)
    @test asymptotic_header == [
        "C",
        "mode_ax",
        "mode_ay",
        "mode_az",
        "log_boundary_slope",
        "min_cells",
        "log_coordinate_jacobian_slope",
        "jacobian_corrected_log_slope",
        "jacobian_corrected_unbounded_below",
        "min_eta",
        "max_eta",
        "asymptotic_unbounded_below",
        "slope_tolerance",
        "stable_fourier_min_k2",
        "mean_field_pressure",
        "dfm_pressure",
        "alpha",
    ]
    @test length(asymptotic_rows) == 2
    slope_idx = findfirst(==("log_boundary_slope"), asymptotic_header)
    unbounded_idx = findfirst(==("asymptotic_unbounded_below"), asymptotic_header)
    min_cells_idx = findfirst(==("min_cells"), asymptotic_header)
    jacobian_slope_idx = findfirst(==("log_coordinate_jacobian_slope"), asymptotic_header)
    corrected_slope_idx = findfirst(==("jacobian_corrected_log_slope"), asymptotic_header)
    corrected_unbounded_idx = findfirst(==("jacobian_corrected_unbounded_below"), asymptotic_header)
    @test all(row[slope_idx] > 0.0 for row in asymptotic_rows)
    @test all(row[unbounded_idx] == 1.0 for row in asymptotic_rows)
    @test all(row[min_cells_idx] > 0.0 for row in asymptotic_rows)
    @test all(row[jacobian_slope_idx] == row[min_cells_idx] for row in asymptotic_rows)
    @test all(row[corrected_slope_idx] ≈ row[slope_idx] - row[min_cells_idx] for row in asymptotic_rows)
    @test all(row[corrected_unbounded_idx] == 1.0 for row in asymptotic_rows)
    asymptotic_text = read(asymptotic_scan.markdown_path, String)
    @test occursin("Direct nonlinear boundary asymptotic scan", asymptotic_text)
    @test occursin("H ~ slope * log(epsilon)", asymptotic_text)
    @test occursin("log-coordinate Jacobian", asymptotic_text)
    @test occursin("unbounded below", asymptotic_text)

    boundary_ray_outdir = mktempdir()
    boundary_ray = run_direct_boundary_ray_mc(boundary_ray_outdir; Nx=4, C=20.0,
        mode=(2, 0, 0), steps=240, burn_fraction=0.25, sample_stride=4,
        epsilon_min=1.0e-4, epsilon_max=0.2, initial_epsilon=0.1,
        log_epsilon_proposal_width=0.8)
    @test isfile(boundary_ray.samples_csv)
    @test isfile(boundary_ray.summary_csv)
    @test isfile(boundary_ray.markdown_path)
    @test isfile(boundary_ray.metadata_path)
    ray_header, ray_rows = DFMMonteCarlo._read_csv_numeric(boundary_ray.summary_csv)
    @test ray_header == [
        "C",
        "mode_ax",
        "mode_ay",
        "mode_az",
        "steps",
        "burn_fraction",
        "sample_stride",
        "stored_samples",
        "acceptance_rate",
        "epsilon_min",
        "epsilon_max",
        "minimum_sample_epsilon",
        "mean_epsilon",
        "epsilon_stderr",
        "epsilon_effective_sample_size",
        "lower_cutoff_hit_fraction",
        "cutoff_dominated",
        "pressure_mean",
        "pressure_stderr",
        "pressure_effective_sample_size",
        "mean_field_pressure",
        "dfm_pressure",
        "cl_pressure_nearest",
        "cl_error_nearest",
        "pressure_minus_dfm",
        "pressure_minus_cl",
        "log_weight_mean",
        "log_weight_stderr",
        "log_weight_effective_sample_size",
        "log_boundary_slope",
        "min_cells",
        "log_coordinate_jacobian_slope",
        "jacobian_corrected_log_slope",
        "jacobian_corrected_unbounded_below",
        "alpha",
    ]
    @test length(ray_rows) == 1
    ray_row = only(ray_rows)
    stored_idx = findfirst(==("stored_samples"), ray_header)
    acceptance_idx = findfirst(==("acceptance_rate"), ray_header)
    eps_min_idx = findfirst(==("epsilon_min"), ray_header)
    min_sample_idx = findfirst(==("minimum_sample_epsilon"), ray_header)
    pressure_idx = findfirst(==("pressure_mean"), ray_header)
    corrected_unbounded_idx = findfirst(==("jacobian_corrected_unbounded_below"), ray_header)
    @test ray_row[stored_idx] > 0.0
    @test 0.0 <= ray_row[acceptance_idx] <= 1.0
    @test ray_row[min_sample_idx] >= ray_row[eps_min_idx]
    @test isfinite(ray_row[pressure_idx])
    @test ray_row[corrected_unbounded_idx] == 1.0
    sample_header, sample_rows = DFMMonteCarlo._read_csv_numeric(boundary_ray.samples_csv)
    @test "log_weight_log_epsilon_measure" in sample_header
    @test length(sample_rows) == Int(ray_row[stored_idx])
    ray_text = read(boundary_ray.markdown_path, String)
    @test occursin("Direct HE-URP boundary-ray Monte Carlo", ray_text)
    @test occursin("direct-measure log-coordinate Jacobian", ray_text)
    ray_metadata = TOML.parsefile(boundary_ray.metadata_path)
    @test ray_metadata["run"]["mode"] == "direct_boundary_ray_mc"
    @test ray_metadata["run"]["target_model"] == "HE-URP"

    audit_with_asymptotic = write_validation_audit(mktempdir(); gaussian_dirs=[audit_gaussian_outdir],
        direct_dirs=[outdir], boundary_asymptotic_dirs=[asymptotic_scan_outdir],
        boundary_ray_dirs=[boundary_ray_outdir])
    asymptotic_audit_text = read(audit_with_asymptotic.markdown_path, String)
    @test occursin("Boundary asymptotic scan", asymptotic_audit_text)
    @test occursin("unbounded below", asymptotic_audit_text)
    @test occursin("Boundary-ray direct MC", asymptotic_audit_text)
    asymptotic_audit_metadata = TOML.parsefile(audit_with_asymptotic.toml_path)
    @test asymptotic_audit_metadata["direct"]["boundary_asymptotic_points"] == 2
    @test asymptotic_audit_metadata["direct"]["boundary_asymptotic_unbounded_points"] == 2
    @test asymptotic_audit_metadata["direct"]["boundary_asymptotic_jacobian_corrected_unbounded_points"] == 2
    @test asymptotic_audit_metadata["direct"]["boundary_ray_points"] == 1
    @test asymptotic_audit_metadata["direct"]["boundary_ray_cutoff_dominated_points"] in (0, 1)
    @test occursin("Jacobian-corrected", asymptotic_audit_text)
    @test asymptotic_audit_metadata["verdict"]["direct_recommendation"] == "reformulate_or_constrain_nonlinear_target"

    accepted_with_asymptotic = write_validation_audit(mktempdir();
        direct_dirs=[accepted_direct_outdir], boundary_asymptotic_dirs=[asymptotic_scan_outdir])
    accepted_with_asymptotic_text = read(accepted_with_asymptotic.markdown_path, String)
    @test occursin("Direct nonlinear Arianna status: accepted", accepted_with_asymptotic_text)
    @test occursin("Target-specific direct status: HE-URP diagnostic, HE-VK1 accepted", accepted_with_asymptotic_text)
    @test occursin("Jacobian-corrected boundary scan", accepted_with_asymptotic_text)

    gate_outdir = mktempdir()
    pressure_csv = DFMMonteCarlo._write_table(
        joinpath(gate_outdir, "pressure_comparison.csv"),
        ["C", "mc_pressure", "mc_error", "dfm_pressure", "cl_pressure_nearest", "cl_error_nearest"],
        [(20.0, 218.0, 0.1, 218.0, 218.0, 0.1)],
    )
    diagnostics_csv = DFMMonteCarlo._write_table(
        joinpath(gate_outdir, "sampling_diagnostics.csv"),
        [
            "C",
            "acceptance_rate",
            "equilibration_steps",
            "sample_stride",
            "trace_samples",
            "integrated_autocorrelation_time",
            "effective_sample_size",
            "pressure_stderr",
            "pressure_stderr_independent",
            "proposal_amplitude",
            "target_acceptance",
            "tuning_rounds",
            "tuning_steps",
            "tuning_acceptance_rate",
            "collective_proposal_amplitude",
            "collective_pair_count",
            "collective_move_weight",
            "fourier_mode_amplitude",
            "fourier_max_mode",
            "fourier_move_weight",
            "ess_target",
            "converged",
            "production_steps",
            "max_production_steps",
            "convergence_attempts",
            "chain_count",
            "chains_converged",
            "chain_rhat",
            "chain_between_stderr",
            "chain_within_stderr",
            "acceptance_min",
            "acceptance_max",
            "acceptance_pass",
            "rhat_threshold",
            "rhat_pass",
            "ess_pass",
            "min_trace_samples",
            "trace_samples_pass",
            "pressure_tolerance",
            "pressure_error_pass",
            "diagnostic_pass",
        ],
        [(20.0, 0.44, 100, 1, 100, 1.0, 5.0, 0.1, 0.1, 0.1, 0.44, 8, 0.44,
            0.01, 4, 0.2, 0.005, 2, 0.1, 10.0, 0.0, 200, 200, 1, 2, 1, 1.0, 0.1, 0.1,
            0.2, 0.8, 1.0, 1.1, 1.0, 0.0, 100, 1.0, 0.05, 0.0, 0.0)],
    )
    gated_text = write_manuscript_text(SweepResult(pressure_csv, diagnostics_csv, run_config(:quick)), gate_outdir)
    @test occursin("Current status: support-metadata diagnostic", read(gated_text, String))

    lower_support_gate_outdir = mktempdir()
    lower_support_pressure_csv = DFMMonteCarlo._write_table(
        joinpath(lower_support_gate_outdir, "pressure_comparison.csv"),
        ["C", "mc_pressure", "mc_error", "dfm_pressure", "cl_pressure_nearest", "cl_error_nearest"],
        [(20.0, 218.0, 0.1, 218.0, 218.0, 0.1)],
    )
    lower_support_diagnostics_csv = DFMMonteCarlo._write_table(
        joinpath(lower_support_gate_outdir, "sampling_diagnostics.csv"),
        ["C", "phi_min", "diagnostic_pass"],
        [(20.0, 0.05, 1.0)],
    )
    lower_support_text_path = write_manuscript_text(
        SweepResult(lower_support_pressure_csv, lower_support_diagnostics_csv, run_config(:quick)),
        lower_support_gate_outdir)
    lower_support_text = read(lower_support_text_path, String)
    @test occursin("lower-bound-regularized target diagnostic", lower_support_text)
    @test !occursin("manuscript-supporting validation", lower_support_text)
end
