# Initial AQCE transfer check for a symmetric ABA triblock

This directory contains an additional two-state architecture-transfer
check. For the seven-state comparison and Figures 7–8, use the
[comprehensive ABA dataset](../bvk2_aba_comprehensive_validation/README.md).
Here `BVK2` identifies AQCE. Its coefficient remains
`c2 = 0.16`; the change from AB to ABA enters through the ideal-chain correlation
and the corresponding quadratic-response coefficients. ABA SCFT results supply
the validation reference.

Model fingerprint: `model=BVK2|architecture=ABA|schema=bvk2-rpa-eta-v1|law_mode=adaptive|discretization=central_eta_forward_edge_exact_v1|coordinate=eta_psi|units=lengths_and_b_same_units|c2=0.16|f=0.5|N=1|b=1|k_star_definition=k_star^4=12*f*(1-f)*A_psi/b^2|k_star=6.4474195909412515`

| χN | BVK2 D/Rg | SCFT D/Rg | signed error | profile RMS | correlation | accepted |
|---:|---:|---:|---:|---:|---:|:---:|
| 30 | 2.768804 | 2.822245 | -1.894% | 0.0070 | 0.99991 | yes |
| 40 | 2.955394 | 3.022881 | -2.233% | 0.0073 | 0.99993 | yes |

Predeclared gates: `|ΔD|/D_SCFT ≤ 5%`, aligned profile RMS `≤ 0.05`, correlation `≥ 0.98`, converged/local-minimum BVK2 fields, and SCFT cell stress `≤ 1e-4`.

Result: **accepted**. Maximum period error is 2.233%, maximum profile RMS is 0.0073, and minimum correlation is 0.99991.

Scope: this establishes quantitative no-refit structural transfer for symmetric linear ABA lamellae at two strong-segregation state points. It does not establish an ABA phase diagram, arbitrary-topology universality, or Frank–Kasper stabilization.

The numerical entry point is
`julia --project=. scripts/write_bvk2_aba_no_refit_validation.jl`, run from the
repository root. It requires the [private SCFT environment](../../SOFTWARE.md#polyorderjl-availability).
Reading the supplied tables does not require that software.
