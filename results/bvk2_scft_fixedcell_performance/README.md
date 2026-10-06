# AQCE and SCFT: evaluation and fixed-cell timings

This auxiliary dataset separates individual energy/force evaluation timings
from the cost of following a sequence of states at a fixed cell size. `BVK2`
denotes AQCE in the records. Cell optimization is excluded; timings for
complete stress-free solves are in the
[separate lamellar-solve benchmark](../bvk2_scft_cost_benchmark/README.md).

## Kernel benchmark

- BVK2 timing unit: one energy plus analytic chemical-potential evaluation.
- SCFT timing unit: one `Polyorder.q!` evaluation, including all OSF MDE propagations, density integration, and force construction.
- Equal full grids, no DCT/CFFT/symmetry compression, one CPU thread, `ds=0.005` (about 200 contour steps).
- 1D grids: `32, 64, 128, 256, 512, 1024`; 3D grids: `12, 16, 24, 32, 40, 48`.

SCFT/BVK2 kernel-time ratios span `3.71--9.38x` in 1D and `5.18--13.12x` in 3D.

## Fixed-cell warm continuation

- AB: anchor `chiN=20`, then `21:30`, fixed `D/Rg=4.2`.
- ABA: anchor `chiN=30`, then `31:40`, fixed `D/Rg=2.9`.
- Anchors and one compilation trajectory are excluded. Every measured point is seeded from the immediately previous accepted point.
- BVK2 uses the previous density profile. Polyorder uses `reset(scft, new_system)`, which carries the converged auxiliary fields and fixed lattice forward.
- BVK2 warm L-BFGS is capped at `50` iterations and then uses its existing residual-Newton polish. This is an accuracy-qualified stop, not a weakened acceptance rule.
- AB cap audit at chiN=21: `|Delta F|=1.036e-15`, profile RMS `4.220e-10`, capped residuals `(1.062e-05, 3.161e-05)` versus the 4000-iteration reference.
- ABA cap audit at chiN=31: `|Delta F|=3.751e-16`, profile RMS `4.961e-11`, capped residuals `(1.959e-06, 7.847e-06)` versus the 4000-iteration reference.
- Reported time is the complete marginal point cost. For SCFT it includes reset/setup plus fixed-cell `solve!`; for BVK2 it is the complete fixed-period minimizer call.

Median SCFT/BVK2 marginal-time ratios are `4.52x` for AB and `3.46x` for ABA.

## Interpretation

The kernel ratio tests the formal removal of the `Ns` MDE factor. The continuation ratio tests whether that lower evaluation cost survives a production-like warm solve. Neither result includes cell optimization, so it must not be compared directly with the earlier stress-free end-to-end benchmark.

## Reproduction

The command below launches new timing calculations and requires the
[SCFT environment](../../SOFTWARE.md#polyorderjl-availability). Run it from the
repository root. The CPU index `12` is the original host's affinity setting;
choose an available logical CPU on your machine.

```bash
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 12 \
  julia --project=. scripts/write_bvk2_scft_fixedcell_performance_benchmark.jl
```
