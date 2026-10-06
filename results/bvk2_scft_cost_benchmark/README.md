# AQCE and SCFT: complete lamellar-solve timings

This auxiliary dataset measures complete lamellar solves, including cell
relaxation. Every measured repetition constructs a fresh seeded field,
converges it, and relaxes a single-period cell. `BVK2` denotes AQCE in the
timing records. For individual evaluations and warm fixed-cell solves, see the
[separate performance dataset](../bvk2_scft_fixedcell_performance/README.md).

## Protocol

- Same host and one pinned logical CPU; Julia, BLAS, and FFTW each use one thread.
- Four accepted symmetric states: AB at `chiN=20,30` and ABA at `chiN=30,40`.
- BVK2 uses the frozen publication setting `c2=0.16`, `nx=128`, and its analytic/shared cell-root workflow.
- Polyorder uses `OSF`, `ds=0.005`, `maxDeltaX=0.04`, residual tolerance `1e-6`, cell-stress tolerance `1e-5`, and the same seeded Anderson/variable-cell workflow as the accepted references.
- Each isolated process performs one compilation warm-up followed by three fresh measured solves using the same fixed RNG seed. A single-period gate rejects tiled multi-period cells. Solver wall time is the median of those three; process wall time and peak RSS include imports, JIT state, warm-up, and measurements.

The spatial discretizations are the accepted publication settings, not artificially identical grid counts. Therefore this measures cost to reproduce the claimed results, not asymptotic cost at equal degrees of freedom.

## Results

| case | BVK2 time (s) | SCFT time (s) | BVK2/SCFT time | BVK2 peak RSS (MiB) | SCFT peak RSS (MiB) | SCFT/BVK2 RSS |
|---|---:|---:|---:|---:|---:|---:|
| AB, chiN=20 | 4.206 | 0.513 | 8.192 | 1077.2 | 1221.0 | 1.134 |
| AB, chiN=30 | 4.084 | 0.427 | 9.565 | 1077.0 | 1193.3 | 1.108 |
| ABA, chiN=30 | 3.501 | 0.210 | 16.685 | 1147.9 | 1164.4 | 1.014 |
| ABA, chiN=40 | 3.725 | 0.230 | 16.165 | 1151.4 | 1202.4 | 1.044 |

Across these four one-dimensional publication states, BVK2 is **12.865x slower** at the median paired post-warm-up solver time, while the median SCFT/BVK2 peak-RSS ratio is **1.076**.

Polyorder is `8.2--16.7x` faster in the measured post-warm-up solver stage. Its isolated process nevertheless takes longer from cold startup because Polyorder's compilation warm-up is about twice as long. Peak RSS is similar, with SCFT `1--13%` higher, whereas SCFT allocates only about `7--9%` as many bytes during the measured solver call.

Interpretation must remain scoped: the result supports an empirical lamellar cost statement on this CPU and these accuracy contracts. It does not establish GPU behavior, three-dimensional morphology cost, or asymptotic complexity. Most importantly, the propagator-free BVK2 formulation does **not** currently deliver an end-to-end speed advantage in this production workflow; Polyorder's mature Anderson/variable-cell algorithm more than offsets its MDE work. This is an implementation-level result, not a reversal of the formal per-evaluation complexity argument.

## Reproduction

These are new timing calculations, not figure-only commands. They require the
[SCFT environment](../../SOFTWARE.md#polyorderjl-availability) and should be run
from the repository root. CPU affinity in the original command is host-specific;
choose an available logical CPU on your system when repeating the measurement.

```bash
julia --project=. scripts/write_bvk2_scft_cost_benchmark.jl --recompute=true
```
