# Diblock weak-response map

This artifact extends the representative Gaussian-kernel panels into a continuous composition scan.  It is the lamellar/saddle analogue of the Uneyama-Doi weak-response and phase-boundary checks: the question is whether the real-density constructions reproduce the full diblock RPA spinodal and weak-segregation wavelength before any finite-amplitude saddle claim is made.

The scan uses `fA <= 0.5` because AB exchange symmetry covers the other half of composition space.

## Aggregate weak-response errors

| model | max abs spinodal error | mean abs spinodal error | max abs period error | mean abs period error | boundary-limited minima | all converged |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Full RPA kernel (BURP-TI Gaussian limit) | 0 | 0 | 0 | 0 | 0 | `true` |
| Literal nonlinear UD Hessian | 0.156 | 0.11507 | 0.045325 | 0.029518 | 0 | `true` |
| Nominal reduced kernel shared by BVK1 and BVK2 | 0.014562 | 0.00692 | 0.045325 | 0.029518 | 0 | `true` |
| OK | 1.170e-15 | 2.882e-16 | 0.045325 | 0.029518 | 0 | `true` |

## Interpretation

1. BURP-TI is exactly the full diblock RPA kernel at Gaussian order for every scanned composition, so its spinodal and weak-period errors are numerically zero.
2. The literal nonlinear Uneyama-Doi functional has a larger local Hessian than the nominal reduced vertex. The shift changes its spinodal but not its preferred weak-segregation wavevector.
3. BVK1 and BVK2 share the nominal reduced kernel. Their finite-amplitude differences therefore arise beyond Gaussian order.
4. Here OK denotes the Ohta--Kawasaki comparator with the published coefficients of Liu et al. Its local quadratic coefficient reproduces the RPA spinodal, but its gradient and nonlocal coefficients leave a composition-dependent error in the preferred weak-segregation period. This distinction helps separate instability-threshold matching from finite-amplitude accuracy.

Outputs:

- `weak_response_map.csv`: minima and relative errors for every model/composition pair.
- `weak_response_summary.csv`: aggregate error metrics by model.
- `diblock_weak_response_map.svg`: spinodal, weak-period, and error curves.

Settings: `N=1.0`, `b=1.0`, `kRg_min=0.05`, `kRg_max=8.0`, `f_count=81`.
