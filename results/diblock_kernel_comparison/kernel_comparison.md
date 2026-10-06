# Diblock Gaussian kernel comparison

This artifact is the Uneyama-Doi-style weak-response companion to the stress-free lamella profile comparison. It compares the chi-independent parts of the second-order composition vertices for an incompressible AB diblock.

The plotted models are:

- `Exact RPA / BURP-TI Gaussian`: `K_RPA(k) = v^T h(k)^{-1} v`, with `v=(1,-1)`. BURP-TI reduces to this kernel because `eta_psi(phi)=phi-f+O((phi-f)^2)`.
- `Literal nonlinear Uneyama-Doi Hessian`: `A_psi/k^2 + C_UD^(2) + b^2 k^2/[12 f(1-f)]`, where `C_UD^(2) = C_AA + C_BB - C_AB/[f(1-f)]`.
- `BVK1/BVK2 reduced Gaussian`: `A_psi/k^2 + C_psi + b^2 k^2/[12 f(1-f)]`. The adaptive stiffness activates beyond second order.
- `OK`: the Ohta--Kawasaki comparator with the published local, gradient, and nonlocal coefficients from Liu et al., Eqs. 15--17. The local coefficient reproduces the RPA spinodal, while the gradient and nonlocal coefficients determine an approximate preferred wavevector.

## Minimum-kernel summary

| fA | model | kRg at min | L/Rg | chiN spinodal | relative spinodal error |
| ---: | --- | ---: | ---: | ---: | ---: |
| 0.1 | Full RPA kernel (BURP-TI Gaussian limit) | 2.35621 | 2.66665 | 74.3314 | 0 |
| 0.1 | Literal nonlinear UD Hessian | 2.40281 | 2.61493 | 85.9268 | 0.155997 |
| 0.1 | Nominal reduced kernel shared by BVK1 and BVK2 | 2.40281 | 2.61493 | 73.2725 | -0.014245 |
| 0.1 | OK | 2.40281 | 2.61493 | 74.3314 | -1.9118e-16 |
| 0.25 | Full RPA kernel (BURP-TI Gaussian limit) | 2.05276 | 3.06084 | 18.1719 | 0 |
| 0.25 | Literal nonlinear UD Hessian | 2 | 3.14159 | 20.2577 | 0.114779 |
| 0.25 | Nominal reduced kernel shared by BVK1 and BVK2 | 2 | 3.14159 | 18.0355 | -0.0075101 |
| 0.25 | OK | 2 | 3.14159 | 18.1719 | 5.8652e-16 |
| 0.3 | Full RPA kernel (BURP-TI Gaussian limit) | 2.00984 | 3.12621 | 14.6349 | 0 |
| 0.3 | Literal nonlinear UD Hessian | 1.94413 | 3.23187 | 16.2283 | 0.108878 |
| 0.3 | Nominal reduced kernel shared by BVK1 and BVK2 | 1.94413 | 3.23187 | 14.5843 | -0.00345615 |
| 0.3 | OK | 1.94413 | 3.23187 | 14.6349 | -1.2138e-16 |
| 0.35 | Full RPA kernel (BURP-TI Gaussian limit) | 1.98008 | 3.1732 | 12.562 | 0 |
| 0.35 | Literal nonlinear UD Hessian | 1.90561 | 3.2972 | 13.8812 | 0.105013 |
| 0.35 | Nominal reduced kernel shared by BVK1 and BVK2 | 1.90561 | 3.2972 | 12.5649 | 2.3113e-04 |
| 0.35 | OK | 1.90561 | 3.2972 | 12.562 | 2.8281e-16 |
| 0.5 | Full RPA kernel (BURP-TI Gaussian limit) | 1.94557 | 3.22949 | 10.4949 | 0 |
| 0.5 | Literal nonlinear UD Hessian | 1.86121 | 3.37586 | 11.5538 | 0.100898 |
| 0.5 | Nominal reduced kernel shared by BVK1 and BVK2 | 1.86121 | 3.37586 | 10.5538 | 0.00561315 |
| 0.5 | OK | 1.86121 | 3.37586 | 10.4949 | 5.0778e-16 |

## Interpretation

1. BURP-TI is not a new weak-amplitude approximation: at Gaussian order it is exactly the full diblock RPA kernel.  Its approximation enters through the nonlinear log-density completion and thermodynamic integration.
2. The literal nonlinear Uneyama-Doi continuation and the nominal reduced kernel have the same nonlocal and gradient coefficients but different local curvatures.  Their preferred weak-segregation wavelength is therefore the same, whereas their spinodals differ.
3. BVK1 and BVK2 retain the nominal reduced kernel by construction. Their finite-amplitude differences arise from the adaptive gradient term, not from a hidden RPA fit.
4. The OK model reproduces the RPA spinodal through its local quadratic coefficient, but its published gradient and nonlocal coefficients do not exactly reproduce the RPA minimizing wavevector. Spinodal agreement alone is therefore not sufficient evidence for accurate finite-amplitude profiles or periods.
5. This motivates a manuscript structure with two evidence layers: first the Gaussian response shown here, then finite-amplitude stress-free lamella profiles and periods against SCFT.

## Outputs

- `kernel_samples.csv`: sampled inverse kernels and normalized susceptibilities.
- `kernel_minima.csv`: minimum kernel, predicted weak-segregation period, and spinodal values.
- `diblock_kernel_comparison.svg`: publication-style comparison figure.

Settings: `N=1.0`, `b=1.0`, `kRg_min=0.05`, `kRg_max=5.0`, `sample_count=420`.
