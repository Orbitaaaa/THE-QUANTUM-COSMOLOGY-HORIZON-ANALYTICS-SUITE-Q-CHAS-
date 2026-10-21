# Numerical methods

All code lives in `CODE/babyuniverse/`; every figure and data file is produced by `CODE/scripts/run_all.py`.
Dependencies: numpy, scipy, sympy, mpmath, matplotlib (`CODE/requirements.txt`).

| task | method | accuracy control | file |
|---|---|---|---|
| symbolic identities | sympy simplification to 0; exact rational arithmetic where sympy cannot simplify (Hamilton–Jacobi check) | identities exact; mpmath quadrature 30 digits for $J(s)$ | `symbolic.py` |
| turning points | roots of the sextic $-\tfrac{L^4}4R^6+\kappa^2R^4+GM(\chi^2-\eta^2-\kappa^2)R^3-G^2M^2=0$, refined by Brent | $\lvert V(R_t)\rvert<10^{-10}$ (tested) | `model.py` |
| critical mass | closed form from $z_*^3=\tfrac12(-\gamma+\sqrt{\gamma^2+8})$ | double root verified | `model.py` |
| itinerary | chart-free classification of monotonic legs cut at horizons and flip radii | exact (no integration) | `branches.py` |
| global trajectories | $\ddot R=-V'/2$ plus null-coordinate equations; DOP853, rtol $10^{-11}$; switch between two algebraically equivalent forms at each step to avoid $0/0$ at horizons; dense output resampled uniformly | finite-difference checks of $u\cdot u=-1$ and $n^r=A$ | `geometry.py` |
| Penrose diagrams | arctan compactification (Kruskal), global conformal coordinates (dS); shading from the integrated wall only | — | `plotting.py` |
| canonical exponent | adaptive quadrature of $P_E$, with the analytic or principal branch; breakpoints at flip radii | $\epsilon_{\rm rel}=10^{-12}$; M = 0 vs closed form $<10^{-14}$ | `tunnelling.py` |
| Euclidean wall | substitution $R=R_a+(R_b-R_a)\sin^2(s/2)$; cancellation-free endpoint distances; exact endpoint limit | Smarr residual $<3\times10^{-8}$ | `tunnelling.py` |
| interior Euclidean volume | Stokes with the 1-form $(r^3-r_h^3)/3\,dt_E$ (regular at the de Sitter horizon) plus axis term | first law to 8 digits | `tunnelling.py` |
| multiplicity (pseudomanifold) | degree theory: winding number of the boundary loop around grid points; window enclosing the whole wall | grid refinement (121→481) leaves the range unchanged | `tunnelling.py` |
| first law | central finite differences, $h=10^{-5}$ relative | agreement to 8 digits | `tunnelling.first_law_check` |
| $M=0$ closed forms | cancellation-free expressions for $1\pm p$, $1\pm q$ | relative $10^{-6}$ at $\kappa/\chi=3.6\times10^{-3}$ (tested) | `tunnelling._stable_cap_cosines` |
| thick-wall initial data | $2rr''=1-r'^2-8\pi Gr^2\rho$ from a regular centre, series start; refined grid across the wall | first-order convergence $O(w)$ verified | `thickwall.py` |
| CdL bounces | bisection on $t=\ln(1/\lvert\phi(0)-v\rvert)$ with a linearised Bessel start; event detection for over/undershoot; analytic dS cap completion | constraint violation $<3\times10^{-10}$; $O(w)$ convergence to thin wall | `thickwall.py` |
| reduced Hamiltonian | damped Newton in $M$ with exact $\partial P/\partial M$, sector-constrained | $dH/dP$ vs $dR/dt$ to $6\times10^{-10}$ | `quantum.py` |

Reproduction: `cd genesis/CODE/scripts && python run_all.py` (about 8 min on one core), then
`cd .. && python -m pytest tests -q`.
