# An exact endpoint energy identity for the multibin scalar Hartree state

At fixed scale, write the gap as u−u_b = λ Σ(n/E)/2, with u=m² and E=sqrt(q²+u). The energy is R=Σ(n E)−(u−u_b)²/(2λ). Let subscripts 0 and 1 mark the two endpoints, let Δ indicate their difference, and let a bar denote their arithmetic mean. Define

```text
A_i = E0_i + E1_i
h_i = (1/E0_i + 1/E1_i)/2
D   = 1 + λ Σ[nbar_i/(E0_i E1_i A_i)]/2
B   = Δu Σ[(n1_i/E1_i − n0_i/E0_i)/A_i²]/4
g_i = A_i/2 + λ B h_i/(2D).
```

Subtracting the two gap equations gives Δu = λ Σ(h_i Δn_i)/(2D). The product and square differences in R then give ΔR = Σ(Ebar_i Δn_i) + Δu B = Σ(g_i Δn_i). Fixed spectators must be included in both gap solutions and in D and B even though their Δn terms vanish. At a repeated endpoint, g is the ordinary Hartree energy gradient. The implementation computes the difference of n/E without subtracting two large nearly equal reciprocals.

An independent finite-state cloud fixture checks 320 positive pairs on two scalar supports, with massless or massive bare inputs and both evolving and fixed spectator populations. The maximum discrete-gradient energy mismatch is 2.04e−20 absolute; plain average and midpoint gradients reach 3.97e−7 and 1.91e−7 for these deliberately finite moves. The secant energies are ordered on all tested supports. A separately integrated line-average entropy gradient agrees with the direct entropy difference to 1.22e−24 absolute. These errors refer to the fixture normalization, not a collision rate.

This identities-and-ledgers check is not a collision operator. If a future stencil changes its remap to conserve this gradient, its event counts, angle/momentum support, positivity and entropy must be verified separately. No evolution or FRW energy statement follows just from the fixed-scale identity.

From `kinetics/`:

```sh
OPENBLAS_NUM_THREADS=1 python early_pole_scalar_hartree_discrete_gradient_reference.py
```

Raw: `results/early-pole-scalar-hartree-discrete-gradient-2026-11-27.json`.
