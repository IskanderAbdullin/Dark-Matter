# Stage 2026-11-23 — exact-angle FD×BE primitive and independent Yukawa source normalization

- **ID:** DM-KINETICS-20261123-THERMAL-EXACT-WEAK-PRIMITIVE.
- **Parent HEAD:** [bbf5fb86](https://github.com/IskanderAbdullin/Dark-Matter/commit/bbf5fb8649fc7de0a3690dbbaa9e3969e07307b6) / existing research branch / draft PR #1.
- **Analytical result:** F/R=exp[E1(1/T0−1/T1)] at T2=Tφ=T0, so Bose/FD inverse thermal angle integral is an explicit log/expm1 endpoint primitive; total φ net source requires only parent 1D integral.
- **Independent evidence:** [standalone code](../../../../kinetics/phi_bose_thermal_exact_weak_reference.py), [3 focused tests](../../../../tests/test_phi_bose_thermal_exact_weak_reference.py), [report](../../../../kinetics/PHI-BOSE-THERMAL-EXACT-WEAK-REFERENCE-AUDIT.md), [JSON](../../../../kinetics/results/phi-bose-thermal-exact-weak-reference-2026-11-23.json). Tests 3 passed 1.32s after correcting independent adaptive quad abs error goal; no external BEST/time evolution asserted.
- **Numerical scope:** pmax2/8/16 analytic-source Rφ −1.295453955495/−3.212015339411/−3.232104327009×10−11; independently normalized positive angle/radial weak source to ≈2×10⁻⁹ relative or better. Physical weak parent/daughter energy sum residual <4×10⁻²⁶, H>0. Qφ and mean Eφ cutoff sensitivity measured.
- **Hard limitations:** fixed thermal FD/BE, no properly weighted finite-capacity strong collision, real Γel/SM/FRW, Bose zero-mode/Ω/CLASS/observations. Public CI for this commit pending.
- **NEXT_ACTION:** derive energy-resolved quantum-positive finite-volume strong operator from this weak reference; maintain same physical source normalization, exact shell and correct graph invariant rank, before time evolution/FRW.

## Extension: independent exact physical Qφ weak-energy moment

Real Li₂ analytic daughter moment yields 1D pmax2/8/16 Qφ −2.670878687857882e−12/−1.3941860432430647e−11/−1.4382367979860587e−11, independently matching 2D physically weighted on-shell event Qφ ≤2e−13 relative. Two extra local focused tests, 5/5 passed (0.58s). Source+JSON/report extended; public CI pending. Not physical finite-capacity strong solver/FRW/Ω/BEC.
