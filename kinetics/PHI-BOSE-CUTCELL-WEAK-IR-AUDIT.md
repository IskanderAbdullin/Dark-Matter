# Quantified physical Bose IR source missed by sparse 12×12 Yukawa angular quadrature (2026-11-23)

## Falsifiable question and decisive result

The [sparse-angle negative audit](PHI-BOSE-SPARSE-ANGLE-FALSE-ZERO-AUDIT.md) established that the first scalar energy cell in the previous discontinuous-P1 entropy solver at N16, parent cutoff P=2 GeV, is **physically on-shell accessible**, despite receiving zero 12×12 Gauss event deposition. We now calculate the **actual signed first-cell number and scalar-energy sources**, without any daughter-cell sampling alias, and independently integrate the entire physical on-shell angular shell.

Use the **same** illustrative physical tree-level Yukawa vertex, masses and thermal state as the earlier [DG-P1 benchmark](PHI-BOSE-CELL-ENTROPY-AUDIT.md): m1=.1 GeV, m2=.02 GeV, mφ=86e−9 GeV, y=.001, g1=4, T1=.9 GeV, T2=Tφ=1 GeV, μ1=μ2=μφ=0, parent p∈[0,2] GeV. This is a fixed thermal *weak source*, not a nonlinear strong time simulation.

The scalar first-cell energy upper edge is Etop=0.00860359519316437 GeV. Exact Lorentz inversion requires p≥0.269991073952339 GeV for any daughter to enter this cell. For each parent momentum p, integrate only the **exact positive-length cut interval**
\[
k\in[k_{\min}(p),\min(k_{\max}(p),E_{\rm top})].
\]
The physical scalar net flux is \(J=f_1(1-f_2)(1+b_\phi)-(1-f_1)f_2b_\phi\). Multiply by the physical parent density-of-states \(g_1p^2dp/(2\pi^2)\), time-dilated tree-level width \(\Gamma m_1/E_1\), and rest-angle Jacobian \(dk/(2p q_\phi^*/m_1)\). No arbitrary daughter capacity or scalar bin point occupancy enters. Compute Qφ with an additional **physical k** factor, not an energy-bin center.

## Independently evaluated adaptive double integrals

| Quantity | Continuous cut-cell physical source | Previous sparse 12×12 event source |
|---|---:|---:|
| First-cell net scalar number Rφ (GeV⁴) | **−2.76214936045×10⁻¹²** | **0** |
| First-cell scalar energy Qφ (GeV⁵) | **−1.22779976195×10⁻¹⁴** | **0** |
| Total scalar net number Rφ (GeV⁴) | **−1.29545395549×10⁻¹¹** | −1.16493012863×10⁻¹¹ |
| First-cell share of physical net scalar number | **21.32%** | 0% |

Thus the old 12×12 quadrature **underestimates the magnitude of the total net scalar destruction by 10.08%** relative to this independent continuous weak reference. This is not just a formally nonempty first cell: the missed physical contribution is substantial. Both integrals are *negative* because the parent is colder than the daughters and inverse absorption dominates at these trial thermal parameters.

The first-cell number is an integral of the true on-shell \(F-R\) over the full continuous cut-cell region. The prior graph's first-cell zero and zero inverse-fugacity tangent remain true **of the sparse sampled graph only**; they are not physical properties of the Yukawa operator. The previous claim of strong solver physical IR convergence at N16 would be unsupported, even if its discrete charge, energy, entropy and backward-Euler residual all passed.

A separate high-order **direct Gauss parent×angle** cross-check (independent of the cut-cell boundary transformation) gives total Rφ ≈−1.2936843×10⁻¹¹ at 48×48, ≈−1.2954453×10⁻¹¹ at 96×96, and ≈−1.2954540×10⁻¹¹ at 192×192. These approach the nested adaptive exact-angular value −1.29545395549×10⁻¹¹, while the original 12×12 is visibly biased. For the first cell, 192×192 Gauss gives approximately −2.76696×10⁻¹², with slower convergence due to the hard cut; the exact cut-cell integral avoids that discontinuity.

## Scientific decision and limitations

**Refuted:** a first-cell zero in the N16/pmax2/12×12 DG-P1 solver reflects a physical kinematic absence or harmless small missed contribution. It is a **quadrature artifact with an ~21% share of the total true number rate**, in the tested state.

**Supported:** the finite-support rapidity theorem, exact angular cut-cell transformation, sign of reverse net source, physical weak number and energy moment at this thermal point, and convergence of independent high-order Gauss total rate to adaptive double integral.

**Not established:** daughter-bin finite-capacity strong operator, nonthermal time evolution, full three-energy cut-cell ledger, convergence under scalar N/parent P/SM+FRW, a condensate zero mode, physical Ω or CLASS. The 10.08% discrepancy applies **only to the stated benchmark**; it is not a universal error for 12×12 Gauss.

[Reproducible source](phi_bose_cutcell_weak.py) and [regression tests](../tests/test_phi_bose_cutcell_weak.py). Commands:

\`\`\`bash
OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_cutcell_weak
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_cutcell_weak.py
\`\`\`

Local independent SciPy adaptive and NumPy high-order Gauss calculations performed. New public CI must be checked separately. **Next:** use exact cut-cell integration to assemble all physical scalar energy cells, prove total R and Q ledger across bins, and then derive shared physical phase-space capacities before nonlinear Pauli/Bose time. Do not use sparse sampled first-cell zero as evidence of freeze-out or Bose condensation.
