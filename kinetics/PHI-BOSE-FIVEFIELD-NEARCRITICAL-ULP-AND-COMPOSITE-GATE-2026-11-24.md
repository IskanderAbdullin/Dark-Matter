# Five-field massive-Bose nearcritical physical quadrature: positive fixed-node integration, but quantum entropy increment below float64 resolution (2026-11-24)

## Scientific reason for this stage

The four-reaction physical Yukawa network evolves FIVE separate species: χ1,antiχ1,χ2,antiχ2 and one real 86-eV φ, with true net Dirac U1, event-level energy and quantum entropy. Our independent analysis already demonstrated that a primitive 18-node p-Gauss grid resolves the total number M0 but underestimates the massive Bose first-cell susceptibility G00 by a factor ~267.5 when μφ=mφ−1e−16 GeV at T=1GeV. The first scalar energy interval is physically *collision disconnected* under hard P=Ppair=2GeV in this network (its top is the exact χ2-pair support floor 99.9975keV), but a future physical SM/FRW network may make it active. Entropy, susceptibility and mode capacity still require physical integration even when collision source is exactly zero.

## Constructive implementation

Added optional phi_IR_composite=True to FullDiracPhysicalMoments, with False as backward-compatible default. In this mode, only the first Bose cell replaces the ordinary physical-p n18 Gauss by a fixed geometry-only 168-node rule: one [0,1e−16GeV] subinterval and 20 logarithmic intervals from pmin=1e−16 to physical pmax, each with positive n8 Gauss weights. The same static physical species DOS grid and Bernstein-P2 basis is used for **M(η), S(η), G(η) and independent FD/BE quantum KL**. It is not a μ-dependent adaptive integration silently omitting basis-weight derivatives. All on-shell χ and scalar event vertices, positive W, hard-P energy supports, true Dirac U1 and shared energy ledger are unchanged.

Independent tests validate:
1. For μφ=mφ−1e−16GeV, new actual five-field phi G00 vs separately transformed physical u-integral is within 2e−7 fractional error, whereas original n18 G00/reference <0.004. Full first moment M0/reference within 3e−10, all five species' integrated p²dp capacities exact to <5e−13.
2. Real 1↔2/2↔2 four-channel microscopic event rates and H match between old/new grids, source in the first strictly inactive scalar cell is mathematically zero.
3. For a generic off-root perturbation, integrated five-species FD/BE quantum KL plus η_new·full residual still matches entropy difference and dt H in the physical DOS discretization.
4. A nearcritical trial at all species T=1GeV, µχ=0, µφ=mφ−1e−16GeV is processed by the fully coupled five-field positive-Bose implicit model, **but** its absolute entropy and charge/energy changes are at the machine precision floor. It does not prove positive entropy production numerically.

See source phi_bose_five_dirac_full_yukawa_implicit.py and targeted test test_phi_bose_five_dirac_ir_composite.py; independently verified IR quadrature phi_bose_critical_ir_susceptibility.py and 50-digit test.

## Genuine negative result: strict global entropy/source-flux tests are physically unresolvable in float64

Initial new nearcritical test failed focused GitHub Scientific gate #45424310325, with 23 successful tests and one rejected strong step. Added actual error diagnostics (#45424582622). They revealed that at nominal dt=5000 GeV^-1:
- physical H_new=**8.428381326601834e−25 GeV4**, hence dt H_new≈**4.2141906633e−21 GeV3**;
- directly subtracted full-entropy certificate ΔS−dt H_new≈**−4.214190663300917e−21 GeV3** (ΔS numerically rounds to ZERO);
- source-scale normalization ~4.9018430903827645e−14 and physical moment root max absolute residual ~9.738050170253925e−19 GeV3;
- when dt=1250, reported full energy defect **+2.7755575615628914e−17 GeV4**, exactly one float64 ULP of a representative total, while original acceptance demanded a much smaller ~1.225e−22 GeV4 tolerance;
- all thirteen attempted dt halvings were rejected, alternating a spurious negative quantum entropy gap and one-ULP energy difference. This was a GLOBAL float64 subtraction problem, distinct from the earlier local Newton DOF 74 <ULP problem. It is not evidence of physical entropy violation, negative Bose occupation or a positive root no-go.

The acceptance logic now computes the actually integrated full before/after physical N/Q/S and adds conservative explicit summation-error budgets:
  entropy_floor = 32 eps_float64 (|S_old|+|S_new|),
  energy_floor = 32 eps_float64 (Σ|Q_old,s|+Σ|Q_new,s|),
  U1_floor = 32 eps_float64 (Σ|CHARGE_s N_old,s|+Σ|CHARGE_s N_new,s|).

It preserves original source-based acceptance tolerance 1e−8*scale, then checks each *raw* observed defect against the maximum of the actual numerical sum-ULP floor and the original source tolerance. It reports all raw defects, each physical roundoff bound, precision_limited_global_invariants, entropy_gain_sign_resolved, entropy_certificate_sign_resolved.

**Important distinction:** a positive *exact* fully implicit FD/BE moment root obeys the physical quantum Bregman entropy identity ΔS−dtH_new=ΣD(old||new)≥0. A floating approximate root has an extra η_new·moment residual. If actual entropy increment or KL certificate is smaller than the physical global double-precision subtraction floor, the software must report **SIGN UNRESOLVED / numerical precision limited**, not claim positive entropy generation. The nearcritical test explicitly asserts this label and uses independent KL/H source checks. Strict scientifically resolved entropy/energy tests remain in normal, nonthermal finite-time P2/P8/P16 tests.

## Reproducibility and remaining gates

New source commits:
- 5cd6160c22ccc7add3761e0df15018fb332e015a — optional 168-node actual physical five-field Bose cell.
- eae4272b18114adf91b7931555e6e23a13e50d0a — independent full physical G/N/capacity/on-shell source/KL/nearcritical tests.
- 0d9ac210d1f8046690501ff2c6ab8caab0625e3b — focused test CI.
- d9f78a5031ed9ee29428ee15575ae821bc1938e6 — quantify genuine negative certificate/energy defect.
- 97c620203e0c2cc695de3ee0283b51c5222a252d — add documented global ULP sum floor.
- 7741d3f64931374784ffb81c92b73d434a92f5c8 — assert unresolved sign explicitly.
Focused five-field CI #45425105944 and #45425264177 **passed**, including the optional nearcritical composite and full physical normal cutoff regression. Full Scientific verification across Python 3.11/3.13 is separately tracked.

**Scope:** only finite vacuum-mass on-shell four-channel Yukawa network and hard P2 first Bose collision-inaccessible cell, and explicitly selected subcritical chemical potentials; not general self-consistent SM+FRW, real 86-eV condensate, Ω or CLASS-ready. The static 168-node rule is only shown accurate across six δ in [1e−18,mφ]GeV, not arbitrarily critical δ or a delta-mode.

**NEXT_ACTION:** independently reconstruct nearcritical discrete entropy change with a compensated or high-precision physical sum and the exact quantum KL residual-work identity, and quantify at which dt/chemical gap the sign becomes measurable. Then scan actual five-field 3×energy/time/angular/parent cutoff meshes. No condensate inference from the G00 divergence.

**Controlled IR-only implementation erratum:** the first optional implementation built n168 Gauss in ALL φ energy cells because rectangular Species arrays share one node count. This confounded first-IR versus global Gauss refinement. In [code correction](phi_bose_five_dirac_full_yukawa_implicit.py) commit 064995d1, other Bose cells are now exactly the ORIGINAL nq Gauss nodes/Bernstein weights, with additional rectangular slots carrying precisely zero physical DOS weights. ONLY the first frozen IR cell uses the fixed positive n168 geometric p quadrature. Targeted test asserts byte-identical non-first physical cell quadrature and on-shell sources unchanged; focused CI #45426158199 passed, including the subsequent explicit equal-weight test. This correction is a scientifically necessary controlled comparison, not a new physical phase source.
