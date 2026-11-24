# Five-field physical Dirac Yukawa nonlinear strong: first time/spectral/UV cutoff and pair radial convergence (2026-11-24 UTC)

## Purpose

An actual *finite-time* fully coupled 5-population χ1,antiχ1,χ2,antiχ2,real φ four-channel Yukawa solver with physically correct massive particle DOS, same positive on-shell event weights, true U1 and total energy, and conditional positive Bose/FD + quantum KL implicit entropy certification was recently established. Green single-step tests DO NOT establish a converged physical continuum nonlinear time evolution or early-Universe freezeout. This stage independently checks 3 scalar spectral grids, 3 time resolutions, 3 hard parent momentum UV cutoffs and source-specific pair radial quadrature.

Fixed illustrative *vacuum* Yukawa masses m1=.1GeV,m2=.02GeV,mφ=86eV, y=.001. The 5-field initial occupations for this test follow the solver's nonthermal/unequal-temperature specified thermal trial, not a realistic SM/FRW history. Event geometry includes both decays and both χχbar→2φ annihilations. Hard cutoff P applies to parent decays and incoming pair momenta. This is NOT yet physical UV extrapolation.

## Real nonlinear full five-field time convergence, independent of weak source

At hard P=2GeV, total nonlinear physical evolution to t=20,000GeV^-1, coupled 4 channels, Fermi/Bose entropy variables positive, three distinct original Bose spectral meshes Nφ=8,12,16, all using nrad4,pair-nrelative8,nCM32,naz4 and full implicit Euler:

| Bose original Nφ | ΔNφ at Δt=10000 (GeV³) | ΔNφ at Δt=5000 | ΔNφ at Δt=2500 | successive time-difference ratio |
|---:|---:|---:|---:|---:|
| 8 | −1.5359482924415e−7 | −1.5376978099169e−7 | −1.5385816552749e−7 | **1.97944** |
| 12 | −1.5352677097680e−7 | −1.5371168797518e−7 | −1.5380524318043e−7 | **1.97655** |
| 16 | −1.5354135002310e−7 | −1.5372435621663e−7 | −1.5381688922567e−7 | **1.97774** |

Observed difference ratios near 2 match the expected **first order** backward-Euler temporal truncation (not evidence of second-order accuracy). All numerical trajectories passed full 4-channel Bose positivity, true U1 and total E ledger, microscopic H and actual finite-step quantum entropy checks, no adaptive rejections in this reference. Per-final-run absolute net U1 differences are O(1e−17)GeV³ and total E differences O(1e−17)GeV⁴, at roundoff of large total state sums.

**Important spatial negative:** At the finest fixed Δt=2500, Nφ8 and Nφ16 ΔNφ differ ~0.0268%, while Nφ12 and Nφ16 differ ~0.00757%. This sequence is **nonmonotonic** with spectral N, and 3 grids do NOT establish an asymptotic spectral convergence rate. We can quote a numerical spread, but do not extrapolate it as a proven independent continuum error bound.

## Hard parent UV cutoff is materially important

For original Nφ12, same t=20000GeV^-1 and Δt=5000GeV^-1:

| P=pair incoming and decay parent cutoff (GeV) | ΔNφ from full 5-field evolution (GeV³) | initial quantum H (GeV⁴) | observed ΔS (GeV³) |
|---:|---:|---:|---:|
| 2 | −1.5371168797518e−7 | 9.7376929959e−13 | +1.9346234281e−8 |
| 8 | −3.9216250392882e−7 | 4.9730577514e−12 | +9.7278999522e−8 |
| 16 | −3.9679229398015e−7 | 5.0931399548e−12 | +1.0009689744e−7 |

The P2 magnitude is **only 38.7%** of P16, despite all discrete conservation and entropy tests passing. P8 differs from P16 by **1.167%**. These are empirical **cutoff+quadrature sensitivity** measurements, not a pure UV truncation error (the discrete incoming pair p,q quadrature also changes with P). The finite-P2 strong test should NEVER be presented as UV-converged thermal physical freezeout.

## Distinct channel radial convergence, and graded quadrature negative

At fixed P16 with nonthermal initial state, pair CM and incoming-relative angular orders held at (nrelative8,nCM32,naz4), increase the ordinary **single-interval** incoming radial Gauss count nrad:

| nrad | χ1+antiχ1→2φ signed event rate (GeV⁴) | χ2+antiχ2→2φ signed event rate (GeV⁴) |
|---:|---:|---:|
| 4 | −1.264253570520e−17 | +8.454062124304e−17 |
| 6 | −1.287991955309e−17 | +8.457093190987e−17 |
| 8 | −1.306543041994e−17 | +8.302357141400e−17 |
| 12 | −1.312448955774e−17 | +8.297828508393e−17 |
| 16 | −1.308993078516e−17 | +8.279782200607e−17 |

While the **total** entropy H is dominated by much larger 1↔2 weak rates and remains ∼5.093138e−12GeV⁴, individual pair channels still vary at the **subpercent to several-percent level**. True U1 and E defects are machine-small even for *underresolved* channels; these algebraic source tests cannot certify convergence.

To independently probe a positive physically normalized radial composite quadrature, added optional **fixed, occupation-independent** incoming momentum breakpoints [0,.25,1,4,P]GeV to the pair on-shell event generator, each with nper=2/3/4 Gaussian nodes. Both incoming p²dp physical capacities are exactly integrated by nper≥2. It preserves on-shell energies, positive event weights and full nonlinear U1/E/entropy/KL at a P16 test step. However it does NOT automatically outperform single-interval Gauss at equal node count:

| Graded radial nodes per segment (total nodes) | χ1 pair rate (GeV⁴) | χ2 pair rate (GeV⁴) |
|---:|---:|---:|
| 2 (8) | −1.154923593738e−17 | +7.710318776209e−17 |
| 3 (12) | −1.288076894351e−17 | +8.284811548553e−17 |
| 4 (16) | −1.311896691418e−17 | +8.291779662266e−17 |

The **8-node graded pair rate is roughly 7–12% below** the refined results, i.e. worse than the ordinary 8-node rule for this example; nper3/4 approach each other, but nper4 and ordinary n16 still differ ~0.15–0.22%. This is an explicit numerical negative, not a failed physical positivity or a BEC effect. Segment boundaries chosen for the benchmark are **not universal**, and no arbitrary nonthermal/high-P collision moment convergence has yet been proved.

## Reproducibility / status / next

- [Full 5-field finite time and cutoff/time audit source](phi_bose_five_dirac_strong_convergence_audit.py)
- [All physical tests](../tests/test_phi_bose_five_dirac_strong_convergence_audit.py)
- [Optional exact-capacity graded on-shell pair radial integration](phi_bose_pair_physical_onshell_quantum.py)
- [Physical 4-channel weak source](phi_bose_five_dirac_full_yukawa_weak.py), [fully coupled implicit](phi_bose_five_dirac_full_yukawa_implicit.py)
- [Machine numerical evidence](results/phi-bose-fivefield-strong-time-uv-radial-2026-11-24.json).

Focused on-branch GitHub CI succeeded for 3-grid/3-time audit #45425641136 and the full P16 graded/radial physics after one **test-only** missing numpy import was corrected (#45427815208). Prior run #45427342217 failed 1/34 due the test harness NameError (numpy missing), NOT a collision physics error; that negative and fix are retained. Additional after-refinement run and full two-Python Scientific verification must be checked separately.

**NEXT_ACTION:** decouple the physical hard-cutoff error from radial and relative/CM angular quadrature by independent high-order pair on-shell source, preferably an alternative positive importance transform or adaptive reference for each channel. Then test simultaneous full five-field finite-time at ≥3 energy grids, ≥3 time grids, ≥3 parent/angle grids and qmin/UV. Only after that derive actual χ1/χ2+SM in-medium scattering and a closed FRW physical fφ(q). This audit is NOT Ω86eV, BEC, CLASS-ready or an observational fit.
