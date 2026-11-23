# Stage 2026-11-23 — constant one-η Bose cell: necessary energy/equilibrium no-go

- **ID:** DM-KINETICS-20261123-BOSE-ONE-DOF-ENERGY-THERMAL-NEGATIVE.
- **Starting branch/HEAD:** research/kinetic-equilibrium-audit-20261121 / [f0db3e9](https://github.com/IskanderAbdullin/Dark-Matter/commit/f0db3e936d738626d6b9b1defcb2ce3fda20cfb8), PR #1 draft.
- **Previous real CI:** [Scientific #45130112042](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45130112042) **success** on f0db3e9 (predecessor DG-P1 negative stage); do not report it as verification of this new stage.
- **Inputs:** exact common-T=1 GeV, μ=0, χ1=.1, χ2=.02, φ=86e−9 GeV, y=.001, on-shell boosted physical Yukawa parent pmax2,8,16, positive log-Eφ angle.
- **Result:** analytic one-η cell particle-number monotonicity proven, but scalar energy U=ĒM cannot reproduce on-shell Eφ for multiple events. At exact thermal equilibrium, fixed cell η matching physical Bose cell N produces false Rφ, H and U. N8→128 false H declines 2.2763e−13→8.6469e−16, not identically zero. This rejects the proposed exact closure, **not** a physical BEC.
- **Independent code:** [source](../../../../kinetics/phi_bose_one_dof_no_go.py), [3 tests](../../../../tests/test_phi_bose_one_dof_no_go.py), [report](../../../../kinetics/PHI-BOSE-ONE-DOF-NO-GO-AUDIT.md), [compact data](../../../../kinetics/results/phi-bose-one-dof-no-go-2026-11-23.json).
- **Observed local checks:** focused 3/3 pytest passed in 0.10 s in isolated Python container; no full external BEST/FRW experiment. Current public Scientific CI needs independent check after Git push.
- **Next:** positive event-energy-resolved or convex constrained two-moment representation, exact microreversible energy/number/H, positivity under physical time step, graph nullspace, p/q/IR angular controls, before any real 86-eV relic/CLASS assertion.
