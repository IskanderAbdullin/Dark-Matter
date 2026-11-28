# Decay equilibrates low moments but does not equilibrate the soft scalar tail

This is an independent local/adiabatic/kinematic diagnostic of the project benchmark, not a kinetic evolution or a prediction of a condensate or relic abundance. It uses the same scalar Yukawa channel as the [upstream draft](https://github.com/IskanderAbdullin/Dark-Matter/blob/90c141722376dab64e4cd10422d545474817787c/kinetics/PHI-BOSE-PORTAL-INITIAL-GAIN-GATE-2026-11-24.md): m1=.1 GeV, m2=.02 GeV, mphi=86 eV, y=.001. No empty-sector injection, Yukawa 2-to-2, scalar self-collisions, real vector dynamics or finite-density propagators enter these controls.

## The decay-only low-moment equilibration gate

Use exact two-body kinematics, time-dilated spin-averaged rest width 2.7501974e-9 GeV, FD fermions and finite-mass BE scalar. The event vector for number and energy/T of (1, 2, phi) is (-1,-E1/T,1,E2/T,1,Ephi/T). The local response at the common neutral equilibrium is a Gram matrix with positive weights given by the equilibrium forward event measure. It has three null vectors: dark-fermion count (n1+n2), decay scalar count (n1+nphi), and total energy. The rates refer to both separately conserved Dirac charge sectors at zero asymmetry. The sum n1+n2 need not be conserved when pair channels are added. Compare the Gram result with the derivative of the actual nonlinear FD/BE event flux and check equilibrium, entropy production, positive projected spectrum, MB Bessel normalization and phase-space refinement.

The slowest nonzero projected mode/H is 3.05e8 at 0.68825 GeV and 2.17e11 at 0.01 GeV, under a common-T SM+dark Friedmann normalization. These clocks test only the six-dimensional thermal moment manifold; they are not a spectral gap for the full momentum collision operator.

## Where a normal thermal closure fails under isolated adiabatic dilution

Condition on rapid decay enforcement of eta1=eta2+eta_phi and common temperature. Keep a^3(n1+n2), a^3(n1+nphi), and a^3S constant. Solve for the dark temperature and both independent chemical potentials with the normal constraint eta_phi<=m_phi/T. With T0=.68825 GeV, eta2,eta_phi=0, ordinary BE saturates already near a=1.0000214. Starting below saturation changes this: the initial gap (mu_phi-m_phi)/T=-.35 gives a=11.587, whereas -2 never hits it in the test through a=40. An extra scalar count beyond the normal distribution's support is only a bookkeeping diagnostic; it neither proves dynamically produced zero modes nor licenses assigning a condensate. SM heat, total-particle-changing dark reactions and any nonthermal evolution are excluded here.

## Momentum-resolved refutation of a uniform fast-decay proxy

At fixed equilibrium fermion distributions, the derivative of the scalar collision term in a scalar bin is proportional to the difference f2-f1 integrated over the allowed parent energies. Separate the exponential tail analytically to avoid mistaking underflow for an exactly zero rate. At p_phi=0 the apparent inverse momentum cancels the width of the allowed parent-energy interval, but the parent must carry about 55.8 TeV. For p_phi/T below the first decay-only damping/H crossing, the calculated local clock is slower than H even while the projected heat/number modes are much faster. The first crossings are 1.66e-4 at 0.68825 GeV, 7.26e-3 at 0.10 GeV and 0.678 at 0.01 GeV in this conditional common-bath comparison. These are local scalar momenta, not a decoupling surface or a trajectory. The extreme IR cannot be swept into an instantaneous normal/condensed equilibrium via this decay mechanism alone. Other interactions must be assessed before any conclusion about the soft tail.

## Reproduction
Use the authors' central [SM EOS table](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat) without redistributing it. The independent adiabatic count test does not require the SM table.
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/three_species_decay_onsager_reference.py EOS_TABLE
python -W error kinetics/decay_only_adiabatic_normal_branch_reference.py
python -W error kinetics/decay_scalar_ir_relaxation_reference.py EOS_TABLE
```
