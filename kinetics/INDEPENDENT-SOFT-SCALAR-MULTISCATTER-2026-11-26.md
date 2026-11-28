# Repeated below-pole elastic transport of soft scalars

This checks whether repeated φχ2 scattering can undo the bottleneck seen in the [single-scattering rest test](INDEPENDENT-SCALAR-REST-ELASTIC-HANDOFF-2026-11-26.md). It linearizes about prescribed FD fermion baths and a normal BE scalar distribution, at fixed temperature, including both fermion charge sectors. It is not a self-consistent expanding background, a coherent zero-mode transport calculation, or a prediction of a scalar spectrum or relic abundance.

## Positive finite-volume projection

Use pair invariant energy u, total lab energy L, incoming and outgoing φ angles along the pair boost, and relative azimuth. The measure times inverse s is proportional to |P| p_star²/u. Both bin boundaries are applied to the exact lab scalar energies. The ordered equilibrium traffic is weighted by f_phi f_chi (1+f_phi') (1−f_chi'); detailed balance implies equality across every forward/reverse bin pair. Use a symmetric s cut below the real χ1 pole and a pair-energy tail cutoff, so this does not count real inverse decay as elastic scattering. The independent rest-angle s/u trace matches the invariant trace. The test deliberately uses different incoming/outgoing angular orders before symmetrization; a linear-angle first draft failed binwise reciprocity in the wider shell. A log-energy angular coordinate restored the sampled convergence of this control.

For the perturbation variable h=δf/[f(1+f)], the weight of a resolved scalar cell is the integral of f(1+f) over that cell. The symmetric elastic conductance has a constant null vector and nonnegative quadratic dissipation. Add the actual local decay clock as a separate positive damping term under the fixed fermion bath, never as an elastic sink. Evaluate first passage to the physical cell boundary where local decay/H first reaches one. This threshold is an absorbing condition **only in the first-passage auxiliary calculation**. The separate reflecting/thermal-boundary checks use continuous local decay within their stated support. Use both the discrete Green's function and Monte Carlo walks with exponential holding times. No inference depends on treating a single-bin eigentime as a full collision spectral gap.

## Controls

At T=.10 GeV for a neutral FD equilibrium fermion bath, the full repeated-transport mean from the lowest finite scalar cell to that threshold is about 9e-5 Hubble times. The narrower-first-cell and alternate support refinements leave this order unchanged. A one-step projection from the lowest cell has an extremely small probability; the earlier p=0 calculation found 0.00960 direct events per H. Even imposing adjacent-bin moves only gives a mean less than one Hubble time on the sampled grids, but its precise value is not a mesh invariant. A 40,000-walk control reproduces the root calculation within sampling uncertainty. At T=.01 GeV, direct handoff is absent on the included thermal support, while sequential handoff gives about 3e-4 Hubble times in this fixed FD bath.

These are **small-perturbation** benchmarks with a **populated** thermal fermion bath. They do not override either exact scalar-number conservation of elastic scattering or the very slow direct decay of the scalar rest mode. Eventual decay and artificial reservoir escape fractions vary with the outer cutoff, especially at low temperature; do not claim an outer cutoff-independent abundance from them. Repeated steps are sufficient to reach the fast-decay window in this frozen test. This does not say that a realistic population reaches it when the bath is depleted or nonthermal.

## Lower populated fermion density changes the conclusion

Keep the shapes FD, shift both fermion fugacities together so the neutral BE scalar still satisfies decay detailed balance, normalize H to the same reference background, and recompute both the elastic graph and the decay-fast threshold. At T=.10 GeV the resulting mean first-passage times in Hubble units are .00345, .189, 1.41 and 10.5 for common fermion μ/T=-4, -8, -10 and -12, respectively. The crossing to slow transport in this sample depends on the assumed fermion population: it is not a cosmological decoupling temperature. Real depleted baths can also be nonthermal or asymmetric, neither of which is modeled here.

### Numerical scale and controls

For the populated .10 GeV benchmark, the 40,000-walk mean is 8.959e-5 Hubble times, standard error 2.65e-7, against the deterministic 8.939e-5. Changing the first-bin cutoff, outer support or below-pole cut changes that first-passage mean by less than 1.3% across the recorded selection. The fraction absorbed by decay rather than the artificial reservoir is less stable: widening the outer domain in the .10 GeV example changes it from .470 to .641. The graph has not established the true loss fraction. Its significance is the contrast between single-event and repeated soft transport, and that contrast can disappear in a depleted fermion bath.

## Reproduce
Use the authors' central [SM EOS table](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat), not redistributed here. The results record the bin counts, quadrature and cutoff choices. Run in the cloud, e.g.:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/soft_scalar_elastic_graph_reference.py EOS_TABLE --T .1 --bins 16 --orders 7 5 4 16
python -W error kinetics/soft_scalar_elastic_graph_reference.py EOS_TABLE --T .01 --bins 16 --orders 7 5 4 16 --limit 3 --tail 32
python -W error kinetics/soft_scalar_elastic_graph_reference.py EOS_TABLE --T .1 --bins 14 --orders 6 4 4 14 --eta-fermion -8
```
