# A shared tadpole suppresses the abundance in more than one momentum shell

The [earlier thin-shell control](INDEPENDENT-EARLY-POLE-SHELL-HARTREE-FEEDBACK-2026-11-26.md)
assigns one occupation by construction.
Here independently evolve all bins of a logarithmic
scalar grid spanning $p/T=10^{-9}$ to 30 in the finer
arm. Other parts of that closure remain deliberately
simple: an externally held FD parent, empty fermion
daughter, no energy closure, no elastic mode
redistribution, and the same prescribed radiation
clock and scalar Yukawa rest-width normalization.
All populated scalar bins contribute to **one**
self-consistent positive Hartree pole at each step.

Let each bin carry uniform occupation and integrate
in $z_i=\ln(1+f_i)$. This avoids overflow in the
stimulated control. For the no-self comparison stop
when the chosen nominal full number is first reached,
rather than evolving the unphysical infinite
parent reservoir indefinitely. The feedback examples
instead reach $x=.05$ below that nominal count.

| Held parent and daughter scenario | Final total yield | Yield in $p/T<10^{-5}$ | $m_\phi^{\rm pole}/T$ |
|---|---:|---:|---:|
| No self feedback; stop early at $x=.00709$ | .005057 | .005057 | $7.1\times10^{-12}$ |
| $\lambda=7.65\times10^{-8}$, massless fermion, $x=.05$ | .00080224 | .00080223 | .0011285 |
| $\lambda=1.085\times10^{-6}$, massless fermion, $x=.05$ | .00005657 | .00005657 | .0011285 |
| Small quartic, fermion pole $b/T=.003$, $x=.05$ | .00097917 | .00097917 | .0012060 |

The narrower early population remains cold by
the **momentum** cut here, while building much less
count. Its stress cannot be inferred from momentum
alone during a self-mass regime; see the separate
[energy and pressure check](INDEPENDENT-EARLY-POLE-HARTREE-ENERGY-PRESSURE-2026-11-26.md).
The 150- and 360-bin arms also change the lower
cutoff. Integral yields in tested bands match to
about $4.1\times10^{-7}$ relative or better at
their respective final comparison points; mean
momentum agrees to within .6%. This does not
prove that the extreme grid-peak occupation is
resolved, or that leading Hartree alone controls
higher-order dynamics in such a population.

No parameter or spectrum is retuned to the
abundance after feedback. The wider grid removes
the imposed single shape, but it does not resolve
backreaction on the parent and daughter, actual
thermal masses or widths, nonlinear field
evolution, entropy production or cosmological
linear perturbations.

~~~sh
python -W error early_pole_multiband_hartree_reference.py
~~~

[Compact](results/early-pole-multiband-hartree-compact-2026-11-26.json.gz),
[fine part 1](results/early-pole-multiband-hartree-fine1-2026-11-26.json.gz)
and [fine part 2](results/early-pole-multiband-hartree-fine2-2026-11-26.json.gz)
retain the per-bin number spectra and the sampled clock histories.
