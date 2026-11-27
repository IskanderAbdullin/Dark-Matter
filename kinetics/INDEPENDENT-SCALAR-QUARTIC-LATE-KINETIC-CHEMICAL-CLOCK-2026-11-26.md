# Quartic elastic and chemical clocks are distinct on the held late shapes

Use the two illustrative values of $\lambda$ in the
[prior scalar quartic bound](INDEPENDENT-SCALAR-QUARTIC-MULTIPLICATION-BOUND-2026-11-26.md):
$7.65175\times10^{-8}$ and $1.08516\times10^{-6}$.
These are not inferred rates or a new halo fit.
Hold the recorded post-2-MeV scalar distributions
and evaluate separate tree-vacuum quartic processes.

For identical elastic scalars the unstimulated bath
kernel is
$\sigma v_{\rm M\o ller}=\lambda^2\sqrt{1-4m_\phi^2/s}/(64\pi E_1 E_2)$.
This follows from the final-state division by two
and the bath-to-CM Møller ratio. Counting the two
particles per collision cancels the identical
incoming division in the gross per-particle clock.
Fold the kernel over both sampled Bose momenta
and the relative angle. This is not a damping
eigenvalue or proof of local kinetic equilibrium.

For the isolated quartic $2\to4$ channel, retain
the prior all-invariant vacuum ceiling and multiply
it by $(1+F)^4$, where $F$ is the fitted held
distribution's maximum occupation. The ratio of
reverse to forward gross rates is exactly
$\exp(2\mu/T)$ for an instantaneous ideal BE bath
with this temperature and chemical potential, by
$f/(1+f)=\exp[(\mu-E)/T]$ and energy conservation.
It is about .105 to .116 across the three held
states. This need not hold on an arbitrary spectrum.

| Visible anchor | Small $\lambda$: elastic rate / photon clock | Small $\lambda$: chemical gain / photon clock, upper | Large $\lambda$: chemical gain / photon clock, upper |
|---:|---:|---:|---:|
| 2 MeV | 314 | $5.36\times10^{-7}$ | .0217 |
| 358 keV | 1,709 | $8.79\times10^{-8}$ | .00355 |
| 80 keV | 6,039 | $9.46\times10^{-9}$ | .000383 |

Normalizing an elastic gross rate by the photon
floor alone is **not** a lower bound on its ratio
to the actual Hubble rate. The instantaneous
chemical gain divided by the unstimulated elastic
rate does not depend on that normalization. Even
for the larger coupling its upper bound drops from
$3.43\times10^{-7}$ at 2 MeV to $3.16\times10^{-10}$
at 80 keV on these three states. Fast encounters
should not be counted as fast net production.
This checks only the selected held shapes and
tree quartic; it does not prove a future freezeout
or bound the coupled nonequilibrium production history.

~~~sh
python -W error scalar_quartic_late_kinetic_chemical_clock_reference.py \
 EOS.late.dat \
 results/warm-prompt-reheating-post-2mev-elastic-2026-11-26.json.gz
~~~

[Held rates and forward/reverse ceilings](results/scalar-quartic-late-kinetic-chemical-clock-2026-11-26.json).
