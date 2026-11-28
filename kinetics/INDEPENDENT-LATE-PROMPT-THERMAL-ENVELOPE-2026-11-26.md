# Full-support bound on late prompt heat for thermal trial tails

The archived momentum mesh excludes retained on-shell production after
b=3.63867. This check does not use that cutoff: it estimates the omitted
on-shell channel over the entire continuum for FD/BE trial distributions
at the same instantaneous separate counts and moment-equivalent
temperatures. It cannot bound arbitrary nonthermal high-momentum tails.

For a parent lab energy E the physical positive inverse angular weight
is A=<f_phi f_chi>. The gross inverse event rate is at most the decay
phase-space measure times A, because the parent Pauli factor is at most
one. In the stationary-parent approximation, total prompt angular
redistribution following each such event transfers at most the full
allowed daughter energy span. This bounds transported scalar kinetic
energy, not formation or annihilation of dark species.

With kinetic fugacities z_phi<1 and z_chi, bound the product by
z_phi z_chi/(1-z_phi) times exp[-(E-m_phi-m_chi)/max(T_phi,T_chi)].
The full-shell event integral then has a K1 Bessel moment; replacing the
daughter interval by E gives a K2 energy bound. An independent positive
angular quadrature gives a tighter numerical envelope within the trial
family. It does not condition on the retained momentum boundary.

For the 48-step radiation-clock trace, at b=3.63867 the analytic heat
envelope is about 0.42 of scalar kinetic energy per Hubble time, with
the explicit allowed-angle envelope about 0.020. At b=4.50 the analytic
envelope drops below 3e-5, and at b=5.10 below 3e-8. Thus a thermal
proxy has a rapidly disappearing missing channel, but the tightness of
the immediate transition and nonthermal tail shape still matter.
This does not rule out a tail-driven effect without an exterior
distribution check. The accompanying 48-step actual retained trace
gives final scalar moment N/N_sat=0.002389 and first-cell fraction
4.154e-8, against 0.002331 and 4.130e-8 at 36 steps.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 python -W error late_prompt_thermal_envelope_reference.py results/12-quantum-48-capacity-shape-2026-11-26.json results/late-archived-radiation-prompt48-2026-11-26.json
```
