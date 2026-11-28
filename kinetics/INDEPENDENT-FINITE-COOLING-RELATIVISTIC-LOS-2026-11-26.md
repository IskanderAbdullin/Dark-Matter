# Memory-bounded relativistic tracer check on recorded potentials

The [fluid-switch test](INDEPENDENT-FINITE-COOLING-TRACER-SUPERPOSITION-2026-11-26.md)
shows a numerical grouping dependence for the 86 eV prompt-plus-recoil
tracer. This quiet check does not evolve the metric or impose a fluid
closure. On the same recorded Newtonian-gauge metric, it integrates
the collisionless distribution along its momentum-dependent geodesics.
It addresses two exact modes, k=4 and 10 h/Mpc at z=20. The dominant
cold species and all other species are prescribed by the recorded
history; no inference for total matter follows.

Let \(G=\delta f\), \(D=q\,df_0/dq\),
\(v=q/\sqrt{q^2+(am/T_{\gamma,0})^2}\), and
\(L(\tau')=\int_{\tau'}^{\tau_f}v\,d\tau\). Then
\[
\dot G + i k\mu vG=-D\,[\dot\phi-i k\mu\psi/v].
\]
The density monopole is the free angular projection of the
initial data plus a \(-D\dot\phi\,j_0(kL)\) and
\(Dk\psi\,j_1(kL)/v\) integral. Integrating the first term by parts
gives boundary terms plus \(Dkv\phi\,j_1(kL)\). Both formulations
are computed independently. Initial angular data are approximated
by the recorded early adiabatic radiation monopole, dipole and
quadrupole, at \(k\tau<0.1\); the sensitivity to dropping the latter
two is retained. The calculation directly integrates \(\delta\rho\)
and \(\delta P\) as separate momentum moments of \(G_0\).

At k=10, the alternate LOS calculation gives
\(\delta_s/\delta_{\rm CDM}\simeq0.759673\) on **either** the
fluid-enabled or fluid-disabled recorded metric. The native raw
time-trajectory probes at the same redshift are 0.522630 with
fluid and 0.759671 without. These raw time-trajectory values
need not equal interpolation of the native transfer output
at its printed redshift; its small difference is separately retained.
At k=4 the independent LOS ratio is 0.881305, versus the native
raw no-fluid value 0.881304. These differences at the sampled
modes are below \(2.2\times10^{-6}\) relative to the native
no-fluid raw trajectory. This is a consistency result, not
all-k convergence or a self-consistent solution with a
dominant scalar component.

The common-metric \(f_0\) split adds with the true density weights.
At k=10 the isolated prompt and recoil LOS ratios are 0.853392
and 0.348962. The computed joined \(\delta P/\delta\rho\) is
\(9.884\times10^{-10}\) at this epoch, compared with the
fluid run's assigned \(3.966\times10^{-9}\); the no-fluid
kinetic moment agrees with the independent calculation.
This measures the pressure-closure issue within the test,
not an error bar on observations or on the production model.

## Portable inputs and reproducibility

The recorded drive was compacted to 1,600 geometric time samples
for each mode, storing only scale factor and two potentials in
float32, with the initial and final native reference values
stored separately in float64. The native source SHA-256 and
compact-file hashes are retained. The fixture can be used without
running large native simulations; it does not regenerate their
metric or constitute raw native output. The LOS implementation
computes the integration in small radial chunks. The largest
portable-vs-raw density change is 1.29×10⁻⁸ absolute;
the largest pressure-ratio change is 4.66×10⁻⁸
relative on the tested tables. Both comparisons are independently
audited, not inferred from agreement of the binary file sizes.

Regenerate unremapped sigma=0.14 mixed and prompt tables with the
previous exporters, and finite_pair_only_ratio5_s14.dat with the
split table generator. Then use:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_relativistic_los_reference.py \
  EOS.dat MIXED.dat PROMPT.dat RECOIL.dat results/los-metric-off/*.npz
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_relativistic_los_reference.py \
  EOS.dat MIXED.dat PROMPT.dat RECOIL.dat results/los-metric-on/*.npz
~~~

To regenerate the compact files from native perturbation tables,
pass the two *_00_perturbations_k*_s.dat inputs rather than the
npz files, adding --pack DIRECTORY --pack-nodes 1600.
