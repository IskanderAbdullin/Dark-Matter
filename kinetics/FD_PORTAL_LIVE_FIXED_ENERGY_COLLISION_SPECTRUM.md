# Fixed-energy collision derivatives along the direct FD histories

Measure the Jacobian of the instantaneous collision field in the physical variables dark number density and dark energy density. For each perturbation, solve for the two temperatures and dark chemical potential while holding total energy fixed. Compare with the equal-temperature, zero-chemical-potential state at the **same total energy**, so the factor H is identical. The matrices act on relative dark number and energy displacements; their positive rate labels are the eigenvalues of minus the collision Jacobian divided by H.

At the original slow reference temperatures the two live slow rates are about 0.9901 and 0.9888 in units of H, about 0.9–1.0% below their same-energy equilibrium comparators. At the original fast reference temperatures, the live fast rates are about 0.7951 (timelike) and 0.8021 (with elastic scattering), compared with 0.9793 and 0.9806 in the corresponding equilibrium states. The negative shifts are about 18.8% and 18.2%. At actual temperature equality the sampled live rates are still 6.2–8.3% below the same-energy equilibrium rates. At the final x=.96 checkpoint neither operator has a sampled rate above H.

The older reference temperatures were defined near equilibrium at the specified visible temperature, whereas the comparator here has the same **total energy** as the live state. Its visible temperature is generally different, and its equilibrium rate at a reference label need not equal exactly H. This accounts for the distinction between about 0.80H and the quoted 18–19% relative shifts.

Two step sizes and two independent quadratures resolve eight live spectra. Step-size differences are below 3.9e-12 relative and quadrature differences below 1.4e-8; at equilibrium the numerical Jacobian agrees with the independently constructed thermodynamic mobility spectrum below 4.7e-12. The two sampled instantaneous eigenvalues are real and positive in both physical channels at all eight states. None of these tests proves stability of the full expanding variational flow or defines a time of thermal or chemical decoupling.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results; e=eos-range.dat
python -W error fd_portal_live_fixed_energy_collision_spectrum.py "$e" \
 "$p/fd-portal-live-states-at-local-clock-temperatures-$d.json" \
 "$p/fd-portal-direct-temperature-equality-with-chemical-lag-$d.json" \
 "$p/fd-portal-direct-expansion-055-096-$d.json"
```

Use the normalized public EOS described in `FD_PORTAL_TABULATED_BATH_ENTROPY_AND_MODES.md`. This remains a vacuum charged-lepton two-bath numerical check, not the full particle production or cosmological transfer problem.
