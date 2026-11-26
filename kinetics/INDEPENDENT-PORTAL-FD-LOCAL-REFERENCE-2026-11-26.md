# Exact Fermi two-bath local number/heat reference, 2026-11-26
Independent numerical response and quadrature controls continuing PR #2. The parent reservoir and Bose dark-matter production are not modeled here. These results neither give an actual freeze-out abundance nor replace the physical program in upstream draft PR #1. Main f77e75a43e99f03c5ae4dc1d732731bbb9ba3fca; preceding PR #2 FV checkpoint 7528cffa0fb6d28430cb316c8a6d1a2a3eb2c40d.

## Scope and conservation
Common benchmark `m_chi2=20 MeV, m_V=10 GeV, g_D=0.1, epsilon=1e-5`, one symmetric Dirac pair species and internally thermal relativistic scalar. The SM bath follows the thermodynamically consistent EOS reconstructed from [Saikawa–Shirai 2020](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat) on its central tabulation. Only direct charged leptons are used in the displayed FD channels. No extra real-vector collision term atop the finite-width effective vacuum process. No hadron spin/crossing model or in-medium propagator. All reactions, numerical tests, and reruns were computed on the cloud, not GitHub Actions.

Reconstruct the ideal FD fugacity from `n(T,eta)` and integrate energy, pressure and entropy directly. At fixed density the partial heat capacity uses energy variance with weight `f(1-f)`, and the thermal scalar adds its independent heat capacity. In 35 controls the density inversion recovered eta to 4.71e-14 absolute, finite-difference heat capacity agreed to 4.02e-10 relative and direct entropy satisfied Gibbs to 4.45e-16 relative.

The timelike event kernel from the preceding angular test is integrated on positive subinterval quadrature in invariant energy and lab event energy. Threshold square-root coordinates and separate finite-width near-pole intervals maintain smoothness. Net forward-minus-reverse is evaluated with `expm1`, with event affinity `2 eta_d + (1/T_SM-1/T_D) E_lab`. Entropy is reconstructed from the positive event affinity times its net rate. Eleven tests (including near-equilibrium reverse cancellation and displaced baths) agree with the separate adaptive angular/energy quadrature within 6.65e-7 relative. The second flow moment uses the same lab energy as the entropy ledger; it is not the mediator pole mass.

## Classical control and model limitations

For comparison, the common-temperature no-fugacity MB collision matrix is derived from the same vacuum event weight. The classical ODE is a separately tested reference; the no-reaction evolution and the opposite energy ledger are controls. The nonlinear FD evolution and its thermal closure are separate follow-up tests, not inferred from local eigenvalues. No continuous hadronic spectrum was assumed. A decay treatment with a real vector requires a consistent separation from the full effective 2-to-2 kernel.

## Full-width spin-correlated FD reference


For a vacuum s-channel vector exchange, the spin-summed CM angular kernel is proportional to (3-v_i^2-v_f^2+v_i^2v_f^2\cos^2\theta_{if}\). Its azimuthal integration leaves products of the two baths' blocking Legendre moments of order 0 and 2. This treatment keeps spin correlation rather than assuming independent angles or separate unpolarized decays. With thermal FD occupations each gain/loss has the exact positive blocking product times an exponential of the relevant bath affinity. Every boosted event obeys detailed balance and nonnegative entropy production. A direct numerical integration over the two polar angles and relative azimuth, with forward and reverse FD occupations evaluated explicitly and extended floating-point precision near cancellation, checked 132 events to 1.98e-11 relative. The earlier double reference was too imprecise at the closest equilibrium states; this is recorded rather than weakening the check.

Then integrate the actual finite-width vacuum e/mu/tau kernel over invariant energy and lab event energy. Recover its independent MB Bessel moments within 4.5e-14; the response to simultaneous chemical and temperature displacements agrees with the equilibrium FD matrix within 2.9e-10. Common-temperature FD/MB number-response ratios are .82067 at T=.1 GeV, .82199 at .3, .99335 at .68825 and .76471 at 2. The energy susceptibility changes differently. Full angular spin corrections are below .04% in these integrated checks, not identically zero. Results are for local e/mu/tau vacuum kinetics with exact bath Fermi statistics, not in-medium propagators, real vector transport, hadron FD, or an evolved FD trajectory. This approach counts the full Breit-Wigner as an effective 2-to-2 process and must not also add real-vector events without a consistent separation.


## Reproduction

Use the authors' central EOS table linked above. Data are not redistributed. From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/mb_number_heat_state_reference.py
python -W error kinetics/fd_number_heat_state_reference.py
python -W error kinetics/finite_width_vector_pauli_angular_reference.py
python -W error kinetics/finite_width_vector_pauli_rate_reference.py
python -W error kinetics/positive_fd_portal_quadrature_reference.py
python -W error kinetics/nonlinear_timelike_mb_trajectory_reference.py EOS_TABLE MB_GRID.npy --build-grid --nodes 321
python -W error kinetics/nonlinear_timelike_mb_trajectory_reference.py EOS_TABLE MB_GRID.npy
```
