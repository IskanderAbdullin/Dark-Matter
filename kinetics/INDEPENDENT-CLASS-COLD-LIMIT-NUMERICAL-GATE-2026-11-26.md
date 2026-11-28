# A cold-limit numerical gate for native CLASS perturbations

Start with an analytic Fermi–Dirac diagnostic with the **same RMS
comoving momentum** as the 20 MeV retained fermion, not the archived
distribution itself. In unmodified [CLASS v3.3.4](https://github.com/lesgourg/class_public/releases/tag/v3.3.4),
replace the retained CDM density by this diagnostic species with
T_ncdm=0.896667207524836, m_ncdm=2e7 eV and omega_ncdm=.11816900374571636.
Match all other standard inputs. Use the same pure CDM reference
density plus the tiny synchronous gauge anchor omega_cdm=1e-8.

With Newtonian-gauge tolerances 1e-6 and 1e-7, the linear P(k)
ratios disagree by up to about 4.0 **in units of the reference P**.
A short Newtonian run with the actual file PSD encountered a
minimum-step integration error. These are numerical failures,
not evidence for structure enhancement or suppression.

With synchronous gauge and the anchor, changing tolerances
1e-6→1e-8 and Fourier sampling 24→45 gives max change about
5.4e-4 in the tested RMS-control ratios. The finer result differs
from matched CDM by at most 3.8e-4 for .01≤k≤2.9 h/Mpc at redshifts
0 and 5. Since the numerical control is larger than the deviation,
no physical suppression is resolved. Reject the earlier Newtonian
power, and validate the actual PSD separately in this gauge.

The comparisons use h=.6736, T_cmb=2.72548 K, omega_b=.02237,
A_s=2.1e-9, n_s=.965, tau_reio=.054 and N_ur=3.044, with linear
mPk output at z=0,5. Use kmax=3 h/Mpc for synchronous and 30
for the Newtonian precision comparison; the reported overlapping
band stops at 2.9. Name the native outputs with roots
`cdm_true_sync_`, `fermion_diagnostic_sync_`, `fd_true_sync_fine_`,
`fermion_fd_rms_` and `fermion_fd_rms_fine_`, respectively.
Run `class_cold_limit_numerical_gate_reference.py <output dir>`.
This is not a solved production perturbation or a transfer result
for the nonthermal PSD.
