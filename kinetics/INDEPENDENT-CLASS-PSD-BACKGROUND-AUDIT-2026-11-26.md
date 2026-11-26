# Native CLASS background audit of the exported low-reheat spectra

The 0.3315 MeV retained snapshot is exported as a CLASS file PSD
for each species, without thermal replacement. The positive
log-momentum remapping preserves the represented number integral,
and the first three momentum moments are checked at two resolutions.
Set T_ncdm=1 and deg_ncdm=1. The dimensionless coordinate is
q=p_comoving/Tγ(today); the second column already includes the
physical spin multiplicity divided by (2π)³. No exterior component
is silently put into either file.

Compiled the official [CLASS v3.3.4 release](https://github.com/lesgourg/class_public/releases/tag/v3.3.4)
in the cloud. The GitHub tag tarball SHA-256 was
35f2599f0eaa680d2f8198ddb7f70099df7c00a6c04f239c142cd8f3c420545c.
See the [release input examples](https://github.com/lesgourg/class_public/blob/v3.3.4/explanatory.ini)
and [Lesgourgues and Tram](https://arxiv.org/abs/1104.2935).
The automatic PSD quadrature did not converge at its requested
tolerance, so the independent audit uses manual finite q grids
at four resolutions. The raw runs specify the two particle
masses, not omega_ncdm, so CLASS must infer density from the
actual file normalization.

The two represented late rest densities are roughly Ωχh²=.118169
and Ωφh²=2.990e-10. This is *retained only*, not the earlier
inclusive .11970. The native energy density and pressure are
compared against independent unsmoothed shell integrals from
today back to the actual snapshot. On the 4160-node grid, raw normalizations differ by 4.2e-7 and 5.6e-7; the maximum pressure-to-energy errors are 7.5e-5 and 2.1e-4. Detailed errors and convergence
are recorded in the output JSON, not assumed from input format.
The phase-space background audit does not establish an
adiabatic source perturbation or matter transfer prediction.
A separate estimate gives kη≈.018 at the snapshot even for
100 h/Mpc; this is a long-wave check and does not implement
perturbations in the earlier production epoch.

## Reproduce

Build unmodified CLASS v3.3.4. From this directory use the
EOS extract and the late wide snapshot:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/low-reheat-quantum-bridge-wide-shape-2026-11-26.json
python -W error export_archived_class_psd_reference.py EOS.late.dat "$a" class_psd_base
python -W error export_archived_class_psd_reference.py EOS.late.dat "$a" class_psd_fine --sigma .006 --nodes 6200
```

For each raw background, create a CLASS ini below, selecting the
base table for 260 nodes and fine table for 520, 1040 and 4160.
Use absolute paths and a different output root per run.

```ini
h = .6736
T_cmb = 2.72548
omega_b = .02237
N_ur = 3.044
omega_cdm = 0
N_ncdm = 2
m_ncdm = 20000000,86
T_ncdm = 1,1
deg_ncdm = 1,1
use_ncdm_psd_files = 1,1
ncdm_psd_filenames = /path/stable.dat,/path/scalar.dat
ncdm_quadrature_strategy = 3,3
ncdm_maximum_q = 20,96
ncdm_N_momentum_bins = 260,260
write background = yes
root = /path/out_
```

Run CLASS, leaving perturbation output unset, and pass the resulting
files `*_background.dat` to
`archived_class_background_audit_reference.py EOS.late.dat <snapshot> <backgrounds ...>`.
