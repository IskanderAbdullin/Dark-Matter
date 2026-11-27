# Replacing the bare-energy target by modeled energy barely changes the initial compensation

Keep the same two aperture-local daughter cells and the same cold
seed. Retune only the 18 outside daughter cells by a common
two-parameter affine occupation correction. Match retained daughter
number plus the full hard-sector energy used by the held model,
including its field term. Contrast this preparation with the
previous matched number and bare daughter energy fit. The same
0.2 coefficient gate applies. The analytic gradient of modeled
energy was independently checked by central differences in the
two fit parameters.

The bare-energy fit leaves a modeled hard energy difference of
-1.0297e-11 in the code's energy units. The modeled-energy
fit makes that residual zero at recorded precision, with
daughter number residual -3.47e-18; its residual bare daughter
energy is +1.0297e-11. Their largest individual moment-component
difference is 1.34e-11. At this precision the instantaneous
local-minus-native cold-rate shifts are identical
(-0.0155110728714033). The initialized local suppression also
remains close to the uncompensated -0.0155105610126819.

This compares two static targets within the specified dressed
model. It neither replaces a check on the finite evolution nor
proves another physical energy model would make the same match.
Use `early_pole_daughter_dressed_energy_retune.py` with the
matched fine raw archive and original cold seed.
