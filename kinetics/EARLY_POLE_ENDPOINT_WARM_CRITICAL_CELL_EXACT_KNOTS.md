# Exact charge-event cuts in the largest warm source-error cell

The independent Volume160 controls select q=[3.96367894345,5.72884111204]
as the cell with the largest absolute net source residual in both
evolved rules. Keep each actual charge occupation and positive scalar
moment fit and the common native optical mass and clock. A charged
parent/daughter breakpoints coincidence occurs at E_s=E_p−E_d. Kinematic
collisions with each charged boundary occur at
E_s=(C E_charge±D p_charge)/m_charge², where C is the corresponding
two-body invariant energy product and D is the half Källén square root.
Use the rationalized lower arm to avoid subtractive cancellation.
The warm volume weight q²/(2π²) is integrated on each resulting q panel.

The 1176 exact scalar panels are evaluated with scalar/radial orders
(4,8) and (6,4). Their net results agree relatively to 5.7e−14 and the
inverse pieces to 4.3e−10. The earlier scalar order 2 is not used for
acceptance: despite close rates, it missed the independently known
count beyond the 3e−12 tolerance. The new orders pass that same gate.

Volume160 underestimates the selected shell by about 8.601e−7 of the
whole reference rate on both histories. The actual shell residuals
against the cut reference are +1.94123e−4 and +1.59722e−4 of the
respective whole references, for G7 and mixed255. The selected shell
still favors the mixed rule; this is not a global quadrature error
bound or a continuum dynamical inference.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
python -W error early_pole_endpoint_warm_critical_cell_exact_knots.py \
 --native "$p-x05-state-$d.json.gz" --states "$p"-{gauss7,mixed255}warm-x05-state-$d.json.gz \
 --control results/early-pole-warm-budgeted-gauss-probe-$d.json \
 --reference results/early-pole-warm-actual-source-higher-volume-check-$d.json \
 --cube results/early-pole-warm-actual-source-factor-cube-$d.json
```
