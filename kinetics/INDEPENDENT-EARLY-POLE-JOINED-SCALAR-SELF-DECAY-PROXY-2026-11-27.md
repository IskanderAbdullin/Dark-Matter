# Scalar redistribution and reciprocal decay with both hard and scalar states live

The earlier [isolated shape diagnostic](INDEPENDENT-EARLY-POLE-LIVE-SCALAR-SHAPE-FEEDBACK-2026-11-27.md) held the charged populations. Here both charged species evolve and feed back on the scalar through the reciprocal decay/inverse kernel. Compare matched histories with and without the closed, scalar-only shape stencil. Its vertex and scattering dispersion are held at their initial values. Scalar entropy and the residual work of the fixed-mass stencil against the changing Hartree energy are integrated independently. The scalar shape stages surround the coupled decay with symmetric half steps.

| Input, 16 scalar bins over exposure 2e−7 | Baseline scalar gain | With scalar self-redistribution |
| --- | ---: | ---: |
| Fine, native steps | +1.58500e−8 | +1.57766e−8 |
| Fine, 97% prepared parent | −1.36436e−9 | −9.88457e−10 |
| Third, native steps | +1.45190e−8 | +1.42637e−8 |
| Third, 97% prepared parent | +1.30336e−8 | +1.04883e−8 |

These are **cumulative** gains, not the final instantaneous source. The included scalar proxy reduces the magnitude of the fine witness's net loss by about 27.6%, without changing its sign. It reduces the third witness's net gain by 19.5%; the native-step fine and third changes are 0.46% and 1.76%. Doubling the fine time steps changes the witness's proxy gain by about 0.10%. Increasing both time and angular resolution on third gives +1.04592e−8 on the witness branch, likewise preserving sign. Neither means all hard phases have been tested.

The live flux integrates its entropy and corrected energy work in separate quadratures on the positive dense shape trajectory. A [finite Hartree secant check](INDEPENDENT-EARLY-POLE-SCALAR-HARTREE-DISCRETE-GRADIENT-2026-11-27.md) on every shape half step of a short fine witness control independently agrees with the actual endpoint scalar energy changes to 4.1e−21. This checks the held graph; it does not alter its weights. The measured number and entropy ledgers close. This comparison includes no screened charged elastic scattering, separate hard photon or flavor processes, moving scalar support, or FRW expansion. The scalar vertex is still a sensitivity proxy, not a controlled N=1 result or a physical correction to other models.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
for g in fine third; do
  python early_pole_joined_scalar_self_proxy_reference.py "$g"
done
python early_pole_joined_scalar_self_proxy_reference.py fine --steps 16
python early_pole_joined_scalar_self_proxy_reference.py third --steps 16 --order 10
python early_pole_joined_scalar_self_proxy_reference.py fine --steps 4 --duration 2e-8 --mixes .97
```

Raw values: `results/early-pole-joined-scalar-self-decay-proxy-2026-11-27.json.gz`.
