# The real-parent pole factors, but the archived parent has not relaxed

This is a frozen single-channel vacuum-matrix-element check for bookkeeping at the pole, not a calculation of the full in-medium parent spectral function. It uses the same two archived distributions and incoming/outgoing daughter momentum support [.45,5] as the [guarded annuli](INDEPENDENT-EARLY-POLE-ARCHIVE-SCALAR-NEAR-POLE-2026-11-27.md). The scalar and daughter occupations are fixed; the on-shell parent is not dynamically evolved here.

For (x=s/M^2-1), the spin-averaged (x^2)-rescaled vacuum squared amplitude at the real pole is angle independent, with value ([((M+m)^2-m_s^2)/M^2]^2). We independently perform the incoming shell integral using its analytic endpoints and the outgoing angular averages by exact overlap of lab-energy intervals. This verifies the fixed-invariant pole density (including the signed source) within 8e−10 relative in both archives and the tested q modes.

In a model with only this daughter–scalar channel, define (A=\langle f_d n_s\rangle) and (B=\langle(1-f_d)(1+n_s)\rangle) at each parent momentum. The same retained outgoing angles enter both. A hypothetical stationary parent has (f_p=A/(A+B)); its damping multiplies the bare width by (A+B=\langle1-f_d+n_s\rangle). Divide each pole-residue packet by this local width in the delta-function narrow-width weighting. The result equals the direct collision term with that stationary parent, including both gross rates and the signed net, within 4e−8 relative at order 32 for q=4.1e−6. Common-temperature FD/BE inputs reproduce the parent FD function locally. The equality is a consistency check on the repeated interaction, not an additional rate to add to the direct calculation with an explicit parent.

| Archive, q=4.1e−6 | Direct net with archived parent / Href | Direct net with stationary replacement / Href | Upper bound on lab relaxation / Href |
| --- | ---: | ---: | ---: |
| Fine | +4.42151 | −2.326e−8 | 158.45 |
| Third | +3.79437 | −3.693e−8 | 162.04 |

The relaxation upper bound uses the integral of the complete positive archived scalar occupation over lab kinetic energy in the angular numerator, plus the vacuum term, divided by the smallest outgoing angular window. It is for this channel, these supports, and the frozen input states. Even the upper bound times the largest existing fixed exposure 7.2e−5 is under 0.012 at the displayed q. The very different stationary source therefore cannot replace the archived source over such a short interval. Other interactions, time-dependent occupations, the dynamic dispersive response, and genuinely in-medium resonance/decay matching remain outside this test. See [Ala-Mattinen et al.](https://doi.org/10.1103/PhysRevD.108.096034) for the general on-shell subtraction problem.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
for g in fine third; do
  for n in 8 14 24 32; do
    python early_pole_archive_scalar_pole_factorization_reference.py "$g" \
      --order "$n" --q 0 2e-6 4.1e-6 8e-6
  done
done
```
Raw: `results/early-pole-archive-scalar-pole-factorization-2026-11-27.json.gz`.
