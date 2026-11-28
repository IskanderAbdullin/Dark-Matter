# Signed bare soft kernel on the cold witness at two frequencies

At q=4.1e−6 use the frozen fine and third spectra, actual archived scalar occupation n_s, and the [bounded parent witness](INDEPENDENT-EARLY-POLE-NATIVE-APERTURE-WITNESS-2026-11-27.md). The witness changes the resolved shape within selected parent momentum cells, keeping three radial moments fixed. The benchmark compares parent mixture fractions 0, .85 and 1 on the **same refined boundaries**, including the unmodified mixture. This removes a photon-angle mesh change that could otherwise be mistaken for nonlinear occupation dependence. This is a vacuum bare real kernel with only the leading soft subtraction, at the same fixed cold scalar mode as the [parent-integrated aperture control](INDEPENDENT-FINITE-MASS-QED-PARENT-INTEGRATED-SIGNED-SOFT-2026-11-27.md).

Let A_+ and B_+ be the exact integrated Pauli factors f_P(1−f_b) and (1−f_P)f_b, with the signed matrix and phase weights for positive outgoing photon energy. Use A_− and B_− for incoming energy, and X=A−B. Remove the same zero-energy factor C = X₀ + A₀/n_s and divide out n_s. The appropriately signed scalar kernel becomes

```text
spontaneous = A_+ − B_− + A_+/n_s − C
stimulus    = X_+ + X_− + (A_+ + A_−)/n_s − 2C
combined    = spontaneous + stimulus / expm1(ω/Tγ),  Tγ=.9
```

The two exact thermal directional forward/inverse factors balance separately before subtraction, in an independent common-temperature check with a shared charged chemical potential and a neutral scalar. This identity does not impose such a relation on the frozen archive. The measured n_s values are 1.06697e16 fine and 1.87234e15 third; the expression still keeps the +1 in scalar emission. The table lists the normalized **combined integrand values prior to the frequency measure dω/ω and the gauge factor**.

| Archive | ω | Baseline parent | .85 witness |
| --- | ---: | ---: | ---: |
| fine | 3e−6 | −.0000297673 | +.0035853224 |
| fine | 1e−5 | −.0001659417 | +.0127874531 |
| third | 3e−6 | +.0000101162 | +.0008432356 |
| third | 1e−5 | +.0000483734 | +.0035497467 |

The final angular/parent mesh difference is at most 4.54e−10 absolute and the checked signs are invariant under that refinement. Direct .85 evaluations agree with the affine reconstruction to 7.98e−16 absolute on the same boundaries. An independent parent-first integral checks the **two signed fine stimulus components** at ±1e−5 and .85 mixture to 3.06e−18 absolute. Common-temperature directional factors close to 8.89e−16 relative.

To regenerate the [receipt](results/finite-mass-qed-signed-cold-witness-2026-11-28.json) from `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for g in fine third; do
  python -W error finite_mass_qed_signed_cold_witness_reference.py "$g" 3e-6 1e-5 \
    --meshes 0 1 2 > "$g.jsonl"
done
python -W error finite_mass_qed_signed_cold_witness_reverse_audit.py > reverse.jsonl
python -W error finite_mass_qed_signed_cold_witness_receipt.py fine.jsonl third.jsonl reverse.jsonl
```

These samples show how the already-subtracted bare kernel responds to the within-cell redistribution. They do not locate a zero of a frequency-integrated correction, give a rate for the screened plasma below its gap, or predict the sign of any matched NLO scalar term.
