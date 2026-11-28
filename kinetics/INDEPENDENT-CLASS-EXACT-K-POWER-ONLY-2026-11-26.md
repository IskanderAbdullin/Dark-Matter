# Avoid writing very large exact-mode histories when only power is needed

CLASS 3.3.4 uses k_output_values both to add exact
modes to the numerical k grid and to request their
perturbation-history files. The saved spectrum still
needs the first behavior. In a separate source tree
add a single opt-in condition to the **output**
stage, leaving the inputs, evolution, and numerical
k samples unchanged. CLASS_EXACTK_POWER_ONLY=1
suppresses only the writing of perturbation time
tables. With the variable unset, normal output is
unchanged. Never use this switch for calculations
that require raw initial times, metric histories,
or phase-space trajectories.

Compared with the unmodified executable on otherwise
identical inputs, the total and baryon+CDM power
files were byte-identical at z=5 and 20 in a two-mode
control (.3 and 2 h/Mpc), and again with seven exact
modes .3 through 10 h/Mpc in the 100 MeV absolute
reheating case with a transformed momentum grid.
The switch avoided writing 26,454,999 and 319,300,255
bytes of perturbation tables in the two controls.
It does **not** lower the array storage required
during evolution or certify equality on other
CLASS versions.

~~~sh
# Use a separate copy of the official CLASS 3.3.4 source.
cd owned-class-copy
python /path/to/class_exactk_power_only_output_patch_reference.py source/output.c output.next.c
mv output.next.c source/output.c
make -j2 class
CLASS_EXACTK_POWER_ONLY=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./class exact.ini
~~~

[Source and table checksums](results/class-exactk-power-only-output-validation-2026-11-26.json).
