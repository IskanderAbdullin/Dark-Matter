# Tightening time tolerances preserves the compensation-only factorial result

Repeat the four held simultaneous preparations (native, local,
compensation only, both) with relative tolerance 3e-10 and base
absolute tolerance 2e-15 rather than 3e-9 and 2e-14. Match each
retune exactly, along with initial data and all reporting times.
Both method settings pass the original numerical closure gates.

The largest paired count-history change across any stored
case/time is 1.04e-17. The final paired local and combined
mean shifts change by 1.04e-11, the outside-only shift is
unchanged, and the factorial interaction is identically
+1.6278645099e-7 at the printed precision. Thus numerical
tolerance tightening does not explain the small nonadditivity
of this particular compensation. This says nothing about a
different physical prescription, hard momentum grid or
elastic kernel.

Use `early_pole_daughter_complement_factorial_compare.py`
with the four nominal receipts and `--tighter` in the same
native, local, compensation-only and combined order.
