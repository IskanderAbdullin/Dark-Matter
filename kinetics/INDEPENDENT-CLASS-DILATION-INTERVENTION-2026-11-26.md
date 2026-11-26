# CLASS derivative intervention: identify a numerical power error

The [unchanged native phase test](INDEPENDENT-CLASS-MOMENTUM-NODE-PHASE-2026-11-26.md)
implicates momentum sampling. An additional **unmodified native**
256-node run, with the same fine tabulated input, matched cold
reference, qmax 20 and the same fluid transition, makes the
maximum excess *larger*, 1.0440% at k = 0.2383 h/Mpc. The absolute
number quadrature error falls to 0.0742%, but the number-dilation
identity residual rises to +1.2486%. Simply adding more nodes at
that endpoint is not a convergence test.

As a causal numerical diagnostic, copy the upstream CLASS
v3.3.4 build and change only the perturbation input
$d\ln f_0/d\ln q$ by one constant factor per species.
Set the raw discrete number dilation moment to the identity
$$
  \sum_i w_i q_i^2 \left(\frac{d\ln f_0}{d\ln q}\right)_i
  = -3\sum_i w_i q_i^2.
$$
Here the native weights $w_i$ already include $f_0$. Leave the
tabulated distribution, cosmology, nodes, background weights and
all other precision controls unchanged. This is an instrumented
experiment, not an upstream fix or a usable solver prescription.

| Fine input, 64 nodes | Native | Perturbation derivative rescaled |
| --- | ---: | ---: |
| Maximum departure from matched CDM at z=0 | 0.00609992 | 1.904e-6 |
| Maximum departure at z=5 | 0.00609940 | 1.497e-6 |

The raw native residual before multiplying by 1.0073723
is −0.007318335. An independent cubic-spline derivative gives
−0.007318266. The two native background arrays and the power
wavenumber grids are exactly identical in the baseline and
instrumented 64-node runs. Since the power departure mostly
disappears when only the discrete derivative normalization is
changed, the departure in this matched case is a numerical
artifact. The tiny diagnostic residual is **not** a prediction of
the physical transfer. This change need not preserve other
moments of the collisionless hierarchy, or the correct initial
conditions before the archived production surface.

## Reproduce

Prepare the configurations and unmodified v3.3.4 cases as in the
pressure-balanced test. Keep the original tree and executable
unchanged. Build a separate copy for the diagnostic:

~~~sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error class_native_dilation_diagnostic_patch_reference.py \
   original-class/source/background.c diagnostic-class/source/background.c
make -C diagnostic-class class
sed 's/psd_raw64_matched_sync_/psd_dilation64_sync_/' \
    native/psd_raw64_matched_sync.ini >native/psd_dilation64_sync.ini
diagnostic-class/class native/psd_dilation64_sync.ini >native/psd_dilation64_sync.log 2>&1
python -W error class_native_dilation_audit_reference.py native fine/stable.dat --nodes 64
~~~
