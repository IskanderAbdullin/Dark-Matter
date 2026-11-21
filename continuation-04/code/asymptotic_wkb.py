"""Near-critical universal λ_Jeans / r_half asymptotics.
Depends on Python mpmath and the existing wkb_consistency.py in this directory.
Run: python asymptotic_wkb.py
"""
import json
from mpmath import mp
from wkb_consistency import half_response_radius
mp.dps=40
nuc=mp.zeta(mp.mpf("1.5"))
C=2*mp.pi*mp.sqrt(nuc)/3
rows=[]
for nu in (2.2,2.5,2.6,2.61,2.612):
    full=half_response_radius(nu)
    asym=C/mp.sqrt(nuc-mp.mpf(str(nu)))
    rows.append(dict(nu=nu,exact_root_ratio=full["lambda_over_half_response_radius"],
       universal_asymptotic_ratio=float(asym),
       asymptotic_to_exact=float(asym/full["lambda_over_half_response_radius"])))
assert abs(rows[-1]["asymptotic_to_exact"]-1)<.0002
assert abs(rows[-2]["asymptotic_to_exact"]-1)<.001
result=dict(nuc=float(nuc),coefficient_2pi_sqrt_nuc_over3=float(C),
   asymptotic="lambda_J/r_half ~ C / sqrt(nuc-nu), as nu approaches nuc from below",
   tests_passed=2,rows=rows,
   warning="The comparison is local homogeneous WKB consistency, not a global stability diagnostic.")
print(json.dumps(result,indent=2))
