"""Critical Bose occupation identity and its singular limit.
Run Python 3 with mpmath: python critical_bose_occupancy.py
No simulations or GitHub Actions.
"""
import json
from mpmath import mp
mp.dps=40
critical=mp.zeta(mp.mpf("1.5"))
rows=[]
for d in ["0.1","0.01","0.0001","0.000001","0.00000001"]:
    delta=mp.mpf(d);z=mp.exp(-delta)
    nu=mp.polylog(mp.mpf("1.5"),z)
    enhancement=mp.polylog(mp.mpf("0.5"),z)/nu
    rows.append(dict(delta=float(delta),nu=float(nu),
       number_weighted_factor=float(enhancement),
       deficit_times_factor=float((critical-nu)*enhancement)))
expect=float(2*mp.pi/critical)
assert abs(rows[-1]["deficit_times_factor"]/expect-1)<2e-5
print(json.dumps(dict(critical_nu=float(critical),
 limit_deficit_factor=expect,rows=rows,
 caution="A number-weighted occupancy average is not the full quantum collision rate."),indent=2))
