"""Compare native single- and double-temperature controls at fixed density."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator,CubicSpline
from warm_direct_vacuum_native_focus_reference import K


def run(standard,exact):
    a=[]
    ks=np.asarray(K)
    for d,names in ((standard,("std-amp3","std-thermal3","std-p12")),(exact,("direct3","thermal3","p12"))):
        reports=[]
        for zi,z in enumerate((5,20),1):
            front=[];alt=[];sha=[]
            for n in names:
                p=d/f"{n}_00_z{zi}_pk.dat"
                x=np.loadtxt(p)
                if d==exact:
                    inds=np.searchsorted(x[:,0],ks)
                    assert np.array_equal(x[inds,0],ks)
                    q=np.log(x[inds,1])
                    h=q
                else:
                    q=PchipInterpolator(np.log(x[:,0]),np.log(x[:,1]))(np.log(ks))
                    h=CubicSpline(np.log(x[:,0]),np.log(x[:,1]))(np.log(ks))
                front.append(q);alt.append(h)
                sha.append(hashlib.sha256(p.read_bytes()).hexdigest())
            r=np.vstack(front)
            s=np.vstack(alt)
            dN=r[1:]-r[0]
            diff=np.max(abs((s[1:]-s[0])-dN))
            assert abs(dN[1,-1])<.17*abs(dN[0,-1])
            reports.append(dict(z=z,raw_power_sha256=sha,
               single_q2_and_double_q1_q2_delta_logP_vs_direct=dN.tolist(),
               residual_double_over_single_at_k1=float(abs(dN[1,-1]/dN[0,-1])),
               max_alternative_interpolation_shift=float(diff)))
        a.append(dict(method="exact early RK" if d==exact else "standard NDF",cases=reports))
    return dict(scope="Native first- and second-moment PSD controls, conditional amplified fraction 1e-3. "
                      "Common reference dark density and 86 eV scalar mass. No thermal scalar microphysics or extrapolation.",
                k_h_per_Mpc=list(K),reports=a)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("standard",type=Path);p.add_argument("exact",type=Path)
    a=p.parse_args();print(json.dumps(run(a.standard,a.exact),indent=2,allow_nan=False))
