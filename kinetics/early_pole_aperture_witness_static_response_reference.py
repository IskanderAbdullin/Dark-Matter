"""Check radial integration by parts for steep witness logits."""
import json
import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_frozen_dynamic_vlasov_reference import C
from early_pole_live_aperture_witness_reprojection_reference import initial_projection
from early_pole_native_aperture_witness_reference import PROFILES


def check(prof):
    M,_,Y0,Y1,*_=initial_projection(prof,8,16)
    fit=[M.M.weights((1-z)*Y0+z*Y1,True)[1] for z in (0.,.50,.85,.97)]
    rows=[]
    for nodes in (12,16,24,32,48,64):
        x,w=roots_legendre(nodes)
        D=[]
        for F in fit:
            summed_direct=summed_endpoint=worst=0.
            for cells in M.M.r[1:]:
                for j in range(cells.start,cells.stop):
                    f=M.M.cells[j]
                    a,b=F[j]
                    h=(f.hi-f.lo)/2
                    p=(f.lo+f.hi)/2+h*x
                    W=h*w*C
                    E=np.hypot(p,f.m)
                    occ=expit(a+b*f.basis(p))
                    v=(-W*b*p*p/f.de)@(occ*(1-occ))
                    v+=C*(f.hi*np.hypot(f.hi,f.m)*expit(a+b)-f.lo*np.hypot(f.lo,f.m)*expit(a-b))
                    u=(W*(E+p*p/E))@occ
                    summed_direct+=u
                    summed_endpoint+=v
                    worst=max(worst,abs(v-u))
            D.append((summed_endpoint,summed_direct,worst))
        D=np.array(D)
        diff=(D[:,0]-D[0,0]-D[:,1]+D[0,1])*.005**2
        rows.append(dict(nodes=nodes,max_absolute_gauss_vs_static_change=float(np.max(abs(diff))),
                         max_unnormalized_cell_error=float(np.max(D[:,2]))))
    assert rows[-1]['max_absolute_gauss_vs_static_change']<1e-15
    return dict(grid=prof['grid'],hard_subdivision=8,rows=rows)


if __name__=='__main__':
    print(json.dumps([check(p) for p in json.loads(PROFILES.read_text())],indent=2,allow_nan=False))
