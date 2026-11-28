"""Finite-band integral with the local leading soft term subtracted first."""
import gzip,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from finite_mass_qed_signed_support_bridge_reference import BridgeSource


def seed(S):
    edge=np.unique(np.clip(np.r_[S.lo,S.hi,S.pe,S.be+S.Es],S.lo,S.hi))
    mid=(edge[1:]+edge[:-1])/2
    p=S.pf[np.searchsorted(S.pe,mid)]
    f=S.bf[np.searchsorted(S.be,mid-S.Es)]
    return S.B*(np.diff(edge)@np.stack([p*(1-f),(1-p)*f],axis=1))


def one(S,n,outside,base,fs):
    X,W=roots_legendre(outside)
    result=[]
    edges=np.array([1e-5,2e-5,4e-5,8e-5,1.6e-4,2e-4],dtype=np.longdouble)
    for a,b in zip(edges[:-1],edges[1:],strict=True):
        x,y=np.log(a),np.log(b)
        k=np.exp((x+y)/2+(y-x)*X/2)
        vals=[]
        for q in k:
            c=S.signed(q,n,9,28)[0,:2]
            d=S.signed(-q,n,9,28)[0,:2]
            vac=np.array([c[0]-base[0],d[1]-base[1]])
            th=(c+d-2*base)/np.expm1(q)
            vals.append(np.r_[vac,th])
        raw=(y-x)/2*(W@np.array(vals,dtype=np.longdouble))
        G,L=raw[:2]+raw[2:]
        Born=(base[0]*(1+1/fs)-base[1])/S.B
        r=(G*(1+1/fs)-L)/Born
        result.append(dict(edges=[float(a),float(b)],vacuum_G_L_then_thermal_G_L=raw.astype(float).tolist(),
                           unscaled_net_ratio_to_Born=float(r),net_ratio_to_Born_at_g_point005=float(r*.005**2/(4*np.pi**2))))
    return result


def run(grids):
    res=[]
    for g in grids:
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        S=BridgeSource(d)
        ref=seed(S)
        err=float(max(abs(S.signed(0.,54,12,40)[0,:2]-ref)))
        if err>5e-12: raise ArithmeticError(('Analytic eikonal baseline did not reproduce',err))
        fs=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],float(S.q))-1]
        a,b=(one(S,n,o,ref,fs) for n,o in ((28,4),(46,8)))
        diff=max(abs(x['unscaled_net_ratio_to_Born']-y['unscaled_net_ratio_to_Born']) for x,y in zip(a,b,strict=True))
        res.append(dict(grid=g,baseline_max_abs_defect=err,window_integrals=b,
                        max_mesh_drift_in_unscaled_net_over_Born=diff,
                        whole_bridge_scaled_net_ratio_to_Born=sum(x['net_ratio_to_Born_at_g_point005'] for x in b)))
    return dict(scope='Paired signed recoil with leading soft reference only over lab photon band [1e-5,2e-4].',data=res)


if __name__=='__main__':
    import sys
    print(json.dumps(run(tuple(sys.argv[1:]) or ('fine','third')),indent=2,allow_nan=False))
