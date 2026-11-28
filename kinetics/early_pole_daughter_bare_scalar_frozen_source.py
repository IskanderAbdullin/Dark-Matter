"""Hold initial hard and cold occupations fixed while probing scalar rest mass."""
import argparse
import gzip
import json
from types import SimpleNamespace

import numpy as np
from scipy.special import expit
from scipy.optimize import brentq

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_daughter_bare_scalar_stress_control import state,frozen_energy_cap
from early_pole_daughter_retained_kinetic_stress import build,r
from early_pole_decay_mode_clock_reference import G,KAPPA
from early_pole_live_pauli_gauge_moment_decay_reference import Tool


def source(C,e,F,N,m):
    M=C.M
    x=.05*np.exp(getattr(M,'eta',0.))
    es=np.hypot(C.q,m)
    pref=C.sw*G*KAPPA*x*((1.06*x)**2-m*m)/(C.q*es)
    out=np.zeros((2,2,len(N)))
    counts=np.zeros((2,len(N)),dtype=int)
    for a,(q,E,occ,A) in enumerate(zip(C.q,es,N/C.sw,pref,strict=True)):
        u=SimpleNamespace(n=M.decay.n,cells=M.decay.cells,q=q,nod=M.z,w=M.w)
        for i,j,z,v,w in Tool.events(u,dict(coeff=e[M.start:],es=float(E))):
            i+=M.start;j+=M.start
            fa=expit(F[i,0]+F[i,1]*z)
            fb=expit(F[j,0]+F[j,1]*v)
            k=M.cells[j]
            layer=0 if (k.lo>=1.4138888888888888-2e-14 and k.hi<=1.4666666666666666+2e-14) else 1
            out[layer,0,a]+=float(np.sum(w*(fa*(1-fb))*A))
            out[layer,1,a]+=float(np.sum(w*(occ*(fa-fb))*A))
            counts[layer,a]+=1
    return out,counts


def run(bundle,angular):
    C,depth=build(bundle,angular)
    with gzip.open(bundle,'rt') as f:
        arms=json.load(f)['branches']
    Y=np.asarray(arms['steps']['snapshots'][0]['hard_moments'])
    V=np.asarray(arms['cell-conservative']['snapshots'][0]['hard_moments'])
    M=C.M
    sel=np.array([i for i in range(M.r[2].start,M.r[2].stop)
                  if M.cells[i].lo>=1.4138888888888888-2e-14 and M.cells[i].hi<=1.4666666666666666+2e-14])
    if len(sel)!=2*depth:
        raise ValueError('Wrong original physical daughter interval.')
    L=Y.copy()
    L[sel]=V[sel]
    N=C.ns
    pole=C.scalar(N)[-1]
    cap=frozen_energy_cap(C,True)[0]
    q=np.r_[C.q,C.spq]
    n=np.r_[N,C.spect]
    states=[M.weights(z,True)[:2] for z in (Y,L)]
    vals=[]
    for ratio in (0.,.1,.2,.23,.245,.252,.254,.257,.27,.3,.4,.667,1.,2.,3.,3.13):
        s=state(q,n,C.lam,pole*ratio)
        m=s['effective_gap']
        a,ac=source(C,*states[0],N,m)
        b,bc=source(C,*states[1],N,m)
        control,_=source(C,states[0][0],states[1][1],N,m)
        if ratio==0.:
            for v,arr in zip((Y,L),(a,b),strict=True):
                ref=Cold.rhs(C,np.r_[v.ravel(),N,0.,0.],False)[C.j0:-2]
                np.testing.assert_allclose(arr.sum(axis=(0,1)),ref,rtol=3e-9,atol=2e-13)
        if m>cap*(1.+1e-6) and np.any(ac[0]):
            raise ValueError('Native local support above its checked analytic ceiling.')
        regular=b-a
        fixed=control-a
        if np.max(abs(fixed[1]))>1e-17:
            raise ValueError('A fixed-geometry local occupation change spilled into the outside source.')
        vals.append(dict(initial_bare_to_massless_gap_ratio=ratio,effective_scalar_gap=m,
                         native_total=float(a.sum()),local_total=float(b.sum()),
                         paired_total=float(regular.sum()),paired_by_window_and_direction=regular.sum(axis=2).tolist(),
                         fixed_geometry_paired=float(fixed.sum()),
                         native_aperture_active_modes=int(np.count_nonzero(ac[0])),
                         local_aperture_active_modes=int(np.count_nonzero(bc[0])),
                         native_total_active_modes=int(np.count_nonzero(ac.sum(axis=0)))))
    crossing=next(i for i,(a,b) in enumerate(zip(vals[:-1],vals[1:],strict=True))
                  if a['paired_total']<0.<b['paired_total'] and b['native_aperture_active_modes']>0)
    lo,hi=(v['initial_bare_to_massless_gap_ratio'] for v in vals[crossing:crossing+2])
    def diff(r):
        m=state(q,n,C.lam,pole*r)['effective_gap']
        A=source(C,*states[0],N,m)[0]
        B=source(C,*states[1],N,m)[0]
        return float((B-A).sum())
    zero=brentq(diff,lo,hi,xtol=1e-10)
    if abs(diff(zero))>2e-8:
        raise ValueError('The first sampled sign bracket had no reproducible root.')
    return dict(first_neg_to_pos_sampled_bracket=[lo,hi],
                first_bracket_sign_root=zero,
                scope='Frozen initial occupations, same inserted bare scalar rest mass in both arms. No scatter/expand continuation.',
                hard_depth=depth,cold_modes=len(C.q),rows=vals)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive')
    p.add_argument('--angular-mesh')
    a=p.parse_args()
    print(json.dumps(run(a.archive,r(a.angular_mesh) if a.angular_mesh else None),indent=2,allow_nan=False))
