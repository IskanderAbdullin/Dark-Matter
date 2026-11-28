"""Independent raw-occupation cold aperture vs projected live initial slope."""
import json
import argparse
import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_snapshot_gauge_ledger_reference import cell_J,photon
from early_pole_decay_mode_clock_reference import G,KAPPA,intervals


def run(grid,n):
    src=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
    h=[]
    for l,m in [('parent',.05),('fermion',.003)]:
        ed=np.array(src[f'endpoint_{l}_momentum_edges'])
        f=np.array(src[f'endpoint_{l}_occupations'])
        h.append((ed,f,cell_J(ed,f,m,.45)))
    Jg=photon(1.,0.,.5)
    M,B=np.sqrt(np.array([.05,.003])**2+.005**2/np.pi**2*(Jg+np.array([a[2] for a in h])))
    out=[]
    for ap in (False,True):
        C=Cold(grid,0,0,n,ap)
        fs=C.ns/C.sw
        direct=[]
        for j,q in enumerate(C.q):
            S,lo,hi,_,_=intervals(M,C.m0,B,q)
            ep=np.hypot(h[0][0],M)
            eb=np.hypot(h[1][0],B)+S
            x=np.unique(np.r_[lo,hi,ep[(ep>lo)&(ep<hi)],eb[(eb>lo)&(eb<hi)]])
            mid=(x[1:]+x[:-1])/2
            a=h[0][1][np.searchsorted(ep,mid)-1]
            b=h[1][1][np.searchsorted(eb,mid)-1]
            pref=C.sw[j]*G*KAPPA*.05*((.05+.003)**2-C.m0*C.m0)/(S*q)
            direct.append(pref*np.diff(x)@(a*(1-b)+fs[j]*(a-b)))
        Z=np.r_[C.M.init.ravel(),C.ns,0.,0.]
        actual=C.rhs(Z,False)[C.j0:C.j0+C.n]
        dd=abs(np.array(direct)-actual).sum()
        out.append(dict(archive=grid,scalar_bins=n,aperture_refined=ap,
                        raw_by_scalar_bin=direct,projected_by_scalar_bin=actual.tolist(),
                        raw_total=sum(direct),projected_total=actual.sum(),
                        normalized_l1=dd/sum(abs(np.array(direct)))))
    return dict(scope='Initial rates in the same finite vacuum-width proxy; raw square-root shift '
                      'versus live variational square-root-order approximation; not physical rates.',
                rows=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',choices=['fine','third'],default='fine')
    p.add_argument('--n',type=int,default=8)
    o=p.parse_args()
    print(json.dumps(run(o.grid,o.n),indent=2,allow_nan=False))
