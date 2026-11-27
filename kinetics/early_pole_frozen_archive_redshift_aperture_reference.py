"""Frozen redshift-only drift of the cold kinematic aperture on actual archives."""
import json
import numpy as np
from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_snapshot_gauge_ledger_reference import cell_J,photon
from early_pole_decay_mode_clock_reference import G,KAPPA,intervals


def one(grid,n):
    d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
    h=[(np.array(d[f'endpoint_{k}_momentum_edges']),np.array(d[f'endpoint_{k}_occupations']))
       for k in ('parent','fermion')]
    old=np.array(d['scalar_momentum_edges'])
    N0=np.array(d['endpoint_scalar_signed_shell_number_over_T3'])
    edges=np.geomspace(1e-6,1e-5,n+1)
    q=np.cbrt((edges[1:]**3+edges[:-1]**3)/2)
    V=np.diff(edges**3)/(6*np.pi**2)
    cap=np.maximum(0.,np.minimum(old[1:,None]**3,edges[1:]**3)-
        np.maximum(old[:-1,None]**3,edges[:-1]**3))/np.diff(old**3)[:,None]
    f=(N0@cap)/V
    ms=d['checkpoints'][-1]['scalar_pole_over_T']
    Jg=photon(1.,0.,.5)
    out=[]
    for t in (0.,8e-5,.0003,.0008,.002,.004,.008,.015,.025,.039):
        x=.05*np.exp(t)
        M,B=np.sqrt(np.array([x,.06*x])**2+.005**2/np.pi**2*(Jg+np.r_[
            cell_J(*h[0],x,.45),cell_J(*h[1],.06*x,.45)]))
        ep=np.hypot(h[0][0],M)
        eb=np.hypot(h[1][0],B)
        gain=[]
        pcent=[]
        aw=[]
        probe=None
        for z,S,v,fs in zip(q,np.hypot(q,ms),V,f,strict=True):
            _,lo,hi,_,_=intervals(M,ms,B,z)
            e=np.unique(np.r_[lo,hi,ep[(ep>lo)&(ep<hi)],(eb+S)[(eb+S>lo)&(eb+S<hi)]])
            mid=(e[1:]+e[:-1])/2
            a=h[0][1][np.searchsorted(ep,mid)-1]
            b=h[1][1][np.searchsorted(eb+S,mid)-1]
            flux=v*G*KAPPA*x*((1.06*x)**2-ms**2)/(S*z)*(np.diff(e)@(a*(1-b)+fs*(a-b)))
            central=(lo+hi)/2
            pp=np.sqrt(central**2-M*M)
            half=(np.sqrt(hi*hi-M*M)-np.sqrt(lo*lo-M*M))/2
            gain.append(flux)
            pcent.append(pp)
            aw.append(half)
        j=np.argmin(abs(np.log(q/(4.1e-6))))
        out.append(dict(eta=t,x=x,vacuum_reference_total_gain=sum(gain),gains=gain,
                         sample_q=q[j],sample_parent_p_center=pcent[j],sample_parent_halfwidth=aw[j]))
    p=out[0]['sample_parent_p_center']
    w=out[0]['sample_parent_halfwidth']
    for r in out:
        r['sample_center_shift_in_initial_halfspans']=(r['sample_parent_p_center']-p)/w
        r['relative_frozen_total_gain']=r['vacuum_reference_total_gain']/out[0]['vacuum_reference_total_gain']
    return dict(grid=grid,n=n,rows=out)


def run():
    return dict(scope='Hold all comoving occupations and scalar conformal gap fixed; redshift changes only '
                      'the hard bare energies and reciprocal moments. Not an evolved spectrum or an actual gate.',
                rows=[one(grid,n) for grid in ('fine','third') for n in (8,16)])


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
