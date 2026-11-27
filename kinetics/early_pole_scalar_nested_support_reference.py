"""Nested closed scalar supports with exactly matched interior decay probes."""
import argparse
import json
from types import SimpleNamespace

import numpy as np
from scipy.integrate import solve_ivp

from early_pole_live_aperture_witness_reprojection_reference import initial_projection
from early_pole_native_aperture_witness_reference import PROFILES
from early_pole_scalar_cold_screened_flux_reference import Screen
from early_pole_scalar_elastic_stencil_reference import Elastic


def run(grid,bins,order,duration):
    if bins<4 or order<2 or not duration>0.:
        raise ValueError("Invalid reference controls.")
    profile=next(p for p in json.loads(PROFILES.read_text()) if p['grid']==grid)
    C,a,Y0,Y1,*_=initial_projection(profile,3,bins,duration)
    central=np.geomspace(1e-6,1e-5,bins+1)
    old=np.asarray(a['scalar_momentum_edges'])**3
    ns=np.maximum(0.,np.asarray(a['endpoint_scalar_signed_shell_number_over_T3']))
    H=np.sqrt(np.pi**2*106.75/90)*(86./a['x_end'])/2.435e18
    out=[]
    for wing in (0,3,6):
        r=10.**(1/bins)
        p=np.r_[central[0]*r**np.arange(-wing,0),central,central[-1]*r**np.arange(1,wing+1)]
        c=p**3
        ov=np.maximum(0.,np.minimum(old[1:,None],c[None,1:])-np.maximum(old[:-1,None],c[None,:-1]))
        W=ov/np.diff(old)[:,None]
        N0=ns@W
        C.n=len(N0);C.ns=N0
        C.spect=ns*(1.-W.sum(axis=1))
        C.q=np.cbrt((c[1:]+c[:-1])/2)
        C.sw=np.diff(c)/(6*np.pi**2)
        C.lambdas=[SimpleNamespace(n=C.M.decay.n,cells=C.M.decay.cells,q=q,nod=C.M.z,w=C.M.w) for q in C.q]
        assert abs(N0.sum()+C.spect.sum()-ns.sum())<1e-15 and min(N0)>0
        _,e0,s0,f,m=C.scalar(N0)
        K=Elastic(C.q,m,p,order,C.lam,Screen(p,f,m,C.lam,order))
        anchor=np.argmax(N0)
        ids=np.flatnonzero(np.arange(C.n)!=anchor)
        number=N0.sum()
        def counts(log):
            n=N0.copy()
            n[ids]*=np.exp(log)
            n[anchor]=number-n[ids].sum()
            if n[anchor]<=0:
                raise ValueError("Nonpositive fixed-count coordinate.")
            return n
        def rate(_,log):
            n=counts(log)
            return (K.flux(n/K.vol)[0]/(H*n))[ids]
        start=np.zeros(len(ids))
        first=min(duration/100,.01/max(1.,np.max(abs(rate(0.,start)))))
        time=np.r_[0.,np.geomspace(duration/128,duration,11)]
        h=solve_ivp(rate,(0.,duration),start,method='Radau',
                    rtol=2e-9,atol=3e-11,t_eval=time,first_step=first,max_step=duration/12)
        assert h.success,h.message
        central_slice=slice(wing,wing+bins)
        rows=[]
        for t,y in zip(time,h.y.T,strict=True):
            n=counts(y)
            es,rho,s,_,ms=C.scalar(n)
            state=[]
            for a0 in (0.,.97):
                Y=(1-a0)*Y0+a0*Y1
                z=np.r_[Y.ravel(),n,0.,0.]
                out_source=C.rhs(z,False)[C.j0:C.j0+C.n]
                state.append(float(out_source[central_slice].sum()))
            assert abs(n.sum()-number)<1e-15 and abs(rho-e0)<1e-13
            rows.append(dict(exposure=t,initial_count_fraction_in_active_range=number/ns.sum(),
                             central_fraction_change=(n[central_slice].sum()-N0[central_slice].sum())/N0[central_slice].sum(),
                             local_shape_l1=float(np.sum(abs(n[central_slice]-N0[central_slice]))/N0[central_slice].sum()),
                             total_scalar_energy_change=rho-e0,entropy_change=s-s0,scalar_pole=ms,
                             matched_central_step_and_witness_source=state))
        assert min(np.diff([r['entropy_change'] for r in rows]))>-1e-27
        out.append(dict(extra_bins_each_side=wing,momentum_support=[p[0],p[-1]],rows=rows,excluded_graph_geometry=K.raw_escape))
    return dict(grid=grid,central_bins=bins,order=order,exposure=duration,cases=out,
                scope="Held scalar-only shape proxy on nested closed supports; hard occupations frozen.")


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--bins',type=int,default=16)
    p.add_argument('--order',type=int,default=6)
    p.add_argument('--duration',type=float,default=2e-7)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
