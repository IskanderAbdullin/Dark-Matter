"""Local forced hard-occupation and chemical feedback diagnostics."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_decay_mode_clock_reference import X
from early_pole_live_parent_energy_reference import state_moments
from early_pole_scalar_angular_chemical_coupled_reference import AngularChemical


def one(s,hot,power):
    nshape={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}
    d=AngularChemical(nshape[tuple(s['cells'])],s['mass_ratio'],1.,power)
    qstar=np.sqrt(3)*s["checkpoints"][-1]["scalar_pole_over_T"]
    s,meta=d.setup_hot_grid(s,hot,(.8*qstar,1.2*qstar))
    f=np.r_[s['endpoint_scalar_occupations'],s['endpoint_fermion_occupations'],
            s['endpoint_parent_occupations']]
    x=s['x_end']
    t=np.log(x/X)
    E,D,P,m=state_moments(d,x,f)
    a,h=d.collider(x,m)
    ix=a.active
    keep=~a.cold
    hard=ix[keep]
    assert min(f[ix])>0
    scale=np.sqrt(d.radiation/(d.radiation+D+P))
    d.strength=0.
    justdecay=d.rhs(t,f)
    d.strength=1.
    justchem=a.flux(f[ix])*scale/h
    # Derivative in physical cell counts or occupations has same eigenvalues.
    F=f[ix]
    lf=np.log(F)
    L=np.logaddexp(0.,-lf)
    common=a.A@lf
    behind=L[a.incoming].sum(axis=1)
    ahead=behind+a.S@L
    R,_=a.events(F)
    st=a.S[:,keep]
    # Hard nodes only on the outgoing side at the tracked snapshots.
    assert np.min(st.data)>=0
    fw=a.W*np.exp(common+ahead)
    B=st.multiply(R[:,None]).multiply((1/F[keep])[None,:])
    Df=st.multiply(fw[:,None]).multiply((-1./(F[keep]*(1+F[keep])))[None,:])
    J=(st.T@(B+Df)).toarray()*scale/(a.vol[keep,None]*h)
    lam=np.linalg.eigvals(J)
    g=float(lam.real.max())
    assert max(abs(lam.imag))<max(g,1)*1e-7
    assert g>0
    v=np.cos(np.arange(len(hard)))*F[keep]
    w=np.zeros(len(F));w[keep]=v
    tstep=2e-5
    df=((a.flux(F+tstep*w)-a.flux(F-tstep*w))[keep]*scale/(2*tstep*a.vol[keep]*h))
    err=max(abs(df-J@v))/max(abs(df))
    assert err<1e-6
    dq=1e-5
    qmean=np.sqrt(3)*m
    en=np.hypot(m,d.d.p)
    n=d.d.w*f[:d.ns]
    coupling=.5*d.d.lam
    mdot=(2*(1e-9*x)**2+coupling*(d.d.w*justdecay[:d.ns])@(1/en))/(2*m+coupling*m*(n@(1/en**3)))
    # Nominal crossing using an analytic full cold angular support width.
    dE=2*dq*np.sqrt(1-m*m/(4*(m*m+dq*dq)))
    lo=2*m-dE;hi=2*np.hypot(m,dq)+dE
    width=np.sqrt((hi-m)*(hi+m))-np.sqrt((lo-m)*(lo+m))
    dwell=width/max(abs(np.sqrt(3)*mdot),1e-99)
    return dict(x=x,mass_ratio=s['mass_ratio'],initial_pole_over_T=m,nominal_hot_q=qmean,
                hard_receiving_nodes=d.d.p[hard].tolist(),hard_occupation=F[keep].tolist(),
                modeled_initial_bare_net_count_loss_per_log_a=float(-justchem.sum()/n.sum()),
                modeled_bare_local_hard_growth_eigenvalue_per_log_a=g,angular_samples=2**power,
                assumed_no_chemical_window_drift_dq_per_log_a=float(np.sqrt(3)*mdot),
                cold_kinematic_support_width=width,nominal_crossing_time_in_log_a=float(dwell),
                formal_multiplier_for_one_hot_growth_efold_in_this_crossing=float(1/(g*dwell)),
                jacobian_directional_relative_error=float(err),**meta)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    a=p.parse_args()
    result={}
    for name in a.spectra:
        file=Path(name)
        raw=gzip.decompress(file.read_bytes()) if file.suffix=='.gz' else file.read_bytes()
        s=json.loads(raw)
        result[file.name]=[one(s,h,power) for h,power in ((0,3),(36,3),(100,3),(300,3),(100,5),(300,5),(300,7))]
    print(json.dumps(result,indent=2,allow_nan=False))
