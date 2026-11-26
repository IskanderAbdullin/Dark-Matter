"""Fixed-sign EdS forcing control for additive nonrelativistic tracer response."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import roots_legendre

from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import mp


def run(eos,files):
    bg=Background(eos)
    Hm=np.sqrt(bg.om)/bg.D0 # Mpc^-1
    K=np.array([.3,1.,2.,4.,8.,10.,20.])
    waves=K*bg.H0/100
    raw=[np.loadtxt(p).T for p in files]
    lo=min(np.log(q[0]) for q,_ in raw)
    hi0=max(q[-1] for q,_ in raw)
    tab=[CubicSpline(q,f,extrapolate=False) for q,f in raw]
    Nfull=np.array([np.trapezoid(f*q*q,q) for q,f in raw])
    # After a_in, Newtonian tracer response to a constant potential in
    # EdS is a weighted j0 integral, normalized to the zero-dispersion
    # response to the same source and with the same zero initial data.
    z=20
    final=1/(1+z)
    allout=[]
    for a0 in (.001,.005,.01):
        # Common physical speed cutoff avoids applying the Newtonian
        # formula to tiny but formally extended lognormal tails.
        cut=min(hi0,.10*mp*a0/bg.T0)
        x=np.exp(np.linspace(lo,np.log(cut),3100))
        w=np.array([np.maximum(0.,np.nan_to_num(F(x),nan=0.))*x*x for F in tab])
        N=np.trapezoid(w,x,axis=1)
        q0=np.linspace(lo,np.log(max(hi0,cut)),5000)
        x0=np.exp(q0)
        tail_use=x0>cut
        v0=x0[tail_use]
        omitted=np.array([np.trapezoid(np.maximum(0,np.nan_to_num(F(v0),nan=0.))*v0*v0,v0) for F in tab])/Nfull
        assert omitted.max()<1e-7
        vmax=x[-1]*bg.T0/(mp*a0)
        assert vmax<=.1000000001
        uv=np.linspace(np.log(a0),np.log(final),25)
        rows=[]
        for n in (7,12):
            t,s=roots_legendre(n)
            a=np.exp(((uv[:-1,None]+uv[1:,None])+np.diff(uv)[:,None]*t)/2).ravel()
            da=a*(np.diff(uv)[:,None]*s/2).ravel()
            travel=(2/Hm)*(1/np.sqrt(a)-1/np.sqrt(final))
            wt=da*(1-np.sqrt(a/final))
            assert np.all(wt>=0.)
            norm=final/3.-a0+2*a0**1.5/(3*np.sqrt(final))
            assert abs(wt.sum()/norm-1.)<2e-13
            v=[]
            for k in waves:
                kern=np.sinc((k*bg.T0/(np.pi*mp))*x[:,None]*travel[None,:])
                I=kern@wt/norm
                v.append(np.trapezoid(w*I,x,axis=1)/N)
            v=np.asarray(v).T
            add=np.average(v[1:],axis=0,weights=N[1:])
            err=np.max(abs(v[0]-add))
            assert err<3e-7 and np.min(v)>-.218 and np.max(v)<=1.+1e-12
            rows.append((v,add,err))
        assert np.max(abs(rows[0][0]-rows[1][0]))<2e-4
        q=rows[-1]
        allout.append(dict(a_start=a0,a_final=final,physical_p_over_m_max=vmax,omitted_number_fraction=omitted.tolist(),k_h_per_Mpc=K.tolist(),components=q[0].tolist(),
             weighted_prompt_recoil=q[1].tolist(),absolute_mix_residual=float(q[2]),
             max_temporal_quadrature_defect=float(np.max(abs(rows[0][0]-rows[1][0])))))
    # Exact q-integrals for the same physical pressure moments; these
    # help diagnose the effect of an effective common sound speed.
    moments=[]
    for a in (1e-6,3e-6,1e-5,1e-4,final):
        Q=[];P=[];S=[]
        for q,F in raw:
            e=np.hypot(q,a*mp/bg.T0)
            Q.append(np.trapezoid(q*q*F*e,q))
            P.append(np.trapezoid(q**4*F/(3*e),q))
            S.append(np.trapezoid(q**6*F/(3*e**3),q))
        Q,P,S=map(np.asarray,(Q,P,S))
        adj=(5*P-S)/(3*(Q+P))
        err=abs(adj[0]/np.average(adj[1:],weights=(Q+P)[1:])-1.)
        assert err<5e-7
        moments.append(dict(a=a,recoil_energy_fraction=float(Q[2]/(Q[1]+Q[2])),
             recoil_pressure_fraction=float(P[2]/(P[1]+P[2])),
             rho_active_adiabatic_speed_squared=adj.tolist(),
             enthalpy_weighted_adiabatic_speed_relative_error=float(err)))
    return dict(components=["joined","prompt","spontaneous_recoil"],cases=allout,background_pressure_moments=moments,horizon_scope="z=20 endpoint; Newtonian subhorizon and nonrelativistic control with EdS clock and prescribed constant same-sign potential, zero initial perturbations and a CDM response with matching zero initial conditions. Gravitationally forced tracer Green function, not the physical CLASS redshift evolution or observational transfer.",
      input_sha256=[hashlib.sha256(p.read_bytes()).hexdigest() for p in (eos,*files)])


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","joined","prompt","recoil"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,[a.joined,a.prompt,a.recoil]),indent=2,allow_nan=False))
