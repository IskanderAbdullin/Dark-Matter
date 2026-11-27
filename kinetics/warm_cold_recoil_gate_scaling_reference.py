"""Resolve the local recoil gate by following its positive angular collar."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import Collision,kin
from joint_visible_heavy_adiabatic_reference import angle_panels,snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp,width
from warm_light_occupation_stage_audit_reference import digest,read
from warm_near_rest_light_inverse_decay_reference import E0,Q0,window
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,source):
    A=read(source)
    S=A["snapshot"]
    assert A["snapshot_sha256"]==digest(S)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    edge=p[1]
    keys=["entropy_variable_"+j for j in ("scalar","fermion")]
    state=[np.asarray(S[k]) for k in keys]
    B,be,F,fe=snapshot_functions(S)
    rho=get_mom(p,state)[2].sum()
    H=np.sqrt((rho+import_restart(S["positive_exterior_restart"]).energy(1.)+
               TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1])/3)/MPL
    threshold=window(edge)[-2]

    def crop(cap):
        upper=np.sqrt((cap-mp)*(cap+mp))
        q=np.r_[p[p<upper],upper]
        d=len(q)-1
        st=[]
        for x,mass in zip(state,(mp,m),strict=True):
            y=x[:,:d].copy()
            y[1,-1]*=(kin(upper,mass)-kin(q[-2],mass))/(kin(p[d],mass)-kin(q[-2],mass))
            st.append(y)
        return q,st

    def react(phi,chi,cap,n):
        q=np.sqrt((cap-mp)*(cap+mp))
        fc=np.hypot(q,m)
        bs=np.r_[be[be<cap],cap]
        fs=np.r_[fe[fe<fc],fc]
        E=(phi+chi).ravel()
        w,weight=angle_panels(E,bs,fs,angle=n)
        y=E[:,None,None]-w
        mask=(w<=cap)&(y<=fc)
        f=B(w); g=F(y)
        a=(weight*f*g*mask).sum(axis=(1,2))
        b=(weight*(1+f)*(1-g)*mask).sum(axis=(1,2))
        h=np.divide(a,a+b,out=np.zeros_like(a),where=a+b>0.).reshape(phi.shape)
        scalar=B(phi); light=F(chi)
        return h*(1-light)*(1+scalar)-(1-h)*light*scalar

    e,_,d,_,_=window(edge)
    slope=Q0/m-(E0/m)*edge/e
    factor=(M/m)**3*m/e*width/H /(2*d)
    # The first shell owns the left limit of the entropy profile.
    # The next cell can have a distinct fitted endpoint value.
    integrand=float(react(np.array([[threshold]]),np.array([[np.nextafter(e,0.)]]),threshold,16)[0,0])
    leading=3/edge*factor/(2*slope)*integrand
    assert np.isfinite(leading) and leading!=0.

    def direct(cap,radial,angular):
        if cap<=threshold:
            return 0.
        root=brentq(lambda k:window(k)[-2]-cap,0.,edge,xtol=1e-18,rtol=1e-14)
        z,w=roots_legendre(radial)
        t=root+(edge-root)*(1+z)/2
        E,_,d,Lo,Hi=window(t)
        z1,w1=roots_legendre(angular)
        span=(np.minimum(cap,Hi)-Lo)/2
        x=Lo[:,None]+span[:,None]*(z1+1)
        f=react(x,E[:,None],cap,angular)
        packed=(M/m)**3*m/E*width/H/(2*d)*span*(f@w1)
        return float(3/edge**3*(edge-root)/2*((t*t*w)@packed))

    out=[]
    for margin in (-5e-7,0.,2e-4,6e-5,2e-5,6e-6,2e-6,6e-7,2e-7,6e-8,2e-8):
        cap=threshold+margin
        q,th=crop(cap)
        angles=[]
        for n in (3,6,12):
            kernel=Collision(q,n,n)
            val=kernel.rate(th)[1,0,0]/kernel.c[1].vol[0]/H
            angles.append([n,int(np.sum(kernel.J==0)),float(val)])
        independent=[direct(cap,n,n) for n in (8,16,24)]
        out.append(dict(E_margin_GeV=margin,cap_GeV=cap,native_orders_events_rate=angles,
                        independent_rates=independent,
                        finest_to_asymptotic=None if margin<=0. else independent[-1]/(leading*margin**2)))
    assert all(not any(z[1] for z in r["native_orders_events_rate"]) for r in out if r["E_margin_GeV"]<0.)
    assert abs(out[-1]["finest_to_asymptotic"]-1.)<5e-5
    assert max(abs(z[2]) for z in out[1]["native_orders_events_rate"])<2e-29
    return dict(scope="Fixed archived shape cropped at a variable common edge; cold-cell geometry only. "
                      "The Hubble normalization and retained profile do not evolve during this test.",
                snapshot_sha256=A["snapshot_sha256"],T_GeV=S["T_visible_GeV"],threshold_GeV=threshold,
                leading_rate_per_Hubble_over_margin_squared_GeV_minus2=leading,cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("source",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.source),indent=2,allow_nan=False))
