"""Audit the finite pad by allowing vacuum exterior daughter emission on frozen phases."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_visible_heavy_adiabatic_reference import angle_panels,snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import Es,M,P,m,mp,width
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def outgoing(S,n):
    B,es,F,ef=snapshot_functions(S)
    lo=min(es[-1]+m,ef[-1]+mp)
    hi=es[-1]+ef[-1]
    def gate(x,d,e):
        p=np.sqrt((x-d)*(x+d))
        u=e*x+p*P
        return np.r_[M*(e*e+p*p)/u,M*u/d**2]
    edge=np.unique(np.clip(np.r_[lo,hi,gate(es,mp,Es),gate(ef,m,M-Es),
                                    (es[:,None]+ef).ravel()],lo,hi))
    w0=np.diff(edge)/2
    x,w=roots_legendre(n)
    E=((edge[1:]+edge[:-1])[:,None]/2+w0[:,None]*x).ravel()
    dE=(w0[:,None]*w).ravel()
    out=np.zeros((2,2))
    for i in range(0,len(E),100):
        e=E[i:i+100]
        a,b=angle_panels(e,es,ef,angle=n)
        f=e[:,None,None]-a
        sb,ff=B(a),F(f)
        u=b*sb*ff
        v=b*(1+sb)*(1-ff)
        D=(u+v).sum(axis=(1,2))
        h=np.divide(u.sum(axis=(1,2)),D,out=np.zeros_like(D),where=D>0.)
        coeff=2/np.pi**2*width*M*dE[i:i+100]*np.sqrt((e-M)*(e+M))*h
        for k,(d,r) in enumerate(((a,es[-1]),(f,ef[-1]))):
            g=v*(d>r)
            out[:,k]+=[coeff@g.sum(axis=(1,2)),coeff@(g*d).sum(axis=(1,2))]
    assert np.all(out>=0.)
    return out


def run(eos,prior,archives):
    root=read(prior)["branches"][-1]
    cut=np.hypot(M,1./root["entropy_expansion_from_original_anchor"])
    result=[]
    for p in archives:
        x=read(p)
        S=x["snapshot"]
        assert x["hash"]==digest(S)
        q=np.asarray(S["lab_momentum_edges_GeV"])
        rho=get_mom(q,[S["entropy_variable_"+z] for z in ("scalar","fermion")])[2].sum()
        H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+rho)/3)/MPL
        rates=[]
        for n in (6,10,14):
            a=outgoing(S,n)
            scale=M*a[0]/(rho*H)
            loss=a[1]/(rho*H)
            b=np.array([integrate(S,cut,n,c)[3] for c in (True,False)])/(rho*H)
            upper_energy=max(np.hypot(m,q[-1]),np.hypot(mp,q[-1]))
            diag=integrate(S,upper_energy,n,False)/(rho*H)
            assert diag[4]==0. and diag[5]<=loss.sum()*(1+2e-12)
            assert abs(diag[5]-diag[3])<1e-12*max(abs(diag[1]),1e-245)
            rates.append(dict(order=n,exterior_count_times_M_over_rhoH=scale.tolist(),
                              exterior_energy_over_rhoH=loss.tolist(),
                              independent_excess_emitted_above_both_pad_energies=float(diag[5]),
                              frozen_cut_flux_closed_and_open=b.tolist(),
                              fractional_change_of_frozen_signed_current=float((b[1]-b[0])/max(abs(b[0]),1e-280))))
        result.append(dict(hard_GeV=x["plan"]["hard_GeV"],time=x["time"],snapshot_sha256=x["hash"],
                           pad_momentum_GeV=float(q[-1]),controls=rates))
    return dict(scope="A frozen comparison on the same saved phase. The open branch changes the "
                      "stationary parent as well as outgoing support. Not independently evolved on a larger pad.",
                cut_GeV=cut,cases=result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("prior",type=Path);p.add_argument("archives",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.archives),indent=2,allow_nan=False))
