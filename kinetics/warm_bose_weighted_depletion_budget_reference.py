"""Turn the conditional fermion floor into a finite-duration Bose-pulse target."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from late_radiation_background_reference import MPL,TabulatedRadiation
from warm_fermion_survival_with_leptons_reference import subpole


def read(p):
    return json.loads(Path(p).read_bytes())


def window_time(eos,T,a,b):
    R=TabulatedRadiation(eos,T)
    l,h=np.log(R.state(b)[0]),np.log(R.state(a)[0])
    e=R.sg.x[(R.sg.x>l)&(R.sg.x<h)]
    u=np.r_[l,e,h]
    s,g=R.sg(u),R.rg(u)
    Hph=np.pi*T*T/(np.sqrt(45)*MPL)
    g0=float(R.sg(R.initial))
    f=np.exp(.5*(np.minimum(g[:-1],g[1:])-np.log(2))+
             (2/3)*(g0-np.maximum(s[:-1],s[1:]))-1e-12)
    t=(R.s0-3*u-s)/3
    assert abs(t[0]-np.log(b))<1e-12 and abs(t[-1]-np.log(a))<1e-12
    return float(np.sum(np.diff(np.exp(-t))/f)/Hph)


def run(eos,prior,leptons,isolated):
    P,L,A=map(read,(prior,leptons,isolated))
    K=L["coefficient"]; goal=P["target"]
    cuts=np.linspace(2,9.9,80)
    cutK=subpole(cuts)
    out=[]
    for p,l,a in zip(P["cases"],L["cases"],A["cases"],strict=True):
        T,N,om = p["anchor_T_GeV"],p["core_start_fermion_number_GeV3"],p["conservative_initial_fermion_Omega"]
        H=p["integration_time_energy_assisted"]
        assert abs(a["T_GeV"]/T-1.)<1e-12 and l["q"]==p["q"]
        loss=2*l["available_energy_upper_GeV4"]/(cuts*N)
        def lower(total_B):
            return float(om*max(0,np.max(1/(1+.5*N*(K*total_B+H*cutK))-loss)))
        target=brentq(lambda w: lower(w*H)-goal,0.,1e5)
        assert abs(np.sqrt(target)-1-l["uniform_F_cap_min_to_permit_012_vacuum_leptons"])<1e-10
        rows=[]
        for lo,hi in ((1.,1.12),(1.,1.5),(1.,3.),(1.5,3.),(3.,10.)):
            J=window_time(eos,T,lo,hi)
            share=J/H
            assert 0<share<1
            f={}
            for outer in (.1,1.):
                ref=(1+outer)**2
                peak=np.sqrt(ref+(target-ref)/share)-1
                assert abs(lower(H*ref+J*((1+peak)**2-ref))-goal)<2e-12
                f[str(outer)]=float(peak)
            rows.append(dict(relative_b=[lo,hi],weighted_dose_share_using_upper_integrals=share,
                             necessary_burst_peak_allowed_at_reference_total=f))
        out.append(dict(q=p["q"],min_weighted_bose_enhancement_relative_to_I_upper=target,
                        windows=rows))
    return dict(scope="Necessary weighted Bose losses in the specified vacuum charged-lepton "
                      "network. Caps outside and inside time windows are conditional, "
                      "not model-derived. Windows give alternatives, not simultaneous budgets.",
                target=goal,cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for s in ("eos","prior","leptons","isolated"):
        p.add_argument(s)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.leptons,a.isolated),indent=2,allow_nan=False))
