"""Add the specified neutral-charge pion-pair visible sink and source."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import gammaincc

from kloe_pion_band_portal_reference import POINT_MUON, table
from late_radiation_background_reference import MPL,TabulatedRadiation
from timelike_portal_onsager_coverage_reference import mpi
from warm_extra_injection_headroom_reference import pion_positive_Q
from warm_late_visible_recycling_envelope_reference import energy_rate_coefficient
from warm_pion_bose_and_threshold_reference import bose,low_strip_per_unit
from warm_signed_fermion_time_comparison_reference import comparison,read,setup
from warm_fermion_survival_with_leptons_reference import subpole


def pion_upper(eos,T,data,n,angles):
    R=TabulatedRadiation(eos,T)
    # Same thermodynamic floor as in the late lepton allowance.
    bstop=500.
    l,h=np.log(R.state(bstop)[0]),R.initial
    knots=R.sg.x[(R.sg.x>l)&(R.sg.x<h)]
    u=np.unique(np.r_[np.linspace(l,h,n+1),knots])[::-1]
    lgs,lgr=R.sg(u),R.rg(u)
    t=(R.s0-3*u-lgs)/3
    assert np.all(np.diff(t)>0) and abs(t[-1]-np.log(bstop))<1e-12
    hrad=np.exp(.5*(np.minimum(lgr[:-1],lgr[1:])-np.log(2))+
                (2/3)*(lgs[0]-np.maximum(lgs[:-1],lgs[1:]))-1e-12)
    Hph=np.pi*T*T/(np.sqrt(45)*MPL)
    dt=np.diff(np.exp(5*t))/(5*Hph*hrad)
    temp=np.exp(u[:-1])
    Q=np.array([(pion_positive_Q(z,data,angles)[0]+4*low_strip_per_unit(z)[0])*bose(z)**2
                for z in temp])
    assert np.min(Q)>=0 and np.all(np.diff(Q)<=1e-200)
    grids,sig,err=data
    # This is a prescribed effective two-massive-pion kernel, with no
    # separate photons. Use a single global R cap for reverse pair loss.
    R_table=max((sig+err)*(grids+.005)/POINT_MUON)
    R_below=(1-4*mpi**2/.10)**1.5  # 4 times the pointlike 1/4 factor.
    Rmax=max(R_table,R_below)*bose(T)**2
    scale=np.exp((lgs[0]-np.log(2))/3+1e-12)*T
    coeff=energy_rate_coefficient(T)[0]*(Rmax/3)
    # Cold continuation: start with the threshold-free T^9 moment,
    # but retain the true massive pion threshold on the Bessel tail.
    # K_nu(u) <= sqrt(pi/198) exp(-100) (100/u)^8 for u>=100.
    uc=2*mpi/(scale/bstop)
    assert uc>100
    # Keep the sharper threshold tail over the total invariant integral:
    # u^6(2K1+uK0) <= A u^-2 (2+u) would have a logarithmic tail,
    # so use a stronger power 16 instead, valid for every u>=100.
    A=np.exp(-100)*np.sqrt(np.pi/198)*100**16
    tail_shape=A*(2*uc**-9/9+uc**-8/8)/3072
    late=coeff*scale**9/(4*Hph*bstop**4)*tail_shape
    e=dt@Q+late
    return float(e),float(Rmax)


def run(eos,kloe,F,L,I,V):
    F,L,I,V=map(read,(F,L,I,V))
    data=table(kloe); cases=[]
    K=L["coefficient"]
    for f,l,i,v in zip(F["cases"],L["cases"],I["cases"],V["cases"],strict=True):
        q=f["q"]
        assert q==l["q"]==i["q"]==v["q"]
        T=f["anchor_T_GeV"]; N=f["core_start_fermion_number_GeV3"]; O=f["conservative_initial_fermion_Omega"]
        E=l["available_energy_upper_GeV4"]
        B=bose(T)**2
        J=4*B*gammaincc(5,mpi/T)*v["elastic_integrating_factor_exponent"]/6
        c=i["bounds"][1]["clock"]
        controls=[]
        for n,ang,nstep in ((400,8,1000),(1000,12,2800)):
            P,r=pion_upper(eos,T,data,n,ang)
            Eplus=(E+P)*np.exp(J)
            W,A,S,H,Z,eta,cut=setup(N,Eplus,K,c,nstep)
            # Actual stated window is below both the mediator pole and
            # every lepton split. Also compare granting this R below
            # every split, deliberately not using the visible support.
            end=max(data[0])+.005
            assert end<min(cut.ravel())**2
            Kpi=subpole(np.sqrt(end))*r/3
            tied=O*comparison((W,A+.5*N*Kpi,S,H,Z,eta,cut),.5)[0]
            A=A+.5*N*subpole(cut)*r/3
            x=O*comparison((W,A,S,H,Z,eta,cut),.5)[0]
            lepton_only=O*comparison(setup(N,E,K,c,nstep),.5)[0]
            assert F["target"]<x<=tied<lepton_only
            # Independent stress tests of the not-yet-included source budget,
            # holding the chosen subpole R allowance fixed.
            stress=[]
            for increase in (.02,.05,.10,.20):
                a,b,d,e,z,g,t=setup(N,Eplus*(1+increase),K,c,nstep)
                b=b+.5*N*subpole(t)*r/3
                val=O*comparison((a,b,d,e,z,g,t),.5)[0]
                stress.append(dict(extra_gross_energy_fraction=increase,lower_Omega_h2=val))
            controls.append(dict(time_cells=nstep,pion_injection_GeV4=P,pion_injection_relative_to_previous_energy_upper=P/E,
                                 pion_extra_heat_exponent=float(J),effective_global_pion_R_cap=r,
                                 new_available_energy_upper_GeV4=float(Eplus),
                                 pion_pair_sigma_v_cap_with_stated_support=float(Kpi),
                                 with_charged_pions_lower_Omega_h2=tied,
                                 weaker_unspecific_subpole_R_lower_Omega_h2=x,
                                 no_pion_control_lower_Omega_h2=lepton_only,
                                 sensitivity_to_extra_energy=stress))
        cases.append(dict(q=q,fermion_per_spin_cap=.5,controls=controls))
    return dict(scope="Specified post-10.442 MeV homogeneous isotropic warm signed network, "
                      "fermion occupancies per spin <= 1/2. Retains vacuum charged leptons, "
                      "adds the effective massive charged pion pair in the selected KLOE "
                      "window plus a bounded threshold strip and capped spacelike heating. "
                      "Not the full hadronic or radiative reverse theory.",cases=cases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for t in ("eos","kloe","fermion","leptons","isolated","visible"):
        p.add_argument(t)
    a=p.parse_args()
    out=run(a.eos,a.kloe,a.fermion,a.leptons,a.isolated,a.visible)
    out["input_sha256"]=[hashlib.sha256(Path(x).read_bytes()).hexdigest() for x in vars(a).values()]
    print(json.dumps(out,indent=2,allow_nan=False))
