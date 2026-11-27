"""Pauli-free heavy regeneration on central, imposed direct-lepton spectra."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_heavy_recycling_gate_reference import angular_integral
from archived_scalar_pair_reaction_reference import amplitude,gross
from empty_lepton_pair_spectrum_reference import spectrum
from sm_dark_entropy_exchange_reference import MPL, StandardModel
from three_species_decay_onsager_reference import M,m,mp,y,Es
from warm_same_charge_heavy_recycling_reference import exact_angular,folded


def cells(edges,density):
    p=(edges[:-1]+edges[1:])/2
    w=np.array(density)
    assert abs(sum(w)-1.)<1e-8
    return p[w>0],w[w>0]


def rates(p,n,inc,angle,want_energy=False):
    e=np.hypot(p,inc)
    X,Z=roots_legendre(angle)
    prod=p[:,None]*p[None,:]
    E=e[:,None]*e[None,:]
    c=np.clip((E+inc*inc-2*M*M)/prod,-1.,1.)
    mu=(c-1)[:,:,None]/2+(c+1)[:,:,None]*X/2
    S=2*inc*inc+2*(E[:,:,None]-prod[:,:,None]*mu)
    k=np.sqrt(np.maximum(S/4-M*M,0.))
    pair=n[:,None]*n[None,:]*(c+1)/4
    base=pair[:,:,None]*Z*k/(np.sqrt(S)*np.pi*E[:,:,None])
    Etotal=e[:,None,None]+e[None,:,None]
    if inc==m:
        # The input is summed over the four stable spin/charge states.
        # Opposite plus both same-charge event rates. Their charge and
        # incoming identical factors give the same 1/4 multiplier each.
        pre=y**4/128
        r=np.stack((base*angular_integral(S),base*exact_angular(S)))
        pairs=pre*np.sum(r,axis=(1,2,3))
        return (pairs,pre*np.sum(r*Etotal,axis=(1,2,3))) if want_energy else pairs
    # Final heavy spins and the initial identical scalar factors are
    # exactly those of archived_scalar_pair_reaction_reference.gross.
    # For the outgoing anisotropy integrate the spin-averaged reverse.
    masked=k>0
    s=S[masked]
    w=base[masked]
    value=0.
    weighted_energy=0.
    q=np.broadcast_to(Etotal,base.shape)[masked]
    for i in range(0,len(s),10000):
        ai=amplitude(s[i:i+10000,None],X,light=M,heavy=m,scalar=mp)@Z
        value+=w[i:i+10000]@ai
        weighted_energy+=(w[i:i+10000]*q[i:i+10000])@ai
    a=np.array([value/16])
    return (a,np.array([weighted_energy/16])) if want_energy else a


def check_normalization(d,T):
    edge=np.array(d["scalar_differential_bins"])*T
    cen=(edge[:-1]+edge[1:])/2
    ws=np.asarray(d["stable_normalized_counts"])
    wf=np.asarray(d["scalar_normalized_counts"])
    def occ(weights,mass,scale):
        f=weights*scale*(2*np.pi**2)/(np.diff(edge)*cen**2)
        def get(E):
            p=np.sqrt(np.maximum(0.,(E-mass)*(E+mass)))
            j=np.clip(np.searchsorted(edge,p)-1,0,len(f)-1)
            return f[j]
        return get
    u,s=cells(edge,ws)
    v,w=cells(edge,wf)
    actual=rates(u,s,m,23)
    target=np.asarray(folded(edge,occ(ws,m,1/4),1,23))
    control_st=float(max(abs(actual/target-1)))
    _,ev,_=gross(edge,mp,M,occ(wf,mp,1),np.zeros_like,True,1,23,30,23,
                  split_final=False,exchange_mass=m)
    control_sc=float(abs(rates(v,w,mp,23)[0]/ev[0]-1))
    assert max(control_st,control_sc)<1e-3 and control_st<5e-13
    return dict(stable_charge_kernel_relative_defect=control_st,
                independent_scalar_angle_relative_defect=control_sc)


def run(eos,cache):
    cache.mkdir(parents=True,exist_ok=True)
    b=StandardModel(eos)
    rows=[]
    for tr,tops in ((.05,(.047,.043,.035,.025,.015)),(.1,(.095,.085,.075,.06,.05,.035,.021,.012))):
        for t in tops:
            control=[]
            for nh,nout,order,a in ((97,17,(9,24,16),14),(171,26,(12,36,24),24)):
                name=cache/f"{tr:.3f}-{t:.3f}-{nh}.json"
                if name.exists():
                    s=json.loads(name.read_text())
                else:
                    s=spectrum(b,tr,t,order=order,nout=nout,nhist=nh)
                    name.write_text(json.dumps(s,allow_nan=False)+"\n")
                rho,_,ent,_=b.state(t)
                H=np.sqrt(rho/3)/MPL
                edge=np.asarray(s["scalar_differential_bins"])*t
                p,u=cells(edge,s["stable_normalized_counts"])
                v,w=cells(edge,s["scalar_normalized_counts"])
                nst=sum(s["stable_yields"])*ent
                nsc=s["stable_yields"][1]*ent
                st,st_e=rates(p,u*nst,m,a,True)
                sc,sc_e=rates(v,w*nsc,mp,a,True)
                fill=np.r_[st,sc]
                # With unpolarized prompt decays each scalar retains the
                # mean energy fraction E*_scalar/M of its parent. Count
                # gain alone does not preserve the incident spectra.
                D=float((Es/M)*st_e.sum()-(1-Es/M)*sc_e.sum())
                rho_l=nst*(u@np.hypot(p,m))
                rho_s=nsc*(w@np.hypot(v,mp))
                assert np.all(fill>=0.)
                control.append(dict(histogram_count=nh,source_time_order=nout,pair_angle_order=a,
                    scalar_yield=nsc/ent,stable_yield=nst/ent,
                    opposite_same_scalar_regeneration_events_per_H_per_entropy=(fill/(H*ent)).tolist(),
                    weighted_incident_energy_per_H_per_entropy=(np.r_[st_e,sc_e]/(H*ent)).tolist(),
                    net_scalar_energy_transfer_per_H_per_entropy_GeV=D/(H*ent),
                    net_energy_transfer_over_passive_scalar_energy_per_H=D/(H*rho_s),
                    input_processing_over_passive_stable_number_per_H=2*st.sum()/(H*nst),
                    input_processing_over_passive_scalar_number_per_H=2*sc.sum()/(H*nsc),
                    relative_scalar_energy_density=rho_s/rho,
                    relative_stable_energy_density=rho_l/rho,
                    extra_Q_per_H_per_entropy=2*fill.sum()/(H*ent)))
            rows.append(dict(TR_GeV=tr,T_GeV=t,controls=control))
    norm=check_normalization(json.loads((cache/"0.100-0.060-97.json").read_text()),.060)
    return dict(normalization_controls=norm,scope="Gross Pauli-free heavy regeneration at labeled epochs on controlled, forced vacuum lepton "
                      "spectra, with all heavy parents decayed promptly. Bin midpoint projection and two quadratures. "
                      "The populations do not respond; not a time-integrated count, not a coupled or upper relic result, "
                      "and does not include hadronic sources or thermal feedback.",rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("cache",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.cache),indent=2,allow_nan=False))
