"""Fold selected omitted slow rates over the actual short fast-return trajectory."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import trapezoid

from archived_scalar_pair_reaction_reference import gross
from joint_prompt_two_moment_evolution_reference import Cells,run
from sm_dark_entropy_exchange_reference import MPL,StandardModel
from three_species_decay_onsager_reference import Es,M,m,mp
from warm_prompt_reheating_heavy_gain_probe_reference import rates


def occ(cell,shape):
    shape=np.asarray(shape)
    def evaluate(e):
        p=np.sqrt(np.maximum((e-cell.mass)*(e+cell.mass),0.))
        ix=np.clip(np.searchsorted(cell.p,p,side="right")-1,0,len(cell.p)-2)
        return np.where(p<cell.p[-1],cell.f(shape,e,ix),0.)
    return evaluate


def sample(edge,shapes,t,b,order):
    phi,chi=from_shapes(edge,shapes)
    rho,_,s,_=b.state(t)
    H=np.sqrt(rho/3)/MPL
    nr,na=order
    J=np.array([0.,0.])
    inc_E=np.zeros(2)
    for i,(mass,g,f) in enumerate(((mp,1,phi),(m,4,chi))):
        x,w=roots_legendre(nr)
        a=(edge[1:]+edge[:-1])/2
        b1=np.diff(edge)/2
        p=a[:,None]+b1[:,None]*x
        weight=b1[:,None]*g*p*p*w*f(np.hypot(p,mass))/(2*np.pi**2)
        keep=weight.ravel()>0
        n=weight.ravel()[keep]
        z=p.ravel()[keep]
        pair,energy=rates(z,n,mass,na,True)
        J[i]=pair.sum()/(H*s)
        inc_E[i]=energy.sum()/(H*s)
    # Full final-state factors for this signed light↔scalar process.
    _,l,_=gross(edge,m,mp,chi,phi,False,nr,na,10,6,False,exchange_mass=M)
    _,q,_=gross(edge,mp,m,phi,chi,True,nr,na,10,6,False,exchange_mass=M)
    net_pair=float((l[1]-q[1])/(H*s))
    pair_heat=float((l[2]-q[2])/(H*s))
    dY=2*(J+np.array([-net_pair,net_pair]))
    heat=pair_heat+Es/M*inc_E[1]-(1-Es/M)*inc_E[0]
    return dict(phi_and_light_heavy_forward_per_H_s=J.tolist(),light_to_scalar_and_inverse_per_H_s=[l[1]/(H*s),q[1]/(H*s)],
                net_pair_only_light_and_scalar_per_H_s=[-2*net_pair,2*net_pair],
                projected_light_and_scalar_omitted_rate_per_H_s=dY.tolist(),
                projected_scalar_energy_gain_per_H_s=heat,opposite_light_energy_gain_per_H_s=-heat)


def from_shapes(edge,z):
    return [occ(Cells(edge,ma,g),p) for (ma,g),p in zip(((mp,1),(m,4)),z,strict=True)]


def main(eos,seed,cache):
    cache.mkdir(parents=True,exist_ok=True)
    b=StandardModel(eos)
    s=json.loads(seed.read_text())
    times=np.array((0.,1e-10,3e-10,7e-10,1.5e-9,3e-9,5e-9,8e-9,1.15e-8,1.5e-8))
    out=[]
    norm=b.state(s["TF"])[2]
    initial=None
    for x in times:
        p=cache/f"{x:.4e}.json"
        if p.exists():
            d=json.loads(p.read_text())
        else:
            d=run(s,b,4,max(x,1e-10),1.25e-10,cap=17.7,record_terminal=True)
            p.write_text(json.dumps(d)+"\n")
        z=d["terminal_state"]
        e=np.asarray(z["momentum_edges_GeV"])
        start=np.asarray(z["initial_species_number_and_kinetic_energy"])
        initial=start[:,0].sum(axis=1)[::-1]/norm
        if x==0:
            f=[Cells(e,ma,g).reconstruct(w)[0] for (ma,g),w in zip(((mp,1),(m,4)),start,strict=True)]
        else:
            f=z["species_entropy_variables"]
        out.append(dict(time_in_Hubble_units=x,controls=[sample(e,f,s["TF"],b,k) for k in ((2,10),(3,15))]))
    integ=[]
    for stride in (2,1):
        keep=np.unique(np.r_[np.arange(0,len(out),stride),len(out)-1])
        X=times[keep]
        rows=[out[i]["controls"][-1] for i in keep]
        names=("phi_and_light_heavy_forward_per_H_s","net_pair_only_light_and_scalar_per_H_s",
               "projected_light_and_scalar_omitted_rate_per_H_s","projected_scalar_energy_gain_per_H_s")
        values={k:trapezoid(np.array([v[k] for v in rows]),x=X,axis=0).tolist() for k in names}
        integ.append(dict(histories=len(keep),integrated=values,
                           projected_count_change_over_initial=(np.array(values[names[2]])/initial).tolist()))
    return dict(scope="Selected omitted rates folded as diagnostics on fixed-background fast-only trajectory. "
                      "Includes Pauli and Bose factors in signed light↔scalar conversion. Heavy production "
                      "is a gross Pauli-free projection with prompt energy sharing; heavy backward reactions "
                      "are not included here. The folded path does not feed back onto the evolved spectra.",
                rows=out,initial_light_scalar_yields=initial.tolist(),time_integrals=integ)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","seed","cache"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(main(a.eos,a.seed,a.cache),indent=2,allow_nan=False))
