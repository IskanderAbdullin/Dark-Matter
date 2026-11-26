"""Kinematic on-shell reachability on two finite momentum supports."""
import argparse,json
from pathlib import Path
import numpy as np

from continuum_dark_normal_capacity_reference import rule
from joint_prompt_two_moment_evolution_reference import Collision
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import M,m,mp,P


def minimum_incident(mass,partner,pmax):
    """Smallest incident p that can reach parent M with a bounded partner."""
    E=np.hypot(partner,pmax)
    C=(M*M-mass*mass-partner*partner)/2
    root=np.sqrt((C-mass*partner)*(C+mass*partner))
    return max(0.,(C*C-(E*mass)**2)/(C*pmax+E*root))


def pair_mask(edges):
    e=np.asarray(edges)
    rp=np.arcsinh(e/mp)
    rc=np.arcsinh(e/m)
    gap=np.maximum(0.,np.maximum(rp[:-1,None]-rc[1:],rc[None,:-1]-rp[1:,None]))
    hi=rp[1:,None]+rc[None,1:]
    C=(M*M-m*m-mp*mp)/(2*m*mp)
    return (np.cosh(hi)>C)&(np.cosh(gap)<C)


def shell_fraction(snap,b,physical_cut,j,order=28):
    p=np.array(snap["lab_momentum_edges_GeV"])
    phi,_,chi,_=snapshot_functions(snap)
    m0,f=(mp,phi) if j==0 else (m,chi)
    x,w=rule(order)
    def integral(a,e):
        e=np.maximum(e,a)
        q=(a[:,None]*(1-x)+e[:,None]*(1+x))/2
        w0=q*q*w*(e-a)[:,None]/2
        E=np.hypot(q,m0)
        E_live=np.hypot(q/b,m0)
        KE=(q/b)**2/(E_live+m0)
        f0=f(E)
        return np.array([np.sum(w0*f0),np.sum(w0*f0*KE)])
    return integral(p[:-1],np.minimum(p[1:],physical_cut*b))/integral(p[:-1],p[1:])


def audit(snap,original_boundaries,extras=()):
    q=np.array(snap["lab_momentum_edges_GeV"])
    out=[]
    for b in (1.,1.5,2.5,3.6,4.,5.,7.,8.):
        e=q/b
        mask=pair_mask(e)
        coll=Collision(e)
        sampled=np.zeros_like(mask)
        sampled[coll.I,coll.J]=True
        assert np.all(mask[sampled]),"quadrature admitted forbidden support"
        fp=minimum_incident(mp,m,e[-1])
        fc=minimum_incident(m,mp,e[-1])
        if e[-1]<P: assert not mask.any() and len(coll.y)==0
        out.append(dict(b=b,represented_maximum_lab_p_GeV=float(e[-1]),
                        scalar_p_floor_GeV=float(fp),fermion_p_floor_GeV=float(fc),
                        scalar_q_over_TR_floor=float(fp*b*original_boundaries),
                        scalar_count_kinetic_fraction_below_floor=shell_fraction(snap,b,fp,0).tolist(),
                        fermion_count_kinetic_fraction_below_floor=shell_fraction(snap,b,fc,1).tolist(),
                        kinematically_open_pairs=int(mask.sum()),populated_quadrature_pairs=int(sampled.sum()),
                        at_least_one_pair_for_first_scalar_shell=bool(mask[0].any())))
    trace=[]
    p0=q[1]
    for U in extras:
        if U==0.: continue
        phys=U/original_boundaries
        threshold=(M*M-m*m-mp*mp)/2
        Ephi=np.hypot(mp,p0)
        root=np.sqrt((threshold-mp*m)*(threshold+mp*m))
        # Required partner momentum at the first scalar bin's upper edge.
        required=(threshold*threshold-(Ephi*m)**2)/(threshold*p0+Ephi*root)
        tr=minimum_incident(mp,m,phys)
        trace.append(dict(upper_q_over_TR=U,floor_q_over_TR=float(tr*original_boundaries),
                          required_q_partner_for_first_scalar_bin_upper_edge=float(required*original_boundaries)))
    return dict(rows=out,tail_requirement_trace_at_snapshot=trace)


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("run");p.add_argument("shape")
    a=p.parse_args()
    r=json.loads(Path(a.run).read_text())["rows"][-1]
    s=json.loads(Path(a.shape).read_text())["snapshots"][-1]
    print(json.dumps(audit(s,r["final_a"]/r["initial_visible_T_GeV"],(16,32,64,128)),indent=2,allow_nan=False))
