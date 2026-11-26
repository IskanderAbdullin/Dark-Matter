"""Frozen tagged f=1/2 prompt and pair tests, not a dynamical cap proof."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_visible_heavy_adiabatic_reference import parent_branch,snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import m,width
from warm_cold_light_finite_pair_reference import tagged
from warm_light_occupation_stage_audit_reference import digest,read
from warm_near_rest_light_inverse_decay_reference import rate,window
from warm_prompt_reentry_shape_audit_reference import get_mom


def frozen(eos,path):
    v=read(path)
    v=v[-1] if isinstance(v,list) else v
    S=v["snapshot"]
    assert digest(S)==v["snapshot_sha256"]
    p=np.asarray(S["lab_momentum_edges_GeV"])
    B,be,F,fe=snapshot_functions(S)
    rad=TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]
    th=[S["entropy_variable_"+z] for z in ("scalar","fermion")]
    H=np.sqrt((rad+get_mom(p,th)[2].sum()+import_restart(S["positive_exterior_restart"]).energy(1.))/3)/MPL
    sampled=np.r_[0.,(p[:-1]+p[1:])/2,np.nextafter(p[-1],0.)]
    orders=((8,6,8,8),(14,10,14,12))
    output=[]
    for (r,i,o,z),n in zip(orders,(10,16),strict=True):
        def parent(E,closed):
            original_shape=E.shape
            E=np.asarray(E).ravel()
            if not closed:
                v=parent_branch(E,B,F,be,fe,angle=n)[:,0]
            else:
                # Preserve original closed geometry conditioning for both
                # directions, not only absorption. At zero overlap h=0.
                from joint_visible_heavy_adiabatic_reference import angle_panels
                v=np.zeros_like(E)
                for j in range(0,len(E),64):
                    A=E[j:j+64]
                    x,w=angle_panels(A,be,fe,angle=n)
                    y=A[:,None,None]-x
                    w=w*((x>=be[0])&(x<be[-1])&(y>=fe[0])&(y<fe[-1]))
                    b,f=B(x),F(y)
                    gain=np.sum(w*b*f,axis=(1,2))
                    loss=np.sum(w*(1+b)*(1-f),axis=(1,2))
                    v[j:j+64]=np.divide(gain,gain+loss,out=np.zeros_like(gain),where=gain+loss>0)
            return v.reshape(original_shape)
        rows=[]
        for pt in sampled:
            Et=np.hypot(m,pt)
            cap=lambda E: np.where(E==Et,.5,F(E))
            # Change the tagged external mode only. An isolated energy
            # coincidence in integration partner nodes is excluded.
            edge=np.unique(np.clip(np.r_[m,fe,2*be-Et,max(fe[-1],2*be[-1]-Et)],m,max(fe[-1],2*be[-1]-Et)))
            from scipy.special import roots_legendre
            x,_=roots_legendre(r)
            E=(edge[:-1,None]+edge[1:,None])/2+np.diff(edge)[:,None]*x/2
            use=F if np.any(E==Et) else cap
            pair=tagged(pt,use,B,fe,be,(r,i,o,z))[0]
            if use is F:
                a=F(np.array([Et]))[0]
                assert 0<a<1
                pair=pair*np.array([.5/(1-a),.5/a])
            P=[]
            for closed in (False,True):
                def heavy(E):
                    v=parent(E,closed)
                    if closed:
                        test=E-Et
                        v=np.where((test>=be[0])&(test<be[-1]),v,0.)
                    return v
                K,*_=rate([pt],lambda e:np.full_like(e,.5),B,heavy,be,order=n)
                P.append([float((K[0,0]+K[1,0])*width/H),float(K[2,0]*width/H)])
            rows.append(dict(p_GeV=float(pt),archived_f=float(F(np.array([Et]))[0]),
                scalar_q_support_gap_GeV=float(window(pt)[-2]-be[-1]),
                pair_gain_loss_over_H=(pair/H).tolist(),open_prompt_gain_loss_over_H=P[0],
                represented_prompt_gain_loss_over_H=P[1]))
        output.append(dict(q_order=[r,i,o,z],prompt_order=n,modes=rows))
    return dict(scope="Frozen bath at tagged light f=1/2; not a continuum or time-dependent cap and does not model elastic.",
                snapshot_sha256=digest(S),file=path.name,T_GeV=S["T_visible_GeV"],controls=output)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("paths",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps([frozen(a.eos,p) for p in a.paths],indent=2,allow_nan=False))
