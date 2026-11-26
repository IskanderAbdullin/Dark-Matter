"""Exact outward on-shell flux immediately after a finite-support snapshot.

Uses the archived FD/BE shell shapes and explicit vacuum outside their
support. For b>1 this is an exact collisionless redshift of the snapshot,
used as a fixed-state probe, not the actual interacting trajectory.
"""
import argparse,json
from pathlib import Path
import numpy as np

from continuum_dark_normal_capacity_reference import from_row,rule
from joint_visible_heavy_adiabatic_reference import angle_panels,snapshot_functions
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import M,m,mp,P,Es,width


def outward(snap,b,order=(40,8),boundary_GeV=None):
    pp=np.array(snap["lab_momentum_edges_GeV"])/b
    phi,_,chi,_=snapshot_functions(snap)
    def transport(f,E,mass):
        k=(E-mass)*(E+mass)
        return f(np.sqrt(np.maximum(0.,mass*mass+b*b*k)))
    ph=lambda e:transport(phi,e,mp)
    ch=lambda e:transport(chi,e,m)
    ep0=np.hypot(mp,pp)
    ec0=np.hypot(m,pp)
    pb=pp[-1] if boundary_GeV is None else boundary_GeV/b
    bp=np.hypot(mp,pb); bc=np.hypot(m,pb)
    hi=ep0[-1]+ec0[-1]
    eg=np.unique(np.r_[ep0,bp]); cg=np.unique(np.r_[ec0,bc])
    if hi<=M:
        return np.zeros(10)
    # Parent lab energies at which either daughter support cuts its rest-frame
    # angular envelope; the whole integral still includes both full tails.
    def gates(E,p,rest,mass):
        upper=rest*E+P*p
        lo=M*(rest*rest+p*p)/upper
        high=M*upper/(mass*mass)
        return np.r_[lo,high]
    edge=np.unique(np.clip(np.r_[M,hi,gates(np.r_[ep0,bp],np.r_[pp,pb],Es,mp),gates(np.r_[ec0,bc],np.r_[pp,pb],M-Es,m)],M,hi))
    u,w=rule(order[0])
    e=(edge[:-1,None]*(1-u)+edge[1:,None]*(1+u)).ravel()/2
    wt=(np.diff(edge)[:,None]*w/2).ravel()
    p=np.sqrt((e-M)*(e+M))
    pred=4*width*M*p*wt/(2*np.pi**2)
    out=np.zeros(10)
    for i in range(0,len(e),160):
        s=slice(i,i+160)
        E,W=angle_panels(e[s],eg,cg,angle=order[1])
        C=e[s,None,None]-E
        B=W*(1+ph(E))*(1-ch(C))
        A=W*ph(E)*ch(C)
        rateA=A.sum(axis=(1,2))
        rateB=B.sum(axis=(1,2))
        F=rateA/(rateA+rateB)
        wrong=(E>ep0[-1])|(C>ec0[-1])
        assert np.max(np.where(wrong,A,0.))==0.
        out[0]+=pred[s]@rateA
        out[1]+=pred[s]@(F*B.sum(axis=(1,2)))
        out[2]+=pred[s]@(F*(B*(E>bp)).sum(axis=(1,2)))
        out[3]+=pred[s]@(F*(B*(C>bc)).sum(axis=(1,2)))
        out[4]+=pred[s]@(F*(B*(E>bp)*(E-mp)).sum(axis=(1,2)))
        out[5]+=pred[s]@(F*(B*(C>bc)*(C-m)).sum(axis=(1,2)))
        out[6]+=pred[s]@((1.-F)*(A*(E>bp)).sum(axis=(1,2)))
        out[7]+=pred[s]@((1.-F)*(A*(C>bc)).sum(axis=(1,2)))
        out[8]+=pred[s]@((1.-F)*(A*(E>bp)*(E-mp)).sum(axis=(1,2)))
        out[9]+=pred[s]@((1.-F)*(A*(C>bc)*(C-m)).sum(axis=(1,2)))
    if boundary_GeV is None:
        assert out[2]+out[3]<=out[1]*(1+1e-11) and np.max(out[6:])==0.
    return out


def audit(eos,radiation,original,snapshot,boundary_GeV=None):
    cap=from_row(original,eos,20)
    ns=np.array(cap["species_counts_GeV3"])
    k=cap["retained_dark_kinetic_energy_GeV4"]
    bg=TabulatedRadiation(radiation,cap["actual_visible_T_GeV"])
    ans=[]
    for b in (1.,1.3,1.8,2.5,3.2,3.55,3.7,5.,7.,8.):
        a=outward(snapshot,b,(24,6),boundary_GeV)
        c=outward(snapshot,b,(40,10),boundary_GeV)
        H=np.sqrt(bg.state(b)[1]/3)/MPL
        # Absolute rates per original reference phase density. For energy use
        # the initial kinetic density scaled by volume, not current kinetic.
        norm=np.r_[ns[0],ns[0],ns,k,k,ns,k,k]*H/b**3
        rel=abs(c-a)/np.maximum(c,1e-300)
        assert np.all(a>=0.) and np.all(c>=0.) and np.all(np.isfinite(c))
        ans.append(dict(b=b,rates_normalized=(c/norm).tolist(),
                        signed_outside_net_normalized=((c[2:6]-c[6:])/norm[2:6]).tolist(),
                        max_relative_order_difference=float(np.max(rel))))
    return ans


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("radiation");p.add_argument("run");p.add_argument("shape")
    p.add_argument("--boundary-q-over-TR",type=float)
    a=p.parse_args()
    r=json.loads(Path(a.run).read_text())["rows"][-1]
    s=json.loads(Path(a.shape).read_text())["snapshots"][-1]
    labels=["gross_inverse", "steady_parent_forward", "exterior_scalar_count",
            "exterior_fermion_count", "exterior_scalar_kinetic", "exterior_fermion_kinetic","reverse_scalar_count","reverse_fermion_count","reverse_scalar_kinetic","reverse_fermion_kinetic"]
    pb=None if a.boundary_q_over_TR is None else a.boundary_q_over_TR*r["initial_visible_T_GeV"]/r["final_a"]
    print(json.dumps(dict(rows=audit(StandardModel(a.eos),a.radiation,r,s,pb),rate_columns=labels,
                          boundary_q_over_TR=a.boundary_q_over_TR),indent=2,allow_nan=False))
