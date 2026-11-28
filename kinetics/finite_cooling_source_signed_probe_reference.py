"""Signed rates on passive prompt spectra at labeled cooling epochs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import roots_legendre
from archived_scalar_pair_reaction_reference import gross
from export_finite_cooling_class_psd_reference import smooth
from finite_cooling_reheat_distribution_reference import deposit
from finite_cooling_reheat_portal_reference import TF,SourceFold
from three_species_decay_onsager_reference import m,mp

def run(eos,reference):
    ref=json.loads(reference.read_text())
    key=hashlib.sha256(eos.read_bytes()).hexdigest()
    assert key==ref["eos_sha256"]
    r0=ref["controls"][-1]["rows"][-1]
    F=SourceFold(eos,(7,18,14),100,4)
    u0,r,h=np.log([TF,r0['T_R_GeV'],r0['Tmax_GeV']])
    sf=F.s(u0); ar=(sf/F.s(r))**(1/3)
    E=np.r_[0.,np.geomspace(.00001,140,53),np.inf]
    X,W=roots_legendre(4)
    Y=np.zeros((2,2,len(E)-1))
    requested=np.sort(np.r_[u0,np.log([.0035,.005,.009,.015,.025]),r])
    grid=np.unique(np.r_[requested,np.clip(F.bg.sg.x,u0,h),np.linspace(u0,r,28),np.linspace(r,h,34)])
    assert grid[0]==u0 and grid[-1]==h

    def coefficients(u):
        if u>r:
            a=ar*np.exp(-(u-r)*8/3)
            H=F.H(r)*np.exp(4*(u-r))
            dt=(8/3)/H
        else:
            a=(sf/F.s(u))**(1/3)
            H=F.H(u)
            dt=(1+F.grad(u)/3)/H
        return a,H,dt

    def make(n,p2,a,mass,g,sig):
        ok=(n>1e-14*n.sum())&(p2>0.)
        assert n[~ok].sum()/n.sum()<1e-11
        pp=TF/a*np.sqrt(p2[ok]/n[ok])
        q,f=smooth(pp,n[ok]*sf/a**3,sig,4200)
        target=n.sum()*sf/a**3
        f*=target/(4*np.pi*np.trapezoid(f*q*q,q))
        Fint=CubicSpline(q,(2*np.pi)**3*f/g,extrapolate=False)
        def fn(e):
            p=np.sqrt(np.maximum(0.,(e-mass)*(e+mass)))
            return np.maximum(0.,np.nan_to_num(Fint(p),nan=0.))
        return q,fn,float(np.max(Fint(q)))

    out=[]
    for H,L in zip(grid[:0:-1],grid[-2::-1]):
        for u,w in zip((H+L+(H-L)*X)/2,(H-L)*W/2):
            a,_,dt=coefficients(u)
            Y+=deposit(F,np.exp(u),E,a,20)*(a**3*dt*w/sf)
        if L not in requested:continue
        a,H,_=coefficients(L)
        for sig,ang in ((.10,8),(.14,12)):
            qs,fs,ps=make(*Y[0],a,mp,1,sig)
            qc,fc,pc=make(*Y[1],a,m,4,sig)
            edges=np.r_[0.,np.geomspace(min(qs[0],qc[0])*.9,max(qs[-1],qc[-1])*1.01,54)]
            na,rF,_=gross(edges,m,mp,fc,fs,False,3,ang,ang,ang,split_final=False)
            ns,rR,_=gross(edges,mp,m,fs,fc,True,3,ang,ang,ang,split_final=False)
            assert abs(na/(Y[1,0].sum()*sf/a**3)-1.)<.035
            assert abs(ns/(Y[0,0].sum()*sf/a**3)-1.)<.035
            row=dict(T_GeV=float(np.exp(L)),sigma=sig,n_scalar_over_stable=ns/na,
                     scalar_max_reconstructed_f=ps,stable_max_reconstructed_f=pc,
                     inverse_over_forward=rR[1]/rF[1],stimulated_forward_over_vacuum=rF[1]/rF[0],
                     fractional_forward_stable_loss_per_H=2*rF[1]/(na*H))
            out.append(row)
    return dict(eos_sha256=key,reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
       rows=out,scope="Signed local fixed-shape scalar-pair rate at labeled epochs of the passive prompt source. The reconstructed scalar distributions omit reaction-produced scalars; this is not an integrated quantum solution or a rigorous error bound.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos",type=Path)
    p.add_argument("reference",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.reference),indent=2,allow_nan=False))
