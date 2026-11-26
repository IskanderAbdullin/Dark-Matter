"""Conditional signed vacuum on evolving prompt source; no feedback on that source."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from archived_scalar_pair_reaction_reference import amplitude
from finite_cooling_reheat_distribution_reference import deposit
from finite_cooling_reheat_portal_reference import SourceFold,TF
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import m,mp

def rate(p,w,rev,n,moment=False):
    I=m if not rev else mp
    O=mp if not rev else m
    i,j=np.tril_indices(len(p))
    x,y=roots_legendre(n)
    E,F=np.hypot(p[i],I),np.hypot(p[j],I)
    EF=E*F
    P=p[i]*p[j]
    cut=np.clip((EF+I*I-2*O*O)/P,-1.,1.) if rev else np.ones_like(P)
    sel=cut>-1.
    cut,EF,P,i,j=[a[sel] for a in (cut,EF,P,i,j)]
    mu=((cut+1)[:,None]*x+(cut-1)[:,None])/2
    W=(cut+1)[:,None]*y/4
    s=2*(I*I+EF[:,None]-P[:,None]*mu)
    weight=w[i]*w[j]*np.where(i==j,1.,2.)/(64*np.pi*EF)
    if moment:
        assert rev
        weight=weight*(p[i]**2+p[j]**2)/2
    scale=np.sqrt(np.maximum(0.,.25-O*O/s))
    C=0.
    for c,y0 in zip(x,y):
        C+=y0*(weight@((W*scale*amplitude(s,c)) .sum(axis=1)))
    return C


def occ(Y,edges,a,k,mode):
    n,p2=Y[k]
    keep=(n>n.sum()*1e-14)&(p2>0.)
    mu=p2[keep]/n[keep]
    if mode==0:
        p=np.sqrt(mu); w=n[keep]
    else:
        cut=keep[:-1]
        lo=edges[:-2][cut];hi=edges[1:-1][cut]
        n=n[:-1][cut]; q2=p2[:-1][cut]
        f=np.clip((q2/n-lo**2)/(hi**2-lo**2),0.,1.)
        p=np.r_[lo,hi];w=np.r_[n*(1-f),n*f]
    sel=(p>0.)&(w>0.)
    return (TF/a)*p[sel],w[sel]

def check():
    from archived_scalar_pair_reaction_reference import gross,points
    E=np.r_[0.,np.geomspace(.0001,.50,42)]
    f=lambda x: np.exp(-x/.018)
    def one(I,O,rev):
        p,_,w=points(E,I,f,3)
        u=rate(p,w,rev,10)*4
        _,r,_=gross(E,I,O,f,np.zeros_like,rev,3,10,16,12,split_final=False)
        return abs(u/r[0]-1.),u
    f,r=one(m,mp,False),one(mp,m,True)
    err=max(f[0],r[0])
    assert err<2e-10
    return dict(max_independent_lab_angle_relative_rate_defect=err,
        finite_support_classical_equilibrium_inverse_over_forward=r[1]/f[1])


def run(eos,reference,recoil_moments=False):
    T=json.loads(reference.read_text())
    hsh=hashlib.sha256(eos.read_bytes()).hexdigest()
    assert T["eos_sha256"]==hsh
    b=SourceFold(eos,(7,18,14),100,4)
    conv=conversion(eos)['Omega_h2_per_GeV_times_yield']*mp
    ref=T['controls'][-1]['rows'][-1]
    L,r,h=np.log([TF,ref['T_R_GeV'],ref['Tmax_GeV']])
    sf=b.s(L);ar=(sf/b.s(r))**(1/3)
    def v(u,pre):
        if pre:return ar*np.exp(-8*(u-r)/3),(8/3)/(b.H(r)*np.exp(4*(u-r)))
        return (sf/b.s(u))**(1/3),(1+b.grad(u)/3)/b.H(u)
    out=[]
    if recoil_moments:
        from finite_cooling_pair_emission_reference import kernel_and_lab_p2
    for nt,nb,na,mode in ((28,48,7,0),(54,62,9,1)):
        edge=np.r_[0.,np.geomspace(.00001,140,nb),np.inf]
        Y=np.zeros((2,2,len(edge)-1));Q=np.zeros((2,2));J=np.zeros((2,2))
        xx,ww=roots_legendre(3)
        maxratio=0.;probes=[]
        for epoch,(hi,lo) in enumerate(((h,r),(r,L))):
            pre=epoch==0
            G=np.unique(np.r_[np.linspace(lo,hi,nt+1),np.clip(np.log([.009,.015,.020,.025,.03]),lo,hi)])
            for U,V in zip(G[:0:-1],G[-2::-1]):
                c=(U+V)/2
                for i,(l,u) in enumerate(((c,U),(V,c))):
                    for x,w in zip((l+u+(u-l)*xx)/2,(u-l)*ww/2):
                        a,dt=v(x,pre)
                        Y+=deposit(b,np.exp(x),edge,a,20)*(a**3*dt*w/sf)
                    if i!=0:continue
                    a,dt=v(c,pre)
                    p,n=occ(Y,edge,a,1,mode);q,mw=occ(Y,edge,a,0,mode)
                    forw=rate(p,n,False,na)/4
                    rev=4*rate(q,mw,True,na)
                    clock=2*(U-V)*dt*sf/a**3
                    Q[epoch]+=np.array([forw,rev])*clock
                    if recoil_moments:
                        F2=kernel_and_lab_p2(p,n,na)[0]/4
                        assert abs(F2[0]/forw-1.)<2e-11
                        I2=4*rate(q,mw,True,na,moment=True)
                        J[epoch]+=np.array([F2[1],I2])*clock*(a/TF)**2
                    maxratio=max(maxratio,rev/forw)
                    probes.append([np.exp(c),rev/forw,forw,rev])
        result=dict(panels=nt,configuration=[nb,na,mode],omega_forward_inverse_by_epoch=(Q*conv).tolist(),
            omega_signed_by_epoch=((Q[:,0]-Q[:,1])*conv).tolist(),max_local_inverse_over_forward=float(maxratio),
            trace_columns=["T_GeV","inverse_over_forward","forward_weighted_yield_squared_kernel","inverse_weighted_yield_squared_kernel"],
            traces=probes)
        if recoil_moments:
            yp,mp2=Y[0].sum(axis=1)
            yf,yr=Q.sum(axis=0)
            jf,jr=J.sum(axis=0)
            assert yp>yr and mp2>jr and np.all(J>0.)
            result["anchor_moments"]=dict(
                 scalar_prompt_Y=float(yp),scalar_prompt_Y_q2=float(mp2),
                 scalar_omega_times_q2_forward_inverse_by_epoch=(J*conv).tolist(),
                 standalone_absorbed_rms_by_epoch=np.sqrt(J[:,1]/Q[:,1]).tolist(),
                 standalone_forward_rms_by_epoch=np.sqrt(J[:,0]/Q[:,0]).tolist(),
                 forward_only_mixed_rms=float(np.sqrt((mp2+jf)/(yp+yf))),
                 signed_passive_mixed_rms=float(np.sqrt((mp2+jf-jr)/(yp+yf-yr))),
                 removed_fraction_of_forward_only_count=float(yr/(yp+yf)),
                 removed_fraction_of_forward_only_q2_yield=float(jr/(mp2+jf)))
        out.append(result)
    return dict(eos_sha256=hsh,reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
        independent_checks=check(),cases=out,scope="Signed spontaneous vacuum rate folded on the forced prompt-only spectra. No evolution due to the rates, no secondary recoil feedback, no Bose stimulation and no Pauli factors in this fold. Near-sampled local quantum rates are checked separately. Not a signed evolving relic calculation.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos",type=Path)
    p.add_argument("reference",type=Path)
    p.add_argument("--recoil-moments",action="store_true")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.reference,a.recoil_moments),separators=(",",":"),allow_nan=False))
