"""Conservative pair collision with both quantum directions on one fermion chart."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import amplitude, gross, m, mp, points
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from scalar_pair_expanding_continuation_reference import e_total
from warm_cold_light_finite_pair_reference import tagged
from warm_light_occupation_stage_audit_reference import digest, read


def rate(p,F,B,order,ext=4,verify_balance=None):
    """Use vacuum occupations above the represented bath, but allow empty final states."""
    n = len(p) - 1
    max_E = max(2*np.hypot(mp,p[-1])-m,m)
    high = np.sqrt((max_E-m)*(max_E+m))
    grid = np.r_[p,np.geomspace(p[-1],high,ext+1)[1:]] if high>p[-1] else p
    nr,na,nz,nf = order
    pp,ee,rad = points(grid,m,np.ones_like,nr)
    ff = F(ee)
    i,j = np.tril_indices(len(pp))
    k = ((ff[i]*ff[j]>0.) | ((ee[i]+ee[j])<2*np.hypot(mp,p[-1])))
    i,j=i[k],j[k]
    u,w=roots_legendre(na)
    t,tw=roots_legendre(nf)
    az=np.cos(np.pi*(np.arange(nz)+.5)/nz)
    eF=np.hypot(m,p); eB=np.hypot(mp,p)
    output=np.zeros((2,2,2,n+2))
    gross_rate=np.zeros((2,2))

    def scatter(array, edge, EE, W, mass):
        a=np.broadcast_to(EE,W.shape).ravel()
        q=W.ravel()
        s=np.searchsorted(edge,a,side='right')
        array[0] += np.bincount(s,q,minlength=n+2)
        array[1] += np.bincount(s,q*np.maximum(0.,a-mass),minlength=n+2)

    for l in range(0,len(i),48):
        I,J=i[l:l+48],j[l:l+48]
        E,F0=ee[I,None],ee[J,None]
        U,V=pp[I,None],pp[J,None]
        s=2*m*m+2*(E*F0-U*V*u)
        q=np.sqrt(np.maximum(0.,s/4-m*m))
        k=np.sqrt(s/4-mp*mp)
        boost=np.sqrt((U-V)**2+2*U*V*(1+u))/np.sqrt(s)
        d=boost*k
        mid=np.broadcast_to((E+F0)/2,s.shape)
        cut=np.clip((eB-mid[...,None])/np.maximum(d[...,None],1e-280),-1.,1.)
        cut=np.sort(np.concatenate([-np.ones_like(cut[...,:1]),cut,-cut,
                                    np.ones_like(cut[...,:1])],axis=-1),axis=-1)
        h=np.diff(cut,axis=-1)/2
        a,b,c=np.where(h>0.)
        z=(cut[a,b,c,None]+cut[a,b,c+1,None])/2+h[a,b,c,None]*t
        c0=np.clip((E-F0)/(2*np.maximum(boost*q,1e-280)),-1.,1.)[a,b,None]
        A0= c0*z
        A1= np.sqrt(np.maximum(1-c0*c0,0.))*np.sqrt(1-z*z)
        cos=np.clip(A0[...,None]+A1[...,None]*az,-1.,1.)
        M=amplitude(s[a,b,None,None],cos).mean(axis=-1)
        weight=4*(rad[I]*rad[J]*np.where(I==J,1.,2.))[:,None]*(w/2)*k/(64*np.pi*E*F0*np.sqrt(s))
        W=weight[a,b,None]*(h[a,b,c,None]*tw)*M
        e0=mid[a,b,None]+d[a,b,None]*z
        e1=mid[a,b,None]-d[a,b,None]*z
        b0,b1=B(e0),B(e1)
        f0,f1=ff[I[a],None],ff[J[a],None]
        G=W*(1.-f0)*(1.-f1)*b0*b1
        L=W*f0*f1*(1.+b0)*(1.+b1)
        if verify_balance is not None:
            interior=(E[a]<eF[-1])&(F0[a]<eF[-1])&(e0<eB[-1])&(e1<eB[-1])
            valid=interior&(G+L>1e-240)
            if np.any(valid):
                diff=np.max(abs(G[valid]-L[valid])/np.maximum(G[valid],L[valid]))
                verify_balance[0]=max(verify_balance[0],float(diff))
                verify_balance[1]+=int(np.sum(valid))
        for n0,Q in enumerate((G,L)):
            output0=output[n0]
            energy=E[a]+F0[a]
            gross_rate[n0]+=[Q.sum(),np.sum(energy*Q)]
            scatter(output0[0],eB,e0,Q,mp)
            scatter(output0[0],eB,e1,Q,mp)
            scatter(output0[1],eF,E[a],Q,m)
            scatter(output0[1],eF,F0[a],Q,m)

    for R in range(2):
        side=output[R]
        ledger=gross_rate[R]
        if ledger[0]>0.:
            count=side[:,0].sum(axis=1)
            energy=side[:,1].sum(axis=1)+np.array([mp,m])*count
            assert max(abs(count/(2*ledger[0])-1.))<1e-12
            assert max(abs(energy/ledger[1]-1.))<1e-12
    delta=(output[0]-output[1])*np.array([-1,1])[:,None,None]
    assert np.all(delta[:,:,[0,-1]]>=0.)
    d=max(abs(delta[:,0].sum())/max(2*gross_rate[:,0].sum(),1e-280),
          abs(e_total(delta))/max(gross_rate[:,1].sum(),1e-280))
    assert d<2e-12
    return delta,output,gross_rate,d,high


def compare(path):
    row=read(path)
    row=row[-1] if isinstance(row,list) else row
    S=row["snapshot"]; assert digest(S)==row["snapshot_sha256"]
    p=np.asarray(S["lab_momentum_edges_GeV"])
    B,be,F,fe=snapshot_functions(S)
    out=[]
    for order in ((3,7,6,2),(4,10,10,4),(6,14,14,6)):
        Z,ret,b,d,lim=rate(p,F,B,order)
        _,Fg,L=gross(p,m,mp,F,B,False,*order,project=True)
        _,Rg,G=gross(p,mp,m,B,F,True,*order,project=True)
        old=np.array([L[1]-G[0],G[1]-L[0]])
        loss=np.array([L[1],L[0]])
        loss_rel=float(abs(ret[1]-loss).sum()/max(abs(loss).sum(),1e-280))
        assert loss_rel<2e-11
        print(f"completed {order}", file=__import__("sys").stderr, flush=True)
        out.append(dict(quadrature=order,direction_gross_number_energy=b.tolist(),
                        inverse_over_original_number_energy=(b[0]/Rg[1:3]).tolist(),
                        forward_match=loss_rel, ledger_closure=d, max_exterior_p_GeV=lim,
                        first_cell_gain_over_loss_original=float(G[1,0,1]/L[0,0,1]),
                        first_cell_gain_over_loss_common=float(ret[0,1,0,1]/ret[1,1,0,1]),
                        normalized_signed_change=(abs(Z-old).sum(axis=2)/abs(old).sum(axis=2)).tolist()))
    x,w=roots_legendre(6); q=p[1]*(1+x)/2
    reference=np.array([tagged(pt,F,B,fe,be,(12,9,12,12))[0] for pt in q]).T@(w*q*q)
    out[-1]["independent_first_cell_gain_over_loss"]=float(reference[0]/reference[1])
    return dict(file=path.name,T_GeV=S["T_visible_GeV"],snapshot_sha256=digest(S),controls=out)


def thermal():
    T, mu = .017, -.03
    p=np.r_[0.,np.geomspace(.00004,.002,7),np.geomspace(.007,1.07,14)]
    energies=[np.hypot(p,m0) for m0 in (mp,m)]
    z=[np.vstack(((u[:-1]-mu)/T,np.diff(u)/T)) for u in energies]
    spec=dict(lab_momentum_edges_GeV=p,entropy_variable_scalar=z[0],entropy_variable_fermion=z[1])
    B,_,F,_=snapshot_functions(spec)
    test=[0.,0]
    net,projected,G,d,lim=rate(p,F,B,(3,7,6,2),verify_balance=test)
    defect=net[1,0,1:8]/projected[1,1,0,1:8]
    assert test[0]<8e-14 and test[1]>=100
    assert np.max(abs(defect))<2e-11
    return dict(max_relative_event_balance_defect=test[0],nonzero_events_checked=test[1],
                first_7_shell_net_over_loss=defect.tolist(),global_conservation=d,
                truncated_domain_GeV=p[-1],max_virtual_GeV=lim)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("snapshots",nargs="*",type=Path)
    p.add_argument("--thermal",action="store_true")
    a=p.parse_args()
    print(json.dumps(dict(thermal=thermal() if a.thermal else None,archives=[compare(path) for path in a.snapshots]),indent=2,allow_nan=False))
