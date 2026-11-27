"""Audit instantaneous cold scalar count currents at the warm archived anchors."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_visible_heavy_adiabatic_reference import angle_panels, snapshot_functions
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import Es, M, P, m, mp, width
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def angular_counts(S, cuts, n, closed):
    f,ep,g,el=snapshot_functions(S)
    A=np.unique(np.r_[ep,cuts])
    B=el
    def wall_parent(w,mass,rest):
        q=np.sqrt((w-mass)*(w+mass))
        d=w*rest+q*P
        return np.r_[M*(rest*rest+q*q)/d,M*d/(mass*mass)]
    peak=ep[-1]+el[-1]
    edges=np.unique(np.clip(np.r_[M,peak,wall_parent(A,mp,Es),wall_parent(B,m,M-Es),
                                  (A[:,None]+B).ravel()],M,peak))
    x,w=roots_legendre(n)
    dx=np.diff(edges)[:,None]/2
    en=(edges[:-1,None]+dx*(1+x)).ravel()
    de=(dx*w).ravel()
    data=np.zeros((len(cuts)+1,3))
    for i in range(0,len(en),96):
        E=en[i:i+96]
        P1,W=angle_panels(E,A,B,angle=n)
        P2=E[:,None,None]-P1
        if closed:
            W=W*((P1>=ep[0])&(P1<=ep[-1])&(P2>=el[0])&(P2<=el[-1]))
        f1,f2=f(P1),g(P2)
        inv=W*f1*f2
        dec=W*(1+f1)*(1-f2)
        s=inv.sum(axis=(1,2))
        q=s+dec.sum(axis=(1,2))
        h=np.divide(s,q,out=np.zeros_like(s),where=q>0.)[:,None,None]
        pos=h*dec
        neg=(1-h)*inv
        K=de[i:i+96]*2/np.pi**2*width*M*np.sqrt((E-M)*(E+M))
        for j,cut in enumerate(np.r_[cuts,np.inf]):
            v=P1<cut
            data[j]+=[K@(pos*v).sum(axis=(1,2)),K@(neg*v).sum(axis=(1,2)),
                      K@((pos-neg)*v).sum(axis=(1,2))]
    return data


def counts(S, pcut, n):
    f,_,_,_=snapshot_functions(S)
    q=np.asarray(S["lab_momentum_edges_GeV"])
    x,w=roots_legendre(n)
    out=[]
    for cut in pcut:
        e=np.unique(np.r_[q[q<cut],cut])
        h=np.diff(e)[:,None]/2
        v=e[:-1,None]+h*(1+x)
        out.append(float(np.sum(h*w*v*v*f(np.hypot(v,mp)))/(2*np.pi**2)))
    return np.asarray(out)


def recoil_gate(q, qmax):
    # Lowest parent energy able to emit a scalar at momentum q. All
    # test q fall below the rest-frame daughter momentum, so this
    # lower threshold decreases through each cumulative cold interval.
    Ec=np.hypot(mp,q)
    assert np.all((q>0.)&(q<P))
    E_min=M*(Es**2+q**2)/(Es*Ec+P*q)
    E_light_cap=np.hypot(m,qmax)
    E_parent_cap=np.hypot(mp,qmax)+E_light_cap
    E_l_min=E_min-Ec
    assert np.all(E_l_min>m)
    q_l_min=np.sqrt((E_l_min-m)*(E_l_min+m))
    return E_min, q_l_min, E_parent_cap, E_light_cap


def run(eos,paths):
    out=[]
    p_over_T=np.array([.01,.03,.1,.3,1.,3.])
    for path in paths:
        x=read(path)
        S=x["snapshot"]
        assert x["snapshot_sha256"]==digest(S)
        T=S["T_visible_GeV"]
        p=np.asarray(S["lab_momentum_edges_GeV"])
        pcut=p_over_T*T
        assert p[0]==0. and pcut[-1]<p[-1]
        N,_,E=get_mom(p,[S["entropy_variable_"+k] for k in ("scalar","fermion")])
        H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+E.sum())/3.)/MPL
        scalar=N[0].sum()
        residence=counts(S,pcut,16)
        np.testing.assert_allclose(counts(S,pcut,24),residence,rtol=4e-10)
        parent_min,light_min,parent_cap,light_cap=recoil_gate(pcut,p[-1])
        controls=[]
        for closed in (True,False):
            A=[]
            for n in (6,10,14):
                a=angular_counts(S,np.hypot(mp,pcut),n,closed)/(H*scalar)
                impossible=parent_min>(np.hypot(mp,pcut)+light_cap if closed else parent_cap)
                assert np.all(a[:-1][impossible]==0.)
                np.testing.assert_allclose((a[:,0]-a[:,1]),a[:,2],atol=1e-12*max(a[-1,:2]))
                assert abs(a[-1,2])<3e-14*max(a[-1,:2])
                A.append(a)
            amounts=np.max(A[-1][:,:2],axis=1)
            zero=amounts==0.
            assert np.all(np.asarray(A)[:,zero,:]==0.)
            assert np.max(abs(A[-1][~zero]-A[-2][~zero])/amounts[~zero,None])<1e-4
            controls.append(dict(closed=closed,orders=[6,10,14],gross_gain_loss_net=A,
                              net_over_resident_per_H=(A[-1][:-1,2]*scalar/residence).tolist()))
        out.append(dict(snapshot=x["snapshot_sha256"],visible_T_GeV=T,pcut_over_T=p_over_T.tolist(),
                        scalar_cumulative_counts_over_represented_N=(residence/scalar).tolist(),
                        N_scalar=float(scalar),represented_rho_GeV4=float(E.sum()),
                        H_GeV=float(H),archive_common_p_cap_GeV=float(p[-1]),
                        open_parent_energy_cap_GeV=float(parent_cap),
                        lowest_required_parent_energy_GeV=parent_min.tolist(),
                        lowest_required_light_partner_momentum_GeV=light_min.tolist(),
                        rate_controls=[
                            {**a,"gross_gain_loss_net":[z.tolist() for z in a["gross_gain_loss_net"]]} for a in controls]))
    return dict(scope="Frozen archived daughters with equilibrium eliminated stationary prompt parent. "
                "Incoming exterior set to vacuum. Outputs are currents per original represented scalar number "
                "and per H time, not an evolution; open geometry counts outgoing scalars even outside the archive.",
                cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("archives",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archives),indent=2,allow_nan=False))
