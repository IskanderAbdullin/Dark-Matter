"""Controlled dilute thermal closure for selected dark reactions after fast return."""
import argparse
import json
from pathlib import Path
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq,root
from scipy.special import kve,roots_legendre

from archived_heavy_recycling_gate_reference import angular_integral
from archived_scalar_pair_reaction_reference import amplitude,gross
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import m,mp,M,y
from warm_same_charge_heavy_recycling_reference import exact_angular
from warm_energy_weighted_recycling_kernels_reference import geom_sq


@lru_cache(None)
def nodes(n):
    return roots_legendre(n)


def pair_mb(T,n):
    """Forward thermal event coefficients at unit fugacity; all GeV^4."""
    z,w=nodes(n)
    x,wx=nodes(n+12)
    out=[]
    for process,mi,mf in (("LH",m,M),("SH",mp,M),("LS",m,mp)):
        Emin=2*max(mi,mf)
        # Resolve the production threshold and retain an exponentially
        # suppressed tail, with fixed dimensionless panel boundaries.
        pieces=np.array((0.,.4,1.5,4.,9.,20.,39.,72.))
        mid=(pieces[:-1,None]+pieces[1:,None])/2
        h=np.diff(pieces)[:,None]/2
        q=Emin+T*(mid+h*z)
        s=q*q
        pa=np.sqrt(s/4-mi*mi)
        pb=np.sqrt(s/4-mf*mf)
        if process=="LH":
            mat=y**4*(angular_integral(s)+exact_angular(s))
            weighted=4*mat/32
        elif process=="SH":
            A=m*m+s/2-M*M-mp*mp
            B=2*pa*pb
            lo=geom_sq(s)/(A+B)
            gap=np.log1p(B/lo)
            u=-gap[:,:,None]*(1-x)/2
            zz=A[:,:,None]*np.exp(u)
            ca=np.clip(-A[:,:,None]*np.expm1(u)/B[:,:,None],0.,1.)
            mat=(zz*amplitude(s[:,:,None],ca,light=M,heavy=m,scalar=mp))@wx*gap/B
            weighted=.5*4*mat/32
        else:
            mat=amplitude(s[:,:,None],x)@wx
            weighted=4*mat/64
        sigma_weighted=weighted*pb/(np.pi*s*pa)
        k=kve(1,q/T)*np.exp(-q/T)
        coeff=T/(32*np.pi**4)*np.sum((2*q*T*h*w)*(s-4*mi*mi)*q*k*sigma_weighted)
        out.append(float(coeff))
    return np.array(out)


def thermo(T):
    masses=np.array((m,mp,M))
    g=np.array((4,1,4))
    x=masses/T
    n=g*masses**2*T*kve(2,x)*np.exp(-x)/(2*np.pi**2)
    e=3*T + masses*kve(1,x)/kve(2,x)
    return n,e


def state(T,L,S):
    n,e=thermo(T)
    c=n[2]/(n[0]*n[1])
    k=1+c*(L+S)
    H=(2*c*L*S)/(k+np.sqrt(k*k-4*c*c*L*S))
    N=np.array((L-H,S-H,H))
    z=N/n
    assert np.all(N>0.) and abs(z[2]/(z[0]*z[1])-1)<4e-12
    return N,z,e


def closed(latent,energy,Tguess):
    L,S=latent
    def test(t):
        n,z,e=state(t,L,S)
        return n@e/energy-1
    t=brentq(test,.15*Tguess,4*Tguess,xtol=2e-14)
    n,z,e=state(t,L,S)
    return t,n,z


def signed(T,z,n):
    c=pair_mb(T,n)
    l,s,h=z
    forward=c*np.array((l*l,s*s,l*l))
    reverse=c*np.array((h*h,h*h,s*s))
    net=forward-reverse
    # Units of physical density per physical time for the two counts
    # invariant under each individual fast decay/inverse decay.
    dl=2*(net[1]-net[2])
    ds=2*(net[0]+net[2])
    assert abs((dl+ds)-2*(net[0]+net[1]))<1e-12*max(forward.max(),reverse.max())
    return np.array((dl,ds)),forward,reverse


def quantum_state(latent,energy,T0,z0,k):
    """Dilute finite-fugacity series on the same unbounded momentum domain."""
    xi=np.arange(1,k+1)
    masses=np.array([m,mp,M])
    sign=np.stack(((-1.)**(xi-1),np.ones(len(xi)),(-1.)**(xi-1)),axis=1)
    g=np.array([4,1,4])
    def state(p):
        T=np.exp(p[0])
        eta=np.array([p[1],p[2],p[1]+p[2]])
        x=xi[:,None]*masses/T
        dens=sign*g*masses**2*(T/xi)[:,None]*kve(2,x)*np.exp(xi[:,None]*eta-x)/(2*np.pi**2)
        e=3*(T/xi)[:,None]+masses*kve(1,x)/kve(2,x)
        return T,eta,dens.sum(axis=0),(dens*e).sum()
    def residual(p):
        _,_,N,E=state(p)
        return np.log(np.r_[(N[:2]+N[2])/latent,E/energy])
    sol=root(residual,np.r_[np.log(T0),np.log(z0[:2])],tol=2e-11)
    assert max(abs(residual(sol.x)))<2e-11
    return state(sol.x)


def enclose_quantum_snapshot(latent,energy,T0,z0,units):
    v=quantum_state(latent,energy,T0,z0,10)
    T,eta,N,E=quantum_state(latent,energy,T0,z0,16)
    np.testing.assert_allclose(np.r_[T,eta,N,E],np.r_[v[0],v[1],v[2],v[3]],rtol=1e-11,atol=0)
    z=np.exp(eta)
    l,s,h=z*np.exp(-np.array([m,mp,M])/T)
    assert max(l,s,h)<.1
    F=pair_mb(T,52)*np.array([z[0]**2,z[1]**2,z[0]**2])/units
    lo=np.array([1/((1+l)*(1+h))**2,1/(1+h)**2,1/(1+l)**2])
    hi=np.array([1.,1/(1-s)**2,1/(1-s)**2])
    # For the genuine Bose/Fermi collision bracket on a common T,
    # reverse/forward is exp(Δη) regardless of the event momenta.
    net_factor=-np.expm1(2*np.array([eta[2]-eta[0],eta[2]-eta[1],eta[1]-eta[0]]))
    net=np.sort(np.stack((F*lo*net_factor,F*hi*net_factor)),axis=0)
    stable=2*np.array([net[0,1]-net[1,2],net[1,1]-net[0,2]])
    scalar=2*(net[:,0]+net[:,2])
    return dict(T_dark_GeV=T,eta_light_scalar_heavy=eta.tolist(),
                heavy_count_fraction=float(N[2]/N.sum()),absolute_moment_relative_defect=float(abs(E/energy-1.)),
                quantum_net_light_scalar_rate_lower_upper=[stable.tolist(),scalar.tolist()],
                quantum_net_channel_rate_lower_upper=net.T.tolist(),
                relative_quantum_event_upper_to_MB=hi.tolist(),
                relative_quantum_event_lower_to_MB=lo.tolist())


def independent_ls_forward_check(T):
    p=np.r_[0.,np.geomspace(.0004,2.7,66)]
    z=pair_mb(T,64)[2]
    defects=[]
    for radial,angle in ((2,12),(3,18)):
        _,a,_=gross(p,m,mp,lambda E:np.exp(-E/T),np.zeros_like,False,radial,angle,
                    12,8,split_final=False,exchange_mass=M)
        defects.append(abs(a[0]/z-1))
    assert max(defects)<4e-5
    return dict(T_dark_GeV=T,relative_defects_two_angle_grids=defects)


def run(eos,cache,ref):
    B=StandardModel(eos)
    index={(z["TR_GeV"],z["T_GeV"]):z for z in json.loads(ref.read_text())["rows"]}
    records=[]
    for tr,t in ((.05,.035),(.1,.085),(.1,.075),(.1,.06),(.1,.05),(.1,.035)):
        d=json.loads((cache/f"{tr:.3f}-{t:.3f}-171.json").read_text())
        vis,_,S,_=B.state(t)
        H=np.sqrt(vis/3)/MPL
        Y=np.array((sum(d["stable_yields"]),d["stable_yields"][1]))
        e=np.array(d["redshifted_means"])
        energy=S*(e[:2]@d["stable_yields"]+e[2]*Y[1])
        T,N,z=closed(Y*S,energy,t)
        c0=pair_mb(T,24);c1=pair_mb(T,52)
        np.testing.assert_allclose(c0,c1,rtol=8e-5,atol=0)
        D,f,r=signed(T,z,52)
        matched=index[(tr,t)]["controls"][-1]
        old=matched["entropy_reference_extra_Q_per_H_per_entropy"]
        new=2*f[:2].sum()/(H*S)
        assert abs(old/new-1.)<.08
        ini=dict(TR_GeV=tr,T_visible_GeV=t,T_dark_MB_GeV=T,
            initial_fugacity_light_scalar_heavy=z.tolist(),
            initial_heavy_number_fraction=float(N[2]/N.sum()),
            raw_chemical_forward_and_reverse_rates_per_entropy_H=[(q/(H*S)).tolist() for q in (f,r)],
            net_fast_invariant_rates_per_entropy_H=(D/(H*S)).tolist(),
            signed_scalar_count_divided_by_gross_heavy_scalar_source=float(D[1]/(2*f[0])),
            gross_heavy_MB_over_independent_quantum_shape=new/old,
            quantum_snapshot=enclose_quantum_snapshot(Y*S,energy,T,z,H*S))
        # A fixed-volume, fixed-energy diagnostic. Do not extrapolate this
        # across an expansion epoch, a new visible source, or a condensation branch.
        if (tr,t)==(.1,.06):
            def rhs(x,yy,q):
                U,_,fug=closed(yy*S,energy,t)
                return signed(U,fug,q)[0]/(H*S)
            path=[]
            for q,err in ((24,2e-7),(52,2e-9)):
                grid=np.array((0.,.0003,.001,.003,.01,.02))
                sol=solve_ivp(lambda x,u:rhs(x,u,q),(0.,.02),Y,method="DOP853",
                              t_eval=grid,rtol=err,atol=err*1e-7,max_step=.002)
                assert sol.success and min(sol.y.ravel())>0
                rows=[]
                for x,v in zip(grid,sol.y.T,strict=True):
                    U,pop,fug=closed(v*S,energy,t)
                    speed,ff,rr=signed(U,fug,q)
                    assert abs(pop@thermo(U)[1]/energy-1.)<4e-12
                    rows.append(dict(t_H=x,light_scalar_fast_invariant_yields=v.tolist(),
                                     actual_light_scalar_heavy_yields=(pop/S).tolist(),
                                     dark_T=U,fugacity=fug.tolist(),net_rates_over_Hs=(speed/(H*S)).tolist()))
                path.append(dict(quadrature=q,tolerance=err,rows=rows))
            a,b=[np.array([r["light_scalar_fast_invariant_yields"] for r in p["rows"]]) for p in path]
            assert np.max(abs(a/b-1))<1e-4
            ini["short_comparison"]=path
        records.append(ini)
    return dict(independent_LS_forward_check=independent_ls_forward_check(.05432665983722525),
                scope="Selected vacuum channels only, common dilute MB T and finite heavy population fixed by the "
                      "fast decay affinity eta_H=eta_L+eta_scalar. Uses full inverse detailed balance in selected "
                      "pair channels, but no visible sources, elastic evolution, expansion, Bose condensation "
                      "or omitted dark processes during the very short fixed-volume test. Not a relic solver.",
                pair_order=["light pair to heavy (opposite plus both same-charge arms)",
                            "scalar pair to heavy","light-antilight to scalar pair"],rows=records)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for s in ("eos","cache","reference"):
        p.add_argument(s,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.cache,a.reference),indent=2,allow_nan=False))
