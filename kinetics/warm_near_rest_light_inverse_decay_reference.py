"""Mode-resolved H <-> light + scalar rate with a stable fixed-light rest limit."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit, roots_legendre

from joint_visible_heavy_adiabatic_reference import parent_branch, snapshot_functions
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M, m, mp, P, width
from warm_light_occupation_stage_audit_reference import digest, read

E0=(M*M-m*m-mp*mp)/(2*m)
Q0=M*P/m


def window(p):
    e=np.hypot(p,m)
    center=e*(E0/m)
    half=p*(Q0/m)
    hi=center+half
    # The rational lower endpoint also holds away from the rest shell.
    lo=(E0**2+(p*mp/m)**2)/hi
    return e,center,half,lo,hi


def rate(p, light, scalar, heavy, gates=(), order=12):
    """Return spontaneous/stimulated/backward occupations averaged over angle."""
    p=np.atleast_1d(p)
    e,c,d,_,_=window(p)
    u=np.divide(np.asarray(gates)[None,:]-c[:,None],d[:,None],
                out=np.zeros((len(p),len(gates))),where=d[:,None]>0.)
    u=np.sort(np.clip(np.c_[-np.ones_like(p),u,np.ones_like(p)],-1.,1.),axis=1)
    x,w=roots_legendre(order)
    dz=np.diff(u,axis=1)/2
    mid=(u[:,1:]+u[:,:-1])/2
    ep=c[:,None,None]+d[:,None,None]*(mid[:,:,None]+dz[:,:,None]*x)
    weight=dz[:,:,None]*w/2
    assert np.max(abs(weight.sum(axis=(1,2))-1))<3e-15
    f=light(e)[:,None,None]
    B=scalar(ep)
    H=heavy(ep+e[:,None,None])
    values=np.array([np.sum(Z*weight,axis=(1,2)) for Z in
                      (H*(1-f),H*(1-f)*B,(1-H)*f*B)])
    factor=(M/m)**3*m/e
    return values*factor, np.max(B,axis=(1,2)), np.max(H,axis=(1,2)), np.sum(B*weight,axis=(1,2))


def verify():
    T=.032
    light=lambda e:expit(-e/T)
    scalar=lambda e:1/np.expm1(e/T)
    heavy=lambda e:expit(-e/T)
    p=np.r_[0.,np.geomspace(1e-8,.35,36)]
    K,*_=rate(p,light,scalar,heavy,order=80)
    mask=K[2]>1e-250
    thermal_err=max(abs((K[0]+K[1])[mask]/K[2,mask]-1))
    assert thermal_err<4e-14
    # Empty daughters: integrate light number and light total energy from
    # the inverse geometry, independently against a parent-coordinate fold.
    zero=lambda e:np.zeros_like(e)
    p=np.r_[0.,np.geomspace(1e-12,.5,31)]
    K,*_=rate(p,zero,zero,heavy,order=48)
    err=abs(K[0,0]/((M/m)**3*heavy(np.array([m+E0]))[0])-1)
    assert err<2e-15
    def by_light(n):
        def f(p):
            K,*_=rate([p],zero,zero,heavy,order=64)
            return p*p*np.hypot(p,m)**n*K[0,0]
        return quad(f,0.,2.,epsabs=1e-17,epsrel=1e-9)[0]
    def by_heavy(n):
        eL=(M*M+m*m-mp*mp)/(2*M)
        return quad(lambda p: p*p*heavy(np.hypot(p,M))*M/np.hypot(p,M)*
                    (np.hypot(p,M)*eL/M)**n,0.,2.,epsabs=1e-17,epsrel=1e-9)[0]
    check=[abs(by_light(i)/by_heavy(i)-1) for i in (0,1)]
    assert max(check)<1e-8
    return dict(thermal_pointwise_defect=thermal_err,rest_limit_defect=err,
                integrated_spontaneous_number_and_energy_relative_defects=check)


def freeze(s,eos,outer,inner,cut=None,crop_heavy_inputs=False):
    B,bg,F,fg=snapshot_functions(s)
    p=s["lab_momentum_edges_GeV"]
    pg=np.linspace(0.,np.nextafter(p[1],0.),21)
    T=s["T_visible_GeV"]
    Hubble=np.sqrt(TabulatedRadiation(eos,T).state(1.)[1]/3)/MPL
    core_B,core_bg=B,bg
    if cut is not None:
        ecut=np.hypot(cut,mp)
        core_B=lambda E:np.where(E<ecut,B(E),0.)
        core_bg=np.r_[bg[bg<ecut],ecut]
    def heavy(E):
        A=parent_branch(E.ravel(),core_B if crop_heavy_inputs else B,F,
                        core_bg if crop_heavy_inputs else bg,fg,angle=inner)
        return A[:,0].reshape(E.shape)
    K,fB,fH,avg_B=rate(pg,F,core_B,heavy,core_bg,order=outer)
    lo,hi=window(pg)[-2:]
    pop=F(window(pg)[0])
    net=K[0]+K[1]-K[2]
    A=np.sum(K[:2,0])/(1-pop[0])
    L=K[2,0]/pop[0]
    assert A+L>0
    return dict(local_scalar_energy_window_GeV=[float(np.min(lo)),float(np.max(hi))],
                widest_known_scalar_E_GeV=core_bg[-1],sample_max_scalar_occupation=float(fB.max()),
                sample_max_frozen_heavy_occupation=float(fH.max()),
                light_f_at_rest=float(pop[0]), cold_shell_edge_GeV=p[1],
                instantaneous_local_fixed_H_and_B_equilibrium_at_rest=A/(A+L),
                light_linear_relaxation_at_rest_over_Hubble=(A+L)*width/Hubble,
                sufficient_fH_uniform_bound_for_half_at_rest=float(avg_B[0]/(1+2*avg_B[0])),
                H_over_rest_decay_width=float(Hubble/width),
                spontaneous_stimulated_and_inverse_at_rest_per_Hubble=(K[:,0]*(width/Hubble)).tolist(),
                net_cold_shell_rate_per_Hubble_range=(np.array([net.min(),net.max()])*width/Hubble).tolist())


def run(eos,D):
    rows=[]
    q16_cuts=[]
    for q in (16,32):
        c=read(D/f"warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz")
        c+=[u["checkpoint"] for u in
            read(D/f"warm-observed-suboperator-q{q}-tape-2026-11-26.json.gz")["blocks"]]
        stages=[]
        for j,u in enumerate(c):
            s=u["snapshot"]
            assert u["snapshot_sha256"]==digest(s)
            a=freeze(s,eos,8,13)
            b=freeze(s,eos,14,21)
            r=np.asarray(a["net_cold_shell_rate_per_Hubble_range"])
            z=np.asarray(b["net_cold_shell_rate_per_Hubble_range"])
            err=float(np.max(abs(r-z))/max(np.max(abs(z)),1e-280))
            assert err<1e-4
            assert (b["sample_max_scalar_occupation"]==0.)==(q==16)
            assert b["spontaneous_stimulated_and_inverse_at_rest_per_Hubble"][0]>0.
            d=dict(T_GeV=s["T_visible_GeV"],snapshot_sha256=u["snapshot_sha256"],
                   fine=b,relative_net_range_order_defect=err)
            if q==16:
                q16_cuts.append((s["T_visible_GeV"],s["lab_momentum_edges_GeV"][-1]))
                assert min(z)>0.
            if q==32:
                TT,cut=q16_cuts[j]
                assert abs(s["T_visible_GeV"]/TT-1)<1e-12
                fixed=freeze(s,eos,14,21,cut=cut)
                reco=freeze(s,eos,14,21,cut=cut,crop_heavy_inputs=True)
                for f in (fixed,reco):
                    assert f["sample_max_scalar_occupation"]==0.
                    assert min(f["net_cold_shell_rate_per_Hubble_range"])>0.
                assert max(z)<0.
                d["same_q32_scalar_support_crop"] = dict(
                    physical_scalar_p_cut_GeV=cut, held_full_support_frozen_H=fixed,
                    recomputed_frozen_H_with_cropped_scalar=reco)
            stages.append(d)
        rows.append(dict(q=q,stages=stages))
    return dict(scope="Frozen local rate on the innermost light shell. Scalar and light fixed "
                      "to archived retained profiles; heavy picked as the stationary balance "
                      "of the represented daughter decay pair; not its actual trajectory or "
                      "a bound on an omitted tail or on other collision operators.",
                EOS_sha256=hashlib.sha256(eos.read_bytes()).hexdigest(),
                per_mode_cold_limit_over_rest_decay_width=(M/m)**3,
                light_rest_scalar_energy_GeV=E0,controls=verify(),cases=rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("eos",type=Path)
    parser.add_argument("--tapes",type=Path,default=Path("results"))
    ns=parser.parse_args()
    print(json.dumps(run(ns.eos,ns.tapes),indent=2,allow_nan=False))
