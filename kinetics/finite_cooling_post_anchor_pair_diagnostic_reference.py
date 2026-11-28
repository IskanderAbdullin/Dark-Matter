"""Vacuum forward-only post-anchor optical depths of frozen prompt distributions."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import amplitude
from conditional_relic_normalization_reference import conversion
from finite_cooling_collisionless_expansion_reference import moments
from finite_cooling_reheat_portal_reference import TF, SourceFold
from three_species_decay_onsager_reference import M, m, mp, y


def grid(case,edge,kind):
    q,x,w=moments(case,edge)
    N=q[1];x=x[1];w=w[1]
    keep=N>1e-15
    assert np.sum(N[~keep])<1e-12
    assert np.sum(N[~keep]*x[~keep])/(N@x)<1e-10
    if kind=="center":
        z,d=np.sqrt(x[keep]),N[keep]
    else:
        assert N[-1]<1e-15
        v=N[:-1][keep[:-1]];ww=w[keep[:-1]]
        z=np.r_[edge[:-1][keep[:-1]],edge[1:][keep[:-1]]]
        d=np.r_[v*(1-ww),v*ww]
    d=d/d.sum()
    assert abs(d@(z*z)/(N@x)-1.)<1e-8
    return z*TF,d


def kernel(p,w,order):
    i,j=np.tril_indices(len(p))
    c,g=roots_legendre(order)
    V=(p*p)/(np.hypot(m,p)+m)
    pp=p[i,None]*p[j,None]
    EF=np.hypot(p[i],m)*np.hypot(p[j],m)
    S=4*m*m+2*(m*(V[i]+V[j])[:,None]+V[i,None]*V[j,None]-pp*c)
    W=w[i]*w[j]*np.where(i==j,1.,2.)
    total=0.
    for xx,gg in zip(c,g):
        F=amplitude(S,xx)*np.sqrt(np.maximum(.25-mp*mp/S,0.))
        total+=gg*float(W@((F@(g/2))/EF))
    return total/(64*np.pi)


def p_wave():
    c,g=roots_legendre(16)
    k2=m*m-mp*mp
    a=mp*mp-m*m-M*M
    ap=4*k2*c*np.sqrt(1-c*c)/(a*a)
    zz=2*(m+M)/a+4*m*k2*c*c/(a*a)
    amp_per_p2=2*y**4*(m*m*ap*ap+zz*zz)
    return float((g@amp_per_p2)*np.sqrt(1-mp*mp/m**2)/(128*np.pi*m*m))


def run(eos,dist):
    data=json.loads(Path(dist).read_text())
    esh=hashlib.sha256(Path(eos).read_bytes()).hexdigest()
    assert esh==data["eos_sha256"]
    K=SourceFold(eos,(4,9,8),65,4)
    factor=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    edge=np.asarray(data["boundary_p_over_T"])
    r=np.log(TF)
    s0=K.s(r)
    C=p_wave()
    rows=[]
    for b in data["cases"]:
        chi=b["species"][1]["Omega_h2"]
        N0=chi/(m*factor)*s0
        Om_prompt=b["species"][0]["Omega_h2"]
        refined=[]
        for kind,order,n in (("center",7,3),("endpoints",9,4)):
            p,w=grid(b,edge,kind)
            p2=w@(p*p)
            def sample(u):
                s=K.s(u)
                a=(s/s0)**(1/3)
                k=kernel(a*p,w,order)
                pp=C*p2*a*a/2
                X=(N0/2)*(s/s0)/K.H(u)
                return np.array([X*k,X*pp,k/pp])
            ratio=[]
            for tm in (TF,.0003,.0001,.00003):
                ratio.append([float(tm),*sample(np.log(tm)).tolist()])
            assert abs(ratio[-1][-1]-1.)<5e-3
            def f(t):
                v=np.stack([sample(z)[:2] for z in t.flat]).reshape(*t.shape,2)
                return (1+K.grad(t[:,:,None])/3)*v
            l=np.log(.00003)
            old=(K.nodes,K.weights)
            K.nodes,K.weights=roots_legendre(n)
            A=K.fold(l,r,np.r_[K.bg.sg.x,np.linspace(l,r,9)],f)
            K.nodes,K.weights=old
            raw=chi*mp/m*A[0]
            refined.append(dict(quadrature=[kind,order,n],instantaneous_rows=ratio,
                 integrated_vacuum_forward_log_expansion_optical_depth=float(A[0]),
                 integrated_p_wave_optical_depth=float(A[1]),
                 perturbative_unblocked_scalar_Omega_h2_equivalent=float(raw),
                 equivalent_new_to_prompt_scalar_rest_count_ratio=float(raw/Om_prompt)))
        rel=abs(refined[0]["integrated_vacuum_forward_log_expansion_optical_depth"]/
                refined[1]["integrated_vacuum_forward_log_expansion_optical_depth"]-1.)
        assert rel<1e-2
        rows.append(dict(T_R_GeV=b["T_R_GeV"],Tmax_GeV=b["Tmax_GeV"],controls=refined,
                        optical_depth_relative_control_change=rel))
    return dict(eos_sha256=esh,archive_sha256=hashlib.sha256(Path(dist).read_bytes()).hexdigest(),
                p_wave_CoM_sigma_v_over_p_cm2_GeV_minus4=C,
                columns_for_instantaneous_rows=["T_GeV","Nchi_sigma_v_over_2H","pwave_proxy","exact_over_pwave"],
                scope="Forward reaction from unmodified frozen stable counts on prescribed post-1MeV radiation clock. Equal particle and antiparticle density, vacuum spontaneous amplitudes, no Bose enhancement, inverse, statistical blocking or depletion/repopulation. Scalar count proxy is first-order optical-depth multiplication, not a prediction or a general lower/upper bound. Excludes all earlier kinetics.",cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("distribution")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.distribution),indent=2,allow_nan=False))
