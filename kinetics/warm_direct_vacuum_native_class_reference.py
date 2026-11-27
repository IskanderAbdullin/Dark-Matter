"""Prepare matched native CLASS power runs for the independent warm scalar PSD."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import roots_legendre,zeta
from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import points
from export_finite_cooling_class_psd_reference import smooth
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import mp
from warm_light_occupation_stage_audit_reference import digest,read


def prepare(eos,archive,source,metric,root):
    root.mkdir(parents=True,exist_ok=True)
    bg=Background(eos);S=read(archive)
    assert S["snapshot_sha256"]==digest(S["snapshot"])
    A=S["snapshot"]
    R=read(source)
    p=np.asarray(A["lab_momentum_edges_GeV"])
    f,*_=snapshot_functions(A)
    q,_,w=points(p,mp,f,24)
    q*=bg.scale(A["T_visible_GeV"])/bg.T0
    w/=w.sum()
    x,d=roots_legendre(4)
    E=np.array([t if t is not None else 36. for t in R["conserved_momentum_edges_q_over_T0"]])
    qq=(E[:-1,None]*(1-x)+E[1:,None]*(1+x))/2
    n=np.asarray(R["controls"][-1]["shell_yields"]).sum(axis=0)
    nn=n[:,None]*(d/2)/n.sum()
    assert n[-1]/n.sum()<2e-9
    shapes={"direct":(qq.ravel(),nn.ravel()),"held":(q,w)}
    D={}
    for key,(q,w) in shapes.items():
        keep=w>1e-16
        assert w[~keep].sum()<1e-12
        q,w=q[keep],w[keep]
        data=[]
        for sigma in (.07,.11):
            a,b=smooth(q,w,sigma,3600)
            b*=w.sum()/(4*np.pi*np.trapezoid(a*a*b,a))
            assert np.min(CubicSpline(a,b)(a[:-1]+.5*np.diff(a)))>0.
            Q=float(4*np.pi*np.trapezoid(a**4*b,a))
            assert abs(Q/(w@(q*q))-1.)<1e-7
            name=f"{key}-sig{int(sigma*100):02d}"
            path=root/(name+".dat")
            np.savetxt(path,np.c_[a,b],fmt="%.16e")
            clips=a<=29.
            cut=float(4*np.pi*np.trapezoid(a[clips]**2*b[clips],a[clips])/w.sum())
            cut2=float(4*np.pi*np.trapezoid(a[clips]**4*b[clips],a[clips])/Q)
            assert 0.99999<cut<1+1e-11 and .99998<cut2<1+1e-11
            data.append(dict(PSD=name,mean_q2=Q,number_fraction_in_Class_q_cap=cut,
                             q2_fraction_in_Class_q_cap=cut2,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        D[key]=data
    fr=read(metric)["held_and_actual_direct_scalar_rest_fractions_with_same_archived_stable_count"]
    spec=[("base","direct",0.),("trace","direct",1e-10),("physical","direct",fr[1]),("amp4","direct",1e-4),("amp3","direct",1e-3),("amp2","direct",.01),("held3","held",1e-3),("thermal3","matched_FD",1e-3)]
    configs=[]
    for tier,eps,nq,lmax,code in (("std",2e-8,70,20,"sig07"),("fine",8e-10,128,32,"sig07"),
                                 ("shape",2e-8,70,20,"sig11")):
        for id,kind,fraction in spec:
            if tier=="shape" and id not in ("base","amp3","held3"):continue
            label=tier+"-"+id
            values={"output":"mPk", "P_k_max_h/Mpc":2.5,"z_pk":"5,20","h":.6736,"T_cmb":2.72548,
                    "omega_b":.02237,"A_s":2.1e-9,"n_s":.965,"tau_reio":.054,"N_ur":3.044,"gauge":"synchronous",
                    "non linear":"none","background_verbose":2,"perturbations_verbose":1,
                    "omega_cdm":f"{.12*(1-fraction):.17g}","write background":"no", "headers":"yes",
                    "k_per_decade_for_pk":13,"tol_perturbations_integration":eps,"root":str(root/(label+"_"))}
            if fraction>0.:
                values.update({"N_ncdm":1,"m_ncdm":mp*1e9,"T_ncdm":1,"omega_ncdm":f"{.12*fraction:.17g}",
                     "use_ncdm_psd_files":1,"ncdm_psd_filenames":str(root/(kind+"-"+code+".dat")),
                     "ncdm_quadrature_strategy":3,"ncdm_maximum_q":29,"ncdm_N_momentum_bins":nq,
                     "l_max_ncdm":lmax,"ncdm_fluid_approximation":3})
                if kind=="matched_FD":
                    ft=float(np.sqrt(D["direct"][0]["mean_q2"]/(15*zeta(5)/zeta(3))))
                    values["T_ncdm"]=ft
                    values.pop("use_ncdm_psd_files")
                    values.pop("ncdm_psd_filenames")
            (root/(label+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in values.items()))
            configs.append(dict(label=label,kind=kind,rest_fraction=fraction,rtol=eps,quadrature=nq,lmax=lmax,smoothing=code))
    return dict(scope="Conditional CLASS: keep omega_cdm+omega_ncdm=.120 with CDM proxy for all other dark matter, "
                      "same radiation/background parameters, fixed post-production scalar distribution. "
                      "No inference of the actual cosmological abundance.",inputs=D,physical_rest_fraction=fr[1],cases=configs)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","archive","source","metric","root"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(prepare(a.eos,a.archive,a.source,a.metric,a.root),indent=2,allow_nan=False))
