"""Passive forward-pair diagnostic during positive prompt production."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from finite_cooling_reheat_distribution_reference import deposit
from finite_cooling_reheat_portal_reference import SourceFold,TF
from finite_cooling_post_anchor_pair_diagnostic_reference import p_wave
from finite_cooling_pair_emission_reference import kernel_and_lab_p2, checks
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import m,mp

def run(eos,reference,collect_shape=False):
    ref0=json.loads(reference.read_text())
    eh=hashlib.sha256(eos.read_bytes()).hexdigest()
    assert eh==ref0["eos_sha256"]
    ref=ref0["controls"][-1]["rows"]
    K=SourceFold(eos,(7,18,14),100,4)
    F=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    thresholds=np.array([.05,.1,.2,.5,1.,2.,3.,4.,6.,8.,12.,16.,24.,40.,80.])
    recoil_edges=np.unique(np.r_[0.,thresholds,np.geomspace(.025,110,30),np.inf])
    recoils=[float(v) for v in recoil_edges[:-1]]
    U=np.log(TF)
    sf=K.s(U)
    cw=p_wave()
    out=[]
    for row in (ref[0],ref[-1]):
        r,h=np.log([row['T_R_GeV'],row['Tmax_GeV']])
        ar=(sf/K.s(r))**(1/3)
        modes=[]
        for n,order,n_bins in ((32,6,41),(64,8,55)):
            E=np.r_[0.,np.geomspace(.001,110,n_bins),np.inf]
            yt=np.zeros((2,2,len(E)-1))
            Q=np.zeros((2,2)); D=np.zeros(2)
            p2out=np.zeros(2)
            B=np.zeros((2,2,len(recoil_edges)-1))
            x,w=roots_legendre(3)
            def vals(u,pre):
                if pre:
                    a=ar*np.exp(-8*(u-r)/3)
                    H=K.H(r)*np.exp(4*(u-r))
                    dt=(8/3)/H
                else:
                    a=(K.s(u)/sf)**(-1/3)
                    H=K.H(u)
                    dt=(1+K.grad(u)/3)/H
                return a,H,dt
            def src(u,pre):
                a,H,dt=vals(u,pre)
                return deposit(K,np.exp(u),E,a,20)*(a**3*dt/sf)
            for i,(hi,lo) in enumerate(((h,r),(r,U))):
                if hi<=lo: continue
                pre=(i==0)
                nodes=np.linspace(hi,lo,n+1)
                for b,a0 in zip(nodes[:-1],nodes[1:]):
                    center=(b+a0)/2
                    for j,(l0,u0) in enumerate(((center,b),(a0,center))):
                        yt+=sum(f*src(u,pre) for u,f in zip((l0+u0)/2+(u0-l0)*x/2,(u0-l0)*w/2))
                        if j!=0:continue
                        ac,H,dt=vals(center,pre)
                        Y=yt[1,0];Z=yt[1,1]
                        keep=Y>Y.sum()*1e-14
                        p=TF/ac*np.sqrt(Z[keep]/Y[keep]); c=Y[keep]/Y[keep].sum()
                        (L,p2L),shell=kernel_and_lab_p2(p,c,order,
                            recoil_edges if collect_shape else None,ac/TF)
                        nphys=sf/ac**3*Y.sum()
                        value=nphys*L/(2*H)
                        D[i]=max(D[i],value)
                        dw=(b-a0)*dt*(nphys*Y.sum()/2)
                        Q[i]+=dw*np.array([L,cw*(c@(p*p))/2])
                        p2out[i]+=dw*p2L*(ac/TF)**2
                        if collect_shape:
                            B[i]+=dw*shell
            sec=Q[:,0].sum()
            qp=float(row["combined_rms_momentum_over_1MeV"][0])
            Yp=yt[0,0].sum()
            r2=np.divide(p2out,Q[:,0],out=np.zeros_like(p2out),where=Q[:,0]>0.)
            shape=None
            if collect_shape:
                assert np.max(abs(B.sum(axis=-1).sum(axis=0)/np.array([sec,p2out.sum()])-1.))<2e-10
                shape=dict(count_fraction_by_epoch=(B[:,0,:]/sec).tolist(),
                     q2_weighted_count_fraction_by_epoch=(B[:,1,:]/sec).tolist(),
                     fraction_below_thresholds=(np.cumsum(B[:,0,:].sum(axis=0))[np.searchsorted(recoil_edges,thresholds)-1]/sec).tolist())
            modes.append(dict(steps=n,incoming_momentum_bins=len(E)-1,
                new_scalar_p_over_TF_rms_by_epoch=np.sqrt(r2).tolist(),
                first_order_mixture_scalar_rms_p_over_TF=float(np.sqrt((Yp*qp*qp+p2out.sum())/(Yp+sec))),
                first_order_pair_fraction_of_scalar_squared_momentum=float(p2out.sum()/(Yp*qp*qp+p2out.sum())),
                new_scalar_recoil_shape=shape,shape_count_omega=F*yt[:,0].sum(axis=1)*np.array([mp,m]),
                        maximum_forward_loss_over_H=D,epoch_order=['pre','post'],new_scalar_omega=F*mp*Q,
                        fraction_of_final_stable_equivalent=Q[:,0]/yt[1,0].sum()))
        for mode in modes:
            assert np.max(np.abs(mode["shape_count_omega"]/np.array([row["scalar_prompt_fixed_count_Omega_h2"],row.get("expected_gross_stable_Omega_h2",.120)])-1.))<6e-5
            mode["new_scalar_vacuum_to_total_prompt_ratio"]=float(mode["new_scalar_omega"][:,0].sum()/mode["shape_count_omega"][0])
        a,b=[v["new_scalar_omega"][:,0] for v in modes]
        j=b>0.
        change=float(np.max(np.abs(a[j]/b[j]-1.)))
        assert change<.006
        out.append(dict(T_R_GeV=row['T_R_GeV'],Tmax_GeV=row['Tmax_GeV'],modes=modes,max_relative_forward_change=change))
    return dict(eos_sha256=eh,angular_checks=checks(),reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
            temperature_panel_quadrature_order=3,recoil_p_over_TF_lower_edges=recoils if collect_shape else None,
        recoil_shape_probe_thresholds=thresholds.tolist() if collect_shape else None,
            episodes=out,scope="First-order spontaneous pair conversion only, evaluated on the evolving unmodified prompt frozen source. No inverse, stimulation, dark rescattering, depletion, external plasma, or heating branch. Not a coupled outcome or a lower or upper bound on the complete model.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos",type=Path)
    p.add_argument("reference",type=Path)
    p.add_argument("--shape",action="store_true")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.reference,a.shape),default=lambda x:x.tolist(),indent=2,allow_nan=False))
