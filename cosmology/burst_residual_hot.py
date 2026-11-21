"""Two-population primordial 86-eV Bose initial PSD: cold burst + residual warm tail.

SINGLE particle species with a sum of positive, number-normalized momentum shapes.
This is an illustrative collisionless CLASS initial condition, NOT the self-consistent
Yin 2023 1<->2 production solution nor a true quantum kinetic simulation.
"""
from __future__ import annotations
import argparse,json,math,tempfile
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from scipy.special import zeta
from thermal_bose_class import thermal_background,params,class_spectra,first_cross

def inv_expm1(x):
    x=np.asarray(x)
    out=np.empty_like(x)
    mask=x<=90
    out[mask]=1/np.expm1(x[mask])
    out[~mask]=np.exp(-x[~mask])
    return out

def write_bimodal(path,alpha,hot_fraction,coarse=False):
    assert 0<alpha<1 and 0<=hot_fraction<=1
    if coarse:
        q=np.unique(np.r_[np.geomspace(1e-6,.005,450),
          np.geomspace(.005,.8,900),np.linspace(.8,45,2300)])
    else:
        q=np.unique(np.r_[np.geomspace(1e-7,.005,800),
          np.geomspace(.005,.8,1600),np.linspace(.8,45,3500)])
    # q=p/T_reference; same species, one REAL spin-zero degree of freedom.
    cold=(1-hot_fraction)*alpha**-3*inv_expm1(q/alpha)
    hot=hot_fraction*inv_expm1(q)
    f=(cold+hot)/(2*math.pi)**3
    np.savetxt(path,np.column_stack([q,f]),fmt="%.17e")
    number=4*math.pi*trapezoid(q*q*f,q)
    relativistic=4*math.pi*trapezoid(q**3*f,q)
    expected_radiation=(1-hot_fraction)*alpha+hot_fraction
    return {"alpha":alpha,"hot_number_fraction":hot_fraction,
      "radiation_fraction_relative_thermal":expected_radiation,
      "hot_fraction_of_early_relativistic_energy":hot_fraction/expected_radiation,
      "number_moment_ratio":float(number/(zeta(3)/math.pi**2)),
      "radiation_moment_ratio":float(relativistic/(expected_radiation*math.pi**2/30)),
      "q_points":int(len(q)),"q_min":float(q[0]),"q_max":float(q[-1])}


def class_component_file(path,fraction,alpha,coarse=False):
    """One smooth PSD in its OWN CLASS q=p/T_ncdm variable, not the mixed-q table."""
    assert 0 < fraction <= 1 and 0 < alpha <= 1
    if coarse:
        q=np.unique(np.r_[np.geomspace(3e-5,.1,320),np.linspace(.1,35.,1000)])
    else:
        q=np.unique(np.r_[np.geomspace(1e-5,.1,620),np.linspace(.1,40.,2200)])
    f=fraction*alpha**-3/((2*math.pi)**3*np.expm1(q))
    np.savetxt(path,np.column_stack((q,f)),fmt="%.17e")
    return {"fraction":fraction,"alpha":alpha,"q_points":int(len(q))}

def class_mixture_params(b,alpha,fhot,folder,coarse=False):
    """Split ONE collisionless particle PSD into two numerical CLASS tracers.

    This is an exact decomposition of linear collisionless phase-space moments,
    not the introduction of two new physical particle species. Do NOT apply it
    when significant cold/hot inter-stream collisions or conversion occur.
    """
    folder=Path(folder)
    suffix="coarse" if coarse else "fine"
    cold=folder/("cold_"+str(alpha)+"_"+str(fhot)+"_"+suffix+".dat")
    class_component_file(cold,1-fhot,alpha,coarse)
    p=params(b,psd=cold)
    p["P_k_max_h/Mpc"]=65.
    p["T_ncdm"]=alpha*b["T_ratio"]
    if fhot>0:
        hot=folder/("hot_"+str(alpha)+"_"+str(fhot)+"_"+suffix+".dat")
        class_component_file(hot,fhot,1.,coarse)
        p.update({"N_ncdm":2,"m_ncdm":"86,86",
          "T_ncdm":"%.15g,%.15g"%(alpha*b["T_ratio"],b["T_ratio"]),
          "deg_ncdm":"1,1","use_ncdm_psd_files":"1,1",
          "ncdm_psd_filenames":str(cold)+","+str(hot)})
    return p

def interpolate(ks,values,k):
    return float(np.interp(math.log(k),np.log(ks),values))

def calculate(analytic_only=False,alpha_only=None):
    b=thermal_background()
    ks=np.geomspace(.002,35,155)
    zs=[5.]
    out={"title":"Two-scale initial PSD for ONE 86-eV scalar",
      "model":"(1-fhot)alpha^-3 fBE(q/alpha)+fhot fBE(q), q=p/T_reference",
      "scope":"Illustrative free-streaming linear perturbations. NOT microscopic stimulated 1<->2 production.",
      "theoretical_motivation":"Yin JHEP05(2023)180 eq (4.12) has burst+cold and residual hot populations, but does not derive this exact PSD",
      "reference_thermal_omega_h2":b["omega_h2"],
      "reference_thermal_deltaNeff":b["delta_Neff_early"],
      "scenarios":[]}
    with tempfile.TemporaryDirectory() as folder:
        if not analytic_only:
            from classy import Class
            pc=params(b,cold_omega=b["omega_h2"]+1e-8)
            pc["P_k_max_h/Mpc"]=65.
            cold,bc=class_spectra(Class,pc,ks,zs)
            cdp=np.array(cold["5.0"])
            assert np.all(cdp>0)
            out["CDM_reference"]=bc
            out["k_1_per_Mpc"]=list(map(float,ks))
        for alpha in ((alpha_only,) if alpha_only is not None else (.02,.01)):
            for fhot in (0.,.01,.05,.10,.20):
                file=Path(folder)/("mix_"+str(alpha)+"_"+str(fhot)+".dat")
                info=write_bimodal(file,alpha,fhot)
                assert abs(info["number_moment_ratio"]-1)<.0001,info
                assert abs(info["radiation_moment_ratio"]-1)<.0001,info
                row={"alpha":alpha,"hot_fraction":fhot,
                   "shape_moments":info,
                   "expected_deltaNeff":b["delta_Neff_early"]*info["radiation_fraction_relative_thermal"]}
                if not analytic_only:
                    p=class_mixture_params(b,alpha,fhot,folder)
                    spectrum,bg=class_spectra(Class,p,ks,zs)
                    assert abs(bg["omega_ncdm_h2"]/b["omega_h2"]-1)<.01,(alpha,fhot,bg)
                    assert abs((bg["CLASS_Neff"]-3.044)-row["expected_deltaNeff"])<.0003,(alpha,fhot,bg)
                    T=np.sqrt(np.array(spectrum["5.0"])/cdp)
                    assert np.all(np.isfinite(T)) and np.all(T>0)
                    row.update({"CLASS_background":bg,
                      "k_first_power_half_1_per_Mpc":first_cross(ks,T,1/math.sqrt(2)),
                      "k_first_amplitude_half_1_per_Mpc":first_cross(ks,T,.5),
                      "P_over_CDM_k1":interpolate(ks,T*T,1.),
                      "P_over_CDM_k5":interpolate(ks,T*T,5.),
                      "P_over_CDM_k10":interpolate(ks,T*T,10.),
                      "P_over_CDM_k20":interpolate(ks,T*T,20.),
                      "transfer_amplitude":list(map(float,T))})
                    if fhot in (0.,.05,.20):
                        co=Path(folder)/("coarse_"+str(alpha)+"_"+str(fhot)+".dat")
                        coinfo=write_bimodal(co,alpha,fhot,coarse=True)
                        pp=class_mixture_params(b,alpha,fhot,folder,coarse=True)
                        alternate,_=class_spectra(Class,pp,ks,zs)
                        Pco=np.array(alternate["5.0"])
                        rel=np.abs(Pco/np.array(spectrum["5.0"])-1)
                        pre=ks<=min(15.,row["k_first_power_half_1_per_Mpc"] or 15.)
                        row["coarse_file_check"]={
                          "moments":coinfo,
                          "max_relative_P_difference_before_first_power_half":float(max(rel[pre])),
                          "relative_P_difference_k1":interpolate(ks,rel,1.),
                          "relative_P_difference_k5":interpolate(ks,rel,5.),
                          "max_relative_P_difference_all_k_not_validated":float(max(rel))}
                    print("alpha",alpha,"f_hot",fhot,"P10",row["P_over_CDM_k10"],"khalf",row["k_first_power_half_1_per_Mpc"],flush=True)
                out["scenarios"].append(row)
    return out

if __name__=="__main__":
    pa=argparse.ArgumentParser()
    pa.add_argument("--analytic-only",action="store_true")
    pa.add_argument("--alpha-only",type=float,choices=(.02,.01))
    pa.add_argument("--output",default="outputs/mixed_bose_primordial_CLASS.json")
    args=pa.parse_args()
    result=calculate(args.analytic_only,args.alpha_only)
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print("saved",path,"scenarios",len(result["scenarios"]))
