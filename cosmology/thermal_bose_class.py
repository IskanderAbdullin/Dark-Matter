"""Linear CLASS benchmark: ONE relativistically decoupled real 86-eV Bose scalar.
No observational likelihood, no nonlinear structure, no galaxy Bose kinetics.
"""
from __future__ import annotations
import argparse,json,math,platform,tempfile
from pathlib import Path
import numpy as np
from scipy.constants import G,Boltzmann,hbar,c,parsec,elementary_charge
from scipy.integrate import trapezoid
from scipy.special import zeta

def thermal_background(m=86.,gdec=55.,gtoday=3.909,T=2.7255):
    tr=(gtoday/gdec)**(1/3)
    ng=2*zeta(3)/math.pi**2*(Boltzmann*T/(hbar*c))**3
    n=ng*tr**3/2
    H100=1e5/(1e6*parsec)
    rhoc=3*H100**2/(8*math.pi*G)
    return dict(m_eV=m,T_ratio=tr,n_scalar_cm3=n/1e6,
        omega_h2=m*elementary_charge/c**2*n/rhoc,
        delta_Neff_early=4/7*(tr/(4/11)**(1/3))**4,
        assumptions="one real BE species, early relativistic freezeout mu=0 and no BEC")

def write_psd(path,qmin=1e-5,qmax=40,nlog=620,nlin=2200):
    q=np.unique(np.r_[np.geomspace(qmin,.1,nlog),np.linspace(.1,qmax,nlin)])
    f=1/((2*math.pi)**3*np.expm1(q))
    np.savetxt(path,np.c_[q,f],fmt="%.17e")  # CLASS fscanf expects ONLY two numeric columns; NO header
    num=4*math.pi*trapezoid(q*q*f,q)
    rad=4*math.pi*trapezoid(q**3*f,q)
    return dict(qmin=qmin,qmax=qmax,npoints=len(q),
        relative_number_moment_error=float(num/(zeta(3)/math.pi**2)-1),
        relative_radiation_moment_error=float(rad/(math.pi**2/30)-1))

def params(b,psd=None,cold_omega=None):
    d={"h":.674,"omega_b":.02237,"T_cmb":2.7255,"N_ur":3.044,
       "A_s":2.1e-9,"n_s":.965,"tau_reio":.0544,
       "gauge":"newtonian","output":"mPk","P_k_max_h/Mpc":35.,
       "z_max_pk":5.1}
    if psd is None:d["omega_cdm"]=float(cold_omega)
    else:d.update({"N_ncdm":1,"m_ncdm":b["m_eV"],"T_ncdm":b["T_ratio"],
       "deg_ncdm":1,"use_ncdm_psd_files":1,
       "ncdm_psd_filenames":str(psd),"omega_cdm":1e-8})
    return d

def class_spectra(Class,p,ks,zs):
    m=Class()
    try:
        m.set(p);m.compute()
        h=float(m.h)
        return ({str(z):[float(m.pk_lin(float(k),float(z))) for k in ks] for z in zs},
           dict(omega_ncdm_h2=float(m.Omega_nu)*h*h,
                omega_m_h2=float(m.Omega_m)*h*h,
                h=h,z_eq=float(m.z_eq),CLASS_Neff=float(m.Neff)))
    finally:
        try:m.struct_cleanup()
        except Exception:pass
        try:m.empty()
        except Exception:pass

def first_cross(k,t,level=.5):
    k=np.asarray(k);t=np.asarray(t)
    ii=np.where((t[:-1]>=level)&(t[1:]<level))[0]
    if not len(ii):return None
    j=int(ii[0]);z=(level-t[j])/(t[j+1]-t[j])
    return float(np.exp(np.log(k[j])*(1-z)+np.log(k[j+1])*z))

def calculate(analytic_only=False):
    b=thermal_background();ks=np.geomspace(.002,12,95)
    with tempfile.TemporaryDirectory() as folder:
        f1=Path(folder)/"be_fine.dat";a=write_psd(f1)
        assert abs(a["relative_number_moment_error"])<1e-4
        assert abs(a["relative_radiation_moment_error"])<1e-4
        out={"benchmark":b,"file_PSD_moments":a,
             "normalization":"CLASS f0(q)=(2pi)^(-3)/(exp(q)-1), ONE real scalar, deg_ncdm=1",
             "assumptions":"cold comparator has same nonbaryonic omega at z=0 but no thermal relic radiation"}
        if not analytic_only:
            from classy import Class
            z=[0.,5.]
            warm,bgw=class_spectra(Class,params(b,f1),ks,z)
            assert abs(bgw["omega_ncdm_h2"]/b["omega_h2"]-1)<.015,(bgw,b)
            cold,bgc=class_spectra(Class,params(b,cold_omega=bgw["omega_ncdm_h2"]+1e-8),ks,z)
            ratios={}
            for zz in z:
                T=np.sqrt(np.asarray(warm[str(zz)])/np.asarray(cold[str(zz)]))
                assert np.all(np.isfinite(T)) and np.all(T>0)
                ratios[str(zz)]={"amplitude":list(map(float,T)),
                    "normalized_shape":list(map(float,T/T[0])),
                    "crossing_T_0p5_1_per_Mpc":first_cross(ks,T),
                    "crossing_shape_0p5_1_per_Mpc":first_cross(ks,T/T[0])}
            f2=Path(folder)/"be_coarse.dat";coarse=write_psd(f2,3e-5,35,320,1000)
            alt,ba=class_spectra(Class,params(b,f2),ks,[5.])
            rel=np.abs(np.asarray(alt["5.0"])/np.asarray(warm["5.0"])-1)
            out.update({"python":platform.python_version(),"k_1_per_Mpc":list(map(float,ks)),
                "warm_background":bgw,"cold_background":bgc,
                "power_amplitude_ratios":ratios,
                "PSD_convergence":{"coarse_file":coarse,
                    "max_relative_P_z5":float(np.max(rel)),
                    "median_relative_P_z5":float(np.median(rel)),
                    "coarse_omega_ncdm_h2":ba["omega_ncdm_h2"],
                    "max_rel_k_below_0p3":float(np.max(rel[ks<=.3])),
                    "max_rel_k_below_1":float(np.max(rel[ks<=1.])),
                    "worst_k":float(ks[int(np.argmax(rel))]),
                    "relative_errors":list(map(float,rel))},
                "limitation":"Linear Einstein-Boltzmann benchmark ONLY; no Lyman-alpha/CMB likelihood, cold stimulated production, boson stars or late halo equilibration"})
            assert np.all(np.isfinite(rel)), "Nonfinite q-grid P convergence"
    return out

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--analytic-only",action="store_true")
    parser.add_argument("--output",default="outputs/thermal_bose_CLASS.json")
    a=parser.parse_args()
    result=calculate(a.analytic_only)
    p=Path(a.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items()
            if k not in ("k_1_per_Mpc","power_amplitude_ratios")},indent=2))
    print("saved",p)
