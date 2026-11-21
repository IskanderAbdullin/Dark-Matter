"""NONTHERMAL illustrative colder Bose phase-space spectra at fixed m=86 eV and N.
f_alpha(q)=alpha^(-3) (2pi)^(-3)/(exp(q)-1), q=p/(alpha*T_thermal).
The number density is unchanged, but early radiation and velocities scale with alpha.
This family is NOT derived from Bose-stimulated production or a particle Lagrangian.
Run: python cosmology/cold_psd_scan.py --analytic-only
Full CLASS: python cosmology/cold_psd_scan.py --output outputs/cold_PSD_CLASS.json
"""
from __future__ import annotations
import argparse,json,math,tempfile
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from scipy.special import zeta
from thermal_bose_class import thermal_background, params, class_spectra, first_cross

def cold_psd_file(path,alpha):
    assert 0<alpha<=1
    q=np.unique(np.r_[np.geomspace(1e-5,.1,620),np.linspace(.1,40.,2200)])
    f=alpha**-3/((2*math.pi)**3*np.expm1(q))
    np.savetxt(path,np.c_[q,f],fmt="%.17e") # CLASS fscanf: NO nonnumeric header!
    number=(alpha**3)*4*math.pi*trapezoid(q*q*f,q)
    radiation=(alpha**4)*4*math.pi*trapezoid(q**3*f,q)
    return dict(alpha=alpha,
       number_ratio=float(number/(zeta(3)/math.pi**2)),
       early_radiation_ratio=float(radiation/(math.pi**2/30)),
       expected_radiation_ratio=alpha,
       q_npoints=len(q))

def calculate(analytic_only=False):
    b=thermal_background()
    ks=np.geomspace(.002,12.,95);zs=[0.,5.]
    out={"model":"Toy number-preserving nonthermal Bose phase-space family; NOT a realistic cold-production model",
         "base_thermal_relic":b,
         "exact_scalings":["N(alpha)=N(1) at fixed m=86eV",
             "rho_rel(alpha)=alpha*rho_rel(1), DeltaNeff(alpha)=alpha*DeltaNeff(1)",
             "mean_momentum(alpha)=alpha*mean_momentum(1)"],
         "scenarios":[]}
    with tempfile.TemporaryDirectory() as folder:
        if not analytic_only:
            from classy import Class
            cold,bgc=class_spectra(Class,params(b,cold_omega=b["omega_h2"]+1e-8),ks,zs)
            out["cold_reference"]=bgc
            out["k_1_per_Mpc"]=list(map(float,ks))
        for alpha in [1.,.5,.3,.1]:
            file=Path(folder)/("Bose_alpha_"+str(alpha)+".dat")
            mom=cold_psd_file(file,alpha)
            assert abs(mom["number_ratio"]-1)<1e-4
            assert abs(mom["early_radiation_ratio"]/alpha-1)<1e-4
            row={"alpha":alpha,"moment_tests":mom,
                "expected_omega_h2":b["omega_h2"],
                "expected_delta_Neff":alpha*b["delta_Neff_early"]}
            if not analytic_only:
                p=params(b,file)
                p["T_ncdm"]=alpha*b["T_ratio"]
                warm,bgw=class_spectra(Class,p,ks,zs)
                assert abs(bgw["omega_ncdm_h2"]/b["omega_h2"]-1)<.02,(alpha,bgw)
                row["CLASS_background"]=bgw
                row["ratio_to_equal_density_CDM"]={}
                for z in zs:
                    T=np.sqrt(np.asarray(warm[str(z)])/np.asarray(cold[str(z)]))
                    assert np.all(np.isfinite(T)) and np.all(T>0)
                    row["ratio_to_equal_density_CDM"][str(z)]={
                        "T":list(map(float,T)),
                        "k_amplitude_0p5_1_per_Mpc":first_cross(ks,T),
                        "k_shape_0p5_1_per_Mpc":first_cross(ks,T/T[0])}
            out["scenarios"].append(row)
        if not analytic_only:
            first=out["scenarios"][0]["ratio_to_equal_density_CDM"]["5.0"]["k_amplitude_0p5_1_per_Mpc"]
            assert first is not None and abs(first/.43275-1)<.03
            out["limitations"]=["Toy PSD shape is not an actual freeze-in/condensation calculation.",
              "Linear Boltzmann solver only, no CMB or Lyman-alpha likelihood or formation dynamics.",
              "No actual kinetic equilibrium or star-formation calculation for a galaxy."]
    return out

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--analytic-only",action="store_true")
    parser.add_argument("--output",default="outputs/cold_PSD_CLASS.json")
    args=parser.parse_args()
    result=calculate(args.analytic_only)
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps([{"alpha":s["alpha"],
             "expected_delta_Neff":s["expected_delta_Neff"],
             "khalf_z5":s.get("ratio_to_equal_density_CDM",{}).get("5.0",{}).get("k_amplitude_0p5_1_per_Mpc")}
          for s in result["scenarios"]],indent=2))
    print("saved",path)
