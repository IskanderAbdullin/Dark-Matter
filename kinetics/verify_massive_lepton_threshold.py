"""Reproduce conditional massive e± Pauli heat and independent electron FD number ratio."""
import argparse,json,math
from pathlib import Path
from scipy.integrate import quad
from microreversible_massive_lepton import thermal_heat_example,ELECTRON_MASS_GEV

def number_ratio(me,T):
    mass=quad(lambda x:x*x/(math.exp(math.hypot(x,me/T))+1),0,45)[0]
    zero=quad(lambda x:x*x/(math.exp(x)+1),0,45)[0]
    return mass/zero

def audit(n=80):
    cases=[]
    for T in (.0002,.0001,.00005):
        a=thermal_heat_example(T=T,m=0.,n=n)
        b=thermal_heat_example(T=T,m=ELECTRON_MASS_GEV,n=n)
        ratio=b["Qchi"]/a["Qchi"]
        assert 0<ratio<1 and b["entropy"]>0
        cases.append({"Te_GeV":T,"me_over_T":ELECTRON_MASS_GEV/T,
            "Q_massless_GeV5":a["Qchi"],"Q_massive_GeV5":b["Qchi"],
            "heat_ratio":ratio,
            "independent_thermal_pair_density_ratio":number_ratio(ELECTRON_MASS_GEV,T)})
    assert cases[2]["heat_ratio"]<cases[1]["heat_ratio"]<cases[0]["heat_ratio"]
    assert cases[1]["heat_ratio"]>cases[1]["independent_thermal_pair_density_ratio"]
    return {"scope":"charge symmetric FD e± only, not physical SM e+ annihilation and baryon electrons",
            "Mchi_GeV":.01,"G_GeVminus2":.001,"Tchi_over_Te":.8,
            "mu_chi":"Mchi-2Te","energy_nodes":n,"cases":cases,
            "no_86eV_Bose_relic_or_CLASS":True}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--n",type=int,default=80)
    p.add_argument("--output",default="")
    args=p.parse_args()
    result=audit(args.n)
    if args.output:
        f=Path(args.output)
        f.parent.mkdir(parents=True,exist_ok=True)
        f.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
