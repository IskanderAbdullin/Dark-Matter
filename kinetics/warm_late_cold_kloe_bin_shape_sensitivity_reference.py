"""Check two continuous within-band spectral choices against held KLOE bins."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.interpolate import PchipInterpolator
from archived_collisionless_streaming_reference import Background
from kloe_pion_band_portal_reference import table
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from warm_late_cold_visible_kloe_recoil_reference import rates
from warm_light_occupation_stage_audit_reference import read


def run(eos,kloe,reference):
    g,cs,e=table(kloe)
    bg=Background(eos)
    radiation=TabulatedRadiation(eos,.05)
    last=read(reference)["controls"][-1]["BE"]
    m=(PchipInterpolator(g,cs),PchipInterpolator(g,np.log(cs)))
    styles={
      "linear_s":lambda q:np.interp(q,g,cs),
      "PCHIP_sigma":lambda q:m[0](np.clip(q,g[0],g[-1])),
      "PCHIP_log_sigma":lambda q:np.exp(m[1](np.clip(q,g[0],g[-1]))),
      "published_piecewise_constant":None}
    results=[]
    for label,F in styles.items():
        vals=[]
        for t,n in ((14,8),(20,12)):
            x,w=roots_legendre(t)
            A,B=np.log([.012,.05]);h=(B-A)/2
            r=[]
            for T in np.exp((A+B)/2+h*x):
                a=bg.scale(T)
                _,rho=radiation.state(a/bg.scale(.05))
                s=2*np.pi**2/45*T**3*np.exp(bg.gs(np.log(T)))
                dt=(3+bg.gs.derivative()(np.log(T)))/3
                r.append(rates(T,.1*bg.T0/a,g,cs,e,n,spectral_profile=F)[1]*dt/(s*np.sqrt(rho/3)/MPL))
            vals.append(dict(kinetic_order=n,time_order=t,BE_yield=float(h*(w@r))))
        assert abs(vals[-1]["BE_yield"]/vals[0]["BE_yield"]-1)<1e-4
        effect=vals[-1]["BE_yield"]/last-1.
        if F is None:
            assert abs(effect)<3e-14
        results.append(dict(shape=label,relative_yield_shift_vs_stored_bin=effect,controls=vals))
    return dict(scope="Common measured invariant window only. Interpolants hold the edge-center value to the "
                      "published boundary. Alternate curves are not error bands or confidence intervals.",cases=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","kloe","reference"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.kloe,a.reference),indent=2,allow_nan=False))
