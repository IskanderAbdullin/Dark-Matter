"""Compute and separate shape and absolute-count effects as the direct source clock changes."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import roots_legendre

from archived_collisionless_streaming_reference import Background
from empty_lepton_pair_spectrum_reference import Source
from export_finite_cooling_class_psd_reference import smooth
from kloe_pion_band_portal_reference import table
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M, mp, Es, P
from warm_direct_vacuum_scalar_psd_reference import pion_parents, scalar_cumulative
from warm_light_occupation_stage_audit_reference import read

TOPS=(.025,.05,.075,.1)
EDGES=np.r_[0.,np.geomspace(.001,40.,100),np.inf]


def integration(bg,measured,order):
    src=Source(((9,28,24),(14,46,32))[order],tail=75.,dark_mass=M)
    nodes,weights=roots_legendre((15,26)[order])
    spectra=np.zeros((len(TOPS),2,len(EDGES)))
    moments=np.zeros((len(TOPS),2))
    for i,(l,r) in enumerate(zip(np.log((.006,*TOPS[:-1])),np.log(TOPS),strict=True)):
        for u,du in zip((r+l)/2+(r-l)*nodes/2,(r-l)*weights/2,strict=True):
            t=np.exp(u)
            a=bg.scale(t)
            Tq=bg.T0/a
            sh=(2*np.pi**2/45)*t**5*np.exp(bg.gs(u)+.5*bg.gr(u))*np.sqrt(np.pi**2/90)/MPL
            F=du*(1+bg.gs.derivative()(u)/3)/sh
            lp=src.particles(t,(26,44)[order])
            pi=pion_parents(t,measured,*((5,24,18),(8,36,26))[order])
            for j,(E,W) in enumerate((lp,pi)):
                C=scalar_cumulative(E,W,EDGES*Tq)
                p2=(E-M)*(E+M)/M**2
                mom=(Es*E/M)**2+P**2*p2/3-mp**2
                assert np.all(mom>=0.) and abs(C[-1]/W.sum()-1)<1e-13
                spectra[i,j]+=F*C
                moments[i,j]+=F*(W@mom)/Tq**2
    return spectra.cumsum(axis=0),moments.cumsum(axis=0)


def configs(template,root,tag,omega,file,case):
    c={k.strip():v.strip() for k,v in (line.split("=",1) for line in template.read_text().splitlines())}
    c.update(omega_cdm=format(.12-omega,'.17g'),omega_ncdm=format(omega,'.17g'),
             root=str((root/(tag+"_"+case+"_")).resolve()),
             ncdm_psd_filenames=str(file.resolve()),ncdm_maximum_q="45")
    (root/(tag+"_"+case+".ini")).write_text("".join(f"{k} = {v}\n" for k,v in c.items()))


def run(eos,kloe,ladder,template,root):
    root.mkdir(parents=True,exist_ok=True)
    bg=Background(eos)
    measured=table(kloe)
    c0,m0=integration(bg,measured,0)
    c,m=integration(bg,measured,1)
    b=read(ladder)
    idx={v["T_reheat_GeV"]:v for v in b["rows"]}
    nodes,weights=roots_legendre(5)
    lo,hi=EDGES[:-1],np.where(np.isfinite(EDGES[1:]),EDGES[1:],55.)
    q=((hi+lo)[:,None]+(hi-lo)[:,None]*nodes)/2
    rows=[]
    for j,t in enumerate(TOPS):
        d0=c0[j]/c0[j,:,-1,None];d=c[j]/c[j,:,-1,None]
        ex=c[j,:,-1]*b["Y_to_Omega_h2"]
        expect=np.asarray(idx[t]["central_lepton_and_KLOE_Omega_h2"])
        assert np.max(abs(ex/expect-1))<2e-4
        numerical=float(np.max(abs(d-d0)))
        assert numerical<.012
        counts=np.maximum(np.diff(c[j],axis=1),0.)
        w=counts.sum(axis=0)[:,None]*(weights/2)/c[j,:,-1].sum()
        assert abs(w.sum()-1)<1e-13
        keep=w.ravel()>1e-17
        x,f=smooth(q.ravel()[keep],w.ravel()[keep],.07,4200)
        f/=4*np.pi*np.trapezoid(x*x*f,x)
        assert np.min(CubicSpline(x,f)(x[:-1]+.5*np.diff(x)))>0
        moment=float(4*np.pi*np.trapezoid(x**4*f,x))
        exact=float(m[j].sum()/c[j,:,-1].sum())
        cap=x<=45.
        assert abs(moment/exact-1.)<.02
        assert 1.-4*np.pi*np.trapezoid(x[cap]**4*f[cap],x[cap])/moment<2e-6
        tag=f"direct{t*1000:.0f}"
        file=root/(tag+".dat")
        np.savetxt(file,np.c_[x,f],fmt="%.16e")
        for case,omega in (("trace",1.2e-11),("absolute",sum(ex)),
                           ("matched50",sum(idx[.05]["central_lepton_and_KLOE_Omega_h2"]))):
            configs(template,root,tag,omega,file,case)
        rows.append(dict(T_reheat_GeV=t,omega_h2_per_source=ex.tolist(),exact_p_rms_over_T0=np.sqrt(exact),
                         rms_over_T0_per_source=np.sqrt(m[j]/c[j,:,-1]).tolist(),
                         shell_lepton_then_pion_yields=counts.tolist(),smoothed_mean_q2=moment,
                         normalized_CDF_convergence=numerical,PSD_sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
    return dict(scope="Post-production prompt decay and positive comoving spectra for central vacuum charged leptons "
                      "and the measured charged-pion window only. Native configurations either retain the direct "
                      "absolute yield or match the .05 GeV scalar abundance while fixing total dark density. "
                      "No stable production and no nonlinear dark kinetics.",
                q_over_T0_edges=[float(x) if np.isfinite(x) else None for x in EDGES],rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("eos","kloe","ladder","template","root"):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.kloe,a.ladder,a.template,a.root),indent=2,allow_nan=False))
