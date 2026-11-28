"""Fixed anchor signed pair probes in two positive reconstructions."""
import argparse
import hashlib
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline

from archived_scalar_pair_reaction_reference import gross,points
from archived_collisionless_streaming_reference import Background
from finite_cooling_post_anchor_pair_diagnostic_reference import kernel
from three_species_decay_onsager_reference import m,mp

def table(root,record,qscale,mass,g):
    path=root/record["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest()==record["sha256"]
    q,F=np.loadtxt(path).T
    sp=CubicSpline(q,(2*np.pi)**3*F/g,extrapolate=False)
    check=np.r_[q,sp.derivative().roots(extrapolate=False)]
    mx=float(np.nanmax(sp(check)))
    def f(e):
        Q=np.sqrt(np.maximum(0.,(e-mass)*(e+mass)))*qscale
        return np.nan_to_num(sp(Q),nan=0.)
    return q[-1],f,mx


def run(eos,root,manifest):
    data=json.loads(manifest.read_text())
    assert hashlib.sha256(eos.read_bytes()).hexdigest()==data["eos_sha256"]
    bg=Background(eos)
    qscale=float(bg.scale(.001))/bg.T0
    def get(tag,name,sig,mass,g):
        a=next(x for x in data["components"] if (x["case"],x["species"],x["sigma_ln_q"])==(tag,name,sig))
        return table(root,a,qscale,mass,g)
    output=[]
    for tag in ('instant','ratio5'):
     for sig in (.10,.14):
      qs,fs,mx=get(tag,'scalar',sig,mp,1)
      qc,fc,cx=get(tag,'stable',sig,m,4)
      r=3
      edges=np.r_[0.,np.geomspace(.000005,1.001*max(qc,qs)/qscale,47)]
      p,_,nf=points(edges,m,fc,r)
      e=kernel(p,nf/nf.sum(),10)*(2*nf.sum())**2
      n0,u,_=gross(edges,m,mp,fc,fs,False,r,10,14,10,split_final=False)
      _,v,_=gross(edges,mp,m,fs,fc,True,r,10,14,10,split_final=False)
      rel=abs(u[0]/e-1.)
      assert abs(n0/(4*nf.sum())-1.)<1e-13 and rel<2e-4
      assert u[1]/u[0]-1.<=2*mx+mx*mx+1e-12
      output.append(dict(history=tag,shape_sigma_ln_q=sig,max_reconstructed_scalar_occupation=float(mx),
       max_reconstructed_fermion_occupation=float(cx),relative_vacuum_event_rate_difference=rel,
       vacuum_event_rate_GeV4=e,forward_over_vacuum=float(u[1]/u[0]),inverse_over_quantum_forward=float(v[1]/u[1])))

    return dict(cases=output,anchor_temperature_GeV=.001,
       scope="Anchor pair rates on smooth synthetic reconstructions of the prompt source only. No time evolution. Tiny inverse counts depend on unresolved tails and are not a general bound.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos",type=Path)
    p.add_argument("psd",type=Path)
    p.add_argument("manifest",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.psd,a.manifest),indent=2,allow_nan=False))
