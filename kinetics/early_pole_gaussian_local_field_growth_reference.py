"""Linear hard-mode probe in independent local random-amplitude quartic patches."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc

from early_pole_coherent_quartic_floquet_reference import background,edges


def evolve(spec,n,step,probe,end=240.):
    r=spec['checkpoints'][-1]
    phi=np.array(spec['endpoint_scalar_signed_shell_number_over_T3'])
    p=np.array(spec['scalar_momentum_nodes'])
    z=np.hypot(r['scalar_pole_over_T'],p)
    n0=np.maximum(phi,0.)
    lam=7.65175e-8
    m0=1e-9*r['x']
    rho_ref=r['integrated_scalar_energy_over_T4']
    a,w,_=background(rho_ref,m0,lam)
    m=m0/w
    g=lam*a*a/(6*w*w)
    u=qmc.Sobol(d=2,scramble=True,seed=414497).random_base2(n)
    f,p=(ndtri(u)*[np.sqrt(n0@(1/z))/a,np.sqrt(n0@z)/(a*w)]).T
    e=.5*p*p+.5*m*m*f*f+.25*g*f**4
    localA,localW,k=background(e,m,g*6)
    gaplo,gaphi=edges(k)
    cond=(probe/localW)**2
    band=(cond>gaplo)&(cond<gaphi)
    q1=np.ones_like(f);q2=np.zeros_like(f)
    v1=np.zeros_like(f);v2=np.full_like(f,probe)
    ts=(40.,80.,120.,180.,240.)
    j0=0
    out=[]
    for j in range(round(end/step)+1):
        if j0<len(ts) and j==round(ts[j0]/step):
            norm=((v1*v1+v2*v2)+probe*probe*(q1*q1+q2*q2))/(2*probe*probe)
            ed=.5*p*p+.5*m*m*f*f+.25*g*f**4
            err=np.max(abs(ed/e-1))
            assert err<.012,err
            out.append(dict(time_in_uniform_Omega_units=j*step,largest_relative_oscillator_energy_drift=float(err),
                            energy_weighted_fraction_above_amplitude_factor_10=float(e[norm>100].sum()/e.sum()),
                            energy_weighted_fraction_above_amplitude_factor_100=float(e[norm>1e4].sum()/e.sum()),
                            max_amplitude_factor=float(np.sqrt(max(norm))),
                            share_of_amplified_patch_energy_outside_local_band=float(e[(norm>100)&~band].sum()/max(e[norm>100].sum(),1e-99))))
            j0+=1
        if j==round(end/step):
            break
        a0=-(m*m+g*f*f)*f
        s0=-(probe*probe+m*m+3*g*f*f)
        fnew=f+step*p+.5*step*step*a0
        q1new=q1+step*v1+.5*step*step*s0*q1
        q2new=q2+step*v2+.5*step*step*s0*q2
        a1=-(m*m+g*fnew*fnew)*fnew
        s1=-(probe*probe+m*m+3*g*fnew*fnew)
        p+=.5*step*(a0+a1)
        v1+=.5*step*(s0*q1+s1*q1new)
        v2+=.5*step*(s0*q2+s1*q2new)
        f,q1,q2=fnew,q1new,q2new
    return dict(sample_count=len(f),dimensionless_time_step=step,
                probe_q_over_uniform_Omega=probe,probe_q_over_T=probe*w,
                conditional_local_band_energy_fraction=float(e[band].sum()/e.sum()),
                conditional_local_band_volume_fraction=float(band.mean()),
                records=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    p.add_argument('--refined',action='store_true')
    a=p.parse_args()
    path=Path(a.spec)
    raw=gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
    d=json.loads(raw)
    out=[]
    for probe in (1.272,1.573):
        for n,dt in (((13,.04),(15,.02)) if a.refined else ((11,.04),(13,.02))):
            out.append(evolve(d,n,dt,probe))
    print(json.dumps(out,indent=2,allow_nan=False))
