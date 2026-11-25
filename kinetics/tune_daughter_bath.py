"""Bracket and tune the physical channel width, never rescale the PSD."""
import argparse
import json
from pathlib import Path
from scipy.optimize import brentq
from daughter_bath_relaxation import DaughterBath


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=120)
    p.add_argument('--fermion-ratio',type=float,required=True);p.add_argument('--kappa2',type=float,default=.03)
    p.add_argument('--bracket',nargs=2,type=float,required=True)
    p.add_argument('--x0',type=float,default=.1);p.add_argument('--xend',type=float,default=80.)
    p.add_argument('--target',type=float,default=.12);p.add_argument('--rtol',type=float,default=2e-6)
    p.add_argument('--qmin',type=float,default=1e-7);p.add_argument('--scalar-ratio',type=float,default=1e-9)
    p.add_argument('--parent-gev',type=float,default=86.);p.add_argument('--output',required=True);a=p.parse_args()
    out=[];cache={}
    def evaluate(k):
        if k not in cache:
            m=DaughterBath(n=a.n,scalar_ratio=a.scalar_ratio,fermion_ratio=a.fermion_ratio,kappa=k,
                           kappa2=a.kappa2,relaxation=100.,qmin=a.qmin)
            r=m.evolve_bath(a.x0,a.xend,a.rtol,1e-15,True,a.parent_gev);out.append(r);cache[k]=r
            s=r['final'];print(json.dumps(dict(n=a.n,kappa=k,kappa2=a.kappa2,b=a.fermion_ratio,
                     omega=s['omega_phi_h2_m86'],hot=s['hot_fraction_p_over_actual_T_gt_0p1'],
                     mean=s['mean_p_over_actual_T'],residual=s['continuity_residual'],
                     fermion_number=s['fermion_number'],matter_fraction_max=max(t['matter_fraction_chi2'] for t in r['history']))),flush=True)
            Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        return cache[k]['final']['omega_phi_h2_m86']-a.target
    lo,hi=a.bracket;flo,fhi=evaluate(lo),evaluate(hi)
    if flo*fhi>0:raise ValueError('No abundance bracket: '+str((lo,flo,hi,fhi)))
    k=brentq(evaluate,lo,hi,xtol=1e-8,rtol=2e-5,maxiter=30)
    err=evaluate(k)
    print('matched',k,err,flush=True)
    assert abs(err)<2e-4,err
    Path(a.output+'.matched.json').write_text(json.dumps([cache[k]],indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
