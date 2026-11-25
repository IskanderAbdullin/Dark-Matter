"""Fit the *physical decay width* at empty-parent initial conditions.

No rescaling of a frozen PSD. Store every evaluated kinetic trajectory,
including the two independent fine-mesh checks. The model assumes kinetic
equilibrium for the mother; this is NOT full SM momentum-dependent production.
"""
import argparse,json
from pathlib import Path
from scipy.optimize import brentq
from yukawa_rate_bath import YukawaRateBath,rate_mapping


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--kappa-sm1',type=float,default=47.)
    p.add_argument('--kappa2',type=float,default=.03)
    p.add_argument('--fermion-ratio',type=float,default=.9)
    p.add_argument('--parent-gev',type=float,default=86000.)
    p.add_argument('--x0',type=float,default=.01)
    p.add_argument('--xend',type=float,default=80.)
    p.add_argument('--n',type=int,default=160)
    p.add_argument('--fine-n',type=int,nargs='+',default=[240,320])
    p.add_argument('--qmin',type=float,default=1e-10)
    p.add_argument('--bracket',nargs=2,type=float,default=[.00224,.008])
    p.add_argument('--target',type=float,default=.12)
    p.add_argument('--rtol',type=float,default=1e-6)
    p.add_argument('--width-xtol',type=float,default=3e-5)
    p.add_argument('--output',required=True)
    a=p.parse_args();out=[];memo={}
    def evaluate(k,n=None):
        n=a.n if n is None else n
        key=(k,n)
        if key in memo:return memo[key]
        m=YukawaRateBath(n=n,scalar_ratio=86./(a.parent_gev*1e9),fermion_ratio=a.fermion_ratio,
            kappa=k,kappa2=a.kappa2,kappa_sm1=a.kappa_sm1,gamma=.01,scatter=True,qmin=a.qmin)
        r=m.run(a.x0,a.xend,a.rtol,1e-15,True,a.parent_gev,'empty')
        r['mapping']=rate_mapping(a.parent_gev,a.fermion_ratio,k,a.kappa2,a.kappa_sm1)
        out.append(r);Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        s=r['final']
        print('tune',json.dumps(dict(n=n,kappa=k,omega=s['omega_phi_h2_m86'],
            hot=s['hot_fraction_p_over_actual_T_gt_0p1'],peak=s['peak_bin_q'],
            continuity=s['continuity_residual'])),flush=True)
        memo[key]=s['omega_phi_h2_m86']-a.target
        return memo[key]
    lo,hi=a.bracket;fa=evaluate(lo);fb=evaluate(hi)
    if fa*fb>=0:raise ValueError(('No sign-changing abundance bracket',fa,fb))
    root=brentq(evaluate,lo,hi,xtol=a.width_xtol,rtol=1e-4)
    print('root',root,flush=True)
    for n in a.fine_n:evaluate(root,n)

if __name__=='__main__':main()
