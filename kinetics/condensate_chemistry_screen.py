"""Conditional 4-zero -> 2-quasiparticle reaction-rate screen.

Vacuum tree amplitude, massive on-shell external/internal lines, replacing
mass by a Hartree gap only as a sensitivity test. Not a controlled
calculation for a coherent condensate. Compare to Arcadi et al. appendix C
(arXiv:1906.07659, eqs 77-78), including the same 4!2! convention.
"""
import argparse,json
from pathlib import Path
from itertools import combinations
import numpy as np
from quartic_consistency import MPL

def amplitude(p,m,lam):
    p=np.asarray(p,float)
    assert p.shape==(6,4)
    # Each unordered 3+3 partition exactly once; include leg 0 on one side.
    den=[]
    for pair in combinations(range(1,6),2):
        v=p[[0,*pair]].sum(axis=0)
        den.append(v[0]**2-v[1:]@v[1:]-m*m)
    return lam*lam*np.sum(1/np.array(den))

def rate(N0,m,lam,T):
    """Total number loss -dN/dt per Tref^4, in an occupied zero mode.
    Forward minus reverse, neglecting the relative 1/N0 finite-volume
    correction to the zero-mode Bose factors. Both zero-mode products N0^4.
    For two gas particles E=2m in saturated gas, f=1/expm1(m/T).
    (1+f)^2-f^2=1+2f, NOT (1+f)^2. Returns vacuum and net loss.
    """
    C=np.sqrt(3)/(12288*np.pi)*lam**4/m**8
    loss=2*C*N0**4
    f=0. if T<=0 or m/T>700 else 1/np.expm1(m/T)
    return loss,loss*(1+2*f),float(f)

def screen(row,limits):
    hist=row['history']
    out=[]
    for case in limits:
        lam=case['lambda_quartic'];states=[]
        for s in case['states']:
            if not s['equilibrium']['saturated']:continue
            x=s['x'];h=min(hist,key=lambda h:abs(h['x']-x))
            b=s['equilibrium'];N=s['initial']['N']
            TH=MPL*x/(1.66*np.sqrt(row['g_bath'])*row['parent_GeV'])*np.sqrt(np.pi**2/30*row['g_bath']/h['total_rho_Tref4'])
            v,net,f=rate(b['N0'],b['m_eff'],lam,b['T'])
            states.append(dict(x=x,number_loss_over_H=net*TH/N,vacuum_loss_over_H=v*TH/N,
                 f_at_2m=f,m_eff=b['m_eff'],condensate_fraction=b['condensate_fraction'],
                 lambda_N0_over_m3=lam*b['N0']/b['m_eff']**3))
        out.append(dict(sigma=case['sigma'],lambda_quartic=lam,states=states))
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--limits',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();r=screen(json.loads(Path(a.input).read_text())[0],json.loads(Path(a.limits).read_text()))
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    for s in r:print('SIGMA',s['sigma'],json.dumps(s['states']))
