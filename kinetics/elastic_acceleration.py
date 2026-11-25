"""Optional compiled evaluation of the existing conservative event graph.

No changed kernel, reduced precision, parallel atomic summation, fast-math,
or physics approximation. Sequential accumulation into compact RHS and J.
The NumPy reference stays available for a numerical comparison.
"""
import sys
from pathlib import Path
DEPS=Path(__file__).resolve().parents[1]/'.runtime'/'elastic-deps'
if DEPS.exists():sys.path.insert(0,str(DEPS))
import numpy as np
from numba import njit

@njit(cache=True,fastmath=False)
def compiled_rates(f,idx,c,stoich,weight,inv,fl,need_jac):
    n=len(f);ne=len(weight)
    ff=np.maximum(f,fl)
    logf=np.log(ff);r=np.log1p(1/ff)
    out=np.zeros(n);J=np.zeros((n,n));D=np.empty(5)
    for t in range(ne):
        i,j,k,l,h=idx[0,t],idx[1,t],idx[2,t],idx[3,t],idx[4,t]
        w=c[4,t];wl=c[3,t]
        prod=logf[i]+logf[j]+logf[k]+wl*logf[l]+w*logf[h]
        rin=r[i]+r[j];rout=r[k]+wl*r[l]+w*r[h]
        delta=rout-rin
        F=np.sign(delta)*weight[t]*np.exp(prod+max(rin,rout))*(-np.expm1(-abs(delta)))
        for a in range(5):
            out[idx[a,t]]+=stoich[a,t]*F
        if need_jac:
            fw=weight[t]*np.exp(prod+rout);bw=weight[t]*np.exp(prod+rin)
            for b in range(5):
                s=idx[b,t];v=ff[s];cb=c[b,t]
                if f[s]<=fl:
                    D[b]=0.
                elif b<2:
                    if v<1e-3:D[b]=cb*(fw/v-bw/(1+v))
                    else:D[b]=cb*(F/v+bw/(v*(1+v)))
                else:
                    if v<1e-3:D[b]=cb*(fw/(1+v)-bw/v)
                    else:D[b]=cb*(F/v-fw/(v*(1+v)))
            for a in range(5):
                for b in range(5):
                    J[idx[a,t],idx[b,t]]+=stoich[a,t]*D[b]
    for a in range(n):
        out[a]*=inv[a]
        if need_jac:
            for b in range(n):J[a,b]*=inv[a]
    return out,J
