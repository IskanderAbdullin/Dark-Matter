"""Independent parent-first integral of the exact signed stimulus kernel."""
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre
from finite_mass_qed_fixed_q_signed_soft_reference import F
from finite_mass_qed_parent_integrated_signed_soft_reference import Source

def independent(S,k):
    M,m,s,Es,q=S.M,S.m,S.s,S.Es,S.q
    roots=[]
    for t in np.r_[0.,np.geomspace(1e-8,1.,20),2.]:
        def fun(ep,sign):
            p=np.sqrt((ep-M)*(ep+M))
            ell=M*M/(ep+p)+p*t
            return (ep-k)*Es-(M*M+s*s-m*m)/2+k*ell-sign*q*np.sqrt((p-k)**2+2*p*k*t)
        roots.extend([brentq(lambda ep:fun(ep,e),float(S.lo)-.12,float(S.hi)+.12,xtol=2e-16,rtol=1e-15) for e in (-1,1)])
    l,h=min(roots),max(roots)
    edges=np.unique(np.r_[roots,S.pe[(S.pe>l)&(S.pe<h)],(S.be+Es+k)[(S.be+Es+k>l)&(S.be+Es+k<h)]])
    x,w=roots_legendre(12)
    ans=0.
    for a,b in zip(edges[:-1],edges[1:],strict=True):
        E=(a+b)/2+(b-a)*x/2
        P=S.pf[np.searchsorted(S.pe,E)-1]
        D=S.bf[np.searchsorted(S.be,E-Es-k)-1]
        val=np.asarray([F(M,m,s,u,q,k,64,24)[0] for u in E])
        ans+=(b-a)/2*(w@((P-D)*val))
    return ans

def run():
    out=[]
    for g in ('fine','third'):
        p=Path(__file__).parent/f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz'
        with gzip.open(p) as f:
            S=Source(json.load(f))
        for k in (-1e-5,-3e-6,3e-6,1e-5):
            a=S.signed(k,96,18,36)[0,2]
            b=independent(S,k)
            delta=float(abs(a-b))
            if delta>3e-15:
                raise ValueError('The exact parent-first and photon-first integrals disagree.')
            out.append(dict(grid=g,signed_photon_energy=k,photon_first=float(a),parent_first=float(b),
                            absolute_defect=delta))
    return out


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
