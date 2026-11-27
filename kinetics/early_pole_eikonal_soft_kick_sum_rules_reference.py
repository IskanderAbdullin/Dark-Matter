"""Eikonal HTL test-particle momentum-kick coefficients, not a two-line rate."""
import argparse
import json
import numpy as np
from scipy.special import roots_legendre


def log_nodes(low,high,n,segments=40):
    a,w=roots_legendre(n)
    e=np.linspace(np.log(low),np.log(high),segments+1)
    h=np.diff(e)[:,None]/2
    y=np.exp((e[1:]+e[:-1])[:,None]/2+h*a)
    return y.ravel(),(h*w*y).ravel()


def kernel(c,n,zcut=None):
    pt,tw=log_nodes(2e-5,c,n)
    pz,zw=log_nodes(2e-15,1e4,n,90)
    pt,pz=pt[:,None],pz[None,:]
    k=np.hypot(pt,pz)
    x=pz/k
    c2=(pt/k)**2
    h=np.arcsinh(pz/pt)
    L=1/((k*k+1-x*h)**2+(np.pi*x/2)**2)
    T=c2*c2/2/((pt*pt+.5*(x*x+x*c2*h))**2+(np.pi*x*c2/4)**2)
    W=(2*pt/k)*(L+T)*tw[:,None]*zw[None,:]
    if zcut is not None:
        W=W*(k<=zcut)
    return W,pt,pz


def run(n):
    rows=[]
    for cut in (2.,5.,10.,20.,40.):
        W,pt,pz=kernel(cut,n)
        long=float(np.sum(W*pz*pz))
        perp=float(np.sum(W*pt*pt))
        la=.5*np.log1p(2*cut*cut)
        ta=np.log1p(cut*cut)
        assert max(abs(long/la-1),abs(perp/ta-1))<2e-4
        rows.append(dict(qperp_cutoff_over_mD=cut,kick_qhat_L=long,kick_qhat_T=perp,
                         longitudinal_sum_rule=la,transverse_sum_rule=ta))
    return dict(scope='Classical equilibrium lightlike Wilson test-particle; transverse cutoff; '
                      'no hard matching, quantum vertex, or cosmological rate.',rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nodes',type=int,default=8)
    print(json.dumps(run(p.parse_args().nodes),indent=2,allow_nan=False))
