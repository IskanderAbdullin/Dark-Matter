"""Exact 1↔2 angle for the reciprocal hard quasiparticle energy gradient."""
import json
import numpy as np
from scipy.optimize import brentq

from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_snapshot_gauge_ledger_reference import cell_J, photon
from early_pole_decay_mode_clock_reference import Q,intervals


def evalcase(g,cut,s):
    x=s['x_end']; M=x
    B=s['mass_ratio']*x
    S=s['checkpoints'][-1]['scalar_pole_over_T']
    Jg=photon(1.,0.,cut)
    P,D=[cell_J(s[f'endpoint_{n}_momentum_edges'],s[f'endpoint_{n}_occupations'],m,cut)
         for n,m in [('parent',M),('fermion',B)]]
    a,b=.5*g*g/np.pi**2*np.array([Jg+P,Jg+D])
    Ws=np.hypot(S,Q)
    def e(v,m,c): 
        v0=np.hypot(v,m)
        return v0+c/v0
    p1=brentq(lambda p:e(p,M,a)-Ws-e(p+Q,B,b),.2,4.,xtol=5e-14)
    p2=brentq(lambda p:e(p,M,a)-Ws-e(abs(p-Q),B,b),.2,4.,xtol=5e-14)
    assert p2>p1
    p=np.linspace(p1,p2,83)[1:-1]
    target=e(p,M,a)-Ws
    dot=(target+np.sqrt(target*target-4*b))/2
    k=np.sqrt((dot-B)*(dot+B))
    mu=(p*p+Q*Q-k*k)/(2*p*Q)
    assert max(abs(mu))<1.
    v=k/dot*(1-b/dot**2)
    jac=k/(p*Q*v)
    def mu_energy(c):
        return e(np.sqrt(p*p+Q*Q-2*p*Q*c),B,b)
    h=4e-4
    fd=-(mu_energy(mu+h)-mu_energy(mu-h))/(2*h)
    check=max(abs(fd*jac-1))
    assert check<2e-6
    # Compare the same squared-mass increment used as an exact square root
    # with the actual first-order variational derivative.
    _,lo,hi,_,_=intervals(np.sqrt(M*M+2*a),S,np.sqrt(B*B+2*b),Q)
    pp=np.sqrt((np.array([lo,hi])-np.sqrt(M*M+2*a))*(np.array([lo,hi])+np.sqrt(M*M+2*a)))
    surrogate=e(p,M,a)-Ws-e(k,B,b)
    assert max(abs(surrogate))<2e-14
    return dict(g=g,cut=cut,parent_p_aperture=[p1,p2],sqrt_mass_shift_aperture=pp.tolist(),
                endpoint_delta_in_halfspan=((np.array([p1,p2])-pp)/(p2-p1)*2).tolist(),
                angular_energy_jacobian_max_fd_error=check,
                word='Angle/momentum only; no spin residue, mean-field feedback, or full transfer spectrum.')


def run():
    s=load('results/early-pole-spectrum-third-m06-2026-11-26.json.gz')
    return dict(scope='Hard isotropic perturbative dispersion; exact event geometry within that model.',
                rows=[evalcase(g,cut,s) for g in (0.,.005,.01,.05) for cut in (0.,.5)])


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
