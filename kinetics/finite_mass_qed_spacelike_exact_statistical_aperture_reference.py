"""Free on-shell finite-recoil statistical ratio for a spacelike transverse vertex."""
import gzip,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import expit,roots_legendre
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Species
from finite_mass_qed_transverse_kinetic_noise_reference import archive


def weight(m,w,k):
    q2=(k-w)*(k+w)
    Emin=(k*np.sqrt(1+4*m*m/q2)-w)/2
    def M(E):
        pz=(w*E-q2/2)/k
        return 4*((E-m)*(E+m)-pz*pz)+2*q2
    return Emin,M


def cell(m,k,w,p,f):
    lo,M=weight(m,w,k)
    edge=np.hypot(m,p)
    hi=edge[-1]
    if hi<=lo:return np.zeros(3)
    a=np.unique(np.clip(np.r_[lo,hi,edge,edge-w],lo,hi))
    x,t=roots_legendre(3)
    E=(a[:-1,None]+a[1:,None])/2+(a[1:,None]-a[:-1,None])*x/2
    y=(a[:-1]+a[1:])/2
    f=np.r_[0.,f,0.]
    p,b=f[np.searchsorted(edge,y)],f[np.searchsorted(edge,y+w)]
    w0=np.diff(a)/2*(M(E)@t)
    return np.array([w0@(p*(1-b)),w0@(b*(1-p)),w0@(p-b)])


def run():
    rows=[]
    thermal=[]
    fluctuation=[]
    for m in (.003,.05):
        for k in (.0002,.00256,.015):
            for u in (.0047,.25,.9):
                w=k*u
                lo,M=weight(m,w,k)
                for T in (.4,1.):
                    G=quad(lambda E:M(E)*expit(-E/T)*expit((E+w)/T),lo,lo+60*T)[0]
                    L=quad(lambda E:M(E)*expit(-(E+w)/T)*expit(E/T),lo,lo+60*T)[0]
                    thermal.append(abs((L/G)/np.exp(-w/T)-1.))
                    D=quad(lambda E:M(E)*(expit(-E/T)-expit(-(E+w)/T)),lo,lo+60*T,
                           epsabs=1e-15,epsrel=5e-12)[0]
                    ratio=w*(G+L)/(2*D)
                    coth=w/(2*np.tanh(w/(2*T)))
                    fluctuation.append(abs(ratio/coth-1.))
    if max(thermal)>3e-13 or max(fluctuation)>8e-10:raise ArithmeticError('Exact statistical equilibrium gate')
    for g in ('fine','third'):
        root=Path(__file__).parent/'results'
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        M=d['x_end'];m=d['mass_ratio']*M
        ch=[Species(q,d[f'endpoint_{a}_momentum_edges'],d[f'endpoint_{a}_occupations']) for q,a in ((M,'parent'),(m,'fermion'))]
        ss=[]
        for k in (.0002,.0008,.00256,.01):
            for u in (.0047,.25,.9,.99):
                w=k*u
                G,L,D=np.sum([cell(s.m,k,w,s.e,s.f) for s in ch],axis=0)
                A,B=np.sum([archive(s,u,36)[:2] for s in ch],axis=0)
                if min(G,L,D)<=0. or abs((G-L)-D)>2e-13*max(G,L):
                    raise ArithmeticError('Exact archived raw sign or Pauli reconstruction')
                x=w*(G+L)/(2*D)
                ss.append(dict(k=k,u=u,exact_symmetrized_over_dissipation_energy=x,soft_limit=A/B,
                               relative_shift=x/(A/B)-1.,exact_cut_positive=D))
        rows.append(dict(grid=g,rows=ss))
    return dict(scope='Free on-shell charged cut statistical ratio; not interacting photon population or matched scalar rate.',thermal_KMS_check=max(thermal),
                thermal_coth_check=max(fluctuation),data=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
