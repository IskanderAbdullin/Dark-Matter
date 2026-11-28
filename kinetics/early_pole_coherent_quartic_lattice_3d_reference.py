"""Small periodic 3D coherent quartic smoke test; not a converged spectrum."""
import argparse
import json

import numpy as np
from scipy.special import ellipk


def run(N,dt,end,seed,stable=False):
    L=4*2*np.pi/1.272
    x=np.arange(N)*(2*np.pi/N)
    k=2*np.pi*np.fft.fftfreq(N,L/N)
    z=k[:N//2+1].copy()
    z[-1]=-z[-1]
    q2=k[:,None,None]**2+k[None,:,None]**2+z[None,None,:]**2
    active=(abs(k[:,None,None])<abs(k[N//4]))&(abs(k[None,:,None])<abs(k[N//4]))&(abs(z[None,None,:])<abs(k[N//4]))
    band=(q2>1.5)&(q2<np.sqrt(3))
    W=np.where((np.arange(len(z))==0)|(np.arange(len(z))==N//2),1.,2.)/N**6
    F=lambda p:np.fft.rfftn(p)
    B=lambda p:np.fft.irfftn(p,s=(N,N,N),axes=(0,1,2))
    mode=2 if stable else 4
    f=1+seed*(np.cos(mode*x[:,None,None])+.7*np.cos(mode*x[None,:,None]+.3)+.9*np.cos(mode*x[None,None,:]+.8))
    v=np.zeros_like(f)
    def kick(a): return B(-q2*F(a)-active*F(a**3))
    a=kick(f)
    j0=0; P=2*ellipk(.5); t=[]; b=[]; e=[]; hom=[]
    for j in range(round(end/dt)+1):
        if j==round(j0*P/dt):
            u,p=F(f),F(v)
            E=.5*W*(abs(p)**2+q2*abs(u)**2)
            total=E.sum()+.25*np.mean(f**4)
            t.append(j*dt);e.append(total);b.append(E[band].sum()/total)
            hom.append(1-(.5*np.mean(v)**2+.25*np.mean(f)**4)/total)
            j0+=1
        if j==round(end/dt):break
        f+=dt*v+.5*dt*dt*a
        aa=kick(f);v+=.5*dt*(a+aa);a=aa
    t,b,e,hom=map(np.asarray,(t,b,e,hom))
    use=(t>=35)&(t<=min(150.,end*.85))
    mu=None if stable else .5*np.polyfit(t[use],np.log(b[use]),1)[0]
    err=np.max(abs(e/e[0]-1))
    assert err<.002
    if stable: assert np.max(hom)<1e-9
    else: assert .034<mu<.038
    return dict(cells_per_axis=N,time_step=dt,requested_end=end,initial_seed_is_off_band=stable,initial_seed_q_over_Omega=float(k[mode]),
                linear_band_amplitude_exponent=None if mu is None else float(mu),
                max_energy_relative_drift=float(err),time=t.tolist(),band_energy_fraction=b.tolist(),
                inhomogeneous_energy_fraction=hom.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells',type=int,default=32)
    p.add_argument('--step',type=float,default=.06)
    p.add_argument('--end',type=float,default=260.)
    p.add_argument('--seed',type=float,default=1e-6)
    p.add_argument('--stable',action='store_true')
    a=p.parse_args()
    print(json.dumps(run(a.cells,a.step,a.end,a.seed,a.stable),indent=2,allow_nan=False))
