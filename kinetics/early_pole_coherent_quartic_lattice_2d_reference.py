"""Conditional periodic classical wave evolution in a two-dimensional plane."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipj, ellipk


def run(N,step,end,seed,stable,seed_rng):
    L=8*2*np.pi/1.272
    ds=L/N
    kx=2*np.pi*np.fft.fftfreq(N,ds)
    ky=2*np.pi*np.fft.rfftfreq(N,ds)
    K=kx[:,None]**2+ky[None,:]**2
    ix=np.fft.fftfreq(N)*N
    jy=np.arange(len(ky))
    active=(abs(ix[:,None])<N//4)&(jy[None,:]<N//4)
    res=(K>1.5)&(K<np.sqrt(3))
    band=(K>.30)&(K<.80) if stable else res
    F=lambda x:np.fft.rfft2(x)
    B=lambda x:np.fft.irfft2(x,s=(N,N),axes=(-2,-1))
    rng=np.random.default_rng(seed_rng)
    phi=np.zeros((N,N))
    coord=np.arange(N)*ds*2*np.pi/L
    # Generate coefficients indexed by the same integer plane waves on
    # both resolutions, not by a resolution-dependent array of noise.
    lo,hi=(.30,.80) if stable else (1.5,np.sqrt(3))
    for iy in range(10):
        for ix1 in range(-9,10):
            if iy==0 and ix1<=0:
                continue
            kk=(2*np.pi/L)**2*(ix1*ix1+iy*iy)
            if lo<kk<hi:
                beta,phase=rng.normal(),rng.uniform(-np.pi,np.pi)
                phi+=beta*np.cos(ix1*coord[:,None]+iy*coord[None,:]+phase)
    phi=1+seed*phi/np.sqrt(np.mean(phi*phi))
    vel=np.zeros_like(phi)
    def force(a):
        return B(-K*F(a)-active*F(a**3))
    a=force(phi)
    initial_seed=F(phi)
    P=2*ellipk(.5)
    j0=0; initial=None; out=[]
    norm=np.broadcast_to(np.where((jy==0)|(jy==N//2),1.,2.),K.shape)/(N**4)
    tot=int(round(end/step))
    for j in range(tot+1):
        if j==round(j0*P/step):
            u,v=F(phi),F(vel)
            kin2=.5*(abs(v)**2+K*abs(u)**2)*norm
            G=float(.5*np.sum(K*abs(u)**2*norm))
            H=kin2.sum()+.25*np.mean(phi**4)
            if initial is None: initial=H
            homogeneous=.5*np.mean(vel)**2+.25*np.mean(phi)**4
            out.append(dict(t_omega=j*step,relative_energy_error=float(H/initial-1),
                            resonant_quadratic_energy_fraction=float(kin2[res].sum()/H),
                            above_point6_quadratic_energy_fraction=float(kin2[K>.36].sum()/H),
                            inactive_uv_fraction=float(kin2[~active].sum()/H),
                            gradient_energy_fraction=G/H,
                            inhomogeneous_energy_fraction=float(1-homogeneous/H)))
            j0+=1
        if j==tot: break
        phi+=step*vel+.5*step*step*a
        aa=force(phi)
        vel+=.5*step*(a+aa)
        a=aa
    t=np.array([q['t_omega'] for q in out])
    iv=np.array([q['inhomogeneous_energy_fraction'] for q in out])
    kr=np.array([q['resonant_quadratic_energy_fraction'] for q in out])
    keep=(t>=35)&(t<=165)
    fit=None if stable else float(.5*np.polyfit(t[keep],np.log(kr[keep]),1)[0])
    predicted=None
    if not stable:
        # Preserve the actual seed weights across the discrete frequencies.
        # A broad seeded band need not grow at its single-mode peak rate.
        ks=K[res]; initial_w=norm[res]*abs(initial_seed[res])**2
        z=len(ks)
        def linear_rhs(t, y):
            u,v=np.split(y,2)
            return np.r_[v,-(ks+3*ellipj(t,.5)[1]**2)*u]
        sol=solve_ivp(linear_rhs,(0,max(t[keep])),np.r_[np.ones(z),np.zeros(z)],
                      t_eval=t[keep],method='DOP853',rtol=1e-10,atol=1e-12)
        u,v=np.split(sol.y,2,axis=0)
        e=initial_w@(v*v+ks[:,None]*u*u)
        predicted=float(.5*np.polyfit(t[keep],np.log(e),1)[0])
        assert abs(predicted-fit)<.0025,(predicted,fit)
    e=max(abs(z['relative_energy_error']) for z in out)
    assert e<.003
    q=np.flatnonzero(iv>.01)
    return dict(cells_per_axis=N,step=step,end=end,box_length=L,seed=seed,
                seed_family="stable" if stable else "resonant",seed_rng=seed_rng,
                seed_band_fourier_count=int(band.sum()),measured_early_amplitude_growth=fit,seeded_linear_pole_prediction=predicted,
                max_abs_energy_drift=e,max_inhomogeneous_fraction=float(iv.max()),
                first_sample_over_1pct_inhomogeneous=None if not len(q) else float(t[q[0]]),
                max_resonant_quadratic_energy_fraction=float(kr.max()),rows=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells',type=int,default=96)
    p.add_argument('--step',type=float,default=.06)
    p.add_argument('--end',type=float,default=340.)
    p.add_argument('--seed',type=float,default=1e-6)
    p.add_argument('--rng',type=int,default=31471)
    p.add_argument('--stable',action='store_true')
    a=p.parse_args()
    print(json.dumps(run(a.cells,a.step,a.end,a.seed,a.stable,a.rng),indent=2,allow_nan=False))
