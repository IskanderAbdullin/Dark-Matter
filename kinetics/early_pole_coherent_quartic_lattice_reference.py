"""A periodic one-dimensional coherent quartic field, not a cosmological spectrum."""
import argparse
import json

import numpy as np
from scipy.special import ellipk


def run(cells,dt,tend,seed,modes,stable=False,modulation=0.):
    peak=1.272
    L=modes*2*np.pi/peak
    k=2*np.pi*np.fft.rfftfreq(cells,L/cells)
    z=np.arange(cells)*(L/cells)
    active=np.arange(len(k))<cells//4
    def spectral(f):
        return np.fft.rfft(f)
    def fx(phi):
        u=spectral(phi)
        v=spectral(phi**3)
        v[~active]=0.
        return np.fft.irfft(-k*k*u-v,n=cells)
    kt=.75 if stable else peak
    kindex=int(np.argmin(abs(k-kt)))
    phi=1.+seed*np.cos(k[kindex]*z)+modulation*np.cos(k[1]*z)
    vel=np.zeros(cells)
    a=fx(phi)
    P=2*ellipk(.5)
    first=0.
    tsave=0
    out=[]
    e0=None
    for j in range(int(round(tend/dt))+1):
        if j==round(tsave*P/dt):
            u=spectral(phi)
            pp=spectral(vel)
            der=np.fft.irfft(1j*k*u,n=cells)
            energy=.5*np.mean(vel**2+der**2)+.25*np.mean(phi**4)
            if e0 is None: e0=energy
            mean=.5*np.mean(vel)**2+.25*np.mean(phi)**4
            mode=2*np.hypot(abs(u[kindex]),abs(pp[kindex])/2)/cells
            grad=.5*np.mean(der**2)
            bin_energy=.5*np.mean((phi-phi.mean())**2)
            # Orthogonal quadratic spectral diagnostics, not a split of
            # the full nonlinear potential energy into particle modes.
            scale=np.ones_like(k)*2
            scale[0]=scale[-1]=1
            spec_energy=scale*(abs(pp)**2+k*k*abs(u)**2)/(2*cells*cells)
            band=(k*k>1.5)&(k*k<np.sqrt(3))
            phard=float(spec_energy[k>=.6].sum()/energy)
            pband=float(spec_energy[band].sum()/energy)
            out.append(dict(t_in_Omega_units=j*dt,mode_norm=float(mode),
                            energy_relative_error=float(energy/e0-1),
                            inhomogeneous_energy_fraction=float(1-mean/energy),
                            gradient_energy_fraction=float(grad/energy),
                            fluctuation_squared_over_two=float(bin_energy),
                            quadratic_resonance_band_energy_fraction=pband,
                            quadratic_k_ge_point6_energy_fraction=phard))
            tsave+=1
        if j>=round(tend/dt): break
        phi=phi+dt*vel+.5*dt*dt*a
        ap=fx(phi)
        vel=vel+.5*dt*(a+ap)
        a=ap
    t=np.array([r['t_in_Omega_units'] for r in out])
    amp=np.array([r['mode_norm'] for r in out])
    if modulation==0 and not stable:
        keep=(t>=30)&(t<=160)
        growth=float(np.polyfit(t[keep],np.log(amp[keep]),1)[0])
        assert .034<growth<.0385
    else:
        growth=None
    dist=np.array([r['inhomogeneous_energy_fraction'] for r in out])
    flag=np.flatnonzero(dist>.01)
    return dict(cells=cells,box_length=L,dt=dt,requested_time=tend,seed=seed,mode_q_over_Omega=float(k[kindex]),
                box_slowest_q_over_Omega=float(k[1]),initial_low_mode_amplitude=modulation,
                pre_backreaction_log_amplitude_slope=growth,maximum_absolute_energy_drift=float(max(abs(r['energy_relative_error']) for r in out)),
                first_sample_over_one_percent_inhomogeneous_energy=None if not len(flag) else float(t[flag[0]]),
                max_inhomogeneous_energy_fraction=float(dist.max()),
                final_mode_norm_over_initial=float(amp[-1]/amp[0]),rows=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells',type=int,default=512)
    p.add_argument('--step',type=float,default=.06)
    p.add_argument('--end',type=float,default=360.)
    p.add_argument('--seed',type=float,default=1e-6)
    p.add_argument('--mode-count',type=int,default=16)
    p.add_argument('--stable',action='store_true')
    p.add_argument('--modulation',type=float,default=0.)
    a=p.parse_args()
    result=run(a.cells,a.step,a.end,a.seed,a.mode_count,a.stable,a.modulation)
    assert result['maximum_absolute_energy_drift']<3e-3
    print(json.dumps(result,indent=2,allow_nan=False))
