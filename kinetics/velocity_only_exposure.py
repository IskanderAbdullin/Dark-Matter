"""Cross-section-independent dilute collision exposure for a fixed v-dependent SIDM.

Assume frozen isotropic comoving scalar PSD and a cold subset p<cut*p_peak.
Use exact relative-angle density of two independent radial momenta.
The cold-cold pair density is UNNORMALIZED by subset, i.e. normalized by
the full number squared. Thus the cosmological event count per total particle
is a conservative positive contribution. This is NOT a net relaxation rate.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import trapezoid
from quartic_consistency import shell_nodes
from late_activation_overlap import H,elapsed

C=2.99792458e10
RHO0=.12*1.8783416e-29
SEC_GYR=1e9*365.25*86400

def pair_intervals(q,w):
    i,j=np.triu_indices(len(q))
    q1,q2=q[i],q[j]
    c=w[i]*w[j]*np.where(i==j,1.,2.)/(2*q1*q2)
    return abs(q1-q2),q1+q2,c

def pair_coefficient(lo,hi,c,x):
    e=np.r_[lo,hi];d=np.r_[c,-c]
    ix=np.argsort(e,kind='stable');e=e[ix];d=d[ix]
    v=np.r_[0.,np.cumsum(d)]
    a=v[np.searchsorted(e,x,side='right')]
    return np.maximum(a,0.)

def halo_pdf(v,vT,escape=None):
    # For individual f(v) proportional exp(-v^2/vT^2),
    # relative speed distribution is Maxwell with 1D std = vT.
    # Use untruncated Maxwell, conservative upper comparison near tail
    # is not guaranteed for all v; no observational fit is implied.
    return np.sqrt(2/np.pi)*v**2/vT**3*np.exp(-v*v/(2*vT*vT))

def calculate(row,cut=6.,order=8,nw=12000,nv=340,zmax=1e7):
    h=row['final'];p0=h['peak_bin_q']
    e=np.asarray(row['frozen_cell_spectrum']['edges']);f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    e2=np.sort(np.unique(np.r_[e,cut*p0]));e2=e2[(e2>=e[0])&(e2<=e[-1])]
    f2=f[np.clip(np.searchsorted(e,(e2[1:]+e2[:-1])/2)-1,0,len(f)-1)]
    q,w=shell_nodes(e2,f2,order);N=w.sum()
    keep=(q<cut*p0)&(w>0)
    q,w=q[keep]/p0,w[keep]/N
    Tref0=2.7255*8.617333262e-5*(3.909/row['g_bath'])**(1/3)/h['theta']
    v0=p0*Tref0/86*C
    # Require each particle, not only the relative velocity, to be safely NR.
    rmax=min(1+zmax,.03*C/(cut*v0))
    a,b,c=pair_intervals(q,w)
    x=np.geomspace(1e-9,2*cut,nw)
    coef=pair_coefficient(a,b,c,x)
    P=x*coef
    # Diagonal radial nodes cause P~x at tiny x rather than the continuum x^2.
    # Its contribution to the integrated exposure is explicitly tested in order.
    pair_mass=float(trapezoid(P,x))
    # Cosmological exposure kernel in relative speed:
    # d tau/dv = rho0 * v * integral_r r/H * P0(v/r) dr,
    # with P0(w)=P(w/v0)/v0. Change variable x=v/(r v0).
    # Divide out v, which cancels against halo kernel.
    vel=np.geomspace(1e2,3e8,nv)
    z_at_v=vel[:,None]/(v0*x[None,:])-1
    valid=(z_at_v>=0)&(z_at_v<=rmax-1)
    z=np.clip(z_at_v,0,rmax-1)
    kernel=RHO0*(vel[:,None]/v0)**2*P[None,:]/(v0*x[None,:]**3*H(z))
    J=trapezoid(np.where(valid,kernel,0.),x,axis=1)
    rows=[]
    for rho,vT in [(1e-24,100e5),(6.700740771077128e-24,np.sqrt(2)*44.1e5),(1e-23,30e5),(1e-22,10e5)]:
        denom=rho*10*SEC_GYR*halo_pdf(vel,vT)
        ratio=J/np.maximum(denom,1e-300)
        # Give explicit supported velocity intervals, not a universal no-go.
        windows=[]
        for vlo,vhi in [(1e3,3*vT),(0.1*vT,3*vT),(.01*vT,6*vT),(1e3,6*vT)]:
            sel=(vel>=vlo)&(vel<=vhi)
            i=np.flatnonzero(sel)[np.argmin(ratio[sel])]
            windows.append(dict(vmin_kms=vlo/1e5,vmax_kms=vhi/1e5,min_exposure_ratio=float(ratio[i]),argmin_kms=float(vel[i]/1e5)))
        rows.append(dict(rho=rho,vT_kms=vT/1e5,halo_age_Gyr=10,windows=windows,ratio=ratio.tolist()))
    return dict(cut_over_peak=cut,order=order,nw=nw,nv=nv,N=N,cold_number_fraction=float(w.sum()),
        pair_probability_expected=float(w.sum()**2),pair_probability_integral=pair_mass,peak_velocity_today_cm_s=v0,
        zmax=float(rmax-1),max_individual_speed_over_c=float(cut*v0*rmax/C),
        velocity_kms=(vel/1e5).tolist(),cosmic_kernel=J.tolist(),halo=rows,
        scope='Cold-cold positive dilute collision count on frozen PSD only. No Bose factors, no net relaxation, no claim of cosmological exclusion.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--cut',type=float,default=6.);p.add_argument('--order',type=int,default=8);p.add_argument('--nw',type=int,default=12000);p.add_argument('--nv',type=int,default=340);a=p.parse_args()
    r=calculate(json.loads(Path(a.input).read_text())[0],a.cut,a.order,a.nw,a.nv)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print('fraction',r['cold_number_fraction'],'pair',r['pair_probability_integral'],r['pair_probability_expected'],'zmax',r['zmax'])
    for h in r['halo']:print(h['rho'],h['vT_kms'],h['windows'])
