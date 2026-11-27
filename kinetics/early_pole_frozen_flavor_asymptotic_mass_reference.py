"""Species-dependent lightlike asymptotic mass: frozen kinematic sensitivity."""
import json
import numpy as np
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import intervals,Q
from early_pole_frozen_dynamic_vlasov_reference import load


def read_species(d,name,reconstruct):
    p=np.asarray(d[f'endpoint_{name}_momentum_edges'])
    f=np.asarray(d[f'endpoint_{name}_occupations'])
    v=p**3/3
    u=(v[1:]+v[:-1])/2
    slope=np.zeros_like(f)
    if reconstruct=='limited':
        z=np.diff(f)/np.diff(u)
        slope[1:-1]=np.where(z[:-1]*z[1:]>0.,np.sign(z[:-1])*np.minimum(abs(z[:-1]),abs(z[1:])),0.)
        lim=2*np.minimum(f,1-f)/np.diff(v)
        slope=np.clip(slope,-lim,lim)
    # Independent Gauss integral in p, retaining the exact piecewise
    # shell number of a linear reconstruction in u=p^3/3.
    y,w=roots_legendre(8)
    h=np.diff(p)[:,None]/2
    rad=(p[1:]+p[:-1])[:,None]/2+h*y
    fsol=f[:,None]+slope[:,None]*(rad**3/3-u[:,None])
    assert np.max(abs(np.sum(h*w*rad**2*fsol,axis=1)-f*np.diff(v)))<1e-12
    moment=np.sum(h*w*rad*fsol)
    return p,f,moment,slope,u


def mean(grid,E,mass):
    b,f,mt,slope,u=grid
    p=np.sqrt(np.maximum((E-mass)*(E+mass),0.))
    i=np.searchsorted(b,p,side='right')-1
    i=np.clip(i,0,len(f)-1)
    return np.where((p>=b[0])&(p<=b[-1]),f[i]+slope[i]*(p**3/3-u[i]),0.)


def window(ms,g,mP,mB,p,b):
    out=intervals(mP,ms,mB,Q)
    if out is None:
        return None
    S,lo,hi,_,_=out
    E=[np.hypot(p[0],mP),np.hypot(b[0],mB)+S]
    breaks=np.unique(np.r_[lo,hi,*[e[(e>lo)&(e<hi)] for e in E]])
    y,w=roots_legendre(8)
    dh=np.diff(breaks)[:,None]/2
    mid=(breaks[1:]+breaks[:-1])[:,None]/2+dh*y
    pa,da=mean(p,mid,mP),mean(b,mid-S,mB)
    seed=np.sum(dh*w*pa*(1-da));stim=np.sum(dh*w*(pa-da))
    return dict(g=g,energy_center=(lo+hi)/2,energy_halfspan=(hi-lo)/2,
                seed_bare_aperture=seed,stimulated_bare_aperture=stim,
                interior_breaks=len(dh)-1, daughter_mass_sq=mB*mB,parent_mass_sq=mP*mP)


def run():
    rows=[]
    snaps={grid:load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
           for grid in ('fine','third')}
    poles={grid:d['checkpoints'][-1]['scalar_pole_over_T'] for grid,d in snaps.items()}
    for grid,d in snaps.items():
        for rec in ('constant','limited'):
            p,b=[read_species(d,k,rec) for k in ('parent','fermion')]
            base=np.array([d['x_end'],d['x_end']*d['mass_ratio']])**2
            pref=np.array([p[2],b[2]])/np.pi**2
            D=pref[0]-pref[1]
            for pole_grid,ms in poles.items():
                vac=window(ms,0.,*np.sqrt(base),p,b)
                for g in (.001,.003,.004,.005,.01,.02):
                    v=[]
                    for label,delta in [('thermal_common',np.r_[.25,.25]),
                                        ('frozen_common_mean',np.repeat(1/6+np.mean(pref),2)),
                                        ('frozen_distinct',1/6+pref),
                                        ('cold_photons_distinct',.65**2/6+pref)]:
                        a=window(ms,g,*np.sqrt(base+g*g*delta),p,b)
                        if a is None:
                            v.append(dict(case=label,closed=True))
                            continue
                        a.update(case=label,relative_seed=a['seed_bare_aperture']/vac['seed_bare_aperture'],
                                 relative_stimulus=a['stimulated_bare_aperture']/vac['stimulated_bare_aperture'],
                                 center_shift_in_baseline_halfspans=(a['energy_center']-vac['energy_center'])/vac['energy_halfspan'])
                        v.append(a)
                    assert abs(v[0]['center_shift_in_baseline_halfspans'])<1e-10
                    assert abs(v[2]['center_shift_in_baseline_halfspans']-v[3]['center_shift_in_baseline_halfspans'])<1e-10
                    rows.append(dict(grid=grid,scalar_pole_grid=pole_grid,reconstruction=rec,g=g,
                                     scalar_pole=ms,asymptotic_mass_sq_difference_over_g2=D,
                                     plasma_fermion_mass_terms_per_g2=pref.tolist(),vacuum=vac,shifts=v))
    return dict(scope='Lightlike massless-loop estimate inserted only as an asymptotic dispersion proxy. '
                      'No matched hard massive self energies, residues, rates, vertex, bath shift, or scalar evolution.',
                rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
