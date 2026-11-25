"""Order-of-magnitude thermal consistency screen; NOT a finite-T collision calculation.

Eq. (18) of Perez-Gonzalez, Sen & Tangarife (2026), arXiv:2606.14441,
gives equilibrium simple-Yukawa delta m_phi^2=g^2 T^2/6 and
delta m_f^2=g^2 T^2/16. These are benchmark coefficients, not
self-energies for a derivative ALP portal or the nonthermal state.
We also compute the vacuum-subtracted scalar variance from the frozen
PSD, as a proxy for the possible occupation-induced medium scale.
No equilibrium Bose enhancement multiplier is applied to our collision
operator, which already includes stimulated emission and inverse decay.
"""
import argparse,json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

from yukawa_rate_bath import rate_mapping


def backward_energy(parent,daughter,scalar,p):
    """Minimum lab energy for fixed parent momentum, stable boost formula."""
    assert parent>daughter+scalar and p>=0
    e=(parent*parent-daughter*daughter+scalar*scalar)/(2*parent)
    ps=np.sqrt(max(0.,(e-scalar)*(e+scalar)))
    ep=np.hypot(parent,p)
    return float((parent*parent*e/(ep+p)+p*scalar*scalar/(e+ps))/parent)


def variance_proxy(spectrum,mass_over_Tref,theta,nq=16):
    # Vacuum-subtracted <phi^2>_particles / T_bath^2 for one real scalar.
    e=np.asarray(spectrum['edges']);f=np.asarray(spectrum['f_phi'])
    assert f.min()>-1e-7
    l,r=e[:-1],e[1:]
    u,w=leggauss(nq)
    q=l[:,None]+(r-l)[:,None]*(u[None,:]+1)/2
    vals=q*q/np.hypot(q,mass_over_Tref)
    var=np.sum(np.maximum(f,0)*(r-l)*.5*(vals@w))/(2*np.pi*np.pi*theta*theta)
    return float(var)


def screen(row):
    M=row['parent_GeV'];r=row['fermion_ratio'];mu=row['scalar_ratio']
    mp=rate_mapping(M,r,row['kappa'],row['kappa2'],row['kappa_sm1'],row['g_bath'])
    g=mp['pseudoscalar_coupling'];y1=mp['y1'];y2=mp['y2'];eta=1-r*r
    massT=g/np.sqrt(6.)
    states=[]
    last=row['final']['number']
    for h in row['history']:
        x=h['x'];theta=h['theta'];xb=x/theta
        if x>1 or h['number']<.01*last:continue
        peak=h['peak_bin_q']/theta
        m1=np.hypot(xb,g/4);m2=np.hypot(r*xb,g/4);mphi=np.hypot(mu*xb,massT)
        vac=backward_energy(xb,r*xb,mu*xb,3.)
        allowed=bool(m1>m2+mphi)
        corrected=backward_energy(m1,m2,mphi,3.) if allowed else None
        var=variance_proxy(row['frozen_cell_spectrum'],mu*x,theta)
        var2=variance_proxy(row['frozen_cell_spectrum'],mu*x,theta,32)
        assert abs(var/var2-1)<1e-7
        states.append(dict(x=x,x_actual=xb,produced_fraction=h['number']/last,peak_p_over_T=peak,
            scalar_thermal_mass_over_peak=float(massT/peak),scalar_thermal_mass_over_backward_energy=float(massT/vac),
            representative_backward_energy_vac_over_T=vac,representative_backward_energy_thermal_over_T=corrected,
            equilibrium_threshold_open=allowed,portal_mass_splitting_proxy=float((y1*y1+y2*y2)/(eta*xb*xb)),
            frozen_PSD_variance_proxy_over_T2=var,
            occupation_induced_mass2_over_vacuum_gap_proxy=float(g*g*var/(eta*xb*xb))))
    return dict(initial_parent=row.get('initial_parent','thermal'),x0=row['x0'],n=row['n'],kappa=row['kappa'],
        kappa_sm1=row['kappa_sm1'],omega=row['final']['omega_phi_h2_m86'],hot=row['final']['hot_fraction_p_over_actual_T_gt_0p1'],
        mapping=mp,scalar_equilibrium_thermal_mass_over_T=float(massT),
        unprotected_vacuum_loop_mass_scale_GeV=float(g*M/(4*np.pi)),
        unprotected_loop_mass_scale_over_mphi=float(g*M/(4*np.pi*86e-9)),
        states=states,limitations='Equilibrium coefficients from simple Yukawa model, p_parent/T=3 representative only; frozen PSD is NOT time-resolved and not a self-energy. No bound, rate correction or physical viability claim.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',nargs='+',required=True);p.add_argument('--output',required=True)
    a=p.parse_args()
    out=[screen(row) for path in a.input for row in json.loads(Path(path).read_text())]
    Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for r in out:
        print(json.dumps({k:r[k] for k in ('initial_parent','n','kappa_sm1','kappa','omega','hot','scalar_equilibrium_thermal_mass_over_T','unprotected_loop_mass_scale_over_mphi')}|
              {'max_thermal_mass_over_peak':max(v['scalar_thermal_mass_over_peak'] for v in r['states']),
               'max_occupation_gap_proxy':max(v['occupation_induced_mass2_over_vacuum_gap_proxy'] for v in r['states'])}),flush=True)

if __name__=='__main__':main()
