"""Conservative fixed-mass scalar shape evolution as a dressed decay-source diagnostic."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp

from early_pole_live_aperture_witness_reprojection_reference import initial_projection
from early_pole_native_aperture_witness_reference import PROFILES
from early_pole_scalar_cold_screened_flux_reference import Screen
from early_pole_scalar_elastic_stencil_reference import Elastic


def run(grid,bins,order,duration,subdivision):
    if bins<4 or order<2 or not duration>0.:
        raise ValueError("invalid scalar mesh, quadrature, or duration")
    p=next(q for q in json.loads(PROFILES.read_text()) if q['grid']==grid)
    C,archive,Y0,Y1,*_=initial_projection(p,subdivision,bins,duration)
    Ns=C.ns.copy()
    _,R0,S0,f0,ms=C.scalar(Ns)
    edges=np.geomspace(1e-6,1e-5,bins+1)
    kernel=Elastic(C.q,ms,edges,order,C.lam,Screen(edges,f0,ms,C.lam,order))
    np.testing.assert_array_equal(C.sw,kernel.vol)
    H=np.sqrt(np.pi**2*106.75/90)*(86./archive['x_end'])/2.435e18
    raw,entropy=kernel.flux(f0)
    init_shape_rate=float(np.sum(abs(raw))/Ns.sum()/H)
    def fun(_,logs):
        N=Ns*np.exp(logs)
        v,_=kernel.flux(N/kernel.vol)
        return v/(H*N)
    times=np.r_[0.,np.geomspace(duration/128,duration,11)]
    first=min(duration/100,0.01/max(1.,np.max(abs(raw/(H*Ns)))))
    h=solve_ivp(fun,(0.,duration),np.zeros_like(Ns),method='Radau',
                t_eval=times,rtol=2e-9,atol=3e-11,max_step=duration/12,first_step=first)
    if not h.success:
        raise RuntimeError(h.message)
    N=Ns[:,None]*np.exp(h.y)
    dn=np.max(abs((N.sum(axis=0)-Ns.sum())/Ns.sum()))
    de=np.max(abs((kernel.K@N-(kernel.K@Ns))/(kernel.K@Ns)))
    out=[]
    for i,exposure in enumerate(times):
        n=N[:,i]
        _,rho,s,f,m=C.scalar(n)
        rates=[]
        for a in (0.,.97):
            X=(1-a)*Y0+a*Y1
            z=np.r_[X.ravel(),n,0.,0.]
            R=C.rhs(z,False)[C.j0:C.j0+C.n]
            rates.append(float(R.sum()))
        out.append(dict(exposure=exposure,signed_number_changes=(n-Ns).tolist(),
                        scalar_L1_shape_fraction=float(np.sum(abs(n-Ns))/Ns.sum()),
                        total_dressed_energy_change=float(rho-R0), entropy_change=float(s-S0),
                        scalar_pole=m,cold_source_step_and_witness=rates))
    assert dn<1e-8 and de<1e-8
    assert min(np.diff([r['entropy_change'] for r in out]))>-2e-28
    return dict(grid=grid,scalar_bins=bins,hard_aperture_subdivision=subdivision,angular_order=order,
                frozen_reference_mass=ms,matched_H_over_T=H,initial_screened_L1_shape_rate=init_shape_rate,
                initial_scalar_entropy_production_per_H=entropy/H,excluded_geometry=kernel.raw_escape,
                number_relative_defect=dn,fixed_dispersion_kinetic_energy_relative_defect=de,rows=out,
                rhs_evaluations=h.nfev,first_solver_step=first,
                scope="Closed finite scalar graph, held initial screening and fixed scattering dispersion; "
                      "read-only dressed cold decay diagnostic with archived charged occupations.")


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('grid',choices=('fine','third'))
    a.add_argument('--bins',type=int,default=8)
    a.add_argument('--order',type=int,default=6)
    a.add_argument('--duration',type=float,default=2e-7)
    a.add_argument('--subdivision',type=int,default=3)
    print(json.dumps(run(**vars(a.parse_args())),indent=2,allow_nan=False))
