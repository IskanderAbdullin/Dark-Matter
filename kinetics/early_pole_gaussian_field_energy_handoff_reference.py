"""An energy-matched zero-mean Gaussian ensemble and local coherent-band diagnostic."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc

from early_pole_coherent_quartic_floquet_reference import background,edges


def coverage(lo,hi,w,q):
    def scan(a):
        ix=np.argsort(a)
        c=np.r_[0.,np.cumsum(w[ix])]
        return c[np.searchsorted(a[ix],q,side='right')]
    return (scan(lo)-scan(hi))/w.sum()


def one(spec,order):
    sr=spec['checkpoints'][-1]
    q=np.array(spec['scalar_momentum_nodes'])
    n=np.maximum(np.array(spec['endpoint_scalar_signed_shell_number_over_T3']),0.)
    lam=7.65175e-8
    m0=1e-9*spec['x_end']
    m=sr['scalar_pole_over_T']
    E=np.hypot(q,m)
    var=n@(1/E)
    varpi=n@E
    grad=n@(q*q/E)
    rho_gauss=.5*(varpi+grad+m0*m0*var)+lam*var*var/8
    rho_ref=sr['integrated_scalar_energy_over_T4']
    g0=abs(rho_gauss/rho_ref-1)
    assert g0<2e-11
    u=qmc.Sobol(d=2,scramble=True,seed=623519).random_base2(order)
    f,p=(ndtri(u)*np.array([np.sqrt(var),np.sqrt(varpi)])).T
    rho=.5*(p*p+m0*m0*f*f)+lam*f**4/24
    energy_overlap=abs(rho.mean()/(rho_ref-.5*grad)-1)
    assert energy_overlap<.001
    _,Ou,ku=background(rho_ref,m0,lam)
    A,O,k=background(rho,m0,lam)
    low,high=edges(k)
    lo,hi=O*np.sqrt(low),O*np.sqrt(high)
    qref=Ou*np.sqrt(sum(edges(ku))/2)
    v=np.geomspace(.08*qref,4*qref,251)
    en=coverage(lo,hi,rho,v)
    vol=coverage(lo,hi,np.ones(len(lo)),v)
    i=int(np.argmax(en))
    return dict(x=spec['x_end'],mass_ratio=spec['mass_ratio'],samples=len(rho),
                coordinate_zero_mean_variance_over_T2=var,velocity_variance_over_T4=varpi,
                gradient_energy_fraction=float(.5*grad/rho_ref),
                gaussian_hartree_energy_identity_error=g0,mc_local_energy_mean_error=float(energy_overlap),
                oscillator_amplitude_absolute_quantiles_over_T=np.quantile(A,[.1,.5,.9]).tolist(),
                local_Omega_quantiles_over_T=np.quantile(O,[.1,.5,.9]).tolist(),
                uniform_energy_matched_resonant_mid_q=qref,
                energy_weighted_cover_at_uniform_mid_q=float(coverage(lo,hi,rho,np.array([qref]))[0]),
                peak_local_band_energy_coverage=float(en[i]),local_band_energy_coverage_peak_q=float(v[i]),
                peak_local_band_volume_coverage=float(max(vol)),
                q_grid=v.tolist(),weighted_coverage=en.tolist(),volume_coverage=vol.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    a=p.parse_args()
    result={}
    for name in a.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        s=json.loads(raw)
        rows=[one(s,i) for i in (16,18)]
        err=abs(rows[0]['peak_local_band_energy_coverage']-rows[1]['peak_local_band_energy_coverage'])
        assert err<.004
        result[f.name]=rows
    print(json.dumps(result,indent=2,allow_nan=False))
