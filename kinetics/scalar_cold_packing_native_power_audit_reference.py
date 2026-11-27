"""Audit mixed-scalar exact native k against a matched pure-CDM spectrum."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson


def power(p,k):
    d=np.loadtxt(p)
    j=np.searchsorted(d[:,0],k)
    np.testing.assert_array_equal(d[j,0],k)
    assert np.all(d[j,1]>0)
    return d[j,1].tolist(),hashlib.sha256(p.read_bytes()).hexdigest()


def run(root):
    spec=json.loads((root/'spec.json').read_text())
    k=np.array(spec['modes_h_per_Mpc'])
    cases=[]
    psd=[]
    for tier in ('compact','dense'):
        f=root/f'{tier}-packed.dat'
        u,d=np.loadtxt(f).T
        g=CubicSpline(u,d)
        um=(u[:-1]+u[1:])/2
        assert np.all(g(um)>=0)
        vals=simpson(d*u*u,x=u)
        first=simpson(d*u**4,x=u)/vals
        first_abs=simpson(d*u**3,x=u)/vals
        shape_checks=[]
        for label,cap in (('F1',1.),('F200000',200000.)):
            cfg={k.strip():v.strip() for k,v in (line.split('=',1) for line in
                 (root/f'{tier}-{label}.ini').read_text().splitlines())}
            scale=float(cfg['T_ncdm'])
            nc=cap*scale**3*vals/(2*np.pi**2)
            dneff=nc*scale*first_abs/(np.pi**2/15)/(7/8*(4/11)**(4/3))
            assert float(cfg['m_ncdm'])==86 and float(cfg['omega_ncdm'])==.06
            shape_checks.append(dict(label=label,scalar_number_over_Tgamma0_cubed=nc,
                                     present_comoving_mean_q_over_Tgamma0=scale*first_abs,
                                     ideal_ultrarelativistic_delta_Neff=dneff))
        a,b=shape_checks
        assert abs(a['scalar_number_over_Tgamma0_cubed']/b['scalar_number_over_Tgamma0_cubed']-1)<1e-12
        assert abs(a['ideal_ultrarelativistic_delta_Neff']/b['ideal_ultrarelativistic_delta_Neff']/200000.**(1/3)-1)<1e-10
        psd.append(dict(tier=tier,integral_u2=vals,second_u=first,
                        independently_computed_psd_moments=shape_checks,
                        sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    for zind,z in ((1,5),(2,20)):
        out={}
        for tier in ('compact','dense'):
            cdm,hs=power(root/f'{tier}-cdm_00_z{zind}_pk.dat',k)
            v={'exact_cdm_power':cdm,'cdm_sha256':hs}
            for label in ('F1','F200000'):
                p,sha=power(root/f'{tier}-{label}_00_z{zind}_pk.dat',k)
                v[label]={'power':p,'ratio':(np.array(p)/cdm).tolist(),'sha256':sha}
            out[tier]=v
        cds=np.array(out['dense']['exact_cdm_power'])/out['compact']['exact_cdm_power']-1
        check=dict(cdm_max_relative_drift=float(np.max(abs(cds))))
        for label in ('F1','F200000'):
            hi=np.array(out['dense'][label]['ratio'])
            lo=np.array(out['compact'][label]['ratio'])
            check[label]={'max_relative_ratio_drift':float(np.max(abs(hi/lo-1))),
                          'max_absolute_ratio_drift':float(np.max(abs(hi-lo)))}
        assert check['cdm_max_relative_drift']<1e-3
        cases.append(dict(z=z,profiles=out,tier_and_edge_joint_check=check))
    return dict(scope="Approximate smooth near-packed shapes, with jointly different "
                      "edge widths, native resolution and hierarchy; a counterfactual fixed-abundance "
                      "input, not the thermal production branch. The coldest ballistic shape is "
                      "not proved to maximize linear power; no observational inference.",
                modes_h_per_Mpc=k.tolist(),psd_checks=psd,results=cases)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.root),indent=2,allow_nan=False))
