"""Replay the physical occupation and positive-energy guards of the two bare poles."""
import gzip,json
from pathlib import Path
import numpy as np
from finite_mass_qed_nonnull_photon_ward_poles_reference import setup,dot


def run():
    root=Path(__file__).parent/'results'
    bag=json.loads((root/'finite-mass-qed-nonnull-photon-ward-poles-2026-11-27.json').read_text())
    data=[]
    for g in bag['data']:
        with gzip.open(root/f"early-pole-spectrum-{g['grid']}-m06-2026-11-26.json.gz") as z:d=json.load(z)
        worst=0.;minprob=1.;count=0
        for row in g['checks']:
            k,w=row['k'],row['omega']
            c={(.0008,.0002):0.,(.00128,.00064):-.5,(.00256,.00064):.5}.get((k,w),0.)
            M,m,q,K,T,geo=setup(d,k,w,c)
            for ph in row['azimuths']:
                roots=ph['found_near_born_intermediate_on_shell_roots']
                if len(roots)!=(2 if row['vacuum_photon_invariant']<0. else 0):
                    raise ArithmeticError('Unexplained root count')
                for z in roots:
                    P,b=geo(z['parent_E'],ph['azimuth'])
                    for vec,label in ((P,'parent'),(b,'fermion')):
                        f=np.r_[0.,d[f'endpoint_{label}_occupations'],0.][np.searchsorted(
                            d[f'endpoint_{label}_momentum_edges'],np.linalg.norm(vec[1:]))]
                        minprob=min(minprob,f,1.-f)
                    inner=P+K if z['arm']==0 else b-K
                    mass=M if z['arm']==0 else m
                    if inner[0]<=mass:raise ArithmeticError('The intermediate state is not future-directed')
                    worst=max(worst,abs(dot(inner,inner)-mass*mass));count+=1
        if count!=24 or worst>2e-15 or minprob<=0.:
            raise ArithmeticError('On-shell occupancy reconstruction failed')
        data.append(dict(grid=g['grid'],future_occupied_sequential_roots=count,
                         minimum_occupied_and_unblocked_probability=float(minprob),
                         maximum_abs_internal_mass_shell_defect=float(worst)))
    return data


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
