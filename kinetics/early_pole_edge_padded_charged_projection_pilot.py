"""Pack thin charged cells against archived steps, conserving energy with two-node events."""
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import G
from early_pole_double_charge_stencil_pilot import case
from early_pole_endpoint_stencil_remap_pilot import LABELS,SPECS,archives,load,parts,remap,MassiveFermion


def pad(e,eta):
    if eta==0:
        return e.copy()
    left,right=e[:-1],e[1:]
    u=right-left
    new=np.sort(np.r_[e,left+eta*u,right-eta*u])
    if not np.array_equal(new[::3],e) or np.any(np.diff(new)<=0):
        raise ArithmeticError('Edge padding moved an archived charge edge')
    return new


def charge_mesh(d,eta_d,eta_p):
    da=pad(d.d.fed,eta_d)
    pa=pad(d.pe,eta_p)
    d.d.fed=da
    d.d.ef=np.hypot(da,d.d.b)
    d.d.wf=4*np.diff(da**3)/(6*np.pi**2)
    q,w=roots_legendre(16)
    step=np.diff(da)[:,None]
    d.daughter_q=da[:-1,None]+step*(q+1)/2
    d.daughter_w=4*step*w*d.daughter_q**2/(4*np.pi**2)
    d.nf=len(da)-1
    d.iP=d.ns+d.nf

    d.pe=pa
    nstep=np.diff(pa)[:,None]
    u,v=roots_legendre(12)
    d.pq=pa[:-1,None]+nstep*(u+1)/2
    d.pqw=G*nstep*v*d.pq**2/(4*np.pi**2)
    d.pw=G*np.diff(pa**3)/(6*np.pi**2)
    d.np=len(pa)-1
    d.n=d.iP+d.np
    d.w=np.r_[d.d.w,d.d.wf,d.pw]


def run():
    logs=[]
    for grid,shape in SPECS.items():
        arc=load(grid)[0]
        x=arc['x_end']
        cells=archives(arc)
        archive=[arc['endpoint_'+s+'_occupations'] for s in LABELS]
        rows=[]
        for eta in (0.,1/8,1/64,1/512,1/4096,1/32768):
            for mode,dpad,ppad in (('both',eta,eta),('daughter',eta,0.),('parent',0.,eta)):
                if eta==0 and mode!='both':
                    continue
                d=MassiveFermion(*shape,arc['mass_ratio'],arc['elastic_rate'],True)
                charge_mesh(d,dpad,ppad)
                remapped=[remap(a,f,b)[0] for a,f,b in zip(cells,archive,parts(d),strict=True)]
                for z,f in ((1,dpad),(2,ppad)):
                    factor=1 if f==0 else 3
                    np.testing.assert_allclose(remapped[z],np.repeat(archive[z],factor),rtol=4e-16,atol=1e-29)
                y=np.concatenate(remapped)
                s=case(d,x,y,orders=(4,8))
                a=s['double_centroid'][-1]
                step=np.asarray(s['exact_step_on_shared_support_gross'])
                geom=np.array([a['gross_forward'],a['gross_inverse']])
                rates=np.r_[geom,geom[0]-geom[1]]
                refs=np.r_[step,step[0]-step[1]]
                rows.append(dict(eta=eta,arm=mode,counts_by_species=[d.ns,d.nf,d.np],
                                 projected_and_explicit_step_rates=[rates.tolist(),refs.tolist()],
                                 rate_ratio_minus_one=(rates/refs-1).tolist(),
                                 quadrature_defect=s['mesh_relative_defect'],
                                 energy_and_thermal_max=max(a['energy_relative_defect'],a['thermal_event_relative_defect']),
                                 omitted_seed=a['omitted_fraction_for_thermal_parent_seed']))
        resolved=next(z for z in rows if z['arm']=='both' and z['eta']==1/32768)
        if max(map(abs,resolved['rate_ratio_minus_one']))>4e-5:
            raise ArithmeticError(('Narrow-edge stencil did not recover archived scalar source',grid,resolved))
        logs.append(dict(grid=grid,rows=rows))
    return logs


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
