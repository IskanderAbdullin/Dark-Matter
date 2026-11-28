"""Explain which invariant-energy regions carry the reciprocal slope changes."""
import json
import numpy as np

from fd_portal_differential_reciprocity_reference import integrate, PositivePortal
from timelike_portal_onsager_coverage_reference import V,width


def measured(q,T,count):
    _,J,_,G=integrate(q,T,0.,0.)
    w=count*width
    bounds=[(0.,V-w),(V-w,V+w),(V+w,np.inf)]
    cells=[integrate(q,T,0.,0.,invariant_window=k) for k in bounds]
    D=np.array([x[1] for x in cells])
    A=np.array([x[3] for x in cells])
    if max(np.max(abs(D.sum(axis=0)/J-1.)),np.max(abs(A.sum(axis=0)/G-1.)))>2e-13:
        raise ArithmeticError('The exact invariant partitions do not reproduce their total')
    M=np.c_[D[:,0,0],D[:,0,1],D[:,1,1],A]
    if np.min(M)<0:
        raise ArithmeticError('Unphysical invariant contribution')
    slopes=np.stack((M[:,1]/M[:,0],M[:,2]/M[:,1],M[:,4]/M[:,3]),axis=1)
    P=np.stack((M[:,0],M[:,1],M[:,3]),axis=1)
    return dict(pole_width_multiple=count,invariant_cut_GeV=[V-w,V+w],
                fractions=M/M.sum(axis=0),conditional_mean_scaled_lab_energy=slopes,
                mean_energy_weights=P/P.sum(axis=0))


def run():
    q=[PositivePortal(x) for x in ((13,36,26),(17,44,32))]
    out=[]
    for T in (.30,.60,.6386086443790474,.64135,.644144559740595,.67,.6882487043):
        lows=[]
        for n in (30,100):
            a,b=[measured(p,T,n) for p in q]
            err=max(np.max(abs(a[x]-b[x])) for x in ('fractions','mean_energy_weights'))
            if err>6e-7:
                raise ArithmeticError('Spectral partition did not resolve on both independent orders')
            lows.append(dict(pole_width_multiple=n,invariant_cut_GeV=b['invariant_cut_GeV'],
                             fractions=b['fractions'].tolist(),conditional_mean_scaled_lab_energy=b['conditional_mean_scaled_lab_energy'].tolist(),
                             mean_energy_weights=b['mean_energy_weights'].tolist(),maximum_paired_fraction_difference=err))
        out.append(dict(T_SM=T,windows=lows))
    return dict(visible_temperature_case=out,partition=['below','near', 'above'],
                fraction_columns=['M0','M1_scaled','M2_scaled2','skew_S0','skew_S1_scaled'],
                moment_mean_columns=['number slope magnitude','heat slope magnitude','reciprocity slope magnitude'],
                scope='Exact partition of the positive quadrature by invariant energy; near-resonance windows deliberately vary. Not independent physical channels or a rate contribution during an evolved history.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
