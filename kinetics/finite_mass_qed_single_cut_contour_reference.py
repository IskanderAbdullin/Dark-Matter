"""Contour algebra and a local absorptive sign control for a scalar-latest triangle."""
import gzip,json,itertools
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import Q,intervals
from finite_mass_qed_charged_complex_guard_reference import chart


def contour_assignment():
    """3-point sum with latest fixed '+' and oriented edges 0->1->2->0."""
    X,R,A,T=[sp.symbols(z+':3') for z in 'XRAT']
    matrices=[sp.Matrix([[R[i]+X[i]+T[i],X[i]+T[i]],
                         [R[i]-A[i]+X[i]+T[i],X[i]-A[i]+T[i]]]) for i in range(3)]
    rows=[]
    for latest in range(3):
        val=0
        for hs in itertools.product((0,1),repeat=3):
            if hs[latest]: continue
            g=(-1)**sum(hs)
            for k in range(3):
                g*=matrices[k][hs[k],hs[(k+1)%3]]
            val+=g
        v=sp.Poly(sp.expand(val),*T)
        if any(sum(p)>1 for p,c in v.terms()) or len(v.terms())!=4:
            raise ArithmeticError('Closed contour did not remove multiple thermal statistical insertions.')
        c=[v.coeff_monomial(T[i]) for i in range(3)]
        expected=[None]*3
        expected[latest]=A[(latest+1)%3]*A[(latest+2)%3]
        expected[(latest+1)%3]=R[latest]*A[(latest+2)%3]
        expected[(latest+2)%3]=R[latest]*R[(latest+1)%3]
        if c!=expected:
            raise ArithmeticError('Uncut causal assignments did not match direct contour expansion.')
        rows.append(dict(latest_vertex=latest,single_thermal_terms=[str(u) for u in c]))
    return rows


def absorptive(d,E,cut,sign):
    M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
    evaluate,_,Kin,C,Z,p,w,mass,A,B=chart(M,m,s,E,cut)
    name='parent' if cut=='parent' else 'fermion'
    bounds=np.hypot(mass,d[f'endpoint_{name}_momentum_edges'])-w
    occupations=np.asarray(d[f'endpoint_{name}_occupations'])
    i=np.searchsorted(bounds,0.)-1
    f=occupations[i]
    L=3e-4
    if min(abs(bounds))<L:
        L=min(abs(bounds))*.6
    expected=sign*np.pi*Z*(2*f)/(2*Q*p*p*np.hypot(A,B))
    X,W=roots_legendre(26)
    def I(eta):
        knots=np.geomspace(1e-22,L,38)
        out=0.
        for l,r in zip(knots[:-1],knots[1:],strict=True):
            a,b=np.log(l),np.log(r)
            x=np.exp((a+b+(b-a)*X)/2)
            z=[t*(np.imag(evaluate(t,eta*sign,True))*Kin(t)[-1]+
                  np.imag(evaluate(-t,eta*sign,True))*Kin(-t)[-1])*f/2 for t in x]
            out+=(b-a)/2*(W@z)
        return float(out)
    a,b=I(1e-13),I(1e-14)
    measured=(b-a)/np.log(10.)
    if abs(measured-expected)>3e-5*max(1.,abs(expected)):
        raise ArithmeticError(f'Assigned positive charged sign does not replay local absorptive slope: {measured}, {expected}.')
    return dict(cut=cut,charged_denominator_eta_sign=sign,occupation=float(f),local_radial_half_width=L,
                analytic_absorptive_log_coefficient=expected,numerical_absorptive_log_coefficient=float(measured),
                absolute_angular_and_radial_defect=float(abs(measured-expected)))


def run():
    grids=[]
    for grid in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        Es,lo,hi,*_=intervals(M,s,m,Q)
        out=[]
        for t in (.22,.5,.78):
            E=lo+t*(hi-lo)
            values=[absorptive(d,E,'parent',1),absorptive(d,E,'daughter',-1)]
            slopes=np.asarray([x['analytic_absorptive_log_coefficient'] for x in values])
            dPb=(M*M+m*m-s*s)/2
            Pq=(M*M-m*m+s*s)/2
            bq=(M*M-m*m-s*s)/2
            H1=np.sqrt((Pq-M*s)*(Pq+M*s))
            H2=np.sqrt((bq-m*s)*(bq+m*s))
            jshape=np.pi*dPb/H1
            analytic=jshape*(values[0]['occupation']-values[1]['occupation'])
            invariant=max(abs(H1/H2-1),abs(analytic-float(sum(slopes))),
                          abs(slopes[0]/values[0]['occupation']-jshape),
                          abs(slopes[1]/values[1]['occupation']+jshape))
            if invariant>2e-12:
                raise ArithmeticError('Absorptive overlap did not map to one shared invariant coefficient.')
            out.append(dict(parent_energy=E,born_aperture_fraction=t,positive_charged_cuts=values,
                            local_coefficient_per_occupation=jshape,
                            signed_coefficient_from_occupation_difference=analytic,invariant_matching_defect=invariant,
                            signed_sum=float(sum(slopes)),sum_of_magnitudes=float(sum(abs(slopes)))))
        grids.append(dict(grid=grid,rows=out))
    return dict(scope='Latest scalar insertion on a directed one-loop triangle and mathematical single-cut absorptive check; no complete retarded medium rate or width inferred.',
                loop_edge_order=['daughter','photon','parent'],
                symbolic_causal_assignment=contour_assignment(),grids=grids)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
