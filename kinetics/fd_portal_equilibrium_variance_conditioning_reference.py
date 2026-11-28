"""Locate the laboratory energy variance that prevents a rank-one heat/number map."""
import json
import numpy as np

from fd_portal_differential_reciprocity_reference import integrate,PositivePortal
from timelike_portal_onsager_coverage_reference import V,width


def summary(D):
    a,b,c=D
    mu=b/a
    var=c/a-mu*mu
    if min(a,b,c,var)<=0.:
        raise ArithmeticError('Expected a positive nontrivial lab energy distribution')
    rho=b/np.sqrt(a*c)
    return dict(mean_scaled_energy=mu,variance_scaled_energy=var,
                heat_schur_over_heat_diagonal=var/(c/a),
                normalized_condition_number=(1+rho)/(1-rho))


def case(p,T,number):
    _,D,_,_=integrate(p,T,0.,0.)
    full=np.array([D[0,0],D[0,1],D[1,1]])
    w=number*width
    bands=((0.,V-w),(V-w,V+w),(V+w,np.inf))
    pieces=[]
    for r in bands:
        _,x,_,_=integrate(p,T,0.,0.,invariant_window=r)
        pieces.append([x[0,0],x[0,1],x[1,1]])
    y=np.asarray(pieces)
    if np.max(abs(y.sum(axis=0)/full-1.))>3e-13:
        raise ArithmeticError('Invariant pieces do not reproduce the unsplit lab moments')
    norm=full[0]
    base=summary(full)
    mu=np.array([summary(k)['mean_scaled_energy'] for k in y])
    vr=np.array([summary(k)['variance_scaled_energy'] for k in y])
    P=y[:,0]/norm
    inside=P*vr
    spread=P*(mu-base['mean_scaled_energy'])**2
    if abs(np.sum(inside+spread)/base['variance_scaled_energy']-1)>5e-11:
        raise ArithmeticError('Independent law of total lab variance failed')
    laws=[]
    for kind,m in (('within_pole_region_mean',mu[1]),('invariant_pole_energy',V/T)):
        a=y.copy()
        a[1]=y[1,0]*np.array([1.,m,m*m])
        laws.append(dict(collapsed_pole_lab_energy=kind,reconstructed_moments=a.sum(axis=0).tolist(),**summary(a.sum(axis=0))))
    return dict(pole_widths=number,actual_lab_moments=full.tolist(),actual=base,
                total_variance_partition_by_invariant_region=[inside.tolist(),spread.tolist()],
                variance_fraction_inside_near_region=float(inside[1]/base['variance_scaled_energy']),
                invariant_region_weights=P.tolist(),collapse_controls=laws)


def run():
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T in (.30,.60,.64135,.6882487043):
        checks=[]
        for z in (30,100):
            a,b=[case(s,T,z) for s in (p,q)]
            vals=lambda x:np.r_[x['actual']['normalized_condition_number'],
                               [c['normalized_condition_number'] for c in x['collapse_controls']]]
            delta=np.max(abs(vals(a)/vals(b)-1))
            if delta>6e-6:
                raise ArithmeticError('Independent actual and collapsed laws do not resolve conditioning')
            checks.append(dict(max_order_condition_defect=delta,**b))
        rows.append(dict(T_SM=T,controls=checks))
    return dict(rows=rows,partition=['below','near','above'],
                scope='Positive equilibrium lab energy moments of a local timelike charged-lepton FD portal. Resonance collapse is deliberately a counterfactual, not a physical model or trajectory.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
