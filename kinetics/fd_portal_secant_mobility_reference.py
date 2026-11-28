"""A symmetric positive nonlinear mobility is not the differential Jacobian."""
import json
import numpy as np

from fd_portal_differential_reciprocity_reference import PositivePortal,integrate


def record(p,T,force,pauli):
    u,v=force
    f,D,skew,_,M=integrate(p,T,u,v,dark_blocking_fraction=pauli,return_secant=True)
    q=np.asarray(force)
    trial=M@q
    e=np.linalg.eigvalsh(M)
    scale=max(np.max(abs(f)),np.max(abs(M))*np.linalg.norm(q),1e-100)
    err=max(np.max(abs(trial-f[:2])),abs(q@trial-f[2]))/scale
    if min(e)<=0. or err>3e-14:
        raise ArithmeticError('Nonlinear mobility or its entropy law did not reconstruct the old rates')
    if u==v==0. and np.max(abs(D/M-1.))>2e-13:
        raise ArithmeticError('At detailed balance the tangent must coincide with the mobility')
    return dict(analytic_mobility=M.tolist(),tangent_jacobian=D.tolist(),
                number_normalized_heat_and_entropy=f.tolist(),reconstruction_error=float(err),
                eigenvalue_ratio=float(e[1]/e[0]),tangent_minus_secant_scaled=float(
                    np.max(abs(D-M))/np.max(abs(M))),actual_tangent_skew_normalized=float(
                    skew/np.sqrt(D[0,1]*D[1,0])))


def run():
    p,q=(PositivePortal(o) for o in ((13,36,26),(17,44,32)))
    rows=[]
    for T,x in ((.05,(0.,0.)),(.05,(.6,.16)),(.05,(-.6,-.16)),
                (.3,(-1.5,.16)),(.3,(1.5,-.16)),
                (.64135,(-.03418787777103929,.002)),
                (.6882487043,(0.,0.)),(.6882487043,(1.8,-.16))):
        modes=[]
        for k in (0.,1.):
            a,b=[record(r,T,x,k) for r in (p,q)]
            d=np.asarray(a['analytic_mobility']); e=np.asarray(b['analytic_mobility'])
            error=np.max(abs(d/e-1.))
            if error>2e-6 or k==0. and b['actual_tangent_skew_normalized']!=0.:
                raise ArithmeticError('Independent mobility orders or blocking control failed')
            modes.append(dict(dark_blocking_fraction=k,paired_mobility_drift=error,**b))
        rows.append(dict(T_SM=T,forces=x,states=modes))
    return dict(rows=rows,scope='Positive state-dependent secant mobility of the fixed kernel, not the nonlinear Jacobian or a stability test for coupled expansion.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
