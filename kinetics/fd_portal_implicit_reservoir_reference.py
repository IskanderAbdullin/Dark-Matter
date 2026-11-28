"""Two-moment closed collision update against a fixed visible reservoir."""
import json
import numpy as np
from scipy.optimize import root

from fd_number_heat_state_reference import evaluate
from fd_portal_differential_reciprocity_reference import PositivePortal,integrate


def dark(T,u,nodes):
    a,b=u
    Td=T/(1-b)
    d=evaluate(Td,a/2,nodes=nodes)
    n,rho=d['n'],d['rho']
    B=d['dn_deta']
    E=d['chemical_energy']
    C=d['cv_fixed_number']
    statejac=np.array([[B/2,B*E/T],[B*E/2,(Td*Td*C+B*E*E)/T]])
    return np.array([n,rho]),float(rho/T-d['entropy']),statejac


def step(p,T,old,nodes,dt):
    y0,q0,_=dark(T,old,nodes)
    def packet(u):
        y,q,A=dark(T,u,nodes)
        f,J,_,_=integrate(p,T,*u)
        # F counts events, so both dark particles in the pair are lost;
        # the second coordinate of the returned flux is heat divided by T.
        v=f[:2]*[2.,T]
        if abs(float(u@f[:2])-f[2])>5e-13*max(abs(f)):
            raise ArithmeticError('Detailed balance entropy ledger lost')
        return (y-y0+dt*v)/y0,(A+dt*np.diag([2.,T])@J)/y0[:,None],y,q,f
    w=root(lambda u:packet(u)[0],old,jac=lambda u:packet(u)[1],tol=1e-11)
    delta,A,y,q,f=packet(w.x)
    scale=max(abs(q),abs(q0),dt*f[2])
    gap=q0-q-dt*f[2]
    if np.max(abs(delta))>4e-13 or min(y)<=0 or f[2]<=0 or gap<-3e-13*scale:
        raise ArithmeticError('Implicit update failed finite local balance or entropy inequality')
    return dict(final_forces=w.x.tolist(),dark_state=y.tolist(),final_step_entropy=f[2]*dt,
                free_energy_decrease=q0-q,strict_remainder=gap,
                normalized_count_heat_residual=np.max(abs(delta)),maximum_absolute_force_update=np.max(abs(w.x-old)))


def run():
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T,a,b in ((.05,.6,.16),(.05,-.6,-.16),(.3,-1.5,.16),(.64135,-.034187877771,.002),(.6882487043,1.8,-.16)):
        u=np.array([a,b])
        y,_,J=dark(T,u,168)
        eps=1e-5
        A=np.stack([(dark(T,u+d,168)[0]-dark(T,u-d,168)[0])/(2*eps) for d in eps*np.eye(2)],axis=1)
        if np.max(abs(A/J-1.))>1e-7:
            raise ArithmeticError('Analytic dark thermodynamic susceptibility differs from direct differences')
        f=integrate(q,T,*u)[0]
        dt=.08/max(abs(f[:2]*[2.,T]/y))
        checks=[]
        for fac in (.5,1.,4.):
            z,v=[step(s,T,u,n,dt*fac) for s,n in ((p,112),(q,168))]
            err=np.max(abs(np.asarray(z['final_forces'])-v['final_forces']))
            if err>2e-7 or v['strict_remainder']<=0.:
                raise ArithmeticError('Unresolved paired implicit solution or nonpositive strict gap')
            checks.append(dict(relative_timestep_factor=fac,force_two_mesh_absolute_defect=err,**v))
        rows.append(dict(T_SM=T,initial_forces=u.tolist(),base_time_in_rate_units=dt,checks=checks))
    return dict(rows=rows,scope='Dark FD and massless-scalar moments against a fixed visible reservoir, on short finite implicit collision steps. Not a cosmological trajectory, an expanding self-consistent visible bath, or a universal timestep bound.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
