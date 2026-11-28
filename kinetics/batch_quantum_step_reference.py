"""Order-independent timestep rejection on an illustrative two-species network."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit


def network(n=24):
    M,m,b=1.,.08,.18
    eb=np.geomspace(.081,6.,n)
    ef=np.geomspace(.181,8.,n)
    pb=np.sqrt((eb-m)*(eb+m))
    pf=np.sqrt((ef-b)*(ef+b))
    K=(M*M-m*m-b*b)/2
    mask=abs(eb[:,None]*ef[None,:]-K)<=pb[:,None]*pf[None,:]
    i,j=np.nonzero(mask)
    parent_energy=eb[i]+ef[j]
    assert np.all(parent_energy>=M)
    capb=np.geomspace(.1,1.,n)
    capf=np.geomspace(.1,1.,n)
    f1=expit(-parent_energy)
    weight=capb[i]*capf[j]*.7
    def rhs(u):
        a,c=u[i]/capb[i],u[n+j]/capf[j]
        flux=weight*(f1*(1-c)*(1+a)-(1-f1)*c*a)
        return np.r_[np.bincount(i,weights=flux,minlength=n),
                     np.bincount(j,weights=flux,minlength=n),-flux@parent_energy]
    bounds=np.r_[np.full(n,np.inf),capf]
    eq=np.r_[capb/np.expm1(eb),capf*expit(-ef),0.]
    assert np.max(abs(rhs(eq)))<1e-15
    start=eq.copy()
    start[:n]*=.15
    start[n:2*n]=capf*np.linspace(.98,.02,n)
    return rhs,bounds,eq,start,eb,ef,len(i)


def heun(rhs,bounds,start,tend,tol,max_step=np.inf):
    """Global adaptive update; reject out-of-domain candidates without clipping."""
    t,h,u=0.,tend,start.copy()
    accepted, rejected_bounds, rejected_tol=0,0,0
    def valid(y):
        return np.all(y[:-1]>=0) and np.all(y[:-1]<=bounds)
    while tend-t > 8*np.spacing(tend):
        h=min(h,tend-t,max_step)
        assert h>0 and t+h>t
        x=u+h*rhs(u)
        if not valid(x):
            h*=.5
            rejected_bounds+=1
            continue
        y=(u+x+h*rhs(x))/2
        if not valid(y):
            h*=.5
            rejected_bounds+=1
            continue
        error=np.max(abs(y-x)/(1+np.maximum(abs(y),abs(u))))
        if error>tol:
            h*=max(.2,.85*np.sqrt(tol/error))
            rejected_tol+=1
            continue
        t,u=t+h,y
        accepted+=1
        h*=min(2.,max(.3,.85*np.sqrt(tol/max(error,1e-300))))
    return u,(accepted,rejected_bounds,rejected_tol)


def run():
    rhs,cap,eq,start,eb,ef,events=network()
    t=5.
    exact=solve_ivp(lambda t,y:rhs(y),(0,t),start,method="DOP853",rtol=2e-13,atol=1e-14)
    assert exact.success
    exact=exact.y[:,-1]
    rows=[]
    for tol in (1e-3,1e-4,1e-5):
        y,counters=heun(rhs,cap,start,t,tol)
        d=y-start
        n=len(eb)
        err_n=abs(d[:n].sum()-d[n:2*n].sum())
        err_e=abs(eb@d[:n]+ef@d[n:2*n]+d[-1])
        rows.append(dict(tol=tol,max_absolute_error_vs_dop853=float(np.max(abs(y-exact))),
                         pair_number_diff_residual=float(err_n),energy_including_bath_residual=float(err_e),
                         accepted=counters[0],rejected_bounds=counters[1],rejected_error=counters[2]))
        assert err_n<2e-12 and err_e<2e-12
        assert np.all(y[:-1]>=0) and np.all(y[:-1]<=cap)
    errors=[r["max_absolute_error_vs_dop853"] for r in rows]
    assert all(a>b for a,b in zip(errors,errors[1:]))
    assert errors[-1]<1e-5
    e,counts=heun(rhs,cap,eq,t,1e-5)
    eqdrift=float(np.max(abs(e-eq)))
    fine_eq,_=heun(rhs,cap,eq,t,1e-5,max_step=.1)
    eqdrift_cap=float(np.max(abs(fine_eq-eq)))
    assert eqdrift<1e-12 and eqdrift_cap<1e-14
    return dict(cells_per_species=len(eb),allowed_node_pairs=events,rows=rows,
                equilibrium_drift_unrestricted_step=eqdrift,
                equilibrium_drift_step_capped_at_p1=eqdrift_cap)

if __name__=="__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
