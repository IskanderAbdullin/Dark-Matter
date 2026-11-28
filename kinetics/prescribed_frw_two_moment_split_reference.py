"""Independent collision/pressure-work check for a prescribed short FRW window."""
import json
import numpy as np
from scipy.special import expit, roots_legendre
from paired_two_moment_fv_reference import Model
from fermion_two_moment_closure_reference import infer as infer_fd
from regular_bose_two_moment_closure_reference import infer as infer_be, occupation
from frw_two_moment_pressure_bounds_reference import k


def streaming(q,m0,b0,m1,b1,y):
    nb,nf,kb,kf=np.split(y[:-1],4)
    cells=len(nb)
    h=np.log(m1/m0)
    assert abs(h-np.log(b1/b0))<4e-15
    if h==0.:
        return y.copy(),0.,0.
    p,w=roots_legendre(112)
    g,z=roots_legendre(10)
    trace=0.
    out=[]
    for mass,new,n,K,bose in ((m0,m1,nb,kb,True),(b0,b1,nf,kf,False)):
        kl,kh=k(q[:-1],mass),k(q[1:],mass)
        lo,hi=q[:-1,None],q[1:,None]
        qp=lo+(hi-lo)*(p+1)/2
        d=(hi-lo)*w/2*qp*qp*(1 if bose else 4)/(2*np.pi**2)
        eps=k(qp,mass)
        Ecut=kh-kl
        f=np.empty_like(qp)
        for j in range(cells):
            fill=n[j]/d[j].sum()
            mean=(K[j]/n[j]-kl[j])/Ecut[j]
            if bose:
                A,B=infer_be(q[j],q[j+1],mass,fill,mean)
                f[j]=occupation(A,B,(eps[j]-kl[j])/Ecut[j])
            else:
                A,B=infer_fd(q[j],q[j+1],mass,fill,mean)
                f[j]=expit(A-B*(eps[j]-kl[j])/Ecut[j])
        population=d*f
        assert np.max(abs(population.sum(1)/n-1))<2e-10
        assert np.max(abs((population*eps).sum(1)/K-1))<2e-10
        out.append((population*k(qp,new)).sum(1))
        for t,wt in zip(g,z):
            mt=mass*np.exp(h*(t+1)/2)
            trace+=h*wt/2*np.sum(population*(mt*mt/np.hypot(qp,mt)))
    old=m0*nb.sum()+b0*nf.sum()+kb.sum()+kf.sum()
    new=m1*nb.sum()+b1*nf.sum()+out[0].sum()+out[1].sum()
    residual=abs((new-old)-trace)
    assert residual<2e-12*max(1.,abs(old),abs(new))
    return np.r_[nb,nf,*out,y[-1]],trace,residual


def integrate(n,steps,method,x0=.1,loga1=.02,collision_time_per_log_a=.05):
    # The time coefficient is prescribed, not a Friedmann closure.
    opt=dict(order=6,thermal="regular")
    first=Model(n,x0,**opt)
    y=first.initial(True)
    initial=first.invariants(y)
    t=0.
    dt=loga1/steps
    total_trace,total_residual=0.,0.
    def redshift(x,p,y):
        return streaming(first.q,.08*x,.18*x,.08*p,.18*p,y)
    rejected=0
    for _ in range(steps):
        xl,xc,xr=(x0*np.exp(v) for v in (t,t+.5*dt,t+dt))
        xm=xc if method=="symmetric" else xl
        y,T,err=redshift(xl,xm,y)
        total_trace+=T
        total_residual+=err
        model=Model(n,xm,**opt)
        assert model.admissible(y)
        rest=collision_time_per_log_a*dt
        while rest>1e-16:
            h=rest
            dy=model.rhs(y)
            while True:
                stage=y+h*dy
                if model.admissible(stage):
                    next_=(y+stage+h*model.rhs(stage))/2
                    if model.admissible(next_): break
                h*=.5
                rejected+=1
                assert h>1e-14 and rejected<100
            y=next_
            rest-=h
        y,T,err=redshift(xm,xr,y)
        total_trace+=T
        total_residual+=err
        t+=dt
    final=Model(n,x0*np.exp(loga1),**opt)
    d=final.invariants(y)-initial
    assert abs(d[0])<5e-13 and abs(d[1]-total_trace)<5e-12
    reference_error=None
    if collision_time_per_log_a==0.:
        direct,_,_=redshift(x0,x0*np.exp(loga1),first.initial(True))
        reference_error=float(np.sum(abs(direct-y)))
    return dict(n=n,steps=steps,method=method,M_initial=x0,log_scale_interval=loga1,
                prescribed_collision_time=collision_time_per_log_a*loga1,
                pure_stream_reference_l1_error=reference_error,
                number_diff_error=float(abs(d[0])),
                energy_minus_pressure_work_error=float(abs(d[1]-total_trace)),
                sum_local_pressure_work_error=float(total_residual),
                integrated_daughter_trace=float(total_trace),last=y.tolist(),rejected=rejected)


def run():
    rows=[integrate(12,n,m) for m in ("lie","symmetric") for n in (3,6,12)]
    for mode in ("lie","symmetric"):
        group=[r for r in rows if r["method"]==mode]
        for i in (0,1):
            group[i]["l1_state_difference_to_next"]=float(
                np.sum(abs(np.array(group[i]["last"])-group[i+1]["last"])))
    s=[r for r in rows if r["method"]=="symmetric"]
    l=[r for r in rows if r["method"]=="lie"]
    assert s[1]["l1_state_difference_to_next"]<s[0]["l1_state_difference_to_next"]/2.5
    assert s[1]["l1_state_difference_to_next"]<l[1]["l1_state_difference_to_next"]/20
    c=[integrate(n,6,"symmetric",m) for m,n in ((.005,12),(.005,24),(.1,24))]
    pure=[integrate(12,n,"lie",.1,.5,0.) for n in (1,3,6)]
    return dict(rows=rows,composition_and_mesh_controls=c,
                isolated_streaming_controls=pure,status="prescribed background and parent; not coupled FRW or relic")
if __name__=="__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
