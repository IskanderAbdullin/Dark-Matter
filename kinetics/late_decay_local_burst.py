"""Local, short-time A<->B+eta Bose kinetic burst (no cosmic redshift).
Exact relativistic daughter energy and time dilation, finite cell
projection; number invariants exact, energy error measured not assumed.
This can only diagnose onset before line redshifts. No B self-scattering,
in-medium effects, coherence, or earlier/later cosmological evolution.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import BDF
from scipy.optimize import brentq
from scipy.special import roots_legendre
from scipy.sparse import coo_matrix
from late_decay_bose_screen import mesh,vgrid
from direct_activation_consistency import T0
from late_activation_overlap import H
from resonant_mediator_inverse import M,SEC_GYR,C,HBAR
from quartic_consistency import lambda_from_sigma
from late_decay_line_gain import fixed_tau

def geometry(row,r,ps,n,ne,cut,order=4):
    p0=T0*(3.909/row['g_bath'])**(1/3)/row['final']['theta']*row['final']['peak_bin_q']
    P=p0*r;ea,eb,f,meta=mesh(row,n,cut,ps/P*1.04)
    va,qa=vgrid(ea);vb,qb=vgrid(eb);va=va*P**3;vb=vb*P**3
    mb=np.sqrt(M*M-2*M*ps);delta=2*M*ps/(M+mb)
    # x=k/ps-1. Exact maximum range in parent support.
    p=ea[-1]*P;gam1=p*p/(M*(np.hypot(p,M)+M));xmax=1.001*(p/M+gam1)
    b=np.geomspace(xmax*1e-6,xmax,max(3,ne//2))
    xe=np.r_[-b[::-1],0.,b]
    kl=1+xe[:-1];kh=1+xe[1:]
    ve=ps**3*np.diff(xe)*(kl*kl+kl*kh+kh*kh)/(6*np.pi**2)
    # k-ps cell mean stable in x.
    xl=xe[:-1];xh=xe[1:];dx=xh-xl;mid=(xh+xl)/2
    xmean=(mid*(1+mid)**2+(2+3*mid)*dx**2/12)/((1+mid)**2+dx**2/12)
    kb=eb**2*P**2/(np.hypot(eb*P,mb)+mb)
    nA=len(va);nB=len(vb);nE=len(ve);i=[];j=[];k=[];w=[]
    xg,wg=roots_legendre(order)
    for a in range(nA):
        lo,hi=ea[a:a+2];qs=(hi+lo)/2+(hi-lo)*xg/2
        wt=wg*(hi-lo)/2*qs**2/((hi**3-lo**3)/3)
        for q,v in zip(qs,wt):
            p=q*P;Ka=p*p/(np.hypot(p,M)+M);g1=Ka/M;gv=p/M
            k0=Ka+(ps*delta/(M+mb))-ps*g1
            # KB=k0-ps*gv*mu, eta x=g1+gv*mu.
            mus=np.unique(np.clip(np.r_[-1.,1.,(xe-g1)/gv,(k0-kb)/(ps*gv)],-1.,1.))
            midu=(mus[:-1]+mus[1:])/2
            B=np.searchsorted(kb,k0-ps*gv*midu,side='right')-1
            E=np.searchsorted(xe,g1+gv*midu,side='right')-1
            good=(B>=0)&(B<nB)&(E>=0)&(E<nE)
            if not np.all(good):raise ValueError(('uncovered event',a,int(np.sum(~good))))
            i.extend([a]*len(B));j.extend(B+nA);k.extend(E+nA+nB)
            w.extend(v*np.diff(mus)/2/(1+g1))
    i,j,k,w=map(np.array,(i,j,k,w));i=i.astype(int);j=j.astype(int);k=k.astype(int)
    vol=np.r_[va,vb,ve];N0=va@f;F=vol/N0;R=1/(ve.sum()/N0);s=F[i]/F[j]
    # Kinetic energy plus eta-minus-splitting. B-eta invariant eliminates
    # need to subtract the huge common rest-mass term.
    E=np.r_[qa**2*P**2/(np.hypot(qa*P,M)+M),
            qb**2*P**2/(np.hypot(qb*P,mb)+mb),
            ps*xmean-ps*delta/(M+mb)]
    return dict(i=i,j=j,k=k,w=w,F=F,R=R,s=s,E=E,N0=N0,
                nA=nA,nB=nB,nE=nE,p0=p0,P=P,ea=ea,eb=eb,xe=xe,f=f,meta=meta)
def run(row,r,ps,lifetime,n,ne,cut,order,wall,rtol,seed=1e-12,elastic_sigma=0.):
    start=time.monotonic();g=geometry(row,r,ps,n,ne,cut,order)
    i,j,k,w,F,R,s,E,N0=[g[x] for x in ['i','j','k','w','F','R','s','E','N0']]
    nA,nB,nE=[g[x] for x in ['nA','nB','nE']];aend=nA;bend=nA+nB
    ndim=len(F);gamma=1/(lifetime*SEC_GYR);Hn=H(r-1)
    y0=np.r_[F[:nA]*g['f'],np.zeros(nB+nE)]
    y0/=y0.sum()
    q=1/(R*F[k]);bcoef=1/(R*F[j])
    el=None;el_fac=0.;el_checks=None
    if elastic_sigma>0:
        from elastic_implicit_trial import ImplicitElastic
        mb=np.sqrt(M*M-2*M*ps)
        el=ImplicitElastic(g['eb'],mb/g['P'])
        assert np.max(abs(el.vol*g['P']**3/(F[aend:bend]*N0)-1))<1e-10
        el_fac=float(lambda_from_sigma(elastic_sigma)**2*g['P']/(2048*np.pi**5*HBAR*gamma*R))
        # Use actual representative kinetic energies for the C22 invariant.
        E[aend:bend]=el.K*g['P']
        el_checks=dict(events=len(el.i),event_energy_error=el.event_energy_error,compiled=el.compiled)
        print('elastic graph',el_checks,'factor',el_fac,flush=True)

    def calc(y):
        a,b,c=y[i],y[j],y[k]
        return w*(a/R+a*b*bcoef+(a-s*b)*c*q)
    def rhs(t,y):
        v=calc(y)
        out=np.bincount(i,weights=-v,minlength=ndim)+np.bincount(j,weights=v,minlength=ndim)+np.bincount(k,weights=v,minlength=ndim)
        if el is not None:out[aend:bend]+=F[aend:bend]*el.rates_jac(y[aend:bend]/F[aend:bend],False)*el_fac
        return out
    ii=np.r_[i,i,i,j,j,j,k,k,k];jj=np.tile(np.r_[i,j,k],3)
    def jac(t,y):
        a,b,c=y[i],y[j],y[k]
        vals=np.r_[w*(1/R+b*bcoef+c*q),w*(a*bcoef-s*c*q),w*((a-s*b)*q)]
        J=coo_matrix((np.r_[-vals,vals,vals],(ii,jj)),shape=(ndim,ndim)).tocsc()
        if el is not None:
            fb=F[aend:bend]
            Jel=el.rates_jac(y[aend:bend]/fb,True)*el_fac*fb[:,None]/fb[None,:]
            a,b=np.indices(Jel.shape)
            J+=coo_matrix((Jel.ravel(),((a+aend).ravel(),(b+aend).ravel())),shape=(ndim,ndim)).tocsc()
        return J
    base=E@y0
    # Independent frozen invariant integrated gain check.
    lin=gamma*np.bincount(k,weights=w*y0[i]/F[k],minlength=ndim)[bend:]
    tau_num=float(lin@np.diff(np.log1p(g['xe']))/Hn)
    tau_analytic=2*np.pi**2*gamma*N0/(Hn*ps**3)
    maxlin=float(max(lin));cross=g['P']/M/Hn
    # Only integrate up to one peak-Doppler crossing time or artificial
    # scaled ceiling; freeze expansion only meaningful at t/cross <<1.
    tend=min(1e4,cross*gamma*R)
    # BDF's implicit step can follow a negative spurious branch if the
    # physical spontaneous seed is below absolute tolerance. Seed with
    # the exact frozen A, B=0 solution, and explicitly resolve growth.
    gscaled=lin/(gamma*R);fe=F[bend:]
    def frozen(t):
        return fe*np.expm1(np.minimum(gscaled*t,650.))
    ts=brentq(lambda t:frozen(t).sum()-seed,0.,min(tend,100/max(gscaled)),xtol=1e-14) if frozen(min(tend,100/max(gscaled))).sum()>seed else min(tend,100/max(gscaled))
    eta_seed=frozen(ts);weight=w*y0[i];norm=np.bincount(k,weights=weight,minlength=ndim)
    events=weight*eta_seed[k-bend]/np.maximum(norm[k],1e-300)
    transfer=np.bincount(i,weights=-events,minlength=ndim)+np.bincount(j,weights=events,minlength=ndim)+np.bincount(k,weights=events,minlength=ndim)
    yy=y0+transfer
    seed_max_loss=float(max(-transfer[:nA]/np.maximum(y0[:nA],1e-300)))
    if min(yy)<-1e-16 or seed_max_loss>1e-3:raise ValueError(('unsafe analytic seed',seed_max_loss))
    sol=BDF(rhs,ts,yy,tend,jac=jac,rtol=rtol,atol=1e-23,max_step=.15/max(gscaled),first_step=min(.01/max(gscaled),(tend-ts)/2))
    def stat(t,y):
        A=y[:aend].sum();B=y[aend:bend].sum();eta=y[bend:].sum()
        sec=t/(gamma*R)
        return dict(t_scaled=float(t),time_s=sec,time_H=sec*Hn,time_over_peak_crossing=sec/cross,
                    A_fraction=float(A),B_fraction=float(B),eta_fraction=float(eta),
                    massive_number_error=float(A+B-1),eta_minus_B=float(eta-B),
                    energy_error_over_initial_K=float((E@y-base)/base),
                    min_number_fraction=float(min(y)),negative_number=float(-np.minimum(y,0.).sum()),
                    eta_max_f=float(max(y[bend:]/F[bend:])))
    hist=[stat(0.,y0),stat(ts,yy)];first={};snapshots={};steps=0;msg='finished';prev=0.
    while sol.status=='running':
        note=sol.step();steps+=1;st=stat(sol.t,sol.y)
        if st['negative_number']>1e-8:msg='negative state';break
        if sol.status=='failed':msg=str(note);break
        for x in [1e-8,1e-6,.0001,.001,.01,.1,.25,.4,.49,.6]:
            if st['B_fraction']>=x and str(x) not in first:
                first[str(x)]=st
                if x in [.001,.01,.1,.4,.49]:snapshots[str(x)]=(sol.y[aend:bend]/F[aend:bend]).tolist()
        if steps%10==0 or sol.status=='finished':hist.append(st)
        if time.monotonic()-start>wall:msg='wall budget';break
        if st['B_fraction']>1e-3:sol.max_step=max(tend/200,.15/max(gscaled))
        if st['B_fraction']>.49 and st['time_over_peak_crossing']<.02:
            msg='reached 49% within local burst approximation';break
    end=stat(sol.t,sol.y);hist.append(end)
    return dict(recoil_eV=ps,r=r,lifetime_Gyr=lifetime,nA=nA,nB=nB,nEta=nE,events=len(w),angular_order=order,cut=cut,
                p0=g['p0'],N0_eV3=N0,tau_analytic=tau_analytic,tau_bin_numeric=tau_num,relative_tau_error=tau_num/tau_analytic-1,
                fastest_linear_growth_s=maxlin,peak_crossing_s=cross,
                status=msg,elastic_sigma_cm2_g=elastic_sigma,elastic_checks=el_checks,elastic_factor=el_fac,seed_fraction=float(eta_seed.sum()),seed_max_parent_cell_loss=seed_max_loss,steps=steps,wall_s=time.monotonic()-start,first_crossings=first,B_snapshots=snapshots,final=end,history=hist,
                edges_A=g['ea'].tolist(),edges_B=g['eb'].tolist(),eta_x_edges=g['xe'].tolist(),initial_f=g['f'].tolist(),
                final_f=(sol.y/F).tolist(),
                scope='Finite cell local burst. Exact two number invariants; energy error monitored, not exactly conserved by cell assignment. Valid only for elapsed time small compared with line crossing and scattering/dispersion times.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--r',type=float,default=1e5);p.add_argument('--recoil',nargs='+',type=float,default=[1e-5,1e-3,.01]);p.add_argument('--n',type=int,default=100);p.add_argument('--ne',type=int,default=160);p.add_argument('--cut',type=float,default=20.);p.add_argument('--order',type=int,default=4);p.add_argument('--rtol',type=float,default=3e-6);p.add_argument('--lifetime',type=float,default=10.);p.add_argument('--wall',type=float,default=100.);p.add_argument('--seed',type=float,default=1e-12);p.add_argument('--elastic-sigma',type=float,default=0.);a=p.parse_args()
    row=json.loads(Path(a.input).read_text())[0];out=[]
    for ps in a.recoil:
        v=run(row,a.r,ps,a.lifetime,a.n,a.ne,a.cut,a.order,a.wall,a.rtol,a.seed,a.elastic_sigma);out.append(v)
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print('FINAL',json.dumps({k:v for k,v in v.items() if k not in ['history','edges_A','edges_B','eta_x_edges','initial_f','final_f','B_snapshots']}),flush=True)
