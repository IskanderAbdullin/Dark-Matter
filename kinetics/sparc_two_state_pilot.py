"""Reproducible one/two-state isothermal Bose halo SPARC pilot.

No condensate, cutoff, covariance, pressure support or selection likelihood.
Use official observations and Q<=2, i>=30, radial interleaved holdout.
The two-state mass gap is a specified EFT scaling, not a fitted parameter.
"""
import argparse,json,time,hashlib,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,minimize_scalar,brentq
from scipy.interpolate import PchipInterpolator
from scipy.special import zeta
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'continuation-02'/'code'))
from sparc_likelihood import read_sparc,nu_from_mass
from sparc_quality import parse
from two_state_halo_equilibrium import pl
from soft_mixing_activation import scales
from resonant_mediator_inverse import M,C

K86=nu_from_mass(86.,1.,1.)

def build(nr=700,fine=False):
    A=np.unique(np.r_[0,np.geomspace(.002,.2,19 if fine else 11),np.linspace(.2,2.,31 if fine else 19),np.linspace(2.,6.,21 if fine else 11)])
    D=np.unique(np.r_[0,np.geomspace(.001,8.,40 if fine else 25)])
    s=np.linspace(np.log(1e-5),np.log(1e4),nr+1);h=s[1]-s[0]
    # Dense, accurate scalar log-polylog table. Interpolate sqrt exponent.
    xx=np.linspace(0,12,24001);yy=np.log(pl(1.5,xx*xx))
    def li(a):
        a=np.maximum(a,0.)
        return np.exp(np.interp(np.sqrt(a),xx,yy)-(np.maximum(a-144.,0.)))
    aa,dd=np.meshgrid(A*A,D*D,indexing='ij')
    n=li(aa)+li(aa+dd)
    Y=np.full_like(aa,np.exp(2*s[0])/6);U=2*Y.copy()
    hist=np.empty((len(A),len(D),nr+1));hist[:,:,0]=U
    def rhs(t,y,u):
        return u,np.exp(2*t)*(li(aa+y)+li(aa+dd+y))/n-u
    for i,t in enumerate(s[:-1]):
        k1,l1=rhs(t,Y,U);k2,l2=rhs(t+h/2,Y+h*k1/2,U+h*l1/2)
        k3,l3=rhs(t+h/2,Y+h*k2/2,U+h*l2/2);k4,l4=rhs(t+h,Y+h*k3,U+h*l3)
        Y+=h*(k1+2*k2+2*k3+k4)/6;U+=h*(l1+2*l2+2*l3+l4)/6
        hist[:,:,i+1]=U
    assert np.all(hist>0)
    return dict(A=A,D=D,s=s,logu=np.log(hist),li_x=xx,li_log=yy)

class Library:
    def __init__(self,table):
        self.A=table['A'];self.D=table['D'];self.s=table['s'];self.u=table['logu'];self.xx=table['li_x'];self.yy=table['li_log']
    def li(self,a):
        return np.exp(np.interp(np.sqrt(a),self.xx,self.yy)-np.maximum(a-144.,0.))
    def nu(self,a,d,two):
        return self.li(a*a)+(self.li(a*a+d) if two else 0.)
    def shape(self,a,d,xi):
        # d=0 has two identical species with same normalized shape as one.
        x=np.sqrt(max(d,0.));i=int(np.clip(np.searchsorted(self.A,a)-1,0,len(self.A)-2))
        j=int(np.clip(np.searchsorted(self.D,x)-1,0,len(self.D)-2))
        w=np.clip((a-self.A[i])/(self.A[i+1]-self.A[i]),0.,1.)
        v=np.clip((x-self.D[j])/(self.D[j+1]-self.D[j]),0.,1.)
        f=(1-w)*((1-v)*self.u[i,j]+v*self.u[i,j+1])+w*((1-v)*self.u[i+1,j]+v*self.u[i+1,j+1])
        return np.where(xi<np.exp(self.s[0]),xi*xi/3,np.exp(np.interp(np.log(xi),self.s,f)))

def validate(lib):
    from scipy.integrate import solve_ivp
    rng=np.random.default_rng(127);out=[]
    for a,d in [(0.,0.),(0.,.02),(.13,.027),(.47,.15),(1.4,1.7),(3.1,20.)]+[(float(rng.uniform(0,2)),float(10**rng.uniform(-3,1))) for _ in range(8)]:
        nu=float(pl(1.5,a*a)+pl(1.5,a*a+d));r=np.geomspace(1e-4,100,100)
        f=lambda x,y:[y[1],np.exp(2*x)*float(pl(1.5,max(a*a+y[0],0))+pl(1.5,max(a*a+d+y[0],0)))/nu-y[1]]
        sol=solve_ivp(f,(np.log(r[0]),np.log(r[-1])),[r[0]**2/6,r[0]**2/3],t_eval=np.log(r),rtol=1e-9,atol=1e-12)
        assert sol.success
        err=float(max(abs(lib.shape(a,d,r)/sol.y[1]-1)))
        out.append(dict(sqrt_alpha=a,d=d,max_relative_V2_error=err))
    return out

def data_arrays(rows):
    rows=sorted(rows,key=lambda x:x.R)
    return np.array([[x.R,x.Vobs,x.sigma,x.Vgas,x.Vdisk,x.Vbul] for x in rows]).T

def fit(data,lib,kind,mass=86.,sigma=.3,prior=False,metadata=None,starts=4):
    r,obs,err,vg,vd,vb=data
    train=np.arange(len(r))%4!=0
    q=np.ones(len(r),dtype=bool)
    two=kind=='two';bose=kind in ('one','two')
    gap=scales(sigma)['mass_gap_eV']*(mass/86.)**2.5 if two else 0.
    if bose:
        # Free sqrt(alpha), log(V0), with R derived from physical density.
        bounds=[(0.,6.),(np.log(3.),np.log(650.))]
    else:bounds=[(np.log(.08),np.log(120.)),(np.log(3.),np.log(650.))]
    if prior:bounds+= [(.05,1.5),(.05,2.),(-3.,3.),(-3.,3.)]
    if prior and metadata is None:raise ValueError('metadata needed')
    def predict(p):
        v=np.exp(p[1]);a=p[0]
        if prior:
            ml,mb,sd,si=p[2:]
            ds=1+sd*metadata['distance_error_Mpc']/metadata['distance_Mpc']
            incl=metadata['inclination_deg']+si*metadata['inclination_error_deg']
            if ds<=.05 or incl<5 or incl>90:return None,None,None
            incl_scale=np.sin(np.deg2rad(incl))/np.sin(np.deg2rad(metadata['inclination_deg']))
        else:ml,mb,ds,incl_scale=.5,.7,1.,1.
        if bose:
            d=gap/(mass*(v*1e5/C)**2)
            nu=lib.nu(a,d,two)
            R=np.sqrt(K86*(86/mass)**4/(v*nu))
            if R<.08 or R>120.:return None,None,None
            h=v*v*lib.shape(a,d,r*ds/R)
        else:
            R=np.exp(a);xi=r*ds/R
            if kind=='classical':h=v*v*lib.shape(6.,0.,xi)
            elif kind=='pis':h=v*v*np.where(xi<.001,xi*xi/3,1-np.arctan(xi)/xi)
            else:h=v*v*np.where(xi<.001,4*xi*xi/3-xi**3,(2*np.log1p(xi)+np.log1p(xi*xi)-2*np.arctan(xi))/xi)
        vv=h+ds*(vg*abs(vg)+ml*vd*vd+mb*vb*vb)
        if np.any(vv<=0):return None,None,None
        pp=0. if not prior else ((ml-.5)/.1)**2+((mb-.7)/.14)**2+sd*sd+si*si
        return np.sqrt(vv)*incl_scale,R,pp
    def fun(p):
        y,R,pp=predict(p)
        if y is None:return 1e12
        return float(np.sum(((y[train]-obs[train])/err[train])**2)+pp)
    vs=np.unique(np.clip(np.r_[np.geomspace(5,500,13),np.median(obs)/np.sqrt(2),np.median(obs)],3,650))
    if bose:aa=np.r_[0.,.025,.1,.25,.5,.8,1.1,1.5,2.,2.5,3.,4.,5.]
    else:aa=np.log(np.geomspace(.1,100,13))
    c=[]
    for a in aa:
        for v in vs:
            p=np.r_[a,np.log(v),[.5,.7,0.,0.] if prior else []]
            c.append((fun(p),p))
    if bose:
        # Independent coarse R,V grid used by the legacy reference, mapped
        # into sqrt(alpha) so the optimizer never crosses the Bose bound.
        for v in np.geomspace(7,450,13):
            d=gap/(mass*(v*1e5/C)**2)
            for R in np.geomspace(.2,70,13):
                target=K86*(86/mass)**4/(v*R*R)
                if not lib.nu(6.,d,two)<=target<=lib.nu(0.,d,two):continue
                a=brentq(lambda a:lib.nu(a,d,two)-target,0.,6.,xtol=1e-8)
                pp=np.r_[a,np.log(v),[.5,.7,0.,0.] if prior else []]
                c.append((fun(pp),pp))
    c.sort(key=lambda x:x[0])
    best=c[0]
    if best[0]>=1e12:
        return dict(success=False,message='No feasible coarse point',train=None,test=None)
    chosen=[]
    for score,p in c:
        if score>=1e12:break
        if all(np.linalg.norm(p[:2]-old[:2])>.35 for old in chosen):chosen.append(p)
        if len(chosen)>=starts:break
    # Ensure two distinct degeneracy regimes get a local start if feasible.
    if bose:
        for lim in [lambda p:p[0]<.15,lambda p:p[0]>1.3]:
            candidate=next((p for score,p in c if score<1e12 and lim(p)),None)
            if candidate is not None and all(np.linalg.norm(candidate[:2]-x[:2])>.2 for x in chosen):chosen.append(candidate)
    for pp in chosen:
        simplex=[pp.copy()]
        steps=np.array([.16,.13]+([.08,.10,.25,.25] if prior else []))
        for i,step in enumerate(steps):
            y=pp.copy();y[i]=pp[i]+step if pp[i]+step<=bounds[i][1] else pp[i]-step
            simplex.append(y)
        opt=minimize(fun,pp,method='Nelder-Mead',bounds=bounds,options={'initial_simplex':np.array(simplex),'maxiter':440 if prior else 200,'xatol':1e-6,'fatol':1e-5})
        if opt.fun<best[0]:best=(float(opt.fun),opt.x)
    # Explicit one-dimensional minimization along the physical critical edge.
    if bose:
        p=best[1].copy();p[0]=0.
        vgrid=np.linspace(bounds[1][0],bounds[1][1],40)
        def fcrit(x):
            t=p.copy();t[1]=x;return fun(t)
        vals=np.array([fcrit(x) for x in vgrid]);j=int(vals.argmin())
        lo,hi=vgrid[max(j-1,0)],vgrid[min(j+1,len(vgrid)-1)]
        opt=minimize_scalar(fcrit,bounds=(lo,hi),method='bounded')
        p[1]=opt.x
        if fun(p)<best[0]:best=(fun(p),p)
    score,p=best;pred,R,pp=predict(p)
    delta=((pred-obs)/err)**2
    return dict(success=True,train=float(delta[train].sum()),test=float(delta[~train].sum()),prior_penalty=float(pp),train_objective=float(score),
        R_kpc=float(R),V0_kms=float(np.exp(p[1])),sqrt_alpha=float(p[0]) if bose else None,
        params=p.tolist(),num_train=int(train.sum()),num_test=int((~train).sum()),fit_boundary=bool(bose and p[0]<1e-4),
        Vpred=pred.tolist())

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--quality',required=True);p.add_argument('--output',required=True)
    p.add_argument('--library',required=True);p.add_argument('--fine',action='store_true');p.add_argument('--masses',nargs='+',type=float,default=[60.,74.,86.,120.,200.])
    p.add_argument('--sample',choices=['quality','all','pilot'],default='quality');p.add_argument('--prior',action='store_true');p.add_argument('--kinds',nargs='+',default=['one','two','pis','burkert']);p.add_argument('--sigma',type=float,default=.3)
    a=p.parse_args();started=time.monotonic();path=Path(a.library)
    if not path.exists():np.savez_compressed(path,**build(1000 if a.fine else 700,a.fine))
    with np.load(path) as f:lib=Library(f);checks=validate(lib)
    print('library error',max(x['max_relative_V2_error'] for x in checks),flush=True)
    data=read_sparc(a.data);qual=parse(a.quality)['per_galaxy']
    assert len(data)==175 and sum(map(len,data.values()))==3391 and set(data)==set(qual)
    arrays={k:data_arrays(v) for k,v in data.items()}
    if a.sample=='quality':names=sorted(k for k in data if qual[k]['quality']<=2 and qual[k]['inclination_deg']>=30)
    elif a.sample=='pilot':
        names=[]
        for k,v in sorted(arrays.items(),key=lambda x:x[0].lower()):
            r,o,e,vg,vd,vb=v;bary=vg*abs(vg)+.5*vd*vd+.7*vb*vb
            if len(r)>=15 and np.sort(1-bary/o**2)[len(r)//2]>.7 and np.all(o**2>=bary):names.append(k)
        names=names[:20]
    else:names=sorted(data)
    out=dict(sample=a.sample,galaxies=names,count=len(names),input_hashes={k:hashlib.sha256(Path(v).read_bytes()).hexdigest() for k,v in [('data',a.data),('quality',a.quality)]},library_checks=checks,prior=a.prior,sigma=a.sigma,cases=[])
    for kind in a.kinds:
        for m in (a.masses if kind in ('one','two') else [None]):
            t=time.monotonic();rows={}
            for name in names:
                rows[name]=fit(arrays[name],lib,kind,86. if m is None else m,a.sigma,a.prior,qual[name])
            ok=[v for v in rows.values() if v['success']]
            rec=dict(kind=kind,mass_eV=m,success_count=len(ok),train_chi2=sum(v['train'] for v in ok),test_chi2=sum(v['test'] for v in ok),
                objective=sum(v['train_objective'] for v in ok),num_train=sum(v['num_train'] for v in ok),num_test=sum(v['num_test'] for v in ok),
                boundary_count=sum(v['fit_boundary'] for v in ok),elapsed_s=time.monotonic()-t,rows=rows)
            out['cases'].append(rec);out['wall_seconds']=time.monotonic()-started
            Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
            print({k:v for k,v in rec.items() if k!='rows'},flush=True)

if __name__=='__main__':main()
