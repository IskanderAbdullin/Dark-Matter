"""Implicit, fixed-geometry evolution of the conservative elastic trial.

Fractional stoichiometry and positive rate weights from Elastic. Differentiate
the net Bose flux stably in log occupation; use the chain rule to scaled f.
Internal Newton iterates can leave f>=0; regularize those evaluations only,
and stop on significant negative accepted occupations. No renormalization.
This is not yet an expanding or chemically reacting calculation.
"""
import argparse,json,time,os
from pathlib import Path
import numpy as np
from scipy.integrate import BDF,Radau
from elastic_conservative_trial import Elastic
from quartic_consistency import shell_nodes,gap,lambda_from_sigma

class ImplicitElastic(Elastic):
    def __init__(self,e,m):
        super().__init__(e,m)
        self.idx=np.stack([self.i,self.j,self.k,self.l,self.l+1])
        self.coef=np.stack([np.ones(len(self.i))]*3+[1-self.w,self.w])
        self.stoich=self.coef.copy();self.stoich[:2]*=-1
        self.vol_inv=1/self.vol
        self.floor=1e-90
        self.compiled=os.environ.get("DM_ELASTIC_COMPILED")=="1"
        if self.compiled:
            from elastic_acceleration import compiled_rates
            self.compiled_fn=compiled_rates
    def rates_jac(self,f,need_jac):
        if self.compiled:
            R,J=self.compiled_fn(f,self.idx,self.coef,self.stoich,self.rate,self.vol_inv,self.floor,need_jac)
            return J if need_jac else R
        return self.rates_jac_numpy(f,need_jac)
    def rates_jac_numpy(self,f,need_jac):
        ff=np.maximum(f,self.floor)
        logf=np.log(ff);r=np.log1p(1/ff)
        i,j,k,l,h=self.idx
        prod=(logf[i]+logf[j]+logf[k]+self.coef[3]*logf[l]+self.coef[4]*logf[h])
        rin=r[i]+r[j];rout=r[k]+self.coef[3]*r[l]+self.coef[4]*r[h]
        diff=rout-rin
        F=self.rate*np.sign(diff)*np.exp(prod+np.maximum(rin,rout))*(-np.expm1(-abs(diff)))
        rhs=np.bincount(self.idx.ravel(),weights=(self.stoich*F).ravel(),minlength=self.n)*self.vol_inv
        if not need_jac:return rhs
        # dF/dlog(f_s) = c_s F +/- c_s A exp(r_side)/(1+f_s).
        # Evaluate derivative in stable grouped form for small f; the
        # near-empty incoming/outgoing derivative of the appropriate term
        # cancels analytically, not as two large floating point terms.
        fw=self.rate*np.exp(prod+rout);bw=self.rate*np.exp(prod+rin)
        G=ff[self.idx]
        D=np.empty_like(self.coef)
        # incoming factors f_i (forward) and 1+f_i (backward)
        D[:2]=self.coef[:2]*(F[None,:]/G[:2]+bw[None,:]/(G[:2]*(1+G[:2])))
        # outgoing factors 1+f_k (forward) and f_k (backward)
        D[2:]=self.coef[2:]*(F[None,:]/G[2:]-fw[None,:]/(G[2:]*(1+G[2:])))
        # Algebraically equivalent direct derivative near zero (avoid
        # cancellation between F and fw/bw when one side dominates).
        near=G<1e-3
        D[:2]=np.where(near[:2],self.coef[:2]*(fw[None,:]/G[:2]-bw[None,:]/(1+G[:2])),D[:2])
        D[2:]=np.where(near[2:],self.coef[2:]*(fw[None,:]/(1+G[2:])-bw[None,:]/G[2:]),D[2:])
        D=np.where(f[self.idx]>self.floor,D,0.)
        J=np.zeros((self.n,self.n))
        for a in range(5):
            for b in range(5):
                J+=np.bincount(self.idx[a]*self.n+self.idx[b],
                    weights=self.stoich[a]*D[b],minlength=self.n*self.n).reshape(self.n,self.n)
        J*=self.vol_inv[:,None]
        return J
    def validate_jac(self):
        f=1/np.expm1(self.K+.3)
        f=f*(1+.2*np.sin(np.arange(self.n)))
        J=self.rates_jac(f,True)
        d=np.random.default_rng(192).normal(size=self.n)*f
        eps=1e-6
        fd=(self.rates_jac(f+eps*d,False)-self.rates_jac(f-eps*d,False))/(2*eps)
        analytic=J@d
        err=np.linalg.norm(fd-analytic)/max(np.linalg.norm(analytic),1e-300)
        assert err<2e-5,err
        eq=1/np.expm1(self.K+.3)
        res=self.rates_jac(eq,False)
        old=self.rhs(f);new=self.rates_jac(f,False)
        err2=np.linalg.norm(old-new)/max(np.linalg.norm(old),1e-300)
        assert err2<1e-12,err2
        return dict(jac_relative_error=float(err),old_rhs_relative_error=float(err2),equilibrium_abs=float(self.vol@abs(res)))

def merge_grid(edges,fp,n,qmin=1e-5,qmax=30.):
    e=np.asarray(edges);f=np.asarray(fp)
    # Retain a subset of original faces; exact number rebin on each interval.
    a=max(0,int(np.searchsorted(e,qmin,side='right')-1))
    b=min(len(f),int(np.searchsorted(e,qmax,side='left')))
    ix=np.unique(np.rint(np.linspace(a,b,min(n,b-a)+1)).astype(int))
    N=np.maximum(f,0)*np.diff(e**3)/(6*np.pi**2)
    g=np.array([N[l:h].sum() for l,h in zip(ix[:-1],ix[1:])])/ (np.diff(e[ix]**3)/(6*np.pi**2))
    return e[ix],g

def evolve(e,f,m,factor,end,rtol=1e-5,wall=180.,method='BDF'):
    start=time.monotonic();model=ImplicitElastic(e,m);check=model.validate_jac()
    v=model.vol;E=model.K+m;N0=float(v@f);E0=float((v*E)@f);M20=float((v*model.p**2)@f)
    # A floor in initially empty cells has a separately measured negligible
    # number and energy cost. Keep it far below integration tolerances.
    ff=np.maximum(f,1e-38)
    seed_number=float(v@(ff-f));seed_energy=float((v*E)@(ff-f))
    norm=np.maximum(ff,1.)
    scale=factor/(2048*np.pi**5)
    def rhs(t,y):return model.rates_jac(y*norm,False)*scale/norm
    def jac(t,y):return model.rates_jac(y*norm,True)*scale*norm[None,:]/norm[:,None]
    cls={'BDF':BDF,'Radau':Radau}[method]
    sol=cls(rhs,0.,ff/norm,end,jac=jac,rtol=rtol,atol=np.full(len(ff),1e-12),max_step=end/20.)
    hist=[];steps=0;bad=0;message='';worst_negative=0.;last=start
    while sol.status=='running':
        tprev=sol.t;yprev=sol.y.copy()
        msg=sol.step();steps+=1
        g=sol.y*norm
        negnum=float(v@np.maximum(-g,0)/N0);negene=float((v*E)@np.maximum(-g,0)/E0)
        worst_negative=max(worst_negative,negnum,negene)
        if negnum>1e-8 or negene>1e-8:
            sol.t=tprev;sol.y=yprev;bad+=1;message='significant negative accepted occupation';break
        if sol.status=='failed':message=str(msg);break
        if steps%10==0 or sol.status=='finished':
            h=dict(t_H=float(sol.t),number_relative_error=float(v@g/N0-1),energy_relative_error=float((v*E)@g/E0-1),
                q2_relative_change=float((v*model.p**2)@g/M20-1),minimum_f=float(g.min()),negative_number_fraction=negnum,
                steps=steps,nfev=sol.nfev,njev=sol.njev,elapsed=time.monotonic()-start)
            hist.append(h)
            if time.monotonic()-last>12 or sol.status=='finished':
                print('progress',json.dumps(h),flush=True);last=time.monotonic()
        if time.monotonic()-start>wall:message='wall budget';break
    g=sol.y*norm
    result=dict(compiled=model.compiled,success=sol.status=='finished' and not bad,message=message or sol.status,n=model.n,method=method,rtol=rtol,
        requested_t_H=end,achieved_t_H=float(sol.t),steps=steps,nfev=sol.nfev,njev=sol.njev,wall_seconds=time.monotonic()-start,
        seed_number_fraction=seed_number/N0,seed_energy_fraction=seed_energy/E0,maximum_negative_fraction=worst_negative,
        initial_number=N0,initial_energy=E0,initial_q2=M20,initial_relative_q2_derivative=float((v*model.p**2)@rhs(0,ff/norm).__mul__(norm)/M20),
        number_relative_error=float(v@g/N0-1),energy_relative_error=float((v*E)@g/E0-1),
        q2_relative_change=float((v*model.p**2)@g/M20-1),checks=check,history=hist,edges=e.tolist(),initial_f=f.tolist(),final_f=g.tolist())
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input',required=True);p.add_argument('--sigma',type=float,default=1e-7)
    p.add_argument('--n',type=int,default=80);p.add_argument('--qmin',type=float,default=1e-5);p.add_argument('--qmax',type=float,default=30.)
    p.add_argument('--end',type=float,default=.01);p.add_argument('--rtol',type=float,default=1e-5);p.add_argument('--wall',type=float,default=180.)
    p.add_argument('--method',choices=['BDF','Radau'],default='BDF');p.add_argument('--output',required=True)
    a=p.parse_args()
    row=json.loads(Path(a.input).read_text())[0];h=row['final']
    e=np.asarray(row['frozen_cell_spectrum']['edges']);f=np.asarray(row['frozen_cell_spectrum']['f_phi'])
    q,w=shell_nodes(e,f,8);lam=float(lambda_from_sigma(a.sigma));m=gap(q,w,row['scalar_ratio']*h['x'],lam)
    fac=lam**2*row['parent_GeV']*h['x']/(1.66*np.sqrt(row['g_bath'])*row['parent_GeV']**2/1.220890e19)*np.sqrt(np.pi**2/30*row['g_bath']/h['total_rho_Tref4'])
    e2,f2=merge_grid(e,f,a.n,a.qmin,a.qmax)
    r=evolve(e2,f2,m,fac,a.end,a.rtol,a.wall,a.method)
    r.update(input=a.input,sigma=a.sigma,mass=m,lambda_quartic=lam,
        original_number=float(w.sum()),original_energy=float(w@np.hypot(q,m)),original_q2=float(w@(q*q)),
        remesh_number_error=r['initial_number']/w.sum()-1,remesh_energy_error=r['initial_energy']/(w@np.hypot(q,m))-1,
        remesh_q2_error=r['initial_q2']/(w@(q*q))-1)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print('FINAL',json.dumps({k:v for k,v in r.items() if k not in ('history','edges','initial_f','final_f')}),flush=True)
if __name__=='__main__':main()
