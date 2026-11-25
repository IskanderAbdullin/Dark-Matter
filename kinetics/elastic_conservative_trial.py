"""Conservative nodal elastic 2<->2 trial operator for a frozen scale factor.

Each unordered incoming representative pair is integrated over the ordered
outgoing energy interval with lumped nonnegative hat weights. The complementary
outgoing kinetic energy is projected onto its two neighbouring representatives.
The reverse contribution uses the same weights and geometric f / (1+f).
Every event conserves discrete particle number and energy, and each discrete
Bose-Einstein state has detailed balance. This is a trial discretization; it
drops out-of-range phase space and approximates angles and occupations.
No condensate mode, 2<->4, expansion, production, or dynamical Hartree mass.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import solve_ivp
from quartic_consistency import shell_nodes,lambda_from_sigma,gap
from quantum_elastic_moments import kinetic,momentum

class Elastic:
    def __init__(self,edges,m):
        self.edges=np.array(edges);self.m=m
        l,h=self.edges[:-1],self.edges[1:];self.n=len(l)
        self.vol=(h**3-l**3)/(6*np.pi**2)
        q,w=shell_nodes(self.edges,np.ones(self.n),8)
        q=q.reshape(self.n,8);w=w.reshape(self.n,8)
        self.K=(w*kinetic(q,m)).sum(axis=1)/self.vol
        self.p=momentum(self.K,m)
        # Nodal weights in kinetic-energy integration, one per original cell.
        self.dK=np.diff(kinetic(self.edges,m))
        I,J=np.triu_indices(self.n)
        S=self.K[I]+self.K[J]
        ii,jj,kk,ll,ww,aa=[],[],[],[],[],[]
        for i,j,s in zip(I,J,S):
            # Ordered outgoing K3 <= K4. Integrate positive linear hats
            # over the allowed energy interval inside represented energies.
            lo=max(self.K[0],s-self.K[-1]);hi=min(s/2,self.K[-1])
            if hi<=lo:continue
            L=np.maximum(self.K[:-1],lo);H=np.minimum(self.K[1:],hi)
            width=np.maximum(H-L,0.);c=(L+H)/2
            t=np.clip((c-self.K[:-1])/np.diff(self.K),0.,1.)
            hats=np.r_[width*(1-t),0.]+np.r_[0.,width*t]
            k=np.flatnonzero(hats>0)
            k4=s-self.K[k]
            valid=(k4>=self.K[0])&(k4<=self.K[-1])
            k=k[valid];k4=k4[valid]
            if len(k)==0:continue
            a=np.clip(np.searchsorted(self.K,k4)-1,0,self.n-2)
            v=np.clip((k4-self.K[a])/(self.K[a+1]-self.K[a]),0.,1.)
            p4=momentum(k4,m);p3=self.p[k]
            angle=np.maximum(np.minimum(self.p[i]+self.p[j],p3+p4)
                -np.maximum(abs(self.p[i]-self.p[j]),abs(p3-p4)),0.)
            weight=2*(1 if i==j else 2)*self.dK[i]*self.dK[j]*hats[k]*angle
            use=weight>0
            k,a,v,weight=k[use],a[use],v[use],weight[use]
            ii.extend(np.full(len(k),i));jj.extend(np.full(len(k),j));kk.extend(k);ll.extend(a)
            ww.extend(v);aa.extend(weight)
        self.i,self.j,self.k,self.l=[np.asarray(a,dtype=np.int32) for a in (ii,jj,kk,ll)]
        self.w=np.array(ww);self.rate=np.array(aa)
        self.indices=np.r_[self.i,self.j,self.k,self.l,self.l+1]
        self.sigma=np.r_[-np.ones(2*len(self.i)),np.ones(len(self.i)),1-self.w,self.w]
        self.event_energy_error=float(np.max(np.abs(self.K[self.i]+self.K[self.j]-self.K[self.k]
                                  -(1-self.w)*self.K[self.l]-self.w*self.K[self.l+1])))
    def flux(self,f):
        f=np.maximum(f,0.)
        # Robust exact cancellation of the quartic term when f>0.
        ff=np.maximum(f,1e-100);lf=np.log(ff);r=np.log1p(1/ff)
        i,j,k,l,w=self.i,self.j,self.k,self.l,self.w
        logprod=lf[i]+lf[j]+lf[k]+(1-w)*lf[l]+w*lf[l+1]
        rin=r[i]+r[j];rout=r[k]+(1-w)*r[l]+w*r[l+1]
        delta=rout-rin
        sign=np.sign(delta)
        # Always exponentiate the smaller, nonpositive log difference.
        event=sign*np.exp(logprod+np.maximum(rin,rout))*(-np.expm1(-np.abs(delta)))
        # Zero occupations are represented only by a tiny computational floor;
        # the resulting error is negligible for our fixed positive test PSD.
        return self.rate*event
    def rhs(self,f):
        event=self.flux(f)
        return np.bincount(self.indices,weights=self.sigma*np.tile(event,5),minlength=self.n)/self.vol
    def positive_local_step(self,f,scale,frac=.01):
        d=self.rhs(f)*scale
        neg=d<0
        if np.any(neg):
            limit=np.min(f[neg]/(-d[neg]))
        else:limit=np.inf
        dt=min(frac*limit,.01)
        if not (dt>0 and np.isfinite(dt)):
            return dict(success=False,reason="No positive finite step")
        g=f+dt*d
        assert np.min(g)>=0
        N=self.vol@f;E=(self.vol*self.K)@f
        return dict(success=True,dt_H=float(dt),min_f=float(np.min(g)),
                    relative_number_change=float(self.vol@g/N-1),
                    relative_energy_change=float((self.vol*self.K)@g/E-1),
                    relative_q2_change=float((self.vol*self.p**2)@g/((self.vol*self.p**2)@f)-1))

    def diagnostics(self,f):
        d=self.rhs(f);N=self.vol@f;M2=self.vol@(self.p**2*f)
        return dict(events=len(self.i),event_energy_error=self.event_energy_error,number=float(N),
            q2_moment=float(M2),number_residual=float(abs(self.vol@d)/max(self.vol@abs(d),1e-300)),
            energy_residual=float(abs((self.vol*self.K)@d)/max((self.vol*self.K)@abs(d),1e-300)),
            q2_raw=float((self.vol*self.p**2)@d),q2_relative_raw=float((self.vol*self.p**2)@d/M2),
            entropy_production_raw=float((self.vol*np.log1p(1/np.maximum(f,1e-100)))@d))

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--sigma',type=float,default=1e-7)
    p.add_argument('--output',required=True);p.add_argument('--evolve',action='store_true')
    a=p.parse_args();row=json.loads(Path(a.input).read_text())[0];h=row['final']
    e=np.array(row['frozen_cell_spectrum']['edges']);f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    q,w=shell_nodes(e,f,8);lam=float(lambda_from_sigma(a.sigma));mu=gap(q,w,row['scalar_ratio']*h['x'],lam)
    t=time.monotonic();m=Elastic(e,mu)
    Href=1.66*np.sqrt(row['g_bath'])*row['parent_GeV']**2/1.220890e19
    Hf=np.sqrt(np.pi**2/30*row['g_bath']/h['total_rho_Tref4'])
    # This event sum uses both forward and reverse orientations, so one half
    # of the forward-only moment prefactor. Includes all incoming orderings.
    scale=lam**2/(2048*np.pi**5)*row['parent_GeV']*h['x']/Href*Hf
    d=m.diagnostics(f)
    eq=1/np.expm1(m.K+.1)
    eqrate=m.rhs(eq);eq_abs=np.bincount(m.indices,weights=np.abs(m.sigma)*np.tile(np.abs(m.rate),5),minlength=m.n)
    diag=dict(n=m.n,mu=mu,lambda_quartic=lam,lambda_f_peak=lam*f.max(),sigma_over_m_cm2_g=a.sigma,
        normalization=scale,relative_q2_change_per_H=d['q2_relative_raw']*scale,
        detailed_balance_weighted_abs=float(m.vol@abs(eqrate)),rate_measure_sum=float(m.rate.sum()),**d)
    print(json.dumps(diag),flush=True)
    result=dict(input=a.input,diagnostics=diag,construction_seconds=time.monotonic()-t,positive_local_step=m.positive_local_step(f,scale))
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    if not a.evolve:return
    # Frozen geometry test only. Scale occupations to reduce solver conditioning.
    norm=np.maximum(f,1.)
    dt=.03/max(abs(diag['relative_q2_change_per_H']),1e-20)
    started=time.monotonic()
    def limited_rhs(t,y):
        if time.monotonic()-started>90:raise TimeoutError("90 s explicit relaxation trial budget")
        return m.rhs(y*norm)*scale/norm
    try:
        sol=solve_ivp(limited_rhs,(0.,dt),f/norm,method='RK23',
                      atol=1e-11,rtol=1e-5,max_step=dt/25.)
    except TimeoutError as err:
        result.update(evolution_success=False,message=str(err))
        Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result),flush=True)
        return
    result.update(evolution_success=bool(sol.success),message=sol.message,dt_in_Hubble_times=dt,
        nfev=sol.nfev,minimum_occupation=float(np.min(sol.y*norm[:,None])),
        final_diagnostics=m.diagnostics(sol.y[:,-1]*norm))
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
