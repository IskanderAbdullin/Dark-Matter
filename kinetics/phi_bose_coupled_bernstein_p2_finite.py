"""Physical-capacity conservative coupled three-species Yukawa finite-step prototype.

Shared positive boosted on-shell Yukawa event measure and nonnegative
Bernstein-P2 energy-cell basis reproducing constants and physical energies.
Mode moments integrate physical g*p²dp/(2π²) density of states.
Fermi/Bose occupations are reconstructed from local entropy-dual eta.

Moment-explicit update + nonlinear local positive inversion, with adaptive
a-posteriori entropy/conservation rejection; NOT an unconditional positivity
or entropy theorem, cosmological time evolution or condensate solver.
"""
from __future__ import annotations
import math
import numpy as np
from dataclasses import dataclass
from numpy.polynomial.legendre import leggauss
from scipy.special import expit
from scipy.sparse import csr_matrix
from kinetics.phi_bose_physical_cutcell_weak import (
    M1,M2,MPHI,G1,KSTAR,ESTAR,WIDTH,shell,scalar_edges,parent_splits)
from kinetics.phi_bose_invariant_parent_support import global_support_floor


def bernstein(x):
    x=np.asarray(x)
    return np.stack(((1-x)**2,2*x*(1-x),x*x),axis=-1)


@dataclass
class Species:
    name:str
    mass:float
    g:int
    kind:str
    edges:np.ndarray
    weights:np.ndarray
    H:np.ndarray
    energy:np.ndarray
    nodal_energies:np.ndarray

    @property
    def ncells(self):return len(self.edges)-1
    @property
    def nvars(self):return 3*self.ncells

    def evaluate(self,eta,with_G=False,with_entropy=False):
        eta=np.reshape(eta,(self.ncells,3))
        local_eta=np.einsum('cqa,ca->cq',self.H,eta)
        if self.kind=='fermion':
            occ=expit(-local_eta)
            var=occ*(1-occ)
        else:
            if np.min(local_eta)<=0:raise ValueError('Bose inverse fugacity <=0')
            occ=1/np.expm1(local_eta)
            var=occ*(1+occ)
        M=np.einsum('cq,cqa,cq->ca',self.weights,self.H,occ)
        ret={'M':M.ravel()}
        if with_G:
            G=np.einsum('cq,cqa,cqb,cq->cab',self.weights,self.H,self.H,var)
            ret['G']=G
        if with_entropy:
            if self.kind=='fermion':
                sdens=local_eta*occ+np.logaddexp(0,-local_eta)
            else:
                sdens=local_eta*occ-np.log(-np.expm1(-local_eta))
            ret['S']=float(np.sum(self.weights*sdens))
        return ret

    def moment_difference(self,eta,eta0):
        # Direct occupation difference avoids subtracting large equilibrium moments.
        shape=(self.ncells,3)
        new=np.einsum('cqa,ca->cq',self.H,np.reshape(eta,shape))
        old=np.einsum('cqa,ca->cq',self.H,np.reshape(eta0,shape))
        if self.kind=='fermion':diff=expit(-new)-expit(-old)
        else:
            if np.min(new)<=0:raise ValueError('Bose eta <=0')
            diff=1/np.expm1(new)-1/np.expm1(old)
        return np.einsum('cq,cqa,cq->ca',self.weights,self.H,diff).ravel()


def build_species(name,mass,g,kind,edges,nq=28):
    edges=np.asarray(edges,float)
    assert np.all(np.diff(edges)>0) and edges[0]>=mass and g>0
    z,w=leggauss(nq)
    pl=np.sqrt(np.maximum(0,(edges[:-1]-mass)*(edges[:-1]+mass)))
    pu=np.sqrt(np.maximum(0,(edges[1:]-mass)*(edges[1:]+mass)))
    p=pl[:,None]+(pu-pl)[:,None]*(z[None,:]+1)/2
    dp=(pu-pl)[:,None]*w[None,:]/2
    E=np.hypot(p,mass)
    h=bernstein((E-edges[:-1,None])/(edges[1:,None]-edges[:-1,None]))
    W=g*dp*p*p/(2*math.pi**2)
    nodes=np.stack((edges[:-1],(edges[1:]+edges[:-1])/2,edges[1:]),axis=1)
    assert np.all(h>=-1e-12)
    assert np.max(abs(np.sum(h,axis=2)-1))<2e-13
    assert np.max(abs(np.sum(h*nodes[:,None,:],axis=2)-E))<2e-12
    return Species(name,mass,g,kind,edges,W,h,E,nodes.ravel())


class Coupled:
    def __init__(self,P=2.,N1=4,N2=4,Nphi=8,nq=28,nrad=8,nangle=12):
        if not (0<P<.99*M1*KSTAR/MPHI):
            raise ValueError('Prototype parent cutoff below scalar-at-rest boost required; rederive cut-cell partitions otherwise')
        self.P=P
        emax=math.hypot(P,M1)
        e2star=(M1*M1+M2*M2-MPHI*MPHI)/(2*M1)
        e2max=(emax*e2star+P*KSTAR)/M1
        fedges=[np.linspace(M1,emax,N1+1),np.linspace(M2,e2max,N2+1)]
        kmin=global_support_floor(MPHI,P)
        kmax=shell(P)[2]
        phi_edges=np.unique(np.r_[scalar_edges(P,Nphi),kmin,kmax])
        self.species=[build_species('chi1',M1,G1,'fermion',fedges[0],nq),
                      build_species('chi2',M2,G1,'fermion',fedges[1],nq),
                      build_species('phi',MPHI,1,'boson',phi_edges,nq)]
        self.offsets=np.cumsum([0]+[s.nvars for s in self.species])
        self.phi_support=(kmin,kmax)
        self._make_events(nrad,nangle)

    def _make_events(self,nrad,nangle):
        e1_edges=self.species[0].edges
        e2_edges=self.species[1].edges
        phi_edges=self.species[2].edges
        cuts=parent_splits(self.P,phi_edges)
        x,w=leggauss(nrad);z,v=leggauss(nangle)
        ids=[[],[],[]]; hats=[[],[],[]]
        weights=[];evalues=[[],[],[]]
        for p0,p1 in zip(cuts[:-1],cuts[1:]):
            for xx,ww in zip(x,w):
                p=p0+(p1-p0)*(xx+1)/2
                dp=(p1-p0)*ww/2
                E1=math.hypot(p,M1)
                lo,hi,B=shell(p)[1:4]
                pre=G1*p*p*dp/(2*math.pi**2)*WIDTH*M1/E1
                for j,(pl,pu) in enumerate(zip(phi_edges[:-1],phi_edges[1:])):
                    l=max(pl,lo);h=min(pu,hi)
                    if h<=l:continue
                    lg=math.log(h/l)
                    K=l*np.exp((z+1)*lg/2)
                    dk=v*K*lg/2
                    E2=E1-K
                    for kk,dkk,e2 in zip(K,dk,E2):
                        k=float(kk);E2v=float(e2)
                        idx1=np.clip(np.searchsorted(e1_edges,E1,side='right')-1,0,len(e1_edges)-2)
                        idx2=np.clip(np.searchsorted(e2_edges,E2v,side='right')-1,0,len(e2_edges)-2)
                        ix=[int(idx1),int(idx2),int(j)]
                        E=[E1,E2v,k]
                        for a in range(3):
                            ed=self.species[a].edges
                            f=(E[a]-ed[ix[a]])/(ed[ix[a]+1]-ed[ix[a]])
                            if f < -5e-10 or f >1+5e-10:
                                raise ArithmeticError(('event outside species energy bins',a,f,E[a]))
                            ids[a].append(ix[a])
                            hats[a].append(bernstein(float(np.clip(f,0,1))).ravel())
                            evalues[a].append(E[a])
                        weights.append(pre*dkk/(2*B))
        self.W=np.array(weights)
        self.E=[np.array(t) for t in evalues]
        n=len(self.W)
        self.A=[]
        for s,i0,H in zip(self.species,ids,hats):
            H=np.array(H)
            row=np.repeat(np.arange(n),3)
            col=(3*np.asarray(i0)[:,None]+np.arange(3)[None,:]).ravel()
            self.A.append(csr_matrix((H.ravel(),(row,col)),shape=(n,s.nvars)))
        assert np.min(self.W)>0 and n>100
        self.event_count=n
        for spec,a,E in zip(self.species,self.A,self.E):
            assert np.max(abs(a@spec.nodal_energies-E))<3e-12

    def initial_etas(self,T=(.9,1.,1.),mu=(0.,0.,0.)):
        return [np.array((s.nodal_energies-mu[i])/T[i]) for i,s in enumerate(self.species)]

    def source(self,etas):
        eta1,eta2,etaphi=[a@eta for a,eta in zip(self.A,etas)]
        if min(etaphi)<=0:raise ValueError('nonpositive boson eta')
        f1=expit(-eta1)
        f2=expit(-eta2)
        b=1/np.expm1(etaphi)
        F=f1*(1-f2)*(1+b)
        R=(1-f1)*f2*b
        J=F-R
        rate=self.W*J
        C=[-self.A[0].T@rate,self.A[1].T@rate,self.A[2].T@rate]
        entropy=float(np.dot(rate,-eta1+eta2+etaphi))
        return dict(C=C,H=entropy,Rphi=float(sum(rate)),
                    Qphi=float(np.dot(rate,self.E[2])),
                    Q1=float(-np.dot(rate,self.E[0])),
                    Q2=float(np.dot(rate,self.E[1])))

    def integrate_spectrum(self,etas):
        vals=[s.evaluate(eta,with_entropy=True) for s,eta in zip(self.species,etas)]
        return dict(N=[float(sum(v['M'])) for v in vals],
                    Q=[float(np.dot(s.nodal_energies,v['M'])) for s,v in zip(self.species,vals)],
                    S=sum(v['S'] for v in vals))

    def step(self,etas,dt,maxiter=16):
        """Explicit moment-Euler step with positive local inverse moment solve.

        This is a conditional numerical construction, not an unconditional
        entropy/positivity theorem. Some large steps fail and must be rejected.
        """
        initial=self.source(etas)
        next_eta=[]; stats=[]
        for species,eta0,C in zip(self.species,etas,initial['C']):
            eta0=np.asarray(eta0).reshape((-1,3))
            C=np.asarray(C).reshape((-1,3))
            current=eta0.copy()
            initial_M=species.evaluate(eta0.ravel())['M'].reshape((-1,3))
            maxres=0
            for j in range(species.ncells):
                eold=eta0[j]
                rhs=dt*C[j]
                resolvable=max(1e-23,32*np.finfo(float).eps*max(abs(initial_M[j])))
                if max(abs(rhs))<resolvable:
                    continue
                e=eold.copy()
                for it in range(maxiter):
                    trial=current.copy();trial[j]=e
                    dM=species.moment_difference(trial.ravel(),eta0.ravel()).reshape((-1,3))[j]
                    residual=dM-rhs
                    relative=max(abs(residual))/max(max(abs(rhs)),1e-23)
                    if relative<5e-8:break
                    G=species.evaluate(trial.ravel(),with_G=True)['G'][j]
                    de=np.linalg.solve(G,residual)
                    alpha=1.
                    while species.kind=='boson' and min(e+alpha*de)<=0 and alpha>1e-8:
                        alpha/=2
                    if alpha<=1e-8:
                        raise ValueError(('Positive-Bose Newton line search failed; no nonrepresentability proof',species.name,j,dt))
                    e=e+alpha*de
                else:
                    raise ArithmeticError(('non-converged cell moment inversion',species.name,j,relative))
                current[j]=e
                maxres=max(maxres,relative)
            next_eta.append(current.ravel())
            stats.append(maxres)
        before=self.integrate_spectrum(etas)
        after=self.integrate_spectrum(next_eta)
        after_source=self.source(next_eta)
        energy_defect=sum(a-b for a,b in zip(after['Q'],before['Q']))
        number_residuals=((after['N'][0]-before['N'][0])+(after['N'][1]-before['N'][1]),
                          (after['N'][0]-before['N'][0])+(after['N'][2]-before['N'][2]))
        return dict(etas=next_eta,source_old=initial,source_new=after_source,
                    before=before,after=after,energy_defect=energy_defect,
                    number_residuals=number_residuals,
                    entropy_increase=after['S']-before['S'],
                    moment_solve_residual=stats,
                    min_bose_eta=float(min(next_eta[2])))


def conservative_adaptive_step(model,eta,requested_dt,max_halvings=16):
    """Reject until an admissible numerical moment step is found; a posteriori only."""
    assert requested_dt>0
    dt=float(requested_dt)
    rejections=[]
    for attempt in range(max_halvings+1):
        try:
            v=model.step(eta,dt)
            fluxscale=max(1e-25,dt*abs(v['source_old']['Rphi']))
            defects=max(abs(v['energy_defect']),*(abs(x) for x in v['number_residuals']))
            if v['entropy_increase'] < -1e-12*fluxscale:
                raise ValueError('Numerical physical entropy decreases')
            if defects>1e-8*fluxscale:
                raise ValueError('Discrete physical charge/energy moment defect too large')
            if v['min_bose_eta']<=0:
                raise ValueError('Nonpositive Bose entropy dual node')
            v['accepted_dt']=dt
            v['rejections']=rejections
            return v
        except (ValueError,ArithmeticError,OverflowError,np.linalg.LinAlgError) as error:
            rejections.append({'dt':dt,'reason':str(error)})
            dt/=2.
    raise ArithmeticError('No positive entropy/conserving numerical step within backtracking budget')


def integrate_to(model,initial_etas,end_time,dt_max):
    assert end_time>0 and dt_max>0
    eta=initial_etas
    t=0.;steps=[]
    while t<end_time-1e-9:
        v=conservative_adaptive_step(model,eta,min(dt_max,end_time-t))
        eta=v['etas'];t+=v['accepted_dt'];steps.append(v)
        if len(steps)>1000:raise ArithmeticError('Excessive numerical backtracking')
    return dict(etas=eta,elapsed=t,accepted_steps=len(steps),
                rejected_attempts=sum(len(v['rejections']) for v in steps),
                last_step=steps[-1],history=steps)


def audit():
    m=Coupled()
    eta=m.initial_etas()
    first=m.source(eta)
    from kinetics.phi_bose_physical_cutcell_weak import analytic_weak
    ref=analytic_weak(2.)
    print('events',m.event_count,'basis',[s.nvars for s in m.species],
          'weak',{k:first[k] for k in ('Rphi','Qphi','Q1','Q2','H')},
          'relative errors',{k:first[k]/ref[k]-1 for k in ('Rphi','Qphi','Q1','Q2')})
    for ix in range(3):
        result=conservative_adaptive_step(m,eta,100000.)
        eta=result['etas']
        print('step',ix,'dt',result['accepted_dt'],
              'entropy increase',result['entropy_increase'],
              'Q defect',result['energy_defect'],
              'charge defects',result['number_residuals'],
              'min Bose eta',result['min_bose_eta'],
              'moment inversion errors',result['moment_solve_residual'])
    return m


if __name__=='__main__':
    audit()
