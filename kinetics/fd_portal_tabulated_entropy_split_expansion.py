"""Finite-volume entropy-split expansion benchmark for the closed FD collision step."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

from fd_number_heat_state_reference import evaluate
from fd_portal_differential_reciprocity_reference import PositivePortal, integrate
from fd_portal_discrete_gradient_reservoir_reference import rule
from fd_portal_tabulated_bath_discrete_entropy import Tabulated
from sm_dark_entropy_exchange_reference import MPL, StandardModel


def hubble(es,ed):
    return np.sqrt((es+ed)/3)/MPL


class Stage:
    def __init__(self,eos,p,nodes):
        self.eos=eos
        self.p=p
        self.nodes=nodes
        self.lo,self.hi=np.exp([eos.lo+2e-14,eos.hi-2e-14])

    def bath(self,u):
        return evaluate(np.exp(u[1]),u[0],nodes=self.nodes)

    def drift_chart(self,a,ds):
        Ts,eta,l=a
        e,pr,s,v=self.eos.state(Ts)
        old=self.bath(a[1:])
        factor=np.exp(-3*ds)
        T=brentq(lambda z:self.eos.state(z)[2]/s-factor,self.lo,self.hi,xtol=4e-15)
        target=np.array([old['n'],old['entropy']])*factor
        def f(u):
            d=self.bath(u)
            return (np.array([d['n'],d['entropy']])-target)/target
        b=root(f,[eta,l-ds],tol=2e-12).x
        if np.max(abs(f(b)))>4e-13:
            raise ArithmeticError('Adiabatic density invariants failed')
        return np.r_[T,b]

    def drift(self,x,a,dt):
        e=self.eos.state(a[0])[0]
        H=hubble(e,self.bath(a[1:])['rho'])
        def f(y):
            m=self.drift_chart(a,y/2)
            return y/dt-hubble(self.eos.state(m[0])[0],self.bath(m[1:])['rho'])
        dx=brentq(f,0.,dt*H*(1.+1e-8),xtol=8e-17)
        b=self.drift_chart(a,dx)
        def work(order):
            t,w=rule(order)
            value=[]
            for u in t:
                r=self.drift_chart(a,dx*u)
                P=self.eos.state(r[0])[1]+self.bath(r[1:])['pressure']
                value.append(3*np.exp(3*(x+dx*u))*P)
            return dx*(w@value)
        old,new=self.energy(x,a),self.energy(x+dx,b)
        w,z=work(5),work(12)
        defect=max(abs(new-old+z),abs(w-z))/old
        if defect>7e-12 or min(w,z)<=0.:
            raise ArithmeticError('Independent pressure-work quadratures do not close the adiabatic energy ledger')
        return x+dx,b,z,defect

    def closed_bath(self,T,u):
        return Tabulated(T,u,self.nodes,self.eos)

    def step(self,x,a,dt):
        w,b,j,k=self.drift(x,a,dt/2)
        c=self.closed_bath(b[0],b[1:])
        u,e,s=c.step(self.p,b[1:],dt)
        p,v,l,m=self.drift(w,np.r_[c.point(u)[3],u],dt/2)
        return p,v,np.exp(3*w)*s,e,j+l,max(k,m)

    def energy(self,x,a):
        return np.exp(3*x)*(self.eos.state(a[0])[0]+self.bath(a[1:])['rho'])

    def entropy(self,x,a):
        return np.exp(3*x)*(self.eos.state(a[0])[2]+self.bath(a[1:])['entropy'])

    def collision(self,Ts,Td,eta):
        f=integrate(self.p,Ts,2*eta,1-Ts/Td)[0]
        return np.array([2.,Ts])*f[:2],f[2]

    def rhs(self,x,a):
        Ts,eta,l=a
        Td=np.exp(l)
        rs,ps,_,cv=self.eos.state(Ts)
        d=self.bath(a[1:])
        rd,pd,n,B,E,C=(d[k] for k in ('rho','pressure','n','dn_deta','chemical_energy','cv_fixed_number'))
        H=hubble(rs,rd)
        k,sigma=self.collision(Ts,Td,eta)
        A=np.array([[B,B*E/Td],[B*E,Td*C+B*E*E/Td]])
        vel=np.linalg.solve(A,-3*H*np.array([n,rd+pd])-k)
        return H,np.r_[(-3*H*(rs+ps)+k[1])/cv,vel],np.exp(3*x)*sigma


def run(path,stage=Stage):
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    ref=stage(eos,q,168)
    secondary=stage(eos,p,112)
    rows=[]
    for T,eta,r in ((.05,.3,1.16),(.3,-.75,1.16),(.6882487043,.9,.86)):
        a=np.array([T,eta,np.log(T*r)])
        iso=ref.closed_bath(T,a[1:])
        y=iso.point(a[1:])[0]
        f=iso.flux(q,a[1:])[0]
        H=ref.rhs(0,a)[0]
        horizon=min(.012/H,.11/np.max(abs(f/y)))
        initial_entropy=ref.entropy(0,a)
        initial_energy=ref.energy(0,a)

        def rhs(t,s):
            g,k,u=ref.rhs(s[0],s[1:4])
            P=eos.state(s[1])[1]+ref.bath(s[2:4])['pressure']
            return horizon*np.r_[g,k,u/initial_entropy,3*g*np.exp(3*s[0])*P/initial_energy]

        s=solve_ivp(rhs,(0.,1.),np.r_[0.,a,0.,0.],method='DOP853',rtol=3e-12,atol=3e-14)
        if not s.success:
            raise ArithmeticError('Independent nonsplit expansion reference missing')
        target=s.y[:,-1]
        if abs(ref.entropy(target[0],target[1:4])/initial_entropy-1.-target[4])>2e-11:
            raise ArithmeticError('Independent continuous entropy budget failed')
        if abs(ref.energy(target[0],target[1:4])/initial_energy-1.+target[5])>2e-11:
            raise ArithmeticError('Independent differential pressure-work check failed')
        mass=np.array([.01,T,1.,1.])
        def follow(g,k):
            x=0.
            z=a.copy()
            area=0.
            work=0.
            work_defects=[]
            scales=[]
            S=g.entropy(x,z)
            E=g.energy(x,z)
            for _ in range(k):
                x,z,u,e,w,d=g.step(x,z,horizon/k)
                work+=w
                work_defects.append(d)
                area+=u
                scales.append(e)
            closure=abs(g.entropy(x,z)/S-1.-area/S)
            error=np.max(abs((np.r_[x,z]-target[:4])/mass))
            energy_closure=abs(g.energy(x,z)/E-1.+work/E)
            if closure>1e-11 or energy_closure>1e-11 or area<=0.:
                raise ArithmeticError('Discrete expansion/comoving-entropy ledger failed')
            return dict(steps=k,log_a=x,state=z.tolist(),relative_comoving_entropy_added=area/S,
                        scaled_independent_error=error,entropy_ledger_defect=closure,max_collision_entropy_defect=max(scales),
                        normalized_pressure_work=work/E,comoving_energy_ledger_defect=energy_closure,max_drift_energy_defect=max(work_defects))
        a1=[follow(ref,k) for k in (2,4,8)]
        z=follow(secondary,4)
        d=np.array([x['scaled_independent_error'] for x in a1])
        ratios=d[:-1]/d[1:]
        independent=np.max(abs(np.asarray(a1[1]['state'])-z['state']))
        if min(ratios)<3.25 or max(ratios)>4.85 or independent>2e-7:
            raise ArithmeticError(('Short expansion time convergences did not resolve',T,ratios.tolist(),independent))
        rows.append(dict(initial_state=a.tolist(),horizon=horizon,independent_terminal_coordinates=target.tolist(),
                         levels=a1,error_halving_ratios=ratios.tolist(),independent_mesh_force_defect=independent))
    return dict(rows=rows,table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Short homogeneous expansion, tabulated visible EOS and restricted timelike vacuum FD collisions only. No cosmological duration or abundance inference.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(run(p.parse_args().eos),indent=2,allow_nan=False))
