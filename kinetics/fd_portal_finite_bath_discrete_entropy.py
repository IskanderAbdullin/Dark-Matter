"""Energy-conserving, entropy-discrete FD collision step with a finite ideal reservoir."""
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from fd_number_heat_state_reference import evaluate
from fd_portal_differential_reciprocity_reference import PositivePortal, integrate
from fd_portal_discrete_gradient_reservoir_reference import rule

C_VIS=10.75*np.pi**2/30.


class Bath:
    def __init__(self,T,u,dn):
        self.dn=dn
        self.E=C_VIS*T**4+self.dark(u)[0][1]

    def dark(self,u):
        eta,l=u
        T=np.exp(l)
        d=evaluate(T,eta,nodes=self.dn)
        n,r,s,B,E,C=(d[k] for k in ('n','rho','entropy','dn_deta','chemical_energy','cv_fixed_number'))
        A=np.array([[B,B*E/T],[B*E,T*C+B*E*E/T]])
        return np.array([n,r]),s,A

    def point(self,u):
        y,s,A=self.dark(u)
        r=self.E-y[1]
        if not r>0.:
            raise ArithmeticError('No energy left in the visible bath')
        T=(r/C_VIS)**.25
        g=np.array([u[0],1./T-np.exp(-u[1])])
        return y,-s-4*C_VIS*T**3/3,g,T,A

    def invert(self,y,guess):
        def fit(u):
            z,_,A=self.dark(u)
            return (z-y)/y,A/y[:,None]
        u=root(lambda z:fit(z)[0],guess,jac=lambda z:fit(z)[1],tol=8e-13).x
        if np.max(abs(fit(u)[0]))>3e-13:
            raise ArithmeticError('Physical moment recovery failed')
        return u

    def path(self,a,b,y,z,n):
        t,w=rule(n)
        return w@np.array([self.point(self.invert(y+s*(z-y),a+s*(b-a)))[2] for s in t])

    def flux(self,p,u):
        _,_,g,T,_=self.point(u)
        f=integrate(p,T,2*g[0],T*g[1])[0]
        v=f[:2]*[2.,T]
        if abs(g@v-f[2])>6e-14*max(abs(f)):
            raise ArithmeticError('Finite reservoir force did not match physical FD entropy')
        return v,f[2]

    def step(self,p,a,dt):
        y,L,_,_,_=self.point(a)

        def packet(b):
            z,M,_,_,_=self.point(b)
            g=self.path(a,b,y,z,5)
            def solve(l):
                zz,s,h,T,A=self.point([g[0],l[0]])
                return [h[1]-g[1]], [[np.exp(-l[0])+A[1,1]/(4*C_VIS*T**5)]]
            l=root(lambda q:solve(q)[0],[(a[1]+b[1])/2],
                   jac=lambda q:solve(q)[1],tol=6e-13).x[0]
            if abs(solve([l])[0][0])>3e-12*max(1.,abs(g[1])):
                raise ArithmeticError('Averaged force does not bind an actual isoenergy state')
            v,s=self.flux(p,[g[0],l])
            return (z-y+dt*v)/y,z,M,g,s

        b=root(lambda x:packet(x)[0],a,tol=3e-12).x
        residual,z,M,g,s=packet(b)
        c=max(abs(L),abs(M),dt*s)
        err=max(abs(M-L+dt*s),abs(g@(z-y)-(M-L)))/c
        if (np.max(abs(residual))>5e-13 or err>5e-13 or s<=0.
                or np.max(abs(g-self.path(a,b,y,z,12)))>4e-11):
            raise ArithmeticError('Independent moment chord or total entropy balance failed')
        return b,err,dt*s


def run():
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T,eta,ratio in ((.05,.3,1.16),(.3,-.75,1.16),(.6882487043,.9,.86)):
        a=np.array([eta,np.log(T*ratio)])
        major=Bath(T,a,168)
        y0,L,_,_,_=major.point(a)
        f,_=major.flux(q,a)
        tau=.10/np.max(abs(f/y0))
        def rhs(t,u):
            y,s,g,v,A=major.point(u)
            return np.linalg.solve(A,-tau*major.flux(q,u)[0])
        ref=solve_ivp(rhs,(0,1),a,method='DOP853',rtol=3e-12,atol=4e-14)
        if not ref.success:
            raise ArithmeticError('Missing independent finite-reservoir differential trajectory')
        target=major.point(ref.y[:,-1])[0]

        def follow(portal,bath,k):
            z=a.copy()
            ent=[]
            defects=[]
            for _ in range(k):
                z,e,s=bath.step(portal,z,tau/k)
                ent.append(s)
                defects.append(e)
            y,M,_,Ts,_=bath.point(z)
            start=bath.point(a)[1]
            if abs(M-start+sum(ent))>2e-12*abs(start):
                raise ArithmeticError('Integrated entropy production did not telescope')
            return dict(steps=k,visible_T=Ts,dark_T=np.exp(z[1]),eta=z[0],
                        normalized_reference_error=np.max(abs(y-target)/y0),total_entropy_produced=sum(ent),
                        max_entropy_balance_error=max(defects))
        out=[follow(q,major,n) for n in (2,4,8)]
        independent=follow(p,Bath(T,a,112),4)
        fields=('visible_T','dark_T','eta')
        agreement=max(abs(out[1][k]-independent[k]) for k in fields)
        x=np.array([r['normalized_reference_error'] for r in out])
        ratios=x[:-1]/x[1:]
        if min(ratios)<3.3 or max(ratios)>4.8 or agreement>2e-7:
            raise ArithmeticError('Second-order agreement or independent finite-bath quadrature failed')
        rows.append(dict(initial_visible_T=T,initial_dark_T=T*ratio,initial_eta=eta,total_energy=major.E,
                         levels=out,halving_ratios=ratios.tolist(),independent_mesh_difference=agreement))
    return dict(rows=rows,scope='Closed homogeneous collision substep with a finite ideal visible reservoir. No expansion, tabulated Standard Model EOS, complete interaction network, or abundance claim.')


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
