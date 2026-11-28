"""Short dynamic Pauli two-moment closure with variational gauge energy."""
import json
import argparse
import numpy as np
from scipy.optimize import brentq, root
from scipy.integrate import solve_ivp
from scipy.special import expit, roots_legendre, xlogy



class Cell:
    def __init__(self,l,h,m,order):
        self.lo,self.hi,self.m=l,h,m
        ee=np.hypot([l,h],m)
        self.ec,self.de=np.mean(ee),np.diff(ee)[0]/2
        x,w=roots_legendre(order)
        p=(l+h)/2+(h-l)*x/2
        self.w=w*(h-l)/2*p*p*2/np.pi**2
        self.energy=np.hypot(p,m)
        self.u=(self.energy-self.ec)/self.de
        self.B=np.stack([np.ones(len(p)),self.u])
        self.I=(self.w[:,None]*self.B.T).T
        self.inv=(np.pi**2/2)/self.energy
        self.W=self.w.sum()

    def initial(self,a,b):
        return self.I@expit(a+b*self.u)

    def fit(self,Y):
        q=Y[0]/self.W
        assert 0.<q<1.,q
        def f(z):return expit(z@self.B)
        def calc(z):
            o=f(z)
            D=o*(1-o)
            return (self.I@o-Y)/self.W, (self.I*D)@self.B.T/self.W
        b0=np.array([np.log(q/(1-q)),0.])
        fit=root(lambda z:calc(z)[0],b0,jac=lambda z:calc(z)[1],tol=2e-10)
        z=fit.x
        assert max(abs(calc(z)[0]))<2e-10,(Y,calc(z)[0])
        F=f(z)
        V=self.W*calc(z)[1]
        J=self.w@(self.inv*F)
        derivative=((self.I*(F*(1-F)*self.inv)).sum(axis=1))@np.linalg.inv(V)
        S=(-self.w)@(xlogy(F,F)+xlogy(1-F,1-F))
        return z,J,derivative,S

    def basis(self,p):
        return (np.hypot(p,self.m)-self.ec)/self.de

    def pfrombare(self,E):
        return np.sqrt(np.maximum((E-self.m)*(E+self.m),0.))


class Tool:
    def __init__(self,n,g,order,hartree):
        self.order=order;self.n=n;self.g=g
        self.hartree=hartree
        self.Js=np.zeros(2)
        self.clocknorm=1.
        self.q=4.1e-6
        self.es=np.hypot(.00086724,self.q)
        self.wscalar=.002
        self.lambda_s=2*.00086724**2*self.es/(100*self.wscalar)
        self.fieldcoin=g*g/np.pi**4
        self.jgamma=np.pi**2/6
        self.nod,self.w=roots_legendre(order)
        self.cells=[]
        for m,l,h in ((.05,1.315,1.555),(.003,1.314,1.554)):
            e=np.linspace(l,h,n+1)
            self.cells +=[Cell(a,b,m,order) for a,b in zip(e[:-1],e[1:],strict=True)]
        self.ec=np.array([c.ec for c in self.cells])
        self.de=np.array([c.de for c in self.cells])
        self.init=np.r_[np.array([c.initial(-1.65-2.4*(c.ec-1.44),-2.4*c.de)
                       if j<n else c.initial(-2.65-1.4*(c.ec-1.44),-1.4*c.de)
                       for j,c in enumerate(self.cells)]).ravel(),100*self.wscalar]

    def state(self,y,mode):
        a=y[:-1].reshape(-1,2)
        f=[c.fit(z) for c,z in zip(self.cells,a,strict=True)]
        J=np.array([q[1] for q in f])
        S=sum(q[3] for q in f)
        x=y[-1]/self.wscalar
        S+=self.wscalar*(np.log1p(x)+x*np.log1p(1/x))
        j=self.Js+np.r_[sum(J[:self.n]),sum(J[self.n:])]
        G=self.fieldcoin*(self.jgamma*sum(j)+.5*j@j)
        E=np.column_stack([self.ec,self.de])
        es,field=self.scalar(y[-1])
        energy=np.sum(E*a)+G+field
        Egrad=E+np.repeat(self.fieldcoin*(self.jgamma+j),self.n)[:,None]*np.array([r[2] for r in f])
        if mode=='bare':
            ee=E
        else:
            ee=Egrad
        return dict(E=energy,entropy=S,grad=Egrad,field=G,coeff=ee,
                    f=np.asarray([r[0] for r in f]),x=x,es=es)

    def scalar(self,ns):
        if self.hartree:
            # Gap m² = lambda N/(2 E) with m0 = 0, so E solves
            # E³ - q² E - lambda N/2 = 0 on the physical branch.
            es=brentq(lambda z:z*(z*z-self.q**2)-.5*self.lambda_s*ns,self.q,.02,xtol=5e-18)
            field=ns*es-(es*es-self.q*self.q)**2/(2*self.lambda_s)
        else:
            es=self.es
            field=es*ns
        return es,field

    def events(self,box):
        energies=box['coeff'];es=box['es']
        for i,A in enumerate(self.cells[:self.n]):
            a0,a1=energies[i]
            for j,B in enumerate(self.cells[self.n:],self.n):
                b0,b1=energies[j]
                assert a1>0. and b1>0.
                ep=lambda p:a0+a1*A.basis(p)
                down=lambda p:B.pfrombare(B.ec+B.de*(ep(p)-es-b0)/b1)
                # Intersection of the two energy supports, then exact spherical triangle.
                el=max(a0-a1,b0-b1+es)
                er=min(a0+a1,b0+b1+es)
                if er<=el: continue
                p0=A.pfrombare(A.ec+A.de*(el-a0)/a1)
                p1=A.pfrombare(A.ec+A.de*(er-a0)/a1)
                # Return disjoint valid intervals, checking both angular boundaries.
                s=np.linspace(p0,p1,12)
                cuts=[p0,p1]
                for shift in (self.q,-self.q):
                    f=lambda p:down(p)-p+shift
                    z=f(s)
                    for k in np.flatnonzero(z[:-1]*z[1:]<0):
                        cuts.append(brentq(f,s[k],s[k+1],xtol=5e-14))
                zz=np.sort(cuts)
                for lo,hi in zip(zz[:-1],zz[1:],strict=True):
                    mid=(lo+hi)/2
                    if not abs(down(mid)-mid)<=self.q:continue
                    q=mid+(hi-lo)*self.nod/2
                    k=down(q)
                    wp=ep(q)
                    wb=wp-es
                    assert max(abs(k-q))<=self.q*(1+1e-6)
                    up=A.basis(q); ub=B.basis(k)
                    Eb=np.hypot(k,B.m)
                    vb=(b1/B.de)*k/Eb
                    Z=(hi-lo)/2*self.w*(q*k/(wp*wb*vb))
                    yield i,j,up,ub,Z

    def rhs(self,t,y,mode):
        B=self.state(y,mode)
        T=np.zeros((2*self.n,2))
        S=0.; Ent=0.
        for i,j,u,v,w in self.events(B):
            a=B['f'][i]@np.stack([np.ones_like(u),u])
            b=B['f'][j]@np.stack([np.ones_like(v),v])
            f,h=expit(a),expit(b)
            gain=(w/self.clocknorm)*(f*(1-h)+B['x']*(f-h))
            p=gain.sum()
            T[i]-=np.r_[p,gain@u]
            T[j]+=np.r_[p,gain@v]
            S+=p
            Ent+=gain@(a-b+np.log1p(1/B['x']))
        assert Ent>=-1e-14
        dot=np.sum(B['grad']*T)+B['es']*S
        if mode=='closed':
            assert abs(dot)<5e-13
        return np.r_[T.ravel(),S],Ent,dot

    def solve(self,mode,tmax):
        def fn(t,z):
            f,ent,W=self.rhs(t,z[:-2],mode)
            return np.r_[f,ent,W]
        ini=self.state(self.init,mode)
        # Confirm the independent scalar and two-local-moment gradients
        # of the full reconstructed free + interaction energy functional.
        v=np.linspace(.1,.5,self.init.size)*self.init
        ex=4e-5
        value=lambda a:self.state(self.init+a*v,mode)['E']
        der=(value(-2*ex)-8*value(-ex)+8*value(ex)-value(2*ex))/(12*ex)
        exact=np.r_[ini['grad'].ravel(),ini['es']]@v
        assert abs(der/exact-1)<5e-8
        sol=solve_ivp(fn,(0.,tmax),np.r_[self.init,0.,0.],method='DOP853',
                      rtol=3e-10,atol=2e-13,t_eval=np.linspace(0.,tmax,11),max_step=tmax/7)
        assert sol.success,sol.message
        n=sol.y[:-2,-1]
        B=self.state(n,mode)
        Y=(n-self.init)[:-1].reshape(-1,2)[:,0]
        S=n[-1]-self.init[-1]
        defects=[sum(Y[:self.n])+S,sum(Y[self.n:])-S]
        assert max(abs(v) for v in defects)<2e-13
        drift=B['E']-ini['E']
        work=sol.y[-1,-1]
        assert abs(drift-work)<2e-12
        entropy=B['entropy']-ini['entropy']
        assert abs(entropy-sol.y[-2,-1])<2e-11
        if mode=='closed': assert abs(drift)<2e-13
        f0,f1=1.,0.
        for x in sol.y[:-2].T:
            fit=self.state(x,mode)
            for cell,theta in zip(self.cells,fit['f'],strict=True):
                endpoints=expit(theta[0]+theta[1]*np.array([-1.,1.]))
                f0=min(f0,np.min(endpoints))
                f1=max(f1,np.max(endpoints))
        return dict(mode=mode,trajectory_time=tmax,scalar_number_gain=S,
                    initial_scalar_energy=ini['es'],final_scalar_energy=B['es'],
                    full_energy_gradient_fd_relative=der/exact-1,
                    energy_drift=drift,independent_work=work,entropy_increase=entropy,
                    entropy_ledger=sol.y[-2,-1],field_energy_change=B['field']-ini['field'],
                    two_number_defects=defects,minmax_occupations=[f0,f1])


def run(n,order,hartree,tmax):
    E=Tool(n,.005,order,hartree)
    return dict(scope='Continuous on-shell 1↔2 angles with two locally closed Pauli moments '
                      'per hard cell, optional live scalar gap and fixed photon moment, dimensionless exposure; '
                      'no physical kinetic rate, scalar freeze-out or FRW.',
                    n=n,order=order,hartree=hartree,rows=[E.solve(m,tmax) for m in ('bare','closed')])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cells',type=int,default=3)
    p.add_argument('--order',type=int,default=12)
    p.add_argument('--fixed-scalar',action='store_true')
    p.add_argument('--exposure',type=float,default=.19)
    a=p.parse_args()
    print(json.dumps(run(a.cells,a.order,not a.fixed_scalar,a.exposure),indent=2,allow_nan=False))
