"""Reusable positive quadrature for symmetric, spin-correlated two-bath FD reactions."""
import json
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre
from timelike_portal_onsager_coverage_reference import mchi,V,width,factor
from timelike_chemical_heat_eigenmodes_reference import leptons
from finite_width_vector_pauli_rate_reference import integrator

@lru_cache
def rule(n):
    return roots_legendre(n)

class PositivePortal:
    def __init__(self,orders=(9,24,16), tail=55.,dark_mass=mchi):
        self.orders=orders
        self.tail=tail
        self.dark_mass=dark_mass

    def geometry(self,T):
        N,M,_=self.orders
        z,w=rule(N)
        y,v=rule(M)
        yy=(y+1)*np.sqrt(self.tail)/2
        hi=max(60*T,V+30*T)
        lo=2*self.dark_mass
        ends=sorted(set((lo,hi,V,*[2*m for m in leptons],*[V+i*width for i in
             (-300,-100,-30,-10,-3,-1,1,3,10,30,100,300)],*[T*j for j in
             (.2,.5,1,2,3.5,5.5,8,12,18,27,39,52,60)])))
        ed=np.array([x for x in ends if lo<=x<=hi])
        # Square-root excess coordinate removes each lower-edge threshold
        # cusp. Constant-width BW features have their own subintervals.
        u=np.sqrt(np.log(ed[1:]/ed[:-1]))[:,None]*(z+1)/2
        E=(ed[:-1,None]*np.exp(u*u)).ravel()
        dE=(w*np.sqrt(np.log(ed[1:]/ed[:-1]))[:,None]*u).ravel()*E
        q=E[:,None]+T*yy*yy
        p=np.sqrt(T*yy*yy*(q+E[:,None]))
        dq=v*T*np.sqrt(self.tail)*yy
        vi2=1-(lo/E)**2
        b0=(factor*(E**2+2*self.dark_mass**2)*E**3*np.sqrt(vi2)*
            (2*E*dE/((E**2-V**2)**2+(V*width)**2)))
        W=b0[:,None]*(dq*p/E[:,None])
        return E,q,p,W,vi2

    def moments(self,Ts,Td,eta,*,common_momentum=False):
        assert Ts>0. and Td>0. and abs(eta)<50.
        hot=max(Ts,Td)
        E,q,p,W,vi2=self.geometry(hot)
        x,w=rule(self.orders[2])
        # Both members carry the same chemical potential. The blocking
        # product is even in x; integrate only x>=0 after remapping.
        x=(x+1)/2; w=w/2
        wx=w*((3*x*x-1)/2)
        def block(T,mu,v2,p_,q_):
            a=q_/(2*T)-mu
            b=p_[:,:,None]*np.sqrt(v2)[:,None,None]*x/(2*T)
            a=a[:,:,None]
            return 1/(1+np.exp(-a+b)+np.exp(-a-b)+np.exp(-2*a))
        b=block(Td,eta,vi2,p,q)
        a0,a2=b@w,b@wx
        B=np.zeros_like(q)
        C=np.zeros_like(q)
        s=E**2
        for m in leptons:
            ok=E>2*m
            vf2=1-4*m*m/s[ok]
            v=vi2[ok]
            c=block(Ts,0.,vf2,p[ok],q[ok])
            f2=((2/3)*v*vf2/(3-v-vf2+v*vf2/3))[:,None]
            amp=((1+2*m*m/s[ok])*np.sqrt(vf2))[:,None]
            B[ok]+=amp*(a0[ok]*(c@w)+f2*a2[ok]*(c@wx))
            C[ok]+=amp
        d=2*eta+q*(Td-Ts)/(Td*Ts)
        a=np.maximum(-q/Ts,2*eta-q/Td)
        I=W*B
        flux=I*np.exp(a)*np.sign(d)*(-np.expm1(-abs(d)))
        eq=np.exp(-q/Ts)*I
        mb=np.exp(-q/Ts)*W*C
        if common_momentum:
            if Ts!=Td or eta!=0.:
                raise ValueError("Common momentum susceptibility requires full FD detailed balance.")
            return np.array([np.sum(eq*q*q)/Ts, np.sum(eq*p*p)/(3*Ts),
                             np.sum(eq*E[:,None]**2)/Ts])
        return np.array([np.sum(flux),np.sum(flux*q),np.sum(flux*d),
                       np.sum(eq),np.sum(eq*q),np.sum(eq*q*q),np.sum(mb)])

    def flow(self,Ts,Td,eta):
        return self.moments(Ts,Td,eta)[:3]

def test():
    p=PositivePortal()
    ref=PositivePortal((13,36,26))
    rows=[]
    for Ts,Td,eta in ((.01,.005,-3.),(.06,.07,.01),(.1,.12,-.43),(.3,.31,-.001),
                      (.45,.41,.05),(.6882487043,.6882487043,0.),(.688,.688001,-.000001),
                      (.69,.88,.3),(.69,.54,-.5),(1.,.8,.7),(2.,1.8,-1.)):
        u=p.moments(Ts,Td,eta)
        v=ref.moments(Ts,Td,eta)
        r=integrator(Td,Ts,eta)
        flu=float(np.max(abs(u[:2]-r[3,:2])/np.maximum(abs(r[3,:2]),1e-290)))
        eq=float(np.max(abs(u[3:6]/r[1]-1)))
        df=float(np.max(abs(u-v)/np.maximum(abs(v),1e-290)))
        assert u[2]>=0. and max(flu,eq,df)<2e-5,(Ts,Td,eta,flu,eq,df)
        rows.append(dict(T_SM=Ts,T_dark=Td,eta_dark=eta,refinement=df,
                    flux_vs_adaptive=flu,reference_weight_vs_adaptive=eq))
    return dict(rows=rows,scope="Positive fixed quadrature for the local vacuum FD vector portal including spin correlations.")

if __name__=="__main__":
    print(json.dumps(test(),indent=2,allow_nan=False))
