"""Coupled medium-aware weights and quantum extents in two events."""
import json
import numpy as np
from scipy.optimize import brentq,least_squares,root
from two_mass_kinetic_offset_deposition_reference import DT,event
from two_mass_discrete_entropy_gradient_reference import averaged_alphas
from two_mass_entropy_extent_reference import entropy


def run():
    e,m,b,P,g=map(DT,(.046,.002,.007,.36,.012))
    M=e+m+b
    energies=[event(e,m,b,P,x)[1]+np.array([m,b],dtype=DT) for x in (-.38,-.31)]
    es,fd=np.array(energies).T
    ep=np.hypot(M,P)
    assert np.max(abs(ep-es-fd))<2e-16*ep
    l,h=min(fd)*DT(.94),max(fd)*DT(1.06)
    bare=np.array([ep,l,h,*es],dtype=DT)
    gamma=g*g/DT(np.pi)**4
    Jext=np.array([.62,.32],dtype=DT)
    B=np.zeros((2,5),dtype=DT); B[0,0]=DT(np.pi)**2/(2*ep)
    B[1,1:3]=DT(np.pi)**2/(2*bare[1:3])
    photon=DT(.27)
    H=gamma*B.T@B
    def energy(z):
        J=Jext+B@z
        return bare@z+gamma*(photon*sum(J)+(J@J)/2)
    def grad(z): return bare+gamma*B.T@(photon+Jext+B@z)
    gcap=np.array([.9,.8,.9,1.,1.],dtype=DT)
    bose=np.array([False,False,False,True,True])
    n=gcap*np.array([.16,.29,.31,.005,1.1],dtype=DT)
    def direction(j,w):
        a=np.array([-1.,1-DT(w),DT(w),0.,0.],dtype=DT); a[3+j]=1.
        return a
    def S(mean,j,v):
        scalar=mean[3+j]
        final=mean[1:3]@v[1:3]
        A=mean@v
        lognot=-np.logaddexp(DT(0),mean[0])-np.logaddexp(DT(0),-final)-np.log(-np.expm1(-scalar))
        return np.sign(A)*np.exp(max(lognot,lognot-A))*(-np.expm1(-abs(A))),A
    res=[]
    for dt in (.03,.1,.25):
        def equations(z,output=False):
            xi=np.array(z[:2],dtype=DT)
            V=np.array([direction(j,w) for j,w in enumerate(z[2:])],dtype=DT).T
            y=n+V@xi
            if min(y)<=0 or max((y/gcap)[~bose])>=1:
                raise ValueError("Invalid test gamut.")
            mean,_=averaged_alphas(n,y,gcap,bose)
            net,aff=np.array([S(mean,j,d) for j,d in enumerate(V.T)]).T
            if output: return y,V,net,aff
            return np.asarray(np.r_[xi-dt*net, (grad((n+y)/2)@V)/(h-l)],float)
        # Independent feasibility box guarantees the relevant occupations.
        x=least_squares(equations,[.005,-.005,.5,.5],bounds=([-.0049,-.035,0,0],[.05,.05,1,1]),
                        xtol=2e-14,gtol=2e-14,ftol=2e-14)
        polished=root(equations,x.x,tol=1e-11)
        vector=polished.x
        y,V,s,a=equations(vector,True)
        xi=vector[:2]
        if (not x.success or max(abs(equations(vector)))>4e-16 or xi[0]<=0 or xi[1]>=0
                or min(vector[2:])<=0. or max(vector[2:])>=1.):
            raise ValueError("Coupled entropy/energy root not accurately resolved.")
        delta=y-n
        number=max(abs(delta[0]+delta[1:3].sum()),abs(delta[0]+delta[3:].sum()))
        if number>2e-16:
            raise ValueError("The coupled medium step did not retain both daughter number ledgers.")
        Ee=abs(energy(y)-energy(n))
        es=abs(entropy(y,gcap,bose)-entropy(n,gcap,bose)-dt*(a@s))
        if min(a*s)<0 or Ee>2e-15 or es>2e-14:
            raise ValueError(f"Energy or entropy ledger failed: {Ee}, {es}, {a*s}.")
        # Match each finite event against the original state in isolation.
        separate=[]
        for j in range(2):
            def f(w):
                v=direction(j,w)
                return float(grad(n+DT(xi[j])*v/2)@v)
            separate.append(direction(j,brentq(f,0.,1.)))
        separate=np.array(separate,dtype=DT)
        summed=n+np.asarray(xi,dtype=DT)@separate
        drift=energy(summed)-energy(n)
        theory=DT(xi[0])*DT(xi[1])*(separate[0]@H@separate[1])
        if abs(drift/theory-1)>5e-7:
            raise ValueError("Quadratic mixed-medium drift did not reproduce.")
        res.append(dict(step=dt,extents=xi.tolist(),weights=vector[2:].tolist(),
                        coupled_energy_defect=float(Ee),entropy_identity_defect=float(es),count_ledger_defect=float(number),
                        separate_weight_energy_drift=float(drift),mixed_term=float(theory)))
    return res


if __name__=="__main__":
    print(json.dumps(dict(scope="Manufactured two-event quadratic medium control; not matched thermal rates.",rows=run()),
                     indent=2,allow_nan=False))
