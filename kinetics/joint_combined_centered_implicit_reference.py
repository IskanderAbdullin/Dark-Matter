"""Joint prompt and exact off-pole elastic implicit test in a fixed finite box."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import expit
from joint_resonant_centered_stiff_reference import advance,Centered
from joint_resonant_kinetic_linear_modes_reference import configuration,metric
from joint_below_pole_elastic_linear_modes_reference import setup
from sm_dark_entropy_exchange_reference import StandardModel

class Combined(Centered):
    def __init__(self,k,H,panels):
        super().__init__(k,H)
        self.panels=panels
    def reference(self,q,guess):
        root,R0,t=super().reference(q,guess)
        self.t=t
        r=self.arrays(root)
        self.eta=[(c.mass+c.low[0])/t-r[j][0,0] for j,c in enumerate(self.k.c)]
        return root,R0,t

    def elastic(self,z):
        # Each panel is symmetric in incoming/outgoing angles, and is
        # sampled on the same finite support in both directions.
        n=2*self.n
        value=np.zeros(self.d)
        saved=[]
        for p,c,L,P,raw in self.panels:
            P1=P[:,:n];P2=P[:,n:]
            z1=P1@z[:n];z2=P2@z[n:]
            aa=p/self.t-self.eta[0]+z1
            bb=c/self.t-self.eta[1]+z2
            if np.min(aa)<=0.:raise ValueError("BE exponent")
            f=np.exp(-aa)/(-np.expm1(-aa));g=expit(-bb)
            block=(1+f)*(1-g)
            K=raw*(block[:,None]*block)*np.exp(sum(self.eta)-L/self.t)/self.H
            K1=K.sum(axis=1)
            u=np.expm1(-z1-z2)
            v=K@u-u*K1
            value+=P.T@v
            saved.append((P,P1,P2,K,K1,u,v,f,g))
        jac=None
        def derivative(d):
            nonlocal jac
            if jac is None:
                jac=np.zeros((self.d,self.d))
                for P,P1,P2,K,K1,u,v,f,g in saved:
                    C=np.concatenate((-f[:,None]*P1,g[:,None]*P2),axis=1)
                    D=-(1+u)[:,None]*P
                    term=v[:,None]*C+K@(u[:,None]*C+D)-u[:,None]*(K@C)-K1[:,None]*D
                    jac+=P.T@term
            return jac@d
        return value,derivative
    def centered(self,z,ref,R0):
        a,da=super().centered(z,ref,R0)
        b,db=self.elastic(z)
        return a+b,lambda v:da(v)+db(v)

def test(r,eos):
    T,H,k,N,t,Q,sh=configuration(r,eos,8,6.5)
    n=len(k.c[0].vol)
    start=np.concatenate([a.T.ravel() for a in sh])
    start[0]-=.09
    res=Centered(k,H);X,_,_=res.moment(start)
    root,R,t=res.reference(start,T)
    Eq=res.moment(root)[0]
    p=[]
    setup(T,k,t,[k.c[i].mass/t-res.arrays(root)[i][0,0] for i in (0,1)],4,3,2,10,30.,.82,cache=p)
    co=Combined(k,H,p)
    ref,R,t=co.reference(start,T)
    z=start-ref
    f,j=co.centered(z,ref,R)
    d=np.cos(np.arange(len(z))*.34)
    e=3e-6
    chk=np.linalg.norm((co.centered(z+d*e,ref,R)[0]-co.centered(z-d*e,ref,R)[0])/(2*e)-j(d))/np.linalg.norm(j(d))
    assert chk<8e-8
    _,_,C=metric(k,sh)
    results=[]
    for name,m in (("resonant",res),("combined",co)):
        y,x=advance(m,start.copy(),.02,T)
        M,_,_=m.moment(y)
        err=np.max(abs((C.T@M)/(C.T@X)-1))
        assert err<3e-10
        keep=float((M[0]-Eq[0])/(X[0]-Eq[0]))
        results.append(dict(channels=name,invariant_defect=err,remaining_cold_cell_excess_relative_to_common_reference=keep,**x))
    return dict(panels=len(p),scalar_first_bin_upper_p_over_visible_T=k.c[0].p[1]/T,
       directional_analytic_derivative_error=chk,steps=results,
       scope="Single fully nonlinear backward-Euler test with actual projected prompt and below-pole elastic scattering, fixed support and manufactured populated initial shapes. Not a physical FRW result or exponential relaxation time.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos');p.add_argument('spectra')
    a=p.parse_args();eos=StandardModel(a.eos)
    r=next(r for r in json.loads(Path(a.spectra).read_text())['rows'] if r['TR']==.05 and abs(r['TF']-.045)<1e-8)
    print(json.dumps(test(r,eos),indent=2,allow_nan=False))
