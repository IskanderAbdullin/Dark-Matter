"""Expansion check for the common hard-moment state with vacuum packet weights."""
import argparse
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit
from early_pole_live_joint_hard_pair_packets_reference import Packets


class Expanding(Packets):
    eta=0.

    def weights(self,Y,closed):
        fits=[]; J=[]; Jgrad=[]; Jred=[]; Efree=[]; Egrad=[]; Ered=[]; S=0.
        a2=np.exp(2*self.eta)
        for k,U in enumerate(Y):
            C=self.cells[k]
            if k<self.r[0].stop:
                a,b,j,grad,s=C.fit(U)
                fits.append(np.r_[a,b])
                J.append(j);Jgrad.append(grad);Jred.append(0.)
                Efree.append(self.shape[k]@U);Egrad.append(self.shape[k])
                Ered.append(0.)
            else:
                ab,_,_,s=C.fit(U)
                v=expit(ab@C.B)
                p2=C.energy*C.energy-C.m*C.m
                bare=np.sqrt(p2+a2*C.m*C.m)
                Jkernel=(np.pi**2/2)/bare
                D=v*(1-v)
                I=C.I*D
                inv=np.linalg.inv(I@C.B.T)
                f=C.w*v
                fits.append(ab)
                J.append(Jkernel@f)
                Jgrad.append((I@Jkernel)@inv)
                Jred.append(-f@(Jkernel*a2*C.m*C.m/(bare*bare)))
                Efree.append(bare@f)
                Egrad.append((I@bare)@inv)
                Ered.append(f@(a2*C.m*C.m/bare))
            S+=s
        background,db,ext,dext=(self.spectators(self.eta) if hasattr(self,'spectators')
                                 else (0.,np.zeros(3),0.,0.))
        j=background+np.array([sum(J[sl]) for sl in self.r])
        vec=np.r_[sum(j[1:]),j[0]+j[1:]]
        C=np.empty(len(J))
        for i,sl in enumerate(self.r):C[sl]=vec[i]
        E=self.c*C[:,None]*np.asarray(Jgrad)+np.asarray(Egrad)
        R=sum(Efree)+ext+.5*self.c*(j@vec)
        self.trace=sum(Ered)+dext+self.c*(C@Jred+vec@db)
        return E if closed else np.asarray(Egrad),np.array(fits),E,R,S

    def go_frw(self,closed,decay,pair,duration=.14,etaend=.028):
        self.pair_scale=pair
        a=np.r_[self.init.ravel(),self.Ns,0.,0.,0.,0.,0.]
        def energy(eta,y):
            self.eta=eta
            r,s=self.weights(y[:-6].reshape(-1,2),closed)[-2:]
            _,rs,ss,_=self.scalar(y[-6])
            return r+rs,s+ss,self.trace
        def rhs(t,y):
            self.eta=t/duration*etaend
            dy=self.rhs_joint(y[:-1],closed,decay)
            return np.r_[dy,self.trace*etaend/duration]
        h=2e-5
        Y=[energy(x*h,a)[0] for x in (-2,-1,1,2)]
        diff=(Y[0]-8*Y[1]+8*Y[2]-Y[3])/(12*h)
        r0,s0,trace=energy(0.,a)
        assert abs(diff-trace)<3e-10
        S=solve_ivp(rhs,(0.,duration),a,method='DOP853',rtol=1e-10,atol=3e-13,max_step=.018)
        assert S.success,S.message
        y=S.y[:,-1]
        r1,s1,_=energy(etaend,y)
        R=r1-r0-y[-1]
        assert abs(R-y[-4])<2e-12
        if closed:assert abs(R)<3e-13
        assert abs(s1-s0-y[-5])<3e-11
        assert abs(sum((y-a)[:-6].reshape(-1,2)[:,0]))<1e-13
        return dict(closed=closed,decay=decay,pairs=pair,log_scale=etaend,
                    scalar_gain=y[-6]-self.Ns,comoving_redshift_trace_integral=y[-1],
                    trace_defect_fixed_state=diff-trace,redshift_subtracted_drift=R,
                    collision_work=y[-4],entropy_gain=s1-s0,remap_rejection=self.reject)


def run(level):
    quad=[(2,3,10,-3),(3,5,16,-4)][level]
    out=[]
    for c,d,p,e in ((True,0.,0.,.028),(True,1.,0.,.028),(True,1.,1.,.028),
                    (True,1.,1.,0.),(False,1.,1.,.028)):
        out.append(Expanding(12,quad,p).go_frw(c,d,p,etaend=e))
    return dict(scope='Tagged comoving hard domains, vacuum packets and arbitrary exposure per log scale; '
                      'no Friedmann fit, calibrated cosmological rates, or matching to hard exterior.',
                quad=quad,rows=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0)
    a=p.parse_args()
    print(json.dumps(run(a.level),indent=2,allow_nan=False))
