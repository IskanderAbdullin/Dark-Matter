"""Phase normalization with independent gamma matrices on the lightlike benchmark."""
import argparse,json
import numpy as np
from scipy.special import roots_legendre
from finite_mass_qed_absorption_transfer_reference import point,transfer
from finite_mass_qed_timelike_pole_angular_reference import angular_knots
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q
from finite_mass_qed_timelike_pole_batched_absorption_reference import load,radial


def guides(M,m,T,s,A,v,lo,hi):
    t=T[0]
    guides=[]
    # Explicitly resolve both nearly collinear charged denominators.
    ladder=np.array([.006,.032,.15,.7,3.,14.,65.,300.,1500.,7000.])
    pz=(t*lo-A/2)/v
    deriv=abs(2*(lo-t*pz/v))
    for mass in (M,m):
        guides.extend(lo+mass*mass/max(deriv,1e-15)*ladder)
    sZ=(t-T[3])*(t+T[3])
    for mass,J,shift in ((M,A/2,0.),(m,A/2-s,t)):
        disc=J*J-mass*mass*sZ
        if disc>0.:
            big=t*J+np.copysign(abs(T[3])*np.sqrt(disc),t*J)
            for E0 in (big/sZ,(J*J+(mass*T[3])**2)/big):
                edge=E0+shift
                if lo-.03*(hi-lo)<edge<hi+.03*(hi-lo):
                    pz=(t*edge-A/2)/v
                    grad=abs(2*(edge-t*pz/v))
                    scale=(mass*mass+2*mass*abs(pz*T[1]/v))/max(grad,1e-15)
                    guides.extend(edge+scale*np.r_[-ladder,0.,ladder])
    return guides

def angle_breaks(ch,k,es,m_s):
    points,_=angular_knots(ch,k,k,es)
    M,m=[x.m for x in ch]
    a=[]
    for c in points[1:-1]:
        E0,_,S,t,V,D,l,h,co,si=transfer(M,m,m_s,k,c)
        p=np.sqrt((l-M)*(l+M))
        eb=l-t
        b=np.sqrt((eb-m)*(eb+m))
        dE=abs(k*Q/(V*(eb/b-l/p)))
        z=(t*l-D/2)/V
        derivative=abs(2*(l-t*z/V))
        scale=(m*m+2*m*abs(z*si))/(derivative*dE)
        a.extend(c+scale*np.r_[-65.,-14.,-3.,-.7,-.15,-.03,.03,.15,.7,3.,14.,65.])
    return np.unique(np.clip(np.r_[points,a,-.985,-.9,-.6,0.,.6,.9,.985],-1.,1.))


def calc(g,ks):
    d,ch,es=load(g)
    for k in ks:
        b=point(d,k,37,58)
        if g=='third' and k==.00256:
            bnd=angle_breaks(ch,k,es,d['checkpoints'][-1]['scalar_pole_over_T'])
            meshes=((10,24,4),(16,48,8))
        else:
            bnd,_=angular_knots(ch,k,k,es)
            meshes=((8,24,10),(16,48,18))
        for n,az,oa in meshes:
            x,w=roots_legendre(oa)
            stat=np.zeros(3)
            for l,h in zip(bnd[:-1],bnd[1:],strict=True):
                for c,u in zip((l+h)/2+(h-l)*x/2,(h-l)*w/4,strict=True):
                    stat+=u*radial(ch,k,k,1/(2*k),es,c,'T',n,az,extra_energy_breaks=guides)[0]
            value=k*k*stat[:2]/(32*np.pi**3*es)
            r=float(max(abs(value/b-1.)))
            if n==16 and r>1e-5:
                raise ArithmeticError(('Lightlike normalization check unresolved',g,k,r))
            print(json.dumps(dict(grid=g,k=k,orders=[n,az,oa],phase_value=value.tolist(),legacy=b.tolist(),rel=r)),
                  flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--k',type=float,choices=(.00128,.00256),action='append')
    o=p.parse_args()
    calc(o.grid,o.k or (.00128,.00256))
