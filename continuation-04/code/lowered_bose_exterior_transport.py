"""Frozen-potential C22 tracer transport including birth radius, angular momentum,
inward/outward travel, and first Kepler excursion. Prescribed initial source held
constant, no subsequent collisions or evolving potential. Tracer ledgers only.
"""
import argparse,json,math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator
from lowered_antonov_hessian import Halo
from lowered_bose_exterior_source import density,exterior_min,period,R,V,GYR,TIME
from lowered_bose_collision_escape import rate

def nodes(n):
    x,w=leggauss(n);return (x+1)/2,w/2

def trajectories(H,r,eps,mu,potential,ng):
    e0=H.Qt/H.rt;rt=H.rt;X=r[:,None]**2
    L2=2*X*(potential(X)-eps)*(1-mu[None,:]**2)
    def D(x):return 2*x*(potential(x)-eps)-L2
    a=np.zeros_like(L2);b=np.broadcast_to(X,L2.shape).copy()
    for _ in range(48):
        m=(a+b)/2;above=D(m)>0;b=np.where(above,m,b);a=np.where(above,a,m)
    p=(a+b)/2;z,w=nodes(ng)
    def time_to(x):
        dx=x-p;xx=p[...,None]+dx[...,None]*z**2
        d=2*xx*(potential(xx)-eps)-L2[...,None]
        if np.any(d<=0):raise RuntimeError('negative radial radicand')
        return TIME*np.sum(w*dx[...,None]*z/np.sqrt(d),axis=-1)
    t_r=time_to(X);t_t=time_to(rt**2)
    delay=np.stack((t_t-t_r,t_t+t_r),axis=-1)
    ecc=np.sqrt(np.maximum(0,1-2*eps*L2/H.Qt**2))
    E=np.arccos(np.clip((1-2*eps*rt/H.Qt)/ecc,-1,1));M=E-ecc*np.sin(E)
    dwell=period(H.Qt,eps)*(1-M/math.pi)
    upper_delay=2*rt*TIME/math.sqrt(e0-2*eps)
    if np.min(delay)<-1e-9 or np.max(delay)>upper_delay*(1+1e-6):raise RuntimeError('delay bound')
    if np.min(dwell)<exterior_min(H.Qt,rt,eps)*(1-1e-10):raise RuntimeError('exterior bound')
    return delay,dwell[...,None],ecc,E,M

def run(H,n=16,nr=48,ne=20,na=10,ng=24,npot=1500):
    e0=H.Qt/H.rt;emax=e0/2;rt=H.rt
    x,w=nodes(nr);r=rt*x*x;ww=w*2*rt*x*r*r;ps=H.psi(r)
    xx=(rt*np.linspace(0,1,npot)**2)**2
    pot=PchipInterpolator(xx,e0+H.psi(np.sqrt(xx)))
    rho=V**2/(4*math.pi*6.6743e-11*R**2);pref=.5*rho*.125*math.sqrt(2)*V*GYR/(H.Ic**2*H.Qt)
    es,we=nodes(ne);mu,wm=nodes(na)
    ages=np.array([1.,5.,10.,13.8]);N=np.zeros((len(ages),4));En=np.zeros_like(N)
    radii=np.array([1.,1.5,2.,3.,5.,10.,20.])*rt;outer=np.zeros((len(ages),len(radii)))
    for ep,weight in zip(emax*es,emax*we):
        mask=ps>e0-ep;rr=r[mask]
        S=pref*ww[mask]*np.array([density(H,p,e0-ep,n) for p in ps[mask]])
        delay,dwell,ecc,E,M=trajectories(H,rr,ep,mu,pot,ng)
        apo=H.Qt/(2*ep)*(1+ecc)
        anom=np.arccos(np.clip((1-2*ep*radii[:,None,None]/H.Qt)/ecc,-1,1))
        tR=period(H.Qt,ep)/(2*math.pi)*(anom-ecc*np.sin(anom)-M)
        duration=np.maximum(0,dwell[...,0]-2*tR)
        duration=np.where(radii[:,None,None]<=apo,duration,0.)
        W=weight*S[:,None,None]*wm[None,:,None]/2
        for i,T in enumerate(ages):
            inside=np.minimum(T,delay);outside=np.clip(T-delay,0,dwell);returned=np.maximum(T-delay-dwell,0)
            ledger=np.array([T*weight*S.sum(),np.sum(W*inside),np.sum(W*outside),np.sum(W*returned)])
            N[i]+=ledger;En[i]-=ep*ledger
            outer[i]+=np.sum(W[None,...]*np.clip(T-delay[None,...]-tR[...,None],0,duration[...,None]),axis=(1,2,3))
    gross=pref*(ww@np.array([rate(H,p,e0-emax,n)[1]-rate(H,p,e0,n)[1] for p in ps]))
    assert np.allclose(outer[:,0],N[:,2],rtol=1e-10,atol=1e-15) and np.all(np.diff(outer,axis=1)<=1e-15)
    rows=[]
    for t,(a,b,c,d),energy,mass in zip(ages,N,En,outer):
        front=(rt**1.5+1.5*math.sqrt(2*H.Qt)*t/TIME)**(2/3)
        assert np.all(mass[radii>front]<=1e-15)
        rows.append(dict(age_Gyr=float(t),born=float(a),inside=float(b),outside=float(c),returned=float(d),ledger_rel=float(abs(a-b-c-d)/a),source_closure=float(abs(a/t/gross-1)),energy_ledger=energy.tolist(),energy_ledger_rel=float(abs(energy[0]-energy[1:].sum())/abs(energy[0])),outside_radii_fraction=mass.tolist(),front_kpc=.6*front))
    return dict(model=H.model,n=n,nr=nr,ne=ne,na=na,ng=ng,npot=npot,source_cumulative=gross,radii_over_rt=(radii/rt).tolist(),rows=rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--n',type=int,default=16);p.add_argument('--nr',type=int,default=48);p.add_argument('--ne',type=int,default=20);p.add_argument('--na',type=int,default=10);p.add_argument('--ng',type=int,default=24);p.add_argument('--npot',type=int,default=1500);a=p.parse_args()
    print(json.dumps(run(Halo(model=a.model),a.n,a.nr,a.ne,a.na,a.ng,a.npot)))
