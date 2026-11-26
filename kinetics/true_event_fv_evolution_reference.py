"""Conditional fixed-temperature true-event evolution from shell numbers.

Uses finite on-shell two-mass cut-cell quadrature and positive timestep
rejection. Parent is prescribed; event and represented energy are different
ledgers. Above-capacity Bose shapes are diagnostics, not condensates.
"""
import argparse
import json
import numpy as np
from scipy.special import expit, logit, roots_legendre
from massive_fv_fixed_bath_reference import bounds, kinematics
from massive_shell_moments import shell_weights


def energy_mesh(n,M,ratios,order):
    m,b=np.array(ratios)*M
    K,D=kinematics(M,m,b)
    q=np.geomspace(2e-5,24.,n+1)
    e,r=np.hypot(q,m),np.hypot(q,b)
    if min(np.diff(e).min(),np.diff(r).min())<=0:
        raise ValueError("collapsed energy cell")
    l,h=np.broadcast_arrays(e[:-1,None],r[None,:-1])
    u,v=np.broadcast_arrays(e[1:,None],r[None,1:])
    z=bounds(h,b,m,K,D)+bounds(v,b,m,K,D)
    pts=np.sort(np.stack([l,u,*[np.clip(a,l,u) for a in z]]),axis=0)
    nodes,w=roots_legendre(order)
    blocks=[]
    for a,c in zip(pts[:-1],pts[1:]):
        i,j=np.nonzero(c>a)
        low,delta=a[i,j],c[i,j]-a[i,j]
        for x,wx in zip(nodes,w):
            y=low+delta*(x+1)/2
            r0,r1=bounds(y,m,b,K,D)
            s,t=np.maximum(h[i,j],r0),np.minimum(v[i,j],r1)
            iw=t>s
            if not np.any(iw):
                continue
            for x2,w2 in zip(nodes,w):
                rr=s[iw]+(t[iw]-s[iw])*(x2+1)/2
                wt=delta[iw]*wx*(t[iw]-s[iw])*w2/4
                blocks.append((i[iw],j[iw],y[iw],rr,wt))
    i,j,y,z,weight=(np.concatenate(v) for v in zip(*blocks))
    weight*=4/(2*np.pi**2)*M/(2*D/M**2)
    return q,m,b,(i,j,y,z,weight)


class Reconstruction:
    def __init__(self,q,m,bose):
        wgt=1 if bose else 4
        N,E,_,_=shell_weights(q,m)
        self.number=N/(2*np.pi**2)*wgt
        self.mean=E/N
        self.mass,self.bose=m,bose
        z,w=roots_legendre(80)
        l,r=q[:-1,None],q[1:,None]
        k=l+(r-l)*(z+1)/2
        self.k=k*k/(np.hypot(k,m)+m)
        self.w=((r-l)/2)*w*k*k/N[:,None]
        self.lo=q[:-1]**2/(np.hypot(q[:-1],m)+m)
        self.hi=q[1:]**2/(np.hypot(q[1:],m)+m)
        self.cap=np.sum(self.w/np.expm1(self.k),axis=1) if bose else np.ones(len(N))

    def average(self,shift):
        d=self.k-np.asarray(shift)[...,None]
        with np.errstate(over="ignore"):
            f=1/np.expm1(d) if self.bose else expit(-d)
        return np.sum(self.w*f,axis=1)

    def solve(self,n):
        target=n/self.number
        if np.any(target<0) or (not self.bose and np.any(target>1+3.e-16)):
            raise ValueError("outside number/occupancy domain")
        # Avoid evaluating infinite root brackets at the exact endpoints.
        active=(target>0)&(target<self.cap)
        ta=target[active]
        k,w=self.k[active],self.w[active]
        if self.bose:
            lo=self.lo[active]-np.log1p(1/ta)
            hi=np.zeros(len(lo))
        else:
            lo=self.lo[active]+logit(ta)
            hi=self.hi[active]+logit(ta)
        for _ in range(45):
            mid=(lo+hi)/2
            d=k-mid[:,None]
            with np.errstate(over="ignore"):
                f=1/np.expm1(d) if self.bose else expit(-d)
            above=(w*f).sum(1)>ta
            lo,hi=np.where(above,lo,mid),np.where(above,mid,hi)
        s=np.zeros_like(target)
        s[active]=(lo+hi)/2
        return target,s

    def at(self,target,s,e,i,mode):
        f=target[i]
        if mode=="constant":
            return f
        below=f<self.cap[i]
        d=(e-self.mass)-s[i]
        if not self.bose:
            return np.where(f==0,0.,np.where(below,expit(-d),1.))
        with np.errstate(over="ignore"):
            v=1/np.expm1(d)
            z=(f/self.cap[i])/np.expm1(e-self.mass)
        return np.where(f==0,0.,np.where(below,v,(f if mode=="flat-high" else z)))

    def represented_energy(self,N,target,s,mode):
        if mode=="constant":
            return float(self.mean@N)
        if self.bose:
            with np.errstate(over="ignore"):
                fit=1/np.expm1(self.k-s[:,None])
                amp=target[:,None]/self.cap[:,None]
            high=target[:,None] if mode=="flat-high" else amp/np.expm1(self.k)
            f=np.where(target[:,None]==0,0.,np.where(target[:,None]<self.cap[:,None],fit,high))
        else:
            f=np.where(target[:,None]==0,0.,np.where(target[:,None]<1,expit(s[:,None]-self.k),1.))
        return float(self.mass*N.sum()+self.number@(self.w*f*self.k).sum(1))


def model(n,M,order,mode):
    q,m,b,data=energy_mesh(n,M,(.08,.18),order)
    i,j,e,r,weight=data
    B,F=Reconstruction(q,m,True),Reconstruction(q,b,False)
    ebath=e+r
    mother=expit(-ebath)
    def averages(y):
        return B.solve(y[:n]),F.solve(y[n:2*n])
    def rhs(y):
        (nb,sb),(nf,sf)=averages(y)
        p,f=B.at(nb,sb,e,i,mode),F.at(nf,sf,r,j,mode)
        g=weight*(mother*(1-f)*(1+p)-(1-mother)*p*f)
        return np.r_[np.bincount(i,weights=g,minlength=n),np.bincount(j,weights=g,minlength=n),
                     e@g,r@g,-ebath@g]
    def physical_energy(y):
        (nb,sb),(nf,sf)=averages(y)
        return B.represented_energy(y[:n],nb,sb,mode)+F.represented_energy(y[n:2*n],nf,sf,mode)
    limits=np.r_[np.full(n,np.inf),F.number]
    eq=np.r_[B.average(-m)*B.number,F.average(-b)*F.number,0.,0.,0.]
    return rhs,physical_energy,limits,eq,B,F,weight.size,float(weight@mother)


def evolve(rhs,upper,initial,end,accuracy,scales):
    y=initial.copy()
    d=np.asarray(scales)
    # Scale driving error in the shell occupations, not by total abundance.
    t,h=0.,end/5
    count=np.zeros(3,int)
    def valid(a):
        return np.all((a[:-3]>=0)&(a[:-3]<=upper))
    while end-t>8*np.spacing(end):
        h=min(h,end-t)
        assert h>0 and t+h>t and count.sum()<20000
        z=y+h*rhs(y)
        if not valid(z):
            count[1]+=1
            h*=.5
            continue
        w=(y+z+h*rhs(z))/2
        if not valid(w):
            count[1]+=1
            h*=.5
            continue
        err=np.max(abs(w-z)/np.maximum(d,abs(y)))
        if err>accuracy:
            count[2]+=1
            h*=max(.2,.85*np.sqrt(accuracy/err))
            continue
        t+=h
        y=w
        count[0]+=1
        h*=min(2,max(.2,.85*np.sqrt(accuracy/max(err,1e-300))))
    return y,count.tolist()


def trial(M,n,order,mode,rate,target,tol):
    rhs,actual,upper,eq,B,F,nodes,vacuum=model(n,M,order,mode)
    eq_error=np.sum(abs(rhs(eq)[:n]))/vacuum
    initial=np.zeros(2*n+3)
    if target:
        initial[:n]=B.number*B.average(-B.mass)
        initial[n:2*n]=F.number*F.average(-F.mass)
        j=np.argmin(abs(B.mean-.08*M-.0002))
        initial[j]=target*B.number[j]*B.cap[j]
    y,cts=evolve(rhs,upper,initial,rate,tol,np.r_[B.number,F.number,1.,1.,1.])
    d=y-initial
    nf=np.sum(d[n:2*n]);nb=np.sum(d[:n])
    assert abs(nb-nf)<1e-10*max(abs(nb),abs(nf),1e-15)+2e-14*max(float(np.sum(initial[:-3])),float(np.sum(y[:-3])),1e-15)
    assert abs(sum(d[-3:]))<1e-10*max(sum(abs(d[-3:])),1e-15)
    if mode!="constant":
        assert eq_error<2e-10
    nr=np.sum(y[:n])
    return dict(M=M,cells=n,order=order,mode=mode,end=rate,accuracy=tol,
                quadrature_nodes=nodes,vacuum_rate=vacuum,eq_abs_scalar_number_source_over_vacuum_rate=float(eq_error),
                net_boson_number=float(nb),total_boson_number=float(nr),
                boson_shell_mean_energy=float(B.mean@y[:n]/nr),
                number_residual=float(abs(nb-nf)),event_energy_residual=float(abs(sum(d[-3:]))),
                represented_energy_defect=float(actual(y)-actual(initial)-sum(d[-3:-1])),
                represented_energy_defect_over_abs_event_gain=float(
                    abs(actual(y)-actual(initial)-sum(d[-3:-1]))/max(sum(abs(d[-3:-1])),1e-300)),
                timesteps=cts,above_capacity_shells=int(np.count_nonzero(y[:n]>B.number*B.cap)))


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--quick",action="store_true")
    a=ap.parse_args()
    grids=(12,24) if a.quick else (24,48,96)
    rows=[trial(.1,n,4,mode,1e-4,0.,1e-6) for n in grids
                          for mode in ("constant","sat-high")]
    stress=[trial(.1,grids[0],4,m,1e-4,10.,1e-7) for m in ("sat-high","flat-high")]
    q_refine=trial(.1,grids[0],6,"sat-high",1e-4,0.,1e-6)
    dt_refine=trial(.1,grids[0],4,"sat-high",1e-4,0.,1e-8)
    early=[trial(.005,n,4,"sat-high",1e-4,0.,1e-6) for n in grids[:2]]
    return dict(rows=rows,stress=stress,finer_cut_quadrature=q_refine,
                smaller_timestep_tolerance=dt_refine,early_mass_ratio_scan=early,
                status="conditional-fixed-bath-reference-not-FRW")
if __name__ == "__main__":
    print(json.dumps(main(),indent=2,allow_nan=False))
