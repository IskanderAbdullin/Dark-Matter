"""Two-mass fixed-bath finite-volume pilot; conditional on prescribed bath."""
import json
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import expit, roots_legendre
from massive_shell_moments import shell_weights
from two_mass_area import allowed_area


def kinematics(M,m,b):
    if not (M>m+b and min(m,b)>0):
        raise ValueError("this reference requires two finite, on-shell masses")
    K=(M*M-m*m-b*b)/2
    D=np.sqrt((M-m-b)*(M+m+b)*(M-m+b)*(M+m-b))/2
    return K,D


def bounds(fixed, mass, other, K, D):
    p=np.sqrt(np.maximum((fixed-mass)*(fixed+mass),0.))
    return other*other*fixed/(K+D)+D/(fixed+p),(K*fixed+D*p)/(mass*mass)


def integrated_area(e,r,M,m,b,order):
    K,D=kinematics(M,m,b)
    lo,hi=e[:-1,None],e[1:,None]
    rlo,rhi=r[None,:-1],r[None,1:]
    a,c=bounds(rlo,b,m,K,D)
    u,v=bounds(rhi,b,m,K,D)
    L,H=np.broadcast_arrays(lo+np.zeros_like(rlo),hi+np.zeros_like(rhi))
    pts=np.sort(np.stack([L,H,*[np.clip(y,L,H) for y in (a,c,u,v)]]),axis=0)
    out,qe,qr,ve,vr=(np.zeros_like(L) for _ in range(5))
    z,w=roots_legendre(order)
    for low,high in zip(pts[:-1],pts[1:]):
        d=high-low
        for x,y in zip(z,w):
            q=low+d*(x+1)/2
            f,g=bounds(q,m,b,K,D)
            l=np.minimum(rhi,np.maximum(rlo,f))
            h=np.maximum(l,np.minimum(rhi,g))
            dr=h-l
            k=np.exp(-q-l)
            wq=d*y/2
            vac=np.log1p(k*(-np.expm1(-dr))/(1+k*np.exp(-dr)))
            out+=wq*vac
            qe+=wq*vac*q
            ve+=wq*dr*q
            vr+=wq*dr*(l+h)/2
            weighted=np.zeros_like(out)
            for zr,wr in zip(z,w):
                rr=l+dr*(zr+1)/2
                weighted+=wr/2*rr*expit(-q-rr)
            qr+=wq*dr*weighted
    return out,qe,qr,ve,vr


def parent_vacuum_gain(e,r,M,m,b,moment=None):
    """Independent isotropic parent momentum and angle-fraction integral."""
    K,D=kinematics(M,m,b)
    rest=(M*M+m*m-b*b)/(2*M)
    k=D/M
    def integrand(p):
        E=np.hypot(p,M)
        mid,span=rest*E/M,k*p/M
        lower=max(e[0],E-r[-1],mid-span)
        upper=min(e[-1],E-r[0],mid+span)
        angle=max(upper-lower,0)/(2*span) if p>0 else 0.
        avg=(lower+upper)/2
        w=avg if moment=="boson" else E-avg if moment=="fermion" else 1.
        return w*4/(2*np.pi**2)*p*p*M/E*expit(-E)*angle
    emax=e[-1]+r[-1]
    pmax=np.sqrt((emax-M)*(emax+M))
    return quad(integrand,0.,pmax,points=[i for i in (.01,.04,.1,.3,1,2,4,8,16,32) if i<pmax],
                epsrel=2e-10,epsabs=1e-15,limit=500)[0]


class FixedBath:
    def __init__(self,n,M=.1,ratios=(.08,.18),qlo=2e-5,qhi=24.,order=12):
        m,b=(x*M for x in ratios)
        self.M,self.m,self.b=M,m,b
        K,D=kinematics(M,m,b)
        eta=2*D/(M*M)
        self.q=q=np.geomspace(qlo,qhi,n+1)
        en,fn=np.hypot(q,m),np.hypot(q,b)
        if min(np.diff(en).min(),np.diff(fn).min())<=0:
            raise ValueError("at least one energy shell collapsed; no silent clipping")
        n0,e0,_,t0=shell_weights(q,m)
        n1,e1,_,t1=shell_weights(q,b)
        self.vb,self.vf=n0/(2*np.pi**2),4*n1/(2*np.pi**2)
        self.eb,self.ef=e0/n0,e1/n1
        self.tb,self.tf=t0/n0,t1/n1
        pref=4/(2*np.pi**2)*M/eta
        ar=allowed_area(en[:-1,None],en[1:,None],fn[None,:-1],fn[None,1:],M,m,b)
        ai,qe,qr,ve,vr=integrated_area(en,fn,M,m,b,order)
        self.A_energy=np.array([pref*qe,pref*qr])
        self.R_energy=np.array([pref*ve,pref*vr])
        self.R=pref*ar
        self.Ai=pref*ai
        self.A=self.R*expit(-self.eb[:,None]-self.ef[None,:])
        self.vacuum=parent_vacuum_gain(en,fn,M,m,b)
        self.vacuum_q=np.array([parent_vacuum_gain(en,fn,M,m,b,k) for k in ("boson","fermion")])
        ix,iy=self.eb[:,None],self.ef[None,:]
        shell_p=np.sqrt((ix-m)*(ix+m))*np.sqrt((iy-b)*(iy+b))
        nonphysical=abs(ix*iy-K)>shell_p
        self.off_center=float(self.Ai[nonphysical].sum()/self.Ai.sum())
        # Relative against the largest matrix element, not an ill-conditioned
        # element-by-element ratio for cells of almost zero measure.
        assert np.max(self.Ai-self.R)<5e-11*np.max(self.R) and np.min(self.Ai)>=0.

    def rhs(self,y,A):
        p,f=np.split(y,2)
        F=A*(1+p[:,None])*(1-f[None,:])-(self.R-A)*p[:,None]*f[None,:]
        return np.r_[F.sum(1)/self.vb,F.sum(0)/self.vf]

    def jac(self,y,A):
        p,f=np.split(y,2)
        U=A-self.R*f[None,:]
        V=-A-self.R*p[:,None]
        return np.block([[np.diag(U.sum(1)/self.vb),V/self.vb[:,None]],
                         [U.T/self.vf[:,None],np.diag(V.sum(0)/self.vf)]])

    def run(self,bath,tend=1e-4,rtol=5e-8):
        A=self.A if bath=="balanced" else self.Ai
        eq=np.r_[1/np.expm1(self.eb),expit(-self.ef)]
        balance=np.abs(self.rhs(eq,A)[:len(self.eb)])@self.vb/self.vacuum
        sol=solve_ivp(lambda t,y:self.rhs(y,A),(0,tend),np.zeros_like(eq),
                      jac=lambda t,y:self.jac(y,A),method="BDF",
                      rtol=rtol,atol=1e-11,max_step=tend/12,t_eval=[tend])
        if not sol.success:
            raise RuntimeError(sol.message)
        p,f=np.split(sol.y[:,-1],2)
        b=self.vb*p
        n=b.sum()
        assert p.min()>=-1e-12 and f.min()>=-1e-12 and f.max()<=1+1e-12
        assert abs(n-self.vf@f)<2e-9*n
        ab=(1+p[:,None]-f[None,:])
        prod=p[:,None]*f[None,:]
        # Physical event-energy sources versus the actual discrete rate of
        # stored shell energy in the integrated (not balanced-center) model.
        physical=np.sum(self.A_energy*ab-self.R_energy*prod,axis=(1,2))
        dy=self.rhs(sol.y[:,-1],A)
        shell=np.array([self.eb@(self.vb*dy[:len(p)]),self.ef@(self.vf*dy[len(p):])])
        return dict(final_shell_energy_source_rel_defect_if_integrated=(
                        None if bath=="balanced" else float(np.sum(abs(shell-physical))/np.sum(abs(physical)))),
                    vacuum_gain_rel_error=float(abs(A.sum()/self.vacuum-1)),
                    equilibrium_net_over_vacuum_gain=float(balance),n=float(n),
                    boson_mean_energy=float(b@self.eb/n),
                    boson_number_difference_rel=float(abs(n-self.vf@f)/n),
                    fermion_max=float(f.max()),boson_max=float(p.max()),
                    daughter_trace_over_energy=float((b@self.tb+(self.vf*f)@self.tf)/
                                                     (b@self.eb+(self.vf*f)@self.ef)),
                    nfev=sol.nfev)


def report():
    out=[]
    settings=[(.1,(.08,.18),(32,64,128,256)),
              (.005,(.08,.18),(64,128,256)),
              (.1,(8.6e-7,.2),(64,128,256))]
    for x,ratios,models in settings:
        for n in models:
            m=FixedBath(n,M=x,ratios=ratios)
            r=m.run("integrated")
            b=m.run("balanced")
            gainq=np.sum(m.A_energy,axis=(1,2))
            qerr=float(np.max(abs(gainq/m.vacuum_q-1)))
            assert r["vacuum_gain_rel_error"]<2e-7 and qerr < 3e-7
            assert b["equilibrium_net_over_vacuum_gain"]<1e-12
            out.append(dict(M=x,daughter_mass_ratios=list(ratios),cells=n,vacuum_reference=m.vacuum,
                            vacuum_physical_energy_gain_by_species=m.vacuum_q.tolist(),
                            vacuum_physical_energy_gain_max_rel_error=qerr,
                            integrated_gain_in_cells_with_forbidden_center=m.off_center,
                            integrated=r,balanced=b))
    return dict(qmin=2e-5,qmax=24.,
                time_in_rest_width_units=1e-4,rows=out)

if __name__ == "__main__":
    print(json.dumps(report(),indent=2,allow_nan=False))
