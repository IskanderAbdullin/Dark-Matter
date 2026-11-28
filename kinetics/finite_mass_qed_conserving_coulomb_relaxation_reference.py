"""Prescribed conserving relaxation of a frozen density response; no physical width."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import root
from scipy.special import roots_legendre

from finite_mass_qed_frozen_transverse_plasmon_gate_reference import PREF, Species


def pi0(ch,z):
    out=0j
    if z==0:
        return -sum(c.D for c in ch)
    for c in ch:
        v=c.v
        x=v/z
        small=abs(x)<.2
        s=x[small]**2
        p=np.empty_like(x,dtype=complex)
        p[small]=sum(s**j/(2*j+1) for j in range(1,17))
        v=v[~small]
        p[~small]=z/(2*v)*(np.log(z+v)-np.log(z-v))-1
        out+=PREF*(c.W@p)
    return out


def pinu(ch,w,k,nu):
    if nu==0:
        return pi0(ch,w/k)
    base=pi0(ch,(w+1j*nu)/k)
    return (w+1j*nu)*base/(w+1j*nu*base/pi0(ch,0))


def direct(ch,w,k,nu,separate=False):
    """Solve mode angular population with an explicit conserved-density target."""
    nodes,weights=roots_legendre(144)
    A=[];B=[];N=[];jA=[];jB=[]
    z=(w+1j*nu)/k
    eta=nu/k
    for c in ch:
        V=c.v[:,None]*nodes
        U=(c.W[:,None]*weights/2)/(z-V)
        A.append(np.sum(U))
        B.append(np.sum(U*V))
        jA.append(np.sum(U*V))
        jB.append(np.sum(U*V*V))
        N.append(sum(c.W))
    A,B,N,jA,jB=map(np.asarray,(A,B,N,jA,jB))
    if separate:
        mu=B/(N-1j*eta*A)
    else:
        mu=sum(B)/sum(N-1j*eta*A)
    density=PREF*sum(B+1j*eta*mu*A)
    current=PREF*sum(jB+1j*eta*mu*jA)
    return density,current


def run(g,d):
    M=d['x_end']
    ch=[Species(m,d[f'endpoint_{tag}_momentum_edges'],d[f'endpoint_{tag}_occupations'])
        for m,tag in ((M,'parent'),(M*d['mass_ratio'],'fermion'))]
    gap=np.sqrt(sum(c.G for c in ch))
    rows=[]
    exact=[]
    continuity=[]; directdef=[]; closuregap=[]
    for k,w,nu in ((.001,.0004,.0002),(.001,.001,.00022),(.002,.0018,.0004)):
        rho,J=direct(ch,w,k,nu)
        P=pinu(ch,w,k,nu)
        directdef.append(abs(rho/P-1))
        continuity.append(abs(k*J/(w*rho)-1))
        R,S=direct(ch,w,k,nu,separate=True)
        closuregap.append(abs(R/P-1))
        directdef.append(abs(R/sum(pinu([c],w,k,nu) for c in ch)-1))
        continuity.append(abs(k*S/(w*R)-1))
    if max(continuity)>2e-12 or max(directdef)>2e-12:
        raise ArithmeticError((continuity,directdef))
    for nu in (.00001,.0001,.0004):
        for k in (1e-7,.0008,.00128,.00256,.004):
            static=pinu(ch,0.,k,nu)
            if abs(static/pi0(ch,0)-1)>5e-13:
                raise ArithmeticError('Conserving relaxation changed static Debye')
            samp=np.linspace(.002,1.8,200)*max(k,gap)
            Y=np.asarray([pinu(ch,w,k,nu) for w in samp])
            cut=-Y.imag
            if min(cut)<-1e-16:
                raise ArithmeticError('Negative spectral damping')
            ref=np.sqrt(gap**2+k*k)-.5j*nu
            def fn(z):
                v=1-pinu(ch,z[0]+1j*z[1],k,nu)/(k*k)
                return [v.real,v.imag]
            x=root(fn,[ref.real,ref.imag],tol=2e-11)
            rr=float(np.max(abs(np.asarray(fn(x.x)))))
            if rr>2e-9:
                raise ArithmeticError(('Controlled pole solve failed',g,nu,k,x.x,rr))
            if not -nu<x.x[1]<0:
                raise ArithmeticError(('Pole outside continued branch',x.x))
            if k==1e-7:
                v=np.sqrt(gap*gap-nu*nu/4)-1j*nu/2
                exact.append(abs(complex(*x.x)/v-1.))
            rows.append(dict(k=k,nu=nu,real=float(x.x[0]),minus_imag=float(-x.x[1]),residual=rr))
    sums=[]
    for r in rows:
        k,nu=r['k'],r['nu']
        if k not in (.0008,.00256,.004):
            continue
        match=r['real']/gap
        h=nu/gap
        joints=np.unique(np.r_[0.,max(k/gap,.05),match,np.maximum(0.,match+np.r_[-16.,-4.,-1.,1.,4.,16.]*h)])
        def weight(t):
            z=pinu(ch,t*gap,k,nu)
            return 2*k*k*(-z.imag)/(np.pi*t*((k*k-z.real)**2+z.imag*z.imag))
        def integral(eta,precision):
            def fn(t):
                return weight(t) if eta==0 else weight(t)*t*t/(t*t+eta*eta)
            bd=np.r_[joints,np.inf]
            return float(sum(quad(fn,a,b,epsabs=precision,epsrel=10*precision,limit=200)[0]
                             for a,b in zip(bd[:-1],bd[1:],strict=True)))
        static=-pi0(ch,0).real
        y=integral(0.,2e-11)
        z=integral(0.,5e-13)
        target=static/(k*k+static)
        err=max(abs(z-target),abs(y-z))
        e=[]
        for eta in (.1,1.,3.):
            R=pinu(ch,1j*eta*gap,k,nu)
            if abs(R.imag)>1e-17:raise ArithmeticError('Nonreal imaginary-axis response')
            e.append(abs(integral(eta,5e-13)-(-R.real/(k*k-R.real))))
        bd=np.r_[joints,np.inf]
        first=float(sum(quad(lambda t:weight(t)*t*t,a,b,epsabs=2e-11,epsrel=1e-10,limit=300)[0]
                        for a,b in zip(bd[:-1],bd[1:],strict=True)))
        if max(err,*e,abs(first-1.))>2e-8:
            raise ArithmeticError(('Conserving spectral normalization',k,nu,err,e,first))
        naive=pi0(ch,1j*nu/k)
        sums.append(dict(k=k,nu=nu,static_weight=z,static_expected=target,
                         static_absolute_defect=abs(z-target),mesh_defect=abs(y-z),imag_axis_defect=max(e),
                         weighted_high_frequency_sum=first,weighted_high_frequency_sum_defect=abs(first-1.),
                         naive_static_relative_error=abs(naive.real/(-static)-1)))
    return dict(grid=g,small_k_drude_defect=max(exact),kinetic_direct_defect=max(directdef),
                sum_checks=sums,
                continuity_defect=max(continuity),separate_vs_total_projection_example_max=max(closuregap),rows=rows)


if __name__=='__main__':
    p=Path(__file__).parent/'results'
    a=[]
    for g in ('fine','third'):
        with gzip.open(p/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:
            a.append(run(g,json.load(f)))
    print(json.dumps(a,indent=2))
