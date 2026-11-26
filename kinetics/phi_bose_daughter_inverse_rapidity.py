"""Physical χ2 daughter-energy-first Yukawa inverse-rapidity mode operator.

For a FIXED massive daughter χ2 energy, the allowed parent momentum is
generally a two-sided finite interval, not [p_threshold,pmax]. Integrating
p/E1 dp = dk (scalar daughter energy) with positive log-k quadrature resolves
Bose IR. χ1 and φ occupations are held fixed for the positive Pauli mode step;
not a full self-consistent three-species strong Boltzmann evolution.
"""
import math
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import expit
from kinetics.phi_bose_physical_cutcell_weak import (
    M1,M2,MPHI,G1,KSTAR,WIDTH,
)

E2STAR=(M1*M1+M2*M2-MPHI*MPHI)/(2*M1)
XI2=math.asinh(KSTAR/M2)


def parent_bounds(e2,pmax):
    """Exact inverse physical daughter support, from daughter rapidity z.

    Parent y=asinh(p/M1), decay daughter rest-frame ξ=asinh(k*/m2),
    lab daughter z=asinh(p2/m2). Allowed |ξ-y|<=z<=ξ+y iff
    |ξ-z|<=y<=ξ+z. Intersect both bounds with 0<=p<=pmax.
    """
    if e2<M2:return None
    z=math.asinh(math.sqrt(max(0,e2*e2-M2*M2))/M2)
    yp=math.asinh(pmax/M1)
    ymin=abs(XI2-z)
    if ymin>=yp:return None
    ymax=min(yp,XI2+z)
    return M1*math.sinh(ymin),M1*math.sinh(ymax)


def daughter_coeff(e2,pmax,nparent=24,T1=.9,Tphi=1.,
                   mu1=0.,muphi=0.,g2=4):
    """G,L,G_E1,L_E1; exact Pauli daughter source dot f2=G(1−f2)−Lf2.

    The physical χ2 energy DOS per volume is g2 E2 sqrt(E2²−m2²)/(2π²).
    The inclusive initial g1 Γrest state counting and independent g2 must
    be fixed consistently with the underlying spin channels.
    """
    assert e2>M2 and g2>0 and T1>0 and Tphi>0
    bounds=parent_bounds(e2,pmax)
    if bounds is None:return 0.,0.,0.,0.
    a,b=bounds
    # Exact monotonic scalar k=E1-e2 turns p dp/E1 into dk. Integrate
    # on log(k) to retain stimulated Bose small-energy contributions.
    kl=math.hypot(a,M1)-e2
    kh=math.hypot(b,M1)-e2
    if kl<MPHI*(1-1e-7):
        raise ArithmeticError("Out of physical daughter energy shell")
    kl=max(kl,MPHI)
    x,w=leggauss(nparent)
    ln=math.log(kh/kl)
    k=kl*np.exp((x+1)*ln/2)
    dk=w*k*ln/2
    e1=e2+k
    f1=expit(-(e1-mu1)/T1)
    bose=1/np.expm1((k-muphi)/Tphi)
    Gparts=f1*(1+bose)
    Lparts=(1-f1)*bose
    p2=math.sqrt((e2-M2)*(e2+M2))
    pref=G1*WIDTH*M1*M1/(2*g2*KSTAR*e2*p2)
    return tuple(float(pref*np.dot(dk,v)) for v in (
        Gparts,Lparts,e1*Gparts,e1*Lparts))


def exact_daughter_step(e2,pmax,dt,f20,nparent=24,
                        T1=.9,Tphi=1.,mu1=0.,muphi=0.):
    """Exact 0<=f2(t)<=1 conditional on fixed positive f1,bφ baths."""
    assert 0<=f20<=1 and dt>=0
    G,L,_,_=daughter_coeff(e2,pmax,nparent,T1,Tphi,mu1,muphi)
    if G+L==0:return f20
    equilibrium=G/(G+L)
    return equilibrium+(f20-equilibrium)*math.exp(-(G+L)*dt)


def _z_support_cuts(pmax):
    yp=math.asinh(pmax/M1)
    zmin=max(0.,XI2-yp)
    zmax=XI2+yp
    return sorted(set([zmin,zmax]+[z for z in
                      (XI2,abs(yp-XI2)) if zmin<z<zmax]))


def _rapidity_weak_integrand(z,pmax,nparent,T1,T2,Tphi,mu1,mu2,muphi):
    e2=M2*math.cosh(z)
    p2=M2*math.sinh(z)
    G,L,G1,L1=daughter_coeff(e2,pmax,nparent,T1,Tphi,mu1,muphi)
    f2=expit(-(e2-mu2)/T2)
    J=G*(1-f2)-L*f2
    J1=G1*(1-f2)-L1*f2
    # p2 physical momentum, de2=p2 dz, DOS2=g2 e2 p2/(2π²).
    weight=4*e2*p2*p2/(2*math.pi**2)
    return weight*np.array([J,J1-e2*J,-J1,e2*J])


def integrate_daughter_weak(pmax,nz=32,nparent=32,T1=.9,T2=1.,
                            Tphi=1.,mu1=0.,mu2=0.,muphi=0.):
    """Daughter-rapidity Gauss first, log scalar energy + parent second."""
    cuts=_z_support_cuts(pmax)
    x,w=leggauss(nz)
    total=np.zeros(4)
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        zs=lo+(hi-lo)*(x+1)/2
        dz=w*(hi-lo)/2
        for z,ww in zip(zs,dz):
            total+=ww*_rapidity_weak_integrand(
                z,pmax,nparent,T1,T2,Tphi,mu1,mu2,muphi)
    return dict(Rphi=float(total[0]),Qphi=float(total[1]),
                Q1=float(total[2]),Q2=float(total[3]),
                energy_defect=float(sum(total[1:])))


def integrate_daughter_adaptive(pmax,nparent=36,T1=.9,T2=1.,
                               Tphi=1.,mu1=0.,mu2=0.,muphi=0.):
    """Independent adaptive outer daughter rapidity near support kinks."""
    cuts=_z_support_cuts(pmax)
    vals=[sum(quad(lambda z:_rapidity_weak_integrand(
                   z,pmax,nparent,T1,T2,Tphi,mu1,mu2,muphi)[i],
                   a,b,epsabs=1e-21,epsrel=2e-9,limit=180)[0]
              for a,b in zip(cuts[:-1],cuts[1:]))
          for i in range(4)]
    return dict(Rphi=vals[0],Qphi=vals[1],Q1=vals[2],
                Q2=vals[3],energy_defect=sum(vals[1:]))


def audit():
    from kinetics.phi_bose_physical_cutcell_weak import analytic_weak
    rows=[]
    for P in (2.,8.,16.):
        ref=analytic_weak(P)
        for nz in (8,16,32,64):
            v=integrate_daughter_weak(P,nz,nz)
            rows.append(dict(P=P,mesh=nz,
                             error=max(abs(v[k]/ref[k]-1)
                               for k in ("Rphi","Qphi","Q1","Q2"))))
        v=integrate_daughter_adaptive(P)
        rows.append(dict(P=P,mesh="adaptive outer",
                         error=max(abs(v[k]/ref[k]-1)
                           for k in ("Rphi","Qphi","Q1","Q2"))))
    return dict(schema_version=1,
                scope="conditional fixed χ1,φ baths and weak 3-energy reference, NOT full coupled strong kinetics",
                rows=rows,class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2))
