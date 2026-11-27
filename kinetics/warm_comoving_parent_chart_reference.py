"""Positive closed prompt events on a fixed comoving parent momentum lattice."""

import numpy as np
from scipy.special import roots_legendre

from joint_prompt_two_moment_evolution_reference import Cells
from three_species_decay_onsager_reference import Es,M,P,m,mp,width


def lab_gates(p):
    """Lab parent endpoints at which daughter-bin geometry changes."""
    e=np.hypot(p,mp)
    f=np.hypot(p,m)
    lim=e[-1]+f[-1]
    assert lim>M
    def cross(x,y,mass,rest):
        v=y*rest+x*P
        return np.r_[M*(rest*rest+x*x)/v,M*v/mass**2]
    return np.unique(np.clip(np.r_[M,lim,cross(p,e,mp,Es),cross(p,f,m,M-Es),
                                  (e[:,None]+f).ravel(),M+np.geomspace(1e-9,lim-M,46)],M,lim))


def common_rule(comoving_daughter_edges, scales, order=4):
    """Union the geometry gates at the locked collision scale factors.

    The resulting quadrature measures comoving number, not lab density.
    Keeping Q, dN_c and f fixed transports without remapping or bin shocks.
    """
    p=np.asarray(comoving_daughter_edges)
    points=[]
    for b in scales:
        walls=lab_gates(p/b)
        points.append(b*np.sqrt((walls-M)*(walls+M)))
    q=np.unique(np.concatenate(points))
    assert q[0]==0. and np.all(np.diff(q)>0.)
    # Nearly identical gates can create zero-width panels in binary64.
    # Coalesce only an absolute-plus-relative roundoff neighborhood.
    tol=8*np.finfo(float).eps*np.maximum(q[1:],1e-3)
    keep=np.r_[np.diff(q)>tol,True]
    keep[0]=True
    q=q[keep]
    z,w=roots_legendre(order)
    d=np.diff(q)[:,None]/2
    Q=(q[:-1,None]+(z+1.)*d).ravel()
    dw=(w*d).ravel()
    assert np.all(Q>0.) and np.all(dw>0.)
    return Q,2/np.pi**2*Q**2*dw


class ComovingCollision:
    """Expose the established closed two-daughter moment chart interface."""
    def __init__(self, q_daughter, scale, Q, dN, angle=4):
        p=np.asarray(q_daughter)/scale
        self.c=[Cells(p,mp,1),Cells(p,m,4)]
        self.E=np.hypot(M,Q/scale)
        self.n=len(Q)
        self.comoving_Q=np.asarray(Q)
        self.comoving_measure=np.asarray(dN)
        self.scale=scale
        e0,e1=[np.hypot(p,c.mass) for c in self.c]
        pp=Q/scale
        span=2.*P*pp/M
        self.w=(dN/scale**3)*(width*M/self.E)/span
        z,w=roots_legendre(angle)
        parts=[[] for _ in range(6)]
        for i in range(0,self.n,160):
            E=self.E[i:i+160]
            dp=pp[i:i+160]
            high=(E*Es+P*dp)/M
            low=(Es**2+(mp*dp/M)**2)/high
            ends=np.concatenate([low[:,None],high[:,None],np.broadcast_to(e0,(len(E),len(e0))),
                                 E[:,None]-e1],axis=1)
            ends=np.sort(np.clip(ends,low[:,None],high[:,None]),axis=1)
            l,r=ends[:,:-1],ends[:,1:]
            mid=l+(r-l)/2
            Em=E[:,None]-mid
            ok=(r>l)&(mid>=e0[0])&(mid<e0[-1])&(Em>=e1[0])&(Em<e1[-1])
            a,b=np.nonzero(ok)
            lab=(l[a,b,None]+(r[a,b,None]-l[a,b,None])*(z+1.)/2).ravel()
            iy=np.repeat(i+a,len(z))
            parts[0].append(iy)
            parts[1].append(np.repeat(np.searchsorted(e0,mid[a,b],side='right')-1,len(z)))
            parts[2].append(np.repeat(np.searchsorted(e1,Em[a,b],side='right')-1,len(z)))
            parts[3].append(lab)
            parts[4].append(self.E[iy]-lab)
            parts[5].append(((r[a,b,None]-l[a,b,None])*w/2).ravel())
        self.y,self.I,self.J,self.ph,self.ch,self.pre=map(np.concatenate,parts)
        self.fac=self.w[self.y]
        self.k1=self.ph-mp
        self.k2=self.ch-m
        assert np.all((0<=self.I)&(self.I<len(p)-1))
        assert np.all((0<=self.J)&(self.J<len(p)-1))
        fullspan=np.bincount(self.y,weights=self.pre,minlength=self.n)
        assert np.all(fullspan<=span+16*np.finfo(float).eps*self.E)
        assert len(self.y)>0 and np.all(self.pre>0.)
        assert np.all(self.E[self.y]<=(e0[-1]+e1[-1])*(1.+5e-15))


def remapless_parent_work(Q,dN,phase,a,b,order=20):
    """Compute a⁴ρ, pressure trace work and invariant comoving number."""
    assert 0.<a<=b and np.all((0.<=phase)&(phase<=1.))
    z,w=roots_legendre(order)
    x=np.log1p((b-a)/a)
    s=a*np.exp(x*(z+1.)/2)
    E=np.hypot(M,Q[:,None]/s)
    n=dN*phase
    E0,E1=np.hypot(M,Q/a),np.hypot(M,Q/b)
    start=a*(n@E0)
    end=b*(n@E1)
    work=x/2*np.sum((s*w)*(n@(M*M/E)))
    # A stable independent endpoint difference. Unlike direct subtraction,
    # this formula resolves tiny expansion intervals and excess energies.
    ea=np.hypot(a*M,Q);eb=np.hypot(b*M,Q)
    stable=(b-a)*(b+a)*M*M*(n@(1./(ea+eb)))
    assert abs((end-start)-work)<3e-13*max(end,1e-275)
    if work>0.:
        assert abs(stable/work-1.)<5e-8
    # Separately check the actual energy per comoving volume, whose
    # redshift work has the opposite sign and a different integrand.
    kinetic_work=-x/2*np.sum(w*(n@((Q[:,None]/s)**2/E)))
    kinetic_exact=np.expm1(-2*x)*(n@((Q/a)**2/(E0+E1)))
    if kinetic_work != 0.:
        assert abs(kinetic_exact/kinetic_work-1.)<2e-13
    return dict(comoving_count=float(n.sum()),a4energy_start=float(start),a4energy_end=float(end),
                pressure_trace_work=float(work),stable_energy_difference=float(stable),
                a3energy_redshift_quadrature=float(kinetic_work),a3energy_redshift_exact=float(kinetic_exact))
