"""Vectorized on-shell two-mass daughter-energy finite-volume area.

The result is a geometric dE dr area, not a rate. Energies/masses share units.
Widths below float64 resolution require a different energy grid; zero-width
cells return zero. No quasiparticle dispersion or off-shell physics is assumed.
"""
import numpy as np


def _bounds(r, b, m, K, D):
    """Energy of the other daughter at fixed r, rationalized lower root."""
    p = np.sqrt(np.maximum((r-b)*(r+b), 0))
    with np.errstate(divide='ignore', invalid='ignore'):
        lo = m*m*r/(K+D) + D/(r+p)
        hi = (K*r+D*p)/(b*b) if b else np.full_like(r, np.inf)
    return lo, hi


def _integrals(lo, hi, b):
    """Integrals of r, 1/(r+p), and p on [lo, hi]; b>0."""
    dr = hi-lo
    pl = np.sqrt(np.maximum((lo-b)*(lo+b), 0))
    ph = np.sqrt(np.maximum((hi-b)*(hi+b), 0))
    dp = np.divide(dr*(hi+lo), ph+pl, out=np.zeros_like(dr), where=ph+pl>0)
    dt = np.log1p((dr+dp)/(lo+pl))
    # Rapidity t: integral dr/(r+p) = t/2 + exp(-2t)/4.
    s = -2*dt
    rem = np.where(abs(s)<1e-3, s*s*(.5+s*(1/6+s*(1/24+s/120))), np.expm1(s)-s)
    e = (b/(lo+pl))**2
    inv = .5*(1-e)*dt+.25*e*rem
    ir = dr*(hi+lo)/2
    ip = ir-b*b*inv
    # Near rest, use a non-cancelling rapidity form for integral p dr.
    near = ph < .1*b
    if np.any(near):
        d = dt[near]; sh = np.sinh(d); ch = np.cosh(d)
        corr = np.where(abs(d)<1e-3, d**3*(2/3+d*d*(2/15+d*d*4/315)), sh*ch-d)
        ip[near] = pl[near]**2*sh*ch + lo[near]*pl[near]*sh**2 + b*b/2*corr
    return ir, inv, ip


def _massless(e0,e1,r0,r1,K):
    # E >= K/(2r). Partial region followed by a fully allowed region.
    with np.errstate(divide='ignore'):
        a = np.maximum(r0, K/(2*e1)); split = K/(2*e0)
    c = np.maximum(a, np.minimum(r1, split)); a = np.minimum(a, r1)
    c = np.minimum(c, r1)
    dr = c-a
    log = np.log1p(np.divide(dr,a,out=np.zeros_like(dr),where=a>0))
    return np.clip(e1*dr-K/2*log+(e1-e0)*np.maximum(r1-np.maximum(r0,split),0),0,(e1-e0)*(r1-r0))


def allowed_area(elo,ehi,rlo,rhi,x,m,b):
    """Broadcastable exact area for decay x -> m + b, with x > m+b.

    Fast full/empty classification; split mixed rectangles at all crossings
    of the two energy boundaries. Integrate analytic primitives in rapidity.
    """
    if min(m,b)<0 or x<=0:
        raise ValueError('invalid masses')
    e0,e1,r0,r1 = np.broadcast_arrays(*[np.asarray(a,dtype=float) for a in (elo,ehi,rlo,rhi)])
    if np.any(e0<m) or np.any(r0<b) or np.any(e1<e0) or np.any(r1<r0):
        raise ValueError('invalid daughter energy cell')
    if x<=m+b:
        return np.zeros(e0.shape)
    # Use the larger mass as the integration variable for numerical stability.
    if m>b:
        return allowed_area(r0,r1,e0,e1,x,b,m)
    K=(x*x-m*m-b*b)/2
    if b==0:
        return _massless(e0,e1,r0,r1,K)
    # Factorized K-m*b near threshold.
    D=np.sqrt(max((x-m-b)*(x+m+b)*(x-m+b)*(x+m-b),0))/2
    lo0,hi0=_bounds(r0,b,m,K,D); lo1,hi1=_bounds(r1,b,m,K,D)
    low=np.minimum(lo0,lo1)
    if m:
        low=np.where((r0<=K/m)&(K/m<=r1),m,low)
    rect=(e1-e0)*(r1-r0)
    empty=(e1<=low)|(e0>=hi1)|(rect==0)
    full=(e0>=np.maximum(lo0,lo1))&(e1<=hi0)&~empty
    out=np.where(full,rect,0.)
    mask=~(full|empty)
    if not np.any(mask):
        return out
    # Extended precision on mixed cells, including heights tiny relative
    # to their energy: avoid cancellation in constant-minus-bound integrals.
    a,c,l,h=[v[mask].astype(np.longdouble) for v in (e0,e1,r0,r1)]
    xx,mm,bb=map(np.longdouble,(x,m,b))
    kk=(xx*xx-mm*mm-bb*bb)/2
    dd=np.sqrt((xx-mm-bb)*(xx+mm+bb)*(xx-mm+bb)*(xx+mm-bb))/2
    u,v=_bounds(a,mm,bb,kk,dd); w,z=_bounds(c,mm,bb,kk,dd)
    pts=np.sort(np.stack([l,h,np.clip(u,l,h),np.clip(v,l,h),np.clip(w,l,h),np.clip(z,l,h)]),axis=0)
    acc=np.zeros_like(a)
    for L,H in zip(pts[:-1],pts[1:]):
        mid=(L+H)/2
        low,high=_bounds(mid,bb,mm,kk,dd)
        live=(H>L)&(np.minimum(c,high)>np.maximum(a,low))
        if not np.any(live):
            continue
        LL,HH=L[live],H[live]; A,C=a[live],c[live]
        dr=HH-LL; ir,inv,ip=_integrals(LL,HH,bb)
        il=mm*mm/(kk+dd)*ir+dd*inv
        ih=(kk*ir+dd*ip)/(bb*bb)
        upper=C<=high[live]; lower=A>=low[live]
        area=np.where(upper,np.where(lower,(C-A)*dr,C*dr-il),np.where(lower,ih-A*dr,2*dd*ip/(bb*bb)))
        acc[live]+=np.clip(area,0,(C-A)*dr)
    out[mask]=np.minimum(np.maximum(acc,0),rect[mask])
    return out
