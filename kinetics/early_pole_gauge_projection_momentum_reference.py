"""Quantify vector momentum error of a hard energy-only radial interpolation."""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_snapshot_gauge_ledger_reference import cell_J,photon


def curve(p,m,A):
    v=np.hypot(p,m)
    return v+A/v


def bracket(ed,target,m,A,mode):
    a,b=ed[:-1],ed[1:]
    p=(a+b)/2 if mode=='linear' else np.cbrt((a**3+b**3)/2)
    q=curve(p,m,A)
    j=np.searchsorted(q,target)-1
    assert j.min()>=0 and j.max()<len(p)-1
    w=(target-q[j])/(q[j+1]-q[j])
    v=p[j]+w*(p[j+1]-p[j])
    e=(1-w)*q[j]+w*q[j+1]
    return v,max(abs(e-target))


def check(d,g,Q,mode,cut=.5):
    x=d['x_end'];M=x;B=d['mass_ratio']*x
    S=d['checkpoints'][-1]['scalar_pole_over_T']
    Jg=photon(1.,0.,cut)
    Jp,Jb=[cell_J(d[f'endpoint_{lab}_momentum_edges'],
                  d[f'endpoint_{lab}_occupations'],m,cut)
           for lab,m in [('parent',M),('fermion',B)]]
    AP,AB=.5*g*g/np.pi**2*(Jg+np.array([Jp,Jb]))
    es=np.hypot(S,Q)
    p1=brentq(lambda p:curve(p,M,AP)-es-curve(p+Q,B,AB),.2,4.,xtol=5e-14)
    p2=brentq(lambda p:curve(p,M,AP)-es-curve(p-Q,B,AB),.2,4.,xtol=5e-14)
    # Split when either nodal moment interpolation switches. High-order
    # angular-aperture fraction computed as union of root-delimited pieces.
    ph=np.asarray(d['endpoint_parent_momentum_edges'])
    bh=np.asarray(d['endpoint_fermion_momentum_edges'])
    pnodes=.5*(ph[:-1]+ph[1:]) if mode=='linear' else np.cbrt((ph[:-1]**3+ph[1:]**3)/2)
    bnodes=.5*(bh[:-1]+bh[1:]) if mode=='linear' else np.cbrt((bh[:-1]**3+bh[1:]**3)/2)
    q=curve(bnodes,B,AB)+es
    D=(q+np.sqrt(q*q-4*AP))/2
    v=np.sqrt(np.maximum((D-M)*(D+M),0.))
    nodes=np.unique(np.r_[p1,pnodes[(pnodes>p1)&(pnodes<p2)],v[(v>p1)&(v<p2)],p2])
    t,w=roots_legendre(36)
    h=np.diff(nodes)[:,None]/2
    p=((nodes[:-1]+nodes[1:])[:,None]/2+h*t).ravel()
    weights=(h*w).ravel()/(p2-p1)
    E=curve(p,M,AP)
    y=E-es
    enc=(y+np.sqrt(y*y-4*AB))/2
    k=np.sqrt((enc-B)*(enc+B))
    mu=(p*p+Q*Q-k*k)/(2*p*Q)
    assert max(abs(mu))<1.
    p1m,ef=bracket(ph,E,M,AP,mode)
    p2m,eb=bracket(bh,y,B,AB,mode)
    bdir=(p-Q*mu)/k
    defect=np.hypot((p1m-p)-(p2m-k)*bdir,(p2m-k)*Q*np.sqrt(1-mu*mu)/k)/Q
    assert max(ef,eb)<3e-14

    def geom(z):
        E=curve(np.array([z]),M,AP)
        X,_=bracket(ph,E,M,AP,mode)
        Y,_=bracket(bh,E-es,B,AB,mode)
        return abs((X-Y)[0])-Q
    fails=0.
    for a,b in zip(nodes[:-1],nodes[1:]):
        # On an interpolation branch X-Y affine in E, so its two signed
        # threshold crossings may be solved independently.
        def linear_diff(z):
            q=curve(np.array([z]),M,AP)
            return (bracket(ph,q,M,AP,mode)[0]-bracket(bh,q-es,B,AB,mode)[0])[0]
        d1,d2=linear_diff(a),linear_diff(b)
        breaks=[a,b]
        for z in [Q,-Q]:
            if (d1-z)*(d2-z)<0:
                breaks.append(brentq(lambda p:linear_diff(p)-z,a,b,xtol=3e-14))
        cuts=sorted(breaks)
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            if geom((lo+hi)/2)>0.: fails+=hi-lo
    return dict(g=g,q=Q,mode=mode,max_energy_defect=max(ef,eb),
                vector_defect_over_q_mean=weights@defect,
                vector_defect_over_q_max=np.max(defect),
                aperture_length_fraction_no_repaired_triangle=fails/(p2-p1))


def run():
    rows=[]
    for grid in ('fine','third'):
        d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
        for mode in ('linear','volume'):
            for g in (0.,.005):
                for q in (1e-6,4.1e-6,1e-5):
                    rows.append(dict(grid=grid,**check(d,g,q,mode)))
    return dict(scope='Actual frozen aperture and representative radial nodes in archived bins; '
                      'not a direct measurement of production solver vector momentum error.',
                rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
