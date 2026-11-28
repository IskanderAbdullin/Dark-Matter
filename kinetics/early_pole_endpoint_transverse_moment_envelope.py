"""Bound the frozen transverse pole using only each charge cell's N and energy."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,root
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import G,Species


def local_bound(m,edges,beta):
    E=np.hypot(m,edges)
    umax=(m/E[:-1])**2
    zmax=1-(m/E[1:])**2
    sums=np.zeros_like(umax)
    cutoff=32
    for n in range(cutoff+1):
        A=2+(2*n+1)*umax
        f=A*zmax**n
        fp=(2*n+1)*zmax**n
        fpp=np.zeros_like(f)
        if n:
            fp+=n*A*zmax**(n-1)
            fpp+=2*n*(2*n+1)*zmax**(n-1)
        if n>1:
            fpp+=n*(n-1)*A*zmax**(n-2)
        coeff=beta**(2*n)/((2*n+1)*(2*n+3))
        sums+=coeff*(2*f+10*umax*fp+4*umax**2*fpp)
    first=128*(cutoff+2)**3*beta**(2*(cutoff+1))
    ratio=beta**2*((cutoff+3)/(cutoff+2))**3
    if ratio>=1:
        raise ArithmeticError('Need a gap-induced bound on beta before summing the tail')
    # The per-cell interpolation error is bounded by curvature/8.
    return (sums+first/(1-ratio))*np.diff(E)**2/(8*E[:-1]**3)


def signed_bound(m,edges,beta,segments=12,derivative=False):
    # Bound the residual relative to the *original* cell secant.
    # Subdivide only in this analytic estimate, without splitting the
    # unknown charge occupation or introducing extra known moments.
    E=np.hypot(m,edges)
    l,h=E[:-1],E[1:]
    t=np.linspace(0.,1.,segments+1)
    ep=l[:,None]+(h-l)[:,None]*t
    u=(m/ep)**2
    z=1-u
    v=u[:,:-1]
    q=z[:,1:]
    small=np.diff(ep,axis=1)**2/8
    ml=np.zeros_like(l)
    mh=np.zeros_like(h)
    cut=32
    for n in range(cut+1):
        if derivative and n==0:
            continue
        p=(2*n+1)
        F=z**n*(2+p*u)/ep
        L=F[:,:1]*(1-t)+F[:,-1:]*t
        r=F-L
        A=2+p*v
        c=A*q**n
        c1=p*q**n
        c2=np.zeros_like(c)
        if n:
            c1+=n*A*q**(n-1)
            c2+=2*n*p*q**(n-1)
        if n>1:
            c2+=n*(n-1)*A*q**(n-2)
        e=(2*c+10*v*c1+4*v*v*c2)*small/ep[:,:-1]**3
        lo=np.min(np.minimum(r[:,:-1],r[:,1:])-e,axis=1)
        hi=np.max(np.maximum(r[:,:-1],r[:,1:])+e,axis=1)
        weight=2*n*beta**(2*n-1) if derivative else beta**(2*n)
        scale=weight/((2*n+1)*(2*n+3))
        if n:
            lo=np.minimum(lo,0.)
            hi=np.maximum(hi,0.)
        ml+=scale*lo
        mh+=scale*hi
    if derivative:
        first=256*(cut+2)**4*beta**(2*(cut+1)-1)
        r=beta**2*((cut+3)/(cut+2))**4
    else:
        first=128*(cut+2)**3*beta**(2*(cut+1))
        r=beta**2*((cut+3)/(cut+2))**3
    if r>=1:
        raise ArithmeticError('Need a positive gap to enclose upper beta')
    tail=first*np.diff(E)**2/(8*l**3*(1-r))
    return ml-tail,mh+tail


def enclosure(d,x,f):
    charged=[(x,d.p_edges,d.qp,d.wp*f[d.iP:]),
             (x*d.ratio,d.d_edges,d.qd,d.wd*f[d.ns:d.iP])]
    def raw(b):
        return sum(charge_field(p,N,m,b,44) for m,e,p,N in charged)
    def error(b):
        return G*G*fsum(
            fsum(local_bound(m,e,b)*N.reshape(-1,2).sum(1)) for m,e,p,N in charged)
    k=.0015
    Pi0=raw(0.)[0]
    e0=error(0.)
    lo=np.sqrt(k*k+Pi0-e0)
    if not lo>k:
        raise ArithmeticError('No strictly positive lower optical gap')
    beta_cap=k/lo
    err=error(beta_cap)
    shape=np.array([G*G*np.array([fsum(row*N.reshape(-1,2).sum(1)) for row in signed_bound(m,e,beta_cap)])
                    for m,e,p,N in charged]).sum(axis=0)
    small=4e-13*Pi0
    delta_low=max(-err,shape[0])-small
    delta_high=min(err,shape[1])+small
    w0=root(raw,k)[1]
    w=lambda x:(x-k)*(x+k)-raw(k/x)[0]
    if w(lo)+err>=0:
        lower=lo
    else:
        lower=brentq(lambda x:w(x)+err,lo,w0,xtol=5e-17)
    hi=.004
    while w(hi)<err:
        hi*=2
    upper=brentq(lambda x:w(x)-err,w0,hi,xtol=5e-17)
    def just(z):
        if w(lo)>=z:
            return lo
        return brentq(lambda x:w(x)-z,lo,hi,xtol=5e-17)
    A,B=just(delta_low),just(delta_high)
    # The derivative of the positive-power response series is itself
    # monotone in beta. Combine its interval and the frequency interval.
    dshape=np.array([G*G*np.array([fsum(z*N.reshape(-1,2).sum(1))
                                  for z in signed_bound(m,e,beta_cap,derivative=True)])
                     for m,e,p,N in charged]).sum(axis=0)
    dlow=max(0.,raw(k/B)[1]+dshape[0]-small)
    dhigh=raw(k/A)[1]+dshape[1]+small
    residue=[1/(2*B+k/A**2*dhigh),1/(2*A+k/B**2*dlow)]

    # Two independent conditional continuum closures are only checks; the
    # enclosure uses their charge moments but not their shapes.
    cats=[]
    fits=[]
    concentrated=[]
    for m,e,p,N in charged:
        V=4*np.diff(e**3)/(6*np.pi**2)
        cats.append(Species(m,e,N.reshape(-1,2).sum(1)/V))
        q,G0,*_=continuum_moments(e,m,N,36)
        fits.append((q,G0,m))
        mass=np.hypot(m,p)
        count=N.reshape(-1,2).sum(1)
        mean=(mass*N).reshape(-1,2).sum(1)/count
        concentrated.append((np.sqrt(np.maximum((mean-m)*(mean+m),0.)),count,m))
    _,moments,moment_residue=root(lambda z:sum(charge_field(*p,z,44) for p in fits),k)
    _,step,step_residue=root(lambda z:sum(np.array(c.field(z)) for c in cats),k)
    _,mean_root,mean_residue=root(lambda z:sum(charge_field(*p,z,44) for p in concentrated),k)
    # The step comparison does not generally preserve energy; a
    # two-moment fit must preserve it and fall in the derived interval.
    if not lower<=moments<=upper:
        raise ArithmeticError(('The moment-preserving control escaped the curvature interval',lower,moments,upper))
    if not A<=moments<=B or B-A>upper-lower or not residue[0]<=moment_residue<=residue[1]:
        raise ArithmeticError(('Signed interval did not enclose two-moment control',A,moments,B))
    if not A<=mean_root<=B or not residue[0]<=mean_residue<=residue[1]:
        raise ArithmeticError('Concentrated mean-energy stress fixture escaped')
    return dict(k=k,beta_upper=beta_cap,frozen_gap_squared_error=e0,
                absolute_field_defect_bound=err,node_frequency=w0,
                number_step_frequency=step,two_moment_frequency=moments,
                frequency_lower_upper=[lower,upper],relative_frequency_width=(upper-lower)/moments,
                tightened_field_defect_interval=[float(delta_low),float(delta_high)],
                tightened_frequency_enclosure=[A,B],tightened_relative_frequency_width=(B-A)/moments,
                residue_enclosure=residue,two_moment_residue=moment_residue,number_step_residue=step_residue,
                concentrated_at_mean_energy_frequency=mean_root,concentrated_at_mean_energy_residue=mean_residue,
                relative_residue_width=(residue[1]-residue[0])/moment_residue)


def run(cards):
    res=[]
    for g in ('fine','third'):
        s=load(g)[0]
        for n in (1,2,4,8):
            c=CellChart(s,SPECS[g],n)
            row=enclosure(c,s['x_end'],c.initial)
            if not row['frequency_lower_upper'][0]<=row['number_step_frequency']<=row['frequency_lower_upper'][1]:
                raise ArithmeticError('Exact original step must satisfy its moment envelope')
            if not row['tightened_frequency_enclosure'][0]<=row['number_step_frequency']<=row['tightened_frequency_enclosure'][1]:
                raise ArithmeticError('Exact original step must satisfy its signed moment envelope')
            if not row['residue_enclosure'][0]<=row['number_step_residue']<=row['residue_enclosure'][1]:
                raise ArithmeticError('Exact original step must satisfy its residue moment envelope')
            res.append(dict(grid=g,kind='fixed_original_step',subcells=n,**row))
    for x in cards:
        g=x['grid']
        s=load(g)[0]
        d=CellChart(s,SPECS[g],x['subcells_per_archived_cell'],reference_x=x['reference_x'])
        y=np.asarray(x['endpoint_occupation'])
        if x['counts']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,):
            raise ArithmeticError('Different endpoint state dimensions')
        res.append(dict(grid=g,kind='evolved_endpoint',subcells=x['subcells_per_archived_cell'],
                        **enclosure(d,x['x_end'],y)))
    return res


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',nargs='*',type=Path)
    args=p.parse_args()
    print(json.dumps(run([json.loads(f.read_text()) for f in args.states]),indent=2,allow_nan=False))
