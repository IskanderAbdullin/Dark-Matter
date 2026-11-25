"""Positive three-scale numerical representation of one collisionless scalar PSD.

The partition is an exact identity, not a physical three-species model.
Additional tracer resolves the intermediate momentum shoulder which is missed
by a two-scale quadrature at extreme mass degeneracy.
"""
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator,CubicSpline
from scipy.integrate import trapezoid


def export_three(row,folder,tag='fine',cold_multiple=8.,middle_factor=1.,points=7000,energy_peak=False):
    raw=row['frozen_cell_spectrum'];edges=np.asarray(raw['edges']);f=np.asarray(raw['f_phi'])
    l,h=edges[:-1],edges[1:];cent=np.sqrt(l*h)
    assert np.min(f)>=-1e-7
    f=np.maximum(f,0.)
    moments=[float(np.sum(f*(h**(k+3)-l**(k+3))/(2*np.pi**2*(k+3)))) for k in (0,1,2)]
    assert abs(moments[0]/row['final']['number']-1)<1e-8
    peak=f*cent**3
    cold=float(cent[np.argmax(peak)])
    mask=cent>30*cold
    assert np.any(mask)
    hot=float(cent[mask][np.argmax((peak*cent**2 if energy_peak else peak)[mask])])
    s1=cold_multiple*cold
    s2=middle_factor*np.sqrt(s1*hot)
    assert s1<s2<hot,(s1,s2,hot)
    # Center each numerical quadrature on the number-weighted peak of its
    # partition, with a distinct scale for the broad middle shoulder.
    def weights(q):
        a=(q/s1)**8;b=(q/s2)**8
        w0=np.exp(-a);r=-np.expm1(-a)
        return [w0,r*np.exp(-b),r*(-np.expm1(-b))]
    middle=float(cent[np.argmax(peak*weights(cent)[1])])
    scales=[cold,middle,hot]
    interp=PchipInterpolator(np.log(cent),np.log(np.maximum(f,1e-100)),extrapolate=False)
    def evaluate(q):
        out=np.full_like(q,1e-100)
        within=(q>=cent[0])&(q<=cent[-1])
        out[within]=np.exp(interp(np.log(q[within])))
        above=q>cent[-1]
        slope=min(-.01,(np.log(max(f[-1],1e-100))-np.log(max(f[-2],1e-100)))/(cent[-1]-cent[-2]))
        out[above]=np.maximum(f[-1],1e-100)*np.exp(slope*(q[above]-cent[-1]))
        return out
    comps=[]
    for i,(label,scale) in enumerate(zip(('cold','middle','hot'),scales)):
        limit=max(60.,(4*s1 if i==0 else 4*s2 if i==1 else edges[-1]+5)/scale)
        u=np.unique(np.r_[np.geomspace(1e-7,.1,points//3),np.geomspace(.1,limit,2*points//3)])
        q=scale*u
        # A strictly positive negligible floor prevents 0/0 in CLASS's
        # logarithmic derivative and extrapolation, including large-u files.
        occupation=np.maximum(evaluate(q)*weights(q)[i],1e-100*np.exp(-u/50.))
        assert np.all(np.isfinite(occupation)) and np.min(occupation)>1e-290
        psd=occupation/(2*np.pi)**3
        file=Path(folder)/(tag+'_'+label+'.dat')
        np.savetxt(file,np.column_stack((u,psd)),fmt='%.17e')
        probe=(u[:-1]+u[1:])/2;cv=CubicSpline(u,psd)(probe)
        neg=float(4*np.pi*scale**3*trapezoid(probe**2*np.maximum(-cv,0.),probe))
        m=[float(4*np.pi*scale**(k+3)*trapezoid(u**(k+2)*psd,u)) for k in (0,1,2)]
        comps.append(dict(label=label,file=str(file),scale=scale,points=len(u),number=m[0],energy=m[1],q2=m[2],spline_negative_number=neg))
    sums=[sum(c[k] for c in comps) for k in ('number','energy','q2')]
    errors=[float(a/b-1) for a,b in zip(sums,moments)]
    assert max(map(abs,errors))<.015,errors
    assert sum(c['spline_negative_number'] for c in comps)<moments[0]*1e-9
    return dict(raw_number=moments[0],raw_energy=moments[1],raw_q2=moments[2],splits=[s1,s2],
                components=comps,moment_relative_errors=errors,
                reconstructed_hot_number_fraction=sum(c['number'] for c in comps[1:])/sums[0],
                energy_peak=energy_peak,method='positive log-PCHIP with three smooth partitions; no abundance rescaling')
