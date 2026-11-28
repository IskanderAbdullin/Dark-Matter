"""Exact area and trace controls for an omitted near-rest strip."""
import json

import numpy as np

from two_mass_kinetic_offset_trace_reference import DT, geometry, width


def sinhm(x):
    if abs(x)<.02:
        y=x*x
        return x*y*(DT(1)/6+y*(DT(1)/120+y*(DT(1)/5040+y*(DT(1)/362880+y/39916800))))
    return np.sinh(x)-x


def atanm(x):
    if abs(x)<.02:
        y=x*x
        return x*y*(DT(1)/3+y*(-DT(1)/5+y*(DT(1)/7+y*(-DT(1)/9+y*(DT(1)/11-y/13)))))
    return x-np.arctan(x)


def strip(e,m,b,p):
    e,m,b,p=map(DT,(e,m,b,p))

    K=m*b+e*(2*(m+b)+e)/2
    D=np.sqrt(e*(2*(m+b)+e)*(2*m+e)*(2*b+e))/2
    Y=2*np.arcsinh(np.sqrt(e*(2*(m+b)+e)/(4*m*b)))
    u=np.arcsinh(p/m)
    if not 0<u<Y:
        raise ValueError("The comparison uses a contained rest-side strip.")
    area=D*sinhm(2*u)/2
    trace=2*D*m*atanm(p/m)
    kinetic_edges=[0.,p]
    # Pad the other daughter interval by half the width to avoid trying
    # to represent exact near-coincident partner endpoints in float64.
    partner=b*np.sinh(np.array([Y-1.5*u,Y+1.5*u]))
    if partner[0]<0 or np.hypot(b,partner[0])==np.hypot(b,partner[1]):
        raise ValueError("Choose representable partner bounds.")
    actual=geometry(e,m,b,kinetic_edges,partner,64)
    rect=width(kinetic_edges,m)*width(partner,b)
    area_ref=actual[0]*rect
    traces=actual[1:3]*rect
    err=float(abs(DT(area_ref)/area-1))
    terr=float(abs(DT(traces[0])/trace-1))
    if max(err,terr)>2e-8:
        raise ValueError("Analytic strip moments did not reproduce the offset quadrature.")
    leading=2*D*(p/m)**3/3
    # Integrate p^2 / sqrt(m^2+p^2) to rigorously trap the leading area.
    lower=leading/np.hypot(DT(1),p/m)
    if not lower<=area*(DT(1)+DT(2e-17)) or not area<=leading*(DT(1)+DT(2e-17)):
        raise ValueError(f"Cubic near-rest area bound was violated at p/m={p/m}: A/L={area/leading}.")
    # The trace weight of the second daughter is bounded by its
    # extremal shell energies on the entire strip.
    tr=(traces[1]/area_ref)
    global_min=b/np.cosh(Y+u)
    global_max=b/np.cosh(Y-u)
    if not global_min<=tr<=global_max:
        raise ValueError("Partner trace escaped its shell bound.")
    return dict(p_over_scalar_mass=float(p/m),area=float(area),leading_cubic_area=float(leading),
                area_vs_leading=float(area/leading),scalar_trace=float(trace),
                partner_mean_trace=tr,limiting_partner_mean_trace=float(m*b*b/K),
                independent_area_relative_error=err,independent_scalar_trace_relative_error=terr,
                naive_scalar_absolute_bin_collapses=bool(float(np.hypot(m,p))==float(m)))


def run():
    grid=(1e-2,1e-4,1e-6,1e-8,1e-10)
    out=[]
    for e,m,b in ((.031,.065,.018),(.08,.024,.102)):
        rows=[strip(e,m,b,r*m) for r in grid]
        if not any(r['naive_scalar_absolute_bin_collapses'] for r in rows):
            raise ValueError("No lost absolute bin exercised.")
        out.append(dict(masses=[m,b],parent_excess=e,results=rows))
    return dict(scope="Exact geometric omitted-strip control; no collision kernel assumed.",
                cases=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
