"""Pauli-safe, number-conservative parent remap in comoving momentum q=p/T.

Piecewise-constant fermion occupations are remapped by exact overlap in
u=q**3, in O(Nold+Nnew). The resulting convex averages preserve number
and Pauli bounds; discrete coarse-grained entropy cannot decrease.
The remap does NOT conserve the massive-particle energy moment exactly,
so it is a Liouville discretization benchmark, NOT a complete FRW solver.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre, expit, xlogy


def edges_checked(values):
    x=np.asarray(values,dtype=float)
    if x.ndim!=1 or len(x)<2 or x[0]<0 or not np.all(np.isfinite(x)) or np.any(np.diff(x)<=0):
        raise ValueError("Momentum edges must increase strictly from a nonnegative finite start")
    return x


def number(q_edges,f,g=4):
    q=edges_checked(q_edges)
    f=np.asarray(f,dtype=float)
    if f.shape!=(len(q)-1,) or not np.all(np.isfinite(f)) or g<=0:
        raise ValueError("Invalid occupation/degeneracy for number moment")
    return float(g/(6*np.pi**2)*np.dot(np.diff(q**3),f))


def rebin(old_edges,new_edges,old_f):
    old=edges_checked(old_edges);new=edges_checked(new_edges)
    f=np.asarray(old_f,dtype=float)
    if f.shape!=(len(old)-1,) or not np.all(np.isfinite(f)) or np.any(f<0) or np.any(f>1):
        raise ValueError("Invalid Pauli occupations")
    if not (np.isclose(old[0],new[0],atol=1e-13,rtol=1e-13)
            and np.isclose(old[-1],new[-1],atol=1e-13,rtol=1e-13)):
        raise ValueError("Shared global momentum domain required")
    u=old**3;v=new**3
    deposition=np.zeros(len(new)-1)
    i=j=0
    while i<len(f) and j<len(deposition):
        low=max(u[i],v[j]); high=min(u[i+1],v[j+1])
        if high>low:deposition[j]+=f[i]*(high-low)
        if u[i+1]<=v[j+1]:i+=1
        else:j+=1
    result=deposition/np.diff(v)
    if abs(float(np.dot(f,np.diff(u))-np.dot(result,np.diff(v))))>3e-12*max(1,float(np.dot(f,np.diff(u)))):
        raise ArithmeticError("Number conservation failed")
    if min(result)<-1e-13 or max(result)>1+1e-13:
        raise ArithmeticError("Pauli violation")
    return result


def entropy(q_edges,f):
    q=edges_checked(q_edges);f=np.asarray(f,dtype=float)
    if f.shape!=(len(q)-1,) or np.any(f<0) or np.any(f>1):
        raise ValueError("Invalid fermionic entropy occupations")
    return float(-np.dot(np.diff(q**3),xlogy(f,f)+xlogy(1-f,1-f))/3)


def energy(q_edges,f,x,gauss_n=16):
    q=edges_checked(q_edges);f=np.asarray(f,dtype=float)
    if f.shape!=(len(q)-1,) or x<=0 or gauss_n<3:
        raise ValueError("Invalid physical-energy inputs")
    z,w=roots_legendre(gauss_n)
    a=q[:-1,None];b=q[1:,None]
    samples=.5*(a+b)+.5*(b-a)*z
    v=.5*(b-a)*w*samples**2
    return float(4/(2*np.pi**2)*np.sum(f[:,None]*v*np.hypot(samples,x)))


def thermal_shell_averages(q_edges,x0,gauss_n=16):
    q=edges_checked(q_edges)
    if x0<=0:raise ValueError("positive initial x required")
    z,w=roots_legendre(gauss_n)
    a=q[:-1,None];b=q[1:,None]
    samples=.5*(a+b)+.5*(b-a)*z
    weighted=.5*(b-a)*w*samples**2
    return np.sum(weighted*expit(-np.hypot(samples,x0)),axis=1)/(np.diff(q**3)/3)


def moving_parent_edges(parent_E,x,qmax):
    E=np.asarray(parent_E,dtype=float)
    if E.ndim!=1 or len(E)<2 or np.any(np.diff(E)<=0):
        raise ValueError("Parent energies must be ordered")
    if x<=0 or x>=E[0] or qmax<=0:
        raise ValueError("Parent center below rest mass: rebuild resonant graph")
    mids=.5*(E[:-1]+E[1:])
    qmid=np.sqrt(mids*mids-x*x)
    if qmid[-1]>=qmax:raise ValueError("Global comoving momentum domain too short")
    return np.r_[0.,qmid,qmax]


def diagnose(old,new,old_f,x):
    mapped=rebin(old,new,old_f)
    n0=number(old,old_f);n1=number(new,mapped)
    E0=energy(old,old_f,x);E1=energy(new,mapped,x)
    s0=entropy(old,old_f);s1=entropy(new,mapped)
    if s1<s0-1e-10:raise ArithmeticError("Numerical entropy decreased")
    return {"number_before":n0,"number_after":n1,
            "absolute_number_error":n1-n0,
            "relative_energy_regrid_error":E1/E0-1,
            "numerical_fermion_entropy_increase":s1-s0,
            "min_fermion_occupation":float(mapped.min()),
            "max_fermion_occupation":float(mapped.max())}


def run():
    from adaptive_resonant_lattice import AdaptiveClosedDecay
    from resonant_graph_repair import bridged_indices
    model=AdaptiveClosedDecay(bridged_indices(max_index=300),M1_T=.02)
    qmax=model.parent_E[-1]+.2
    old=moving_parent_edges(model.parent_E,.02,qmax)
    occupations=thermal_shell_averages(old,.02)
    rows=[]
    for x in (.02005,.0201,.0202):
        new=moving_parent_edges(model.parent_E,x,qmax)
        rows.append({"initial_x":.02,"later_x":x,**diagnose(old,new,occupations,x)})
    return {"scope":"Pauli- and number-conservative q^3 overlap, not energy-conserving",
            "model":"moving parent energy bins, shared comoving momentum domain",
            "rows":rows,"CLASS_ready":False,
            "warning":"Artificial coarse-graining produces entropy and physical energy error; no coupled FRW collisions or frozen primordial PSD."}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/parent_q_remap.json")
    option=p.parse_args()
    out=run()
    path=Path(option.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))
