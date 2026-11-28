"""Measure instantaneous collision derivatives in fixed-total-energy physical densities."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, root

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_combined_tabulated_modes import sample
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_tabulated_entropy_split_expansion import Stage, hubble
from sm_dark_entropy_exchange_reference import StandardModel


def fixed_energy(s, center):
    Ts,eta,l = center
    d0 = s.bath([eta,l])
    scale = np.array([d0['n'],d0['rho']])
    rho = s.eos.state(Ts)[0]+scale[1]
    H = hubble(rho,0.)

    def restore(v):
        n,r=scale*v
        if min(n,r,rho-r)<=0.:
            raise ArithmeticError('Perturbed physical densities must stay positive')
        ts=brentq(lambda t:s.eos.state(t)[0]/(rho-r)-1.,s.lo,s.hi,xtol=5e-15)
        def f(a):
            d=s.bath(a)
            B,E,C=(d[k] for k in ('dn_deta','chemical_energy','cv_fixed_number'))
            td=np.exp(a[1])
            M=np.array([[B,B*E/td],[B*E,td*C+B*E*E/td]]) / np.array([n,r])[:,None]
            return np.array([d['n']/n-1.,d['rho']/r-1.]),M
        q=root(lambda a:f(a)[0],[eta,l],jac=lambda a:f(a)[1],tol=2e-12).x
        if max(abs(f(q)[0]))>4e-13 or abs((s.eos.state(ts)[0]+s.bath(q)['rho'])/rho-1.)>4e-13:
            raise ArithmeticError('Perturbation leaves the fixed-energy contour')
        return np.r_[ts,q]
    if max(abs(restore(np.ones(2))-center))>4e-13:
        raise ArithmeticError('Contour inversion does not preserve its exact center')
    def vel(v):
        ts,e,x=restore(v)
        q,_=s.collision(ts,np.exp(x),e)
        return -q/(H*scale)
    return vel,rho,H


def jac(f,h):
    v=np.ones(2)
    return np.column_stack([(f(v-2*d)-8*f(v-d)+8*f(v+d)-f(v+2*d))/(12*h) for d in h*np.eye(2)])


def eigen(M):
    q=np.linalg.eigvals(-M)
    return np.array(sorted(q,key=lambda k:(k.real,k.imag)))


def measure(s, state, mesh, tag):
    f,rho,H=fixed_energy(s,state)
    Te=brentq(lambda t:(s.eos.state(t)[0]+s.bath([0.,np.log(t)])['rho'])/rho-1.,
              s.lo,s.hi,xtol=5e-15)
    g,R,J=fixed_energy(s,[Te,0.,np.log(Te)])
    if abs(R/rho-1.)>4e-13 or abs(J/H-1.)>4e-13:
        raise ArithmeticError('Equilibrium comparator must have the same total energy and H')
    F,G=[(jac(u,6e-4),jac(u,3e-4)) for u in (f,g)]
    ef,eg=eigen(F[1]),eigen(G[1])
    conv=max(np.max(abs(eigen(y[0])-eigen(y[1]))/np.maximum(abs(v),1e-20)) for y,v in ((F,ef),(G,eg)))
    eq=np.array(sample(s.p,mesh,Te,s.eos.state(Te)[3])[tag+'_modes'])/H
    eq_error=np.max(abs(eg/eq-1.))
    if conv>8e-6 or eq_error>8e-6:
        raise ArithmeticError(('Jacobian not resolved or exact-equilibrium spectrum mismatch',conv,eq_error))
    return Te,ef,eg,F[1],G[1],conv,eq_error


def pair_set(a,b,c):
    for r in a['rows']:
        T,e,beta=[r[k] for k in ('local_crossing_T','actual_eta','actual_beta')]
        yield r['operator'],r['mode']+'_reference_temperature',r['log_a'],[T,e,np.log(T/(1.-beta))]
    for r in b['rows']:
        T,e=r['T_SM'],r['eta_dark']
        yield r['operator'],'equal_temperatures',r['x'],[T,e,np.log(T)]
    for r in c['rows']:
        p=r['points'][-1]
        yield r['operator'],'direct_x096',p[0],p[1:4]


def run(eos_path, clock, equal, direct):
    b=[p.read_bytes() for p in (clock,equal,direct)]
    A,B,C=map(json.loads,b)
    h=hashlib.sha256(eos_path.read_bytes()).hexdigest()
    tail=hashlib.sha256(b[2]).hexdigest()
    if (any(d['table_input_sha256']!=h for d in (A,B,C)) or B['input_sha256']!=tail
            or A['prior_continuation_and_clocks_sha256'][1]!=tail
            or [(r['operator'],r['mode']) for r in A['rows']]!=[(o,m) for o in ('timelike','combined') for m in ('slow','fast')]
            or any([r['operator'] for r in d['rows']]!=['timelike','combined'] for d in (B,C))):
        raise ArithmeticError('The exact direct archives and common EOS must match')
    e=StandardModel(eos_path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for tag,label,x,u in pair_set(A,B,C):
        K={'timelike':Stage,'combined':CombinedExpansion}[tag]
        a,b0=[measure(K(e,i,n),u,m,tag) for i,n,m in ((p,112,(18,29.)),(q,168,(25,36.)))]
        mesh=max(max(abs(a[i]-b0[i])/np.maximum(abs(b0[i]),1e-20)) for i in (1,2))
        if mesh>3e-5:
            raise ArithmeticError('Distinct phase-space meshes do not resolve the physical-density spectrum')
        T,F,G,M,N,n,err=b0
        rows.append(dict(operator=tag,label=label,log_a=x,state=list(u),same_energy_equilibrium_T=T,
                         live_rate_over_H_real=F.real.tolist(),live_rate_over_H_imaginary=F.imag.tolist(),
                         equilibrium_rate_over_H=G.real.tolist(),
                         live_real_over_equilibrium_minus_one=(F.real/G.real-1.).tolist(),live_jacobian=M.tolist(),
                         equilibrium_jacobian=N.tolist(),relative_defects=[n,mesh,err]))
    return dict(rows=rows,relative_defect_fields=['step_difference','phase_space_mesh','independent_equilibrium_spectrum'],
                variables=['relative_dark_number_density','relative_dark_energy_density'],
                input_sha256=[hashlib.sha256(v).hexdigest() for v in b],table_input_sha256=h,
                scope='Instantaneous Jacobian of the collision-only field at fixed total energy in physical-density coordinates; not the expanding variational propagator or a decoupling criterion.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','clock','equal','direct'):
        ap.add_argument(k,type=Path)
    print(json.dumps(run(*vars(ap.parse_args()).values()),indent=2,allow_nan=False))
