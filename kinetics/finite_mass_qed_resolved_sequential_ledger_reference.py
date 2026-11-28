"""Exactly solvable resolved two-step occupation, not a matched physical rate."""

import gzip
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import expit, logit

from finite_mass_qed_nonnull_photon_ward_poles_reference import setup
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q
from finite_mass_qed_spacelike_exact_statistical_aperture_reference import cell


def case(a,b,c,d,f0,t):
    """R1=a(1-f)-bf and R2=cf-d(1-f), for an internal fermionic mode."""
    H=a+b+c+d
    steady=(a+d)/H
    J=(a*c-b*d)/H
    scale=-np.expm1(-H*t)/H
    state=steady+(f0-steady)*np.exp(-H*t)
    first=J*t-(a+b)*(f0-steady)*scale
    second=J*t+(c+d)*(f0-steady)*scale
    return state,first,second,J


def audit(a,b,c,d,f0,t):
    def rhs(_, y):
        f=y[0]
        R1=a*(1-f)-b*f
        R2=c*f-d*(1-f)
        return [R1-R2,R1,R2]
    direct=solve_ivp(rhs,[0,t],[f0,0.,0.],method='DOP853',rtol=4e-13,atol=8e-16)
    f,L,R,J=case(a,b,c,d,f0,t)
    err=np.max(abs(np.array([f,L,R])-direct.y[:,-1]))
    affinity=np.log(a*c/(b*d))
    entropy=lambda f: -f*np.log(f)-(1-f)*np.log1p(-f)
    def local(x):
        f=case(a,b,c,d,f0,x)[0]
        r,s=rhs(x,[f])[1:]
        return (r*np.log(a*(1-f)/(b*f))+s*np.log(c*f/(d*(1-f))))
    budget=L*np.log(a/b)+R*np.log(c/d)+entropy(f)-entropy(f0)
    ent=abs(quad(local,0,t,epsabs=4e-13,epsrel=4e-13)[0]-budget)
    return float(err),float(ent),float(abs((L-R)-(f-f0))),float(J*affinity)


def rates(p,b,n,g,l,arm,rate_ratio):
    """Directional rates for the same forward cycle and bath data in either order.

    The coefficient sums are fixed to 1 for the charge-bath step and
    rate_ratio for the scalar step. Time units here are illustrative.
    """
    if not (min(p,1-p,b,1-b,n,g,l,rate_ratio)>0 and arm in (0,1)):
        raise ValueError("A paired diagnostic requires open forward and inverse support")
    def scale(x,y,s):
        h=s/(x+y)
        return h*x,h*y
    if arm==0:
        a,B=scale(p*l,(1-p)*g,1.)
        c,d=scale((1-b)*(1+n),b*n,rate_ratio)
    else:
        a,B=scale(p*(1+n),(1-p)*n,rate_ratio)
        c,d=scale((1-b)*l,b*g,1.)
    return a,B,c,d


def run(roots_dir):
    samples=[audit(a,b,c,d,f,t) for a,b,c,d in [(1.,.3,.7,.2),(.1,1.,.3,.7),(.35,.7,.8,.4),(.004,10.,.5,.23)]
             for f in (.03,.49,.97) for t in (.01,.2,3.)]
    x=np.array(samples)
    if np.any(x[:,:3] > 3e-12) or np.any(x[:,3]< -1e-15):
        raise ArithmeticError(x)
    receipt=json.loads((roots_dir/'finite-mass-qed-nonnull-photon-ward-poles-2026-11-27.json').read_text())
    shells=[]; thermal=[]; energy=[]; affinities=[]; memory=[]; identity=[]; conserved=[]; entropy=[]
    for source in receipt['data']:
        grid=source['grid']
        with gzip.open(roots_dir/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as fd:
            st=json.load(fd)
        masses=st['x_end'],st['x_end']*st['mass_ratio']
        def take(name,mag):
            edges=st['scalar_momentum_edges' if name=='scalar' else f'endpoint_{name}_momentum_edges']
            values=st[f'endpoint_{name}_occupations']
            return np.r_[0.,values,0.][np.searchsorted(edges,mag)]
        rows=[]
        for geo in source['checks']:
            if geo['vacuum_photon_invariant']>=0:
                continue
            k,w=geo['k'],geo['omega']
            c0={(.0008,.0002):0.,(.00128,.00064):-.5,(.00256,.00064):.5}[k,w]
            _,_,q,K,_,kin=setup(st,k,w,c0)
            G,L,D=np.sum([cell(m,k,w,np.asarray(st[f'endpoint_{tag}_momentum_edges']),
                               np.asarray(st[f'endpoint_{tag}_occupations']))
                          for m,tag in zip(masses,('parent','fermion'),strict=True)],axis=0)
            if min(G,L,D)<=0 or abs(G-L-D)>3e-13*max(G,L):
                raise ArithmeticError('Free cut does not supply the expected positive pair')
            n=take('scalar',Q)
            for az in geo['azimuths']:
                for pole in az['found_near_born_intermediate_on_shell_roots']:
                    P,C=kin(pole['parent_E'],az['azimuth'])
                    arm=pole['arm']
                    internal,ns=(P+K,'parent') if arm==0 else (C-K,'fermion')
                    vertex_a=(internal-P-K) if arm==0 else (P-internal-q)
                    vertex_b=(internal-C-q) if arm==0 else (internal+K-C)
                    energy.append(float(max(np.max(abs(vertex_a)),np.max(abs(vertex_b)))))
                    P0,C0,F0=(take(tag,np.linalg.norm(v[1:])) for tag,v in
                              (('parent',P),('fermion',C),(ns,internal)))
                    affinity=logit(P0)-logit(C0)+np.log(L/G)+np.log1p(1/n)
                    values=[]
                    for scale in (.01,1.,100.):
                        a,b,c,d=rates(P0,C0,n,G,L,arm,scale)
                        A=np.log(a/b)+np.log(c/d)
                        affinities.append(abs(A-affinity))
                        F,Lt,Rt,J=case(a,b,c,d,F0,1/(a+b+c+d))
                        identity.append(abs((Lt-Rt)-(F-F0)))
                        if J*affinity< -1e-15 or min(F,1-F)<=0.:
                            raise ArithmeticError('Sequential resolved sign or occupation failed')
                        scalar,current=(Rt,Lt) if arm==0 else (Lt,Rt)
                        memory.append(abs(scalar-J/(a+b+c+d)))
                        charged=np.array([-Lt,Rt,F-F0])
                        energies=np.array([P[0],C[0],internal[0]])
                        conserved.append(abs(energies @ charged + q[0]*scalar - w*current))
                        entropy.append(audit(a,b,c,d,F0,1/(a+b+c+d))[1])
                        values.append(dict(scalar_to_current_sum=scale,held_fraction=(a+d)/(a+b+c+d),
                                           stationary_scalar_current=J,
                                           scalar_count_over_one_held_relaxation=scalar,
                                           cut_current_over_one_held_relaxation=current,
                                           stored_internal_occupation_change=F-F0))
                    for T in (.4,1.,1.3):
                        mug=.12
                        fp,fb,fi=(expit((mug-E)/T) for E in (P[0],C[0],internal[0]))
                        n0=1/np.expm1(q[0]/T)
                        aa,bb,cc,dd=rates(fp,fb,n0,1.,np.exp(-w/T),arm,1.7)
                        thermal.extend((abs(aa/(aa+bb)-fi),abs(dd/(cc+dd)-fi),
                                        abs((aa*cc-bb*dd)/(aa+bb+cc+dd))))
                    rows.append(dict(k=k,omega=w,azimuth=az['azimuth'],arm=arm,parent_energy=P[0],
                                     actual_internal_occupation=F0,affinity=affinity,
                                     scalar_occupation=n,current_lesser_over_greater=L/G,diagnostics=values))
        if len(rows)!=24:
            raise ArithmeticError('A distinct occupied root disappeared')
        shells.append(dict(grid=grid,rows=rows))
    if max(thermal)>3e-14 or max(identity)>3e-14 or max(energy)>3e-14 or max(affinities)>3e-14 or max(conserved)>3e-14 or max(entropy)>3e-13:
        raise ArithmeticError((max(thermal),max(identity),max(energy),max(affinities),max(conserved),max(entropy)))
    return dict(scope="Normalized time diagnostic for a separately resolved fermion; not a sequential subtraction, integrated rate or physical regulator.",
                toy_cases=len(x),toy_state_integral_absolute_defect=x[:,0].max(),
                toy_entropy_absolute_defect=x[:,1].max(),toy_number_absolute_defect=x[:,2].max(),
                occupied_root_four_momentum_defect=max(energy),occupied_root_storage_defect=max(identity),
                occupied_root_energy_ledger_defect=max(conserved),occupied_root_integrated_entropy_defect=max(entropy),
                occupied_root_common_temperature_defect=max(thermal),cycle_affinity_defect=max(affinities),
                largest_memory_difference_at_one_held_relaxation=max(memory),data=shells)


if __name__=='__main__':
    data=run(Path(__file__).parent/'results')
    grids=data.pop('data')
    print(json.dumps(data,indent=2,allow_nan=False)[:-2]+',')
    print('  "data": [')
    for i,g in enumerate(grids):
        print('    {"grid": '+json.dumps(g['grid'])+', "rows": [')
        for j,row in enumerate(g['rows']):
            print('      '+json.dumps(row,allow_nan=False,separators=(',',':'))+
                  (',' if j+1<len(g['rows']) else ''))
        print('    ]}'+(',' if i+1<len(grids) else ''))
    print('  ]\n}')
