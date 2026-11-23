"""Non-overlap DG P1 has an intrinsic *within-cell* Bose η tangent obstruction.

Even with both local scalar source moments negative, the inverse physical
positive 2x2 susceptibility can drive eta_rest toward 0: C=-[u,v],
J=-[[a,b],[b,c]], eta_dot_rest=(c*u-b*v)/(ac-b²). No false claim
of physical BEC or impossibility of every finite-step solution.
"""
import math
import numpy as np
from kinetics.phi_bose_cell_entropy import fixture,source,mode_moments,implicit_step
from kinetics.phi_bose_daughter_log_angle import physical_events


def local_spd_counterexample(a=3.,b=1.,c=2.,u=.05,v=1.):
    det=a*c-b*b
    if not(a>0 and b>=0 and c>0 and det>0 and u>=0 and v>=0):
        raise ValueError('strict SPD Gram, nonpositive moment sources')
    G=np.array([[a,b],[b,c]])
    C=-np.array([u,v])
    rate=np.linalg.solve(-G,C)
    explicit=(c*u-b*v)/det
    return {'G':G.tolist(),'C':C.tolist(),
            'eta_rate':rate.tolist(),'analytical_first_eta_rate':explicit,
            'rest_eta_decreases_despite_both_local_sources_nonpositive':bool(rate[0]<0)}


def physical_case(n=32,nangle=48,nrad=24,scheme='log_ephi',*,run_time_tests=False):
    modes,eta,_=fixture(n=n,pmax=2.)
    events,_=physical_events(pmax=2.,nr=nrad,nangle=nangle,scheme=scheme)
    C,_,H=source(modes,eta,events)
    k=modes[2].size
    cp=C[-k:]
    J=mode_moments(modes[2],eta[2],order=24)['J']
    rate=np.linalg.solve(J,cp)
    G=-J[:2,:2]
    a,b,c=G[0,0],G[0,1],G[1,1]
    u,v=-cp[0],-cp[1]
    analytic=(c*u-b*v)/(a*c-b*b)
    out={'N':n,'nangle':nangle,'nrad':nrad,'angular_scheme':scheme,
         'first_cell_C':cp[:2].tolist(),'first_cell_G':G.tolist(),
         'first_rest_eta':float(eta[2][0]),'first_eta_derivative':float(rate[0]),
         'second_eta_derivative':float(rate[1]),
         'analytic_eta_derivative':float(analytic),
         'boundary_linear_dt':float(eta[2][0]/(-rate[0])) if rate[0]<0 else math.inf,
         'source_Rphi':float(sum(cp)),'H':H,
         'all_scalar_source_negative':bool(np.all(cp<=0)),
         'not_physical_BEC':True}
    if run_time_tests:
        trials={}
        for dt in (20.,30.,1000.):
            try:
                q=implicit_step(modes,eta,events,dt)
                trials[str(dt)]={'accepted':True,'min_phi_eta':float(min(q['eta'][2])),
                    'Nphi_change':q['delta']['Nphi'],'S_gain':q['delta']['S'],
                    'max_charge_energy_error':max(abs(q['delta'][key]) for key in ('N12','N1phi','U'))}
            except ArithmeticError as exc:
                trials[str(dt)]={'accepted':False,'failure':str(exc)}
        out['implicit_trials']=trials
    return out


def audit():
    return {'schema_version':1,'analytic_two_local_cell':local_spd_counterexample(),
            'physical_cases':[physical_case(n,angle,run_time_tests=angle in (12,48))
                              for n in (16,32) for angle in (12,24,48,96)],
            'uniform_original_n32_control':physical_case(32,12,scheme='uniform_cos',run_time_tests=True),
            'scope':'synthetic on-shell Yukawa quantum Bose IR tangent: numerical feasible domain, not physical BEC or 86eV relic',
            'class_ready':False}


if __name__=='__main__':
    import json
    print(json.dumps(audit(),indent=2))
