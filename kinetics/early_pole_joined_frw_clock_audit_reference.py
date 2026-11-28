"""Compare matched expansion and Hubble clocks without swapping collision controls."""
import argparse
import json

import numpy as np
from scipy.interpolate import PchipInterpolator


def compare(fixed, reference, full, optical=None, shorter=None):
    n=full['steps']
    T=full['duration']
    if n<2 or not T>0.:
        raise ValueError('At least two macro steps at the same positive exposure are required.')
    flags=((False,False),(True,False),(True,True))
    arms=(fixed,reference,full)
    for a,f in zip(arms,flags,strict=True):
        if (a['expands'],a['full_H'])!=f:
            raise ValueError('Supply fixed, expanding radiation, and expanding conditional full H, in that order.')
        for k in ('grid','steps','duration','projection','initial_dark_visible'):
            if a[k]!=full[k]:
                raise ValueError(f'Mismatched {k}.')
    c,r,d=(np.array(a['state_history'],dtype=float) for a in arms)
    if not (len(c)==len(r)==len(d)==n+1 and c[0]==r[0]==d[0]):
        raise ValueError('Histories must be complete and share a common initial count.')
    if not np.all(np.isfinite([c,r,d])):
        raise ValueError('Nonfinite clock history.')
    t=np.linspace(0,T,n+1)
    c,r,d=(v-v[0] for v in (c,r,d))
    z=t/np.sqrt(1+full['initial_dark_visible'])
    # These are postprocessed guide curves, not an independently evolved
    # shorter expansion. Report the interpolation sensitivity explicitly.
    clock=PchipInterpolator(t,r)(z)
    linear=np.interp(z,t,r)
    extras={}
    for name,extra in [('same_eta_constant_clock',optical),('shorter_expansion_reference',shorter)]:
        if extra is None:
            continue
        if not extra['expands'] or extra['full_H']:
            raise ValueError(f"Unexpected clock flag in {name}.")
        if any(extra[k]!=full[k] for k in ('grid','steps','projection','initial_dark_visible')):
            raise ValueError(f"Unexpected initial state in {name}.")
        hist=np.array(extra['state_history'])
        if len(hist)!=n+1 or hist[0]!=full['state_history'][0]:
            raise ValueError(f"Unexpected history in {name}.")
        expected_dt=T if name=='same_eta_constant_clock' else z[-1]
        expected_scale=z[-1]/T if name=='same_eta_constant_clock' else 1.
        if not np.isclose(extra['duration'],expected_dt,rtol=1e-14,atol=0.) or not np.isclose(
            extra.get('constant_clock_scale',1.),expected_scale,rtol=1e-14,atol=0.):
            raise ValueError(f"Unexpected duration or rate scale in {name}.")
        g=hist-hist[0]
        extras[name]=dict(final_gain=float(g[-1]),final_minus_full=float(g[-1]-d[-1]),
                          max_index_aligned_deviation=float(np.max(abs(g-d))),
                          norm_note='Identical scale path' if name=='same_eta_constant_clock' else
                                    'Different scale path; matching fractions of their own total durations.')
    vis=np.pi**2*106.75/30
    initial_energy=full['initial_dark_visible']*vis
    final_rel=(initial_energy+full['trace_integral']+full['energy_ledger'])/vis
    accepted_bisections=[sum('elastic_substeps' in i for i in a.get('implicit',[])) for a in arms]
    return dict(grid=full['grid'],steps=n,duration=T,projection=full['projection'],
                initial_dark_visible=full['initial_dark_visible'],
                reference_clock_rescale=float(z[-1]/T),
                final_gains=dict(fixed=float(c[-1]),expanding_radiation=float(r[-1]),expanding_full=float(d[-1])),
                expansion_relative_to_fixed=float((r[-1]-c[-1])/abs(c[-1])) if c[-1] else None,
                full_clock_relative_to_radiation=float((d[-1]-r[-1])/abs(r[-1])) if r[-1] else None,
                final_full_minus_interpolated_initial_clock=float(d[-1]-clock[-1]),
                final_clock_interpolation_spread=float(abs(clock[-1]-linear[-1])),
                maximum_clock_interpolation_spread=float(np.max(abs(clock-linear))),
                final_interval_rates=[float((x[-1]-x[-2])/(T/n)) for x in (c,r,d)],
                energy_ledgers=[a['energy_ledger'] for a in arms],
                entropy_ledgers=[a['entropy_ledger'] for a in arms],
                trace_integrals=[a['trace_integral'] for a in arms],
                macro_steps_with_elastic_bisection=accepted_bisections,
                final_full_clock_rate_scale_from_energy=float(1/np.sqrt(1+final_rel)),
                independent_controls=extras,
                reference_clock_scope='Interpolation diagnostic only; not an independent shorter reference trajectory.')


def fixture():
    a=dict(grid='fine',steps=4,duration=.04,projection='steps',initial_dark_visible=3.,
           energy_ledger=0.,entropy_ledger=0.,trace_integral=0.)
    s=np.linspace(0.,1.,5)
    f=dict(a,expands=False,full_H=False,state_history=.01+s)
    r=dict(a,expands=True,full_H=False,state_history=.01+2*s)
    d=dict(a,expands=True,full_H=True,state_history=.01+s)
    answer=compare(f,r,d)
    assert answer['reference_clock_rescale']==.5
    assert abs(answer['final_full_minus_interpolated_initial_clock'])<1e-14
    assert abs(answer['full_clock_relative_to_radiation']+.5)<1e-15


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('fixed')
    p.add_argument('radiation')
    p.add_argument('full')
    p.add_argument('--optical')
    p.add_argument('--shorter')
    a=p.parse_args()
    fixture()
    def load(path):
        if path is None:
            return None
        with open(path) as f:
            return json.load(f)
    print(json.dumps(compare(load(a.fixed),load(a.radiation),load(a.full),load(a.optical),load(a.shorter)),indent=2,allow_nan=False))
