"""Preserve positive exterior cohorts across restart epochs."""
import json

import numpy as np

from frw_exterior_cohort_reference import Exterior

SCHEMA = "positive-exterior-restart-v1"
FIELDS = ("a_injection_over_restart", "mass_GeV", "source_channel_id",
          "high_tail_flag", "injection_boundary_energy_GeV",
          "number_in_restart_a3_units", "mean_injection_energy_GeV")


def export_restart(exterior):
    """Serialize in a local a=1 frame at the current endpoint."""
    a=float(exterior.a)
    assert np.isfinite(a) and a>=1.
    exterior.result()
    return dict(schema=SCHEMA,fields=list(FIELDS),
                cohorts=[[float(ai/a),float(m),int(j),int(side),float(bound),float(n/a**3),float(e)]
                         for ai,m,j,side,bound,n,e in exterior.cohorts],
                injected_in_restart_a4_units=float(exterior.injected/a**4),
                trace_in_restart_a4_units=float(exterior.trace/a**4))


def import_restart(state):
    """Use only cohort data saved by the positive projection reference."""
    if state["schema"] != SCHEMA or state["fields"] != list(FIELDS):
        raise ValueError("Unknown exterior cohort format.")
    res=Exterior()
    x=np.asarray(state["cohorts"],dtype=float).reshape(-1,7)
    ai,m,j,side,bd,n,e=x.T
    assert np.all(np.isfinite(x)) and np.all((0.<ai)&(ai<=1.))
    assert np.all((m>0.)&(bd>=m)&(n>0.)&(e>=m))
    assert np.all((0<=j)&(j<=4)&(j==np.floor(j))) and np.all((side==0)|(side==1))
    assert np.all(np.where(side==0,e<=bd*(1+2e-12),e>=bd*(1-2e-12)))
    res.cohorts=[(float(a),float(b),int(c),int(d),float(f),float(g),float(h)) for a,b,c,d,f,g,h in x]
    res.trace=float(state["trace_in_restart_a4_units"])
    res.injected=float(state["injected_in_restart_a4_units"])
    assert np.isfinite(res.trace+res.injected) and res.trace>=0. and res.injected>=0.
    assert abs(np.sum(ai*n*e)-res.injected)<3e-12*max(res.injected,1e-280)
    res.result()
    return res


def test():
    f=Exterior()
    for a in (1.,1.2,1.3,1.9):
        f.advance(a)
        for i,m in enumerate((8.6e-8,.02)):
            for s,p in enumerate((.01,.04)):
                n=1/(i+2)*(a+s)
                K=p*p/(np.hypot(p,m)+m)
                f.add(a,m,i,s,.02,n,a*n*K)
    A=f.a
    g=import_restart(json.loads(json.dumps(export_restart(f),allow_nan=False)))
    energy0,trace0=g.energy(1.),g.trace
    for b in (1.01,1.5,3.,4.):
        f.advance(A*b)
        g.advance(b)
        F=np.array([f.energy(A*b),f.trace,f.upper_bound(A*b)])/A**4
        G=np.array([g.energy(b),g.trace,g.upper_bound(b)])
        assert np.max(abs(G/F-1.))<2e-14
    assert abs(g.energy(g.a)-energy0-(g.trace-trace0))<2e-12
    print("direct and portable continuation give identical exterior gravity and trace.")
    zero=import_restart(export_restart(Exterior()))
    assert zero.result()["exterior_cohorts"]==0


if __name__=="__main__":
    test()
