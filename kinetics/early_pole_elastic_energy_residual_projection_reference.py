"""Independently project a screened elastic solve residual onto the discrete energy gradient."""
import argparse
import json

import numpy as np

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_joined_temporal_refinement_reference import select
from early_pole_screened_active_response_update_reference import ActiveResponse
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels
from early_pole_screened_elastic_discrete_gradient_reference import step


def score(model,Y0,Y,h):
    M=model.M
    d=Y-Y0
    E0=M.weights(Y0,True)[3]
    _,F,_,E1,_=M.weights(Y,True)
    _,_,G,*_=M.weights((Y+Y0)/2,True)
    k=M.r[1].start
    def currents(Y):
        j0=sum(M.cells[j].fit(Y[j])[2] for j in range(k))
        v=[sum(M.cells[j].fit(Y[j])[1] for j in range(a.start,a.stop)) for a in M.r[1:]]
        return M.background+np.r_[j0,v]
    def polar(J):
        return np.r_[sum(J[1:]),J[0]+J[1:]]
    v=np.sum(d*d)
    J0,J1=currents(Y0),currents(Y)
    exact=np.sum(M.shape*d)+.5*M.c*(J1+J0)@polar(J1-J0)
    if v>1e-26:
        G=G+(exact-np.sum(G*d))/v*d
    D,_,W,_=model.collision(G,F,G)
    R=d-h*D
    sc=np.repeat([c.W for c in M.cells[k:]],2).reshape(-1,2)
    r_inf=float(np.max(abs(R[k:]/sc)))
    ledger=float(E1-E0)
    r_E=float(np.sum(G*R))
    work=float(h*np.sum(G*D))
    assert abs(ledger-r_E-work)<2e-15
    assert abs(work)<5e-17 and abs(W*h)<5e-17
    assert abs(r_E) <= r_inf*float(np.sum(abs(G[k:]*sc)))+1e-16
    return dict(endpoint_energy_change=ledger,residual_energy_projection=r_E,work_integral=work,
                stable_total_increment=float(exact),residual_infinity_norm=r_inf,
                rescaled_energy_norm_L1=float(np.sum(abs(G[k:]*sc))),
                maximum_charged_moment_change=float(np.max(abs(d[k:]))))


def run(path,projection):
    run=select(path,projection)
    model=Screened(run['grid'],1,len(run['snapshots'][0]['cold_bins']),run['duration'],(1,1,2),.08,
                   1.,run['aperture_grid'],projection)
    soft=ElasticChannels(model,run['grid'],1.)
    active=ActiveResponse(model,.12)
    for a,b in zip(model.M.edge_sets,run['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(a,b)
    h=run['duration']/432
    rows=[]
    for r in run['snapshots']:
        if r['step'] not in {0,run['steps']//2,3*run['steps']//4,run['steps']}:
            continue
        y=np.asarray(r['hard_moments'])
        active.refresh(y,'full')
        ds=[]
        for budget in (None,5e-15):
            v,d=step(soft,y,h,energy_budget=budget)
            info=score(soft,y,v,h)
            assert d['rho_drift']==info['endpoint_energy_change']
            ds.append((v,dict(energy_budget=budget,iterations=len(d['history']),**info)))
        rows.append(dict(snapshot_step=r['step'],snapshot_exposure=r['exposure'],
                         moment_difference_from_new_acceptance=float(np.max(abs(ds[1][0]-ds[0][0]))),
                         default=ds[0][1],tight=ds[1][1]))
    return dict(projection=projection,elastic_width=h,rows=rows,
                scope='One elastic step on fixed snapshots; no new long time evolution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    a=p.parse_args()
    print(json.dumps(run(a.bundle,a.projection),indent=2,allow_nan=False))
