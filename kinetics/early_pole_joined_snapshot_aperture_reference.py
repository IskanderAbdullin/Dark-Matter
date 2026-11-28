"""Replay the dressed 1↔2 kernel and decompose archived live snapshot sources."""
import argparse
import json
import numpy as np
from scipy.special import expit

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_decay_mode_clock_reference import G, KAPPA
from early_pole_live_pauli_gauge_moment_decay_reference import Tool
from early_pole_joined_history_memory_reference import load


def replay(a):
    grid=a['grid']
    C=Cold(grid,1,0,len(a['snapshots'][0]['cold_bins']),True,exposure=a['duration'],
           reference_aperture_grid=a['aperture_grid'],archive_projection=a['projection'],
           aperture_subdivision=a.get('aperture_subdivision',1))
    for model,replay in zip(C.M.edge_sets,a['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(model,replay)
    rows=[]
    for s in a['snapshots']:
        Y=np.asarray(s['hard_moments'])
        e,f,*_=C.M.weights(Y,True)
        np.testing.assert_allclose(f,s['hard_logits'],rtol=2e-13,atol=2e-13)
        N=np.asarray(s['cold_bins'])
        E,_,_,_,m=C.scalar(N)
        factor=G*KAPPA*.05*((1.06*.05)**2-m*m)/(C.q*E)
        I=np.zeros((len(C.q),3))
        net_stim=np.zeros_like(N)
        for j,es in enumerate(E):
            for i,k,u,v,w in Tool.events(C.lambdas[j],dict(coeff=e[C.M.start:],es=es)):
                i+=C.M.start
                k+=C.M.start
                P=expit(f[i,0]+f[i,1]*u)
                B=expit(f[k,0]+f[k,1]*v)
                I[j] += [w@P,w@B,w@(P*(1-B))]
                net_stim[j] += w@(P-B)
        split=np.stack([N*factor*I[:,0],N*factor*I[:,1],C.sw*factor*I[:,2]],axis=1)
        net=N*factor*net_stim+split[:,2]
        saved=np.asarray(s['instantaneous_source_by_q'])
        defect=np.max(np.abs(net-saved))
        assert defect<5e-13
        rows.append(dict(step=s['step'],exposure=s['exposure'],scalar_pole=m,
                         parent_daughter_seed_integrals=I.tolist(),
                         parent_stimulus=float(split[:,0].sum()),
                         daughter_stimulus=float(split[:,1].sum()),
                         spontaneous=float(split[:,2].sum()),
                         scalar_rate=float(net.sum()),
                         maximum_absolute_source_per_mode_replay_defect=float(defect)))
    return dict(grid=grid,projection=a['projection'],q=C.q.tolist(),rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs',nargs='+')
    a=p.parse_args()
    print(json.dumps([replay(load(x)) for x in a.inputs],indent=2,allow_nan=False))
