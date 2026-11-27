"""Audit the bare real photon cutoff ladder against the independent Born soft scale."""
import gzip
import json
from pathlib import Path

import numpy as np
from early_pole_decay_mode_clock_reference import Q, intervals


def theoretical(d):
    M=d['x_end']; m=M*d['mass_ratio']
    s=d['checkpoints'][-1]['scalar_pole_over_T']
    Es,a,b,*_=intervals(M,s,m,Q)
    pe=np.hypot(d['endpoint_parent_momentum_edges'],M)
    be=np.hypot(d['endpoint_fermion_momentum_edges'],m)+Es
    ps=np.r_[0.,d['endpoint_parent_occupations'],0.]
    bs=np.r_[0.,d['endpoint_fermion_occupations'],0.]
    edge=np.unique(np.clip(np.r_[a,b,pe,be],a,b))
    y=(edge[1:]+edge[:-1])/2
    p=ps[np.searchsorted(pe,y)]
    f=bs[np.searchsorted(be,y)]
    fs=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1]
    pref=2*((M+m)**2-s*s)/(16*np.pi*Es*Q)
    G,L=pref*np.diff(edge)@np.stack([p*(1-f),f*(1-p)],axis=1)
    chi=np.arccosh((M*M+m*m-s*s)/(2*M*m))
    A=2*(chi/np.tanh(chi)-1)/(4*np.pi**2)
    return A*np.array([G,L]),float(fs),float(G+fs*(G-L))


def run(root,grids):
    result=[]
    for g in grids:
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        s=json.loads((root/f'finite-mass-qed-hard-photo-angular-cut-{g}-2026-11-27.json').read_text())['data'][0]
        scale,fs,B=theoretical(d)
        raw=np.array([[x['fine'][k] for k in ('gain','loss')] for x in s['rows']])
        cut=np.array([x['lab_photon_minimum_energy'] for x in s['rows']])
        frac=raw*cut[:,None]/scale
        compare=np.array([x['net_relative_to_Born_at_g_point005'] for x in s['rows']])
        reconstructed=.005**2*(raw[:,0]+fs*(raw[:,0]-raw[:,1]))/B
        np.testing.assert_allclose(reconstructed,compare,rtol=6e-14,atol=1e-16)
        shape=raw[:-1]-raw[1:]
        if shape.min()<=0.:
            raise ArithmeticError('Finite bins of a positive real gain or loss are not positive')
        max_drift=max(x['refinement_relative_drift'] for x in s['rows'])
        if max_drift>3e-3 or max(abs(frac[0]-1.))>.08:
            raise ArithmeticError((g,'refinement or asymptotic comparison',max_drift,frac[0]))
        result.append(dict(grid=g,coefficient_theoretical=scale.tolist(),cut_times_raw_over_theoretical=frac.tolist(),
                            cut_ladder=cut.tolist(),max_mesh_relative_drift=max_drift,
                            minimum_gain_and_loss_into_finite_bins=shape.min(axis=0).tolist()))
    return dict(scope='Born eikonal normalization and positivity gate only, no inclusive QED rate.',result=result)


if __name__=='__main__':
    import sys
    print(json.dumps(run(Path(__file__).parent/'results',tuple(sys.argv[1:]) or ('fine',)),indent=2,allow_nan=False))
