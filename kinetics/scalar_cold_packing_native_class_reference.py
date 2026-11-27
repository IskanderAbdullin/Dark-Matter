"""Build native CLASS controls near the coldest finite-occupation scalar PSD."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import expit

MODES=(.05,.3,1.,2.,3.,5.)
H=.6736


def run(capacity,root):
    d=json.loads(capacity.read_text())
    row=d['fixed_targets'][0]
    assert row['target_omega']==.06
    q1=next(r['q_edge_over_Tgamma0'] for r in row['cap_scan'] if r['cap']==1.)
    root.mkdir(parents=True,exist_ok=True)
    cases=[]
    for tier,eps,n,l,tol in (('compact',.045,96,28,2e-8),('dense',.0225,192,44,1.2e-9)):
        # Smooth the artificial edge to keep the CLASS derivative regular.
        # Correct the q0 scale to preserve the target count at the specified cap.
        u=np.linspace(0.,2.9,1801)
        f=expit((1-u)/eps)
        I=lambda power:quad(lambda v:v**power*expit((1-v)/eps),0.,2.9,
                            epsabs=3e-12)[0]
        i2=I(2)
        scale=(3*i2)**(-1/3)
        moment=scale*np.sqrt(I(4)/i2)/np.sqrt(3/5)
        assert 1.<moment<1.04
        dat=root/f'{tier}-packed.dat'
        np.savetxt(dat,np.c_[u,f],fmt='%.16e')
        for cap in (0.,1.,2e5):
            label=f'{tier}-'+('cdm' if cap==0 else f'F{cap:g}')
            cfg=dict(output='mPk',z_pk='5,20',h=H,T_cmb=2.72548,
                     omega_b=.02237,A_s=2.1e-9,n_s=.965,tau_reio=.054,
                     N_ur=3.044,gauge='synchronous',omega_cdm=.12 if cap==0 else .06,
                     headers='yes',root=str(root.resolve()/(label+'_')),evolver=0,
                     start_small_k_at_tau_c_over_tau_h=2e-7,tol_perturbations_integration=tol,
                     perturbations_verbose=1,background_verbose=2)
            cfg.update({'P_k_max_h/Mpc':5.12,'k_per_decade_for_pk':8,
                        'k_output_values':','.join(f'{k*H:.13g}' for k in MODES)})
            if cap:
                cfg.update(dict(N_ncdm=1,m_ncdm=86.,T_ncdm=q1*cap**(-1/3)*scale,
                                omega_ncdm=.06,use_ncdm_psd_files=1,ncdm_psd_filenames=str(dat.resolve()),
                                ncdm_quadrature_strategy=3,ncdm_maximum_q=2.9,ncdm_N_momentum_bins=n,
                                l_max_ncdm=l,ncdm_fluid_approximation=3))
            (root/(label+'.ini')).write_text(''.join(f'{k} = {v}\n' for k,v in cfg.items()))
            cases.append(dict(label=label,cap=cap,fraction=.5 if cap else 0.,edge_width=eps,
                              momentum_nodes=n,hierarchy=l,
                              rms_ratio_to_sharp_edge=moment,psd_sha256=hashlib.sha256(dat.read_bytes()).hexdigest()))
    return dict(scope="Near-packed smooth shapes considered only as counterfactual initial data. "
                      "Fixed total dark rest density; remaining dark matter is an imposed CDM complement.",
                modes_h_per_Mpc=MODES,cases=cases)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('capacity',type=Path)
    p.add_argument('root',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.capacity,a.root),indent=2,allow_nan=False))
