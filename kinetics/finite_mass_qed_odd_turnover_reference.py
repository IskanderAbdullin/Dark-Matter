"""Separately audit even net and odd turnover in the leading signed subtraction."""
import gzip
import json
from pathlib import Path

import numpy as np
from finite_mass_qed_signed_support_bridge_reference import BridgeSource
from finite_mass_qed_integrated_signed_bridge_reference import seed


def grouped(k,T,N0,Np,Nm,Pp,Pm):
    q=np.longdouble(k)
    return ((Np+Nm-2*N0)/(2*np.tanh(q/(2*T)))+(Pp-Pm)/2)/q


def run():
    root=Path(__file__).parent/'results'
    rows=[]
    for r in json.loads((root/'finite-mass-qed-signed-support-bridge-2026-11-27.json').read_text())['data']:
        g=r['grid']
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        S=BridgeSource(d)
        fs=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],float(S.q))-1]
        factor=1+np.longdouble(1)/fs
        ref=seed(S)*[factor,1.]
        N0=ref[0]-ref[1]
        rr=[]
        for p in r['rows']:
            q=np.longdouble(p['k'])
            Fp,Rp=np.array(p['bare_emission_G_L'],dtype=np.longdouble)*[factor,1.]
            Fm,Rm=np.array(p['bare_absorption_G_L'],dtype=np.longdouble)*[factor,1.]
            Np,Nm=Fp-Rp,Fm-Rm
            Pp,Pm=Fp+Rp,Fm+Rm
            coldG,coldL=np.array(p['vacuum_subtracted_G_L'],dtype=np.longdouble)
            old0=Fp-q*coldG*factor-(Rm-q*coldL)
            y=grouped(q,1.,old0,Np,Nm,Pp,Pm)/N0
            analytic=grouped(q,1.,N0,Np,Nm,Pp,Pm)/N0
            ref0=p['local_net_derivative_over_Born_eik_scale']
            wrong=grouped(q,1.,old0,Np,Nm,Np,Nm)/N0
            if abs(analytic/y-1)>3e-7:raise ArithmeticError(('Baseline shift too large',g,q,analytic,y))
            rr.append(dict(k=float(q),grouped=float(y),previous=float(ref0),analytic_baseline_relative_change=float(analytic/y-1),
                           wrong_odd_net_relative_shift=float((wrong-y)/y)))
            if abs(y/ref0-1)>4e-9:raise ArithmeticError((g,q,y,ref0))
        rows.append(dict(grid=g,controls=rr))
    q=np.array([2e-7,1e-5,4e-5,.001,.03],dtype=np.longdouble)
    T=np.longdouble('.7')
    Pp=1.5+.021*q
    Pm=1.5-.021*q
    h=np.tanh(q/(2*T))
    null=grouped(q,T,0,-h*Pp,h*Pm,Pp,Pm)
    N0=np.longdouble('.15'); n1=np.longdouble('.13');n2=np.longdouble('-.16');p1=np.longdouble('.06')
    d=q[:3]
    smooth=grouped(d,1.,N0,N0+n1*d+n2*d*d/2,N0-n1*d+n2*d*d/2,1.4+p1*d,1.4-p1*d)
    if max(abs(null))>8e-12 or max(abs(smooth-(n2+p1)))>2e-6:
        raise ArithmeticError('Artificial smooth/equilibrium gate failed')
    return dict(scope='Algebra and manufactured smooth/thermal checks; does not impose smoothness on archive edges.',
                thermal_max_defect=float(max(abs(null))), smooth_limit_expected=float(n2+p1),
                smooth_values=smooth.astype(float).tolist(),archive=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
