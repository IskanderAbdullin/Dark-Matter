"""Gates for a finite massive absorption chart across the neutral null transfer."""
import gzip,json
from pathlib import Path
import numpy as np
from finite_mass_qed_absorption_transfer_reference import transfer,angle_pilot,Q
from finite_mass_qed_signed_support_bridge_reference import BridgeSource


def run():
    p=Path(__file__).parent/'results'
    samples=json.loads((p/'finite-mass-qed-absorption-transfer-sampled-2026-11-27.json').read_text())['data']
    data=[]
    for g in ('fine','third'):
        with gzip.open(p/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        S=BridgeSource(d)
        M,m,s=map(float,(S.M,S.m,S.s))
        Es=float(S.Es)
        null=[]
        for c in (-.97,-.5,.2,.7,.999):
            k0=s*s/(2*(Es-Q*c))
            Es0,J,U,e,V,D,lo,hi,co,si=transfer(M,m,s,k0,c)
            analytic=(D*D+4*M*M*e*e)/(4*e*D)
            l1=transfer(M,m,s,k0*(1-1e-9),c)
            l2=transfer(M,m,s,k0*(1+1e-9),c)
            if l1[2]<=0 or l2[2]>=0 or not np.isfinite(l1[7]) or not np.isinf(l2[7]):
                raise ArithmeticError('Crossing did not change only the upper branch')
            drift=max(abs(lo/analytic-1),abs(l1[6]/analytic-1),abs(l2[6]/analytic-1))
            if drift>7e-9:raise ArithmeticError('Lower endpoint did not continue at null transfer')
            null.append(dict(cos=c,zero_transfer_k=k0,lower_limit_relative_defect=drift))
        pref=2*((S.M+S.m)**2-S.s**2)/(64*np.pi**3*S.q*S.Es)
        comp=[]
        for q in samples[g]:
            if q['k']>.0002:continue
            k=q['k']
            A=np.array(q['gain_loss_before_Planck'])
            B=np.asarray(pref*S.signed(-k,64,13,44)[0,:2]/k,dtype=float)
            err=np.abs(A/B-1.)
            if max(err)>1.5e-4:
                raise ArithmeticError((g,k,A,B,err))
            comp.append(dict(k=k,independent_parent_signed_chart_relative_errors=err.tolist()))
        if any(min(z['gain_loss_before_Planck'])<=0 or z['two_mesh_relative_drift']>3e-5 for z in samples[g]):
            raise ArithmeticError('Sampled absorption quadrature gate')
        data.append(dict(grid=g,angle_and_direct_trace=angle_pilot(M,m,s),null_continuation=null,shared_chart=comp))
    return dict(scope='Positive bare real absorption; no finite virtual interpretation.',data=data)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
