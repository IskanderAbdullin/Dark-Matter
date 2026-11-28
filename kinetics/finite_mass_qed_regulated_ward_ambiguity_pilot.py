"""Test whether imposing K.J=0 on a width-regulated diagram fixes its current."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_hard_yukawa_photoproduction_reference import G,I,dot,slash,spinors
from finite_mass_qed_nonnull_photon_ward_poles_reference import setup


def shot(archive,row,ds):
    k,w=row['physical_k'],row['physical_omega']
    h=row['separation_in_widths']
    eps=row['width']
    M,m,q,K,T,geo=setup(archive,k,w,.25)
    ph=row['azimuth']
    ref=row['edge_parent_energy']
    K2=dot(K,K)
    e=brentq(lambda e:((K2+2*dot(geo(e,ph)[0],K))/eps)-ds,
             ref-.003,ref+.003,xtol=3e-16)
    p,b=geo(e,ph)
    raw=(K2+2*dot(p,K))/eps
    # Use the exact geometric invariant sum, not the separately rounded
    # second denominator divided by the tiny diagnostic epsilon.
    d1,d2=raw,2*dot(q,K)/eps-raw
    if abs(h-d1-d2)>3e-11 or abs(raw-ds)>3e-6:
        raise ArithmeticError(('Archived pole chart was not reproduced',h,d1,d2,ds))
    spin=spinors(p,M)
    bar=spinors(b,m).conj().T@G[0]
    A=slash(p)+slash(K)+M*I
    B=slash(b)-slash(K)+m*I
    S=bar@spin
    # J below represents epsilon times the ordinary finite-width current.
    # Derive epsilon K.J analytically to avoid cancellation in the gauge test.
    mask=np.array([1.,-1.,-1.,-1.])
    a=np.array([bar@A@slash(v)@spin for v in np.eye(4)])*mask[:,None,None]
    bcur=np.array([bar@slash(v)@B@spin for v in np.eye(4)])*mask[:,None,None]
    left=np.array([0.,1.,0.,0.])
    ell=np.array([k,0.,0.,w])/np.sqrt(-K2)
    ref_vectors=np.array([K,[1.,0.,0.,0.],[1.,.7,0.,0.]])
    basenorm=np.sum(abs(a[1])**2)+np.sum(abs(a[2])**2)
    result=[]
    for sign in (1,-1):
        j=a/(d1+1j)+bcur/(d2+sign*1j)
        defect=eps*S*(d1/(d1+1j)-d2/(d2+sign*1j))
        independent=np.einsum('a,aij->ij',mask*K,j)
        current=np.sum(abs(j)**2)
        # This checks the width-induced Ward identity, including cases in
        # which the analytic defect vanishes, against explicit spinors.
        defect_delta=np.linalg.norm(independent-defect)/max(eps*np.linalg.norm(S),1e-30)
        if defect_delta>2e-6:
            raise ArithmeticError(('Contracted spinor identity',row,ds,sign,defect_delta))
        jd=j[None,:]-np.array([r/dot(r,K) for r in ref_vectors])[:,:,None,None]*defect[None,None,:]
        Ward=np.einsum('a,baij->bij',mask*K,jd)
        ward_remains=np.max(abs(Ward))/max(np.sqrt(current),1e-30)
        T=np.einsum('a,baij->bij',mask*left,jd)
        L=np.einsum('a,baij->bij',mask*ell,jd)
        rate=np.array([np.sum(abs(z)**2,axis=(1,2)) for z in (T,L)])/basenorm
        if ward_remains>2e-13:
            raise ArithmeticError(('Projected Ward identity',ward_remains))
        if rate[0,0]!=rate[0,1] or rate[1,1]!=rate[1,2]:
            raise ArithmeticError('A reference shifted an exactly orthogonal polarization')
        if h==0. and sign==-1 and np.linalg.norm(defect)>3e-12*eps*np.linalg.norm(S):
            raise ArithmeticError('Equal-width opposite orientation failed at exact coalescence')
        result.append(dict(imaginary_second_sign=sign,analytic_defect_l2=float(np.linalg.norm(defect)),
                           spinor_contraction_discrepancy=float(defect_delta),completed_ward_max=ward_remains,
                           transverse_and_longitudinal_proxies_by_reference=rate.tolist()))
    return dict(grid=row['grid'],separation=h,scaled_parent_denominator=ds,physical_width=eps,
                k=k,omega=w,theoretical_references=['parallel K','bath time','tilted timelike'],data=result)


def run():
    p=Path(__file__).parent/'results'
    raw=json.loads((p/'finite-mass-qed-coalescing-histogram-edge-2026-11-28.json').read_text())['data']
    out=[]
    for grid in ('fine','third'):
        with gzip.open(p/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:
            d=json.load(f)
        for r in raw:
            if r['grid']==grid and r['edge_displacement_in_widths']==0. and r['separation_in_widths'] in (0.,2.,16.):
                for x in (-2.,-.5,.5,2.5):
                    out.append(shot(d,r,x))
    diff=[np.ptp(y['transverse_and_longitudinal_proxies_by_reference'],axis=1)
          for r in out for y in r['data']]
    peak=np.max(diff,axis=0)
    if len(out)!=48 or np.min(peak)<1e-5:
        raise ArithmeticError(('Failed to verify transverse-current ambiguity',len(out),peak))
    return dict(scope='Divergence repair of a diagnostic bare amplitude, not a dressed current or rate.',
                maximum_absolute_reference_effect_T_and_L=peak.tolist(),data=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
