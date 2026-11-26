"""Direct light channel: FD blocking and inverse reaction on one geometry.

The visible bath is thermal (zero lepton chemical potential). Both final
dark legs are drawn from the same isotropic stable fermion distribution.
This check does not evolve the dark bath or model the heavy production channel.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from empty_lepton_pair_spectrum_reference import Source
from sm_dark_entropy_exchange_reference import StandardModel
from three_species_decay_onsager_reference import m


def evaluate(T, theta, order=(8,22,16), angle=32):
    src=Source(order,dark_mass=m)
    q,p,v,N,A=src.event_moments(T)
    x,w=roots_legendre(angle)
    E1=(q[:,:,None]+p[:,:,None]*np.sqrt(v)[:,None,None]*x)/2
    E2=q[:,:,None]-E1
    G=(N[:,:,None]+A[:,:,None]*(3*x*x-1)/2)*w
    assert np.all(G>=0.) and np.min(E1)>m and np.min(E2)>m
    a,b=theta(E1),theta(E2)
    logB=-np.logaddexp(0.,-a)-np.logaddexp(0.,-b)
    aff=q[:,:,None]/T-a-b
    block=np.exp(logB)
    rev=np.exp(logB+aff)
    # B-exp(logB+aff) with midpoint scaling for very different occupations
    # and exact expm1 cancellation close to detailed balance.
    signed=-np.sign(aff)*np.exp(logB+np.maximum(aff,0.))*(-np.expm1(-np.abs(aff)))
    normalize=np.array([G.sum(),np.sum(G*E1)])
    def moments(z):
        c=G*z
        return np.array([c.sum(),np.sum(c*E1)])/normalize
    result=dict(net=moments(signed).tolist(),forward_pauli=moments(block).tolist(),inverse=moments(rev).tolist())
    assert np.max(abs(np.array(result['net'])-np.array(result['forward_pauli'])+np.array(result['inverse'])))<6e-14
    result["net_partner_exchange_L1"]=float(np.sum(G*abs(signed-signed[:,:,::-1]))/G.sum())
    assert result["net_partner_exchange_L1"]<1e-12
    return result


def test(eos,row):
    T=row['TF']
    edge=T*np.array(row['scalar_differential_bins'])
    n=eos.state(T)[2]*sum(row['stable_yields'])
    vol=4*np.diff(edge**3)/(6*np.pi**2)
    f=n*np.array(row['stable_normalized_counts'])/vol
    assert np.min(f)>=0. and np.max(f)<1.
    def raw(E):
        p=np.sqrt((E-m)*(E+m))
        j=np.searchsorted(edge,p,side="right")-1
        ff=f[np.clip(j,0,len(f)-1)]
        inside=(0<=j)&(j<len(f))&(ff>0.)
        return np.where(inside,np.log1p(-ff)-np.log(np.where(ff>0.,ff,1.)),1000.)
    def matched(E):
        return E/(row["FD_moment_matched_T_over_SM"]*T)-row["FD_moment_matched_eta"]
    models=dict(vacuum=lambda E:np.full_like(E,1000.),
                detailed_balance=lambda E:E/T,
                overpopulated_equal_T=lambda E:E/T-.2,
                raw_empty_production_snapshot=raw,
                FD_moment_matched_snapshot=matched)
    controls=[]
    for order,angle in (((8,22,16),32),((12,32,20),48)):
        r={name:evaluate(T,fn,order,angle) for name,fn in models.items()}
        assert max(abs(np.array(r["detailed_balance"]["net"])))<1e-14
        assert r["vacuum"]["net"]==[1.,1.]
        assert max(r["overpopulated_equal_T"]["net"])<0.
        controls.append(dict(order=order,angle=angle,results=r))
    return dict(T_GeV=T,maximum_raw_bin_occupation=float(f.max()),controls=controls,
                moments="number and total energy, as fractions of their separate vacuum production rates")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("spectra")
    a=p.parse_args()
    rows=json.loads(Path(a.spectra).read_text())["rows"]
    row=next(r for r in rows if r['TR']==.05 and abs(r['TF']-.045)<1e-8)
    print(json.dumps(test(StandardModel(a.eos),row),indent=2,allow_nan=False))
