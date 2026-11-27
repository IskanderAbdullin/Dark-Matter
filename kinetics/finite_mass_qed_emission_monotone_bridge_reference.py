"""Bare one-charge emission in finite bands, independent of a large raw cutoff integral."""
import gzip,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from finite_mass_qed_signed_support_bridge_reference import BridgeSource

base=Path(__file__).parent/'results'


class EmissionBridge(BridgeSource):
    def __init__(self,data):
        super().__init__(data,emission_only_extension=True)


def integrate(src,lo,hi,order,angular):
    x,w=roots_legendre(order)
    a,b=np.log(np.longdouble(lo)),np.log(np.longdouble(hi))
    values=[]
    pref=2*((src.M+src.m)**2-src.s**2)/(64*np.pi**3*src.q*src.Es)
    for u in x:
        k=np.exp((a+b)/2+(b-a)*u/2)
        Born=src.signed(k,angular,10,32)[0,:2]
        n=1/np.expm1(k)
        values.append(pref*Born*np.array([1+n,n]))
    return np.asarray((b-a)/2*(w@np.asarray(values)),dtype=float)


def run(grids):
    from finite_mass_qed_signed_support_bridge_reference import BridgeSource
    agg=[]
    for g in grids:
        with gzip.open(base/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        S=EmissionBridge(d)
        check=BridgeSource(d)
        defect=float(max(abs(S.signed(2e-4,38,10,28).ravel()-check.signed(2e-4,38,10,28).ravel())))
        if defect>1e-12:raise ArithmeticError(('shared endpoint defect',defect))
        cuts=[4e-5,8e-5,1.6e-4,3.2e-4,6.4e-4,.00128,.00256]
        windows=[]
        for a,b in zip(cuts[:-1],cuts[1:],strict=True):
            eq0=integrate(S,a,b,5,32)
            eq1=integrate(S,a,b,10,56)
            drift=max(abs(eq0/eq1-1.))
            windows.append(dict(edges=[a,b],gain_loss_one_charge=eq1.tolist(),max_two_rule_relative_drift=float(drift)))
        agg.append(dict(grid=g,bridge_endpoint_abs_defect=defect,windows=windows))
    return dict(scope='One-charge bare emission in finite bands. Upper threshold is the certified limit of the monotone parent chart, not a physical separation.',data=agg)


if __name__=='__main__':
    import sys
    print(json.dumps(run(tuple(sys.argv[1:]) or ('fine','third')),indent=2,allow_nan=False))
