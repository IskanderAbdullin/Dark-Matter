"""Match a spacelike absorption aperture to the independently monotone emission."""
import gzip,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from finite_mass_qed_absorption_transfer_reference import point
from finite_mass_qed_emission_monotone_bridge_reference import EmissionBridge
from finite_mass_qed_signed_support_bridge_reference import BridgeSource
from finite_mass_qed_integrated_signed_bridge_reference import seed


def setup(d):
    S=EmissionBridge(d)
    unit=2*((S.M+S.m)**2-S.s**2)/(64*np.pi**3*S.q*S.Es)
    fs=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],float(S.q))-1]
    f=1+np.longdouble(1)/fs
    v=seed(S)*[f,1.]
    N0=v[0]-v[1]
    return S,unit,f,N0


def controls(d,setup0):
    S,unit,F,N0=setup0
    A=BridgeSource(d)
    matches=[]
    for k in (.00016,.00019,.0002):
        x=A.signed(-k,64,13,44)[0,:2]
        y=point(d,k,42,58)*k/unit
        rel=np.abs(y/x-1.)
        if max(rel)>1.5e-4:
            raise ArithmeticError(('Independent absorption chart does not match controlled shared band',k,rel))
        matches.append(dict(k=k,absolute_defect=np.abs(x-y).astype(float).tolist(),relative_defect=rel.astype(float).tolist()))
    return matches


def integrals(d,setup0,mesh):
    S,unit,F,N0=setup0
    o,ac,ep,em=mesh
    X,W=roots_legendre(o)
    lo=[.0002,.00032,.00064,.00128,.00256]
    for v in (-1,1):
        h=float(S.s**2/(2*(S.Es+v*S.q)))
        if lo[0]<h<lo[-1]:
            lo.append(h)
    lo=np.sort(lo)
    out=[]
    g2=.005**2/(4*np.pi**2)
    B=N0/S.B
    for a,b in zip(lo[:-1],lo[1:],strict=True):
        l,h=np.log(a),np.log(b)
        kn=np.exp((l+h)/2+(h-l)*X/2)
        rates=[]
        for k in kn:
            Fp,Rp=S.signed(k,em,11,32)[0,:2]*[F,1.]
            Fm,Rm=np.asarray(point(d,k,ac,ep),dtype=np.longdouble)*[F,1.]*k/unit
            Np,Nm=Fp-Rp,Fm-Rm
            Pp,Pm=Fp+Rp,Fm+Rm
            nb=1/np.expm1(np.longdouble(k))
            net=(.5+nb)*(Np+Nm-2*N0)+(Pp-Pm)/2
            real_heat=(.5+nb)*(Np-Nm)+(Pp+Pm)/2
            rates.append([net,real_heat*k])
        Qn,Qh=(h-l)/2*(W@np.array(rates))
        out.append(dict(k_band=[float(a),float(b)],
                        local_leading_subtracted_scalar_ratio=float(g2*Qn/B),
                        real_photon_heat_over_Born_scalar_energy=float(g2*Qh/(B*S.Es))))
    return out


def run(grids):
    root=Path(__file__).parent/'results'
    rows=[]
    for g in grids:
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        x=setup(d)
        chk=controls(d,x)
        a,b=(integrals(d,x,m) for m in ((5,20,40,34),(9,38,58,58)))
        err=max(abs(s[t]-r[t]) for s,r in zip(a,b,strict=True) for t in
                ('local_leading_subtracted_scalar_ratio','real_photon_heat_over_Born_scalar_energy'))
        rows.append(dict(grid=g,overlapping_unweighted_sign_controls=chk,refined=b,
                         max_scaled_band_absolute_change=err,gain_and_heat_totals=[
                             sum(r[t] for r in b) for t in ('local_leading_subtracted_scalar_ratio','real_photon_heat_over_Born_scalar_energy')]))
    return dict(scope='Finite leading soft subtraction and real photon heat in an extended support; finite virtual terms not included.',data=rows)


if __name__=='__main__':
    import sys
    print(json.dumps(run(tuple(sys.argv[1:]) or ('fine','third')),indent=2,allow_nan=False))
