"""At the bare poles the scalar Yukawa residue resolves into sequential spin sums."""
import gzip,json
from pathlib import Path
import numpy as np

from early_pole_hard_yukawa_photoproduction_reference import I,G,slash,spinors
from finite_mass_qed_nonnull_photon_ward_poles_reference import setup


def study(data,row,ph,root):
    key=(row['k'],row['omega'])
    c={(.0008,.0002):0.,(.00128,.00064):-.5,(.00256,.00064):.5}[key]
    M,m,q,K,transfer,geo=setup(data,*key,c)
    P,b=geo(root['parent_E'],ph)
    out=spinors(b,m).conj().T@G[0]
    incoming=spinors(P,M)
    which=root['arm']
    L,mass=(P+K,M) if which==0 else (b-K,m)
    ui=spinors(L,mass)
    ui_bar=ui.conj().T@G[0]
    scal=(out@ui) if which==0 else (ui_bar@incoming)
    gauge=(lambda p:ui_bar@slash(p)@incoming) if which==0 else (lambda p:out@slash(p)@ui)
    radi=[gauge(p) for p in np.eye(4)[1:3]]
    D=scal.conj().T@scal if which==0 else scal@scal.conj().T
    unpol=np.trace(D).real/2
    rel_D=np.max(abs(D-unpol*np.eye(2)))/unpol
    total=ref=0.
    for eps,V in zip(np.eye(4)[1:3],radi,strict=True):
        X=scal@V if which==0 else V@scal
        mat=slash(L)+mass*I
        direct=out@(mat@slash(eps) if which==0 else slash(eps)@mat)@incoming
        ref+=np.sum(abs(direct)**2)
        total+=np.sum(abs(X)**2)
    traced=unpol*sum(np.sum(abs(x)**2) for x in radi)
    ward=np.sum(abs(gauge(K))**2)/sum(np.sum(abs(x)**2) for x in radi)
    errors=[abs(total/ref-1.),abs(traced/ref-1.),rel_D,ward]
    if max(errors)>3e-8:
        raise ArithmeticError((row['k'],which,errors))
    return dict(k=row['k'],omega=row['omega'],angle_azimuth=ph,pole_arm=which,
                parent_energy=root['parent_E'],scalar_subprocess_spin_anisotropy=float(rel_D),
                factored_amplitude_weight_relative_defect=float(errors[0]),
                spin_summed_product_relative_defect=float(errors[1]),
                normalized_on_shell_current_ward_squared=float(ward))


def run():
    p=Path(__file__).parent/'results'
    objects=json.loads((p/'finite-mass-qed-nonnull-photon-ward-poles-2026-11-27.json').read_text())['data']
    result=[]
    for g in objects:
        with gzip.open(p/f"early-pole-spectrum-{g['grid']}-m06-2026-11-26.json.gz") as f:d=json.load(f)
        z=[study(d,row,az['azimuth'],pole) for row in g['checks'] for az in row['azimuths']
           for pole in az['found_near_born_intermediate_on_shell_roots']]
        if len(z)!=24:raise ArithmeticError('Gate lost one of the occupied roots')
        result.append(dict(grid=g['grid'],rows=z))
    return dict(scope='Static on-shell numerator and spin density only. No width, sequential collision probabilities, statistical weighting or in-medium rate.',data=result)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
