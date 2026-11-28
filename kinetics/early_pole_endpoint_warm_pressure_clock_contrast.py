"""Compare coarse/fine pressure contrasts at one fixed final measurement mass."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_spectrum_audit import read


def values(card,kind):
    if (card['x_start'],card['x_end'],card['steps'])!=(.04,.05,[1024,2048]):
        raise ArithmeticError('Expected paired last-segment clocks')
    if kind in ('Gauss','Gauss3','Gauss4','Gauss5','Gauss6'):
        if card.get('warm_gauss_order',2)!={'Gauss':2,'Gauss3':3,'Gauss4':4,'Gauss5':5,'Gauss6':6}[kind]:
            raise ArithmeticError('Invalid matched warm Gaussian chart')
        x=np.asarray([card['coarse_count_and_pressures'],card['fine_count_and_pressures']])
        declared=np.asarray(card['coarse_over_fine_minus_one'])
        trace=card['maximum_trace_absolute_difference']
    else:
        if card['warm_div']!=kind:
            raise ArithmeticError('Wrong geometric warm chart')
        x=np.asarray(card['coarse_fine_count_pressure_and_warm_pressure'])
        declared=np.asarray(card['coarse_fine_relative_change'])
        trace=card['max_archived_trace_difference']
    if (x.shape!=(2,3) or declared.shape!=(3,) or min(x.ravel())<=0 or not np.all(np.isfinite(x))
            or not np.all(np.isfinite(declared)) or not 0<=trace<=5e-13):
        raise ArithmeticError('Invalid paired moments or trace replay')
    if np.max(abs(x[0]/x[1]-1-declared))>1e-15:
        raise ArithmeticError('Declared coarse-to-fine shift disagrees with saved rates')
    return x


def run(w1,w2,g,w4=None,gauss3=None,gauss4=None,gauss5=None,gauss6=None):
    mass=w1['fixed_pressure_mass']
    if not np.isfinite(mass) or mass<=0 or any(card['fixed_pressure_mass']!=mass for card in (w2,g) if card is not None):
        raise ArithmeticError('Not all moments share one fixed pressure mass')
    a,b,c=[values(r,k) for r,k in ((w1,1),(w2,2),(g,'Gauss'))]
    pairs=dict(W2_over_W1=b/a-1,Gauss_over_W1=c/a-1,Gauss_over_W2=c/b-1)
    reference=(a,b,c)
    if w4 is not None:
        if w4['fixed_pressure_mass']!=mass:
            raise ArithmeticError('W4 measurement mass changed')
        d=values(w4,4)
        pairs.update(W4_over_W1=d/a-1,W4_over_W2=d/b-1,Gauss_over_W4=c/d-1)
        reference=(*reference,d)
    scope='Three same-restart last-segment sensitivities, with each chart using its own archived .04 state. No full-history timestep error bound.'
    if gauss3 is not None:
        if w4 is None or gauss3['fixed_pressure_mass']!=mass:
            raise ArithmeticError('Need all four reference charts and the same mass for Gauss3')
        h=values(gauss3,'Gauss3')
        pairs.update(Gauss3_over_W1=h/a-1,Gauss3_over_W2=h/b-1,Gauss3_over_Gauss2=h/c-1,Gauss3_over_W4=h/d-1)
        reference=(*reference,h)
        scope='Last-segment time-step contrast on five separately archived warm chart restarts; no full-history error bound.'
    if gauss4 is not None:
        if gauss3 is None or gauss4['fixed_pressure_mass']!=mass:
            raise ArithmeticError('Need all five reference charts and the same mass for Gauss4')
        j=values(gauss4,'Gauss4')
        pairs.update(Gauss4_over_W1=j/a-1,Gauss4_over_W2=j/b-1,Gauss4_over_Gauss2=j/c-1,
                     Gauss4_over_W4=j/d-1,Gauss4_over_Gauss3=j/h-1)
        reference=(*reference,j)
        scope='Last-segment time-step contrast on six separately archived warm chart restarts; no full-history error bound.'
    if gauss5 is not None:
        if gauss4 is None or gauss5['fixed_pressure_mass']!=mass:
            raise ArithmeticError('Need all six reference charts and their common mass for Gauss5')
        k=values(gauss5,'Gauss5')
        pairs.update(Gauss5_over_W1=k/a-1,Gauss5_over_W2=k/b-1,Gauss5_over_Gauss2=k/c-1,
                     Gauss5_over_W4=k/d-1,Gauss5_over_Gauss3=k/h-1,Gauss5_over_Gauss4=k/j-1)
        reference=(*reference,k)
        scope='Last-segment time-step contrast on seven separately archived warm chart restarts; no full-history error bound.'
    if gauss6 is not None:
        if gauss5 is None or gauss6['fixed_pressure_mass']!=mass:
            raise ArithmeticError('Need all seven reference charts and their common mass for Gauss6')
        u=values(gauss6,'Gauss6')
        pairs.update(Gauss6_over_W1=u/a-1,Gauss6_over_W2=u/b-1,Gauss6_over_Gauss2=u/c-1,
                     Gauss6_over_W4=u/d-1,Gauss6_over_Gauss3=u/h-1,Gauss6_over_Gauss4=u/j-1,Gauss6_over_Gauss5=u/k-1)
        reference=(*reference,u)
        scope='Last-segment time-step contrast on eight separately archived warm chart restarts; no full-history error bound.'
    return dict(fixed_mass=mass,columns=('number','pressure','warm_pressure'),steps=(1024,2048),
                individual_coarse_over_fine_minus_one=[(r[0]/r[1]-1).tolist() for r in reference],
                same_clock_relative_contrasts={k:v.tolist() for k,v in pairs.items()},
                fine_minus_coarse_contrast={k:(v[1]-v[0]).tolist() for k,v in pairs.items()},
                scope=scope)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('w1',type=Path)
    p.add_argument('w2',type=Path)
    p.add_argument('gauss',type=Path)
    p.add_argument('--w4',type=Path)
    p.add_argument('--gauss3',type=Path)
    p.add_argument('--gauss4',type=Path)
    p.add_argument('--gauss5',type=Path)
    p.add_argument('--gauss6',type=Path)
    a=p.parse_args()
    print(json.dumps(run(*map(read,(a.w1,a.w2,a.gauss)),w4=read(a.w4) if a.w4 else None,gauss3=read(a.gauss3) if a.gauss3 else None,gauss4=read(a.gauss4) if a.gauss4 else None,gauss5=read(a.gauss5) if a.gauss5 else None,gauss6=read(a.gauss6) if a.gauss6 else None),indent=2,allow_nan=False))
