"""Diagnostic, number-matched occupations beyond thermal capacity and at vacuum."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from thermal_cell_inversion_reference import infer_shift, shell_average


def cell_shape(low,high,mass,temp,target,statistics):
    """Return occupation function and closure tag; never claim a Bose condensate."""
    if not 0 <= low < high or mass <= 0 or temp <= 0 or not np.isfinite(target):
        raise ValueError("invalid state")
    if statistics not in ("FD","BE") or target<0 or (statistics=="FD" and target>1):
        raise ValueError("invalid cell occupation")
    if target==0 or (statistics=="FD" and target==1):
        return lambda q: np.full_like(np.asarray(q,float),target),"exact-boundary"
    scale=1.
    kind="thermal-inversion"
    if statistics=="BE":
        cap=shell_average(low,high,mass,temp,0.,"BE")
        if target>cap:
            scale=target/cap
            shift=0.
            kind="artificial-saturated-shape-rescaled"
        else:
            shift=infer_shift(low,high,mass,temp,target,statistics)
    else:
        shift=infer_shift(low,high,mass,temp,target,statistics)
    def shape(q):
        q=np.asarray(q,float)
        if np.any((q<=0)|(q<low)|(q>high)):
            raise ValueError("evaluate finite positive samples inside the cell only")
        kinetic=q*q/(np.hypot(mass,q)+mass)
        if statistics=="FD":
            return expit((shift-kinetic)/temp)
        with np.errstate(over="ignore"):
            return scale/np.expm1((kinetic-shift)/temp)
    return shape,kind


def run():
    rows=[]
    for low,high in [(0.,1e-7),(1e-5,1e-4),(1e-4,.1),(.1,3.)]:
        mass,T=.0004,1.
        cap=shell_average(low,high,mass,T,0.,"BE")
        vf=(high*high+low*high+low*low)/3
        for stat, values in [("BE",(0.,.5*cap,cap,np.nextafter(cap,np.inf),10*cap)),
                             ("FD",(0.,1e-13,.2,1-1e-12,1.))]:
            for value in values:
                f,l=cell_shape(low,high,mass,T,value,stat)
                q=low+(high-low)*np.linspace(1e-5,.99999,100)
                y=f(q)
                assert y.min()>=0 and (stat=="BE" or y.max()<=1.)
                def g(s):
                    q=low+(high-low)*s
                    return q*q*f(q)/vf
                recovered=quad(g,0.,1.,epsabs=0.,epsrel=2e-10)[0]
                err=abs(recovered-value)/value if value else abs(recovered)
                assert err<2e-9
                rows.append(dict(statistics=stat,qlo=low,qhi=high,kind=l,target=value,
                                 relative_number_residual=float(err)))
    # Two extreme diagnostic shapes with the same N need not have the same E.
    hi,mass,T=1e-7,.0004,1.
    value=10*shell_average(0.,hi,mass,T,0.,"BE")
    sat,_=cell_shape(0.,hi,mass,T,value,"BE")
    def energy(shape):
        return quad(lambda z: z*z*float(shape(hi*z))*
                    (hi*z)**2/(np.hypot(hi*z,mass)+mass),
                    0.,1.,epsabs=0.,epsrel=2e-10)[0]
    ratio=energy(sat)/energy(lambda _:value)
    assert abs(ratio-5/9)<5e-8
    return dict(equal_number_kinetic_energy_ratio_saturated_over_flat=ratio,cases=len(rows),max_relative_number_residual=max(x["relative_number_residual"] for x in rows),
                samples=rows)

if __name__ == "__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
