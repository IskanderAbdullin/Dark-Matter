"""Reference finite-mass decay phase-space rectangle integration (not an evolution solver)."""
import math
import numpy as np
from scipy.integrate import quad


def bounds(E,m,b,K,D):
    p=math.sqrt(max((E-m)*(E+m),0.))
    # Rationalized lower root avoids subtraction of almost equal terms.
    lo=(b*b*E*E+D*D)/(K*E+D*p)
    hi=(K*E+D*p)/(m*m) if m else math.inf
    return lo,hi


def cell_area(elo,ehi,rlo,rhi,x,m,b):
    """Allowed area in daughter energies E,r; masses and energies in Tref units."""
    if x<=m+b:return 0.
    if elo<m or rlo<b:raise ValueError('cell extends below rest mass')
    K=(x*x-m*m-b*b)/2;D=math.sqrt((K-m*b)*(K+m*b))
    # Split wherever a kinematic bound intersects an r cell edge.
    points=[]
    for r in (rlo,rhi):
        if r>0:
            for e in bounds(r,b,m,K,D):
                if elo<e<ehi:points.append(e)
    def fun(E):
        lo,hi=bounds(E,m,b,K,D)
        return max(0.,min(rhi,hi)-max(rlo,lo))
    return quad(fun,elo,ehi,points=sorted(set(points)),epsabs=1e-15,epsrel=1e-10,limit=100)[0]


def verify():
    from massive_daughter_bath import allowed_area
    rng=np.random.default_rng(7321);max_sym=0.;max_old=0.;nonzero=0
    for x,m,b in [(0.1,0.,0.),(.1,.003,0.),(.1,.003,.02),(1.,.2,.4),(1.,.49,.49)]:
        for _ in range(300):
            elo,ehi=np.sort(m+np.exp(rng.uniform(-8,2,2)))
            rlo,rhi=np.sort(b+np.exp(rng.uniform(-8,2,2)))
            a=cell_area(elo,ehi,rlo,rhi,x,m,b)
            c=cell_area(rlo,rhi,elo,ehi,x,b,m)
            scale=(ehi-elo)*(rhi-rlo)
            assert -1e-14<=a<=scale+1e-14
            nonzero+=a>0
            max_sym=max(max_sym,abs(a-c)/scale)
            if b==0:
                d=allowed_area(elo,ehi,rlo,rhi,x,m)
                max_old=max(max_old,abs(a-d)/scale)
    assert max_sym<1e-6 and max_old<1e-6
    assert cell_area(1.,2.,1.,2.,1.,.5,.5)==0
    return dict(cells=1500,nonzero_cells=nonzero,max_symmetry_scaled_error=max_sym,max_previous_kernel_scaled_error=max_old)


if __name__=='__main__':
    import json
    print(json.dumps(verify()))
