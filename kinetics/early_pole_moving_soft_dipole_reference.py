"""Doppler-sampled longitudinal and transverse parts of a local pair correlator."""
import json
import numpy as np
from scipy.special import roots_legendre
from early_pole_static_longitudinal_dipole_reference import static_shape
from early_pole_differential_phase_separation_reference import geometry


def weights(v,y,n=150,zmax=8000.):
    # k measured in m_D. Frequency is omega=k.v, not zero.
    # Coulomb-gauge longitudinal retarded denominator only.
    mu,w=roots_legendre(n)
    mu=(mu+1)/2
    x=v*mu
    R=1-x*np.arctanh(x)
    I=np.pi*x/2
    grid=np.r_[0.,np.geomspace(.003,zmax,210)]
    zz,ww=roots_legendre(32)
    k=(grid[1:]+grid[:-1])[:,None]/2+np.diff(grid)[:,None]*zz/2
    wk=np.diff(grid)[:,None]*ww/2
    k=k.ravel()[:,None]
    dw=wk.reshape(-1,1)*w/2
    W=2*k*dw/((k*k+R)**2+I**2)
    # Parallel separation and arithmetic average over growing same-axis length.
    c=1-np.sinc(k*y*mu/np.pi)
    reT=.5*(x*x+x*(1-x*x)*np.arctanh(x))
    imT=np.pi*x*(1-x*x)/4
    numT=.5*(1-x*x)*v*v*(1-mu*mu)
    WT=2*k*dw*numT/((k*k*(1-x*x)+reT)**2+imT**2)
    return np.array([np.sum(W),np.sum(W*c),np.sum(WT*c)])


def run():
    out=[]
    for g in (.001,.003,.004,.005,.01):
        dv,a=geometry(g)
        y=g*np.sqrt(2/3)*dv/a
        s=static_shape(y,True,100)[0]
        base=weights(0.,y)
        assert abs(base[1]-s)<4e-7 and abs(base[0]-1)<4e-8
        for v in (.5,.9,.999):
            p=weights(v,y)
            error=max(abs(p-weights(v,y,200)))
            assert error<2e-6
            out.append(dict(g=g,v=v,length=y,static_average=s,
                            doppler_relative_longitudinal_independent=p[0],
                            doppler_average_relative_static_independent=p[1],
                            doppler_average_relative_own_independent=p[1]/p[0],
                            exploratory_transverse_paired_relative_static_independent=p[2],
                            quadrature_absolute_change=error))
    return dict(scope='Separate HTL longitudinal and transverse pair terms; local straight common-velocity approximation, '
                      'no endpoint connectors, no physical width; transverse independent IR integral deliberately omitted.',rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
