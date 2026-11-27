"""A common soft kick model: compare its characteristic exponent to diffusion."""
import argparse
import json
import numpy as np
from scipy.special import roots_legendre

from early_pole_eikonal_soft_kick_sum_rules_reference import kernel,log_nodes
from early_pole_matched_wilson_geometry_reference import inputs
from early_pole_moving_soft_dipole_reference import weights


def sphere(c,n,v):
    k,kw=log_nodes(2e-5,c,n,60)
    a,w=roots_legendre(n)
    # Log panels resolve the quasi-static transverse component for soft k.
    u,uw=log_nodes(2e-14,.25,n,70)
    left=np.linspace(.25,.98,31)
    h=np.diff(left)[:,None]/2
    um=(left[1:]+left[:-1])[:,None]/2+h*a
    we=(h*w).ravel()
    e,ew=log_nodes(1e-10,.02,n,18)
    u,uw=np.r_[u,um.ravel(),1-e],np.r_[uw,we,ew]
    k=k[:,None];u=u[None,:]
    x=v*u
    d=1-x*x
    h=np.arctanh(x)
    R=1-x*h
    P=.5*(x*x+x*d*h)
    L=1/((k*k+R)**2+(np.pi*x/2)**2)
    T=.5*d*v*v*(1-u*u)/((k*k*d+P)**2+(np.pi*x*d/4)**2)
    W=2*k*(L+T)*kw[:,None]*uw[None,:]
    return W,k*u,k*np.sqrt(1-u*u)


def evalset(kind,c,g,n):
    D=inputs(g)
    v=D['vP'];a=D['angular_halfspan'];mD=g*np.sqrt(2/3)
    W,kz,kp=kernel(c,n) if kind=='transverse' else sphere(c,n,v)
    M2=np.sum(W*kz*kz)
    MT=np.sum(W*kp*kp)
    out=[]
    for u in (.025,.1,.5,1.,2.,4.):
        y=u*D['r']
        arg=y*kz
        ph=np.sum(W*(1-np.sinc(arg/np.pi)))
        diff=y*y*M2/6
        assert 0.<ph<=diff
        pre=g*g*u/(4*np.pi*a)
        lcut=kz>=1/D['r']
        exact_frac=np.sum(W*(1-np.sinc(arg/np.pi))*lcut)/ph
        m2w=W*kz*kz
        moment_frac=np.sum(m2w*lcut)/M2
        QL=g*g*mD*mD*M2/(4*np.pi)
        QT=g*g*mD*mD*MT/(4*np.pi)
        dp=np.sqrt(QL*u/a)
        point=1.4358
        eta=QL/(2*point)+(QT-2*QL)/(2*point*point)
        H_over_T=np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18
        if kind=='sphere':
            # Algebraically the same Doppler-support integral as the local
            # straight-path soft phase experiment. Do not add them.
            ref=sum(weights(v,y,zmax=c)[1:])
            assert abs(ref/ph-1)<.001,(ref,ph,g,u,c)
        else:
            ref=None
        out.append(dict(kind=kind,g=g,cut=c,u=u,characteristic_chi=pre*ph,
                        diffusion_chi=pre*diff,diffusion_over_kick=diff/ph,
                        local_phase_shape_normalized=ref,
                        radial_rms_in_T=dp,radial_halfspan_crossings=dp/(mD/D['r']),
                        transverse_rms_angle=np.sqrt(QT*u/a)/point,
                        one_line_hat_q_longitudinal_over_T3=QL,
                        one_line_hat_q_transverse_over_T3=QT,
                        reference_drag_over_H=eta/H_over_T,
                        reference_drag_over_angular_halfspan=eta/a,
                        fraction_M2_beyond_one_angle_detuning=moment_frac,
                        fraction_phase_beyond_one_angle_detuning=exact_frac))
    return out


def run(n):
    rows=[]
    for g in (.001,.003,.004,.005):
        for c in (5.,10.,20.):
            for kind in ('sphere','transverse'):
                rows+=evalset(kind,c,g,n)
    return dict(scope='Classical independent symmetric kicks without drag, recoil, hard matching or '
                      'nonthermal effects. In the spherical run this is the same local soft field '
                      'as the existing Doppler phase calculation, not an additional damping input.',
                rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nodes',type=int,default=8)
    print(json.dumps(run(p.parse_args().nodes),indent=2,allow_nan=False))
