"""Forward projection of small halo-core contrasts at the sampled F583-4 scales."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre
from halo_bose_primordial_cap_core_reference import integrate


def halo(d,shift=0.,physical_critical_normalizer=None):
    yg=np.r_[0.,np.geomspace(1e-12,18.,7800),d]
    yg=np.unique(yg)
    n=np.array([integrate(y+shift,d) for y in yg])
    dn=n[0] if physical_critical_normalizer is None else physical_critical_normalizer
    rho=PchipInterpolator(yg,n/dn)
    u0=-18.
    z=np.exp(2*u0)
    s=solve_ivp(lambda u,a:[a[1],np.exp(2*u)*rho(a[0])-a[1]],(u0,np.log(50.)),
                [z*n[0]/(6*dn),z*n[0]/(3*dn)],rtol=4e-11,atol=1e-13,max_step=.023,dense_output=True)
    assert s.success
    return lambda r:np.sqrt(s.sol(np.log(np.maximum(r.ravel(),np.exp(u0))))[1]).reshape(r.shape),float(n[0])


def make_projection(Rs, inc, step_radius, emission_scale, order):
    """A filled circular aperture on the sky, not the source instrument pipeline."""
    if step_radius==0.:
        return Rs[:,None],np.ones((len(Rs),1))
    x,w=roots_legendre(order)
    u=(x+1)/2
    ang=2*np.pi*(np.arange(4*order)+.5)/(4*order)
    Dx=(step_radius*u[:,None]*np.cos(ang)).ravel()
    Dy=(step_radius*u[:,None]*np.sin(ang)).ravel()
    wg=np.broadcast_to((u*w)[:,None],(order,4*order)).ravel()
    sky=Rs[:,None]+Dx[None,:]
    rad=np.hypot(sky,Dy[None,:]/np.cos(inc))
    tracer=np.exp(-rad/emission_scale)*wg
    W=tracer/np.sum(tracer,axis=1,keepdims=True)
    return rad,W*sky/np.maximum(rad,1e-30)


def run(path):
    cat=json.loads(path.read_text())
    row=np.asarray(cat['table'],float)
    Rs,er=row[:,0],row[:,2]
    D=cat['distance_Mpc']
    ik=np.deg2rad(cat['inclination_deg'])
    acs=D*1e3/206264.8062471
    for i in [0,1]:
        assert row[i,1]>0 and er[i]>0
    # These are fiducial source benchmark scales, not fits to the selected table.
    RR=.6*D/49.
    vv=44.1
    vg=PchipInterpolator(np.r_[0.,Rs],np.r_[0.,row[:,3]])
    vd=PchipInterpolator(np.r_[0.,Rs],np.r_[0.,row[:,4]])
    vals,w=make_projection(Rs,ik,1.5*acs,1e12,24)
    assert max(abs(np.sum(vals*w,axis=1)-Rs))<1e-10
    C,n0=halo(0.)
    cases=[]
    for F in (100.,5000.,7000.,200000.):
        d=np.log1p(1/F)
        P,n1=halo(d)
        R1=RR*np.sqrt(n0/n1)
        case=[]
        scaling_check=None
        if F==7000.:
            verify,_=halo(d,physical_critical_normalizer=n0)
            scaling_check=float(np.max(abs(verify(Rs/RR)-P(Rs/R1)))*vv)
            assert scaling_check<2e-7,scaling_check
        for r_ap,light in ((0.,1.93),(1.5,.5),(1.5,1.93),(1.5,4.)):
            r,w=make_projection(Rs,ik,r_ap*acs,light,32)
            g=vg(np.minimum(r,Rs[-1]))**2
            b=vd(np.minimum(r,Rs[-1]))**2
            def obs(fun,rr,vel=44.1,ups=.5):
                return np.sum(np.sqrt(vel**2*fun(r/rr)**2+g+ups*b)*w,axis=1)
            if F==7000. and r_ap:
                re,we=make_projection(Rs,ik,r_ap*acs,light,48)
                ge=vg(np.minimum(re,Rs[-1]))**2
                be=vd(np.minimum(re,Rs[-1]))**2
                xc=np.sum(np.sqrt(vv**2*C(re/RR)**2+ge+.5*be)*we,axis=1)
                yp=np.sum(np.sqrt(vv**2*P(re/R1)**2+ge+.5*be)*we,axis=1)
                assert np.max(abs((yp-xc)-(obs(P,R1)-obs(C,RR))))<5e-6
            x=obs(C,RR)
            y=obs(P,R1)
            s=y-x
            aux=[(obs(C,RR*np.exp(1e-4))-obs(C,RR*np.exp(-1e-4)))/2e-4,
                 (obs(C,RR,vv*np.exp(1e-4))-obs(C,RR,vv*np.exp(-1e-4)))/2e-4,
                 (obs(C,RR,ups=.5001)-obs(C,RR,ups=.4999))/2e-4]
            mat=np.array(aux).T/er[:,None]
            drift=np.linalg.lstsq(mat,s/er,rcond=None)[0]
            orth=s/er-mat@drift
            rot=obs(P,RR)-x
            case.append(dict(aperture_radius_arcsec=r_ap,illustrative_emission_scale_kpc=light,
                             point_radius_kpc=Rs.tolist(),
                             model_velocity_difference_same_mass_and_Vscale_km_s=s.tolist(),
                             diagnostic_same_dimensionless_R_difference_km_s=rot.tolist(),
                             max_abs_model_difference_km_s=float(np.max(abs(s))),
                             optimistic_diagonal_contrast_norm=float(np.linalg.norm(s/er)),
                             remaining_contrast_norm_after_nuisance_projection=float(np.linalg.norm(orth))))
        cases.append(dict(F=F,rounding_scale_pc=1e3*np.sqrt(6*d)*R1,
                          rounding_scale_arcsec=np.sqrt(6*d)*R1/acs,
                          physical_core_scale_kpc=R1,
                          independent_fixed_physical_radius_form_km_s=scaling_check,cases=case))
    return dict(scope="Detector forward control at public sampling and error scales, not a fit or "
                      "reanalysis of the catalog velocities. Synthetic uniform circular aperture "
                      "and idealized exponential emission; no PSF, fiber layout, gas pressure, "
                      "inclination covariance or line profiles. Point errors treated independently "
                      "only to set a scale. Fixed source benchmark velocities/radii and illustrative M/L.",
                galaxy=cat['galaxy'],physical_arcsec_kpc=acs,weakest_sample_radius_kpc=Rs[0],
                benchmark_R_at_table_D_kpc=RR,V0_km_s=vv,outputs=cases)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('data',type=Path)
    print(json.dumps(run(p.parse_args().data),indent=2,allow_nan=False))
