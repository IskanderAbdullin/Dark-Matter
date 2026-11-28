"""Same physical soft photon resolution in the plasma and parent charts."""
import gzip,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
from finite_mass_qed_signed_real_soft_recoil_reference import setup,angular_scaled

D=Path(__file__).parent/'results'
coeffs=json.loads((D/'finite-mass-qed-signed-real-soft-recoil-2026-11-27.json').read_text())['rows']
bases=json.loads((D/'finite-mass-qed-soft-subtraction-2026-11-27.json').read_text())['angular']
rows=[]
for data,b in zip(coeffs,bases,strict=True):
    g=b['grid']
    with gzip.open(D/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:h=json.load(f)
    M=h['x_end'];m=M*h['mass_ratio'];s=h['checkpoints'][-1]['scalar_pole_over_T']
    c=np.asarray(data['fine']['normalized_soft_coefficients'])
    u=b['plasma_parent_momentum']/M
    gam=np.hypot(u,1)
    beta=u/gam
    x,w,p,E,B=setup(M,m,s,52)
    norm=(p/M)*4*M*(E+m)
    const=float(w@angular_scaled(0.,x,M,m,s,np.sqrt)/norm)
    def G(r):
        return (angular_scaled(r[:,None],x[None,:],M,m,s,np.sqrt)@w)/norm
    def angle_rate(z,system,n):
        v,k=roots_legendre(n)
        if system=='rest':
            L=-np.arcsinh(u)
            H=-L
            ell=(L+H)/2+(H-L)*v/2
            d=np.exp(ell)
            weight=(H-L)/4*k*d/u
        else:
            L=np.log((1/gam**2)/(1+beta)); H=np.log1p(beta)
            ell=(L+H)/2+(H-L)*v/2
            yy=np.exp(ell)
            d=1/(gam*yy)
            weight=(H-L)/4*k*yy/beta*d*d
        r=z/(M*d)
        small=r<=1e-3
        rr=r[small]
        h1=np.zeros(n)
        h2=np.zeros(n)
        h1[small]=sum(c[j]*rr**j for j in range(1,5))
        h2[small]=2*(c[2]*rr**2+c[4]*rr**4)
        h1[~small]=G(r[~small])-const
        h2[~small]=(G(r[~small])+G(-r[~small])-2*const)
        return np.array([weight@h1,(weight@h2)/np.expm1(z/.9)])/z
    runs=[]
    for cut in (2.5e-6,5e-6,1e-5):
        a=np.array([quad(lambda z:angle_rate(z,'rest',56)[i],0,cut,epsabs=8e-12)[0] for i in range(2)])
        r=np.array([quad(lambda z:angle_rate(z,'lab',76)[i],0,cut,epsabs=8e-12)[0] for i in range(2)])
        err=max(abs(a-r))
        if err>2e-11:
            raise ValueError('Photon bath/cut charts disagree.')
        # Isotropic algebraic leading moments after holding fixed the same lab cutoff.
        rapid=np.arcsinh(u)
        opposite=np.exp(rapid)
        first=rapid/u
        second=(opposite-1/opposite)/(2*u)
        if abs(second-1)>2e-14:
            raise ValueError('Doppler squared inverse average lost its exact normalization.')
        leading=[c[1]*first*cut/M,2*.9*c[2]*cut/M**2]
        runs.append(dict(lab_photon_cap=cut,rest_angle=a.tolist(),plasma_angle=r.tolist(),
                         chart_abs_defect=float(err),small_cap_leading=leading))
    rows.append(dict(grid=g,parent_gamma=float(gam),inverse_doppler_mean=float(first),
                     inverse_doppler_squared_mean=float(second),comparisons=runs))
print(json.dumps(dict(scope='Isotropic lab photon bath, all scalar angles integrated, no finite virtual remainder.',
                      columns=['vacuum_recoil','even_thermal_recoil'],rows=rows),
                 indent=2,allow_nan=False))
