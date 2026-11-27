"""Leading soft transverse noise from an isotropic Fermi archive, not a photon KMS assumption."""
import gzip,json
from pathlib import Path
import numpy as np
from scipy.special import expit
from scipy.integrate import quad
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import PREF,Species,gl
from finite_mass_qed_frozen_transverse_cut_reference import response


def archive(ch,u,n):
    m=ch.m
    p0=m*u/np.sqrt((1-u)*(1+u))
    lo=np.minimum(np.maximum(ch.e[:-1],p0),ch.e[1:])
    hi=ch.e[1:]
    L,H=np.hypot(lo,m),np.hypot(hi,m)
    fundamental=(H-L)*((1-u*u)*(H*H+H*L+L*L)/3-m*m)
    f=ch.f
    noise=(f*(1-f))@fundamental
    dissip=(1-u*u)*(f@((hi-lo)*(hi+lo)))
    x,w=gl(n)
    mid=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*x/2
    E=np.hypot(mid,m)
    independent=((f*(1-f))*(hi-lo)/2)@((mid/E*((1-u*u)*mid**2-u*u*m*m))@w)
    defect=abs(noise-independent)/max(noise,1e-30)
    return noise,dissip,defect


def check_thermal(u,T,m):
    a=m*u/np.sqrt((1-u)*(1+u))
    def e(p):return np.hypot(p,m)
    def p0(p):return expit(-e(p)/T)
    def W(p):return p/e(p)*(p*p*(1-u*u)-m*m*u*u)
    noise=quad(lambda p:W(p)*p0(p)*(1-p0(p)),a,a+60*T,epsabs=1e-14,epsrel=3e-12)[0]
    dissip=2*(1-u*u)*quad(lambda p:p*p0(p),a,a+60*T,epsabs=1e-14,epsrel=3e-12)[0]
    return noise/dissip/T-1.


def run():
    out=[]
    alpha=(.001,.0047,.025,.1,.25,.5,.75,.92,.99,.9996)
    ref=max(abs(check_thermal(u,T,m)) for u in alpha for T in (.4,1.2) for m in (.003,.05,.2))
    if ref>2e-10:raise ArithmeticError(('Massive Fermi equilibrium noise limit',ref))
    for g in ('fine','third'):
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        x=d['x_end'];m=d['mass_ratio']*x
        ch=[Species(ma,d[f'endpoint_{n}_momentum_edges'],d[f'endpoint_{n}_occupations'])
            for ma,n in ((x,'parent'),(m,'fermion'))]
        rows=[]
        for u in alpha:
            diag=np.array([archive(c,u,40) for c in ch])
            N,D=diag[:,:2].sum(axis=0)
            _,imag,_=response(ch,u)
            corr=abs(PREF*np.pi*u/4*D/imag-1.)
            if corr>5e-13 or max(diag[:,2])>5e-11 or min(N,D)<=0:
                raise ArithmeticError(('Landau noise independent primitive',g,u,corr,diag))
            rows.append(dict(phase_speed=u,leading_noise_over_response_energy=float(N/D),
                             separate_species_energy_ratios=(diag[:,0]/diag[:,1]).tolist(),
                             relative_defect_in_retarded_imaginary_normalization=float(corr),
                             cell_primitives_max_rel_defect=float(max(diag[:,2]))))
        out.append(dict(grid=g,rows=rows))
    return dict(scope='Leading collisionless transverse soft noise; not an imposed KMS law, finite-frequency Wightman completion, or dressed scalar rate.',
                common_FD_temperature_control_max_relative_defect=ref,data=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
