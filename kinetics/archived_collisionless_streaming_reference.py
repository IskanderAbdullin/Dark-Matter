"""Conditional single-particle streaming of represented archived distributions."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import constants
from scipy.interpolate import PchipInterpolator
from scipy.special import expit, roots_legendre, zeta

from archived_scalar_pair_reaction_reference import points,snapshot_functions
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import m,mp


class Background:
    def __init__(self,path,T_kelvin=2.72548,H0=67.36,omega_m=.3153):
        t,s,g=np.loadtxt(path).T
        assert t[0]==0 and s[0]==s[1] and g[0]==g[1]
        self.T0=constants.k*T_kelvin/(1e9*constants.e)
        self.gs0=s[0]
        self.gr0=g[0]
        self.H0=H0
        self.om=omega_m
        h_0=H0*1000/(1e6*constants.parsec)
        rho_crit=3*h_0*h_0/(8*np.pi*constants.G)*constants.c**2
        self.ogamma=4*constants.sigma*T_kelvin**4/(constants.c*rho_crit)
        self.or0=self.ogamma*g[0]/2
        self.ol=1-self.om-self.or0
        self.aeq=self.or0/self.om
        self.D0=constants.c/1000/H0
        self.u=np.log(t[1:])
        self.v=np.log(self.T0)+np.log(self.gs0/s[1:])/3-self.u
        self.v_to_u=PchipInterpolator(self.v[::-1],self.u[::-1],extrapolate=False)
        self.gs=PchipInterpolator(self.u,np.log(s[1:]))
        self.gr=PchipInterpolator(self.u,np.log(g[1:]))

    def scale(self,T):
        u=np.log(T)
        gs=self.gs(u) if u>=self.u[0] else np.log(self.gs0)
        return self.T0/T*np.exp((np.log(self.gs0)-gs)/3)

    def weighted_radiation(self,a):
        v=np.log(a)
        lv=np.clip(v,self.v[-1],self.v[0])
        u=self.v_to_u(lv)
        res=self.ogamma*np.exp(self.gr(u)+4*(np.log(self.gs0)-self.gs(u))/3)/2
        return np.where(v>self.v[0],self.or0,res)

    def streams(self,q,mass,ai,af,n):
        assert ai<af and np.all(q>=0.)
        vv=np.log([ai,af])
        # Resolve EOS knots even if the species traverses the MeV period
        # in a small fraction of its eventual comoving path.
        bounds=np.unique(np.r_[vv,np.clip(self.v,*vv),
                                np.linspace(*vv,50)])
        dx,dw=roots_legendre(n)
        l,r=bounds[:-1,None],bounds[1:,None]
        u=((l+r)+(r-l)*dx)/2
        a=np.exp(u.ravel())
        weights=((r-l)/2*dw).ravel()*a
        den=np.sqrt(self.weighted_radiation(a)+self.om*a+self.ol*a**4)
        vel=q[:,None]/np.hypot(q[:,None],mass*a)
        result=self.D0*(vel@(weights/den))
        rmin=self.weighted_radiation(a).min()
        ultra_rad=self.D0*(q/mass)*(np.arcsinh(af*mass/q)-np.arcsinh(ai*mass/q))/np.sqrt(rmin)
        assert np.all(result<=ultra_rad*(1+1e-10))
        return result


def summarize(bg,p,w,ai,mass,order):
    n=w.sum()
    w=w/n
    q=p*ai
    assert bg.scale(bg.T0)<=1.+1e-12
    ae=bg.aeq
    distances=bg.streams(q,mass,ai,ae,order)
    v=q/np.hypot(q,mass*ae)
    cf=(np.cumsum(w)-w/2)
    assert np.all(np.diff(distances)>=0.)
    sa=np.array([.1,.5,.9])
    return dict(initial_number_density_GeV3=n,
        comoving_particle_momentum_over_present_photon_T_mean=float(w@q/bg.T0),
        comoving_particle_momentum_over_present_photon_T_rms=float(np.sqrt(w@(q*q))/bg.T0),
        mean_comoving_distance_to_equality_Mpc=float(w@distances),
        characteristic_quantiles_probability=sa.tolist(),
        free_stream_distance_quantiles_Mpc=np.interp(sa,cf,distances).tolist(),
        rms_speed_at_equality_over_c=float(np.sqrt(w@(v*v))),
        rho_over_rest_density_at_equality=float(w@np.hypot(1.,q/(mass*ae))),
        kinetic_pressure_over_energy_at_equality=float((w@((v*v)*np.hypot(1.,q/(mass*ae))))/(3*(w@np.hypot(1.,q/(mass*ae))))))


def archive(bg,paths):
    out=[]
    for path in paths:
        s=json.loads(Path(path).read_text())["snapshots"][-1]
        T=s["T_visible_GeV"]
        ai=float(bg.scale(T))
        ed=np.array(s["lab_momentum_edges_GeV"])
        fp,_,fc,_=snapshot_functions(s)
        rows=[]
        for rad,ti in ((12,4),(32,7)):
            zz=[]
            for mass,f,g in ((mp,fp,1),(m,fc,4)):
                p,_,w=points(ed,mass,f,rad)
                zz.append(summarize(bg,p,w*g,ai,mass,ti))
            n=np.array([z["initial_number_density_GeV3"] for z in zz])
            ref=np.array(s["represented_species_count_a3"])/s["scale_factor"]**3
            assert max(abs(n/ref-1))<3e-8
            check_T=.001
            if T>=check_T:
                ab=bg.scale(check_T)
                den=np.pi**2/30*check_T**4*np.exp(bg.gr(np.log(check_T)))
                ene=[]
                for mass,f,g in ((mp,fp,1),(m,fc,4)):
                    nodes,_,ww=points(ed,mass,f,rad)
                    ene.append(float((ai/ab)**3*g*(ww@np.hypot(mass,nodes*ai/ab))/den))
            else:
                ene=None
            rows.append(dict(momentum_order=rad,expansion_order=ti,scalar=zz[0],stable=zz[1],
                              collisionless_represented_energy_fraction_of_visible_radiation_at_1_MeV=ene))
        for key in ("scalar","stable"):
            fld="mean_comoving_distance_to_equality_Mpc"
            assert abs(rows[0][key][fld]/rows[1][key][fld]-1.)<3e-7
        out.append(dict(snapshot=Path(path).name,T_visible_GeV=T,a_snapshot=ai,controls=rows))
    return out


def thermal_control(bg,table,mass=5.7e-6):
    fac=conversion(table)["Omega_h2_per_GeV_times_yield"]
    Y=.120/(mass*fac)
    entropy=2*np.pi**2/45*bg.gs0*bg.T0**3
    temperature=(entropy*Y*(4*np.pi**2)/(3*zeta(3,1)*2))**(1/3)
    ai=bg.scale(.05)
    p,_,w=points(temperature/ai*np.r_[0.,np.geomspace(.01,80.,40)],mass,
            lambda e:expit(-np.sqrt(np.maximum((e-mass)*(e+mass),0.))*ai/temperature),14)
    r=summarize(bg,p,w*2,ai,mass,7)
    assert abs(r["initial_number_density_GeV3"]*ai**3/(entropy*Y)-1.)<1e-6
    return dict(species_mass_GeV=mass,present_reference_distribution_temperature_GeV=temperature,
                assumed_full_thermal_Omega_h2=.120,**r)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table");p.add_argument("paths",nargs="+")
    args=p.parse_args()
    bg=Background(args.table)
    out=dict(cosmology=dict(H0_km_per_s_Mpc=bg.H0,matter_fraction=bg.om,
              present_effective_massless_radiation_fraction=bg.or0,equality_scale_factor=bg.aeq),
             controls=archive(bg,args.paths),thermal_reference=thermal_control(bg,args.table),
             scope="Exact geodesic streaming after each finite-support snapshot in prescribed reference radiation/matter/constant background; collisionless, no perturbations or transfer function, no power/halo or Lyman-alpha exclusion. Number-weighted for retained modes only.")
    print(json.dumps(out,indent=2,allow_nan=False))
