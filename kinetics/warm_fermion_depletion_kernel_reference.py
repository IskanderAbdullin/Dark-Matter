"""Uniform forward vacuum fermion annihilation bounds for both masses."""
import json
import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import amplitude
from three_species_decay_onsager_reference import m,M,mp,y


def geom(s,initial,exchange):
    return (exchange**4+exchange**2*(s-2*initial**2-2*mp**2)
            +(initial**2-mp**2)**2)


def finite(initial,exchange,n_s,n_x,end):
    s=np.geomspace(4*initial**2,end,n_s+1)
    x=np.linspace(0.,1.,n_x+1)
    zl,zh=x[:-1][None,:]**2,x[1:][None,:]**2
    z=np.clip(.5,zl,zh)
    lo,hi=s[:-1,None],s[1:,None]
    al=exchange**2+lo/2-initial**2-mp**2
    ah=exchange**2+hi/2-initial**2-mp**2
    pl,ph=np.maximum(0.,lo/4-initial**2),hi/4-initial**2
    kl,kh=lo/4-mp**2,hi/4-mp**2
    Ll=(initial+exchange)*al-2*initial*kl*zh
    Lh=(initial+exchange)*ah-2*initial*kh*zl
    den=geom(lo,initial,exchange)+4*pl*kl*(1-zh)
    upper=8*ph*y**4*(hi*kh**2*z*(1-z)+np.maximum(Ll*Ll,Lh*Lh))/den**2
    sigma=(upper @ (2*np.diff(x)))/(32*np.pi*lo[:,0])*(1+1e-11)
    j=np.argmax(sigma)
    return float(sigma[j]),[float(s[j]),float(s[j+1])]


def run():
    rows=[]
    for d,e in ((m,M),(M,m)):
        s=np.geomspace(4*d*d*(1+1e-8),20.,1500)
        x,w=roots_legendre(240)
        exact=(amplitude(s[:,None],x,light=d,heavy=e)@w)/(32*np.pi*s)*np.sqrt(1-4*mp**2/s)
        end=2.
        tail=y**4/(8*np.pi*np.sqrt(geom(end,d,e)))*(1+1e-11)
        r=[]
        for ncell in (220,440,880):
            q,at=finite(d,e,ncell,2*ncell,end)
            r.append(dict(s_cells=ncell,angle_cells=2*ncell,upper_on_finite_s=q,dominant_s_cell=at))
        final=max(tail,r[-1]["upper_on_finite_s"])
        assert max(exact)<final
        rows.append(dict(initial_mass_GeV=d,exchange_mass_GeV=e,split_s_GeV2=end,
                         sampled_peak=float(max(exact)),analytic_tail_upper=float(tail),
                         final_uniform_K_GeVminus2=final,controls=r))
    return dict(scope="Vacuum fermion forward sigma times relative velocity, "
                      "unpolarized, with the identical scalar-pair symmetry. "
                      "Multiply any outgoing Bose factors separately.",
                max_pair_K_GeVminus2=max(x["final_uniform_K_GeVminus2"] for x in rows),cases=rows)


if __name__=="__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
