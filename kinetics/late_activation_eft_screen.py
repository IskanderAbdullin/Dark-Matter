"""Illustrative VEV/cutoff/uncancelled vacuum energy inequalities."""
import argparse,json,math
from pathlib import Path

def screen(z=3.,m=86.):
    Tgamma=2.7255*8.617333262e-5*(1+z)
    rho=3*(2.1927e-18*6.582119569e-16)**2*(2.435e27)**2*(.315*(1+z)**3+.685)
    rows=[]
    for sigma in [.1,1.25]:
        lam=math.sqrt(128*math.pi*m**3*sigma*1.7826619216278976e-33/(1.973269804e-5)**2)
        for ratio in [1.,10.,100.]:
            cutoff=ratio*m;vev=math.sqrt(lam)*cutoff
            rows.append(dict(sigma=sigma,lambda_late=lam,cutoff_eV=cutoff,vev_eV=vev,vev_over_Tgamma=vev/Tgamma,
                max_chi_quartic_for_one_percent_total_vacuum=4*.01*rho/vev**4))
    return dict(z=z,Tgamma_eV=Tgamma,total_density_eV4=rho,
        assumptions='illustrative coefficient lambda_eff=VEV^2/Lambda^2, vacuum scale lambda_chi*VEV^4/4, no cancellation, no specified finite-temperature potential or UV completion',rows=rows)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    Path(a.output).write_text(json.dumps(screen(),indent=2)+'\n',encoding='utf-8')
