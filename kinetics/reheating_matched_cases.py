"""Reproduce the selected near-abundance-matched early-start cases."""
import argparse,json
from pathlib import Path
from energy_closed_bath import EnergyClosedBath

CASES=[(.05,320,1e-7,.09521813777831052),(.02,320,1e-7,.0553),
       (.01,320,1e-8,.0353),(.005,320,2.5e-9,.0217),
       (.005,480,2.5e-10,.0217),(.002,320,4e-10,.0109)]

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',nargs='+',type=int,default=list(range(len(CASES))))
    p.add_argument('--output',required=True);a=p.parse_args();out=[]
    for i in a.case:
        x,n,q,k=CASES[i]
        m=EnergyClosedBath(n=n,qmin=q,kappa=k,relaxation=100.)
        r=m.evolve_closed(x,40.,5e-7,True)
        s=r['final'];assert abs(s['omega_phi_h2_m86']/.12-1)<.01
        print(i,json.dumps(s),flush=True);out.append(r)
    Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
