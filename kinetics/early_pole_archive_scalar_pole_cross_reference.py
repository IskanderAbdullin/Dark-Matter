"""Cross-check guarded Yukawa integrals by reversing the two invariant integrals."""
import argparse,gzip,json
import numpy as np
from scipy.special import roots_legendre
from early_pole_archive_scalar_below_pole_finite_q_reference import integrate
from early_pole_native_aperture_witness_reference import GRID

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('grid',choices=('fine','third'));p.add_argument('source');p.add_argument('--order',type=int,default=12)
 args=p.parse_args()
 with gzip.open(GRID/f'early-pole-spectrum-{args.grid}-m06-2026-11-26.json.gz','rt') as f:a=json.load(f)
 with open(args.source) as f:data=json.load(f)
 nodes,weights=roots_legendre(8);out=[]
 for case in data['cases']:
  hi,lo=case['pole_guard'],case['lower_invariant_fraction']
  if (lo,hi) not in ((.9999,.99999),(1.00001,1.0001)):continue
  vhi,vlo=1/(hi**2-1),1/(lo**2-1)
  h=(vlo-vhi)/2
  v=(vhi+vlo)/2+h*nodes
  y=np.zeros((2,3))
  for vi,w in zip(v,h*weights):
   z=integrate(a,case['q'],args.order,.90,probe_x=1/vi)
   y+=w*np.c_[z['gross_density_per_reference_H'],z['net_density_per_reference_H']]
  ref=np.c_[case['gross_by_band_per_reference_H'],case['net_by_band_per_reference_H']]
  relative=abs(ref-y)/np.maximum(abs(ref),1e-50)
  assert np.max(relative[1])<6e-6 and np.max(relative[:,:2])<.003
  out.append(dict(q=case['q'],annulus=[lo,hi],direct_annulus=ref.tolist(),invariant_slice_quadrature=y.tolist(),
                  relative_difference=relative.tolist()))
 assert len(out)==4, 'Use the matched q=0 and q=4.1e-6 order-14 annuli output.'
 print(json.dumps(dict(grid=args.grid,results=out),indent=2,allow_nan=False))
