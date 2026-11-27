"""Inspect the parity and small-offset regularity of an archived pole density."""
import argparse,gzip,json
import numpy as np
from early_pole_native_aperture_witness_reference import GRID
from early_pole_archive_scalar_pole_core_reference import density,conserved_ramp

def regularity_fixture():
 from scipy.integrate import quad
 R=.002
 for r in (2.5e-11,2.24e-10,2e-4):
  def bw(x):return 1/(x*x+r*r)
  def pv(x):return (x*x-r*r)/(x*x+r*r)**2
  def integral(kernel,k):
   # Integrate in a logarithmic coordinate without a narrow peak.
   hi=np.log1p(R/r)
   def integrand(z):
    x=r*np.expm1(z)
    return x**k*kernel(x)*r*np.exp(z)
   return quad(integrand,0.,hi,epsabs=2e-13,epsrel=2e-13)[0]
  cb=.5*np.log1p((R/r)**2)
  cp=cb+r*r/(R*R+r*r)-1
  # Local constant subtraction uses only the area in the window. It is
  # distinct from subtracting the global delta weight; their difference
  # is minus the exterior constant tail.
  tail=-2*np.arctan(r/R)/r
  numeric_tail=-2*quad(lambda u:1/(1+(r*u)**2),0.,1/R,epsabs=1e-12)[0]
  np.testing.assert_allclose(tail,numeric_tail,rtol=2e-14)
  sb=R-r*np.arctan(R/r)
  sp=sb-r*(np.arctan(R/r)-r*R/(R*R+r*r))
  np.testing.assert_allclose([integral(bw,1),integral(pv,1),integral(bw,2),integral(pv,2)],
                             [cb,cp,sb,sp],rtol=2e-12,atol=2e-15)

def run(g):
 with gzip.open(GRID/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as f:
  orig=json.load(f)
 out=[]
 for n in (0,8,16,32,64):
  a=conserved_ramp(orig,g,.03,n)[0] if n else orig
  for o in (18,32):
   q=4.1e-6
   x0=density(a,q,o,0)[0]
   rows=[]
   for x in (2e-3,1e-3,2e-4,1e-4,2e-5,1e-5,2e-6,1e-6,2e-7,2e-8,2e-9):
    plus,minus=(density(a,q,o,s*x)[0] for s in (1,-1))
    even=(plus-x0)+(minus-x0)
    rows.append(dict(x=x,g_over_x=(even/x).tolist(),g_over_x2=(even/(x*x)).tolist(),
                     pos_minus_center=(plus-x0).tolist(),neg_minus_center=(minus-x0).tolist()))
   out.append(dict(ramp_subcells=n,order=o,pole=x0.tolist(),samples=rows))
 return dict(grid=g,tagged_q=4.1e-6,width_convention='Density only; no width inserted.',cases=out)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('grid',choices=('fine','third'));o=p.parse_args()
 regularity_fixture()
 print(json.dumps(run(o.grid),indent=2,allow_nan=False))
