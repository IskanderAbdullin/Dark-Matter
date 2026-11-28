import argparse,json
import numpy as np
from scipy.special import expit
from scipy.linalg import block_diag
from early_pole_archived_cold_screened_elastic_reference import Screened


def study(grid,level):
 m=Screened(grid,1,8,8e-8,[(1,1,2),(1,2,3)][level],.08,1.)
 M=m.M
 Y=M.init.copy()
 r0=M.r[1].start
 for _ in range(32):
  z=Y.copy()
  _,_,E,_,_=M.weights(Y,True)
  for j in range(r0,len(Y)):
   mu=.1 if j<M.r[2].start else -.4
   c=M.cells[j]
   Y[j]=c.I@expit(-E[j,0]-E[j,1]*c.u+mu)
  if np.max(abs(Y-z))<3e-15:break
 e,F,E,_,_=M.weights(Y,True)
 ix,u,pb,Q,W=m.graph
 d=M.r[2]
 val=e[ix,0]+e[ix,1]*u
 ec=val[0]+val[1]-val[2]
 j=d.start+np.searchsorted(e[d,0]+e[d,1],ec)
 bounded=j<d.stop
 j=np.minimum(j,d.stop-1)
 v=(ec-e[j,0])/e[j,1]
 q=M.cells[d.start].pfrombare(M.shape[j,0]+M.shape[j,1]*v)
 ok=bounded&(abs(v)<=1.)&(abs(pb-q)<=Q)&(Q<=pb+q)
 v=np.where(ok,v,0.)
 a,b,c=expit(F[ix,0]+F[ix,1]*u)
 z=expit(F[j,0]+F[j,1]*v)
 wh=m.factor*W*ok*np.sqrt(a*b*(1-c)*(1-z)*c*z*(1-a)*(1-b))
 idx=np.stack([*ix,j])-r0
 mi=np.stack([*u,v])
 sin=np.array([-1.,-1.,1.,1.])[:,None]
 B=np.stack([sin*np.ones_like(mi),sin*mi],axis=1).reshape(8,-1)
 ind=(2*idx[:,None,:]+np.array([0,1])[None,:,None]).reshape(8,-1)
 size=(len(Y)-r0)*2
 rows=[]
 H=[]
 for k,c in enumerate(M.cells[r0:],start=r0):
  p=expit(F[k]@c.B)
  val=(c.I*(p*(1-p)))@c.B.T
  eig,U=np.linalg.eigh(val)
  H.append((U/np.sqrt(eig))@U.T)
 H=block_diag(*H)
 C=np.zeros((size,3));C[:,:2]=np.repeat(np.eye(2),[(M.r[1].stop-r0)*2,(M.r[2].stop-M.r[2].start)*2],axis=0)
 C[1::2,:2]=0.
 C[:,2]=E[r0:].ravel()
 for mode in ('all','screened_part','vacuum_part'):
  w=wh*(1 if mode=='all' else Q<.08 if mode=='screened_part' else Q>=.08)
  K=np.zeros((size,size))
  for i in range(8):
   for j in range(i+1):
    k=np.bincount(ind[i]*size+ind[j],weights=w*B[i]*B[j],minlength=size*size).reshape(size,size)
    K+=k if i==j else k+k.T
  # Pauli susceptibility weighted relaxation, with gauge remapping held fixed.
  norm=np.linalg.norm(K,2)
  residual=np.linalg.norm(K@C,axis=0)/np.linalg.norm(C,axis=0)/norm
  assert np.max(residual)<2e-12
  rate=np.linalg.eigvalsh(H@K@H)
  assert min(rate)>-2e-12*max(rate)
  rows.append(dict(mode=mode,largest_rate=rate[-1],relative_invariant_residual=residual.tolist(),
                   slowest_positive_after_nulls=rate[3],null_values=rate[:3].tolist()))
 return dict(grid=grid,level=level,scope='Common-T frozen geometry linearization; not non-equilibrium eigenvalues.',rows=rows)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('grid');p.add_argument('--level',type=int,default=0)
 a=p.parse_args();print(json.dumps(study(a.grid,a.level),indent=2,allow_nan=False))
