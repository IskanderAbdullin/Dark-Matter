"""Conditional fixed-bath two-moment evolution with alternative local BE shapes."""
import json
import numpy as np
from scipy.special import expit
from true_event_fv_evolution_reference import energy_mesh, Reconstruction
from cell_number_energy_realizability_reference import fd_bounds
from fermion_two_moment_closure_reference import infer as infer_fd
from bose_two_moment_shape_sensitivity_reference import infer as infer_be
from regular_bose_two_moment_closure_reference import infer as infer_regular_be, occupation


class Model:
    def __init__(self,n,M,order,thermal):
        q,m,b,(i,j,e,r,w)=energy_mesh(n,M,(.08,.18),order)
        self.q,self.n=q,n
        self.m,self.b=m,b
        self.i,self.j,self.e,self.r,self.w=i,j,e,r,w
        self.eb,self.ef=e-m,r-b
        self.B,self.F=Reconstruction(q,m,True),Reconstruction(q,b,False)
        self.thermal=thermal
        self.gf=4/(2*np.pi**2)
        self.Bindex=[np.flatnonzero(i==c) for c in range(n)]
        self.Findex=[np.flatnonzero(j==c) for c in range(n)]
        self.parent=expit(-e-r)
        self.neq=self.w@self.parent
    def split(self,y):
        return np.split(y[:-1],4)
    def admissible(self,y):
        nb,nf,kb,kf=self.split(y)
        if np.any(nb<0) or np.any(nf<0) or np.any(nf>self.F.number):
            return False
        if np.any(kb<self.B.lo*nb) or np.any(kb>self.B.hi*nb):
            return False
        for l,t in enumerate(nf/self.F.number):
            lo,hi=fd_bounds(self.q[l],self.q[l+1],self.b,t)
            if not lo*self.gf<=kf[l]<=hi*self.gf:
                return False
        return True
    def initial(self,perturbed):
        B,F=self.B,self.F
        phi=1/np.expm1(B.k+self.m)
        if perturbed:
            phi=phi*(1+.3*np.sin(np.arange(self.n)))[:,None]*.3
        fd=expit(-F.k-self.b)
        if perturbed:
            fd*=.73
        def moments(mesh,f):
            return ((mesh.w*f).sum(1)*mesh.number,
                    (mesh.w*f*mesh.k).sum(1)*mesh.number)
        nb,kb=moments(B,phi)
        nf,kf=moments(F,fd)
        y=np.r_[nb,nf,kb,kf,0.]
        assert self.admissible(y)
        return y
    def reconstruct(self,y):
        n=self.n
        nb,nf,kb,kf=self.split(y)
        p=np.zeros_like(self.e)
        f=np.zeros_like(self.r)
        thermo=np.zeros((4,n))
        for c in range(n):
            a,z=self.q[c:c+2]
            ib,iff=self.Bindex[c],self.Findex[c]
            if nb[c] and ib.size:
                target=(kb[c]/nb[c]-self.B.lo[c])/(self.B.hi[c]-self.B.lo[c])
                if self.thermal=="regular":
                    A,C=infer_regular_be(a,z,self.m,nb[c]/self.B.number[c],target,order=96)
                    slope=C/(self.B.hi[c]-self.B.lo[c])
                    thermo[0,c]=A-slope*self.B.lo[c]
                    thermo[2,c]=slope
                    s=(self.eb[ib]-self.B.lo[c])/(self.B.hi[c]-self.B.lo[c])
                    p[ib]=occupation(A,C,s)
                else:
                    shape,beta=infer_be(a,z,self.m,nb[c]/self.B.number[c],target,self.thermal,80)
                    p[ib]=shape(self.e[ib],self.eb[ib])
            if nf[c] and iff.size:
                if nf[c]>=self.F.number[c]:
                    f[iff]=1.
                else:
                    target=(kf[c]/nf[c]-self.F.lo[c])/(self.F.hi[c]-self.F.lo[c])
                    A,B=infer_fd(a,z,self.b,nf[c]/self.F.number[c],target,order=80)
                    s=(self.ef[iff]-self.F.lo[c])/(self.F.hi[c]-self.F.lo[c])
                    f[iff]=expit(A-B*s)
                    slope=B/(self.F.hi[c]-self.F.lo[c])
                    thermo[1,c]=-A-slope*self.F.lo[c]
                    thermo[3,c]=slope
        return p,f,thermo

    def flux(self,p,f):
        w,J1=self.w,self.parent
        G=w*J1*(1-f)*(1+p)
        L=w*(1-J1)*f*p
        return G,L

    def project(self,J):
        n=self.n
        return np.r_[np.bincount(self.i,weights=J,minlength=n),
                     np.bincount(self.j,weights=J,minlength=n),
                     np.bincount(self.i,weights=J*self.eb,minlength=n),
                     np.bincount(self.j,weights=J*self.ef,minlength=n),
                     -(self.e+self.r)@J]
    def rhs(self,y):
        p,f,_=self.reconstruct(y)
        G,L=self.flux(p,f)
        return self.project(G-L)

    def entropy_check(self,y):
        p,f,thermo=self.reconstruct(y)
        g,l=self.flux(p,f)
        v=g-l
        e=float(v@(np.log(g)-np.log(l)))
        ydot=self.project(v)
        moment=float(thermo.ravel()@ydot[:-1]+ydot[-1])
        err=abs(e-moment)/max(abs(e),1e-300)
        assert e>0 and err<2e-9
        return dict(entropy_production_event=e,entropy_production_moment=moment,
                    relative_identity_error=err)

    def invariants(self,y):
        nb,nf,kb,kf=self.split(y)
        return np.array([(nb-nf).sum(),
                         self.m*nb.sum()+self.b*nf.sum()+kb.sum()+kf.sum()+y[-1]])


def solve(n,M,order,thermal,steps,dtot=1e-3):
    m=Model(n,M,order,thermal)
    eq=m.initial(False)
    eql=np.sum(abs(m.rhs(eq)[:n]))/m.neq
    if thermal:
        assert eql<2e-10
    y=m.initial(True)
    y0=y.copy()
    t,reject=0.,0
    hmax=dtot/steps
    for _ in range(8*steps):
        if dtot-t<=8*np.spacing(dtot):
            break
        dy=m.rhs(y)
        h=min(hmax,dtot-t)
        while True:
            v=y+h*dy
            if m.admissible(v):
                z=(y+v+h*m.rhs(v))/2
                if m.admissible(z):
                    break
            h*=.5
            reject+=1
            assert reject<300 and h>0 and t+h>t
        y=z
        t+=h
    assert abs(t-dtot)<1e-13
    n0,f0,k0,l0=m.split(y0)
    nb,nf,kb,kf=m.split(y)
    d=m.invariants(y)-m.invariants(y0)
    assert np.max(abs(d))<3e-13
    return dict(n=n,M=M,order=order,hmax_steps=steps,thermal_prior=thermal,t_end=t,
                regular_BE_entropy_identity=m.entropy_check(y) if thermal=="regular" else None,
                all_rejections=reject,equilibrium_scalar_source_l1_over_vacuum=eql,
                net_boson=float((nb-n0).sum()),net_boson_kinetic_energy=float((kb-k0).sum()),
                boson_N_l1_change=float(abs(nb-n0).sum()),parent_energy_debit=float(y[-1]),
                max_conservation_drift=float(np.max(abs(d))),
                final_state=y.tolist())

def run():
    rows=[solve(n,.1,order,thermal,s) for n in (12,24)
          for order,s in ((4,6),(4,12),(6,12))
          for thermal in (True,False)]
    rows += [solve(n,M,6,"regular",12) for M,n in ((.1,12),(.1,24),(.005,12))]
    for r in rows: assert r["boson_N_l1_change"]>0.
    return dict(rows=rows,status="conditional fixed-bath evolution; not FRW or a physical relic")
if __name__=="__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
