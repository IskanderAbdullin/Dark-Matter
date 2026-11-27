import argparse
import json
import numpy as np
from scipy.special import expit
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_archived_cold_screened_elastic_reference import Screened

def one(grid,mesh,level,projection,cells=False,subdivision=1,shape_amplitude=0.,shape_degree=3):
    quad=[(1,1,2),(1,2,3)][level]
    M=Screened(grid,1,8,8e-8,quad,.08,1.,mesh,projection,subdivision,shape_amplitude,shape_degree)
    N=M.M
    e,f,G,_,_=N.weights(N.init,True)
    z=np.r_[N.init.ravel(),M.ns,0.,0.]
    select=slice(M.j0,M.j0+M.n)
    source=lambda x: Cold.rhs(M,x,False)[select]
    C=source(z)
    x=M.collision(e,f,G,split_affinity=True)
    rows={}
    for k in ('total','heat','shape'):
        ray=np.r_[x[k].ravel(),np.zeros(M.n+2)]
        values=[]
        for h in (1e-14,5e-14,2e-13):
            y=(source(z+h*ray)-C)/h
            values.append(y.tolist())
        rows[k]=values
    # Independent symmetric differentiation in Fermi logit coordinates.
    # An ordinary backward displacement of moments can leave the Pauli cone,
    # although there is a perfectly regular infinitesimal tangent.
    base=N.r[1].start
    occ0=[expit(f[j]@B.B) for j,B in enumerate(N.cells[base:],start=base)]
    chi=[(c.I*(o*(1-o)))@c.B.T for c,o in zip(N.cells[base:],occ0,strict=True)]
    y0=N.init
    projected={'total':x['total'],'heat':x['heat'],'shape':x['shape']}
    for label,sl in zip(('shape_parent','shape_daughter'),N.r[1:],strict=True):
        J=np.zeros_like(x['shape'])
        J[sl]=x['shape'][sl]
        projected[label]=J
    centered={}
    for key,ray in projected.items():
        coeff=np.array([np.linalg.solve(a,b) for a,b in zip(chi,ray[base:],strict=True)])
        cmax=max(np.max(abs(coeff)),1.)
        qrows=[]
        for size in (2e-5,1e-4,5e-4):
            h=size/cmax
            def only_f(t):
                R=y0.copy()
                for j,(B,o,v) in enumerate(zip(N.cells[base:],occ0,coeff,strict=True),start=base):
                    R[j]+=B.I@(expit((f[j]+t*v)@B.B)-o)
                return source(np.r_[R.ravel(),M.ns,0.,0.])
            qrows.append(((only_f(h)-only_f(-h))/(2*h)).tolist())
        centered[key]=qrows
    discrepancy=np.max(abs(np.asarray(centered['total'])-rows['total']))/max(np.max(np.abs(rows['total'])),1.)
    assert discrepancy<1e-2, (grid,mesh,level,discrepancy)
    cell_source=[]
    if cells:
        for i,(cell,O,cov) in enumerate(zip(N.cells[base:],occ0,chi,strict=True),start=base):
            coeff=np.linalg.solve(cov,x['shape'][i])
            h=1e-4/max(np.max(abs(coeff)),1.)
            def perturb(sign):
                Y=y0.copy()
                Y[i]+=cell.I@(expit((f[i]+sign*h*coeff)@cell.B)-O)
                return source(np.r_[Y.ravel(),M.ns,0.,0.])
            by_q=(perturb(1.)-perturb(-1.))/(2*h)
            cell_source.append(dict(flavor='parent' if i<N.r[2].start else 'daughter',
                                    momentum_interval=[cell.lo,cell.hi],
                                    scalar_source_derivative=by_q.tolist()))
        s1=np.sum([j['scalar_source_derivative'] for j in cell_source],axis=0)
        s2=np.asarray(centered['shape'][1])
        closure=np.max(abs(s2-s1))/max(np.max(abs(s2)),1.)
        assert closure<2e-2,(grid,mesh,level,closure)
    target=f.copy()
    elasticity=[]
    for delta in (.01,.3,1.2):
        for r,b,mu in zip(N.r[1:],(1.,1.+delta),(.1,-.4),strict=True):
            target[r]=-b*e[r]
            target[r,0]+=mu
        y=M.collision(e,target,G,split_affinity=True)
        assert abs(y['heat_entropy_fraction']-1.)<2e-10
        assert abs(y['heat_flux']/delta/y['conductance']-1.)<2e-10
        elasticity.append([delta,y['heat_flux'],y['conductance'],y['heat_entropy_fraction']])
    rel=np.max(np.abs(np.asarray(rows['total'])-rows['heat']-rows['shape']))/max(np.max(abs(np.asarray(rows['total']))),1.)
    assert rel<1e-3
    return dict(grid=grid,mesh=mesh,quad=quad,projection=projection,subdivision=subdivision,shape_amplitude=shape_amplitude,shape_degree=shape_degree,archive_scalar_source=C.tolist(),
                archived={k:v for k,v in x.items() if k not in rows},
                frozen_two_temperature=elasticity,source_derivatives=rows,centered_logit_derivatives=centered,
                coordinate_derivative_relative_difference=float(discrepancy),relative_tangent_closure=rel,
                heat_free_cell_tangent=cell_source)

if __name__=='__main__':
    p=argparse.ArgumentParser(description='Separate daughter heat from its zero-heat elastic residual locally.')
    p.add_argument('--level',type=int,choices=(0,1),default=0)
    p.add_argument("--projection",choices=("steps","cell-conservative"),default="steps")
    p.add_argument("--cells",action="store_true",help="attribute the heat-free source derivative to individual cells")
    p.add_argument("--subdivision",type=int,choices=range(1,7),default=1)
    p.add_argument("--grid",choices=("fine","third"))
    p.add_argument("--mesh",choices=("fine","third"))
    p.add_argument("--shape-amplitude",type=float,default=0.)
    p.add_argument("--shape-degree",type=int,choices=(3,4,5),default=3)
    a=p.parse_args()
    print(json.dumps([one(g,m,a.level,a.projection,a.cells,a.subdivision,a.shape_amplitude,a.shape_degree)
                      for g in ((a.grid,) if a.grid else ('fine','third'))
                      for m in ((a.mesh,) if a.mesh else ('fine','third'))],indent=2,allow_nan=False))
