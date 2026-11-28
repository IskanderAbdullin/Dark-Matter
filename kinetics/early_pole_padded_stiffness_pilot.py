"""Fixed-chart derivative and stiffness controls for thin charged face cells."""
import json

import numpy as np
from scipy.sparse import coo_matrix

from early_pole_double_charge_stencil_pilot import chart,event
from early_pole_edge_padded_charged_projection_pilot import charge_mesh
from early_pole_endpoint_stencil_remap_pilot import LABELS,SPECS,archives,load,parts,remap,MassiveFermion
from early_pole_live_parent_energy_reference import state_moments


def derivative(ch,d,y):
    es,ed,ep,i,j,k,a,b,W=ch
    fs,fd,fp=np.split(y,(d.ns,d.iP))
    pg,pc,pu,pv,cu,cv=d.parent_products(fp[k-1],fp[k],1-b,b,True)
    dg,dc,du,dv,bu,bv=d.parent_products(fd[j-1],fd[j],1-a,a,True)
    u=(1+fs[i]);v=fs[i]
    A=(i,d.ns+j-1,d.ns+j,d.iP+k-1,d.iP+k)
    weights=(np.ones_like(a),1-a,a,-(1-b),-b)
    partials=(
        W*(pg*dc-pc*dg),
        W*(pg*u*bu-pc*v*du),
        W*(pg*u*bv-pc*v*dv),
        W*(pu*u*dc-cu*v*dg),
        W*(pv*u*dc-cv*v*dg))
    rows=[];cols=[];entries=[]
    for row,w in zip(A,weights,strict=True):
        for column,part in zip(A,partials,strict=True):
            rows.append(row)
            cols.append(column)
            entries.append(w*part/d.w[row])
    J=coo_matrix((np.concatenate(entries),(np.concatenate(rows),np.concatenate(cols))),
                 shape=(d.n,d.n)).tocsc()
    return J


def run():
    output=[]
    for g,shape in SPECS.items():
        arc=load(g)[0]
        x=arc['x_end']
        data=[arc['endpoint_'+lab+'_occupations'] for lab in LABELS]
        base=archives(arc)
        cases=[]
        for f in (0.,1/8,1/64,1/512,1/4096):
            d=MassiveFermion(*shape,arc['mass_ratio'],arc['elastic_rate'],True)
            charge_mesh(d,f,f)
            y=np.concatenate([remap(a,h,b)[0] for a,h,b in zip(base,data,parts(d),strict=True)])
            X,ra,rb,m=state_moments(d,x,y)
            H=np.sqrt(d.radiation/(d.radiation+ra+rb))
            ch,_=chart(d,x,m,4)
            J=derivative(ch,d,y)*H
            v=np.sin(.31+np.arange(len(y)))*y
            dz=1e-5
            estimate=H*(event(ch,y+dz*v,d)[0]-event(ch,y-dz*v,d)[0])/(2*dz*d.w)
            exact=J@v
            dev=float(np.max(abs(estimate-exact))/np.max(abs(exact)))
            if dev>5e-7:
                raise ArithmeticError(('Fixed-chart analytic derivative not reproduced',g,f,dev))
            basis=np.zeros((3,d.n))
            basis[0,:d.ns]=1.
            basis[0,d.ns:d.iP]=-1.
            basis[1,:d.ns]=1.
            basis[1,d.iP:]=1.
            basis[2]=X
            shadow=d.w[None,:]*basis
            null_defect=float(np.max(abs(J.T@shadow.T))/max(np.max(abs(J).T@abs(shadow).T),1e-300))
            if null_defect>8e-13:
                raise ArithmeticError('Frozen stencil derivative changed an invariant')
            diag=J.diagonal()
            idx=d.ns+int(np.argmin(diag[d.ns:]))
            lo,hi=(d.d.fed[idx-d.ns:idx-d.ns+2] if idx<d.iP else d.pe[idx-d.iP:idx-d.iP+2])
            cases.append(dict(eta=f,dimension=d.n,jacobian_nonzeros=J.nnz,
                             largest_daughter_or_parent_diagonal_magnitude=float(abs(diag[idx])),
                             inverse_diagonal_scale_log_a=float(1/abs(diag[idx])),
                             smallest_scale_cell_species='daughter' if idx<d.iP else 'parent',
                             smallest_scale_cell_edges=[lo,hi],finite_difference_relative_defect=dev,
                             relative_weighted_jacobian_left_null_defect=null_defect))
        output.append(dict(grid=g,rows=cases))
    return output


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
