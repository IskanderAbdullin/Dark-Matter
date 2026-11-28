"""Two real crossed processes: number and photon-energy Onsager matrix at common equilibrium."""
import gzip,json
from pathlib import Path
import numpy as np
from finite_mass_qed_crossed_annihilation_source_reference import value,Q
from finite_mass_qed_scalar_absorption_reference import audit


def matrix(data,T,mesh):
    pair=[]; abs0=[]
    value(data,T,68.,mesh[0],0,mesh[1],mesh[2],equilibrium=True,thermal_measure=pair)
    audit(data,T,56.,*mesh,equilibrium=True,thermal_measure=abs0)
    Es=np.hypot(Q,data['checkpoints'][-1]['scalar_pole_over_T'])
    results=[]
    b1=1/T
    for name,sign,data0 in [('pair',1,pair),('scalar_absorption',-1,abs0)]:
        k,A,B=np.concatenate(data0,axis=0).T
        q=np.vstack([np.ones_like(k),sign*k])
        def stable_nb(x):return np.exp(-x)/(-np.expm1(-x))
        nb0=stable_nb(b1*k); fs0=stable_nb(b1*Es)
        G0=A*(1+nb0 if sign==1 else nb0)*(1+fs0)
        L0=B*(nb0 if sign==1 else 1+nb0)*fs0
        weights=(G0+L0)/2
        C=(q*weights)@q.T
        def rhs(X,Y):
            fs=stable_nb(b1*Es+X)
            nb=stable_nb((b1+Y)*k)
            G=A*(1+nb if sign==1 else nb)*(1+fs)
            L=B*(nb if sign==1 else 1+nb)*fs
            return q@(G-L)
        eps=1e-7*Es/T
        D=np.column_stack([(rhs(eps,0.)-rhs(-eps,0.))/(2*eps),
                           (rhs(0.,1e-5)-rhs(0.,-1e-5))/(2e-5)])
        err=np.max(abs(D-C)/np.maximum(abs(C),1e-30))
        thermo=abs(rhs(0.,0.)).max()/max(C[0,0],1e-30)
        if err>2e-5 or thermo>3e-13 or min(np.linalg.eigvalsh(C))<=0:
            raise ArithmeticError((name,T,err,thermo,C,D))
        results.append(dict(channel=name,exchange_matrix=C.tolist(),independent_fd_relative_defect=float(err),
                            equilibrium_defect=float(thermo)))
    C=sum(np.array(x['exchange_matrix']) for x in results)
    return dict(channels=results,combined_matrix=C.tolist(),
                determinant=float(np.linalg.det(C)),number_heat_correlation=float(C[0,1]/np.sqrt(C[0,0]*C[1,1])))


def run(grids):
    rows=[]
    for g in grids:
       with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
       for T in (.4,1.):
           a,b=(matrix(d,T,mesh) for mesh in ((4,4,12),(7,7,22)))
           w=np.array(a['combined_matrix']);t=np.array(b['combined_matrix'])
           rows.append(dict(grid=g,T=T,refined=b,maximum_matrix_mesh_relative_difference=float(np.max(abs(w/t-1.)))))
    return dict(scope='FD/BE equilibrium matrix of two separate infrared-finite real processes; held scalar q and masses, no chemical or thermal history.',rows=rows)


if __name__=='__main__':
    import sys
    print(json.dumps(run(tuple(sys.argv[1:]) or ('fine','third')),indent=2,allow_nan=False))
