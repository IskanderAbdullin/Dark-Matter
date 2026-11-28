"""Independent assembled-vs-directional benchmark on a wide physical step."""
import argparse,json,time
from pathlib import Path
import numpy as np
from compact_prompt_jacobian_reference import compile_prompt,local_prompt
from continuum_dark_normal_capacity_reference import from_row
from joint_resonant_centered_stiff_reference import Centered
from joint_resonant_frw_split_reference import pack,redshift
from joint_prompt_two_moment_evolution_reference import Collision
from late_archived_shape_continuation_reference import Late
from late_radiation_background_reference import TabulatedRadiation
from late_soft_two_moment_test_reference import totals
from sm_dark_entropy_exchange_reference import StandardModel,MPL


def test(eos,table,row,checkpoint):
    cap=from_row(row,eos,20)
    bg=TabulatedRadiation(table,cap["actual_visible_T_GeV"])
    b0=checkpoint["b"]
    b=b0*np.exp(np.log(1.5)/8)
    q=np.array(checkpoint["mesh_q"])
    th=[np.array(x) for x in checkpoint["theta"]]
    th,*_=redshift(q,b0,b,th)
    k=Collision(q/b)
    state=pack(th)
    m0=Late(k,1.,[],False)
    U=m0.moment(state)[0]
    rho=totals(k,U,b)[1].sum()/b**4
    H=np.sqrt((bg.state(b)[1]+rho)/3)/MPL
    op=Late(k,H,[],True)
    ref,R,T=op.reference(state,cap["common_T_GeV"]/b)
    z=state-ref
    ts=time.perf_counter()
    panels=compile_prompt(op)
    setup_s=time.perf_counter()-ts
    dirs=np.eye(len(state))
    def diff(z):
        ts=time.perf_counter()
        F,a=Centered.centered(op,z,ref,R)
        X=np.stack([a(p) for p in dirs],axis=1)
        dense=time.perf_counter()-ts
        ts=time.perf_counter()
        G,b=local_prompt(op,z,ref,R,panels)
        Y=np.stack([b(p) for p in dirs],axis=1)
        local=time.perf_counter()-ts
        derr=np.linalg.norm(X-Y)/np.linalg.norm(X)
        assert np.array_equal(F,G)
        assert derr<2e-12
        v=np.cos(np.arange(len(z))*.23)
        h=4e-6
        pp=local_prompt(op,z+h*v,ref,R,panels)[0]
        qq=local_prompt(op,z-h*v,ref,R,panels)[0]
        dv=(pp-qq)/(2*h)
        j=Y@v
        fd=float(np.linalg.norm(j-dv)/np.linalg.norm(j))
        assert fd<8e-7
        n=len(op.k.c[0].vol)
        V=np.zeros((len(z),3))
        for j,c in enumerate(op.k.c):
            rr=slice(2*j*n,2*(j+1)*n)
            V[rr,2]=np.stack([c.mass+c.low,c.dk],axis=1).ravel()
            V[rr.start:rr.stop:2,j]=1.
        V/=np.linalg.norm(V,axis=0)
        inv=np.linalg.norm(V.T@Y)/np.linalg.norm(Y)
        assert inv<3e-13
        return dict(J_relative_defect=derr,JVP_finite_difference_error=fd,invariant_null_error=inv,
                    directional_assembly_seconds=dense,parent_local_assembly_seconds=local)
    r=diff(z)
    s=diff(z*1e-6)
    return dict(b0=b0,midpoint_b=b,reference_T_GeV=T,prompt_nodes=len(k.y),
                parent_groups=len(panels),mean_active_columns=float(np.mean([len(x[1]) for x in panels])),
                global_columns=len(z),compilation_seconds=setup_s,archived_state=r,near_reference=s)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("table");p.add_argument("run");p.add_argument("checkpoint")
    a=p.parse_args()
    row=json.loads(Path(a.run).read_text())["rows"][-1]
    c=json.loads(Path(a.checkpoint).read_text())
    print(json.dumps(test(StandardModel(a.eos),a.table,row,c),indent=2,allow_nan=False))
