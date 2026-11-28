"""Fixed reported comoving-momentum nodes for forward-only depletion."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from archived_scalar_pair_reaction_reference import amplitude
from conditional_relic_normalization_reference import conversion
from finite_cooling_reheat_distribution_reference import deposit
from finite_cooling_reheat_portal_reference import SourceFold,TF
from three_species_decay_onsager_reference import m,mp

def matrix(p,n):
    eps=np.hypot(m,p)
    V=p*p/(eps+m)
    E,F=eps[:,None,None],eps[None,:,None]
    pp=p[:,None,None]*p[None,:,None]
    a,g=roots_legendre(n)
    S=4*m*m+2*(m*(V[:,None,None]+V[None,:,None])+V[:,None,None]*V[None,:,None]-pp*a)
    K=np.zeros((2,len(p),len(p)))
    pf=np.sqrt(np.maximum(.25-mp*mp/S,0.))*g/2
    boost2=np.maximum(0.,p[:,None,None]**2+p[None,:,None]**2+2*pp*a)/S
    kin=boost2*np.maximum(0.,S/4-mp*mp)
    ci2=np.clip(((E-F)/2)**2/np.maximum(1e-300,boost2*np.maximum(0.,S/4-m*m)),0.,1.)
    rest=((E+F)/2)**2-mp*mp
    for x,gx in zip(a,g):
        weight=gx*amplitude(S,x)*pf
        K[0]+=weight.sum(axis=-1)
        K[1]+=(weight*(rest+kin*(ci2*x*x+(1-ci2)*(1-x*x)/2))).sum(axis=-1)
    assert np.min(K)>=0.
    return K/(64*np.pi*eps[:,None]*eps[None,:])


def matrix_check():
    from finite_cooling_pair_emission_reference import kernel_and_lab_p2
    p=np.r_[0.,np.geomspace(.00001,.25,45)]
    w=np.exp(-p/.009)
    w/=w.sum()
    M=matrix(p,8)
    k=kernel_and_lab_p2(p,w,8)[0]
    q=w@M@w
    err=np.max(abs(q/k-1.))
    assert err<2e-13 and np.max(abs(M-np.transpose(M,(0,2,1))))<1e-25
    return float(err)


def exact_substep_optical_bracket(K,lo,frozen):
    # Linear no-source lower trajectory obtained by holding both partner
    # densities at the full passive upper solution. Integrate its
    # positive bilinear emission exactly over this split substep.
    depth=K[0]@frozen
    A=depth[:,None]+depth[None,:]
    factor=np.divide(-np.expm1(-A),A,out=np.ones_like(A),where=A>0.)
    weighted=lo[:,None]*lo[None,:]*factor
    low=(K*weighted).sum(axis=(1,2))
    high=frozen@K@frozen
    return low,high,depth


def run(eos,reference,high=False):
    config=json.loads(reference.read_text())
    eh=hashlib.sha256(eos.read_bytes()).hexdigest()
    assert config["eos_sha256"]==eh
    B=SourceFold(eos,(7,18,14),100,4)
    C=conversion(eos)['Omega_h2_per_GeV_times_yield']
    z=np.log(TF)
    sf=B.s(z)
    out=[]
    for r0 in (config['controls'][-1]['rows'][0], config['controls'][-1]['rows'][-1]):
        r,h=np.log([r0['T_R_GeV'],r0['Tmax_GeV']])
        ar=(sf/B.s(r))**(1/3)
        for nt,ng,nang in (((96,95,10),) if high else ((36,48,6),(68,64,8))):
            p=np.r_[0.,np.geomspace(.001,130,ng)]
            edge=np.r_[p,np.inf]
            frozen=np.zeros(len(p))
            actual=frozen.copy()
            Bopt=frozen.copy()
            lower_state=frozen.copy()
            up=np.zeros((2,2)); lower=np.zeros((2,2)); net=np.zeros((2,2))
            source_p2=0.
            gx,gw=roots_legendre(3)
            for epoch,(hi,lo) in enumerate(((h,r),(r,z))):
                if hi<=lo:continue
                def get(u):
                    if epoch==0:
                        a=ar*np.exp(-8*(u-r)/3)
                        H=B.H(r)*np.exp(4*(u-r))
                        return a,(8/3)/H
                    a=(sf/B.s(u))**(1/3)
                    return a,(1+B.grad(u)/3)/B.H(u)
                def inject(l,u):
                    ans=np.zeros((2,len(p)))
                    for v,w in zip((l+u+(u-l)*gx)/2,(u-l)*gw/2):
                        a,dt=get(v)
                        ans+=deposit(B,np.exp(v),edge,a,20)[1]*(a**3*dt*w/sf)
                    N,P2=ans
                    assert N[-1]<1e-12*N.sum()
                    f=np.clip((P2[:-1]-N[:-1]*p[:-1]**2)/np.diff(p*p),0.,N[:-1])
                    S=N.copy()
                    S[:-1]-=f
                    S[1:]+=f
                    assert abs((S@(p*p))/P2.sum()-1.)<2e-8
                    return S,float(P2.sum())
                u=np.linspace(hi,lo,nt+1)
                for t1,t0 in zip(u[:-1],u[1:]):
                    tm=(t0+t1)/2
                    s,sp=inject(tm,t1)
                    frozen+=s;actual+=s;lower_state+=s;source_p2+=sp
                    a,dt=get(tm)
                    K=matrix(p*(TF/a),nang)*(dt*(t1-t0)*sf/(2*a**3))
                    K[1]*=(a/TF)**2
                    lb,ub,step=exact_substep_optical_bracket(K,lower_state,frozen)
                    Bopt+=step
                    lower_state*=np.exp(-step)
                    up[epoch]+=ub
                    lower[epoch]+=lb
                    M=K[0]
                    mid=actual*np.exp(-(M@actual)/2)
                    d=-actual*np.expm1(-(M@mid))
                    actual-=d
                    ww=mid@K@mid
                    net[epoch]+=d.sum()*ww/ww[0]
                    assert np.all(actual<=frozen*(1.+1e-14))
                    assert np.all(actual>=lower_state*(1.-1e-14))
                    s,sp=inject(t0,tm)
                    frozen+=s;actual+=s;lower_state+=s;source_p2+=sp
            assert abs(net[:,0].sum()/(frozen.sum()-actual.sum())-1.)<2e-9
            val=dict(T_R=r0['T_R_GeV'],T_max=r0['Tmax_GeV'],panels=nt,nodes=len(p),
                maximum_accumulated_node_optical_depth=float(Bopt.max()),
                spontaneous_frozen_omega_scalar=(C*mp*up[:,0]).tolist(),
                feedback_midpoint_converted_omega_scalar=(C*mp*net[:,0]).tolist(),
                absorbing_linear_lower_omega_scalar=(C*mp*lower[:,0]).tolist(),
                emitted_scalar_p_over_TF2_moment_yield_lower_feedback_upper=np.stack([lower[:,1],net[:,1],up[:,1]]).tolist(),
                emitted_scalar_rms_p_over_TF_feedback=np.sqrt(np.divide(net[:,1],net[:,0],out=np.zeros(2),where=net[:,0]>0.)).tolist(),
                conservative_global_suppression_bound_for_split_operator=float(-np.expm1(-2*Bopt.max())),
                stable_omega_before_after=(C*m*np.array([frozen.sum(),actual.sum()])).tolist(),
                p2_projection_defect=float((frozen@(p*p))/source_p2-1.))
            jj=up>0.
            assert np.all(lower[jj]<=net[jj]*(1+2e-9)) and np.all(net[jj]<=up[jj]*(1+2e-9))
            out.append(val)
    return dict(eos_sha256=eh,independent_matrix_mean_relative_defect=matrix_check(),reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),cases=out,
      scope="Forward-only positive depletion on fixed comoving nodes with externally injected prompt production. Exact analytic bounds refer to the split fixed-kernel dynamics, not integration with continuous thermal source between panels or a continuum physical bound. Feedback column uses a positive exponential midpoint approximation. Omits inverse, Bose response, scalar feedback and unrepresented physics.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos",type=Path)
    p.add_argument("reference",type=Path)
    p.add_argument("--high",action="store_true")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.reference,a.high),indent=2,allow_nan=False))
