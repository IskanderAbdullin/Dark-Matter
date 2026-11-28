"""Number-conserving on-shell inverse decay followed by prompt scalar return."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre,expit
from tagged_scalar_nonthermal_bath_reference import Bath
from decay_scalar_ir_relaxation_reference import limits
from three_species_decay_onsager_reference import M,mp,P,Es,width
from sm_dark_entropy_exchange_reference import StandardModel,MPL

def kernel(p,bath,edges,n=8,blocked=True,occupation=None):
    """Outward rates from a single physical scalar momentum; last bin external."""
    e,a,d=limits(p)
    sz=len(edges)
    if p==0.:return np.zeros(sz),0.
    A=max(a,bath.q[0]+e,M)
    B=min(a+d,bath.q[-1]+e)
    if not A<B:return np.zeros(sz),0.
    occupation=bath if occupation is None else occupation
    location,wei=roots_legendre(n)
    ref=np.asarray(edges)
    ee=np.hypot(edges,mp)
    lefts,rightss=[],[]
    for j in ref:
        _,L,dL=limits(j)
        lefts.append(L);rightss.append(L+dL)
    # Parent-energy panels include return-bin openings plus broad coverage
    # over the injected partner support, including its low-energy turn-on.
    edges1=np.unique(np.r_[A,B,np.clip(lefts,A,B),np.clip(rightss,A,B),
          np.clip(e+np.geomspace(bath.q[0],bath.q[-1],24),A,B)])
    l=np.log(edges1[:-1])
    h=np.log(edges1[1:])
    es=np.exp((l[:,None]*(1-location)+h[:,None]*(1+location))/2).ravel()
    dE=((h-l)[:,None]*wei/2).ravel()*es
    fresh=dE*occupation(es-e)
    if fresh.sum()==0.:return np.zeros(sz),0.
    p1=np.sqrt(np.maximum(0.,(es-M)*(es+M)))
    upper=(es*Es+p1*P)/M
    lower=(Es*Es+(p1*mp/M)**2)/upper
    lo=np.maximum(lower[:,None],ee)
    hi=np.minimum(upper[:,None],np.r_[ee[1:],np.inf])
    v=np.maximum(0.,hi-lo)
    if blocked:
        x,w=roots_legendre(4)
        eF=lo[:,:,None]+v[:,:,None]*(1+x)/2
        flow=np.sum((1-occupation(es[:,None,None]-eF))*(w/2),axis=-1)*v
    else:
        flow=v
    prob=flow/flow.sum(axis=1)[:,None]
    scale=4*width*M*M/(2*P*p*e)
    Bins=scale*(fresh@prob)
    initial=scale*fresh.sum()
    assert np.all(Bins>=0.) and abs(Bins.sum()/initial-1)<2e-12
    return Bins,initial

def matrix(bath,edges,occupation=None,npin=10,n=8):
    v,w=roots_legendre(npin)
    cubed=((edges[:-1,None]**3*(1-v)+edges[1:,None]**3*(1+v))/2)
    q=np.cbrt(cubed)
    R=np.zeros((len(q),len(edges)))
    counts=np.zeros(len(q))
    for j,row in enumerate(q):
        for s,a in zip(row,w/2):
            g,rate=kernel(s,bath,edges,n,occupation=occupation)
            R[j]+=g*a
            counts[j]+=rate*a
    err=np.max(abs(R.sum(axis=1)-counts))/max(counts.max(),1e-280)
    assert err<3e-13
    return R,counts

def check(eos,data):
    out=[]
    for r in data["rows"]:
        if not (r["TR"] in (.02,.05,.1) and abs(r["TF"]/r["TR"]-.9)<1e-6):continue
        T=r["TF"]
        bath=Bath(r,eos)
        H=np.sqrt(eos.state(T)[0]/3)/MPL
        q=T*np.array([0.,.02,.05,.125,.3,.5,1.,2.,4.])
        for j in (.012,.03,.08,.2,.8):
            a,v=kernel(j*T,bath,q)
            b,u=kernel(j*T,bath,q,n=13)
            c,_=kernel(j*T,bath,q,n=13,blocked=False)
            if u/H<1e-4:continue
            refined=np.sum(abs(a-b))/u
            old=abs(bath.rate(j*T)/u-1)
            assert refined<.025 and old<.12,(j,T,refined,old)
            out.append(dict(TR=r["TR"],TF=T,input_p_over_T=j,rate_over_H=u/H,
                            absolute_distribution_refine=refined,total_rate_to_cumulative_relative_error=old,
                            Pauli_angle_effect_l1=float(np.sum(abs(c-b))/u),
                            probability_return_below_half_T=b[:5].sum()/u,
                            probability_outside_four_T=b[-1]/u,
                            width_dilation_floor_over_H=width*M/(bath.q[-1]+j*T)/H))
    return dict(rows=out,scope="Prompt 1↔2 sequential resonant scalar redistribution tagged limit, no coherent or stimulated scalar effects, frozen nonthermal fermions; do not add full resonant elastic cross-section in parallel.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("spectra")
    a=p.parse_args()
    print(json.dumps(check(StandardModel(a.eos),json.loads(Path(a.spectra).read_text())),indent=2,allow_nan=False))
