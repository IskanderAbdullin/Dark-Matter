"""Exact stationary neutral correlator for two species sharing a toy bath."""
import json
import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp


def one(n,temp,shift,bath,pbias,rotate):
    # N_total=1 sector with distinct conserved species labels. The orbital
    # two-level and oscillator degrees of freedom are identical in dimension.
    a=np.diag(np.sqrt(np.arange(1.,n)),1)
    o=.43*np.diag(np.arange(n))
    I,D,X=np.eye(n),np.diag([1.,-1.]),np.array([[0.,1.],[1.,0.]])
    H=lambda orb,g,e: np.kron(np.eye(2),o+e*I)+.5*orb*np.kron(D,I)+g*np.kron(X,a+a.T)
    ep,P=eigh(H(.31+shift,bath,.73))
    eb,B=eigh(H(.31-shift,bath,0.))
    M=B.T@np.kron(np.cos(rotate)*np.eye(2)+1j*np.sin(rotate)*X,I)@P
    m=np.abs(M)**2
    E=ep[None,:]-eb[:,None]
    # Chemical potential difference (parent minus daughter).
    lp=(pbias-ep)/temp
    lb=-eb/temp
    Z=logsumexp(np.r_[lp,lb])
    emit=m*np.exp(lp-Z)[None,:]
    reverse=m*np.exp(lb-Z)[:,None]
    mask=(emit>1e-16)&(reverse>1e-16)
    bath_defect=np.max(abs(np.log(reverse[mask])-np.log(emit[mask])-(E[mask]-pbias)/temp))
    assert bath_defect<2e-13
    assert abs(emit.sum()/np.exp(lp-Z).sum()-1)<2e-13
    assert abs(reverse.sum()/np.exp(lb-Z).sum()-1)<2e-13
    # Compare directly with the KMS ratio in both stimulated directions.
    err=np.max(abs(reverse-emit*np.exp((E-pbias)/temp)))
    assert err<2e-13
    center=(E*emit).sum()/emit.sum()
    std=np.sqrt((emit*(E-center)**2).sum()/emit.sum())
    off=emit[abs(E-.73)>.02].sum()/emit.sum()
    reverse_off=reverse[abs(E-.73)>.02].sum()/reverse.sum()
    # Convolve with one positive Gaussian in frequency, first naively,
    # then by broadening only the symmetric spectral core in the thermal split.
    omega=pbias+np.array([-.015,0.,.015])
    sigma=.035
    smooth=[]
    for x in omega:
        k=np.exp(-.5*((E-x)/sigma)**2)/(sigma*np.sqrt(2*np.pi))
        e0=float(np.sum(k*emit)); b0=float(np.sum(k*reverse))
        e1=float(np.sum(k*np.sqrt(emit*reverse))*np.exp(-(x-pbias)/(2*temp)))
        b1=e1*np.exp((x-pbias)/temp)
        smooth.append(dict(omega=x,naive_emit_minus_inverse=e0-b0,
                           kms_emit_minus_inverse=e1-b1))
    assert smooth[0]['kms_emit_minus_inverse']>0
    assert smooth[1]['kms_emit_minus_inverse']==0.
    assert smooth[2]['kms_emit_minus_inverse']<0
    # Spectral weights are nonnegative by their Lehmann construction.
    return dict(n=n,T=temp,orbital_halfshift=shift,bath=bath,mu_gap=pbias,
                probe_rotation=rotate,mean_emitted_energy=center,spread_emitted_energy=std,
                emit_weight_off_nominal=off,reverse_weight_off_nominal=reverse_off,coarse_smearing=smooth,
                max_log_kms_defect=bath_defect,
                max_unitary_sum_rule_defect=float(max(np.max(abs(m.sum(axis=0)-1)),np.max(abs(m.sum(axis=1)-1)))))


def run():
    rows=[]
    for n in (20,36,52):
        for shift,rot in [(0.,0.),(.035,0.),(.14,0.),(0.,.3),(.035,.3)]:
            r=one(n,.8,shift,.12,.68,rot)
            rows.append(r)
            if shift==rot==0:
                assert r['spread_emitted_energy']<1e-10 and r['emit_weight_off_nominal']<1e-16
    for a,b in zip(rows[5:10],rows[10:]):
        assert abs(a['spread_emitted_energy']-b['spread_emitted_energy'])<1e-7
        assert abs(a['emit_weight_off_nominal']-b['emit_weight_off_nominal'])<1e-7
    return dict(scope='exact single-particle two-flavor oscillator witness; not an Abelian momentum-space rate or cosmology',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
