"""All parent energies: exact fixed-q signed real recoil, no virtual finite part."""
import gzip,json
from pathlib import Path
import functools
import numpy as np
from scipy.special import roots_legendre
from scipy.optimize import brentq
from early_pole_decay_mode_clock_reference import intervals,Q

base=Path(__file__).parent/'results'
@functools.lru_cache(maxsize=8)
def gl(n):return roots_legendre(n)


class BridgeSource:
    def __init__(self,d,emission_only_extension=False):
        self.emission_only_extension=bool(emission_only_extension)
        self.M=np.longdouble(str(d['x_end']));self.m=self.M*np.longdouble(str(d['mass_ratio']))
        self.s=np.longdouble(str(d['checkpoints'][-1]['scalar_pole_over_T'])); self.q=np.longdouble(str(Q))
        self.Es,self.lo,self.hi,*_=intervals(self.M,self.s,self.m,Q)
        self.pe=np.hypot(np.asarray(d['endpoint_parent_momentum_edges'],dtype=np.longdouble),self.M)
        self.be=np.hypot(np.asarray(d['endpoint_fermion_momentum_edges'],dtype=np.longdouble),self.m)
        self.pf=np.r_[0.,d['endpoint_parent_occupations'],0.]
        self.bf=np.r_[0.,d['endpoint_fermion_occupations'],0.]
        dy=np.unique(np.r_[self.pe,self.be+self.Es])
        xx=np.r_[self.lo,dy[(dy>self.lo)&(dy<self.hi)],self.hi]
        mid=(xx[1:]+xx[:-1])/2
        self.stim=np.diff(xx)@(self.pf[np.searchsorted(self.pe,mid)]-self.bf[np.searchsorted(self.be,mid-self.Es)])
        t=np.arccosh((self.M*self.M+self.m*self.m-self.s*self.s)/(2*self.M*self.m))
        self.B=2*(t/np.tanh(t)-1)

    def signed(self,k,n=30,o=8,az=20):
        k=np.longdouble(k)
        M,m,s,q,Es=self.M,self.m,self.s,self.q,self.Es
        cut=np.unique(np.r_[self.pe,self.be+Es+k,self.lo,self.hi])
        cut=cut[(cut>M)&(cut-k-Es>m)]
        # This bridge is limited to a sector in which all parent
        # aperture equations are strictly monotone, including translated
        # archive edges far outside the zero-photon aperture.
        if self.emission_only_extension:
            if k<0. or k>np.longdouble('2.561e-3'):
                raise ValueError('Outside verified monotone emission geometry')
        elif abs(k)>np.longdouble('2.01e-4'):
            raise ValueError('Outside verified monotone signed geometry')
        breaks=[0.,2.]
        if k:
            for ep in cut:
                p=np.sqrt((ep-M)*(ep+M))
                eb=ep-k-Es
                kb=np.sqrt((eb-m)*(eb+m))
                A=(ep-k)*Es-(M*M+s*s-m*m)/2
                low=M*M/(ep+p)
                for v in (q*q-A-q*kb,q*q-A+q*kb):
                    t=(v/k-low)/p
                    if 0.<t<2.:
                        breaks.append(t)
        shift=(m/((self.lo+self.hi)/2))**2/2
        X,W=np.asarray(gl(n),dtype=np.longdouble);x,w=np.asarray(gl(o),dtype=np.longdouble)
        phi=np.cos((np.arange(az,dtype=np.longdouble)+.5)*2*np.pi/az)
        out=np.zeros((2,3),dtype=np.longdouble)
        for thl,thr in zip(sorted(set(breaks))[:-1],sorted(set(breaks))[1:],strict=True):
            phase=np.log([thl+shift,thr+shift])
            zz=(phase[0]+phase[1])/2+np.diff(phase)[0]*X/2
            for t,sup in zip(np.exp(zz)-shift,W*np.diff(phase)[0]*np.exp(zz)/4,strict=True):
                def bounds(E,sign):
                    p=np.sqrt((E-M)*(E+M))
                    el=M*M/(E+p)+p*t
                    zR=np.sqrt((p-k)**2+2*p*k*t)
                    formula=((E-k)*Es-(M*M+s*s-m*m)/2+k*el)-sign*q*zR
                    dP=E/p
                    dL=-(M*M/(E+p)**2)*(1+dP)+t*dP
                    dR=dP*(p-k+k*t)/zR
                    return formula,Es+k*dL-sign*q*dR
                def endpoint(E,sign):
                    test_lo=np.longdouble('.18' if self.emission_only_extension else '.22');test_hi=np.longdouble('6.')
                    f_lo,d_lo=bounds(test_lo,sign)
                    f_hi,d_hi=bounds(test_hi,sign)
                    if not f_lo<0.<f_hi or min(d_lo,d_hi)<np.longdouble('.00043'):
                        raise ArithmeticError(('Broad signed bracket invalid',k,t,sign))
                    for z in range(10):
                        a,b=bounds(E,sign)
                        if abs(a)<np.longdouble('3e-20'):
                            break
                        if a>0: test_hi=E
                        else: test_lo=E
                        v=E-a/b
                        E=v if test_lo<v<test_hi else (test_hi+test_lo)/2
                    if abs(bounds(E,sign)[0])>6e-20:
                        raise ValueError('Broad exact support boundary did not close.')
                    return E
                left=endpoint(self.lo,-1.)
                right=endpoint(self.hi,1.)
                xx=np.r_[left,cut[(cut>left)&(cut<right)],right]
                hs=np.diff(xx)/2
                ep=((xx[1:]+xx[:-1])/2)[:,None]+hs[:,None]*x
                p=np.sqrt((ep-M)*(ep+M))
                el=M*M/(ep+p)+p*t
                Qr=np.sqrt((p-k)**2+2*p*k*t)
                cs=((ep-k)*Es-(M*M+s*s-m*m)/2+k*el)/(q*Qr)
                if max(abs(cs).ravel())>1+1e-9: raise ValueError('Unexpected parent support.')
                z=1-t
                zR=(p*z-k)/Qr
                si=p*np.sqrt(t*(2-t))/Qr
                Bn=el[...,None]-Es+q*((cs*zR)[...,None]+(si*np.sqrt(np.maximum(0,1-cs*cs)))[...,None]*phi)
                a,b=-el[...,None],-Bn
                d=(M*M+m*m-s*s)/2+k*(Bn-el[...,None])
                V=(1/a-1/b)/2
                H=M*M/(a*a)+m*m/(b*b)-2*d/(a*b)
                U=b*M*M/a+a*m*m/b-2*d
                A=-4*(d+M*m)*H-8*k*V*U+16*k*k*V*V*a*b
                Wtot=(hs[:,None]*w)*A.mean(axis=-1)*(p/Qr)/(2*((M+m)**2-s*s))*sup
                fp=self.pf[np.searchsorted(self.pe,ep)]
                fb=self.bf[np.searchsorted(self.be,ep-k-Es)]
                factors=np.asarray([fp*(1-fb),(1-fp)*fb,fp-fb])
                out[0]+=np.sum(factors*Wtot,axis=(1,2))
                outside=(ep<self.lo)|(ep>self.hi)
                out[1]+=np.sum(factors*(Wtot*outside),axis=(1,2))
        return out



def run():
    from finite_mass_qed_parent_integrated_signed_soft_reference import Source
    out=[]
    for g in ('fine','third'):
        with gzip.open(base/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        X=BridgeSource(d);old=Source(d)
        old_width=max(np.max(abs(X.signed(z,32,8,24)-old.signed(z,32,8,24))) for z in (-1e-5,0.,1e-5))
        if old_width>4e-13:
            raise ArithmeticError(('Original thin-support normalization lost',old_width))
        v0=X.signed(0.,32,8,24)[0,:2]
        fs=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],float(X.q))-1]
        results=[]
        for k in (1e-5,2e-5,4e-5,8e-5,1.6e-4,2e-4):
            def one(n,o,p):
                A=X.signed(k,n,o,p)[0,:2]
                B=X.signed(-k,n,o,p)[0,:2]
                C=X.signed(0.,n,o,p)[0,:2]
                q=np.longdouble(k)
                T=q/np.expm1(q)
                # Pair signed spectra before integrating over k; retain
                # the vacuum and thermal event coefficients separately.
                cold=np.array([A[0]-C[0],B[1]-C[1]])/q
                warm=T*(A+B-2*C)/(q*q)
                return np.concatenate([A,B,cold,warm])
            a=one(38,10,28); b=one(62,14,48)
            base_net=v0[0]/fs+v0[0]-v0[1]
            dnet=(b[4]+b[6])/fs+(b[4]+b[6])-(b[5]+b[7])
            results.append(dict(k=k,bare_emission_G_L=b[:2].astype(float).tolist(),
                                bare_absorption_G_L=b[2:4].astype(float).tolist(),
                                vacuum_subtracted_G_L=b[4:6].astype(float).tolist(),
                                thermal_subtracted_G_L=b[6:].astype(float).tolist(),
                                max_mesh_abs_residual=float(max(abs(a[4:]-b[4:]))),
                                local_net_derivative_over_Born_eik_scale=float(dnet/base_net)))
        out.append(dict(grid=g,original_matching_max_abs_defect=float(old_width),rows=results))
    return dict(scope='Leading soft counterterm only, exact moving signed Born apertures on a finite bridge. No finite virtual matching or physical width.',data=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
