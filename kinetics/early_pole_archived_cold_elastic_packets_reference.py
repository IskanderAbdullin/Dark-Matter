"""Cold-window reference with a reciprocal strict hard distinct-fermion elastic subset."""
import argparse
import json

import numpy as np
from scipy.special import expit, roots_legendre

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_live_joint_hard_pair_packets_reference import gauss_each


class Elastic(Cold):
    """Vacuum massless t-channel weights, dressed conservative discrete energy remap."""

    def __init__(self, grid, wide, count, exposure, order, cut, strength):
        super().__init__(grid, wide, 0, count, True, exposure=exposure)
        self.cut = cut
        self.strength = strength
        self.reject_E = 0.
        self.uncut_s_weight = 0.
        self.graph = self.packets(order)
        H = np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18
        self.factor = .005**4/H

    def packets(self, order):
        nr, ni, nt, na = order
        M = self.M
        z, w = roots_legendre(ni)
        edges = [np.array([M.cells[span.start].lo, *[c.hi for c in M.cells[span]]])
                 for span in M.r[1:]]
        (pa, da), (pb, db) = [gauss_each(e, nr) for e in edges]
        moment = lambda q, dq: 2*dq*q*q/np.pi**2
        a, b = np.meshgrid(pa, pb, indexing='ij')
        inc = np.outer(moment(pa, da), moment(pb, db))[..., None]
        lim = 4*a*b
        left = np.minimum(1., 1/lim)
        u = left[..., None]+(1-left[..., None])*(1+z)/2
        W = inc*(1-left[..., None])*u*w/2
        # Distinct initial flavors. Factor half compensates for integrating
        # forward minus reverse over the full reciprocal domain.
        ok = W > 0.
        s = (lim[..., None]*u)[ok]
        P, B = (np.broadcast_to(p[..., None], W.shape)[ok] for p in (a, b))
        W = W[ok]
        E = P+B
        boost = np.sqrt(np.maximum(E*E-s, 0.))
        A = np.clip(np.divide(P-B, boost, out=np.zeros_like(boost), where=boost>0), -1., 1.)
        az, azw = roots_legendre(na)
        out = []
        high = [e[-1] for e in edges]
        low = [e[0] for e in edges]
        # Split at both outgoing momentum grid edges before dressing the energy.
        mass = [M.cells[span.start].m for span in M.r[1:]]
        angles = np.array([.15,.25,.4,.7,1.])
        tau, tw = gauss_each(np.r_[self.cut, angles[angles>self.cut]], nt)
        for t, wt in zip(tau, tw, strict=True):
            spin = (1+(1-t)**2)/(16*np.pi*s*t*t)
            ct = 1-2*t
            v = (E+boost*A*ct)/2
            r = boost*np.sqrt((1-ct*ct)*(1-A*A))/2
            totw = 2*wt*W*spin
            self.uncut_s_weight += totw.sum()
            aa = np.concatenate([(edges[0]-v[:,None])/np.maximum(r[:,None],1e-30),
                                 (E[:,None]-edges[1]-v[:,None])/np.maximum(r[:,None],1e-30),
                                 np.broadcast_to([-1.,1.],(len(E),2))],axis=1)
            azedge = np.sort(np.arccos(np.clip(aa,-1.,1.)),axis=1)
            ln = np.diff(azedge,axis=1)/2
            i,j = np.nonzero(ln > 1e-14)
            h = ln[i,j,None]
            m = (azedge[i,j]+azedge[i,j+1])[:,None]/2
            p = v[i,None]+r[i,None]*np.cos(m+h*az)
            b = E[i,None]-p
            w = totw[i,None]*h*azw/np.pi
            take = (p>=low[0])&(p<=high[0])&(b>=low[1])&(b<=high[1])
            out.append(np.array([np.broadcast_to(z[i,None],p.shape)[take] for z in
                                 (P,B,boost)]+[p[take],b[take],w[take]]))
        antiderivative = lambda t: -2/t - 2*np.log(t) + t
        exact = (W/(8*np.pi*s)).sum()*(antiderivative(1.)-antiderivative(self.cut))
        self.angular_normalization_error = abs(self.uncut_s_weight/exact - 1.)
        assert self.angular_normalization_error < 5e-5
        P, B, boost, Pf, Bf, W = np.concatenate(out,axis=1)
        ix = np.stack([np.searchsorted(e[1:],q)+f.start for q,e,f in
              zip((P,B,Pf,Bf),edges+edges,M.r[1:]*2,strict=True)])
        val = np.hypot([P,B,Pf,Bf],np.asarray(mass*2)[:,None])
        uv = (val-M.shape[ix,0])/M.shape[ix,1]
        assert np.max(abs(uv))<1.+2e-12 and min(W)>0.
        return ix[:3], uv[:3], Pf, boost, W

    def collision(self, e, f, grad):
        M = self.M
        ix,u,p,boost,W = self.graph
        d = M.r[2]
        val = e[ix,0]+e[ix,1]*u
        ec = val[0]+val[1]-val[2]
        hi = e[d,0]+e[d,1]
        assert min(np.diff(hi))>0.
        j = d.start+np.searchsorted(hi,ec)
        bounded = j<d.stop
        j = np.minimum(j,d.stop-1)
        v = (ec-e[j,0])/e[j,1]
        bare = M.shape[j,0]+M.shape[j,1]*v
        m = M.cells[d.start].m
        q = np.sqrt(np.maximum((bare-m)*(bare+m),0.))
        ok = bounded&(abs(v)<=1.)&(abs(p-q)<=boost)&(boost<=p+q)
        self.reject_E = max(self.reject_E, float(W[~ok].sum()/W.sum()))
        v = np.where(ok,v,0.)
        eta = f[ix,0]+f[ix,1]*u
        ef = f[j,0]+f[j,1]*v
        aff = eta[0]+eta[1]-eta[2]-ef
        A,B,C,D = expit(np.vstack((eta,ef)))
        R = self.factor*self.strength*W*ok*(A*B*(1-C)*(1-D)-C*D*(1-A)*(1-B))
        out = np.zeros_like(e)
        for k, w, n in zip([*ix,j],[*u,v],[-R,-R,R,R],strict=True):
            out[:,0] += np.bincount(k,weights=n,minlength=len(out))
            out[:,1] += np.bincount(k,weights=n*w,minlength=len(out))
        norm = max(1.,np.sum(abs(out)))
        energy = np.sum(grad*out)
        S = float(R@aff)
        assert abs(energy)<2e-12*norm and S>=-1e-13
        assert max(abs(np.sum(out[d,0])),abs(np.sum(out[M.r[1],0])))<2e-12*norm
        return out, S, energy, float(np.sum(grad[d]*out[d]))

    def rhs(self, y, pair):
        ans = super().rhs(y, pair)
        if self.strength:
            Y = y[:self.j0].reshape(-1,2)
            e,f,grad,_,_ = self.M.weights(Y,True)
            D,S,E,_ = self.collision(e,f,grad)
            ans[:self.j0] += D.ravel()
            ans[-2:] += [S,E]
        return ans

    def static(self):
        M = self.M
        r = []
        for mode in ('archive','common_T','parent_hotter'):
            Y = M.init.copy()
            if mode != 'archive':
                for _ in range(48):
                    old = Y.copy()
                    e = M.weights(Y,True)[2]
                    for j in range(M.r[1].start,len(Y)):
                        b = 1. if j<M.r[2].start or mode=='common_T' else 1/.65
                        mu = .1 if j<M.r[2].start else -.40/b
                        C = M.cells[j]
                        Y[j] = C.I@expit(-b*(e[j,0]+e[j,1]*C.u-mu))
                    if np.max(abs(Y-old))<3e-15:
                        break
            e,g,E,_,_ = M.weights(Y,True)
            D,S,work,heat = self.collision(e,g,E)
            if mode=='common_T':
                assert max(abs(heat),np.max(abs(D)))<2e-9
            if mode=='parent_hotter' and self.strength:
                assert S>0. and heat>0.
            r.append(dict(mode=mode,entropy=S,energy=work,daughter_heat=heat,rhs_inf=np.max(abs(D))))
        return r


def run(grid, wide, n, exposure, level, cut, strength, static):
    q = [(1,3,4,1),(2,4,5,2),(3,6,7,2)][level]
    M = Elastic(grid, wide, n, exposure, q, cut, strength)
    rows = M.static()
    rejection = M.reject_E
    out = dict(scope='Strict vacuum-hard subset with conservative dressed moment remap. '
               'Frozen projections and omitted soft/other hard gauge channels preclude a physical rate.',
               order=q, cut=cut, strength=strength, static=rows, grid=grid, wide=wide,
               packets=len(M.graph[-1]), angular_normalization_relative_error=M.angular_normalization_error,
               bare_geometry_retained=M.graph[-1].sum()/M.uncut_s_weight,
               static_maximum_dressed_remap_rejection=rejection)
    M.reject_E = 0.
    if not static:
        out['evolution'] = M.run()
        out['maximum_dressed_remap_rejection'] = M.reject_E
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',default='fine',choices=('fine','third'))
    p.add_argument('--wide',type=int,default=1,choices=(0,1,2))
    p.add_argument('--n',type=int,default=8)
    p.add_argument('--exposure',type=float,default=.0003)
    p.add_argument('--level',type=int,default=0,choices=range(3))
    p.add_argument('--cut',type=float,default=.15,choices=(.10,.15,.20,.25))
    p.add_argument('--strength',type=float,default=1.)
    p.add_argument('--static',action='store_true')
    args=p.parse_args()
    print(json.dumps(run(**vars(args)),indent=2,allow_nan=False))
