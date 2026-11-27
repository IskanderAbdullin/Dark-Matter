"""Positive, locally balanced nodal energy quadrature of early paired decay."""

import argparse
import json

import numpy as np
from scipy.special import expit

from early_pole_decay_mode_clock_reference import G, KAPPA
from early_pole_paired_hartree_matrix_reference import DirectDriver


class NodalDriver(DirectDriver):
    """Massless fermion; integrate geometric hat weights before sampling FD."""
    def __init__(self, ns, nf, lam, qlow=1e-9):
        super().__init__(ns, nf, lam, 0., qlow)
        # Interior nodes closely match old cell centroids; edge nodes allow
        # non-negative barycentric weights over the entire finite support.
        old = self.fed.copy()
        self.nodes = np.r_[old[0], (old[:-1] + old[1:]) / 2, old[-1]]
        self.fed = np.r_[old[0], (self.nodes[:-1] + self.nodes[1:]) / 2, old[-1]]
        self.ef = self.fed
        self.wf = 4 * np.diff(self.fed**3) / (6 * np.pi**2)

    def geometric(self, x, m):
        p = self.p
        E = np.hypot(p, m)
        if x <= m:
            return E, np.zeros((len(E), len(self.nodes)))
        A0, D = (x*x + m*m) / 2, (x*x - m*m) / 2
        lower = (x*x*E*E + D*D) / (A0*E+p*D) - E
        span = 2*p*D/(m*m)
        # For each adjacent node interval, the integral of its linear hats
        # is its intersected width multiplied by the weights at its midpoint.
        left, right = self.nodes[:-1], self.nodes[1:]
        lo = np.maximum(lower[:, None], left)
        hi = np.minimum((lower + span)[:, None], right)
        width = np.minimum(np.maximum(hi-lo, 0.), span[:, None])
        theta = np.clip((lo + width/2 - left)/(right-left), 0., 1.)
        q = np.zeros((len(E), len(self.nodes)))
        q[:, :-1] = width*(1-theta)
        q[:, 1:] += width*theta
        return E, q

    def kernels(self, x, m):
        E, w = self.geometric(x, m)
        D = w * (KAPPA * G * x * max(x*x-m*m, 0.) / (E*self.p))[:, None]
        # Effective parent energy per *discrete* channel; not a claim that
        # each such point inside the represented cell is kinematically exact.
        return D * expit(-E[:, None]-self.nodes), D

    def matrix_check(self):
        x, m = .006, 5e-5
        E, w = self.geometric(x, m)
        p = self.p
        Delta = (x*x-m*m)/2
        lower = (x*x*E*E+Delta*Delta)/(((x*x+m*m)/2)*E+p*Delta)-E
        upper = lower+2*p*Delta/(m*m)
        left = np.maximum(lower, self.nodes[0])
        right = np.minimum(upper, self.nodes[-1])
        width = np.maximum(right-left, 0)
        first = width * (left + width/2)
        assert np.max(abs(w.sum(axis=1)-width)) < 4e-14
        assert np.max(abs(w@self.nodes-first)) < 5e-13
        A, D = self.kernels(x, m)
        f, b = 1/np.expm1(E), expit(-self.nodes)
        forward = A*(1+f[:, None])*(1-b)
        reverse = (D-A)*f[:, None]*b
        keep = forward > np.max(forward)*1e-13
        balance = np.max(abs(reverse[keep]/forward[keep]-1))
        assert balance < 1e-12
        # The operator is inward at every boundary and conserves shared
        # count. At a generic interior state each channel deposits E_i+E_j.
        u = 3*np.exp(-self.p)
        b = .2*expit(-self.nodes)
        R = A*(1-b)+u[:, None]*(A-D*b)
        v = self.rhs_at_m(x,u,b,m)
        primal, daughter = self.w*v[:self.n], self.wf*v[self.n:]
        lhs = primal@E + daughter@self.nodes
        rhs = np.sum(self.w[:, None]*R*(E[:, None]+self.nodes))
        assert abs(lhs-rhs) <= 1e-13*abs(rhs)
        assert abs(primal.sum()-daughter.sum()) <= 1e-13*abs(primal.sum())
        assert np.min(self.rhs_at_m(x,u,np.zeros_like(b),m)[self.n:]) >= 0
        assert np.max(self.rhs_at_m(x,u,np.ones_like(b),m)[self.n:]) <= 0
        old = DirectDriver(self.n-2, len(self.nodes), self.lam, 0., self.ed[0])
        old.p = self.p
        old.ef = np.r_[self.nodes[0], self.nodes[-1]]
        exact = old.kernels(x,m)[0][:,0]
        active = exact > np.max(exact)*1e-7
        return dict(hat_width_max_abs=float(np.max(abs(w.sum(axis=1)-width))),
                    hat_energy_max_abs=float(np.max(abs(w@self.nodes-first))),
                    thermal_detailed_balance_max_relative=float(balance),
                    signed_projected_energy_identity_relative=float(abs(lhs/rhs-1)),
                    empty_daughter_seed_relative_error_max=float(
                        np.max(abs(A.sum(axis=1)[active]/exact[active]-1))))

    def run(self, xend):
        out = super().run(xend)
        f = np.array(out['final_daughter_per_spin_and_charge_occupations'])
        n = f*self.wf
        # Override inherited cell-centroid diagnostic with the nodal energy.
        for r in out['checkpoints']:
            del r['fermion_mean_p_over_T']
        out['final_fermion_mean_p_over_T'] = float(n@self.nodes/n.sum())
        out['fermion_energy_nodes_over_T'] = self.nodes.tolist()
        out['quadrature_checks'] = out.pop('empty_daughter_row_sum_control')
        return out


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar', type=int, default=65)
    p.add_argument('--fermion', type=int, default=90)
    p.add_argument('--xend', type=float, default=.05)
    p.add_argument('--static-only', action='store_true')
    a = p.parse_args()
    d = NodalDriver(a.scalar, a.fermion, 7.65175e-8)
    print(json.dumps(d.matrix_check() if a.static_only else d.run(a.xend),
                     indent=2, allow_nan=False))
