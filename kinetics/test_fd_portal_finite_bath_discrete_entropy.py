"""Check analytic fixed-energy thermodynamics against direct perturbations."""
import json
from pathlib import Path
import unittest

import numpy as np

from fd_portal_finite_bath_discrete_entropy import Bath, C_VIS


class FiniteBathThermodynamics(unittest.TestCase):
    def test_direct_fixed_energy_entropy_gradient(self):
        x=np.array([-.65,np.log(.34)])
        a=Bath(.3,x,168)
        y,L,g,_,A=a.point(x)
        t=2e-6*np.eye(2)
        z=np.array([(a.point(x+v)[1]-a.point(x-v)[1])/(4e-6) for v in t])
        np.testing.assert_allclose(z,g@A,rtol=2e-8)

    def test_stored_states_preserve_total_energy(self):
        r=Path(__file__).parent/'results/fd-portal-finite-bath-discrete-entropy-2026-11-28.json'
        for z in json.loads(r.read_text())['rows']:
            u=[z['initial_eta'],np.log(z['initial_dark_T'])]
            a=Bath(z['initial_visible_T'],u,168)
            for c in z['levels']:
                u=[c['eta'],np.log(c['dark_T'])]
                E=C_VIS*c['visible_T']**4+a.dark(u)[0][1]
                self.assertAlmostEqual(E/z['total_energy'],1.,places=13)
