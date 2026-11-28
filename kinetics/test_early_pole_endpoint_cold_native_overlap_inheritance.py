"""Test inherited overlap loss separately from an actual refined optimum."""
import unittest
import numpy as np
from early_pole_endpoint_cold_native_overlap_inheritance import inherit,idx

class InheritedOverlap(unittest.TestCase):
    def test_refinement_can_lose_overlap_even_with_exact_parent_mass(self):
        U=np.array([0.,1.])
        F=np.array([0.,.5,1.])
        r=inherit(U,U,np.array([1.]),np.array([1.]),
                  F,F,np.array([.99,.01]),np.array([.01,.99]))
        for q in r['lost_from_original_primal_bounds']:
            self.assertAlmostEqual(q,.98)
        self.assertAlmostEqual(r['alternative_refined_common_mass'],.02)
        self.assertAlmostEqual(r['lifted_parent_dual_common_mass_upper'],1.)

    def test_geometry_requires_all_old_edges(self):
        with self.assertRaisesRegex(ArithmeticError,'refinement'):
            idx(np.array([0.,.5,1.]),np.array([0.,.6,1.]))
