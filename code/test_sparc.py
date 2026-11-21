"""Synthetic tests only; no real galaxy data loaded."""
import tempfile, unittest
from pathlib import Path
import numpy as np
from sparc_likelihood import Datum,read_sparc,nu_from_mass,predict_model,diagonal_loglike
class TestSPARC(unittest.TestCase):
    def test_catalog(self):
        src="Title: SPARC\nCamB 3.36 0.16 1.99 1.50 1.86 3.75 0.00 30.32 0.00\nCamB 3.36 0.41 4.84 1.50 4.24 9.47 0.00 23.77 0.00\n"
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"data.mrt";p.write_text(src)
            d=read_sparc(p);self.assertEqual(len(d["CamB"]),2)
            self.assertAlmostEqual(d["CamB"][0].Vgas,1.86)
    def test_shared_mass_relation(self):
        self.assertAlmostEqual(nu_from_mass(86.22456551472253,0.6,44.1),2.61,places=7)
    def test_signed_gas_contribution(self):
        a=Datum("test",10,.6,45,2,-12,0,0,1,0)
        b=Datum("test",10,.6,45,2,+12,0,0,1,0)
        neg=predict_model([a],86.22456551472253,.6,44.1)[0]
        pos=predict_model([b],86.22456551472253,.6,44.1)[0]
        self.assertAlmostEqual(pos**2-neg**2,288.,places=6)
    def test_finite_log_likelihood(self):
        g=[Datum("g",10,.6,45,2,0,5,0,0,0)]
        self.assertTrue(np.isfinite(diagonal_loglike(g,86.22456551472253,.6,44.1)))
if __name__=="__main__":unittest.main()
