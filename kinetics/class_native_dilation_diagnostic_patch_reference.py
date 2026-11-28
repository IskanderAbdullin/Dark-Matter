"""Diagnostic CLASS v3.3.4 perturbation-only rescaling, not a production solver."""
import argparse
from pathlib import Path

ANCHOR="    pba->factor_ncdm[k]=pba->deg_ncdm[k]*4*_PI_*pow"
PATCH=r'''    /* Diagnostic only: set weighted number dilation to its continuum identity.
       Deliberately leave the native nodes, weights, f0, and background alone. */
    if ((pba->got_files!=NULL) && (pba->got_files[k]==_TRUE_) &&
        pba->ncdm_quadrature_strategy[k]==qm_trapz) {
      double a0=0.0, d0=0.0, weight, scale;
      for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
        q = pba->q_ncdm[k][index_q];
        weight = q*q*pba->w_ncdm[k][index_q];
        a0 += weight;
        d0 += weight*pba->dlnf0_dlnq_ncdm[k][index_q];
      }
      class_test(!(a0>0.0 && d0<0.0),pba->error_message,
                 "Invalid diagnostic ncdm derivative number moment");
      scale = -3.0*a0/d0;
      printf("DILATION_DIAGNOSTIC n=%d raw_residual=%+.14e factor=%.14e\n",
             pba->q_size_ncdm[k],d0/(-3.0*a0)-1.0,scale);
      for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
        pba->dlnf0_dlnq_ncdm[k][index_q] *= scale;
      }
    }

'''


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("original_background_c",type=Path)
    p.add_argument("output_background_c",type=Path)
    a=p.parse_args()
    assert a.original_background_c.resolve()!=a.output_background_c.resolve()
    s=a.original_background_c.read_text()
    assert s.count(ANCHOR)==1 and "DILATION_DIAGNOSTIC" not in s
    a.output_background_c.write_text(s.replace(ANCHOR,PATCH+ANCHOR))
