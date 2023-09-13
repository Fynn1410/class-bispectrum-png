#include "ext_storage.h"

int ext_save(struct ext_storage * pext,
            struct background * pba,
            struct thermodynamics * pth,
            struct perturbations * ppt,
            struct primordial * ppm,
            struct fourier * pfo,
            struct transfer * ptr,
            struct harmonic * phr,
            struct lensing * ple,
            struct distortions * psd) {

  int i, j;

  /** EFT loop matrices */
  // if (pfo->eft_hp.use_time_independent_kernels) {
  //   for (j = 0; j < pfo->peft->index_num; j++) {
  //     pext->loop_matrices[j] = pfo->peft->loop_matrices[j];
  //     pfo->peft->loop_matrices[j] = NULL;
  //   }
  // }
  // else {   not yet implemented
  //   for (i = 0; i < pfo->tau_size; i++) {
      
  //   }
  // }

  return _SUCCESS_;
}

int ext_insert_eft(struct ext_storage * pext,
                   struct eft * peft,
                   const int index,
                   ErrorMsg errmsg) {


  /** - check matching configuration (size & symmetry & period) */
                   
  return _SUCCESS_;
}