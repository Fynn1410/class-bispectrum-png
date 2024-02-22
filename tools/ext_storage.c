/** @file ext_storage.c
 * 
 * author: Christian Radermacher, 2024
 * 
 * this module offers memory storage for pointers that are reused in repeated CLASS evaluations
*/

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
                   const int num_matrices,
                   ErrorMsg errmsg) {

  if (pext) {
    if (!(pext->loop_matrices_stored) || index >= pext->eft_size || pext->eft_index_num < num_matrices) {
      class_test_message(errmsg, "!pext || !(pext->loop_matrices_stored) || index >= pext->eft_size || pext->eft_index_num < num_matrices", \
                        "Error in ext_insert_eft: pext = %p, stored = %d, stored number of indices = %d, stored number of matrices = %d", \
                        pext, pext->loop_matrices_stored, pext->eft_size, pext->eft_index_num);
      /** - query an update of the storage, flag will be reset at the update */
      pext->update_required = _TRUE_;

      return _FAILURE_; /** - insertion is impossible */
    }
  }
  else {
    return _FAILURE_; /** - not considered an error, therefore no message */
  }
                   
  return _SUCCESS_;
}