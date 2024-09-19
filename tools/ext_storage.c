/** @file ext_storage.c
 *
 * author: Christian Radermacher, 2024
 *
 * this module offers memory storage for pointers that are reused in repeated CLASS evaluations
*/

#include "ext_storage.h"
#include "header.h"

int ext_init(struct ext_storage * pext) {
  pext->eft_size = 0;
  pext->loop_matrices_stored = _FALSE_;

  return _SUCCESS_;
}

int ext_cleanup(struct ext_storage * pext) {
  int i, j;

  if (pext->loop_matrices_stored) {
    for (i = 0; i < pext->eft_size; i++) {
      if (pext->loop_matrices[i]) { /** - if memory block is still linked, delete it */
        for (j = 0; j < pext->eft_index_num[i]; j++) {
          free(pext->loop_matrices[i][j]);
        }
        free(pext->loop_matrices[i]);
      }
      if (pext->loop_matrices_size[i]) { free(pext->loop_matrices_size[i]); }
      if (pext->symmetry[i]) { free(pext->symmetry[i]); }
      if (pext->use_tracer[i]) { free(pext->use_tracer[i]); }
      if (pext->spectra_contributions_dimension[i]) { free(pext->spectra_contributions_dimension[i]); }

      free(pext->period[i]);
    }
    free(pext->loop_matrices);
    free(pext->loop_matrices_size);
    free(pext->symmetry);
    free(pext->use_tracer);
    free(pext->spectra_contributions_dimension);

    free(pext->period);
    free(pext->fourier_coeff_size);
    free(pext->eft_index_num);
  }

  pext->loop_matrices_stored = _FALSE_;
  pext->eft_size = 0;

  // pext->update_required = _TRUE_;

  return _SUCCESS_;
}

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

  /** - clean the storage: this will delete memory blocks that are still associated to the pointers in the ext-structure */
  ext_cleanup(pext);

  /**----------------- EFT loop matrices -----------------*/
  if (pfo && pfo->method == nl_oneloopPT) {
    /** - count master eft-structures */
    for (i = 0; i < pfo->eft_size; i++) {
      if ((pfo->peft + i)->role == eft_master) { pext->eft_size += 1; }
    }

    /** - allocate pointer arrays */
    class_alloc(pext->loop_matrices, pext->eft_size*sizeof(class_complex **), pext->error_message);
    class_alloc(pext->loop_matrices_size, pext->eft_size*sizeof(int *), pext->error_message);
    class_alloc(pext->symmetry, pext->eft_size*sizeof(short *), pext->error_message);
    class_alloc(pext->use_tracer, pext->eft_size*sizeof(short *), pext->error_message);
    class_alloc(pext->spectra_contributions_dimension, pext->eft_size*sizeof(short *), pext->error_message);

    class_alloc(pext->period, pext->eft_size*sizeof(double *), pext->error_message);
    class_alloc(pext->fourier_coeff_size, pext->eft_size*sizeof(int), pext->error_message);
    for (i = 0; i < pext->eft_size; i++) {
      class_alloc(pext->period[i], eft_tracer_num*sizeof(double), pext->error_message);
    }
    class_alloc(pext->eft_index_num, pext->eft_size*sizeof(int), pext->error_message);
    for (i = 0; i < pfo->eft_size; i++) {
      if ((pfo->peft + i)->role == eft_master) { pext->eft_index_num[i] = (pfo->peft + i)->index_num; }
    }

    /** - copy pointers from eft-structures */
    for (i = 0; i < pfo->eft_size; i++) {
      if ((pfo->peft + i)->role == eft_master) {
        pext->loop_matrices[i] = (pfo->peft + i)->loop_matrices;
        pext->loop_matrices_size[i] = (pfo->peft + i)->loop_matrices_size;
        pext->symmetry[i] = (pfo->peft + i)->symmetry;
        pext->use_tracer[i] = (pfo->peft + i)->use_tracer;
        pext->spectra_contributions_dimension[i] = (pfo->peft + i)->spectra_contributions_dimension;

        for (j = 0; j < eft_tracer_num; j++)
          pext->period[i][j] = (pfo->peft + i)->hp->period[j];  /** data copied */

        pext->fourier_coeff_size[i] = (pfo->peft + i)->hp->fourier_coeff_size;  /** data copied */
      }
    }

    /** - replace with NULL at the original place */
    for (i = 0; i < pfo->eft_size; i++) {
      if ((pfo->peft + i)->role == eft_master) {
        (pfo->peft + i)->loop_matrices = NULL;
        (pfo->peft + i)->loop_matrices_size = NULL;
        (pfo->peft + i)->symmetry = NULL;
        (pfo->peft + i)->use_tracer = NULL;
        (pfo->peft + i)->spectra_contributions_dimension = NULL;
      }
    }

    pext->loop_matrices_stored = _TRUE_;
  }
  /**-----------------------------------------------------*/


  return _SUCCESS_;
}

int ext_insert_eft(struct ext_storage * pext,
                   struct eft * peft,
                   const int index_eft,
                   const int num_matrices,
                   ErrorMsg errmsg) {

  int index_moment, index_tracer;

  if (pext) {
    if (!(pext->loop_matrices_stored)) {  /** - if the matrices are not stored, then return immediately (not considered an error) */
      class_sprintf(errmsg, "Nothing is stored in pext = %p", pext);
      return _FAILURE_;
    }
    else if (index_eft >= pext->eft_size || pext->eft_index_num[index_eft] < num_matrices) {  /** - if they are stored, but don't conform to the request, return an error message */
      class_test_message(errmsg, "index_eft >= pext->eft_size || pext->eft_index_num[index_eft] < num_matrices", \
                        "Error in ext_insert_eft: pext = %p, stored = %d, stored number of indices = %d, stored number of matrices = %d", \
                        pext, pext->loop_matrices_stored, pext->eft_size, pext->eft_index_num[index_eft]);
      // /** - query an update of the storage, flag will be reset at the update */
      // pext->update_required = _TRUE_;

      return _FAILURE_; /** - insertion is impossible */
    }
    else {  /** - attempt insertion */
      for (index_tracer = 0; (index_tracer < eft_tracer_num)  \
                            && (pext->period[index_eft][index_tracer] == peft->hp->period[index_tracer])  \
                            && (pext->fourier_coeff_size[index_eft] == peft->hp->fourier_coeff_size); index_tracer++);
      if (index_tracer == eft_tracer_num) { /** - stored period and size is equal to the one requested in peft->hp */
        /** - copy the pointers back over to the eft-structure */
        peft->loop_matrices = pext->loop_matrices[index_eft];
        peft->loop_matrices_size = pext->loop_matrices_size[index_eft];
        peft->symmetry = pext->symmetry[index_eft];
        peft->use_tracer = pext->use_tracer[index_eft];
        peft->spectra_contributions_dimension = pext->spectra_contributions_dimension[index_eft];

        /** - disassociate the pointers in the ext-structure */
        pext->loop_matrices[index_eft] = NULL;
        pext->loop_matrices_size[index_eft] = NULL;
        pext->symmetry[index_eft] = NULL;
        pext->use_tracer[index_eft] = NULL;
        pext->spectra_contributions_dimension[index_eft] = NULL;

        /** - write the number of matrices copied */
        peft->moments_allocated = pext->eft_index_num[index_eft];
        pext->eft_index_num[index_eft] = 0;

        if (peft->hp->eft_verbose > 0) {
          printf("Kernel matrices inserted from external storage for index = %d. \n", index_eft);
        }
      }
      else {  /** - mismatch btw. periods */
        /** TODO: cleanup (dump memory) */
        /** guess: cleanup will be done automatically when calling ext_save */
        class_protect_sprintf(errmsg, "Mismatch of kernel matrix specifications for index = %d: Requested period in ln(k)-space and size is (%.4f, %d), external storage has (%.4f, %d).", \
                                      index_eft, peft->hp->period[index_tracer], peft->hp->fourier_coeff_size, pext->period[index_eft][index_tracer], pext->fourier_coeff_size[index_eft]);
        return _FAILURE_;
      }

      return _SUCCESS_;
    }
  }
  else {
    return _FAILURE_; /** - not considered an error, therefore no message */
  }

  return _SUCCESS_;
}
