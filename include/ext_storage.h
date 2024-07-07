#include "common.h"

#ifndef __EXTERNAL_STORAGE__
#define __EXTERNAL_STORAGE__

// extern struct background;
// extern struct thermodynamics;
// extern struct perturbations;
// extern struct primordial;
// extern struct fourier;
// extern struct transfer;
// extern struct harmonic;
// extern struct lensing;
// extern struct distortions;
// extern struct eft;

#include "fourier.h"
#include "carray.h"
#include "common_header.h"

struct ext_storage
{
  /** - EFT loop matrix pointers */
  class_complex *** loop_matrices; /**< pointers to loop_matrices; indexed as loop_matrices[index_eft][index_moment] */
  int ** loop_matrices_size;  /**< pointers to loop_matrices_size; indexed as loop_matrices_size[index_eft] */
  short ** symmetry;  /**< pointers to symmetry type; indexed as symmetry[index_eft] */
  short ** use_tracer;  /**< pointers to tracer type; indexed as use_tracer[index_eft] */
  short ** spectra_contributions_dimension;  /**< pointers to dimensions of spectra contributions; indexed as spectra_contributions_dimension[index_eft] */

  double ** period; /**< array of logarithmic period of tracer types; indexed as period[index_eft][index_tracer] */
  int loop_matrices_stored; // = _FALSE_;
  int * eft_index_num;  /**< array of index_num for each eft-structure */
  int eft_size; // = 0;

  // short update_required;  // = _TRUE_;

  ErrorMsg error_message;
};

int ext_init(struct ext_storage * pext);

int ext_cleanup(struct ext_storage * pext);

int ext_save(struct ext_storage * pext,
            struct background * pba,
            struct thermodynamics * pth,
            struct perturbations * ppt,
            struct primordial * ppm,
            struct fourier * pfo,
            struct transfer * ptr,
            struct harmonic * phr,
            struct lensing * ple,
            struct distortions * psd);

int ext_insert_eft(struct ext_storage * pext,
                   struct eft * peft,
                   const int index_eft,
                   const int num_matrices,
                   ErrorMsg errmsg);

#endif
