#include "common.h"
#include <complex.h>
#undef I

#ifndef __EXTERNAL_STORAGE__
#define __EXTERNAL_STORAGE__

extern struct background;
extern struct thermodynamics;
extern struct perturbations;
extern struct primordial;
extern struct fourier;
extern struct transfer;
extern struct harmonic;
extern struct lensing;
extern struct distortions;
extern struct eft;

struct ext_storage
{
  /** - EFT loop matrix pointers */
  double complex *** loop_matrices;
  int ** loop_matrices_size;
  short ** symmetry;
  int loop_matrices_stored; // = _FALSE_;
  int eft_index_num; // = 0;
  int eft_size; // = 0;

  short update_required;  // = _TRUE_;
};

//int ext_alloc(struct ext_storage * pext);

//int ext_cleanup(struct ext_storage * pext);

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
                   const int index,
                   ErrorMsg errmsg);


#endif