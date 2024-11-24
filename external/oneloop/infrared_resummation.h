
#include "header.h"

#ifndef __EFT_IR_RESUMMATION__
#define __EFT_IR_RESUMMATION__

struct indexed_real_arg
{
  int index;
  double ln_k;
};

struct indexed_rsd_arg
{
  int index;
  double ln_k;
  double mu;
};

double eft_ir_sigma2(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo, 
            double z, 
            double k_split,
            double k_bao
            );

double eft_ir_dsigma2(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo, 
            double z, 
            double k_split,
            double k_bao
            );

int eft_real_argument_list_rect(
            double * ln_kvec,
            int k_size,
            int mu_size,
            struct indexed_real_arg ** vec,
            ErrorMsg errmsg
            );

int eft_ir_pk_lo(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            enum linear_or_logarithmic mode,
            struct indexed_real_arg * vec,
            int vec_size,
            double z,
            int linear_spectrum_index,
            double sigma2_ir_at_z,
            double * out_pk
            );

int eft_ir_pk_nlo(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            enum linear_or_logarithmic mode,
            struct indexed_real_arg * vec,
            int vec_size,
            double z,
            int linear_spectrum_index,
            double sigma2_ir_at_z,
            double * out_pk
            );

int eft_rsd_argument_list_rect(
            double * ln_kvec,
            int k_size,
            double * muvec,
            int mu_size,
            struct indexed_rsd_arg ** vec,
            ErrorMsg errmsg
            );

int eft_rsd_argument_list(
            double * ln_kvec,
            double * muvec,
            int size,
            struct indexed_rsd_arg ** vec,
            ErrorMsg errmsg
            );

int eft_ir_pk_rsd_lo(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            enum linear_or_logarithmic mode,
            struct indexed_rsd_arg * vec,
            int vec_size,
            double z,
            double f,
            int linear_spectrum_index,
            double sigma2_ir_at_z,
            double dsigma2_ir_at_z,
            double * out_pk
            );

int eft_ir_pk_rsd_nlo(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            enum linear_or_logarithmic mode,
            struct indexed_rsd_arg * vec,
            int vec_size,
            double z,
            double f,
            int linear_spectrum_index,
            double sigma2_ir_at_z,
            double dsigma2_ir_at_z,
            double * out_pk
            );

// int eft_compute_pk_moments(
//             struct precision * ppr,
//             struct background * pba,
//             struct primordial * ppm,
//             struct fourier * pfo, 
//             struct eft_hyper_parameters hp,
//             struct eft * peft,
//             double * biases,
//             int biases_size,
//             double z, 
//             double sigma2_ir_at_z,
//             enum eft_pk_type use_pk_type,
//             double * moments,
//             int * moments_n,
//             int * moments_size
//             );

#endif
