
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
            const double z, 
            const double k_split,
            const double k_bao
            );

double eft_ir_dsigma2(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo, 
            const double z, 
            const double k_split,
            const double k_bao
            );

int eft_real_argument_list_rect(
            const double * const ln_kvec,
            const int k_size,
            const int mu_size,
            struct indexed_real_arg ** vec,
            ErrorMsg errmsg
            );

int eft_ir_pk_lo(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            enum linear_or_logarithmic mode,
            const struct indexed_real_arg * const vec,
            const int vec_size,
            const double z,
            const int linear_spectrum_index,
            const double sigma2_ir_at_z,
            double * out_pk
            );

int eft_ir_pk_nlo(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            enum linear_or_logarithmic mode,
            const struct indexed_real_arg * const vec,
            const int vec_size,
            const double z,
            const int linear_spectrum_index,
            const double sigma2_ir_at_z,
            double * out_pk
            );

int eft_rsd_argument_list_rect(
            const double * const ln_kvec,
            const int k_size,
            const double * const muvec,
            const int mu_size,
            struct indexed_rsd_arg ** vec,
            ErrorMsg errmsg
            );

int eft_rsd_argument_list(
            const double * const ln_kvec,
            const double * const muvec,
            const int size,
            struct indexed_rsd_arg ** vec,
            ErrorMsg errmsg
            );

int eft_ir_pk_rsd_lo(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            enum linear_or_logarithmic mode,
            const struct indexed_rsd_arg * const vec,
            const int vec_size,
            const double z,
            const double f,
            const int linear_spectrum_index,
            const double sigma2_ir_at_z,
            const double dsigma2_ir_at_z,
            double * out_pk
            );

int eft_ir_pk_rsd_nlo(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            enum linear_or_logarithmic mode,
            const struct indexed_rsd_arg * const vec,
            const int vec_size,
            const double z,
            const double f,
            const int linear_spectrum_index,
            const double sigma2_ir_at_z,
            const double dsigma2_ir_at_z,
            double * out_pk
            );

// int eft_compute_pk_moments(
//             struct precision * ppr,
//             struct background * pba,
//             struct primordial * ppm,
//             struct fourier * pfo, 
//             const struct eft_hyper_parameters hp,
//             struct eft * peft,
//             double * biases,
//             int biases_size,
//             const double z, 
//             const double sigma2_ir_at_z,
//             enum eft_pk_type use_pk_type,
//             double * moments,
//             int * moments_n,
//             int * moments_size
//             );

#endif
