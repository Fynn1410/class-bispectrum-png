#include "header.h"
#include "../../include/common.h"

#ifndef __EFT_POWER_SPECTRUM__
#define __EFT_POWER_SPECTRUM__

int eft_necessary_spectra_contributions(
                      struct eft * peft,
                      enum eft_pk_out_type pk_out_type,
                      int ** list_spectra_contributions,
                      int * list_spectra_contributions_size
                      );

int eft_necessary_pk_types_loops(
                      struct eft * peft,
                      enum eft_pk_out_type pk_out_type,
                      int ** list_pk_types,
                      int * list_pk_types_size
                      );

int eft_set_sampling_points_all(
                      struct eft * peft0,
                      int eft_size,
                      double * kvec_Mpc,
                      double * muvec,
                      int k_size,
                      int mu_size
                      );

int eft_set_sampling_points(
                      struct eft * peft,
                      double * kvec_Mpc,
                      double * muvec,
                      int k_size,
                      int mu_size
                      );

int eft_set_sampling_points_mu_only(
                      struct eft * peft,
                      double * muvec,
                      int mu_size
                      );

int eft_get_sampling_points(
                      struct eft * peft,
                      double * kvec_Mpc,
                      double * muvec
                      );

int eft_get_sampling_grid_size(
                      struct eft * peft,
                      int * k_size,
                      int * mu_size
                      );

int eft_allocate_spectra_contributions(
                      struct eft * peft,
                      int size,
                      int * moment_list,
                      int moment_list_size
                      );

int eft_compute_spectra_contributions(
                      struct eft * peft,
                      int * moment_list,
                      int moment_list_size
                      );

int eft_load_linear_spectra(
                      struct background * pba,
                      struct fourier * pfo,
                      struct primordial * ppm,
                      struct eft * peft,
                      double z,
                      double f,
                      double D,
                      int * pk_types,
                      int pk_types_size,
                      double * muvec,
                      int mu_size
                      );

int eft_compute_divergences(
                      struct eft * peft,
                      int * moment_list,
                      int moment_list_size
                      );

int eft_build_nonlinear_power_spectrum_wedges(
                      struct eft * peft,
                      struct background * pba,
                      struct primordial * ppm,
                      struct fourier * pfo,
                      short index_pk_out_type,
                      double z,
                      double D_z,
                      double f_z,
                      double * muvec,
                      int mu_size,
                      double As_ratio,
                      struct eft_input_parameters eft_ip,
                      double * pkmu
                      );

double sigma_sq(        \
                      struct eft * peft, 
                      short n, 
                      enum eft_pk_type pk_type
                      );

// int eft_build_nonlinear_power_spectrum_wedges_multiple(
//                       struct eft * peft,
//                       struct background * pba,
//                       struct primordial * ppm,
//                       struct fourier * pfo,
//                       short index_pk_out_type,
//                       double * z,
//                       double * D_z,
//                       double * f_z,
//                       int z_size,
//                       double * muvec,
//                       int mu_size,
//                       struct eft_input_parameters eft_ip,
//                       double ** pkmu
//                       );

/** TODO: add free spectra_contributions */

#endif
