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

int eft_set_sampling_points(
                      struct eft * peft,
                      const double * const kvec_Mpc,
                      const double * const muvec,
                      const int k_size,
                      const int mu_size
                      );

int eft_set_sampling_points_mu_only(
                      struct eft * peft,
                      const double * const muvec,
                      const int mu_size
                      );

int eft_allocate_spectra_contributions(
                      struct eft * peft,
                      const int size,
                      const int * const moment_list,
                      const int moment_list_size
                      );

int eft_compute_spectra_contributions(
                      struct eft * peft,
                      const int * const moment_list,
                      const int moment_list_size
                      );

int eft_load_linear_spectra(
                      struct background * pba,
                      struct fourier * pfo,
                      struct primordial * ppm,
                      struct eft * peft,
                      const double z,
                      const double f,
                      const double D,
                      const int * const pk_types,
                      const int pk_types_size,
                      const double * const muvec,
                      const int mu_size
                      );

int eft_compute_divergences(
                      struct eft * peft,
                      const int * const moment_list,
                      const int moment_list_size
                      );

int eft_build_nonlinear_power_spectrum_wedges(
                      struct eft * peft,
                      struct background * pba,
                      struct primordial * ppm,
                      struct fourier * pfo,
                      const short index_pk_out_type,
                      const double z,
                      const double D_z,
                      const double f_z,
                      const double * muvec,
                      int mu_size,
                      struct eft_input_parameters eft_ip,
                      double * pkmu
                      );

int eft_build_nonlinear_power_spectrum_wedges_multiple(
                      struct eft * peft,
                      struct background * pba,
                      struct primordial * ppm,
                      struct fourier * pfo,
                      const short index_pk_out_type,
                      const double * const z,
                      const double * const D_z,
                      const double * const f_z,
                      const int z_size,
                      const double * muvec,
                      int mu_size,
                      struct eft_input_parameters eft_ip,
                      double ** pkmu
                      );

/** TODO: add free spectra_contributions */

#endif
