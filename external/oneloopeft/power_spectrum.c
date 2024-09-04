/** @file power_spectrum.c
 *
 * author: Christian Radermacher, 2024
 * based on prototype version of Azadeh Moradinezhad Dizgah & Dennis Linde
 *
 * contains non-linear power spectrum definitions along with routines to compute the individual loop terms
*/

#include "power_spectrum.h"

int eft_necessary_spectra_contributions(struct eft * peft,
                                        enum eft_pk_out_type pk_out_type,
                                        int ** list_spectra_contributions,
                                        int * list_spectra_contributions_size) {

  int index_pk_type, index_moment, num_moments, index_list, it = 0;

  switch (pk_out_type)
  {
  case Pdd_mm_real:
    *list_spectra_contributions_size = 2;
    class_alloc(*list_spectra_contributions, (*list_spectra_contributions_size)*sizeof(int), peft->error_message);
    for (index_moment = peft->index_I2200; index_moment <= peft->index_I1300; index_moment++) {
      (*list_spectra_contributions)[index_moment - peft->index_I2200] = pk_ir_resummed_lo * peft->index_num + index_moment;
    }
    break;

  case Pdd_mm_real_no_IR_resum:
    *list_spectra_contributions_size = 2;
    class_alloc(*list_spectra_contributions, (*list_spectra_contributions_size)*sizeof(int), peft->error_message);
    for (index_moment = peft->index_I2200; index_moment <= peft->index_I1300; index_moment++) {
      (*list_spectra_contributions)[index_moment - peft->index_I2200] = pk_lin * peft->index_num + index_moment;
    }
    break;

  case Pdd_mm_22:
    *list_spectra_contributions_size = 1;
    class_alloc(*list_spectra_contributions, (*list_spectra_contributions_size)*sizeof(int), peft->error_message);
    (*list_spectra_contributions)[0] = pk_ir_resummed_lo * peft->index_num + peft->index_I2200;
    break;

  case Pdd_mm_22_no_IR_resum:
    *list_spectra_contributions_size = 1;
    class_alloc(*list_spectra_contributions, (*list_spectra_contributions_size)*sizeof(int), peft->error_message);
    (*list_spectra_contributions)[0] = pk_lin * peft->index_num + peft->index_I2200;
    break;

  case Pdd_mm_13:
    *list_spectra_contributions_size = 1;
    class_alloc(*list_spectra_contributions, (*list_spectra_contributions_size)*sizeof(int), peft->error_message);
    (*list_spectra_contributions)[0] = pk_ir_resummed_lo * peft->index_num + peft->index_I1300;
    break;

  case Pdd_mm_13_no_IR_resum:
    *list_spectra_contributions_size = 1;
    class_alloc(*list_spectra_contributions, (*list_spectra_contributions_size)*sizeof(int), peft->error_message);
    (*list_spectra_contributions)[0] = pk_lin * peft->index_num + peft->index_I1300;
    break;

  case Pdd_mm_rsd:
    /** TODO: */
    *list_spectra_contributions_size = 0;
    class_stop(peft->error_message, "Not implemented yet.");
    break;

  case Pdd_hh_real:
    *list_spectra_contributions_size = 8;
    class_alloc(*list_spectra_contributions, (*list_spectra_contributions_size)*sizeof(int), peft->error_message);
    for (index_moment = peft->index_I2200; index_moment <= peft->index_FG200; index_moment++) {
      (*list_spectra_contributions)[index_moment - peft->index_I2200] = pk_ir_resummed_lo * peft->index_num + index_moment;
    }
    break;

  case Pdd_hh_rsd:
    if (peft->hp->integration_mode == fftlog) {
      num_moments = peft->index_N22z - peft->index_I2200 + 1;
      *list_spectra_contributions_size = 2 * num_moments;
      class_alloc(*list_spectra_contributions, (*list_spectra_contributions_size)*sizeof(int), peft->error_message);
      for (index_pk_type = pk_lin; index_pk_type <= pk_nowiggle; index_pk_type++) {
        for (index_moment = peft->index_I2200; index_moment <= peft->index_N22z; index_moment++) {
          (*list_spectra_contributions)[(index_pk_type - pk_lin)*num_moments + (index_moment - peft->index_I2200)] = index_pk_type * peft->index_num + index_moment;
        }
      }
    }
    else {
      // const int exclusion_list_size = 10;
      // const int exclusion_list[10] = {peft->index_J12102y, peft->index_J21102y, peft->index_Jdelta202y, peft->index_JG202y, peft->index_N11y, peft->index_J21112y, peft->index_J12112y, peft->index_N12y, peft->index_N22y, peft->index_N22z};
      num_moments = peft->index_sigmav_mu - peft->index_I2200 + 1; // - exclusion_list_size;
      *list_spectra_contributions_size = num_moments;
      class_alloc(*list_spectra_contributions, (*list_spectra_contributions_size)*sizeof(int), peft->error_message);
      for (index_moment = peft->index_I2200; (index_moment <= peft->index_sigmav_mu); index_moment++) {
        // for (index_list = 0; (index_list < exclusion_list_size) && (index_moment != exclusion_list[index_list]); index_list++);
        // if (index_list >= exclusion_list_size) {
        (*list_spectra_contributions)[it++] = pkmu_rsd_ir_resummed_lo * peft->index_num + index_moment;
        // }
      }
    }
    break;

  default:
    *list_spectra_contributions_size = 0;
    class_stop(peft->error_message, "Invalid pk_out_type = %d given.", pk_out_type);
    break;
  }

  /** if this list contains an RSD spectrum, the computation of the spectra has to be changed in eft_build_nonlinear_power_spectrum_wedges */

  return _SUCCESS_;
}

int eft_necessary_pk_types_loops(struct eft * peft,
                                 enum eft_pk_out_type pk_out_type,
                                 int ** list_pk_types,
                                 int * list_pk_types_size) {

  switch (pk_out_type)
  {
  case Pdd_mm_real:
  case Pdd_mm_22:
  case Pdd_mm_13:
    *list_pk_types_size = 1;
    class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
    (*list_pk_types)[0] = pk_ir_resummed_lo;
    break;

  case Pdd_mm_real_no_IR_resum:
  case Pdd_mm_22_no_IR_resum:
  case Pdd_mm_13_no_IR_resum:
    *list_pk_types_size = 1;
    class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
    (*list_pk_types)[0] = pk_lin;
    break;

  case Pdd_mm_rsd:
    /** TODO: */
    *list_pk_types_size = 0;
    class_stop(peft->error_message, "Not implemented yet.");
    break;

  case Pdd_hh_real:
    *list_pk_types_size = 1;
    class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
    (*list_pk_types)[0] = pk_ir_resummed_lo;
    break;

  case Pdd_hh_rsd:
    if (peft->hp->integration_mode == fftlog) {
      *list_pk_types_size = 2;
      class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
      (*list_pk_types)[0] = pk_lin;
      (*list_pk_types)[1] = pk_nowiggle;
    }
    else {
      *list_pk_types_size = 1;
      class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
      (*list_pk_types)[0] = pkmu_rsd_ir_resummed_lo;
    }
    break;

  default:
    *list_pk_types_size = 0;
    class_stop(peft->error_message, "Invalid pk_out_type = %d given.", pk_out_type);
    break;
  }

  /** if this list contains an RSD spectrum, the computation of the spectra has to be changed in eft_build_nonlinear_power_spectrum_wedges */

  return _SUCCESS_;
}

int eft_necessary_pk_types_total(struct eft * peft,
                                 const short pk_out_type,
                                 int ** list_pk_types,
                                 int * list_pk_types_size) {

  switch (pk_out_type)
  {
  case Pdd_mm_real:
    *list_pk_types_size = 2;
    class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
    (*list_pk_types)[0] = pk_ir_resummed_lo;
    (*list_pk_types)[1] = pk_ir_resummed_nlo;
    break;

  case Pdd_mm_22:
  case Pdd_mm_13:
    *list_pk_types_size = 1;
    class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
    (*list_pk_types)[0] = pk_ir_resummed_lo;
    break;

  case Pdd_mm_real_no_IR_resum:
  case Pdd_mm_22_no_IR_resum:
  case Pdd_mm_13_no_IR_resum:
    *list_pk_types_size = 1;
    class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
    (*list_pk_types)[0] = pk_lin;
    break;

  case Pdd_mm_rsd:
    /** TODO: */
    *list_pk_types_size = 0;
    class_stop(peft->error_message, "Not implemented yet.");
    break;

  case Pdd_hh_real:
    *list_pk_types_size = 2;
    class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
    (*list_pk_types)[0] = pk_ir_resummed_lo;
    (*list_pk_types)[1] = pk_ir_resummed_nlo;
    break;

  case Pdd_hh_rsd:
    if (peft->hp->integration_mode == fftlog) {
      *list_pk_types_size = 3;
      class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
      (*list_pk_types)[0] = pk_lin;
      (*list_pk_types)[1] = pk_nowiggle;
      (*list_pk_types)[2] = pkmu_rsd_ir_resummed_nlo;
    }
    else {
      *list_pk_types_size = 2;
      class_alloc(*list_pk_types, (*list_pk_types_size)*sizeof(int), peft->error_message);
      (*list_pk_types)[0] = pkmu_rsd_ir_resummed_lo;
      (*list_pk_types)[1] = pkmu_rsd_ir_resummed_nlo;
    }
    break;

  default:
    *list_pk_types_size = 0;
    class_stop(peft->error_message, "Invalid pk_out_type = %d given.", pk_out_type);
    break;
  }

  return _SUCCESS_;
}

int eft_set_sampling_points_all(struct eft * peft0,
                                const int eft_size,
                                const double * const kvec_Mpc,
                                const double * const muvec,
                                const int k_size,
                                const int mu_size) {

    int index_eft;

    for (index_eft = 0; index_eft < eft_size; index_eft++) {
      class_call(eft_set_sampling_points(peft0 + index_eft,
                                         kvec_Mpc,
                                         muvec,
                                         k_size,
                                         mu_size),
                  (peft0 + index_eft)->error_message, (peft0 + index_eft)->error_message);
    }
    
    return _SUCCESS_;
}

int eft_set_sampling_points(struct eft * peft,
                            const double * const kvec_Mpc,
                            const double * const muvec,
                            const int k_size,
                            const int mu_size) {

  int it;

  /** - set the k(mu) and mu arrays where spectra_contributions (and the nonlinear spectrum) shall be evaluated */
  peft->k_size = k_size;
  switch (peft->hp->integration_mode)
    {
    case fftlog:
      peft->mu_size = 1;
      break;
    case direct_integration:
      peft->mu_size = mu_size;
      break;
    }

  /** - allocate (non)linear spectra and wavenumber arrays and copy from input */
  class_realloc(peft->ln_k, peft->ln_k, peft->mu_size*peft->k_size*sizeof(double), peft->error_message);
  for (it = 0; it < peft->mu_size*peft->k_size; it++) {
    peft->ln_k[it] = log( kvec_Mpc[it] );
  }
  if (peft->hp->integration_mode == direct_integration) {
    class_realloc(peft->mu, peft->mu, peft->mu_size*sizeof(double), peft->error_message);
    class_protect_memcpy(peft->mu, muvec, peft->mu_size*sizeof(double));
  }
  else {
    class_realloc(peft->mu, peft->mu, 1*sizeof(double), peft->error_message);
    peft->mu[0] = 0.;
  }

  /** - mark all spectra_contributions for re-computation after changing the sampling points */
  for (it = 0; it < pk_type_num*peft->index_num; it++) {
    peft->spectra_contributions_size[it] = 0;
  }

  return _SUCCESS_;
}

int eft_set_sampling_points_mu_only(struct eft * peft,
                                    const double * const muvec,
                                    const int mu_size) {

  int it;

  /** - set the mu arrays where spectra_contributions (and the nonlinear spectrum) shall be evaluated */
  switch (peft->hp->integration_mode)
    {
    case fftlog:
      peft->mu_size=1;
      break;
    case direct_integration:
      peft->mu_size = mu_size;
      break;
    }

  /** - allocate mu arrays and copy from input */
  if (peft->hp->integration_mode == direct_integration) {
    class_realloc(peft->mu, peft->mu, peft->mu_size*sizeof(double), peft->error_message);
    class_protect_memcpy(peft->mu, muvec, peft->mu_size*sizeof(double));
  }
  else {
    class_realloc(peft->mu, peft->mu, 1*sizeof(double), peft->error_message);
    peft->mu[0] = 0.;
  }

  /** - mark all spectra_contributions for re-computation after changing the sampling points */
  for (it = 0; it < pk_type_num*peft->index_num; it++) {
    peft->spectra_contributions_size[it] = 0;
  }

  return _SUCCESS_;
}

/**
 * @brief Allocates the spectra_contributions array for a list of (pk_type*peft->index_num + index_moment) with a given size;
 *        guarantees that the corresponding kernel matrices and Fourier coefficients are loaded
 * @param peft
 * @param size
 * @param moment_list
 * @param moment_list_size
 *
 * @return the error status
 */
int eft_allocate_spectra_contributions(struct eft * peft,
                                       const int size,
                                       const int * const moment_list,
                                       const int moment_list_size) {

  int index_list, index_pk_type, index_moment, index_part;
  div_t list_elem;

  for (index_list = 0; index_list < moment_list_size; index_list++) {
    list_elem = div(moment_list[index_list], peft->index_num);
    index_pk_type = list_elem.quot;
    index_moment = list_elem.rem;

    class_test(index_pk_type >= pk_type_num, peft->error_message, "Invalid pk_type = %d given.", index_pk_type);
    if (peft->hp->integration_mode == fftlog) {
      class_test(!peft->pk_type_loaded[index_pk_type], peft->error_message, "No Fourier coefficients available for pk_type = %d", index_pk_type);
      class_test((peft->symmetry[index_moment] > no_finite_part) && (peft->loop_matrices_size[index_moment] < 1), peft->error_message, "No kernel matrix available for index_moment = %d", index_moment);
    }

    if (peft->spectra_contributions_size[moment_list[index_list]] != size) {  /** - different size allocated: realloc */
      peft->spectra_contributions_size[moment_list[index_list]] = size;

      for (index_part = 0; index_part < eft_spectra_contribution_num; index_part++) {
        /** - realloc is okay here, since every spectra_contributions[pk_type][index_moment*eft_spectra_contribution_num + index_part] is initialized as NULL */
        class_realloc(peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + index_part],
                      peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + index_part],
                      peft->spectra_contributions_size[moment_list[index_list]]*sizeof(double), peft->error_message);
      }
    } /** - else: correct size allocated (no-op) */
  }

  return _SUCCESS_;
}


int eft_compute_spectra_contributions(struct eft * peft,
                                      const int * const moment_list,
                                      const int moment_list_size) {

  int abort = _FALSE_, index_k, index_list, index_pk_type, index_tracer, index_freq1, index_freq2, index_moment, index_k_loaded = -1, pk_type_loaded = -1;
  div_t list_elem;
  double complex * vec[eft_tracer_num]; /**< Fourier basis elements vec[index_tracer][index_freq1] */
  // double complex * vec_IR;
  register double complex sum1, sum2;

  if (moment_list_size < 1) { return _SUCCESS_; }

  #pragma omp parallel shared(peft, abort, moment_list, moment_list_size), default(none)   \
                       private(vec, index_k, index_tracer, index_freq1, index_freq2, index_pk_type, index_moment, list_elem, sum1, sum2, index_k_loaded, pk_type_loaded)
  {
  for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
    class_alloc_parallel(vec[index_tracer], peft->hp->fourier_coeff_size * sizeof(double complex), peft->error_message);
  }

  index_k_loaded = -1; pk_type_loaded = -1;

  if (!abort) {
    /** - main loop over momenta: parallelized */
    #pragma omp for schedule(static), collapse(2)
    for (index_k = 0; index_k < peft->mu_size*peft->k_size; index_k++) {
      for (index_list = 0; index_list < moment_list_size; index_list++) {
        list_elem = div(moment_list[index_list], peft->index_num);
        index_pk_type = list_elem.quot;
        index_moment = list_elem.rem;

        /** - compute Fourier basis elements */
        if (index_k_loaded != index_k || pk_type_loaded != index_pk_type) {
          index_k_loaded = index_k; pk_type_loaded = index_pk_type;
          for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
            for (index_freq1 = 0; index_freq1 < peft->hp->fourier_coeff_size; index_freq1++) {
              vec[index_tracer][index_freq1] = peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_freq1]    \
                                              * cexp((peft->hp->bias[index_tracer] + _Complex_I * peft->fourier_frequencies[index_tracer][index_freq1]) * peft->ln_k[index_k]);
            }
          }
        }

        /** - multiply with matrices dependent on symmetry type */
        sum1 = 0.;

        switch (peft->symmetry[index_moment])
        {
        case no_finite_part:
          /** - sum1 = 0.; */
          break;

        case sym_vec:
          for (index_freq1 = 0; index_freq1 < peft->hp->fourier_coeff_size; index_freq1++) {
            sum1 += peft->loop_matrices[index_moment][index_freq1] * vec[peft->use_tracer[index_moment]][index_freq1];
          }
          break;

        case sym_mat_none:
          for (index_freq2 = 0; index_freq2 < peft->hp->fourier_coeff_size; index_freq2++) {
            sum2 = 0.;
            for (index_freq1 = 0; index_freq1 < peft->hp->fourier_coeff_size; index_freq1++) {
              sum2 += peft->loop_matrices[index_moment][index_freq1 + index_freq2*peft->hp->fourier_coeff_size]     \
                      * vec[peft->use_tracer[index_moment]][index_freq1];
            }
            sum1 += sum2 * vec[peft->use_tracer[index_moment]][index_freq2];
          }
          break;

        case sym_mat_symmetric:
          /** - off-diagonal */
          for (index_freq2 = 0; index_freq2 < peft->hp->fourier_coeff_size; index_freq2++) {
            sum2 = 0.;
            for (index_freq1 = 0; index_freq1 < index_freq2; index_freq1++) {
              sum2 += peft->loop_matrices[index_moment][index_freq1 + index_freq2*(index_freq2+1)/2]                \
                        * vec[peft->use_tracer[index_moment]][index_freq1];
            }
            sum1 += 2. * sum2 * vec[peft->use_tracer[index_moment]][index_freq2];
          }
          /** - diagonal */
          for (index_freq1 = 0; index_freq1 < peft->hp->fourier_coeff_size; index_freq1++) {
            sum1 += peft->loop_matrices[index_moment][index_freq1*(index_freq1+3)/2]                                \
                      * vec[peft->use_tracer[index_moment]][index_freq1] * vec[peft->use_tracer[index_moment]][index_freq1];
          }
          break;

        default:
          if (!abort) {
            class_protect_sprintf(peft->error_message, "Invalid symmetry type for index_moment = %d at index_k = %d", index_moment, index_k);
          }
          abort = _TRUE_;
          break;
        }

        /** write the sum into the spectra_contributions array and apply the wavenumber dimension, this throws away any residual imaginary part */
        peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + finite_part][index_k]    \
            = exp(peft->spectra_contributions_dimension[index_moment] * peft->ln_k[index_k]) * creal(sum1);
      } /** - end of loop over moment_list */
    } /** - end of k-loop */
  }

  for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
    free(vec[index_tracer]);
  }

  } /** - end of parallel region */

  return abort;
}

int eft_load_linear_spectra(struct background * pba,
                            struct fourier * pfo,
                            struct primordial * ppm,
                            struct eft * peft,
                            const double z,
                            const double f,
                            const double D,
                            const int * const pk_types,
                            const int pk_types_size,
                            const double * const muvec,
                            const int mu_size) {

  int index_list, index_mu, index_pk_type, size, rsd_indicator, abort = _FALSE_;
  double * ln_k;

  /** - allocate arrays for the linear spectra in this list */
  for (index_list = 0; index_list < pk_types_size; index_list++) {
    index_pk_type = pk_types[index_list];
    class_realloc(peft->pk_l[index_pk_type], peft->pk_l[index_pk_type], mu_size*peft->k_size*sizeof(double), peft->error_message);
    class_realloc(peft->ddpk_l[index_pk_type], peft->ddpk_l[index_pk_type], mu_size*peft->k_size*sizeof(double), peft->error_message);
  }

  /** - if we want to use interpolation, repeat the ln_k array mu_size times */
  if (peft->hp->use_interpolation) {
    class_alloc(ln_k, mu_size*peft->k_size*sizeof(double), peft->error_message);
    for (index_mu = 0; index_mu < mu_size; index_mu++) {
      class_protect_memcpy(ln_k + index_mu*peft->k_size, peft->ln_k, peft->k_size*sizeof(double));
    }
  }
  else {
    ln_k = peft->ln_k;
  }

  /** - load the spectra into the arrays */
  #pragma omp parallel for schedule(static, 1), shared(peft, pba, pfo, ppm, z, f, D, pk_types, pk_types_size, ln_k, muvec, mu_size, abort),  \
                           private(index_list, index_mu, index_pk_type, rsd_indicator), default(none)
  for (index_list = 0; index_list < pk_types_size; index_list++) {
    index_pk_type = pk_types[index_list];
    rsd_indicator = eft_rsd_indicator(index_pk_type);
    if (rsd_indicator) {
      class_call_parallel(eft_linear_spectrum_rsd(pba, ppm, pfo, peft,
                                                  linear,
                                                  ln_k,
                                                  peft->k_size,
                                                  muvec,
                                                  mu_size,
                                                  cartesian_product,
                                                  z,
                                                  f,
                                                  D,
                                                  index_pk_type,
                                                  peft->pk_l[index_pk_type]),
                          peft->error_message, peft->error_message);
    }
    else {
      class_call_parallel(eft_linear_spectrum_real(pba, ppm, pfo, peft,
                                                  linear,
                                                  ln_k,
                                                  peft->k_size,
                                                  mu_size,
                                                  z,
                                                  f,
                                                  D,
                                                  index_pk_type,
                                                  peft->pk_l[index_pk_type]),
                          peft->error_message, peft->error_message);
    }

    for (index_mu = 0; index_mu < (rsd_indicator ? mu_size : 1); index_mu++) {
      class_call_parallel(array_spline_table_columns(ln_k + index_mu*peft->k_size,
                                                     peft->k_size,
                                                     peft->pk_l[index_pk_type] + index_mu*peft->k_size,
                                                     1,
                                                     peft->ddpk_l[index_pk_type] + index_mu*peft->k_size,
                                                     _SPLINE_EST_DERIV_,
                                                     peft->error_message),
                          peft->error_message, peft->error_message);
    }
  }

  if (peft->hp->use_interpolation) {
    free(ln_k);
  }

  return _SUCCESS_;
}


double sigma_sq(struct eft * peft, const short n, enum eft_pk_type pk_type) {
  short found = _FALSE_;
  int index_moment;
  double complex value;

  /** - search for the right moment in the storage */
  for (index_moment = 0; (index_moment < EFT_DISPERSION_SIZE) && (peft->dispersion[pk_type][index_moment].index_bias > SHRT_MIN); index_moment++) {
    if (peft->dispersion[pk_type][index_moment].index_bias == n) { found = _TRUE_; break; }
  }

  if (found) {
    return peft->dispersion[pk_type][index_moment].moment;
  }
  else if (index_moment >= EFT_DISPERSION_SIZE)
  {
    class_protect_sprintf(peft->error_message, "%s(L:%d) : not enough storage capacity in dispersion for index_moment = %d", __func__, __LINE__, index_moment);
    return 0.;
  }
  else {
    /** - compute the moment and store at index_moment */
    if (array_integrate_all_spline_table_lines_exponential(
                                  peft->ln_k_moments,
                                  peft->hp->k_size_moments,
                                  peft->pk_l_moments[pk_type],
                                  1,
                                  peft->ddpk_l_moments[pk_type],
                                  class_complex(-(2.*n + 3), 0.),
                                  &value,
                                  peft->error_message) == _FAILURE_) {
      value = 0.;
    }

    peft->dispersion[pk_type][index_moment].moment = value / (2.*_PI_*_PI_);
    peft->dispersion[pk_type][index_moment].index_bias = n;
    return peft->dispersion[pk_type][index_moment].moment;
  }
}

double shot_sym_sq(struct eft * peft, const short n, const short m, enum eft_pk_type pk_type) {
  short found = _FALSE_;
  int index_moment;
  double value;

  /** - search for the right moment in the storage */
  for (index_moment = 0; (index_moment < EFT_UV_CORRECTIONS_SIZE)     \
                      && (peft->ps_uv_shot_noise_corrections_underlying[pk_type][index_moment].index_bias > SHRT_MIN)    \
                      && (peft->ps_uv_shot_noise_corrections_underlying[pk_type][index_moment].index_derivative > SHRT_MIN); index_moment++) {
    if ((peft->ps_uv_shot_noise_corrections_underlying[pk_type][index_moment].index_bias == n)   \
        && (peft->ps_uv_shot_noise_corrections_underlying[pk_type][index_moment].index_derivative == m)) { found = _TRUE_; break; }
  }

  if (found) {
    return peft->ps_uv_shot_noise_corrections_underlying[pk_type][index_moment].moment;
  }
  else if (index_moment >= EFT_UV_CORRECTIONS_UNDERLYING_SIZE)
  {
    class_protect_sprintf(peft->error_message, "%s(L:%d) : not enough storage capacity in dispersion for index_moment = %d", __func__, __LINE__, index_moment);
    return 0.;
  }
  else {
    /** - compute the moment and store at index_moment */
    if (array_square_integrate_exponential_all_spline_table_lines(
                                        peft->ln_k_moments,
                                        peft->hp->k_size_moments,
                                        peft->pk_l_moments[pk_type],
                                        1,
                                        peft->ddpk_l_moments[pk_type],
                                        (double)n + 1.5,
                                        m,
                                        &value,
                                        peft->error_message) == _FAILURE_) {
      value = 0.;
    }

    peft->ps_uv_shot_noise_corrections_underlying[pk_type][index_moment].moment = value / (2.*_PI_*_PI_);
    peft->ps_uv_shot_noise_corrections_underlying[pk_type][index_moment].index_bias = n;
    peft->ps_uv_shot_noise_corrections_underlying[pk_type][index_moment].index_derivative = m;
    return peft->ps_uv_shot_noise_corrections_underlying[pk_type][index_moment].moment;
  }
}

double plin_of_lnk_uv(struct eft * peft, const short derivative_order, enum eft_pk_type pk_type) {
  double value;
  int inf = peft->hp->k_size_moments-2;

  if (array_interpolate_spline_derivative_closeby(
                                peft->ln_k_moments,
                                peft->hp->k_size_moments,
                                peft->pk_l_moments[pk_type],
                                peft->ddpk_l_moments[pk_type],
                                1,
                                log(peft->hp->k_UV_cutoff) - _EPSILON_,
                                derivative_order,
                                &inf,
                                &value,
                                1,
                                peft->error_message) == _FAILURE_) {
    value = 0.;
  }

  return value;
}

double plin_of_lnk_ir(struct eft * peft, const short derivative_order, enum eft_pk_type pk_type) {
  double value;
  int inf = 0;

  if (array_interpolate_spline_derivative_closeby(
                                peft->ln_k_moments,
                                peft->hp->k_size_moments,
                                peft->pk_l_moments[pk_type],
                                peft->ddpk_l_moments[pk_type],
                                1,
                                log(peft->hp->k_IR_cutoff) + _EPSILON_,
                                derivative_order,
                                &inf,
                                &value,
                                1,
                                peft->error_message) == _FAILURE_) {
    value = 0.;
  }

  return value;
}

double shot_sq(struct eft * peft, const short n, const short m, enum eft_pk_type pk_type) {
  short found = _FALSE_;
  int index_moment, code;
  double value = 0., uv_boundary = 0., ir_boundary = 0.;

  /** - search for the right moment in the storage */
  for (index_moment = 0; (index_moment < EFT_UV_CORRECTIONS_SIZE)     \
                      && (peft->ps_uv_shot_noise_corrections[pk_type][index_moment].index_bias > SHRT_MIN)    \
                      && (peft->ps_uv_shot_noise_corrections[pk_type][index_moment].index_derivative > SHRT_MIN); index_moment++) {
    if ((peft->ps_uv_shot_noise_corrections[pk_type][index_moment].index_bias == n)   \
        && (peft->ps_uv_shot_noise_corrections[pk_type][index_moment].index_derivative == m)) { found = _TRUE_; break; }
  }

  if (found) {
    return peft->ps_uv_shot_noise_corrections[pk_type][index_moment].moment;
  }
  else if (index_moment >= EFT_UV_CORRECTIONS_SIZE)
  {
    class_protect_sprintf(peft->error_message, "%s(L:%d) : not enough storage capacity in dispersion for index_moment = %d", __func__, __LINE__, index_moment);
    return 0.;
  }
  else {
    /** - compute the moment and store at index_moment
     *    some depend on each other through partial integration
     *    which leads to recursive calls of this function (max. order of 10 calls) */
    code = (m << 16) | n;
    switch (code)
    {
    /** ------------ bias >= -3/2 ------------ */
    case -2: /** S_0,-2 */
      value = shot_sym_sq(peft, 0, 0, pk_type);
      break;

    /** ------------ bias >= -1/2 ------------ */
    case 0: /** S_0,0 */
      value = shot_sym_sq(peft, -1, 0, pk_type);
      break;

    case 65534: /** S_1,-2 */
      value = -0.5 * shot_sq(peft, 0, 0, pk_type);
      uv_boundary = 0.5 * peft->hp->k_UV_cutoff * pow(plin_of_lnk_uv(peft, 0, pk_type), 2);
      ir_boundary = 0.5 * peft->hp->k_IR_cutoff * pow(plin_of_lnk_ir(peft, 0, pk_type), 2);
      break;

    case 131068:  /** S_2,-4 */
      value = -2. * shot_sq(peft, -2, 1, pk_type) - shot_sym_sq(peft, -1, 1, pk_type);
      uv_boundary = peft->hp->k_UV_cutoff * plin_of_lnk_uv(peft, 0, pk_type) * plin_of_lnk_uv(peft, 1, pk_type);
      ir_boundary = peft->hp->k_IR_cutoff * plin_of_lnk_ir(peft, 0, pk_type) * plin_of_lnk_ir(peft, 1, pk_type);

    /** ------------ bias >= 1/2 ------------ */
    case 2: /** S_0,2 */
      value = shot_sym_sq(peft, -2, 0, pk_type);
      break;

    case 65536: /** S_1,0 */
      value = 0.5 * shot_sq(peft, 2, 0, pk_type);
      uv_boundary = 0.5 * pow(plin_of_lnk_uv(peft, 0, pk_type), 2) / peft->hp->k_UV_cutoff;
      ir_boundary = 0.5 * pow(plin_of_lnk_ir(peft, 0, pk_type), 2) / peft->hp->k_IR_cutoff;
      break;

    case 131070:  /** S_2,-2 */
      value = -shot_sym_sq(peft, -2, 1, pk_type);
      uv_boundary = plin_of_lnk_uv(peft, 0, pk_type) * plin_of_lnk_uv(peft, 1, pk_type) / peft->hp->k_UV_cutoff;
      ir_boundary = plin_of_lnk_ir(peft, 0, pk_type) * plin_of_lnk_ir(peft, 1, pk_type) / peft->hp->k_IR_cutoff;
      break;

    case 196604:  /** S_3,-4 */
      value = -1.5 * shot_sq(peft, -2, 2, pk_type);
      uv_boundary = 0.5 * (plin_of_lnk_uv(peft, 0, pk_type) * (2.*plin_of_lnk_uv(peft, 2, pk_type) - plin_of_lnk_uv(peft, 1, pk_type))    \
                           - pow(plin_of_lnk_uv(peft, 1, pk_type), 2)) / peft->hp->k_UV_cutoff;
      ir_boundary = 0.5 * (plin_of_lnk_ir(peft, 0, pk_type) * (2.*plin_of_lnk_ir(peft, 2, pk_type) - plin_of_lnk_ir(peft, 1, pk_type))    \
                           - pow(plin_of_lnk_ir(peft, 1, pk_type), 2)) / peft->hp->k_IR_cutoff;
      break;

    case 262138:  /** S_4,-6 */
      value = 2. * shot_sq(peft, -2, 2, pk_type) - 2. * shot_sym_sq(peft, -2, 1, pk_type) + shot_sym_sq(peft, -2, 2, pk_type);
      uv_boundary = (plin_of_lnk_uv(peft, 0, pk_type) * (plin_of_lnk_uv(peft, 3, pk_type) - 5.*plin_of_lnk_uv(peft, 2, pk_type) + 4.*plin_of_lnk_uv(peft, 1, pk_type))    \
                    - plin_of_lnk_uv(peft, 1, pk_type) * (plin_of_lnk_uv(peft, 2, pk_type) - 2.*plin_of_lnk_uv(peft, 1, pk_type))) / peft->hp->k_UV_cutoff;
      ir_boundary = (plin_of_lnk_ir(peft, 0, pk_type) * (plin_of_lnk_ir(peft, 3, pk_type) - 5.*plin_of_lnk_ir(peft, 2, pk_type) + 4.*plin_of_lnk_ir(peft, 1, pk_type))    \
                    - plin_of_lnk_ir(peft, 1, pk_type) * (plin_of_lnk_ir(peft, 2, pk_type) - 2.*plin_of_lnk_ir(peft, 1, pk_type))) / peft->hp->k_IR_cutoff;
      break;

    /** for bias >~ 1 the method of extracting the divergence breaks down (exactly at bias > n_s) */

    default:
      break;
    }

    peft->ps_uv_shot_noise_corrections[pk_type][index_moment].moment = value + (uv_boundary - ir_boundary)/(2.*_PI_*_PI_);
    peft->ps_uv_shot_noise_corrections[pk_type][index_moment].index_bias = n;
    peft->ps_uv_shot_noise_corrections[pk_type][index_moment].index_derivative = m;
    return peft->ps_uv_shot_noise_corrections[pk_type][index_moment].moment;
  }
}

// double plin_of_lnk_out(struct eft * peft, double ln_k, int * last_index, const short derivative_order, enum eft_pk_type pk_type) {
//   double value;

//   if (array_interpolate_spline_derivative_closeby(
//                                 peft->ln_k,
//                                 peft->k_size,
//                                 peft->pk_l[pk_type],
//                                 peft->ddpk_l[pk_type],
//                                 1,
//                                 ln_k,
//                                 derivative_order,
//                                 last_index,
//                                 &value,
//                                 1,
//                                 peft->error_message) == _FAILURE_) {
//     value = 0.;
//   }

//   return value;
// }

int add_plin_contribution(struct eft * peft,
                          const int k_power,
                          const int mu_power,
                          const short derivative_order,
                          const double prefactor,
                          double * out_pkmu0,
                          const int pk_type_and_moment) {

  int it_k, it_mu, pk_type, mu_size, last_index = 0;
  double plin;

  mu_size = peft->spectra_contributions_size[pk_type_and_moment] / peft->k_size;
  pk_type = pk_type_and_moment / peft->index_num;

  for (it_mu = 0; it_mu < mu_size; it_mu++) {
    for (it_k = 0; it_k < peft->k_size; it_k++) {
      /** - get the (derivative of) linear spectrum at ln_k (evaluate spline at the sample points) */
      if (array_interpolate_spline_derivative_closeby(
                                peft->ln_k + it_mu*peft->k_size,
                                peft->k_size,
                                peft->pk_l[pk_type] + it_mu*peft->k_size,
                                peft->ddpk_l[pk_type] + it_mu*peft->k_size,
                                1,
                                peft->ln_k[it_mu*peft->k_size + it_k],
                                derivative_order,
                                &last_index,
                                &plin,
                                1,
                                peft->error_message) == _FAILURE_) {
        plin = 0.;
      }
      /** - prefactor * k^(k_power) * P^(m)_lin(ln(k)) */
      out_pkmu0[it_mu*peft->k_size + it_k] = prefactor * exp(k_power * peft->ln_k[it_mu*peft->k_size + it_k]) * plin;
    }
  }

  /** - multiply by the mu-factor if the mu-dependencies are not added later */
  if (peft->hp->integration_mode == direct_integration) {
    for (it_mu = 0; it_mu < peft->mu_size; it_mu++) {
      for (it_k = 0; it_k < peft->k_size; it_k++) {
        /** - prefactor * k^(k_power) * mu^(mu_power) * P^(m)_lin(ln(k)) */
        out_pkmu0[it_mu*peft->k_size + it_k] *= pow(peft->mu[it_mu], mu_power);
      }
    }
  }

  return _SUCCESS_;
}

int add_shot_contribution(struct eft * peft,
                          const int k_power,
                          const int mu_power,
                          const double prefactor,
                          double * out_pkmu0,
                          const int pk_type_and_moment) {

  int it_k, it_mu, mu_size;

  mu_size = peft->spectra_contributions_size[pk_type_and_moment] / peft->k_size;

  for (it_mu = 0; it_mu < mu_size; it_mu++) {
    for (it_k = 0; it_k < peft->k_size; it_k++) {
      /** - prefactor * k^(k_power) */
      out_pkmu0[it_mu*peft->k_size + it_k] = prefactor * exp(k_power * peft->ln_k[it_mu*peft->k_size + it_k]);
    }
  }

  /** - multiply by the mu-factor if the mu-dependencies are not added later */
  if (peft->hp->integration_mode == direct_integration) {
    for (it_mu = 0; it_mu < mu_size; it_mu++) {
      for (it_k = 0; it_k < peft->k_size; it_k++) {
        /** - prefactor * k^(k_power) * mu^(mu_power) */
        out_pkmu0[it_mu*peft->k_size + it_k] *= pow(peft->mu[it_mu], mu_power);
      }
    }
  }

  return _SUCCESS_;
}


int eft_compute_divergences(struct eft * peft,
                            const int * const moment_list,
                            const int moment_list_size) {

  int index_pk_type, index_list, index_moment, index_part, it, max_moments = 0;
  div_t list_elem;
  double bias;
  short abort = _FALSE_;

  if (moment_list_size < 1) { return _SUCCESS_; }

  #pragma omp parallel for shared(peft, moment_list, moment_list_size, abort), private(index_pk_type, index_moment, index_list, list_elem, index_part, it, bias), default(none),  \
                           schedule(static), collapse(2)
  for (index_list = 0; index_list < moment_list_size; index_list++) {
    for (index_part = 1; index_part < eft_spectra_contribution_num; index_part++) {
      list_elem = div(moment_list[index_list], peft->index_num);
      index_pk_type = list_elem.quot;
      index_moment = list_elem.rem;
      if (abort) { continue; }
      /** - initialize the spectra contribution array */
      for (it = 0; it < peft->mu_size*peft->k_size; it++) {
        peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + index_part][it] = 0.;
      }
      bias = peft->hp->bias[peft->use_tracer[index_moment]];

      switch (index_part)
      {
      /** -------------------------- UV divergences (q->oo) -------------------------- */
      case uv_divergence:
        if (index_moment == peft->index_I2200) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 4, 0, 9./196.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_I1300) { // * P_lin(k)
          if (bias >= -1.) {
            add_shot_contribution(peft, 2, 0, -61./1890.*sigma_sq(peft, -1, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Idelta200) {
          if (bias >= -0.5) {
            add_shot_contribution(peft, 2, 0, -1./42.*shot_sq(peft, 0, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          if (bias >= 0.5) {
            add_shot_contribution(peft, 4, 0, (2.*shot_sq(peft, 2, 0, index_pk_type) - shot_sq(peft, -2, 2, index_pk_type))/28.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_IG200) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 4, 0, -1./21.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Idelta2delta200) {
          if (bias < -1.5) { /** TODO: decide on shot noise removal */
            add_shot_contribution(peft, 0, 0, -shot_sq(peft, -2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          else if (bias >= -1.5) {
            /* add_shot_contribution(peft, 0, 0, shot_sq(peft, -2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]); */
            /** - cancel with the contribution before */
          }
          if (bias >= -0.5) {
            add_shot_contribution(peft, 2, 0, (2.*shot_sq(peft, -2, 1, index_pk_type) + shot_sq(peft, -4, 2, index_pk_type))/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          if (bias >= 0.5) {
            add_shot_contribution(peft, 4, 0, (4.*shot_sq(peft, -4, 3, index_pk_type) + shot_sq(peft, -6, 4, index_pk_type))/120.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_IG2G200) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 4, 0, 8./15.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Idelta2G200) {
          if (bias >= -0.5) {
            add_shot_contribution(peft, 2, 0, -2./3.*shot_sq(peft, 0, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          if (bias >= 0.5) {
            add_shot_contribution(peft, 4, 0, (2.*shot_sq(peft, 2, 0, index_pk_type) - shot_sq(peft, -2, 2, index_pk_type))/15.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_FG200) { // * P_lin(k)
          if (bias >= -1.) {
            add_shot_contribution(peft, 2, 0, -8./21.*sigma_sq(peft, -1, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        /** ------------ 1-st moment ------------ */
        else if (index_moment == peft->index_I2201) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 4, 0, 19./588.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_I1301p3101) {  // * P_lin(k)
          if (bias >= -1.) {
            add_shot_contribution(peft, 2, 0, -25./189.*sigma_sq(peft, -1, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Idelta201) {
          if (bias >= -0.5) {
            add_shot_contribution(peft, 2, 0, -3./14.*shot_sq(peft, 0, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          if (bias >= 0.5) {
            add_shot_contribution(peft, 4, 0, 23./420.*(2.*shot_sq(peft, 2, 0, index_pk_type) - shot_sq(peft, -2, 2, index_pk_type)),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_IG201) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 4, 0, 11./105.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_FG201) { // * P_lin(k)
          if (bias >= -1.) {
            add_shot_contribution(peft, 2, 0, -8./21.*sigma_sq(peft, -1, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J12101) {  // * P_lin(k) mu
          if (bias >= -3.) {
            add_shot_contribution(peft, -1, 1, sigma_sq(peft, 0, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          if (bias >= -1.) {
            add_shot_contribution(peft, 1, 1, 6./35.*sigma_sq(peft, -1, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J11201) {  // * P_lin(k) mu
          if (bias >= -3.) {
            add_shot_contribution(peft, -1, 1, -sigma_sq(peft, 0, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          if (bias >= -1.) {
            add_shot_contribution(peft, 1, 1, -sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J21101) {  // * mu
          if (bias >= 0.5) {
            add_shot_contribution(peft, 3, 1, (shot_sq(peft, 2, 0, index_pk_type) + 3.*shot_sq(peft, 0, 1, index_pk_type))/42.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Jdelta201) { // * mu
          if (bias >= -0.5) {
            add_shot_contribution(peft, 1, 1, -shot_sq(peft, -2, 1, index_pk_type)/3.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          if (bias >= 0.5) {
            add_shot_contribution(peft, 3, 1, (2.*shot_sq(peft, 0, 1, index_pk_type) - 2.*shot_sq(peft, -2, 2, index_pk_type) - shot_sq(peft, -4, 3, index_pk_type))/30.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_JG201) { // * mu
          if (bias >= 0.5) {
            add_shot_contribution(peft, 3, 1, 2./15.*(-2.*shot_sq(peft, 2, 0, index_pk_type) + shot_sq(peft, 0, 1, index_pk_type)),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        /** ------------ 2-nd moment ------------ */
        else if (index_moment == peft->index_J12102x) { // * P_lin(k)
          if (bias >= -1.) {
            add_shot_contribution(peft, 0, 0, -4./35.*sigma_sq(peft, -1, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J12102y) { // * P_lin(k) mu^2
          if (bias >= -1.) {
            add_shot_contribution(peft, 0, 2, -23./210.*sigma_sq(peft, -1, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J21102x) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 2, 0, -shot_sq(peft, 2, 0, index_pk_type)/42.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J21102y) { // * mu^2
          if (bias >= 0.5) {
            add_shot_contribution(peft, 2, 2, 2./21.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Jdelta202x) {
          if (bias >= -0.5) {
            add_shot_contribution(peft, 0, 0, -1./3.*shot_sq(peft, 0, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          if (bias >= 0.5) {
            add_shot_contribution(peft, 2, 0, (2.*shot_sq(peft, 2, 0, index_pk_type) - shot_sq(peft, -2, 2, index_pk_type))/30.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Jdelta202y) {  // * mu^2
          if (bias >= 0.5) {
            add_shot_contribution(peft, 2, 2, (2.*shot_sq(peft, 2, 0, index_pk_type) - shot_sq(peft, -2, 2, index_pk_type))/15.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_JG202x) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 2, 0, 4./15.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_JG202y) {  // * mu^2
          if (bias >= 0.5) {
            add_shot_contribution(peft, 2, 2, -2./15.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_I2211) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 4, 0, 61./980.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_I1311) { // * P_lin(k)
          if (bias >= -1.) {
            add_shot_contribution(peft, 2, 0, -sigma_sq(peft, -1, index_pk_type)/10.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J12111) {  // * P_lin(k) mu
          if (bias >= -3.) {
            add_shot_contribution(peft, -1, 1, sigma_sq(peft, 0, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          if (bias >= -1.) {
            add_shot_contribution(peft, 1, 1, 6./35.*sigma_sq(peft, -1, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J11211) {  // * P_lin(k) mu
          if (bias >= -3.) {
            add_shot_contribution(peft, -1, 1, -sigma_sq(peft, 0, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
          if (bias >= -1.) {
            add_shot_contribution(peft, 1, 1, -sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J21111) {  // * mu
          if (bias >= 0.5) {
            add_shot_contribution(peft, 3, 1, (-11.*shot_sq(peft, 2, 0, index_pk_type) + 23.*shot_sq(peft, 0, 1, index_pk_type))/210.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N11x) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 2, 0, (shot_sq(peft, 2, 0, index_pk_type) + 2.*shot_sq(peft, 0, 1, index_pk_type))/15.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N11y) {  // * mu^2
          if (bias >= 0.5) {
            add_shot_contribution(peft, 2, 2, (2.*shot_sq(peft, 2, 0, index_pk_type) - shot_sq(peft, 0, 1, index_pk_type))/15.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        /** ------------ 3-rd moment ------------ */
        else if (index_moment == peft->index_J21112x) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 2, 0, 11./210.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J21112y) { // * mu^2
          if (bias >= 0.5) {
            add_shot_contribution(peft, 2, 2, 2./35.*shot_sq(peft, 2, 0, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J12112x) { // * P_lin(k)
          if (bias >= -1.) {
            add_shot_contribution(peft, 0, 0, -4./35.*sigma_sq(peft, -1, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J12112y) { // * P_lin(k) mu^2
          if (bias >= -1.) {
            add_shot_contribution(peft, 0, 2, -23./210.*sigma_sq(peft, -1, index_pk_type),
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N12x) {  // * mu
          if (bias >= 0.5) {
            add_shot_contribution(peft, 1, 1, (-shot_sq(peft, 2, 0, index_pk_type) + 3.*shot_sq(peft, 0, 1, index_pk_type))/15.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N12y) {  // * mu^3
          /** - nothing for bias < 1 */
        }
        /** ------------ 4-th moment ------------ */
        else if (index_moment == peft->index_N22x) {
          if (bias >= 0.5) {
            add_shot_contribution(peft, 0, 0, shot_sq(peft, 2, 0, index_pk_type)/5.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N22y) {  // * mu^2
          /** - nothing for bias < 1 */
        }
        else if (index_moment == peft->index_N22z) {  // * mu^4
          /** - nothing for bias < 1 */
        }
        else {
          abort = _TRUE_;
          ErrorMsg errmsg;
          class_protect_sprintf(errmsg, "index_moment = %d is out of range for uv_divergence part.", index_moment);
          class_build_error_string(peft->error_message, "error; %s", errmsg);
        }
        break;

      /** -------------------------- IR divergences (q->0) --------------------------- */
      case ir_divergence:
        if (index_moment == peft->index_I2200) {
          if (bias <= -1.) {
            add_plin_contribution(peft, 2, 0, 0, sigma_sq(peft, -1, index_pk_type)/12.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_I1300) { // * P_lin(k)
          if (bias <= -1.) {
            add_shot_contribution(peft, 2, 0, -sigma_sq(peft, -1, index_pk_type)/18.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Idelta200) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_IG200) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_Idelta2delta200) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_IG2G200) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_Idelta2G200) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_FG200) { // * P_lin(k)
          /** - nothing for bias > -3 */
        }
        /** ------------ 1-st moment ------------ */
        else if (index_moment == peft->index_I2201) {
          if (bias <= -1.) {
            add_plin_contribution(peft, 2, 0, 0, sigma_sq(peft, -1, index_pk_type)/12.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_I1301p3101) {  // * P_lin(k)
          if (bias <= -1.) {
            add_shot_contribution(peft, 2, 0, -sigma_sq(peft, -1, index_pk_type)/9.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Idelta201) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_IG201) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_FG201) { // * P_lin(k)
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_J12101) {  // * P_lin(k) mu
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_J11201) {  // * P_lin(k) mu
          if (bias < -1.) {
            add_shot_contribution(peft, 1, 1, -sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J21101) {  // * mu
          if (bias <= -1.) {
            add_plin_contribution(peft, 1, 1, 0, sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Jdelta201) { // * mu
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_JG201) { // * mu
          /** - nothing for bias > -3 */
        }
        /** ------------ 2-nd moment ------------ */
        else if (index_moment == peft->index_J12102x) { // * P_lin(k)
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_J12102y) { // * P_lin(k) mu^2
          if (bias <= -1.) {
            add_shot_contribution(peft, 0, 2, -sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J21102x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_J21102y) { // * mu^2
          if (bias <= -1.) {
            add_plin_contribution(peft, 0, 2, 0, sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Jdelta202x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_Jdelta202y) {  // * mu^2
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_JG202x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_JG202y) {  // * mu^2
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_I2211) {
          if (bias <= -1.) {
            add_plin_contribution(peft, 2, 0, 0, sigma_sq(peft, -1, index_pk_type)/12.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_I1311) { // * P_lin(k)
          if (bias < -1.) {
            add_shot_contribution(peft, 2, 0, -sigma_sq(peft, -1, index_pk_type)/18.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J12111) {  // * P_lin(k) mu
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_J11211) {  // * P_lin(k) mu
          if (bias < -1.) {
            add_shot_contribution(peft, 1, 1, -sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J21111) {  // * mu
          if (bias <= -1.) {
            add_plin_contribution(peft, 1, 1, 0, sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N11x) {
          if (bias <= -1.) {
            add_plin_contribution(peft, 0, 0, 0, sigma_sq(peft, -1, index_pk_type)/3.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N11y) {  // * mu^2
          /** - nothing for bias > -3 */
        }
        /** ------------ 3-rd moment ------------ */
        else if (index_moment == peft->index_J21112x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_J21112y) { // * mu^2
          if (bias <= -1.) {
            add_plin_contribution(peft, 0, 2, 0, sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J12112x) { // * P_lin(k)
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_J12112y) { // * P_lin(k) mu^2
          if (bias <= -1.) {
            add_shot_contribution(peft, 0, 2, -sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N12x) {  // * mu
          if (bias <= -1.) {
            add_plin_contribution(peft, -1, 1, 0, sigma_sq(peft, -1, index_pk_type)/3.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N12y) {  // * mu^3
          /** - nothing for bias > -3 */
        }
        /** ------------ 4-th moment ------------ */
        else if (index_moment == peft->index_N22x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_N22y) {  // * mu^2
          if (bias <= -1.) {
            add_plin_contribution(peft, -2, 2, 0, sigma_sq(peft, -1, index_pk_type)/3.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N22z) {  // * mu^4
          /** - nothing for bias > -3 */
        }
        else {
          abort = _TRUE_;
          ErrorMsg errmsg;
          class_protect_sprintf(errmsg, "index_moment = %d is out of range for ir_divergence part.", index_moment);
          class_build_error_string(peft->error_message, "error; %s", errmsg);
        }
        break;

      /** ------------------------- Pole divergences (q->k) -------------------------- */
      case pole_divergence:
        if (index_moment == peft->index_I2200) {
          if (bias <= -1.) {
            add_plin_contribution(peft, 2, 0, 0, sigma_sq(peft, -1, index_pk_type)/12.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + pole_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_I1300) { // * P_lin(k)
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_Idelta200) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_IG200) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_Idelta2delta200) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_IG2G200) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_Idelta2G200) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_FG200) { // * P_lin(k)
          // TODO
        }
        /** ------------ 1-st moment ------------ */
        else if (index_moment == peft->index_I2201) {
          if (bias <= -1.) {
            add_plin_contribution(peft, 2, 0, 0, sigma_sq(peft, -1, index_pk_type)/12.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + pole_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_I1301p3101) {  // * P_lin(k)
          // TODO
        }
        else if (index_moment == peft->index_Idelta201) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_IG201) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_FG201) { // * P_lin(k)
          // TODO
        }
        else if (index_moment == peft->index_J12101) {  // * P_lin(k) mu
          // TODO
        }
        else if (index_moment == peft->index_J11201) {  // * P_lin(k) mu
          // TODO
        }
        else if (index_moment == peft->index_J21101) {  // * mu
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_Jdelta201) { // * mu
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_JG201) { // * mu
          /** - nothing for bias > -3 */
        }
        /** ------------ 2-nd moment ------------ */
        else if (index_moment == peft->index_J12102x) { // * P_lin(k)
          // TODO
        }
        else if (index_moment == peft->index_J12102y) { // * P_lin(k) mu^2
          // TODO
        }
        else if (index_moment == peft->index_J21102x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_J21102y) { // * mu^2
          if (bias <= -1.) {
            add_plin_contribution(peft, 0, 2, 0, sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + pole_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_Jdelta202x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_Jdelta202y) {  // * mu^2
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_JG202x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_JG202y) {  // * mu^2
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_I2211) {
          if (bias <= -1.) {
            add_plin_contribution(peft, 2, 0, 0, sigma_sq(peft, -1, index_pk_type)/12.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + pole_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_I1311) { // * P_lin(k)
          // TODO
        }
        else if (index_moment == peft->index_J12111) {  // * P_lin(k) mu
          // TODO
        }
        else if (index_moment == peft->index_J11211) {  // * P_lin(k) mu
          // TODO
        }
        else if (index_moment == peft->index_J21111) {  // * mu
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_N11x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_N11y) {  // * mu^2
          /** - nothing for bias > -3 */
        }
        /** ------------ 3-rd moment ------------ */
        else if (index_moment == peft->index_J21112x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_J21112y) { // * mu^2
          if (bias <= -1.) {
            add_plin_contribution(peft, 0, 2, 0, sigma_sq(peft, -1, index_pk_type)/6.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + pole_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_J12112x) { // * P_lin(k)
          // TODO
        }
        else if (index_moment == peft->index_J12112y) { // * P_lin(k) mu^2
          // TODO
        }
        else if (index_moment == peft->index_N12x) {  // * mu
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_N12y) {  // * mu^3
          /** - nothing for bias > -3 */
        }
        /** ------------ 4-th moment ------------ */
        else if (index_moment == peft->index_N22x) {
          /** - nothing for bias > -3 */
        }
        else if (index_moment == peft->index_N22y) {  // * mu^2
          if (bias <= -1.) {
            add_plin_contribution(peft, -2, 2, 0, sigma_sq(peft, -1, index_pk_type)/3.,
                                  peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + pole_divergence],
                                  moment_list[index_list]);
          }
        }
        else if (index_moment == peft->index_N22z) {  // * mu^4
          /** - nothing for bias > -3 */
        }
        else {
          abort = _TRUE_;
          ErrorMsg errmsg;
          class_protect_sprintf(errmsg, "index_moment = %d is out of range for pole_divergence part.", index_moment);
          class_build_error_string(peft->error_message, "error; %s", errmsg);
        }
        break;

      default:
        abort = _TRUE_;
        ErrorMsg errmsg;
        class_protect_sprintf(errmsg, "divergence identifier = %d not recognized.", index_part);
        class_build_error_string(peft->error_message, "error; %s", errmsg);
        break;
      }
    }
  }


  return _SUCCESS_;
}

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
                                  double * pkmu) {  /**< pkmu[index_mu*peft->k_size + index_k] */

  int index_pk_type, index_list, index_moment, index_part, index_z, index_mu, index_mu_k, index_k, it, size;
  int * pk_types_loops, pk_types_loops_size = 0, * pk_types, pk_types_size = 0, * pk_types_outside_loops, pk_types_outside_loops_size = 0;
  double * pkmu_loop[pk_type_num*eft_spectra_contribution_num]; /**< nonlinear power spectrum still with separated finite and divergent parts */
  double * mu_real = NULL;
  short abort = _FALSE_;
  double  k, mu, D2, D4, Plin, I2200, I1300, Idelta200, IG200, Idelta2delta200, IG2G200, Idelta2G200, FG200,  \
          I2201, Idelta201, IG201, J21101, Jdelta201, JG201, FG201, I1301p3101, J12101, J11201,       \
          J21102, Jdelta202, JG202, I2211, J21111, N11, J12102, I1311, J12111, J11211,                \
          J21112, N12, J12112,                                                                        \
          N22, sigmav_mu,                                                                             \
          Preal_loop, Prsd0_loop, Prsd1_loop, Prsd2_loop, Prsd3_loop, Prsd4_loop, sum, sigma2_tot_z0;  /** RSD moments multiplied by i^n/n! (k mu)^n only containing loop terms */
  /** sigma2_ir_at_z = D2 * sigma2_ir */

  if (!peft->hp->use_interpolation) {
    muvec = peft->mu; mu_size = peft->mu_size;
    /** TODO: implement */
    //class_stop(peft->error_message, "Not implemented yet.")
  }

  /** - if a real space spectrum is requested, ignore the mu input and let muvec point to a valid address with a constant value */
  if (!eft_rsd_out_indicator(index_pk_out_type)) {
    class_alloc(mu_real, mu_size*sizeof(double), peft->error_message);
    for (index_mu = 0; index_mu < mu_size; index_mu++) { mu_real[it] = 0.; }
    muvec = mu_real;
  }

  /** - get a list of necessary pk_types for the chosen output spectrum of which loop corrections shall be calculated */
  class_call(eft_necessary_pk_types_loops(peft, index_pk_out_type, &pk_types_loops, &pk_types_loops_size),
             peft->error_message, peft->error_message);
  /** - get a list of all necessary pk_types for the chosen output spectrum */
  class_call(eft_necessary_pk_types_total(peft, index_pk_out_type, &pk_types, &pk_types_size),
             peft->error_message, peft->error_message);

  /** - find the spectra in this list that have not yet been loaded in FFTLog mode */
  pk_types_outside_loops_size = pk_types_size - pk_types_loops_size;
  class_alloc(pk_types_outside_loops, pk_types_outside_loops_size*sizeof(int), peft->error_message);
  it = 0;
  for (index_list = 0; index_list < pk_types_size; index_list++) {
    for (index_part = 0; index_part < pk_types_loops_size && pk_types[index_list] != pk_types_loops[index_part]; index_part++);
    if (index_part == pk_types_loops_size) {
      pk_types_outside_loops[it++] = pk_types[index_list];
    }
  }
  /** - and load them into peft->pk_l (will be rescaled by D2 later) */
  class_call(eft_load_linear_spectra(pba, pfo, ppm,
                                     peft,
                                     z,
                                     f_z,
                                     D_z,
                                     (peft->hp->integration_mode == fftlog) ? pk_types_outside_loops : pk_types,
                                     (peft->hp->integration_mode == fftlog) ? pk_types_outside_loops_size : pk_types_size,
                                     muvec,
                                     mu_size),
              peft->error_message, peft->error_message);

  /** - allocate arrays for the loop corrections */
  for (index_list = 0; index_list < pk_types_loops_size; index_list++) {
    for (index_part = 0; index_part < eft_spectra_contribution_num; index_part++) {
      index_pk_type = pk_types_loops[index_list];
      class_alloc(pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part], mu_size*peft->k_size*sizeof(double), peft->error_message);
    }
  }

  /** - growth function rescale factors w.r.t. peft->z0 for tree and 1-loop contributions */
  D2 = pow(D_z/peft->D_z0, 2.);
  D4 = pow(D_z/peft->D_z0, 4.);

  #pragma omp parallel shared(peft, f_z, D2, D4, pk_types_loops, pk_types_loops_size, muvec, mu_size, index_pk_out_type, eft_ip, pkmu_loop, pkmu),   \
                       private(index_list, index_part, index_k, index_mu, index_mu_k, index_pk_type,  \
                               k, mu, Plin, I2200, I1300, Idelta200, IG200, Idelta2delta200, IG2G200, Idelta2G200, FG200,  \
                               I2201, Idelta201, IG201, J21101, Jdelta201, JG201, FG201, I1301p3101, J12101, J11201,       \
                               J21102, Jdelta202, JG202, I2211, J21111, N11, J12102, I1311, J12111, J11211,                \
                               J21112, N12, J12112,                                                                        \
                               N22, sigmav_mu,                                                                             \
                               Preal_loop, Prsd0_loop, Prsd1_loop, Prsd2_loop, Prsd3_loop, Prsd4_loop, sum, sigma2_tot_z0), default(none)
  {
  /** - compute the loop power spectra of every internal pk_type in the list for every k, mu, still separated into finite, UV, IR and pole parts */
  #pragma omp for schedule(static), collapse(4)
  for (index_list = 0; index_list < pk_types_loops_size; index_list++) {
    for (index_part = 0; index_part < eft_spectra_contribution_num; index_part++) {
      for (index_k = 0; index_k < peft->k_size; index_k++) {
        for (index_mu = 0; index_mu < mu_size; index_mu++) {
          index_pk_type = pk_types_loops[index_list];

          if (peft->hp->use_interpolation) { index_mu_k = index_k; }
          else { index_mu_k = index_mu*peft->k_size + index_k; }

          k = exp(peft->ln_k[index_mu_k]);
          mu = muvec[index_mu];
          Plin = peft->pk_l[index_pk_type][index_mu_k];

          switch (index_pk_out_type)
          {
          case Pdd_mm_real:
          case Pdd_mm_real_no_IR_resum:
            I2200           =        peft->spectra_contributions[index_pk_type][peft->index_I2200*eft_spectra_contribution_num + index_part][index_mu_k];
            I1300           = Plin * peft->spectra_contributions[index_pk_type][peft->index_I1300*eft_spectra_contribution_num + index_part][index_mu_k];

            Preal_loop = D4 * 2.*(I2200 + 3.*I1300);
            if (index_part == uv_divergence) {
              Preal_loop += -D2 * 2.*k*k * eft_ip.cs2 * Plin;
            }

            pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part][index_mu*peft->k_size + index_k] = Preal_loop;
            break;

          case Pdd_mm_22:
          case Pdd_mm_22_no_IR_resum:
            I2200           =        peft->spectra_contributions[index_pk_type][peft->index_I2200*eft_spectra_contribution_num + index_part][index_mu_k];
            
            Preal_loop = D4 * I2200;

            pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part][index_mu*peft->k_size + index_k] = Preal_loop;
            break;

          case Pdd_mm_13:
          case Pdd_mm_13_no_IR_resum:
            I1300           = Plin * peft->spectra_contributions[index_pk_type][peft->index_I1300*eft_spectra_contribution_num + index_part][index_mu_k];

            Preal_loop = D4 * I1300;

            pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part][index_mu*peft->k_size + index_k] = Preal_loop;
            break;

          case Pdd_mm_rsd:
            // TODO
            break;

          case Pdd_hh_real:
            I2200           =        peft->spectra_contributions[index_pk_type][peft->index_I2200*eft_spectra_contribution_num + index_part][index_mu_k];
            I1300           = Plin * peft->spectra_contributions[index_pk_type][peft->index_I1300*eft_spectra_contribution_num + index_part][index_mu_k];
            Idelta200       =        peft->spectra_contributions[index_pk_type][peft->index_Idelta200*eft_spectra_contribution_num + index_part][index_mu_k];
            IG200           =        peft->spectra_contributions[index_pk_type][peft->index_IG200*eft_spectra_contribution_num + index_part][index_mu_k];
            Idelta2delta200 =        peft->spectra_contributions[index_pk_type][peft->index_Idelta2delta200*eft_spectra_contribution_num + index_part][index_mu_k];
            IG2G200         =        peft->spectra_contributions[index_pk_type][peft->index_IG2G200*eft_spectra_contribution_num + index_part][index_mu_k];
            Idelta2G200     =        peft->spectra_contributions[index_pk_type][peft->index_Idelta2G200*eft_spectra_contribution_num + index_part][index_mu_k];
            FG200           = Plin * peft->spectra_contributions[index_pk_type][peft->index_FG200*eft_spectra_contribution_num + index_part][index_mu_k];

            Prsd0_loop = D4 * (2.*eft_ip.b1*eft_ip.b1 * (I2200 + 3.*I1300) + 2.*eft_ip.b1*eft_ip.b2 * Idelta200 + 2.*eft_ip.bG2*eft_ip.bG2 * IG2G200 + 4.*eft_ip.b1*eft_ip.bG2 * IG200        \
                              + 0.5*eft_ip.b2*eft_ip.b2 * Idelta2delta200 + 2.*eft_ip.b2*eft_ip.bG2 * Idelta2G200 + 8.*eft_ip.b1*(eft_ip.bG2 + 0.4*eft_ip.btd) * FG200);
            if (index_part == uv_divergence) {
              Prsd0_loop += -D2 * 2.*k*k * (eft_ip.b1 * (eft_ip.R2 + eft_ip.b1 * eft_ip.cs2)) * Plin;
            }

            pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part][index_mu*peft->k_size + index_k] = Prsd0_loop;
            break;

          case Pdd_hh_rsd:
            /** 0-th RSD moment */
            I2200           =        peft->spectra_contributions[index_pk_type][peft->index_I2200*eft_spectra_contribution_num + index_part][index_mu_k];
            I1300           = Plin * peft->spectra_contributions[index_pk_type][peft->index_I1300*eft_spectra_contribution_num + index_part][index_mu_k];
            Idelta200       =        peft->spectra_contributions[index_pk_type][peft->index_Idelta200*eft_spectra_contribution_num + index_part][index_mu_k];
            IG200           =        peft->spectra_contributions[index_pk_type][peft->index_IG200*eft_spectra_contribution_num + index_part][index_mu_k];
            Idelta2delta200 =        peft->spectra_contributions[index_pk_type][peft->index_Idelta2delta200*eft_spectra_contribution_num + index_part][index_mu_k];
            IG2G200         =        peft->spectra_contributions[index_pk_type][peft->index_IG2G200*eft_spectra_contribution_num + index_part][index_mu_k];
            Idelta2G200     =        peft->spectra_contributions[index_pk_type][peft->index_Idelta2G200*eft_spectra_contribution_num + index_part][index_mu_k];
            FG200           = Plin * peft->spectra_contributions[index_pk_type][peft->index_FG200*eft_spectra_contribution_num + index_part][index_mu_k];
            /** 1-st RSD moment */
            I2201      =             peft->spectra_contributions[index_pk_type][peft->index_I2201     *eft_spectra_contribution_num + index_part][index_mu_k];
            I1301p3101 = Plin *      peft->spectra_contributions[index_pk_type][peft->index_I1301p3101*eft_spectra_contribution_num + index_part][index_mu_k];
            Idelta201  =             peft->spectra_contributions[index_pk_type][peft->index_Idelta201 *eft_spectra_contribution_num + index_part][index_mu_k];
            IG201      =             peft->spectra_contributions[index_pk_type][peft->index_IG201     *eft_spectra_contribution_num + index_part][index_mu_k];
            FG201      = Plin *      peft->spectra_contributions[index_pk_type][peft->index_FG201     *eft_spectra_contribution_num + index_part][index_mu_k];
            J12101     = Plin * mu * peft->spectra_contributions[index_pk_type][peft->index_J12101    *eft_spectra_contribution_num + index_part][index_mu_k];
            J11201     = Plin * mu * peft->spectra_contributions[index_pk_type][peft->index_J11201    *eft_spectra_contribution_num + index_part][index_mu_k];
            J21101     = mu *        peft->spectra_contributions[index_pk_type][peft->index_J21101    *eft_spectra_contribution_num + index_part][index_mu_k];
            Jdelta201  = mu *        peft->spectra_contributions[index_pk_type][peft->index_Jdelta201 *eft_spectra_contribution_num + index_part][index_mu_k];
            JG201      = mu *        peft->spectra_contributions[index_pk_type][peft->index_JG201     *eft_spectra_contribution_num + index_part][index_mu_k];
            /** 2-nd RSD moment */
            J12102    = Plin * (peft->spectra_contributions[index_pk_type][peft->index_J12102x   *eft_spectra_contribution_num + index_part][index_mu_k]   \
                      + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_J12102y   *eft_spectra_contribution_num + index_part][index_mu_k]);
            J21102    =         peft->spectra_contributions[index_pk_type][peft->index_J21102x   *eft_spectra_contribution_num + index_part][index_mu_k]   \
                      + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_J21102y   *eft_spectra_contribution_num + index_part][index_mu_k];
            Jdelta202 =         peft->spectra_contributions[index_pk_type][peft->index_Jdelta202x*eft_spectra_contribution_num + index_part][index_mu_k]   \
                      + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_Jdelta202y*eft_spectra_contribution_num + index_part][index_mu_k];
            JG202     =         peft->spectra_contributions[index_pk_type][peft->index_JG202x    *eft_spectra_contribution_num + index_part][index_mu_k]   \
                      + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_JG202y    *eft_spectra_contribution_num + index_part][index_mu_k];
            I2211     =         peft->spectra_contributions[index_pk_type][peft->index_I2211     *eft_spectra_contribution_num + index_part][index_mu_k];
            I1311     = Plin *  peft->spectra_contributions[index_pk_type][peft->index_I1311     *eft_spectra_contribution_num + index_part][index_mu_k];
            J12111    = Plin * mu * peft->spectra_contributions[index_pk_type][peft->index_J12111*eft_spectra_contribution_num + index_part][index_mu_k];
            J11211    = Plin * mu * peft->spectra_contributions[index_pk_type][peft->index_J11211*eft_spectra_contribution_num + index_part][index_mu_k];
            J21111    = mu *    peft->spectra_contributions[index_pk_type][peft->index_J21111    *eft_spectra_contribution_num + index_part][index_mu_k];
            N11       =         peft->spectra_contributions[index_pk_type][peft->index_N11x      *eft_spectra_contribution_num + index_part][index_mu_k]   \
                      + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_N11y      *eft_spectra_contribution_num + index_part][index_mu_k];
            /** 3-rd RSD moment */
            J21112 =         peft->spectra_contributions[index_pk_type][peft->index_J21112x*eft_spectra_contribution_num + index_part][index_mu_k]   \
                   + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_J21112y*eft_spectra_contribution_num + index_part][index_mu_k];
            J12112 = Plin * (peft->spectra_contributions[index_pk_type][peft->index_J12112x*eft_spectra_contribution_num + index_part][index_mu_k]   \
                   + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_J12112y*eft_spectra_contribution_num + index_part][index_mu_k]);
            N12    = mu *   (peft->spectra_contributions[index_pk_type][peft->index_N12x   *eft_spectra_contribution_num + index_part][index_mu_k]   \
                   + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_N12y   *eft_spectra_contribution_num + index_part][index_mu_k]);
            /** 4-th RSD moment */
            N22 =          peft->spectra_contributions[index_pk_type][peft->index_N22x*eft_spectra_contribution_num + index_part][index_mu_k]   \
                + mu*mu * (peft->spectra_contributions[index_pk_type][peft->index_N22y*eft_spectra_contribution_num + index_part][index_mu_k]   \
                + mu*mu *  peft->spectra_contributions[index_pk_type][peft->index_N22z*eft_spectra_contribution_num + index_part][index_mu_k]);


            /** - assemble the power spectrum contributions from the RSD moments only from FFTLog loop terms */
            Prsd0_loop = D4 * (2.*eft_ip.b1*eft_ip.b1 * (I2200 + 3.*I1300) + 2.*eft_ip.b1*eft_ip.b2 * Idelta200 + 2.*eft_ip.bG2*eft_ip.bG2 * IG2G200 + 4.*eft_ip.b1*eft_ip.bG2 * IG200        \
                              + 0.5*eft_ip.b2*eft_ip.b2 * Idelta2delta200 + 2.*eft_ip.b2*eft_ip.bG2 * Idelta2G200 + 8.*eft_ip.b1*(eft_ip.bG2 + 0.4*eft_ip.btd) * FG200);
            Prsd1_loop = D4*f_z * (2.*mu*mu * (2.*eft_ip.b1 * I2201 + 3.*eft_ip.b1 * I1301p3101 + eft_ip.b2 * Idelta201 + 2.*eft_ip.bG2 * IG201 + 4.*(eft_ip.bG2 + 0.4*eft_ip.btd) * FG201)   \
                                  + 4.*mu*k * (eft_ip.b1*eft_ip.b1 * (J12101 + J11201 + J21101) + 0.5*eft_ip.b1*eft_ip.b2 * Jdelta201 + eft_ip.b1*eft_ip.bG2 * JG201));
            Prsd2_loop = D4*f_z*f_z * (mu*mu*k*k * (2.*eft_ip.b1 * (2.*J12102 + J21102) + eft_ip.b2 * Jdelta202 + 2.*eft_ip.bG2 * JG202)
                                      + 2.*mu*mu*mu*mu * (I2211 + 3.*I1311) + 4.*mu*mu*mu*k * eft_ip.b1 * (J12111 + J11211 + J21111) + mu*mu*k*k * eft_ip.b1*eft_ip.b1 * N11);
            Prsd3_loop = D4*f_z*f_z*f_z * (2.*mu*mu*mu*mu*k*k * (J21112 + 2.*J12112) + 2.*mu*mu*mu*k*k*k * eft_ip.b1 * N12);
            Prsd4_loop = D4*f_z*f_z*f_z*f_z * 0.5*mu*mu*mu*mu*k*k*k*k * N22;

            /** - add counterterm contributions to UV-divergent parts */
            if (index_part == uv_divergence) {
              Prsd0_loop += D2 * eft_ip.c00 * k*k * Plin;
              Prsd1_loop += D2*f_z * eft_ip.c10 * mu*mu*k*k * Plin;
              Prsd2_loop += D2*f_z*f_z * (eft_ip.c20 + eft_ip.c22 * mu*mu) * mu*mu*k*k * Plin;
              Prsd3_loop += D2*f_z*f_z*f_z * (eft_ip.c30 + eft_ip.c32 * mu*mu) * mu*mu*mu*mu*k*k * Plin;
              Prsd4_loop += D2*f_z*f_z*f_z*f_z * eft_ip.c42 * mu*mu*mu*mu*mu*mu*k*k * Plin;
            }

            if (index_part == finite_part) {
              if (peft->hp->integration_mode == fftlog) {
                Prsd2_loop += -D4*f_z*f_z * eft_ip.b1*eft_ip.b1 * sigma_sq(peft, -1, index_pk_type)/3. * mu*mu*k*k * Plin;
                Prsd3_loop += -2.*D4*f_z*f_z*f_z * eft_ip.b1 * sigma_sq(peft, -1, index_pk_type)/3. * mu*mu*mu*mu*k*k * Plin;
                Prsd4_loop += -D4*f_z*f_z*f_z*f_z * sigma_sq(peft, -1, index_pk_type)/3. * mu*mu*mu*mu*mu*mu*k*k * Plin;
              }
              else {
                sigmav_mu = peft->spectra_contributions[index_pk_type][peft->index_sigmav_mu*eft_spectra_contribution_num + index_part][index_mu_k];
                Prsd2_loop += -D4*f_z*f_z * eft_ip.b1*eft_ip.b1 * sigmav_mu * mu*mu*k*k * Plin;
                Prsd3_loop += -2.*D4*f_z*f_z*f_z * eft_ip.b1 * sigmav_mu * mu*mu*mu*mu*k*k * Plin;
                Prsd4_loop += -D4*f_z*f_z*f_z*f_z * sigmav_mu * mu*mu*mu*mu*mu*mu*k*k * Plin;
              }
            }

            pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part][index_mu*peft->k_size + index_k] = Prsd0_loop + Prsd1_loop + Prsd2_loop + Prsd3_loop + Prsd4_loop;
            break;

          default:
            break;
          }
        }
      }
    }
  }

  /** - sum the parts of the power spectrum and accumulate in the last defined part index,
   *    add the finite part last to utilize cancellations between divergent parts */
  #pragma omp for schedule(static), collapse(3)
  for (index_list = 0; index_list < pk_types_loops_size; index_list++) {
    for (index_mu = 0; index_mu < mu_size; index_mu++) {
      for (index_k = 0; index_k < peft->k_size; index_k++) {
        index_pk_type = pk_types_loops[index_list];

        for (index_part = eft_spectra_contribution_num-2; index_part >= 0; index_part--) {
          pkmu_loop[index_pk_type*eft_spectra_contribution_num + eft_spectra_contribution_num-1][index_mu*peft->k_size + index_k]   \
            += pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part][index_mu*peft->k_size + index_k];
        }
      }
    }
  }

  /** - combine with the linear spectra to produce the full spectrum with loop corrections */
  #pragma omp for schedule(static), collapse(2)
  for (index_mu = 0; index_mu < mu_size; index_mu++) {
    for (index_k = 0; index_k < peft->k_size; index_k++) {
      if (peft->hp->use_interpolation) { index_mu_k = index_k; }
      else { index_mu_k = index_mu*peft->k_size + index_k; }

      k = exp(peft->ln_k[index_mu_k]);
      mu = muvec[index_mu];

      switch (index_pk_out_type)
      {
      case Pdd_mm_real:
        pkmu[index_mu*peft->k_size + index_k] = peft->pk_l[pk_ir_resummed_nlo][index_mu_k]   \
                        + pkmu_loop[pk_ir_resummed_lo*eft_spectra_contribution_num + eft_spectra_contribution_num-1][index_mu*peft->k_size + index_k];
        break;

      case Pdd_mm_real_no_IR_resum:
        pkmu[index_mu*peft->k_size + index_k] = peft->pk_l[pk_lin][index_mu_k]   \
                        + pkmu_loop[pk_lin*eft_spectra_contribution_num + eft_spectra_contribution_num-1][index_mu*peft->k_size + index_k];
        break;
      
      case Pdd_mm_22:
      case Pdd_mm_13:
        pkmu[index_mu*peft->k_size + index_k] = pkmu_loop[pk_ir_resummed_lo*eft_spectra_contribution_num + eft_spectra_contribution_num-1][index_mu*peft->k_size + index_k];
        break;

      case Pdd_mm_22_no_IR_resum:
      case Pdd_mm_13_no_IR_resum:
        pkmu[index_mu*peft->k_size + index_k] = pkmu_loop[pk_lin*eft_spectra_contribution_num + eft_spectra_contribution_num-1][index_mu*peft->k_size + index_k];
        break;

      case Pdd_mm_rsd:
        // TODO
        break;

      case Pdd_hh_real:
        pkmu[index_mu*peft->k_size + index_k] = eft_ip.b1*eft_ip.b1 * peft->pk_l[pk_ir_resummed_nlo][index_mu_k]   \
                        + pkmu_loop[pk_ir_resummed_lo*eft_spectra_contribution_num + eft_spectra_contribution_num-1][index_mu*peft->k_size + index_k];
        break;

      case Pdd_hh_rsd:
        if (peft->hp->integration_mode == fftlog) {
          sigma2_tot_z0 = (1. + f_z*mu*mu*(2. + f_z)) * peft->Sigma2_ir + f_z*f_z*mu*mu*(mu*mu - 1.) * peft->dSigma2_ir;  /** D2 * sigma2_tot_z0 = sigma2_tot at z */
          pkmu[index_mu*peft->k_size + index_k] = pow(eft_ip.b1 + f_z*mu*mu, 2.) * peft->pk_l[pkmu_rsd_ir_resummed_nlo][index_mu*peft->k_size + index_k]    \
                                                  + pkmu_loop[pk_nowiggle*eft_spectra_contribution_num + eft_spectra_contribution_num-1][index_mu*peft->k_size + index_k]   \
                                          + exp(-k*k * D2 * sigma2_tot_z0) * (pkmu_loop[pk_lin*eft_spectra_contribution_num + eft_spectra_contribution_num-1][index_mu*peft->k_size + index_k]   \
                                                                            - pkmu_loop[pk_nowiggle*eft_spectra_contribution_num + eft_spectra_contribution_num-1][index_mu*peft->k_size + index_k]);
        }
        else {
          pkmu[index_mu*peft->k_size + index_k] = pow(eft_ip.b1 + f_z*mu*mu, 2.) * peft->pk_l[pkmu_rsd_ir_resummed_nlo][index_mu*peft->k_size + index_k]    \
                                                  + pkmu_loop[pkmu_rsd_ir_resummed_lo*eft_spectra_contribution_num + eft_spectra_contribution_num-1][index_mu*peft->k_size + index_k];
        }
        break;

      default:
        break;
      }
    }
  }

  } /** - end of parallel region */


  for (index_list = 0; index_list < pk_types_loops_size; index_list++) {
    for (index_part = 0; index_part < eft_spectra_contribution_num; index_part++) {
      index_pk_type = pk_types_loops[index_list];
      free(pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part]);
    }
  }
  free(pk_types_loops);
  free(pk_types);
  free(pk_types_outside_loops);
  free(mu_real);

  return abort;
}

//assumes EdS scaling
/** TODO: deprecated, to be removed */
// int eft_build_nonlinear_power_spectrum_wedges_multiple(
//                                   struct eft * peft,
//                                   struct background * pba,
//                                   struct primordial * ppm,
//                                   struct fourier * pfo,
//                                   const short index_pk_out_type,
//                                   const double * const z,
//                                   const double * const D_z,
//                                   const double * const f_z,
//                                   const int z_size,
//                                   const double * muvec,
//                                   int mu_size,
//                                   struct eft_input_parameters eft_ip,
//                                   double ** pkmu) {  /**< pkmu[index_z][index_mu*peft->k_size + index_k] */

//   int index_pk_type, index_list, index_moment, index_part, index_z, index_mu, index_mu_k, index_k, it, size;
//   int * pk_types_loops, pk_types_loops_size = 0, * pk_types, pk_types_size = 0, * pk_types_rsd, pk_types_outside_loops_size = 0;
//   double * pkmu_loop[pk_type_num*eft_spectra_contribution_num]; /**< nonlinear power spectrum still with separated finite and divergent parts */
//   double * mu_real = NULL, D2[z_size], D4[z_size];
//   short abort = _FALSE_;
//   double  k, mu, Plin, I2200, I1300, Idelta200, IG200, Idelta2delta200, IG2G200, Idelta2G200, FG200,  \
//           I2201, Idelta201, IG201, J21101, Jdelta201, JG201, FG201, I1301p3101, J12101, J11201,       \
//           J21102, Jdelta202, JG202, I2211, J21111, N11, J12102, I1311, J12111, J11211,                \
//           J21112, N12, J12112,                                                                        \
//           N22,                                                                                        \
//           Preal_loop, Prsd0_loop, Prsd1_loop, Prsd2_loop, Prsd3_loop, Prsd4_loop, sum, sigma2_tot_z0;  /** RSD moments multiplied by i^n/n! (k mu)^n only containing loop terms */
//   /** sigma2_ir_at_z = D2 * sigma2_ir */

//   if (!peft->hp->use_interpolation) {
//     muvec = peft->mu; mu_size = peft->mu_size;
//     // TODO
//     class_stop(peft->error_message, "Not implemented yet.")
//   }

//   /** - if a real space spectrum is requested, ignore the mu input and let muvec point to a valid address with a constant value */
//   if (!eft_rsd_out_indicator(index_pk_out_type)) {
//     class_alloc(mu_real, mu_size*sizeof(double), peft->error_message);
//     for (index_mu = 0; index_mu < mu_size; index_mu++) { mu_real[it] = 0.; }
//     muvec = mu_real;
//   }

//   /** - get a list of necessary pk_types for the chosen output spectrum of which loop corrections shall be calculated */
//   class_call(eft_necessary_pk_types_loops(peft, index_pk_out_type, &pk_types_loops, &pk_types_loops_size),
//              peft->error_message, peft->error_message);
//   /** - get a list of all necessary pk_types for the chosen output spectrum */
//   class_call(eft_necessary_pk_types_total(peft, index_pk_out_type, &pk_types, &pk_types_size),
//              peft->error_message, peft->error_message);

//   /** - find the RSD / mu-dependent spectra in this list (since the mu-independent spectra have to be loaded before) */
//   for (index_list = 0; index_list < pk_types_size; index_list++) {
//     if (eft_rsd_indicator(pk_types[index_list])) { pk_types_outside_loops_size++; }
//   }
//   class_alloc(pk_types_rsd, pk_types_outside_loops_size*sizeof(int), peft->error_message);
//   it = 0;
//   for (index_list = 0; index_list < pk_types_size; index_list++) {
//     if (eft_rsd_indicator(pk_types[index_list])) { pk_types_rsd[it++] = pk_types[index_list]; }
//   }
//   /** - and load them into peft->pk_l (will be rescaled by D2 later) */
//   class_call(eft_load_linear_spectra(pba, pfo, ppm,
//                                      peft,
//                                      peft->z0,
//                                      peft->f_z0,
//                                      peft->D_z0,
//                                      pk_types_rsd,
//                                      pk_types_outside_loops_size,
//                                      muvec,
//                                      mu_size),
//               peft->error_message, peft->error_message);

//   /** - allocate arrays for the loop corrections */
//   for (index_list = 0; index_list < pk_types_loops_size; index_list++) {
//     for (index_part = 0; index_part < eft_spectra_contribution_num; index_part++) {
//       index_pk_type = pk_types_loops[index_list];
//       class_alloc(pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part], z_size*mu_size*peft->k_size*sizeof(double), peft->error_message);
//     }
//   }

//   /** - growth function rescale factors w.r.t. peft->z0 for tree and 1-loop contributions */
//   for (index_z = 0; index_z < z_size; index_z++) {
//     D2[index_z] = pow(D_z[index_z]/peft->D_z0, 2.);
//     D4[index_z] = pow(D_z[index_z]/peft->D_z0, 4.);
//   }

//   #pragma omp parallel shared(peft, f_z, D2, D4, pk_types_loops, pk_types_loops_size, muvec, mu_size, z_size, index_pk_out_type, eft_ip, pkmu_loop, pkmu),   \
//                        private(index_list, index_part, index_k, index_mu, index_mu_k, index_z, index_pk_type,  \
//                                k, mu, Plin, I2200, I1300, Idelta200, IG200, Idelta2delta200, IG2G200, Idelta2G200, FG200,  \
//                                I2201, Idelta201, IG201, J21101, Jdelta201, JG201, FG201, I1301p3101, J12101, J11201,       \
//                                J21102, Jdelta202, JG202, I2211, J21111, N11, J12102, I1311, J12111, J11211,                \
//                                J21112, N12, J12112,                                                                        \
//                                N22,                                                                                        \
//                                Preal_loop, Prsd0_loop, Prsd1_loop, Prsd2_loop, Prsd3_loop, Prsd4_loop, sum, sigma2_tot_z0), default(none)
//   {
//   /** - compute the loop power spectra of every internal pk_type in the list for every k, mu, still separated into finite, UV, IR and pole parts */
//   #pragma omp for schedule(static), collapse(4)
//   for (index_list = 0; index_list < pk_types_loops_size; index_list++) {
//     for (index_part = 0; index_part < eft_spectra_contribution_num; index_part++) {
//       for (index_k = 0; index_k < peft->k_size; index_k++) {
//         for (index_mu = 0; index_mu < mu_size; index_mu++) {
//           index_pk_type = pk_types_loops[index_list];

//           if (peft->hp->use_interpolation) { index_mu_k = index_k; }
//           else { index_mu_k = index_mu*peft->k_size + index_k; }

//           k = exp(peft->ln_k[index_mu_k]);
//           mu = muvec[index_mu];
//           Plin = peft->pk_l[index_pk_type][index_mu_k];

//           switch (index_pk_out_type)
//           {
//           case Pdd_mm_real:
//             I2200           =        peft->spectra_contributions[index_pk_type][peft->index_I2200*eft_spectra_contribution_num + index_part][index_mu_k];
//             I1300           = Plin * peft->spectra_contributions[index_pk_type][peft->index_I1300*eft_spectra_contribution_num + index_part][index_mu_k];

//             for (index_z = 0; index_z < z_size; index_z++) {
//               Preal_loop = D4[index_z] * 2.*(I2200 + 3.*I1300);
//               if (index_part == uv_divergence) {
//                 Preal_loop += -D2[index_z] * 2.*k*k * eft_ip.cs2 * Plin;
//               }

//               pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part][(index_z*mu_size + index_mu)*peft->k_size + index_k] = Preal_loop;
//             }
//             break;

//           case Pdd_mm_rsd:
//             // TODO
//             break;

//           case Pdd_hh_real:
//             I2200           =        peft->spectra_contributions[index_pk_type][peft->index_I2200*eft_spectra_contribution_num + index_part][index_mu_k];
//             I1300           = Plin * peft->spectra_contributions[index_pk_type][peft->index_I1300*eft_spectra_contribution_num + index_part][index_mu_k];
//             Idelta200       =        peft->spectra_contributions[index_pk_type][peft->index_Idelta200*eft_spectra_contribution_num + index_part][index_mu_k];
//             IG200           =        peft->spectra_contributions[index_pk_type][peft->index_IG200*eft_spectra_contribution_num + index_part][index_mu_k];
//             Idelta2delta200 =        peft->spectra_contributions[index_pk_type][peft->index_Idelta2delta200*eft_spectra_contribution_num + index_part][index_mu_k];
//             IG2G200         =        peft->spectra_contributions[index_pk_type][peft->index_IG2G200*eft_spectra_contribution_num + index_part][index_mu_k];
//             Idelta2G200     =        peft->spectra_contributions[index_pk_type][peft->index_Idelta2G200*eft_spectra_contribution_num + index_part][index_mu_k];
//             FG200           = Plin * peft->spectra_contributions[index_pk_type][peft->index_FG200*eft_spectra_contribution_num + index_part][index_mu_k];

//             for (index_z = 0; index_z < z_size; index_z++) {
//               Prsd0_loop = D4[index_z] * (2.*eft_ip.b1*eft_ip.b1 * (I2200 + 3.*I1300) + 2.*eft_ip.b1*eft_ip.b2 * Idelta200 + 2.*eft_ip.bG2*eft_ip.bG2 * IG2G200 + 4.*eft_ip.b1*eft_ip.bG2 * IG200        \
//                                 + 0.5*eft_ip.b2*eft_ip.b2 * Idelta2delta200 + 2.*eft_ip.b2*eft_ip.bG2 * Idelta2G200 + 8.*eft_ip.b1*(eft_ip.bG2 + 0.4*eft_ip.btd) * FG200);
//               if (index_part == uv_divergence) {
//                 Prsd0_loop += -D2[index_z] * 2.*k*k * (eft_ip.b1 * (eft_ip.R2 + eft_ip.b1 * eft_ip.cs2)) * Plin;
//               }

//               pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part][(index_z*mu_size + index_mu)*peft->k_size + index_k] = Prsd0_loop;
//             }
//             break;

//           case Pdd_hh_rsd:
//             /** 0-th RSD moment */
//             I2200           =        peft->spectra_contributions[index_pk_type][peft->index_I2200*eft_spectra_contribution_num + index_part][index_mu_k];
//             I1300           = Plin * peft->spectra_contributions[index_pk_type][peft->index_I1300*eft_spectra_contribution_num + index_part][index_mu_k];
//             Idelta200       =        peft->spectra_contributions[index_pk_type][peft->index_Idelta200*eft_spectra_contribution_num + index_part][index_mu_k];
//             IG200           =        peft->spectra_contributions[index_pk_type][peft->index_IG200*eft_spectra_contribution_num + index_part][index_mu_k];
//             Idelta2delta200 =        peft->spectra_contributions[index_pk_type][peft->index_Idelta2delta200*eft_spectra_contribution_num + index_part][index_mu_k];
//             IG2G200         =        peft->spectra_contributions[index_pk_type][peft->index_IG2G200*eft_spectra_contribution_num + index_part][index_mu_k];
//             Idelta2G200     =        peft->spectra_contributions[index_pk_type][peft->index_Idelta2G200*eft_spectra_contribution_num + index_part][index_mu_k];
//             FG200           = Plin * peft->spectra_contributions[index_pk_type][peft->index_FG200*eft_spectra_contribution_num + index_part][index_mu_k];
//             /** 1-st RSD moment */
//             I2201      =             peft->spectra_contributions[index_pk_type][peft->index_I2201     *eft_spectra_contribution_num + index_part][index_mu_k];
//             I1301p3101 = Plin *      peft->spectra_contributions[index_pk_type][peft->index_I1301p3101*eft_spectra_contribution_num + index_part][index_mu_k];
//             Idelta201  =             peft->spectra_contributions[index_pk_type][peft->index_Idelta201 *eft_spectra_contribution_num + index_part][index_mu_k];
//             IG201      =             peft->spectra_contributions[index_pk_type][peft->index_IG201     *eft_spectra_contribution_num + index_part][index_mu_k];
//             FG201      = Plin *      peft->spectra_contributions[index_pk_type][peft->index_FG201     *eft_spectra_contribution_num + index_part][index_mu_k];
//             J12101     = Plin * mu * peft->spectra_contributions[index_pk_type][peft->index_J12101    *eft_spectra_contribution_num + index_part][index_mu_k];
//             J11201     = Plin * mu * peft->spectra_contributions[index_pk_type][peft->index_J11201    *eft_spectra_contribution_num + index_part][index_mu_k];
//             J21101     = mu *        peft->spectra_contributions[index_pk_type][peft->index_J21101    *eft_spectra_contribution_num + index_part][index_mu_k];
//             Jdelta201  = mu *        peft->spectra_contributions[index_pk_type][peft->index_Jdelta201 *eft_spectra_contribution_num + index_part][index_mu_k];
//             JG201      = mu *        peft->spectra_contributions[index_pk_type][peft->index_JG201     *eft_spectra_contribution_num + index_part][index_mu_k];
//             /** 2-nd RSD moment */
//             J12102    = Plin * (peft->spectra_contributions[index_pk_type][peft->index_J12102x   *eft_spectra_contribution_num + index_part][index_mu_k]   \
//                       + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_J12102y   *eft_spectra_contribution_num + index_part][index_mu_k]);
//             J21102    =         peft->spectra_contributions[index_pk_type][peft->index_J21102x   *eft_spectra_contribution_num + index_part][index_mu_k]   \
//                       + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_J21102y   *eft_spectra_contribution_num + index_part][index_mu_k];
//             Jdelta202 =         peft->spectra_contributions[index_pk_type][peft->index_Jdelta202x*eft_spectra_contribution_num + index_part][index_mu_k]   \
//                       + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_Jdelta202y*eft_spectra_contribution_num + index_part][index_mu_k];
//             JG202     =         peft->spectra_contributions[index_pk_type][peft->index_JG202x    *eft_spectra_contribution_num + index_part][index_mu_k]   \
//                       + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_JG202y    *eft_spectra_contribution_num + index_part][index_mu_k];
//             I2211     =         peft->spectra_contributions[index_pk_type][peft->index_I2211     *eft_spectra_contribution_num + index_part][index_mu_k];
//             I1311     = Plin *  peft->spectra_contributions[index_pk_type][peft->index_I1311     *eft_spectra_contribution_num + index_part][index_mu_k];
//             J12111    = Plin * mu * peft->spectra_contributions[index_pk_type][peft->index_J12111*eft_spectra_contribution_num + index_part][index_mu_k];
//             J11211    = Plin * mu * peft->spectra_contributions[index_pk_type][peft->index_J11211*eft_spectra_contribution_num + index_part][index_mu_k];
//             J21111    = mu *    peft->spectra_contributions[index_pk_type][peft->index_J21111    *eft_spectra_contribution_num + index_part][index_mu_k];
//             N11       =         peft->spectra_contributions[index_pk_type][peft->index_N11x      *eft_spectra_contribution_num + index_part][index_mu_k]   \
//                       + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_N11y      *eft_spectra_contribution_num + index_part][index_mu_k];
//             /** 3-rd RSD moment */
//             J21112 =         peft->spectra_contributions[index_pk_type][peft->index_J21112x*eft_spectra_contribution_num + index_part][index_mu_k]   \
//                    + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_J21112y*eft_spectra_contribution_num + index_part][index_mu_k];
//             J12112 = Plin * (peft->spectra_contributions[index_pk_type][peft->index_J12112x*eft_spectra_contribution_num + index_part][index_mu_k]   \
//                    + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_J12112y*eft_spectra_contribution_num + index_part][index_mu_k]);
//             N12    = mu *   (peft->spectra_contributions[index_pk_type][peft->index_N12x   *eft_spectra_contribution_num + index_part][index_mu_k]   \
//                    + mu*mu * peft->spectra_contributions[index_pk_type][peft->index_N12y   *eft_spectra_contribution_num + index_part][index_mu_k]);
//             /** 4-th RSD moment */
//             N22 =          peft->spectra_contributions[index_pk_type][peft->index_N22x*eft_spectra_contribution_num + index_part][index_mu_k]   \
//                 + mu*mu * (peft->spectra_contributions[index_pk_type][peft->index_N22y*eft_spectra_contribution_num + index_part][index_mu_k]   \
//                 + mu*mu *  peft->spectra_contributions[index_pk_type][peft->index_N22z*eft_spectra_contribution_num + index_part][index_mu_k]);

//             /** - assemble the power spectrum controbutions from the RSD moments only from FFTLog loop terms */
//             for (index_z = 0; index_z < z_size; index_z++) {
//               Prsd0_loop = D4[index_z] * (2.*eft_ip.b1*eft_ip.b1 * (I2200 + 3.*I1300) + 2.*eft_ip.b1*eft_ip.b2 * Idelta200 + 2.*eft_ip.bG2*eft_ip.bG2 * IG2G200 + 4.*eft_ip.b1*eft_ip.bG2 * IG200        \
//                                 + 0.5*eft_ip.b2*eft_ip.b2 * Idelta2delta200 + 2.*eft_ip.b2*eft_ip.bG2 * Idelta2G200 + 8.*eft_ip.b1*(eft_ip.bG2 + 0.4*eft_ip.btd) * FG200);
//               Prsd1_loop = D4[index_z]*f_z[index_z] * (2.*mu*mu * (2.*eft_ip.b1 * I2201 + 3.*eft_ip.b1 * I1301p3101 + eft_ip.b2 * Idelta201 + 2.*eft_ip.bG2 * IG201 + 4.*(eft_ip.bG2 + 0.4*eft_ip.btd) * FG201)   \
//                                     + 4.*mu*k * (eft_ip.b1*eft_ip.b1 * (J12101 + J11201 + J21101) + 0.5*eft_ip.b1*eft_ip.b2 * Jdelta201 + eft_ip.b1*eft_ip.bG2 * JG201));
//               Prsd2_loop = D4[index_z]*f_z[index_z]*f_z[index_z] * (mu*mu*k*k * (2.*eft_ip.b1 * (2.*J12102 + J21102) + eft_ip.b2 * Jdelta202 + 2.*eft_ip.bG2 * JG202)
//                                         + 2.*mu*mu*mu*mu * (I2211 + 3.*I1311) + 4.*mu*mu*mu*k * eft_ip.b1 * (J12111 + J11211 + J21111) + mu*mu*k*k * eft_ip.b1*eft_ip.b1 * N11);
//               Prsd3_loop = D4[index_z]*f_z[index_z]*f_z[index_z]*f_z[index_z] * (2.*mu*mu*mu*mu*k*k * (J21112 + 2.*J12112) + 2.*mu*mu*mu*k*k*k * eft_ip.b1 * N12);
//               Prsd4_loop = D4[index_z]*f_z[index_z]*f_z[index_z]*f_z[index_z]*f_z[index_z] * 0.5*mu*mu*mu*mu*k*k*k*k * N22;

//               /** - add counterterm contributions to UV-divergent parts */
//               if (index_part == uv_divergence) {
//                 Prsd0_loop += D2[index_z] * eft_ip.c00 * k*k * Plin;
//                 Prsd1_loop += D2[index_z]*f_z[index_z] * eft_ip.c10 * mu*mu*k*k * Plin;
//                 Prsd2_loop += D2[index_z]*f_z[index_z]*f_z[index_z] * ((eft_ip.c20 + eft_ip.c22 * mu*mu) - D2[index_z] * eft_ip.b1*eft_ip.b1 * sigma_sq(peft, -1, index_pk_type)/3.) * mu*mu*k*k * Plin;
//                 Prsd3_loop += D2[index_z]*f_z[index_z]*f_z[index_z]*f_z[index_z] * ((eft_ip.c30 + eft_ip.c32 * mu*mu) - 2.*D2[index_z] * eft_ip.b1 * sigma_sq(peft, -1, index_pk_type)/3.) * mu*mu*mu*mu*k*k * Plin;
//                 Prsd4_loop += D2[index_z]*f_z[index_z]*f_z[index_z]*f_z[index_z]*f_z[index_z] * (eft_ip.c42 - D2[index_z] * sigma_sq(peft, -1, index_pk_type)/3.) * mu*mu*mu*mu*mu*mu*k*k * Plin;
//               }

//               pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part][(index_z*mu_size + index_mu)*peft->k_size + index_k] = Prsd0_loop + Prsd1_loop + Prsd2_loop + Prsd3_loop + Prsd4_loop;
//             }
//             break;

//           default:
//             break;
//           }
//         }
//       }
//     }
//   }

//   /** - sum the parts of the power spectrum and accumulate in the last defined part index,
//    *    add the finite part last to utilize cancellations between divergent parts */
//   #pragma omp for schedule(static), collapse(4)
//   for (index_list = 0; index_list < pk_types_loops_size; index_list++) {
//     for (index_z = 0; index_z < z_size; index_z++) {
//       for (index_mu = 0; index_mu < mu_size; index_mu++) {
//         for (index_k = 0; index_k < peft->k_size; index_k++) {
//           index_pk_type = pk_types_loops[index_list];

//           for (index_part = eft_spectra_contribution_num-2; index_part >= 0; index_part--) {
//             pkmu_loop[index_pk_type*eft_spectra_contribution_num + eft_spectra_contribution_num-1][(index_z*mu_size + index_mu)*peft->k_size + index_k]   \
//               += pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part][(index_z*mu_size + index_mu)*peft->k_size + index_k];
//           }
//         }
//       }
//     }
//   }

//   /** - combine with the linear spectra to produce the full spectrum with loop corrections */
//   #pragma omp for schedule(static), collapse(3)
//   for (index_z = 0; index_z < z_size; index_z++) {
//     for (index_mu = 0; index_mu < mu_size; index_mu++) {
//       for (index_k = 0; index_k < peft->k_size; index_k++) {
//         if (peft->hp->use_interpolation) { index_mu_k = index_k; }
//         else { index_mu_k = index_mu*peft->k_size + index_k; }

//         k = exp(peft->ln_k[index_mu_k]);
//         mu = muvec[index_mu];

//         switch (index_pk_out_type)
//         {
//         case Pdd_mm_real:
//           pkmu[index_z][index_mu*peft->k_size + index_k] = D2[index_z] * peft->pk_l[pk_ir_resummed_nlo][index_mu_k]   \
//                           + pkmu_loop[pk_ir_resummed_lo*eft_spectra_contribution_num + eft_spectra_contribution_num-1][(index_z*mu_size + index_mu)*peft->k_size + index_k];
//           break;

//         case Pdd_mm_rsd:
//           // TODO
//           break;

//         case Pdd_hh_real:
//           pkmu[index_z][index_mu*peft->k_size + index_k] = D2[index_z] * eft_ip.b1*eft_ip.b1 * peft->pk_l[pk_ir_resummed_nlo][index_mu_k]   \
//                           + pkmu_loop[pk_ir_resummed_lo*eft_spectra_contribution_num + eft_spectra_contribution_num-1][(index_z*mu_size + index_mu)*peft->k_size + index_k];
//           break;

//         case Pdd_hh_rsd:
//           sigma2_tot_z0 = (1. + f_z[index_z]*mu*mu*(2. + f_z[index_z])) * peft->Sigma2_ir + f_z[index_z]*f_z[index_z]*mu*mu*(mu*mu - 1.) * peft->dSigma2_ir;  /** D2 * sigma2_tot_z0 = sigma2_tot at z */
//           pkmu[index_z][index_mu*peft->k_size + index_k] = D2[index_z] * pow(eft_ip.b1 + f_z[index_z]*mu*mu, 2.) * peft->pk_l[pkmu_rsd_ir_resummed_nlo][index_mu*peft->k_size + index_k]    \
//                                                   + pkmu_loop[pk_nowiggle*eft_spectra_contribution_num + eft_spectra_contribution_num-1][(index_z*mu_size + index_mu)*peft->k_size + index_k]   \
//                                           + exp(-k*k * D2[index_z] * sigma2_tot_z0) * (pkmu_loop[pk_lin*eft_spectra_contribution_num + eft_spectra_contribution_num-1][(index_z*mu_size + index_mu)*peft->k_size + index_k]   \
//                                                                             - pkmu_loop[pk_nowiggle*eft_spectra_contribution_num + eft_spectra_contribution_num-1][(index_z*mu_size + index_mu)*peft->k_size + index_k]);
//           break;

//         default:
//           break;
//         }
//       }
//     }
//   }

//   } /** - end of parallel region */


//   for (index_list = 0; index_list < pk_types_loops_size; index_list++) {
//     for (index_part = 0; index_part < eft_spectra_contribution_num; index_part++) {
//       index_pk_type = pk_types_loops[index_list];
//       free(pkmu_loop[index_pk_type*eft_spectra_contribution_num + index_part]);
//     }
//   }
//   free(pk_types_loops);
//   free(pk_types);
//   free(pk_types_rsd);
//   free(mu_real);

//   return abort;
// }
