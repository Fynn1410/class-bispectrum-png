/** @file IR_res.c 
 * 
 * author: Christian Radermacher, 2024
 * based on prototype version of Azadeh Moradinezhad Dizgah & Dennis Linde
 * 
 * this module can be used by others without ever referencing eft structures
 * contains infrared-resummed linear matter power spectra while the wiggle-nowiggle separation is performed in wnw_split.c
*/

#include "IR_res.h"

/** Integration measure for Sigma^2(z):  1 - j_0(q_osc) + 2 j_2(q_osc) 
 * error is always below 1e-10
*/
inline double ir_sigma2_measure(const double q_osc) {
  return (q_osc > 0.01) ? (1. - sin(q_osc) / q_osc + 2. * ((3. / (q_osc*q_osc) - 1.) * sin(q_osc) / q_osc - 3. * cos(q_osc) / (q_osc*q_osc) ))    \
                        : (3.*q_osc*q_osc/10. - q_osc*q_osc*q_osc*q_osc/56.);
}

/**
 * @brief Compute the suppression factor Sigma^2(z) in real space
 *        d(Sigma^2(z)) = 1/(6 pi^2) dln(q) q P_nw(z, q) (1 - j_0(q/k_bao) + 2 j_2(q/k_bao))
 * @param pba       Input: pointer to background structure
 * @param ppm       Input: pointer to primordial structure
 * @param pfo       Input: pointer to fourier structure
 * @param z         Input: redshift
 * @param k_split   Input: long/short mode separating scale (end of integration region)
 * @param k_bao     Input: wavenumber of BAO scale (default h/(110. Mpc))
 *
 * @return value of IR resummation suppression factor at z
 */
double eft_ir_sigma2(
              struct background * pba,
              struct primordial * ppm,
              struct fourier * pfo, 
              const double z, 
              const double k_split,
              const double k_bao) {

  int split_index = 0, index_y, index_ddy, index_num = 0, it_q;
  double q, q_osc, result;
  double *ln_q, *intg_splines, *pk_nw;

  /** - search for k_split in pfo->ln_k */
  class_call(array_hunt_ascending(pfo->ln_k,
                                  pfo->k_size,
                                  log(k_split),
                                  &split_index,
                                  pfo->error_message),
            pfo->error_message,
            pfo->error_message);
  
  class_define_index(index_y, _TRUE_, index_num, 1);
  class_define_index(index_ddy, _TRUE_, index_num, 1);

  class_alloc(ln_q, (split_index+1)*sizeof(double),
              pfo->error_message);
  class_alloc(pk_nw, (split_index+1)*sizeof(double),
              pfo->error_message);
  class_alloc(intg_splines, 2*(split_index+1)*sizeof(double),
              pfo->error_message);

  /** - copy pfo->ln_k to ln_q and overwrite the last value with ln(k_split) */
  memcpy(ln_q, pfo->ln_k, split_index*sizeof(double));
  ln_q[split_index] = log(k_split);

  /** - get the nowiggle spectrum at ln_q and z */
  class_call(fourier_pk_nw_at_kvec_and_z(pba, ppm, pfo,
                                        linear,
                                        ln_q,
                                        split_index+1,
                                        z,
                                        pk_nw),
            pfo->error_message,
            pfo->error_message);

  /** - prepare the integrand */
  for (it_q = 0; it_q < split_index+1; it_q++)
  {
    q = exp(ln_q[it_q]);
    q_osc = q / k_bao;
    intg_splines[it_q*index_num + index_y] = q * pk_nw[it_q] * ir_sigma2_measure(q_osc);
                                /** q * P_nw(z, q) * (1 - j_0(q/k_BAO) + 2 j_2(q/k_bao)) */
  }

  /** - spline the integrand */
  class_call(array_spline_table_line_to_line(ln_q,
                                            split_index+1,
                                            intg_splines,
                                            index_num,
                                            index_y,
                                            index_ddy,
                                            _SPLINE_EST_DERIV_,
                                            pfo->error_message),
            pfo->error_message,
            pfo->error_message);

  /** - integrate the spline */
  class_call(array_integrate_all_spline_table_line_to_line(ln_q,
                                                          split_index+1,
                                                          intg_splines,
                                                          index_num,
                                                          index_y,
                                                          index_ddy,
                                                          &result,
                                                          pfo->error_message),
            pfo->error_message,
            pfo->error_message);

  // for (it_q = 0; it_q < split_index+1; it_q++)
  //   fprintf(stderr, "%.15e  %.15e  %.15e \n", exp(ln_q[it_q]), intg_splines[it_q*index_num + index_y], intg_splines[it_q*index_num + index_ddy]);


  /** - multiply with prefactor */
  result *= 1./(6.*_PI_*_PI_);

  free(intg_splines);
  free(ln_q);
  free(pk_nw);

  return result;
}

/** Integration measure for Sigma^2(z):  j_2(q_osc)
 * error is always below 1e-10, except at the roots
*/
inline double ir_dsigma2_measure(const double q_osc) {
  return (q_osc > 0.01) ? ((3. / (q_osc*q_osc) - 1.) * sin(q_osc) / q_osc - 3. * cos(q_osc) / (q_osc*q_osc))    \
                        : (q_osc*q_osc/15. - q_osc*q_osc*q_osc*q_osc/210.);
}

/**
 * @brief Compute the suppression factor IR_sigma2 in redshift space (CLASS-PT: 2004.106007v2 p. 15)
 *        d(DSigma^2(z)) = 1/(2 pi^2) dln(q) q P_nw(z, q) j_2(q/k_bao)
 *
 * @param pba       Input: pointer to background structure
 * @param ppm       Input: pointer to primordial structure
 * @param pfo       Input: pointer to fourier structure
 * @param z         Input: redshift
 * @param k_split   Input: long/short mode separating scale (end of integration region)
 * @param k_bao     Input: wavenumber of BAO scale (default h/(110. Mpc))
 *
 * @return value of IR resummation suppression factor at z
 */
double eft_ir_dsigma2(
              struct background * pba,
              struct primordial * ppm,
              struct fourier * pfo, 
              const double z, 
              const double k_split,
              const double k_bao) {

  int split_index = 0, index_y, index_ddy, index_num = 0, it_q;
  double q, q_osc, result;
  double *ln_q, *intg_splines, *pk_nw;

  /** - search for k_split in pfo->ln_k */
  class_call(array_hunt_ascending(pfo->ln_k,
                                  pfo->k_size,
                                  log(k_split),
                                  &split_index,
                                  pfo->error_message),
            pfo->error_message,
            pfo->error_message);
  
  class_define_index(index_y, _TRUE_, index_num, 1);
  class_define_index(index_ddy, _TRUE_, index_num, 1);

  class_alloc(ln_q, (split_index+1)*sizeof(double),
              pfo->error_message);
  class_alloc(pk_nw, (split_index+1)*sizeof(double),
              pfo->error_message);
  class_alloc(intg_splines, 2*(split_index+1)*sizeof(double),
              pfo->error_message);

  /** - copy pfo->ln_k to ln_q and overwrite the last value with ln(k_split) */
  memcpy(ln_q, pfo->ln_k, split_index*sizeof(double));
  ln_q[split_index] = log(k_split);

  /** - get the nowiggle spectrum at ln_q and z */
  class_call(fourier_pk_nw_at_kvec_and_z(pba, ppm, pfo,
                                        linear,
                                        ln_q,
                                        split_index+1,
                                        z,
                                        pk_nw),
            pfo->error_message,
            pfo->error_message);

  /** - prepare the integrand */
  for (it_q = 0; it_q < split_index+1; it_q++)
  {
    q = exp(ln_q[it_q]);
    q_osc = q / k_bao;
    intg_splines[it_q*index_num + index_y] = q * pk_nw[it_q] * ir_dsigma2_measure(q_osc);
                                              /** q * P_nw(z, q) * j_2(q/k_bao) */
  }

  /** - spline the integrand */
  class_call(array_spline_table_line_to_line(ln_q,
                                            split_index+1,
                                            intg_splines,
                                            index_num,
                                            index_y,
                                            index_ddy,
                                            _SPLINE_EST_DERIV_,
                                            pfo->error_message),
            pfo->error_message,
            pfo->error_message);

  /** - integrate the spline */
  class_call(array_integrate_all_spline_table_line_to_line(ln_q,
                                                          split_index+1,
                                                          intg_splines,
                                                          index_num,
                                                          index_y,
                                                          index_ddy,
                                                          &result,
                                                          pfo->error_message),
            pfo->error_message,
            pfo->error_message);

  // for (it_q = 0; it_q < split_index+1; it_q++)
  //   fprintf(stderr, "%.15e  %.15e  %.15e \n", exp(ln_q[it_q]), intg_splines[it_q*index_num + index_y], intg_splines[it_q*index_num + index_ddy]);

  /** - multiply with prefactor */
  result *= 1./(2.*_PI_*_PI_);

  free(intg_splines);
  free(ln_q);
  free(pk_nw);

  return result;
}

static int indexed_real_arg_cmp_k(const void * a, const void * b) {
  struct indexed_real_arg * a_ = (struct indexed_real_arg *)a;
  struct indexed_real_arg * b_ = (struct indexed_real_arg *)b;
  if ((*a_).ln_k < (*b_).ln_k)
      return -1;
  else if ((*a_).ln_k > (*b_).ln_k)
      return 1;
  else
      return 0;
}

int eft_real_argument_list_rect(const double * const ln_kvec,
                                const int k_size,
                                const int n_columns,
                                struct indexed_real_arg ** vec,
                                ErrorMsg errmsg) {

  int it_k, it_mu;
  
  /** - allocate the output list */
  class_alloc(*vec, n_columns*k_size*sizeof(struct indexed_real_arg), errmsg);

  /** - transpose ln_kvec array to presort */
  for (it_k = 0; it_k < k_size; it_k++) {
    for (it_mu = 0; it_mu < n_columns; it_mu++) {
      (*vec)[it_k*n_columns + it_mu].index = it_mu*k_size + it_k;  /**< index in ln_kvec */
      (*vec)[it_k*n_columns + it_mu].ln_k  = ln_kvec[it_mu*k_size + it_k]; /**< corresponding ln(k) value */
    }
  }
  /** - sort using quicksort (fourier_pk functions require ln_kvec to be sorted in ascending order) */
  qsort(*vec, n_columns*k_size, sizeof(struct indexed_real_arg), indexed_real_arg_cmp_k);

  return _SUCCESS_;
}

/**
 * @brief Compute the leading-order IR-resummed matter power spectrum, ala Ivanovic et al.
 *
 * @param pba       Input: pointer to background structure
 * @param ppm       Input: pointer to primordial structure
 * @param pfo       Input: pointer to fourier structure
 * @param mode      Input: linear or logarithmic
 * @param ln_kvec   Input: array of logarithmic wavenumbers in ascending order (in 1/Mpc)
 * @param kvec_size Input: size of array of wavenumbers
 * @param z         Input: redshift
 * @param sigma2    Input: infrared-suppression factor in the exponent at z
 * @param out_pk    Output: the leading-order IR resummed power spectrum
 *
 * @return the error status
 */
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
            double * out_pk) {
  
  int it;
  double pk_w, k;
  double *ln_kvec_sorted, *pk_lin, *pk_nw;

  class_alloc(ln_kvec_sorted, vec_size*sizeof(double), pfo->error_message);
  class_alloc(pk_lin, vec_size*sizeof(double), pfo->error_message);
  class_alloc(pk_nw,  vec_size*sizeof(double), pfo->error_message);

  for (it = 0; it < vec_size; it++) {
    ln_kvec_sorted[it] = vec[it].ln_k;
  }

  /** - get the linear power spectrum at z */
  class_call(fourier_pk_at_kvec_and_z(pba, ppm, pfo,
                                      linear,
                                      pk_linear,
                                      ln_kvec_sorted,
                                      vec_size,
                                      z,
                                      linear_spectrum_index,
                                      pk_lin),
              pfo->error_message,
              pfo->error_message);

  /** - get the nowiggle power spectrum at z */
  class_call(fourier_pk_nw_at_kvec_and_z(pba, ppm, pfo,
                                        linear,
                                        ln_kvec_sorted,
                                        vec_size,
                                        z,
                                        pk_nw),
             pfo->error_message,
             pfo->error_message);
  
  for (it = 0; it < vec_size; it++)
  {
    /** - compute the IR resummed spectrum */
    k = exp( vec[it].ln_k );
    pk_w = pk_lin[it] - pk_nw[it];
    out_pk[vec[it].index] = pk_nw[it] + exp(-k*k * sigma2_ir_at_z) * pk_w;
  }

  free(ln_kvec_sorted);
  free(pk_lin);
  free(pk_nw);

  /** - convert to logarithmic output if needed */
  if (mode == logarithmic) {
    for (it = 0; it < vec_size; it++)
      out_pk[it] = log( out_pk[it] );
  }

  return _SUCCESS_;
}


/**
 * @brief Compute the next-to-leading-order IR-resummed matter power spectrum, ala Ivanovic et al.
 *
 * @param pba       Input: pointer to background structure
 * @param ppm       Input: pointer to primordial structure
 * @param pfo       Input: pointer to fourier structure
 * @param mode      Input: linear or logarithmic
 * @param ln_kvec   Input: array of logarithmic wavenumbers in ascending order (in 1/Mpc)
 * @param kvec_size Input: size of array of wavenumbers
 * @param z         Input: redshift
 * @param sigma2    Input: infrared-suppression factor in the exponent at z
 * @param out_pk    Output: the leading-order IR resummed power spectrum
 *
 * @return the error status
 */
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
            double * out_pk) {
  
  int it;
  double pk_w, k;
  double *ln_kvec_sorted, *pk_lin, *pk_nw;

  class_alloc(ln_kvec_sorted, vec_size*sizeof(double), pfo->error_message);
  class_alloc(pk_lin, vec_size*sizeof(double), pfo->error_message);
  class_alloc(pk_nw,  vec_size*sizeof(double), pfo->error_message);

  for (it = 0; it < vec_size; it++) {
    ln_kvec_sorted[it] = vec[it].ln_k;
  }

  /** - get the linear power spectrum at z */
  class_call(fourier_pk_at_kvec_and_z(pba, ppm, pfo,
                                      linear,
                                      pk_linear,
                                      ln_kvec_sorted,
                                      vec_size,
                                      z,
                                      linear_spectrum_index,
                                      pk_lin),
              pfo->error_message,
              pfo->error_message);

  /** - get the nowiggle power spectrum at z */
  class_call(fourier_pk_nw_at_kvec_and_z(pba, ppm, pfo,
                                        linear,
                                        ln_kvec_sorted,
                                        vec_size,
                                        z,
                                        pk_nw),
             pfo->error_message,
             pfo->error_message);
  
  for (it = 0; it < vec_size; it++)
  {
    /** - compute the IR resummed spectrum */
    k = exp( vec[it].ln_k );
    pk_w = pk_lin[it] - pk_nw[it];
    out_pk[vec[it].index] = pk_nw[it] + exp(-k*k * sigma2_ir_at_z) * pk_w * (1. + k*k * sigma2_ir_at_z);
  }

  free(ln_kvec_sorted);
  free(pk_lin);
  free(pk_nw);

  /** - convert to logarithmic output if needed */
  if (mode == logarithmic) {
    for (it = 0; it < vec_size; it++)
      out_pk[it] = log( out_pk[it] );
  }

  return _SUCCESS_;
}

static int indexed_rsd_arg_cmp_k(const void * a, const void * b) {
  struct indexed_rsd_arg * a_ = (struct indexed_rsd_arg *)a;
  struct indexed_rsd_arg * b_ = (struct indexed_rsd_arg *)b;
  if ((*a_).ln_k < (*b_).ln_k)
      return -1;
  else if ((*a_).ln_k > (*b_).ln_k)
      return 1;
  else
      return 0;
}

/**
 * @brief Compiles a list of sorted and indexed arguments for IR-resummed linear spectra in RSD space.
 *        The argument list will contain the tensor product of the input k- and mu-lists.
 * @param ln_kvec   Input: array of logarithmic wavenumbers (in 1/Mpc)
 * @param k_size    Input: size of the wavenumber array
 * @param muvec     Input: array of line-of-sight angles w.r.t. the RSD direction (cos(theta))
 * @param mu_size   Input: size of the line-of-sight angle array
 * @param vec       Output: argument list
 * @param errmsg 
 * 
 * @return the error status
 */
int eft_rsd_argument_list_rect(const double * const ln_kvec,
                               const int k_size,
                               const double * const muvec,
                               const int mu_size,
                               struct indexed_rsd_arg ** vec,
                               ErrorMsg errmsg) {

  int it_k, it_mu;
  
  /** - allocate the output list */
  class_alloc(*vec, mu_size*k_size*sizeof(struct indexed_rsd_arg), errmsg);

  /** - transpose ln_kvec array to presort */
  for (it_k = 0; it_k < k_size; it_k++) {
    for (it_mu = 0; it_mu < mu_size; it_mu++) {
      (*vec)[it_k*mu_size + it_mu].index = it_mu*k_size + it_k;  /**< index in ln_kvec */
      (*vec)[it_k*mu_size + it_mu].ln_k  = ln_kvec[it_mu*k_size + it_k]; /**< corresponding ln(k) value */
      (*vec)[it_k*mu_size + it_mu].mu    = muvec[it_mu]; /**< corresponding mu value */
    }
  }
  /** - sort using quicksort (fourier_pk functions require ln_kvec to be sorted in ascending order) */
  qsort(*vec, mu_size*k_size, sizeof(struct indexed_rsd_arg), indexed_rsd_arg_cmp_k);

  return _SUCCESS_;
}

int eft_rsd_argument_list(const double * const ln_kvec,
                          const double * const muvec,
                          const int size,
                          struct indexed_rsd_arg ** vec,
                          ErrorMsg errmsg) {

  int it;
  
  /** - allocate the output list */
  class_alloc(*vec, size*sizeof(struct indexed_rsd_arg), errmsg);

  /** - transpose ln_kvec array to presort */
  for (it = 0; it < size; it++) {
    (*vec)[it].index = it;  /**< index in ln_kvec */
    (*vec)[it].ln_k  = ln_kvec[it]; /**< corresponding ln(k) value */
    (*vec)[it].mu    = muvec[it]; /**< corresponding mu value */
  }
  /** - sort using quicksort (fourier_pk functions require ln_kvec to be sorted in ascending order) */
  qsort((void *)*vec, size, sizeof(struct indexed_rsd_arg), indexed_rsd_arg_cmp_k);

  return _SUCCESS_;
}

/**
 * @brief Compute the leading-order IR-resummed matter power spectrum in RSD space, ala Ivanovic et al.;
 *        assumes ln_kvec is sorted in ascending order
 *
 * @param pba       Input: pointer to background structure
 * @param ppm       Input: pointer to primordial structure
 * @param pfo       Input: pointer to fourier structure
 * @param mode      Input: linear or logarithmic
 * @param ln_kvec   Input: array of logarithmic wavenumbers in ascending order (in 1/Mpc)
 * @param muvec     Input: array of line-of-sight angles w.r.t. the RSD direction (cos(theta))
 * @param vec_size  Input: size of array of wavenumbers and angles
 * @param z         Input: redshift
 * @param f         Input: growth factor f(z)
 * @param sigma2_ir_at_z    Input: infrared-suppression factor in the exponent at z
 * @param dsigma2_ir_at_z   Input: next order of infrared-suppression factor
 * @param out_pk    Output: the leading-order IR-resummed RSD power spectrum
 *                          indexed as out_pk[index_point]
 *
 * @return the error status
 */
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
            double * out_pk) {
  
  int it;
  double pk_w, k, mu, sigma2_tot;
  double *ln_kvec_sorted, *pk_lin, *pk_nw;

  class_alloc(ln_kvec_sorted, vec_size*sizeof(double), pfo->error_message);
  class_alloc(pk_lin, vec_size*sizeof(double), pfo->error_message);
  class_alloc(pk_nw,  vec_size*sizeof(double), pfo->error_message);

  for (it = 0; it < vec_size; it++) {
    ln_kvec_sorted[it] = vec[it].ln_k;
  }

  /** - get the linear power spectrum at z */
  class_call(fourier_pk_at_kvec_and_z(pba, ppm, pfo,
                                      linear,
                                      pk_linear,
                                      ln_kvec_sorted,
                                      vec_size,
                                      z,
                                      linear_spectrum_index,
                                      pk_lin),
              pfo->error_message,
              pfo->error_message);

  /** - get the nowiggle power spectrum at z */
  class_call(fourier_pk_nw_at_kvec_and_z(pba, ppm, pfo,
                                        linear,
                                        ln_kvec_sorted,
                                        vec_size,
                                        z,
                                        pk_nw),
            pfo->error_message,
            pfo->error_message);
  
  for (it = 0; it < vec_size; it++) {
    /** - compute the total mu-dependent suppression factor */
    k = exp( vec[it].ln_k );
    mu = vec[it].mu;
    sigma2_tot = (1. + f*mu*mu*(2. + f)) * sigma2_ir_at_z     \
                + f*f*mu*mu*(mu*mu - 1.) * dsigma2_ir_at_z;
    /** - compute the IR resummed spectrum */
    pk_w = pk_lin[it] - pk_nw[it];
    out_pk[vec[it].index] = pk_nw[it] + exp(-k*k * sigma2_tot) * pk_w;
  }

  free(ln_kvec_sorted);
  free(pk_lin);
  free(pk_nw);

  /** - convert to logarithmic output if needed */
  if (mode == logarithmic) {
    for (it = 0; it < vec_size; it++) {
      out_pk[it] = log( out_pk[it] );
    }
  }

  return _SUCCESS_;
}

/**
 * @brief Compute the next-to-leading-order IR-resummed matter power spectrum in RSD space, ala Ivanovic et al.
 *
 * @param pba       Input: pointer to background structure
 * @param ppm       Input: pointer to primordial structure
 * @param pfo       Input: pointer to fourier structure
 * @param mode      Input: linear or logarithmic
 * @param ln_kvec   Input: array of logarithmic wavenumbers in ascending order (in 1/Mpc)
 * @param kvec_size Input: size of array of wavenumbers
 * @param z         Input: redshift
 * @param f         Input: growth factor f(z)
 * @param muvec     Input: array of line-of-sight angles w.r.t. the RSD direction (cos(theta))
 * @param mu_size   Input: size of array of angles
 * @param sigma2_ir_at_z    Input: infrared-suppression factor in the exponent at z
 * @param dsigma2_ir_at_z   Input: next order of infrared-suppression factor
 * @param out_pk    Output: the next-to-leading-order IR-resummed RSD power spectrum
 *                          indexed as out_pk[index_mu][index_k]
 *
 * @return the error status
 */
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
            double * out_pk) {
  
  int it;
  double pk_w, k, mu, sigma2_tot;
  double *ln_kvec_sorted, *pk_lin, *pk_nw;

  class_alloc(ln_kvec_sorted, vec_size*sizeof(double), pfo->error_message);
  class_alloc(pk_lin, vec_size*sizeof(double), pfo->error_message);
  class_alloc(pk_nw,  vec_size*sizeof(double), pfo->error_message);

  for (it = 0; it < vec_size; it++) {
    ln_kvec_sorted[it] = vec[it].ln_k;
  }

  /** - get the linear power spectrum at z */
  class_call(fourier_pk_at_kvec_and_z(pba, ppm, pfo,
                                      linear,
                                      pk_linear,
                                      ln_kvec_sorted,
                                      vec_size,
                                      z,
                                      linear_spectrum_index,
                                      pk_lin),
              pfo->error_message,
              pfo->error_message);

  /** - get the nowiggle power spectrum at z */
  class_call(fourier_pk_nw_at_kvec_and_z(pba, ppm, pfo,
                                        linear,
                                        ln_kvec_sorted,
                                        vec_size,
                                        z,
                                        pk_nw),
            pfo->error_message,
            pfo->error_message);
  
  for (it = 0; it < vec_size; it++) {
    /** - compute the total mu-dependent suppression factor */
    k = exp( vec[it].ln_k );
    mu = vec[it].mu;
    sigma2_tot = (1. + f*mu*mu*(2. + f)) * sigma2_ir_at_z     \
                + f*f*mu*mu*(mu*mu - 1.) * dsigma2_ir_at_z;
    /** - compute the IR resummed spectrum */
    pk_w = pk_lin[it] - pk_nw[it];
    out_pk[vec[it].index] = pk_nw[it] + exp(-k*k * sigma2_tot) * pk_w * (1. + k*k * sigma2_tot);
  }

  free(ln_kvec_sorted);
  free(pk_lin);
  free(pk_nw);

  /** - convert to logarithmic output if needed */
  if (mode == logarithmic) {
    for (it = 0; it < vec_size; it++) {
      out_pk[it] = log( out_pk[it] );
    }
  }

  return _SUCCESS_;
}

// TODO: not needed anymore
// int eft_necessary_pk_moments(
//               double * biases,
//               int biases_size,
//               int * moments_n,
//               int * moments_size,
//               ErrorMsg errmsg) {

//   int j;
//   double min_bias, max_bias;

//   /** - find the minimum bias */
//   min_bias = biases[0];
//   for (j = 0; j < biases_size; j++) {
//     if (biases[j] < min_bias) { min_bias = biases[j]; }
//   }
//   /** - find the maximum bias */
//   max_bias = biases[0];
//   for (j = 0; j < biases_size; j++) {
//     if (biases[j] > max_bias) { max_bias = biases[j]; }
//   }

//   /** TODO: find the function momentum range dependent on min- and max-bias 
//    *  currently just fixed for the preset biases: sigma_0^2 and sigma_(-1)^2 = 3 sigma_v^2
//   */
  
//   *moments_size = 2;
//   class_alloc(moments_n, (*moments_size)*sizeof(int), errmsg);
//   moments_n[0] = -1;
//   moments_n[1] = 0;

//   return _SUCCESS_;
// }

// /**
//  * @brief Compute moments of the different power spectra as
//  *        d(s_n^2(z)) = 1/(2 pi^2) dln(q) q^(2(n+1)+1) P_xx(q)
//  *        and compute the shot-noise contribution of Idelta2delta200 if peft is given
//  *
//  * @param ppr         Input: pointer to precision structure
//  * @param pba         Input: pointer to background structure
//  * @param ppm         Input: pointer to primordial structure
//  * @param pfo         Input: pointer to fourier structure
//  * @param hp          Input: hyperparameters necessary for the computation
//  * @param peft        Input: pointer to eft structure (can be set to NULL if only the moments are wanted)
//  * @param biases      Input: array of FFTLog biases (only read if peft is given)
//  * @param biases_size Input: size of the bias array
//  * @param z           Input: redshift
//  * @param sigma2_ir_at_z  Input: exponential IR suppression factor (only read if use_pk_type is an IR-resummed spectrum)
//  * @param use_pk_type Input: sets the type of linear power spectrum {pk_lin, pk_nowiggle, pk_ir_resummed_lo, pk_ir_resummed_nlo}
//  * @param moments     Output: array of computed moments of the linear power spectrum (if peft is not given, has to be pre-allocated)
//  * @param moments_n   In/Output: array of n of the computed moments (will be set if peft is given, otherwise has to be filled)
//  * @param moments_size  In/Output: size of the moments array (will be set if peft is given, otherwise has to be set before)
//  *
//  * @return the error status
//  */
// int eft_compute_pk_moments(
//               struct precision * ppr,
//               struct background * pba,
//               struct primordial * ppm,
//               struct fourier * pfo, 
//               const struct eft_hyper_parameters hp,
//               struct eft * peft,
//               double * biases,
//               int biases_size,
//               const double z, 
//               const int linear_spectrum_index,
//               const double sigma2_ir_at_z,
//               enum eft_pk_type use_pk_type,
//               double * moments,
//               int * moments_n,
//               int * moments_size) {

//   int index_y, index_ddy, intg_size = 0, it_q, j;
//   double q, result;
//   double *ln_k, *pk, *intg_splines;

//   if (peft) {
//     class_call(eft_necessary_pk_moments(biases, biases_size, moments_n, moments_size, pfo->error_message),
//                 pfo->error_message, pfo->error_message);
//     class_alloc(moments, (*moments_size)*sizeof(double), pfo->error_message);
//   }
  
//   class_define_index(index_y, _TRUE_, intg_size, ppr->eft_pk_moments_points);
//   class_define_index(index_ddy, _TRUE_, intg_size, ppr->eft_pk_moments_points);

//   class_alloc(ln_k, ppr->eft_pk_moments_points*sizeof(double),
//               pfo->error_message);
//   class_alloc(pk,   ppr->eft_pk_moments_points*sizeof(double),
//               pfo->error_message);
//   class_alloc(intg_splines, intg_size*sizeof(double),
//               pfo->error_message);

//   /** - get sample points for the spline of P_lin against ln(k) in the range of [k_IR, k_UV] */
//   eft_spline_sample_points_nonuniform(hp.k_IR_cutoff, exp(ppr->k_bao_center), hp.k_UV_cutoff,         \
//                                       hp.bao_oversampling, hp.ln_k_oversampling_width,                \
//                                       ppr->eft_pk_moments_points, ln_k, pfo->fourier_verbose);

//   /** - get the chosen spectrum at ln_k and z */
//   switch (use_pk_type)
//   {
//   case pk_linear:
//     class_call(fourier_pk_at_kvec_and_z(pba, ppm, pfo,
//                                         linear,
//                                         pk_linear,
//                                         ln_k,
//                                         ppr->eft_pk_moments_points,
//                                         z,
//                                         linear_spectrum_index,
//                                         pk),
//                 pfo->error_message,
//                 pfo->error_message);
//     break;
  
//   case pk_nowiggle:
//     class_call(fourier_pk_nw_at_kvec_and_z(pba, ppm, pfo,
//                                           linear,
//                                           ln_k,
//                                           ppr->eft_pk_moments_points,
//                                           z,
//                                           pk),
//                 pfo->error_message,
//                 pfo->error_message);
//     break;

//   case pk_ir_resummed_lo:
//     class_call(eft_ir_pk_lo(pba, ppm, pfo,
//                             linear,
//                             ln_k,
//                             ppr->eft_pk_moments_points,
//                             z,
//                             sigma2_ir_at_z,
//                             pk),
//                 pfo->error_message,
//                 pfo->error_message);
//     break;

//   case pk_ir_resummed_nlo:
//     class_call(eft_ir_pk_nlo(pba, ppm, pfo,
//                             linear,
//                             ln_k,
//                             ppr->eft_pk_moments_points,
//                             z,
//                             sigma2_ir_at_z,
//                             pk),
//                 pfo->error_message,
//                 pfo->error_message);
//     break;

//   default:
//     free(intg_splines);
//     free(pk);
//     free(ln_k);
//     class_stop(pfo->error_message, "Instructed to use pk_type = %d which is not supported in eft_pk_moments", use_pk_type);
//     break;
//   }

//   /** TODO: could be parallelized */
//   for (j = 0; j < *moments_size; j++) {
//     /** - prepare the integrand */
//     for (it_q = 0; it_q < ppr->eft_pk_moments_points; it_q++)
//     {
//       intg_splines[index_y + it_q] = exp((2.*((double)moments_n[j] + 1.) + 1.) * ln_k[it_q]) * pk[it_q];
//     }

//     /** - spline the integrand */
//     class_call(array_spline_table_one_column(ln_k,
//                                              ppr->eft_pk_moments_points,
//                                              intg_splines + index_y,
//                                              0,
//                                              intg_splines + index_ddy,
//                                              _SPLINE_EST_DERIV_,
//                                              pfo->error_message),
//                 pfo->error_message,
//                 pfo->error_message);

//     /** - integrate the spline */
//     class_call(array_integrate_all_spline_table_lines(ln_k,
//                                                       ppr->eft_pk_moments_points,
//                                                       intg_splines + index_y,
//                                                       1,
//                                                       intg_splines + index_ddy,
//                                                       &result,
//                                                       pfo->error_message),
//                 pfo->error_message,
//                 pfo->error_message);

//     /** - multiply with prefactor */
//     moments[j] = 1./(2.*_PI_*_PI_) * result;
//   }

//   /** - if eft struct is given, also compute shot-noise part of Idelta2delta200, 
//    *    assumes that spectra_contributions is allocated for the given pk_type */
//   if (peft) {
//     /** - prepare the integrand */
//     for (it_q = 0; it_q < ppr->eft_pk_moments_points; it_q++)
//     {
//       intg_splines[index_y + it_q] = exp(3./2. * ln_k[it_q]) * pk[it_q];
//     }

//     /** - spline the integrand */
//     class_call(array_spline_table_one_column(ln_k,
//                                              ppr->eft_pk_moments_points,
//                                              intg_splines + index_y,
//                                              0,
//                                              intg_splines + index_ddy,
//                                              _SPLINE_EST_DERIV_,
//                                              pfo->error_message),
//                 pfo->error_message,
//                 pfo->error_message);

//     /** - integrate the spline */
//     class_call(array_square_integrate_all_spline_table_lines(ln_k,
//                                                              ppr->eft_pk_moments_points,
//                                                              intg_splines + index_y,
//                                                              1,
//                                                              intg_splines + index_ddy,
//                                                              &result,
//                                                              pfo->error_message),
//                 pfo->error_message,
//                 pfo->error_message);

//     /** - multiply with prefactor and subtract from Idelta2delta200 */
//     result *= 4.*_PI_;
//     for (it_q = 0; it_q < peft->k_size; it_q++) {
//       peft->spectra_contributions[use_pk_type][peft->index_Idelta2delta200][it_q] -= result;
//     }
//   }

//   free(intg_splines);
//   free(pk);
//   free(ln_k);

//   return _SUCCESS_;
// }





/**
 * Compute the leading-order IR-resummed matter power spectrum, ala Ivanovic et al
 *
 * @param Cx           Input: pointer to cosmology structure
 * @param k            Input: wavenumber in unit of 1/Mpc.
 * @param z            Input: redshift
 * @param SPLIT        Input: switch to set the method of wiggle-nowiggle split
 *
 * @return value of leading IR-ressumed power spectrum
 */
double pm_IR_LO(struct background * pba,
                struct primordial * ppm,
                struct fourier * pfo,
                double k,
                double z,
                long SPLIT)
{
    double kf0 = 1.e-4;
    double kmin = 2.e-4;
    double kmax = 15.; // ??
    static double sig2_LO = - 1.;

    if(sig2_LO == - 1.){
          sig2_LO = IR_Sigma2(pba, ppm, pfo, z, kf0, SPLIT);
    }
    //fprintf(stderr,"Call pm_nowiggle from pm_IR_LO with k=%e, z=%e, cleanup=%d\n",k,z,cleanup);
    double p_nowiggle = pm_nowiggle(pba, ppm, pfo, k, z, kf0, 0, SPLIT);
    double plin       = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    double p_wiggle   = plin - p_nowiggle;
    double sup        = exp(-k * k * sig2_LO);
    double f;

    //fprintf(stderr, "%.15e  %.15e  %.15e \n", kmin, kmax, sig2_LO);
    //fprintf(stderr, "%.15e  %.15e \n", k, p_nowiggle);

    // if (k < kmax && k > kmin){
    //   f = p_nowiggle + sup * p_wiggle;
    // }
    // else{
    //   f = p_nowiggle + p_wiggle;
    // }

    f = p_nowiggle + sup * p_wiggle;

    // FILE *fpa;
    // char file_name[50];
    // sprintf(file_name, "NOIRvsWIR.txt");
    // fpa = fopen(file_name, "a");
    // fprintf(fpa, "%e %e %e %e %e\n", k, p_nowiggle, p_wiggle, plin, f);
    // fclose(fpa);

    return f;
}

/**
 * Compute the next-to-leading-order IR-resummed matter power spectrum, ala Ivanovic et al
 *
 * @param Cx           Input: pointer to cosmology structure
 * @param k            Input: wavenumber in unit of 1/Mpc.
 * @param z            Input: redshift
 * @param SPLIT        Input: switch to set the method of wiggle-nowiggle split
 *
 * @return value of Tree IR-ressumed power spectrum
 */
double pm_IR_NLO(struct background * pba,
                struct primordial * ppm,
                struct fourier * pfo,
                double k,
                double z,
                long SPLIT)
{
    double kf0 = 1.e-4;
    double kmin = 2.e-4;
    double kmax = 15.;
    static double sig2_LO = - 1.;

    if(sig2_LO == - 1.){
          sig2_LO = IR_Sigma2(pba, ppm, pfo, z, kf0, SPLIT);
    }
    //fprintf(stderr,"Call pm_nowiggle from pm_IR_LO with k=%e, z=%e, cleanup=%d\n",k,z,cleanup);
    double p_nowiggle = pm_nowiggle(pba, ppm, pfo, k, z, kf0, 0, SPLIT);
    double plin       = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    double p_wiggle   = plin - p_nowiggle;
    double sup        = exp(-k * k * sig2_LO);

    double f     = p_nowiggle + sup * p_wiggle * (1. + k * k * sig2_LO);
    //fprintf(stderr, "%e %e ",p22_IR, p13_IR);
    return f;
}



/**
 * Integrand to compute the suppression factor IR_sigma2
 *
 * @param x            Input: integration variable, k-values
 * @param par          Input: integration parameters
 *
 * @return integrand to be used in IR_sigma2() function
 */
double IR_Sigma2_integrand(double x, void *par)
{
    double result = 0;

    struct integrand_parameters2 pij;
    pij = *((struct integrand_parameters2 *)par);
    struct background *pba = pij.pba;
    struct primordial *ppm = pij.ppm;
    struct fourier *pfo    = pij.pfo;
    double bao_scale       = pij.p4;
    double z               = pij.p5;
    double kf0             = pij.p6;
    long   SPLIT           = pij.p13;

    double q_osc = x * bao_scale;  // = q/k_osc , BAO_scale = 110. Mpc/h.
    //fprintf(stderr,"Call pm_nowiggle from IR_Sigma2_integrand with x=%e\n",x);
    result = 1./(6.*M_PI*M_PI) * pm_nowiggle(pba, ppm, pfo, x, z, kf0, 0, SPLIT) 
            * (1. - sin(q_osc) / q_osc 
                + 2. * ((3. / (q_osc*q_osc) - 1.) * sin(q_osc) / q_osc - 3. * cos(q_osc) / (q_osc*q_osc) ) );
              // 1 - j_0(q/k_osc) + 2 j_2(q/k_osc)
    return result;

}


/**
 * Compute the suppression factor IR_sigma2 in real space
 *
 * @param Cx           Input: pointer to cosmology structure
 * @param z            Input: redshift
 * @param kf0          Input: first element of the k-array, used in normalization of EH no-wiggle spectrum
 * @param SPLIT        Input: switch to set the method of wiggle-nowiggle split
 *
 * @return value of IR resummation suppression factor
 */
double IR_Sigma2(struct background * pba,
                struct primordial * ppm,
                struct fourier * pfo,
                double z,
                double kf0,
                long SPLIT)
{
    //extern struct globals gb;
    double result=0., error=0.;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);

    struct integrand_parameters2 par;

    double kmin = 1.e-6;
    double kmax = 0.2;

    gsl_function F;
    F.function = &IR_Sigma2_integrand;
    F.params = &par;

    par.ppm = ppm;
    par.pba = pba;
    par.pfo = pfo;
    par.p4  = 110. / pba->h;
    par.p5  = z;
    par.p6  = kf0;
    par.p13 = SPLIT;
    gsl_integration_qags(&F,kmin,kmax,0.0,1.0e-3,1000000,w,&result,&error);
    gsl_integration_workspace_free(w);

    return result;
}

/**
 * Integrand to compute the suppression factor delta_IR_sigma2 for RSD
 *
 * @param x            Input: integration variable, k-values
 * @param par          Input: integration parameters
 *
 * @return integrand to be used in IR_sigma2() function
 */
double IR_del_Sigma2_integrand(double x, void *par)
{
    double result = 0;

    struct integrand_parameters2 pij;
    pij = *((struct integrand_parameters2 *)par);
    struct background *pba = pij.pba;
    struct primordial *ppm = pij.ppm;
    struct fourier *pfo    = pij.pfo;
    double bao_scale       = pij.p4;
    double z               = pij.p5;
    double kf0             = pij.p6;
    long   SPLIT           = pij.p13;

    double q_osc = x * bao_scale;  // = q/k_osc , BAO_scale = 110. Mpc/h.
    //fprintf(stderr,"Call pm_nowiggle from IR_del_Sigma2_integrand with x=%e\n",x);
    result = 1./(2.*M_PI*M_PI) * pm_nowiggle(pba, ppm, pfo, x, z, kf0, 0, SPLIT) 
            * ((3. / (q_osc*q_osc) - 1.) * sin(q_osc) / q_osc - 3. * cos(q_osc) / (q_osc*q_osc) );
            // j_2(q/k_osc)
    return result;

}

/**
 * Compute the suppression factor IR_sigma2 in redshift space (CLASS-PT: 2004.106007v2 p. 15)
 *
 * @param Cx           Input: pointer to cosmology structure
 * @param z            Input: redshift
 * @param kf0          Input: first element of the k-array, used in normalization of EH no-wiggle spectrum
 * @param SPLIT        Input: switch to set the method of wiggle-nowiggle split
 *
 * @return value of IR resummation suppression factor
 */
double IR_del_Sigma2(struct background * pba,
                struct primordial * ppm,
                struct fourier * pfo,
                double z,
                double kf0,
                long SPLIT)
{
    //extern struct globals gb;
    double del_sigma_2=0., error=0. ;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);

    struct integrand_parameters2 par;

    double kmin = 1.e-6;
    double kmax = 0.2;

    double f = 1.; // Needs to be substituted with a CLASS call obtaining the logarithmic growth

    gsl_function F;
    F.function = &IR_del_Sigma2_integrand;
    F.params = &par;

    par.ppm = ppm;
    par.pba = pba;
    par.pfo = pfo;
    par.p4  = 110. / pba->h;
    par.p5  = z;
    par.p6  = kf0;
    par.p13 = SPLIT;
    gsl_integration_qags(&F,kmin,kmax,0.0,1.0e-3,1000000,w,&del_sigma_2,&error);
    gsl_integration_workspace_free(w);

    // double sigma_IR = IR_Sigma2(pba, ppm, pfo, z, kf0, SPLIT);

    // double result = (1 + f * pow(mu,2.) * (2 + f)) * sigma_IR + pow(f*mu, 2.) * (pow(mu,2.) - 1) * del_sigma_2;

    return del_sigma_2;
}

/**
 * Compute the no-wiggle componenet of the matter power spectrum
 *
 * @param Cx           Input: pointer to cosmology structure
 * @param k            Input: wavenumber in unit of h/Mpc.
 * @param z            Input: redshift
 * @param kf0          Input: first element of the k-array, used in normalization of EH no-wiggle spectrum
 * @param cleanup      Input: switch to set whether to free the memory allocated to no-wiggle interpolators
 * @param SPLIT        Input: switch to set the method of wiggle-nowiggle split
 *
 * @return double value of no-wiggle power spectrum
 */
double pm_nowiggle(struct background * pba,
                struct primordial * ppm,
                struct fourier * pfo,
                double k,
                double z,
                double kf0, //TODO: unused
                int cleanup,
                long SPLIT)
{
  double pm_nw = 0.;

  if(SPLIT == DST)
    pm_nw = pm_nowiggle_dst(pba, ppm, pfo, k, z, cleanup);
  else if(SPLIT == GFILTER)
    pm_nw = pm_nowiggle_gfilter(pba, ppm, pfo, k, z, cleanup);
  else if(SPLIT == BSPLINE){
    pm_nw = pm_nowiggle_bspline(pba, ppm, pfo, k, z, cleanup);
  }

  return pm_nw;
}


/**
 * Compute the no-wiggle componenet of the matter power spectrum, reading in and interpolating the output of apython code which computed the broadband by fitting families of Bsplines (see Vlah et al 2015)
 *
 * @param Cx           Input: pointer to cosmology structure
 * @param k            Input: wavenumber in unit of h/Mpc.
 * @param z            Input: redshift
 * @param cleanup      Input: switch to set whether to free the memory allocated to no-wiggle interpolators
 *
 * @return double value of no-wiggle power spectrum
 */
double pm_nowiggle_bspline(struct background * pba,
                           struct primordial * ppm,
                           struct fourier * pfo,
                           double k,
                           double z,
                           int cleanup)
{

  static gsl_interp_accel   *pknw_accel_ptr;
  static gsl_spline         *pknw_spline_ptr;
  FILE *pksplit_file;

  int i,j;
  static double kmin=0., kmax=0.;
  static int first = 1;
  if(first == 1){
    char pksplit_filename[FILENAME_MAX];
    sprintf(pksplit_filename,"Input/wnw_split/pk_Bspline.txt");

    int nlines = 0;
    nlines     = count_lines_in_file(pksplit_filename);

    double *k_in, *pk_nw, *log_k, *log_pknw;
    k_in     = make_1Darray(nlines);
    log_k    = make_1Darray(nlines);
    pk_nw    = make_1Darray(nlines);
    log_pknw = make_1Darray(nlines);

    pksplit_file = fopen(pksplit_filename,"r");

    if(pksplit_file==NULL){
      printf("Failed to open the file with k-values");
      exit(1);
    }

    char line[MAXL];
    int err;
    while(fgets(line, sizeof line, pksplit_file) != NULL )
    {
        if(*line == '#')  continue;
          for(i=0;i<nlines;i++){
            err         = fscanf(pksplit_file,"%lg %lg \n",&k_in[i],&pk_nw[i]);
            log_k[i]    = log(k_in[i]);
            log_pknw[i] = log(pk_nw[i]/pow(2.*M_PI,3.));
          }
    }
    fclose(pksplit_file);

    pknw_accel_ptr  = gsl_interp_accel_alloc();
    pknw_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,nlines);

    gsl_spline_init(pknw_spline_ptr,log_k,log_pknw,nlines);

    kmin = k_in[0];
    kmax = k_in[nlines-1];

    free(k_in);
    free(log_k);
    free(pk_nw);
    free(log_pknw);

    first = 0;
  }


  double pknw, logpknw, logk;

  if(k<kmin || k >kmax) {
      pknw = 0.;
  }
  else{
      logk    = log(k);
      logpknw = gsl_spline_eval(pknw_spline_ptr,logk,pknw_accel_ptr);
      pknw    = exp(logpknw);
  }



  double growth2 = pow(growth_D(pba, z),2.);
  double pknw_z = growth2 * pknw;


  if (cleanup == 1){
        gsl_interp_accel_free(pknw_accel_ptr);
        gsl_spline_free(pknw_spline_ptr);
  }


  return pknw_z;
}

/**
 * Compute the no-wiggle componenet of the matter power spectrum, using Gaussian filter (see Vlah et al 2015)
 *
 * @param Cx           Input: pointer to cosmology structure
 * @param k            Input: wavenumber in unit of h/Mpc.
 * @param z            Input: redshift
 * @param cleanup      Input: switch to set whether to free the memory allocated to no-wiggle interpolators
 *
 * @return double value of no-wiggle power spectrum
 */
double pm_nowiggle_gfilter(struct background *pba, struct primordial *ppm, struct fourier *pfo, double k, double z, int cleanup)
{
  double pm_nw = 0.;
  static gsl_interp_accel   *pknw_accel_ptr;
  static gsl_spline         *pknw_spline_ptr;
  FILE *pksplit_file;

  int i,j;
  static double kmin=0., kmax=0.;
  static int first = 1;

  if(first == 1){
    //fprintf(stderr,"Initialises nowiggle interpolator with first=%d\n",first);
    int nlines       = 400;
    double *k_in     = loginit_1Darray(nlines, 1.e-6, 21.); // was 40.
    double *log_k    = make_1Darray(nlines);
    double *pk_nw    = make_1Darray(nlines);
    double *log_pknw = make_1Darray(nlines);

    for(i=0;i<nlines;i++){
        pk_nw[i]    = pk_Gfilter_nw(pba, ppm, pfo, k_in[i], k_in[0], z);
        log_k[i]    = log(k_in[i]);
        log_pknw[i] = log(pk_nw[i]);
        // printf("nw %d %12.6e %12.6e \n",i, k_in[i],pk_nw[i]);
    }


    pknw_accel_ptr  = gsl_interp_accel_alloc();
    pknw_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,nlines);

    gsl_spline_init(pknw_spline_ptr,log_k,log_pknw,nlines);

    kmin = k_in[0];
    kmax = k_in[nlines-1];

    free(k_in);
    free(log_k);
    free(pk_nw);
    free(log_pknw);

    first = 0;
  }

  double pknw, logpknw, logk;

  if(k<kmin || k >kmax) {
      pknw = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
  }
  else{
      logk    = log(k);
      logpknw = gsl_spline_eval(pknw_spline_ptr,logk,pknw_accel_ptr);
      pknw    = exp(logpknw);
  }


  double growth2 = pow(growth_D(pba, z),2.);
  double pknw_z = growth2 * pknw;

  if (cleanup == 1){
    //fprintf(stderr,"Cleaning nowiggle interpolator with cleanup=1\n");
        gsl_interp_accel_free(pknw_accel_ptr);
        gsl_spline_free(pknw_spline_ptr);
        first = 1;
  }


  return pknw_z;
}


/**
 * Compute the no-wiggle componenet of the matter power spectrum, reading in and interpolating the output of apython code which computed the broadband by discrete sin-transform, See Hamann et al 2010.
 *
 * @param Cx           Input: pointer to cosmology structure
 * @param k            Input: wavenumber in unit of h/Mpc.
 * @param z            Input: redshift
 * @param cleanup      Input: switch to set whether to free the memory allocated to no-wiggle interpolators
 *
 * @return double value of no-wiggle power spectrum
 */
double pm_nowiggle_dst(struct background * pba,
                       struct primordial * ppm,
                       struct fourier * pfo,
                       double k,
                       double z,
                       int cleanup)
{
  double pm_nw = 0.;
    static gsl_interp_accel *pknw_accel_ptr;
  static gsl_spline         *pknw_spline_ptr;
  FILE *pksplit_file;

  int i,j;
  static double kmin=0., kmax=0.;
  static int first = 1;
  if(first == 1){
    char pksplit_filename[FILENAME_MAX];
    sprintf(pksplit_filename,"Input/wnw_split/pk_dst.txt");

    int nlines = 0;
    nlines     = count_lines_in_file(pksplit_filename);

    double *k_in, *pk_nw, *log_k, *log_pknw;
    k_in     = make_1Darray(nlines);
    log_k    = make_1Darray(nlines);
    pk_nw    = make_1Darray(nlines);
    log_pknw = make_1Darray(nlines);

    pksplit_file = fopen(pksplit_filename,"r");

    if(pksplit_file==NULL){
      printf("Failed to open the file with k-values");
      exit(1);
    }

    char line[MAXL];
    int err;
    while(fgets(line, sizeof line, pksplit_file) != NULL )
    {
        if(*line == '#')  continue;
          for(i=0;i<nlines;i++){
            err         = fscanf(pksplit_file,"%lg %lg \n",&k_in[i],&pk_nw[i]);
            log_k[i]    = log(k_in[i]);
            log_pknw[i] = log(pk_nw[i]/pow(2.*M_PI,3.));
            //printf("%12.6e %12.6e \n",log_k[i], log_pknw[i]);
          }
    }
    fclose(pksplit_file);

    pknw_accel_ptr  = gsl_interp_accel_alloc();
    pknw_spline_ptr = gsl_spline_alloc(gsl_interp_cspline,nlines);

    gsl_spline_init(pknw_spline_ptr,log_k,log_pknw,nlines);

    kmin = k_in[0];
    kmax = k_in[nlines-1];

    free(k_in);
    free(log_k);
    free(pk_nw);
    free(log_pknw);

    first = 0;
  }


  double pknw, logpknw, logk;

  if(k<kmin || k >kmax) {
      pknw = 0.;
  }
  else{
      logk    = log(k);
      logpknw = gsl_spline_eval(pknw_spline_ptr,logk,pknw_accel_ptr);
      pknw    = exp(logpknw);
  }



  double growth2 = pow(growth_D(pba, z),2.);
  double pknw_z = growth2 * pknw;

  if (cleanup == 1){
        gsl_interp_accel_free(pknw_accel_ptr);
        gsl_spline_free(pknw_spline_ptr);
  }

  return pknw_z;
}

      /**
 * Compute various moments of the matter power spectrum. For all moments, we need an integrand function and an integrator function
 */

double sigman2_integrand(double x, void *p)
{
	double f=0;
	double q = exp(x);

	struct integrand_parameters2 pij;
	pij = *((struct integrand_parameters2 *)p);

  struct background *pba = pij.pba;
  struct primordial *ppm = pij.ppm;
  struct fourier *pfo    = pij.pfo;
  double z 		           = pij.p4;
  long 	 n 		           = pij.p14;
  long   SPLIT           = pij.p15;
  double window  = 1;

  double pm = Pk_dlnPk(pba, ppm, pfo, q, z, LPOWER);
  // double pm = pm_IR_LO(pba, ppm, pfo, q, z, SPLIT);

  f = 1./(2. * pow(M_PI, 2.)) * pow(q, 2. * (n + 1.) + 1.) * pow(window, 2.) * pm;

	return f;
}


double sigman(struct background * pba, struct primordial * ppm, struct fourier * pfo, double z, double k_min, double k_max, long n, long SPLIT)
{
  double result, error;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);


  struct integrand_parameters2 par;
  double logqmin = log(k_min);
  double logqmax = log(k_max);

  gsl_function F;
  F.function = &sigman2_integrand;
  F.params = &par;

  par.ppm = ppm;
  par.pba = pba;
  par.pfo = pfo;
  par.p4  = z;
  par.p14 = n;
  par.p15 = SPLIT;

	gsl_integration_qags(&F,logqmin,logqmax,0.0,1.e-4,1000000,w,&result,&error);
	gsl_integration_workspace_free (w);

	return result;
}
