#include "header.h"

int eft_spline_sample_points_nonuniform(
                    const double k_min,
                    const double k_feature,
                    const double k_max,
                    const double rel_amplitude,
                    const double width,
                    const int num_points,
                    double * ln_k_sample,
                    const short verbosity) {

  int N;
  double density_outer;
  double ln_k_feature, tau_min, tau_max;
  int N_bound1, N_bound2;
  double sqrt_ampl_1, arctan_sqrt_ampl_1;

  ln_k_feature = log(k_feature);
  tau_min = log(k_min) - ln_k_feature; tau_max = log(k_max) - ln_k_feature;
  sqrt_ampl_1 = sqrt(rel_amplitude - 1.); arctan_sqrt_ampl_1 = atan(sqrt_ampl_1);

  if (tau_min >= -width*sqrt_ampl_1) {
    /** - left and right border is out of bounds */
    if (tau_max <= width*sqrt_ampl_1) {
      density_outer = (double)(num_points - 1) / (rel_amplitude*width*(atan(tau_max/width) - atan(tau_min/width)));
      for (N = 0; N < num_points; N++)
        ln_k_sample[N] = ln_k_feature + width * tan( atan(tau_min/width) + (double)N/(density_outer*rel_amplitude*width) );
      
      if (verbosity > 1) {
        printf("Generated %d sample points for the linear power spectrum with\na minimal point density of %.1f and of %.1f per decade at the feature. \n",     \
                num_points, log(10.)*rel_amplitude*density_outer/(1. + pow(fmax(fabs(tau_min), fabs(tau_max)), 2)), log(10.)*rel_amplitude*density_outer);
      }
      return _SUCCESS_;
    }
    /** - left border is out of bounds */
    else {
      density_outer = (double)(num_points - 1) / (tau_max - width*sqrt_ampl_1 + rel_amplitude*width*(arctan_sqrt_ampl_1 - atan(tau_min/width)));
      N_bound2 = (int)ceil(density_outer*rel_amplitude*width*(arctan_sqrt_ampl_1 - atan(tau_min/width)));

      for (N = 0; N < num_points && N < N_bound2; N++)
        ln_k_sample[N] = ln_k_feature + width * tan( atan(tau_min/width) + (double)N/(density_outer*rel_amplitude*width) );
      for (N; N < num_points; N++)
        ln_k_sample[N] = ln_k_feature + width*sqrt_ampl_1 + rel_amplitude*width*(atan(tau_min/width) - arctan_sqrt_ampl_1) + (double)N/density_outer;
    }
  }
  /** - right border is out of bounds */
  else if (tau_max <= width*sqrt_ampl_1) {
    density_outer = (double)(num_points - 1) / (-tau_min - width*sqrt_ampl_1 + rel_amplitude*width*(atan(tau_max/width) + arctan_sqrt_ampl_1));
    N_bound1 = (int)ceil(-density_outer*(tau_min + width*sqrt_ampl_1));

    for (N = 0; N < num_points && N < N_bound1; N++)
      ln_k_sample[N] = ln_k_feature + tau_min + (double)N/density_outer;
    for (N; N < num_points; N++)
      ln_k_sample[N] = ln_k_feature + width * tan( (tau_min + width*sqrt_ampl_1 + (double)N/density_outer)/(rel_amplitude*width) - arctan_sqrt_ampl_1 );
  }
  /** - both borders are inside the region [k_min, k_max] */
  else {
    density_outer = (double)(num_points - 1) / (tau_max - tau_min - 2.*width*sqrt_ampl_1 + 2.*rel_amplitude*width*arctan_sqrt_ampl_1);
    N_bound1 = (int)ceil(density_outer*(-width*sqrt_ampl_1 - tau_min));
    N_bound2 = (int)ceil(density_outer*(-width*sqrt_ampl_1 - tau_min + 2.*rel_amplitude*width*arctan_sqrt_ampl_1));
    
    for (N = 0; N < num_points && N < N_bound1; N++)
      ln_k_sample[N] = ln_k_feature + tau_min + (double)N/density_outer;
    for (N; N < num_points && N < N_bound2; N++)
      ln_k_sample[N] = ln_k_feature + width * tan( (tau_min + width*sqrt_ampl_1 + (double)N/density_outer)/(rel_amplitude*width) - arctan_sqrt_ampl_1 );
    for (N; N < num_points; N++)
      ln_k_sample[N] = ln_k_feature + tau_min + 2.*width*sqrt_ampl_1 - 2.*rel_amplitude*width*arctan_sqrt_ampl_1 + (double)N/density_outer;
  }
  
  if (verbosity > 1) {
    printf("Generated %d sample points for the linear power spectrum with\na base point density of %.1f and of %.1f per decade at the feature. \n",     \
            num_points, log(10.)*density_outer, log(10.)*rel_amplitude*density_outer);
  }

  return _SUCCESS_;
}


int eft_indices(struct eft * peft,
                struct ext_storage * pext,
                const short eft_role,
                const int index) {

  int i = 0, j;

  if (eft_role == eft_master) { /** - assign indices and allocate all shared quantities */
    /** - define all indices for the different spectra contributions (guaranteed to be consecutive) */
    /** Real space moments */
    class_define_index(peft->index_I2200,           _TRUE_, i, 1);
    class_define_index(peft->index_I1300,           _TRUE_, i, 1);
    class_define_index(peft->index_Idelta200,       _TRUE_, i, 1);
    class_define_index(peft->index_IG200,           _TRUE_, i, 1);
    class_define_index(peft->index_Idelta2delta200, _TRUE_, i, 1);
    class_define_index(peft->index_IG2G200,         _TRUE_, i, 1);
    class_define_index(peft->index_Idelta2G200,     _TRUE_, i, 1);
    class_define_index(peft->index_FG200,           _TRUE_, i, 1);
    
    if (peft->hp->has_rsd) {
      /** - 1-st moment of RSD expansion */
      class_define_index(peft->index_I2201,     _TRUE_, i, 1);
      class_define_index(peft->index_Idelta201, _TRUE_, i, 1);
      class_define_index(peft->index_IG201,     _TRUE_, i, 1);
      class_define_index(peft->index_J21101,    _TRUE_, i, 1);
      class_define_index(peft->index_Jdelta201, _TRUE_, i, 1);
      class_define_index(peft->index_JG201,     _TRUE_, i, 1);
      class_define_index(peft->index_FG201,     _TRUE_, i, 1);
      class_define_index(peft->index_I1301p3101,_TRUE_, i, 1);
      class_define_index(peft->index_J12101,    _TRUE_, i, 1);

      /** - 2-nd moment of RSD expansion */
      class_define_index(peft->index_J21102x,   _TRUE_, i, 1);
      class_define_index(peft->index_J21102y,   _TRUE_, i, 1);
      class_define_index(peft->index_Jdelta202x,_TRUE_, i, 1);
      class_define_index(peft->index_Jdelta202y,_TRUE_, i, 1);
      class_define_index(peft->index_JG202x,    _TRUE_, i, 1);
      class_define_index(peft->index_JG202y,    _TRUE_, i, 1);
      class_define_index(peft->index_I2211,     _TRUE_, i, 1);
      class_define_index(peft->index_J21111,    _TRUE_, i, 1);
      class_define_index(peft->index_N11x,      _TRUE_, i, 1);
      class_define_index(peft->index_N11y,      _TRUE_, i, 1);
      class_define_index(peft->index_J12102x,   _TRUE_, i, 1);
      class_define_index(peft->index_J12102y,   _TRUE_, i, 1);
      class_define_index(peft->index_I1311,     _TRUE_, i, 1);
      class_define_index(peft->index_J12111p11211,_TRUE_, i, 1);
      //class_define_index(peft->index_J11211,    _TRUE_, i, 1);

      /** - 3-rd moment of RSD expansion */
      class_define_index(peft->index_J21112x,   _TRUE_, i, 1);
      class_define_index(peft->index_J21112y,   _TRUE_, i, 1);
      class_define_index(peft->index_N12x,      _TRUE_, i, 1);
      class_define_index(peft->index_N12y,      _TRUE_, i, 1);
      class_define_index(peft->index_J12112x,   _TRUE_, i, 1);
      class_define_index(peft->index_J12112y,   _TRUE_, i, 1);

      /** - 4-th moment of RSD expansion */
      class_define_index(peft->index_N22x,      _TRUE_, i, 1);
      class_define_index(peft->index_N22y,      _TRUE_, i, 1);
      class_define_index(peft->index_N22z,      _TRUE_, i, 1);
    }

    peft->index_num = i;

    if (!pext || !(pext->loop_matrices_stored) || index >= pext->eft_size || pext->eft_index_num != peft->index_num) {
      /** - allocate loop matrices for computation or loading them from files later */
      /** - store symmetry information */
      class_alloc(peft->symmetry, peft->index_num*sizeof(short), peft->error_message);

      peft->symmetry[peft->index_I2200] = mat_symmetric;
      peft->symmetry[peft->index_I1300] = vec;

      peft->symmetry[peft->index_Idelta200]       = mat_symmetric;
      peft->symmetry[peft->index_IG200]           = mat_symmetric;
      peft->symmetry[peft->index_Idelta2delta200] = mat_symmetric;
      peft->symmetry[peft->index_IG2G200]         = mat_symmetric;
      peft->symmetry[peft->index_Idelta2G200]     = mat_symmetric;
      peft->symmetry[peft->index_FG200]           = vec;

      if (peft->hp->has_rsd) {
        peft->symmetry[peft->index_I2201]     = mat_symmetric;
        peft->symmetry[peft->index_I1301p3101]= mat_none;
        peft->symmetry[peft->index_Idelta201] = mat_symmetric;
        peft->symmetry[peft->index_IG201]     = mat_symmetric;
        peft->symmetry[peft->index_FG201]     = vec;
        peft->symmetry[peft->index_J12101]    = vec;
        peft->symmetry[peft->index_J21101]    = mat_none;
        peft->symmetry[peft->index_Jdelta201] = mat_none;
        peft->symmetry[peft->index_JG201]     = mat_none;

        peft->symmetry[peft->index_J12102x]     = vec;
        peft->symmetry[peft->index_J12102y]     = vec;
        peft->symmetry[peft->index_J21102x]     = mat_symmetric;
        peft->symmetry[peft->index_J21102y]     = mat_symmetric;
        peft->symmetry[peft->index_Jdelta202x]  = mat_symmetric;
        peft->symmetry[peft->index_Jdelta202y]  = mat_symmetric;
        peft->symmetry[peft->index_JG202x]      = mat_symmetric;
        peft->symmetry[peft->index_JG202y]      = mat_symmetric;
        peft->symmetry[peft->index_I2211]       = mat_symmetric;
        peft->symmetry[peft->index_I1311]       = vec;
        peft->symmetry[peft->index_J12111p11211]= vec;
        peft->symmetry[peft->index_J21111]      = mat_none;
        peft->symmetry[peft->index_N11x]        = mat_none;
        peft->symmetry[peft->index_N11y]        = mat_none;
        //peft->symmetry[peft->index_J11211]      = vec;

        peft->symmetry[peft->index_J21112x] = mat_symmetric;
        peft->symmetry[peft->index_J21112y] = mat_symmetric;
        peft->symmetry[peft->index_J12112x] = vec;
        peft->symmetry[peft->index_J12112y] = vec;
        peft->symmetry[peft->index_N12x]    = mat_none;
        peft->symmetry[peft->index_N12y]    = mat_none;

        peft->symmetry[peft->index_N22x] = mat_none;
        peft->symmetry[peft->index_N22y] = mat_none;
        peft->symmetry[peft->index_N22z] = mat_none;
      }
    
      /** - loop kernel matrices */
      class_alloc(peft->loop_matrices,      peft->index_num*sizeof(double complex *), peft->error_message);
      class_alloc(peft->loop_matrices_size, peft->index_num*sizeof(int), peft->error_message);
      for (i = 0; i < peft->index_num; i++) { /** - i = index_moment */
        switch (peft->symmetry[i]) {
        case vec:
          peft->loop_matrices_size[i] = peft->hp->fourier_coeff_size; break;
        case mat_none:
          peft->loop_matrices_size[i] = peft->hp->fourier_coeff_size * peft->hp->fourier_coeff_size; break;
        case mat_symmetric:
          peft->loop_matrices_size[i] = peft->hp->fourier_coeff_size * (peft->hp->fourier_coeff_size + 1) / 2; break;
        default:
          peft->loop_matrices_size[i] = 0; break;
        }

        class_alloc(peft->loop_matrices[i], peft->loop_matrices_size[i]*sizeof(double complex), peft->error_message);
      }
    }

    /** - k-grid and frequencies for Fourier transform */
    for (i = 0; i < eft_tracer_num; i++) {  /** - i = index_tracer */
      class_alloc(peft->ln_k_fourier[i],        peft->hp->k_size_fourier*sizeof(double), peft->error_message);
      class_alloc(peft->fourier_frequencies[i], peft->hp->fourier_coeff_size*sizeof(double), peft->error_message);
    }
  }

  /** TODO: ln_k and pk_l */
  /** - non-linear corrections computed from different pk_types: allocate for P_lin and P_nw if using the approximated mu-dependence
   *    spectra with exact mu-dependence are computed on-demand */
  if (peft->hp->compute_mu_approximation) {
    for (j = 0; j <= pk_nowiggle; j++) {
      class_alloc(peft->spectra_contributions[j], peft->index_num*sizeof(double *), peft->error_message);
      for (i = 0; i < peft->index_num; i++) {
        class_alloc(peft->spectra_contributions[j][i], peft->k_size, peft->error_message);
      }
    }

    /** - allocate arrays for the Fourier series inputs, coefficients and frequencies */
    for (j = 0; j <= pk_nowiggle; j++) {  /** - j = index_pk_type */
      for (i = 0; i < eft_tracer_num; i++) {  /** - i = index_tracer */
        class_alloc(peft->pk_l_biased[j*eft_tracer_num + i],           peft->hp->k_size_fourier*sizeof(double), peft->error_message);
        class_alloc(peft->ddpk_l_biased[j*eft_tracer_num + i],         peft->hp->k_size_fourier*sizeof(double), peft->error_message);
        class_alloc(peft->fourier_coeff[j*eft_tracer_num + i],         peft->hp->fourier_coeff_size*sizeof(double complex), peft->error_message);
        class_calloc(peft->fourier_condition_num[j*eft_tracer_num + i], peft->hp->fourier_coeff_size, sizeof(double complex), peft->error_message);
      }
    }
  }
  else {
    /** TODO: request specific mu's (maybe unnecessary) */
  }

  // if (peft->hp->fourier_mode == fourier_mode_fft) {
  //   FFT_planner_init(peft->hp->fourier_coeff_size - 1, &(peft->fft_plan));
  // }

  /** - allocate arrays for different P_lin types, pk_index_num is the last entry in enum eft_pk_type */
  // class_alloc(peft->pk_l, pk_index_num*sizeof(double *), peft->error_message);
  // for (i = 0; i < pk_index_num; i++) {  /** - i = index_pk_type */
  //   class_alloc(peft->pk_l[i], peft->k_size*sizeof(double), peft->error_message); /** TODO: all different spectra types needed? */
  // }

  

  return _SUCCESS_;
}

int eft_free(struct eft * peft) {
  int i, j;

  if (peft->role == eft_master) {
    free(peft->symmetry);

    for (i = 0; i < peft->index_num; i++) { /** - i = index_moment */
      free(peft->loop_matrices[i]);
    }
    free(peft->loop_matrices);
    free(peft->loop_matrices_size);

    for (i = 0; i < eft_tracer_num; i++) {  /** - i = index_tracer */
      free(peft->ln_k_fourier[i]);
      free(peft->fourier_frequencies[i]);
    }

    if (peft->hp->fourier_mode == fourier_mode_fft) {
      FFT_planner_free(&(peft->fft_plan));
    }
  }

  if (peft->hp->compute_mu_approximation) {
    for (j = 0; j <= pk_nowiggle; j++) {
      for (i = 0; i < peft->index_num; i++) {
        free(peft->spectra_contributions[j][i]);
      }
      free(peft->spectra_contributions[j]);
    }

    // for (i = 0; i < pk_index_num; i++) {  /** - i = index_pk_type */
    //   free(peft->pk_l[i]);
    // }
    // free(peft->pk_l);

    for (j = 0; j <= pk_nowiggle; j++) {  /** - j = index_pk_type */
      for (i = 0; i < eft_tracer_num; i++) {  /** - i = index_tracer */
        free(peft->pk_l_biased[j*eft_tracer_num + i]);
        free(peft->ddpk_l_biased[j*eft_tracer_num + i]);
        free(peft->fourier_coeff[j*eft_tracer_num + i]);
        free(peft->fourier_condition_num[j*eft_tracer_num + i]);
      }
    }
  }

  return _SUCCESS_;
}

int eft_init(struct precision * ppr,
            struct background * pba,
            struct eft * peft,
            struct eft_hyper_parameters * eft_hp,
            struct eft_input_parameters * eft_ip,
            struct ext_storage * pext,
            const short eft_role,
            const int index) {

  /** - register pointers to the EFT parameter structs */
  peft->hp = eft_hp;
  peft->ip = eft_ip;
  peft->role = eft_role;

  class_call(eft_indices(peft, pext, eft_role, index), peft->error_message, peft->error_message);

  if (eft_role == eft_master) {
    for (int index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
      /** Generate the sampling points in k; the BAO feature position is given by a precision parameter */
      eft_spline_sample_points_nonuniform(peft->hp->kmin_lin[index_tracer], exp(ppr->k_bao_center), peft->hp->kmax_lin[index_tracer],    \
                                          peft->hp->bao_oversampling, peft->hp->ln_k_oversampling_width,                \
                                          peft->hp->k_size_fourier-1, peft->ln_k_fourier[index_tracer], peft->hp->eft_verbose);
      /** - last point in the array is exactly one period removed from k_min */
      peft->ln_k_fourier[index_tracer][peft->hp->k_size_fourier-1] = log(peft->hp->kmin_lin[index_tracer]) + peft->hp->period[index_tracer];
      /** Generate the Fourier frequencies */
      eft_get_fourier_frequencies(peft->hp->period[index_tracer], peft->hp->fourier_coeff_size, peft->fourier_frequencies[index_tracer]);
    }

    if (peft->hp->fourier_mode == fourier_mode_fft) {
      FFT_planner_init(peft->hp->fourier_coeff_size-1, &(peft->fft_plan));
    }

  }

  return _SUCCESS_;
}

/**
 * @brief Get the frequencies at which the spectrum is sampled by FFT in standard order. 
 *        Returns f = [0, 1, ..., ceil(N/2)-1, -floor(N/2), -floor(N/2)+1, ..., -1] / T.
 *        Definition is compatible with numpy.fft
 *         
 * @param period          Input: Period of the input signal
 * @param num_frequencies Input: number of Fourier components
 * @param frequencies     Input: pointer to frequency array of size num_frequencies
 * 
 * @return the error status (always succeeds)
 */
int eft_get_fourier_frequencies(
                  const double period,
                  const int num_frequencies,
                  double * frequencies) {

  int m;
  for (m = 0; m < (int)ceil(num_frequencies/2.); m++) {
    frequencies[m] = _TWOPI_ * m/period;
  }
  for (m = -(int)floor(num_frequencies/2.); m < 0; m++) {
    frequencies[m + num_frequencies] = _TWOPI_ * m/period;
  }
  return _SUCCESS_;
}


int eft_load_linear_spectra(struct precision * ppr,
                            struct background * pba,
                            struct fourier * pfo,
                            struct primordial * ppm,
                            struct eft * peft,
                            double z,
                            double f,
                            double D,
                            const short use_mu_approximation) {

  int it, index_mu, index_tracer, index_pk_type, abort = _FALSE_;
  int mu_size, index_pk_start, index_pk_stop;
  const int num_independent_coefficients = peft->hp->num_positive_fourier_freq + 1;
  double ** pk_l_biased_p;
  double spline_2nd_der_L1_norm = 0., spline_L2_norm, series_l2_norm = 0.;

  peft->z0 = z;
  peft->f_z0 = f; /**< growth rate may be supplied externally to evaluate at some fiducial */
  peft->D_z0 = D;
  
  /** Compute the suppression factor for IR-Resummation */
  peft->Sigma2_ir  =  eft_ir_sigma2(pba, ppm, pfo, peft->z0, pfo->wnw_k_split*pba->h, pfo->wnw_k_feature*pba->h);  /** k_split = 0.2 h/Mpc, k_feature = 1/110 h/Mpc */
  peft->dSigma2_ir = eft_ir_dsigma2(pba, ppm, pfo, peft->z0, pfo->wnw_k_split*pba->h, pfo->wnw_k_feature*pba->h);  /** according to arXiV:1804.05080 by Ivanov & Sibiryakov */

  #pragma omp parallel shared(ppr, pba, ppm, pfo, peft, abort, stderr, use_mu_approximation, num_independent_coefficients, index_pk_start, index_pk_stop, spline_2nd_der_L1_norm, spline_L2_norm, series_l2_norm), \
                       private(it, index_mu, index_tracer, index_pk_type, mu_size, pk_l_biased_p), \
                       default(none)
  { /** , last_index, ln_k_fft, pk_l_biased_fft, fft_plan, fourier_coeff_real, fourier_coeff_imag */

  /** Retrieve the power spectra */
  if (use_mu_approximation) {
    /** ---------------------- Load pk_lin and pk_nowiggle -------------------- */
    #pragma omp for schedule(static, 1)
    for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
      class_call_parallel(fourier_pk_at_kvec_and_z(pba, ppm, pfo, linear, pk_linear,
                                                  peft->ln_k_fourier[index_tracer],
                                                  peft->hp->k_size_fourier-1,
                                                  peft->z0,
                                                  pfo->index_pk_cluster,
                                                  peft->pk_l_biased[pk_lin*eft_tracer_num + index_tracer]),
                          pfo->error_message,
                          peft->error_message);
      
      class_call_parallel(fourier_pk_nw_at_kvec_and_z(pba, ppm, pfo, linear, 
                                                      peft->ln_k_fourier[index_tracer],
                                                      peft->hp->k_size_fourier-1,
                                                      peft->z0,
                                                      peft->pk_l_biased[pk_nowiggle*eft_tracer_num + index_tracer]),
                          pfo->error_message,
                          peft->error_message);
    }

    #pragma omp master
    {
    /** - used range of pk_type */
    index_pk_start = pk_lin;
    index_pk_stop = pk_nowiggle;
    }
    #pragma omp barrier

    if (!abort) {
      #pragma omp for schedule(static, 1) collapse(2)
      for (index_pk_type = index_pk_start; index_pk_type <= index_pk_stop; index_pk_type++) {
        for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
          for (it = 0; it < peft->hp->k_size_fourier-1; it++) {
            /** Apply the bias to the power spectra */
            peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][it]      \
                  *= exp(-peft->hp->bias[index_tracer] * peft->ln_k_fourier[index_tracer][it]);
          }
          /** Add the last point to make the spectra periodic: P_bias(ln(k_min) + period) = P_bias(ln(k_min)) */
          peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][peft->hp->k_size_fourier-1]    \
                = peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][0];
        
          /** Spline the biased power spectra periodically */
          class_call_parallel(array_spline_table_columns(peft->ln_k_fourier[index_tracer],
                                                        peft->hp->k_size_fourier,
                                                        peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer],
                                                        1,
                                                        peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer],
                                                        _SPLINE_PERIODIC_,
                                                        peft->error_message),
                              peft->error_message,
                              peft->error_message);
          
        }
      }
    }

  }
  else {
    /** ---------------- Load pk_ir_resummed_lo for different mu -------------- */
    /** mu_size, mu are set before & dd/pk_l_biased[pk_ir_resummed_lo*eft_tracer_num + index_tracer] needs to be allocated with size peft->mu_size*peft->hp->k_size_fourier */
    /** fourier_coeff/condition_num[pk_ir_resummed_lo*eft_tracer_num + index_tracer] needs to be allocated with size peft->mu_size*peft->hp->fourier_coeff_size */
    class_alloc_parallel(pk_l_biased_p, peft->mu_size*sizeof(double **), peft->error_message);
    
    #pragma omp for schedule(static, 1)
    for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
      for (index_mu = 0; index_mu < peft->mu_size; index_mu++) {
        pk_l_biased_p[index_mu] = peft->pk_l_biased[pk_ir_resummed_lo*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier;
      }
      class_call_parallel(eft_ir_pk_rsd_lo(pba, ppm, pfo, linear, 
                                          peft->ln_k_fourier[index_tracer],
                                          peft->hp->k_size_fourier-1,
                                          peft->z0,
                                          peft->f_z0,
                                          peft->mu,
                                          peft->mu_size,
                                          peft->Sigma2_ir,
                                          peft->dSigma2_ir,
                                          pk_l_biased_p),
                          peft->error_message,
                          peft->error_message);

      // for (index_mu = 0; index_mu < peft->mu_size; index_mu++) {
      //   pk_l_biased_p[index_mu] = peft->pk_l_biased[pk_ir_resummed_nlo*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier;
      // }
      // class_call(eft_ir_pk_rsd_nlo(pba, ppm, pfo, linear, 
      //                             peft->ln_k_fourier[index_tracer],
      //                             peft->hp->k_size_fourier-1,
      //                             peft->z0,
      //                             peft->f_z0,
      //                             peft->mu,
      //                             peft->mu_size,
      //                             peft->Sigma2_ir,
      //                             peft->dSigma2_ir,
      //                             pk_l_biased_p),
      //             peft->error_message,
      //             peft->error_message);
    }
    free(pk_l_biased_p);

    #pragma omp master
    {
    /** - used range of pk_type */
    index_pk_start = pk_ir_resummed_lo;
    index_pk_stop = pk_ir_resummed_lo;
    }
    #pragma omp barrier

    if (!abort) {
      #pragma omp for schedule(static, 1) collapse(3)
      for (index_pk_type = index_pk_start; index_pk_type <= index_pk_stop; index_pk_type++) {
        for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
          for (index_mu = 0; index_mu < peft->mu_size; index_mu++) {
            for (it = 0; it < peft->hp->k_size_fourier-1; it++) {
              /** Apply the bias to the power spectra */
              peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->k_size_fourier + it]      \
                    *= exp(-peft->hp->bias[index_tracer] * peft->ln_k_fourier[index_tracer][it]);
            }
            /** Add the last point to make the spectra periodic: P_bias(ln(k_min) + period) = P_bias(ln(k_min)) */
            peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->k_size_fourier + peft->hp->k_size_fourier-1]    \
                  = peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][0];
          
            /** Spline the biased power spectra periodically */
            class_call_parallel(array_spline_table_columns(peft->ln_k_fourier[index_tracer],
                                                          peft->hp->k_size_fourier,
                                                          peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
                                                          1,
                                                          peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
                                                          _SPLINE_PERIODIC_,
                                                          peft->error_message),
                                peft->error_message,
                                peft->error_message);
          
          }
        }
      }
    }

  }

  #pragma omp barrier

  /* Compute the Fourier coefficients of the power spectra */
  switch (peft->hp->fourier_mode)
  {
  case fourier_mode_spline:
    for (index_pk_type = index_pk_start; index_pk_type <= index_pk_stop; index_pk_type++) {
      mu_size = (index_pk_type <= pk_nowiggle) ? 1 : peft->mu_size;
      for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
        for (index_mu = 0; index_mu < mu_size; index_mu++) {

          #pragma omp master
          {
          peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + 0] = CMPLX(0., 0.);
          /** The DC component has to be treated differently */
          class_call_parallel(array_integrate_all_spline_table_lines_compensated(
                                      peft->ln_k_fourier[index_tracer],
                                      peft->hp->k_size_fourier,
                                      peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
                                      1,
                                      peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
                                      (double *)(peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->fourier_coeff_size),
                                      ppr->eft_fourier_condition_num_threshold,
                                      (double *)(peft->fourier_condition_num[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->fourier_coeff_size),
                                      peft->error_message),
                              peft->error_message,
                              peft->error_message);
          peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + 0] /= peft->hp->period[index_tracer];
          }

          #pragma omp for schedule(static)
          for (it = 1; it < num_independent_coefficients; it++) {
            array_integrate_all_spline_table_lines_fourier_compensated(
                              peft->ln_k_fourier[index_tracer],
                              peft->hp->k_size_fourier,
                              peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
                              1,
                              peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
                              peft->fourier_frequencies[index_tracer][it],
                              peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->fourier_coeff_size + it,
                              ppr->eft_fourier_condition_num_threshold,
                              peft->fourier_condition_num[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->fourier_coeff_size + it,
                              peft->error_message);
            
            peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + it] /= peft->hp->period[index_tracer];
          }

          #pragma omp for schedule(static)
          /** - negative frequency components in standard order are fully dependent on the others */
          for (it = -(num_independent_coefficients-1); it < 0; it++) {
            peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + peft->hp->fourier_coeff_size + it]     \
                = conj( peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size - it] );
          }

          /** - debug information */
          if (peft->hp->eft_verbose > 2) {
            #pragma omp master
            {
            spline_2nd_der_L1_norm = 0.;
            series_l2_norm = 0.;
            }
            double h;
            double * M0 = peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier;
            #pragma omp barrier
            /** - compute the L1 norm of the second derivative: since the (biased) spline is a C^2 function,
             *    the magnitude of its Fourier coefficients is bounded by || S^(2) ||_L1 / frequency^2 */
            #pragma omp for schedule(static) reduction(+:spline_2nd_der_L1_norm)
            for (it = 1; it < peft->hp->k_size_fourier; it++) {
              h = peft->ln_k_fourier[index_tracer][it] - peft->ln_k_fourier[index_tracer][it-1];
              if (M0[it-1] * M0[it] >= 0.) {
                spline_2nd_der_L1_norm += 0.5*h * fabs(M0[it-1] + M0[it]);
              }
              else {
                spline_2nd_der_L1_norm += 0.5*h * (M0[it-1]*M0[it-1] + M0[it]*M0[it]) / fabs(M0[it-1] - M0[it]);
              }
            }
            /** - compute the l2 norm of its Fourier coefficients: Because of Parseval's theorem, 
             *    these norms must be equal when taking into account infinite Fourier components.
             *    By truncating the series we generate a loss whose norm we can evaluate. */
            #pragma omp for schedule(static) reduction(+:series_l2_norm)
            for (it = 0; it < peft->hp->fourier_coeff_size; it++) {
              series_l2_norm += pow(cabs(peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + it]), 2);
            }
            /** - compute the L2 norm of the biased power spectra */
            #pragma omp master
            {
            class_call_parallel(array_square_integrate_all_spline_table_lines(
                                                                      peft->ln_k_fourier[index_tracer],
                                                                      peft->hp->k_size_fourier,
                                                                      peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
                                                                      1,
                                                                      peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
                                                                      &spline_L2_norm,
                                                                      peft->error_message),
                                peft->error_message,
                                peft->error_message);
            spline_L2_norm /= peft->hp->period[index_tracer];
            double rel_loss_norm = (spline_L2_norm - series_l2_norm) / spline_L2_norm;
            printf("index_pk_type = %d, index_tracer = %d, index_mu = %d: Performed Spline Fourier with bias %.2f and relative power loss of %.3e at frequencies higher than %.3e \n", index_pk_type, index_tracer, index_mu, peft->hp->bias[index_tracer], rel_loss_norm, _TWOPI_*(num_independent_coefficients-1)/peft->hp->period[index_tracer]);
            printf("                                                      Fourier coefficients are bounded by %.3e / frequency^2 \n", spline_2nd_der_L1_norm);
            }
          }

        }
      }
    }
    
    break;
  
  default:
    #pragma omp master
    {
    class_protect_fprintf(stderr, "Fourier mode setting is invalid. Defaulting to FFT! \n");
    }
    /** -------------- fall-through -----------------------*/
  case fourier_mode_fft:
    /** - sample fourier_coeff_size-1 equidistant points for FFT */
    int last_index = 0;
    double ln_k_fft;
    double * pk_l_biased_fft[eft_tracer_num];
    struct FFT_plan * fft_plan;
    double * fourier_coeff_real[2], * fourier_coeff_imag[2];
    
    class_alloc_parallel(pk_l_biased_fft[eft_matter], (peft->hp->fourier_coeff_size-1)*sizeof(double), peft->error_message);
    class_alloc_parallel(pk_l_biased_fft[eft_halo], (peft->hp->fourier_coeff_size-1)*sizeof(double), peft->error_message);
    class_alloc_parallel(fourier_coeff_real[0], (peft->hp->fourier_coeff_size-1)*sizeof(double), peft->error_message);
    class_alloc_parallel(fourier_coeff_imag[0], (peft->hp->fourier_coeff_size-1)*sizeof(double), peft->error_message);
    class_alloc_parallel(fourier_coeff_real[1], (peft->hp->fourier_coeff_size-1)*sizeof(double), peft->error_message);
    class_alloc_parallel(fourier_coeff_imag[1], (peft->hp->fourier_coeff_size-1)*sizeof(double), peft->error_message);

    /** - copy the master FFT plan */
    FFT_planner_alloc(peft->fft_plan->N, &fft_plan);
    class_protect_memcpy(fft_plan->cos_vals, peft->fft_plan->cos_vals, peft->fft_plan->N * sizeof(double));
    class_protect_memcpy(fft_plan->sin_vals, peft->fft_plan->sin_vals, peft->fft_plan->N * sizeof(double));

    mu_size = use_mu_approximation ? 1 : peft->mu_size;

    #pragma omp for schedule(static, 1) collapse(2)
    for (index_pk_type = index_pk_start; index_pk_type <= index_pk_stop; index_pk_type++) {
      for (index_mu = 0; index_mu < mu_size; index_mu++) {
        for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
          for (it = 0; it < peft->hp->fourier_coeff_size-1; it++) {
            /** - equidistant sample points; last point ln(k_min) + period is left out since it is known by periodicity  */
            ln_k_fft = log(peft->hp->kmin_lin[index_tracer]) + peft->hp->period[index_tracer] * it/(double)(peft->hp->fourier_coeff_size-1);
            class_call_parallel(array_interpolate_spline_growing_closeby(peft->ln_k_fourier[index_tracer],
                                                                        peft->hp->k_size_fourier,
                                                                        peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
                                                                        peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
                                                                        1,
                                                                        ln_k_fft,
                                                                        &last_index,
                                                                        pk_l_biased_fft[index_tracer] + it,
                                                                        1,
                                                                        peft->error_message),
                                peft->error_message,
                                peft->error_message);
          }
        }

        if (!abort) {
          /** - compute the FFT of the spectra for both tracers simultaneously */
          if (eft_tracer_num == 2) {
            FFT_real_planned(pk_l_biased_fft[eft_matter], pk_l_biased_fft[eft_halo], 
                            fourier_coeff_real[0], fourier_coeff_imag[0], 
                            fourier_coeff_real[1], fourier_coeff_imag[1], 
                            fft_plan);

            for (index_tracer = 0; index_tracer < 2; index_tracer++) {
              /** - DC coefficient is halved to keep the original power while symmetrizing the spectrum btw. -N/2...N/2 */
              peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + 0]      \
                    = 0.5 * CMPLX(fourier_coeff_real[index_tracer][0], fourier_coeff_imag[index_tracer][0]) / (double)(peft->hp->fourier_coeff_size-1); 
              /** - positive frequency components in standard order, last iteration extends the spectrum symmetrically */
              for (it = 1; it < num_independent_coefficients; it++) {
                peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + it]      \
                    = CMPLX(fourier_coeff_real[index_tracer][it], fourier_coeff_imag[index_tracer][it]) * cpow(peft->hp->kmin_lin[index_tracer], -_Complex_I*peft->fourier_frequencies[index_tracer][it]) / (double)(peft->hp->fourier_coeff_size-1);
              }
              /** - negative frequency components in standard order are fully dependent on the others */
              for (it = -(num_independent_coefficients-1); it < 0; it++) {
                peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + peft->hp->fourier_coeff_size + it]     \
                    = conj( peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size - it] );
              }
            }
          }
          else {
            abort = _TRUE_;
            ErrorMsg errmsg;
            class_protect_sprintf(errmsg, "FFT implementation requires 2 tracers for performance reasons, but the eft_tracer contains %d entries. If you need more tracers, implement the FFT schedule here.", eft_tracer_num);
            class_build_error_string(peft->error_message, "error; %s", errmsg);
          }
        }
      }
    }


    free(pk_l_biased_fft[eft_matter]); free(pk_l_biased_fft[eft_halo]);
    free(fourier_coeff_real[0]); free(fourier_coeff_imag[0]);
    free(fourier_coeff_real[1]); free(fourier_coeff_imag[1]);
    FFT_planner_free(&fft_plan);
    break;
  }

  } /** - end of parallel region */

  if (use_mu_approximation && (pfo->fourier_verbose > 2)) {
    FILE *ffourier = fopen("output/halo_biased_pk_lin_samples.dat", "w");

    fprintf(ffourier, "# Biased samples of the linear power spectrum at z=%.3f \n", peft->z0);
    fprintf(ffourier, "# for k=%.4e to %.4e \n", exp(peft->ln_k_fourier[eft_halo][0]), exp(peft->ln_k_fourier[eft_halo][peft->hp->k_size_fourier-1]));
    fprintf(ffourier, "# number of wavenumbers equal to %d \n", peft->hp->k_size_fourier);
    fprintf(ffourier, "#    1:k (1/Mpc)            2:P_bias (Mpc)^3        3:d^2P_bias/dln(k)^2 (Mpc)^3 \n");
    for (int i = 0; i < peft->hp->k_size_fourier; i++)
      fprintf(ffourier, "  %.16e       %.16e       %+.16e \n", \
              exp(peft->ln_k_fourier[eft_halo][i]), peft->pk_l_biased[pk_lin*eft_tracer_num + eft_halo][i], peft->ddpk_l_biased[pk_lin*eft_tracer_num + eft_halo][i]);

    fclose(ffourier);

    ffourier = fopen("output/halo_pk_lin_fourier_coefficients.dat", "w");

    fprintf(ffourier, "# Fourier coefficients of the linear power spectrum at z=%.3f \n", peft->z0);
    fprintf(ffourier, "# for omega=%.3f to %.3f \n", peft->fourier_frequencies[eft_halo][num_independent_coefficients], peft->fourier_frequencies[eft_halo][num_independent_coefficients-1]);
    fprintf(ffourier, "# number of frequencies equal to %d \n", peft->hp->fourier_coeff_size);
    fprintf(ffourier, "#    1:omega                  2:Re c (Mpc)^3                  3:Im c (Mpc)^3                  4:Cond[Re c]                  5:Cond[Im c] \n");
    for (int i = 0; i < peft->hp->fourier_coeff_size; i++)
      fprintf(ffourier, "  %.16e       %+.16e       %+.16e       %+.16e       %+.16e \n", \
              peft->fourier_frequencies[eft_halo][i], creal(peft->fourier_coeff[pk_lin*eft_tracer_num + eft_halo][i]), cimag(peft->fourier_coeff[pk_lin*eft_tracer_num + eft_halo][i]),     \
              creal(peft->fourier_condition_num[pk_lin*eft_tracer_num + eft_halo][i]), cimag(peft->fourier_condition_num[pk_lin*eft_tracer_num + eft_halo][i]));

    fclose(ffourier);
  }


  return abort;
}

int eft_save_matrix_to_file(const double complex * matrix,
                            const int size,
                            const short symmetry,
                            const short tracer,
                            const double period,
                            FileName filename,
                            const short permit_overwrite,
                            ErrorMsg errmsg,
                            const short print_message) {

  FILE * file;
  size_t write_count = 0, read_count = 0;
  const long info = (size << 32) | (symmetry << 16) | tracer;
  long read_info;
  double read_period;
  short write;

  /** - check if the file exists  */
  file = fopen(filename, "rb");
  if (file) {
    /** - file already exists: compare metadata (tests for changes in size, symmetry, tracer species or period) */
    read_count += fread(&read_info, sizeof(long), 1, file);
    read_count += fread(&read_period, sizeof(double), 1, file);
    if (read_count == 2) {
      /** - file is intact: overwrite if metadata changed and it is permitted */
      write = ((read_info == info) && (fabs(period - read_period)/period < _EPSILON_)) ? _FALSE_ : permit_overwrite;
    }
    else {
      /** - file is invalid: overwrite if permitted, otherwise fail */
      write = permit_overwrite;
      if (!write) { class_stop(errmsg, "Kernel matrix file %s exists, but is invalid and overwriting is forbidden.", filename); }
    }

    if (write) { freopen(NULL, "wb+", file); }
    else { fclose(file); }
  }
  else {
    /** - file does not exist: write is always permitted */
    write = _TRUE_;
    file = fopen(filename, "wb+");
  }

  /** - file is already open if write is set */
  if (write) {
    /** Loop matrix file should contain 
     * i. (long) info = (size | symmetry | tracer)
     * ii.(double) period in ln(k)-space used to generate the matrix
     * iii.(double complex) matrix entries according to LAPACK storage schemes
    */
    write_count += fwrite(&info, sizeof(long), 1, file);
    write_count += fwrite(&period, sizeof(double), 1, file);
    write_count += fwrite(matrix, sizeof(double complex), size, file);
    fclose(file);
    
    if (ferror(file) || write_count < size+2) {
      /** - I/O error */
      class_stop(errmsg, "Writing kernel matrix in %s failed after %d elements.", filename, write_count);
    }
    else if (print_message) {
      printf("Saved kernel matrix from %s with size %s, symmetry type %d, tracer species %d and period %.3e \n", filename, size, symmetry, tracer, period);
    }
  }

  return _SUCCESS_;
}

int eft_read_matrix_from_file(double complex * matrix,
                            const int size,
                            const short symmetry,
                            const short tracer,
                            const double period,
                            FileName filename,
                            ErrorMsg errmsg,
                            short print_message) {

  FILE * file;
  size_t read_count = 0;
  long read_info;
  int read_size;
  short read_symmetry, read_tracer;
  double read_period;

  class_open(file, filename, "rb", errmsg);

  /** Loop matrix file should contain 
   * i. (long) info = (size | symmetry | tracer)
   * ii.(double) period in ln(k)-space used to generate the matrix
   * iii.(double complex) matrix entries according to LAPACK storage schemes
  */
  read_count += fread(&read_info, sizeof(long), 1, file);
  read_count += fread(&read_period, sizeof(double), 1, file);
  read_count += fread(matrix, sizeof(double complex), size, file);
  fclose(file);
  
  if (ferror(file) || read_count < size+2) {
    /** - I/O error */
    class_stop(errmsg, "Reading kernel matrix in %s failed after %d elements.", filename, write_count);
  }
  else {
    /** - recover metadata and compare to expected values */
    read_size = read_info >> 32; read_symmetry = read_info >> 16; read_tracer = read_info;
    if (read_size != size || read_symmetry != symmetry || read_tracer != tracer || (fabs(period - read_period)/period > _EPSILON_)) {
      class_stop(errmsg, "Stored kernel matrix in %s has size %d, symmetry type %d, tracer species %d and period %.3e; expected size %d, symmetry type %d, tracer species %d and period %.3e", \
                  filename, read_size, read_symmetry, read_tracer, read_period, size, symmetry, tracer, period);
    }

    if (print_message) {
      printf("Read kernel matrix from %s with size %s, symmetry type %d, tracer species %d and period %.3e \n", filename, read_size, read_symmetry, read_tracer, read_period);
    }
  }

  return _SUCCESS_;
}

enum eft_loop_decision {decision_compute, decision_load_file, decision_load_ext};

int eft_get_loop_matrices(struct eft * peft,
                          struct ext_storage * pext,
                          const int index) {

  int decision, done;  /** - loop matrices for the current eft struct may need to be recomputed independently of the hyperparameters, if stored matrices are invalid */

  switch (peft->role)
  {
  case eft_master:
    /** master structure is responsible for writing the loop matrices to the arrays, while the slaves only read them */
    /** decision tree */
    if (!compute_loop_matrices) {
      if (pext && pext->loop_matrices_stored) {
        /** - external storage holds some loop matrices, check if they are compatible to the current configuration */
        if (index < pext->eft_size && pext->eft_index_num >= peft->index_num) {
          /** - insert loop matrices from external storage */
          class_call(ext_insert_eft(pext, peft, index, peft->error_message),
                    peft->error_message,
                    peft->error_message);
        }
        else {
          /** - recompute the matrices if possible, arrays must have been allocated */
          if (peft->hp->use_time_independent_kernels) {
            compute_loop_matrices = _TRUE_; /** - compute matrices for the current eft struct */
          }



        }
      }
      
    }
    break;

  case eft_slave:
    /** no-op */
    break;
  default:
    break;
  }

  return _SUCCESS_;
}

