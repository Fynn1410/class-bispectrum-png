/** @file eft_main.c
 *
 * author: Christian Radermacher, 2024
 * based on prototype version of Azadeh Moradinezhad Dizgah & Dennis Linde
 *
 * contains many routines for the EFT workflow and handles all input/output operations
*/

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

int eft_allocate_loop_matrices(struct eft * peft) {
  int i;

  /** - allocate loop matrices for computation or loading them from files later */
  /** - store symmetry information */
  class_alloc(peft->symmetry,   peft->index_num*sizeof(short), peft->error_message);
  class_alloc(peft->use_tracer, peft->index_num*sizeof(short), peft->error_message);
  class_alloc(peft->spectra_contributions_dimension, peft->index_num*sizeof(short), peft->error_message);

  for (i = 0; i < peft->index_num; i++) { /** i = index_moment*/
    peft->spectra_contributions_dimension[i] = 3; /** default wavenumber dimension is 3 for any integral without the prefactors (from integration measure, integrands are mostly dimensionless) */
  }

  peft->symmetry[peft->index_I2200] = sym_mat_symmetric; peft->use_tracer[peft->index_I2200] = eft_matter;
  peft->symmetry[peft->index_I1300] = sym_vec;           peft->use_tracer[peft->index_I1300] = eft_matter;

  peft->symmetry[peft->index_Idelta200]       = sym_mat_symmetric; peft->use_tracer[peft->index_Idelta200]       = eft_halo;
  peft->symmetry[peft->index_IG200]           = sym_mat_symmetric; peft->use_tracer[peft->index_IG200]           = eft_halo;
  peft->symmetry[peft->index_Idelta2delta200] = sym_mat_symmetric; peft->use_tracer[peft->index_Idelta2delta200] = eft_halo;
  peft->symmetry[peft->index_IG2G200]         = sym_mat_symmetric; peft->use_tracer[peft->index_IG2G200]         = eft_halo;
  peft->symmetry[peft->index_Idelta2G200]     = sym_mat_symmetric; peft->use_tracer[peft->index_Idelta2G200]     = eft_halo;
  peft->symmetry[peft->index_FG200]           = sym_vec;           peft->use_tracer[peft->index_FG200]           = eft_halo;

  if (peft->hp->has_rsd) {
    peft->symmetry[peft->index_I2201]     = sym_mat_symmetric; peft->use_tracer[peft->index_I2201]     = eft_matter;
    peft->symmetry[peft->index_I1301p3101]= sym_vec;           peft->use_tracer[peft->index_I1301p3101]= eft_matter;
    peft->symmetry[peft->index_Idelta201] = sym_mat_symmetric; peft->use_tracer[peft->index_Idelta201] = eft_halo;
    peft->symmetry[peft->index_IG201]     = sym_mat_symmetric; peft->use_tracer[peft->index_IG201]     = eft_halo;
    peft->symmetry[peft->index_FG201]     = sym_vec;           peft->use_tracer[peft->index_FG201]     = eft_halo;
    peft->symmetry[peft->index_J12101]    = sym_vec;           peft->use_tracer[peft->index_J12101]    = eft_halo;   peft->spectra_contributions_dimension[peft->index_J12101]    = 2;
    peft->symmetry[peft->index_J21101]    = sym_mat_none;      peft->use_tracer[peft->index_J21101]    = eft_matter; peft->spectra_contributions_dimension[peft->index_J21101]    = 2;
    peft->symmetry[peft->index_Jdelta201] = sym_mat_none;      peft->use_tracer[peft->index_Jdelta201] = eft_halo;   peft->spectra_contributions_dimension[peft->index_Jdelta201] = 2;
    peft->symmetry[peft->index_JG201]     = sym_mat_none;      peft->use_tracer[peft->index_JG201]     = eft_halo;   peft->spectra_contributions_dimension[peft->index_JG201]     = 2;
    peft->symmetry[peft->index_J11201]    = no_finite_part;    peft->use_tracer[peft->index_J11201]    = eft_halo;   peft->spectra_contributions_dimension[peft->index_J11201]    = 2;

    peft->symmetry[peft->index_J12102x]     = sym_vec;           peft->use_tracer[peft->index_J12102x]     = eft_halo;   peft->spectra_contributions_dimension[peft->index_J12102x]    = 1;
    peft->symmetry[peft->index_J12102y]     = sym_vec;           peft->use_tracer[peft->index_J12102y]     = eft_halo;   peft->spectra_contributions_dimension[peft->index_J12102y]    = 1;
    peft->symmetry[peft->index_J21102x]     = sym_mat_symmetric; peft->use_tracer[peft->index_J21102x]     = eft_matter; peft->spectra_contributions_dimension[peft->index_J21102x]    = 1;
    peft->symmetry[peft->index_J21102y]     = sym_mat_symmetric; peft->use_tracer[peft->index_J21102y]     = eft_matter; peft->spectra_contributions_dimension[peft->index_J21102y]    = 1;
    peft->symmetry[peft->index_Jdelta202x]  = sym_mat_symmetric; peft->use_tracer[peft->index_Jdelta202x]  = eft_halo;   peft->spectra_contributions_dimension[peft->index_Jdelta202x] = 1;
    peft->symmetry[peft->index_Jdelta202y]  = sym_mat_symmetric; peft->use_tracer[peft->index_Jdelta202y]  = eft_halo;   peft->spectra_contributions_dimension[peft->index_Jdelta202y] = 1;
    peft->symmetry[peft->index_JG202x]      = sym_mat_symmetric; peft->use_tracer[peft->index_JG202x]      = eft_halo;   peft->spectra_contributions_dimension[peft->index_JG202x]     = 1;
    peft->symmetry[peft->index_JG202y]      = sym_mat_symmetric; peft->use_tracer[peft->index_JG202y]      = eft_halo;   peft->spectra_contributions_dimension[peft->index_JG202y]     = 1;
    peft->symmetry[peft->index_I2211]       = sym_mat_symmetric; peft->use_tracer[peft->index_I2211]       = eft_matter;
    peft->symmetry[peft->index_I1311]       = sym_vec;           peft->use_tracer[peft->index_I1311]       = eft_matter;
    peft->symmetry[peft->index_J12111]      = sym_vec;           peft->use_tracer[peft->index_J12111]      = eft_halo;   peft->spectra_contributions_dimension[peft->index_J12111]     = 2;
    peft->symmetry[peft->index_J21111]      = sym_mat_none;      peft->use_tracer[peft->index_J21111]      = eft_matter; peft->spectra_contributions_dimension[peft->index_J21111]     = 2;
    peft->symmetry[peft->index_N11x]        = sym_mat_none;      peft->use_tracer[peft->index_N11x]        = eft_matter; peft->spectra_contributions_dimension[peft->index_N11x]       = 1;
    peft->symmetry[peft->index_N11y]        = sym_mat_none;      peft->use_tracer[peft->index_N11y]        = eft_matter; peft->spectra_contributions_dimension[peft->index_N11y]       = 1;
    peft->symmetry[peft->index_J11211]      = no_finite_part;    peft->use_tracer[peft->index_J11211]      = eft_halo;   peft->spectra_contributions_dimension[peft->index_J11211]     = 2;

    peft->symmetry[peft->index_J21112x] = sym_mat_symmetric; peft->use_tracer[peft->index_J21112x] = eft_matter; peft->spectra_contributions_dimension[peft->index_J21112x] = 1;
    peft->symmetry[peft->index_J21112y] = sym_mat_symmetric; peft->use_tracer[peft->index_J21112y] = eft_matter; peft->spectra_contributions_dimension[peft->index_J21112y] = 1;
    peft->symmetry[peft->index_J12112x] = sym_vec;           peft->use_tracer[peft->index_J12112x] = eft_halo;   peft->spectra_contributions_dimension[peft->index_J12112x] = 1;
    peft->symmetry[peft->index_J12112y] = sym_vec;           peft->use_tracer[peft->index_J12112y] = eft_halo;   peft->spectra_contributions_dimension[peft->index_J12112y] = 1;
    peft->symmetry[peft->index_N12x]    = sym_mat_none;      peft->use_tracer[peft->index_N12x]    = eft_matter; peft->spectra_contributions_dimension[peft->index_N12x]    = 0;
    peft->symmetry[peft->index_N12y]    = sym_mat_none;      peft->use_tracer[peft->index_N12y]    = eft_matter; peft->spectra_contributions_dimension[peft->index_N12y]    = 0;

    peft->symmetry[peft->index_N22x] = sym_mat_none; peft->use_tracer[peft->index_N22x] = eft_matter; peft->spectra_contributions_dimension[peft->index_N22x] = -1;
    peft->symmetry[peft->index_N22y] = sym_mat_none; peft->use_tracer[peft->index_N22y] = eft_matter; peft->spectra_contributions_dimension[peft->index_N22y] = -1;
    peft->symmetry[peft->index_N22z] = sym_mat_none; peft->use_tracer[peft->index_N22z] = eft_matter; peft->spectra_contributions_dimension[peft->index_N22z] = -1;

    peft->symmetry[peft->index_sigmav_mu] = no_finite_part; peft->use_tracer[peft->index_sigmav_mu] = eft_matter; peft->spectra_contributions_dimension[peft->index_sigmav_mu] = -2;
  }

  /** - loop kernel matrices */
  class_alloc(peft->loop_matrices,      peft->index_num*sizeof(double complex *), peft->error_message);
  class_alloc(peft->loop_matrices_size, peft->index_num*sizeof(int), peft->error_message);
  for (i = 0; i < peft->index_num; i++) { /** - i = index_moment */
    switch (peft->symmetry[i]) {
    case sym_vec:
      peft->loop_matrices_size[i] = peft->hp->fourier_coeff_size; break;
    case sym_mat_none:
      peft->loop_matrices_size[i] = peft->hp->fourier_coeff_size * peft->hp->fourier_coeff_size; break;
    case sym_mat_symmetric:
      peft->loop_matrices_size[i] = peft->hp->fourier_coeff_size * (peft->hp->fourier_coeff_size + 1) / 2; break;
    default:  /** also no_finite_part */
      peft->loop_matrices_size[i] = 0; break;
    }

    class_alloc(peft->loop_matrices[i], peft->loop_matrices_size[i]*sizeof(double complex), peft->error_message);
  }

  peft->moments_allocated = peft->index_num;

  return _SUCCESS_;
}


int eft_indices(struct eft * peft,
                const short eft_role) {

  int i = 0, index_pk_type, index_moment, index_part, index_tracer;

  peft->moments_allocated = 0;

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
      class_define_index(peft->index_J11201,    _TRUE_, i, 1);

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
      class_define_index(peft->index_J12111,    _TRUE_, i, 1);
      class_define_index(peft->index_J11211,    _TRUE_, i, 1);

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

      /** - additional RSD moments */
      class_define_index(peft->index_sigmav_mu, _TRUE_, i, 1);
    }

    peft->index_num = i;

    /** - k-grid and frequencies for Fourier transform */
    for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
      class_alloc(peft->ln_k_fourier[index_tracer],        peft->hp->k_size_fourier*sizeof(double), peft->error_message);
      class_alloc(peft->fourier_frequencies[index_tracer], peft->hp->fourier_coeff_size*sizeof(double), peft->error_message);
    }

    class_alloc(peft->ln_k_moments, peft->hp->k_size_moments*sizeof(double), peft->error_message);

    /** - initialize the output sampling grid */
    if (peft->hp->use_interpolation) {
      peft->k_size = peft->hp->k_size_nl;
      class_alloc(peft->ln_k, peft->k_size*sizeof(double), peft->error_message);
    }
    else {
      peft->ln_k = NULL;
    }

    peft->mu_size = 0;
  }

  /** - allocate and initialize spectra contributions and Fourier coefficient indicators */
  class_calloc(peft->spectra_contributions_size, pk_type_num*peft->index_num, sizeof(int), peft->error_message);
  for (index_pk_type = 0; index_pk_type < pk_type_num; index_pk_type++) {
    peft->pk_type_loaded[index_pk_type] = _FALSE_;
    class_alloc(peft->spectra_contributions[index_pk_type], peft->index_num*eft_spectra_contribution_num*sizeof(double *), peft->error_message);
    for (index_moment = 0; index_moment < peft->index_num; index_moment++) {
      for (index_part = 0; index_part < eft_spectra_contribution_num; index_part++) {
        peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + index_part] = NULL;
      }
    }
    for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
      peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer] = NULL;
      peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer] = NULL;
      peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer] = NULL;
      peft->fourier_condition_num[index_pk_type*eft_tracer_num + index_tracer] = NULL;
    }
    peft->pk_l[index_pk_type] = NULL;
    peft->ddpk_l[index_pk_type] = NULL;
    peft->pk_l_moments[index_pk_type] = NULL;
    peft->ddpk_l_moments[index_pk_type] = NULL;
  }

  /** - initialize the power spectrum moments */
  for (index_pk_type = 0; index_pk_type < pk_type_num; index_pk_type++) {
    for (index_moment = 0; index_moment < EFT_DISPERSION_SIZE; index_moment++) {
      peft->dispersion[index_pk_type][index_moment].index_bias = SHRT_MIN;
    }
    for (index_moment = 0; index_moment < EFT_UV_CORRECTIONS_SIZE; index_moment++) {
      peft->ps_uv_shot_noise_corrections[index_pk_type][index_moment].index_bias = SHRT_MIN;
      peft->ps_uv_shot_noise_corrections[index_pk_type][index_moment].index_derivative = SHRT_MIN;
    }
    for (index_moment = 0; index_moment < EFT_UV_CORRECTIONS_UNDERLYING_SIZE; index_moment++) {
      peft->ps_uv_shot_noise_corrections_underlying[index_pk_type][index_moment].index_bias = SHRT_MIN;
      peft->ps_uv_shot_noise_corrections_underlying[index_pk_type][index_moment].index_derivative = SHRT_MIN;
    }
  }


  return _SUCCESS_;
}

int eft_free(struct eft * peft) {
  int index_pk_type, index_tracer, index_moment, index_part;

  if (peft->role == eft_master) {
    free(peft->symmetry);
    free(peft->use_tracer);
    free(peft->spectra_contributions_dimension);

    if (peft->loop_matrices) {
      for (index_moment = 0; index_moment < peft->index_num; index_moment++) {
        free(peft->loop_matrices[index_moment]);
      }
    }
    free(peft->loop_matrices);
    free(peft->loop_matrices_size);

    for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
      free(peft->ln_k_fourier[index_tracer]);
      free(peft->fourier_frequencies[index_tracer]);
    }
    free(peft->ln_k_moments);
    free(peft->mu);
    free(peft->ln_k);

    if (peft->hp->fourier_mode == fourier_mode_fft) {
      FFT_planner_free(&(peft->fft_plan));
    }
  }

  for (index_pk_type = 0; index_pk_type < pk_type_num; index_pk_type++) {
    if (peft->pk_type_loaded[index_pk_type]) {
      for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
        free(peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer]);
        free(peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer]);
        free(peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer]);
        free(peft->fourier_condition_num[index_pk_type*eft_tracer_num + index_tracer]);
      }
      free(peft->pk_l_moments[index_pk_type]);
      free(peft->ddpk_l_moments[index_pk_type]);
    }

    for (index_moment = 0; index_moment < peft->index_num; index_moment++) {
      for (index_part = 0; index_part < eft_spectra_contribution_num; index_part++) {
        free(peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + index_part]);
      }
    }
    free(peft->spectra_contributions[index_pk_type]);

    free(peft->pk_l[index_pk_type]);
    free(peft->ddpk_l[index_pk_type]);
  }
  free(peft->spectra_contributions_size);

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

  int index_tracer, index_moment, index_pk_type;

  /** - register pointers to the EFT parameter structs */
  peft->hp = eft_hp;
  peft->ip = eft_ip;
  peft->role = eft_role;

  class_call(eft_indices(peft, eft_role), peft->error_message, peft->error_message);

  if (eft_role == eft_master) {
    for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
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

    /** Generate the sampling points for moment compuations */
    eft_spline_sample_points_nonuniform(peft->hp->k_IR_cutoff, exp(ppr->k_bao_center), peft->hp->k_UV_cutoff,   \
                                        peft->hp->bao_oversampling, peft->hp->ln_k_oversampling_width,          \
                                        peft->hp->k_size_moments, peft->ln_k_moments, peft->hp->eft_verbose);

    /** Generate the initial sampling points for the output spectrum (might be overwritten) */
    if (peft->hp->use_interpolation) {
      eft_spline_sample_points_nonuniform(peft->hp->kmin_nl, peft->hp->k_feature_nl, peft->hp->kmax_nl,       \
                                          peft->hp->bao_oversampling, peft->hp->ln_k_oversampling_width_nl,   \
                                          peft->hp->k_size_nl, peft->ln_k, peft->hp->eft_verbose);
    }
    /** - initialize the mu-grid to real-space, but will not be read if use_mu_approximation = _TRUE_ */
    peft->mu_size = 1;
    class_alloc(peft->mu, 1*sizeof(double), peft->error_message);
    peft->mu[0] = 0.;
  }


  return _SUCCESS_;
}

/** TODO: inline this */
int eft_linear_spectrum_real(
                  struct background * pba,
                  struct primordial * ppm,
                  struct fourier * pfo,
                  struct eft * peft,
                  enum linear_or_logarithmic mode,
                  const double * const ln_kvec,
                  const int kvec_size,
                  const int n_columns,
                  const double z,
                  const double f_z,
                  const double D_z,
                  const int index_pk_type,
                  double * out_pk) {

  int it, abort = _FALSE_;
  double * ln_kvec_sorted, * pk_l;
  struct indexed_real_arg * vec;
  const double D2 = pow(D_z/peft->D_z0, 2.);

  class_call_parallel(eft_real_argument_list_rect(ln_kvec, kvec_size,
                                                  n_columns,
                                                  &vec,
                                                  peft->error_message),
                      peft->error_message, peft->error_message);

  switch (index_pk_type)
  {
  case pk_lin:
    class_alloc_parallel(ln_kvec_sorted, kvec_size*n_columns*sizeof(double), peft->error_message);
    class_alloc_parallel(pk_l, kvec_size*n_columns*sizeof(double), peft->error_message);
    for (it = 0; it < kvec_size*n_columns; it++) {
      ln_kvec_sorted[it] = vec[it].ln_k;
    }
    class_call_parallel(fourier_pk_at_kvec_and_z(pba, ppm, pfo, mode, pk_linear,
                                                 ln_kvec_sorted,
                                                 kvec_size*n_columns,
                                                 z,
                                                 peft->hp->linear_spectrum_index,
                                                 pk_l),
                        pfo->error_message,
                        peft->error_message);
    for (it = 0; it < kvec_size*n_columns; it++) {
      out_pk[vec[it].index] = pk_l[it];
    }
    free(ln_kvec_sorted);
    free(pk_l);
    break;

  case pk_nowiggle:
    class_alloc_parallel(ln_kvec_sorted, kvec_size*n_columns*sizeof(double), peft->error_message);
    class_alloc_parallel(pk_l, kvec_size*n_columns*sizeof(double), peft->error_message);
    for (it = 0; it < kvec_size*n_columns; it++) {
      ln_kvec_sorted[it] = vec[it].ln_k;
    }
    class_call_parallel(fourier_pk_nw_at_kvec_and_z(pba, ppm, pfo, mode,
                                                    ln_kvec_sorted,
                                                    kvec_size*n_columns,
                                                    z,
                                                    pk_l),
                        pfo->error_message,
                        peft->error_message);
    for (it = 0; it < kvec_size*n_columns; it++) {
      out_pk[vec[it].index] = pk_l[it];
    }
    free(ln_kvec_sorted);
    free(pk_l);
    break;

  case pk_ir_resummed_lo:
    class_call_parallel(eft_ir_pk_lo(pba, ppm, pfo, mode,
                                     vec,
                                     kvec_size*n_columns,
                                     z,
                                     peft->hp->linear_spectrum_index,
                                     peft->Sigma2_ir * D2,
                                     out_pk),
                        pfo->error_message,
                        peft->error_message);
    break;

  case pk_ir_resummed_nlo:
    class_call_parallel(eft_ir_pk_nlo(pba, ppm, pfo, mode,
                                      vec,
                                      kvec_size*n_columns,
                                      z,
                                      peft->hp->linear_spectrum_index,
                                      peft->Sigma2_ir * D2,
                                      out_pk),
                        pfo->error_message,
                        peft->error_message);
    break;

  default:
    class_stop(pfo->error_message, "index_pk_type = %d not recognized.", index_pk_type);
    break;
  }

  free(vec);

  return abort;
}

// TODO: inline this
int eft_linear_spectrum_rsd(
                  struct background * pba,
                  struct primordial * ppm,
                  struct fourier * pfo,
                  struct eft * peft,
                  enum linear_or_logarithmic mode,
                  const double * const ln_kvec,
                  const int kvec_size,
                  const double * const muvec,
                  const int muvec_size,
                  const double z,
                  const double f_z,
                  const double D_z,
                  const int index_pk_type,
                  double * out_pk) {

  int abort = _FALSE_;
  const double D2 = pow(D_z/peft->D_z0, 2.);
  struct indexed_rsd_arg * vec;

  class_call_parallel(eft_rsd_argument_list_rect(ln_kvec, kvec_size,
                                                 muvec, muvec_size,
                                                 &vec,
                                                 peft->error_message),
                      peft->error_message, peft->error_message);

  switch (index_pk_type)
  {
  case pkmu_rsd_ir_resummed_lo:
    class_call_parallel(eft_ir_pk_rsd_lo(pba, ppm, pfo, mode,
                                         vec,
                                         muvec_size*kvec_size,
                                         z,
                                         f_z,
                                         peft->hp->linear_spectrum_index,
                                         peft->Sigma2_ir * D2,
                                         peft->dSigma2_ir * D2,
                                         out_pk),
                        pfo->error_message,
                        peft->error_message);
    break;

  case pkmu_rsd_ir_resummed_nlo:
    class_call_parallel(eft_ir_pk_rsd_nlo(pba, ppm, pfo, mode,
                                          vec,
                                          muvec_size*kvec_size,
                                          z,
                                          f_z,
                                          peft->hp->linear_spectrum_index,
                                          peft->Sigma2_ir * D2,
                                          peft->dSigma2_ir * D2,
                                          out_pk),
                        pfo->error_message,
                        peft->error_message);
    break;

  default:
    class_stop(pfo->error_message, "index_pk_type = %d not recognized.", index_pk_type);
    break;
  }

  free(vec);

  return abort;
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



int eft_fourier_transform_linear_spectra(
                            struct precision * ppr,
                            struct background * pba,
                            struct fourier * pfo,
                            struct primordial * ppm,
                            struct eft * peft,
                            const double z,
                            const double f,
                            const double D,
                            const int * const index_pk_types,
                            const int index_pk_types_size) {

  int it, index_mu, index_tracer, index_pk_type, index_list, abort = _FALSE_;
  int mu_size;
  // int * index_pk_types, index_pk_types_size = 0;
  const int num_independent_coefficients = peft->hp->num_positive_fourier_freq + 1;
  double ** pk_l_biased_p;
  double spline_lipschitz_const = 0., spline_gradient_change, spline_L2_norm, series_l2_norm = 0.;

  /** if the list is empty there is nothing to do */
  if (index_pk_types_size < 1) { return _SUCCESS_; }

  peft->z0 = z;
  peft->f_z0 = f; /**< growth rate may be supplied externally */
  peft->D_z0 = D;

  /** Compute the suppression factor for IR-Resummation */
  peft->Sigma2_ir  =  eft_ir_sigma2(pba, ppm, pfo, peft->z0, peft->hp->ir_resummation_k_split*pba->h, peft->hp->ir_resummation_k_feature*pba->h);  /** k_split = 0.2 h/Mpc, k_feature = 1/110 h/Mpc */
  peft->dSigma2_ir = eft_ir_dsigma2(pba, ppm, pfo, peft->z0, peft->hp->ir_resummation_k_split*pba->h, peft->hp->ir_resummation_k_feature*pba->h);  /** according to arXiV:1804.05080 by Ivanov & Sibiryakov */

  /** prepare a list of indices for the pk_types to compute */
  // for (index_pk_type = 0; index_pk_type < pk_type_num; index_pk_type++) {
  //   if (use_pk_types[index_pk_type] && !(peft->pk_type_loaded[index_pk_type])) {
  //     index_pk_types_size++;
  //   }
  // }
  // if (index_pk_types_size == 0) { return _SUCCESS_; /** - nothing to do */}
  // class_alloc(index_pk_types, index_pk_types_size*sizeof(int), peft->error_message);
  // index_pk_types_size = 0;
  // for (index_pk_type = 0; index_pk_type < pk_type_num; index_pk_type++) {
  //   if (use_pk_types[index_pk_type] && !(peft->pk_type_loaded[index_pk_type])) {
  //     index_pk_types[index_pk_types_size++] = index_pk_type;
  //   }
  // }
  /** list of indices for the pk_types to compute is prepared externally */


  #pragma omp parallel shared(ppr, pba, ppm, pfo, peft, abort, stderr, num_independent_coefficients, index_pk_types, index_pk_types_size, spline_lipschitz_const, spline_gradient_change, spline_L2_norm, series_l2_norm), \
                       private(it, index_mu, index_tracer, index_pk_type, index_list, mu_size, pk_l_biased_p), \
                       default(none)
  { /** , last_index, ln_k_fft, pk_l_biased_fft, fft_plan, fourier_coeff_real, fourier_coeff_imag */

  /** Retrieve the power spectra */
  if (peft->hp->use_mu_approximation) {
    /** ----------------- Load requested linear spectra w/o mu dependence ---------------------- */
    /** - for the Fourier components */
    #pragma omp for schedule(static, 1), collapse(2)
    for (index_list = 0; index_list < index_pk_types_size; index_list++) {
      for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
        index_pk_type = index_pk_types[index_list];
        /** - allocate spectra and fourier coefficient arrays */
        class_realloc_parallel(peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer], peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer],      \
                               peft->hp->k_size_fourier*sizeof(double), peft->error_message);
        class_realloc_parallel(peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer], peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer],  \
                               peft->hp->k_size_fourier*sizeof(double), peft->error_message);
        class_realloc_parallel(peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer], peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer],  \
                               peft->hp->fourier_coeff_size*sizeof(double complex), peft->error_message);
        class_realloc_parallel(peft->fourier_condition_num[index_pk_type*eft_tracer_num + index_tracer], peft->fourier_condition_num[index_pk_type*eft_tracer_num + index_tracer],  \
                               peft->hp->fourier_coeff_size*sizeof(double complex), peft->error_message);
        memset(peft->fourier_condition_num[index_pk_type*eft_tracer_num + index_tracer], 0, peft->hp->fourier_coeff_size*sizeof(double complex));
        /** - get the linear spectra values */
        class_call_parallel(eft_linear_spectrum_real(pba, ppm, pfo, peft, linear,
                                                     peft->ln_k_fourier[index_tracer],
                                                     peft->hp->k_size_fourier-1,
                                                     1,
                                                     peft->z0,
                                                     peft->f_z0,
                                                     peft->D_z0,
                                                     index_pk_type,
                                                     peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer]),
                            peft->error_message,
                            peft->error_message);
      }
    }
    /** - for the moments of the lin. spectra */
    #pragma omp for schedule(static, 1)
    for (index_list = 0; index_list < index_pk_types_size; index_list++) {
      index_pk_type = index_pk_types[index_list];
      /** - allocate spectra for moments */
      class_realloc_parallel(peft->pk_l_moments[index_pk_type], peft->pk_l_moments[index_pk_type],      \
                             peft->hp->k_size_moments*sizeof(double), peft->error_message);
      class_realloc_parallel(peft->ddpk_l_moments[index_pk_type], peft->ddpk_l_moments[index_pk_type],  \
                             peft->hp->k_size_moments*sizeof(double), peft->error_message);
      /** - get the linear spectra values */
      class_call_parallel(eft_linear_spectrum_real(pba, ppm, pfo, peft, linear,
                                                    peft->ln_k_moments,
                                                    peft->hp->k_size_moments,
                                                    1,
                                                    peft->z0,
                                                    peft->f_z0,
                                                    peft->D_z0,
                                                    index_pk_type,
                                                    peft->pk_l_moments[index_pk_type]),
                          peft->error_message,
                          peft->error_message);
    }


    if (!abort) {
      #pragma omp for schedule(static, 1), collapse(2)
      for (index_list = 0; index_list < index_pk_types_size; index_list++) {
        for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
          index_pk_type = index_pk_types[index_list];
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

      #pragma omp for schedule(static, 1)
      for (index_list = 0; index_list < index_pk_types_size; index_list++) {
        index_pk_type = index_pk_types[index_list];
        /** Spline the power spectra for moments with estimated boundary derivative */
        class_call_parallel(array_spline_table_columns(peft->ln_k_moments,
                                                      peft->hp->k_size_moments,
                                                      peft->pk_l_moments[index_pk_type],
                                                      1,
                                                      peft->ddpk_l_moments[index_pk_type],
                                                      _SPLINE_EST_DERIV_,
                                                      peft->error_message),
                            peft->error_message,
                            peft->error_message);

      }
    }

  }
  else {
    // TODO: is this needed?
    abort = _TRUE_;
    // /** ---------------- Load pk_ir_resummed_lo for different mu -------------- */
    // /** mu_size, mu are set before & dd/pk_l_biased[pk_ir_resummed_lo*eft_tracer_num + index_tracer] needs to be allocated with size peft->mu_size*peft->hp->k_size_fourier */
    // /** fourier_coeff/condition_num[pk_ir_resummed_lo*eft_tracer_num + index_tracer] needs to be allocated with size peft->mu_size*peft->hp->fourier_coeff_size */
    // class_alloc_parallel(pk_l_biased_p, peft->mu_size*sizeof(double **), peft->error_message);

    // #pragma omp for schedule(static, 1), collapse(2)
    // for (index_list = 0; index_list < index_pk_types_size; index_list++) {
    //   for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
    //     index_pk_type = index_pk_types[index_list];
    //     for (index_mu = 0; index_mu < peft->mu_size; index_mu++) {
    //       pk_l_biased_p[index_mu] = peft->pk_l_biased[pk_ir_resummed_lo*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier;
    //     }
    //     class_call_parallel(eft_linear_spectrum_rsd(pba, ppm, pfo, peft, linear,
    //                                                 peft->ln_k_fourier[index_tracer],
    //                                                 peft->hp->k_size_fourier-1,
    //                                                 peft->mu,
    //                                                 peft->mu_size,
    //                                                 peft->z0,
    //                                                 peft->f_z0,
    //                                                 peft->D_z0,
    //                                                 index_pk_type,
    //                                                 pk_l_biased_p),
    //                         peft->error_message,
    //                         peft->error_message);
    //   }
    // }
    // free(pk_l_biased_p);


    // if (!abort) {
    //   #pragma omp for schedule(static, 1), collapse(3)
    //   for (index_list = 0; index_list < index_pk_types_size; index_list++) {
    //     for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
    //       for (index_mu = 0; index_mu < peft->mu_size; index_mu++) {
    //         index_pk_type = index_pk_types[index_list];
    //         for (it = 0; it < peft->hp->k_size_fourier-1; it++) {
    //           /** Apply the bias to the power spectra */
    //           peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->k_size_fourier + it]      \
    //                 *= exp(-peft->hp->bias[index_tracer] * peft->ln_k_fourier[index_tracer][it]);
    //         }
    //         /** Add the last point to make the spectra periodic: P_bias(ln(k_min) + period) = P_bias(ln(k_min)) */
    //         peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->k_size_fourier + peft->hp->k_size_fourier-1]    \
    //               = peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][0];

    //         /** Spline the biased power spectra periodically */
    //         class_call_parallel(array_spline_table_columns(peft->ln_k_fourier[index_tracer],
    //                                                       peft->hp->k_size_fourier,
    //                                                       peft->pk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
    //                                                       1,
    //                                                       peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier,
    //                                                       _SPLINE_PERIODIC_,
    //                                                       peft->error_message),
    //                             peft->error_message,
    //                             peft->error_message);

    //       }
    //     }
    //   }
    // }

  }

  #pragma omp barrier

  /* Compute the Fourier coefficients of the power spectra */
  switch (peft->hp->fourier_mode)
  {
  case fourier_mode_spline:
    for (index_list = 0; index_list < index_pk_types_size; index_list++) {
      index_pk_type = index_pk_types[index_list];
      mu_size = (peft->hp->use_mu_approximation) ? 1 : peft->mu_size;
      for (index_tracer = 0; index_tracer < eft_tracer_num; index_tracer++) {
        for (index_mu = 0; index_mu < mu_size; index_mu++) {

          #pragma omp master
          {
          peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + 0] = class_complex(0., 0.);
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
            spline_lipschitz_const = 0.;
            series_l2_norm = 0.;
            }
            double h1, h2;
            double * M0 = peft->ddpk_l_biased[index_pk_type*eft_tracer_num + index_tracer] + index_mu*peft->hp->k_size_fourier;
            #pragma omp barrier
            /** - compute the Lipschitz constant L of the second derivative: since the (biased) spline is a C^2 Lipschitz-continuous function,
             *    the magnitude of its Fourier coefficients is bounded by pi/2 * L / frequency^3 */
            #pragma omp for schedule(static) reduction(max:spline_lipschitz_const)
            for (it = 0; it <= peft->hp->k_size_fourier; it++) {
              h1 = peft->ln_k_fourier[index_tracer][it % peft->hp->k_size_fourier] - peft->ln_k_fourier[index_tracer][(it-1) % peft->hp->k_size_fourier];
              h2 = peft->ln_k_fourier[index_tracer][(it+1) % peft->hp->k_size_fourier] - peft->ln_k_fourier[index_tracer][it % peft->hp->k_size_fourier];
              spline_gradient_change = fabs((M0[(it+1) % peft->hp->k_size_fourier] - M0[it % peft->hp->k_size_fourier])*h1    \
                                          - (M0[it % peft->hp->k_size_fourier] - M0[(it-1) % peft->hp->k_size_fourier])*h2) / (h1*h2);
              if (spline_gradient_change > spline_lipschitz_const) { spline_lipschitz_const = spline_gradient_change; }
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
            printf("                                                      Fourier coefficients are bounded by %.3e / frequency^3 \n", _PI_/2.*spline_lipschitz_const);
            }
          }

        }
      }
      if (!abort) { peft->pk_type_loaded[index_pk_type] = _TRUE_; }
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

    mu_size = (peft->hp->use_mu_approximation) ? 1 : peft->mu_size;

    #pragma omp for schedule(static, 1), collapse(2)
    for (index_list = 0; index_list < index_pk_types_size; index_list++) {
      for (index_mu = 0; index_mu < mu_size; index_mu++) {
        index_pk_type = index_pk_types[index_list];
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
              /** - TODO: Not really consistent! DC coefficient is halved to keep the original power while symmetrizing the spectrum btw. -N/2...N/2 */
              peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + 0]      \
                    = class_complex(fourier_coeff_real[index_tracer][0], fourier_coeff_imag[index_tracer][0]) / (double)(peft->hp->fourier_coeff_size-1);
              /** - positive frequency components in standard order, last iteration extends the spectrum symmetrically */
              for (it = 1; it < num_independent_coefficients; it++) {
                peft->fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*peft->hp->fourier_coeff_size + it]      \
                    = class_complex(fourier_coeff_real[index_tracer][it], fourier_coeff_imag[index_tracer][it]) * cpow(peft->hp->kmin_lin[index_tracer], -_Complex_I*peft->fourier_frequencies[index_tracer][it]) / (double)(peft->hp->fourier_coeff_size-1);
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
            class_protect_sprintf(errmsg, "FFT implementation requires 2 tracers for performance reasons, but the eft_tracer enum contains %d entries. If you need more tracers, implement the FFT schedule here.", eft_tracer_num);
            class_build_error_string(peft->error_message, "error; %s", errmsg);
          }
        }

        if (!abort) { peft->pk_type_loaded[index_pk_type] = _TRUE_; }
      }
    }


    free(pk_l_biased_fft[eft_matter]); free(pk_l_biased_fft[eft_halo]);
    free(fourier_coeff_real[0]); free(fourier_coeff_imag[0]);
    free(fourier_coeff_real[1]); free(fourier_coeff_imag[1]);
    FFT_planner_free(&fft_plan);
    break;
  }

  } /** - end of parallel region */

  if (peft->hp->use_mu_approximation && (pfo->fourier_verbose > 3)) {
    FILE *ffourier = fopen("output/halo_biased_pk_lin_samples.dat", "w");

    if (ffourier) {
      fprintf(ffourier, "# Biased samples of the linear power spectrum at z=%.3f \n", peft->z0);
      fprintf(ffourier, "# for k=%.4e to %.4e \n", exp(peft->ln_k_fourier[eft_halo][0]), exp(peft->ln_k_fourier[eft_halo][peft->hp->k_size_fourier-1]));
      fprintf(ffourier, "# number of wavenumbers equal to %d \n", peft->hp->k_size_fourier);
      fprintf(ffourier, "#    1:k (1/Mpc)            2:P_bias (Mpc)^3        3:d^2P_bias/dln(k)^2 (Mpc)^3 \n");
      for (int i = 0; i < peft->hp->k_size_fourier; i++)
        fprintf(ffourier, "  %.16e       %.16e       %+.16e \n", \
                exp(peft->ln_k_fourier[eft_halo][i]), peft->pk_l_biased[pk_lin*eft_tracer_num + eft_halo][i], peft->ddpk_l_biased[pk_lin*eft_tracer_num + eft_halo][i]);

      fclose(ffourier);
    }

    ffourier = fopen("output/halo_pk_lin_fourier_coefficients.dat", "w");

    if (ffourier) {
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
  const long info = ((long)size << 32) | (symmetry << 16) | tracer;
  long read_info;
  double read_period;
  short write;

  /** - check if the file exists  */
  file = fopen(filename, "rb");
  if (file) {
    /** - file already exists: compare metadata (tests for changes in size, symmetry, tracer species or period) */
    // read_count += fread(&read_info, sizeof(long), 1, file);
    // read_count += fread(&read_period, sizeof(double), 1, file);
    // if (read_count == 2) {
    //   /** - file is intact: overwrite if metadata changed and it is permitted */
    //   write = ((read_info == info) && (fabs(period - read_period)/period < _EPSILON_)) ? _FALSE_ : permit_overwrite;
    // }
    // else {
    //   /** - file is invalid: overwrite if permitted, otherwise fail */
    //   write = permit_overwrite;
    //   if (!write) { class_stop(errmsg, "Kernel matrix file %s exists, but is invalid and overwriting is forbidden.", filename); }
    // }
    write = permit_overwrite;

    if (write) {
      freopen(NULL, "wb+", file);
      if (!file) { class_stop(errmsg, "Kernel matrix file %s could not be created.", filename); }
    }
    else {
      fclose(file);
    }
  }
  else {
    /** - file does not exist: write is always permitted */
    write = _TRUE_;
    file = fopen(filename, "wb+");
    if (!file) { class_stop(errmsg, "Kernel matrix file %s could not be created.", filename); }
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

    if (ferror(file) || write_count < size+2) {
      /** - I/O error */
      fclose(file);
      class_stop(errmsg, "Writing kernel matrix in %s failed after %d of %d elements.", filename, write_count, size+2);
    }
    fclose(file);
    // else if (print_message) {
    //   printf("Saved kernel matrix in %s with size %s, symmetry type %d, tracer species %d and period %.3e \n", filename, size, symmetry, tracer, period);
    // }
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
                            const short print_message) {

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

  if (ferror(file) || read_count < size+2) {
    /** - I/O error */
    fclose(file);
    class_stop(errmsg, "Reading kernel matrix in %s failed after %d elements.", filename, read_count);
  }
  else {
    fclose(file);
    /** - recover metadata and compare to expected values */
    read_size = read_info >> 32; read_symmetry = read_info >> 16; read_tracer = read_info;
    if (read_size != size || read_symmetry != symmetry || read_tracer != tracer || (fabs(period - read_period)/period > _EPSILON_)) {
      class_stop(errmsg, "Stored kernel matrix in %s has size %d, symmetry type %d, tracer species %d and period %.3e; expected size %d, symmetry type %d, tracer species %d and period %.3e", \
                  filename, read_size, read_symmetry, read_tracer, read_period, size, symmetry, tracer, period);
    }

    // if (print_message) {
    //   printf("Read kernel matrix from %s with size %s, symmetry type %d, tracer species %d and period %.3e \n", filename, read_size, read_symmetry, read_tracer, read_period);
    // }
  }

  return _SUCCESS_;
}

int eft_get_loop_matrices(struct eft * peft,
                          struct ext_storage * pext,
                          const int index) {

  int index_moment, counter = 0, abort = _FALSE_;  /** - loop matrices for the current eft struct may need to be recomputed independently of the hyperparameters, if stored matrices are invalid */
  FileName filename;
  ErrorMsg errmsg;
  #ifdef _OPENMP
  int num_parallel_files = MIN(omp_get_max_threads(), FOPEN_MAX); /** - number of threads is either bounded by OMP_NUM_THREADS or max. concurrently opened files */
  #endif

  switch (peft->role)
  {
  case eft_master:
    /** master structure is responsible for writing the loop matrices to the arrays, while the slaves only read them */
    /** - external storage takes precedent: insert loop matrices from external storage */
    /** - condition here is negated: if insertion fails, allocate and compute matrices / load them from files */
    if (ext_insert_eft(pext, peft, index, peft->moments_allocated, errmsg) != _SUCCESS_) {
      if (pext && peft->hp->eft_verbose > 1) { printf("Could not insert loop matrices from external storage for index = %d: \n => %s \n", index, errmsg); }

      class_call(eft_allocate_loop_matrices(peft), peft->error_message, peft->error_message);
      /** - if external storage is not available/compatible, either */
      if (peft->hp->compute_loop_matrices) {  /** - compute the loop matrices or */
        class_call(eft_compute_loop_matrices(peft), peft->error_message, peft->error_message);

        if (peft->hp->write_loop_matrices == _TRUE_) {  /** - and save them if flag is set */
          #pragma omp parallel for schedule(dynamic), shared(peft, index, abort), private(index_moment, filename), default(none), num_threads(num_parallel_files)
          for (index_moment = 0; index_moment < peft->moments_allocated; index_moment++) {
            if (peft->loop_matrices_size[index_moment] == 0) { continue; }
            sprintf(filename, "%s_%03d.mat", peft->hp->eft_loop_matrix_files[index_moment], index);
            class_call_parallel(eft_save_matrix_to_file(peft->loop_matrices[index_moment],
                                                        peft->loop_matrices_size[index_moment],
                                                        peft->symmetry[index_moment],
                                                        peft->use_tracer[index_moment],
                                                        peft->hp->period[peft->use_tracer[index_moment]],
                                                        filename,
                                                        _TRUE_,
                                                        peft->error_message,
                                                        (peft->hp->eft_verbose > 2) ? _TRUE_ : _FALSE_),
                                peft->error_message,
                                peft->error_message);
          }
          if (!abort && peft->hp->eft_verbose > 1) {
            printf("Wrote %d loop matrix files to '%s' \n", peft->moments_allocated, peft->hp->eft_loop_matrix_directory);
          }
        }
      }
      else {  /** - load the matrices from files without raising errors at this point if flag is set */
        #pragma omp parallel for schedule(dynamic), shared(peft, index, abort), private(index_moment, filename), default(none), num_threads(num_parallel_files), reduction(+:counter)
        for (index_moment = 0; index_moment < peft->moments_allocated; index_moment++) {
          if (peft->loop_matrices_size[index_moment] == 0) { continue; }
          sprintf(filename, "%s_%03d.mat", peft->hp->eft_loop_matrix_files[index_moment], index);
          // class_call_parallel(eft_read_matrix_from_file(peft->loop_matrices[index_moment],
          //                                               peft->loop_matrices_size[index_moment],
          //                                               peft->symmetry[index_moment],
          //                                               peft->use_tracer[index_moment],
          //                                               peft->hp->period[peft->use_tracer[index_moment]],
          //                                               filename,
          //                                               peft->error_message,
          //                                               (peft->hp->eft_verbose > 2) ? _TRUE_ : _FALSE_),
          //                     peft->error_message,
          //                     peft->error_message);

          if (!abort) {
            if (eft_read_matrix_from_file(peft->loop_matrices[index_moment],
                                          peft->loop_matrices_size[index_moment],
                                          peft->symmetry[index_moment],
                                          peft->use_tracer[index_moment],
                                          peft->hp->period[peft->use_tracer[index_moment]],
                                          filename,
                                          peft->error_message,
                                          (peft->hp->eft_verbose > 2) ? _TRUE_ : _FALSE_) == _FAILURE_) {
              if (!peft->hp->ignore_missing_files) {
                class_call_message(peft->error_message, "eft_get_loop_matrices", peft->error_message);
                abort = _TRUE_;
              }
              else {
                peft->loop_matrices_size[index_moment] = 0;
              }
            }
            else {
              counter += 1;
            }
          }
        }
        /** - update the number of allocated moments */
        peft->moments_allocated = counter;
        if (!abort && peft->hp->eft_verbose > 2) {
          printf("Loaded %d loop matrix files out of %d defined indices from '%s' \n", counter, peft->index_num, peft->hp->eft_loop_matrix_directory);
        }
      }

    }

    break;

  case eft_slave:
    /** no-op */
    break;
  default:
    class_stop(peft->error_message, "Role of EFT struct was invalid: %d for index %d", peft->role, index);
    break;
  }

  return abort;
}


int eft_apply_ap_effect_in_place(double ** kvec,
                                 const int * k_sizevec,
                                 double ** muvec,
                                 const int * mu_sizevec,
                                 const int z_size,
                                 const double * ap_parallel,
                                 const double * ap_perpendicular) {

  int index_z, index_mu, index_k;
  double ap_ratio, sqrt_factor;

  for (index_z = 0; index_z < z_size; index_z++) {
    ap_ratio = ap_perpendicular[index_z] / ap_parallel[index_z];
    for (index_mu = 0; index_mu < mu_sizevec[index_z]; index_mu++) {
      sqrt_factor = sqrt(1. + (ap_ratio*ap_ratio - 1.) * muvec[index_z][index_mu]); // here muvec must still be fiducial
      for (index_k = 0; index_k < k_sizevec[index_z]; index_k++) {
        kvec[index_z][index_mu*k_sizevec[index_z] + index_k] *= 1./ap_perpendicular[index_z] * sqrt_factor;
      }
      muvec[index_z][index_mu] *= ap_ratio / sqrt_factor;
    }
  }

  return _SUCCESS_;
}


double eft_temporal_distance(struct background * pba,
                             const double z,
                             const double z_eft) {

  /** - just linear distance in redshift */
  return fabs(z - z_eft);
}

static int indexed_double_cmp_inc(const void * a, const void * b) {
  struct indexed_double * a_ = (struct indexed_double *)a;
  struct indexed_double * b_ = (struct indexed_double *)b;
  if ((*a_).value < (*b_).value)
      return -1;
  else if ((*a_).value > (*b_).value)
      return 1;
  else
      return 0;
}



/**
 * @brief Computes power-spectrum wedges at given redshifts, wavenumbers and l.o.s. angles
 *        using a given growth rate consistently.
 *
 * @param peft0       Input: pointer to the first eft-structure
 * @param peft_size   Input: number of eft-structures
 * @param f_z_pk_eft  Input: logarithmic growth rate at z_pk of each eft-structure
 * @param D_z_pk_eft  Input: growth rate D(z) at z_pk of each eft-structure
 * @param pba         Input: pointer to the background structure
 * @param pfo         Input: pointer to the fourier structure
 * @param ppm         Input: pointer to the primordial structure
 * @param ppr         Input: pointer to the precisions structure
 * @param pk_out_type Input: type of power-spectrum to compute
 * @param zvec        Input: redshifts at which to compute it
 * @param f_zvec      Input: logarithmic growth rate at each z in zvec
 * @param D_zvec      Input: growth rate D(z) at each z in zvec
 * @param peft_ip     Input: Input parameter structures for each eft-structure
 * @param z_size      Input: size of zvec
 * @param kvec        Input: wavenumbers [in 1/Mpc]; indexed as kvec[index_z][index_mu*k_sizevec[index_z] + index_k]
 * @param k_sizevec   Input: size of the wavenumber array for each mu at index_z; indexed as k_sizevec[index_z]
 * @param muvec       Input: cosine of line-of-sight angles; indexed as muvec[index_z][index_mu]
 * @param mu_sizevec  Input: size of muvec at each index_z; indexed as mu_sizevec[index_z]
 * @param out_pkmu    Output: power-spectrum at zvec,muvec,kvec; indexed as out_pkmu[index_z][index_mu*k_sizevec[index_z] + index_k];
 *                            needs to be pre-allocated with size k_sizevec[index_z]*mu_sizevec[index_z] at each index_z
 *
 * @return the error status
 */
int eft_job_powerspectrum_wedges_ext_growth_rate(
                                struct eft * peft0,
                                const int peft_size,
                                const double * const f_z_pk_eft,
                                const double * const D_z_pk_eft,
                                struct background * pba,
                                struct fourier * pfo,
                                struct primordial * ppm,
                                struct precision * ppr,
                                enum eft_pk_out_type pk_out_type,
                                const double * const zvec,
                                const double * const f_zvec,
                                const double * const D_zvec,
                                const struct eft_input_parameters * peft_ip,
                                const int z_size,
                                double ** kvec,
                                const int * const k_sizevec,
                                double ** muvec,
                                const int * const mu_sizevec,
                                double ** out_pkmu
                                ) {

  int it, index_z, index_z_sort, index_mu, index_k, index_eft, index_eft_min_dist[z_size], k_size, mu_size, sub_index_zvec[z_size], sub_index_z_size = 0, sub_z_size, sub_indexvec[z_size], sub_zvec[z_size], sub_f_zvec[z_size], sub_D_zvec[z_size];
  int * list_spectra_contributions, list_spectra_contributions_size, * list_pk_types_loops, list_pk_types_loops_size;
  int * list_pk_types_loops_not_loaded, list_pk_types_loops_not_loaded_size, * list_spectra_contributions_not_loaded, list_spectra_contributions_not_loaded_size;
  int index_pk_type, index_spectra_contribution, index_list, last_index;
  double z, f_z, D_z, min_distance, distance, mu_real = 0., * pkmu_nl = NULL, * ddpkmu_nl = NULL, * pkmu_out;
  struct eft * peft;
  struct eft_hyper_parameters * eft_hp;
  struct indexed_double sort_arr[z_size];
  int sorted_indexvec[z_size];
  double sorted_zvec[z_size], sorted_f_zvec[z_size], sorted_D_zvec[z_size];

  /** - sort the redshift array in increasing order */
  for (index_z = 0; index_z < z_size; index_z++) {
    sort_arr[index_z].index = index_z;
    sort_arr[index_z].value = zvec[index_z];
  }
  qsort(sort_arr, z_size, sizeof(struct indexed_double), indexed_double_cmp_inc);
  for (index_z = 0; index_z < z_size; index_z++) {
    sorted_indexvec[index_z] = sort_arr[index_z].index;
    sorted_zvec[index_z]   = sort_arr[index_z].value;
    sorted_f_zvec[index_z] = f_zvec[sort_arr[index_z].index];
    sorted_D_zvec[index_z] = D_zvec[sort_arr[index_z].index];
  }

  /** - find the neasrest eft structures for each redshift */
  for (index_z_sort = 0; index_z_sort < z_size; index_z_sort++) {
    z = sorted_zvec[index_z_sort];

    /** - find nearest eft structure */
    index_eft_min_dist[index_z_sort] = 0;
    min_distance = eft_temporal_distance(pba, z, pfo->z_pk_eft[0]);
    for (index_eft = 1; index_eft < pfo->z_pk_eft_num; index_eft++) {
      distance = eft_temporal_distance(pba, z, pfo->z_pk_eft[index_eft]);
      if (distance < min_distance) { min_distance = distance; index_eft_min_dist[index_z_sort] = index_eft; }
    }
  }

  for (index_z_sort = z_size-1; index_z_sort >= 0; index_z_sort--) {  /** - main z-loop */
    index_eft = index_eft_min_dist[index_z_sort];
    // if (eft_hp->use_EdS_time_scaling) { /** if we assume EdS time scaling, different redshifts assigned to the same struct can be evaluated concurrently */
    //   if (index_z_sort > 0 && index_eft == index_eft_min_dist[index_z_sort-1]) {
    //     sub_index_z_sortvec[sub_index_z_sort_size++] = index_z_sort; /** - add this redshift to the list */
    //     continue;
    //   }
    //   else {
    //     sub_z_size = sub_index_z_size; sub_index_z_size = 0;
    //     for (it = 0; it < sub_z_size; it++) { /** - build the list of redshifts for this eft struct */
    //       sub_indexvec[it] = sorted_indexvec[sub_index_zvec[it]];
    //       sub_zvec[it]   = sorted_zvec[sub_index_zvec[it]];
    //       sub_f_zvec[it] = sorted_f_zvec[sub_index_zvec[it]];
    //       sub_D_zvec[it] = sorted_D_zvec[sub_index_zvec[it]];
    //     }
    //   }
    // }
    // else {  /** else, everything will be reloaded */
    //   sub_z_size = 1;
    //   sub_indexvec[0] = sorted_indexvec[index_z_sort];
    //   sub_zvec[0]   = sorted_zvec[index_z_sort];
    //   sub_f_zvec[0] = sorted_f_zvec[index_z_sort];
    //   sub_D_zvec[0] = sorted_D_zvec[index_z_sort];
    // }

    peft = peft0 + index_eft;
    eft_hp = peft->hp;
    index_z = sorted_indexvec[index_z_sort];

    /** - if interpolation is deactivated, set the output sampling points for each redshift */

    if (!eft_hp->use_interpolation) {
      class_call(eft_set_sampling_points(peft,
                                         kvec[index_z],
                                         muvec[index_z],
                                         k_sizevec[index_z],
                                         mu_sizevec[index_z]),
                  peft->error_message, peft->error_message);
    }
    else if (!eft_hp->use_mu_approximation) {
      class_call(eft_set_sampling_points_mu_only(peft,
                                                 muvec[index_z],
                                                 mu_sizevec[index_z]),
                  peft->error_message, peft->error_message);
    }

    class_call(eft_necessary_spectra_contributions(peft,
                                                   pk_out_type,
                                                   &list_spectra_contributions,
                                                   &list_spectra_contributions_size),
               peft->error_message, peft->error_message);

    class_call(eft_necessary_pk_types_loops(peft,
                                            pk_out_type,
                                            &list_pk_types_loops,
                                            &list_pk_types_loops_size),
               peft->error_message, peft->error_message);

    /** - compile a list of pk_types of which to compute the Fourier transform,
     *    if use_interpolation = FALSE then this list will contain all of list_pk_types_loops */

    class_alloc(list_pk_types_loops_not_loaded, list_pk_types_loops_size*sizeof(int), peft->error_message);
    list_pk_types_loops_not_loaded_size = 0;
    for (index_list = 0; index_list < list_pk_types_loops_size; index_list++) {
      index_pk_type = list_pk_types_loops[index_list];
      if (!(eft_hp->use_interpolation && peft->pk_type_loaded[index_pk_type])) {
        list_pk_types_loops_not_loaded[list_pk_types_loops_not_loaded_size++] = index_pk_type;
      }
    }

    class_alloc(list_spectra_contributions_not_loaded, list_spectra_contributions_size*sizeof(int), peft->error_message);
    list_spectra_contributions_not_loaded_size = 0;
    for (index_list = 0; index_list < list_spectra_contributions_size; index_list++) {
      index_spectra_contribution = list_spectra_contributions[index_list];
      if (!(eft_hp->use_interpolation && peft->spectra_contributions_size[index_spectra_contribution] > 0)) {
        list_spectra_contributions_not_loaded[list_spectra_contributions_not_loaded_size++] = index_spectra_contribution;
      }
    }

    /* if the kernels have been precomputed for a specific redshift, we have to load spectra at this exact redshift */
    if (!eft_hp->use_time_independent_kernels) {
      z = pfo->z_pk_eft[index_eft];
      f_z = f_z_pk_eft[index_eft];
      D_z = D_z_pk_eft[index_eft];
    }
    /* else, load the spectra at the latest redshift in sub_zvec */
    else {
      z = zvec[index_z];
      f_z = f_zvec[index_z];
      D_z = D_zvec[index_z];
    }

    if (list_pk_types_loops_not_loaded_size > 0 && peft->hp->integration_mode == fftlog) {
      /** - load the real-space linear spectra */
      class_call(eft_load_linear_spectra(pba, pfo, ppm,
                                         peft,
                                         z,
                                         f_z,
                                         D_z,
                                         list_pk_types_loops_not_loaded,
                                         list_pk_types_loops_not_loaded_size,
                                         &mu_real,
                                         1),
                  peft->error_message, peft->error_message);

      /** - compute the Fourier transforms */
      class_call(eft_fourier_transform_linear_spectra(ppr, pba, pfo, ppm,
                                                      peft,
                                                      z,
                                                      f_z,
                                                      D_z,
                                                      list_pk_types_loops_not_loaded,
                                                      list_pk_types_loops_not_loaded_size),
                  peft->error_message, peft->error_message);
    }

    /** - allocate arrays for the necessary spectra_contributions */
    class_call(eft_allocate_spectra_contributions(peft,
                                                  (eft_hp->use_interpolation && eft_hp->integration_mode == fftlog) ? peft->k_size : k_sizevec[index_z]*mu_sizevec[index_z],
                                                  list_spectra_contributions,
                                                  list_spectra_contributions_size),
                peft->error_message, peft->error_message);

    if (eft_hp->integration_mode == fftlog) {
      /** - then compute the spectra_contributions that are not yet loaded */
      class_call(eft_compute_spectra_contributions(peft,
                                                  list_spectra_contributions_not_loaded,
                                                  list_spectra_contributions_not_loaded_size),
                  peft->error_message, peft->error_message);

      /** - compute the divergent parts */
      class_call(eft_compute_divergences(peft,
                                        list_spectra_contributions_not_loaded,
                                        list_spectra_contributions_not_loaded_size),
                  peft->error_message, peft->error_message);
    }
    else if (eft_hp->integration_mode == direct_integration) {
      #ifdef DIRECT_INTEGRATION
      class_call(eft_di_compute_spectra_contributions(peft,
                                                      pba,
                                                      pfo,
                                                      ppm,
                                                      ppr,
                                                      list_spectra_contributions_not_loaded,
                                                      list_spectra_contributions_not_loaded_size,
                                                      z,
                                                      f_z,
                                                      D_z),
                  peft->error_message, peft->error_message);
      #else
      class_stop(peft->error_message, "You have requested direct integration, but the associated module was not compiled!");
      #endif
    }
    else { class_stop(peft->error_message, "EFT integration mode not recognized, was %d", eft_hp->integration_mode); }

    if (eft_hp->use_interpolation) {
      class_realloc(pkmu_nl, pkmu_nl, peft->k_size*mu_sizevec[index_z]*sizeof(double),
                    peft->error_message);
      class_realloc(ddpkmu_nl, ddpkmu_nl, peft->k_size*mu_sizevec[index_z]*sizeof(double),
                    peft->error_message);
      pkmu_out = pkmu_nl;
    }
    else {
      pkmu_out = out_pkmu[index_z];
    }

    /** - build the power spectrum; this will always be done */
    class_call(eft_build_nonlinear_power_spectrum_wedges(peft, pba, ppm, pfo,
                                                         pk_out_type,
                                                         sorted_zvec[index_z_sort],
                                                         sorted_D_zvec[index_z_sort],
                                                         sorted_f_zvec[index_z_sort],
                                                         muvec[index_z],
                                                         mu_sizevec[index_z],
                                                         peft_ip[index_z],
                                                         pkmu_out),
                peft->error_message, peft->error_message);

    /** - spline the output power spectrum and interpolate */
    if (eft_hp->use_interpolation) {
      class_call(array_spline_table_columns_parallel(peft->ln_k,
                                                     peft->k_size,
                                                     pkmu_nl,
                                                     mu_sizevec[index_z],
                                                     ddpkmu_nl,
                                                     _SPLINE_EST_DERIV_,
                                                     peft->error_message),
                  peft->error_message, peft->error_message);

      #pragma omp parallel for schedule(static), shared(peft, pkmu_nl, ddpkmu_nl, out_pkmu, kvec, k_sizevec, mu_sizevec, index_z), private(index_mu, index_k, last_index), default(none)
      for (index_mu = 0; index_mu < mu_sizevec[index_z]; index_mu++) {
        last_index = 0;
        for (index_k = 0; index_k < k_sizevec[index_z]; index_k++) {
          array_interpolate_spline_growing_closeby(peft->ln_k,
                                                   peft->k_size,
                                                   pkmu_nl + index_mu*peft->k_size,
                                                   ddpkmu_nl + index_mu*peft->k_size,
                                                   1,
                                                   log(kvec[index_z][index_mu*k_sizevec[index_z] + index_k]),
                                                   &last_index,
                                                   out_pkmu[index_z] + index_mu*k_sizevec[index_z] + index_k,
                                                   1,
                                                   peft->error_message);
        }
      }
    }

    free(list_pk_types_loops);
    free(list_pk_types_loops_not_loaded);
    free(list_spectra_contributions);
    free(list_spectra_contributions_not_loaded);

  } /** - end of z-loop */

  free(pkmu_nl);
  free(ddpkmu_nl);

  return _SUCCESS_;
}


/**
 * @brief Computes power-spectrum wedges at given redshifts, wavenumbers and l.o.s. angles
 *        using internal CLASS growth rates.
 *
 * @param peft0       Input: pointer to the first eft-structure
 * @param peft_size   Input: number of eft-structures
 * @param pba         Input: pointer to the background structure
 * @param pfo         Input: pointer to the fourier structure
 * @param ppm         Input: pointer to the primordial structure
 * @param ppr         Input: pointer to the precisions structure
 * @param pk_out_type Input: type of power-spectrum to compute
 * @param zvec        Input: redshifts at which to compute it
 * @param peft_ip     Input: Input parameter structures for each eft-structure
 * @param z_size      Input: size of zvec
 * @param kvec        Input: wavenumbers [in 1/Mpc]; indexed as kvec[index_z][index_mu*k_sizevec[index_z] + index_k]
 * @param k_sizevec   Input: size of the wavenumber array for each mu at index_z; indexed as k_sizevec[index_z]
 * @param muvec       Input: cosine of line-of-sight angles; indexed as muvec[index_z][index_mu]
 * @param mu_sizevec  Input: size of muvec at each index_z; indexed as mu_sizevec[index_z]
 * @param out_pkmu    Output: power-spectrum at zvec,muvec,kvec; indexed as out_pkmu[index_z][index_mu*k_sizevec[index_z] + index_k];
 *                            needs to be pre-allocated with size k_sizevec[index_z]*mu_sizevec[index_z] at each index_z
 *
 * @return the error status
 */
int eft_job_powerspectrum_wedges(struct eft * peft0,
                                 const int peft_size,
                                 struct background * pba,
                                 struct fourier * pfo,
                                 struct primordial * ppm,
                                 struct precision * ppr,
                                 enum eft_pk_out_type pk_out_type,
                                 const double * const zvec,
                                 const struct eft_input_parameters * peft_ip,
                                 const int z_size,
                                 double ** kvec,
                                 const int * const k_sizevec,
                                 double ** muvec,
                                 const int * const mu_sizevec,
                                 double ** out_pkmu
                                 ) {

  int index_z, last_index = pba->bt_size-1;
  double * pvecback;
  double f_z_pk_eft[pfo->z_pk_eft_num], D_z_pk_eft[pfo->z_pk_eft_num];
  double f_zvec[z_size], D_zvec[z_size];

  class_alloc(pvecback, pba->bg_size*sizeof(double), pfo->error_message);

  for (index_z = pfo->z_pk_eft_num-1; index_z >= 0; index_z--) {
    class_call(background_at_z(pba,
                               pfo->z_pk_eft[index_z],
                               long_info,
                               inter_normal,
                               &last_index,
                               pvecback),
                pba->error_message,
                pfo->error_message);

    f_z_pk_eft[index_z] = pvecback[pba->index_bg_f];
    D_z_pk_eft[index_z] = pvecback[pba->index_bg_D];
  }

  for (index_z = z_size-1; index_z >= 0; index_z--) {
    class_call(background_at_z(pba,
                               zvec[index_z],
                               long_info,
                               inter_normal,
                               &last_index,
                               pvecback),
                pba->error_message,
                pfo->error_message);

    f_zvec[index_z] = pvecback[pba->index_bg_f];
    D_zvec[index_z] = pvecback[pba->index_bg_D];
  }

  free(pvecback);

  class_call(eft_job_powerspectrum_wedges_ext_growth_rate(peft0,
                                                          peft_size,
                                                          f_z_pk_eft,
                                                          D_z_pk_eft,
                                                          pba,
                                                          pfo,
                                                          ppm,
                                                          ppr,
                                                          pk_out_type,
                                                          zvec,
                                                          f_zvec,
                                                          D_zvec,
                                                          peft_ip,
                                                          z_size,
                                                          kvec,
                                                          k_sizevec,
                                                          muvec,
                                                          mu_sizevec,
                                                          out_pkmu),
              peft0->error_message, peft0->error_message);

  return _SUCCESS_;
}

/* same as previous function with different input/output format (flattened arrays with lower rank pointers, assuming that the number of k and mu values is independent of z) */

int eft_job_powerspectrum_wedges_grid(struct eft * peft0,
                                      const int peft_size,
                                      struct background * pba,
                                      struct fourier * pfo,
                                      struct primordial * ppm,
                                      struct precision * ppr,
                                      enum eft_pk_out_type pk_out_type,
                                      const double * const z, // indexed as z[index_z]
                                      const struct eft_input_parameters * peft_ip,
                                      const int z_size,
                                      double * k, // indexed as k[index_z + z_size*(index_mu + mu_size*index_k)]
                                      int k_size,
                                      double * mu, // indexed as mu[index_z + z_size*index_mu]
                                      int mu_size,
                                      double * out_pkmuz // indexed as out_pkmu[index_z + z_size*(index_mu + mu_size*index_k)], already allocated
                                      ) {

  double ** kvec;
  double ** muvec;
  int * k_sizevec;
  int * mu_sizevec;
  int index_z,index_k,index_mu;
  double ** out_pkmu;

  class_alloc(kvec,z_size*sizeof(double*),peft0->error_message);
  class_alloc(k_sizevec,z_size*sizeof(int),peft0->error_message);
  class_alloc(muvec,z_size*sizeof(double*),peft0->error_message);
  class_alloc(mu_sizevec,z_size*sizeof(int),peft0->error_message);

  for (index_z=0;index_z<z_size;index_z++) {
    k_sizevec[index_z]=k_size;
    class_alloc(kvec[index_z],k_size*mu_size*sizeof(double),peft0->error_message);
    for (index_k=0;index_k<k_size;index_k++) {
      for (index_mu=0;index_mu<mu_size;index_mu++) {
        kvec[index_z][index_mu*k_sizevec[index_z]+index_k] = k[index_z + z_size*(index_mu + mu_size*index_k)];
      }
    }

    mu_sizevec[index_z]=mu_size;
    class_alloc(muvec[index_z],mu_size*sizeof(double),peft0->error_message);
    for (index_mu=0;index_mu<mu_size;index_mu++) {
      muvec[index_z][index_mu] = mu[index_z + z_size*index_mu];
    }
  }

  class_alloc(out_pkmu,z_size*sizeof(double*),peft0->error_message);
  for (index_z=0;index_z<z_size;index_z++) {
    class_alloc(out_pkmu[index_z],k_size*mu_size*sizeof(double),peft0->error_message);
  }

  class_call(eft_job_powerspectrum_wedges(peft0,
                                          peft_size,
                                          pba,
                                          pfo,
                                          ppm,
                                          ppr,
                                          pk_out_type,
                                          z,
                                          peft_ip,
                                          z_size,
                                          kvec,
                                          k_sizevec,
                                          muvec,
                                          mu_sizevec,
                                          out_pkmu
                                          ),
             peft0->error_message, peft0->error_message);

  for (index_z=0;index_z<z_size;index_z++) {
    for (index_k=0;index_k<k_size;index_k++) {
      for (index_mu=0;index_mu<mu_size;index_mu++) {
        out_pkmuz[index_z + z_size*(index_mu + mu_size*index_k)]=out_pkmu[index_z][k_sizevec[index_z]*index_mu+index_k];
      }
    }
  }
  for (index_z=0;index_z<z_size;index_z++) {
      free(out_pkmu[index_z]);
  }
  free(out_pkmu);

  return _SUCCESS_;
}

#define MULTIPOLE_SIZE 3
// static const int gl_rule_size = 9;
// static const double gl_abscissa[] = {-1., -0.89975799541146015731, -0.67718627951073775345, -0.36311746382617815871, 0., 0.36311746382617815871, 0.67718627951073775345, 0.89975799541146015731, 1.};
// static const double gl_weights[] = {0.027777777777777777778, -0.067801868533905777326, 0.087327403186094540303, -0.098096963223617334406, 0.10158730158730158730, -0.098096963223617334406, 0.087327403186094540303, -0.067801868533905777326, 0.027777777777777777778};
static const int gl_rule_sym_size = 5;
static const double gl_sym_abscissa[] = {0., 0.36311746382617815871, 0.67718627951073775345, 0.89975799541146015731, 1.};
static const double gl_sym_weights[] = {0.10158730158730158730, -0.19619392644723466881, 0.17465480637218908061, -0.13560373706781155465, 0.055555555555555555556};
static const double lg_measure[][MULTIPOLE_SIZE] = { {1., -0.5, 0.375}, {1., -0.30221856119666629454, -0.043391796245607083376}, {1., 0.18787188573639255829, -0.42463134814493003484}, {1., 0.71434667546027373625, 0.20648468285207557975}, {1., 1., 1.} };

/**
 * @brief Computes power-spectrum multipoles at given redshifts and fiducial wavenumbers
 *        using a given growth rate consistently.
 *
 * @param peft0       Input: pointer to the first eft-structure
 * @param peft_size   Input: number of eft-structures
 * @param f_z_pk_eft  Input: logarithmic growth rate at z_pk of each eft-structure
 * @param D_z_pk_eft  Input: growth rate D(z) at z_pk of each eft-structure
 * @param pba         Input: pointer to the background structure
 * @param pfo         Input: pointer to the fourier structure
 * @param ppm         Input: pointer to the primordial structure
 * @param ppr         Input: pointer to the precisions structure
 * @param pk_out_type Input: type of power-spectrum to compute
 * @param zvec        Input: redshifts at which to compute it
 * @param f_zvec      Input: logarithmic growth rate at each z in zvec
 * @param D_zvec      Input: growth rate D(z) at each z in zvec
 * @param peft_ip     Input: Input parameter structures for each eft-structure
 * @param z_size      Input: size of zvec
 * @param kvec        Input: wavenumbers [in 1/Mpc]; indexed as kvec[index_z][index_k]
 * @param k_sizevec   Input: size of the wavenumber array for each mu at index_z; indexed as k_sizevec[index_z]
 * @param ap_parallel Input: parallel AP-effect ratio at each z; defined as H^fid(z)/H^true(z)
 * @param ap_perpendicular  Input: perpendicular AP-effect ratio at each z; defined as D_A^true(z)/D_A^fid(z)
 * @param out_pkl     Output: power-spectrum multipoles of order 0,2,4,...,2*(MULTIPOLE_SIZE-1) at zvec,kvec; indexed as out_pkl[index_z][(l/2)*k_sizevec[index_z] + index_k];
 *                            needs to be pre-allocated with size k_sizevec[index_z]*MULTIPOLE_SIZE at each index_z
 *
 * @return the error status
 */
int eft_job_powerspectrum_multipoles_ext_growth_rate(
                                struct eft * peft0,
                                const int peft_size,
                                const double * const f_z_pk_eft,
                                const double * const D_z_pk_eft,
                                struct background * pba,
                                struct fourier * pfo,
                                struct primordial * ppm,
                                struct precision * ppr,
                                enum eft_pk_out_type pk_out_type,
                                const double * const zvec,
                                const double * const f_zvec,
                                const double * const D_zvec,
                                const struct eft_input_parameters * peft_ip,
                                const int z_size,
                                double ** kvec,
                                const int * const k_sizevec,
                                const double * ap_parallel,
                                const double * ap_perpendicular,
                                double ** out_pkl
                                ) {

  int index_z, index_l, index_mu, index_k, mu_sizevec[z_size];
  double ap_ratio, * kvec_true[z_size], * muvec_true[z_size], * out_pkmu[z_size];

  /** - allocate additional input/output arrays for power-spectrum wedges */
  for (index_z = 0; index_z < z_size; index_z++) {
    mu_sizevec[index_z] = gl_rule_sym_size;
    class_alloc(muvec_true[index_z], mu_sizevec[index_z]*sizeof(double), peft0->error_message);
    class_alloc(kvec_true[index_z], k_sizevec[index_z]*mu_sizevec[index_z]*sizeof(double), peft0->error_message);
    class_alloc(out_pkmu[index_z], k_sizevec[index_z]*mu_sizevec[index_z]*sizeof(double), peft0->error_message);

    /** - copy the fiducial wavenumbers and l.o.s. angles */
    for (index_mu = 0; index_mu < mu_sizevec[index_z]; index_mu++) {
      class_protect_memcpy(kvec_true[index_z] + index_mu*k_sizevec[index_z], kvec[index_z], k_sizevec[index_z]*sizeof(double));
    }
    class_protect_memcpy(muvec_true[index_z], gl_sym_abscissa, mu_sizevec[index_z]*sizeof(double));
  }

  /** - apply the AP effect */
  eft_apply_ap_effect_in_place(kvec_true,
                               k_sizevec,
                               muvec_true,
                               mu_sizevec,
                               z_size,
                               ap_parallel,
                               ap_perpendicular);

  /** - get the power-spectrum wedges */
  class_call(eft_job_powerspectrum_wedges_ext_growth_rate(peft0,
                                                          peft_size,
                                                          f_z_pk_eft,
                                                          D_z_pk_eft,
                                                          pba,
                                                          pfo,
                                                          ppm,
                                                          ppr,
                                                          pk_out_type,
                                                          zvec,
                                                          f_zvec,
                                                          D_zvec,
                                                          peft_ip,
                                                          z_size,
                                                          kvec_true,
                                                          k_sizevec,
                                                          muvec_true,
                                                          mu_sizevec,
                                                          out_pkmu),
             peft0->error_message, peft0->error_message);

  for (index_z = 0; index_z < z_size; index_z++) {
    free(kvec_true[index_z]); free(muvec_true[index_z]);
  }

  /** - compute the multipoles with fixed Gauss-Lobatto quadrature */
  #pragma omp parallel for schedule(static), collapse(2), shared(z_size, k_sizevec, mu_sizevec, out_pkl, out_pkmu, ap_parallel, ap_perpendicular),   \
                           private(index_z, index_l, index_k, index_mu), firstprivate(gl_sym_weights, lg_measure), default(none)
  for (index_z = 0; index_z < z_size; index_z++) {
    for (index_l = 0; index_l < MULTIPOLE_SIZE; index_l++) {
      for (index_k = 0; index_k < k_sizevec[index_z]; index_k++) {
        out_pkl[index_z][index_l*k_sizevec[index_z] + index_k] = 0.;
        for (index_mu = 0; index_mu < mu_sizevec[index_z]; index_mu++) {
          out_pkl[index_z][index_l*k_sizevec[index_z] + index_k] += gl_sym_weights[index_mu] * lg_measure[index_mu][index_l]    \
                                                                    * out_pkmu[index_z][index_mu*k_sizevec[index_z] + index_k]  \
                                                                    / (ap_parallel[index_z] * ap_perpendicular[index_z] * ap_perpendicular[index_z]);
        }
        out_pkl[index_z][index_l*k_sizevec[index_z] + index_k] *= 0.5*(4*index_l + 1);
      }
    }
  }

  for (index_z = 0; index_z < z_size; index_z++)
    free(out_pkmu[index_z]);

  return _SUCCESS_;
}


/**
 * @brief Computes power-spectrum multipoles at given redshifts and fiducial wavenumbers
 *        using a given growth rate consistently.
 *
 * @param peft0       Input: pointer to the first eft-structure
 * @param peft_size   Input: number of eft-structures
 * @param f_z_pk_eft  Input: logarithmic growth rate at z_pk of each eft-structure
 * @param D_z_pk_eft  Input: growth rate D(z) at z_pk of each eft-structure
 * @param pba         Input: pointer to the background structure
 * @param pfo         Input: pointer to the fourier structure
 * @param ppm         Input: pointer to the primordial structure
 * @param ppr         Input: pointer to the precisions structure
 * @param pk_out_type Input: type of power-spectrum to compute
 * @param zvec        Input: redshifts at which to compute it
 * @param f_zvec      Input: logarithmic growth rate at each z in zvec
 * @param D_zvec      Input: growth rate D(z) at each z in zvec
 * @param peft_ip     Input: Input parameter structures for each eft-structure
 * @param z_size      Input: size of zvec
 * @param kvec        Input: wavenumbers [in 1/Mpc]; indexed as kvec[index_z][index_k]
 * @param k_sizevec   Input: size of the wavenumber array for each mu at index_z; indexed as k_sizevec[index_z]
 * @param ap_parallel Input: parallel AP-effect ratio at each z; defined as H^fid(z)/H^true(z)
 * @param ap_perpendicular  Input: perpendicular AP-effect ratio at each z; defined as D_A^true(z)/D_A^fid(z)
 * @param out_pkl     Output: power-spectrum multipoles of order 0,2,4,...,2*(MULTIPOLE_SIZE-1) at zvec,kvec; indexed as out_pkl[index_z][(l/2)*k_sizevec[index_z] + index_k];
 *                            needs to be pre-allocated with size k_sizevec[index_z]*MULTIPOLE_SIZE at each index_z
 *
 * @return the error status
 */
int eft_job_powerspectrum_multipoles(
                                struct eft * peft0,
                                const int peft_size,
                                struct background * pba,
                                struct fourier * pfo,
                                struct primordial * ppm,
                                struct precision * ppr,
                                enum eft_pk_out_type pk_out_type,
                                const double * const zvec,
                                const struct eft_input_parameters * peft_ip,
                                const int z_size,
                                double ** kvec,
                                const int * const k_sizevec,
                                const double * ap_parallel,
                                const double * ap_perpendicular,
                                double ** out_pkl
                                ) {

  int index_z, index_l, index_mu, index_k, mu_sizevec[z_size];
  double ap_ratio, * kvec_true[z_size], * muvec_true[z_size], * out_pkmu[z_size];

  /** - allocate additional input/output arrays for power-spectrum wedges */
  for (index_z = 0; index_z < z_size; index_z++) {
    mu_sizevec[index_z] = gl_rule_sym_size;
    class_alloc(muvec_true[index_z], mu_sizevec[index_z]*sizeof(double), peft0->error_message);
    class_alloc(kvec_true[index_z], k_sizevec[index_z]*mu_sizevec[index_z]*sizeof(double), peft0->error_message);
    class_alloc(out_pkmu[index_z], k_sizevec[index_z]*mu_sizevec[index_z]*sizeof(double), peft0->error_message);

    /** - copy the fiducial wavenumbers and l.o.s. angles */
    for (index_mu = 0; index_mu < mu_sizevec[index_z]; index_mu++) {
      class_protect_memcpy(kvec_true[index_z] + index_mu*k_sizevec[index_z], kvec[index_z], k_sizevec[index_z]*sizeof(double));
    }
    class_protect_memcpy(muvec_true[index_z], gl_sym_abscissa, mu_sizevec[index_z]*sizeof(double));
  }

  /** - apply the AP effect */
  eft_apply_ap_effect_in_place(kvec_true,
                               k_sizevec,
                               muvec_true,
                               mu_sizevec,
                               z_size,
                               ap_parallel,
                               ap_perpendicular);

  /** - get the power-spectrum wedges */
  class_call(eft_job_powerspectrum_wedges(peft0,
                                          peft_size,
                                          pba,
                                          pfo,
                                          ppm,
                                          ppr,
                                          pk_out_type,
                                          zvec,
                                          peft_ip,
                                          z_size,
                                          kvec_true,
                                          k_sizevec,
                                          muvec_true,
                                          mu_sizevec,
                                          out_pkmu),
             peft0->error_message, peft0->error_message);

  for (index_z = 0; index_z < z_size; index_z++) {
    free(kvec_true[index_z]); free(muvec_true[index_z]);
  }

  /** - compute the multipoles with fixed Gauss-Lobatto quadrature */
  #pragma omp parallel for schedule(static), collapse(2), shared(z_size, k_sizevec, mu_sizevec, out_pkl, out_pkmu, ap_parallel, ap_perpendicular),   \
                           private(index_z, index_l, index_k, index_mu), firstprivate(gl_sym_weights, lg_measure), default(none)
  for (index_z = 0; index_z < z_size; index_z++) {
    for (index_l = 0; index_l < MULTIPOLE_SIZE; index_l++) {
      for (index_k = 0; index_k < k_sizevec[index_z]; index_k++) {
        out_pkl[index_z][index_l*k_sizevec[index_z] + index_k] = 0.;
        for (index_mu = 0; index_mu < mu_sizevec[index_z]; index_mu++) {
          out_pkl[index_z][index_l*k_sizevec[index_z] + index_k] += gl_sym_weights[index_mu] * lg_measure[index_mu][index_l]    \
                                                                    * out_pkmu[index_z][index_mu*k_sizevec[index_z] + index_k]  \
                                                                    / (ap_parallel[index_z] * ap_perpendicular[index_z] * ap_perpendicular[index_z]);
        }
        out_pkl[index_z][index_l*k_sizevec[index_z] + index_k] *= 0.5*(4*index_l + 1);
      }
    }
  }

  for (index_z = 0; index_z < z_size; index_z++)
    free(out_pkmu[index_z]);

  return _SUCCESS_;
}
