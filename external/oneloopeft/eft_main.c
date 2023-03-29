#include "header.h"

int eft_spline_sample_points_nonuniform(
                    const double k_min,
                    const double k_feature,
                    const double k_max,
                    const double rel_amplitude,
                    const double width,
                    const int num_points,
                    double * ln_k_sample) {

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
  
  return _SUCCESS_;
}


int eft_indices(struct eft * peft) {

  int i = 0;
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
  
  if (peft->has_rsd) {
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

  if (peft->has_rsd) {
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


  /** - allocate TODO */
  class_alloc(peft->spectra_contributions, peft->index_num*sizeof(double *), peft->error_message);
  class_alloc(peft->loop_matrices, peft->index_num*sizeof(double complex *), peft->error_message);
  class_alloc(peft->loop_matrices_size, peft->index_num*sizeof(int), peft->error_message);


  for (i = 0; i < peft->index_num; i++) {
    class_alloc(peft->spectra_contributions[i], peft->k_size, peft->error_message);

    switch (peft->symmetry[i]) {
    case vec:
      peft->loop_matrices_size[i] = peft->fourier_coeff_size; break;
    case mat_none:
      peft->loop_matrices_size[i] = peft->fourier_coeff_size * peft->fourier_coeff_size; break;
    case mat_symmetric:
      peft->loop_matrices_size[i] = peft->fourier_coeff_size * (peft->fourier_coeff_size + 1) / 2; break;
    default:
      peft->loop_matrices_size[i] = 0; break;
    }

    class_alloc(peft->loop_matrices[i], peft->loop_matrices_size[i]*sizeof(double complex), peft->error_message);
  }

  /** - allocate arrays for different P_lin types, pk_index_num is the last entry in enum eft_pk_type */
  class_alloc(peft->pk_l, pk_index_num*sizeof(double *), peft->error_message);
  class_alloc(peft->pk_l_biased, pk_index_num*sizeof(double *), peft->error_message);
  class_alloc(peft->ddpk_l_biased, pk_index_num*sizeof(double *), peft->error_message);

  for (i = 0; i < pk_index_num; i++) {
    class_alloc(peft->pk_l[i], peft->k_size*sizeof(double), peft->error_message);
    class_alloc(peft->pk_l_biased[i], peft->k_size_fourier*sizeof(double), peft->error_message);
    class_alloc(peft->ddpk_l_biased[i], peft->k_size_fourier*sizeof(double), peft->error_message);
  }

  class_alloc(peft->ln_k_matter_fourier, peft->k_size_fourier*sizeof(double), peft->error_message);
  class_alloc(peft->ln_k_halo_fourier, peft->k_size_fourier*sizeof(double), peft->error_message);

  /** - allocate arrays for the Foruier series coefficients and frequencies */
  class_alloc(peft->fourier_frequencies, peft->fourier_coeff_size*sizeof(double), peft->error_message);
  class_alloc(peft->fourier_coeff, pk_index_num*sizeof(double complex *), peft->error_message);

  for (i = 0; i < pk_index_num; i++) {
    class_alloc(peft->fourier_coeff[i], peft->fourier_coeff_size*sizeof(double complex), peft->error_message);
  }

  return _SUCCESS_;
}

int eft_init(struct background * pba,
            struct fourier * pfo,
            struct primordial * ppm,
            struct eft * peft,
            double z,
            short compute_rsd_spectrum,
            double k_min,
            double fourier_period,
            int num_positive_coefficients,
            int num_sample_points) {

  int it;
  double *ln_kvec;
  double *out_pk_at_z;

  peft->z0 = z;
  peft->has_rsd = compute_rsd_spectrum;
  peft->fourier_coeff_size = 2*num_positive_coefficients - 1;
  peft->k_size_fourier = num_sample_points + 1;
  peft->period = fourier_period;
  peft->bias_halo = -1.55;
  peft->bias_matter = -0.3;
  
  
  class_call(eft_indices(peft), peft->error_message, peft->error_message);

  /* Generate the sampling points in k; the first BAO feature has a wavenumber of 0.07 h/Mpc */
  eft_spline_sample_points_nonuniform(k_min, 0.07*pba->h, 5.e5, 7.5, 5., num_sample_points, peft->ln_k_halo_fourier);

  /* Retrieve the power spectra */
  class_call(fourier_pk_at_kvec_and_z(pba, ppm, pfo, linear, pk_linear,
                                      peft->ln_k_halo_fourier,
                                      num_sample_points,
                                      peft->z0,
                                      pfo->index_pk_cluster,
                                      peft->pk_l_biased[pk_lin]),
              pfo->error_message,
              peft->error_message);

  /** Apply the bias to the power spectra */
  peft->ln_k_halo_fourier[peft->k_size_fourier-1] = log(k_min) + peft->period;
  for (it = 0; it < peft->k_size_fourier; it++) {
    peft->pk_l_biased[pk_lin][it] *= exp(-peft->bias_halo * peft->ln_k_halo_fourier[it]);
  }
  /** Add the last point to make the spectra periodic: P_bias(ln(k_min) + period) = P_bias(ln(k_min)) */
  peft->pk_l_biased[pk_lin][peft->k_size_fourier-1] = peft->pk_l_biased[pk_lin][0];

  /** Spline the biased power spectra periodically */
  class_call(array_spline_table_lines(peft->ln_k_halo_fourier,
                                      peft->k_size_fourier,
                                      peft->pk_l_biased[pk_lin],
                                      1,
                                      peft->ddpk_l_biased[pk_lin],
                                      _SPLINE_PERIODIC_,
                                      peft->error_message),
              peft->error_message,
              peft->error_message);

  /* Compute the Fourier coefficients of the power spectra */
  eft_get_fourier_frequencies(peft->period, peft->fourier_coeff_size, peft->fourier_frequencies);

  /** The DC component has to be treated differently */
  class_call(array_integrate_all_spline_table_lines(
                      peft->ln_k_halo_fourier,
                      peft->k_size_fourier,
                      peft->pk_l_biased[pk_lin],
                      1,
                      peft->ddpk_l_biased[pk_lin],
                      peft->fourier_coeff[pk_lin],
                      peft->error_message),
              peft->error_message,
              peft->error_message);
  peft->fourier_coeff[pk_lin][0] /= peft->period;

  for (it = 1; it < num_positive_coefficients; it++) {
    array_integrate_all_spline_table_lines_fourier(
                      peft->ln_k_halo_fourier,
                      peft->k_size_fourier,
                      peft->pk_l_biased[pk_lin],
                      1,
                      peft->ddpk_l_biased[pk_lin],
                      peft->fourier_frequencies[it],
                      peft->fourier_coeff[pk_lin] + it,
                      peft->error_message);
    
    peft->fourier_coeff[pk_lin][it] /= peft->period;
  }

  for (it = 0; it < num_positive_coefficients-1; it++) {
    peft->fourier_coeff[pk_lin][num_positive_coefficients + it] = conj( peft->fourier_coeff[pk_lin][num_positive_coefficients - 1 - it] );
  }

  if (pfo->fourier_verbose > 2) {
    FILE *ffourier = fopen("output/biased_pk_samples.dat", "w");

    fprintf(ffourier, "# Biased samples of the linear power spectrum at z=%.3f \n", peft->z0);
    fprintf(ffourier, "# for k=%.4e to %.4e \n", exp(peft->ln_k_halo_fourier[0]), exp(peft->ln_k_halo_fourier[peft->k_size_fourier-1]));
    fprintf(ffourier, "# number of wavenumbers equal to %d \n", peft->k_size_fourier);
    fprintf(ffourier, "#    1:k (1/Mpc)            2:P_bias (Mpc)^3        3:d^2P_bias/dln(k)^2 (Mpc)^3 \n");
    for (int i = 0; i < peft->k_size_fourier; i++)
      fprintf(ffourier, "  %.16e       %.16e       %+.16e \n", \
              exp(peft->ln_k_halo_fourier[i]), peft->pk_l_biased[pk_lin][i], peft->ddpk_l_biased[pk_lin][i]);

    fclose(ffourier);

    ffourier = fopen("output/fourier_coefficients.dat", "w");

    fprintf(ffourier, "# Fourier coefficients of the linear power spectrum at z=%.3f \n", peft->z0);
    fprintf(ffourier, "# for omega=%.3f to %.3f \n", exp(peft->fourier_frequencies[num_positive_coefficients]), exp(peft->fourier_frequencies[num_positive_coefficients-1]));
    fprintf(ffourier, "# number of frequencies equal to %d \n", peft->fourier_coeff_size);
    fprintf(ffourier, "#    1:omega                  2:c (Mpc)^3 \n");
    for (int i = 0; i < peft->fourier_coeff_size; i++)
      fprintf(ffourier, "  %.16e       %+.16e       %+.16e \n", \
              peft->fourier_frequencies[i], creal(peft->fourier_coeff[pk_lin][i]), cimag(peft->fourier_coeff[pk_lin][i]));

    fclose(ffourier);
  }

  return _SUCCESS_;
}

int eft_get_fourier_frequencies(
                  const double period,
                  const int num_frequencies,
                  double * frequencies) {

  int m;
  for (m = 0; m < (int)ceil(num_frequencies/2.); m++) {
    frequencies[m] = 2.*_PI_ * m/period;
  }
  for (m = -(int)floor(num_frequencies/2.); m < 0; m++) {
    frequencies[m + num_frequencies] = 2.*_PI_ * m/period;
  }
  return _SUCCESS_;
}
