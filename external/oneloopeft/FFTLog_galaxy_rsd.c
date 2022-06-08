/** @file FFTLog_galaxy_real.c Documented FFT-Log based 1loop integrals of galaxy/halo cpow spectrum in perturbation theory 
 * 
 * Azadeh Moradinezhad Dizgah, June 16th 2021
 *
 * This module performs fast computation of the integrals appearing in the expression of 1loop galaxy/halo cpow sprtcurm.
 * The computation can be performed either in real or redshift-space. IR-resummation and EFT counter terms are included. 
 * The integrals are computed uccsing FFTLog techniques.
 *
 * The algorithm closely follows Ref. arXiv:1708.08130 by Simonovic et al. After computing the FFT coefficents of matter cpow spectrum, 
 * sampled in logarithmic scale, the algorithm involves re-casting the integrals into a form that is analytically calculable (Matrices M_ij, 
 * which can be written in terms of ratios of Gamma functions) and finally performing vactor-matrix-vector or matrix-vector multiplications. 
 *
 * An imporcctant feauture of fast computation of loop integrals is that all the cosmology-dependance of the loop integrals is captured by the FFT coeffcients
 * which have ~NlogN complexity. The matrices involving the Gamma functions, are computed only once, and at each of MCMC varying cosmological parmaters, 
 * these coeffcients are evaluated for all k-values at once and then vector-matrix-vector multiplication are computed.
 *
 * The FFT coeffcients are computed uccsing FFTW package, while the vector-matrix-vector computations are performed uccsing blas library implemented in gsl.
 * The analytic formulas for M_ij matrices are computed in Mathematica uccsing a modified version of the publicaly available notebook by Simpnovic. 
 * 
 * The choice of number of points for FFT decomposition of the cpow spectrum is very imporcctant, in terms of accuracy and execusion time. 
 * For halo loops, if chooccsing 512 points *which takes ~ 0.45 seconds, I do get subpercent descrepencies between direct numeric integration of the loops.
 * The best value seem to be around 600 points, which takes ~ 0.6 seconds and is in exquisit agreement with direct integration. 
 *
 * In summary, the following functions can be called from other modules:
 * -# pgloops_nonpropag()
 * -# pgloops_propag()
 * -# M_cIF2()
 * -# M_cIF2S2()
 * -# M_cIdelta2delta2()
 * -# M_cIcG2cG2()
 * -# M_cIdelta2cG2()
 * -# M_cFcG2()
 */

#include "header.h"

/*
      _________________ 0-th Moment _______________________
 */

/**
 * Compute the non-propagator type loop contribution to non-linear galaxy cpow spectrum given the FFTLog coefficents and frequencies (Eq. 2.38, 2.39, 2.41, 2.42, 2.43 of Simonovic 2017)
 * 
 * @param fft_input    Input: structure containing fft coefficents and params
 * @param k             Input: wavenumber in unit of h/Mpc. 
 * @param z             Input: redshift
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_0_FFTLog(struct fft_struct *fft_input, double k, double *loops)
{
      int Nmax     = fft_input -> nfft;

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      double complex *vec_min = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, fft_input->kmin_fft, vec);

      // FFTLog matrices (non-propagator)
      double complex **M22_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(M22, fft_input, k, M22_mat);
      double complex **IF2_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IF2, fft_input, k, IF2_mat);
      double complex **IF2S2_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IF2S2, fft_input, k, IF2S2_mat);
      double complex **IPP_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IPP, fft_input, k, IPP_mat);
      double complex **IS2_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IS2, fft_input, k, IS2_mat);
      double complex **IS2S2_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IS2S2, fft_input, k, IS2S2_mat);

      // FFTLog matrices (propagator)
      double complex *M13_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(M13, fft_input, k, M13_mat);
      double complex *I_SPsi_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(I_SPsi, fft_input, k, M13_mat);

      // non-propagator calculations
      c_nonprop(vec, M22_mat,   vec, Nmax+1, &loops[0]);
      c_nonprop(vec, IF2_mat,   vec, Nmax+1, &loops[2]);
      c_nonprop(vec, IF2S2_mat, vec, Nmax+1, &loops[3]);
      c_nonprop(vec, IPP_mat,   vec, Nmax+1, &loops[4]);
      c_nonprop(vec, IS2_mat,   vec, Nmax+1, &loops[5]);
      c_nonprop(vec, IS2S2_mat, vec, Nmax+1, &loops[6]);

      // d
      // loops[4] -= 

      // propagator calculations
      c_dot(vec, M13_mat, Nmax+1, &loops[1]);
      c_dot(vec, I_SPsi_mat, Nmax+1, &loops[7]);
}

/*
      _________________ 1-st Moment _______________________
 */

/**
 * Compute the non-propagator type loop contribution to non-linear galaxy cpow spectrum given the FFTLog coefficents and frequencies (Eq. 2.38, 2.39, 2.41, 2.42, 2.43 of Simonovic 2017)
 * 
 * @param fft_input    Input: structure containing fft coefficents and params
 * @param k             Input: wavenumber in unit of h/Mpc. 
 * @param z             Input: redshift
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_1_FFTLog(struct fft_struct *fft_input, double k, double *loops)
{
      int Nmax     = fft_input -> nfft;

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      // FFTLog matrices (non-propagator)
      double complex **IF2G2_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IF2G2, fft_input, k, IF2G2_mat);
      double complex **IG2_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IG2, fft_input, k, IG2_mat);
      double complex **IS2G2_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IS2G2, fft_input, k, IS2G2_mat);
      double complex **IF2p_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IF2p, fft_input, k, IF2p_mat);
      double complex **IPPp_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IPPp, fft_input, k, IPPp_mat);
      double complex **IS2p_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IS2p, fft_input, k, IS2p_mat);

      // FFTLog matrices (propagator)
      double complex *IF3G3_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(IF3G3, fft_input, k, IF3G3_mat);
      double complex *IF2p2_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(IF2p2, fft_input, k, IF2p2_mat);
      double complex *IG2p_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(IG2p, fft_input, k, IG2p_mat);
      
      // non-propagator calculations
      c_nonprop(vec, IF2G2_mat, vec, Nmax+1, &loops[0]);
      c_nonprop(vec, IG2_mat,   vec, Nmax+1, &loops[2]);
      c_nonprop(vec, IS2G2_mat, vec, Nmax+1, &loops[3]);
      c_nonprop(vec, IF2p_mat,  vec, Nmax+1, &loops[5]);
      c_nonprop(vec, IPPp_mat,  vec, Nmax+1, &loops[8]);
      c_nonprop(vec, IS2p_mat,  vec, Nmax+1, &loops[9]);

      // propagator calculations
      c_dot(vec, IF3G3_mat, Nmax+1, &loops[1]);
      c_dot(vec, IF2p2_mat, Nmax+1, &loops[6]);
      c_dot(vec, IG2p_mat,  Nmax+1, &loops[7]);
      loops[4] = 0. ;// psi
}

/*
      _________________ 2-nd Moment _______________________
 */

/**
 * Compute the non-propagator type loop contribution to non-linear galaxy cpow spectrum given the FFTLog coefficents and frequencies (Eq. 2.38, 2.39, 2.41, 2.42, 2.43 of Simonovic 2017)
 * 
 * @param fft_input    Input: structure containing fft coefficents and params
 * @param k             Input: wavenumber in unit of h/Mpc. 
 * @param z             Input: redshift
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_2_FFTLog(struct fft_struct *fft_input, double k, double mu, double *loops)
{
      int Nmax     = fft_input -> nfft;

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      // FFTLog matrices (non-propagator)
      double complex **IG2G2_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IG2G2, fft_input, k, IG2G2_mat);
      double complex **IG2pm_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IG2pm, fft_input, k, IG2pm_mat);
      double complex **IPPmm2_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IPPmm2, fft_input, k, IPPmm2_mat);

      // FFTLog matrices (non-propagator) mu dependent
      double complex **IF2pm3_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IF2pm3, fft_input, k, mu, IF2pm3_mat);
      double complex **IG2pm2_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IG2pm2, fft_input, k, mu, IG2pm2_mat);
      double complex **IPPpm_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IPPpm, fft_input, k, mu, IPPpm_mat);
      double complex **IS2pm_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IS2pm, fft_input, k, mu, IS2pm_mat);
      double complex **IPPmm1_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IPPmm1, fft_input, k, mu, IPPmm1_mat);
      double complex **IPPmm3_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IPPmm3, fft_input, k, mu, IPPmm3_mat);

      // FFTLog matrices (propagator)
      double complex *IG3_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(IG3, fft_input, k, IG3_mat);
      double complex *IF2pm_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(IF2pm, fft_input, k, IF2pm_mat);

      // FFTLog matrices (propagator) mu dependent
      double complex *IG2pmm_mat = make_1D_c_array(Nmax+1);
      p_mu_mat_fill(IG2pmm, fft_input, k, mu, IG2pmm_mat);

      // non-propagator calculations
      c_nonprop(vec, IG2G2_mat,  vec, Nmax+1, &loops[4]);
      c_nonprop(vec, IG2pm_mat,  vec, Nmax+1, &loops[6]);
      c_nonprop(vec, IPPmm2_mat, vec, Nmax+1, &loops[10]);
      c_nonprop(vec, IF2pm3_mat, vec, Nmax+1, &loops[0]);
      c_nonprop(vec, IG2pm2_mat, vec, Nmax+1, &loops[1]);
      c_nonprop(vec, IPPpm_mat,  vec, Nmax+1, &loops[2]);
      c_nonprop(vec, IS2pm_mat,  vec, Nmax+1, &loops[3]);
      c_nonprop(vec, IPPmm1_mat, vec, Nmax+1, &loops[9]);
      c_nonprop(vec, IPPmm3_mat, vec, Nmax+1, &loops[11]);

      // propagator calculations
      c_dot(vec, IG3_mat,    Nmax+1, &loops[5]);
      c_dot(vec, IF2pm_mat,  Nmax+1, &loops[8]);
      c_dot(vec, IG2pmm_mat, Nmax+1, &loops[7]);
}

/*
      _________________ 3-rd Moment _______________________
 */

/**
 * Compute the non-propagator type loop contribution to non-linear galaxy cpow spectrum given the FFTLog coefficents and frequencies (Eq. 2.38, 2.39, 2.41, 2.42, 2.43 of Simonovic 2017)
 * 
 * @param fft_input    Input: structure containing fft coefficents and params
 * @param k             Input: wavenumber in unit of h/Mpc. 
 * @param z             Input: redshift
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_3_FFTLog(struct fft_struct *fft_input, double k, double mu, double *loops)
{
      int Nmax     = fft_input -> nfft;

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      // FFTLog matrices (non-propagator) mu dependent
      double complex **IG2pm3_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IG2pm3, fft_input, k, mu, IG2pm3_mat);
      double complex **IPPpm31_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IPPpm31, fft_input, k, mu, IPPpm31_mat);
      double complex **IPPpm32_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IPPpm32, fft_input, k, mu, IPPpm32_mat);

      // FFTLog matrices (propagator) mu dependent
      double complex *IG2pmm3_mat = make_1D_c_array(Nmax+1);
      p_mu_mat_fill(IG2pmm3, fft_input, k, mu, IG2pmm3_mat);

      // non-propagator calculations
      c_nonprop(vec, IG2pm3_mat,  vec, Nmax+1, &loops[0]);
      c_nonprop(vec, IPPpm31_mat,  vec, Nmax+1, &loops[2]);
      c_nonprop(vec, IPPpm32_mat, vec, Nmax+1, &loops[3]);

      // propagator calculations
      c_dot(vec, IG2pmm3_mat, Nmax+1, &loops[1]);
}

/*
      _________________ 4-th Moment _______________________
 */

/**
 * Compute the non-propagator type loop contribution to non-linear galaxy cpow spectrum given the FFTLog coefficents and frequencies (Eq. 2.38, 2.39, 2.41, 2.42, 2.43 of Simonovic 2017)
 * 
 * @param fft_input    Input: structure containing fft coefficents and params
 * @param k             Input: wavenumber in unit of h/Mpc. 
 * @param z             Input: redshift
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_4_FFTLog(struct fft_struct *fft_input, double k, double mu, double *loops)
{
      int Nmax     = fft_input -> nfft;

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      // FFTLog matrices (non-propagator) mu dependent
      double complex **IPPm41_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IPPm41, fft_input, k, mu, IPPm41_mat);
      double complex **IPPm42_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IPPm42, fft_input, k, mu, IPPm42_mat);
      double complex **IPPm43_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mu_mat_fill(IPPm43, fft_input, k, mu, IPPm43_mat);

      // non-propagator calculations
      c_nonprop(vec, IPPm41_mat, vec, Nmax+1, &loops[0]);
      c_nonprop(vec, IPPm42_mat, vec, Nmax+1, &loops[1]);
      c_nonprop(vec, IPPm43_mat, vec, Nmax+1, &loops[2]);
}

double P22_new(struct fft_struct *fft_input, double k, double z, int cleanup)
{
      double result;

      int Nmax     = fft_input -> nfft;
      double complex **M22_mat = make_2D_c_array(Nmax+1, Nmax+1);
      double complex *vec1 = make_1D_c_array(Nmax+1);

      vec_fill(fft_input, k, vec1);
      np_mat_fill(M22, fft_input, k, M22_mat);

      c_nonprop(vec1, M22_mat, vec1, Nmax+1, &result);

      return result * pow(k,3.);
}

double P13_new(struct fft_struct *fft_input, double k, double z, int cleanup)
{
      double result;

      int Nmax     = fft_input -> nfft;
      double complex *M13_mat = make_1D_c_array(Nmax+1);
      double complex *vec1 = make_1D_c_array(Nmax+1);

      vec_fill(fft_input, k, vec1);
      p_mat_fill(M13, fft_input, k, M13_mat);

      c_dot(M13_mat, vec1, Nmax+1, &result);

      return result * pow(k,3.);
}


/*
      _________________ 0-th Moment _______________________
 */

double complex IF2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-4 + 7*n1 + 7*n2));
      double complex denominator = (14.*n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex IF2S2(double complex n1, double complex n2)
{
      double complex numerator   = ((3 - 2*n1 - 2*n2)*(-1 + 2*n1 + 2*n2)*(6 + 7*n1 + 7*n2));
      double complex denominator = (28.*n1*(1 + n1)*n2*(1 + n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex IPP(double complex n1, double complex n2)
{
      double complex out         = 2. * J(n1, n2);
      return out;
}

double complex IS2(double complex n1, double complex n2)
{
      double complex numerator   = (3 - 2*n1 - 2*n2);
      double complex denominator = (n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex IS2S2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2));
      double complex denominator = (n1*(1 + n1)*n2*(1 + n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex I_SPsi(double complex nu1)
{
      double complex numerator    = 15. * ctan(nu1 * M_PI);
      double complex denominator  = 28. * M_PI * (nu1 + 1.) * nu1 * (nu1 - 1.) * (nu1 - 2.) * (nu1 - 3.);
      double complex out          = - numerator/denominator;

      return out;

}

/*
      _________________ 1-st Moment _______________________
 */

double complex IF2G2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(46 + 98*cpow(n1,3)*n2 + (13 - 63*n2)*n2 + 
       7*cpow(n1,2)*(-9 + 2*n2*(-5 + 14*n2)) + n1*(13 + 2*n2*(-69 + 7*n2*(-5 + 7*n2)))));
      double complex denominator = (196.*n1*(1 + n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex IF3G3(double complex n1)
{
      double complex numerator   = ((-7 + 9*n1)*ctan(n1*M_PI));
      double complex denominator = (112.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out         = numerator/denominator;

      return out;
} 

double complex IG2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-8 + 7*n1 + 7*n2));
      double complex denominator = (28.*n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex IS2G2(double complex n1, double complex n2)
{
      double complex numerator   = -((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-2 + 7*n1 + 7*n2));
      double complex denominator = (28*n1*(1 + n1)*n2*(1 + n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex IF2p(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-5 + n1*(-4 + 7*n1 + 7*n2)));
      double complex denominator = (14.*n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2);
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex IF2p2(double complex n1)
{
      return 0;
} 

double complex IG2p(double complex n1)
{
      double complex numerator   = (3*(9 - 4*n1)*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6 + 11*n1 - 6*cpow(n1,2) + cpow(n1,3))*M_PI);
      double complex out         = numerator/denominator;

      return out;
} 

double complex IPPp(double complex n1, double complex n2)
{
      double complex numerator   = -(18 - 21*n1 + 6*cpow(n1,2) - 18*n2 + 8*n1*n2 + 4*cpow(n2,2));
      double complex denominator = (6*n1 - 4*cpow(n1,2));
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex IS2p(double complex n1, double complex n2)
{
      double complex numerator   = -0.5*((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2));
      double complex denominator = (n1*(1 + n1)*(-3 + 2*n1)*n2);
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

/*
      _________________ 2-nd Moment _______________________
 */

double complex IF2pm3(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*
      (-((-1 + 2*n1)*(-1 + 2*n2)*(6 + 7*n1 + 7*n2)) + 
        (34 + n1 + n2 + 56*cpow(n1,3)*n2 - 54*cpow(n2,2) + 2*cpow(n1,2)*(-27 - 2*n2 + 56*cpow(n2,2)) + 
           4*n1*n2*(-21 + n2*(-1 + 14*n2)))*cpow(mu,2))*Gamma(-2*n1)*Gamma(-2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*
      csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (14*(1 + n1)*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex IG2pm2(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + n1 + mu)*(-3 + 2*n1 + 2*mu)*(-1 + 2*n1 + 2*mu)*(-3 + n1*(-8 + 7*n1 + 7*mu)));
      double complex denominator = (14.*n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*mu);
      double complex out = numerator/denominator;

      return out;
} 

double complex IPPpm(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + n1*(2 - 4*n2) + 2*n2 + 
       2*(-1 + 2*n1 + 2*n2)*(-10 + 2*cpow(n1,2) + (31 - 16*n2)*n2 + n1*(-5 + 2*n2))*cpow(mu,2))*Gamma(3 - 2*n1)*
     Gamma(2 - 2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (8.*(-1 + n1)*n1*(-1 + 2*n1)*n2*(-1 + 2*n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex IS2pm(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-1 + (1 + n1 + n2)*cpow(mu,2))*Gamma(2 - 2*n1)*Gamma(2 - 2*n2)*
     Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (2.*n1*(1 + n1)*n2*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex IG2G2(double complex n1, double complex n2)
{
      double complex numerator   =((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(50 + 98*cpow(n1,3)*n2 - n2*(9 + 35*n2) + 
       7*cpow(n1,2)*(-5 + 2*n2*(-9 + 14*n2)) + n1*(-9 + 2*n2*(-33 + 7*n2*(-9 + 7*n2)))));
      double complex denominator = (196.*n1*(1 + n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2));
      double complex out = numerator/denominator * J(n1,n2);

      return out;
} 

double complex IG3(double complex n1)
{
      double complex numerator   = (3*(-5 + 3*n1)*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex IG2pm(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-3 + n1*(-8 + 7*n1 + 7*n2)));
      double complex denominator = (14.*n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2);
      double complex out = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex IG2pmm(double complex n1, double mu)
{
      double complex numerator   = (-9*(5 - 9*n1 + 2*cpow(n1,2))*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex IF2pm(double complex n1)
{
      return 0.;
} 

double complex IPPmm1(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -4 - (8*(-1 + n1)*(-4 + n1 + n2)*(-2 + n1 + 8*(-3 + n1 + n2)*(-5 + 2*n2)*cpow(mu,2)));
      double complex denominator = ((-9 + 2*n1 + 2*n2)*(-7 + 2*n1 + 2*n2)*(-3 + 2*n1 + 8*(-5 + 2*n2)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)))*ctan(n1*M_PI)*ctan((n1 + n2)*M_PI);
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

double complex IPPmm2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(n1*(-3 + 2*n1) + n2*(-1 + 2*n2)));
      double complex denominator = (n1*(-3 + 2*n1)*n2*(-1 + 2*n2));
      double complex out = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex IPPmm3(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -(n1*(-3 + 2*n1)*(-4 + n1 + n2)*(-1 + 2*n2)*(-5 + 2*n1 + 2*n2) + 
       n2*(-4 + n1 + n2)*(-3 + 2*n2)*(-1 + 2*n2)*(-5 + 2*n1 + 2*n2) + 
       n1*(-5 + 2*n1)*(-3 + 2*n1)*(-4 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*cpow(mu,2) + 
       (-3 + 2*n1)*n2*(-4 + n1 + n2)*(-1 + 2*n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*cpow(mu,2) - 
       n1*(-3 + 2*n1)*n2*(-1 + 2*n2)*(-3 + 2*n2 + (-5 + 2*n1)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator =  ((3 - 2*n1)*n1*(1 - 2*n2)*n2*(3 - 2*n2 - (-5 + 2*n1)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

/*
      _________________ 3-rd Moment _______________________
 */

double complex IG2pm3(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-((-1 + 2*n1)*(-1 + 2*n2)*(-2 + 7*n1 + 7*n2)) + 
       (26 + 56*cpow(n1,3)*n2 + (9 - 38*n2)*n2 + 2*cpow(n1,2)*(-19 + 2*n2*(-9 + 28*n2)) + 
          n1*(9 + 4*n2*(-21 + n2*(-9 + 14*n2))))*cpow(mu,2))*Gamma(-2*n1)*Gamma(-2*n2)*Gamma(2*(-2 + n1 + n2))*
     csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (14.*(1 + n1)*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex IG2pmm3(double complex n1, double mu)
{
      double complex numerator   = (-3*(3 + (-69 - 7*n1 + 2*cpow(n1,2))*cpow(mu,2))*ctan(n1*M_PI));
      double complex denominator = (56.*(-3 + n1)*n1*(2 - n1 - 2*cpow(n1,2) + cpow(n1,3))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 


double complex IPPpm31(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (2*(-4 + n1 + n2)*(-3 + n1 + n2)*(-2 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*
     (-1 + 2*n2 + (-1 + 2*n1)*(1 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = (n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2*(-1 + 2*n2)*(-3 + 2*n2 + (-5 + 2*n1)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

double complex IPPpm32(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (-2*(-5 + n1 + n2)*(-4 + n1 + n2)*(-3 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*
     (-3 + 6*n2 + (-3 + 2*n1)*(1 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = n1*(1 + n1)*(-5 + 2*n1)*(-3 + 2*n1)*n2*(-1 + 2*n2)*(-9 + 6*n2 + (-7 + 2*n1)*(-5 + 2*n1 + 2*n2)*cpow(mu,2));
      double complex out = numerator/denominator * M3(n1, n2, mu);

      return out;
}

/*
      _________________ 4-th Moment _______________________
 */

double complex IPPm41(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (2*(-4 + n1 + n2)*(-3 + n1 + n2)*(-2 + n1 + n2)*(-1 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*
     (-1 + 2*n1 + 2*n2)*(1 + 2*n1 + 2*n2)*(1 + 2*n1 + 8*(-1 + 2*n2)*(3 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = (n1*(1 + n1)*(-1 + 2*n1)*(1 + 2*n1)*n2*(1 + n2)*(-3 + 2*n2)*(-1 + 2*n2)*
     (-3 + 2*n1 + 8*(-5 + 2*n2)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

double complex IPPm42(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (2*(-6 + n1 + n2)*(-5 + n1 + n2)*(-4 + n1 + n2)*(-3 + n1 + n2)*(-7 + 2*n1 + 2*n2)*(-5 + 2*n1 + 2*n2)*
     (-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-3 + 12*cpow(n2,2) + 
       6*(-3 + 2*n1)*(1 + 2*n2)*(1 + 2*n1 + 2*n2)*cpow(mu,2) + 
       (-5 + 2*n1)*(-3 + 2*n1)*(1 + 2*n1 + 2*n2)*(3 + 2*n1 + 2*n2)*cpow(mu,4)));
      double complex denominator = (n1*(1 + n1)*(-5 + 2*n1)*(-3 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2)*(1 + 2*n2)*
     (3*(-5 + 2*n2)*(-3 + 2*n2) + 6*(-7 + 2*n1)*(-3 + 2*n2)*(-7 + 2*n1 + 2*n2)*cpow(mu,2) + 
       (-9 + 2*n1)*(-7 + 2*n1)*(-7 + 2*n1 + 2*n2)*(-5 + 2*n1 + 2*n2)*cpow(mu,4)));
      double complex out = numerator/denominator * M4(n1, n2, mu);

      return out;
} 

double complex IPPm43(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (-4*(-5 + n1 + n2)*(-4 + n1 + n2)*(-3 + n1 + n2)*(-2 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*
     (-1 + 2*n1 + 2*n2)*(1 + 2*n1 + 2*n2)*(3 + 6*n2 + (-3 + 2*n1)*(3 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = (n1*(1 + n1)*(-5 + 2*n1)*(-3 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2)*(1 + 2*n2)*
     (-9 + 6*n2 + (-7 + 2*n1)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M3(n1, n2, mu);

      return out;
} 