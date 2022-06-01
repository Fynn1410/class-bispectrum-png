/** @file FFTLog_galaxy_real.c Documented FFT-Log based 1loop integrals of galaxy/halo power spectrum in perturbation theory 
 * 
 * Azadeh Moradinezhad Dizgah, June 16th 2021
 *
 * This module performs fast computation of the integrals appearing in the expression of 1loop galaxy/halo power sprtcurm.
 * The computation can be performed either in real or redshift-space. IR-resummation and EFT counter terms are included. 
 * The integrals are computed using FFTLog techniques.
 *
 * The algorithm closely follows Ref. arXiv:1708.08130 by Simonovic et al. After computing the FFT coefficents of matter power spectrum, 
 * sampled in logarithmic scale, the algorithm involves re-casting the integrals into a form that is analytically calculable (Matrices M_ij, 
 * which can be written in terms of ratios of Gamma functions) and finally performing vactor-matrix-vector or matrix-vector multiplications. 
 *
 * An important feauture of fast computation of loop integrals is that all the cosmology-dependance of the loop integrals is captured by the FFT coeffcients
 * which have ~NlogN complexity. The matrices involving the Gamma functions, are computed only once, and at each of MCMC varying cosmological parmaters, 
 * these coeffcients are evaluated for all k-values at once and then vector-matrix-vector multiplication are computed.
 *
 * The FFT coeffcients are computed using FFTW package, while the vector-matrix-vector computations are performed using blas library implemented in gsl.
 * The analytic formulas for M_ij matrices are computed in Mathematica using a modified version of the publicaly available notebook by Simpnovic. 
 * 
 * The choice of number of points for FFT decomposition of the power spectrum is very important, in terms of accuracy and execusion time. 
 * For halo loops, if choosing 512 points *which takes ~ 0.45 seconds, I do get subpercent descrepencies between direct numeric integration of the loops.
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
 * Compute the non-propagator type loop contribution to non-linear galaxy power spectrum given the FFTLog coefficents and frequencies (Eq. 2.38, 2.39, 2.41, 2.42, 2.43 of Simonovic 2017)
 * 
 * @param fft_input    Input: structure containing fft coefficents and params
 * @param k             Input: wavenumber in unit of h/Mpc. 
 * @param z             Input: redshift
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the pipeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_0_FFTLog(struct fft_struct *fft_input, double k, double *loops)
{
      int Nmax     = fft_input -> nfft;

      // Linear Power Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      // FFTLog matrices (non-propagator)
      double complex **P22_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(P22, fft_input, k, P22_mat);
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
      double complex *P13_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(M13, fft_input, k, M13_mat);
      // ISPsi is missing
      
}

double P22_new(struct fft_struct *fft_input, double k, double z, int cleanup)
{
      double complex result;

      int Nmax     = fft_input -> nfft;
      double complex **M22_mat = make_2D_c_array(Nmax+1, Nmax+1);
      double complex *vec1 = make_1D_c_array(Nmax+1);

      vec_fill(fft_input, k, vec1);
      np_mat_fill(M22, fft_input, k, M22_mat);

      c_nonprop(vec1, M22_mat, vec1, Nmax+1, &result);

      return creal(result) * pow(k,3.);
}

double P13_new(struct fft_struct *fft_input, double k, double z, int cleanup)
{
      double complex result;

      int Nmax     = fft_input -> nfft;
      double complex *M13_mat = make_1D_c_array(Nmax+1);
      double complex *vec1 = make_1D_c_array(Nmax+1);

      vec_fill(fft_input, k, vec1);
      p_mat_fill(M13, fft_input, k, M13_mat);

      c_dot(M13_mat, vec1, Nmax+1, &result);

      return creal(result) * pow(k,3.);
}


/*
      _________________ 0-th Moment _______________________
 */

double complex IF2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3. + 2.*n1 + 2.*n2)*(-4. + 7.*n1 + 7.*n2));
      double complex denominator = (14.*n1*n2);
      double complex out         = numerator/denominator * Ifunc(n1, n2);

      return out;
}

double complex IF2S2(double complex n1, double complex n2)
{
      double complex numerator   = ((3. - 2.*n1 - 2.*n2)*(-1. + 2.*n1 + 2.*n2)*(6. + 7.*n1 + 7.*n2));
      double complex denominator = (14.*n1*(1 + n1)*n2*(1 + n2);
      double complex out         = numerator/denominator * Ifunc(n1, n2);

      return out;
}

double complex IPP(double complex n1, double complex n2)
{
      double complex out         = 2. * Ifunc(n1, n2);
      return out;
}

double complex IS2(double complex n1, double complex n2)
{
      double complex numerator   = (3. - 2.*n1 - 2.*n2);
      double complex denominator = (n1*n2);
      double complex out         = numerator/denominator * Ifunc(n1, n2);

      return out;
}

double complex IS2S2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3. + 2.*n1 + 2.*n2)*(-1. + 2.*n1 + 2.*n2));
      double complex denominator = (n1*(1. + n1)*n2*(1. + n2));
      double complex out         = numerator/denominator * Ifunc(n1, n2);

      return out;
}

/*
      _________________ 1-st Moment _______________________
 */

double complex IF2G2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3. + 2.*n1 + 2.*n2)*(-1. + 2.*n1 + 2.*n2)*(46. + 98.*Power(n1,3.)*n2 + (13. - 63.*n2)*n2 + 
       7.*Power(n1,2.)*(-9. + 2.*n2*(-5. + 14.*n2)) + n1*(13. + 2.*n2*(-69. + 7.*n2*(-5. + 7.*n2)))));
      double complex denominator = (196.*n1*(1. + n1)*(-1. + 2.*n1)*n2*(1. + n2)*(-1. + 2.*n2));
      double complex out         = numerator/denominator * Ifunc(n1, n2);

      return out;
} 

double complex IF3G3(double complex n1)
{
      double complex numerator   = ((-7. + 9.*n1)*Tan(n1*Pi));
      double complex denominator = (112.*n1*(-6. + 5.*n1 + 5.*Power(n1,2.) - 5.*Power(n1,3.) + Power(n1,4.))*Pi);
      double complex out         = numerator/denominator;

      return out;
} 

double complex IG2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3. + 2.*n1 + 2.*n2)*(-8. + 7.*n1 + 7.*n2));
      double complex denominator = (196.*n1*(1. + n1)*(-1. + 2.*n1)*n2*(1. + n2)*(-1. + 2.*n2));
      double complex out         = numerator/denominator * Ifunc(n1, n2);

      return out;
} 

double complex IG2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3. + 2.*n1 + 2.*n2)*(-8. + 7.*n1 + 7.*n2));
      double complex denominator = (196.*n1*(1. + n1)*(-1. + 2.*n1)*n2*(1. + n2)*(-1. + 2.*n2));
      double complex out         = numerator/denominator * Ifunc(n1, n2);

      return out;
} 

double complex IS2G2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3. + 2.*n1 + 2.*n2)*(-8. + 7.*n1 + 7.*n2));
      double complex denominator = (196.*n1*(1. + n1)*(-1. + 2.*n1)*n2*(1. + n2)*(-1. + 2.*n2));
      double complex out         = numerator/denominator * Ifunc(n1, n2);

      return out;
} 