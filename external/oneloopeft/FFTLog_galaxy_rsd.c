/** @file FFTLog_galaxy_real.c Documented FFT-Log based 1loop integrals of galaxy/halo cpow spectrum in perturbation theory 
 * 
 * Azadeh Moradinezhad Dizgah, June 16th 2021
 *
 * This module performs fast computation of the integrals appearing in the expression of 1loop galaxy/halo cpow sprtcurm.
 * The computation can be performed either in real or redshift-space. IR-resummation and EFT counter terms are included. 
 * The integrals are computed ucccsing FFTLog techniques.
 *
 * The algorithm closely follows Ref. arXiv:1708.08130 by Simonovic et al. After computing the FFT coefficents of matter cpow spectrum, 
 * sampled in logarithmic scale, the algorithm involves re-casting the integrals into a form that is analytically calculable (Matrices M_ij, 
 * which can be written in terms of ratios of Gamma functions) and finally performing vactor-matrix-vector or matrix-vector multiplications. 
 *
 * An imporccctant feauture of fast computation of loop integrals is that all the cosmology-dependance of the loop integrals is captured by the FFT coeffcients
 * which have ~NlogN complexity. The matrices involving the Gamma functions, are computed only once, and at each of MCMC varying cosmological parmaters, 
 * these coeffcients are evaluated for all k-values at once and then vector-matrix-vector multiplication are computed.
 *
 * The FFT coeffcients are computed ucccsing FFTW package, while the vector-matrix-vector computations are performed ucccsing blas library implemented in gsl.
 * The analytic formulas for M_ij matrices are computed in Mathematica ucccsing a modified version of the publicaly available notebook by Simpnovic. 
 * 
 * The choice of number of points for FFT decomposition of the cpow spectrum is very imporccctant, in terms of accuracy and execusion time. 
 * For halo loops, if choocccsing 512 points *which takes ~ 0.45 seconds, I do get subpercent descrepencies between direct numeric integration of the loops.
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
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_M_M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_0_FFTLog(struct fft_struct *fft_input, double k, double *np_loops, double *p_loops)
{
      int Nmax     = fft_input -> nfft;

      double *np = make_1Darray(6);
      double *p  = make_1Darray(2);

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      double complex *vec_min = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, fft_input->kmin_fft, vec_min);

      // FFTLog matrices non-propagator
      double complex **I2200_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(I2200, fft_input, k, 0., I2200_mat);
      double complex **Idelta200_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(Idelta200, fft_input, k, 0., Idelta200_mat);
      double complex **IG200_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IG200, fft_input, k, 0., IG200_mat);
      double complex **Idelta2delta200_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(Idelta2delta200, fft_input, k, 0., Idelta2delta200_mat);
      double complex **IG2G200_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IG2G200, fft_input, k, 0., IG2G200_mat);
      double complex **Idelta2G200_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(Idelta2G200, fft_input, k, 0., Idelta2G200_mat);

      // FFTLog matrices propagator
      double complex *I1300_mat   = make_1D_c_array(Nmax+1);
      p_mat_fill(I1300, fft_input, k, 0., I1300_mat);
      double complex *FG200_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(FG200, fft_input, k, 0., FG200_mat);

      // non-propagator calculations
      c_nonprop(vec, I2200_mat,           vec, Nmax+1, &np[0]);
      c_nonprop(vec, Idelta200_mat,       vec, Nmax+1, &np[1]);
      c_nonprop(vec, IG200_mat,           vec, Nmax+1, &np[2]);
      c_nonprop(vec, Idelta2delta200_mat, vec, Nmax+1, &np[3]);
      c_nonprop(vec, IG2G200_mat,         vec, Nmax+1, &np[4]);
      c_nonprop(vec, Idelta2G200_mat,     vec, Nmax+1, &np[5]);

      double Idelta2delta200_const;
      c_nonprop(vec_min, Idelta2delta200_mat, vec_min, Nmax+1, &Idelta2delta200_const);

      // propagator calculations
      c_dot(vec, I1300_mat, Nmax+1, &p[0]);
      c_dot(vec, FG200_mat, Nmax+1, &p[1]);

      // adding factored out k and mu dependencies
      np_loops[0] = pow(k, 3.) * np[0];
      np_loops[1] = pow(k, 3.) * np[1];
      np_loops[2] = pow(k, 3.) * np[2];
      np_loops[3] = pow(k, 3.) * np[3] - pow(fft_input->kmin_fft, 3.) * Idelta2delta200_const ;
      np_loops[4] = pow(k, 3.) * np[4];
      np_loops[5] = pow(k, 3.) * np[5];

      p_loops[0] = pow(k, 3.) * p[0];
      p_loops[1] = pow(k, 3.) * p[1];
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
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_M_M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_1_FFTLog(struct fft_struct *fft_input, double k, double mu, double *np_loops, double *p_loops)
{
      int Nmax     = fft_input -> nfft;

      double *np = make_1Darray(7);
      double *p  = make_1Darray(3);

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      // FFTLog matrices (non-propagator)
      double complex **I2201_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(I2201, fft_input, k, mu, I2201_mat);
      double complex **Idelta201_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(Idelta201, fft_input, k, mu, Idelta201_mat);
      double complex **IG201_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(IG201, fft_input, k, mu, IG201_mat);
      double complex **FG201_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(FG201, fft_input, k, mu, FG201_mat);
      double complex **J21101_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(J21101, fft_input, k, mu, J21101_mat);
      double complex **Jdelta201_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(Jdelta201, fft_input, k, mu, Jdelta201_mat);
      double complex **JG201_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(JG201, fft_input, k, mu, JG201_mat);

      // FFTLog matrices (propagator)
      double complex *I1301_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(I1301, fft_input, k, mu, I1301_mat);
      double complex *J12101_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(J12101, fft_input, k, mu, J12101_mat);
      double complex *J11201_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(J11201, fft_input, k, mu, J11201_mat);
      
      // non-propagator calculations
      c_nonprop(vec, I2201_mat,     vec, Nmax+1, &np[0]);
      c_nonprop(vec, Idelta201_mat, vec, Nmax+1, &np[1]);
      c_nonprop(vec, IG201_mat,     vec, Nmax+1, &np[2]);
      c_nonprop(vec, FG201_mat,     vec, Nmax+1, &np[3]);
      c_nonprop(vec, J21101_mat,    vec, Nmax+1, &np[4]);
      c_nonprop(vec, Jdelta201_mat, vec, Nmax+1, &np[5]);
      c_nonprop(vec, JG201_mat,     vec, Nmax+1, &np[6]);

      // propagator calculations
      c_dot(vec, I1301_mat,  Nmax+1, &p[0]);
      c_dot(vec, J12101_mat, Nmax+1, &p[1]);
      c_dot(vec, J11201_mat, Nmax+1, &p[2]);

      // adding factored out k and mu dependencies
      np_loops[0] = pow(k, 3.) * np[0];
      np_loops[1] = pow(k, 3.) * np[1];
      np_loops[2] = pow(k, 3.) * np[2];
      np_loops[3] = pow(k, 3.) * np[3];
      np_loops[4] = pow(k, 3.) * np[4] * k * mu;
      np_loops[5] = pow(k, 3.) * np[5] * k * mu;
      np_loops[6] = pow(k, 3.) * np[6] * k * mu;

      p_loops[0] = pow(k, 3.) * p[0];
      p_loops[1] = pow(k, 3.) * p[1] * k * mu;
      p_loops[2] = pow(k, 3.) * p[2] * k * mu;
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
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_M_M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_2_FFTLog(struct fft_struct *fft_input, double k, double mu, double *np_loops, double *p_loops)
{
      int Nmax     = fft_input -> nfft;

      double *np = make_1Darray(8);
      double *p  = make_1Darray(4);

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      // FFTLog matrices (non-propagator)
      double complex **J21102_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(J21102, fft_input, k, mu, J21102_mat);
      double complex **Jdelta202_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(Jdelta202, fft_input, k, mu, Jdelta202_mat);
      double complex **JG202_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(JG202, fft_input, k, mu, JG202_mat);
      double complex **I2211_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(I2211, fft_input, k, mu, I2211_mat);
      double complex **J21111_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(J21111, fft_input, k, mu, J21111_mat);
      double complex **N11a_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(N11a, fft_input, k, mu, N11a_mat);
      double complex **N11b_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(N11b, fft_input, k, mu, N11b_mat);
      double complex **N11c_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(N11c, fft_input, k, mu, N11c_mat);

      // FFTLog matrices (propagator)
      double complex *J12102_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(J12102, fft_input, k, mu, J12102_mat);
      double complex *I1311_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(I1311, fft_input, k, mu, I1311_mat);
      double complex *J12111_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(J12111, fft_input, k, mu, J12111_mat);
      double complex *J11211_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(J11211, fft_input, k, mu, J11211_mat);

      // non-propagator calculations
      c_nonprop(vec, J21102_mat,    vec, Nmax+1, &np[0]);
      c_nonprop(vec, Jdelta202_mat, vec, Nmax+1, &np[1]);
      c_nonprop(vec, JG202_mat,     vec, Nmax+1, &np[2]);
      c_nonprop(vec, I2211_mat,     vec, Nmax+1, &np[3]);
      c_nonprop(vec, J21111_mat,    vec, Nmax+1, &np[4]);
      c_nonprop(vec, N11a_mat,      vec, Nmax+1, &np[5]);
      c_nonprop(vec, N11b_mat,      vec, Nmax+1, &np[6]);
      c_nonprop(vec, N11c_mat,      vec, Nmax+1, &np[7]);

      // propagator calculations
      c_dot(vec, J12102_mat, Nmax+1, &p[0]);
      c_dot(vec, I1311_mat,  Nmax+1, &p[1]);
      c_dot(vec, J12111_mat, Nmax+1, &p[2]);
      c_dot(vec, J11211_mat, Nmax+1, &p[3]);

      // adding factored out k and mu dependencies
      np_loops[0] = pow(k, 3.) * np[0] * pow(k, 2.);
      np_loops[1] = pow(k, 3.) * np[1] * pow(k, 2.);
      np_loops[2] = pow(k, 3.) * np[2] * pow(k, 2.);
      np_loops[3] = pow(k, 3.) * np[3];
      np_loops[4] = pow(k, 3.) * np[4] * k * mu;
      np_loops[5] = pow(k, 3.) * (np[5] * pow(k, 2.) + np[6] * pow(k*mu, 2.) + np[7] * pow(k, 2.));

      p_loops[0] = pow(k, 3.) * p[0] * pow(k, 2.);
      p_loops[1] = pow(k, 3.) * p[1];
      p_loops[2] = pow(k, 3.) * p[2] * k * mu;
      p_loops[3] = pow(k, 3.) * p[3] * k * mu;
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
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_M_M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_3_FFTLog(struct fft_struct *fft_input, double k, double mu, double *np_loops, double *p_loops)
{
      int Nmax     = fft_input -> nfft;

      double *np = make_1Darray(3);
      double *p  = make_1Darray(1);

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      // FFTLog matrices (non-propagator)
      double complex **J21112_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(J21112, fft_input, k, mu, J21112_mat);
      double complex **N12a_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(N12a, fft_input, k, mu, N12a_mat);
      double complex **N12b_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(N12b, fft_input, k, mu, N12b_mat);

      // FFTLog matrices (propagator)
      double complex *J12112_mat = make_1D_c_array(Nmax+1);
      p_mat_fill(J12112, fft_input, k, mu, J12112_mat);

      // non-propagator calculations
      c_nonprop(vec, J21112_mat, vec, Nmax+1, &np[0]);
      c_nonprop(vec, N12a_mat,   vec, Nmax+1, &np[1]);
      c_nonprop(vec, N12b_mat,   vec, Nmax+1, &np[2]);

      // propagator calculations
      c_dot(vec, J12112_mat, Nmax+1, &p[0]);

      // adding factored out k and mu dependencies
      np_loops[0] = pow(k, 3.) * np[0] * pow(k, 2.);
      np_loops[1] = pow(k, 3.) * (np[1] + np[2]) * pow(k, 3.) * mu;

      p_loops[0] = pow(k, 3.) * p[0] * pow(k, 2.);
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
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_M_M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_4_FFTLog(struct fft_struct *fft_input, double k, double mu, double *np_loops)
{
      int Nmax     = fft_input -> nfft;

      double *np = make_1Darray(3);

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

      // FFTLog matrices (non-propagator)
      double complex **N22a_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(N22a, fft_input, k, mu, N22a_mat);
      double complex **N22b_mat   = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(N22b, fft_input, k, mu, N22b_mat);
      double complex **N22c_mat = make_2D_c_array(Nmax+1, Nmax+1);
      np_mat_fill(N22c, fft_input, k, mu, N22c_mat);

      // non-propagator calculations
      c_nonprop(vec, N22a_mat, vec, Nmax+1, &np[0]);
      c_nonprop(vec, N22b_mat, vec, Nmax+1, &np[1]);
      c_nonprop(vec, N22c_mat, vec, Nmax+1, &np[2]);

      // adding factored out k and mu dependencies
      np_loops[0] = pow(k, 3.) * pow(k, 4.) * (np[0] * pow(mu, 2.) +  np[1] + np[2] * pow(mu, 2.));
}

double P22_new(struct fft_struct *fft_input, double k, double z, int cleanup)
{
      double result;

      int Nmax     = fft_input -> nfft;
      double complex **M22_mat = make_2D_c_array(Nmax+1, Nmax+1);
      double complex *vec1 = make_1D_c_array(Nmax+1);

      vec_fill(fft_input, k, vec1);
      np_mat_fill(I2200, fft_input, k, 0., M22_mat);

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
      p_mat_fill(I1300, fft_input, k, 0., M13_mat);

      c_dot(M13_mat, vec1, Nmax+1, &result);

      return result * pow(k,3.);
}


/*
      _________________ 0-th Moment _______________________
 */

double complex I2200(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(58 + 98*cpow(n1,3)*n2 + (3 - 91*n2)*n2 + 
       7*cpow(n1,2)*(-13 - 2*n2 + 28*cpow(n2,2)) + n1*(3 + 2*n2*(-73 + 7*n2*(-1 + 7*n2)))));
      double complex denominator = (392.*n1*(1 + n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex I1300(double complex n1, double mu)
{
      double complex numerator   = ((1 + 9*n1)*ctan(n1*M_PI));
      double complex denominator = (672.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out         = numerator/denominator;

      return out;
}

double complex Idelta200(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-4 + 7*n1 + 7*n2));
      double complex denominator = (28.*n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex IG200(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(6 + 7*n1 + 7*n2));
      double complex denominator = (56*n1*(1 + n1)*n2*(1 + n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex Idelta2delta200(double complex n1, double complex n2, double mu)
{
      double complex out         = J(n1, n2);
      return out;
}

double complex IG2G200(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2));
      double complex denominator = (2.*n1*(1 + n1)*n2*(1 + n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex Idelta2G200(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (3 - 2*n1 - 2*n2);
      double complex denominator = (2.*n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex FG200(double complex n1, double mu)
{
      double complex numerator    = (-15*ctan(n1*M_PI));
      double complex denominator  = (112.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out          = numerator/denominator;

      return out;

}

/*
      _________________ 1-st Moment _______________________
 */

double complex I2201(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(46 + 98*cpow(n1,3)*n2 + (13 - 63*n2)*n2 + 
       7*cpow(n1,2)*(-9 + 2*n2*(-5 + 14*n2)) + n1*(13 + 2*n2*(-69 + 7*n2*(-5 + 7*n2)))));
      double complex denominator = (392.*n1*(1 + n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex I1301(double complex n1, double mu)
{
      double complex numerator   = ((-7 + 9*n1)*ctan(n1*M_PI));
      double complex denominator = (336.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out         = numerator/denominator;

      return out;
} 

double complex Idelta201(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-8 + 7*n1 + 7*n2));
      double complex denominator = (28.*n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex IG201(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-2 + 7*n1 + 7*n2));
      double complex denominator = (56*n1*(1 + n1)*n2*(1 + n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex FG201(double complex n1, double complex n2, double mu)
{
      double complex numerator   =  ((-6 + 5*n1 + 14*cpow(n1,2))*ctan(n1*M_PI));
      double complex denominator = (448.*n1*(2 - n1 - 2*cpow(n1,2) + cpow(n1,3))*M_PI);
      double complex out         = numerator/denominator;

      return out;
} 

double complex J12101(double complex n1, double mu)
{
      double complex numerator   = (9*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6 + 11*n1 - 6*cpow(n1,2) + cpow(n1,3))*M_PI);
      double complex out         = numerator/denominator;

      return out;
} 

double complex J11201(double complex n1, double mu)
{
      return 0;
} 

double complex J21101(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-5 + n1*(-4 + 7*n1 + 7*n2)));
      double complex denominator = (14.*n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2);
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex Jdelta201(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2));
      double complex denominator = (n1*(-3 + 2*n1));
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex JG201(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -0.5*((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2));
      double complex denominator = (n1*(1 + n1)*(-3 + 2*n1)*n2);
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

/*
      _________________ 2-nd Moment _______________________
 */

double complex J12102(double complex n1, double mu)
{
      double complex numerator   = (9*(-1 + n1*cpow(mu,2))*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex J21102(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-((-1 + 2*n1)*(-1 + 2*n2)*(6 + 7*n1 + 7*n2)) + 
       (34 + n1 + n2 + 56*cpow(n1,3)*n2 - 54*cpow(n2,2) + 2*cpow(n1,2)*(-27 - 2*n2 + 56*cpow(n2,2)) +4*n1*n2*(-21 + n2*(-1 + 14*n2)))*cpow(mu,2))*Gamma(-2*n1)*Gamma(-2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*
     csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (28.*(1 + n1)*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex Jdelta202(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + (-1 + 2*n1 + 2*n2)*cpow(mu,2))*Gamma(2 - 2*n1)*Gamma(2 - 2*n2)*Gamma(2*(-2 + n1 + n2))*
     csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (4.*n1*n2*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex JG202(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -0.25*((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-1 + (1 + n1 + n2)*cpow(mu,2))*Gamma(2 - 2*n1)*
      Gamma(2 - 2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (n1*(1 + n1)*n2*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex I2211(double complex n1, double complex n2, double mu)
{
      double complex numerator   =((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(50 + 98*cpow(n1,3)*n2 - n2*(9 + 35*n2) + 
       7*cpow(n1,2)*(-5 + 2*n2*(-9 + 14*n2)) + n1*(-9 + 2*n2*(-33 + 7*n2*(-9 + 7*n2)))));
      double complex denominator = (392.*n1*(1 + n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2));
      double complex out = numerator/denominator * J(n1,n2);

      return out;
} 

double complex I1311(double complex n1, double mu)
{
      double complex numerator   = ((-5 + 3*n1)*ctan(n1*M_PI));
      double complex denominator = (448.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex J12111(double complex n1, double mu)
{
      double complex numerator   = (9*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6 + 11*n1 - 6*cpow(n1,2) + cpow(n1,3))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex J11211(double complex n1, double mu)
{
      return 0.;
} 

double complex J21111(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-3 + n1*(-8 + 7*n1 + 7*n2)));
      double complex denominator = (14.*n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2);
      double complex out = numerator/denominator * M1(n1,n2);

      return out;
} 

double complex N11a(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -4 - (8*(-2 + n1)*(-1 + n1)*(-4 + n1 + n2)*(-3 + 2*n2 + 4*(-3 + n1)*(-3 + n1 + n2)*cpow(mu,2)));
      double complex denominator = ((-3 + 2*n1)*(-9 + 2*n1 + 2*n2)*(-7 + 2*n1 + 2*n2)*(-3 + 2*n2 + (-5 + 2*n1)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)))*ctan(n1*M_PI) * ctan((n1 + n2)*M_PI);
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

double complex N11b(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(n1*(-3 + 2*n1) + n2*(-1 + 2*n2)));
      double complex denominator = (n1*(-3 + 2*n1)*n2*(-1 + 2*n2));
      double complex out = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex N11c(double complex n1, double complex n2, double mu)
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

double complex J21112(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-((-1 + 2*n1)*(-1 + 2*n2)*(-2 + 7*n1 + 7*n2)) + 
       (26 + 56*cpow(n1,3)*n2 + (9 - 38*n2)*n2 + 2*cpow(n1,2)*(-19 + 2*n2*(-9 + 28*n2)) + 
          n1*(9 + 4*n2*(-21 + n2*(-9 + 14*n2))))*cpow(mu,2))*Gamma(-2*n1)*Gamma(-2*n2)*Gamma(2*(-2 + n1 + n2))*
     csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (28.*(1 + n1)*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex J12112(double complex n1, double mu)
{
      double complex numerator   = (3*(3 - (-2 + n1)*(-3 + 4*n1)*cpow(mu,2))*ctan(n1*M_PI));
      double complex denominator = (224.*(-3 + n1)*(-2 + n1)*(-1 + n1)*n1*(1 + n1)*M_PI);
      double complex out = numerator/denominator;

      return out;
} 


double complex N12a(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-4 + n1 + n2)*(-3 + n1 + n2)*(-2 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*
     (-1 + 2*n2 + (-1 + 2*n1)*(1 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = (n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2*(-1 + 2*n2)*(-3 + 2*n2 + (-5 + 2*n1)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

double complex N12b(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -((-5 + n1 + n2)*(-4 + n1 + n2)*(-3 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*
       (-3 + 6*n2 + (-3 + 2*n1)*(1 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = (n1*(1 + n1)*(-5 + 2*n1)*(-3 + 2*n1)*n2*(-1 + 2*n2)*(-9 + 6*n2 + (-7 + 2*n1)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M3(n1, n2, mu);

      return out;
}

/*
      _________________ 4-th Moment _______________________
 */

double complex N22a(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-4 + n1 + n2)*(-3 + n1 + n2)*(-2 + n1 + n2)*(-1 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*
     (1 + 2*n1 + 2*n2)*(1 + 2*n2 + (-1 + 2*n1)*(3 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = (n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2)*(1 + 2*n2)*
     (-3 + 2*n2 + (-5 + 2*n1)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

double complex N22b(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-6 + n1 + n2)*(-5 + n1 + n2)*(-4 + n1 + n2)*(-3 + n1 + n2)*(-7 + 2*n1 + 2*n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*
     (-1 + 2*n1 + 2*n2)*(-3 + 12*cpow(n2,2) + 6*(-3 + 2*n1)*(1 + 2*n2)*(1 + 2*n1 + 2*n2)*cpow(mu,2) + 
       (-5 + 2*n1)*(-3 + 2*n1)*(1 + 2*n1 + 2*n2)*(3 + 2*n1 + 2*n2)*cpow(mu,4)));
      double complex denominator = (n1*(1 + n1)*(-5 + 2*n1)*(-3 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2)*(1 + 2*n2)*
     (3*(-5 + 2*n2)*(-3 + 2*n2) + 6*(-7 + 2*n1)*(-3 + 2*n2)*(-7 + 2*n1 + 2*n2)*cpow(mu,2) + 
       (-9 + 2*n1)*(-7 + 2*n1)*(-7 + 2*n1 + 2*n2)*(-5 + 2*n1 + 2*n2)*cpow(mu,4)));
      double complex out = numerator/denominator * M4(n1, n2, mu);

      return out;
} 

double complex N22c(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (-4*(-5 + n1 + n2)*(-4 + n1 + n2)*(-3 + n1 + n2)*(-2 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*
     (-1 + 2*n1 + 2*n2)*(1 + 2*n1 + 2*n2)*(3 + 6*n2 + (-3 + 2*n1)*(3 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = (n1*(1 + n1)*(-5 + 2*n1)*(-3 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2)*(1 + 2*n2)*
     (-9 + 6*n2 + (-7 + 2*n1)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M3(n1, n2, mu);

      return out;
} 