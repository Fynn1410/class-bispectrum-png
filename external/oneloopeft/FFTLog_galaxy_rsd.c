/** @file FFTLog_galaxy_real.c Documented FFT-Log based 1loop integrals of galaxy/halo cpow spectrum in perturbation theory 
 * 
 * Azadeh Moradinezhad Dizgah, June 16th 2021
 *
 * This module performs fast computation of the integrals appearing in the expression of 1loop galaxy/halo cpow sprtcurm.
 * The computation can be performed either in real or redshift-space. IR-resummation and EFT counter terms are included. 
 * The integrals are computed uccccsing FFTLog techniques.
 *
 * The algorithm closely follows Ref. arXiv:1708.08130 by Simonovic et al. After computing the FFT coefficents of matter cpow spectrum, 
 * sampled in logarithmic scale, the algorithm involves re-casting the integrals into a form that is analytically calculable (Matrices M_ij, 
 * which can be written in terms of ratios of Gamma functions) and finally performing vactor-matrix-vector or matrix-vector multiplications. 
 *
 * An imporccctant feauture of fast computation of loop integrals is that all the ccosmology-dependance of the loop integrals is captured by the FFT coeffcients
 * which have ~NlogN complexity. The matrices involving the Gamma functions, are computed only once, and at each of MCMC varying ccosmological parmaters, 
 * these coeffcients are evaluated for all k-values at once and then vector-matrix-vector multiplication are computed.
 *
 * The FFT coeffcients are computed uccccsing FFTW package, while the vector-matrix-vector computations are performed uccccsing blas library implemented in gsl.
 * The analytic formulas for M_ij matrices are computed in Mathematica uccccsing a modified version of the publicaly available notebook by Simpnovic. 
 * 
 * The choice of number of points for FFT decomposition of the cpow spectrum is very imporccctant, in terms of accuracy and execusion time. 
 * For halo loops, if chooccccsing 512 points *which takes ~ 0.45 seconds, I do get subpercent descrepencies between direct numeric integration of the loops.
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

void rsd_0_FFTLog(struct oneloop_fftlog_workspace *fft_ws, double k, double *np_loops, double *p_loops)
{
      int Nmax = fft_ws -> fft_input -> nfft;

      double *np = make_1Darray(6);
      double *p  = make_1Darray(2);
      
      // Linear cpow Spectrum vector for halos
      double complex *vec_h = make_1D_c_array(Nmax+1);
      vec_fill(fft_ws -> fft_input, k, HALO, vec_h);

      double complex *vec_h_min = make_1D_c_array(Nmax+1);
      vec_fill(fft_ws -> fft_input, fft_ws -> fft_input->kmin_fft_g, HALO, vec_h_min);

      // Linear cpow Spectrum vector for matter
      double complex *vec_m = make_1D_c_array(Nmax+1);
      vec_fill(fft_ws -> fft_input, k, MATTER, vec_m);
      
      // non-propagator calculations
      c_nonprop(vec_m, fft_ws -> fft_matrix -> I2200_mat,           vec_m, Nmax+1, &np[0]);
      c_nonprop(vec_h, fft_ws -> fft_matrix -> Idelta200_mat,       vec_h, Nmax+1, &np[1]);
      c_nonprop(vec_h, fft_ws -> fft_matrix -> IG200_mat,           vec_h, Nmax+1, &np[2]);
      c_nonprop(vec_h, fft_ws -> fft_matrix -> Idelta2delta200_mat, vec_h, Nmax+1, &np[3]);
      c_nonprop(vec_h, fft_ws -> fft_matrix -> IG2G200_mat,         vec_h, Nmax+1, &np[4]);
      c_nonprop(vec_h, fft_ws -> fft_matrix -> Idelta2G200_mat,     vec_h, Nmax+1, &np[5]);
      
      double Idelta2delta200_const;
      c_nonprop(vec_h_min, fft_ws -> fft_matrix -> Idelta2delta200_mat, vec_h_min, Nmax+1, &Idelta2delta200_const);

      // propagator calculations
      c_dot(vec_m, fft_ws -> fft_matrix -> I1300_mat, Nmax+1, &p[0]);
      c_dot(vec_h, fft_ws -> fft_matrix -> FG200_mat, Nmax+1, &p[1]);
           
      // adding factored out k and mu dependencies
      np_loops[0] = cpow(k, 3.) * np[0];
      np_loops[1] = cpow(k, 3.) * np[1];
      np_loops[2] = cpow(k, 3.) * np[2];
      np_loops[3] = cpow(k, 3.) * np[3] - cpow(fft_ws -> fft_input->kmin_fft_g, 3.) * Idelta2delta200_const ;
      np_loops[4] = cpow(k, 3.) * np[4];
      np_loops[5] = cpow(k, 3.) * np[5];

      p_loops[0] = cpow(k, 3.) * p[0];
      p_loops[1] = cpow(k, 3.) * p[1];
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

void rsd_1_FFTLog(struct oneloop_fftlog_workspace *fft_ws, double k, double mu, double *np_loops, double *p_loops)
{
      int Nmax     = fft_ws -> fft_input -> nfft;

      double *np = make_1Darray(7);
      double *p  = make_1Darray(3);

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_ws -> fft_input, k, HALO, vec);
      
      // non-propagator calculations
      c_nonprop(vec, fft_ws -> fft_matrix -> I2201_mat,     vec, Nmax+1, &np[0]);
      c_nonprop(vec, fft_ws -> fft_matrix -> Idelta201_mat, vec, Nmax+1, &np[1]);
      c_nonprop(vec, fft_ws -> fft_matrix -> IG201_mat,     vec, Nmax+1, &np[2]);
      c_nonprop(vec, fft_ws -> fft_matrix -> FG201_mat,     vec, Nmax+1, &np[3]);
      c_nonprop(vec, fft_ws -> fft_matrix -> J21101_mat,    vec, Nmax+1, &np[4]);
      c_nonprop(vec, fft_ws -> fft_matrix -> Jdelta201_mat, vec, Nmax+1, &np[5]);
      c_nonprop(vec, fft_ws -> fft_matrix -> JG201_mat,     vec, Nmax+1, &np[6]);

      // propagator calculations
      c_dot(vec, fft_ws -> fft_matrix -> I1301_mat,  Nmax+1, &p[0]);
      c_dot(vec, fft_ws -> fft_matrix -> J12101_mat, Nmax+1, &p[1]);
      c_dot(vec, fft_ws -> fft_matrix -> J11201_mat, Nmax+1, &p[2]);

      // adding factored out k and mu dependencies
      np_loops[0] = cpow(k, 3.) * np[0];
      np_loops[1] = cpow(k, 3.) * np[1];
      np_loops[2] = cpow(k, 3.) * np[2];
      np_loops[3] = cpow(k, 3.) * np[3];
      np_loops[4] = cpow(k, 3.) * np[4] * k * mu;
      np_loops[5] = cpow(k, 3.) * np[5] * k * mu;
      np_loops[6] = cpow(k, 3.) * np[6] * k * mu;

      p_loops[0] = cpow(k, 3.) * p[0];
      p_loops[1] = cpow(k, 3.) * p[1] * k * mu;
      p_loops[2] = cpow(k, 3.) * p[2] * k * mu;
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

void rsd_2_FFTLog(struct oneloop_fftlog_workspace *fft_ws, double k, double mu, double *np_loops, double *p_loops)
{
      int Nmax     = fft_ws -> fft_input -> nfft;

      double *np = make_1Darray(10);
      double *p  = make_1Darray(5);

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_ws -> fft_input, k, HALO, vec);

      // non-propagator calculations
      c_nonprop(vec, fft_ws -> fft_matrix -> J21102x_mat,    vec, Nmax+1, &np[0]);
      c_nonprop(vec, fft_ws -> fft_matrix -> J21102y_mat,    vec, Nmax+1, &np[1]);
      c_nonprop(vec, fft_ws -> fft_matrix -> Jdelta202x_mat, vec, Nmax+1, &np[2]);
      c_nonprop(vec, fft_ws -> fft_matrix -> Jdelta202y_mat, vec, Nmax+1, &np[3]);
      c_nonprop(vec, fft_ws -> fft_matrix -> JG202x_mat,     vec, Nmax+1, &np[4]);
      c_nonprop(vec, fft_ws -> fft_matrix -> JG202y_mat,     vec, Nmax+1, &np[5]);
      c_nonprop(vec, fft_ws -> fft_matrix -> I2211_mat,      vec, Nmax+1, &np[6]);
      c_nonprop(vec, fft_ws -> fft_matrix -> J21111_mat,     vec, Nmax+1, &np[7]);
      c_nonprop(vec, fft_ws -> fft_matrix -> N11x_mat,       vec, Nmax+1, &np[8]);
      c_nonprop(vec, fft_ws -> fft_matrix -> N11y_mat,       vec, Nmax+1, &np[9]);

      // propagator calculations
      c_dot(vec, fft_ws -> fft_matrix -> J12102x_mat, Nmax+1, &p[0]);
      c_dot(vec, fft_ws -> fft_matrix -> J12102y_mat, Nmax+1, &p[1]);
      c_dot(vec, fft_ws -> fft_matrix -> I1311_mat,   Nmax+1, &p[2]);
      c_dot(vec, fft_ws -> fft_matrix -> J12111_mat,  Nmax+1, &p[3]);
      c_dot(vec, fft_ws -> fft_matrix -> J11211_mat,  Nmax+1, &p[4]);

      // adding factored out k and mu dependencies
      np_loops[0] = cpow(k, 3.) * (np[0] + np[1] * cpow(mu,2.)) * cpow(k, 2.);
      np_loops[1] = cpow(k, 3.) * (np[2] + np[3] * cpow(mu,2.)) * cpow(k, 2.);
      np_loops[2] = cpow(k, 3.) * (np[4] + np[5] * cpow(mu,2.)) * cpow(k, 2.);
      np_loops[3] = cpow(k, 3.) * np[6];
      np_loops[4] = cpow(k, 3.) * np[7] * k * mu;
      np_loops[5] = cpow(k, 3.) * (np[8] + np[9] * cpow(mu,2.)) * cpow(k, 2.) ;

      p_loops[0] = cpow(k, 3.) * (p[0] + p[1] * cpow(mu,2.)) * cpow(k, 2.);
      p_loops[1] = cpow(k, 3.) * p[2];
      p_loops[2] = cpow(k, 3.) * p[3] * k * mu;
      p_loops[3] = cpow(k, 3.) * p[4] * k * mu;
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

void rsd_3_FFTLog(struct oneloop_fftlog_workspace *fft_ws, double k, double mu, double *np_loops, double *p_loops)
{
      int Nmax     = fft_ws -> fft_input -> nfft;

      double *np = make_1Darray(4);
      double *p  = make_1Darray(2);

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_ws -> fft_input, k, HALO, vec);

      // non-propagator calculations
      c_nonprop(vec, fft_ws -> fft_matrix -> J21112x_mat, vec, Nmax+1, &np[0]);
      c_nonprop(vec, fft_ws -> fft_matrix -> J21112y_mat, vec, Nmax+1, &np[1]);
      c_nonprop(vec, fft_ws -> fft_matrix -> N12x_mat,    vec, Nmax+1, &np[2]);
      c_nonprop(vec, fft_ws -> fft_matrix -> N12y_mat,    vec, Nmax+1, &np[3]);

      // propagator calculations
      c_dot(vec, fft_ws -> fft_matrix -> J12112x_mat, Nmax+1, &p[0]);
      c_dot(vec, fft_ws -> fft_matrix -> J12112y_mat, Nmax+1, &p[1]);

      // adding factored out k and mu dependencies
      np_loops[0] = cpow(k, 3.) * (np[0] + np[1] * cpow(mu,2.)) * cpow(k, 2.);
      np_loops[1] = cpow(k, 3.) * (np[0] * mu + np[1] * cpow(mu,3.)) * cpow(k, 3.) * mu;

      p_loops[0] = cpow(k, 3.) * (p[0] + p[1] * cpow(mu,2.)) * cpow(k, 2.);
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

void rsd_4_FFTLog(struct oneloop_fftlog_workspace *fft_ws, double k, double mu, double *np_loops)
{
      int Nmax     = fft_ws -> fft_input -> nfft;

      double *np = make_1Darray(3);

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_ws -> fft_input, k, HALO, vec);

      // non-propagator calculations
      c_nonprop(vec, fft_ws -> fft_matrix -> N22x_mat, vec, Nmax+1, &np[0]);
      c_nonprop(vec, fft_ws -> fft_matrix -> N22y_mat, vec, Nmax+1, &np[1]);
      c_nonprop(vec, fft_ws -> fft_matrix -> N22z_mat, vec, Nmax+1, &np[2]);

      // adding factored out k and mu dependencies
      np_loops[0] = cpow(k, 3.) * cpow(k, 4.) * (np[0] +  np[1] * cpow(mu, 2.) + np[2] * cpow(mu, 4.));
}

double P22_new(struct fft_struct *fft_input, double k, double z, int cleanup)
{
      double result;

      int Nmax     = fft_input -> nfft;
      double complex **M22_mat = make_2D_c_array(Nmax+1, Nmax+1);
      double complex *vec1 = make_1D_c_array(Nmax+1);

      vec_fill(fft_input, k, MATTER, vec1);
      np_mat_fill(I2200, fft_input, 0., MATTER, M22_mat);

      c_nonprop(vec1, M22_mat, vec1, Nmax+1, &result);

      return result * cpow(k,3.);
}

double P13_new(struct fft_struct *fft_input, double k, double z, int cleanup)
{
      double result;

      int Nmax     = fft_input -> nfft;
      double complex *M13_mat = make_1D_c_array(Nmax+1);
      double complex *vec1 = make_1D_c_array(Nmax+1);

      vec_fill(fft_input, k, MATTER, vec1);
      p_mat_fill(I1300, fft_input, 0., MATTER, M13_mat);

      c_dot(M13_mat, vec1, Nmax+1, &result);

      return result * cpow(k,3.);
}


/*
      _________________ 0-th Moment _______________________
 */

double complex I2200(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(58 + 98*cpow(n1,3)*n2 + (3 - 91*n2)*n2 + 
       7*cpow(n1,2)*(-13 - 2*n2 + 28*cpow(n2,2)) + n1*(3 + 2*n2*(-73 + 7*n2*(-1 + 7*n2)))));
      double complex denominator = (392.*n1*(1 + n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex I1300(double complex n1)
{
      double complex numerator   = ((1 + 9*n1)*ctan(n1*M_PI));
      double complex denominator = (672.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out         = numerator/denominator;

      return out;
}

double complex Idelta200(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-4 + 7*n1 + 7*n2));
      double complex denominator = (28.*n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex IG200(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(6 + 7*n1 + 7*n2));
      double complex denominator = (56*n1*(1 + n1)*n2*(1 + n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex Idelta2delta200(double complex n1, double complex n2)
{
      double complex out         = J(n1, n2);
      return out;
}

double complex IG2G200(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2));
      double complex denominator = (2.*n1*(1 + n1)*n2*(1 + n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex Idelta2G200(double complex n1, double complex n2)
{
      double complex numerator   = (3 - 2*n1 - 2*n2);
      double complex denominator = (2.*n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex FG200(double complex n1)
{
      double complex numerator    = (-15*ctan(n1*M_PI));
      double complex denominator  = (112.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out          = numerator/denominator;

      return out;

}

/*
      _________________ 1-st Moment _______________________
 */

double complex I2201(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(46 + 98*cpow(n1,3)*n2 + (13 - 63*n2)*n2 + 
       7*cpow(n1,2)*(-9 + 2*n2*(-5 + 14*n2)) + n1*(13 + 2*n2*(-69 + 7*n2*(-5 + 7*n2)))));
      double complex denominator = (392.*n1*(1 + n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex I1301(double complex n1)
{
      double complex numerator   = ((-7 + 9*n1)*ctan(n1*M_PI));
      double complex denominator = (336.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out         = numerator/denominator;

      return out;
} 

double complex Idelta201(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-8 + 7*n1 + 7*n2));
      double complex denominator = (28.*n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex IG201(double complex n1, double complex n2)
{
      double complex numerator   = -((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-2 + 7*n1 + 7*n2));
      double complex denominator = (56*n1*(1 + n1)*n2*(1 + n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex FG201(double complex n1, double complex n2)
{
      double complex numerator   =  ((-6 + 5*n1 + 14*cpow(n1,2))*ctan(n1*M_PI));
      double complex denominator = (448.*n1*(2 - n1 - 2*cpow(n1,2) + cpow(n1,3))*M_PI);
      double complex out         = numerator/denominator;

      return out;
} 

double complex J12101(double complex n1)
{
      double complex numerator   = (9*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6 + 11*n1 - 6*cpow(n1,2) + cpow(n1,3))*M_PI);
      double complex out         = numerator/denominator;

      return out;
} 

double complex J11201(double complex n1)
{
      return 0;
} 

double complex J21101(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-5 + n1*(-4 + 7*n1 + 7*n2)));
      double complex denominator = (14.*n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2);
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex Jdelta201(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2));
      double complex denominator = (n1*(-3 + 2*n1));
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex JG201(double complex n1, double complex n2)
{
      double complex numerator   = -0.5*((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2));
      double complex denominator = (n1*(1 + n1)*(-3 + 2*n1)*n2);
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

/*
      _________________ 2-nd Moment _______________________
 */

double complex J12102x(double complex n1)
{
      double complex numerator   = (-9*ctan(n1*M_PI));
      double complex denominator = (224.*(-3 + n1)*(-2 + n1)*(-1 + n1)*n1*(1 + n1)*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex J12102y(double complex n1)
{
      double complex numerator   = (3*(3 + 2*(-2 + n1)*n1)*ctan(n1*M_PI))/(224.*(-3 + n1)*(-2 + n1)*(-1 + n1)*n1*(1 + n1)*M_PI) + (3*(1 - 2*n1)*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(2 - n1 - 2*cpow(n1,2) + cpow(n1,3))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex J21102x(double complex n1, double complex n2)
{
      double complex numerator   = ((-1 + 2*n1)*(-1 + 2*n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(6 + 7*n1 + 7*n2)*Gamma(-2*n1)*
      Gamma(-2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (28*(1 + n1)*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex J21102y(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(34 + n1 + n2 + 56*cpow(n1,3)*n2 - 54*cpow(n2,2) + 
       2*cpow(n1,2)*(-27 - 2*n2 + 56*cpow(n2,2)) + 4*n1*n2*(-21 + n2*(-1 + 14*n2)))*Gamma(-2*n1)*Gamma(-2*n2)*
       Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (28.*(1 + n1)*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex Jdelta202x(double complex n1, double complex n2)
{
      double complex numerator   = -0.25*((-3 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*Gamma(2 - 2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (n1*n2*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex Jdelta202y(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*Gamma(2 - 2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (4.*n1*n2*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex JG202x(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*Gamma(2 - 2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (4.*n1*(1 + n1)*n2*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex JG202y(double complex n1, double complex n2)
{
      double complex numerator   = -0.25*((1 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*Gamma(2 - 2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (n1*(1 + n1)*n2*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex I2211(double complex n1, double complex n2)
{
      double complex numerator   =((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(50 + 98*cpow(n1,3)*n2 - n2*(9 + 35*n2) + 
       7*cpow(n1,2)*(-5 + 2*n2*(-9 + 14*n2)) + n1*(-9 + 2*n2*(-33 + 7*n2*(-9 + 7*n2)))));
      double complex denominator = (392.*n1*(1 + n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2));
      double complex out = numerator/denominator * J(n1,n2);

      return out;
} 

double complex I1311(double complex n1)
{
      double complex numerator   = ((-5 + 3*n1)*ctan(n1*M_PI));
      double complex denominator = (448.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex J12111(double complex n1)
{
      double complex numerator   = (9*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6 + 11*n1 - 6*cpow(n1,2) + cpow(n1,3))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex J11211(double complex n1)
{
      return 0.;
} 

double complex J21111(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-3 + n1*(-8 + 7*n1 + 7*n2)));
      double complex denominator = (14.*n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2);
      double complex out = numerator/denominator * M1(n1,n2);

      return out;
} 

double complex N11x(double complex n1, double complex n2)
{
      double complex numerator   = (Gamma(2 - 2*n2)*csin(n2*M_PI)*((-12*(-2 + n1)*ccos(n1*M_PI)*ccos((n1 + n2)*M_PI)*Gamma(3 - 2*n1)*
          Gamma(-7 + 2*n1 + 2*n2))/(-9 + 2*n1 + 2*n2) + 
       (8*(-2 + n1)*n2*ccos(n1*M_PI)*ccos((n1 + n2)*M_PI)*Gamma(3 - 2*n1)*Gamma(-7 + 2*n1 + 2*n2))/
        (-9 + 2*n1 + 2*n2) - ((-3 + 2*n1)*(-3 + 2*n2)*Gamma(2 - 2*n1)*Gamma(2*(-3 + n1 + n2))*csin(n1*M_PI)*
          csin((n1 + n2)*M_PI))/(-4 + n1 + n2) + 
       ((-3 + 2*n1)*(-5 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*Gamma(2*(-3 + n1 + n2))*csin(n1*M_PI)*
          csin((n1 + n2)*M_PI))/n2 + ((-3 + 2*n2)*(-5 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*Gamma(2*(-3 + n1 + n2))*
          csin(n1*M_PI)*csin((n1 + n2)*M_PI))/n1 - 
       (6*(-3 + 2*n1)*(-7 + 2*n1 + 2*n2)*Gamma(3 - 2*n1)*Gamma(-7 + 2*n1 + 2*n2)*csin(n1*M_PI)*
          csin((n1 + n2)*M_PI))/((-1 + n1)*(-4 + n1 + n2)) + 
       (4*(-3 + 2*n1)*n2*(-7 + 2*n1 + 2*n2)*Gamma(3 - 2*n1)*Gamma(-7 + 2*n1 + 2*n2)*csin(n1*M_PI)*
          csin((n1 + n2)*M_PI))/((-1 + n1)*(-4 + n1 + n2))));
      double complex denominator = (2.*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex N11y(double complex n1, double complex n2)
{
      double complex numerator   = (2*(-3 + n1 + n2)*(-3 + 2*n1 + 2*n2)*(n1*(-3 + 2*n1) + n2*(-1 + 2*n2))*M1(n1,n2))/
      (n1*(-3 + 2*n1)*n2*(-1 + 2*n2)) + (16*(-3 + n1)*(-2 + n1)*(-3 + n1 + n2)*ccos(n1*M_PI)*
        ccos((n1 + n2)*M_PI)*Gamma(3 - 2*n1)*Gamma(2 - 2*n2)*Gamma(-7 + 2*n1 + 2*n2)*csin(n2*M_PI))/
      ((-9 + 2*n1 + 2*n2)*cpow(M_PI,3)) - ((-5 + 2*n1)*(-3 + 2*n1)*(-5 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*
        Gamma(2 - 2*n2)*Gamma(2*(-3 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI))/
      ((-4 + n1 + n2)*cpow(M_PI,3)) + ((-3 + 2*n1)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*
        Gamma(2 - 2*n2)*Gamma(2*(-3 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI))/(n1*cpow(M_PI,3))\
      + ((-5 + 2*n1)*(-3 + 2*n1)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*Gamma(2 - 2*n2)*
        Gamma(2*(-3 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI))/(n2*(-1 + 2*n2)*cpow(M_PI,3)) + 
      (2*(-5 + 2*n1)*(-3 + 2*n1)*(-7 + 2*n1 + 2*n2)*(-5 + 2*n1 + 2*n2)*Gamma(3 - 2*n1)*Gamma(2 - 2*n2)*
        Gamma(-7 + 2*n1 + 2*n2)*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI))/
      ((-1 + n1)*(-4 + n1 + n2)*cpow(M_PI,3));
      double complex denominator = 2.;
      double complex out = numerator/denominator;

      return out;
} 


/*
      _________________ 3-rd Moment _______________________
 */

double complex J21112x(double complex n1, double complex n2)
{
      double complex numerator   = -((-1 + 2*n1)*(-1 + 2*n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(-2 + 7*n1 + 7*n2)*Gamma(-2*n1)*
       Gamma(-2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = 28*(1 + n1)*(1 + n2)*cpow(M_PI,3);
      double complex out = numerator/denominator;

      return out;
} 

double complex J21112y(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(26 + 56*cpow(n1,3)*n2 + (9 - 38*n2)*n2 + 
       2*cpow(n1,2)*(-19 + 2*n2*(-9 + 28*n2)) + n1*(9 + 4*n2*(-21 + n2*(-9 + 14*n2))))*Gamma(-2*n1)*Gamma(-2*n2)*
       Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*cpow(csin(n2*M_PI),2)*cpow(csin((n1 + n2)*M_PI),2));
      double complex denominator = (28.*(1 + n1)*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex J12112x(double complex n1)
{
      double complex numerator   = (-9*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex J12112y(double complex n1)
{
      double complex numerator   = (9*ctan(n1*M_PI));
      double complex denominator = (224.*(-6 + 5*n1 + 5*cpow(n1,2) - 5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 


double complex N12x(double complex n1, double complex n2)
{
      double complex numerator1   = (3*(-1 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*Gamma(-2*n2)*Gamma(2*(-1 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator1 = (8.*n1*(1 + n1)*(-2 + n1 + n2)*cpow(M_PI,3));
      
      double complex numerator2   = -(3*n2*(-1 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*Gamma(-2*n2)*Gamma(2*(-1 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator2 = (4.*n1*(1 + n1)*(-2 + n1 + n2)*cpow(M_PI,3));
      
      double complex numerator3   = (Gamma(-2*n1)*Gamma(2 - 2*n2)*Gamma(2*(n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator3 = (8.*(1 + n1)*n2*(-1 + n1 + n2)*cpow(M_PI,3));

      double complex out = numerator1/denominator1 + numerator2/denominator2 + numerator3/denominator3;

      return out;
} 

double complex N12y(double complex n1, double complex n2)
{
      double complex numerator1   = -0.125*((-3 + 2*n1)*(-1 + 2*n1 + 2*n2)*(1 + 2*n1 + 2*n2)*Gamma(2 - 2*n1)*Gamma(-2*n2)*Gamma(2*(-1 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator1 = (n1*(1 + n1)*(-2 + n1 + n2)*cpow(M_PI,3));

      double complex numerator2   = ((-1 + 2*n1)*(1 + 2*n1 + 2*n2)*Gamma(-2*n1)*Gamma(-2*n2)*Gamma(2*(n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator2 = (4.*(1 + n1)*(-1 + n1 + n2)*cpow(M_PI,3));

      double complex out = numerator1/denominator1 + numerator2/denominator2;

      return out;
}

/*
      _________________ 4-th Moment _______________________
 */

double complex N22x(double complex n1, double complex n2)
{
      double complex numerator1   = (-3*(-1 + 2*n1)*(-1 + 2*n1 + 2*n2)*Gamma(-2*n1)*Gamma(-2*(1 + n2))*Gamma(2*(-1 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator1 = (4.*(1 + n1)*(-2 + n1 + n2)*cpow(M_PI,3));

      double complex numerator2   = (3*(-1 + 2*n1)*cpow(n2,2)*(-1 + 2*n1 + 2*n2)*Gamma(-2*n1)*Gamma(-2*(1 + n2))*Gamma(2*(-1 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator2 = ((1 + n1)*(-2 + n1 + n2)*cpow(M_PI,3));

      double complex out = numerator1/denominator1 + numerator2/denominator2;

      return out;
} 

double complex N22y(double complex n1, double complex n2)
{
      double complex numerator1   = (3*(-3 + 2*n1)*(-1 + 2*n1)*(1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(1 + 2*n1 + 2*n2)*Gamma(-2*n1)*Gamma(-2*(1 + n2))*
      Gamma(2*(-1 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator1 = (2.*(1 + n1)*(-2 + n1 + n2)*cpow(M_PI,3));

      double complex numerator2   = - (3*(-1 + 2*n1)*(1 + 2*n1 + 2*n2)*Gamma(-2*n1)*Gamma(-2*(1 + n2))*Gamma(2*(n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*
      csin((n1 + n2)*M_PI));
      double complex denominator2 = (2.*(1 + n1)*(-1 + n1 + n2)*cpow(M_PI,3));

      double complex numerator3   = - (3*(-1 + 2*n1)*n2*(1 + 2*n1 + 2*n2)*Gamma(-2*n1)*Gamma(-2*(1 + n2))*Gamma(2*(n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*
      csin((n1 + n2)*M_PI));
      double complex denominator3 = ((1 + n1)*(-1 + n1 + n2)*cpow(M_PI,3));

      double complex numerator4   = (Gamma(-2*n1)*Gamma(-2*(1 + n2))*Gamma(2*(1 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator4 = (4.*(1 + n1)*(n1 + n2)*cpow(M_PI,3));

      double complex numerator5   = (n2*Gamma(-2*n1)*Gamma(-2*(1 + n2))*Gamma(2*(1 + n1 + n2))*csin(n1*M_PI)*
      csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator5 = (2.*(1 + n1)*(n1 + n2)*cpow(M_PI,3));

      double complex out = numerator1/denominator1 + numerator2/denominator2 + numerator3/denominator3 + numerator4/denominator4 + numerator5/denominator5;

      return out;
} 

double complex N22z(double complex n1, double complex n2)
{
      double complex numerator1   = ((-5 + 2*n1)*(-3 + 2*n1)*(-1 + 2*n1)*(-1 + 2*n1 + 2*n2)*(1 + 2*n1 + 2*n2)*(3 + 2*n1 + 2*n2)*Gamma(-2*n1)*
      Gamma(-2*(1 + n2))*Gamma(2*(-1 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator1 = (4.*(1 + n1)*(-2 + n1 + n2)*cpow(M_PI,3));

      double complex numerator2   = - ((-3 + 2*n1)*(-1 + 2*n1)*(1 + 2*n1 + 2*n2)*(3 + 2*n1 + 2*n2)*
      Gamma(-2*n1)*Gamma(-2*(1 + n2))*Gamma(2*(n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator2 = (2.*(1 + n1)*(-1 + n1 + n2)*cpow(M_PI,3));

      double complex numerator3   = ((-1 + 2*n1)*(3 + 2*n1 + 2*n2)*Gamma(-2*n1)*Gamma(-2*(1 + n2))*
      Gamma(2*(1 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator3 = (4.*(1 + n1)*(n1 + n2)*cpow(M_PI,3));

      double complex out = numerator1/denominator1 + numerator2/denominator2 + numerator3/denominator3;

      return out;
} 