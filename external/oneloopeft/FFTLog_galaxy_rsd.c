/** @file FFTLog_galaxy_real.c Documented FFT-Log based 1loop integrals of galaxy/halo cpow spectrum in perturbation theory 
 * 
 * Azadeh Moradinezhad Dizgah, June 16th 2021
 *
 * This module performs fast computation of the integrals appearing in the expression of 1loop galaxy/halo cpow sprtcurm.
 * The computation can be performed either in real or redshift-space. IR-resummation and EFT counter terms are included. 
 * The integrals are computed ucsing FFTLog techniques.
 *
 * The algorithm closely follows Ref. arXiv:1708.08130 by Simonovic et al. After computing the FFT coefficents of matter cpow spectrum, 
 * sampled in logarithmic scale, the algorithm involves re-casting the integrals into a form that is analytically calculable (Matrices M_ij, 
 * which can be written in terms of ratios of Gamma functions) and finally performing vactor-matrix-vector or matrix-vector multiplications. 
 *
 * An imporctant feauture of fast computation of loop integrals is that all the cosmology-dependance of the loop integrals is captured by the FFT coeffcients
 * which have ~NlogN complexity. The matrices involving the Gamma functions, are computed only once, and at each of MCMC varying cosmological parmaters, 
 * these coeffcients are evaluated for all k-values at once and then vector-matrix-vector multiplication are computed.
 *
 * The FFT coeffcients are computed ucsing FFTW package, while the vector-matrix-vector computations are performed ucsing blas library implemented in gsl.
 * The analytic formulas for M_ij matrices are computed in Mathematica ucsing a modified version of the publicaly available notebook by Simpnovic. 
 * 
 * The choice of number of points for FFT decomposition of the cpow spectrum is very imporctant, in terms of accuracy and execusion time. 
 * For halo loops, if choocsing 512 points *which takes ~ 0.45 seconds, I do get subpercent descrepencies between direct numeric integration of the loops.
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
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the M_PIpeline
 * @param pg_loops      Output: an array containing the values of galaxy loop corrections
 * @return void           
 */

void rsd_0_FFTLog(struct fft_struct *fft_input, double k, double *loops)
{
      int Nmax     = fft_input -> nfft;

      // Linear cpow Spectrum vector
      double complex *vec = make_1D_c_array(Nmax+1);
      vec_fill(fft_input, k, vec);

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
      // ISPsi is miscsing
      
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
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex IF2S2(double complex n1, double complex n2)
{
      double complex numerator   = ((3. - 2.*n1 - 2.*n2)*(-1. + 2.*n1 + 2.*n2)*(6. + 7.*n1 + 7.*n2));
      double complex denominator = 14.*n1*(1 + n1)*n2*(1 + n2);
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
      double complex numerator   = (3. - 2.*n1 - 2.*n2);
      double complex denominator = (n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

double complex IS2S2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3. + 2.*n1 + 2.*n2)*(-1. + 2.*n1 + 2.*n2));
      double complex denominator = (n1*(1. + n1)*n2*(1. + n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
}

/*
      _________________ 1-st Moment _______________________
 */

double complex IF2G2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3. + 2.*n1 + 2.*n2)*(-1. + 2.*n1 + 2.*n2)*(46. + 98.*cpow(n1,3.)*n2 + (13. - 63.*n2)*n2 + 
       7.*cpow(n1,2.)*(-9. + 2.*n2*(-5. + 14.*n2)) + n1*(13. + 2.*n2*(-69. + 7.*n2*(-5. + 7.*n2)))));
      double complex denominator = (196.*n1*(1. + n1)*(-1. + 2.*n1)*n2*(1. + n2)*(-1. + 2.*n2));
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex IF3G3(double complex n1)
{
      double complex numerator   = ((-7. + 9.*n1)*ctan(n1*M_PI));
      double complex denominator = (112.*n1*(-6. + 5.*n1 + 5.*cpow(n1,2.) - 5.*cpow(n1,3.) + cpow(n1,4.))*M_PI);
      double complex out         = numerator/denominator;

      return out;
} 

double complex IG2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3. + 2.*n1 + 2.*n2)*(-8. + 7.*n1 + 7.*n2));
      double complex denominator = (28.*n1*n2);
      double complex out         = numerator/denominator * J(n1, n2);

      return out;
} 

double complex IS2G2(double complex n1, double complex n2)
{
      double complex numerator   = -(-3. + 2.*n1 + 2.*n2)*(-1. + 2.*n1 + 2.*n2)*(-2. + 7.*n1 + 7.*n2);
      double complex denominator = (28.*n1*(1. + n1)*n2*(1. + n2));
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

double complex IF2p2(double complex n1, double complex n2)
{
      return 0;
} 

double complex IG2p(double complex n1)
{
      double complex numerator   = (3.*(9. - 4.*n1)*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6. + 11.*n1 - 6.*cpow(n1,2.) + cpow(n1,3.))*M_PI);
      double complex out         = numerator/denominator;

      return out;
} 

double complex IPPp(double complex n1, double complex n2)
{
      double complex numerator   = -(18. - 21.*n1 + 6.*cpow(n1,2.) - 18.*n2 + 8.*n1*n2 + 4.*cpow(n2,2.));
      double complex denominator = (6.*n1 - 4.*cpow(n1,2.));
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex IS2p(double complex n1, double complex n2)
{
      double complex numerator   = -0.5*((-3. + n1 + n2)*(-3. + 2.*n1 + 2.*n2)*(-1. + 2.*n1 + 2.*n2));
      double complex denominator = n1*(1. + n1)*(-3. + 2.*n1)*n2;
      double complex out         = numerator/denominator * M1(n1, n2);

      return out;
} 

/*
      _________________ 2-nd Moment _______________________
 */

double complex IF2p21(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -((-4. + n1 + n2)*(-3. + n1 + n2)*(-5. + 2.*n1 + 2.*n2)*(-3. + 2.*n1 + 2.*n2)*
      (-1. + 2.*n1 + 2.*n2)*((-1. + 2.*n1)*(6. + 7.*n1 + 7.*n2) + 
        8.*(-6. + 14.*cpow(n1,2.)*(-1. + 2.*n2) + n1*(-59. - 30.*n2 + 56.*cpow(n2,2.)) + 
           n2*(-47. + 4.*n2*(-4. + 7.*n2)))*cpow(mu,2.)));
      double complex denominator = (14.*n1*(1 + n1)*(-1. + 2.*n1)*n2*(1. + n2)*(-3. + 2.*n2)*
      (-3. + 2.*n1 + 8.*(-5. + 2.*n2)*(-5. + 2.*n1 + 2.*n2)*cpow(mu,2.)));
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

double complex IF2p22(double complex n1, double complex n2)
{
      double complex numerator   = -((-3. + n1 + n2)*(-2. + n1 + n2)*(-3. + 2.*n1 + 2.*n2)*(-1. + 2.*n1 + 2.*n2)*
      (1. + 2.*n1 + 2.*n2)*(-5. + n1*(3. + 7.*n1 + 7.*n2)));
      double complex denominator = (7.*n1*(1. + n1)*(-3. + 2.*n1)*(-1. + 2.*n1)*n2*(1. + n2)*(-1. + 2.*n2));
      double complex out = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex IG2pm1(double complex n1, double complex n2)
{
      double complex numerator   = -((-3 + n1 + n2)*(-2 + n1 + n2)*
       (-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2));
      double complex denominator = (56.*n1*(2 - n1 - 2*cpow(n1,2) + 
       cpow(n1,3))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex IG2pm2(double complex n1, double mu)
{
      double complex numerator   = (9*(1 - 24*cpow(mu,2))*ctan(n1*M_PI));
      double complex denominator = (56.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 
       5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex IPPpm1(double complex n1, double complex n2)
{
      double complex numerator   = -((-3 + n1 + n2)*(-2 + n1 + n2)*
       (-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2));
      double complex denominator = (n1*(-3 + 2*n1)*n2*(-1 + 2*n2));
      double complex out = numerator/denominator * M1(n1, n2);

      return out;
} 

double complex IPPpm2(double complex n1, double complex n2, double mu)
{
      double complex numerator   = ((-4 + n1 + n2)*(-3 + n1 + n2)*
     (-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*
     (-1 + 2*n1 + 8*(-3 + 2*n2)*
        (-1 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = (n1*(-1 + 2*n1)*n2*(-3 + 2*n2)*
     (-3 + 2*n1 + 8*(-5 + 2*n2)*
        (-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

double complex IS2pm(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (2*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*
     (-1 + n1*(2 - 4*n2) + 2*n2 + (-6 + 4*cpow(n1,3) + 8*cpow(n1,2)*(-1 + n2) + n1*(-33 + 14*(5 - 2*n2)*n2) + 
          n2*(3 + 46*n2 - 32*cpow(n2,2)))*cpow(mu,2))*Gamma(-2*n1)*Gamma(-2*n2)*Gamma(2*(-2 + n1 + n2))*csin(n1*M_PI)*
     csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = ((1 + n1)*(1 + n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex IG2G2(double complex n1, double complex n2)
{
      double complex numerator   = ((-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(50 + 98*cpow(n1,3)*n2 - n2*(9 + 35*n2) + 
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
      double complex numerator   = (3*(-3 + 2*n1 + 6*mu - 2*n1*mu)*ctan(n1*M_PI));
      double complex denominator = (224.*n1*(-6 + 11*n1 - 6*cpow(n1,2) + 
       cpow(n1,3))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex IF2pm(double complex n1)
{
      return 0.;
} 

double complex IPPmm1(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -4 - 8*(-1 + n1)*(-4 + n1 + n2)*(-2 + n1 + 8*(-3 + n1 + n2)*(-5 + 2*n2)*cpow(mu,2))/ctan(n1*M_PI);
      double complex denominator = ((-9 + 2*n1 + 2*n2)*(-7 + 2*n1 + 2*n2)*(-3 + 2*n1 + 8*(-5 + 2*n2)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
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
      double complex numerator   = (((-1 + 2*n1)*(4*cpow(n1,4) + 8*cpow(n1,3)*(-4 + n2) + (-4 + n2)*n2*(-5 + 2*n2)*(-3 + 2*n2) + 
          cpow(n1,2)*(79 - 38*n2 + 4*cpow(n2,2)) + n1*(-5 + 2*n2)*(12 + n2*(-9 + 4*n2))) + 
       8*(-3 + 2*n2)*(-5 + 2*n1 + 2*n2)*((-4 + n1)*n1*(-3 + 2*n1)*(-1 + 2*n1) + 
          (-3 + 2*n1)*(20 + n1*(-7 + 4*n1))*n2 + (79 - 42*n1 + 4*cpow(n1,2))*cpow(n2,2) + 
          8*(-4 + n1)*cpow(n2,3) + 4*cpow(n2,4))*cpow(mu,2)));
      double complex denominator = n1*(-1 + 2*n1)*n2*(-3 + 2*n2)*(-3 + 2*n1 + 8*(-5 + 2*n2)*(-5 + 2*n1 + 2*n2)*cpow(mu,2));
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

/*
      _________________ 3-rd Moment _______________________
 */

double complex IG2pm3(double complex n1, double complex n2, double mu)
{
      double complex numerator   = -((-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*mu*
      ((-1 + 2*n1)*(1 + n2)*(-3 + 2*n2)*(-1 + 2*n2)*(-12 + 7*n2) + 
        2*(n1*(855 + (-2 + n1)*n1*(58 + n1*(-7 + 4*n1*(-22 + 7*n1)))) + 
           n1*(-2753 + n1*(204 + n1*(495 + 4*n1*(-144 + 35*n1))))*n2 + 
           (-4665 + n1*(1256 + n1*(471 + 8*n1*(-108 + 35*n1))))*cpow(n2,2) + 
           (-292 + n1*(2229 + 8*n1*(-72 + 35*n1)))*cpow(n2,3) + 4*(1109 + n1*(-524 + 35*n1))*cpow(n2,4) + 
           4*(-656 + 119*n1)*cpow(n2,5) + 448*cpow(n2,6) + 69*(-6 + 41*n2))*cpow(mu,2))*Gamma(3 - 2*n1)*
      Gamma(2 - 2*n2)*Gamma(2*(-3 + n1 + n2))*csin(n1*M_PI)*csin(n2*M_PI)*csin((n1 + n2)*M_PI));
      double complex denominator = (28*(-1 + n1)*n1*(1 + n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2)*cpow(M_PI,3));
      double complex out = numerator/denominator;

      return out;
} 

double complex IG2pmm1(double complex n1, double mu)
{
      double complex numerator   = (3*(1 - 2*n1)*mu*ctan(n1*M_PI));
      double complex denominator = (56.*n1*(2 - n1 - 2*cpow(n1,2) + 
       cpow(n1,3))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex IG2pmm2(double complex n1, double mu)
{
      double complex numerator   = (9*(-1 + 24*cpow(mu,2))*
     ctan(n1*M_PI));
      double complex denominator = (56.*n1*(-6 + 5*n1 + 5*cpow(n1,2) - 
       5*cpow(n1,3) + cpow(n1,4))*M_PI);
      double complex out = numerator/denominator;

      return out;
} 

double complex IPPpm31(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (2*(-4 + n1 + n2)*(-3 + n1 + n2)*(-2 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*
     (1 + 2*n1 + 8*(-3 + 2*n2)*(1 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = (n1*(1 + n1)*(-1 + 2*n1)*(1 + 2*n1)*n2*(-3 + 2*n2)*(-3 + 2*n1 + 8*(-5 + 2*n2)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

double complex IPPpm32(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (-2*(-5 + n1 + n2)*(-4 + n1 + n2)*(-3 + n1 + n2)*(-7 + 2*n1 + 2*n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*
     (6*(-1 + 2*n1)*(-1 + 2*n2)*(-1 + 2*n1 + 2*n2) + 
       (32*cpow(n1,4) + 64*cpow(n1,3)*(-1 + n2) + cpow(n1,2)*(16 + n2*(-143 + 62*n2)) + 
          6*(-1 + 2*n2)*(1 + n2*(12 + 5*(-3 + n2)*n2)) + n1*(16 + n2*(108 + n2*(-229 + 90*n2))))*cpow(mu,2)));
      double complex denominator = (n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2*(-1 + 2*n2)*
     (6*(-5 + 2*n1)*(-3 + 2*n2)*(-7 + 2*n1 + 2*n2) + 
       (n1*(-3855 + n1*(1757 + 32*(-12 + n1)*n1)) + 50*(76 - 93*n2) + n1*(2672 + n1*(-651 + 64*n1))*n2 + 
          (2470 + n1*(-747 + 62*n1))*cpow(n2,2) + 90*(-7 + n1)*cpow(n2,3) + 60*cpow(n2,4))*cpow(mu,2)));
      double complex out = numerator/denominator * M3(n1, n2, mu);

      return out;
}

/*
      _________________ 4-th Moment _______________________
 */

double complex IPPm41(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (2*(-4 + n1 + n2)*(-3 + n1 + n2)*(-2 + n1 + n2)*(-1 + n1 + n2)*(-5 + 2*n1 + 2*n2)*(-3 + 2*n1 + 2*n2)*
     (-1 + 2*n1 + 2*n2)*(1 + 2*n1 + 2*n2)*cpow(mu,2)*(1 + 2*n1 + 8*(-1 + 2*n2)*(3 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex denominator = (n1*(1 + n1)*(-1 + 2*n1)*(1 + 2*n1)*n2*(1 + n2)*(-3 + 2*n2)*(-1 + 2*n2)*
     (-3 + 2*n1 + 8*(-5 + 2*n2)*(-5 + 2*n1 + 2*n2)*cpow(mu,2)));
      double complex out = numerator/denominator * M2(n1, n2, mu);

      return out;
} 

double complex IPPm42(double complex n1, double complex n2, double mu)
{
      double complex numerator   = (-2*(-6 + n1 + n2)*(-5 + n1 + n2)*(-4 + n1 + n2)*(-3 + n1 + n2)*(-7 + 2*n1 + 2*n2)*(-5 + 2*n1 + 2*n2)*
     (-3 + 2*n1 + 2*n2)*(-1 + 2*n1 + 2*n2)*(3*(-1 + 2*n1)*(-1 + 4*cpow(n2,2)) - 
       6*(1 + 2*n2)*(3 + 3*(-3 + n2)*n2*(1 + n2) + cpow(n1,3)*(-5 + 3*n2) + 
          cpow(n1,2)*(6 - 8*n2 + 6*cpow(n2,2)) + n1*(-1 + 10*n2 + 3*cpow(n2,3)))*cpow(mu,2) + 
       (-5 + 2*n1)*(-3 + 2*n1)*(-1 + 2*n1)*(1 + 2*n1 + 2*n2)*(3 + 2*n1 + 2*n2)*cpow(mu,4)));
      double complex denominator = (n1*(1 + n1)*(-3 + 2*n1)*(-1 + 2*n1)*n2*(1 + n2)*(-1 + 2*n2)*(1 + 2*n2)*
     (3*(-5 + 2*n1)*(-5 + 2*n2)*(-3 + 2*n2) - 6*(-3 + 2*n2)*
        (335 + n1*(-361 + (112 - 11*n1)*n1) - 193*n2 + n1*(210 + n1*(-50 + 3*n1))*n2 + 
          6*(-6 + n1)*(-1 + n1)*cpow(n2,2) + 3*(-1 + n1)*cpow(n2,3))*cpow(mu,2) + 
       (-9 + 2*n1)*(-7 + 2*n1)*(-5 + 2*n1)*(-7 + 2*n1 + 2*n2)*(-5 + 2*n1 + 2*n2)*cpow(mu,4)));
      double complex out = numerator/denominator * M4(n1, n2, mu);

      return out;
} 