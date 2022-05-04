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
      _________________ 1-st Moment _______________________
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

void nprop_rsd_1_FFTLog(struct fft_struct  *fft_input, double k, double z, int cleanup, double *pg_loops)
{
      int Nmax     = fft_input -> nfft;

      static double P_b2b2_ls = 0.;
      static gsl_matrix_complex *M_cIF2G2_mat;
      static gsl_matrix_complex *M_cIG2_mat;
      static gsl_matrix_complex *M_cIS2G2_mat;

      gsl_complex out_cIF2G2_complex;
      gsl_complex out_cIG2_complex;
      gsl_complex out_cIS2G2_complex;

      gsl_vector_complex *matvec_cIF2G2 = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *matvec_cIG2   = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *matvec_cIS2G2 = gsl_vector_complex_alloc(Nmax+1);

      double complex nu1, nu2;
      gsl_complex cmk_gsl;
      gsl_complex cmkmin_gsl;

      gsl_vector_complex *cmk_vec     = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *cmkmin_vec  = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *matvec_cIPP = gsl_vector_complex_alloc(Nmax+1);

      int i, j;
      static int first = 1;
      if(first==1)
      {
         gsl_complex M_cIF2G2_gsl;
         gsl_complex M_cIG2_gsl;
         gsl_complex M_cIS2G2_gsl;

         M_cIF2G2_mat = gsl_matrix_complex_alloc(Nmax+1,Nmax+1);
         M_cIG2_mat   = gsl_matrix_complex_alloc(Nmax+1,Nmax+1);
         M_cIS2G2_mat = gsl_matrix_complex_alloc(Nmax+1,Nmax+1);


         for (i=0; i<Nmax+1; i++){
               for (j=i; j<Nmax+1; j++){                        
                     nu1 = -0.5 * fft_input->etam_g[i];
                     nu2 = -0.5 * fft_input->etam_g[j];                     

                     GSL_SET_COMPLEX(&M_cIF2G2_gsl,creal(M_cIF2G2(nu1,nu2)),cimag(M_cIF2G2(nu1,nu2)));
                     gsl_matrix_complex_set(M_cIF2G2_mat, i, j, M_cIF2G2_gsl);
                     gsl_matrix_complex_set(M_cIF2G2_mat, j, i, M_cIF2G2_gsl); 

                     GSL_SET_COMPLEX(&M_cIG2_gsl,creal(M_cIG2(nu1,nu2)),cimag(M_cIG2(nu1,nu2)));
                     gsl_matrix_complex_set(M_cIG2_mat, i, j, M_cIG2_gsl);
                     gsl_matrix_complex_set(M_cIG2_mat, j, i, M_cIG2_gsl); 

                     GSL_SET_COMPLEX(&M_cIS2G2_gsl,creal(M_cIS2G2(nu1,nu2)),cimag(M_cIS2G2(nu1,nu2)));
                     gsl_matrix_complex_set(M_cIS2G2_mat, i, j, M_cIS2G2_gsl);
                     gsl_matrix_complex_set(M_cIS2G2_mat, j, i, M_cIS2G2_gsl); 

               }
         }
         first =0;
      }


      for (i=0; i<Nmax+1; i++){ 
         nu1 = -0.5 * fft_input->etam_g[i];                       
         GSL_SET_COMPLEX(&cmkmin_gsl, creal(fft_input->cmsym_g[i] * cpow(fft_input->kmin_fft,-2.*nu1)), cimag(fft_input->cmsym_g[i] * cpow(fft_input->kmin_fft,-2.*nu1)));
         gsl_vector_complex_set(cmkmin_vec, i, cmkmin_gsl);
      }

      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M_cIPP_mat, cmkmin_vec, GSL_COMPLEX_ZERO, matvec_cIPP);
      gsl_blas_zdotu(cmkmin_vec, matvec_cIPP, &out_cIdelta2delta20_complex);

      gsl_vector_complex_free(cmkmin_vec);
      gsl_vector_complex_free(matvec_cIPP);


      for (i=0; i<Nmax+1; i++){ 
            nu1 = -0.5 * fft_input->etam_g[i];                       
            GSL_SET_COMPLEX(&cmk_gsl, creal(fft_input->cmsym_g[i] * cpow(k,-2.*nu1)), cimag(fft_input->cmsym_g[i] * cpow(k,-2.*nu1)));
            gsl_vector_complex_set(cmk_vec, i, cmk_gsl);
      }

      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M_cIF2G2_mat, cmk_vec,    GSL_COMPLEX_ZERO, matvec_cIF2G2);
      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M_cIG2_mat,   cmk_vec,    GSL_COMPLEX_ZERO, matvec_cIG2);
      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M_cIS2G2_mat, cmk_vec,    GSL_COMPLEX_ZERO, matvec_cIS2G2);

      gsl_blas_zdotu(cmk_vec, matvec_cIF2G2, &out_cIF2G2_complex);
      gsl_blas_zdotu(cmk_vec, matvec_cIG2,   &out_cIG2_complex);
      gsl_blas_zdotu(cmk_vec, matvec_cIS2G2, &out_cIS2G2_complex);

      pg_loops[0] = pow(k,3.) * GSL_REAL(out_cIF2G2_complex);
      pg_loops[1] = pow(k,3.) * GSL_REAL(out_cIG2_complex);
      pg_loops[2] = pow(k,3.) * GSL_REAL(out_cIS2G2_complex);

      gsl_vector_complex_free(cmk_vec);

      gsl_vector_complex_free(matvec_cIF2G2);
      gsl_vector_complex_free(matvec_cIF3G3);
      gsl_vector_complex_free(matvec_cIG2);
      gsl_vector_complex_free(matvec_cIS2G2);

      if(cleanup == 1){
            gsl_matrix_complex_free(M_cIF2G2_mat);
            gsl_matrix_complex_free(M_cIG2_mat);
            gsl_matrix_complex_free(M_cIS2G2_mat);

      }
      // printf("suceessfully computed galaxy loops in real space using FFTLog\n");

      return;
}

/**
 Compute propagator type loop contribution to non-linear galaxy power spectrum given the FFTLog coefficents and frequencies (Eq. 2.40 of Simonovic 2017)
 *
 * @param fft_input    Input: structure containing fft coefficents and params
 * @param k             Input: wavenumber in unit of h/Mpc. 
 * @param z             Input: redshift
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the pipeline
 * @return value of IF3G3 and IG2m term in unit of (Mpc/h)^3 as a function of redshift and wavenumber
 */

void prop_rsd_1_FFTLog(struct fft_struct  *fft_input, double k, double z, int cleanup, double *pg_loops)
{
      fftw_complex McIF3G3_elems;
      fftw_complex McIG2m_elems;

      int Nmax     = fft_input -> nfft;

      static gsl_vector_complex *M_cIF3G3_vec;
      static gsl_vector_complex *M_cIG2m_vec;
      gsl_vector_complex *cmk_vec  = gsl_vector_complex_alloc(Nmax+1);

      gsl_complex cmk_gsl;
      double complex nu1;

      /*
      * The M_ij matrices are cosmology-independant and need to be computed only once. So we use the first switch
      */
      static int first = 1;
      if(first == 1)
      {
         gsl_complex M_cIF3G3_gsl, M_cIG2m_gsl;
         M_cIF3G3_vec = gsl_vector_complex_alloc(Nmax+1);
         M_cIG2m_vec  = gsl_vector_complex_alloc(Nmax+1);

         for (int i=0; i<Nmax+1; i++){
               nu1 = -0.5 * fft_input->etam_g[i];
               GSL_SET_COMPLEX(&M_cIF3G3_gsl,creal(M_cIF3G3(nu1)),cimag(M_cIF3G3(nu1)));
               GSL_SET_COMPLEX(&M_cIG2m_gsl,creal( M_cIG2m(nu1)), cimag(M_cIG2m(nu1)));
               // printf("M_cIF3G3_gsl %d %12.6e %12.6e %12.6e  \n", i, k,  GSL_REAL(M_cIF3G3_gsl), GSL_IMAG(M_cIF3G3_gsl));

               gsl_vector_complex_set(M_cIF3G3_vec, i, M_cIF3G3_gsl);
               gsl_vector_complex_set(M_cIG2m_vec, i, M_cIG2m_gsl);
         }
         first = 0;
      }         
      for (int i=0; i<Nmax+1; i++){
             nu1 = -0.5 * fft_input->etam_g[i];
             GSL_SET_COMPLEX(&cmk_gsl, creal(fft_input->cmsym_g[i] * cpow(k,-2.*nu1)), cimag(fft_input->cmsym_g[i] * cpow(k,-2.*nu1)));
             // printf("cmk13 %d %12.6e %12.6e %12.6e  \n", i, k,  GSL_REAL(cmk_gsl), GSL_IMAG(cmk_gsl));

            gsl_vector_complex_set(cmk_vec, i, cmk_gsl);
      }
      gsl_blas_zdotu(cmk_vec, M_cIF3G3_vec, &out_cIF3G3_complex);
      gsl_blas_zdotu(cmk_vec, M_cIG2m_vec, &out_cIG2m_complex);
      pg_loops[0] = pow(k,3.) * GSL_REAL(out_cIF3G3_complex);
      pg_loops[1] = pow(k,3.) * GSL_REAL(out_cIG2m_complex);

      gsl_vector_complex_free(cmk_vec);

      if(cleanup == 1)
         gsl_vector_complex_free(M_cIF3G3_vec);
         gsl_vector_complex_free(M_cIG2m_vec);

      // printf("suceessfully computed pgloops_propag using FFTLog\n");

      return;
}


/**
 * Analytic expression of M_cIF2G2 matrix
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cIF2G2
 */

double complex M_cIF2G2(double complex nu1, double complex nu2)
{
      double complex numerator1  = (-3. + 2. * (nu1 + nu2));
      double complex numerator2  = (-1. + 2. * (nu1 + nu2));
      double complex numerator3  = (46. + 98. * pow(nu1, 3) * nu2 + (13. - 63.*nu2)*nu2 + 7.*pow(n1, 2)*(-9. + 2. * nu2*(-5. + 14.*n2)) + nu1*(13 + 2*nu2*(-69. + 7.*nu2*(-5. + 7.*n2)))
      double complex denominator = 196. * nu1 * (1. + nu1) * (-1. +2.*nu1) * nu2 * (1. + nu2) * (-1. + 2.*nu2);
      double complex out         = numerator1 * numerator2 * numerator3 /denominator * Ifunc(nu1, nu2);

      return out;
}

/**
 * Analytic expression of M_cIG2 matrix
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cIG2
 */

double complex M_cIG2(double complex nu1, double complex nu2)
{
      double complex numerator   = (-3. + 2.* (nu1 + nu2)) * (-8. + 7.* (nu1 + nu2));
      double complex denominator = 28. * nu1 * nu2;
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}

/**
 * Analytic expression of M_cIS2G2
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cIS2G2
 */

double complex M_cIS2G2(double complex nu1, double complex nu2)
{
      double complex numerator   = (-3. + 2. * (nu1 + nu2)) *  (1. - 2. * (nu1 + nu2)) * (-2. + 7.* (nu1 + nu2));
      double complex denominator = 28. * nu1 * (1. + nu1) * nu2 * (1. + nu2);
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}

/**
 * Analytic expression of M_cIF3G3 matrix
 *
 * @param nu1         Input: FFTLog coeffcients
 * @return value of M_cIF3G3
 */

double complex M_cIF3G3(double complex nu1)
{
      double complex numerator   = (-7. + 9.* nu1) * ctan(nu1 * M_PI);
      double complex denominator = 112. * nu1 * (-6. + 5.*nu1 + 5.*pow(nu1, 2) - 5.*pow(nu1, 3) + pow(nu1, 4));
      double complex out         = numerator/denominator;

      return out;
}

/**
 * Analytic expression of M_cIG2m matrix
 *
 * @param nu1         Input: FFTLog coeffcients
 * @return value of M_cIG2m
 */

double complex M_cIG2m(double complex nu1)
{
      double complex numerator   = -3. * ctan(nu1 * M_PI);
      double complex denominator = 112. * nu1 * (2. - 3.*nu1 + pow(nu1, 2) * M_PI;
      double complex out         = numerator/denominator;

      return out;
}

/*
      _________________ 2-nd Moment _______________________
 */

/**
 * Analytic expression of M_cIF2m
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_ccIF2m
 */

double complex M_cIF2m(double complex nu1, double complex nu2)
{
      double complex numerator   = (-2. + nu1 + nu2) *  (-3. + 2. * (nu1 + nu2)) * (-1. + 2. * (nu1 + nu2)) * (3. + 7. * (nu1 + nu2));
      double complex denominator = 42. * nu1 * nu2 * (1. + nu2) * (-1. + 2.*nu2);
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}
 
 /**
 * Analytic expression of M_cIG2m2 matrix
 *
 * @param nu1         Input: FFTLog coeffcients
 * @return value of M_cIG2m
 */

double complex M_cIG2m2(double complex nu1)
{
      double complex numerator   = ctan(nu1 * M_PI);
      double complex denominator = 28. * nu1 * (2. - 3.*nu1 + pow(nu1, 2) * M_PI;
      double complex out         = - numerator/denominator;

      return out;
}

/**
 * Analytic expression of M_cPPm
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cPPm
 */

double complex M_cIPPm(double complex nu1, double complex nu2)
{
      double complex numerator   = (-2. + nu1 + nu2) *  (-3. + 2. * (nu1 + nu2)) ;
      double complex denominator = 3. * nu2 * (-1. + 2.*nu2);
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}

/**
 * Analytic expression of M_cSm
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cSm
 */

double complex M_cISm(double complex nu1, double complex nu2)
{
      double complex numerator   = (-2. + nu1 + nu2) *  (-3. + 2. * (nu1 + nu2)) * (-1. + 2. *(nu1 + nu2)) ;
      double complex denominator = 3. * nu1 * nu2 * (1. + nu2) * (-1. + 2.*nu2);
      double complex out         = - numerator/denominator * Ifunc(nu1, nu2);

      return out;
}
 
 /**
 * Analytic expression of M_cIG2G2 matrix
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cIG2G2
 */

double complex M_cIG2G2(double complex nu1, double complex nu2)
{
      double complex numerator1  = (-3. + 2. * (nu1 + nu2));
      double complex numerator2  = (-1. + 2. * (nu1 + nu2));
      double complex numerator3  = (50. + 98. * pow(nu1, 3) * nu2 - (9. + 35.*nu2) * nu2 + 7.*pow(n1, 2)*(-5. + 2. * nu2*(-9. + 14.*n2)) + nu1*(-9. + 2.*nu2*(-33. + 7.*nu2*(-9. + 7.*n2)))
      double complex denominator = 196. * nu1 * (1. + nu1) * (-1. + 2.*nu1) * nu2 * (1. + nu2) * (-1. + 2.*nu2);
      double complex out         = numerator1 * numerator2 * numerator3 /denominator * Ifunc(nu1, nu2);

      return out;
}

 /**
 * Analytic expression of M_cIG2m3 matrix
 *
 * @param nu1         Input: FFTLog coeffcients
 * @return value of M_cIG2m3
 */

double complex M_cIG2m3(double complex nu1)
{
      double complex numerator   = -3. * ctan(nu1 * M_PI);
      double complex denominator = 112. * nu1 * (2. - 3.*nu1 + pow(nu1, 2) * M_PI;
      double complex out         = numerator/denominator;

      return out;
}

/**
 * Analytic expression of M_cIPPm2
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cPPm2
 */

double complex M_cIPPm2(double complex nu1, double complex nu2)
{
      double complex numerator   = (nu1 - nu2) *  (-2. + nu1 + nu2) * (-3. + 2. *(nu1 + nu2)) * (-1. + 2. * (nu1 + nu2));
      double complex denominator = 3. * nu1 * nu2 * (-1. + 2.*nu1) * (-1. + 2.*nu2);
      double complex out         = - numerator/denominator * Ifunc(nu1, nu2);

      return out;
}

/*
      _________________ 3-rd Moment _______________________
 */

 /**
 * Analytic expression of M_cIG2m4
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cG2m4
 */

double complex M_cIG2m4(double complex nu1, double complex nu2)
{
      double complex numerator   = (-2. + nu1 + nu2) *  (-3. + 2. * (nu1 + nu2)) * (-1. + 2. *(nu1 + nu2)) * (-1. + 7. *(nu1 + nu2));
      double complex denominator = 42. * nu1 * nu2 * (1. + nu2) * (-1. + 2.*nu2);
      double complex out         = - numerator/denominator * Ifunc(nu1, nu2);

      return out;
}

 /**
 * Analytic expression of M_cIPPm3
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cPPm3
 */

double complex M_cIPPm3(double complex nu1, double complex nu2)
{
      double complex numerator   = 2. * (-2. + nu1 + nu2) *  (-3. + 2. * (nu1 + nu2)) * (-1. nu1 + nu2) * (-1. + 2. *(nu1 + nu2));
      double complex denominator = 3. * nu1 * nu2 * (-1. + 2.*nu2) * (-1. + 2.*nu2);
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}

/*
      _________________ 4-th Moment _______________________
 */

 /**
 * Analytic expression of M_cIPPm41
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cPPm41
 */

double complex M_cIPPm41(double complex nu1, double complex nu2)
{
      double complex numerator   = 2. * (-2. + nu1 + nu2) *  (-3. + 2. * (nu1 + nu2)) * (-1. nu1 + nu2) * (-1. + 2. *(nu1 + nu2));
      double complex denominator = 3. * nu1 * nu2 * (-1. + 2.*nu2) * (-1. + 2.*nu2);
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}

 /**
 * Analytic expression of M_cIPPm42
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cPPm42
 */

double complex M_cIPPm42(double complex nu1, double complex nu2)
{
      double complex numerator   = 2. * (-2. + nu1 + nu2) *  (-3. + 2. * (nu1 + nu2));
      double complex denominator = 5. * nu2 * (-1. + 2.*nu2);
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}