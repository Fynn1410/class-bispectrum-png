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
 * -# M_cIdelta2()
 * -# M_cIcG2()
 * -# M_cIdelta2delta2()
 * -# M_cIcG2cG2()
 * -# M_cIdelta2cG2()
 * -# M_cFcG2()
 */

#include "header.h"


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

void pgloops_nonpropag(struct fft_struct  *fft_input, double k, double z, int cleanup, double *pg_loops)
{
      int Nmax     = fft_input -> nfft;

      static double P_b2b2_ls = 0.;
      static gsl_matrix_complex *M_cIdelta2_mat;
      static gsl_matrix_complex *M_cIcG2_mat;
      static gsl_matrix_complex *M_cIdelta2delta2_mat;
      static gsl_matrix_complex *M_cIcG2cG2_mat;
      static gsl_matrix_complex *M_cIdelta2cG2_mat;

      gsl_complex out_cIdelta2_complex;
      gsl_complex out_cIcG2_complex;
      gsl_complex out_cIdelta2delta2_complex;
      gsl_complex out_cIcG2cG2_complex;
      gsl_complex out_cIdelta2cG2_complex;

      gsl_vector_complex *matvec_cIdelta2        = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *matvec_cIcG2           = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *matvec_cIdelta2delta2  = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *matvec_cIcG2cG2        = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *matvec_cIdelta2cG2     = gsl_vector_complex_alloc(Nmax+1);

      double complex nu1, nu2;
      gsl_complex cmk_gsl;
      gsl_complex cmkmin_gsl;

      gsl_vector_complex *cmk_vec                = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *cmkmin_vec             = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *matvec_cIdelta2delta20 = gsl_vector_complex_alloc(Nmax+1);
      gsl_complex out_cIdelta2delta20_complex;

      int i, j;
      static int first = 1;
      if(first==1)
      {
         gsl_complex M_cIdelta2_gsl;
         gsl_complex M_cIcG2_gsl;
         gsl_complex M_cIdelta2delta2_gsl;
         gsl_complex M_cIcG2cG2_gsl;
         gsl_complex M_cIdelta2cG2_gsl;

         M_cIdelta2_mat       = gsl_matrix_complex_alloc(Nmax+1,Nmax+1);
         M_cIcG2_mat          = gsl_matrix_complex_alloc(Nmax+1,Nmax+1);
         M_cIdelta2delta2_mat = gsl_matrix_complex_alloc(Nmax+1,Nmax+1);
         M_cIcG2cG2_mat       = gsl_matrix_complex_alloc(Nmax+1,Nmax+1);
         M_cIdelta2cG2_mat    = gsl_matrix_complex_alloc(Nmax+1,Nmax+1);


         for (i=0; i<Nmax+1; i++){
               for (j=i; j<Nmax+1; j++){                        
                     nu1 = -0.5 * fft_input->etam_g[i];
                     nu2 = -0.5 * fft_input->etam_g[j];                     

                     GSL_SET_COMPLEX(&M_cIdelta2_gsl,creal(M_cIdelta2(nu1,nu2)),cimag(M_cIdelta2(nu1,nu2)));
                     gsl_matrix_complex_set(M_cIdelta2_mat, i, j, M_cIdelta2_gsl);
                     gsl_matrix_complex_set(M_cIdelta2_mat, j, i, M_cIdelta2_gsl); 

                     GSL_SET_COMPLEX(&M_cIcG2_gsl,creal(M_cIcG2(nu1,nu2)),cimag(M_cIcG2(nu1,nu2)));
                     gsl_matrix_complex_set(M_cIcG2_mat, i, j, M_cIcG2_gsl);
                     gsl_matrix_complex_set(M_cIcG2_mat, j, i, M_cIcG2_gsl); 

                     GSL_SET_COMPLEX(&M_cIdelta2delta2_gsl,creal(M_cIdelta2delta2(nu1,nu2)),cimag(M_cIdelta2delta2(nu1,nu2)));
                     gsl_matrix_complex_set(M_cIdelta2delta2_mat, i, j, M_cIdelta2delta2_gsl);
                     gsl_matrix_complex_set(M_cIdelta2delta2_mat, j, i, M_cIdelta2delta2_gsl); 

                     GSL_SET_COMPLEX(&M_cIcG2cG2_gsl,creal(M_cIcG2cG2(nu1,nu2)),cimag(M_cIcG2cG2(nu1,nu2)));
                     gsl_matrix_complex_set(M_cIcG2cG2_mat, i, j, M_cIcG2cG2_gsl);
                     gsl_matrix_complex_set(M_cIcG2cG2_mat, j, i, M_cIcG2cG2_gsl); 

                     GSL_SET_COMPLEX(&M_cIdelta2cG2_gsl,creal(M_cIdelta2cG2(nu1,nu2)),cimag(M_cIdelta2cG2(nu1,nu2)));
                     gsl_matrix_complex_set(M_cIdelta2cG2_mat, i, j, M_cIdelta2cG2_gsl);
                     gsl_matrix_complex_set(M_cIdelta2cG2_mat, j, i, M_cIdelta2cG2_gsl); 

               }
         }
         first =0;
      }


      for (i=0; i<Nmax+1; i++){ 
         nu1 = -0.5 * fft_input->etam_g[i];                       
         GSL_SET_COMPLEX(&cmkmin_gsl, creal(fft_input->cmsym_g[i] * cpow(fft_input->kmin_fft,-2.*nu1)), cimag(fft_input->cmsym_g[i] * cpow(fft_input->kmin_fft,-2.*nu1)));
         gsl_vector_complex_set(cmkmin_vec, i, cmkmin_gsl);
      }

      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M_cIdelta2delta2_mat, cmkmin_vec, GSL_COMPLEX_ZERO, matvec_cIdelta2delta20);
      gsl_blas_zdotu(cmkmin_vec, matvec_cIdelta2delta20, &out_cIdelta2delta20_complex);
      
      P_b2b2_ls = pow(fft_input->kmin_fft,3.) * GSL_REAL(out_cIdelta2delta20_complex);
      // printf("Pb2b2 %12.6e \n", P_b2b2_ls);

      gsl_vector_complex_free(cmkmin_vec);
      gsl_vector_complex_free(matvec_cIdelta2delta20);


      for (i=0; i<Nmax+1; i++){ 
            nu1 = -0.5 * fft_input->etam_g[i];                       
            GSL_SET_COMPLEX(&cmk_gsl, creal(fft_input->cmsym_g[i] * cpow(k,-2.*nu1)), cimag(fft_input->cmsym_g[i] * cpow(k,-2.*nu1)));
            gsl_vector_complex_set(cmk_vec, i, cmk_gsl);
      }

      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M_cIdelta2_mat,       cmk_vec,    GSL_COMPLEX_ZERO, matvec_cIdelta2);
      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M_cIcG2_mat,          cmk_vec,    GSL_COMPLEX_ZERO, matvec_cIcG2);
      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M_cIdelta2delta2_mat, cmk_vec,    GSL_COMPLEX_ZERO, matvec_cIdelta2delta2);
      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M_cIcG2cG2_mat,       cmk_vec,    GSL_COMPLEX_ZERO, matvec_cIcG2cG2);
      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M_cIdelta2cG2_mat,    cmk_vec,    GSL_COMPLEX_ZERO, matvec_cIdelta2cG2);

      gsl_blas_zdotu(cmk_vec, matvec_cIdelta2,        &out_cIdelta2_complex);
      gsl_blas_zdotu(cmk_vec, matvec_cIcG2,           &out_cIcG2_complex);
      gsl_blas_zdotu(cmk_vec, matvec_cIdelta2delta2,  &out_cIdelta2delta2_complex);
      gsl_blas_zdotu(cmk_vec, matvec_cIcG2cG2,        &out_cIcG2cG2_complex);
      gsl_blas_zdotu(cmk_vec, matvec_cIdelta2cG2,     &out_cIdelta2cG2_complex);

      pg_loops[0] = pow(k,3.) * GSL_REAL(out_cIdelta2_complex);
      pg_loops[1] = pow(k,3.) * GSL_REAL(out_cIcG2_complex);
      pg_loops[2] = pow(k,3.) * GSL_REAL(out_cIdelta2delta2_complex) - P_b2b2_ls;
      pg_loops[3] = pow(k,3.) * GSL_REAL(out_cIcG2cG2_complex);
      pg_loops[4] = pow(k,3.) * GSL_REAL(out_cIdelta2cG2_complex);

      gsl_vector_complex_free(cmk_vec);

      gsl_vector_complex_free(matvec_cIdelta2);
      gsl_vector_complex_free(matvec_cIcG2);
      gsl_vector_complex_free(matvec_cIdelta2delta2);
      gsl_vector_complex_free(matvec_cIcG2cG2);
      gsl_vector_complex_free(matvec_cIdelta2cG2);

      if(cleanup == 1){
            gsl_matrix_complex_free(M_cIdelta2_mat);
            gsl_matrix_complex_free(M_cIcG2_mat);
            gsl_matrix_complex_free(M_cIdelta2delta2_mat);
            gsl_matrix_complex_free(M_cIcG2cG2_mat);
            gsl_matrix_complex_free(M_cIdelta2cG2_mat);

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
 * @return value of P13 term in unit of (Mpc/h)^3 as a function of redshift and wavnumber
 */

double pgloops_propag(struct fft_struct  *fft_input, double k, double z, int cleanup)
{
      fftw_complex McFcG2_elems;

      int Nmax     = fft_input -> nfft;

      static gsl_vector_complex *M_cFcG2_vec;
      gsl_vector_complex *cmk_vec  = gsl_vector_complex_alloc(Nmax+1);

      gsl_complex cmk_gsl, out_cFcG2_complex;
      double complex nu1;

      /*
      * The M_ij matrices are cosmology-independant and need to be computed only once. So we use the first switch
      */
      static int first = 1;
      if(first == 1)
      {
         gsl_complex M_cFcG2_gsl;
         M_cFcG2_vec  = gsl_vector_complex_alloc(Nmax+1);

         for (int i=0; i<Nmax+1; i++){
               nu1 = -0.5 * fft_input->etam_g[i];
               GSL_SET_COMPLEX(&M_cFcG2_gsl,creal(M_cFcG2(nu1)),cimag(M_cFcG2(nu1)));
               // printf("M_cFcG2_gsl %d %12.6e %12.6e %12.6e  \n", i, k,  GSL_REAL(M_cFcG2_gsl), GSL_IMAG(M_cFcG2_gsl));

               gsl_vector_complex_set(M_cFcG2_vec, i, M_cFcG2_gsl);
         }
         first = 0;
      }         
      for (int i=0; i<Nmax+1; i++){
             nu1 = -0.5 * fft_input->etam_g[i];
             GSL_SET_COMPLEX(&cmk_gsl, creal(fft_input->cmsym_g[i] * cpow(k,-2.*nu1)), cimag(fft_input->cmsym_g[i] * cpow(k,-2.*nu1)));
             // printf("cmk13 %d %12.6e %12.6e %12.6e  \n", i, k,  GSL_REAL(cmk_gsl), GSL_IMAG(cmk_gsl));

            gsl_vector_complex_set(cmk_vec, i, cmk_gsl);
      }
      gsl_blas_zdotu(cmk_vec, M_cFcG2_vec, &out_cFcG2_complex);
      double out_cFcG2 = pow(k,3.) * GSL_REAL(out_cFcG2_complex);

      gsl_vector_complex_free(cmk_vec);

      if(cleanup == 1)
         gsl_vector_complex_free(M_cFcG2_vec);

      // printf("suceessfully computed pgloops_propag using FFTLog\n");

      return out_cFcG2;
}



/**
 * Analytic expression of M_cIdelta2 matrix. (Eq. 2.44 of Simonovic 2017)
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cIdelta2
 */

double complex M_cIdelta2(double complex nu1, double complex nu2)
{
      double complex numerator   = (3. - 2. * (nu1 + nu2)) * (4. - 7. * (nu1 + nu2));
      //double complex denominator = 17. * nu1 * nu2;
      double complex denominator = 14. * nu1 * nu2;
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}



/**
 * Analytic expression of M_cIcG2 matrix (Eq. 2.45 of Simonovic 2017))
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cIcG2
 */

double complex M_cIcG2(double complex nu1, double complex nu2)
{
      double complex numerator   = (3. - 2. * (nu1 + nu2)) *  (1. - 2. * (nu1 + nu2)) *  (6. + 7. * (nu1 + nu2));
      double complex denominator = 28. * nu1 * (1. + nu1) * nu2 * (1. + nu2);
      double complex out         = - numerator/denominator * Ifunc(nu1, nu2);

      return out;
}


/**
 * Analytic expression of M_cIdelta2delta2 matrix (Eq. 2.47 of Simonovic 2017)
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cIdelta2delta2
 */

double complex M_cIdelta2delta2(double complex nu1, double complex nu2)
{
      double complex out  = 2. * Ifunc(nu1, nu2);

      return out;
}

/**
 * Analytic expression of M_cIcG2cG2  (Eq. 2.48 of Simonovic 2017)
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cIcG2cG2
 */

double complex M_cIcG2cG2(double complex nu1, double complex nu2)
{
      double complex numerator   = (3. - 2. * (nu1 + nu2)) *  (1. - 2. * (nu1 + nu2));
      double complex denominator = nu1 * (1. + nu1) * nu2 * (1. +nu2);
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}


/**
 * Analytic expression of M_cIdelta2cG2  (Eq. 2.49 of Simonovic 2017)
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M_cIdelta2cG2
 */

double complex M_cIdelta2cG2(double complex nu1, double complex nu2)
{
      double complex numerator   = (3. - 2. * (nu1 + nu2));
      double complex denominator = nu1 * nu2 ;
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}

/**
 * Analytic expression of M_cFcG2 matrix (Eq. 2.46 of Simonovic 2017))
 *
 * @param nu1         Input: power of q
 * @return value of M_cFcG2
 */

double complex M_cFcG2(double complex nu1)
{
      // double complex numerator    = 15. * ctan(nu1 * M_PI);
      // double complex denominator  = 28. * M_PI * (nu1 + 1.) * nu1 * (nu1 - 1.) * (nu1 - 2.) * (nu1 - 3.);
      double complex numerator    = 30. *  pow(M_PI,2.) * ctan(nu1 * M_PI);
      double complex denominator  = 7. * (nu1 + 1.) * nu1 * (nu1 - 1.) * (nu1 - 2.) * (nu1 - 3.);
      double complex out          = - numerator/denominator;

      return out;

}
