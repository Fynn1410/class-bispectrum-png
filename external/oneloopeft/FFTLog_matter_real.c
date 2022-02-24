/** @file FFTLog_matter_real.c Documented FFT-Log based 1loop integrals of matter power spectrum in perturbation theory 
 * 
 * Azadeh Moradinezhad Dizgah, June 16th 2021
 *
 * This module performs fast computation of the integrals appearing in the expression of 1loop matter and galaxy power sprtcurm.
 * The computation can be performed either in real or redshift-space. IR-resummation and EFT counter terms are included. 
 * The integrals are computed using FFTLog techniques.
 *
 * The algorithm closely follows Ref. arXiv:1708.08130 by Simonovic et al. After computing the FFT coefficents of matter power spectrum, 
 * sampled in logarithmic scale, the algorithm involves re-casting the integrals into a form that is analytically calculable (Matrices M_ij), 
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
 * For matter loops, if choosing 512 points *which takes ~ 0.45 seconds, I do get subpercent descrepencies between direct numeric integration of the loops.
 * The best value seem to be around 600 points, which takes ~ 0.6 seconds and is in exquisit agreement with direct integration. 
 *
 * In summary, the following functions can be called from other modules:
 * -# P22()
 * -# P13()
 * -# M22()
 * -# M13()
 */

#include "header.h"


/**
 * Compute the P22 contribution to non-linear matter power spectrum given the FFTLog coefficents and frequencies
 *
 * @param fft_struct    Input: structure containing fft coefficents and params
 * @param k             Input: wavenumber in unit of h/Mpc. 
 * @param z             Input: redshift
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the pipeline
 * @return value of P22 term in unit of (Mpc/h)^3 as a function of redshift and wavnumber
 */

double P22(struct fft_struct *fft_input, double k, double z, int cleanup)
{
      int Nmax     = fft_input -> nfft;

      static gsl_matrix_complex *M22_mat;
      double complex nu1, nu2;
      gsl_complex cmk_gsl, out_complex;

      gsl_vector_complex *cmk_vec   = gsl_vector_complex_alloc(Nmax+1);
      gsl_vector_complex *matvec    = gsl_vector_complex_alloc(Nmax+1);

      int i, j;
      static int first = 1;
      if(first==1)
      {
         gsl_complex M22_gsl;

         M22_mat = gsl_matrix_complex_alloc(Nmax+1,Nmax+1);

         for (i=0; i<Nmax+1; i++){
               for (j=i; j<Nmax+1; j++){                        
                     nu1 = -0.5 * fft_input->etam_m[i];
                     nu2 = -0.5 * fft_input->etam_m[j];
                     //printf("M22 %d %d %12.6e %12.6e %12.6e  \n", i, j, k,  creal(M22(nu1,nu2)), cimag(M22(nu1,nu2)));
                     
                     GSL_SET_COMPLEX(&M22_gsl,creal(M22(nu1,nu2)),cimag(M22(nu1,nu2)));
                     gsl_matrix_complex_set(M22_mat, i, j, M22_gsl);
                     gsl_matrix_complex_set(M22_mat, j, i, M22_gsl); 
               }
         }
         first=0;
      }   

      for (i=0; i<Nmax+1; i++){ 
            nu1 = -0.5 * fft_input->etam_m[i];                       
            GSL_SET_COMPLEX(&cmk_gsl, creal(fft_input->cmsym_m[i] * cpow(k,-2.*nu1)), cimag(fft_input->cmsym_m[i] * cpow(k,-2.*nu1)));
            gsl_vector_complex_set(cmk_vec, i, cmk_gsl);

      }

      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, M22_mat, cmk_vec, GSL_COMPLEX_ZERO, matvec);
           
      gsl_blas_zdotu(cmk_vec, matvec, &out_complex);
      double out = pow(k,3.) * GSL_REAL(out_complex);

      gsl_vector_complex_free(cmk_vec);
      gsl_vector_complex_free(matvec);

      if(cleanup == 1)
            gsl_matrix_complex_free(M22_mat);
      

      // printf("suceessfully computed P22 using FFTLog\n");

      return out;

}


/**
 * Compute the P13 contribution to non-linear matter power spectrum given the FFTLog coefficents and frequencies
 *
 * @param fft_struct    Input: structure containing fft coefficents and params
 * @param k             Input: wavenumber in unit of h/Mpc. 
 * @param z             Input: redshift
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the pipeline
 * @return value of P13 term in unit of (Mpc/h)^3 as a function of redshift and wavnumber
 */

double P13(struct fft_struct *fft_input, double k, double z, int cleanup)
{
      fftw_complex M22_elems;

      int Nmax     = fft_input -> nfft;

      static gsl_vector_complex *M13_vec;
      gsl_vector_complex *cmk_vec  = gsl_vector_complex_alloc(Nmax+1);

      gsl_complex cmk_gsl, out_complex;
      double complex nu1;

      static int first = 1;
      if(first==1)
      {
         gsl_complex M13_gsl;
         M13_vec  = gsl_vector_complex_alloc(Nmax+1);

         for (int i=0; i<Nmax+1; i++){
               nu1 = -0.5 * fft_input->etam_m[i];
               GSL_SET_COMPLEX(&M13_gsl,creal(M13(nu1)),cimag(M13(nu1)));
               //printf("M13 %d %12.6e %12.6e %12.6e  %12.6e %12.6e \n", i, k, creal(nu1), cimag(nu1), GSL_REAL(M13_gsl), GSL_IMAG(M13_gsl));

               gsl_vector_complex_set(M13_vec, i, M13_gsl);
         }
         first=0;
      }         
      for (int i=0; i<Nmax+1; i++){
             nu1 = -0.5 * fft_input->etam_m[i];
             GSL_SET_COMPLEX(&cmk_gsl, creal(fft_input->cmsym_m[i] * cpow(k,-2.*nu1)), cimag(fft_input->cmsym_m[i] * cpow(k,-2.*nu1)));
             // printf("cmk13 %d %12.6e %12.6e %12.6e  \n", i, k,  GSL_REAL(cmk_gsl), GSL_IMAG(cmk_gsl));

            gsl_vector_complex_set(cmk_vec, i, cmk_gsl);
      }
      gsl_blas_zdotu(cmk_vec, M13_vec, &out_complex);
      double out = pow(k,3.) * GSL_REAL(out_complex);


      gsl_vector_complex_free(cmk_vec);

      if(cleanup == 1)
         gsl_vector_complex_free(M13_vec);

      // printf("suceessfully computed P13 using FFTLog\n");

      return out;
}



/**
 * Analytic expression of M22 matrix
 *
 * @param nu1         Input: FFTLog coeffcients
 * @param nu2         Input: FFTLog frequency exponents
 * @return value of M22
 */

double complex M22(double complex nu1, double complex nu2)
{
      double complex numerator   = ((-3. + 2. * nu1 + 2. * nu2)*(-1. + 2. * nu1 + 2. * nu2) * (58. + 98. * cpow(nu1, 3.) * nu2 \
                                 + (3. - 91. * nu2) * nu2 + 7. * cpow(nu1, 2.) * (- 13. - 2. * nu2 + 28. * cpow(nu2, 2.))
                                 + nu1 * (3. + 2. * nu2 * (- 73. + 7. * nu2 * (- 1. + 7. * nu2)))));
      double complex denominator = (196. * nu1 * (1. + nu1) * (-1. + 2. * nu1) * nu2 * (1. + nu2) * (- 1. + 2. * nu2));
      double complex out         = numerator/denominator * Ifunc(nu1, nu2);

      return out;
}



/**
 * Analytic expression of M13 matrix
 *
 * @param nu1         Input: power of q
 * @return value of M13
 */

double complex M13(double complex nu1)
{
      double complex numerator   = (1. + 9. * nu1) * ctan(nu1 * M_PI);
      double complex denominator = (112. *  M_PI * nu1 * (-6. + 5. * nu1 + 5. * cpow(nu1, 2.) - 5. * cpow(nu1, 3.) + cpow(nu1, 4.)));
      // double complex numerator   = (1. + 9. * nu1) * pow(M_PI,2.) * ctan(nu1 * M_PI);
      // double complex denominator = (14. * nu1 * (-6. + 5. * nu1 + 5. * cpow(nu1, 2.) - 5. * cpow(nu1, 3.) + cpow(nu1, 4.)));
      double complex out          = numerator/denominator;

      return out;

}

