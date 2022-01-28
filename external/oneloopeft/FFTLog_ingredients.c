

/** @file FFTLog_ingredients.c Documented ingredients for FFT-Log based 1loop integrals of matter and galaxy/halo power spectra in perturbation theory 
 * 
 * Azadeh Moradinezhad Dizgah, June 16th 2021
 *
 * This module computes the three functions that are needed for evaluating the real-space loop integrals using FFTLog technique
 * These functions are called by FFTLog_matter_real.c and and FFTLog_galaxy_real.c modules
 *
 * The FFTLog algorithm for loop integrals closely follows Ref. arXiv:1708.08130 by Simonovic et al. After computing the FFT coefficents of matter power spectrum, 
 * sampled in logarithmic scale, the algorithm involves re-casting the integrals into a form that is analytically calculable (Matrices M_ij, 
 * which can be written in terms of ratios of Gamma functions) and finally performing vactor-matrix-vector or matrix-vector multiplications, which are performed
 * in FFTLog_matter_real.c and and FFTLog_galaxy_real.c modules.
 *
 * An important feauture of fast computation of loop integrals is that all the cosmology-dependance of the loop integrals is captured by the FFT coeffcients
 * which have ~NlogN complexity. The matrices involving the Gamma functions, are computed only once, and at each of MCMC varying cosmological parmaters, 
 * these coeffcients are evaluated for all k-values at once and then vector-matrix-vector multiplication are computed.
 *
 * The FFT coeffcients are computed using FFTW package, while the vector-matrix-vector computations are performed using blas library implemented in gsl.
 * The analytic formulas for M_ij matrices are computed in Mathematica using a modified version of the publicaly available notebook by Simpnovic. 
 *
 * In summary, the following functions can be called from other modules:
 * -# Ifunc()
 * -# Gamma()
 * -# FFT_compute_coeff()
 */


#include "header.h"


/**
 * Analytic expression of the integrals of the form
 * I(nu1,nu2) = k^(-3 + 2 nu12) * int_q 1/(q^(2 nu1) kmq^(2 nu2) = 
 *
 * @param nu1         Input: power of q 
 * @param nu2         Input: power of kmq
 * @return value of I(nu1,nu2)
 */

double complex Ifunc(double complex nu1, double complex nu2)
{
      double complex nu12 = nu1 + nu2;
      double complex numerator    = Gamma(3./2.-nu1) * Gamma(3./2.-nu2) * Gamma(nu12-3./2.);
      double complex denominator  = Gamma(nu1) * Gamma(nu2) * Gamma(3.-nu12);
      // double complex out          = 1./(8. * pow(M_PI, 3./2.)) * numerator/denominator; 
      double complex out          = pow(M_PI, 3./2.) * numerator/denominator; 

      return out;
}



/**
 * Gamma function
 *
 * @param z         Input: argument of Gamma function
 * @return value of Gamma(z)
 */

double complex Gamma(double complex z)
{
      double complex  out=0.;
      gsl_sf_result arg;
      gsl_sf_result lnr;

      gsl_sf_lngamma_complex_e(creal(z), cimag(z), &lnr, &arg);
      out = exp(lnr.val) * cexp(_Complex_I*arg.val);

      return out;
}

/**
 * Compute the FFT coefficients and frequencies of the matter power spectrum 
 * 
 * First the matter power spectrum computed by class is extrapolated for smaller values of k to make 
 * sure the fft is well behaved at the largest scales of interest by extending the edges to scales larger 
 * than the kmin we are interested in MCMC.
 * 
 * Computing fourier coeffcients of the power spectrum using fftw package requires three steps: 
 * creating a ffw_plan, executing it with fftw_execute(plan) and destroying it with fftw_destroy_plan(plan)
 *
 * @param block        Input: cosmosis datablock
 * @param z            Input: redshift of the pk_m
 * @param IR_switch    Input: WIR or NOIR. This is a switch for whether or not using IR resummed pk_m
 * @param SPLIT        Input: GFILTER or BSPLINE or DST. Method for performing wiggle nowiggle split of pk_m.
 * @param biased_etam  Output: biased_etam = b + i etam output from fft decomposition of pk_m
 * @param cmsym        Output: fourier coeeficents of pk_m
 * @return void
 */

void FFT_compute_coeff(struct background * pba,
                       struct primordial * ppm,
                       struct fourier * pfo,
                       double z, 
                       struct fft_struct *fft_input,
                       long SPLIT, 
                       long hm_switch)
{
      int i, j, m;  
      int Nmax        = *fft_input->nfft;
      int Nmaxf       = *fft_input->nfft + 1;
      double Nmaxd    = (double)Nmax;

      // choosing a given bias for Matter or Galaxy/Halo calculations
      if(hm_switch == MATTER){
            double fft_bias = *fft_input->fft_bias_m;
      }
      else{
            double fft_bias = *fft_input->fft_bias_g;
      }

      double k_min    = *fft_input->kmin_fft;
      double kmax_fft = FFT_kmax_Brent_solver(pba, ppm, pfo, z, kmin_fft, fft_bias);
      double Delta    = log(kmax_fft/kmin_fft)/(Nmaxd-1); 

      double *k      = make_1Darray(Nmax);
      double *pkz    = make_1Darray(Nmax);
      double *pk_bin = make_1Darray(Nmax);
      double *window = make_1Darray(Nmax);

      int mleft     = - 0.75 * Nmax/2;
      int mright    =  0.75 * Nmax/2;
      double kleft  = k[mleft];
      double kright = k[mright];

      fftw_complex *biased_etam;
      fftw_complex *cmsym;

      for(int n=0;n<Nmax;n++){
            k[n] = kmin_fft * exp(Delta * n);
            window[n] = 1;
            // window[n] = FFT_window(k[n], kmin_fft, kmax_fft, kleft, kright);
      }
      

      for(int i=0;i<Nmax;i++){
            pk_bin[i] = pm_IR_LO(pba, ppm, pfo, k[i], z, SPLIT);
            // printf("all pk %12.6e %12.6e \n", k[i],pk_bin[i]);
      }      


      /*
       * Now perform fftw to compute the Fourier coeffcients
       */

      fftw_complex *input   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nmax);
      fftw_complex *output  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nmax);

      /*
      * etam = m * 2 pi/(Nmax Delta) with m= -N/2, -N/2+1, ...,N/2, N/2+1
      */
      double js[Nmax+1];
      double etam[Nmax+1];
      int index_c  = 0;
      int index_kd = 0;

      for (index_c =0; index_c < Nmax + 1 ; index_c ++){
         js[index_c]          = (double)index_c - Nmaxd/2.;
         etam[index_c]        = 2. * M_PI *  js[index_c] / (Nmaxd * Delta);
         biased_etam[index_c] = fft_bias + _Complex_I * etam[index_c];
      } 

      for (index_kd = 0; index_kd < Nmax; index_kd++)
         input[index_kd] = pk_bin[index_kd] * exp(- (double)index_kd * fft_bias * Delta);

      fftw_plan my_plan = fftw_plan_dft_1d(Nmax, input, output, FFTW_FORWARD, FFTW_ESTIMATE);

      fftw_execute(my_plan);


      /*
       * construct the rescaled fourier coeeficents
       */

      for (index_c=0; index_c < Nmax+1; index_c++){
         if (index_c < Nmax/2) {
             cmsym[index_c]= window[index_c] * cpow(kmin_fft,-biased_etam[index_c]) *(creal(output[Nmax/2 - index_c]) - _Complex_I * cimag(output[Nmax/2 - index_c]))/Nmaxd;
         }
         else{
             cmsym[index_c]= window[index_c] * cpow(kmin_fft,-biased_etam[index_c]) *(creal(output[index_c - Nmax/2]) + _Complex_I * cimag(output[index_c - Nmax/2]))/Nmaxd;
         }
      }

      cmsym[0]    = cmsym[0]/2.;
      cmsym[Nmax] = cmsym[Nmax]/2.;

      // for(i=0; i < Nmax+1; i++)
      //       printf("%d %12.6e %12.6e %12.6e %12.6e \n",i , creal(cmsym[i]), cimag(cmsym[i]), creal(biased_etam[i]), cimag(biased_etam[i]));
      
      // saving results in fft_struct for a giving Matter or Galaxy/Halo calculation
      if(hm_switch == MATTER){
            fft_input->etam_m  = biased_etam;
            fft_input->cmsym_m = cmsym;
      }
      else{
            fft_input->etam_g  = biased_etam;
            fft_input->cmsym_g = cmsym;
      }

      fftw_destroy_plan(my_plan);
      fftw_free(input);
      fftw_free(output);

      free(k);
      free(pkz);
      free(pk_bin);
      free(window);

      // printf("suceessfully computed FFT coeffcients\n");

      return ;
}

/**
 * Tapering window function to reduce the edge effect giving rise to ringing (See section of C of FastPT
 * @param k        Input: kvalue
 * return value of window function 
 */
double FFT_window(double k, double kmin, double kmax, double kleft, double kright)
{

  double out= 0.;
  if(k<kleft)
      out = (k - kmin)/(kleft-kmin) - 1./(2.*M_PI) * sin(2.*M_PI * (k - kmin)/(kleft-kmin));
  else if ( k>kleft && k<kright)
      out = 1.;
  else if (k>kright)
      out = (kmax - k)/(kmax-kright) - 1./(2.*M_PI) * sin(2.*M_PI * (kmax - k)/(kmax-kright));
  return out;
}

/**
 * The following two functions are used to compute the kmax value to be used in computaing the FFT coefficients. The kmax is set such that the low and high end 
 * of the matter power spectrum have the same value. 
 * 
 * The kmax value is determined by solving an algebraic equations using gsl_Brent solver. This routine requires a function for the equation 
 * This is important for loweing the amount of aliasing.FFT_kmax_Brent() and a function that solves the equation using Brent solver FFT_kmax_Brent_solver()
 *
 * parameters of FFT_kmax_Brent
 * @param kmax         Input: value of kmax to be determined by Brent solver
 * @param params       Input: integration parameters
 * @ return double    (the equation to be avaluated)    
 *
 * parameters of FFT_kmax_Brent_solver
 * @param block        Input: cosmosis data block
 * @param z            Input: redshift at which to evaluate the kmax
 * @param kmax         Input: kmin of fft grid
 * @param fft_bias     Input: number of points for fft grid
 * @return value of kmax     
 */
double FFT_kmax_Brent(double kmax, void *params)
{
      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)params);

      struct background *pba = pij.pba;
      struct primordial *ppm = pij.ppm;
      struct fourier *pfo    = pij.pfo;
      double z               = pij.p4;
      double kmin_fft        = pij.p5;
      double bias_fft        = pij.p6;

      double f = pow(kmax,-bias_fft) * Pk_dlnPk(pba, ppm, pfo, kmax, z, LPOWER) - pow(kmin_fft, - bias_fft) * Pk_dlnPk(pba, ppm, pfo, kmin_fft, z, LPOWER);
      
      return f;   
}

double FFT_kmax_Brent_solver(struct background * pba,
                             struct primordial * ppm,
                             struct fourier * pfo, 
                             double z, 
                             double kmin_fft, 
                             double fft_bias)
{
      int status;

      int iter = 0, max_iter = 100000;
      const gsl_root_fsolver_type *T;
      gsl_root_fsolver *s;
      double r = 0;
      double k_lo = 0.01, k_hi = 4000.;
      gsl_function F;

      struct integrand_parameters2 par; 

      par.ppm = ppm;
      par.pba = pba;
      par.pfo = pfo;
      par.p4  = z;
      par.p5  = kmin_fft;
      par.p6  = fft_bias; 

      F.function = &FFT_kmax_Brent;
      F.params   = &par;

      T = gsl_root_fsolver_brent;
      s = gsl_root_fsolver_alloc (T);
      gsl_root_fsolver_set (s, &F, k_lo, k_hi);

      do
      {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r      = gsl_root_fsolver_root (s);
        k_lo   = gsl_root_fsolver_x_lower (s);
        k_hi   = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (k_lo, k_hi, 0., 1e-5);
      }
      while (status == GSL_CONTINUE && iter < max_iter);

      gsl_root_fsolver_free (s);
      
      return r;

      }

