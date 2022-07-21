

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
 * J(nu1,nu2) = k^(-3 + 2 nu12) * int_q 1/(q^(2 nu1) kmq^(2 nu2) = 
 *
 * @param nu1         Input: power of q 
 * @param nu2         Input: power of kmq
 * @return value of I(nu1,nu2)
 */

double complex J(double complex nu1, double complex nu2)
{
      double complex nu12 = nu1 + nu2;
      double complex numerator    = Gamma(3./2.-nu1) * Gamma(3./2.-nu2) * Gamma(nu12-3./2.);
      double complex denominator  = Gamma(nu1) * Gamma(nu2) * Gamma(3.-nu12);
      double complex out          = 1./(8. * pow(M_PI, 3./2.)) * numerator/denominator; 

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
 * Analytic expression of the integrals of the form
 * M1(nu1,nu2) = k^(-3 + 2 nu12) * int_q q|| /(q^(2 nu1) kmq^(2 nu2)
 *
 * @param nu1         Input: power of q 
 * @param nu2         Input: power of kmq
 * @return value of I(nu1,nu2)
 */

double complex M1(double complex nu1, double complex nu2)
{
      double complex out = (J(-1 + nu1,nu2) - J(nu1,-1 + nu2) + J(nu1,nu2))/2.; 

      return out;
}

/**
 * Analytic expression of the integrals of the form
 * M2 = k^(-3 + 2 nu12) * int_q q|| q|| /(q^(2 nu1) kmq^(2 nu2) 
 *
 * @param nu1         Input: power of q 
 * @param nu2         Input: power of kmq
 * @return value of I(nu1,nu2)
 */

double complex M2(double complex nu1, double complex nu2, double mu)
{
      double complex A2 = (-J(-2 + nu1,nu2) + 2*J(-1 + nu1,-1 + nu2) + 
     2*J(-1 + nu1,nu2) - J(nu1,-2 + nu2) + 2*J(nu1,-1 + nu2) - 
     J(nu1,nu2))/8.;
      double complex B2 = (3*J(-2 + nu1,nu2) - 6*J(-1 + nu1,-1 + nu2) + 
     2*J(-1 + nu1,nu2) + 3*J(nu1,-2 + nu2) - 
     6*J(nu1,-1 + nu2) + 3*J(nu1,nu2))/8.;
      
      double complex out = A2 + pow(mu, 2.) * B2; 

      return out;
}

/**
 * Analytic expression of the integrals of the form
 * M3 = k^(-3 + 2 nu12) * int_q q|| q|| q|| /(q^(2 nu1) kmq^(2 nu2) 
 *
 * @param nu1         Input: power of q 
 * @param nu2         Input: power of kmq
 * @return value of I(nu1,nu2)
 */

double complex M3(double complex nu1, double complex nu2, double mu)
{
      double complex A3 = (-3*(J(-3 + nu1,nu2) - 3*J(-2 + nu1,-1 + nu2) - 
       J(-2 + nu1,nu2) + 3*J(-1 + nu1,-2 + nu2) - 
       2*J(-1 + nu1,-1 + nu2) - J(-1 + nu1,nu2) - 
       J(nu1,-3 + nu2) + 3*J(nu1,-2 + nu2) - 
       3*J(nu1,-1 + nu2) + J(nu1,nu2)))/16.;
      double complex B3 = (5*J(-3 + nu1,nu2) - 15*J(-2 + nu1,-1 + nu2) + 
     3*J(-2 + nu1,nu2) + 15*J(-1 + nu1,-2 + nu2) - 
     18*J(-1 + nu1,-1 + nu2) + 3*J(-1 + nu1,nu2) - 
     5*J(nu1,-3 + nu2) + 15*J(nu1,-2 + nu2) - 
     15*J(nu1,-1 + nu2) + 5*J(nu1,nu2))/16.;
      
      double complex out = A3 + pow(mu, 2.) * B3; 

      return out;
}

/**
 * Analytic expression of the integrals of the form
 * M4 = k^(-3 + 2 nu12) * int_q q|| q|| q|| q|| /(q^(2 nu1) kmq^(2 nu2) 
 *
 * @param nu1         Input: power of q 
 * @param nu2         Input: power of kmq
 * @return value of I(nu1,nu2)
 */

double complex M4(double complex nu1, double complex nu2, double mu)
{

      double complex A4 = (3*(J(-4 + nu1,nu2) - 4*J(-3 + nu1,-1 + nu2) - 4*J(-3 + nu1,nu2) + 
       6*J(-2 + nu1,-2 + nu2) + 4*J(-2 + nu1,-1 + nu2) + 
       6*J(-2 + nu1,nu2) - 4*J(-1 + nu1,-3 + nu2) + 
       4*J(-1 + nu1,-2 + nu2) + 4*J(-1 + nu1,-1 + nu2) - 
       4*J(-1 + nu1,nu2) + J(nu1,-4 + nu2) - 4*J(nu1,-3 + nu2) + 
       6*J(nu1,-2 + nu2) - 4*J(nu1,-1 + nu2) + J(nu1,nu2)))/128.;
      double complex B4 = (-3*(5*J(-4 + nu1,nu2) - 20*J(-3 + nu1,-1 + nu2) - 4*J(-3 + nu1,nu2) + 
       30*J(-2 + nu1,-2 + nu2) - 12*J(-2 + nu1,-1 + nu2) - 
       2*J(-2 + nu1,nu2) - 20*J(-1 + nu1,-3 + nu2) + 
       36*J(-1 + nu1,-2 + nu2) - 12*J(-1 + nu1,-1 + nu2) - 
       4*J(-1 + nu1,nu2) + 5*J(nu1,-4 + nu2) - 20*J(nu1,-3 + nu2) + 
       30*J(nu1,-2 + nu2) - 20*J(nu1,-1 + nu2) + 5*J(nu1,nu2)))/64.;
      double complex C4 = (35*J(-4 + nu1,nu2) - 140*J(-3 + nu1,-1 + nu2) + 20*J(-3 + nu1,nu2) + 
      210*J(-2 + nu1,-2 + nu2) - 180*J(-2 + nu1,-1 + nu2) + 
      18*J(-2 + nu1,nu2) - 140*J(-1 + nu1,-3 + nu2) + 
      300*J(-1 + nu1,-2 + nu2) - 180*J(-1 + nu1,-1 + nu2) + 
      20*J(-1 + nu1,nu2) + 35*J(nu1,-4 + nu2) - 140*J(nu1,-3 + nu2) + 
      210*J(nu1,-2 + nu2) - 140*J(nu1,-1 + nu2) + 35*J(nu1,nu2))/128.;
      
      double complex out = A4 + pow(mu, 2.) * B4 + pow(mu, 4.) * C4; 

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
      int Nmax        = fft_input->nfft;
      int Nmaxf       = fft_input->nfft + 1;
      double Nmaxd    = (double)Nmax;

      // choosing a given bias for Matter or Galaxy/Halo calculations
      double fft_bias = 0.0;
      double kmin_fft = 0.0;
      // double kmax_fft = 0.0;
      if(hm_switch == MATTER){
            fft_bias = fft_input->fft_bias_m;
            kmin_fft = fft_input->kmin_fft_m;
            // kmax_fft = 113.4721; // Model 2 parameters for Model 1 testing
            // kmax_fft = 130.2380; // Model 1 parameters for Model 2 testing
      }
      else{
            fft_bias = fft_input->fft_bias_g;
            kmin_fft = fft_input->kmin_fft_g;
            // kmax_fft = 7854.064; // Model 2 parameters for Model 1 testing
            // kmax_fft = 1062.590; // Model 1 parameters for Model 2 testing
      }
      // fprintf(stderr, "Hi %e %e\n", fft_bias, kmin_fft);
      double kmax_fft = FFT_kmax_Brent_solver(pba, ppm, pfo, z, kmin_fft, fft_bias);
      // fprintf(stderr, "no Hi\n");
      double Delta   = log(kmax_fft/kmin_fft)/(Nmaxd-1); 
      double *k      = make_1Darray(Nmax);
      double *pkz    = make_1Darray(Nmax);
      double *pk_bin = make_1Darray(Nmax);
      double *window = make_1Darray(Nmax);
      
      int mleft     = - 0.75 * Nmax/2;
      int mright    =  0.75 * Nmax/2;
      double kleft  = k[mleft];
      double kright = k[mright];

      // if(hm_switch == MATTER){
      //       fprintf(stderr, "kmax_m = %12.6e \n", kmax_fft);
      // }
      // else{
      //       fprintf(stderr, "kmax_g = %12.6e \n", kmax_fft);
      // }

      fftw_complex *biased_etam = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(Nmax+1));
      fftw_complex *cmsym       = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(Nmax+1));

      for(int n=0;n<Nmax;n++){
            k[n] = kmin_fft * exp(Delta * n);
            window[n] = 1;
            // window[n] = FFT_window(k[n], kmin_fft, kmax_fft, kleft, kright);
      }
      

      for(int i=0;i<Nmax;i++){
            pk_bin[i] = pm_IR_LO(pba, ppm, pfo, k[i], z, SPLIT);
            // pk_bin[i] = Pk_dlnPk(pba, ppm, pfo, k[i], z, LPOWER);
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

      // FILE *fpa;
      // char file_name[50];
      // if(hm_switch == MATTER){
      //       sprintf(file_name, "FFTcoeff_m2_128.txt");
      // }
      // else{
      //       sprintf(file_name, "FFTcoeff_g2_128.txt");
      // }

      // fpa = fopen(file_name, "w");
      // fprintf(fpa, "kmin: %12.6e; kmax: %12.6e; fftbias: %12.6e; nfft: %d \n \n",\
      //               kmin_fft, kmax_fft, fft_bias, Nmax+1);
      // for (index_c=0; index_c < Nmax+1; index_c++){
      //       fprintf(fpa, "%e %e %e %e\n", creal(biased_etam[index_c]), cimag(biased_etam[index_c]), creal(cmsym[index_c]), cimag(cmsym[index_c]));
      // }
      // fclose(fpa);

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

      // fprintf(stderr, "k = %e, bias = %e, P(k) = %e, k_min**bias * P(k_min) = %e, k**bias * P(k) = %e\n", kmax, bias_fft, Pk_dlnPk(pba, ppm, pfo, kmax, z, LPOWER), Pk_dlnPk(pba, ppm, pfo, kmin_fft, z, LPOWER)*pow(kmin_fft, -bias_fft), Pk_dlnPk(pba, ppm, pfo, kmax, z, LPOWER)*pow(kmax, -bias_fft));
      double f = pow(kmax, -bias_fft) * Pk_dlnPk(pba, ppm, pfo, kmax, z, LPOWER) - pow(kmin_fft, -bias_fft) * Pk_dlnPk(pba, ppm, pfo, kmin_fft, z, LPOWER);
      
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
      double k_lo = 0.01, k_hi = 20000.;
      gsl_function F;

      struct integrand_parameters2 par; 

      par.ppm = ppm;
      par.pba = pba;
      par.pfo = pfo;
      par.p4  = z;
      par.p5  = kmin_fft;
      par.p6  = fft_bias; 
      // fprintf(stderr, "z = %e, k_min = %e, bias = %e, P(k_min) = %e, k_min**bias * P(k_min) = %e\n", z, kmin_fft, fft_bias, Pk_dlnPk(pba, ppm, pfo, kmin_fft, z, LPOWER), Pk_dlnPk(pba, ppm, pfo, kmin_fft, z, LPOWER)*pow(kmin_fft, -fft_bias));

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

int FFTLog_rsd_init(struct background *pba, struct primordial *ppm, struct fourier *pfo, double z){

      pfo -> pk_halo_nl     = (struct oneloop_fftlog_halo_real *)malloc(sizeof(struct oneloop_fftlog_halo_real));
      pfo -> pk_halo_rsd_nl = (struct oneloop_fftlog_halo_rsd *)malloc(sizeof(struct oneloop_fftlog_halo_rsd));

      pfo -> pk_halo_nl -> pmm = make_1Darray(pfo->k_size);
      pfo -> pk_halo_nl -> pb1b2 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_nl -> pb1bg2 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_nl -> pb22 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_nl -> pbg22 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_nl -> pb2bg2 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_nl -> pb1b3nl = make_1Darray(pfo->k_size);    
    
      pfo -> pk_halo_rsd_nl -> P_mm = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> I2200 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> Idelta200 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> IG200 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> Idelta2delta200 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> IG2G200 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> Idelta2G200 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> I1300 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> FG200 = make_1Darray(pfo->k_size);

      pfo -> pk_halo_rsd_nl -> I2201 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> Idelta201 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> IG201 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> FG201 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> J21101 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> Jdelta201 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> JG201 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> I1301 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> J12101 = make_1Darray(pfo->k_size);

      pfo -> pk_halo_rsd_nl -> J21102x = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> J21102y = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> Jdelta202x = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> Jdelta202y = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> JG202x = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> JG202y = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> I2211 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> J21111 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> N11x = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> N11y = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> J12102x = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> J12102y = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> I1311 = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> J12111 = make_1Darray(pfo->k_size);

      pfo -> pk_halo_rsd_nl -> J21112x = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> J21112y = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> N12x = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> N12y = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> J12112x = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> J12112y = make_1Darray(pfo->k_size);
      
      pfo -> pk_halo_rsd_nl -> N22x = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> N22y = make_1Darray(pfo->k_size);
      pfo -> pk_halo_rsd_nl -> N22z = make_1Darray(pfo->k_size);


      pfo -> fft_ws = (struct oneloop_fftlog_workspace *)malloc(sizeof(struct oneloop_fftlog_workspace));
      pfo -> fft_ws -> fft_input  = (struct fft_struct *)malloc(sizeof(struct fft_struct));
      pfo -> fft_ws -> fft_matrix = (struct fft_matrices *)malloc(sizeof(struct fft_matrices));


      /* Important values for the calculation */
      pfo -> fft_ws -> sigma_v0 = 3. * sigman(pba, ppm, pfo, z, 1.e-5,  1.e3,  0, 142L); // density variance
      pfo -> fft_ws -> sigma_v2 = sigman(pba, ppm, pfo, z, 1.e-5,  1.e3, -1, 142L); // Linear displacement field (velocity) variance
      pfo -> fft_ws -> sigma_2_IR = IR_Sigma2(pba, ppm, pfo, z, 1e-4, 142L); // IR-Ressumation supression exponent

      // fprintf(stderr, "sigma_v0 = %e\nsigma_v2 = %e\n", pfo -> fft_ws -> sigma_v0, pfo -> fft_ws -> sigma_v2);

      /* FFTLog parameters */
      int    N_FFTLog   = 256; // fast -> 128, precise -> 256
      double kmin_fft_m = 1.e-8;
      double kmin_fft_g = 1.e-4;

      /* Setting the FFTLog parameters and calculating the etam and cmsym */
      pfo -> fft_ws -> fft_input -> nfft 	   = N_FFTLog;
      pfo -> fft_ws -> fft_input -> fft_bias_m = - 0.3;  
      pfo -> fft_ws -> fft_input -> kmin_fft_m = kmin_fft_m;    
      pfo -> fft_ws -> fft_input -> fft_bias_g = - 1.6; 
      pfo -> fft_ws -> fft_input -> kmin_fft_g = kmin_fft_g;

      pfo -> fft_ws -> fft_input -> etam_m  = make_1D_c_array(N_FFTLog + 1);
      pfo -> fft_ws -> fft_input -> cmsym_m = make_1D_c_array(N_FFTLog + 1);
      pfo -> fft_ws -> fft_input -> etam_g  = make_1D_c_array(N_FFTLog + 1);
      pfo -> fft_ws -> fft_input -> cmsym_g = make_1D_c_array(N_FFTLog + 1);

      FFT_compute_coeff(pba, ppm, pfo, z, pfo -> fft_ws -> fft_input, 142L, HALO);
      FFT_compute_coeff(pba, ppm, pfo, z, pfo -> fft_ws -> fft_input, 142L, MATTER);
      
      /* Setting the matrices */

      /* 0-th moment */
      // FFTLog matrices non-propagator
      pfo -> fft_ws -> fft_matrix -> I2200_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(I2200, pfo -> fft_ws -> fft_input, MATTER, pfo -> fft_ws -> fft_matrix ->  I2200_mat);
      pfo -> fft_ws -> fft_matrix -> Idelta200_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(Idelta200, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> Idelta200_mat);
      pfo -> fft_ws -> fft_matrix -> IG200_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(IG200, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> IG200_mat);
      pfo -> fft_ws -> fft_matrix -> Idelta2delta200_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(Idelta2delta200, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> Idelta2delta200_mat);
      pfo -> fft_ws -> fft_matrix -> IG2G200_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(IG2G200, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> IG2G200_mat);
      pfo -> fft_ws -> fft_matrix -> Idelta2G200_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(Idelta2G200, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> Idelta2G200_mat);

      // FFTLog matrices propagator
      pfo -> fft_ws -> fft_matrix -> I1300_mat   = make_1D_c_array(N_FFTLog+1);
      p_mat_fill(I1300, pfo -> fft_ws -> fft_input, MATTER, pfo -> fft_ws -> fft_matrix ->  I1300_mat);
      pfo -> fft_ws -> fft_matrix -> FG200_mat = make_1D_c_array(N_FFTLog+1);
      p_mat_fill(FG200, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> FG200_mat);

      /* 1-st moment */
      // FFTLog matrices (non-propagator)
      pfo -> fft_ws -> fft_matrix -> I2201_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(I2201, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> I2201_mat);
      pfo -> fft_ws -> fft_matrix -> Idelta201_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(Idelta201, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> Idelta201_mat);
      pfo -> fft_ws -> fft_matrix -> IG201_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(IG201, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> IG201_mat);
      pfo -> fft_ws -> fft_matrix -> FG201_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(FG201, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> FG201_mat);
      pfo -> fft_ws -> fft_matrix -> J21101_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(J21101, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J21101_mat);
      pfo -> fft_ws -> fft_matrix -> Jdelta201_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(Jdelta201, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> Jdelta201_mat);
      pfo -> fft_ws -> fft_matrix -> JG201_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(JG201, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> JG201_mat);

      // FFTLog matrices (propagator)
      pfo -> fft_ws -> fft_matrix -> I1301_mat = make_1D_c_array(N_FFTLog+1);
      p_mat_fill(I1301, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> I1301_mat);
      pfo -> fft_ws -> fft_matrix -> J12101_mat = make_1D_c_array(N_FFTLog+1);
      p_mat_fill(J12101, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J12101_mat);

      /* 2-nd moment */
      // FFTLog matrices (non-propagator)
      pfo -> fft_ws -> fft_matrix -> J21102x_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(J21102x, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J21102x_mat);
      pfo -> fft_ws -> fft_matrix -> J21102y_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(J21102y, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J21102y_mat);
      pfo -> fft_ws -> fft_matrix -> Jdelta202x_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(Jdelta202x, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> Jdelta202x_mat);
      pfo -> fft_ws -> fft_matrix -> Jdelta202y_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(Jdelta202y, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> Jdelta202y_mat);
      pfo -> fft_ws -> fft_matrix -> JG202x_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(JG202x, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> JG202x_mat);
      pfo -> fft_ws -> fft_matrix -> JG202y_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(JG202y, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> JG202y_mat);
      pfo -> fft_ws -> fft_matrix -> I2211_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(I2211, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> I2211_mat);
      pfo -> fft_ws -> fft_matrix -> J21111_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(J21111, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J21111_mat);
      pfo -> fft_ws -> fft_matrix -> N11x_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(N11x, pfo -> fft_ws -> fft_input, HALO,  pfo -> fft_ws -> fft_matrix -> N11x_mat);
      pfo -> fft_ws -> fft_matrix -> N11y_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(N11y, pfo -> fft_ws -> fft_input, HALO,  pfo -> fft_ws -> fft_matrix -> N11y_mat);

      // FFTLog matrices (propagator)
      pfo -> fft_ws -> fft_matrix -> J12102x_mat = make_1D_c_array(N_FFTLog+1);
      p_mat_fill(J12102x, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J12102x_mat);
      pfo -> fft_ws -> fft_matrix -> J12102y_mat = make_1D_c_array(N_FFTLog+1);
      p_mat_fill(J12102y, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J12102y_mat);
      pfo -> fft_ws -> fft_matrix -> I1311_mat = make_1D_c_array(N_FFTLog+1);
      p_mat_fill(I1311, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> I1311_mat);
      pfo -> fft_ws -> fft_matrix -> J12111_mat = make_1D_c_array(N_FFTLog+1);
      p_mat_fill(J12111, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J12111_mat);

      /* 3-rd moment */
      // FFTLog matrices (non-propagator)
      pfo -> fft_ws -> fft_matrix -> J21112x_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(J21112x, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J21112x_mat);
      pfo -> fft_ws -> fft_matrix -> J21112y_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(J21112y, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J21112y_mat);
      pfo -> fft_ws -> fft_matrix -> N12x_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(N12x, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> N12x_mat);
      pfo -> fft_ws -> fft_matrix -> N12y_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(N12y, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> N12y_mat);
      
      // FFTLog matrices (propagator)
      pfo -> fft_ws -> fft_matrix -> J12112x_mat = make_1D_c_array(N_FFTLog+1);
      p_mat_fill(J12112x, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J12112x_mat);
      pfo -> fft_ws -> fft_matrix -> J12112y_mat = make_1D_c_array(N_FFTLog+1);
      p_mat_fill(J12112y, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> J12112y_mat);

      /* 4-th moment */
      // FFTLog matrices (non-propagator)
      pfo -> fft_ws -> fft_matrix -> N22x_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(N22x, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> N22x_mat);
      pfo -> fft_ws -> fft_matrix -> N22y_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(N22y, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> N22y_mat);
      pfo -> fft_ws -> fft_matrix -> N22z_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
      np_mat_fill(N22z, pfo -> fft_ws -> fft_input, HALO, pfo -> fft_ws -> fft_matrix -> N22z_mat);

      return _SUCCESS_;
}