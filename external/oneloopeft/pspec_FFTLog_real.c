/** @file pspec_FFTLog.c Documented FFT-Log based 1loop matter and galaxy power spectrum 
 * 
 * Azadeh Moradinezhad Dizgah,June 16th 2021
 *
 *
 * In summary, the following functions can be called from other modules:
 * -# pm_IR_FFTLog()
 * -# pg_IR_FFTLog()
 *
 */
 
#include "header.h"


/**
 * Compute the IR-resummed 1loop matter power spectrum, using FFTLog algorithm for the loop integrals
 *
 * @param cosmo        Input: cosmo structure cotaining the PK and TK interpolators
 * @param block        Input: cosmosis datablock
 * @param config_in    Input: cosmosis configuration structure constructed from the .ini file
 * @param ptrs         Input: pointers to interpolatir ibjects for nw power spectrum
 * @param fft_struct   Input: structure containing fft coefficents and params
 * @param k            Input: wavenumber 
 * @param z            Input: redshift of the pk_m
* @return value of 1loop pk_m
 */

int pm_IR_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                 double k,  double z, long SPLIT, double * pk_nl)
{
    static int cleanup_mloops  = 0;

    // setting fft_parameters
    struct fft_struct *fft_input;
	fft_input = (struct fft_struct *)malloc(sizeof(struct fft_struct));

	fft_input -> nfft 	    = 200;
	fft_input -> kmin_fft   = 1.e-8;
	fft_input -> fft_bias_g = - 1.6;  //for halos
	fft_input -> fft_bias_m = - 0.3; //for matter

    fft_input -> fft_first   = 1;
	fft_input -> etam_m  	 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> cmsym_m     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> etam_g      = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> cmsym_g     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));

    FFT_compute_coeff(pba, ppm, pfo, z, fft_input, SPLIT, MATTER);

    double cs2 =  0.2;

    double k0         = 1.e-5;
    double k_max      = 1.e3;
    double plin       = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    double pm_lin_IR  = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
    // grwoth factor, needs to be found within CLASS -> workaround: working with the time evolution of pm_nonwiggle
    // double growth2    = pow(gsl_spline_eval(cosmo->DZ, z, NULL)/gsl_spline_eval(cosmo->DZ, 0., NULL), 2.);
    double p_nowiggle = pm_nowiggle(pba, ppm, pfo, k, z, k0, 0, SPLIT);
    double p_wiggle   = plin - p_nowiggle;
    double sigma2     = IR_Sigma2(pba, ppm, pfo, z, k0, SPLIT);
    double sup        = exp(-k * k * sigma2);

    double sigmav2    = sigman(pba, ppm, pfo, z, k0, k_max, -1, SPLIT);
    double P22_IR     = P22(fft_input, k, z, cleanup_mloops);
    double P13_IR     = pm_lin_IR * P13(fft_input, k, z, cleanup_mloops);
    double P13_uv     = - 61./105. * pm_lin_IR * pow(k, 2.) * sigmav2;
    double P13_IR_tot = P13_IR + P13_uv;
    
    double P22_dif = P22(fft_input, k, z, cleanup_mloops) - P22_new(fft_input, k, z, cleanup_mloops);
    double P13_dif = P13(fft_input, k, z, cleanup_mloops) - P13_new(fft_input, k, z, cleanup_mloops);
    
    FILE *fpa;
    char file_name[50];
    sprintf(file_name, "FFTLog_new.txt");
    fpa = fopen(file_name, "a");
    fprintf(fpa, "%e %e %e %e %e\n", k, P22(fft_input, k, z, cleanup_mloops), P22_new(fft_input, k, z, cleanup_mloops), P13(fft_input, k, z, cleanup_mloops), P13_new(fft_input, k, z, cleanup_mloops));
    fclose(fpa);

    fprintf(stderr, "P22 Diff: %e\n", P22_dif);
    fprintf(stderr, "P13 Diff: %e\n", P13_dif);

    /* 
     * Compute the EFT counter-term contribution
     */
    double pm_ct   = - 2. * cs2 * pow(k, 2.) * pm_lin_IR;

    double ph_tot = (p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2) + P22_IR + P13_IR_tot) + pm_ct; 
    
    // FILE *fpa;
    // char file_name[50];
    // sprintf(file_name, "NOIRvsWIR.txt");
    // fpa = fopen(file_name, "a");
    // fprintf(fpa, "%e %e %e %e %e\n",k, p_nowiggle, p_wiggle, plin, pm_lin_IR);
    // fclose(fpa);
    
    //fprintf(stderr,"sigmav2 = %e, Plin = %e, P13_IR = %e, P13_UV = %e, P22__IR = %e, Plin_IR = %e, P_tot = %e\n",sigmav2, plin, P13_IR, P13_uv, P22_IR, p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2), ph_tot);
    
    // FILE *fpa;
    // char file_name[50];
    // sprintf(file_name, "pm_FFTLog.txt");
    // fpa = fopen(file_name, "a");
    // fprintf(fpa,"%e %e %e %e %e %e %e %e %e\n", k, pm_lin_IR, p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2), P22_IR, P13_IR, P13_uv, P13_IR_tot, pm_ct, ph_tot);
    // fclose(fpa);

    *pk_nl = ph_tot;
    return _SUCCESS_;
}


/**
 * Compute the IR-resummed 1loop galaxy power spectrum, using FFTLog algorithm for the loop integrals
 *
 * @param cosmo        Input: cosmo structure cotaining the PK and TK interpolators
 * @param block        Input: cosmosis datablock
 * @param config_in    Input: cosmosis configuration structure constructed from the .ini file
 * @param ptrs         Input: pointers to interpolatir ibjects for nw power spectrum
 * @param fft_struct   Input: structure containing fft coefficents and params
 * @param k            Input: wavenumber 
 * @param z            Input: redshift of the pk_m
 * @return value of 1loop pk_g
 */

int pg_IR_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    double k,  double z, long SPLIT, double * pk_nl)

{ 
    static int cleanup_gloops = 0;

    // setting fft_parameters
    struct fft_struct *fft_input;
	fft_input = (struct fft_struct *)malloc(sizeof(struct fft_struct));

	fft_input -> nfft 	    = 200;
	fft_input -> kmin_fft   = 5.e-4;
	fft_input -> fft_bias_g = - 1.6;  //for halos
	fft_input -> fft_bias_m = - 0.3; //for matter

    fft_input -> fft_first   = 1;
	fft_input -> etam_m  	 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> cmsym_m     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> etam_g      = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> cmsym_g     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));

    FFT_compute_coeff(pba, ppm, pfo, z, fft_input, SPLIT, HALO);

      //DL CLASS-PT values page 30
      double b1  =  2.0;
      double b2  = -1.0;
      double bG2 =  0.1;
      double btd = -0.1;
      double cs2 =  0.2;
      double R2  =  5.0;

    double pm_1loop_IR;
    pm_IR_FFTLog(pba, ppm, pfo, k, z, SPLIT, &pm_1loop_IR);
    double pm_lin_IR   = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
    double pm_lin      = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);

    /* 
     * Compute the 1loop IR-resummed loops of galaxy power spectrum
     */
    double prop_term = pgloops_propag(fft_input, k, z, cleanup_gloops);
    double *ps_hloops = make_1Darray(5);
    pgloops_nonpropag(fft_input, k, z, cleanup_gloops, ps_hloops);

    double pb1b2      = b1 * b2  * ps_hloops[0];
    double pb1bg2     = 2. * b1 * bG2 * ps_hloops[1];
    double pb22       = 0.25 * pow(b2, 2.) * ps_hloops[2];
    double pbg22      = pow(bG2, 2.)  * ps_hloops[3];
    double pb2bg2     = b2 * bG2 * ps_hloops[4];
    double pb1b3nl    = 2. * b1 * (bG2 + 2./5. * btd) * pm_lin_IR * prop_term;
    double ph_loops   = pb1b2 + pb1bg2 + pb22+ pbg22 + pb2bg2 + pb1b3nl;

    /* 
     * Compute the EFT counter-term contribution
     */
    double pm_ct   = - 2. * b1 * (R2 + cs2 * b1) * pow(k, 2.) * pm_lin_IR;

    /* 
     * Compute the EFT counter-term contribution
     */
    double ph_tot  = pow(b1, 2.) * pm_1loop_IR + pm_ct + ph_loops;


    FILE *fpa;
    char file_name[50];
    sprintf(file_name, "pg_FFTLog.txt");
    fpa = fopen(file_name, "a");
    fprintf(fpa, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",\
                k, pm_lin, pow(b1, 2.) * pm_1loop_IR, pm_ct, pb1b2, pb1bg2, pb22, pbg22, pb2bg2, pb1b3nl, ph_loops, ph_tot);
    fclose(fpa);

    *pk_nl = ph_tot;
    return _SUCCESS_;
}


int rsd_oneloop_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                 double k,  double z, double mu, double f, long SPLIT, double * pk_nl)
{
    // getting the linear power spectrum
    double k0         = 1.e-5;
    double k_max      = 1.e3;
    double pm_lin_IR = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
    double pm_lin    = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    double sigmav2   = sigman(pba, ppm, pfo, z, k0, k_max, -1, SPLIT);

    // setting fft_parameters
    struct fft_struct *fft_input;
	fft_input = (struct fft_struct *)malloc(sizeof(struct fft_struct));

	fft_input -> nfft 	    = 200;
	fft_input -> kmin_fft   = 1.e-8;
	fft_input -> fft_bias_g = - 1.6;  //for halos
	fft_input -> fft_bias_m = - 0.3; //for matter

    fft_input -> fft_first   = 1;
	fft_input -> etam_m  	 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> cmsym_m     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> etam_g      = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> cmsym_g     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));

    FFT_compute_coeff(pba, ppm, pfo, z, fft_input, SPLIT, HALO);

    // setting bias vectors, first translating Eulerian into Lagrangian biases and then distributing them onto the loops
    //DL CLASS-PT values page 30
    double b1  =  2.0;
    double b2  = -1.0;
    double bG2 =  0.1;
    double btd = -0.1;
    double cs2 =  0.2;
    double R2  =  5.0;

    // EPT to LPT correspondance Chen page 14
    double c1 = 1 + b1;
    double c2 = b2 + (8./21) * b1;
    double c3 = btd + (1./3.) * b1; // should be b3, not sure of the prefactor
    double cs = bG2 - (2./7.) *b1; // should be bs, but no clue which one this should be

    // Distributing vector for the biases and LoS factors on the loops
    double *bias_vec0 = make_1Darray(9);
    bias_vec0[8] = pow(c1, 2.) * pm_lin_IR;
    bias_vec0[0] = pow(c1, 2.);
    bias_vec0[1] = pow(c1, 2.);
    bias_vec0[2] = c1 * c2;
    bias_vec0[3] = c1 * cs;
    bias_vec0[4] = pow(c2, 2.);
    bias_vec0[5] = c2 * cs;
    bias_vec0[6] = pow(c2, 2.);
    bias_vec0[7] = c1 * c3 * pm_lin_IR; // psi

    double *bias_vec1 = make_1Darray(10);
    bias_vec1[9] = c1 * pm_lin_IR * (2. * pow(mu, 2.) *f);
    bias_vec1[0] = c1 * (2. * pow(mu, 2.) *f);
    bias_vec1[1] = c1 * pm_lin_IR * (2. * pow(mu, 2.) *f);
    bias_vec1[2] = c2 * (2. * pow(mu, 2.) *f);
    bias_vec1[3] = cs * (2. * pow(mu, 2.) *f);
    bias_vec1[4] = c3 * pm_lin_IR * (2. * pow(mu, 2.) *f); // psi

    bias_vec1[5] = pow(c1, 2.) * (2. * 2. * f * mu * k) * (mu * k);
    bias_vec1[6] = pow(c1, 2.) * pm_lin_IR * (2. * 2. * f * mu * k) * (mu * k);
    bias_vec1[7] = c1 * c2 * (2. * 2. * f * mu * k) * (mu * k);
    bias_vec1[8] = c1 * cs * (2. * 2. * f * mu * k) * (mu * k);

    double *bias_vec2 = make_1Darray(14);
    bias_vec2[0]  = c1 * (-1. * pow(f*mu*k, 2.)) * pow(k, 2.);
    bias_vec2[1]  = c1 * pm_lin_IR * (-1. * pow(f*mu*k, 2.)) * pow(k, 2.);
    bias_vec2[2]  = c2 * (-1. * pow(f*mu*k, 2.)) * pow(k, 2.);
    bias_vec2[3]  = cs * (-1. * pow(f*mu*k, 2.)) * pow(k, 2.);
    bias_vec2[12] = pow(c1, 2.) * pm_lin_IR * sigmav2 * (-1. * pow(f*mu*k, 2.));

    bias_vec2[13] = pm_lin_IR * (-1. * pow(f, 2.) * pow(mu, 4.));
    bias_vec2[4]  = (-1. * pow(f, 2.) * pow(mu, 4.));
    bias_vec2[5]  = pm_lin_IR * (-1. * pow(f, 2.) * pow(mu, 4.));
    bias_vec2[6]  = 4. * c1 * (-1. * k * pow(f, 2.) * pow(mu, 3.)) * (mu * k);
    bias_vec2[7]  = 4. * c1 * pm_lin_IR * (-1. * k * pow(f, 2.) * pow(mu, 3.)) * (mu * k);
    bias_vec2[8]  = 4. * c1 * pm_lin_IR * (-1. * k * pow(f, 2.) * pow(mu, 3.)) * (mu * k);
    bias_vec2[9]  = pow(c1, 2.) * (-1. * pow(f*mu*k, 2.)) * pow(k, 2.);
    bias_vec2[10] = pow(c1, 2.) * (-1. * pow(f*mu*k, 2.)) * pow(mu*k, 2.);
    bias_vec2[11] = pow(c1, 2.) * (-1. * pow(f*mu*k, 2.)) * pow(k, 2.);

    double *bias_vec3 = make_1Darray(6);
    bias_vec3[4]  = c1 * pm_lin_IR * sigmav2 * (pow(f, 3.) * pow(mu, 4.) * pow(k, 2.));
    bias_vec3[0]  = pow(mu, 4.) * pow(k, 2.) * pow(f, 3.) * pow(k, 2.);
    bias_vec3[1]  = pm_lin_IR * pow(mu, 4.) * pow(k, 2.) * pow(f, 3.) * pow(k, 2.);
    bias_vec3[2]  = c1 * pow(mu, 3.) * pow(k, 3.) * pow(f, 3.) * (pow(k, 3.) * mu);
    bias_vec3[3]  = c1 * pow(mu, 3.) * pow(k, 3.) * pow(f, 3.) * (pow(k, 3.) * mu);
    bias_vec3[5]  = c1 * pm_lin_IR * sigmav2 * (pow(f, 3.) * pow(mu, 4.) * pow(k, 2.));

    double *bias_vec4 = make_1Darray(6);
    bias_vec4[3]  = 3. * pow(c1*sigmav2, 2.) * pm_lin_IR * pow(f*mu*k, 4) / 12.;
    bias_vec4[4]  = -1. *sigmav2 * pm_lin_IR * pow(f, 4) * pow(mu, 6) * pow(k, 2);
    bias_vec4[0]  = 3. * pow(f*mu*k, 4) / 12. * (pow(k, 4.) * pow(mu, 2.));
    bias_vec4[1]  = 3. * pow(f*mu*k, 4) / 12. * pow(k, 4.);
    bias_vec4[2]  = 3. * pow(f*mu*k, 4) / 12. * (pow(k, 4.) * pow(mu, 2.));
    bias_vec4[5]  = 3. * pow(c1*sigmav2, 2.) * pm_lin_IR * pow(f*mu*k, 4) / 12.;

    double *rsd_0 = make_1Darray(9);
    rsd_0_FFTLog(fft_input, k, rsd_0);
    rsd_0[8] = 1.0; // Adding the linear PS through the bias vector

    double *rsd_1 = make_1Darray(10);
    rsd_1_FFTLog(fft_input, k, rsd_1);
    rsd_1[9] = 1.0; // Adding the linear PS through the bias vector

    double *rsd_2 = make_1Darray(14);
    rsd_2_FFTLog(fft_input, k, mu, rsd_2);
    rsd_2[12] = 1.0; // Adding the linear PS through the bias vector
    rsd_2[13] = 1.0; // Adding the linear PS through the bias vector

    double *rsd_3 = make_1Darray(6);
    rsd_3_FFTLog(fft_input, k, mu, rsd_3);
    rsd_3[4] = 1.0; // Adding the linear PS through the bias vector
    rsd_3[5] = 1.0; // Adding the linear PS through the bias vector

    double *rsd_4 = make_1Darray(6);
    rsd_4_FFTLog(fft_input, k, mu, rsd_4);
    rsd_4[3] = 1.0; // Adding the linear PS through the bias vector
    rsd_4[4] = 1.0; // Adding the linear PS through the bias vector
    rsd_4[5] = 1.0; // Adding the linear PS through the bias vector

    double mom_0 = 0.;
    dot(rsd_0, bias_vec0, 9, &mom_0);
    mom_0 = mom_0 * pow(k, 3.);

    double mom_1 = 0.;
    dot(rsd_1, bias_vec0, 10, &mom_1);
    mom_1 = mom_1 * pow(k, 3.);

    double mom_2 = 0.;
    dot(rsd_2, bias_vec0, 14, &mom_2);
    mom_2 = mom_2 * pow(k, 3.);

    double mom_3 = 0.;
    dot(rsd_3, bias_vec0, 6, &mom_3);
    mom_3 = mom_3 * pow(k, 3.);

    double mom_4 = 0.;
    dot(rsd_4, bias_vec0, 6, &mom_4);
    mom_4 = mom_4 * pow(k, 3.);

    double p_tot = mom_0 + mom_1 + mom_2 + mom_3 + mom_4;

    FILE *fpa;
    char file_name[50];
    sprintf(file_name, "FFTLog_rsd.txt");
    fpa = fopen(file_name, "a");
    fprintf(fpa, "%e %e %e %e %e %e %e %e\n", mu, k, mom_0, mom_1, mom_2, mom_3, mom_4, p_tot);
    fclose(fpa);

    *pk_nl = p_tot;
    return _SUCCESS_;
}

