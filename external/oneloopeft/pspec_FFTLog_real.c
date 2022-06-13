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

	fft_input -> nfft 	    = 256;
    fft_input -> kmin_fft_g = 1.e-4;
	fft_input -> fft_bias_g = - 1.6;  //for halos
	fft_input -> kmin_fft_m = 1.e-8;
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
    
    // double P22_dif = P22(fft_input, k, z, cleanup_mloops) - P22_new(fft_input, k, z, cleanup_mloops);
    // double P13_dif = P13(fft_input, k, z, cleanup_mloops) - P13_new(fft_input, k, z, cleanup_mloops);
    
    // FILE *fpa;
    // char file_name[50];
    // sprintf(file_name, "FFTLog_new.txt");
    // fpa = fopen(file_name, "a");
    // fprintf(fpa, "%e %e %e %e %e\n", k, P22(fft_input, k, z, cleanup_mloops), P22_new(fft_input, k, z, cleanup_mloops), P13(fft_input, k, z, cleanup_mloops), P13_new(fft_input, k, z, cleanup_mloops));
    // fclose(fpa);

    // fprintf(stderr, "P22 Diff: %e\n", P22_dif);
    // fprintf(stderr, "P13 Diff: %e\n", P13_dif);

    /* 
     * Compute the EFT counter-term contribution
     */
    double pm_ct   = - 2. * cs2 * pow(k, 2.) * pm_lin_IR;

    // double ph_tot = (p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2) + P22_IR + P13_IR_tot) + pm_ct; 
    double ph_tot = (pm_lin_IR + P22_IR + P13_IR_tot) + pm_ct; 
    // FILE *fpa;
    // char file_name[50];
    // sprintf(file_name, "NOIRvsWIR.txt");
    // fpa = fopen(file_name, "a");
    // fprintf(fpa, "%e %e %e %e %e\n",k, p_nowiggle, p_wiggle, plin, pm_lin_IR);
    // fclose(fpa);
    
    //fprintf(stderr,"sigmav2 = %e, Plin = %e, P13_IR = %e, P13_UV = %e, P22__IR = %e, Plin_IR = %e, P_tot = %e\n",sigmav2, plin, P13_IR, P13_uv, P22_IR, p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2), ph_tot);
    
    FILE *fpa;
    char file_name[50];
    sprintf(file_name, "pm_FFTLog.txt");
    fpa = fopen(file_name, "a");
    fprintf(fpa,"%e %e %e %e %e %e %e %e %e\n", k, pm_lin_IR, p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2), P22_IR, P13_IR, P13_uv, P13_IR_tot, pm_ct, ph_tot);
    fclose(fpa);

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

	fft_input -> nfft 	    = 256;
	fft_input -> kmin_fft_g = 1.e-4;
	fft_input -> fft_bias_g = - 1.6;  //for halos
    fft_input -> kmin_fft_m = 1.e-8;
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

	fft_input -> nfft 	    = 256;
    fft_input -> kmin_fft_g = 1.e-4;
	fft_input -> fft_bias_g = - 1.6;  //for halos
	fft_input -> kmin_fft_m = 1.e-8;
	fft_input -> fft_bias_m = - 0.3; //for matter

    fft_input -> fft_first   = 1;
	fft_input -> etam_m  	 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> cmsym_m     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> etam_g      = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));
    fft_input -> cmsym_g     = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(fft_input->nfft+1));

    FFT_compute_coeff(pba, ppm, pfo, z, fft_input, SPLIT, HALO);
    FFT_compute_coeff(pba, ppm, pfo, z, fft_input, SPLIT, MATTER);

    // setting bias vectors, first translating Eulerian into Lagrangian biases and then distributing them onto the loops
    //DL CLASS-PT values page 30
    double b1  =  2.0;
    double b2  = -1.0;
    double bG2 =  0.1;
    double btd = -0.1;
    double cs2 =  0.2;
    double R2  =  5.0;

    // Distributing vectors for the biases and moment expansion factors, which will be applied on the loops
    
    // 0-th moment
    double *np_bias_vec0 = make_1Darray(6);
    double *p_bias_vec0  = make_1Darray(2);
    np_bias_vec0[0] = 2.* pow(b1, 2.);
    np_bias_vec0[1] = 2.* b1 * b2;
    np_bias_vec0[2] = 4.* b1 * bG2;
    np_bias_vec0[3] =  0.5 * pow(b2, 2.);
    np_bias_vec0[4] = 2.* pow(bG2, 2.);
    np_bias_vec0[5] = 2.* b2 * bG2;
    p_bias_vec0[0] = 6.* pow(b1, 2.);
    p_bias_vec0[1] = 8. * b1 * (bG2 + 2./5. * btd);

    // 1-st moment
    double *np_bias_vec1 = make_1Darray(7);
    double *p_bias_vec1  = make_1Darray(3);
    np_bias_vec1[0] = 2. * b1;
    np_bias_vec1[1] = b2;
    np_bias_vec1[2] = 2. * bG2;
    np_bias_vec1[3] = 4. * (bG2 + 2./5. * btd);
    np_bias_vec1[4] = pow(b1, 2.);
    np_bias_vec1[5] = 0.5 * b1 * b2;
    np_bias_vec1[6] = b1 * bG2;
    p_bias_vec1[0] = 3. * b1;
    p_bias_vec1[1] = pow(b1, 2.);
    p_bias_vec1[2] = pow(b1, 2.);
        
    double *np_expans_vec1 = make_1Darray(7);
    double *p_expans_vec1  = make_1Darray(3);
    np_expans_vec1[0] = f * pow(mu, 2.);
    np_expans_vec1[1] = f * pow(mu, 2.);
    np_expans_vec1[2] = f * pow(mu, 2.);
    np_expans_vec1[3] = f * pow(mu, 2.);
    np_expans_vec1[4] = 2. * f * mu * k;
    np_expans_vec1[5] = 2. * f * mu * k;
    np_expans_vec1[6] = 2. * f * mu * k;
    p_expans_vec1[0] = f * pow(mu, 2.);
    p_expans_vec1[1] = 2. * f * mu * k;
    p_expans_vec1[2] = 2. * f * mu * k;

    // 2-nd moment
    double *np_bias_vec2 = make_1Darray(6);
    double *p_bias_vec2  = make_1Darray(4);
    np_bias_vec2[0] = 2. * b1;
    np_bias_vec2[1] = b2;
    np_bias_vec2[2] = 2. * bG2;
    np_bias_vec2[3] = 2. ;
    np_bias_vec2[4] = b1;
    np_bias_vec2[5] = pow(b1, 2.);
    p_bias_vec2[0] = 4. * b1;
    p_bias_vec2[1] = 6.;
    p_bias_vec2[2] = b1;
    p_bias_vec2[3] = b1;

    double *np_expans_vec2 = make_1Darray(6);
    double *p_expans_vec2  = make_1Darray(4);
    np_expans_vec2[0] = 0.5 * pow(f * k * mu, 2.);
    np_expans_vec2[1] = 0.5 * pow(f * k * mu, 2.);
    np_expans_vec2[2] = 0.5 * pow(f * k * mu, 2.);
    np_expans_vec2[3] = - pow(f, 2.) * pow(mu, 4.);
    np_expans_vec2[4] = -4. * pow(f, 2.) * pow(mu, 3.) * k;
    np_expans_vec2[5] = - pow(f, 2.) * pow(k * mu, 4.);
    p_expans_vec2[0] = 0.5 * pow(f * k * mu, 2.);
    p_expans_vec2[1] = - pow(f, 2.) * pow(mu, 4.);
    p_expans_vec2[2] = -4. * pow(f, 2.) * pow(mu, 3.) * k;
    p_expans_vec2[3] = -4. * pow(f, 2.) * pow(mu, 3.) * k;

    // 3-rd moment
    double *np_bias_vec3 = make_1Darray(2);
    double *p_bias_vec3  = make_1Darray(1);
    np_bias_vec3[0] = 1.;
    np_bias_vec3[1] = b1;
    p_bias_vec3[0] = 2.;

    double *np_expans_vec3 = make_1Darray(2);
    double *p_expans_vec3  = make_1Darray(1);
    np_expans_vec3[0] = - pow(f, 3.) * pow(mu, 4.) * pow(k, 2.) / 3.;
    np_expans_vec3[1] = - pow(f * k * mu, 3.) / 3.;
    p_expans_vec3[0] = - pow(f, 3.) * pow(mu, 4.) * pow(k, 2.) / 3.;

    // 4-th moment
    double *np_bias_vec4 = make_1Darray(1);
    np_bias_vec4[0] = 1.;

    double *np_expans_vec4 = make_1Darray(1);
    np_expans_vec4[0] = pow(f * k * mu, 4.) / 12.;


    // Integrals calculation
    double *np_rsd_0 = make_1Darray(6);
    double *p_rsd_0  = make_1Darray(2);
    rsd_0_FFTLog(fft_input, k, np_rsd_0, p_rsd_0);

    double *np_rsd_1 = make_1Darray(7);
    double *p_rsd_1  = make_1Darray(3);
    // rsd_1_FFTLog(fft_input, k, mu, np_rsd_1, p_rsd_1);

    double *np_rsd_2 = make_1Darray(6);
    double *p_rsd_2  = make_1Darray(4);
    // rsd_2_FFTLog(fft_input, k, mu, np_rsd_2, p_rsd_2);

    double *np_rsd_3 = make_1Darray(2);
    double *p_rsd_3  = make_1Darray(1);
    // rsd_3_FFTLog(fft_input, k, mu, np_rsd_3, p_rsd_3);

    double *np_rsd_4 = make_1Darray(1);
    // rsd_4_FFTLog(fft_input, k, mu, np_rsd_4);

    // Bringing everything together
    double np_0  = 0.;
    double p_0   = 0.;
    dot(np_rsd_0, np_bias_vec0, 6, &np_0);
    dot(p_rsd_0,  p_bias_vec0,  2, &p_0);
    double ct_0   = - 2. * b1 * (R2 + cs2 * b1) * pow(k, 2.) * pm_lin_IR; // 0-th order counter term
    double P13_uv = - 61./105. * pm_lin_IR * pow(k, 2.) * sigmav2;
    double mom_0  = np_0 + pm_lin_IR * p_0 + pow(b1, 2.) * (pm_lin_IR + P13_uv + ct_0);
    
    FILE *fpa;
    char file_name[50];
    sprintf(file_name, "rsd_0_elements.txt");
    fpa = fopen(file_name, "a");
    double *np = make_1Darray(6);
    double *p = make_1Darray(2);
    vecmult(np_rsd_0, np_bias_vec0, 6, np);
    vecmult(p_rsd_0,  p_bias_vec0,  2, p);

    fprintf(fpa, "%12.6e ", k);
    fprintf(fpa, "%12.6e ", pm_lin);

    // fprintf(fpa, "%12.6e ", 2. * np_rsd_0[0]);
    // fprintf(fpa, "%12.6e ", 6. * pm_lin_IR * p_rsd_0[0]);
    // fprintf(fpa, "%12.6e ", P13_uv);
    // fprintf(fpa, "%12.6e ", ct_0);
    // fprintf(fpa, "%12.6e ", (pm_lin + 2. * np_rsd_0[0] + 6. * p_rsd_0[0] + P13_uv + ct_0));
    
    fprintf(fpa, "%12.6e ", pow(b1, 2.) * (np[0] + p[0] * pm_lin_IR + P13_uv + ct_0 + pm_lin_IR));
    fprintf(fpa, "%12.6e ", ct_0);
    fprintf(fpa, "%12.6e ", np[1]);
    fprintf(fpa, "%12.6e ", np[2]);
    fprintf(fpa, "%12.6e ", np[3]);
    fprintf(fpa, "%12.6e ", np[4]);
    fprintf(fpa, "%12.6e ", np[5]);
    fprintf(fpa, "%12.6e ", p[1] * pm_lin_IR);
    fprintf(fpa, "%12.6e ", np_0 + pm_lin_IR * p_0);
    fprintf(fpa, "%12.6e", mom_0);  
    fprintf(fpa, "\n");
    fclose(fpa);

    // double *np_prefactors_1 = make_1Darray(7);
    // double *p_prefactors_1  = make_1Darray(3);
    // vecmult(np_bias_vec1, np_expans_vec1, 7, np_prefactors_1);
    // vecmult(p_bias_vec1,  p_expans_vec1,  3, p_prefactors_1);
    // double np_1  = 0.;
    // double p_1   = 0.;
    // dot(np_rsd_1, np_prefactors_1, 7, &np_1);
    // dot(p_rsd_1,  p_prefactors_1,  3, &p_1);
    // double mom_1 = np_1 + pm_lin_IR * p_1 + b1 * pm_lin_IR * f * pow(mu,2.);

    // double *np_prefactors_2 = make_1Darray(6);
    // double *p_prefactors_2  = make_1Darray(4);
    // vecmult(np_bias_vec2, np_expans_vec2, 6, np_prefactors_2);
    // vecmult(p_bias_vec2,  p_expans_vec2,  4, p_prefactors_2);
    // double np_2  = 0.;
    // double p_2   = 0.;
    // dot(np_rsd_2, np_prefactors_2, 6, &np_2);
    // dot(p_rsd_2,  p_prefactors_2,  4, &p_2);
    // double mom_2 = np_2 + pm_lin_IR * p_2 + 2. * b1 * pm_lin_IR * pow(f * k * mu, 2.);

    // double *np_prefactors_3 = make_1Darray(2);
    // double *p_prefactors_3  = make_1Darray(1);
    // vecmult(np_bias_vec3, np_expans_vec3, 2, np_prefactors_3);
    // vecmult(p_bias_vec3,  p_expans_vec3,  1, p_prefactors_3);
    // double np_3  = 0.;
    // double p_3   = 0.;
    // dot(np_rsd_3, np_prefactors_3, 2, &np_3);
    // dot(p_rsd_3,  p_prefactors_3,  1, &p_3);
    // double mom_3 = np_3 + pm_lin_IR * p_3;

    // double *np_prefactors_4 = make_1Darray(1);
    // vecmult(np_bias_vec4, np_expans_vec4, 1, np_prefactors_4);
    // double np_4  = 0.;
    // dot(np_rsd_4, np_prefactors_4, 1, &np_4);
    // double mom_4 = np_4 ;

    // Summing over all moments
    // double p_tot = mom_0 + mom_1 + mom_2 + mom_3 + mom_4;
    double p_tot = mom_0;

    // FILE *fpa2;
    // char file_name2[50];
    // sprintf(file_name2, "FFTLog_rsd.txt");
    // fpa2 = fopen(file_name2, "a");
    // fprintf(fpa2, "%e %e %e %e %e %e %e %e\n", mu, k, mom_0, mom_1, mom_2, mom_3, mom_4, p_tot);
    // fclose(fpa2);

    *pk_nl = p_tot;
    return _SUCCESS_;
}

