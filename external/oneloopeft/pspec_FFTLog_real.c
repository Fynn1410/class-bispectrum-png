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

    double cs2 =  0.2;

    double k0         = 1.e-5;
    double k_max      = 1.e3;
    double plin       = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    double pm_lin_IR  = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
    double p_nowiggle = pm_nowiggle(pba, ppm, pfo, k, z, k0, 0, SPLIT);
    double p_wiggle   = plin - p_nowiggle;
    double sigma2     = IR_Sigma2(pba, ppm, pfo, z, k0, SPLIT);
    double sup        = exp(-k * k * sigma2);

    double sigmav2    = sigman(pba, ppm, pfo, z, k0, k_max, -1, SPLIT);
    double P22_IR     = P22_new(pfo -> fft_ws -> fft_input, k, z, cleanup_mloops);
    double P13_IR     = pm_lin_IR * P13_new(pfo -> fft_ws -> fft_input, k, z, cleanup_mloops);
    double P13_uv     = - 61./105. * pm_lin_IR * pow(k, 2.) * sigmav2;
    double P13_IR_tot = P13_IR + P13_uv;

    /* 
     * Compute the EFT counter-term contribution
     */
    double pm_ct   = - 2. * cs2 * pow(k, 2.) * pm_lin_IR;

    double ph_tot = (p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2) + P22_IR + P13_IR_tot) + pm_ct; 

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
                    double k,  double z, long SPLIT)

{ 
    static int cleanup_gloops = 0;

    double pm_1loop_IR;
    pm_IR_FFTLog(pba, ppm, pfo, k, z, SPLIT, &pm_1loop_IR);
    double pm_lin_IR   = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
    double pm_lin      = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);

    /* 
     * Compute the 1loop IR-resummed loops of galaxy power spectrum
     */
    double prop_term = pgloops_propag(pfo -> fft_ws -> fft_input, k, z, cleanup_gloops);
    double *ps_hloops = make_1Darray(5);
    pgloops_nonpropag(pfo -> fft_ws -> fft_input, k, z, cleanup_gloops, ps_hloops);

    pfo -> oneloop_fftloh_halo_real -> pmm     = pm_1loop_IR;
    pfo -> oneloop_fftloh_halo_real -> pb1b2   = ps_hloops[0];
    pfo -> oneloop_fftloh_halo_real -> pb1bg2  = ps_hloops[1];
    pfo -> oneloop_fftloh_halo_real -> pb22    = ps_hloops[2];
    pfo -> oneloop_fftloh_halo_real -> pbg22   = ps_hloops[3];
    pfo -> oneloop_fftloh_halo_real -> pb2bg2  = ps_hloops[4];
    pfo -> oneloop_fftloh_halo_real -> pb1b3nl = pm_lin_IR * prop_term;
    pfo -> oneloop_fftloh_halo_real -> pm_ct   = pow(k, 2.) * pm_lin_IR;

    return _SUCCESS_;
}

