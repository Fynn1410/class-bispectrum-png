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

double pm_IR_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    struct fft_struct *fft_input, double k,  double z)
{
    static int cleanup_mloops  = 0;

    double k0         = 1.e-4;
    double plin       = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    // grwoth factor, needs to be found within CLASS
    double growth2    = pow(gsl_spline_eval(cosmo->DZ, z, NULL)/gsl_spline_eval(cosmo->DZ, 0., NULL), 2.);
    double p_nowiggle = growth2 * pm_nowiggle(pba, ppm, pfo, k, z);
    double p_wiggle   = plin - p_nowiggle;
    double sup        = exp(-k * k * cosmo->sig2_IR);

    double sigmav2    = cosmo->sig2;
    double P22_IR     = P22(fft_input, k, z, cleanup_mloops);
    double P13_IR     = pm_IR_LO(cosmo, block, config, ptrs, k, z) 
                      * P13(fft_input, k, z, cleanup_mloops);
    double P13_uv     = - 61./315. * pm_IR_LO(cosmo, block, config, ptrs,  k, z) * pow(k, 2.) * sigmav2;
    double P13_IR_tot = P13_IR + P13_uv;
    
    double out = (p_nowiggle + sup * p_wiggle * (1. + k * k * cosmo->sig2_IR) + P22_IR + P13_IR_tot); 
    
    printf("%12.6e %12.6e %12.6e %12.6e %12.6e \n",  k, pk_lin(cosmo->PK, k, z), P22_IR, P13_IR_tot, out);

    return out;
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

double pg_IR_FFTLog(struct Cosmology *cosmo, c_datablock * block, void *config_in, 
            struct interp_ptrs *ptrs, struct fft_struct  *fft_input, double k,  double z)

{    
    int Nmax  = fft_input->nfft;

    static int cleanup_gloops = 0;
    options_config * config = (options_config*) config_in;

    double b1 = 0., cs2 = 0., b2 = 0., bG2 = 0., btd = 0.;

    DATABLOCK_STATUS status = 0;
    
    double As;
    status |= c_datablock_get_double_default(block, "cosmological_parameters", "A_s", 3.e-9, &As);
    status |= c_datablock_get_double_default(block, "galaxy_bias", "b1", 0., &b1);
    status |= c_datablock_get_double_default(block, "galaxy_bias", "b2", 0., &b2);
    status |= c_datablock_get_double_default(block, "galaxy_bias", "bG2", 0., &bG2);
    status |= c_datablock_get_double_default(block, "galaxy_bias", "btd", 0., &btd);
    status |= c_datablock_get_double_default(block, "cosmological_parameters", "cs2", 0., &cs2);


    double pm_1loop_IR = pm_IR_FFTLog(cosmo, block, config, ptrs, fft_input, k, z);
    double pm_lin_IR   = pm_IR_LO(cosmo, block, config, ptrs, k, z);
    double pm_lin      = pk_lin(cosmo->PK, k, z);
    
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
    double pm_ct   = - 2. * cs2 * pow(k, 2.) * pm_lin_IR;

    /* 
     * Compute the EFT counter-term contribution
     */
    double ph_tot  = (pow(b1, 2.) * (pm_1loop_IR + pm_ct) + ph_loops);

    printf("%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",\
        k, pm_lin, pm_1loop_IR, pm_ct, pb1b2, pb1bg2, pb22, pbg22, pb2bg2, pb1b3nl,ph_tot);

    return ph_tot;
}


