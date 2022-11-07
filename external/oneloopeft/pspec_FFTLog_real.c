/** @file pspec_FFTLog_real.c Documented FFT-Log based 1loop matter and galaxy power spectrum 
 * 
 * Dennis Linde, September 19th 2022
 * credits to: Azadeh Moradinezhad Dizgah
 *
 *
 * In summary, the following functions can be called from other modules:
 * -# Real_Oneloop_FFTLog()
 * -# Real_Matter_IR_Resummed()
 * -# Real_Galaxy_IR_Resummed()
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

int Real_Oneloop_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    int index_k, double z, long SPLIT)
{
    double k = pfo->k[index_k];

    // Storing IR-resummed Plin in pfo
    double pm_lin_IR  = pm_IR_LO(pba, ppm, pfo, k, 0., SPLIT); // Power Spectrum at redshift z=0., loops get scaled individually

    // Calculation of the Loop-Integrals for the real-space IR-Resummed linear power spectrum
    P_mm_FFTLog(pfo, index_k, pm_lin_IR);
    P_gg_FFTLog(pfo, index_k, pm_lin_IR);

    return _SUCCESS_;
}

int Real_Matter_IR_Resummed(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                 int index_k,  double z, long SPLIT, double *pk)
{
    double k = pfo->k[index_k];

    double cs2 = pfo->cs2;

    double D  = growth_D(pba, z);
    double D2 = pow(D,2.);
    double D4 = pow(D,4.);

    double pm_lin_IR  = pfo -> pk_matter_real_nl -> Plin_IR[index_k]; // Power Spectrum at redshift z=0., loops get scaled individually
    double plin       = Pk_dlnPk(pba, ppm, pfo, k, 0., LPOWER);
    double p_nowiggle = pm_nowiggle(pba, ppm, pfo, k, 0., 1.e-4, 0, SPLIT);
    double p_wiggle   = plin - p_nowiggle;
    double sigma2     = pfo -> fft_ws -> sigma_2_IR * D2;
    double sup        = exp(-k * k * sigma2);

    double sigmav2    = pfo -> fft_ws -> sigma_v2;

    double P22  = 2. * pfo -> pk_matter_real_nl -> I2200[index_k];
    double P13  = 6. * pfo -> pk_matter_real_nl -> I1300[index_k];

    /* 
     * Compute the EFT counter-term contribution
     */
    double pm_ct   = - 2. * cs2 * pow(k, 2.) * pm_lin_IR;
    pfo -> pk_matter_real_nl -> P_mm[index_k] = (p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2) + pm_ct)*D2 + (P22 + P13)*D4; 

    *pk = pfo -> pk_matter_real_nl -> P_mm[index_k];

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

int Real_Galaxy_IR_Resummed_default(struct fourier *pfo, struct background *pba, struct primordial *ppm,
                    int index_k, double z, long SPLIT, double *pk)
{
    double b1  = pfo->b1;
    double b2  = pfo->b2; 
    double bG2 = pfo->bG2;
    double btd = pfo->btd; 
    double R2  = pfo->R2; 
    double cs2 = pfo->cs2;

    class_call(Real_Galaxy_IR_Resummed(pfo, pba, ppm, index_k, z, b1, b2, bG2, btd, R2, cs2, SPLIT, pk),
                pfo->error_message,
                pfo->error_message);

    return _SUCCESS_;
}


int Real_Galaxy_IR_Resummed(struct fourier *pfo, struct background *pba, struct primordial *ppm,
                    int index_k, double z,
                    double b1, double b2, double bG2, double btd, double R2,
                    double cs2,
                    long SPLIT,
                    double *pk)
{ 
    double k = pfo->k[index_k];

    double D  = growth_D(pba, z);
    double D2 = pow(D,2.);
    double D4 = pow(D,4.);

    double pm_lin_IR  = pfo -> pk_halo_real_nl -> Plin_IR[index_k]; // Power Spectrum at redshift z=0., loops get scaled individually
    double plin       = Pk_dlnPk(pba, ppm, pfo, k, 0., LPOWER);
    double p_nowiggle = pm_nowiggle(pba, ppm, pfo, k, 0., 1.e-4, 0, SPLIT);
    double p_wiggle   = plin - p_nowiggle;
    double sigma2     = pfo -> fft_ws -> sigma_2_IR *D2;
    double sup        = exp(-k * k * sigma2);

    double sigmav2    = pfo -> fft_ws -> sigma_v2;

    double P22  = 2. * pfo -> pk_halo_real_nl -> I2200[index_k];
    double P13  = 6. * pfo -> pk_halo_real_nl -> I1300[index_k];
    double P_mm = pow(b1,2.) * ((p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2))*D2 + (P22 + P13)*D4); 
    pfo -> pk_halo_real_nl -> P_mm[index_k] = P_mm/pow(b1,2.);

    double p_r2     = -2. * pow(k,2.) * (b1 * R2 + pow(b1,2.) *cs2) * pm_lin_IR;
    double pb1b2    = 2. * b1 * b2  *  pfo -> pk_halo_real_nl -> Idelta200[index_k];
    double pb1bg2   = 4. * b1 * bG2 * pfo -> pk_halo_real_nl -> IG200[index_k];
    double pb22     = 0.5 * pow(b2, 2.) * pfo -> pk_halo_real_nl -> Idelta2delta200[index_k];
    double pbg22    = 2. * pow(bG2, 2.)  * pfo -> pk_halo_real_nl -> IG2G200[index_k];
    double pb2bg2   = 2. * b2 * bG2 * pfo -> pk_halo_real_nl -> Idelta2G200[index_k];
    double pb1b3nl  = 8. * b1 * (bG2 + 2./5. * btd) * pfo -> pk_halo_real_nl -> FG200[index_k];

    double ph_loops = p_r2*D2 + (pb1b2 + pb1bg2 + pb22 + pbg22 + pb2bg2 + pb1b3nl)*D4 + P_mm;
    pfo -> pk_halo_real_nl -> P_hh[index_k] = ph_loops;

    *pk = ph_loops;

    return _SUCCESS_;
}

