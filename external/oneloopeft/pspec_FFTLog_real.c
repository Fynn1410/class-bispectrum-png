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
                 int index_k,  double z, long SPLIT)
{
    double k = pfo->k[index_k];

    double plin       = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    double pm_lin_IR  = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
    double p_nowiggle = pm_nowiggle(pba, ppm, pfo, k, z, 1.e-4, 0, SPLIT);
    double p_wiggle   = plin - p_nowiggle;
    double sigma2     = pfo -> fft_ws -> sigma_2_IR;
    double sup        = exp(-k * k * sigma2);

    double sigmav2    = pfo -> fft_ws -> sigma_v2;
    
    P_mm_FFTLog(pfo, index_k, pm_lin_IR);

    double P22  = 2. * pfo -> pk_matter_real_nl -> I2200[index_k];
    double P13  = 6. * pfo -> pk_matter_real_nl -> I1300[index_k];

    /* 
     * Compute the EFT counter-term contribution
     */
    double cs2 =  0.2;
    double pm_ct   = - 2. * cs2 * pow(k, 2.) * pm_lin_IR;

    pfo -> pk_matter_real_nl -> P_mm = (p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2) + P22 + P13) + pm_ct; 

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
                    int index_k, double z, long SPLIT)

{ 
    double k = pfo->k[index_k];

    double b1  = 2.0;
    double b2  = -1.0;
    double bG2 = 0.1;
    double btd = -0.1;
    double R2  =  5.0;
    double cs2 = 0.2;

    double plin       = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    double pm_lin_IR  = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
    double p_nowiggle = pm_nowiggle(pba, ppm, pfo, k, z, 1.e-4, 0, SPLIT);
    double p_wiggle   = plin - p_nowiggle;
    double sigma2     = pfo -> fft_ws -> sigma_2_IR;
    double sup        = exp(-k * k * sigma2);

    double sigmav2    = pfo -> fft_ws -> sigma_v2;

    P_hh_FFTLog(pfo, index_k, pm_lin_IR);

    /* 
     * Compute the EFT counter-term contribution 
     */
    double pm_ct   = - 2. * cs2 * pow(k, 2.) * pm_lin_IR;

    double P22  = 2. * pfo -> pk_halo_real_nl -> I2200;
    double P13  = 6. * pfo -> pk_halo_real_nl -> I1300;
    double P_mm = pow(b1,2.) * ((p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2) + P22 + P13) + pm_ct); 

    double p_r2     = b1 * R2 * pfo -> pk_halo_real_nl[real_ir] -> IR2[index_k];
    double pb1b2    = 2. * b1 * b2  *  pfo -> pk_halo_real_nl[real_ir] -> Idelta200[index_k];
    double pb1bg2   = 4. * b1 * bG2 * pfo -> pk_halo_real_nl[real_ir] -> IG200[index_k];
    double pb22     = 0.5 * pow(b2, 2.) * pfo -> pk_halo_real_nl[real_ir] -> Idelta2delta200[index_k];
    double pbg22    = 2. * pow(bG2, 2.)  * pfo -> pk_halo_real_nl[real_ir] -> IG2G200[index_k];
    double pb2bg2   = 2. * b2 * bG2 * pfo -> pk_halo_real_nl[real_ir] -> Idelta2G200[index_k];
    double pb1b3nl  = 8. * b1 * (bG2 + 2./5. * btd) * pfo -> pk_halo_real_nl[real_ir] -> FG200[index_k];

    double ph_loops = p_r2 + pb1b2 + pb1bg2 + pb22 + pbg22 + pb2bg2 + pb1b3nl + P_mm;
    pfo -> pk_halo_real_nl -> P_hh[index_k] = ph_loops;


    // FILE *fpa;
    // char file_name[50];
    // sprintf(file_name, "data/pg_FFT_const.txt");
    // fpa = fopen(file_name, "a");
    // fprintf(fpa, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",\
    //     k, pm_lin_IR, pow(b1, 2.) * pm_1loop_IR, pm_ct, pb1b2, pb1bg2, pb22, pbg22, pb2bg2, pb1b3nl, ph_loops, ph_tot);
    // fclose(fpa);

    return _SUCCESS_;
}

