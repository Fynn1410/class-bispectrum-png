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

int rsd_oneloop_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    int index_k, double z, long SPLIT)
{
    // Storing P_mm in pfo
    double pm_1loop_IR;
    pm_IR_FFTLog(pba, ppm, pfo, k, z, SPLIT, &pm_1loop_IR);
    pfo -> oneloop_fftlog_halo_rsd -> N22x[index_k] = pm_1loop_IR;

    // Calculation of the Loop-Integrals
    rsd_0_FFTLog(pfo, index_k);
    rsd_1_FFTLog(pfo, index_k);
    rsd_2_FFTLog(pfo, index_k);
    rsd_3_FFTLog(pfo, index_k);
    rsd_4_FFTLog(pfo, index_k);

    return _SUCCESS_;
}



int FFTLog_fill_bias_vector(struct fourier *pfo, double b1, double b2, double bG2, double btd, double cs2, double R2){
    pfo -> fft_ws -> b1  = b1;
    pfo -> fft_ws -> b2  = b2;
    pfo -> fft_ws -> bG2 = bG2;
    pfo -> fft_ws -> btd = btd;
    pfo -> fft_ws -> cs2 = cs2;
    pfo -> fft_ws -> R2  = R2;

    // 0-th moment
    pfo -> fft_ws -> bias -> np_bias_vec0 = make_1Darray(6);
    pfo -> fft_ws -> bias -> p_bias_vec0  = make_1Darray(2);
    pfo -> fft_ws -> bias -> np_bias_vec0[0] = 2.* pow(b1, 2.);
    pfo -> fft_ws -> bias -> np_bias_vec0[1] = 2.* b1 * b2;
    pfo -> fft_ws -> bias -> np_bias_vec0[2] = 4.* b1 * bG2;
    pfo -> fft_ws -> bias -> np_bias_vec0[3] = 0.5 * pow(b2, 2.);
    pfo -> fft_ws -> bias -> np_bias_vec0[4] = 2.* pow(bG2, 2.);
    pfo -> fft_ws -> bias -> np_bias_vec0[5] = 2.* b2 * bG2;
    pfo -> fft_ws -> bias -> p_bias_vec0[0] = 6.* pow(b1, 2.);
    pfo -> fft_ws -> bias -> p_bias_vec0[1] = 8. * b1 * (bG2 + 2./5. * btd);

    // 1-st moment
    pfo -> fft_ws -> bias -> np_bias_vec1 = make_1Darray(7);
    pfo -> fft_ws -> bias -> p_bias_vec1  = make_1Darray(2);
    pfo -> fft_ws -> bias -> np_bias_vec1[0] = 2. * b1;
    pfo -> fft_ws -> bias -> np_bias_vec1[1] = b2;
    pfo -> fft_ws -> bias -> np_bias_vec1[2] = 2. * bG2;
    pfo -> fft_ws -> bias -> np_bias_vec1[3] = 4. * (bG2 + 2./5. * btd);
    pfo -> fft_ws -> bias -> np_bias_vec1[4] = pow(b1, 2.);
    pfo -> fft_ws -> bias -> np_bias_vec1[5] = 0.5 * b1 * b2;
    pfo -> fft_ws -> bias -> np_bias_vec1[6] = b1 * bG2;
    pfo -> fft_ws -> bias -> p_bias_vec1[0] = 3. * b1;
    pfo -> fft_ws -> bias -> p_bias_vec1[1] = pow(b1, 2.);

    // 2-nd moment
    pfo -> fft_ws -> bias -> np_bias_vec2 = make_1Darray(6);
    pfo -> fft_ws -> bias -> p_bias_vec2  = make_1Darray(3);
    pfo -> fft_ws -> bias -> np_bias_vec2[0] = 2. * b1;
    pfo -> fft_ws -> bias -> np_bias_vec2[1] = b2;
    pfo -> fft_ws -> bias -> np_bias_vec2[2] = 2. * bG2;
    pfo -> fft_ws -> bias -> np_bias_vec2[3] = 2.;
    pfo -> fft_ws -> bias -> np_bias_vec2[4] = b1;
    pfo -> fft_ws -> bias -> np_bias_vec2[5] = pow(b1, 2.);
    pfo -> fft_ws -> bias -> p_bias_vec2[0] = 4. * b1;
    pfo -> fft_ws -> bias -> p_bias_vec2[1] = 6.;
    pfo -> fft_ws -> bias -> p_bias_vec2[2] = b1;

    // 3-rd moment
    pfo -> fft_ws -> bias -> np_bias_vec3 = make_1Darray(2);
    pfo -> fft_ws -> bias -> p_bias_vec3  = make_1Darray(1);
    pfo -> fft_ws -> bias -> np_bias_vec3[0] = 1.;
    pfo -> fft_ws -> bias -> np_bias_vec3[1] = b1;
    pfo -> fft_ws -> bias -> p_bias_vec3[0] = 2.;

    // 4-th moment
    pfo -> fft_ws -> bias -> np_bias_vec4 = make_1Darray(1);
    pfo -> fft_ws -> bias -> np_bias_vec4[0] = 1.;

    return _SUCCESS_;
}