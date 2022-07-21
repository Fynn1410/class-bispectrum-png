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

    double k = pfo->k[index_k];

    // getting the linear power spectrum
    double pm_lin_IR = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
    double pm_lin    = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    double sigmav0   = pfo -> fft_ws -> sigma_v0;
    double sigmav2   = pfo -> fft_ws -> sigma_v2;

    // Integrals calculation
    double *np_rsd_0 = make_1Darray(6);
    double *p_rsd_0  = make_1Darray(2);
    rsd_0_FFTLog(pfo -> fft_ws, k, np_rsd_0, p_rsd_0);

    double *np_rsd_1 = make_1Darray(7);
    double *p_rsd_1  = make_1Darray(2);
    rsd_1_FFTLog(pfo -> fft_ws, k, mu, np_rsd_1, p_rsd_1);

    double *np_rsd_2 = make_1Darray(6);
    double *p_rsd_2  = make_1Darray(3);
    rsd_2_FFTLog(pfo -> fft_ws, k, mu, np_rsd_2, p_rsd_2);

    double *np_rsd_3 = make_1Darray(2);
    double *p_rsd_3  = make_1Darray(1);
    rsd_3_FFTLog(pfo -> fft_ws, k, mu, np_rsd_3, p_rsd_3);

    double *np_rsd_4 = make_1Darray(1);
    rsd_4_FFTLog(pfo -> fft_ws, k, mu, np_rsd_4);

    // Bringing everything together
    double np_0  = 0.;
    double p_0   = 0.;
    dot(np_rsd_0, pfo -> fft_ws -> bias -> np_bias_vec0, 6, &np_0);
    dot(p_rsd_0,  pfo -> fft_ws -> bias -> p_bias_vec0,  2, &p_0);
    double ct_0   = - 2. * pfo -> fft_ws -> b1 * (pfo -> fft_ws -> R2 + pfo -> fft_ws -> cs2 * pfo -> fft_ws -> b1) * pow(k, 2.) * pm_lin_IR; // 0-th order counter term
    double P13_uv = - 61./105. * pm_lin_IR * pow(k, 2.) * sigmav2;
    double p_nowiggle = pm_nowiggle(pba, ppm, pfo, k, z, 1.e-5, 0, SPLIT);
    double sup        = exp(-k * k * pfo -> fft_ws -> sigma_2_IR);
    double p_wiggle   = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER) - p_nowiggle;
    double mom_0  = np_0 + pm_lin_IR * p_0 + pow(pfo -> fft_ws -> b1, 2.) * (p_nowiggle + sup * p_wiggle * (1. + k * k * pfo -> fft_ws -> sigma_2_IR) + P13_uv) + ct_0;
    
    // FILE *fpa;
    // char file_name[50];
    // sprintf(file_name, "data/rsd_m_elements.txt");
    // fpa = fopen(file_name, "a");

    // fprintf(fpa, "%12.6e ", k);
    // fprintf(fpa, "%12.6e ", p_nowiggle + sup * p_wiggle * (1. + k * k * pfo -> fft_ws -> sigma_2_IR));

    // fprintf(fpa, "%12.6e ", 2. * np_rsd_0[0]);
    // fprintf(fpa, "%12.6e ", 6. * pm_lin_IR * p_rsd_0[0]);
    // fprintf(fpa, "%12.6e ", P13_uv);
    // fprintf(fpa, "%12.6e ", P13_uv + 6. * pm_lin_IR * p_rsd_0[0]);
    // fprintf(fpa, "%12.6e ", - 2. * pfo -> fft_ws -> cs2 * pow(k, 2.) * pm_lin_IR);
    // fprintf(fpa, "%12.6e ", (p_nowiggle + sup * p_wiggle * (1. + k * k * pfo -> fft_ws -> sigma_2_IR) + 2. * np_rsd_0[0] + 6. * pm_lin_IR * p_rsd_0[0] + P13_uv - 2. * pfo -> fft_ws -> cs2 * pow(k, 2.) * pm_lin_IR));
    // fprintf(fpa, "\n");
    // fclose(fpa);


    // FILE *fpa0;
    // char file_name0[50];
    // sprintf(file_name0, "data/rsd_0.txt");
    // fpa0 = fopen(file_name0, "a");
    // double *np0 = make_1Darray(6);
    // double *p0 = make_1Darray(2);
    // vecmult(np_rsd_0, pfo -> fft_ws -> bias -> np_bias_vec0, 6, np0);
    // vecmult(p_rsd_0,  pfo -> fft_ws -> bias -> p_bias_vec0,  2, p0);
    
    // fprintf(fpa0, "%12.6e ", k);
    // fprintf(fpa0, "%12.6e ", pm_lin);
    // // fprintf(fpa0, "%12.6e ", pow(pfo -> fft_ws -> b1, 2.) *(pm_lin_IR + 2. * np_rsd_0[0] + 6. * pm_lin_IR * p_rsd_0[0] + P13_uv - 2. * pfo -> fft_ws -> cs2 * pow(k, 2.) * pm_lin_IR));
    // fprintf(fpa0, "%12.6e ", pow(pfo -> fft_ws -> b1, 2.) *(p_nowiggle + sup * p_wiggle * (1. + k * k * pfo -> fft_ws -> sigma_2_IR)) + np0[0] + p0[0] *  pm_lin_IR + pow(pfo -> fft_ws -> b1, 2.) * P13_uv + ct_0);
    // fprintf(fpa0, "%12.6e ", ct_0);
    // fprintf(fpa0, "%12.6e ", pow(pfo -> fft_ws -> b1, 2.) * P13_uv);
    // fprintf(fpa0, "%12.6e ", np0[0]);
    // fprintf(fpa0, "%12.6e ", np0[1]);
    // fprintf(fpa0, "%12.6e ", np0[2]);
    // fprintf(fpa0, "%12.6e ", np0[3]);
    // fprintf(fpa0, "%12.6e ", np0[4]);
    // fprintf(fpa0, "%12.6e ", np0[5]);
    // fprintf(fpa0, "%12.6e ", p0[0] * pm_lin_IR);
    // fprintf(fpa0, "%12.6e ", p0[1] * pm_lin_IR);
    // fprintf(fpa0, "%12.6e ", np_0 + pm_lin_IR * p_0);
    // fprintf(fpa0, "%12.6e", mom_0);  
    // fprintf(fpa0, "\n");
    // fclose(fpa0);

    double *np_prefactors_1 = make_1Darray(7);
    double *p_prefactors_1  = make_1Darray(2);
    vecmult(pfo -> fft_ws -> bias -> np_bias_vec1, np_expans_vec1, 7, np_prefactors_1);
    vecmult(pfo -> fft_ws -> bias -> p_bias_vec1,  p_expans_vec1,  2, p_prefactors_1);
    double np_1  = 0.;
    double p_1   = 0.;
    dot(np_rsd_1, np_prefactors_1, 7, &np_1);
    dot(p_rsd_1,  p_prefactors_1,  2, &p_1);
    double J11201    = -0.5 * k * mu * (sigmav2 + sigmav0 / (3. * k * k)) * pm_lin_IR;
    double I1301_UV  = -25./63. * pow(k, 2.) * sigmav2 * pm_lin_IR; 
    double J12101_UV = 1./6. * sigmav0 * pm_lin_IR / k;
    // fprintf(stderr, "J11201 = %e\nI1301_UV = %e\nJ12101_UV = %e\n", J11201, I1301_UV, J12101_UV);
    double mom_1 = np_1 + pm_lin_IR * p_1 + pfo -> fft_ws -> b1 * f * pow(mu,2.) * pm_lin_IR + 2. * 2. * f * mu * k * pow(pfo -> fft_ws -> b1, 2.) * J11201 + 2.* f * pow(mu, 2.) * 3. * pfo -> fft_ws -> b1 * I1301_UV + 2. * 2. * f * mu * k * pow(pfo -> fft_ws -> b1, 2.) * J12101_UV;

    // FILE *fpa1;
    // char file_name1[50];
    // sprintf(file_name1, "data/rsd_1.txt");
    // fpa1 = fopen(file_name1, "a");

    // double *np1 = make_1Darray(7);
    // double *p1 = make_1Darray(2);
    // vecmult(np_rsd_1, np_prefactors_1, 7, np1);
    // vecmult(p_rsd_1,  p_prefactors_1,  2, p1);

    // fprintf(fpa1, "%12.6e ", k);
    // fprintf(fpa1, "%12.6e ", pm_lin);
    // fprintf(fpa1, "%12.6e ", np1[0]);
    // fprintf(fpa1, "%12.6e ", np1[1]);
    // fprintf(fpa1, "%12.6e ", np1[2]);
    // fprintf(fpa1, "%12.6e ", np1[3]);
    // fprintf(fpa1, "%12.6e ", np1[4]);
    // fprintf(fpa1, "%12.6e ", np1[5]);
    // fprintf(fpa1, "%12.6e ", np1[6]);
    // fprintf(fpa1, "%12.6e ", p1[0] * pm_lin_IR);
    // fprintf(fpa1, "%12.6e ", p1[1] * pm_lin_IR);
    // fprintf(fpa1, "%12.6e ", 2. * 2. * f * mu * k * pow(pfo -> fft_ws -> b1, 2.) * J11201);
    // fprintf(fpa1, "%12.6e ", 2.* f * pow(mu, 2.) * 3. * pfo -> fft_ws -> b1 * I1301_UV);
    // fprintf(fpa1, "%12.6e ", 2. * 2. * f * mu * k * pow(pfo -> fft_ws -> b1, 2.) * J12101_UV);
    // fprintf(fpa1, "%12.6e", mom_1);  
    // fprintf(fpa1, "\n");
    // fclose(fpa1);

    // // Velociraptor
    // fprintf(fpa1, "%12.6e ", k);
    // fprintf(fpa1, "%12.6e ", pm_lin);
    // fprintf(fpa1, "%12.6e ", np1[0] + p1[0] * pm_lin_IR);
    // fprintf(fpa1, "%12.6e ", np1[4] + p1[1]  * pm_lin_IR);
    // fprintf(fpa1, "%12.6e ", np1[1]);
    // fprintf(fpa1, "%12.6e ", np1[5]);
    // fprintf(fpa1, "%12.6e ", np1[2]);
    // fprintf(fpa1, "%12.6e ", np1[6]);
    // fprintf(fpa1, "%12.6e ", np1[3]);
    // fprintf(fpa1, "%12.6e", mom_1);  
    // fprintf(fpa1, "\n");
    // fclose(fpa1);

    double *np_prefactors_2 = make_1Darray(6);
    double *p_prefactors_2  = make_1Darray(3);
    vecmult(pfo -> fft_ws -> bias -> np_bias_vec2, np_expans_vec2, 6, np_prefactors_2);
    vecmult(pfo -> fft_ws -> bias -> p_bias_vec2,  p_expans_vec2,  3, p_prefactors_2);
    double np_2  = 0.;
    double p_2   = 0.;
    dot(np_rsd_2, np_prefactors_2, 6, &np_2);
    dot(p_rsd_2,  p_prefactors_2,  3, &p_2);
    double J12102_UV = -0.5 * mu * sigmav2 * pm_lin_IR;
    double I3111_UV  = -0.3 * k * k * sigmav2 * pm_lin_IR;
    double J12111_UV = 1./6. * sigmav0 * pm_lin_IR / k;
    double J11211    = -0.5 * k * mu * (sigmav2 + sigmav0 / (3. * k * k)) * pm_lin_IR;
    double mom_2 = np_2 + pm_lin_IR * p_2 + pm_lin_IR * pow(f, 2.) * pow(mu, 4.) + pm_lin_IR * pow(f * mu * k * pfo -> fft_ws -> b1, 2.) * sigmav2 + (4. * pow(f, 2.) * pow(mu, 3.) * k) * pfo -> fft_ws -> b1 * (J11211+J12111_UV) + pow(f*k*mu, 2.) * J12102_UV + 6. * pow(mu * mu * f, 2.) * I3111_UV;

    // FILE *fpa2;
    // char file_name2[50];
    // sprintf(file_name2, "data/rsd_2.txt");
    // fpa2 = fopen(file_name2, "a");

    // double *np2 = make_1Darray(6);
    // double *p2 = make_1Darray(3);
    // vecmult(np_rsd_2, np_prefactors_2, 6, np2);
    // vecmult(p_rsd_2,  p_prefactors_2,  3, p2);

    // fprintf(fpa2, "%12.6e ", k);
    // fprintf(fpa2, "%12.6e ", pm_lin);
    // fprintf(fpa2, "%12.6e ", np2[0]);
    // fprintf(fpa2, "%12.6e ", np2[1]);
    // fprintf(fpa2, "%12.6e ", np2[2]);
    // fprintf(fpa2, "%12.6e ", np2[3]);
    // fprintf(fpa2, "%12.6e ", np2[4]);
    // fprintf(fpa2, "%12.6e ", np2[5]);
    // fprintf(fpa2, "%12.6e ", p2[0] * pm_lin_IR);
    // fprintf(fpa2, "%12.6e ", p2[1] * pm_lin_IR);
    // fprintf(fpa2, "%12.6e ", p2[2] * pm_lin_IR);
    // fprintf(fpa2, "%12.6e ", 4. * pow(f, 2.) * pow(mu, 3.) * k * pfo -> fft_ws -> b1 * J11211);
    // fprintf(fpa2, "%12.6e ", 4. * pow(f, 2.) * pow(mu, 3.) * k * pfo -> fft_ws -> b1 * J12111_UV);  
    // fprintf(fpa2, "%12.6e ", pow(f*k*mu, 2.) * J12102_UV);  
    // fprintf(fpa2, "%12.6e ", 6. * pow(mu * mu * f, 2.) * I3111_UV);    
    // fprintf(fpa2, "%12.6e", mom_2);  
    // fprintf(fpa2, "\n");
    // fclose(fpa2);

    // // Velociraptor
    // fprintf(fpa2, "%12.6e ", k);
    // fprintf(fpa2, "%12.6e ", pm_lin);
    // fprintf(fpa2, "%12.6e ", pm_lin_IR * pow(f, 2.) * pow(mu, 4.));
    // fprintf(fpa2, "%12.6e ", np2[0] + p2[0] * pm_lin_IR);
    // fprintf(fpa2, "%12.6e ", np2[5] + pm_lin_IR * pow(f * mu * k * pfo -> fft_ws -> b1, 2.) * sigmav2);
    // fprintf(fpa2, "%12.6e ", np2[1]);
    // fprintf(fpa2, "%12.6e ", 0.);
    // fprintf(fpa2, "%12.6e ", np2[3] + p2[1] * pm_lin_IR);
    // fprintf(fpa2, "%12.6e ", np2[4] + p2[2] * pm_lin_IR);
    // fprintf(fpa2, "%12.6e ", 0.);
    // fprintf(fpa2, "%12.6e ", 0.);
    // fprintf(fpa2, "%12.6e ", np2[2]);
    // fprintf(fpa2, "%12.6e", mom_2);  
    // fprintf(fpa2, "\n");
    // fclose(fpa2);

    double *np_prefactors_3 = make_1Darray(2);
    double *p_prefactors_3  = make_1Darray(1);
    vecmult(pfo -> fft_ws -> bias -> np_bias_vec3, np_expans_vec3, 2, np_prefactors_3);
    vecmult(pfo -> fft_ws -> bias -> p_bias_vec3,  p_expans_vec3,  1, p_prefactors_3);
    double np_3  = 0.;
    double p_3   = 0.;
    dot(np_rsd_3, np_prefactors_3, 2, &np_3);
    dot(p_rsd_3,  p_prefactors_3,  1, &p_3);
    double J12112_IR = -0.5 * sigmav2 * pm_lin_IR / mu;
    double mom_3 = np_3 + pm_lin_IR * p_3 - 2. * pm_lin_IR * pow(f, 3.) * pow(mu, 4.) * pow(k, 2.) * pfo -> fft_ws -> b1 * sigmav2 + 4. * pow(f, 3.) * pow(mu, 4.) * pow(k, 2.) * J12112_IR;

    // FILE *fpa3;
    // char file_name3[50];
    // sprintf(file_name3, "data/rsd_3.txt");
    // fpa3 = fopen(file_name3, "a");

    // double *np3 = make_1Darray(2);
    // double *p3 = make_1Darray(1);
    // vecmult(np_rsd_3, np_prefactors_3, 2, np3);
    // vecmult(p_rsd_3,  p_prefactors_3,  1, p3);

    // fprintf(fpa3, "%12.6e ", k);
    // fprintf(fpa3, "%12.6e ", pm_lin);
    // fprintf(fpa3, "%12.6e ", np3[0]);
    // fprintf(fpa3, "%12.6e ", np3[1]);
    // fprintf(fpa3, "%12.6e ", p3[0] * pm_lin_IR);
    // fprintf(fpa3, "%12.6e ", 4. * pow(f, 3.) * pow(mu, 4.) * pow(k, 2.) * J12112_IR);
    // fprintf(fpa3, "%12.6e", mom_3);  
    // fprintf(fpa3, "\n");
    // fclose(fpa3);

    // Velociraptor
    // fprintf(fpa3, "%12.6e ", k);
    // fprintf(fpa3, "%12.6e ", pm_lin);
    // fprintf(fpa3, "%12.6e ", 0.);
    // fprintf(fpa3, "%12.6e ", np3[1]);
    // fprintf(fpa3, "%12.6e ", np3[0] + p3[0] * pm_lin_IR);
    // fprintf(fpa3, "%12.6e ", - 2. * pm_lin_IR * pow(f, 3.) * pow(mu, 4.) * pow(k, 2.) * pfo -> fft_ws -> b1 * sigmav2);
    // fprintf(fpa3, "%12.6e", mom_3);  
    // fprintf(fpa3, "\n");
    // fclose(fpa3);

    double *np_prefactors_4 = make_1Darray(1);
    vecmult(pfo -> fft_ws -> bias -> np_bias_vec4, np_expans_vec4, 1, np_prefactors_4);
    double np_4  = 0.;
    dot(np_rsd_4, np_prefactors_4, 1, &np_4);
    double mom_4 = np_4 - pm_lin_IR * pow(f, 4.) * pow(mu, 6.) * pow(k, 2.) * sigmav2;

    // FILE *fpa4;
    // char file_name4[50];
    // sprintf(file_name4, "data/rsd_4.txt");
    // fpa4 = fopen(file_name4, "a");

    // double *np4 = make_1Darray(1);
    // vecmult(np_rsd_4, np_prefactors_4, 1, np4);

    // fprintf(fpa4, "%12.6e ", k);
    // fprintf(fpa4, "%12.6e ", pm_lin);
    // fprintf(fpa4, "%12.6e ", np4[0]);
    // fprintf(fpa4, "%12.6e", mom_4);  
    // fprintf(fpa4, "\n");
    // fclose(fpa4);

    // Summing over all moments
    double p_tot = mom_0 + mom_1 + mom_2 + mom_3 + mom_4;
    
    // double p_tot = mom_0;

    // FILE *fpa5;
    // char file_name5[50];
    // sprintf(file_name5, "data/FFTLog_rsd.txt");
    // fpa5 = fopen(file_name5, "a");
    // fprintf(fpa5, "%e %e %e %e %e %e %e %e\n", k, mom_0, mom_1, mom_2, mom_3, mom_4, p_tot, mu);
    // fclose(fpa5);

    *pk_nl = p_tot;
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