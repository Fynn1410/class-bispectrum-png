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
                 double k, struct oneloop_fftlog_workspace *fft_ws, long SPLIT, double * pk_nl)
{

    // getting the linear power spectrum
    double pm_lin_IR = pm_IR_LO(pba, ppm, pfo, k, fft_ws -> z, SPLIT);
    double pm_lin    = Pk_dlnPk(pba, ppm, pfo, k, fft_ws -> z, LPOWER);
    double sigmav2   = fft_ws -> sigma_v2;

    // Distributing vectors for the biases and moment expansion factors, which will be applied on the loops

    double mu = fft_ws -> mu;
    double f  = fft_ws -> f;

    // 1-st moment      
    double *np_expans_vec1 = make_1Darray(7);
    double *p_expans_vec1  = make_1Darray(3);
    np_expans_vec1[0] = 2.* f * pow(mu, 2.);
    np_expans_vec1[1] = 2.* f * pow(mu, 2.);
    np_expans_vec1[2] = 2.* f * pow(mu, 2.);
    np_expans_vec1[3] = 2.* f * pow(mu, 2.);
    np_expans_vec1[4] = 2.* 2. * f * mu * k;
    np_expans_vec1[5] = 2.* 2. * f * mu * k;
    np_expans_vec1[6] = 2.* 2. * f * mu * k;
    p_expans_vec1[0] = 2.* f * pow(mu, 2.);
    p_expans_vec1[1] = 2.* 2. * f * mu * k;
    p_expans_vec1[2] = 2.* 2. * f * mu * k;

    // 2-nd moment
    double *np_expans_vec2 = make_1Darray(6);
    double *p_expans_vec2  = make_1Darray(4);
    np_expans_vec2[0] = pow(f * k * mu, 2.);
    np_expans_vec2[1] = pow(f * k * mu, 2.);
    np_expans_vec2[2] = pow(f * k * mu, 2.);
    np_expans_vec2[3] = pow(f, 2.) * pow(mu, 4.);
    np_expans_vec2[4] = pow(f, 2.) * pow(mu, 3.) * k;
    np_expans_vec2[5] = pow(f, 2.) * pow(k * mu, 4.);
    p_expans_vec2[0] = pow(f * k * mu, 2.);
    p_expans_vec2[1] = pow(f, 2.) * pow(mu, 4.);
    p_expans_vec2[2] = 4. * pow(f, 2.) * pow(mu, 3.) * k;
    p_expans_vec2[3] = 4. * pow(f, 2.) * pow(mu, 3.) * k;

    // 3-rd moment
    double *np_expans_vec3 = make_1Darray(2);
    double *p_expans_vec3  = make_1Darray(1);
    np_expans_vec3[0] = 2. * pow(f, 3.) * pow(mu, 4.) * pow(k, 2.);
    np_expans_vec3[1] = 2. * pow(f * k * mu, 3.);
    p_expans_vec3[0] = 2. * pow(f, 3.) * pow(mu, 4.) * pow(k, 2.);

    // 4-th moment
    double *np_expans_vec4 = make_1Darray(1);
    np_expans_vec4[0] = pow(f * k * mu, 4.) / 2.;


    // Integrals calculation
    double *np_rsd_0 = make_1Darray(6);
    double *p_rsd_0  = make_1Darray(2);
    rsd_0_FFTLog(fft_ws, k, np_rsd_0, p_rsd_0);

    double *np_rsd_1 = make_1Darray(7);
    double *p_rsd_1  = make_1Darray(3);
    rsd_1_FFTLog(fft_ws, k, mu, np_rsd_1, p_rsd_1);

    double *np_rsd_2 = make_1Darray(6);
    double *p_rsd_2  = make_1Darray(4);
    rsd_2_FFTLog(fft_ws, k, mu, np_rsd_2, p_rsd_2);

    double *np_rsd_3 = make_1Darray(2);
    double *p_rsd_3  = make_1Darray(1);
    rsd_3_FFTLog(fft_ws, k, mu, np_rsd_3, p_rsd_3);

    double *np_rsd_4 = make_1Darray(1);
    rsd_4_FFTLog(fft_ws, k, mu, np_rsd_4);

    // Bringing everything together
    double np_0  = 0.;
    double p_0   = 0.;
    dot(np_rsd_0, fft_ws -> bias -> np_bias_vec0, 6, &np_0);
    dot(p_rsd_0,  fft_ws -> bias -> p_bias_vec0,  2, &p_0);
    double ct_0   = - 2. * fft_ws -> b1 * (fft_ws -> R2 + fft_ws -> cs2 * fft_ws -> b1) * pow(k, 2.) * pm_lin_IR; // 0-th order counter term
    double P13_uv = - 61./105. * pm_lin_IR * pow(k, 2.) * sigmav2;
    double mom_0  = np_0 + pm_lin_IR * p_0 + pow(fft_ws -> b1, 2.) * (pm_lin_IR + P13_uv ) + ct_0;
    
    // FILE *fpa;
    // char file_name[50];
    // sprintf(file_name, "rsd_m_elements.txt");
    // fpa = fopen(file_name, "a");

    // fprintf(fpa, "%12.6e ", k);
    // fprintf(fpa, "%12.6e ", pm_lin);

    // fprintf(fpa, "%12.6e ", 2. * np_rsd_0[0]);
    // fprintf(fpa, "%12.6e ", 6. * pm_lin_IR * p_rsd_0[0]);
    // fprintf(fpa, "%12.6e ", P13_uv);
    // fprintf(fpa, "%12.6e ", - 2. * fft_ws -> cs2 * pow(k, 2.) * pm_lin_IR);
    // fprintf(fpa, "%12.6e ", (pm_lin_IR + 2. * np_rsd_0[0] + 6. * pm_lin_IR * p_rsd_0[0] + P13_uv - 2. * fft_ws -> cs2 * pow(k, 2.) * pm_lin_IR));
    // fprintf(fpa, "\n");
    // fclose(fpa);


    // FILE *fpa;
    // char file_name[50];
    // sprintf(file_name, "rsd_0_elements.txt");
    // fpa = fopen(file_name, "a");

    // fprintf(fpa, "%12.6e ", k);
    // fprintf(fpa, "%12.6e ", pm_lin);

    // double *np = make_1Darray(6);
    // double *p = make_1Darray(2);
    // vecmult(np_rsd_0, fft_ws -> bias -> np_bias_vec0, 6, np);
    // vecmult(p_rsd_0,  fft_ws -> bias -> p_bias_vec0,  2, p);
    // // fprintf(fpa, "%12.6e ", pow(fft_ws -> b1, 2.) *(pm_lin_IR + 2. * np_rsd_0[0] + 6. * pm_lin_IR * p_rsd_0[0] + P13_uv - 2. * fft_ws -> cs2 * pow(k, 2.) * pm_lin_IR));
    // fprintf(fpa, "%12.6e ", pow(fft_ws -> b1, 2.) * pm_lin_IR + np[0] + p[0] *  pm_lin_IR + P13_uv - 2. * fft_ws -> cs2 * pow(k, 2.) * pm_lin_IR);
    // fprintf(fpa, "%12.6e ", ct_0);
    // fprintf(fpa, "%12.6e ", np[1]);
    // fprintf(fpa, "%12.6e ", np[2]);
    // fprintf(fpa, "%12.6e ", np[3]);
    // fprintf(fpa, "%12.6e ", np[4]);
    // fprintf(fpa, "%12.6e ", np[5]);
    // fprintf(fpa, "%12.6e ", p[1] * pm_lin_IR);
    // fprintf(fpa, "%12.6e ", np_0 + pm_lin_IR * p_0);
    // fprintf(fpa, "%12.6e", mom_0);  
    // fprintf(fpa, "\n");
    // fclose(fpa);

    double *np_prefactors_1 = make_1Darray(7);
    double *p_prefactors_1  = make_1Darray(3);
    vecmult(fft_ws -> bias -> np_bias_vec1, np_expans_vec1, 7, np_prefactors_1);
    vecmult(fft_ws -> bias -> p_bias_vec1,  p_expans_vec1,  3, p_prefactors_1);
    double np_1  = 0.;
    double p_1   = 0.;
    dot(np_rsd_1, np_prefactors_1, 7, &np_1);
    dot(p_rsd_1,  p_prefactors_1,  3, &p_1);
    double mom_1 = np_1 + pm_lin_IR * p_1 + fft_ws -> b1 * pm_lin_IR * f * pow(mu,2.);

    double *np_prefactors_2 = make_1Darray(6);
    double *p_prefactors_2  = make_1Darray(4);
    vecmult(fft_ws -> bias -> np_bias_vec2, np_expans_vec2, 6, np_prefactors_2);
    vecmult(fft_ws -> bias -> p_bias_vec2,  p_expans_vec2,  4, p_prefactors_2);
    double np_2  = 0.;
    double p_2   = 0.;
    dot(np_rsd_2, np_prefactors_2, 6, &np_2);
    dot(p_rsd_2,  p_prefactors_2,  4, &p_2);
    double mom_2 = np_2 + pm_lin_IR * p_2 + pm_lin_IR * pow(f, 2.) * pow(mu, 4.) - pm_lin_IR * pow(f * mu * k * fft_ws -> b1, 2.) * sigmav2;

    double *np_prefactors_3 = make_1Darray(2);
    double *p_prefactors_3  = make_1Darray(1);
    vecmult(fft_ws -> bias -> np_bias_vec3, np_expans_vec3, 2, np_prefactors_3);
    vecmult(fft_ws -> bias -> p_bias_vec3,  p_expans_vec3,  1, p_prefactors_3);
    double np_3  = 0.;
    double p_3   = 0.;
    dot(np_rsd_3, np_prefactors_3, 2, &np_3);
    dot(p_rsd_3,  p_prefactors_3,  1, &p_3);
    double mom_3 = np_3 + pm_lin_IR * p_3 - 2 * pm_lin_IR * pow(f, 3.) * pow(mu, 4.) * pow(k, 2.) * fft_ws -> b1 * sigmav2;

    double *np_prefactors_4 = make_1Darray(1);
    vecmult(fft_ws -> bias -> np_bias_vec4, np_expans_vec4, 1, np_prefactors_4);
    double np_4  = 0.;
    dot(np_rsd_4, np_prefactors_4, 1, &np_4);
    double mom_4 = np_4 - pm_lin_IR * pow(f, 4.) * pow(mu, 6.) * pow(k, 2.) * sigmav2;

    // Summing over all moments
    double p_tot = mom_0 + mom_1 + mom_2 + mom_3 + mom_4;
    
    // double p_tot = mom_0;

    FILE *fpa2;
    char file_name2[50];
    sprintf(file_name2, "FFTLog_rsd.txt");
    fpa2 = fopen(file_name2, "a");
    fprintf(fpa2, "%e %e %e %e %e %e %e %e\n", mu, k, mom_0, mom_1, mom_2, mom_3, mom_4, p_tot);
    fclose(fpa2);

    *pk_nl = p_tot;
    return _SUCCESS_;
}

void FFTLog_rsd_init(struct background *pba, struct primordial *ppm, struct fourier *pfo, struct oneloop_fftlog_workspace *fft_ws){
    fft_ws -> fft_input  = (struct fft_struct *)malloc(sizeof(struct fft_struct));
    fft_ws -> fft_matrix = (struct fft_matrices *)malloc(sizeof(struct fft_matrices));

    int    N_FFTLog   = 256; // for a more precise calculation, gor for 256
    double kmin_fft_m = 1.e-8;
    double kmin_fft_g = 1.e-4;

    /* Setting the FFTLog parameters and calculating the etam and cmsym */
    fft_ws -> fft_input -> nfft 	  = N_FFTLog;
    fft_ws -> fft_input -> fft_bias_m = - 0.3;  
    fft_ws -> fft_input -> kmin_fft_m = kmin_fft_m;    
    fft_ws -> fft_input -> fft_bias_g = - 1.6; 
    fft_ws -> fft_input -> kmin_fft_g = kmin_fft_g;

    fft_ws -> fft_input -> etam_m  = make_1D_c_array(N_FFTLog + 1);
    fft_ws -> fft_input -> cmsym_m = make_1D_c_array(N_FFTLog + 1);
    fft_ws -> fft_input -> etam_g  = make_1D_c_array(N_FFTLog + 1);
    fft_ws -> fft_input -> cmsym_g = make_1D_c_array(N_FFTLog + 1);

    FFT_compute_coeff(pba, ppm, pfo, fft_ws -> z, fft_ws -> fft_input, 142L, HALO);
    FFT_compute_coeff(pba, ppm, pfo, fft_ws -> z, fft_ws -> fft_input, 142L, MATTER);

    /* Setting the bias vectors */

    double b1  = fft_ws -> b1;
    double b2  = fft_ws -> b2;
    double bG2 = fft_ws -> bG2;
    double btd = fft_ws -> btd;
    double cs2 = fft_ws -> cs2;
    double R2  = fft_ws -> R2;

    // 0-th moment
    fft_ws -> bias -> np_bias_vec0 = make_1Darray(6);
    fft_ws -> bias -> p_bias_vec0  = make_1Darray(2);
    fft_ws -> bias -> np_bias_vec0[0] = 2.* pow(b1, 2.);
    fft_ws -> bias -> np_bias_vec0[1] = 2.* b1 * b2;
    fft_ws -> bias -> np_bias_vec0[2] = 4.* b1 * bG2;
    fft_ws -> bias -> np_bias_vec0[3] =  0.5 * pow(b2, 2.);
    fft_ws -> bias -> np_bias_vec0[4] = 2.* pow(bG2, 2.);
    fft_ws -> bias -> np_bias_vec0[5] = 2.* b2 * bG2;
    fft_ws -> bias -> p_bias_vec0[0] = 6.* pow(b1, 2.);
    fft_ws -> bias -> p_bias_vec0[1] = 8. * b1 * (bG2 + 2./5. * btd);

    // 1-st moment
    fft_ws -> bias -> np_bias_vec1 = make_1Darray(7);
    fft_ws -> bias -> p_bias_vec1  = make_1Darray(3);
    fft_ws -> bias -> np_bias_vec1[0] = 2. * b1;
    fft_ws -> bias -> np_bias_vec1[1] = b2;
    fft_ws -> bias -> np_bias_vec1[2] = 2. * bG2;
    fft_ws -> bias -> np_bias_vec1[3] = 4. * (bG2 + 2./5. * btd);
    fft_ws -> bias -> np_bias_vec1[4] = pow(b1, 2.);
    fft_ws -> bias -> np_bias_vec1[5] = 0.5 * b1 * b2;
    fft_ws -> bias -> np_bias_vec1[6] = b1 * bG2;
    fft_ws -> bias -> p_bias_vec1[0] = 3. * b1;
    fft_ws -> bias -> p_bias_vec1[1] = pow(b1, 2.);
    fft_ws -> bias -> p_bias_vec1[2] = pow(b1, 2.);

    // 2-nd moment
    fft_ws -> bias -> np_bias_vec2 = make_1Darray(6);
    fft_ws -> bias -> p_bias_vec2  = make_1Darray(4);
    fft_ws -> bias -> np_bias_vec2[0] = 2. * b1;
    fft_ws -> bias -> np_bias_vec2[1] = b2;
    fft_ws -> bias -> np_bias_vec2[2] = 2. * bG2;
    fft_ws -> bias -> np_bias_vec2[3] = 2. ;
    fft_ws -> bias -> np_bias_vec2[4] = b1;
    fft_ws -> bias -> np_bias_vec2[5] = pow(b1, 2.);
    fft_ws -> bias -> p_bias_vec2[0] = 4. * b1;
    fft_ws -> bias -> p_bias_vec2[1] = 6.;
    fft_ws -> bias -> p_bias_vec2[2] = b1;
    fft_ws -> bias -> p_bias_vec2[3] = b1;

    // 3-rd moment
    fft_ws -> bias -> np_bias_vec3 = make_1Darray(2);
    fft_ws -> bias -> p_bias_vec3  = make_1Darray(1);
    fft_ws -> bias -> np_bias_vec3[0] = 1.;
    fft_ws -> bias -> np_bias_vec3[1] = b1;
    fft_ws -> bias -> p_bias_vec3[0] = 2.;

    // 4-th moment
    fft_ws -> bias -> np_bias_vec4 = make_1Darray(1);
    fft_ws -> bias -> np_bias_vec4[0] = 1.;

    /* Setting the matrices */

    /* 0-th moment */
    // FFTLog matrices non-propagator
    fft_ws -> fft_matrix -> I2200_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(I2200, fft_ws -> fft_input, 0., MATTER, fft_ws -> fft_matrix ->  I2200_mat);
    fft_ws -> fft_matrix -> Idelta200_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(Idelta200, fft_ws -> fft_input, 0., HALO, fft_ws -> fft_matrix -> Idelta200_mat);
    fft_ws -> fft_matrix -> IG200_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(IG200, fft_ws -> fft_input, 0., HALO, fft_ws -> fft_matrix -> IG200_mat);
    fft_ws -> fft_matrix -> Idelta2delta200_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(Idelta2delta200, fft_ws -> fft_input, 0., HALO, fft_ws -> fft_matrix -> Idelta2delta200_mat);
    fft_ws -> fft_matrix -> IG2G200_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(IG2G200, fft_ws -> fft_input, 0., HALO, fft_ws -> fft_matrix -> IG2G200_mat);
    fft_ws -> fft_matrix -> Idelta2G200_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(Idelta2G200, fft_ws -> fft_input, 0., HALO, fft_ws -> fft_matrix -> Idelta2G200_mat);

    // FFTLog matrices propagator
    fft_ws -> fft_matrix -> I1300_mat   = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(I1300, fft_ws -> fft_input, 0., MATTER, fft_ws -> fft_matrix ->  I1300_mat);
    fft_ws -> fft_matrix -> FG200_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(FG200, fft_ws -> fft_input, 0., HALO, fft_ws -> fft_matrix -> FG200_mat);

    /* 1-st moment */
    // FFTLog matrices (non-propagator)
    fft_ws -> fft_matrix -> I2201_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(I2201, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> I2201_mat);
    fft_ws -> fft_matrix -> Idelta201_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(Idelta201, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> Idelta201_mat);
    fft_ws -> fft_matrix -> IG201_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(IG201, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> IG201_mat);
    fft_ws -> fft_matrix -> FG201_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(FG201, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> FG201_mat);
    fft_ws -> fft_matrix -> J21101_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(J21101, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J21101_mat);
    fft_ws -> fft_matrix -> Jdelta201_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(Jdelta201, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> Jdelta201_mat);
    fft_ws -> fft_matrix -> JG201_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(JG201, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> JG201_mat);

    // FFTLog matrices (propagator)
    fft_ws -> fft_matrix -> I1301_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(I1301, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> I1301_mat);
    fft_ws -> fft_matrix -> J12101_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(J12101, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J12101_mat);
    fft_ws -> fft_matrix -> J11201_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(J11201, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J11201_mat);

    /* 2-nd moment */
    // FFTLog matrices (non-propagator)
    fft_ws -> fft_matrix -> J21102x_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(J21102x, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J21102x_mat);
    fft_ws -> fft_matrix -> J21102y_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(J21102y, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J21102y_mat);
    fft_ws -> fft_matrix -> Jdelta202x_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(Jdelta202x, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> Jdelta202x_mat);
    fft_ws -> fft_matrix -> Jdelta202y_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(Jdelta202y, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> Jdelta202y_mat);
    fft_ws -> fft_matrix -> JG202x_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(JG202x, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> JG202x_mat);
    fft_ws -> fft_matrix -> JG202y_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(JG202y, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> JG202y_mat);
    fft_ws -> fft_matrix -> I2211_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(I2211, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> I2211_mat);
    fft_ws -> fft_matrix -> J21111_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(J21111, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J21111_mat);
    fft_ws -> fft_matrix -> N11x_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(N11x, fft_ws -> fft_input, fft_ws -> mu, HALO,  fft_ws -> fft_matrix -> N11x_mat);
    fft_ws -> fft_matrix -> N11y_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(N11y, fft_ws -> fft_input, fft_ws -> mu, HALO,  fft_ws -> fft_matrix -> N11y_mat);

    // FFTLog matrices (propagator)
    fft_ws -> fft_matrix -> J12102x_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(J12102x, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J12102x_mat);
    fft_ws -> fft_matrix -> J12102y_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(J12102y, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J12102y_mat);
    fft_ws -> fft_matrix -> I1311_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(I1311, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> I1311_mat);
    fft_ws -> fft_matrix -> J12111_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(J12111, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J12111_mat);
    fft_ws -> fft_matrix -> J11211_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(J11211, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J11211_mat);

    /* 3-rd moment */
    // FFTLog matrices (non-propagator)
    fft_ws -> fft_matrix -> J21112x_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(J21112x, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J21112x_mat);
    fft_ws -> fft_matrix -> J21112y_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(J21112y, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J21112y_mat);
    fft_ws -> fft_matrix -> N12x_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(N12x, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> N12x_mat);
    fft_ws -> fft_matrix -> N12y_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(N12y, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> N12y_mat);
    
    // FFTLog matrices (propagator)
    fft_ws -> fft_matrix -> J12112x_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(J12112x, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J12112x_mat);
    fft_ws -> fft_matrix -> J12112y_mat = make_1D_c_array(N_FFTLog+1);
    p_mat_fill(J12112y, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> J12112y_mat);

    /* 4-th moment */
    // FFTLog matrices (non-propagator)
    fft_ws -> fft_matrix -> N22x_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(N22x, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> N22x_mat);
    fft_ws -> fft_matrix -> N22y_mat   = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(N22y, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> N22y_mat);
    fft_ws -> fft_matrix -> N22z_mat = make_2D_c_array(N_FFTLog+1, N_FFTLog+1);
    np_mat_fill(N22z, fft_ws -> fft_input, fft_ws -> mu, HALO, fft_ws -> fft_matrix -> N22z_mat);
}
