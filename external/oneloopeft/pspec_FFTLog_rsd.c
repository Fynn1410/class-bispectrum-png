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

    // Storing IR-resummed Plin in pfo
    double pm_lin = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    pfo -> pk_halo_rsd_nl[linear]    -> Plin[index_k] = pm_lin;
    double pm_lin_nw = pm_nowiggle(pba, ppm, pfo, k, z, 1.e-4, 0, SPLIT);
    pfo -> pk_halo_rsd_nl[no_wiggle] -> Plin[index_k] = pm_lin_nw;

    // Calculation of the Loop-Integrals for the full linear power spectrum and the no-wiggle linear power spectrum
    double plin;
    // for (int idx = lin; idx <= lin; idx++){
    for (int idx = lin; idx <= no_wiggle; idx++){
        rsd_0_FFTLog(pfo, idx, index_k, pfo -> pk_halo_rsd_nl[idx] -> Plin[index_k]);
        rsd_1_FFTLog(pfo, idx, index_k, pfo -> pk_halo_rsd_nl[idx] -> Plin[index_k]);
        rsd_2_FFTLog(pfo, idx, index_k, pfo -> pk_halo_rsd_nl[idx] -> Plin[index_k]);
        rsd_3_FFTLog(pfo, idx, index_k, pfo -> pk_halo_rsd_nl[idx] -> Plin[index_k]);
        rsd_4_FFTLog(pfo, idx, index_k, pfo -> pk_halo_rsd_nl[idx] -> Plin[index_k]);
    }

    return _SUCCESS_;
}

int RSD_IR_Ressummed(struct fourier *pfo, struct background *pba, int index_k, double z, double mu, double * result)
{
    double k = pfo->k[index_k];

    // Biases
    double b1  =  2.0;
    double b2  = -1.0;
    double bG2 =  0.1;
    double btd = -0.1;
    double R2  =  5.0;

    // Counter terms
    double cs2 =  0.2;
    double c00  = 0.;
    double c10  = 0.;
    double c20  = 0.;
    double c22  = 0.;
    double c30  = 0.;
    double c32  = 0.;
    double c42  = 0.;

    // double * pvecback;
    // int last_index;
    // class_alloc(pvecback,pba->bg_size*sizeof(double),pfo->error_message);
    // class_call(background_at_z(pba, z, short_info, inter_normal, &last_index, pvecback),
    //                             pba->error_message,
    //                             pfo->error_message);

    // double D = pvecback[pba->index_bg_D]; // Growth factor
    // double f = pvecback[pba->index_bg_f]; // logarithmic growth
    // // fprintf(stderr, "RSD at z = %g, D(z) = %g, f(z) = %g\n",z,D,f);

    double D = 1.;
    double f = pow(0.3,5./9.);

    double D2 = pow(D,2.); 
    double D4 = pow(D,4.); 


    double sigma_v2 = pfo -> fft_ws -> sigma_v2;
    double sigma_2_IR = pfo -> fft_ws -> sigma_2_IR;
    double del_sigma_2_IR = pfo -> fft_ws -> del_sigma_2_IR;

    double sigma_tot = (1. + f*pow(mu,2.)*(2 + f))*sigma_2_IR + pow(f*mu,2.)*(pow(mu,2.) - 1.)*del_sigma_2_IR;
    double suppression = exp(-pow(k,2.)*sigma_tot);

    //0-th moment
    double * Plin = make_1Darray(2);
    double * P_mm = make_1Darray(2);
    double * I2200 = make_1Darray(2);
    double * Idelta200 = make_1Darray(2);
    double * IG200 = make_1Darray(2);
    double * Idelta2delta200 = make_1Darray(2);
    double * IG2G200 = make_1Darray(2);
    double * Idelta2G200 = make_1Darray(2);
    double * I1300 = make_1Darray(2);
    double * FG200 = make_1Darray(2);
    double * IR2 = make_1Darray(2);
    double * Moment_0 = make_1Darray(2);

    //1-st moment
    double * I2201 = make_1Darray(2);
    double * Idelta201 = make_1Darray(2);
    double * IG201 = make_1Darray(2);
    double * FG201 = make_1Darray(2);
    double * J21101 = make_1Darray(2);
    double * Jdelta201 = make_1Darray(2);
    double * JG201 = make_1Darray(2);
    double * I1301p3101 = make_1Darray(2);
    double * J12101 = make_1Darray(2);
    double * J11201 = make_1Darray(2);
    double * Moment_1 = make_1Darray(2);


    //2-nd moment
    double * J21102x = make_1Darray(2);
    double * J21102y = make_1Darray(2);
    double * J21102 = make_1Darray(2);
    double * Jdelta202x = make_1Darray(2);
    double * Jdelta202y = make_1Darray(2);
    double * Jdelta202 = make_1Darray(2);
    double * JG202x = make_1Darray(2);
    double * JG202y = make_1Darray(2);
    double * JG202 = make_1Darray(2);
    double * I2211 = make_1Darray(2);
    double * J21111 = make_1Darray(2);
    double * N11x = make_1Darray(2);
    double * N11y = make_1Darray(2);
    double * N11 = make_1Darray(2);
    double * J12102x = make_1Darray(2);
    double * J12102y = make_1Darray(2);
    double * J12102 = make_1Darray(2);
    double * I1311 = make_1Darray(2);
    double * J12111 = make_1Darray(2);
    double * J11211 = make_1Darray(2);
    double * Moment_2 = make_1Darray(2);


    //3-rd moment
    double * J21112x = make_1Darray(2);
    double * J21112y = make_1Darray(2);
    double * J21112 = make_1Darray(2);
    double * N12x = make_1Darray(2);
    double * N12y = make_1Darray(2);
    double * N12 = make_1Darray(2);
    double * J12112x = make_1Darray(2);
    double * J12112y = make_1Darray(2);
    double * J12112 = make_1Darray(2);
    double * Moment_3 = make_1Darray(2);

    //4-th moment
    double * N22x = make_1Darray(2);
    double * N22y = make_1Darray(2);
    double * N22z = make_1Darray(2);
    double * N22  = make_1Darray(2);
    double * Moment_4 = make_1Darray(2);

    //RSD expansion
    double * RSD = make_1Darray(2);

    // Append the prefactors in terms of mu dependencies from the integrand
    for (int idx = lin; idx <= no_wiggle; idx++){
        //0-th moment
        Plin[idx]            = pfo -> pk_halo_rsd_nl[idx] -> Plin[index_k];
        I2200[idx]           = pfo -> pk_halo_rsd_nl[idx] -> I2200[index_k];
        Idelta200[idx]       = pfo -> pk_halo_rsd_nl[idx] -> Idelta200[index_k];
        IG200[idx]           = pfo -> pk_halo_rsd_nl[idx] -> IG200[index_k];
        Idelta2delta200[idx] = pfo -> pk_halo_rsd_nl[idx] -> Idelta2delta200[index_k];
        IG2G200[idx]         = pfo -> pk_halo_rsd_nl[idx] -> IG2G200[index_k];
        Idelta2G200[idx]     = pfo -> pk_halo_rsd_nl[idx] -> Idelta2G200[index_k];
        I1300[idx]           = pfo -> pk_halo_rsd_nl[idx] -> I1300[index_k];
        FG200[idx]           = pfo -> pk_halo_rsd_nl[idx] -> FG200[index_k];
        IR2[idx]             = pfo -> pk_halo_rsd_nl[idx] -> IR2[index_k];
        
        Moment_0[idx] = (pow(b1,2.)*Plin[idx] - 2.*cs2*pow(k,2.)*Plin[idx] + b1*R2*IR2[idx]) *D2 \
                        + (pow(b1,2.)*(2.*I2200[idx] + 6.*I1300[idx]) + 2.*b1*b2*Idelta200[idx] + 4.*b1*bG2*IG200[idx] \
                        + 0.5*pow(b2,2.)*Idelta2delta200[idx] + 2.*pow(bG2,2.)*IG2G200[idx] + 8.*b1*(bG2 + 0.4*btd)*FG200[idx]) *D4;

        //1-st moment
        I2201[idx]      = pfo -> pk_halo_rsd_nl[idx] -> I2201[index_k];
        Idelta201[idx]  = pfo -> pk_halo_rsd_nl[idx] -> Idelta201[index_k];
        IG201[idx]      = pfo -> pk_halo_rsd_nl[idx] -> IG201[index_k];
        FG201[idx]      = pfo -> pk_halo_rsd_nl[idx] -> FG201[index_k];
        J21101[idx]     = pfo -> pk_halo_rsd_nl[idx] -> J21101[index_k] * mu;
        Jdelta201[idx]  = pfo -> pk_halo_rsd_nl[idx] -> Jdelta201[index_k] * mu;
        JG201[idx]      = pfo -> pk_halo_rsd_nl[idx] -> JG201[index_k] * mu;
        I1301p3101[idx] = pfo -> pk_halo_rsd_nl[idx] -> I1301p3101[index_k];
        J12101[idx]     = pfo -> pk_halo_rsd_nl[idx] -> J12101[index_k] * mu;
        J11201[idx]     = pfo -> pk_halo_rsd_nl[idx] -> J11201[index_k] * mu;

        Moment_1[idx] =   2. * (f*mu/k) *b1*Plin[idx] *D2\
                        + 2. * (f*mu/k) * (2.*b1*I2201[idx] + 3.*b1*I1301p3101[idx] + b2*Idelta201[idx] + 2.*bG2*IG201[idx] + 4.*(bG2 + 0.4*btd)*FG201[idx]) *D4\
                        + 2. * (2.*f)    * (pow(b1,2.)*(J12101[idx] + J11201[idx] + J21101[idx]) + 0.5*b1*b2*Jdelta201[idx] + b1*bG2*JG201[idx]) *D4\
                        + c10*f*mu*k*Plin[idx] *D2;

        //2-nd moment
        J21102x[idx]    = pfo -> pk_halo_rsd_nl[idx] -> J21102x[index_k];
        J21102y[idx]    = pfo -> pk_halo_rsd_nl[idx] -> J21102y[index_k];
        J21102[idx]     = (J21102x[idx] + J21102y[idx] * pow(mu,2.));
        Jdelta202x[idx] = pfo -> pk_halo_rsd_nl[idx] -> Jdelta202x[index_k];
        Jdelta202y[idx] = pfo -> pk_halo_rsd_nl[idx] -> Jdelta202y[index_k];
        Jdelta202[idx]  = (Jdelta202x[idx] + Jdelta202y[idx] * pow(mu,2.));
        JG202x[idx]     = pfo -> pk_halo_rsd_nl[idx] -> JG202x[index_k];
        JG202y[idx]     = pfo -> pk_halo_rsd_nl[idx] -> JG202y[index_k];
        JG202[idx]      = (JG202x[idx] + JG202y[idx] * pow(mu,2.));
        I2211[idx]      = pfo -> pk_halo_rsd_nl[idx] -> I2211[index_k];
        J21111[idx]     = pfo -> pk_halo_rsd_nl[idx] -> J21111[index_k] * mu;
        N11x[idx]       = pfo -> pk_halo_rsd_nl[idx] -> N11x[index_k];
        N11y[idx]       = pfo -> pk_halo_rsd_nl[idx] -> N11y[index_k];
        N11[idx]        = (N11x[idx] + N11y[idx] * pow(mu,2.));
        J12102x[idx]    = pfo -> pk_halo_rsd_nl[idx] -> J12102x[index_k];
        J12102y[idx]    = pfo -> pk_halo_rsd_nl[idx] -> J12102y[index_k];
        J12102[idx]     = (J12102x[idx] + J12102y[idx] * pow(mu,2.));
        I1311[idx]      = pfo -> pk_halo_rsd_nl[idx] -> I1311[index_k];
        J12111[idx]     = pfo -> pk_halo_rsd_nl[idx] -> J12111[index_k] * mu;
        J11211[idx]     = pfo -> pk_halo_rsd_nl[idx] -> J11211[index_k] * mu;

        Moment_2[idx] =   2. * pow(f*mu/k,2.) * Plin[idx] *D2\
                        + 2. * pow(f*mu/k,2.) * (2.*I2211[idx] + 6.*I1311[idx]) *D4\
                        + 8. * pow(f,2.)*(mu/k) * (b1*(J12111[idx] + J11211[idx] + J21111[idx])) *D4\
                        + 2. * pow(f,2.) * (pow(b1,2.)*N11[idx]) *D4\
                        + 2. * pow(f,2.) * (4.*b1*J12102[idx] + 2.*b1*J21102[idx] + b2*Jdelta202[idx] + 2.*bG2*JG202[idx] + pow(b1,2.)*Plin[idx]*sigma_v2) *D4\
                        - 2. * pow(f,2.) * (c20 + c22 * pow(mu,2.))*Plin[idx] *D2;

        //3-rd moment
        J21112x[idx] = pfo -> pk_halo_rsd_nl[idx] -> J21112x[index_k];
        J21112y[idx] = pfo -> pk_halo_rsd_nl[idx] -> J21112y[index_k];
        J21112[idx]  = (J21112x[idx] + J21112y[idx] * pow(mu,2.));
        N12x[idx]    = pfo -> pk_halo_rsd_nl[idx] -> N12x[index_k];
        N12y[idx]    = pfo -> pk_halo_rsd_nl[idx] -> N12y[index_k];
        N12[idx]     = (N12x[idx] * mu + N12y[idx] * pow(mu,3.));
        J12112x[idx] = pfo -> pk_halo_rsd_nl[idx] -> J12112x[index_k];
        J12112y[idx] = pfo -> pk_halo_rsd_nl[idx] -> J12112y[index_k];
        J12112[idx]  = (J12112x[idx] + J12112y[idx] * pow(mu,2.));

        Moment_3[idx] = - 6. * pow(f,3.)*(mu/k) * b1*Plin[idx]*sigma_v2 *D4\
                        + 12.* pow(f,3.)*(mu/k) * (J21112[idx] + 2.*J12112[idx]) *D4\
                        - 6. * pow(f,3.)*(mu/k) * b1*Plin[idx]*sigma_v2 *D4\
                        + 12.* pow(f,3.) * b1*N12[idx] *D4\
                        + 6. * pow(f,3.)*(mu/k) * (c30 + c32*pow(mu,2.))*Plin[idx] *D2;

        //4-th moment
        N22x[idx] = pfo -> pk_halo_rsd_nl[idx] -> N22x[index_k];
        N22y[idx] = pfo -> pk_halo_rsd_nl[idx] -> N22y[index_k];
        N22z[idx] = pfo -> pk_halo_rsd_nl[idx] -> N22z[index_k];
        N22[idx]  = (N22x[idx] + N22y[idx] * pow(mu,2.) + N22z[idx] * pow(mu,4.));

        Moment_4[idx] = -24. * pow(f,4.)*pow(mu/k,2.) * Plin[idx]*sigma_v2 *D4\
                        +12. * pow(f,4.) * N22[idx] *D4\
                        +24. * pow(f,4.)*pow(mu/k,2.) * c42*Plin[idx] *D2;

        //RSD expansion
        RSD[idx] = Moment_0[idx] + k*mu * Moment_1[idx] + (1./2.) * pow(k*mu,2.) * Moment_2[idx] \
                 + (1./6.) * pow(k*mu,3.) * Moment_3[idx] + (1./24.) * pow(k*mu,4.) * Moment_4[idx];
    }

    double P_wiggle = Plin[lin] - Plin[no_wiggle];
    double RSD_IR_Ressummed = pow(b1 + f*pow(mu,2.), 2.) * (Plin[no_wiggle] + suppression*P_wiggle*(1. + pow(k,2.)*sigma_tot))\
                            + RSD[no_wiggle] + suppression * (RSD[lin] - RSD[no_wiggle]);
    
    *result = RSD_IR_Ressummed;

    // FILE *fpa9;
    // char file_name9[50];
    // sprintf(file_name9, "data/RSD_lin_%g.txt",mu);
    // fpa9 = fopen(file_name9, "a");
    // fprintf(fpa9, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n", k, Moment_0[lin], Moment_1[lin], Moment_2[lin], Moment_3[lin], Moment_4[lin], RSD[lin]);
    // fclose(fpa9);

    // FILE *fpa10;
    // char file_name10[50];
    // sprintf(file_name10, "data/RSD_nw_%g.txt",mu);
    // fpa10 = fopen(file_name10, "a");
    // fprintf(fpa10, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n", k, Moment_0[no_wiggle], Moment_1[no_wiggle], Moment_2[no_wiggle], Moment_3[no_wiggle], Moment_4[no_wiggle], RSD[no_wiggle]);
    // fclose(fpa10);

    return _SUCCESS_;
}

int RSD_Multipole(struct fourier *pfo, struct background *pba, int index_k, double z, int l, double * result)
{
    //extern struct globals gb;
    double integ=0., error=0.;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);

    struct integrand_parameters2 par; 

    double mu_min = -1.;
    double mu_max =  1.;

    gsl_function F;
    F.function = &RSD_Multipole_integrand;
    F.params = &par;

    par.pfo = pfo;
    par.pba = pba;
    par.p5  = z;
    par.p19 = l;
    par.p23 = index_k;
    gsl_integration_qags(&F,mu_min,mu_max,0.0,1.0e-5,1000000,w,&integ,&error);
    gsl_integration_workspace_free(w);

    *result = integ;

    double k = pfo->k[index_k];

    FILE *fpa10;
    char file_name10[50];
    sprintf(file_name10, "data/%d_Multipol_%g.txt",l,z);
    fpa10 = fopen(file_name10, "a");
    fprintf(fpa10, "%12.6e %12.6e\n",k,integ);
    fclose(fpa10);

    return _SUCCESS_;
}

double RSD_Multipole_integrand(double x, void *par)
{
    double result = 0;
    
    struct integrand_parameters2 pij;
    pij = *((struct integrand_parameters2 *)par);

    struct fourier *pfo    = pij.pfo;
    struct background *pba = pij.pba;
    double z               = pij.p5;
    int    l               = pij.p19;
    int    index_k         = pij.p23;

    double rsd;
    RSD_IR_Ressummed(pfo, pba, index_k, z, x, &rsd);

    result = rsd * Legendre_Polynomial(l,x)* (2.*l + 1.)/2.;

    return result;
}

double Legendre_Polynomial(int l, double mu)
{
    double result = 0.;
    if(l==0){
        result = 1.;
    }
    else if(l==1){
        result = mu;
    }
    else if(l==2){
        result = 0.5*(3.*pow(mu,2.)-1.);
    }
    else if(l==3){
        result = 0.5*(5.*pow(mu,3.)-3.*mu);
    }
    else if(l==4){
        result = (35.*pow(mu,4.)-30.*pow(mu,2.)+3.)/8.;
    }
    else{
        result = 0.;
    }
    return result;
}

// int FFTLog_fill_bias_vector(struct fourier *pfo, double b1, double b2, double bG2, double btd, double cs2, double R2){
//     pfo -> fft_ws -> b1  = b1;
//     pfo -> fft_ws -> b2  = b2;
//     pfo -> fft_ws -> bG2 = bG2;
//     pfo -> fft_ws -> btd = btd;
//     pfo -> fft_ws -> cs2 = cs2;
//     pfo -> fft_ws -> R2  = R2;

//     // 0-th moment
//     pfo -> fft_ws -> bias -> np_bias_vec0 = make_1Darray(6);
//     pfo -> fft_ws -> bias -> p_bias_vec0  = make_1Darray(2);
//     pfo -> fft_ws -> bias -> np_bias_vec0[0] = 2.* pow(b1, 2.);
//     pfo -> fft_ws -> bias -> np_bias_vec0[1] = 2.* b1 * b2;
//     pfo -> fft_ws -> bias -> np_bias_vec0[2] = 4.* b1 * bG2;
//     pfo -> fft_ws -> bias -> np_bias_vec0[3] = 0.5 * pow(b2, 2.);
//     pfo -> fft_ws -> bias -> np_bias_vec0[4] = 2.* pow(bG2, 2.);
//     pfo -> fft_ws -> bias -> np_bias_vec0[5] = 2.* b2 * bG2;
//     pfo -> fft_ws -> bias -> p_bias_vec0[0] = 6.* pow(b1, 2.);
//     pfo -> fft_ws -> bias -> p_bias_vec0[1] = 8. * b1 * (bG2 + 2./5. * btd);

//     // 1-st moment
//     pfo -> fft_ws -> bias -> np_bias_vec1 = make_1Darray(7);
//     pfo -> fft_ws -> bias -> p_bias_vec1  = make_1Darray(2);
//     pfo -> fft_ws -> bias -> np_bias_vec1[0] = 2. * b1;
//     pfo -> fft_ws -> bias -> np_bias_vec1[1] = b2;
//     pfo -> fft_ws -> bias -> np_bias_vec1[2] = 2. * bG2;
//     pfo -> fft_ws -> bias -> np_bias_vec1[3] = 4. * (bG2 + 2./5. * btd);
//     pfo -> fft_ws -> bias -> np_bias_vec1[4] = pow(b1, 2.);
//     pfo -> fft_ws -> bias -> np_bias_vec1[5] = 0.5 * b1 * b2;
//     pfo -> fft_ws -> bias -> np_bias_vec1[6] = b1 * bG2;
//     pfo -> fft_ws -> bias -> p_bias_vec1[0] = 3. * b1;
//     pfo -> fft_ws -> bias -> p_bias_vec1[1] = pow(b1, 2.);

//     // 2-nd moment
//     pfo -> fft_ws -> bias -> np_bias_vec2 = make_1Darray(6);
//     pfo -> fft_ws -> bias -> p_bias_vec2  = make_1Darray(3);
//     pfo -> fft_ws -> bias -> np_bias_vec2[0] = 2. * b1;
//     pfo -> fft_ws -> bias -> np_bias_vec2[1] = b2;
//     pfo -> fft_ws -> bias -> np_bias_vec2[2] = 2. * bG2;
//     pfo -> fft_ws -> bias -> np_bias_vec2[3] = 2.;
//     pfo -> fft_ws -> bias -> np_bias_vec2[4] = b1;
//     pfo -> fft_ws -> bias -> np_bias_vec2[5] = pow(b1, 2.);
//     pfo -> fft_ws -> bias -> p_bias_vec2[0] = 4. * b1;
//     pfo -> fft_ws -> bias -> p_bias_vec2[1] = 6.;
//     pfo -> fft_ws -> bias -> p_bias_vec2[2] = b1;

//     // 3-rd moment
//     pfo -> fft_ws -> bias -> np_bias_vec3 = make_1Darray(2);
//     pfo -> fft_ws -> bias -> p_bias_vec3  = make_1Darray(1);
//     pfo -> fft_ws -> bias -> np_bias_vec3[0] = 1.;
//     pfo -> fft_ws -> bias -> np_bias_vec3[1] = b1;
//     pfo -> fft_ws -> bias -> p_bias_vec3[0] = 2.;

//     // 4-th moment
//     pfo -> fft_ws -> bias -> np_bias_vec4 = make_1Darray(1);
//     pfo -> fft_ws -> bias -> np_bias_vec4[0] = 1.;

//     return _SUCCESS_;
// }