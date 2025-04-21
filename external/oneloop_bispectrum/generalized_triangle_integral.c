/** @file oneloop.c
 *
 * author: Fynn Janssen, 2025
 * strongly based on https://arxiv.org/abs/2212.07421, https://github.com/dbraganca/python-integer-powers
*/

#include "generalized_triangle_integral.h"





// ************************************************************************************************
// Master Integrals
// ************************************************************************************************


/** Solves triangle master integral: integral over 1/(q2+M1)/(k1mq2 + M2)/(k2pq2 + M3)
 * @param k12           Input: wavenumber squared k1^2
 * @param k22           Input: wavenumber squared k2^2
 * @param k32           Input: wavenumber squared k3^2
 * @param M1            Input: complex mass of propagator
 * @param M2            Input: complex mass of propagator
 * @param M3            Input: complex mass of propagator
 * @param T_out         Output: pointer to value for the triangle master integral
 * @return the error status
 */

 int T_master(
              struct gen_tri_integral *pti,
              double k12,
              double k22,
              double k32,
              class_complex M1,
              class_complex M2,
              class_complex M3,
              class_complex *T_out
              ){
    
    // special case in which all masses are zero
    class_complex M0 = class_complex(0., 0.);
    if (M1==M0 && M2==M0 && M3==M0){
        *T_out = _PI_*_PI_*_PI_/sqrt(k12)/sqrt(k22)/sqrt(k32);
    } else {
        class_complex T_master_contr0, T_master_contr1;
        class_call(util_Tmaster_contr(pti,
                                      0, 
                                      k12,
                                      k22,
                                      k32,
                                      M1,
                                      M2,
                                      M3,
                                      &T_master_contr0),
                    pti->error_message,
                    pti->error_message);

        class_call(util_Tmaster_contr(pti,
                                      1,
                                      k12,
                                      k22,
                                      k32,
                                      M1,
                                      M2,
                                      M3,
                                      &T_master_contr1),
                    pti->error_message, pti->error_message);

        *T_out = T_master_contr1 - T_master_contr0;
    }

    return _SUCCESS_;
}



/** Solves bubble master integral: integral over 1/(q2+M1)/(kmq2 + M2), note: is the same as integral over 1/(q2+M2)/(kpq2 + M1)
 * @param k2            Input: wavenumber squared k^2
 * @param M1            Input: complex mass of propagator
 * @param M2            Input: complex mass of propagator
 * @param B_out         Output: pointer to value for the bubble master integral
 * @return the error status
 */

int B_master(
             struct gen_tri_integral *pti,
             double k2,
             class_complex M1,
             class_complex M2,
             class_complex *B_out
             ){
    
    class_complex m1 = M1/k2;
    class_complex m2 = M2/k2;
    class_complex I = class_complex(0., 1.);
    int sign = 0;

    class_complex arglog0 = I*(m1 - m2 - 1.) + 2.*sqrt(m1);
    class_complex arglog1 = I*(m1 - m2 + 1.) + 2.*sqrt(m2);

    if (cimag(arglog0) > 0 && cimag(arglog1) < 0){
        sign = 1;
    } else {
        sign = 0;
    }

    *B_out = _PI_*_PI_/sqrt(k2)*(I*(log(arglog0)-log(arglog1))+2.*_PI_*sign);
    return _SUCCESS_;
}


/** Solves tadpole master integral: integral over q^2n/(q2+M)^d
 * @param n             Input: power of numerator q^2
 * @param d             Input: power of denominator q^2+M
 * @param M             Input: complex mass of propagator
 * @param Tad_out       Output: pointer to value for the triangle master integral
 * @return the error status
 */

int Tad_master(
               struct gen_tri_integral *pti,
               int n,
               int d,
               class_complex M,
               class_complex *Tad_out
               ){
    
    class_complex temp;


    // class_test(n<-1 || d<1,
    //     pti->error_message,
    //     "Divergent integral: Cannot call Tad_master for n<-1 or d<1. You have n=%d, d=%d", n, d);

    // class_test(d<n+1,
    //     pti->error_message,
    //    "Divergent integral: Cannot call Tad_master for d<n+1. You have n=%d, d=%d", n, d);

    if (M==class_complex(0., 0.)){
        *Tad_out = class_complex(0., 0.); // TODO: dont understand why zero -> dimensional argument?
    } else {
        *Tad_out = _TWOPI_*sqrt(M)* pow(M, n+1-d) *tgamma(n+1.5)*tgamma(d-n-1.5)/tgamma(d);               
    }
    return _SUCCESS_;
}




// ************************************************************************************************
// Recursion Steps
// ************************************************************************************************

/** Solves general trianlge integral: L(n1, d1, n2, d2, n3, d3) = int k1mq^2n1 q^2n2 k2pq^2n3 / (k1mq^2+M1)^d1 / (q^2+M2)^d2 / (k2pq^2+M3)^d3
 * @param n1                Input: power of numerator (k1_vec-q_vec)^2
 * @param d1                Input: power of denominator (k1_vec-q_vec)^2 + M1
 * @param n2                Input: power of numerator q^2
 * @param d2                Input: power of denominator q^2 + M2
 * @param n3                Input: power of numerator (k2_vec+q_vec)^2
 * @param d3                Input: power of denominator (k2_vec+q_vec)^2 + M3
 * @param k12               Input: wavenumber squared k1^2
 * @param k22               Input: wavenumber squared k2^2
 * @param k32               Input: wavenumber squared k3^2
 * @param M1                Input: complex mass of propagator
 * @param M2                Input: complex mass of propagator
 * @param M3                Input: complex mass of propagator
 * @param L_out             Output: pointer to value for the generalized triangle integral
 * @return the error status
 */

int L_recursion(
                struct gen_tri_integral *pti,
                int n1,
                int d1,
                int n2,
                int d2,
                int n3,
                int d3,
                double k12,
                double k22,
                double k32,
                class_complex M1,
                class_complex M2,
                class_complex M3,
                class_complex *L_out
                ){

    class_complex M0, Ms[3];
    doublets_L idxs[3];
    int doublet_idx;

    idxs[0].n = n1; idxs[1].n = n2; idxs[2].n = n3;
    idxs[0].d = d1; idxs[1].d = d2; idxs[2].d = d3; 
    Ms[0] = M1;     Ms[1] = M2;     Ms[2] = M3;

    
    // recursion terminates when n1=n2=n3=0 -> move on with T recursion
    if (n1==0 and n2==0 and n3==0){
        class_call(T_recursion(pti,
                               d1, 
                               d2,
                               d3,
                               k12,
                               k22,
                               k32,
                               M1,
                               M2,
                               M3,
                               L_out),
                    pti->error_message,
                    pti->error_message);
        // printf("\nd1=%d, d2=%d, d3=%d", d1, d2, d3);
        
        return _SUCCESS_;
    }

    // This function can be split in 3 doublet (n_i, d_i) which are saved in idxs
    // The procedure is the same for each doublet -> util_L_step, and has to be done until each n_i=0
    if (n1!=0){
        doublet_idx = 0;
    } else if (n2!=0){
        doublet_idx = 1;
    } else if (n3!=0){
        doublet_idx = 2;
    }

    class_call(util_L_step(pti,
                           idxs,
                           doublet_idx,
                           k12,
                           k22,
                           k32,
                           Ms,
                           L_out),
                pti->error_message,
                pti->error_message);

    return _SUCCESS_;
}


/** Solves T integral: T(d1, d2, d3) = L(0, d1, 0, d2, 0, d3) = int 1 / (k1mq^2+M1)^d1 / (q^2+M2)^d2 / (k2pq^2+M3)^d3
 * @param d1                Input: power of denominator (k1_vec-q_vec)^2 + M1
 * @param d2                Input: power of denominator q^2+M2
 * @param d3                Input: power of denominator (k2_vec+q_vec)^2 + M3
 * @param k12               Input: wavenumber squared k1^2
 * @param k22               Input: wavenumber squared k2^2
 * @param k32               Input: wavenumber squared k3^2
 * @param M1                Input: complex mass of propagator
 * @param M2                Input: complex mass of propagator
 * @param M3                Input: complex mass of propagator
 * @param L_out             Output: pointer to value for the T integral
 * @return the error status
 */

int T_recursion(
                struct gen_tri_integral *pti,
                int d1,
                int d2,
                int d3,
                double k12,
                double k22,
                double k32,
                class_complex M1,
                class_complex M2,
                class_complex M3,
                class_complex *T_out
                ){
    double cos12, cos23, cos31;

    cos12 = (k32-k12-k22)/(2.*sqrt(k12*k22));
    cos23 = (k12-k22-k32)/(2.*sqrt(k22*k32));
    cos31 = (k22-k32-k12)/(2.*sqrt(k32*k12));

    if (d1+d2+d3<1){
        printf("Error: T_recursion is divergent. d1+d2+d3=%d", d1+d2+d3);
    }
    if (d1<-4){
        printf("Warning: not allowed exponent in T_recursion: d1=%d", d1);
    }
    if (d2<-4){
        printf("Warning: not allowed exponent in T_recursion: d2=%d", d2);
    }
    if (d2<-4){
        printf("Warning: not allowed exponent in T_recursion: d2=%d", d2);
    }

    // Consider call conditions that could terminate the recursion:
    
    // d1=d2=d3=1:
    if (d1==1 && d2==1 && d3==1){
        class_call(T_master(pti,
                            k12,
                            k22,
                            k32,
                            M1,
                            M2,
                            M3,
                            T_out),
                    pti->error_message,
                    pti->error_message);
        return _SUCCESS_;
    }

    // cases where one of the exponents is zero: redefine the integration variable and identify the B_recursion integral
    if (d1==0){
        class_call(B_recursion(pti,
                               d2,
                               d3,
                               k22,
                               M2,
                               M3,
                               T_out),
                    pti->error_message,
                    pti->error_message);
        return _SUCCESS_;
    }
    if (d2==0){
        class_call(B_recursion(pti,
                               d3,
                               d1,
                               k32,
                               M3,
                               M1,
                               T_out),
                    pti->error_message,
                    pti->error_message);
        return _SUCCESS_;
    }
    if (d3==0){
        class_call(B_recursion(pti,
                               d1,
                               d2,
                               k12,
                               M1,
                               M2,
                               T_out),
                    pti->error_message,
                    pti->error_message);
        return _SUCCESS_;
    }

    // cases with negative exponents -> tensor reduction
    // two cases: one exponent is less than 0 -> tensor_red_one, two exponents are less than 0 -> tensor_red_two
    if (d1<0){
        if (d2>0 && d3>0){
            class_call(tensor_red_one(pti,
                                      -d1,
                                      d2,
                                      d3,
                                      k12,
                                      k22,
                                      cos12,
                                      M2,
                                      M3,
                                      T_out),
                        pti->error_message,
                        pti->error_message);
            return _SUCCESS_;

        } else if (d2>0 && d3<0){
            class_call(tensor_red_two(pti,
                                      -d1,
                                      -d3,
                                      d2,
                                      k12,
                                      k22,
                                      cos12,
                                      M2,
                                      T_out),
                        pti->error_message,
                        pti->error_message);
            return _SUCCESS_;

        } else if (d2<0 && d3>0){
            class_call(tensor_red_two(pti,
                                      -d2,
                                      -d1,
                                      d3,
                                      k22,
                                      k32,
                                      cos23,
                                      M3,
                                      T_out),
                        pti->error_message,
                        pti->error_message);
            return _SUCCESS_;

        } else {
            // last case is d2<0, d3<0 which automatically violates d1+d2+d3>1
            printf("Undefined case in T_recursion: d1, d2, d3 = %d, %d, %d < 0", d1, d2, d3);
            return _FAILURE_;
        }
    } else if (d2<0){
        if (d1>0 && d3>0){
            class_call(tensor_red_one(pti,
                                      -d2,
                                      d1,
                                      d3,
                                      k12,
                                      k32,
                                      cos31,
                                      M1,
                                      M3,
                                      T_out),
                        pti->error_message,
                        pti->error_message);
            return _SUCCESS_;

        } else if (d1>0 && d3<0){
            class_call(tensor_red_two(pti,
                                      -d3,
                                      -d2,
                                      d1,
                                      k32,
                                      k12,
                                      cos31,
                                      M1,
                                      T_out),
                        pti->error_message,
                        pti->error_message);
            return _SUCCESS_;

        }
    } else if (d3<0){
        // only case that is left: d1, d2 > 0
        class_call(tensor_red_one(pti,
                                  -d3,
                                  d1,
                                  d2,
                                  k32,
                                  k12,
                                  cos12,
                                  M1,
                                  M2,
                                  T_out),
                    pti->error_message,
                    pti->error_message);
        return _SUCCESS_;

    }

    // define utility variables
    class_complex cpm0, cmp0, cm0p, cp0m, c0pm, c0mp, c000, k1s, k2s, k3s, jac, ks11, ks12, ks22, ks23, ks31, ks33, k1s2, k2s2, k3s2, nu1_c, nu2_c, nu3_c, Ndim, T_recursion_temp[7];
    int i, idx1, idx2, idx3, nu1, nu2, nu3;
    
    k1s = k12 + M1 + M2;
    k2s = k22 + M2 + M3;
    k3s = k32 + M3 + M1;

    k1s2 = k1s*k1s;
    k2s2 = k2s*k2s;
    k3s2 = k3s*k3s;

    jac = -8.*M1*M2*M3 + 2.*k1s2*M3 + 2.*k2s2*M1 + 2.*k3s2*M2 - 2.*k1s*k2s*k3s;
    
    ks11 = (-4.*M1*M2 + k1s2)/jac;
    ks12 = (-2.*k3s*M2 + k1s*k2s)/jac;
    ks22 = (-4.*M2*M3 + k2s2)/jac;
    ks23 = (-2.*k1s*M3 + k2s*k3s)/jac;
    ks31 = (-2.*k2s*M1+k1s*k3s)/jac;
    ks33 = (-4.*M1*M3+k3s2)/jac;


    // TODO: maybe this is not the most efficient way to handle this
    // especially this complex devision
    if (d1>1){
        nu1 = d1 - 1;
        nu2 = d2;
        nu3 = d3;

        nu1_c = class_complex(nu1, 0.);
        nu2_c = class_complex(nu2, 0.);
        nu3_c = class_complex(nu3, 0.);

        Ndim = class_complex(3. - nu1 - nu2 - nu3, 0.);

        cpm0 = -ks23;
        cmp0 = (ks22*nu2_c)/nu1_c;
        cm0p = (ks22*nu3_c)/nu1_c;
        cp0m = -ks12;
        c0pm = -(ks12*nu2_c)/nu1_c;
        c0mp = -(ks23*nu3_c)/nu1_c;
        c000 = (-nu3_c+Ndim)*ks12/nu1_c - (-nu1_c+Ndim)*ks22/nu1_c + (-nu2_c+Ndim)*ks23/nu1_c;

    } else if (d2>1){
        nu1 = d1;
        nu2 = d2 - 1;
        nu3 = d3;

        nu1_c = class_complex(nu1, 0.);
        nu2_c = class_complex(nu2, 0.);
        nu3_c = class_complex(nu3, 0.);

        Ndim = class_complex(3. - nu1 - nu2 - nu3, 0.);

        cpm0 = (ks33*nu1_c)/nu2_c;
        cmp0 = -ks23;
        cm0p = -(ks23*nu3_c)/nu2_c;
        cp0m = -(ks31*nu1_c)/nu2_c;
        c0pm = -ks31;
        c0mp = (ks33*nu3_c)/nu2_c;
        c000 = (-nu1_c + Ndim)*ks23/nu2_c + (-nu3_c + Ndim)*ks31/nu2_c - (-nu2_c + Ndim)*ks33/nu2_c;

    } else if (d3>1){
        nu1 = d1;
        nu2 = d2;
        nu3 = d3 - 1;

        nu1_c = class_complex(nu1, 0.);
        nu2_c = class_complex(nu2, 0.);
        nu3_c = class_complex(nu3, 0.);

        Ndim = class_complex(3. - nu1 - nu2 - nu3, 0.);

        cpm0 = -(ks31*nu1_c)/nu3_c;
        cmp0 = -(ks12*nu2_c)/nu3_c;
        cm0p = -ks12;
        cp0m = (ks11*nu1_c)/nu3_c;
        c0pm = (ks11*nu2_c)/nu3_c;
        c0mp = -ks31;
        c000 = -(-nu3_c + Ndim)*ks11/nu3_c + (-nu1_c + Ndim)*ks12/nu3_c + (-nu2_c + Ndim)*ks31/nu3_c;
    }

    // iterate through all relevant class_calls to avoid clutter
    // TODO: maybe create global cache
    for (i=0; i<7; i++){
        if (i==0){
            idx1 = nu1;
            idx2 = nu2;
            idx3 = nu3;
        } else if(i==1){
            idx1 = nu1;
            idx2 = nu2-1;
            idx3 = nu3+1;
        } else if(i==2){
            idx1 = nu1;
            idx2 = nu2+1;
            idx3 = nu3-1;
        } else if(i==3){
            idx1 = nu1-1;
            idx2 = nu2;
            idx3 = nu3+1;
        } else if(i==4){
            idx1 = nu1+1;
            idx2 = nu2;
            idx3 = nu3-1;
        } else if(i==5){
            idx1 = nu1-1;
            idx2 = nu2+1;
            idx3 = nu3;
        } else if(i==6){
            idx1 = nu1+1;
            idx2 = nu2-1;
            idx3 = nu3;
        }

        class_call(T_recursion(pti,
                               idx1, 
                               idx2,
                               idx3,
                               k12,
                               k22,
                               k32,
                               M1,
                               M2,
                               M3,
                               &T_recursion_temp[i]),
                    pti->error_message,
                    pti->error_message);
    }



    *T_out = (c000*T_recursion_temp[0]
            + c0mp*T_recursion_temp[1]
            + c0pm*T_recursion_temp[2]
            + cm0p*T_recursion_temp[3]
            + cp0m*T_recursion_temp[4]
            + cmp0*T_recursion_temp[5]
            + cpm0*T_recursion_temp[6]);


    return _SUCCESS_;
}


/** Solves general bubble integral: B(d1, d2) = int 1 / (kmq^2+M1)^d1 / (q^2+M2)^d2
 * @param d1                Input: power of denominator (k_vec-q_vec)^2 + M1
 * @param d2                Input: power of denominator q^2 + M2
 * @param k2                Input: wavenumber squared k^2
 * @param M1                Input: complex mass of propagator
 * @param M2                Input: complex mass of propagator
 * @param B_out             Output: pointer to value for the general bubble integral
 * @return the error status
 */

int B_recursion(
                struct gen_tri_integral *pti,
                int d1,
                int d2,
                double k2,
                class_complex M1,
                class_complex M2,
                class_complex *B_out
                ){
    
    // TODO: any special cases if the masses are zero?
    int nu1, nu2;
    class_complex B_temp1, B_temp2, B_temp3, c1, c2, c3, jac, k1s, Ndim, nu1_c, nu2_c;
    // terminate the recursion for the following cases:
    if (d1==1 && d2==1){
        // Bubble master integral
        class_call(B_master(pti,
                            k2,
                            M1,
                            M2,
                            B_out),
                    pti->error_message,
                    pti->error_message);
    
    } else if (d1==0){
        if (d2==1){
            printf("Error in B_recursion: d1=0, d2=1, divergent integral\n");
            *B_out = class_complex(0., 0.);
        } else {
            // Tad pole master integral with numerator n=1
            class_call(Tad_master(pti,
                                  0,
                                  d2,
                                  M2,
                                  B_out),
                        pti->error_message,
                        pti->error_message);
        }

    } else if (d2==0){
        if (d1==1){
            printf("Error in B_recursion: d1=1, d2=0, divergent integral\n");
            *B_out = class_complex(0., 0.);
        } else{
            // Tad pole master integral with numerator n=1
            class_call(Tad_master(pti,
                                  0,
                                  d1,
                                  M1,
                                  B_out),
                        pti->error_message,
                        pti->error_message);
        }

    } else if (d1<0 && d2<0){
        // TODO: this integral is divergent, why set it to zero? maybe because it cannot appear? If so -> class_test
        *B_out = class_complex(0., 0.);
    } else if (d1<0 && d2>0){
        if (d1+d2<1){ 
            printf("Error in B_recursion: d1=%d, d2=%d, divergent integral\n", d1, d2);
            *B_out = class_complex(0., 0.);
        } else{
            // massive numerator integral
            class_call(massive_num(pti, 
                                   -d1,
                                   d2,
                                   k2,
                                   M1,
                                   M2,
                                   B_out),
                        pti->error_message,
                        pti->error_message);
        }
    } else if (d1>0 && d2<0){
        if (d1+d2<1){
            printf("Error in B_recursion: d1=%d, d2=%d, divergent integral\n", d1, d2);
            *B_out = class_complex(0., 0.);
        } else{
            // TODO: what it 2+|d2|>d1 (2 coming from the jacobian) then the integral is divergent. return 0?
            // massive numerator integral
            class_call(massive_num(pti, 
                                   -d2,
                                   d1,
                                   k2,
                                   M2,
                                   M1,
                                   B_out),
                        pti->error_message,
                        pti->error_message);
        }

    } else if (d1>=1 && d2>=1){
        // first apply the recursion to reduce d1 until d1=1, the do the same for d2
        // The order is arbitrary
        Ndim = class_complex(4. - d1 - d2, 0.);
        k1s = k2 + M1 + M2;
        jac = k1s*k1s - 4.*M1*M2;
        if (d1>1){
            nu1 = d1-1;
            nu2 = d2;
            nu1_c = class_complex(nu1, 0.);
            nu2_c = class_complex(nu2, 0.);
            c1 = ((2.*M2-k1s)/nu1_c*Ndim - 2.*M2 + (k1s*nu2_c)/nu1_c)/jac;
            c2 = k1s/jac;
            c3 = (-2.*M2/nu1_c*nu2_c)/jac;
        } else if (d2>1){
            nu1 = d1;
            nu2 = d2-1;
            nu1_c = class_complex(nu1, 0.);
            nu2_c = class_complex(nu2, 0.);
            c1 = ((2.*M1 - k1s)/nu2_c*Ndim + (k1s*nu1_c)/nu2_c - 2.*M1)/jac;
            c2 = (-2.*M1/nu2_c*nu1_c)/jac;
            c3 = k1s/jac;
        }
        class_call(B_recursion(pti,
                               nu1,
                               nu2,
                               k2,
                               M1,
                               M2,
                               &B_temp1),
                    pti->error_message,
                    pti->error_message);
        
        class_call(B_recursion(pti,
                               nu1+1,
                               nu2-1,
                               k2,
                               M1,
                               M2,
                               &B_temp2),
                    pti->error_message,
                    pti->error_message);

        class_call(B_recursion(pti,
                               nu1-1,
                               nu2+1,
                               k2,
                               M1,
                               M2,
                               &B_temp3),
                    pti->error_message,
                    pti->error_message);

        *B_out = c1*B_temp1 + c2*B_temp2 + c3*B_temp3;
    } else {
        printf("Error in B_recursion: case not considered");
    }
        
    return _SUCCESS_;
}


/** Intermediate integral of tensor reduction case 1: solves integrals of the form (2*k_vec.q_vec)^m * q^2n / (q^2+M1)^d1 / (kpq^2+M2)^d2
 * @param m                 Input: power of numerator scalar product k2_vec.q_vec
 * @param n                 Input: power of numerator a^2
 * @param d1                Input: power of denominator q^2+M1
 * @param d2                Input: power of denominator (k2_vec+q_vec)^2 + M2
 * @param k2                Input: wavenumber squared k^2
 * @param M1                Input: complex mass of propagator
 * @param M2                Input: complex mass of propagator
 * @param I_out             Output: pointer to value for the scalar product integral (case 1)
 * @return the error status
 */

 int scalar_prod_one(
                     struct gen_tri_integral *pti,
                     int m,
                     int n,
                     int d1,
                     int d2,
                     double k2,
                     class_complex M1,
                     class_complex M2,
                     class_complex *I_out
                    ) {
    
    class_complex M1_pow, result=class_complex(0., 0.), B_cache[m+n+1][m+1];
    double binomial_k, binomial_i, binomial_j;
    int k, i, j, d1_eff, d2_eff, d1_idx, d2_idx;
    class_complex B_temp;

    // Initially set all cache values to zero
    for (i=0; i<=m+n; i++){
        for (j=0; j<=m; j++){
            B_cache[i][j] = class_complex(0., 0.);
        }
    }


    for (k=0; k<=n; k++){
        // binomial coefficient for the q^2n term
        class_call(util_binomial(pti,
                                 n,
                                 k, 
                                 &binomial_k),
                    pti->error_message,
                    pti->error_message);

        M1_pow = pow(M1, n-k);
        for (i=0; i<=m; i++){
            // first binomial coefficient for the k.q^m term
            class_call(util_binomial(pti,
                                    m,
                                    i,
                                    &binomial_i),
                        pti->error_message,
                        pti->error_message);

            for (j=0; j<=m-i; j++){
                // second binomial coefficient for the k.q^m term
                class_call(util_binomial(pti,
                                         m-i,
                                         j,
                                         &binomial_j),
                            pti->error_message,
                            pti->error_message);
                            
                d1_idx = i+k;
                d2_idx = j;
                B_temp = B_cache[d1_idx][d2_idx];
                // Check if this contribution of B_recurions has been computed already
                if (B_temp == class_complex(0., 0.)){
                    class_call(B_recursion(pti,
                                           d2-d2_idx,
                                           d1-d1_idx,
                                           k2,
                                           M2,
                                           M1,
                                           &B_temp),
                               pti->error_message,
                               pti->error_message);

                    // Save this value in the cache
                    B_cache[d1_idx][d2_idx] = B_temp;
                }


                result += binomial_k*binomial_i*binomial_j * pow(-1., n+i-k)*M1_pow * pow(M1-M2-k2, m-i-j)*B_temp;
            }
        }
    }
    *I_out = result;

    return _SUCCESS_;
}


/** Tensor reduction case 1: solves integrals of the form k1mq^2n1 / (q^2+M1)^d1 / (k2pq^2+M2)^d2
 * @param n                 Input: power of numerator k1_vec-q_vec
 * @param d1                Input: power of denominator q^2+M1
 * @param d2                Input: power of denominator (k2_vec+q_vec)^2 + M2
 * @param k12               Input: wavenumber squared k1^2
 * @param k22               Input: wavenumber squared k2^2
 * @param cos12             Input: angle between k1_vec and k2_vec
 * @param M1                Input: complex mass of propagator
 * @param M2                Input: complex mass of propagator
 * @param I_out             Output: pointer to value for the tensor reduction integral (case 1)
 * @return the error status
 */

int tensor_red_one(
                   struct gen_tri_integral *pti,
                   int n,
                   int d1,
                   int d2,
                   double k12,
                   double k22,
                   double cos12,
                   class_complex M1,
                   class_complex M2,
                   class_complex *I_out
                   ){
    
    class_complex I_temp1, I_temp2, I_temp3, A1, A2, A3, result, result_int;
    double k14, k24, cos12_sq, sin12_sq, binomial_i, binomial_j;
    int i, j;

    cos12_sq = cos12*cos12;

    // Iterate over the main loop
    // Since different number of scalar product lead to different forms of the solution
    // there are many if statements in this loop
    for (i=0; i<=n; i++){
        class_call(util_binomial(pti,
                                 n,
                                 i,
                                 &binomial_i),
                   pti->error_message,
                   pti->error_message);

        for (j=0; j<=n-i; j++){
            if (d1+d2<1+j+i) {
                printf("Problem in tensor_red_one, d1=%d, d2=%d, j=%d, i=%d\n", d1, d2, j, i);
            }
            class_call(util_binomial(pti,
                                     n-i,
                                     j, 
                                     &binomial_j),
                       pti->error_message,
                       pti->error_message);
            
            if (j==0){
                class_call(scalar_prod_one(pti,
                                           j,
                                           i,
                                           d1,
                                           d2,
                                           k22,
                                           M1,
                                           M2,
                                           &I_temp1),
                           pti->error_message,
                           pti->error_message);
                
                result_int = I_temp1;

            } else if (j==1){
                class_call(scalar_prod_one(pti,
                                           j,
                                           i,
                                           d1,
                                           d2,
                                           k22,
                                           M1,
                                           M2,
                                           &I_temp1),
                           pti->error_message,
                           pti->error_message);

                result_int = sqrt(k12/k22)*cos12 * I_temp1;

            } else if (j==2) {
                class_call(scalar_prod_one(pti,
                                           j,
                                           i,
                                           d1,
                                           d2,
                                           k22,
                                           M1,
                                           M2,
                                           &I_temp1),
                           pti->error_message,
                           pti->error_message);
            
                class_call(scalar_prod_one(pti,
                                           j-2,
                                           i+1,
                                           d1,
                                           d2,
                                           k22,
                                           M1,
                                           M2,
                                           &I_temp2),
                           pti->error_message,
                           pti->error_message);
            
                A1 = 0.5 * (3.*I_temp1/k22 - I_temp2);
                A2 = 0.5 * (-I_temp1/k22 + I_temp2);

                result_int = k12*(cos12_sq*A1 + A2);

            } else if (j==3){
                class_call(scalar_prod_one(pti,
                                           j,
                                           i,
                                           d1,
                                           d2,
                                           k22,
                                           M1,
                                           M2,
                                           &I_temp1),
                           pti->error_message,
                           pti->error_message);

                class_call(scalar_prod_one(pti,
                                            j-2,
                                            i+1,
                                            d1,
                                            d2,
                                            k22,
                                            M1,
                                            M2,
                                            &I_temp2),
                                pti->error_message,
                                pti->error_message);
                
                A1 = 0.5 * (5.*I_temp1 - 3.*I_temp2*k22);
                A2 = 0.5 * (-I_temp1 + I_temp2*k22);

                result_int = sqrt(k12/k22)*k12/k22*cos12*(cos12_sq * A1 + 3.*A2);

            } else if (j==4) {
                class_call(scalar_prod_one(pti,
                                           j,
                                           i,
                                           d1,
                                           d2,
                                           k22,
                                           M1,
                                           M2,
                                           &I_temp1),
                           pti->error_message,
                           pti->error_message);

                class_call(scalar_prod_one(pti,
                                           j-2,
                                           i+1,
                                           d1,
                                           d2,
                                           k22,
                                           M1,
                                           M2,
                                           &I_temp2),
                           pti->error_message,
                           pti->error_message);

                class_call(scalar_prod_one(pti,
                                           j-4,
                                           i+2,
                                           d1,
                                           d2,
                                           k22,
                                           M1,
                                           M2,
                                           &I_temp3),
                           pti->error_message,
                           pti->error_message);
                
                k24 = k22*k22;
                k14 = k12*k12;
                A1 = 0.125 * (35.*I_temp1 - 30.*I_temp2*k22   + 3.*I_temp3*k24);
                A2 = 0.125 * (-5.*I_temp1 + 6.*I_temp2*k22    - I_temp3*k24);
                A3 = 0.125 * (I_temp1     - 2.*I_temp2*k22    + I_temp3*k24);

                result_int = k14/k24 * (cos12_sq*cos12_sq*A1 + 6.*cos12_sq*A2 + 3.*A3);
            }


            result += binomial_i*binomial_j*pow(k12, n-i-j)*pow(-1., j)*result_int;
            result_int = class_complex(0., 0.);
        }
    }
    *I_out = result;

    return _SUCCESS_;
}


/** Tensor reduction case 2: solves integrals of the form k1mq^2n1 * k2pq^2n2 / (q^2+M)^d1
 * @param n1                Input: power of numerator k1_vec-q_vec
 * @param n2                Input: power of numerator k2_vec+q_vec
 * @param d1                Input: power of denominator
 * @param k12               Input: wavenumber squared k1^2
 * @param k22               Input: wavenumber squared k2^2
 * @param cos12             Input: cos alpha12, alpha12 is the angle between k1 and k2
 * @param M                 Input: complex mass of propagator
 * @param I_out             Output: pointer to value for the tensor reduction integral (case 2)
 * @return the error status
 */

int tensor_red_two(
                   struct gen_tri_integral *pti,
                   int n1,
                   int n2,
                   int d,
                   double k12,
                   double k22,
                   double cos12,
                   class_complex M,
                   class_complex *I_out
                   ){

    class_complex tad_master_precompute[n1+n2+1], result_temp, result;
    double binomial_i1, binomial_i2, binomial_j1, binomial_j2, k1, k2, k13, k23, k14, k24, cos12_sq, cos12_sqsq;
    int i1, i2, i, j1, j2, j, n;


    double binomial_i1c, binomial_i2c, binomial_j1c, binomial_j2c;


    // initialize some double that are used multiple times
    k1 = sqrt(k12);
    k2 = sqrt(k22);
    k13 = k1*k12;
    k23 = k2*k22;
    k14 = k12*k12;
    k24 = k22*k22;
    cos12_sq = cos12*cos12;
    cos12_sqsq = cos12_sq*cos12_sq;

    result = class_complex(0., 0.);

    // check for illegal cases, TODO: Make these class calls
    if (n1>4 || n2>4){
        printf("Problem in tensor_red_two. Only values n1<=4, n2<=4 are allowed. You have n1=%d, n2=%d", n1, n2);
    }
    if (d<=1+n1+n2){
        printf("Problem in tensor_red_two: Integral is divergent. Requirement for convergence is d1>1+n1+n2. You have n1=%d, n2=%d, d=%d", n1, n2, d);
    }

    // Precompute all relevant Tad_master integrals
    for (j=0; j<=n1+n2; j++){
        class_call(Tad_master(pti,
                              j,
                              d,
                              M,
                              &tad_master_precompute[j]),
                   pti->error_message,
                   pti->error_message);
    }

    // this is the main loop
    // Since different number of scalar product lead to different forms of the solution
    // there are many if statements in this loop
    for (i1=0; i1<=n1; i1++){
        class_call(util_binomial(pti,
                                 n1,
                                 i1,
                                 &binomial_i1),
                   pti->error_message,
                   pti->error_message);

        for (i2=0; i2<=n2; i2++){
            class_call(util_binomial(pti,
                                     n2,
                                     i2,
                                     &binomial_i2),
                       pti->error_message,
                       pti->error_message);

            for (j1=0; j1<=n1-i1; j1++){
                class_call(util_binomial(pti,
                                         n1-i1,
                                         j1,
                                         &binomial_j1),
                           pti->error_message,
                           pti->error_message);

                for (j2=0; j2<=n2-i2; j2++){
                    class_call(util_binomial(pti,
                                             n2-i2,
                                             j2,
                                             &binomial_j2),
                               pti->error_message,
                               pti->error_message);
                    
                    // i is the power of q^2
                    // j is the number of the q_vec in the integral
                    i = i1+i2;
                    j = j1+j2;

                    if (j==0){                    
                        result_temp = tad_master_precompute[i];

                    } else if (j==2){
                        if (j1==0) {
                            result_temp = k22/3. * tad_master_precompute[i+1];
                        } else if (j1==1) {
                            result_temp = k1*k2*cos12/3. * tad_master_precompute[i+1];
                        } else if (j1==2) {
                            result_temp = k12/3. * tad_master_precompute[i+1];
                        }

                    } else if (j==4) {
                        if (j1==0) {
                            result_temp = k24/5. * tad_master_precompute[i+2];
                        } else if (j1==1) {
                            result_temp = k1*k23*cos12/5. * tad_master_precompute[i+2];
                        } else if (j1==2) {
                            result_temp = k12*k22/15. * (2.*cos12_sq + 1.) * tad_master_precompute[i+2];
                        } else if (j1==3) {
                            result_temp = k13*k2*cos12/5. * tad_master_precompute[i+2];
                        } else if (j1==4) {
                            result_temp = k14/5. * tad_master_precompute[i+2];
                        }

                    } else if (j==6) {
                        // since j1, j2 <= 4, the lowest contribution is j1=2, j2=4
                        if (j1==2) {
                            result_temp = k12*k24/35. * (1. + 4.*cos12_sq) * tad_master_precompute[i+3];
                        } else if (j1==3) {
                            result_temp = k13*k23*cos12/35. * (3. + 2.*cos12_sq) * tad_master_precompute[i+3];
                        } else if (j1==4) {
                            result_temp = k14*k22/35. * (1. + 4.*cos12_sq) * tad_master_precompute[i+3];
                        }
                    
                    } else if (j==8) {
                        // since j1, j2 <= 4, the only contribution is j1=j2=4
                        if (j1==4) {
                            result_temp = k14*k14/315. * (3. + 24.*cos12_sq + 8.*cos12_sqsq) * tad_master_precompute[i+4];
                        }
                    } else {
                        // the integral vanishes if j is odd because the integrant is odd in q
                        result_temp = class_complex(0., 0.);
                    }

                    result += binomial_i1*binomial_i2*binomial_j1*binomial_j2 * pow(k12, n1-i1-j1)*pow(k22, n2-i2-j2) * pow(-2., j1) * pow(2., j2) * result_temp;
                }
            }
        }
    }
    *I_out = result;

    return _SUCCESS_;
}


/** Solves a variation of the Tad pole master integral: solves integrals of the form kmq^2n / (q^2+M)^d
 * @param n                 Input: power of numerator (k_vec-q_vec)^2
 * @param d                 Input: power of denominator q^2 + M
 * @param k2                Input: wavenumber squared k^2
 * @param M                 Input: complex mass of propagator
 * @param I_out             Output: pointer to value for the integral
 * @return the error status
 */

 int Tad_var(
             struct gen_tri_integral *pti,
             int n,
             int d,
             double k2,
             class_complex M,
             class_complex *I_out
             ){
    
    class_complex Tad_master_precompute[n+1], result=0;  // in this array we will compute all the relevant Tad_master calls
    double binomial_i, binomial_j;
    int i, j, temp_idx, check_precompute[n+1];

    // fill the arrays with zeros such that the if statement does not have unexpected behaviour
    for (i = 0; i <= n; i++) {
        check_precompute[i] = 0;
        Tad_master_precompute[i] = class_complex(0., 0.);
    }

    
    // use the binomial theorem twice to rewrite (k2 - 2k.q + q2)^n as a double sum over k2^j*(-2)^(n-i-j)*q2^i
    for (i=0; i<=n; i++){
        class_call(util_binomial(pti,
                                 n,
                                 i,
                                 &binomial_i),
                    pti->error_message,
                    pti->error_message);

        for (j=0; j<=n-i; j++){
            if ((n-i-j)%2==0){
                class_call(util_binomial(pti,
                                         n-i,
                                         j,
                                         &binomial_j),
                            pti->error_message,
                            pti->error_message);

                temp_idx = (n+i-j)/2;
                // only compute Tad_master if it hasn't been already computed
                // check_precompute[idx]==1: already has been computed
                // check_precompute[idx]==0: has not been computed yet
                // TODO: this can be a class_Test
                if (temp_idx < 0 || temp_idx > n) {
                    printf("Error: temp_idx out of bounds (%d)\n", temp_idx);
                }
                if (check_precompute[temp_idx]==0){
                    class_call(Tad_master(pti,
                                          temp_idx,
                                          d,
                                          M,
                                          &Tad_master_precompute[temp_idx]),
                                pti->error_message,
                                pti->error_message);
                }
                check_precompute[temp_idx]=1;
                
                // TODO: dont use pow() since it is slow
                result += 1./(n-i-j+1)* pow(k2, double(n-i+j)/2) * binomial_i * binomial_j * pow(2., double(n-i-j)) * Tad_master_precompute[temp_idx];
            }
        }
    }

    *I_out = result;

    return _SUCCESS_;
}

/** Solves an integral with a massive numerator: integrate over (kmq2+M1)^n/(q2+M2)^d
 * @param n                 Input: power of numerator (k_vec-q_vec)^2+M1
 * @param d                 Input: power of denominator q^2 + M
 * @param k2                Input: wavenumber squared k^2
 * @param M1                Input: complex mass of numerator
 * @param M2                Input: complex mass of propagator (denominator)
 * @param I_out             Output: pointer to value for the integral
 * @return the error status
 */

 int massive_num(
                 struct gen_tri_integral *pti,
                 int n,
                 int d,
                 double k2,
                 class_complex M1,
                 class_complex M2,
                 class_complex *I_out
                 ){
    
                    class_complex Tad_var_temp, result = class_complex(0., 0.);
    double binomial;
    int i;
    class_test(n<0 || d<0,
        pti->error_message,
        "Cannot call massive_num for negative values of n or d. You have n=%d, d=%d", n, d);
    for (i=0; i<=n; i++){
        class_call(util_binomial(pti,
                                 n,
                                 i,
                                 &binomial),
                    pti->error_message,
                    pti->error_message);

        class_call(Tad_var(pti,
                           i,
                           d,
                           k2,
                           M2,
                           &Tad_var_temp),
                    pti->error_message,
                    pti->error_message);
        
        result += binomial*pow(M1, n-i)*Tad_var_temp;
    }
    *I_out = result;
    return _SUCCESS_;
}


// ************************************************************************************************
// Utility Fucntions
// ************************************************************************************************

/** Calculate binomial coefficients n over k
 * @param n                     Input: integer
 * @param k                     Input: integer
 * @param n_over_k              Output: pointer to output value
 * @return the error status
 */

int util_binomial(
                  struct gen_tri_integral *pti,
                  int n,
                  int k,
                  double *n_over_k
                  ){
    
    int i, j;
    int C[n + 1][k + 1];

    class_test(n<0 || k<0,
               pti->error_message,
               "Cannot call util_binomial for negative values. You have n=%d, k=%d", n, k);

    class_test(k>n,
               pti->error_message,
               "Cannot call util_binomial for k>n. You have n=%d, k=%d", n, k);

    for (i = 0; i <= n; i++) {
        for (j = 0; j <= k && j <= i; j++) {
            if (j == 0 || j == i)
                C[i][j] = 1;
            else
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
        }
    }

    *n_over_k = C[n][k]; 
    return _SUCCESS_;
}


/** Calculate trinomial coefficients (n over k) * (n-k over j)
 * @param n                     Input: integer
 * @param k                     Input: integer
 * @param j                     Input: integer
 * @param out                   Output: pointer to output value
 * @return the error status
 */

int util_trinomial(
                   struct gen_tri_integral *pti,
                   int n,
                   int k,
                   int j,
                   double *out
                   ){

    return _SUCCESS_;
}


/** Calculate arctan*prefactor, TODO: reference equation in thesis or paper
 * @param x                     Input: integration variable, evaluated at x=0 and x=1
 * @param y1                    Input: complex variable z_+
 * @param y2                    Input: complex variable z_-
 * @param x0                    Input: complex variable, evaluated for x_- or x_+
 * @param out                   Output: pointer to output value
 * @return the error status
 */

int util_antideriv(
                   struct gen_tri_integral *pti,
                   double x_in,
                   class_complex y1,
                   class_complex y2,
                   class_complex x0,
                   class_complex *out
                   ){
    double CHOP_TOL = 1.e-14;
    class_complex x = class_complex(x_in, 0.);

    if (cabs(y1-x0)<CHOP_TOL || cabs(x0-y2)<CHOP_TOL){
        printf("\nERROR in util_antideriv: division by zero. either x0_re_c=y1 or x0_re_c=y2.\n");
    }

    // case where x0 = y2 = 0 or 1
    if (cabs(y2 - x0)<CHOP_TOL){
        if (cabs(y2-x)<CHOP_TOL){
            printf("\nspecial case: antideriv = 0");
            *out = class_complex(0., 0.);
            return _SUCCESS_;
        }
        *out = 2.*sqrt(x-y1)/(-x0+y1)/sqrt(x-y2);
        return _SUCCESS_;
    }

    // case where x0 = y2 = 0 or 1
    if (cabs(x0-y1)<CHOP_TOL){
        printf("WARNING: switching var in Antideriv");
        class_call(util_antideriv(pti,
                                  x_in,
                                  y2,
                                  y1,
                                  x0,
                                  out),
                    pti->error_message,
                    pti->error_message);
        return _SUCCESS_;
    }

    class_complex prefac = 2./(sqrt(-x0+y1)*sqrt(x0-y2));
    class_complex temp = sqrt(x-y1)*sqrt(x0-y2)/sqrt(-x0+y1);
    class_complex LimArcTan;

    if (x_in==1. && cabs(y2-class_complex(1., 0.))<CHOP_TOL ){
        LimArcTan = class_complex(0., 1.)*sqrt(-temp*temp) * _PI_/(2.*temp);
        *out = prefac * LimArcTan;
        return _SUCCESS_;
    }
    if (x_in==0. && cabs(y2)<CHOP_TOL){
        LimArcTan = sqrt(temp*temp) * _PI_/(2.*temp);
        *out = prefac * LimArcTan;
        return _SUCCESS_;
    }

    class_complex atan_ = atan(temp/sqrt(x-y2));
    *out = prefac * atan_;

    return _SUCCESS_;
}


/** Calculate TODO: not completely sure what this is exactly
 * I thought this is s(z+, z-), but I am confused about the 1/sqrt(a) since a<0
 * @param a                     Input: integer
 * @param y1                    Input: complex variable z_+
 * @param y2                    Input: complex variable z_-
 * @param out                   Output: pointer to output value
 * @return the error status
 */

int util_prefactor(
                   struct gen_tri_integral *pti,
                   double a_in,
                   class_complex y1,
                   class_complex y2,
                   class_complex *out
                   ){
    class_complex y1_re, y1_im, y2_re, y2_im, m_y1_re, m_y2_re, a, temp_ay1y2;
    double CHOP_TOL;

    CHOP_TOL = 1.e-14;

    a = class_complex(a_in, 0.);

    y1_re = class_complex(creal(y1),    0.);
    y1_im = class_complex(0,            cimag(y1));
    y2_re = class_complex(creal(y2),    0.);
    y2_im = class_complex(0.,           cimag(y2));

    // the m_ is for "minus ..."
    // sqrt(-real - I*0.) = - I*sqrt(real)          -> incorrect since imag part is exactly zero
    // sqrt(-real + I*0.) = + I*sqrt(real)          -> correct
    m_y1_re = class_complex(-creal(y1),    0.);
    m_y2_re = class_complex(-creal(y2),    0.);



    // compute s(y1, y2) = sqrt(-y1)*sqrt(-y2)/sqrt(a*y1*y2)
    // special cases if the real or imaginary part of y1 and/or y2 is close to zero

    // The first case requires a special treatment since both imaginary parts are zero
    // It is important that we don't calculate sqrt(real_number - I*0.) as explained above
    if (cabs(y2_im) < CHOP_TOL && cabs(y1_im) < CHOP_TOL){
        if (cabs(y1_re) >= CHOP_TOL && cabs(y2_re) >= CHOP_TOL){
            temp_ay1y2 = class_complex(a_in*creal(y1)*creal(y2), 0.);
            *out = sqrt(m_y1_re)*sqrt(m_y2_re)/sqrt(temp_ay1y2);
            return _SUCCESS_;
        }
        if (cabs(y1_re) < CHOP_TOL && cabs(y2_re) >= CHOP_TOL){
            temp_ay1y2 = class_complex(-a_in*creal(y2), 0.);
            *out = sqrt(m_y2_re)/sqrt(temp_ay1y2);
            return _SUCCESS_;
        }
        if (cabs(y1_re) >= CHOP_TOL && cabs(y2_re) < CHOP_TOL){
            temp_ay1y2 = class_complex(-a_in*creal(y1), 0.);
            *out = sqrt(m_y1_re)/sqrt(temp_ay1y2);
            return _SUCCESS_;
        }
        if (cabs(y1_re) < CHOP_TOL && cabs(y2_re) < CHOP_TOL){
            temp_ay1y2 = class_complex(a_in, 0.);
            *out = 1./sqrt(temp_ay1y2);
            return _SUCCESS_;
        }
    } else if (cabs(y2_im) >= CHOP_TOL && cabs(y1_im) < CHOP_TOL){
        if (cabs(y1_re) >= CHOP_TOL){
            *out = sqrt(m_y1_re)*sqrt(-y2)/sqrt(a*y1_re*y2);
            return _SUCCESS_;
        }
        if (cabs(y1_re) < CHOP_TOL){
            *out = sqrt(-y2)/sqrt(-a*y2);
            return _SUCCESS_;
        }
    } else if (cabs(y2_im) < CHOP_TOL && cabs(y1_im) >= CHOP_TOL){
        if (cabs(y2_re) > CHOP_TOL){
            *out = sqrt(-y1)*sqrt(m_y2_re)/sqrt(a*y1*y2_re);
            return _SUCCESS_;
        }
        if (cabs(y2_re) < CHOP_TOL){
            *out = sqrt(-y1)/sqrt(-a*y1);
            return _SUCCESS_;
        }
    } else {
        // case where cabs(y2_im) >= CHOP_TOL && cabs(y1_im) >= CHOP_TOL
        *out = sqrt(-y1)*sqrt(-y2)/sqrt(a*y1*y2);
        return _SUCCESS_;
    }
}


/** Calculate F_int
 * @param R2                    Input: negative double
 * @param y1                    Input: complex double z_+
 * @param y2                    Input: complex double z_-
 * @param x0                    Input: complex double
 * @param out                   Output: pointer to output value
 * @return the error status
 */

int util_F_int(
               struct gen_tri_integral *pti,
               double R2,
               class_complex y1,
               class_complex y2,
               class_complex x0,
               class_complex *out
               ){
    class_complex cutx0, cut, derivx0, deriv, x0_re_c, atan_arg, prefactor, antideriv0, antideriv1;
    double CHOP_TOL, y1_re, y1_im, y2_re, y2_im, x0_re, x0_im, y1_re2, y1_im2, y2_re2, y2_im2, x0_re2, x0_im2, a, b, c;
    double xcross_temp, xsol[2], xbranch[2]; // at most two branch cuts, so length of list is two
    int idx, signx0, sign, n_branchpoints;

    CHOP_TOL = 1.e-10;

    // contribution of a branch cut crossing at the boundary <=> B=0 (only pi/2 gap)
    cutx0   = class_complex(0., 0.);

    // regular contribution from branch cut crossing (pi gap)
    cut     = class_complex(0., 0.);

    // initialize values for the lists that information on the branch cuts.
    // branch cut crossings can only occur between 0 and 1, so the value -1 indicates no branch cut crossing
    for (idx = 0; idx<2; idx++){
        xsol[idx] = -1.;
        xbranch[idx] = -1.;
    }

    // initialize number of branch cut crossings
    n_branchpoints = 0;

    // split y1, y2, and x0 into real and imag parts
    y1_re = creal(y1);
    y1_im = cimag(y1);
    y2_re = creal(y2);
    y2_im = cimag(y2);
    x0_re = creal(x0);
    x0_im = cimag(x0);

    // squared
    y1_re2 = y1_re*y1_re;
    y1_im2 = y1_im*y1_im;
    y2_re2 = y2_re*y2_re;
    y2_im2 = y2_im*y2_im;
    x0_re2 = x0_re*x0_re;
    x0_im2 = x0_im*x0_im;


    // create special cases if the imaginary part is close to zero
    if (abs(cimag(y1)) < CHOP_TOL*fmax(1.0, cabs(y1))){
        y1 = class_complex(creal(y1), 0.);
    }
    if (abs(cimag(y2)) < CHOP_TOL*fmax(1.0, cabs(y2))){
        y2 = class_complex(creal(y2), 0.);
    }
    if (abs(cimag(x0)) < CHOP_TOL*fmax(1.0, cabs(x0))){
        x0 = class_complex(creal(x0), 0.);
    }

    a = y1_im*x0_re-y2_im*x0_re-x0_im*y1_re+y2_im*y1_re+x0_im*y2_re-y1_im*y2_re;
    b = -x0_im2*y1_im + x0_im*y1_im2+x0_im2*y2_im-y1_im2*y2_im-x0_im*y2_im2+y1_im*y2_im2-y1_im*x0_re2+y2_im*x0_re2+x0_im*y1_re2-y2_im*y1_re2-x0_im*y2_re2+y1_im*y2_re2;
    c = y1_im2*y2_im*x0_re - y1_im*y2_im2*x0_re-x0_im2*y2_im*y1_re + x0_im*y2_im2*y1_re-y2_im*x0_re2*y1_re + y2_im*x0_re*y1_re2 + x0_im2*y1_im*y2_re-x0_im*y1_im2*y2_re+y1_im*x0_re2*y2_re-x0_im*y1_re2*y2_re-y1_im*x0_re*y2_re2+x0_im*y1_re*y2_re2;

    // case B=0 <=> 0 < x0_re < 1 and x0_im=0 <=> branch cut crossing at the boundary of the branch cut
    // determine the sign of the extra contribution by taking the derivative of the argument of arctan
    // w.r.t. x
    if (0. < x0_re && x0_re < 1. && abs(x0_im)<CHOP_TOL*fmax(1.0, cabs(x0))){
        // printf("\nBranch cut crossing at B=0!!!");
        x0_re_c = class_complex(x0_re, 0.);
        if (cabs(y1-x0_re_c)<CHOP_TOL || cabs(x0_re_c-y2)<CHOP_TOL){
            // printf("\nERROR in util_F_int: division by zero. either x0_re_c=y1 or x0_re_c=y2.\n");
        }
        derivx0 = (y1 - y2)/2./sqrt(-(x0_re_c-y1)*(x0_re_c-y1))/(x0_re_c-y2);
        if (creal(derivx0) < 0){
            signx0 = 1;
        } else {
            signx0 = -1;
        }
        cutx0 = signx0*_PI_/(sqrt(-x0_re_c+y1)*sqrt(x0_re_c-y2));
        // printf("\ncutx0 = %.10f + %.10fj", creal(cutx0), cimag(cutx0));
    }

    // other possible branch cut crossings:
    // x_1,2 = 1/(2a) * (-b +- sqrt( b^2 -4ac ) ) -> needs to be between 0 and 1

    // special case if a = 0
    if (abs(a)<CHOP_TOL){
        if (abs(b)>CHOP_TOL){
            if ( 0<-c/b && -c/b<1){
                xsol[0] = -c/b;
            }
        }
    } else {
        // check that the sqrt is not imagniary
        if (b*b-4*a*c > CHOP_TOL){
            xcross_temp = (-b + sqrt(b*b-4*a*c))/(2.*a);           // first potential branch cut crossing
            // check that the solutions are between 0 and 1
            // and that x≠x0 since this corresponds to B=0
            if ((CHOP_TOL < xcross_temp && xcross_temp < 1.-CHOP_TOL) && !(abs(xcross_temp-x0_re) < CHOP_TOL && x0_im < CHOP_TOL)){
                xsol[0] = xcross_temp;
            }
            xcross_temp = (-b - sqrt(b*b-4*a*c))/(2.*a);           // second potential branch cut crossing
            if ((CHOP_TOL < xcross_temp && xcross_temp < 1.-CHOP_TOL) && !(abs(xcross_temp-x0_re) < CHOP_TOL && x0_im < CHOP_TOL)){
                xsol[1] = xcross_temp;
            }
        }
    }
    for (idx=0; idx<2; idx++){
        if (xsol[idx] != -1.){
            if (cabs(y1-x0)<CHOP_TOL || cabs(x0-y2)<CHOP_TOL){
                // printf("\nERROR in util_F_int: division by zero. either x0=y1 or x0=y2\n you have y1-x0=%.10f + %.10fj, y2-x0=%.10f + %.10fj, cabs(y1-x0)=%.10f , cabs(y2-x0)=%.10f , CHOP_TOL = %.10f.\n", creal(x0-y1), cimag(x0-y1), creal(x0-y2), cimag(x0-y2), cabs(x0-y1), cabs(x0-y2), CHOP_TOL);
            }
            atan_arg = sqrt(xsol[idx]-y1)*sqrt(x0-y2)/(sqrt(-x0+y1)*sqrt(xsol[idx]-y2));
            // critereon for branch cut crossing: argument of atan lies on the imaginary axis either between -i*infty to -i or between i to i*infty
            if (cabs(atan_arg)>1. && abs(creal(atan_arg))<CHOP_TOL){
                n_branchpoints += 1;
            }
        }
    }

    // calculate the contributions from the branch cut crossings
    // only relevant contribution if there is exactly one branch cut crossing
    // printf("\nn_branchpoints=%d",n_branchpoints);
    if (n_branchpoints==1){
        // printf("\nbranch cut crossing!!!\n");
        // printf("xsol[0] = %.10f, xsol[1] = %.10f", xsol[0], xsol[1]);
        for (idx=0; idx<2; idx++){
            if (xsol[idx] != -1.){
                if (cabs(y1-xsol[idx])<CHOP_TOL || cabs(xsol[idx]-y2)<CHOP_TOL){
                    // printf("\nERROR in util_F_int: division by zero. either xsol[idx]=y1 or xsol[idx]=y2.\n");
                }
                deriv = sqrt(x0-y2)/sqrt(-x0+y1)*(1./(2.*sqrt(xsol[idx]-y1)*sqrt(xsol[idx]-y2)) - sqrt(xsol[idx]-y1)/(2.*(xsol[idx]-y2)*sqrt(xsol[idx]-y2)));
            }
        }
        
        // determine the sign of the extra contribution by taking the derivative of the argument of arctan
        // w.r.t. x
        if (creal(deriv)<0.){
            sign = 1.;
        } else {
            sign = -1.;
        }
        if (cabs(y1-x0)<CHOP_TOL || cabs(x0-y2)<CHOP_TOL){
            // printf("\nERROR in util_F_int: division by zero. either x0=y1 or x0=y2.\n");
        }
        cut = sign*_PI_*2./(sqrt(-x0+y1)*sqrt(x0-y2));

        // printf("\ncut = %.10f + %.10fj", creal(cut), cimag(cut));

        // check for mistakes
        if (xsol[0] == -1. && xsol[1] == -1.){
            // printf("\n ERROR: number of branch points is supposed to be equal to one, but both entries of xsol are -1\n");
            // printf("xsol[0] = %.10f, xsol[1] = %.10f", xsol[0], xsol[1]);
        }
            
        
    } else {
        // I think unnecessary since we initialized it as zero
        cut = class_complex(0., 0.);
    }
    class_call(util_prefactor(pti, 
                              R2,
                              y1,
                              y2,
                              &prefactor),
                pti->error_message,
                pti->error_message);


    class_call(util_antideriv(pti,
                              0,
                              y1,
                              y2,
                              x0,
                              &antideriv0),
                pti->error_message,
                pti->error_message);

    class_call(util_antideriv(pti,
                              1,
                              y1,
                              y2,
                              x0,
                              &antideriv1),
                pti->error_message,
                pti->error_message);

    *out = prefactor * _PI_*_PI_/2. * (cut + cutx0 + antideriv1 - antideriv0);
    
    return _SUCCESS_;
}


/** Calculate the T_master contributions
 * @param y                     Input: evaluation of the integral at 0 and 1
 * @param k12                   Input: wavenumber k1^2
 * @param k22                   Input: wavenumber k2^2
 * @param k32                   Input: wavenumber k3^2
 * @param M1                    Input: complex mass
 * @param M1                    Input: complex mass
 * @param M1                    Input: complex mass
 * @param out                   Output: pointer to output value
 * @return the error status
 */
int util_Tmaster_contr(
                       struct gen_tri_integral *pti,
                       int y,
                       double k12,
                       double k22,
                       double k32,
                       class_complex M1,
                       class_complex M2,
                       class_complex M3,
                       class_complex *out
                       ){
    class_complex Num0, R1, R0, S1, S1_overk22, S0, DiscS_sqrt, DiscR_sqrt, solS1, solS2, solR1, solR2, c1, c2, Fint_temp1, Fint_temp2;
    double Num1, R2, S2, k1, k2, k3, k14, k24, k34, CHOP_TOL;

    class_complex R1_overk32, R1_overk12, tmp_disc1, tmp_disc1_overk22, DiscR_sqrt_overk32, DiscR_sqrt_overk12, DiscS_sqrt_overk22, Num0_overk22;
    double S2_overk22, Num1_overk22;


    // antiderivative of the master integral integrand in y 
    
    CHOP_TOL = 1.e-14;


    k1 = sqrt(k12);
    k2 = sqrt(k22);
    k3 = sqrt(k32);
    k14 = k12*k12;
    k24 = k22*k22;
    k34 = k32*k32;

    // ---------------- Eq. 5.79 ----------------
    // N_1, N_0
    // Num1 = 4.*k22*y+2.*k12-2.*k22-2.*k32;
    // Num0 = -4.*k22*y+2.*M2-2.*M3+2.*k22;

    // R_2, R_1, R_0
    // R2 = -k12*y+k32*y-k32;
    // R1 = -M2*y+M3*y+k12*y-k32*y+M1-M3+k32;
    // R0 = M2*y-M3*y+M3;

    // implement these formulas separately to avoid minimize erros
    if (y==0){
        R2 = -k32;
        R1 = M1-M3+k32;
        R0 = M3;

        DiscR_sqrt_overk32 = sqrt((1. + (M1 - M3)/k32)*(1. + (M1 - M3)/k32) + 4.*M3/k32);
        R1_overk32 = (M1-M3)/k32 + 1.;
        solR1 = -0.5*(-R1_overk32 + DiscR_sqrt_overk32);                 // z_+
        solR2 = -0.5*(-R1_overk32 - DiscR_sqrt_overk32);                 // z_-

        Num1 = 2.*(k12-k22-k32);
        Num0 = 2.*(M2-M3+k22);

        Num1_overk22 = 2.*(k12/k22 - 1. - k32/k22);
        Num0_overk22 = 2.*((M2-M3)/k22 + 1.);
    } else if (y==1){
        R2 = -k12;
        R1 = M1-M2+k12;
        R0 = M2;
        DiscR_sqrt_overk12 = sqrt((1. + (M1 - M2)/k12)*(1. + (M1 - M2)/k12) + 4.*M2/k12);
        R1_overk12 = (M1-M2)/k12 + 1.;
        solR1 = -0.5*(-R1_overk12 + DiscR_sqrt_overk12);                 // z_+
        solR2 = -0.5*(-R1_overk12 - DiscR_sqrt_overk12);                 // z_-

        Num1 = 2.*(k12+k22-k32);
        Num0 = 2.*(M2-M3-k22);

        Num1_overk22 = 2.*(k12/k22 + 1. - k32/k22);
        Num0_overk22 = 2.*((M2-M3)/k22 - 1.);
    } else {
        printf("Error in util_Tmaster_contr. Case not considered. y=%d",y);
    }


    // S_2, S_1, S_0
    S2 = -k14+2.*k12*k22+2.*k12*k32-k24+2.*k22*k32-k34;
    S1 = -4.*M1*k22-2.*M2*k12+2.*M2*k22+2.*M2*k32+2.*M3*k12+2.*M3*k22-2.*M3*k32-2.*k12*k22+2.*k24-2.*k22*k32;
    S0 = -M2*M2+2.*M2*M3-2.*M2*k22-M3*M3-2.*M3*k22-k24;

    S2_overk22 = -k12/k22 *k12+2.*k12+2.*k12/k22*k32-k22+2.*k32-k32/k22 *k32;
    S1_overk22 = -4.*M1-2.*M2*k12/k22+2.*M2+2.*M2*k32/k22+2.*M3*k12/k22+2.*M3-2.*M3*k32/k22-2.*k12+2.*k22-2.*k32;

    // Diakr(a,b,c)= b^2 - 4*a*c -> defined after ---------------- Eq. 5.94 ----------------
    // Sqrt of Diakr

    // DiscS_sqrt = sqrt(S1*S1 - 4.*S2*S0); // Old implementation 

    // New implementation -> does the same thing, hoping that is has better numerical precision
    tmp_disc1 = -k24 + k22*(k32 + 2.*M1 - M2 - M3) + k12*(k22 + M2 - M3) + k32*(-M2 + M3);
    tmp_disc1_overk22 = -k22 + k32 + k12 + 2.*M1 - M2 - M3 + k12/k22*(M2 - M3) + k32/k22*(-M2 + M3);
    DiscS_sqrt = 2.*sqrt(tmp_disc1*tmp_disc1 - (k1 - k2 - k3)*(k1 + k2 - k3)*(k1 - k2 + k3)*(k1 + k2 + k3)*(k24 + (M2 - M3)*(M2 - M3) + 2.*k22*(M2 + M3)));
    DiscS_sqrt_overk22 = 2.*sqrt(tmp_disc1_overk22*tmp_disc1_overk22 - (k1/k2 - 1. - k3/k2) * (k1/k2 + 1. - k3/k2) * (k1/k2 - 1. + k3/k2) * (k1/k2 + 1. + k3/k2)*(k24 + (M2 - M3)*(M2 - M3) + 2.*k22*(M2 + M3)));


    if (cabs(DiscS_sqrt)<CHOP_TOL){
        printf("\nERROR in util_Tmaster_contr: DiscS_sqrt=%.10f + %.10fj\n", creal(DiscS_sqrt), cimag(DiscS_sqrt));
    }
    if (abs(S2)<CHOP_TOL){
        printf("\nERROR in util_Tmaster_contr: S2 = 4*k12*k22*(cos12^2 - 1)=%.10f\n", S2);
        printf("\nk12 = %.10f, k22 = %.10f, cos12 = %.10f", k12, k22, (k32-k12-k12)/(2.*sqrt(k12)*sqrt(k22)));
    }

    //a*x^2 + b*x + c = 0 has the solutions
    // x_1,2 = (-b +- sqrt(b^2 - 4*a*c))/(2*a)

    // ---------------- Eq. 5.80 ----------------
    // solS1 and solS2 are the solution of 
    // DeltaS2*x^2 + DeltaS1*x + DeltaS0 = 0
    solS1 = 0.5*(-S1_overk22+DiscS_sqrt_overk22)/S2_overk22;                // x_+
    solS2 = 0.5*(-S1_overk22-DiscS_sqrt_overk22)/S2_overk22;                  // x_-

    c2 = -(Num1_overk22*solS2 + Num0_overk22)/DiscS_sqrt_overk22;                    // c1
    c1 =  (Num1_overk22*solS1 + Num0_overk22)/DiscS_sqrt_overk22;                        // c2
        
    DiscR_sqrt = sqrt(R1*R1 - 4.*R2*R0);    // sqrt( DeltaR1^2 - 4*DeltaR2*DeltaR0 )

    if (cabs(DiscR_sqrt)<CHOP_TOL){
        printf("\nERROR in util_Tmaster_contr: DiscR_sqrt=0\n");
    }
    if (abs(R2)<CHOP_TOL){
        printf("\nERROR in util_Tmaster_contr: R2=0\n");
    }


    // solR1 = 0.5*(-R1+DiscR_sqrt)/R2;                 // z_+
    // solR2 = 0.5*(-R1-DiscR_sqrt)/R2;                 // z_-

    // ---------------- Eq. 5.82 ----------------
    if (cabs(c1) < CHOP_TOL){
        class_call(util_F_int(pti,
                              R2,
                              solR1,
                              solR2,
                              solS2,
                              &Fint_temp2),
                    pti->error_message,
                    pti->error_message);
        
        *out = c2*Fint_temp2;
        return _SUCCESS_;

    } else if (cabs(c2) < CHOP_TOL){
        class_call(util_F_int(pti,
                              R2,
                              solR1,
                              solR2,
                              solS1,
                              &Fint_temp1),
                    pti->error_message,
                    pti->error_message);

        *out = c1*Fint_temp1;
        return _SUCCESS_;

    } else {
        class_call(util_F_int(pti,
                              R2,
                              solR1,
                              solR2,
                              solS1,
                              &Fint_temp1),
                    pti->error_message,
                    pti->error_message);

        class_call(util_F_int(pti,
                              R2,
                              solR1,
                              solR2,
                              solS2,
                              &Fint_temp2),
                    pti->error_message,
                    pti->error_message);
                    

        *out = c1*Fint_temp1 + c2*Fint_temp2;
        return _SUCCESS_;

    }
}

/** Solves general trianlge integral: L(n1, d1, n2, d2, n3, d3) = int k1mq^2n1 q^2n2 k2pq^2n3 / (k1mq^2+M1)^d1 / (q^2+M2)^d2 / (k2pq^2+M3)^d3
 * @param n1                Input: power of numerator (k1_vec-q_vec)^2
 * @param d1                Input: power of denominator (k1_vec-q_vec)^2 + M1
 * @param n2                Input: power of numerator q^2
 * @param d2                Input: power of denominator q^2 + M2
 * @param n3                Input: power of numerator (k2_vec+q_vec)^2
 * @param d3                Input: power of denominator (k2_vec+q_vec)^2 + M3
 * @param k12               Input: wavenumber squared k1^2
 * @param k22               Input: wavenumber squared k2^2
 * @param k32               Input: wavenumber squared k3^2
 * @param M1                Input: complex mass of propagator
 * @param M2                Input: complex mass of propagator
 * @param M3                Input: complex mass of propagator
 * @param L_out             Output: pointer to value for the generalized triangle integral
 * @return the error status
 */

int util_L_step(
                struct gen_tri_integral *pti,
                doublets_L idxs[],
                int idx,
                double k12,
                double k22,
                double k32,
                class_complex Ms[],
                class_complex *L_out
                ){

    class_complex M0, tmp_Ms[3], L_temp1, L_temp2;
    doublets_L tmp_idxs[3];
    int i, op_idx[3];

    M0 = class_complex(0., 0.);

    // tmp_idxs are modified indices that are used if one of the d_i==0 while n_i!=0. Here they are initialized to the default values.
    // op_idxs are operations of the indices, which are used in the recursion relations. Here they are initialized to zero
    for (i=0; i<3; i++){
        tmp_idxs[i].n = idxs[i].n;
        tmp_idxs[i].d = idxs[i].d;
        tmp_Ms[i] = Ms[i];
        op_idx[i] = 0;
    }

    // if one of the denominator exponents is zero, redefine the numerator exponent (n_i -> -d_i) and set the corresponding mass to zero (M0)
    if (idxs[idx].d==0 && idxs[idx].n!=0){
        tmp_idxs[idx].n = 0;
        tmp_idxs[idx].d = -idxs[idx].n;
        tmp_Ms[idx] = M0;

        class_call(L_recursion(pti,
                               tmp_idxs[0].n,
                               tmp_idxs[0].d,
                               tmp_idxs[1].n,
                               tmp_idxs[1].d,
                               tmp_idxs[2].n,
                               tmp_idxs[2].d,
                               k12,
                               k22,
                               k32,
                               tmp_Ms[0],
                               tmp_Ms[1],
                               tmp_Ms[2],
                               L_out),
                    pti->error_message,
                    pti->error_message);
        return _SUCCESS_;

    } else if (idxs[idx].n>0){
        op_idx[idx] = 1;

        // recursion relation of n_i>0
        class_call(L_recursion(pti,
                               tmp_idxs[0].n - op_idx[0],
                               tmp_idxs[0].d - op_idx[0],
                               tmp_idxs[1].n - op_idx[1],
                               tmp_idxs[1].d - op_idx[1],
                               tmp_idxs[2].n - op_idx[2],
                               tmp_idxs[2].d - op_idx[2],
                               k12,
                               k22,
                               k32,
                               Ms[0],
                               Ms[1],
                               Ms[2],
                               &L_temp1),
                    pti->error_message,
                    pti->error_message);

        class_call(L_recursion(pti,
                               tmp_idxs[0].n - op_idx[0],
                               tmp_idxs[0].d,
                               tmp_idxs[1].n - op_idx[1],
                               tmp_idxs[1].d,
                               tmp_idxs[2].n - op_idx[2],
                               tmp_idxs[2].d,
                               k12,
                               k22,
                               k32,
                               Ms[0],
                               Ms[1],
                               Ms[2],
                               &L_temp2),
                    pti->error_message,
                    pti->error_message);

        *L_out = L_temp1 - Ms[idx]*L_temp2;
        return _SUCCESS_;

    } else if (idxs[idx].n<0){
        op_idx[idx] = 1;

        // recursion relation of n_i<0
        class_call(L_recursion(pti,
                               tmp_idxs[0].n,
                               tmp_idxs[0].d - op_idx[0],
                               tmp_idxs[1].n,
                               tmp_idxs[1].d - op_idx[1],
                               tmp_idxs[2].n,
                               tmp_idxs[2].d - op_idx[2],
                               k12,
                               k22,
                               k32,
                               Ms[0],
                               Ms[1],
                               Ms[2],
                               &L_temp1),
                    pti->error_message,
                    pti->error_message);

        class_call(L_recursion(pti,
                               tmp_idxs[0].n + op_idx[0],
                               tmp_idxs[0].d,
                               tmp_idxs[1].n + op_idx[1],
                               tmp_idxs[1].d,
                               tmp_idxs[2].n + op_idx[2],
                               tmp_idxs[2].d,
                               k12,
                               k22,
                               k32,
                               Ms[0],
                               Ms[1],
                               Ms[2],
                               &L_temp2),
                    pti->error_message,
                    pti->error_message);

        *L_out = (L_temp1 - L_temp2)/Ms[idx];
        return _SUCCESS_;
    } else {
        printf("Error in util_L_step: case not considered, n1=%d, n2=%d, n3=%d, d1=%d, d2=%d, d3=%d", tmp_idxs[0].n, tmp_idxs[1].n, tmp_idxs[2].n, tmp_idxs[0].d, tmp_idxs[1].d, tmp_idxs[2].d);
        return _FAILURE_;
    }
}