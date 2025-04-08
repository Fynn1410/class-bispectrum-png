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
    if (M1==M0 and M2==M0 and M3==M0){
        *T_out = _PI_*_PI_*_PI_/sqrt(k12)/sqrt(k22)/sqrt(k32);
    } else {
        class_complex TrMxy0, TrMxy1;
        //util_tri_master(0., k12, k22, k32, M1, M2, M3, &TrMxy0);
        //util_tri_master(1., k12, k22, k32, M1, M2, M3, &TrMxy1);
        //*T_out = TrMxy1 - TrMxy0;
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

    if (cimag(arglog0) > 0 and cimag(arglog1) < 0){
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
            nu1_c = class_complex(d1-1, 0.);
            nu2_c = class_complex(d2, 0.);
            c1 = ((2.*M2-k1s)/nu1_c*Ndim - 2.*M2 + (k1s*nu2_c)/nu1_c)/jac;
            c2 = k1s/jac;
            c3 = (-2.*M2/nu1_c*nu2_c)/jac;
        } else if (d2>1){
            nu1 = d1;
            nu2 = d2-1;
            nu1_c = class_complex(d1, 0.);
            nu2_c = class_complex(d2-1, 0.);
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
                   int m,
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
