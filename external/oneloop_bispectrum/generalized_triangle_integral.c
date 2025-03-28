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

 int T_master(double k12, double k22, double k32, class_complex M1, class_complex M2, class_complex M3, class_complex *T_out){
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

int B_master(double k2, class_complex M1, class_complex M2, class_complex *B_out){
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

int Tad_master(int n, int d, class_complex M, class_complex *Tad_out){
    class_complex temp;

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

int L_recursion(int n1, int d1, int n2, int d2, int n3, int d3, double k12, double k22, double k32, class_complex M1, class_complex M2, class_complex M3, class_complex *L_out){

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

int T_recursion(int d1, int d2, int d3, double k12, double k22, double k32, class_complex M1, class_complex M2, class_complex M3, class_complex *T_out){

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

int B_recursion(int d1, int d2, double k2, class_complex M1, class_complex M2, class_complex *B_out){
    // TODO: any special cases if the masses are zero?
    int nu1, nu2;
    class_complex B_temp1, B_temp2, B_temp3, c1, c2, c3, jac, k1s, Ndim, nu1_c, nu2_c;
    // terminate the recursion for the following cases:
    if (d1==1 && d2==1){
        // Bubble master integral
        B_master(k2, M1, M2, B_out);
    } else if (d1==0){
        if (d2==1){
            printf("Error in B_recursion: d1=0, d2=1, divergent integral");
        }
        // Tad pole master integral with numerator n=1
        Tad_master(0, d2, M2, B_out);
    } else if (d2==0){
        if (d1==1){
            printf("Error in B_recursion: d1=1, d2=0, divergent integral");
        }
        // Tad pole master integral with numerator n=1
        Tad_master(0, d1, M1, B_out);
    } else if (d1<0 && d2<0){
        // TODO: this integral is divergent, why set it to zero? maybe because it cannot appear? If so -> class_test
        *B_out = class_complex(0., 0.);
    } else if (d1<0 && d2>0){
        if (2-d1>d2){
            printf("Error in B_recursion: d1=%d, d2=%d, divergent integral", d1, d2);
        }
        // massive numerator integral
        massive_num(-d1, d2, k2, M1, M2, B_out);
    } else if (d1>0 && d2<0){
        if (2-d2>d1){
            printf("Error in B_recursion: d1=%d, d2=%d, divergent integral", d1, d2);
        }
        // TODO: what it 2+|d2|>d1 (2 coming from the jacobian) then the integral is divergent. return 0?
        // massive numerator integral
        massive_num(-d2, d1, k2, M2, M1, B_out);
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
        B_recursion(nu1, nu2, k2, M1, M2, &B_temp1);
        B_recursion(nu1+1, nu2-1, k2, M1, M2, &B_temp2);
        B_recursion(nu1-1, nu2+1, k2, M1, M2, &B_temp3);
        *B_out = c1*B_temp1 + c2*B_temp2 + c3*B_temp3;
    } else {
        printf("Error in B_recursion: case not considered");
    }
        
    return _SUCCESS_;
}


/** Tensor reduction case 1: solves integrals of the form k1mq^2n1 / (q^2+M1)^d1 / (k2pq^2+M2)^d2
 * @param n                 Input: power of numerator k1_vec-q_vec
 * @param d1                Input: power of denominator q^2+M1
 * @param d2                Input: power of denominator (q2_vec+q_vec)^2 + M2
 * @param k12               Input: wavenumber squared k1^2
 * @param k22               Input: wavenumber squared k2^2
 * @param M1                Input: complex mass of propagator
 * @param M2                Input: complex mass of propagator
 * @param I_out             Output: pointer to value for the tensor reduction integral (case 1)
 * @return the error status
 */

int tensor_red_one(int n, int d1, int d2, double k12, double k22, class_complex M1, class_complex M2, class_complex *I_out){

    return _SUCCESS_;
}


/** Tensor reduction case 2: solves integrals of the form k1mq^2n1 * k2pq^2n2 / (q^2+M)^d1
 * @param n1                Input: power of numerator k1_vec-q_vec
 * @param n2                Input: power of numerator k2_vec+q_vec
 * @param d1                Input: power of denominator
 * @param k12               Input: wavenumber squared k1^2
 * @param k22               Input: wavenumber squared k2^2
 * @param M                 Input: complex mass of propagator
 * @param I_out             Output: pointer to value for the tensor reduction integral (case 2)
 * @return the error status
 */

int tensor_red_two(int n1, int n2, int d1, double k12, double k22, class_complex M, class_complex *I_out){

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

 int Tad_var(int n, int d, double k2, class_complex M, class_complex *I_out){
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
        util_binomial(n, i, &binomial_i);
        for (j=0; j<=n-i; j++){
            if ((n-i-j)%2==0){
                util_binomial(n-i, j, &binomial_j);
                temp_idx = (n+i-j)/2;
                // only compute Tad_master if it hasn't been already computed
                // check_precompute[idx]==1: already has been computed
                // check_precompute[idx]==0: has not been computed yet
                // TODO: this can be a class_Test
                if (temp_idx < 0 || temp_idx > n) {
                    printf("Error: temp_idx out of bounds (%d)\n", temp_idx);
                }
                if (check_precompute[temp_idx]==0){
                    Tad_master(temp_idx, d, M, &Tad_master_precompute[temp_idx]); 
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

 int massive_num(int n, int d, double k2, class_complex M1, class_complex M2, class_complex *I_out){
    int i;
    double binomial;
    class_complex Tad_var_temp, result = class_complex(0., 0.);
    for (i=0; i<=n; i++){
        util_binomial(n, i, &binomial);
        Tad_var(i, d, k2, M2, &Tad_var_temp);
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

int util_binomial(int n, int k, double *n_over_k){
    int i, j;
    int C[n + 1][k + 1];

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

int util_trinomial(int n, int k, int j, double *out){

    return _SUCCESS_;
}
