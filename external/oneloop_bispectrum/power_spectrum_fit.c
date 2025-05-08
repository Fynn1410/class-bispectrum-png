/** @file power_spectrum_fit.c
 *
 * author: Fynn Janssen, 2025
 * strongly based on https://arxiv.org/abs/2212.07421, https://github.com/dbraganca/python-integer-powers
*/

#include "power_spectrum_fit.h"


// ************************************************************************************************
// Global variables
// ************************************************************************************************

// parameters used in the power spectrum fit, following: https://arxiv.org/abs/2212.07421
static const double k2_UV[5] = {1e-4, 6.9e-2, 8.2e-3, 1.3e-3, 1.35e-5};         // in units [h/Mpc]^2
static const double k2_peak[5] = {0., -3.4e-2, -1e-3, -7.6e-5, -1.56e-5};       // in units [h/Mpc]^2


// Indices that are used to compute the power spectrum expansion in massive propagators
static const int D = 17;
static const int K_size = 5;

static const int index_UV[D] = {0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4};
static const int index_Kj[D] = {1, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 5, 1, 2, 3, 4};
static const int index_alpha[D] = {0, 0, 0, 0, 4, 4, 4, 4, 8, 8, 8, 8, 8, 12, 12, 12, 12};
static const int index_exp_k[D] = {1, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 4, 5, 1, 2, 3, 4};
static const int index_Ki[D] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1};
static const int index_loop_lower[D] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 0, 1, 2, 3, 0, 1, 2, 3};

// Indices that are used to compute the power spectrum expansion in fitting functions f_n
static const int N = 16; // number of terms in the fitting function of the power spectrum

static const int index_exp_num[N] = {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
static const int index_exp_den[N] = {1, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 1, 2, 3, 4};
static const int group_idx[N]   = {0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};


// ************************************************************************************************
// Compute the fitted linear matter power spectrum
// ************************************************************************************************

/** Pfit(k) = sum^D_d C_d * k^2i_d / (k^2 + M_d)^j_d
 * @param ppf           Input: pointer to ps_fit structure
 * @param pti           Input: pointer to gen_tri_integral structure
 * @param pba           Input: pointer to background structure
 * @param ppm           Input: pointer to primordial structure
 * @param pfo           Input: pointer to fourier structure
 * @param N_fit         Input: number of fitting points
 * @param k_min         Input: lower limit of fitting region
 * @param k_max         Input: upper limit of fitting region
 * @param z             Input: redshift z
 * @param C_out         Output: double pointer (matrix) which contains all the 
 *                            information for the linear matter power spectrum fit
 * @return the error status
*/

int pfit_massive_propagator(
                            struct ps_fit *ppf,
                            struct gen_tri_integral *pti,
                            struct background *pba,
                            struct primordial * ppm,
                            struct fourier *pfo,
                            int N_fit,
                            double kmin,
                            double kmax, //class_complex **out     // needs to be preallocated. Size (4, 17)
                            double z,
                            class_complex *C_out
                            ){
    // Important note:
    // this expansion has 33 terms: 
    //  - first term is real
    //  - terms 2-17 are complex
    //  - terms 18-33 are the complex conjugate of terms 2-17
    // -> only store information on the first 17 (= D) indices 
    class_complex C_d[D], term;
    double h, h2, *alpha;
    int idx_alpha, idx_exp, idx_UV, idx_Ki, idx_Kj, idx_ll, i, j;       // declare utility indices
    h = pba->h;
    h2 = h*h;
    

    // allocate memory
    class_alloc(alpha, N*sizeof(double),
                ppf->error_message);

    // compute the initial expansion coeffs alpha_n in P_fit(k^2) = sum_n alpha_n f_n(k^2)
    class_call(pfit_coeffs(ppf,
                           pba,
                           ppm,
                           pfo,
                           N_fit,
                           kmin,
                           kmax,
                           z,
                           alpha),
                ppf->error_message,
                ppf->error_message);

    // procompute the coeffs K[n][l]
    class_call(util_compute_K(ppf),
                pti->error_message,
                ppf->error_message);


    // special case
    C_out[0] = alpha[0] * k2_UV[0]*h2;

    // compute C_d coeffs for all other cases
    for (i=1; i<D; i++){
        term = class_complex(0., 0.);
        idx_UV = index_UV[i];
        idx_exp = index_exp_k[i];
        idx_ll = index_loop_lower[i];

        for (j=idx_ll; j<=3; j++){
            idx_alpha = index_alpha[i] + j;
            idx_Ki = j + index_Ki[i] - 1;
            idx_Kj = index_Kj[i] - 1;
            term += alpha[idx_alpha] * ppf->K_coeff[idx_Ki][idx_Kj] * pow(k2_UV[idx_UV]*h2, idx_exp);
        }
        C_out[i] = term;
    }

    free(alpha);
    
    return _SUCCESS_;
}


/** Pfit(k) = sum^16_n a_n * f_n(k^2)
 * @param ppf           Input: pointer to ps_fit structure
 * @param pba           Input: pointer to background structure
 * @param ppm           Input: pointer to primordial structure
 * @param pfo           Input: pointer to fourier structure
 * @param N_fit         Input: number of fitting points
 * @param kmin          Input: lower limit of fitting region
 * @param kmax          Input: upper limit of fitting region
 * @param z             Input: redshift z
 * @param alpha         Output: array of fittings coeffs
 * @return the error status
*/

int pfit_coeffs(
                struct ps_fit *ppf,
                struct background *pba,
                struct primordial * ppm,
                struct fourier *pfo,
                int N_fit,
                double kmin,
                double kmax,
                double z,
                double *alpha     // needs to be preallocated
                ){
    
    double *X, *W, *Pk, *A, *y, k, dk, h, h2; // X_ij = f_i(k2_j)
    int i, j, m;
    // reduces Hubble constant
    h = pba->h;
    h2 = h*h;

    // flattened 2D array -> X[i][j] -> X[i*N + j] (i goes from 0 to N_fit-1 and j goes from 0 to N-1)
    class_alloc(X, N*N_fit*sizeof(double),
                ppf->error_message);

    // data vector for the linear matter power spectrum
    class_alloc(Pk, N_fit*sizeof(double),
                ppf->error_message);

    // diagonal weight matrix stored as a 1D array
    class_alloc(W, N_fit*sizeof(double),
                ppf->error_message);

    // auxiliary matrix: A = X^T * W * X
    class_alloc(A, N*N*sizeof(double),
                ppf->error_message);

    // auxiliary vector: A = X^T * W * P_lin
    class_alloc(y, N*sizeof(double),
                ppf->error_message);



    // initialize the matrix A and the vector y:
    for (i=0; i<N; i++){
        y[i] = 0.;
        for (j=0; j<N; j++){
            A[i*N + j] = 0.;
        }
    }


    // define log k-grid
    dk = (log10(kmax) - log10(kmin))/(N_fit - 1);

    // fill the X matrix, the diagonal weight matrix W (store as 1D array), and the linear power spectrum array Pk
    for (i=0; i<N_fit; i++){
        k = pow(10., log10(kmin) + dk*i);

        // compute the fitting function
        class_call(util_f_n(ppf,
                            h2,
                            k*k,
                            &X[i*N]),
                    ppf->error_message, 
                    ppf->error_message);

        // compute P_lin
        class_call(fourier_pk_at_k_and_z(pba,
                                         ppm,
                                         pfo,
                                         pk_linear,
                                         k,
                                         z,
                                         pfo->index_pk_m,
                                         &Pk[i], 
                                         NULL),
                    pfo->error_message,
                    ppf->error_message);
        
        // fill the array for the diagnonal weight matrix
        W[i] = 1./(Pk[i] * Pk[i]);
    }

    
    // compute the fitting coeffs a = (X^T * W * X)^-1 * X^T * W * P_lin     ->      matrix-vector products
    // compute the auxialiry matrix A = X^T * W * X, which is symmetric
    for (i=0; i<N_fit; i++){
        for (j=0; j<N; j++){
            // y = X^T * W * P_lin, note: W*P_lin = 1/P_lin (element wise)
            y[j] += X[i*N + j] * W[i] * Pk[i];

            // A = X^T * W * X
            for (m=0; m<N; m++){
                A[j*N + m] += X[i*N + m] * W[i] * X[i*N + j];
            }
        }
    }

    class_call(util_matrix_vector_prod(ppf,
                                       A,
                                       y,
                                       N,
                                       alpha),
                ppf->error_message,
                ppf->error_message);

    // free the memory
    free(A);
    free(W);
    free(y);
    free(X);
    free(Pk);
    
    return _SUCCESS_;
}



// ************************************************************************************************
// Utility functions
// ************************************************************************************************


/** utility function f_n(k^2), which appears in Pfit(k) = sum^N_n a_n * f_n(k^2)
 * @param ppf           Input: pointer to ps_fit structure
 * @param h2            Input: reduced Hubble constant squared
 * @param k2            Input: wavenumber squared
 * @param out           Output: array of fitting functions f_n
 * @return the error status
*/

int util_f_n(
             struct ps_fit *ppf,
             double h2,
             double k2,
             double *out         // needs to be pre-allocated
             ){
    double ratio, deltak2, denominator, num_pow; // tmp auxiliary variables
    int i, j;
    

    for (j=0; j<N; j++){
        i = group_idx[j];
        ratio = k2/(k2_UV[i]*h2);
    
        if (j==0){
            out[j] = 1./(1. + ratio);
        } else {
            deltak2     = k2 - k2_peak[i]*h2;
            denominator = 1. + (deltak2 * deltak2) / (k2_UV[i] * k2_UV[i]* h2*h2);

            // To avoid calculating powers of zero
            if (index_exp_num[j]==0){
                num_pow = 1.;
            } else {
                num_pow = pow(ratio, index_exp_num[j]);
            }

            out[j] = num_pow / pow(denominator, index_exp_den[j]);
        }
    }

    return _SUCCESS_;
}



/** utility function to compute the invers of a NxN matrix
 * @param ppf           Input: pointer to ps_fit structure
 * @param A             Input: matrix A = X^T * W * X
 * @param y             Input: vector y = X^T * W * P_lin
 * @param N             Input: matrix dimension
 * @param alpha         Output: coeff vector alpha_n = A^-1 * b
 * @return the error status
*/

int util_matrix_vector_prod(
                            struct ps_fit *ppf,
                            double *A,       // input matrix (symmetric positive-definite)
                            double *b,       // right-hand side
                            int N,
                            double *alpha // output
                            ){
    int i, j, m;
    double **L, *y, sum;

    // memory allocation of the lower triangle matrix L, defined through A = L*L^T
    class_alloc(L, N*sizeof(double*), ppf->error_message);
    for (i=0; i<N; i++){
        class_alloc(L[i], (i+1)*sizeof(double), ppf->error_message);
    }

    // allocate memory for vector y, defined throught L*y = b 
    class_alloc(y, N*sizeof(double), ppf->error_message);


    // Cholesky decomposition: A = L * L^T
    for (i = 0; i < N; i++) {
        for (j = 0; j < i+1; j++) {
            sum = 0.;
            for (m = 0; m < j; m++) {
                sum += L[j][m] * L[i][m];
            }

            if (i == j) {
                // formula for diagonal elements
                L[i][i] = sqrt(A[i*N+i] - sum);
            } else {
                // formula for off-diagonal elements
                L[i][j] = (A[i*N + j] - sum) / L[j][j];
            }
            // TODO: class_test
            if (L[j][j]==0.){
                printf("\ndiv. by zero 1");
            }
        }
    }

    // Forward substitution: L * y = b
    for (i = 0; i < N; i++) {
        sum = 0.;
        for (j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }

        y[i] = (b[i] - sum) / L[i][i];
    }

    // Backward substitution: L^T * alpha = y
    for (i = N-1; i >= 0; i--) {
        sum = 0.;
        for (j = i+1; j < N; j++) {
            sum += L[j][i] * alpha[j];
        }

        alpha[i] = (y[i] - sum) / L[i][i];
    }


    // free the memory
    for (i=0; i<N; i++){
        free(L[i]);
    }
    free(L);
    free(y);
    return _SUCCESS_;
}


/** Calculate binomial coefficients n over k
 * @param n                     Input: integer
 * @param k                     Input: integer
 * @param K_out                 Output: pointer to output value
 * @return the error status
 */

 int util_compute_K(struct ps_fit *ppf) {
    class_complex K0;
    double sign, tmp_bin;
    int i, j;
    
    K0 = class_complex(0., 0.5);

    // Allocate top-level array
    class_alloc(ppf->K_coeff, K_size * sizeof(class_complex *),
                ppf->error_message);

    // Allocate each row (note: i+1 elements)
    for (i = 0; i < K_size; i++) {
        class_alloc(ppf->K_coeff[i], (i + 1) * sizeof(class_complex),
                    ppf->error_message);
    }

    // Fill the triangular K matrix
    // K = binomial(2*i - j - 1, i - j) * K0^i * conj(K0)^(i-j)
    // since K0 is purely imaginary simply to K = binomial(2*i - j - 1, i - j) * K0^(2i-j) * (-1)^(i-j))
    for (i = 1; i <= K_size; i++) {
        for (j = 1; j <= i; j++) {
            // get binomial coeff
            // TODO: this should be a class_call using util_binomial from the generalized_triangle_integral.c file
            class_call(util_binomial_tmp(ppf, 2*i - j - 1, i - j, &tmp_bin),
                       ppf->error_message,
                       ppf->error_message);
            
            sign = ((i - j) % 2 == 0) ? 1. : -1.;

            ppf->K_coeff[i - 1][j - 1] = tmp_bin * pow(K0, 2*i - j) * sign;
        }
    }

    return _SUCCESS_;
}



// TODO: remove this

/** Calculate binomial coefficients n over k
 * @param n                     Input: integer
 * @param k                     Input: integer
 * @param n_over_k              Output: pointer to output value
 * @return the error status
 */

int util_binomial_tmp(
                    struct ps_fit *ppf,
                    int n,
                    int k,
                    double *n_over_k
                    ){

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
