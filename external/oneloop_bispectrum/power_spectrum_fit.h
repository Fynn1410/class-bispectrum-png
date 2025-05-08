#include <stdio.h>
#include <math.h>

#include "../../include/carray.h"
#include "../../include/common.h"
#include "../../include/fourier.h"
#include "../../include/background.h"
#include "generalized_triangle_integral.h"


#ifndef __PS_FIT__
#define __PS_FIT__

struct ps_fit
{
    ErrorMsg error_message;
    class_complex **K_coeff;
};


#ifdef __cplusplus
extern "C" {
#endif

    int pfit_massive_propagator(struct ps_fit *ppf,
                                struct gen_tri_integral *pti,
                                struct background *pba,
                                struct primordial * ppm,
                                struct fourier *pfo,
                                int N_fit,
                                double k_min,
                                double k_max,
                                double z,
                                class_complex *C_out);

    int pfit_coeffs(struct ps_fit *ppf,
                    struct background *pba,
                    struct primordial * ppm,
                    struct fourier *pfo,
                    int N_fit,
                    double kmin,
                    double kmax,
                    double z,
                    double *alpha);

    int util_f_n(struct ps_fit *ppf,
                 double h2,
                 double k2,
                 double *out);

    int util_matrix_vector_prod(struct ps_fit *ppf,
                                double *A,       // input matrix (symmetric positive-definite)
                                double *b,       // right-hand side
                                int N,
                                double *alpha);
    
    int util_compute_K(struct ps_fit *ppf);

    int util_binomial_tmp(struct ps_fit *ppf,
                          int n,
                          int k,
                          double *n_over_k);

#ifdef __cplusplus
}
#endif

#endif

