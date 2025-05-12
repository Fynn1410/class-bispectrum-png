#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "../../include/carray.h"
#include "../../include/common.h"

#ifndef __LOOP_INTEGRALS__
#define __LOOP_INTEGRALS__

struct gen_tri_integral
{
    ErrorMsg error_message;
};

struct doublets_L
{
    int n;
    int d;
};


#ifdef __cplusplus
extern "C" {
#endif

    int T_master(struct gen_tri_integral *pti,
                 double k12,
                 double k22,
                 double k32,
                 class_complex M1,
                 class_complex M2,
                 class_complex M3,
                 class_complex *T_out);

    int B_master(struct gen_tri_integral *pti,
                 double k2,
                 class_complex M1,
                 class_complex M2,
                 class_complex *B_out);

    int Tad_master(struct gen_tri_integral *pti,
                   int n,
                   int d,
                   class_complex M,
                   class_complex *Tad_out);


    int L_recursion(struct gen_tri_integral *pti,
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
                    class_complex *L_out);

    int T_recursion(struct gen_tri_integral *pti,
                    int d1,
                    int d2,
                    int d3,
                    double k12,
                    double k22,
                    double k32,
                    class_complex M1,
                    class_complex M2,
                    class_complex M3,
                    class_complex *T_out);

    int B_recursion(struct gen_tri_integral *pti,
                    int d1,
                    int d2,
                    double k2,
                    class_complex M1,
                    class_complex M2,
                    class_complex *B_out);

    int scalar_prod_one(struct gen_tri_integral *pti,
                       int m, 
                       int n,
                       int d1,
                       int d2,
                       double k22,
                       class_complex M1,
                       class_complex M2,
                       class_complex *I_out);

    int tensor_red_one(struct gen_tri_integral *pti,
                       int n,
                       int d1,
                       int d2,
                       double k12,
                       double k22,
                       double cos12,
                       class_complex M1,
                       class_complex M2,
                       class_complex *I_out);

    int tensor_red_two(struct gen_tri_integral *pti,
                       int n1,
                       int n2,
                       int d,
                       double k12,
                       double k22,
                       double cos12,
                       class_complex M,
                       class_complex *I_out);

    int Tad_var(struct gen_tri_integral *pti,
                int n,
                int d,
                double k2,
                class_complex M,
                class_complex *I_out);

    int massive_num(struct gen_tri_integral *pti,
                    int n, 
                    int d,
                    double k2,
                    class_complex M1, 
                    class_complex M2, 
                    class_complex *I_out);

    int util_binomial(struct gen_tri_integral *pti,
                      int n,
                      int k,
                      double *n_over_k);

    int util_trinomial(struct gen_tri_integral *pti,
                       int n,
                       int k,
                       int j,
                       double *out);

    int util_antideriv(struct gen_tri_integral *pti,
                       double x,
                       class_complex y1,
                       class_complex y2,
                       class_complex x0,
                       class_complex *out);

    int util_prefactor(struct gen_tri_integral *pti,
                       double a,
                       class_complex y1,
                       class_complex y2,
                       class_complex *out);

    int util_F_int(struct gen_tri_integral *pti,
                   double R2,
                   class_complex y1,
                   class_complex y2,
                   class_complex x0,
                   class_complex *out);

    int util_Tmaster_contr(struct gen_tri_integral *pti,
                           int y,
                           double k12,
                           double k22,
                           double k32,
                           class_complex M1,
                           class_complex M2,
                           class_complex M3,
                           class_complex *out);

    int util_L_step(struct gen_tri_integral *pti,
                    doublets_L idxs[],
                    int idx,
                    double k12,
                    double k22,
                    double k32,
                    class_complex Ms[],
                    class_complex *L_out);

//     int util_tri_master(double y,
//                         double k12,
//                         double k22,
//                         double k32,
//                         class_complex M1,
//                         class_complex M2,
//                         class_complex M3,
//                         class_complex *out);

//    int util_Fint(double aa, 
//                  class_complex Y1, 
//                  class_complex Y2, 
//                  class_complex X0,
//                  class_complex *out);

#ifdef __cplusplus
}
#endif

#endif
