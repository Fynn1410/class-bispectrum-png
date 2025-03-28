// TODO: maybe remove
#include <stdio.h>
#include <math.h>

#include "../../include/carray.h"
#include "../../include/common.h"

#ifndef __LOOP_INTEGRALS__
#define __LOOP_INTEGRALS__

#ifdef __cplusplus
extern "C" {
#endif

    int T_master(double k12,
                 double k22,
                 double k32,
                 class_complex M1,
                 class_complex M2,
                 class_complex M3,
                 class_complex *T_out);

    int B_master(double k2,
                 class_complex M1,
                 class_complex M2,
                 class_complex *B_out);

    int Tad_master(int n,
                   int d,
                   class_complex M,
                   class_complex *Tad_out);


    int L_recursion(int n1,
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

    int T_recursion(int d1,
                    int d2,
                    int d3,
                    double k12,
                    double k22,
                    double k32,
                    class_complex M1,
                    class_complex M2,
                    class_complex M3,
                    class_complex *T_out);

    int B_recursion(int d1,
                    int d2,
                    double k2,
                    class_complex M1,
                    class_complex M2,
                    class_complex *B_out);

    int tensor_red_one(int n,
                       int d1,
                       int d2,
                       double k12,
                       double k22,
                       class_complex M1,
                       class_complex M2,
                       class_complex *I_out);

    int tensor_red_two(int n1,
                       int n2,
                       int d1,
                       double k12,
                       double k22,
                       class_complex M,
                       class_complex *I_out);

    int Tad_var(int n,
                int d,
                double k2,
                class_complex M,
                class_complex *I_out);

    int massive_num(int n, 
                    int d,
                    double k2,
                    class_complex M1, 
                    class_complex M2, 
                    class_complex *I_out);

    int util_binomial(int n,
                      int k,
                      double *n_over_k);

    int util_trinomial(int n,
                       int k,
                       int j,
                       double *out);

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
