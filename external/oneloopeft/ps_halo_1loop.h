#ifndef LINE_PS_PT_H
#define LINE_PS_PT_H
#include "cuba.h"

double PS_hh_G(
               //struct Cosmology *Cx,
               struct precision * ppr,
               struct background * pba,
               struct perturbations * ppt,
               struct primordial * ppm,
               struct fourier * pfo,
               int index_k,
               double z,
               double M,
               short has_loop,
               short has_ir,
               long SPLIT,
               long mode_mf);

int PS_mm_G(
            struct precision * ppr,
            struct background * pba,
            struct perturbations * ppt,
            struct primordial * ppm,
            struct fourier * pfo,
            double k,
            double z,
            short has_loop,
            short has_ir,
            long SPLIT,
            double * pk_nl);

int Compute_G_loops(
                    //struct Cosmology *Cx,
                    struct background * pba,
                    struct primordial * ppm,
                    struct fourier * pfo,
                    double k,
                    double z,
                    short has_ir,
                    long hm_switch,
                    long SPLIT,
                    double *result);


static int G_loop_integrands(
                             const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p);

double F2_s(double k1,double k2,double mu);
double S2_s(double k1,double k2,double mu);
double F3_s(double k,double q, double mu);
double S2(double mu);
double F2(double k1,double k2,double mu);


#endif
