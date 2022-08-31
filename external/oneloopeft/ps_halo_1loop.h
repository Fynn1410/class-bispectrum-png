#ifndef LINE_PS_PT_H
#define LINE_PS_PT_H
#include "/home/dennis/Software/class/external/oneloopeft/library/Cuba-4.2.1/cuba.h"

int PS_hh_G(
               struct background * pba,
               struct primordial * ppm,
               struct fourier * pfo,
               double k,
               double z,
               short has_loop,
               short has_ir,
               long SPLIT, 
               double * pk_nl);

int PS_mm_G(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            double k,
            double z,
            short has_loop,
            short has_ir,
            long SPLIT,
            double * pk_nl);

int PS_hh_1(   struct background * pba,
               struct primordial * ppm,
               struct fourier * pfo,
               double k,
               double z,
               double mu,
               long SPLIT);

int PS_hh_2(   struct background * pba,
               struct primordial * ppm,
               struct fourier * pfo,
               double k,
               double z,
               double mu,
               long SPLIT);

int PS_hh_3(   struct background * pba,
               struct primordial * ppm,
               struct fourier * pfo,
               double k,
               double z,
               double mu,
               long SPLIT);

int PS_hh_4(   struct background * pba,
               struct primordial * ppm,
               struct fourier * pfo,
               double k,
               double z,
               double mu,
               long SPLIT);

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

int Compute_1_loops(
                    struct background * pba,
                    struct primordial * ppm,
                    struct fourier * pfo,
                    double k,
                    double z,
                    double mu,
                    long SPLIT,
                    double *result);
                    
static int G1_loop_integrands(
                             const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p);

int Compute_2_loops(
                    struct background * pba,
                    struct primordial * ppm,
                    struct fourier * pfo,
                    double k,
                    double z,
                    double mu,
                    long SPLIT,
                    double *result);
                    
static int G2_loop_integrands(
                             const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p);

int Compute_3_loops(
                    struct background * pba,
                    struct primordial * ppm,
                    struct fourier * pfo,
                    double k,
                    double z,
                    double mu,
                    long SPLIT,
                    double *result);
                    
static int G3_loop_integrands(
                             const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p);

int Compute_4_loops(
                    struct background * pba,
                    struct primordial * ppm,
                    struct fourier * pfo,
                    double k,
                    double z,
                    double mu,
                    long SPLIT,
                    double *result);
                    
static int G4_loop_integrands(
                             const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p);

double F2_s(double k1,double k2,double mu);
double G2_s(double q,double k,double mu);
double S2_s(double k1,double k2,double mu);
double F3_s(double k,double q, double mu);
double G3_s(double k,double q, double mu);
double S2(double mu);
double F2(double k1,double k2,double mu);
double G2(double k1,double k2,double mu);

double LoS1(double q, double k, double cos, double mu);
double LoS2(double q, double k, double cos, double mu);
double LoS3(double q, double k, double cos, double mu);
double LoS4(double q, double k, double cos, double mu);
double LoS5(double q, double k, double cos, double mu);
double LoS6(double q, double k, double cos, double mu);

#endif
