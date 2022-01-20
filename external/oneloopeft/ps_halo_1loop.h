#ifndef LINE_PS_PT_H
#define LINE_PS_PT_H

double PS_hh_G(
               //struct Cosmology *Cx,
               struct precision * ppr,
               struct background * pba,
               struct perturbations * ppt,
               struct primordial * ppm,
               struct fourier * pfo,
               double k, double z, double M, long mode, long IR_switch, long SPLIT, long mode_mf);

void Compute_G_loops(
                     //struct Cosmology *Cx,
                     struct precision * ppr,
                     struct background * pba,
                     struct perturbations * ppt,
                     struct fourier * pfo,
                     double k, double z, long IR_switch, long hm_switch, long SPLIT, double *result);


static int G_loop_integrands(
                             struct primordial * ppm,
                             struct background * pba,
                             struct fourier * pfo,
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
