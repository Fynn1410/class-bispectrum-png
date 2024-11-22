#include "header.h"

#ifndef __NOWIGGLE_FILTER__
#define __NOWIGGLE_FILTER__

int eft_ln_pk_nw_gfilter(struct precision *ppr, struct background *pba, struct primordial *ppm, struct fourier *pfo, int index_pk, int index_k0, int index_kmin, int k_size, double *ln_pknw_array);
int eft_ln_pk_nw_gfilter_3d(struct precision *ppr, struct background *pba, struct primordial *ppm, struct fourier *pfo, int index_pk, int index_k0, int index_kmin, int k_size, double *ln_pknw_array);
double eft_pk_nw_eisenstein_hu_factor(struct background *pba, struct primordial *ppm, struct fourier *pfo, double k, double k0);
double eft_pk_w_eisenstein_hu_factor(struct background *pba, struct primordial *ppm, struct fourier *pfo, double k, double k0);

double T0( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k);
double T( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k);
double Tt0( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double x1, double x2);

#endif
