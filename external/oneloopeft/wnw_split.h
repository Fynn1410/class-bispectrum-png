#include "header.h"

#ifndef __NOWIGGLE_FILTER__
#define __NOWIGGLE_FILTER__

static inline double eft_gfilter_smoothing_scale(const double ln_k) {
  return 0.6907755279 * exp( -pow((ln_k - (-3.4538776395)) / 2.752115578658, 2) ) + 0.06907755279;
  // return ppr->nowiggle_filter_amplitude * exp( -pow((ln_k - ppr->nowiggle_filter_ln_k_center) / ppr->nowiggle_filter_ln_k_width, 2) ) + ppr->nowiggle_filter_const;
}

int eft_ln_pk_nw_gfilter(struct precision *ppr, struct background *pba, struct primordial *ppm, struct fourier *pfo, const int index_pk, const int index_k0, const int index_kmin, const int k_size, double *ln_pknw_array);
int eft_ln_pk_nw_gfilter_parallel(struct precision *ppr, struct background *pba, struct primordial *ppm, struct fourier *pfo, const int index_pk, const int index_k0, const int index_kmin, const int k_size, double *ln_pknw_array);
int eft_ln_pk_nw_gfilter_3d(struct precision *ppr, struct background *pba, struct primordial *ppm, struct fourier *pfo, const int index_pk, const int index_k0, const int index_kmin, const int k_size, double *ln_pknw_array);
double eft_pk_nw_eisenstein_hu_factor(struct background *pba, struct primordial *ppm, struct fourier *pfo, const double k, const double k0);
double eft_pk_w_eisenstein_hu_factor(struct background *pba, struct primordial *ppm, struct fourier *pfo, const double k, const double k0);

double T0( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k);
double T( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k);
double Tt0( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double x1, double x2);

#endif



