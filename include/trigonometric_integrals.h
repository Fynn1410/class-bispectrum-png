/**
 * definitions for module trigonometric_integrals.c
 */

#ifndef __TRIGONOMETRIC_INTEGRALS__
#define __TRIGONOMETRIC_INTEGRALS__

#include "common.h"
#include "arrays.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int cosine_integral(
          double x,
          double *Ci,
          ErrorMsg error_message
          );

  int sine_integral(
          double x,
          double *Si,
          ErrorMsg error_message
          );

  int spherical_bessel_j(
          const int n, 
          double x,
          double * result,
          ErrorMsg errmsg
          );

#ifdef __cplusplus
}
#endif

#endif
