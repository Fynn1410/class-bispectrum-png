#include "header.h"
#include "utilities.h"
#include "../../include/common.h"

#ifndef __EFT_KERNEL_MATRICES__
#define __EFT_KERNEL_MATRICES__

#define _PI_CUBED_ 31.006276680299820175476315067101395202225288565885
#define J_PREFACTOR 2.2448390265645820211135247953461594410273242929353e-2  /**< 1/(8*pi^3/2) */

const static double complex (*loop_mat[NUM_MOMENTS])(const double complex * const n);

struct eft;
int eft_compute_loop_matrices(struct eft * peft);

#endif
