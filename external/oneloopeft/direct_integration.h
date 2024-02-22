/** #if DIRECT_INTEGRATION : make will only compile this file if DIRECT_INTEGRATION=yes is given */

#include "header.h"
#include "../../include/common.h"

#ifndef __EFT_DIRECT_INTEGRATION__
#define __EFT_DIRECT_INTEGRATION__

#include <stdlib.h>
#include "library/Cuba-4.2.2/cuba.h"

#undef _FAILURE_
#define _FAILURE_ -999

struct direct_integration_parameters
{
  struct eft * peft;
  struct background * pba;
  struct primordial * ppm;
  struct fourier * pfo;

  struct eft_hyper_parameters * eft_hp;

  const int * const moment_list;
  const int moment_list_size;

  const double z;
  const double f_z;
  double ln_k;  /**< external logarithmic wavenumber (in 1/Mpc) */
  double mu;    /**< angle w.r.t. line-of-sight mu = e_k . n */
};


int eft_di_compute_spectra_contributions(
                    struct eft * peft,
                    struct background * pba,
                    struct fourier * pfo,
                    struct primordial * ppm,
                    struct precision * ppr,
                    const int * const moment_list,
                    const int moment_list_size,
                    const double z,
                    const double f_z
                    );


#endif