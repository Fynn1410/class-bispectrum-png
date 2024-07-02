/** @file direct_integration.c
 *
 * author: Christian Radermacher, 2024
 * based on prototype version of Azadeh Moradinezhad Dizgah & Dennis Linde
 *
 * contains kernel definitions for direct integration of loop terms
 * only compiled if DIRECT_INTEGRATION=yes is given
*/

#include "direct_integration.h"

/** All kernels can be expressed in terms of k^2, q^2 and |k-q|^2 */

/** F2_s(q, k-q) */
static double eft_di_kernel_F2s_q_kmq(const double k_sq, const double q_sq, const double kmq_sq) {
  return (2.*k_sq*k_sq + 3.*k_sq*(kmq_sq + q_sq) - 5.*(kmq_sq - q_sq)*(kmq_sq - q_sq)) / (28.*kmq_sq*q_sq);
}

/** F2_s(q, -k) */
static double eft_di_kernel_F2s_q_mk(const double k_sq, const double q_sq, const double kmq_sq) {
  return (-5.*k_sq*k_sq + k_sq*(3.*kmq_sq + 10.*q_sq) + (2.*kmq_sq*kmq_sq + 3.*kmq_sq*q_sq - 5.*q_sq*q_sq)) / (28.*k_sq*q_sq);
}

/** G2_s(q, k-q) */
static double eft_di_kernel_G2s_q_kmq(const double k_sq, const double q_sq, const double kmq_sq) {
  return (4.*k_sq*k_sq - k_sq*(kmq_sq + q_sq) - 3.*(kmq_sq - q_sq)*(kmq_sq - q_sq)) / (28.*kmq_sq*q_sq);
}

/** G2_s(q, -k) */
static double eft_di_kernel_G2s_q_mk(const double k_sq, const double q_sq, const double kmq_sq) {
  return -(3.*k_sq*k_sq + k_sq*(kmq_sq - 6.*q_sq) + (-4.*kmq_sq*kmq_sq + kmq_sq*q_sq + 3.*q_sq*q_sq)) / (28.*k_sq*q_sq);
}

/** F3_s(k, q, -q) */
static double eft_di_kernel_F3s(const double k_sq, const double q_sq, const double kmq_sq) {
  return -(6.*k_sq*k_sq*k_sq*k_sq                                                                                 \
          +k_sq*k_sq*k_sq*(31.*kmq_sq + 3.*q_sq)                                                                  \
          -k_sq*k_sq*(66.*kmq_sq*kmq_sq - 5.*kmq_sq*q_sq + 45.*q_sq*q_sq)                                         \
          +k_sq*(15.*kmq_sq*kmq_sq*kmq_sq + 13.*kmq_sq*kmq_sq*q_sq - 85.*kmq_sq*q_sq*q_sq + 57.*q_sq*q_sq*q_sq)   \
          +7.*(kmq_sq - q_sq)*(kmq_sq - q_sq)*(kmq_sq - q_sq)*(2.*kmq_sq + 3.*q_sq)) / (1512.*k_sq*kmq_sq*q_sq*q_sq);
}

/** G3_s(k, q, -q) */
static double eft_di_kernel_G3s(const double k_sq, const double q_sq, const double kmq_sq) {
  return -(6.*k_sq*k_sq*k_sq*k_sq                                                                               \
          +k_sq*k_sq*k_sq*(kmq_sq - 15.*q_sq)                                                                   \
          -k_sq*k_sq*(18.*kmq_sq*kmq_sq - 11.*kmq_sq*q_sq - 9.*q_sq*q_sq)                                       \
          +k_sq*(9.*kmq_sq*kmq_sq*kmq_sq + 7.*kmq_sq*kmq_sq*q_sq - 19.*kmq_sq*q_sq*q_sq + 3.*q_sq*q_sq*q_sq)    \
          +(kmq_sq - q_sq)*(kmq_sq - q_sq)*(kmq_sq - q_sq)*(2.*kmq_sq + 3.*q_sq)) / (504.*k_sq*kmq_sq*q_sq*q_sq);
}

/** S2(q, k-q) = (q.(k-q))^2/(|q|^2*|k-q|^2) - 1 = (k mu_q - q)^2/|k-q|^2 - 1 */
static double eft_di_kernel_S2(const double k_sq, const double q_sq, const double kmq_sq) {
  double k, q, kmq;
  k = sqrt(k_sq); q = sqrt(q_sq); kmq = sqrt(kmq_sq);
  return (k + kmq + q)*(k + kmq - q)*(-k + kmq + q)*(-k + kmq - q) / (4.*kmq_sq*q_sq);
}

/** (k-q)_par / |k-q|^2 */
static double eft_di_los_kernel_kmq(const double k_sq, const double q_sq, const double kmq_sq, const double cos_thetakmq) {
  double kmq;
  kmq = sqrt(kmq_sq);
  return cos_thetakmq / kmq;
}

/** (q)_par / |q|^2 */
static double eft_di_los_kernel_q(const double k_sq, const double q_sq, const double kmq_sq, const double cos_thetaq, const double cos_phiq, const double cos_theta) {
  double q;
  q = sqrt(q_sq);
  return (cos_theta*cos_thetaq - sqrt(1. - cos_theta*cos_theta)*sqrt(1. - cos_thetaq*cos_thetaq)*cos_phiq) / q;
}


/**
 * @brief
 * @param ndim      Input: number of variables to integrate over
 * @param x         Input: array of arguments, indexed as x[index_vec][index_dim]
 * @param ncomp     Input: number of components of the integrand
 * @param f         Output: integrand values at x, indexed as f[index_vec][index_comp]
 * @param userdata  Input: additional constant parameters
 * @param nvec      Input: number of points to evaluate at the same time
 * @param core      Input: id of the core (< 0 for accelerators and > 0 for cores)
 *
 * @return the error status (_CUBA_ERROR_ for immediate abortion)
 */
int eft_di_integrands(const int * ndim,
                      const double * x,
                      const int * ncomp,
                      double * f,
                      void * userdata,
                      const int * nvec,
                      const int * core) {

  int it, index_list, index_moment, index_comp = 0;
  struct direct_integration_parameters * params = userdata;
  struct eft_hyper_parameters * hp = (struct eft_hyper_parameters *)params->eft_hp;
  const double k = exp(params->ln_k), cos_theta = params->mu, sin_theta = sqrt(1. - params->mu * params->mu), k_sq = exp(2.*params->ln_k);
  double * ln_q, * q_sq, * cos_thetaq, * ln_kmq, *kmq_sq, * cos_thetakmq, * cos_phiq, * measure;
  double * plin_rsd_q, * plin_rsd_kmq, q, kmq, sin_thetaq;
  short * remove_point;
  struct indexed_rsd_arg * argqvec, * argkmqvec;
  ErrorMsg errmsg;

  class_alloc(ln_q,         *nvec*sizeof(double), errmsg);
  class_alloc(q_sq,         *nvec*sizeof(double), errmsg);
  class_alloc(cos_thetaq,   *nvec*sizeof(double), errmsg);
  class_alloc(ln_kmq,       *nvec*sizeof(double), errmsg);
  class_alloc(kmq_sq,       *nvec*sizeof(double), errmsg);
  class_alloc(cos_thetakmq, *nvec*sizeof(double), errmsg);
  class_alloc(cos_phiq,     *nvec*sizeof(double), errmsg);
  class_alloc(measure,      *nvec*sizeof(double), errmsg);
  class_alloc(plin_rsd_q,   *nvec*sizeof(double), errmsg);
  class_alloc(plin_rsd_kmq, *nvec*sizeof(double), errmsg);
  class_calloc(remove_point, *nvec, sizeof(short), errmsg);

  /** - extract the arguments into (ln(q), mu_q, phi_q) from x */
  for (it = 0; it < *nvec; it++) {
    ln_q[it] = log(hp->k_IR_cutoff) + x[it*(*ndim) + 0] * log(hp->k_UV_cutoff/hp->k_IR_cutoff);
    q = exp(ln_q[it]);
    q_sq[it] = q*q;
    cos_thetaq[it]  = 2.*x[it*(*ndim) + 1] - 1.; /** angle w.r.t. k, not line-of-sight!; mu_q in (-1, 1), alternatively (0, 1) with additional symmetry factor */
    sin_thetaq = sqrt(1. - cos_thetaq[it]*cos_thetaq[it]);
    cos_phiq[it] = cos(2.*_PI_*x[it*(*ndim) + 2]);
    ln_kmq[it] = 0.5 * log( k_sq + q*q - 2.*k*q*cos_thetaq[it] ); /** values outside the stored interpolation range for the power spectrum will be removed from the integral (essentially pole cutoff = IR cutoff) */
    kmq = exp(ln_kmq[it]);
    kmq_sq[it] = kmq*kmq;
    cos_thetakmq[it] = ((k - q*cos_thetaq[it]) * cos_theta + q*sin_thetaq*cos_phiq[it] * sin_theta) / kmq;
    measure[it] = log(hp->k_UV_cutoff/hp->k_IR_cutoff) * q*q*q / (2.*_PI_*_PI_);  /**< integral measure after transforming to x_i in [0,1], starting from d^3q/(2 pi)^3 */

    /** - screen for points close to the pole q=k or the boundaries */
    if ((kmq < hp->k_pole_cutoff) || (kmq > hp->k_UV_cutoff)) {
      ln_kmq[it] = -4.; /** just an arbitrary default value inside the interpolation range */
      kmq_sq[it] = 3.35462628e-4;
      cos_thetakmq[it] = 0.;
      remove_point[it] = 1;
    }
  }

  /** - create the q argument list */
  class_call(eft_rsd_argument_list(ln_q,
                                   cos_thetaq,
                                   *nvec,
                                   &argqvec,
                                   errmsg),
             errmsg, errmsg);

  /** - create the k-q argument list */
  class_call(eft_rsd_argument_list(ln_kmq,
                                   cos_thetakmq,
                                   *nvec,
                                   &argkmqvec,
                                   errmsg),
             errmsg, errmsg);

  /** - get the (IR-resummed) linear power spectrum */
  class_call(eft_ir_pk_rsd_lo(params->pba,
                              params->ppm,
                              params->pfo,
                              linear,
                              argqvec,
                              *nvec,
                              params->z,
                              params->f_z,
                              hp->linear_spectrum_index,
                              params->peft->Sigma2_ir,
                              params->peft->dSigma2_ir,
                              plin_rsd_q),
             errmsg, errmsg);

  class_call(eft_ir_pk_rsd_lo(params->pba,
                              params->ppm,
                              params->pfo,
                              linear,
                              argkmqvec,
                              *nvec,
                              params->z,
                              params->f_z,
                              hp->linear_spectrum_index,
                              params->peft->Sigma2_ir,
                              params->peft->dSigma2_ir,
                              plin_rsd_kmq),
             errmsg, errmsg);

  free(argqvec); free(argkmqvec);

  for (it = 0; it < *nvec; it++) {
    if (remove_point[it])
      plin_rsd_kmq[it] = 0.;
  }

  for (index_list = 0; index_list < params->moment_list_size; index_list++) {
    index_moment = (params->moment_list[index_list] % params->peft->index_num);

    /** --------------------- Real-space moments ------------------------------ */
    if (index_moment == params->peft->index_I2200) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * pow(eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]), 2);
    }
    else if (index_moment == params->peft->index_I1300) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * eft_di_kernel_F3s(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_Idelta200) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_IG200) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]) * eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_Idelta2delta200) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * (plin_rsd_kmq[it] - plin_rsd_q[it]);
    }
    else if (index_moment == params->peft->index_IG2G200) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * pow(eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]), 2);
    }
    else if (index_moment == params->peft->index_Idelta2G200) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_FG200) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]) * eft_di_kernel_F2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    /** --------------------- 1-st RSD moments -------------------------------- */
    else if (index_moment == params->peft->index_I2201) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]) * eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_I1301p3101) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * (eft_di_kernel_F3s(k_sq, q_sq[it], kmq_sq[it]) + eft_di_kernel_G3s(k_sq, q_sq[it], kmq_sq[it]));
    }
    else if (index_moment == params->peft->index_Idelta201) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_IG201) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]) * eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_FG201) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]) * eft_di_kernel_F2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_J12101) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_G2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_J11201) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_kernel_F2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_J21101) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_Jdelta201) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta);
    }
    else if (index_moment == params->peft->index_JG201) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]);
    }
    /** --------------------- 2-nd RSD moments -------------------------------- */
    else if (index_moment == params->peft->index_J12102x) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_G2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_J12102y) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = 0.;
    }
    else if (index_moment == params->peft->index_J21102x) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_J21102y) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = 0.;
    }
    else if (index_moment == params->peft->index_Jdelta202x) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]);
    }
    else if (index_moment == params->peft->index_Jdelta202y) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = 0.;
    }
    else if (index_moment == params->peft->index_JG202x) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_JG202y) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = 0.;
    }
    else if (index_moment == params->peft->index_I2211) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * pow(eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]), 2);
    }
    else if (index_moment == params->peft->index_I1311) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * eft_di_kernel_G3s(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_J12111) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_G2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_J11211) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_kernel_F2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_J21111) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_N11x) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * (eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) + eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]));
    }
    else if (index_moment == params->peft->index_N11y) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = 0.;
    }
    /** --------------------- 3-rd RSD moments -------------------------------- */
    else if (index_moment == params->peft->index_J21112x) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_J21112y) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = 0.;
    }
    else if (index_moment == params->peft->index_J12112x) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_G2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == params->peft->index_J12112y) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = 0.;
    }
    else if (index_moment == params->peft->index_N12x) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * pow(eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta), 2) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]);
    }
    else if (index_moment == params->peft->index_N12y) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = 0.;
    }
    /** --------------------- 4-th RSD moments -------------------------------- */
    else if (index_moment == params->peft->index_N22x) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * pow(eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]), 2);
    }
    else if (index_moment == params->peft->index_N22y) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = 0.;
    }
    else if (index_moment == params->peft->index_N22z) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = 0.;
    }
    else if (index_moment == params->peft->index_sigmav_mu) {
      for (it = 0; it < *nvec; it++)
        f[it*(*ncomp) + index_comp] = log(hp->k_UV_cutoff/hp->k_IR_cutoff) * exp(ln_q[it]) / (2.*_PI_*_PI_) * pow(cos_thetaq[it], 2) * plin_rsd_q[it];
    }
    else {
      /** - abort immediately */
      return _FAILURE_;
    }

    index_comp++;
  }

  free(ln_q); free(q_sq); free(cos_thetaq);
  free(ln_kmq); free(kmq_sq); free(cos_thetakmq); free(cos_phiq);
  free(measure); free(plin_rsd_q); free(plin_rsd_kmq);
  free(remove_point);

  return _SUCCESS_;
}


int eft_di_compute_spectra_contributions(struct eft * peft,
                                         struct background * pba,
                                         struct fourier * pfo,
                                         struct primordial * ppm,
                                         struct precision * ppr,
                                         const int * const moment_list,
                                         const int moment_list_size,
                                         const double z,
                                         const double f,
                                         const double D) {

  int index_k, index_mu, index_list, index_pk_type, index_moment, nregions, neval, fail;
  int ncores = 0, npoints = ppr->eft_di_vecsize, naccel = 0, npointsaccel = 0, max_error_index;
  double * result, * error, * prob, max_rel_error;
  div_t list_elem;
  struct direct_integration_parameters di_params = { .peft = peft, .eft_hp = peft->hp,    \
                                                     .pba = pba, .pfo = pfo, .ppm = ppm,  \
                                                     .moment_list = moment_list, .moment_list_size = moment_list_size,  \
                                                     .z = z, .f_z = f, .ln_k = 0., .mu = 0.};

  peft->z0 = z;
  peft->f_z0 = f; /**< growth rate may be supplied externally */
  peft->D_z0 = D;

  /** Compute the suppression factor for IR-Resummation */
  peft->Sigma2_ir  =  eft_ir_sigma2(pba, ppm, pfo, peft->z0, peft->hp->ir_resummation_k_split*pba->h, peft->hp->ir_resummation_k_feature*pba->h);  /** k_split = 0.2 h/Mpc, k_feature = 1/110 h/Mpc */
  peft->dSigma2_ir = eft_ir_dsigma2(pba, ppm, pfo, peft->z0, peft->hp->ir_resummation_k_split*pba->h, peft->hp->ir_resummation_k_feature*pba->h);  /** according to arXiV:1804.05080 by Ivanov & Sibiryakov */

  cubacores(&ncores, &npoints);
  cubaaccel(&naccel, &npointsaccel);

  printf("Starting direct integration.\n");

  class_alloc_parallel(result, moment_list_size*sizeof(double), peft->error_message);
  class_alloc_parallel(error,  moment_list_size*sizeof(double), peft->error_message);
  class_alloc_parallel(prob,   moment_list_size*sizeof(double), peft->error_message);

  for (index_mu = 0; index_mu < peft->mu_size; index_mu++) {
    for (index_k = 0; index_k < peft->k_size; index_k++) {
      /** - load the sampling point */
      di_params.ln_k = peft->ln_k[index_mu*peft->k_size + index_k];
      di_params.mu   = peft->mu[index_mu];

      /** - compute the integral using cubature rules */
      Cuhre(3,
            moment_list_size,
            eft_di_integrands,
            &di_params,
            ppr->eft_di_vecsize,
            ppr->eft_di_epsrel,
            ppr->eft_di_epsabs,
            ppr->eft_di_flags,
            ppr->eft_di_mineval,
            ppr->eft_di_maxeval,
            ppr->eft_di_key,
            NULL,
            NULL,
            &nregions,
            &neval,
            &fail,
            result,
            error,
            prob);

      // Suave(3,
      //       moment_list_size,
      //       eft_di_integrands,
      //       &di_params,
      //       ppr->eft_di_vecsize,
      //       ppr->eft_di_epsrel,
      //       ppr->eft_di_epsabs,
      //       ppr->eft_di_flags,
      //       0,
      //       ppr->eft_di_mineval,
      //       ppr->eft_di_maxeval,
      //       )

      /** - copy the result to the spectra_contributions array */
      for (index_list = 0; index_list < moment_list_size; index_list++) {
        list_elem = div(moment_list[index_list], peft->index_num);
        index_pk_type = list_elem.quot;
        index_moment = list_elem.rem;
        peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + finite_part][index_mu*peft->k_size + index_k]     = result[index_list];
        peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + uv_divergence][index_mu*peft->k_size + index_k]   = 0.;
        peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + ir_divergence][index_mu*peft->k_size + index_k]   = 0.;
        peft->spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + pole_divergence][index_mu*peft->k_size + index_k] = 0.;
      }

      if (peft->hp->eft_verbose > 1) {
        if (fail == 0) {
          printf("Direct integration (mu = %6.3f, ln(k) = %6.3f): finished in %.2e evaluations \n", di_params.mu, di_params.ln_k, (double)neval);
        }
        else if (fail > 0) {
          max_rel_error = 0.;
          for (index_list = 0; index_list < moment_list_size; index_list++) {
            if ((result[index_list] != 0.) && (fabs(error[index_list]/result[index_list]) > max_rel_error)) {
              max_rel_error = fabs(error[index_list]/result[index_list]);
              max_error_index = index_list;
            }
          }
          index_moment = (moment_list[max_error_index] % peft->index_num);
          printf("Direct integration (mu = %6.3f, ln(k) = %6.3f): accuracy goal not met in %.2e evaluations; highest rel. error is %.2e for index_moment = %d \n", di_params.mu, di_params.ln_k, (double)neval, max_rel_error, index_moment);
        }
      }
    }
  }

  free(result);
  free(error);
  free(prob);

  return _SUCCESS_;
}
