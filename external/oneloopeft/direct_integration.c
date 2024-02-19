/** #if DIRECT_INTEGRATION : make will only compile this file if DIRECT_INTEGRATION=yes is given */

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
int integrands(const int * ndim,
               const double x[],
               const int * ncomp,
               double f[],
               struct direct_integration_parameters * userdata,
               const int * nvec,
               const int * core) {

  int it, index_list, index_comp = 0;
  struct eft_hyper_parameters * hp = (struct eft_hyper_parameters *)userdata->eft_hp;
  const double k = exp(userdata->ln_k), cos_theta = userdata->mu, sin_theta = sqrt(1. - userdata->mu * userdata->mu), k_sq;
  double * ln_q, * q_sq, * cos_thetaq, * ln_kmq, *kmq_sq, * cos_thetakmq, * measure;
  double * plin_rsd_q, * plin_rsd_kmq, q, kmq, sin_thetaq, cos_phiq;
  struct indexed_rsd_arg ** argqvec, ** argkmqvec;
  ErrorMsg errmsg;
  k_sq = k*k;

  class_alloc(ln_q,         *nvec*sizeof(double), errmsg);
  class_alloc(q_sq,         *nvec*sizeof(double), errmsg);
  class_alloc(cos_thetaq,   *nvec*sizeof(double), errmsg);
  class_alloc(ln_kmq,       *nvec*sizeof(double), errmsg);
  class_alloc(kmq_sq,       *nvec*sizeof(double), errmsg);
  class_alloc(cos_thetakmq, *nvec*sizeof(double), errmsg);
  class_alloc(measure,      *nvec*sizeof(double), errmsg);
  class_alloc(plin_rsd_q,   *nvec*sizeof(double), errmsg);
  class_alloc(plin_rsd_kmq, *nvec*sizeof(double), errmsg);

  /** - extract the arguments into (ln(q), mu_q, phi_q) from x */
  for (it = 0; it < *nvec; it++) {
    ln_q[it] = log(hp->k_IR_cutoff) + x[it][0] * log(hp->k_UV_cutoff/hp->k_IR_cutoff);
    q = exp(ln_q[it]);
    q_sq[it] = q*q;
    cos_thetaq[it]  = 2.*x[it][1] - 1.; /** angle w.r.t. k, not line-of-sight!; mu_q in (-1, 1), alternatively (0, 1) with additional symmetry factor */
    sin_thetaq = sqrt(1. - cos_thetaq[it]*cos_thetaq[it]);
    cos_phiq = cos(2.*_PI_*x[it][2]);
    ln_kmq[it] = 0.5 * log( k_sq + q*q - 2.*k*q*cos_thetaq[it] ); /** values outside the stored interpolation range for the power spectrum will be removed from the integral (essentially pole cutoff = IR cutoff) */
    kmq = exp(ln_kmq[it]);
    kmq_sq[it] = kmq*kmq;
    cos_thetakmq[it] = ((k - q*cos_thetaq[it]) * cos_theta[it] + q*sin_thetaq[it]*cos_phiq[it] * sin_theta[it]) / kmq;
    measure[it] = log(hp->k_UV_cutoff/hp->k_IR_cutoff) * q*q*q / (2.*_PI_*_PI_);  /**< integral measure after transforming to x_i in [0,1], starting from d^3q/(2 pi)^3 */
  }

  /** - create the q argument list */
  class_call(eft_rsd_argument_list(ln_q,
                                   cos_thetaq,
                                   *nvec,
                                   argqvec,
                                   errmsg),
             errmsg, errmsg);

  /** - create the k-q argument list */
  class_call(eft_rsd_argument_list(ln_kmq,
                                   cos_thetakmq,
                                   *nvec,
                                   argkmqvec,
                                   errmsg),
             errmsg, errmsg);
  
  /** - get the (IR-resummed) linear power spectrum */
  class_call(eft_ir_pk_rsd_lo(userdata->pba,
                              userdata->ppm,
                              userdata->pfo,
                              linear,
                              argqvec,
                              *nvec,
                              userdata->z,
                              userdata->f_z,
                              hp->linear_spectrum_index,
                              userdata->peft->Sigma2_ir,
                              userdata->peft->dSigma2_ir,
                              plin_rsd_q),
             errmsg, errmsg);

  class_call(eft_ir_pk_rsd_lo(userdata->pba,
                              userdata->ppm,
                              userdata->pfo,
                              linear,
                              argkmqvec,
                              *nvec,
                              userdata->z,
                              userdata->f_z,
                              hp->linear_spectrum_index,
                              userdata->peft->Sigma2_ir,
                              userdata->peft->dSigma2_ir,
                              plin_rsd_kmq),
             errmsg, errmsg);

  for (index_list = 0; index_list < userdata->moment_list_size; index_list++) {
    index_moment = (userdata->moment_list[index_list] % peft->index_num);

    /** --------------------- Real-space moments ------------------------------ */
    if (index_moment == userdata->peft->index_I2200) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * pow(eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]), 2);
    }
    else if (index_moment == userdata->peft->index_I1300) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * eft_di_kernel_F3s(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_Idelta200) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_IG200) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]) * eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_Idelta2delta200) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * (plin_rsd_kmq[it] - plin_rsd_q[it]);
    }
    else if (index_moment == userdata->peft->index_IG2G200) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * pow(eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]), 2);
    }
    else if (index_moment == userdata->peft->index_Idelta2G200) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_FG200) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]) * eft_di_kernel_F2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    /** --------------------- 1-st RSD moments -------------------------------- */
    else if (index_moment == userdata->peft->index_I2201) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]) * eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_I1301p3101) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * (eft_di_kernel_F3s(k_sq, q_sq[it], kmq_sq[it]) + eft_di_kernel_G3s(k_sq, q_sq[it], kmq_sq[it]));
    }
    else if (index_moment == userdata->peft->index_Idelta201) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_IG201) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]) * eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_FG201) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]) * eft_di_kernel_F2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_J12101) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_G2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_J11201) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_kernel_F2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_J21101) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_Jdelta201) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta);
    }
    else if (index_moment == userdata->peft->index_JG201) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]);
    }
    /** --------------------- 2-nd RSD moments -------------------------------- */
    else if (index_moment == userdata->peft->index_J12102x) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_G2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_J12102y) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = 0.;
    }
    else if (index_moment == userdata->peft->index_J21102x) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_F2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_J21102y) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = 0.;
    }
    else if (index_moment == userdata->peft->index_Jdelta202x) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]);
    }
    else if (index_moment == userdata->peft->index_Jdelta202y) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = 0.;
    }
    else if (index_moment == userdata->peft->index_JG202x) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_S2(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_JG202y) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = 0.;
    }
    else if (index_moment == userdata->peft->index_I2211) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * pow(eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]), 2);
    }
    else if (index_moment == userdata->peft->index_I1311) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * eft_di_kernel_G3s(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_J12111) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_G2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_J11211) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_kernel_F2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_J21111) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_N11x) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * (eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) + eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]));
    }
    else if (index_moment == userdata->peft->index_N11y) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = 0.;
    }
    /** --------------------- 3-rd RSD moments -------------------------------- */
    else if (index_moment == userdata->peft->index_J21112x) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_G2s_q_kmq(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_J21112y) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = 0.;
    }
    else if (index_moment == userdata->peft->index_J12112x) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]) * eft_di_kernel_G2s_q_mk(k_sq, q_sq[it], kmq_sq[it]);
    }
    else if (index_moment == userdata->peft->index_J12112y) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = 0.;
    }
    else if (index_moment == userdata->peft->index_N12x) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * pow(eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta), 2) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]);
    }
    else if (index_moment == userdata->peft->index_N12y) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = 0.;
    }
    /** --------------------- 4-th RSD moments -------------------------------- */
    else if (index_moment == userdata->peft->index_N22x) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = measure[it] * plin_rsd_q[it] * plin_rsd_kmq[it] * pow(eft_di_los_kernel_q(k_sq, q_sq[it], kmq_sq[it], cos_thetaq[it], cos_phiq[it], cos_theta) * eft_di_los_kernel_kmq(k_sq, q_sq[it], kmq_sq[it], cos_thetakmq[it]), 2);
    }
    else if (index_moment == userdata->peft->index_N22y) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = 0.;
    }
    else if (index_moment == userdata->peft->index_N22z) {
      for (it = 0; it < *nvec; it++)
        f[it][index_comp] = 0.;
    }

    index_comp++;
  }

  return _SUCCESS_;
}