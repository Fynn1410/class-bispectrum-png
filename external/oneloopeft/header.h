/** @file header.h
 *
 * contains the eft module and parameter struct definitions and includes all other source files needed for CLASS-OneLoop
 */

#ifndef __EFFECTIVE_FIELD_THEORY__
#define __EFFECTIVE_FIELD_THEORY__


#include "../../include/common.h" //Use here ONLY the things required for defining the struct (i.e. common.h for the ErrorMsg)
#include "../../include/background.h"
#include "../../include/primordial.h"
#include "../../include/fft.h"
#include "common_header.h"
#include "parallel.h"

#include <limits.h>


enum eft_pk_type {pk_lin, pk_nowiggle, pk_ir_resummed_lo, pkmu_rsd_ir_resummed_lo, pk_ir_resummed_nlo, pkmu_rsd_ir_resummed_nlo, pk_type_num};
enum eft_pk_out_type {Pdd_mm_real, Pdd_mm_rsd, Pdd_hh_real, Pdd_hh_rsd, Pdd_mm_real_no_IR_resum, Pdd_mm_22, Pdd_mm_13, Pdd_mm_22_no_IR_resum, Pdd_mm_13_no_IR_resum, pk_out_type_num};
enum eft_arg_type {points, cartesian_product};
enum eft_spectra_contribution {finite_part, uv_divergence, ir_divergence, pole_divergence, eft_spectra_contribution_num};
enum sym_type {no_finite_part, sym_vec, sym_mat_none, sym_mat_symmetric};


static inline int eft_rsd_indicator(enum eft_pk_type pk_type) {
  switch (pk_type)
  {
  case pk_lin:
  case pk_nowiggle:
  case pk_ir_resummed_lo:
  case pk_ir_resummed_nlo:
    return _FALSE_;
    break;

  case pkmu_rsd_ir_resummed_lo:
  case pkmu_rsd_ir_resummed_nlo:
    return _TRUE_;
    break;

  default:
    return -1;
    break;
  }
}

static inline int eft_rsd_out_indicator(enum eft_pk_out_type pk_out_type) {
  switch (pk_out_type)
  {
  case Pdd_mm_real:
  case Pdd_mm_real_no_IR_resum:
  case Pdd_hh_real:
  case Pdd_mm_22:
  case Pdd_mm_13:
  case Pdd_mm_22_no_IR_resum:
  case Pdd_mm_13_no_IR_resum:
    return _FALSE_;
    break;

  case Pdd_mm_rsd:
  case Pdd_hh_rsd:
    return _TRUE_;
    break;

  default:
    return -1;
    break;
  }
}

#include "wnw_split.h"
#include "kernel_matrices.h"

struct indexed_double
{
  int index;
  double value;
};



struct eft_moment_single  /**< structure for lin. power spectrum moments with 1 index */
{
  double moment;
  short index_bias;
};

struct eft_moment_double  /**< structure for lin.power spectrum moments with 2 indices */
{
  double moment;
  short index_bias;
  short index_derivative;
};

#define EFT_DISPERSION_SIZE 2
#define EFT_UV_CORRECTIONS_SIZE 9
#define EFT_UV_CORRECTIONS_UNDERLYING_SIZE 6

struct eft
{
  double z0;
  double f_z0;
  double D_z0;
  double * ln_k;  /**< ln(k * Mpc) for every mu; indexed as ln_k[index_mu*k_size + index_k] -> size = k_size * mu_size */
  int k_size;
  double * mu;  /**< cos(theta) to line-of-sight */
  int mu_size;
  double * pk_l[pk_type_num];
  double * ddpk_l[pk_type_num];

  int * spectra_contributions_size;  /**< size of every part of each moment in spectra_contributions[index_pk_type*index_num + index_moment] */
  double ** spectra_contributions[pk_type_num];  /**< spectra_contributions[index_pk_type][index_moment*eft_spectra_contribution_num + index_part][index_mu*k_size + index_k] */
  short * spectra_contributions_dimension;  /**< spectra_contributions_dimension[index_moment]: physical wavenumber dimensions of the contribution without any prefactors */

  double * ln_k_moments;
  double * pk_l_moments[pk_type_num];  /**< samples of the linear power spectra for moment computation */
  double * ddpk_l_moments[pk_type_num];  /**< second derivative of the linear power spectra for moment computation */
  struct eft_moment_single dispersion[pk_type_num][EFT_DISPERSION_SIZE];  /**< integrals of dq/(2 pi^2) q^2(n+1) P_lin(q) from k_IR_cutoff to k_UV_cutoff; n = index_bias; strongly UV-cutoff dependent for n >= 0 and strongly IR-cutoff dependent for n <= -2 */
  struct eft_moment_double ps_uv_shot_noise_corrections[pk_type_num][EFT_UV_CORRECTIONS_SIZE];  /**< integrals of dq/(2 pi^2) q^-(n+m) P_lin(q) (d^m P_lin(q) / dq^m) from k_IR_cutoff to k_UV_cutoff; n = index_bias, m = index_derivative; strongly UV-cutoff dependent for 2m + n <= -5 and strongly IR-cutoff dependent for 2m + n >= 3 */
  struct eft_moment_double ps_uv_shot_noise_corrections_underlying[pk_type_num][EFT_UV_CORRECTIONS_UNDERLYING_SIZE];  /**< integrals of dq/(2 pi^2) q^2(n+1) (d^m P_lin(q) / dln(q)^m)^2 from k_IR_cutoff to k_UV_cutoff; n = index_bias, m = index_derivative; strongly UV-cutoff dependent for n >= 3/2 and strongly IR-cutoff dependent for n <= -5/2 */

  double Sigma2_ir;
  double dSigma2_ir;

  double * pk_nl[pk_out_type_num];  /**< final nonlinear power spectrum pk_nl[index_pk_out_type][index_mu*k_size + index_k] */
  double * ddpk_nl[pk_out_type_num];  /**< second derivative of the final nonlinear power spectrum pk_nl[index_pk_out_type][index_mu*k_size + index_k]; spline of P_nl(ln(k)) */

  /** @name - indices of spectra contributions (guaranteed to be consecutive in {0,..., index_num}) */
  //@{

  int index_num;

  /** - Matter terms */
  int index_I2200;
  int index_I1300;

  /** - Halo terms / 0-th moment of RSD expansion */
  int index_Idelta200;
  int index_IG200;
  int index_Idelta2delta200;
  int index_IG2G200;
  int index_Idelta2G200;
  int index_FG200;

  /** ------ has_rsd = TRUE -------- */
  /** - 1-st moment of RSD expansion */
  int index_I2201;
  int index_Idelta201;
  int index_IG201;
  int index_J21101;
  int index_Jdelta201;
  int index_JG201;
  int index_FG201;
  int index_I1301p3101;
  int index_J12101;
  int index_J11201;

  /** - 2-nd moment of RSD expansion */
  int index_J21102x;
  int index_J21102y;
  int index_Jdelta202x;
  int index_Jdelta202y;
  int index_JG202x;
  int index_JG202y;
  int index_I2211;
  int index_J21111;
  int index_N11x;
  int index_N11y;
  int index_J12102x;
  int index_J12102y;
  int index_I1311;
  int index_J12111;
  int index_J11211;

  /** - 3-rd moment of RSD expansion */
  int index_J21112x;
  int index_J21112y;
  int index_N12x;
  int index_N12y;
  int index_J12112x;
  int index_J12112y;

  /** - 4-th moment of RSD expansion */
  int index_N22x;
  int index_N22y;
  int index_N22z;

  /** - additional RSD moments if the mu-approximation is not used */
  int index_sigmav_mu;

  //@}

  double * ln_k_fourier[eft_tracer_num];  /**< sample points for Fourier transform; ln_k_fourier[index_tracer][index_k] */
  double * pk_l_biased[eft_tracer_num*pk_type_num]; /**< biased linear spectra for Fourier transform; pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][index_mu*k_size + index_k], but mu_size is always 1 for pk_lin, pk_nowiggle and pk_ir_resummed_... */
  double * ddpk_l_biased[eft_tracer_num*pk_type_num];

  int moments_allocated;
  int * loop_matrices_size; /**< loop_matrices_size[index_moment] */
  class_complex ** loop_matrices;  /**< loop_matrices[index_moment][index_matrix] */
  short * symmetry; /**< symmetry[index_moment] */
  short * use_tracer; /**< use_tracer[index_moment] */
  short pk_type_loaded[pk_type_num];  /**< keeps track of which power spectra the Fourier coefficients were computed */
  class_complex * fourier_coeff[eft_tracer_num*pk_type_num];  /**< fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*fourier_coeff_size + index_freq] */
  class_complex * fourier_condition_num[eft_tracer_num*pk_type_num];  /**< fourier_condition_num[index_pk_type*eft_tracer_num + index_tracer][index_mu*fourier_coeff_size + index_freq] */
  double * fourier_frequencies[eft_tracer_num]; /**< fourier_frequencies[index_tracer][index_freq] */

  struct FFT_plan * fft_plan;

  struct eft_hyper_parameters * hp;
  struct eft_input_parameters * ip;

  short role;
  ErrorMsg error_message;
};


#include "../../include/fourier.h"


/** - forward declaration of ext_storage methods */
struct ext_storage;
int ext_insert_eft(struct ext_storage * pext,
                   struct eft * peft,
                   const int index_eft,
                   const int num_matrices,
                   ErrorMsg errmsg);
/** -------------------------------------------- */


int eft_init(
        struct precision * ppr,
        struct background * pba,
        struct eft * peft,
        struct eft_hyper_parameters * eft_hp,
        struct eft_input_parameters * eft_ip,
        struct ext_storage * pext,
        const short eft_role,
        const int index
        );

int eft_spline_sample_points_nonuniform(
        const double k_min,
        const double k_feature,
        const double k_max,
        const double rel_amplitude,
        const double width,
        const int num_points,
        double * ln_k_sample,
        const short verbosity
        );

int eft_nearest_structure_in_time(
        struct eft * peft0,
        const int peft_size,
        struct background * pba,
        struct fourier * pfo,
        const double z,
        int * index_eft_min_dist,
        struct eft * peft_min_dist,
        ErrorMsg errmsg_out
        );

int eft_fourier_transform_linear_spectra(
        struct precision * ppr,
        struct background * pba,
        struct fourier * pfo,
        struct primordial * ppm,
        struct eft * peft,
        const double z,
        const double f,
        const double D,
        const int * const index_pk_types,
        const int index_pk_types_size
        );

int eft_get_fourier_frequencies(
        const double period,
        const int num_frequencies,
        double * frequencies
        );

int eft_get_loop_matrices(struct eft * peft,
                          struct ext_storage * pext,
                          const int index);

int eft_free(struct eft * peft);

int eft_linear_spectrum_real(
        struct background * pba,
        struct primordial * ppm,
        struct fourier * pfo,
        struct eft * peft,
        enum linear_or_logarithmic mode,
        const double * const ln_kvec,
        const int kvec_size,
        const int n_columns,
        const double z,
        const double f_z,
        const double D_z,
        const int index_pk_type,
        double * out_pk
        );

int eft_linear_spectrum_rsd(
        struct background * pba,
        struct primordial * ppm,
        struct fourier * pfo,
        struct eft * peft,
        enum linear_or_logarithmic mode,
        const double * const ln_kvec,
        const int kvec_size,
        const double * const muvec,
        const int muvec_size,
        enum eft_arg_type arg_type,
        const double z,
        const double f_z,
        const double D_z,
        const int index_pk_type,
        double * out_pk
        );

int eft_apply_ap_effect_in_place(
        double ** kvec,
        int * k_sizevec,
        double ** muvec,
        const int * mu_sizevec,
        const int z_size,
        const double * ap_parallel,
        const double * ap_perpendicular
        );

int eft_job_powerspectrum_wedges_ext_growth_rate(
        struct eft * peft0,
        const int peft_size,
        const double * const f_z_pk_eft,
        const double * const D_z_pk_eft,
        struct background * pba,
        struct fourier * pfo,
        struct primordial * ppm,
        struct precision * ppr,
        enum eft_pk_out_type pk_out_type,
        const double * const zvec,
        const double * const f_zvec,
        const double * const D_zvec,
        const double As_correction,
        const struct eft_input_parameters * peft_ip,
        const int z_size,
        double ** kvec,
        int * const k_sizevec,
        double ** muvec,
        const int * const mu_sizevec,
        double ** out_pkmu,
        double ** ddout_pkmu
        );

int eft_job_powerspectrum_wedges(
        struct eft * peft0,
        const int peft_size,
        struct background * pba,
        struct fourier * pfo,
        struct primordial * ppm,
        struct precision * ppr,
        enum eft_pk_out_type pk_out_type,
        const double * const zvec,
        const double As_correction,
        const struct eft_input_parameters * peft_ip,
        const int z_size,
        double ** kvec,
        int * const k_sizevec,
        double ** muvec,
        const int * const mu_sizevec,
        double ** out_pkmu,
        double ** ddout_pkmu
        );

int eft_job_powerspectrum_wedges_grid(struct eft * peft0,
                                      const int peft_size,
                                      struct background * pba,
                                      struct fourier * pfo,
                                      struct primordial * ppm,
                                      struct precision * ppr,
                                      enum eft_pk_out_type pk_out_type,
                                      const double * const zvec,
                                      const struct eft_input_parameters * peft_ip,
                                      const int z_size,
                                      double * k, // indexed as k[index_z + z_size*(index_mu + mu_size*index_k)]
                                      int k_size,
                                      double * mu, // indexed as mu[index_z + z_size*index_mu]
                                      int mu_size,
                                      double * out_pkmu // indexed as out_pkmu[index_z + z_size*(index_mu + mu_size*index_k)]
                                      );

int eft_job_powerspectrum_multipoles_ext_growth_rate(
        struct eft * peft0,
        const int peft_size,
        const double * const f_z_pk_eft,
        const double * const D_z_pk_eft,
        struct background * pba,
        struct fourier * pfo,
        struct primordial * ppm,
        struct precision * ppr,
        enum eft_pk_out_type pk_out_type,
        const double * const zvec,
        const double * const f_zvec,
        const double * const D_zvec,
        const double As_correction,
        const struct eft_input_parameters * peft_ip,
        const int z_size,
        double ** kvec,
        int * const k_sizevec,
        const double * ap_parallel,
        const double * ap_perpendicular,
        double ** out_pkl
        );

int eft_job_powerspectrum_multipoles(
        struct eft * peft0,
        const int peft_size,
        struct background * pba,
        struct fourier * pfo,
        struct primordial * ppm,
        struct precision * ppr,
        enum eft_pk_out_type pk_out_type,
        const double * const zvec,
        const double As_correction,
        const struct eft_input_parameters * peft_ip,
        const int z_size,
        double ** kvec,
        int * const k_sizevec,
        const double * ap_parallel,
        const double * ap_perpendicular,
        double ** out_pkl
        );

#include "infrared_resummation.h"
#include "power_spectrum.h"
//#include "direct_integration.h"

#endif
