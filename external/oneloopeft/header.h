

/** @file header.h
 *
 */
#ifndef __EFFECTIVE_FIELD_THEORY__
#define __EFFECTIVE_FIELD_THEORY__


#include "../../include/common.h" //Use here ONLY the things required for defining the struct (i.e. common.h for the ErrorMsg)
#include "../../include/background.h"
#include "../../include/primordial.h"
#include "../../include/fft.h"

//TODO: get rid of forward declarations
struct ext_storage
{
  /** - EFT loop matrix pointers */
  double complex *** loop_matrices;
  int ** loop_matrices_size;
  short ** symmetry;
  int loop_matrices_stored; // = _FALSE_;
  int eft_index_num; // = 0;
  int eft_size; // = 0;
};

//#include <time.h>
//#include <unistd.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <float.h>
//#include <string.h>
//#include <omp.h>
//#include <mpi.h>

#include <complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_odeiv2.h>  /// For solving ODER
#include <gsl/gsl_roots.h>	 // For finding the root of algebraic equation
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <fftw3.h>	// for DST
//#include <ctype.h>
//#include "../Class/include/class.h"
#include "library/Cuba-4.2.1/cuba.h"


#define PSC  		101L
#define ST 		    	102L
#define TR 		    	103L

#define GROWTH		104L
#define DERGROWTH		105L

#define NONLINEAR	   	106L
#define LINEAR		107L

#define GAUSSIAN		114L
#define NONGAUSSIAN	115L

#define INIT  		116L
#define LOCAL		117L
#define EQUILATERAL	118L
#define ORTHOGONAL	119L
#define QSF			120L
#define HS			121L
#define NGLOOP		122L
#define derNGLOOP		123L

#define QUADRATIC		124L
#define TIDE		125L
#define GAMMA		126L

#define LPOWER		127L
#define NLPOWER		128L
#define TRANS		129L
#define DER 		130L

#define CO10  		131L
#define CO21  		132L
#define CO32  		133L
#define CO43  		134L
#define CO54  		135L
#define CO65  		136L
#define CII 		137L

#define MATTER		138L
#define LINEMATTER	139L
#define LINE 		140L

#define DST			141L
#define GFILTER		142L
#define BSPLINE		143L


#define TREE		144L
#define LOOP		145L
#define WIR			146L
#define NOIR		147L

#define HALO 		148L


#define CDM 		90L
#define BA 			91L
#define TOT 		92L


#define PS_KMIN        1.0e-7
#define PS_KMAX        1.0e4

#define CLEANUP        1

#define DO_NOT_EVALUATE -1.0

#define MAXL 2000



#define NUM_MOMENTS 40

enum eft_struct_role {eft_master, eft_slave};
enum eft_tracer {eft_matter, eft_halo, eft_tracer_num};
enum eft_pk_type {pk_lin, pk_nowiggle, pk_ir_resummed_lo, pk_ir_resummed_nlo, pk_type_num};
enum sym_type {vec, mat_none, mat_symmetric};
enum fourier_mode {fourier_mode_fft, fourier_mode_spline};


const static FileName eft_loop_matrix_files_default[NUM_MOMENTS] = {"I2200.mat", "I1300.mat", "Idelta200.mat", "IG200.mat", "Idelta2delta200.mat", "IG2G200.mat", "Idelta2G200.mat", "FG200.mat",         \
                                                            "I2201.mat", "Idelta201.mat", "IG201.mat", "J21101.mat", "Jdelta201.mat", "JG201.mat", "FG201.mat", "I1301p3101.mat", "J12101.mat",   \
                                                            "J21102x.mat", "J21102y.mat", "Jdelta202x.mat", "Jdelta202y.mat", "JG202x.mat", "JG202y.mat", "I2211.mat", "J21111.mat", "N11x.mat", "N11y.mat", "J12102x.mat", "J12102y.mat", "I1311.mat", "J12111p11211.mat",   \
                                                            "J21112x.mat", "J21112y.mat", "N12x.mat", "N12y.mat", "J12112x.mat", "J12112y.mat",                                                   \
                                                            "N22x.mat", "N22y.mat", "N22z.mat"};  /**< default filenames for the loop matrices */

/**
 * List of limHaloPT header files
 */
#include "cosmology.h"
#include "wnw_split.h"
#include "IR_res.h"
#include "kernel_matrices.h"
// #include "ps_halo_1loop.h"


struct eft_input_parameters
{
  double b1;
  double b2;
  double bG2;
  double btd;

  double cs2;
  double R2;

  short has_rsd;
  double c00;
  double c10;
  double c20;
  double c22;
  double c30;
  double c32;
  double c42;
};

struct eft_hyper_parameters
{
  /** - linear spectrum sampling settings */
  double kmin_lin[eft_tracer_num];  /**< begin of sampling range */
  double kmax_lin[eft_tracer_num];  /**< end of sampling range */
  double period[eft_tracer_num];  /**< period in logarithmic space, must be larger than ln(kmax/kmin) */
  int k_size_fourier; /**< number of sample points + 1 of the linear spectra */
  double bao_oversampling;  /**< sample point density increases by this factor at the BAO */
  double ln_k_oversampling_width; /**< width of the region in which to oversample: n(ln(k)) = oversampling/(1 + (ln(k/k_feature)/width)^2) */

  /** - fourier transform settings */
  short fourier_mode; /**< use either FFT (=1) or exact Spline Fourier transform (=0) */
  double bias[eft_tracer_num];  /**< FFTLog bias parameter */
  int num_positive_fourier_freq;  /**< number of positive frequencies to sample */
  int fourier_coeff_size; /**< number of Fourier components = 2*num_positive_fourier_freq + 1 */

  short has_rsd; /**< if set, all RSD moments will be allocated */
  short use_EdS_time_scaling; /**< if set, use EdS scaling relations for the time-dependence of the non-linear spectra, otherwise interpolate btw. z_pk_eft */
  short use_time_independent_kernels; /**< if set, use EdS approximation for the loop matrices [time-dependence is not yet implemented] */
  short compute_loop_matrices;  /**< if set, compute all loop matrices from scratch, otherwise load them from file */
  short compute_mu_approximation; /**< if set, compute the halo power spectrum corrections on P_lin(k, z0) and P_nw(k, z0) only, else on P_mm^RSD,LO(k, z0, mu) for different mu */

  FileName eft_loop_matrix_files[NUM_MOMENTS];

  short eft_verbose;
};

struct eft
{
  double z0;
  double f_z0;
  double D_z0;
  double * ln_k;
  int k_size;
  double * mu;
  int mu_size;
  double ** pk_l;

  double ** spectra_contributions[pk_type_num];  /**< spectra_contributions[index_pk_type][index_moment][index_mu*k_size + index_k] */

  double sigma_v0;
  double sigma_v2;
  double Sigma2_ir;
  double dSigma2_ir;

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
  int index_J12111p11211;
  //int index_J11211;

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

  //@}

  double * ln_k_fourier[eft_tracer_num];  /**< sample points for Fourier transform; ln_k_fourier[index_tracer][index_k] */
  double * pk_l_biased[eft_tracer_num*pk_type_num]; /**< biased linear spectra for Fourier transform; pk_l_biased[index_pk_type*eft_tracer_num + index_tracer][index_mu*k_size + index_k], but mu_size is always 1 for pk_lin and pk_nowiggle */
  double * ddpk_l_biased[eft_tracer_num*pk_type_num];

  int * loop_matrices_size; /**< loop_matrices_size[index_moment]*/
  double complex ** loop_matrices;  /**< loop_matrices[index_moment][index_matrix] */
  short * symmetry; /**< symmetry[index_moment]*/
  double complex * fourier_coeff[eft_tracer_num*pk_type_num];  /**< fourier_coeff[index_pk_type*eft_tracer_num + index_tracer][index_mu*fourier_coeff_size + index_freq] */
  double complex * fourier_condition_num[eft_tracer_num*pk_type_num];  /**< fourier_condition_num[index_pk_type*eft_tracer_num + index_tracer][index_mu*fourier_coeff_size + index_freq] */
  double * fourier_frequencies[eft_tracer_num]; /**< fourier_frequencies[index_tracer][index_freq] */

  struct FFT_plan * fft_plan;

  struct eft_hyper_parameters * hp;
  struct eft_input_parameters * ip;

  short role;
  ErrorMsg error_message;
};


#include "../../include/fourier.h"




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

int eft_load_linear_spectra(
        struct precision * ppr,
        struct background * pba,
        struct fourier * pfo,
        struct primordial * ppm,
        struct eft * peft,
        double z,
        double f,
        double D,
        const short use_mu_approximation
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


/**
* A structure passed to the integrators to hold the parameters fixed in the integration
*/
struct integrand_parameters
{
	double p1;
	double p2;
	double p3;
	double p4;
	double p5;
	double p6;
	double p7;
	double p8;
	double p9;
	double p10;
	double p11;
	long 	p12;
	long 	p13;
};

/**
* Another structure passed to the integrators to hold the parameters fixed in the integration
*/
struct integrand_parameters2
{

    struct background *pba;
	struct primordial *ppm;
    struct fourier *pfo;

	double p4;
	double p5;
	double p6;
	double p7;
	double p8;
	double p9;
	double p10;
	double p11;
	double p12;

	double p24;
	double p25;
	double p26;
	double p27;

	long p13;
	long p14;
	long p15;
	long p16;
	long p17;
	long p18;

	int p19;
	int p23;

	double *p20;
	size_t p22;
};

#include "utilities.h"
#include "FFTLog_ingredients.h"
#include "FFTLog_real.h"
#include "FFTLog_rsd.h"
#include "pspec_FFTLog_real.h"
#include "pspec_FFTLog_rsd.h"

#endif
