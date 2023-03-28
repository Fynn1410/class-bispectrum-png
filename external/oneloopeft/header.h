

/** @file header.h
 *
 */
#ifndef HEADER_H_
#define HEADER_H_


#include "../../include/common.h" //Use here ONLY the things required for defining the struct (i.e. common.h for the ErrorMsg)
#include "../../include/background.h"
#include "../../include/fourier.h"
#include "../../include/primordial.h"

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

enum eft_pk_type {pk_lin, pk_nowiggle, pk_ir_resummed_lo, pk_ir_resummed_nlo, pk_index_num};
enum sym_type {vec, mat_none, mat_symmetric};

/**
 * List of limHaloPT header files
 */
#include "cosmology.h"
#include "wnw_split.h"
#include "IR_res.h"
// #include "ps_halo_1loop.h"



struct eft
{
  double z0;
  double * ln_k;
  int k_size;
  double ** pk_l;

  double ** spectra_contributions;

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

  short has_rsd; /**< if set, all moments will be allocated */

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
  
  double bias_matter;
  double bias_halo;
  double * ln_k_matter_fourier;
  double * ln_k_halo_fourier;
  int k_size_fourier;
  double ** pk_l_biased;
  double ** ddpk_l_biased;

  int * loop_matrices_size;
  double complex ** loop_matrices;
  short * symmetry;
  int fourier_coeff_size;
  double period;
  double complex ** fourier_coeff;
  double * fourier_frequencies;

  ErrorMsg error_message;
};


int eft_init(struct background * pba,
            struct fourier * pfo,
            struct primordial * ppm,
            struct eft * peft,
            double z,
            short compute_rsd_spectrum,
            double k_min,
            double fourier_period,
            int num_positive_coefficients,
            int num_sample_points);

int eft_get_fourier_frequencies(
                  const double period,
                  const int num_frequencies,
                  double * frequencies);


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
