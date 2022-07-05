

/** @file header.h
 *
 */
#ifndef HEADER_H_
#define HEADER_H_

#define _GNU_SOURCE

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

/**
 * List of limHaloPT header files
 */
#include "cosmology.h"
#include "wnw_split.h"
#include "IR_res.h"
#include "ps_halo_1loop.h"



/**
 * Function declarations of main.c module
 */
void  initialize();
void 	cleanup();


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

	long p13;
	long p14;
	long p15;
	long p16;
	long p17;
	long p18;

	int p19;

	double *p20;
	size_t p22;
};

typedef struct fft_struct
{
	int nfft;
	int fft_first;
	double kmin_fft_m;
	double kmin_fft_g;
	double fft_bias_m;
	double fft_bias_g;
	double complex * etam_m;  
  double complex * cmsym_m; 
  double complex * etam_g;  
  double complex * cmsym_g;
       
}fft_struct;

typedef struct fft_matrices
{
  // Matter loops
  double complex ** I2200_mat; 
  double complex *  I1300_mat;

  // Realspace (0-th moment in RSD) Biased loops
  double complex ** Idelta200_mat; 
  double complex ** IG200_mat; 
  double complex ** Idelta2delta200_mat; 
  double complex ** IG2G200_mat; 
  double complex ** Idelta2G200_mat; 
  double complex *  FG200_mat; 

  // 1-st moment in RSD Biased loops
  double complex ** I2201_mat; 
  double complex ** Idelta201_mat; 
  double complex ** IG201_mat; 
  double complex ** FG201_mat; 
  double complex ** J21101_mat; 
  double complex ** Jdelta201_mat; 
  double complex ** JG201_mat; 
  double complex *  I1301_mat; 
  double complex *  J12101_mat; 
  
  // 2-nd moment in RSD Biased loops
  double complex ** J21102x_mat;
  double complex ** J21102y_mat;
  double complex ** Jdelta202x_mat;
  double complex ** Jdelta202y_mat;
  double complex ** JG202x_mat;
  double complex ** JG202y_mat;
  double complex ** I2211_mat;
  double complex ** J21111_mat;
  double complex ** N11x_mat;
  double complex ** N11y_mat;
  double complex *  J12102x_mat;
  double complex *  J12102y_mat;
  double complex *  I1311_mat;
  double complex *  J12111_mat;
  double complex *  J11211_mat;

  // 3-rd moment in RSD Biased loops
  double complex ** J21112x_mat;
  double complex ** J21112y_mat;
  double complex ** N12x_mat;
  double complex ** N12y_mat;
  double complex *  J12112x_mat;
  double complex *  J12112y_mat;

  // 4-th moment in RSD Biased loops
  double complex ** N22x_mat;
  double complex ** N22y_mat;
  double complex ** N22z_mat;

}fft_matrices;

typedef struct rsd_bias
{
  double *np_bias_vec0;
  double *p_bias_vec0;

  double *np_bias_vec1;
  double *p_bias_vec1;

  double *np_bias_vec2;
  double *p_bias_vec2;

  double *np_bias_vec3;
  double *p_bias_vec3;

  double *np_bias_vec4;
  double *p_bias_vec4;

}rsd_bias;

/**
 * Structure containing variables, calculated in fourier and used only in nl_oneloopPT by various functions.
 *
 */

struct oneloop_fftlog_workspace {

  /** @name - quantitites used by nl_oneloopPT */

  //@{

  double z; /** Analyzed red-shift **/

  double mu; /** Scalar-product between LoS-direction and direction of the analyzed wavenumber (normalized) **/

  double f;  /** Velocity divergence growth rate **/

  double sigma_v2; /** Value of the integrated linear power spectrum (for the UV- / IR-divergences of the integrals) **/

  double sigma_2_IR; /** Value of the supression factor of the wiggle part for the IR-resummation **/

  struct fft_struct *fft_input; /** Containing the details of the FFTLog and the FFTLog transform of the IR-resummed power spectrum **/

  // FFTLog matrices for the Loop-Integrals
  struct fft_matrices *fft_matrix;  /** Containing the complex matrices used for the FFTLog loop calculations **/

  // Eulerian Biases
  double b1;
  double b2;
  double bG2;
  double btd;
  double cs2;
  double R2;

  // RSD Bias vectors
  struct rsd_bias *bias; /** Bias vector for distributing in the end on the loops **/

  //@}

};

#include "utilities.h"
#include "FFTLog_ingredients.h"
#include "FFTLog_matter_real.h"
#include "FFTLog_galaxy_real.h"
#include "FFTLog_galaxy_rsd.h"
#include "pspec_FFTLog_real.h"
#include "pspec_FFTLog_rsd.h"

#endif
