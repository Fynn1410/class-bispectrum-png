
#ifndef __EFFECTIVE_FIELD_THEORY_COMMON__
#define __EFFECTIVE_FIELD_THEORY_COMMON__

#define gfilter_smoothing_scale(ln_k) (0.6907755279 * exp( -pow((ln_k - (-3.4538776395)) / 2.752115578658, 2) ) + 0.06907755279)

//Using an enum to define the length of the enum is discouraged
enum eft_tracer {eft_matter, eft_halo, eft_tracer_num}; /**< nothing more than the number of distinct bias choices we want to use */
enum eft_fourier_mode {fourier_mode_fft, fourier_mode_spline};
enum eft_integration_mode {fftlog, direct_integration};

enum eft_struct_role {eft_master, eft_slave};

#define NUM_MOMENTS 43



const static FileName eft_loop_matrix_files_default[NUM_MOMENTS] = {"I2200", "I1300", "Idelta200", "IG200", "Idelta2delta200", "IG2G200", "Idelta2G200", "FG200",                 \
                                                                    "I2201", "Idelta201", "IG201", "J21101", "Jdelta201", "JG201", "FG201", "I1301p3101", "J12101", "J11201",     \
                                                                    "J21102x", "J21102y", "Jdelta202x", "Jdelta202y", "JG202x", "JG202y", "I2211", "J21111", "N11x", "N11y", "J12102x", "J12102y", "I1311", "J12111", "J11211",   \
                                                                    "J21112x", "J21112y", "N12x", "N12y", "J12112x", "J12112y",                                                   \
                                                                    "N22x", "N22y", "N22z"};  /**< default filenames for the loop matrices */


int eft_ln_pk_nw_gfilter_parallel(struct precision *ppr, struct background *pba, struct primordial *ppm, struct fourier *pfo, const int index_pk, const int index_k0, const int index_kmin, const int k_size, double *ln_pknw_array);

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
  int linear_spectrum_index;  /**< index of the linear spectrum in CLASS on which to operate (e.g. pfo->index_pk_cluster)*/

  /** - Infrared resummation settings */
  double ir_resummation_k_split;  /**< upper bound of integration region for the IR resummation suppression factors in the exponential (in h/Mpc); default 0.2 h/Mpc */
  double ir_resummation_k_feature;  /**< position of the BAO feature to be resummed (in h/Mpc); default 1/110. h/Mpc */

  /** - fourier transform settings */
  enum eft_fourier_mode fourier_mode; /**< use either FFT (=1) or exact Spline Fourier transform (=0) */
  double bias[eft_tracer_num];  /**< FFTLog bias parameter */
  int num_positive_fourier_freq;  /**< number of positive frequencies to sample */
  int fourier_coeff_size; /**< number of Fourier components = 2*num_positive_fourier_freq + 1 */

  /** - divergence / direct integration settings */
  double k_UV_cutoff;
  double k_IR_cutoff;
  double k_pole_cutoff; /**< |k-q| > pole_cutoff (only relevant for direct integration) */
  int k_size_moments;

  /** - output sampling settings: in effect if use_interpolation = _TRUE_ */
  short use_interpolation;  /**< if set, precompute the spectra and interpolate */
  int k_size_nl;
  double kmin_nl;
  double kmax_nl;
  double k_feature_nl;
  double ln_k_oversampling_width_nl;

  enum eft_integration_mode integration_mode; /**< either FFTLog or Direct integration using Cuba */
  short has_rsd; /**< if set, all RSD moments will be allocated */
  short use_EdS_time_scaling; /**< if set, use EdS scaling relations for the time-dependence of the non-linear spectra, otherwise interpolate btw. z_pk_eft */
  short compute_loop_matrices;  /**< if set, compute all loop matrices from scratch, otherwise load them from file */
  short reload_linear_spectra;  /**< if set, reload the linear spectra at each evaluation, not compatible with use_interpolation */

  short write_loop_matrices;
  short ignore_missing_files; /**< ignore errors when loading loop matrices from files */
  FileName eft_loop_matrix_directory;
  FileName eft_loop_matrix_files[NUM_MOMENTS];

  short eft_verbose;
};

#endif
