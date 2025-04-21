# Bunch of declarations from C to python. The idea here is to define
# only the quantities that will be used, for input, output or
# intermediate manipulation, by the python wrapper. For instance, in
# the precision structure, the only item used here are the error
# message and one parameter used for an error message. That is why
# nothing more is defined from this structure. The rest is internal in
# Class.  If, for whatever reason, you need an other, existing
# parameter from Class, remember to add it inside this cdef.

cdef extern from "class.h":

    cdef char[10] _VERSION_

    ctypedef char FileArg[40]

    ctypedef char* ErrorMsg

    ctypedef char FileName[256]

    cdef enum interpolation_method:
        inter_normal
        inter_growing_closeby

    cdef enum vecback_format:
        short_info
        normal_info
        long_info

    cdef enum linear_or_logarithmic:
        linear
        logarithmic

    cdef enum file_format:
        class_format
        camb_format

    cdef enum non_linear_method:
        nl_none
        nl_halofit
        nl_HMcode
        nl_oneloopPT

    cdef enum pk_outputs:
        pk_linear
        pk_nonlinear

    cdef enum out_sigmas:
        out_sigma
        out_sigma_prime
        out_sigma_disp

    cdef enum eft_fourier_mode:
        fourier_mode_fft
        fourier_mode_spline

    cdef enum eft_integration_mode:
        fftlog
        direct_integration

    cdef struct FFT_plan:
        int N;
        double* cos_vals;
        double* sin_vals;

        double* temp_real;
        double* temp_imag;

    cdef struct ext_storage:
        double complex *** loop_matrices;
        int ** loop_matrices_size;
        short ** symmetry;
        short ** use_tracer;
        short ** spectra_contributions_dimension;

        double ** period;
        int loop_matrices_stored;
        int * eft_index_num;
        int eft_size;

        ErrorMsg error_message;

    cdef struct precision:
        double nonlinear_min_k_max
        ErrorMsg error_message

    cdef struct background:
        short is_allocated
        ErrorMsg error_message
        int bg_size
        int index_bg_ang_distance
        int index_bg_lum_distance
        int index_bg_conf_distance
        int index_bg_a
        int index_bg_H
        int index_bg_D
        int index_bg_f
        int index_bg_Omega_m
        int index_bg_rho_b
        int index_bg_rho_cdm
        int index_bg_rho_ncdm1
        int index_bg_rho_crit
        short has_cdm
        short  has_ncdm
        int N_ncdm
        double T_cmb
        double h
        double H0
        double age
        double conformal_age
        double K
        double * m_ncdm_in_eV
        double Neff
        double Omega0_g
        double Omega0_b
        double Omega0_idr
        double T_idr
        double Omega0_cdm
        double Omega0_idm
        double Omega0_dcdm
        double Omega0_ncdm_tot
        double Omega0_lambda
        double Omega0_fld
        double w0_fld
        double wa_fld
        double cs2_fld
        double Omega0_ur
        double Omega0_dcdmdr
        double Omega0_dr
        double Omega0_scf
        double Omega0_k
        int bt_size
        double Omega0_m
        double Omega0_r
        double Omega0_de
        double a_eq
        double H_eq
        double z_eq
        double tau_eq

    cdef struct thermodynamics:
        short is_allocated
        ErrorMsg error_message
        int th_size
        int index_th_xe
        int index_th_Tb
        double tau_reio
        double z_reio
        double z_rec
        double tau_rec
        double rs_rec
        double ds_rec
        double da_rec
        double z_star
        double tau_star
        double rs_star
        double ds_star
        double ra_star
        double da_star
        double rd_star
        double z_d
        double tau_d
        double ds_d
        double rs_d
        double YHe
        double n_e
        double a_idm_dr
        double b_idr
        double nindex_idm_dr
        double m_idm
        double cross_idm_g
        double u_idm_g
        double cross_idm_b
        double n_index_idm_b
        int tt_size

    cdef struct perturbations:
        short is_allocated
        ErrorMsg error_message
        short has_scalars
        short has_vectors
        short has_tensors

        short has_density_transfers
        short has_velocity_transfers

        int has_pk_matter
        int l_lss_max

        int store_perturbations
        int k_output_values_num
        double k_output_values[30]
        double k_max_for_pk
        int index_k_output_values[30]
        char scalar_titles[8000]
        char vector_titles[8000]
        char tensor_titles[8000]
        int number_of_scalar_titles
        int number_of_vector_titles
        int number_of_tensor_titles
        int index_md_scalars

        double * scalar_perturbations_data[30]
        double * vector_perturbations_data[30]
        double * tensor_perturbations_data[30]
        int size_scalar_perturbation_data[30]
        int size_vector_perturbation_data[30]
        int size_tensor_perturbation_data[30]

        double * alpha_idm_dr
        double * beta_idr

        # add source functions for comparison
        short has_source_t
        short has_source_p
        short has_source_delta_m
        short has_source_delta_cb
        short has_source_delta_tot
        short has_source_delta_g
        short has_source_delta_b
        short has_source_delta_cdm
        short has_source_delta_idm
        short has_source_delta_idr
        short has_source_delta_dcdm
        short has_source_delta_fld
        short has_source_delta_scf
        short has_source_delta_dr
        short has_source_delta_ur
        short has_source_delta_ncdm
        short has_source_theta_m
        short has_source_theta_cb
        short has_source_theta_tot
        short has_source_theta_g
        short has_source_theta_b
        short has_source_theta_cdm
        short has_source_theta_idm
        short has_source_theta_idr
        short has_source_theta_dcdm
        short has_source_theta_fld
        short has_source_theta_scf
        short has_source_theta_dr
        short has_source_theta_ur
        short has_source_theta_ncdm
        short has_source_phi
        short has_source_phi_prime
        short has_source_phi_plus_psi
        short has_source_psi
        short has_source_h
        short has_source_h_prime
        short has_source_eta
        short has_source_eta_prime
        short has_source_H_T_Nb_prime
        short has_source_k2gamma_Nb

        int index_tp_t0
        int index_tp_t1
        int index_tp_t2
        int index_tp_p
        int index_tp_delta_m
        int index_tp_delta_cb
        int index_tp_delta_tot
        int index_tp_delta_g
        int index_tp_delta_b
        int index_tp_delta_cdm
        int index_tp_delta_idm
        int index_tp_delta_dcdm
        int index_tp_delta_fld
        int index_tp_delta_scf
        int index_tp_delta_dr
        int index_tp_delta_ur
        int index_tp_delta_idr
        int index_tp_delta_ncdm1

        int index_tp_theta_m
        int index_tp_theta_cb
        int index_tp_theta_tot
        int index_tp_theta_g
        int index_tp_theta_b
        int index_tp_theta_cdm
        int index_tp_theta_dcdm
        int index_tp_theta_fld
        int index_tp_theta_scf
        int index_tp_theta_ur
        int index_tp_theta_idr
        int index_tp_theta_idm
        int index_tp_theta_dr
        int index_tp_theta_ncdm1

        int index_tp_phi
        int index_tp_phi_prime
        int index_tp_phi_plus_psi
        int index_tp_psi
        int index_tp_h
        int index_tp_h_prime
        int index_tp_eta
        int index_tp_eta_prime
        int index_tp_H_T_Nb_prime
        int index_tp_k2gamma_Nb


        double *** sources
        double * tau_sampling
        int tau_size
        int k_size_pk
        int * k_size
        double ** k
        int * ic_size
        int index_ic_ad
        int md_size
        int * tp_size
        double * ln_tau
        int ln_tau_size

    cdef struct transfer:
        short is_allocated
        ErrorMsg error_message

    cdef struct primordial:
        short is_allocated
        ErrorMsg error_message
        double k_pivot
        double A_s
        double n_s
        double alpha_s
        double beta_s
        double r
        double n_t
        double alpha_t
        double V0
        double V1
        double V2
        double V3
        double V4
        double f_cdi
        double n_cdi
        double c_ad_cdi
        double n_ad_cdi
        double f_nid
        double n_nid
        double c_ad_nid
        double n_ad_nid
        double f_niv
        double n_niv
        double c_ad_niv
        double n_ad_niv
        double phi_min
        double phi_max
        int lnk_size

    cdef struct harmonic:
        short is_allocated
        ErrorMsg error_message
        int has_tt
        int has_te
        int has_ee
        int has_bb
        int has_pp
        int has_tp
        int has_dd
        int has_td
        int has_ll
        int has_dl
        int has_tl
        int l_max_tot
        int ** l_max_ct
        int ct_size
        int * ic_size
        int * ic_ic_size
        int md_size
        int d_size
        int non_diag
        int index_ct_tt
        int index_ct_te
        int index_ct_ee
        int index_ct_bb
        int index_ct_pp
        int index_ct_tp
        int index_ct_dd
        int index_ct_td
        int index_ct_pd
        int index_ct_ll
        int index_ct_dl
        int index_ct_tl
        int * l_size
        int index_md_scalars

    cdef struct output:
        ErrorMsg error_message

    cdef struct distortions:
        short is_allocated
        double * sd_parameter_table
        int index_type_g
        int index_type_mu
        int index_type_y
        int index_type_PCA
        int type_size
        double * DI
        double * x
        double DI_units
        double x_to_nu
        int has_distortions
        int x_size
        ErrorMsg error_message

    cdef struct lensing:
        short is_allocated
        int has_tt
        int has_ee
        int has_te
        int has_bb
        int has_pp
        int has_tp
        int has_dd
        int has_td
        int has_ll
        int has_dl
        int has_tl
        int index_lt_tt
        int index_lt_te
        int index_lt_ee
        int index_lt_bb
        int index_lt_pp
        int index_lt_tp
        int index_lt_dd
        int index_lt_td
        int index_lt_ll
        int index_lt_dl
        int index_lt_tl
        int * l_max_lt
        int lt_size
        int has_lensed_cls
        int l_lensed_max
        int l_unlensed_max
        ErrorMsg error_message

    struct eft_input_parameters:
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

    cdef struct eft_hyper_parameters:
        double kmin_lin[2];
        double kmax_lin[2];
        double period[2];
        int k_size_fourier;
        double bao_oversampling;
        double ln_k_oversampling_width;
        int linear_spectrum_index;

        double ir_resummation_k_split;
        double ir_resummation_k_feature;

        int fourier_mode;
        double bias[2];
        int num_positive_fourier_freq;
        int fourier_coeff_size;

        double k_UV_cutoff;
        double k_IR_cutoff;
        double k_pole_cutoff;
        int k_size_moments;

        short use_interpolation;
        int k_size_nl;
        double kmin_nl;
        double kmax_nl;
        double k_feature_nl;
        double ln_k_oversampling_width_nl;

        int integration_mode;
        short has_rsd;
        short use_EdS_time_scaling;
        short use_time_independent_kernels;
        short compute_loop_matrices;
        short use_mu_approximation;
        short reload_linear_spectra;

        short write_loop_matrices;
        short ignore_missing_files;
        FileName eft_loop_matrix_directory;
        FileName eft_loop_matrix_files[43];

        short eft_verbose;

    cdef struct eft_moment_single:
        double moment;
        short index_bias;

    cdef struct eft_moment_double:
        double moment;
        short index_bias;
        short index_derivative;

    cdef struct eft:
        double z0;
        double f_z0;
        double D_z0;
        double * ln_k;
        int k_size;
        double * mu;
        int mu_size;
        double * pk_l[10];
        double * ddpk_l[10];

        int * spectra_contributions_size;
        double ** spectra_contributions[10];
        short * spectra_contributions_dimension;

        double * ln_k_moments;
        double * pk_l_moments[6];
        double * ddpk_l_moments[6];
        eft_moment_single dispersion[6][2];
        eft_moment_double ps_uv_shot_noise_corrections[6][9];
        eft_moment_double ps_uv_shot_noise_corrections_underlying[6][6];

        double Sigma2_ir;
        double dSigma2_ir;

        int index_num;

        int index_I2200;
        int index_I1300;

        int index_Idelta200;
        int index_IG200;
        int index_Idelta2delta200;
        int index_IG2G200;
        int index_Idelta2G200;
        int index_FG200;

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

        int index_J21112x;
        int index_J21112y;
        int index_N12x;
        int index_N12y;
        int index_J12112x;
        int index_J12112y;

        int index_N22x;
        int index_N22y;
        int index_N22z;

        int index_sigmav_mu;

        double * ln_k_fourier[2];
        double * pk_l_biased[2*6];
        double * ddpk_l_biased[2*6];

        int moments_allocated;
        int * loop_matrices_size;
        double complex ** loop_matrices;
        short * symmetry;
        short * use_tracer;
        short pk_type_loaded[6];
        double complex * fourier_coeff[2*6];
        double complex * fourier_condition_num[2*6];
        double * fourier_frequencies[2];

        FFT_plan * fft_plan;

        eft_hyper_parameters * hp;
        eft_input_parameters * ip;

        short role;
        ErrorMsg error_message;

    cdef struct fourier:
        short is_allocated
        short has_pk_matter
        int method
        int ic_size
        int ic_ic_size
        int k_size
        int k_size_pk
        int ln_tau_size
        int index_ln_tau_pk
        int tau_size
        int index_tau_min_nl
        double * k
        double * ln_tau
        double * tau
        double ** ln_pk_l
        double ** ln_pk_nl
        eft * peft;
        int eft_size;
        double * z_pk_eft;
        int z_pk_eft_num;
        eft_input_parameters * eft_ip;
        eft_hyper_parameters eft_hp;
        double * sigma8
        int has_pk_m
        int has_pk_cb
        int index_pk_m
        int index_pk_cb
        int index_pk_total
        int index_pk_cluster
        ErrorMsg error_message

    cdef struct file_content:
        char * filename
        int size
        FileArg * name
        FileArg * value
        short * read

    void lensing_free(void*)
    void harmonic_free(void*)
    void transfer_free(void*)
    void primordial_free(void*)
    void perturbations_free(void*)
    void thermodynamics_free(void*)
    void background_free(void*)
    void fourier_free(void*)
    void distortions_free(void*)

    cdef int _FAILURE_
    cdef int _FALSE_
    cdef int _TRUE_

    int input_read_from_file(void*, void*, void*, void*, void*, void*, void*, void*, void*,
        void*, void*, void*, char*)
    int background_init(void*,void*)
    int thermodynamics_init(void*,void*,void*)
    int perturbations_init(void*,void*,void*,void*)
    int primordial_init(void*,void*,void*)
    int fourier_init(void*,void*,void*,void*,void*,void*,void*)
    int transfer_init(void*,void*,void*,void*,void*,void*)
    int harmonic_init(void*,void*,void*,void*,void*,void*,void*)
    int lensing_init(void*,void*,void*,void*,void*)
    int distortions_init(void*,void*,void*,void*,void*,void*)

    int background_tau_of_z(void* pba, double z,double* tau)
    int background_z_of_tau(void* pba, double tau,double* z)
    int background_at_z(void* pba, double z, int return_format, int inter_mode, int * last_index, double *pvecback)
    int background_at_tau(void* pba, double tau, int return_format, int inter_mode, int * last_index, double *pvecback)
    int background_output_titles(void * pba, char titles[8000])
    int background_output_data(void *pba, int number_of_titles, double *data)

    int thermodynamics_at_z(void * pba, void * pth, double z, int inter_mode, int * last_index, double *pvecback, double *pvecthermo)
    int thermodynamics_output_titles(void * pba, void *pth, char titles[8000])
    int thermodynamics_output_data(void *pba, void *pth, int number_of_titles, double *data)

    int perturbations_output_data_at_z(void *pba,void *ppt, file_format output_format, double z, int number_of_titles, double *data)
    int perturbations_output_data_at_index_tau(void *pba,void *ppt, file_format output_format, int ondex_tau, int number_of_titles, double *data)
    int perturbations_output_data(void *pba,void *ppt, file_format output_format, double * tkfull, int number_of_titles, double *data)
    int perturbations_output_firstline_and_ic_suffix(void *ppt, int index_ic, char first_line[1024], FileName ic_suffix)
    int perturbations_output_titles(void *pba, void *ppt,  file_format output_format, char titles[8000])

    int primordial_output_titles(void * ppt, void *ppm, char titles[8000])
    int primordial_output_data(void *ppt, void *ppm, int number_of_titles, double *data)

    int harmonic_cl_at_l(void* phr,double l,double * cl,double * * cl_md,double * * cl_md_ic)
    int lensing_cl_at_l(void * ple,int l,double * cl_lensed)

    int harmonic_pk_at_z(
        void * pba,
        void * phr,
        int mode,
        double z,
        double * output_tot,
        double * output_ic,
        double * output_cb_tot,
        double * output_cb_ic
        )
    int fourier_pk_at_z(
        void * pba,
        void *pfo,
        int mode,
        int pk_output,
        double z,
        int index_pk,
        double * out_pk,
        double * out_pk_ic
        )

    int harmonic_pk_at_k_and_z(
        void* pba,
        void * ppm,
        void * phr,
        double k,
        double z,
        double * pk,
        double * pk_ic,
        double * pk_cb,
        double * pk_cb_ic)

    int harmonic_pk_nl_at_k_and_z(
        void* pba,
        void * ppm,
        void * phr,
        double k,
        double z,
        double * pk,
        double * pk_cb)

    int harmonic_pk_nl_at_z(
        void * pba,
        void * phr,
        int mode,
        double z,
        double * output_tot,
        double * output_cb_tot)

    int fourier_pk_at_k_and_z(
        void * pba,
        void * ppm,
        void * pfo,
        int pk_output,
        double k,
        double z,
        int index_pk,
        double * out_pk,
        double * out_pk_ic)

    int fourier_pk_tilt_at_k_and_z(
        void * pba,
        void * ppm,
        void * pfo,
        int pk_output,
        double k,
        double z,
        int index_pk,
        double * pk_tilt)

    int fourier_sigmas_at_z(
        void * ppr,
        void * pba,
        void * pfo,
        double R,
        double z,
        int index_pk,
        int sigma_output,
        double * result)

    int fourier_pks_at_kvec_and_zvec(
        void * pba,
        void * pfo,
        int pk_output,
        double * kvec,
        int kvec_size,
        double * zvec,
        int zvec_size,
        double * out_pk,
        double * out_pk_cb)

    int fourier_hmcode_sigma8_at_z(void* pba, void* pfo, double z, double* sigma_8, double* sigma_8_cb)
    int fourier_hmcode_sigmadisp_at_z(void* pba, void* pfo, double z, double* sigma_disp, double* sigma_disp_cb)
    int fourier_hmcode_sigmadisp100_at_z(void* pba, void* pfo, double z, double* sigma_disp_100, double* sigma_disp_100_cb)
    int fourier_hmcode_sigmaprime_at_z(void* pba, void* pfo, double z, double* sigma_prime, double* sigma_prime_cb)
    int fourier_hmcode_window_nfw(void* pfo, double k, double rv, double c, double* window_nfw)

    int fourier_k_nl_at_z(void* pba, void* pfo, double z, double* k_nl, double* k_nl_cb)

    int fourier_B_tree_at_k_and_z(
                                void * pfo,
                                linear_or_logarithmic mode,
                                double f,
                                double b1,
                                double b2,
                                double bG2,
                                double d1,
                                double d2,
                                double d3,
                                double P_eps,
                                double k1,
                                double k2,
                                double k3,
                                double mu1,
                                double mu2,
                                double Pk1,
                                double Pk2,
                                double Pk3,
                                double z,
                                double * bispectrum_treelevel
                                )

    int fourier_B_ell_tree_at_k_and_z(
                                    void *pba,
                                    void * ppm,
                                    void * pfo,
                                    linear_or_logarithmic mode,
                                    int use_IR_resum,
                                    int index_pk,
                                    double b1,
                                    double b2,
                                    double bG2,
                                    double d1,
                                    double d2,
                                    double d3,
                                    double P_eps,
                                    double k1,
                                    double k2,
                                    double k3,
                                    int l,
                                    double z,
                                    double * B_l
                                    )

    int fourier_pk_l_nw_extra_at_kvec_and_z(
                                            void * pba,
                                            void * ppm,
                                            void * pfo,
                                            linear_or_logarithmic mode,
                                            double * ln_kvec,
                                            int kvec_size,
                                            double z,
                                            double * out_pk)

    int harmonic_firstline_and_ic_suffix(void *ppt, int index_ic, char first_line[1024], FileName ic_suffix)

    int harmonic_fast_pk_at_kvec_and_zvec(
                  void * pba,
                  void * phr,
                  double * kvec,
                  int kvec_size,
                  double * zvec,
                  int zvec_size,
                  double * pk_tot_out,
                  double * pk_cb_tot_out,
                  int nonlinear)

    int array_convert_spline_table_columns_to_local_power_basis(
        double * x,
        int x_size,
        double * y_array,
        int y_size,
        double * ddy_array,
        double * coefficients,
        double * breakpoints,
        ErrorMsg errmsg
        )

cdef extern from "header.h":

    cdef enum eft_struct_role:
        eft_master
        eft_slave

    cdef enum eft_tracer:
        eft_matter
        eft_halo
        eft_tracer_num

    cdef enum eft_pk_type:
        pk_lin
        pk_nowiggle
        pk_ir_resummed_lo
        pkmu_rsd_ir_resummed_lo
        pk_ir_resummed_nlo
        pkmu_rsd_ir_resummed_nlo
        pk_type_num

    cdef enum eft_pk_out_type:
        Pdd_mm_real
        Pdd_mm_rsd
        Pdd_hh_real
        Pdd_hh_rsd
        Pdd_hm_real
        Pdd_hm_rsd
        Pdd_mm_real_no_IR_resum
        Pdd_mm_rsd_no_IR_resum
        Pdd_hh_real_no_IR_resum
        Pdd_hh_rsd_no_IR_resum
        Pdd_hm_real_no_IR_resum
        Pdd_hm_rsd_no_IR_resum
        pk_out_type_num

    cdef enum eft_arg_type:
        points
        cartesian_product

    cdef enum eft_spectra_contribution:
        finite_part
        uv_divergence
        ir_divergence
        pole_divergence
        eft_spectra_contribution_num

    int oneloop_nearest_structure_in_time(
        void * peft0,
        int peft_size,
        void * pba,
        void * pfo,
        double z,
        int * index_eft_min_dist,
        void * peft_min_dist,
        ErrorMsg errmsg
        )

    int eft_set_sampling_points_all(
        void * peft0,
        int eft_size,
        double * kvec_Mpc,
        double * muvec,
        int k_size,
        int mu_size
        )

    int eft_set_sampling_points(
        void * peft,
        double * kvec_Mpc,
        double * muvec,
        int k_size,
        int mu_size
        )

    int eft_get_sampling_points(
        void * peft,
        double * kvec_Mpc,
        double * muvec
        )

    int eft_get_sampling_grid_size(
        void * peft,
        int * k_size,
        int * mu_size
        )

    double sigma_sq(
        void * peft,
        short n,
        int index_pk_type
        )

    int oneloop_linear_spectrum_real(
        void * pba,
        void * ppm,
        void * pfo,
        void * peft,
        linear_or_logarithmic mode,
        double * ln_kvec,
        int kvec_size,
        int n_columns,
        double z,
        double f_z,
        double D_z,
        int index_pk_type,
        double * out_pk
        )

    int oneloop_linear_spectrum_rsd(
        void * pba,
        void * ppm,
        void * pfo,
        void * peft,
        linear_or_logarithmic mode,
        double * ln_kvec,
        int kvec_size,
        double * muvec,
        int muvec_size,
        eft_arg_type arg_type,
        double z,
        double f_z,
        double D_z,
        int index_pk_type,
        double * out_pk
        )

    int oneloop_job_powerspectrum_wedges_grid(
        void * peft0,
        int peft_size,
        void * pba,
        void * pfo,
        void * ppm,
        void * ppr,
        eft_pk_out_type pk_out_type,
        double * z,
        void * peft_ip,
        int z_size,
        double * k,
        int k_size,
        double * mu,
        int mu_size,
        double * out_pkmuz
        )

    int oneloop_job_powerspectrum_wedges(
        void * peft0,
        int peft_size,
        void * pba,
        void * pfo,
        void * ppm,
        void * ppr,
        eft_pk_out_type pk_out_type,
        double * zvec,
        double As_correction,
        void * peft_ip,
        int z_size,
        double ** kvec,
        int * k_sizevec,
        double ** muvec,
        int * mu_sizevec,
        double ** out_pkmu,
        double ** ddout_pkmu
        )

    int oneloop_job_powerspectrum_multipoles(
        void * peft0,
        int peft_size,
        void * pba,
        void * pfo,
        void * ppm,
        void * ppr,
        eft_pk_out_type pk_out_type,
        double * zvec,
        double As_correction,
        void * peft_ip,
        int z_size,
        double ** kvec,
        int * k_sizevec,
        double * ap_parallel,
        double * ap_perpendicular,
        double ** out_pkl
        )

    int oneloop_spectra_contributions_output(
        void * peft,
        int pk_type,
        int index_moment,
        double * zout,
        double * muvec_out,
        double * kvec_Mpc_out,
        double ** out_pkmu
        )

cdef extern from "ext_storage.h":
    cdef struct ext_storage:
        double complex *** loop_matrices
        int ** loop_matrices_size
        short ** symmetry
        short ** use_tracer
        short ** spectra_contributions_dimension
        int loop_matrices_stored
        int * eft_index_num
        int eft_size
        double ** period
        ErrorMsg error_message

    int ext_init(void * pext)

    int ext_cleanup(void * pext)

    int ext_save( void * pext,
                  void * pba,
                  void * pth,
                  void * ppt,
                  void * ppm,
                  void * pfo,
                  void * ptr,
                  void * phr,
                  void * ple,
                  void * psd)

    int ext_insert_eft(
                  void * pext,
                  void * peft,
                  const int index_eft,
                  const int num_matrices,
                  ErrorMsg errmsg)

cdef extern from "generalized_triangle_integral.h":
    cdef struct gen_tri_integral:
        ErrorMsg error_message

    int T_master(void *pti,
                 double k12,
                 double k22,
                 double k32,
                 double complex M1,
                 double complex M2,
                 double complex M3,
                 double complex *T_out)

    int B_master(void *pti,
                 double k2,
                 double complex M1,
                 double complex M2,
                 double complex *B_out)

    int Tad_master(void *pti,
                  int n,
                  int d,
                  double complex M,
                  double complex *Tad_out)

    int Tad_var(void *pti,
                int n,
                int d,
                double k2,
                double complex M,
                double complex *I_out)

    int massive_num(void *pti,
                    int n, 
                    int d,
                    double k2,
                    double complex M1, 
                    double complex M2, 
                    double complex *I_out)

    int L_recursion(void *pti,
                    int n1,
                    int d1,
                    int n2,
                    int d2,
                    int n3,
                    int d3,
                    double k12,
                    double k22,
                    double k32,
                    double complex M1,
                    double complex M2,
                    double complex M3,
                    double complex *L_out)

    int T_recursion(void *pti,
                    int d1,
                    int d2,
                    int d3,
                    double k12,
                    double k22,
                    double k32,
                    double complex M1,
                    double complex M2,
                    double complex M3,
                    double complex *T_out)

    int B_recursion(void *pti,
                    int d1,
                    int d2,
                    double k2,
                    double complex M1,
                    double complex M2,
                    double complex *B_out)

    int scalar_prod_one(void *pti,
                        int m, 
                        int n,
                        int d1,
                        int d2,
                        double k22,
                        double complex M1,
                        double complex M2,
                        double complex *I_out)

    int tensor_red_one(void *pti,
                       int n,
                       int d1,
                       int d2,
                       double k12,
                       double k22,
                       double cos12,
                       double complex M1,
                       double complex M2,
                       double complex *I_out);

    int tensor_red_two(void *pti,
                       int n1,
                       int n2,
                       int d,
                       double k12,
                       double k22,
                       double cos12,
                       double complex M,
                       double complex *I_out);

    int util_binomial(void *pti,
                      int n,
                      int k,
                      double *n_over_k)

    int util_antideriv(void *pti,
                       double x,
                       double complex y1,
                       double complex y2,
                       double complex x0,
                       double complex *out)

    int util_prefactor(void *pti,
                       double a,
                       double complex y1,
                       double complex y2,
                       double complex *out)

    int util_F_int(void *pti,
                   double R2,
                   double complex y1,
                   double complex y2,
                   double complex x0,
                   double complex *out)

    int util_Tmaster_contr(void *pti,
                           double y,
                           double k12,
                           double k22,
                           double k32,
                           double complex M1,
                           double complex M2,
                           double complex M3,
                           double complex *out)


