/** @file FFTLog_matter_real.c Documented FFT-Log based 1loop integrals of matter power spectrum in perturbation theory 
 * 
 * Azadeh Moradinezhad Dizgah, June 16th 2021
 *
 * This module performs fast computation of the integrals appearing in the expression of 1loop matter and galaxy power sprtcurm.
 * The computation can be performed either in real or redshift-space. IR-resummation and EFT counter terms are included. 
 * The integrals are computed using FFTLog techniques.
 *
 * The algorithm closely follows Ref. arXiv:1708.08130 by Simonovic et al. After computing the FFT coefficents of matter power spectrum, 
 * sampled in logarithmic scale, the algorithm involves re-casting the integrals into a form that is analytically calculable (Matrices M_ij), 
 * which can be written in terms of ratios of Gamma functions) and finally performing vactor-matrix-vector or matrix-vector multiplications. 
 *
 * An important feauture of fast computation of loop integrals is that all the cosmology-dependance of the loop integrals is captured by the FFT coeffcients
 * which have ~NlogN complexity. The matrices involving the Gamma functions, are computed only once, and at each of MCMC varying cosmological parmaters, 
 * these coeffcients are evaluated for all k-values at once and then vector-matrix-vector multiplication are computed.
 *
 * The FFT coeffcients are computed using FFTW package, while the vector-matrix-vector computations are performed using blas library implemented in gsl.
 * The analytic formulas for M_ij matrices are computed in Mathematica using a modified version of the publicaly available notebook by Simpnovic. 
 *
 * The choice of number of points for FFT decomposition of the power spectrum is very important, in terms of accuracy and execusion time. 
 * For matter loops, if choosing 512 points *which takes ~ 0.45 seconds, I do get subpercent descrepencies between direct numeric integration of the loops.
 * The best value seem to be around 600 points, which takes ~ 0.6 seconds and is in exquisit agreement with direct integration. 
 *
 * In summary, the following functions can be called from other modules:
 * -# P22()
 * -# P13()
 * -# M22()
 * -# M13()
 */

#include "header.h"


/**
 * Compute the P22 contribution to non-linear matter power spectrum given the FFTLog coefficents and frequencies
 *
 * @param fft_struct    Input: structure containing fft coefficents and params
 * @param k             Input: wavenumber in unit of h/Mpc. 
 * @param z             Input: redshift
 * @param cleanup       Input: switch whether to free the M_ij matrix. Only freed at the end of the pipeline
 * @return value of P22 term in unit of (Mpc/h)^3 as a function of redshift and wavnumber
 */

void P_mm_FFTLog(struct fourier *pfo, int index_k, double Plin)
{
      int Nmax = pfo -> fft_ws -> fft_input[real_ir] -> nfft;
      double k = pfo->k[index_k];

      double *np = make_1Darray(1);
      double *p  = make_1Darray(1);

      // Linear cpow Spectrum vector
      double complex *vec_m = make_1D_c_array(Nmax+1);
      vec_fill(pfo -> fft_ws -> fft_input[real_ir], k, MATTER, vec_m);

      // non-propagator calculations
      c_nonprop(vec_m, pfo -> fft_ws -> fft_matrix[real_ir] -> I2200_mat, vec_m, Nmax+1, &np[0]);

      c_dot(vec_m, pfo -> fft_ws -> fft_matrix[real_ir] -> I1300_mat, Nmax+1, &p[0]);

      // adding factored out k and mu dependencies
      pfo -> pk_matter_real_nl -> I2200[index_k] = pow(k, 3.) * np[0];
      pfo -> pk_matter_real_nl -> I1300[index_k] = pow(k, 3.) * Plin * p[0] - 61./630. * Plin * pow(k, 2.) * pfo->fft_ws->sigma_v2;
}

void P_gg_FFTLog(struct fourier *pfo, int index_k, double Plin)
{
      int Nmax = pfo -> fft_ws -> fft_input[real_ir] -> nfft;
      double k = pfo->k[index_k];

      double *np = make_1Darray(6);
      double *p  = make_1Darray(2);
      
      // Linear cpow Spectrum vector for halos
      double complex *vec_h = make_1D_c_array(Nmax+1);
      vec_fill(pfo -> fft_ws -> fft_input[real_ir], k, HALO, vec_h);

      double complex *vec_h_min = make_1D_c_array(Nmax+1);
      vec_fill(pfo -> fft_ws -> fft_input[real_ir], pfo -> fft_ws -> fft_input[real_ir]->kmin_fft_g, HALO, vec_h_min);

      // Linear cpow Spectrum vector for matter
      double complex *vec_m = make_1D_c_array(Nmax+1);
      vec_fill(pfo -> fft_ws -> fft_input[real_ir], k, MATTER, vec_m);
      
      // non-propagator calculations
      c_nonprop(vec_m, pfo -> fft_ws -> fft_matrix[real_ir] -> I2200_mat,           vec_m, Nmax+1, &np[0]);
      c_nonprop(vec_h, pfo -> fft_ws -> fft_matrix[real_ir] -> Idelta200_mat,       vec_h, Nmax+1, &np[1]);
      c_nonprop(vec_h, pfo -> fft_ws -> fft_matrix[real_ir] -> IG200_mat,           vec_h, Nmax+1, &np[2]);
      c_nonprop(vec_h, pfo -> fft_ws -> fft_matrix[real_ir] -> Idelta2delta200_mat, vec_h, Nmax+1, &np[3]);
      c_nonprop(vec_h, pfo -> fft_ws -> fft_matrix[real_ir] -> IG2G200_mat,         vec_h, Nmax+1, &np[4]);
      c_nonprop(vec_h, pfo -> fft_ws -> fft_matrix[real_ir] -> Idelta2G200_mat,     vec_h, Nmax+1, &np[5]);
      
      double Idelta2delta200_const;
      c_nonprop(vec_h_min, pfo -> fft_ws -> fft_matrix[real_ir] -> Idelta2delta200_mat, vec_h_min, Nmax+1, &Idelta2delta200_const);

      // propagator calculations
      c_dot(vec_m, pfo -> fft_ws -> fft_matrix[real_ir] -> I1300_mat, Nmax+1, &p[0]);
      c_dot(vec_h, pfo -> fft_ws -> fft_matrix[real_ir] -> FG200_mat, Nmax+1, &p[1]);
           
      // adding factored out k and mu dependencies
      pfo -> pk_halo_real_nl -> I2200[index_k]           = pow(k, 3.) * np[0];
      pfo -> pk_halo_real_nl -> Idelta200[index_k]       = pow(k, 3.) * np[1];
      pfo -> pk_halo_real_nl -> IG200[index_k]           = pow(k, 3.) * np[2];
      pfo -> pk_halo_real_nl -> Idelta2delta200[index_k] = pow(k, 3.) * np[3] - pow(pfo->fft_ws->fft_input[real_ir]->kmin_fft_g, 3.) * Idelta2delta200_const ;
      pfo -> pk_halo_real_nl -> IG2G200[index_k]         = pow(k, 3.) * np[4];
      pfo -> pk_halo_real_nl -> Idelta2G200[index_k]     = pow(k, 3.) * np[5];

      pfo -> pk_halo_real_nl -> I1300[index_k] = pow(k, 3.) * Plin * p[0] - 61./630. * Plin * pow(k, 2.) * pfo->fft_ws->sigma_v2;
      pfo -> pk_halo_real_nl -> FG200[index_k] = pow(k, 3.) * Plin * p[1];

      pfo -> pk_halo_real_nl -> IR2[index_k] = - 2. * pow(k, 2.) * Plin;
}
