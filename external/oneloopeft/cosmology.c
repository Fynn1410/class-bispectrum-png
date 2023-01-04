/** @file cosmology.c  Documented cosmology module 
 * 
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *   
 * The first routine of this module initalizes the Cosmology structure, which is the main building block of this entire code. 
 * This structure includes two sub-structures: the CLASS cosmology structure and line structure. Once the CLASS cosmology is initialized,
 * various useful functions can be directly called from CLASS, example to compute matter power spectrum and transfer function, 
 * angular and comoving radii, growth factor and growth rate, variance of matter fluctuations and its derivative. 
 * Lastly, the module also includes various window functions and their derivatives. 
 * 
 * In summary, the following functions can be called from other modules:
 * 
 * -# Cosmology_init()              allocates memory to and initializes cosmology structure
 * -# Cosmology_free()              frees the memory allocated to cosmology structure
 * -# CL_Cosmology_initilize()      initializes the class cosmology structure
 * -# CL_Cosmology_free()           frees the class cosmology structure
 * -# PS()                          computes matter power spectrum calling class function
 * -# Transfer()                    computes matter transfer function calling class function
 * -# growth_D()                    computes the scale-dep growth factor
 * -# growth_f()                    computes the scale-dep growth rate dlnD(k,a)/dlna
 * -# scale_indep_growth_D()        computes the scale-indep growth factor using directly CLASS functions
 * -# scale_indep_growth_f()        computes the scale-indep growth rate dlnD(k,a)/dlna  using directly CLASS functions
 * -# Hubble()                      computes hubbble parameter using directly CLASS functions
 * -# angular_distance()            computes angular diamtere distance  using directly CLASS functions
 * -# comoving_radial_distance()    computes radial distance  using directly CLASS functions
 * -# sig_sq()                      computes variance of smoothed matter fluctuations
 * -# der_sig_sq()                  computes derivative of the variance of smoothed matter fluctuations w.r.t. smoothing scale
 * -# sigma0_sq()                   computes variance of unsmoothed matter fluctuations
 * -# rhoc()                        computes the critical density of the universe 
 * -# R_scale()                     computes the size of a spherical halo corresponding to a given mass at z=0
 * -# R_scale_wrong()               computes the size of a spherical halo corresponding to a given mass at a given redshift
 * -# window_rth()                  computes top-hat filter in real space
 * -# window_g()                    computes Gaussian window
 * -# window_kth()                  computes top-hat filter in Fourier space
 * -# derR_window_rth()             computes derivative of top-hat filter in real space w.r.t. smoothing scale 
 * -# derR_logwindow_g()            computes derivative of top-hat filter in Fourier space w.r.t. smoothing scale
 */

 #include "header.h"

/**
 * Compute the matter power spectra (in unit of (Mpc)^3) as a function of k (in unit of 1/Mpc) and z, 
 * Setting the switch "mode", to LINEAR or NONLINEAR, we can compute the linear or nonlinear spectrum respectively.
 * 
 * The CLASS spectra_pk_at_k_and_z() and spectra_pk_nl_at_k_and_z, evaluate the matter power spectrum 
 * at a given value of k and z by interpolating in a table of all P(k)'s computed at this z 
 * by spectra_pk_at_z() (when kmin <= k <= kmax), or eventually by using directly the primordial 
 * spectrum (when 0 <= k < kmin): the latter case is an approximation, valid when kmin << comoving Hubble scale today.
 * Returns zero when k=0. Returns an error when k<0 or k > kmax.
 * 
 * @param Cx                Input: pointer to Cosmology structure
 * @param k                 Input: wavenumbber in unit of 1/Mpc
 * @param z                 Input: redshift to compute the spectrum
 * @param modes             Input: switch to decide whether to compute linear or nonlinear spectrum
 *                              It can  be set to sheth-Tormen (ST), Tinker (TR) or Press-Schecter (PSC)
 * 
 * @return the double value of matter power spectrum
 */
double Pk_dlnPk(struct background * pba,
                struct primordial * ppm,
                struct fourier * pfo, 
                double k, 
                double z, 
                int mode)
{

  double result;
  double pk_cb;
  int index_pk;
  
  if (pfo->has_pk_cb == _TRUE_)
    index_pk = pfo->index_pk_cb;
  else if (pfo->has_pk_m == _TRUE_)
    index_pk = pfo->index_pk_m;
  else
    fprintf(stderr,"Temporary debugging message in Pk_dlnPk of comsology.c");
  
  if (mode == LPOWER){
    fourier_pk_at_k_and_z(pba, ppm, pfo, pk_linear, k, z, index_pk, &pk_cb, NULL);
  }
  else if (mode == NLPOWER){
    fourier_pk_at_k_and_z(pba, ppm, pfo, pk_nonlinear, k, z, index_pk, &pk_cb, NULL);
  }
  
  result = pk_cb;

  return result;
} 


// /**
//  * Compute the growth factor D(k,z) which is scale-indep if mode_nu = NUM, and scale-dep if mode_nu = MASS
//  * The scale-dep growth is calculated by taking the ratio of the transfer function at redshift z and zero.
//  * The scale-indep growth is computed by CLASS directly
//  * The switch "mode" can be set to CDM, BA, TOT to return the growth factor of 
//  * cdm, baryon and total matter. 
//  * 
//  * @param Cx                Input: pointer to Cosmology structure
//  * @param k                 Input: wavenumbber in unit of 1/Mpc
//  * @param z                 Input: redshift to compute the spectrum
//  * @return the growth factor, can be k-dep (ex. with nonzero neutrino mass)
//  */

double growth_D(struct background * pba,
                double z)
{ 

  double tau;
  int last_index; ///junk
  double * pvecback;

  pvecback = (double*) calloc(pba->bg_size,sizeof(double));

  background_tau_of_z(pba,z,&tau);
  background_at_tau(pba,tau,long_info,inter_normal,&last_index,pvecback);

  // DL so perturbations should not be provided
  // class_call(background_tau_of_z(pba,z,&tau),
  //           pba.error_message,Cx->ccs.pt.error_message);
  // class_call(background_at_tau(&Cx -> pba,tau,long_info,inter_normal,&last_index,pvecback),
  //           Cx->ccs.ba.error_message,Cx->ccs.pt.error_message);

  double Dz = pvecback[pba->index_bg_D];

  free(pvecback);

  return Dz;
}



/**
 * Compute the scale-dependant linear growth rate f(k,z) (i.e the velocity growth factor)
 * by taking numerical derivative of the scale_dep_growth_D() function f(k,a) = d ln D(k,a)/d ln a.
 * The switch "mode" can be set to CDM, BA, TOT to return the growth factor of the corresponding matter component.
 *
 * This is a useful function when constraining physics that induces scale-dependant growth 
 * such as massive neutrinos. 
 * 
 * @param Cx                Input: pointer to Cosmology structure
 * @param k                 Input: wavenumbber in unit of 1/Mpc
 * @param z                 Input: redshift to compute the spectrum
 * @return the growth rate, can be k-dep (ex. with nonzero neutrino mass)
 */

double growth_f(struct background * pba, double z)
{ 

  double tau;
  int last_index; ///junk
  double * pvecback;

  pvecback = (double*) calloc(pba->bg_size,sizeof(double));

  background_tau_of_z(pba,z,&tau);
  background_at_tau(pba,tau,long_info,inter_normal,&last_index,pvecback);

  double fz = pvecback[pba->index_bg_f];


  free(pvecback);

  return fz;
}







