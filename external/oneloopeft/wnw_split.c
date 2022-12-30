
/** @file wnw_split.c Documented wiggle-nowiggle split based on 1d Gaussian filter in logarithmic k, 
 * and using the Eisentein-Hu wiggle-nowiggle template arXiv:astro-ph/9709112 
 * 
 * Azadeh Moradinezhad Dizgah, June 16th 2021
 * Christian Radermacher, December 2022
 * 
 * The algorithm closely follows Ref. arXiv:1509.02120 by Vlah et al. (described in Appendix A)
 *
 * The following function will be called from other modules:
 * -# pk_Gfilter_nw()
 */



#include "header.h"


/**
 * @brief Compute the nowiggle component of linear matter power spectrum using 1d logarithmic Gaussian filter.
 * @param pba         Input: pointer to background structure
 * @param ppm         Input: pointer to primordial structure
 * @param pfo         Input: pointer to fourier structure
 * @param k0          Input: fixing scale k in pfo->k, i.e. the largest usable scale
 * @param ln_pknw_array  Output: Dewiggled power spectrum for every pfo->k and pfo->tau
 *                            in units of (Mpc)^3 [size = pfo->k_size_extra * pfo->tau_size].
 *                            Indexed as ln_pknw_array[index_tau * pfo->k_size_extra + index_k]
 * @return the error status
 */
int eft_ln_pk_nw_gfilter(
                      struct background *pba,
                      struct primordial *ppm,
                      struct fourier *pfo,
                      const int index_pk,
                      const int index_k0,
                      const int index_kmin,
                      const int k_size,
                      double *ln_pknw_array) {
  
  //double a = 0.3 * exp(-0.7 * pow(log10(k)+1.5, 2.)) + 0.03;
  //const double smoothing_scale = 0.25 * log(10); // lambda * (1 Mpc)
  int it_k = 0, it_q = 0, it_tau, index_x, index_y, index_ddy, index_num, last_index;
  double ln_pk0_z, k0, smoothing_scale;
  double *pk_approx_f, *intg_splines, *intg_result;

  /** - define indices for the spline array */
  class_define_index(index_x, _TRUE_, it_q, 1);
  class_define_index(index_y, _TRUE_, it_q, 1);
  class_define_index(index_ddy, _TRUE_, it_q, 1);
  index_num = it_q;

  class_alloc(pk_approx_f, pfo->k_size_extra * sizeof(double), pfo->error_message);
  class_alloc(intg_splines, pfo->k_size_extra * index_num * sizeof(double), pfo->error_message);
  class_alloc(intg_result, k_size * sizeof(double), pfo->error_message);

  k0 = pfo->k[index_k0];

  /** - compute the Eisenstein-Hu approximation to the nowiggle power spectrum */
  for (it_q = 0; it_q < pfo->k_size_extra; it_q++) {
    pk_approx_f[it_q] = eft_pk_nw_eisenstein_hu_factor(pba, ppm, pfo, pfo->k[it_q], k0);
    
    /** - also write the array of x-values for splining */
    intg_splines[it_q*index_num + index_x] = pfo->ln_k[it_q];
  }

  
  /** - compute the gaussian window integral at every tau */
  for (it_tau = 0; it_tau < pfo->ln_tau_size; it_tau++)
  {
    /** - power spectrum at fixing scale at tau */
    ln_pk0_z = pfo->ln_pk_l_extra[index_pk][it_tau*pfo->k_size_extra + index_k0];

    /** - compute the integrand at the control points once */
    for (it_q = 0; it_q < pfo->k_size_extra; it_q++)
      intg_splines[it_q*index_num + index_y] = exp( pfo->ln_pk_l_extra[index_pk][it_tau*pfo->k_size_extra + it_q] - ln_pk0_z) \
                                                  / pk_approx_f[it_q];

    /** - spline the integrand function without exponential once */
    class_call(array_spline(intg_splines,
                              index_num,
                              pfo->k_size_extra,
                              index_x,
                              index_y,
                              index_ddy,
                              _SPLINE_EST_DERIV_,
                              pfo->error_message),
                  pfo->error_message,
                  pfo->error_message);

    /** - integrate the same function with different windows */
    for (it_k = 0; it_k < k_size; it_k++)
    {
      /** - compute the running smoothing scale */
      smoothing_scale = 0.691 * exp(-0.132 * pow(pfo->ln_k[index_kmin + it_k] + 3.454, 2.)) + 0.0691;

      /** - integrate the spline */
      class_call(array_integrate_all_spline_gaussian_window(intg_splines,
                                                            index_num,
                                                            pfo->k_size_extra,
                                                            index_x,
                                                            index_y,
                                                            index_ddy,
                                                            pfo->ln_k[index_kmin + it_k],
                                                            smoothing_scale,
                                                            intg_result + it_k,
                                                            pfo->error_message),
                  pfo->error_message,
                  pfo->error_message);

      //fprintf(stderr, "%.15e  %.15e  %.15e  %.15e \n", pfo->k[index_kmin + it_k], intg_splines[(index_kmin+it_k)*index_num + index_y], intg_splines[(index_kmin+it_k)*index_num + index_ddy], intg_result[it_k]);

      /** - multiply with prefactors and take log */
      intg_result[it_k] = log( pk_approx_f[index_kmin + it_k] * intg_result[it_k] ) + ln_pk0_z;
    }

    /** - write to output array */
    memcpy(ln_pknw_array + it_tau*pfo->k_size_extra + index_kmin, 
            intg_result,
            k_size * sizeof(double));
  }

  free(pk_approx_f);
  free(intg_splines);

  return _SUCCESS_;
}


// /**
//  * @brief Compute the nowiggle component of linear matter power spectrum using 1d logarithmic Gaussian filter.
//  * @param pba         Input: pointer to background structure
//  * @param ppm         Input: pointer to primordial structure
//  * @param pfo         Input: pointer to fourier structure
//  * @param k0          Input: fixing scale k in pfo->k, i.e. the largest usable scale
//  * @param ln_pknw_array  Output: Dewiggled power spectrum for every pfo->k and pfo->tau
//  *                            in units of (Mpc)^3 [size = pfo->k_size_extra * pfo->tau_size].
//  *                            Indexed as ln_pknw_array[index_tau * pfo->k_size_extra + index_k]
//  * @return the error status
//  */
// int eft_ln_pk_nw_gfilter_eds(
//                       struct background *pba,
//                       struct primordial *ppm,
//                       struct fourier *pfo,
//                       const int index_k0,
//                       const int index_pk,
//                       double *ln_pknw_array) {
  
//   int it_k, it_tau, last_index;
//   double pk0_z, pk0, growthD, *bg_vec;

//   class_alloc(bg_vec, pba->bg_size * sizeof(double), pfo->error_message);

//   /** - get the time-independent integrals */
//   class_call(eft_pk_nw_gfilter_time_indep(pba, ppm, pfo, index_k0, index_pk,
//                                             ln_pknw_array + (pfo->ln_tau_size-1)*pfo->k_size_extra),
//               pfo->error_message,
//               pfo->error_message);
  
//   /** - power spectrum at fixing scale today */
//   pk0 = pfo->ln_pk_l_extra[index_pk][(pfo->ln_tau_size-1)*pfo->k_size_extra + index_k0];

//   for (it_tau = pfo->ln_tau_size - 2; it_tau >= 0; it_tau--)
//   {
//     /** - get the linear growth factor */
//     class_call(background_at_tau(pba, 
//                                 exp(pfo->ln_tau[it_tau]), 
//                                 long_info, 
//                                 inter_normal, 
//                                 &last_index, 
//                                 bg_vec),
//                 pba->error_message,
//                 pfo->error_message)

//     /** - power spectrum at fixing scale at tau */
//     pk0_z = pfo->ln_pk_l_extra[index_pk][it_tau*pfo->k_size_extra + index_k0];
    
//     /** - copy array to the former time point */
//     memcpy(ln_pknw_array + it_tau*pfo->k_size_extra, 
//             ln_pknw_array + (pfo->ln_tau_size-1)*pfo->k_size_extra, 
//             pfo->k_size_extra * sizeof(double));

//     /** - multiply with the correct time dependence */
//     for (it_k = 0; it_k < pfo->k_size_extra; it_k++)
//     {
//       ln_pknw_array[it_tau*pfo->k_size_extra + it_k] += log( pow(bg_vec[pba->index_bg_D], 2.) * pk0_z / pk0 );
//     }
//   }

//   return _SUCCESS_;
// }



// /**
//  * @brief Compute the nowiggle component of linear matter power spectrum using 1d logarithmic Gaussian filter.
//  * Time dependence (D(z)^2 * Plin(k0, z)/Plin(k0, 0)) has been factorised and needs to be multiplied separately. 
//  * Computing the nowiggle component involves calculating an integral (Eq. A3 of Vlah et al)
//  * 
//  * @param pba         Input: pointer to background structure
//  * @param ppm         Input: pointer to primordial structure
//  * @param pfo         Input: pointer to fourier structure
//  * @param k0          Input: fixing scale k in pfo->k, i.e. the largest usable scale
//  * @param ln_pknw_array  Output: Dewiggled power spectrum for every pfo->k without the time dependence (D(z)^2 * Plin(k0, z)/Plin(k0, 0)) 
//  *                            in units of (Mpc)^3 [size = pfo->k_size_extra]
//  * @return the error status
//  */
// int eft_ln+_pk_nw_gfilter_time_indep(
//                                 struct background *pba, 
//                                 struct primordial *ppm, 
//                                 struct fourier *pfo, 
//                                 const int index_k0,
//                                 const int index_pk,
//                                 double *ln_pknw_array) {

//   const double smoothing_scale = 0.25 * log(10); // lambda * (1 Mpc)
  
//   double *pk_approx_f, *intg_splines, k0;
//   int index_x, index_y, index_ddy;
//   int it_k = 0, it_q = 0, index_num;

//   /** - define indices for the spline array */
//   class_define_index(index_x, _TRUE_, it_q, 1);
//   class_define_index(index_y, _TRUE_, it_q, 1);
//   class_define_index(index_ddy, _TRUE_, it_q, 1);
//   index_num = it_q;

//   class_alloc(pk_approx_f, pfo->k_size_extra * sizeof(double), pfo->error_message);
//   class_alloc(intg_splines, pfo->k_size_extra * index_num * sizeof(double), pfo->error_message);

//   k0 = pfo->k[index_k0];

//   /** - compute the Eisenstein-Hu approximation to the nowiggle power spectrum */
//   for (int it_k = 0; it_k < pfo->k_size_extra; it_k++)
//   {
//     pk_approx_f[it_k] = eft_pk_nw_eisenstein_hu_factor(pba, ppm, pfo, exp(pfo->ln_k[it_k]), k0);
//   }
  
//   /** Prepare the integrand for every ln(k) */
//   for (it_k = 0; it_k < pfo->k_size_extra; it_k++)
//   {
//     for (it_q = 0; it_q < pfo->k_size_extra; it_q++)
//     {
      
//       intg_splines[it_q*index_num + index_x] = pfo->ln_k[it_q];
//       intg_splines[it_q*index_num + index_y] = exp( pfo->ln_pk_l_extra[index_pk][it_q] - pow((pfo->ln_k[it_k] - pfo->ln_k[it_q])/smoothing_scale, 2.)/2. ) / pk_approx_f[it_q];
//     }

//     /** - spline the integrand function */
//     class_call(array_spline(intg_splines,
//                             index_num,
//                             pfo->k_size_extra,
//                             index_x,
//                             index_y,
//                             index_ddy,
//                             _SPLINE_NATURAL_,
//                             pfo->error_message),
//                 pfo->error_message,
//                 pfo->error_message);

//     /** - integrate the spline */
//     class_call(array_integrate_all_spline(intg_splines,
//                                           index_num,
//                                           pfo->k_size_extra,
//                                           index_x,
//                                           index_y,
//                                           index_ddy,
//                                           ln_pknw_array+it_k,
//                                           pfo->error_message),
//                 pfo->error_message,
//                 pfo->error_message);

//     /** - multiply with prefactors */
//     *(ln_pknw_array + it_k) *= pk_approx_f[it_k] / (sqrt(2*_PI_) * smoothing_scale);
//     *(ln_pknw_array + it_k) = log( *(ln_pknw_array + it_k) );
//   }

//   free(pk_approx_f);
//   free(intg_splines);

//   return _SUCCESS_;
// }


/**
 * @brief Compute the Eisentein-Hu approximate nowiggle factor of the linear matter power spectrum
 *
 * @param pba     Input: pointer to background structure
 * @param ppm     Input: pointer to primordial structure
 * @param pfo     Input: pointer to fourier structure
 * @param k       Input: wavenumber in units of 1/Mpc
 * @param k0      Input: smallest value of k, i.e. the largest scale
 * @return Approximate dewiggled power spectrum at z divided by the original linear spectrum
 */
double eft_pk_nw_eisenstein_hu_factor(
                                      struct background *pba, 
                                      struct primordial *ppm, 
                                      struct fourier *pfo, 
                                      const double k, 
                                      const double k0) {
  
    return pow(k/k0, ppm->n_s) * pow(T0(pba, ppm, pfo, k)/T(pba, ppm, pfo, k0), 2.);
}


/**
 * @brief Compute the Eisentein-Hu approximate wiggle factor of the linear matter power spectrum
 *
 * @param pba     Input: pointer to background structure
 * @param ppm     Input: pointer to primordial structure
 * @param pfo     Input: pointer to fourier structure
 * @param k       Input: wavenumber in units of 1/Mpc
 * @param k0      Input: smallest value of k, i.e. the largest scale
 * @return Approximate wiggle power spectrum at z divided by the original linear spectrum
 */
double eft_pk_w_eisenstein_hu_factor(
                                      struct background *pba, 
                                      struct primordial *ppm, 
                                      struct fourier *pfo, 
                                      const double k, 
                                      const double k0) {
  
    return pow(k/k0, ppm->n_s) * pow(T(pba, ppm, pfo, k)/T(pba, ppm, pfo, k0), 2.);
}


/**
 * Compute the nowiggle component of linear matter power spectrum using 3d (NO!) Gaussian filter
 * Computing the nowiggle component involves calculating an integral (Eq. A3 of Vlah et al)
 * Below, pk_nw_integrand()is the corresponding integrand and pk_Gfilter_nw() is the integrator 
 * 
 * @param Cx            Input: pointer to Cosmology structure
 * @param k             Input: wavenumber in unit of 1/Mpc
 * @param k0            Input: smallest value of k, i.e. the largest scale
 * @return broadband component in unit of (Mpc)^3
 */
double pk_Gfilter_nw(struct background *pba, struct primordial *ppm, struct fourier *pfo, double k, double k0, double z)    
{ 
    struct integrand_parameters2 par; 
    par.ppm = ppm;
    par.pba = pba;
    par.pfo = pfo;
    par.p4  = k;
    par.p5  = k0;

    double a       = 0.3 * exp(-0.7 * pow(log10(k) + 1.5, 2.)) + 0.03; 
    double logqmin =  log10(k) - 4.* a;
    double logqmax =  log10(k) + 4.* a;

    double result, error;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);

    gsl_function F;
    F.function = &pk_nw_integrand;
    F.params = &par;
    gsl_integration_qags(&F,logqmin, logqmax,0.0,1.0e-3,1000000,w,&result,&error);
    gsl_integration_workspace_free (w);
   
    double pk0     = Pk_dlnPk(pba, ppm, pfo, k0, z, LPOWER);
    double ratio   = 1./(sqrt(2.*M_PI)* a) * result; 
    double out     =  EH_PS_nw(pba, ppm, pfo, k, k0, pk0) * ratio;

    // printf("%12.6e %12.6e %12.6e \n", k, out, Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER));

    return out;
}


/**
 * Integrand to compute the nowiggle matter power spectrum
 * 
 * @param x            Input: integration variable, k-values
 * @param par          Input: integration parameters
 * @return integrand to be used in pk_Gfilter_nw() function         
 */

double pk_nw_integrand(double x, void *par)  ///integration variable x = logq
{  
    struct integrand_parameters2 pij;
    pij = *((struct integrand_parameters2 *)par);

    struct background *pba = pij.pba;
    struct primordial *ppm = pij.ppm;
    struct fourier *pfo    = pij.pfo;
    double k               = pij.p4;
    double kf0             = pij.p5;

    double logq  = x;
    double q     = pow(10.,logq);
    double logk  = log10(k);
    double a     = 0.3 * exp(-0.7 * pow(logk + 1.5, 2.)) + 0.03; 
    double pkf0  = Pk_dlnPk(pba, ppm, pfo, kf0, 0., LPOWER);
    double ratio = Pk_dlnPk(pba, ppm, pfo, q, 0., LPOWER)/EH_PS_nw(pba, ppm, pfo, q, kf0, pkf0);  
    double out   = ratio * exp(-1./(2.*pow(a,2.))*pow(logk-logq,2.));

    return out;
}        


/**
 * Compute the Eisentein-Hu approximate wiggle component of linear matter power spectrum
 *
 * @param Cx            Input: pointer to Cosmology structure
 * @param k             Input: wavenumber in unit of 1/Mpc
 * @param k0            Input: smallest value of k, i.e. the largest scale
 * @param pk0           Input: value of the power spectrum at the largest scale
 * @return P_w(k) in unit of (Mpc)^3
 */

double EH_PS_w(struct background *pba, struct primordial *ppm, struct fourier *pfo, double k, double k0, double pk0)
{
    double h    = pba->h;
    double ns   = ppm->n_s;
    double kh   = h*k;
    double kh0  = h*k0;
    // double norm = pk0/(pow(k0,ns)*pow(T(kh0,Cx),2.));
    // double out  = pow(k,ns) * pow(T(kh,Cx),2.) * norm;
    double norm = pk0/(pow(k0,ns)*pow(T(pba, ppm, pfo, k0),2.));
    double out  = pow(k,ns) * pow(T(pba, ppm, pfo, k),2.) * norm;

    return out;
}    


/**
 * Compute the Eisentein-Hu approximate nowiggle component of the linear matter power spectrum
 *
 * @param Cx            Input: pointer to Cosmology structure
 * @param k           Input: wavenumber in unit of 1/Mpc
 * @param k0          Input: smallest value of k, i.e. the largest scale
 * @param pk0         Input: value of the power spectrum at the largest scale
 * @return P_nw(k) in unit of (Mpc)^3
 */

double EH_PS_nw(struct background *pba, struct primordial *ppm, struct fourier *pfo, double k,double k0,double pk0) 
{
    double h    = pba->h;
    double ns   = ppm->n_s;
    double kh   = h*k;
    double kh0  = h*k0;
    // double norm = p0/(pow(k0,ns)*pow(T(kh0,Cx),2.));
    // double out  = pow(k,ns) * pow(T0(kh,Cx),2.) * norm;
    double norm = pk0/(pow(k0,ns)*pow(T(pba, ppm, pfo, k0),2.));
    double out  = pow(k,ns) * pow(T0(pba, ppm, pfo, k),2.) * norm;

    return out;
}


/**
 * Compute the no-baryon transfer function given in Eq. 29 of EH ref
 *
 * @param pba     Input: pointer to background structure
 * @param ppm     Input: pointer to primordial structure
 * @param pfo     Input: pointer to fourier structure
 * @param k           Input: wavenumber in unit of 1/Mpc
 * @return value of nor-baryon transfer fit
 */
double T0(
          struct background *pba, 
          struct primordial *ppm, 
          struct fourier *pfo, 
          double k) {

  double h     = pba->h;
  double ombh2 = pow(h,2.) * pba->Omega0_b; 
  double omch2 = pow(h,2.) * pba->Omega0_cdm; 
  double om0h2 = ombh2 + omch2;
  double om0   = om0h2/pow(h,2.);
  double theta = 2.728/2.7;    //OBBE-FIRAS value 

  double s     = 44.5 * log(9.83/om0h2)/sqrt(1.+10.*pow(ombh2,3./4.));    ///approximate sound speed given in Eq. (26) of EH
  double AG    = 1. - 0.328*log(431.*om0h2)*ombh2/om0h2 + 0.38*log(22.3*om0h2)*pow(ombh2/om0h2,2.);                                            
  double Gamma = om0 * h * (AG + (1. - AG)/(1.+pow(0.43*k*s,4.)));
  double q     = k/h *pow(theta,2.)/Gamma ;
  double L0    = log(2.*exp(1.) + 1.8 * q);
  double C0    = 14.2 + 731./(1.+62.5*q);
  double out   = L0/(L0+C0*pow(q,2.));

  return out;
}

/** 
 * Compute the total baryon+CDM transfer function
 *
 * @param pba     Input: pointer to background structure
 * @param ppm     Input: pointer to primordial structure
 * @param pfo     Input: pointer to fourier structure
 * @param k       Input: wavenumber in units of 1/Mpc
 * @return value of baryon+cdm transfer function
 */
double T(
          struct background *pba, 
          struct primordial *ppm, 
          struct fourier *pfo, 
          double k) {

  double h     = pba->h;
  double ombh2 = pow(h,2.) * pba->Omega0_b;
  double omch2 = pow(h,2.) * pba->Omega0_cdm;
  double om0h2 = ombh2 + omch2;
  double om0   = om0h2/pow(h,2.);
  double theta = 2.728/2.7;    //OBBE-FIRAS value 

  double HH0     = 1.e3*1.e2*h/299792458.;   ///H0 value devided by the speed of light
  double zeq     = 2.5e4 * om0h2 * pow(theta,-4.);
  double keq     = sqrt(2.*om0*pow(HH0,2.)*zeq);
  double k_silk  = 1.6*pow(ombh2,0.52)*pow(om0h2,0.73)*(1.+pow(10.4*om0h2,-0.95));  ////in 1/Mpc
  
  double beta_node = 8.41*pow(om0h2,0.435);
  double s         = 44.5 * log(9.83/om0h2)/sqrt(1.+10.*pow(ombh2,3./4.));    ///approximate sound speed given in Eq. (26) of EH
  double st        = s/pow(1.+ pow(beta_node/(k*s),3.),1./3.);

  double bb1   = 0.313 * pow(om0h2,-0.419) * (1.+0.607*pow(om0h2,0.674));
  double bb2   = 0.238 * pow(om0h2,0.223);
  double zd    = 1291. * pow(om0h2,0.251)/(1.+0.659*pow(om0h2,0.828))*(1.+bb1*pow(ombh2,bb2));
  double Rd    = 31.5 * ombh2*pow(theta,-4.)*pow(zd/1.e3,-1.);
  double y     = (1.+zeq)/(1.+zd);
  double G     = y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)));

  double a1     = pow(46.9*om0h2,0.670)*(1.+pow(32.1*om0h2,-0.532));
  double a2     = pow(12.0*om0h2,0.424)*(1.+pow(45.*om0h2,-0.582));
  double alphac = pow(a1,-ombh2/om0h2) * pow(a2,-pow(ombh2/om0h2,3.));
  double b1     = 0.944*pow(1.+pow(458.*om0h2,-0.708),-1.);
  double b2     = pow(0.395 * om0h2,-0.0266);
  double betac  = pow(1. + b1*(pow(omch2/om0h2,b2)-1.),-1.);

  double alphab = 2.07*keq*s*pow(1.+Rd,-3./4.)*G;
  double betab  = 0.5 + ombh2/om0h2 + (3.-2.*ombh2/om0h2) * sqrt(pow(17.2*om0h2,2.)+1.); 
  double f      = 1./(1.+pow(k*s/5.4, 4.));

  double Tb = (Tt0(pba,ppm,pfo,k,1.,1.)/(1.+pow(k*s/5.2,2.)) + alphab/(1.+pow(betab/(k*s),3.)) * exp(-pow(k/k_silk,1.4)))* sin(k*st) / (k*st); ///Eq. 21 of EH ref.
  double Tc = f*Tt0(pba,ppm,pfo,k,1.,betac) + (1.-f) *Tt0(pba,ppm,pfo,k,alphac,betac);       ///Eq. 17 of EH ref                                    

  double out = ombh2/om0h2 * Tb + omch2/om0h2 * Tc;

  return out;
}    


/**
 * Compute the function defined in Eq. 19 of EH ref, which will be used to compute the fit for CDM transfer function in Eq. 17.
 *
 * @param pba     Input: pointer to background structure
 * @param ppm     Input: pointer to primordial structure
 * @param pfo     Input: pointer to fourier structure
 * @param k       Input: wavenumber in unit of 1/Mpc. 
 * @param x1      Input: alpha_c
 * @param x2      Input: beta_c
 * @return value of the function
 */
double Tt0(
            struct background *pba, 
            struct primordial *ppm, 
            struct fourier *pfo, 
            double k, 
            double x1, 
            double x2) {

  double h     = pba->h;
  double ombh2 = pow(h,2.) * pba->Omega0_b;
  double omch2 = pow(h,2.) * pba->Omega0_cdm;
  double om0h2 = ombh2 + omch2;
  double theta = 2.728/2.7;    //OBBE-FIRAS value 

  double qq  = k*pow(theta,2.)*pow(om0h2,-1.);
  double C   = 14.2/x1+386./(1.+69.9*pow(qq,1.08));
  double out = log(exp(1.)+1.8*x2*qq)/(log(exp(1.)+1.8*x2*qq)+C*pow(qq,2.));

  return out;
}     

