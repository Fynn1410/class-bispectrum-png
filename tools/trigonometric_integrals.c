/**
 * Module with tools for trigonometric integrals
 * Samuel Brieden, 2018
 */

#include "trigonometric_integrals.h"

/** this is the Cosine Integral function Ci(x) */
int cosine_integral(
                    double x,
                    double *Ci,
                    ErrorMsg error_message
                    ){

  double x2, y, f, g, ci8;
  double em_const = 0.577215664901532861e0;

  if (fabs(x)<=4.){
    x2=x*x;

    ci8=em_const+log(x)+x2*(-0.25e0+x2*(7.51851524438898291e-3+x2*(-1.27528342240267686e-4
            +x2*(1.05297363846239184e-6+x2*(-4.68889508144848019e-9+x2*(1.06480802891189243e-11
            +x2*(-9.93728488857585407e-15)))))))/ (1.+x2*(1.1592605689110735e-2+
            x2*(6.72126800814254432e-5+x2*(2.55533277086129636e-7+x2*(6.97071295760958946e-10+
            x2*(1.38536352772778619e-12+x2*(1.89106054713059759e-15+x2*(1.39759616731376855e-18))))))));

    *Ci=ci8;
  }
  else {
    y=1./(x*x);

    f = (1.e0 + y*(7.44437068161936700618e2 + y*(1.96396372895146869801e5 +
            y*(2.37750310125431834034e7 +y*(1.43073403821274636888e9 + y*(4.33736238870432522765e10
            + y*(6.40533830574022022911e11 + y*(4.20968180571076940208e12 + y*(1.00795182980368574617e13
            + y*(4.94816688199951963482e12 +y*(-4.94701168645415959931e11)))))))))))/
            (x*(1. +y*(7.46437068161927678031e2 +y*(1.97865247031583951450e5 +
            y*(2.41535670165126845144e7 + y*(1.47478952192985464958e9 +
            y*(4.58595115847765779830e10 +y*(7.08501308149515401563e11 + y*(5.06084464593475076774e12
            + y*(1.43468549171581016479e13 + y*(1.11535493509914254097e13)))))))))));

    g = y*(1.e0 + y*(8.1359520115168615e2 + y*(2.35239181626478200e5 + y*(3.12557570795778731e7
            + y*(2.06297595146763354e9 + y*(6.83052205423625007e10 +
            y*(1.09049528450362786e12 + y*(7.57664583257834349e12 +
            y*(1.81004487464664575e13 + y*(6.43291613143049485e12 +y*(-1.36517137670871689e12)))))))))))
            / (1. + y*(8.19595201151451564e2 +y*(2.40036752835578777e5 +
            y*(3.26026661647090822e7 + y*(2.23355543278099360e9 + y*(7.87465017341829930e10
            + y*(1.39866710696414565e12 + y*(1.17164723371736605e13 + y*(4.01839087307656620e13 +y*(3.99653257887490811e13))))))))));
    *Ci=f*sin(x)-g*cos(x);
  }
  return _SUCCESS_;
}

/** this is the Sine Integral function Si(x) */
int sine_integral(
                  double x,
                  double *Si,
                  ErrorMsg error_message
                  ){

  double x2, y, f, g, si8;
  double pi8=3.1415926535897932384626433;

  if (fabs(x)<=4.){
    x2=x*x;

    si8 = x*(1.e0+x2*(-4.54393409816329991e-2+x2*(1.15457225751016682e-3
            +x2*(-1.41018536821330254e-5+x2*(9.43280809438713025e-8+x2*(-3.53201978997168357e-10
            +x2*(7.08240282274875911e-13+x2*(-6.05338212010422477e-16))))))))/
            (1.+x2*(1.01162145739225565e-2 +x2*(4.99175116169755106e-5+
            x2*(1.55654986308745614e-7+x2*(3.28067571055789734e-10+x2*(4.5049097575386581e-13
            +x2*(3.21107051193712168e-16)))))));

    *Si=si8;
  }
  else {
    y=1./(x*x);

    f = (1.e0 + y*(7.44437068161936700618e2 + y*(1.96396372895146869801e5 +
            y*(2.37750310125431834034e7 +y*(1.43073403821274636888e9 + y*(4.33736238870432522765e10
            + y*(6.40533830574022022911e11 + y*(4.20968180571076940208e12 +
            y*(1.00795182980368574617e13 + y*(4.94816688199951963482e12 +
            y*(-4.94701168645415959931e11)))))))))))/ (x*(1. +y*(7.46437068161927678031e2 +
            y*(1.97865247031583951450e5 +y*(2.41535670165126845144e7 +
            y*(1.47478952192985464958e9 + y*(4.58595115847765779830e10 +
            y*(7.08501308149515401563e11 + y*(5.06084464593475076774e12 +
            y*(1.43468549171581016479e13 + y*(1.11535493509914254097e13)))))))))));


    g = y*(1.e0 + y*(8.1359520115168615e2 + y*(2.35239181626478200e5 +
            y*(3.12557570795778731e7 + y*(2.06297595146763354e9 + y*(6.83052205423625007e10 +
            y*(1.09049528450362786e12 + y*(7.57664583257834349e12 +y*(1.81004487464664575e13 +
            y*(6.43291613143049485e12 +y*(-1.36517137670871689e12)))))))))))/
            (1. + y*(8.19595201151451564e2 +y*(2.40036752835578777e5 + y*(3.26026661647090822e7
            + y*(2.23355543278099360e9 + y*(7.87465017341829930e10 + y*(1.39866710696414565e12
            + y*(1.17164723371736605e13 + y*(4.01839087307656620e13 +y*(3.99653257887490811e13))))))))));

    *Si=pi8/2.-f*cos(x)-g*sin(x);
  }
  return _SUCCESS_;
}


/**
 * @brief Computes the spherical Bessel functions of the first kind up to order n
 * 
 * @param n       Input: maximum order
 * @param x       Input: argument
 * @param result  Output: array of output values j_{n}(x), must be initialised previously
 * @param errmsg
 * 
 * @returns the error status
*/
int spherical_bessel_j(const int n, 
                      const double x,
                      double * result,
                      ErrorMsg errmsg) {
  int m;
  const double coeff_border_ab[4] = {9.9123380e-1, 2.1567720e-2, -3.1745158e-5, -3.6693021e-8};
  const double coeff_border_bc[4] = {3.7588265, 1.0238752e-2, 8.7822303e-6, -9.9921351e-8};
  const double ln_x = log(x);

  if (ln_x < array_polynomial(coeff_border_ab, 4, 0, (double)n)) {
    /** - series expansion */
    int k, ksup;
    double prefactor = 1., term, sum;
    double x2_half = x*x/2;

    /** - write the prefactors to the series */
    *(result + 0) = 1.;
    for (m = 1; m <= n; m++) {
      prefactor *= x/(2*m+1);
      *(result + m) = prefactor;
    }

    /** - compute the expansion for every order */
    /** m = 0 */
    term = sum = 1.;
    for (k = 1; fabs(term) > __DBL_EPSILON__; k++) {
      term *= -x2_half / (k * (2*k + 1));
      sum += term;
    }
    *(result + 0) *= sum;
    ksup = k;
    /** m > 0 with the same number of terms */
    for (m = 1; m <= n; m++) {
      term = sum = 1.;
      for (k = 1; k < ksup; k++) {
        term *= -x2_half / (k * (2*k + 2*m + 1));
        sum += term;
      }
    
      *(result + m) *= sum;
    }
  }
  else if (ln_x < array_polynomial(coeff_border_bc, 4, 0, (double)n)) {
    /** - backward recursion */
    const int nmax = n + 100;
    double jn_fid[nmax+2], factor;
    
    /** - starting values set at nmax */
    jn_fid[nmax+1] = 0.; jn_fid[nmax] = 1.;

    for (m = nmax; m > 0; m--) {
      jn_fid[m-1] = (2*m + 1) * jn_fid[m] / x - jn_fid[m+1];
    }
    /** - proportionality factor */
    factor = sin(x)/x / jn_fid[0];
    for (m = 0; m <= n; m++) {
      *(result + m) = factor * jn_fid[m];
    }
  }
  else {
    /** - forward recursion */
    for (m = 0; m <= n && m <= 0; m++) {
      /** j_0(x) */
      *(result + m) = sin(x) / x;
    }
    for (m; m <= n && m <= 1; m++) {
      /** j_1(x) */
      *(result + m) = (sin(x)/x - cos(x))/x;
    }
    for (m; m <= n; m++) {
      /** j_n(x) = (2n-1)/x * j_{n-1}(x) - j_{n-2}(x) */
      *(result + m) = (2*m - 1) * *(result + m-1) / x - *(result + m-2);
    }
  }
  
  return _SUCCESS_;
}