
/** @file utilities.c
 *
 * author: Christian Radermacher, 2023
 *
 * contains the complex Gamma function using the Lanczos approximation
 */


#include "header.h"


int _errno_util = 0;
class_complex _err_last_arg = class_complex(0., 0.);

/**
 * @brief Retrieves the error status of the EFT tools module.
 *
 * @param last_argument   Input: allocated double[2] array for the input arguments that lead to an error
 * @return the error status
*/
int get_error_status(double * last_argument)
{
  last_argument[0] = creal(_err_last_arg);
  last_argument[1] = cimag(_err_last_arg);
  return _errno_util;
}


#define LANCZOS_GAMMA_N 13
/**
 * @brief Complex Gamma function approximated using the
 *        Lanczos method with g ~ 6 and N = 13.
 *
 * @param z     Input: Complex argument
 * @return      value of Gamma(z) / NaN at a singularity
*/
class_complex cGamma(class_complex z)
{
  const double g = 6.0246800407767295;
  const double co[LANCZOS_GAMMA_N] = {2.506628274631000,
      589.5106040667278, -888.0253935502019,
      395.8387847487511, -53.21395931507683,
      1.277182848200161, -4.046172558017667e-4,
      -7.347584327845915e-6, 8.208805790146655e-6,
      -5.159543403041225e-6, 2.319631454949221e-6,
      -6.671246136975432e-7, 9.060393467651553e-8};
  class_complex Lg;
  class_complex t;
  double r;

  r = creal(z);
  if (r < 0.5)
  {
    /** - check the distance to the singularity */
    if (cabs(z - rint(r)) < _EPSILON_) {  /** - uses current rounding-mode, which defaults to round-to-nearest */
      _errno_util |= _ERR_RES_OUT_OF_RANGE_;
      _err_last_arg = z;
      #ifdef NAN
      return class_complex(nan(""), nan(""));   /** Gamma(z) is undefined */
      #else
      return class_complex(0., 0.);
      #endif
    }
    /** - use mirror identity */
    return _PI_ / (csin(_PI_ * z) * cGamma(1.0 - z));
  }

  Lg = co[0];
  for (int k = 1; k < LANCZOS_GAMMA_N; k++)
  { Lg += co[k] / (z + (k - 1.0)); }  /** Partial fraction sum */
  t = z + g - 0.5;
  return cpow(t, z - 0.5) * cexp(-t) * Lg;
}


// Why is this not just tgamma???

/**
 * @brief Real Gamma function approximated using the
 *        Lanczos method with g ~ 6 and N = 13.
 *
 * @param x     Input: Real argument
 * @return      value of Gamma(x) / NaN at a singularity
*/
double rGamma(double x)
{
  const double g = 6.0246800407767295;
  const double co[LANCZOS_GAMMA_N] = {2.506628274631000,
      589.5106040667278, -888.0253935502019,
      395.8387847487511, -53.21395931507683,
      1.277182848200161, -4.046172558017667e-4,
      -7.347584327845915e-6, 8.208805790146655e-6,
      -5.159543403041225e-6, 2.319631454949221e-6,
      -6.671246136975432e-7, 9.060393467651553e-8};
  double Lg;
  double t;

  if (x < 0.5)
  {
    /** - check the distance to the singularity */
    if ((fabs(x - rint(x)) < _EPSILON_)) {  /** - uses current rounding-mode, which defaults to round-to-nearest */
      _errno_util = _ERR_RES_OUT_OF_RANGE_;
      _err_last_arg = x;
      #ifdef NAN
      return nan("");   /** Gamma(z) is undefined */
      #else
      return 0.;
      #endif
    }
    /** - use mirror identity */
    return _PI_ / (sin(_PI_ * x) * rGamma(1 - x));
  }

  Lg = co[0];
  for (int k = 1; k < LANCZOS_GAMMA_N; k++)
  { Lg += co[k] / (x + k - 1); }  /** Partial fraction sum */
  t = x + g - 0.5;
  return pow(t, x - 0.5) * exp(-t) * Lg;
}
