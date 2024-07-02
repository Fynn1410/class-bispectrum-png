/** @file kernel_matrices.c
 *
 * author: Christian Radermacher, 2023
 * based on prototype version of Azadeh Moradinezhad Dizgah & Dennis Linde
 *
 * contains analytic FFTLog kernel definitions using dimensional regularisation with d = 3
*/

#include "kernel_matrices.h"

/*
      --------------- General functions -------------------
*/

/**
 * @brief Analytic expression of the integrals of the form
 *        J(nu1,nu2) = k^(-3 + 2 nu12) * int_q  1/( q^(2 nu1) |k - q|^(2 nu2) ).
 *        Expression is symmetric
 *
 * @param nu1         Input: power of q
 * @param nu2         Input: power of |k - q|
 *
 * @return value of I(nu1,nu2)
 */
static double complex eft_mat_J(const double complex nu1, const double complex nu2)
{
  double complex nu12 = nu1 + nu2;
  double complex numerator    = cGamma(3./2.-nu1) * cGamma(3./2.-nu2) * cGamma(nu12-3./2.);
  double complex denominator  = cGamma(nu1) * cGamma(nu2) * cGamma(3.-nu12);
  double complex out          = J_PREFACTOR * numerator/denominator;

  return out;
}

/**
 * @brief Analytic expression of the integrals of the form
 *        M1(nu1,nu2) = k^(-3 + 2 nu12) * int_q  q|| /( q^(2 nu1) |k-q|^(2 nu2) ).
 *
 * @param nu1         Input: power of q^2
 * @param nu2         Input: power of |k - q|^2
 *
 * @return value of M1(nu1,nu2)
 */
static double complex eft_mat_M1(const double complex nu1, const double complex nu2)
{
  double complex out = 0.5*(eft_mat_J(-1. + nu1,nu2) - eft_mat_J(nu1,-1. + nu2) + eft_mat_J(nu1,nu2));

  return out;
}

/**
 * Analytic expression of the integrals of the form
 * M2 = k^(-3 + 2 nu12) * int_q  q|| q|| /(q^(2 nu1) kmq^(2 nu2)
 *
 * @param nu1   Input: power of q^2
 * @param nu2   Input: power of |k - q|^2
 * @param mu    Input: cos( angle between k and z )
 *
 * @return value of M2(nu1,nu2,mu)
 */
static double complex eft_mat_M2(const double complex nu1, const double complex nu2, const double mu)
{
  double complex A2 = 0.125*(-eft_mat_J(-2 + nu1,nu2) + 2*eft_mat_J(-1 + nu1,-1 + nu2) +
                            2*eft_mat_J(-1 + nu1,nu2) - eft_mat_J(nu1,-2 + nu2) + 2*eft_mat_J(nu1,-1 + nu2) -
                            eft_mat_J(nu1,nu2));
  double complex B2 = 0.125*(3*eft_mat_J(-2 + nu1,nu2) - 6*eft_mat_J(-1 + nu1,-1 + nu2) +
                            2*eft_mat_J(-1 + nu1,nu2) + 3*eft_mat_J(nu1,-2 + nu2) -
                            6*eft_mat_J(nu1,-1 + nu2) + 3*eft_mat_J(nu1,nu2));

  double complex out = A2 + mu*mu * B2;

  return out;
}

/**
 * Analytic expression of the integrals of the form
 * M3 = k^(-3 + 2 nu12) * int_q  q|| q|| q|| /(q^(2 nu1) kmq^(2 nu2)
 *
 * @param nu1   Input: power of q^2
 * @param nu2   Input: power of |k - q|^2
 * @param mu    Input: cos( angle between k and z )
 *
 * @return value of M3(nu1,nu2,mu)
 */
static double complex eft_mat_M3(const double complex nu1, const double complex nu2, const double mu)
{
  double complex A3 = -0.1875*(eft_mat_J(-3 + nu1,nu2) - 3*eft_mat_J(-2 + nu1,-1 + nu2) -
                              eft_mat_J(-2 + nu1,nu2) + 3*eft_mat_J(-1 + nu1,-2 + nu2) -
                              2*eft_mat_J(-1 + nu1,-1 + nu2) - eft_mat_J(-1 + nu1,nu2) -
                              eft_mat_J(nu1,-3 + nu2) + 3*eft_mat_J(nu1,-2 + nu2) -
                              3*eft_mat_J(nu1,-1 + nu2) + eft_mat_J(nu1,nu2));
  double complex B3 = 0.0625*(5*eft_mat_J(-3 + nu1,nu2) - 15*eft_mat_J(-2 + nu1,-1 + nu2) +
                              3*eft_mat_J(-2 + nu1,nu2) + 15*eft_mat_J(-1 + nu1,-2 + nu2) -
                              18*eft_mat_J(-1 + nu1,-1 + nu2) + 3*eft_mat_J(-1 + nu1,nu2) -
                              5*eft_mat_J(nu1,-3 + nu2) + 15*eft_mat_J(nu1,-2 + nu2) -
                              15*eft_mat_J(nu1,-1 + nu2) + 5*eft_mat_J(nu1,nu2));

  double complex out = A3 + mu*mu * B3;

  return out;
}

/**
 * Analytic expression of the integrals of the form
 * M4 = k^(-3 + 2 nu12) * int_q q|| q|| q|| q|| /(q^(2 nu1) kmq^(2 nu2)
 *
 * @param nu1   Input: power of q^2
 * @param nu2   Input: power of |k - q|^2
 * @param mu    Input: cos( angle between k and z )
 *
 * @return value of M4(nu1,nu2,mu)
 */
static double complex eft_mat_M4(const double complex nu1, const double complex nu2, const double mu)
{
  double complex A4 = (3/128.)*(eft_mat_J(-4 + nu1,nu2) - 4*eft_mat_J(-3 + nu1,-1 + nu2) - 4*eft_mat_J(-3 + nu1,nu2) +
                                6*eft_mat_J(-2 + nu1,-2 + nu2) + 4*eft_mat_J(-2 + nu1,-1 + nu2) +
                                6*eft_mat_J(-2 + nu1,nu2) - 4*eft_mat_J(-1 + nu1,-3 + nu2) +
                                4*eft_mat_J(-1 + nu1,-2 + nu2) + 4*eft_mat_J(-1 + nu1,-1 + nu2) -
                                4*eft_mat_J(-1 + nu1,nu2) + eft_mat_J(nu1,-4 + nu2) - 4*eft_mat_J(nu1,-3 + nu2) +
                                6*eft_mat_J(nu1,-2 + nu2) - 4*eft_mat_J(nu1,-1 + nu2) + eft_mat_J(nu1,nu2));
  double complex B4 = (-3/64.)*(5*eft_mat_J(-4 + nu1,nu2) - 20*eft_mat_J(-3 + nu1,-1 + nu2) - 4*eft_mat_J(-3 + nu1,nu2) +
                                30*eft_mat_J(-2 + nu1,-2 + nu2) - 12*eft_mat_J(-2 + nu1,-1 + nu2) -
                                2*eft_mat_J(-2 + nu1,nu2) - 20*eft_mat_J(-1 + nu1,-3 + nu2) +
                                36*eft_mat_J(-1 + nu1,-2 + nu2) - 12*eft_mat_J(-1 + nu1,-1 + nu2) -
                                4*eft_mat_J(-1 + nu1,nu2) + 5*eft_mat_J(nu1,-4 + nu2) - 20*eft_mat_J(nu1,-3 + nu2) +
                                30*eft_mat_J(nu1,-2 + nu2) - 20*eft_mat_J(nu1,-1 + nu2) + 5*eft_mat_J(nu1,nu2));
  double complex C4 = (1/128.)*(35*eft_mat_J(-4 + nu1,nu2) - 140*eft_mat_J(-3 + nu1,-1 + nu2) + 20*eft_mat_J(-3 + nu1,nu2) +
                                210*eft_mat_J(-2 + nu1,-2 + nu2) - 180*eft_mat_J(-2 + nu1,-1 + nu2) +
                                18*eft_mat_J(-2 + nu1,nu2) - 140*eft_mat_J(-1 + nu1,-3 + nu2) +
                                300*eft_mat_J(-1 + nu1,-2 + nu2) - 180*eft_mat_J(-1 + nu1,-1 + nu2) +
                                20*eft_mat_J(-1 + nu1,nu2) + 35*eft_mat_J(nu1,-4 + nu2) - 140*eft_mat_J(nu1,-3 + nu2) +
                                210*eft_mat_J(nu1,-2 + nu2) - 140*eft_mat_J(nu1,-1 + nu2) + 35*eft_mat_J(nu1,nu2));
  double mu_sq = mu*mu;
  double complex out = A4 + mu_sq * (B4 + mu_sq * C4);

  return out;
}

/*
      _________________ 0-th Moment _______________________
 */

static double complex eft_mat_I2200(const double complex * const n)
{
  double complex numerator   = ((-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(58. + 98.*n[0]*n[0]*n[0]*n[1] + (3. - 91.*n[1])*n[1] +
    7.*n[0]*n[0]*(-13. - 2.*n[1] + 28.*n[1]*n[1]) + n[0]*(3. + 2.*n[1]*(-73. + 7.*n[1]*(-1. + 7.*n[1])))));
  double complex denominator = (392.*n[0]*(1. + n[0])*(-1. + 2.*n[0])*n[1]*(1. + n[1])*(-1. + 2.*n[1]));
  double complex out         = numerator/denominator * eft_mat_J(n[0], n[1]);

  return out;
}

static double complex eft_mat_I1300(const double complex * const n)
{
  double complex numerator   = ((1. + 9.*n[0])*ctan(n[0]*_PI_));
  double complex denominator = (672.*n[0]*(n[0] + 1.)*(n[0] - 1.)*(n[0] - 2.)*(n[0] - 3.)*_PI_);
  double complex out         = numerator/denominator;

  return out;
}

static double complex eft_mat_Idelta200(const double complex * const n)
{
  double complex numerator   = ((-3. + 2.*n[0] + 2.*n[1])*(-4. + 7.*n[0] + 7.*n[1]));
  double complex denominator = (28.*n[0]*n[1]);
  double complex out         = numerator/denominator * eft_mat_J(n[0], n[1]);

  return out;
}

static double complex eft_mat_IG200(const double complex * const n)
{
  double complex numerator   = ((-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(6. + 7.*n[0] + 7.*n[1]));
  double complex denominator = (56.*n[0]*(1. + n[0])*n[1]*(1. + n[1]));
  double complex out         = - numerator/denominator * eft_mat_J(n[0], n[1]);

  return out;
}

static double complex eft_mat_Idelta2delta200(const double complex * const n)
{
  double complex out         = eft_mat_J(n[0], n[1]);
  return out;
}

static double complex eft_mat_IG2G200(const double complex * const n)
{
  double complex numerator   = ((-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1]));
  double complex denominator = (2.*n[0]*(1. + n[0])*n[1]*(1. + n[1]));
  double complex out         = numerator/denominator * eft_mat_J(n[0], n[1]);

  return out;
}

static double complex eft_mat_Idelta2G200(const double complex * const n)
{
  double complex numerator   = (3. - 2.*n[0] - 2.*n[1]);
  double complex denominator = (2.*n[0]*n[1]);
  double complex out         = numerator/denominator * eft_mat_J(n[0], n[1]);

  return out;
}

static double complex eft_mat_FG200(const double complex * const n)
{
  double complex numerator    = (-15.*ctan(n[0]*_PI_));
  double complex denominator  = (112.*n[0]*(n[0] + 1.)*(n[0] - 1.)*(n[0] - 2.)*(n[0] - 3.)*_PI_);
  double complex out          = numerator/denominator;

  return out;
}

/*
      _________________ 1-st Moment _______________________
 */

static double complex eft_mat_I2201(const double complex * const n)
{
  double complex numerator   = ((-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(46. + 98.*n[0]*n[0]*n[0]*n[1] + (13. - 63.*n[1])*n[1] +
    7.*n[0]*n[0]*(-9. + 2.*n[1]*(-5. + 14.*n[1])) + n[0]*(13. + 2.*n[1]*(-69. + 7.*n[1]*(-5. + 7.*n[1])))));
  double complex denominator = (392.*n[0]*(1. + n[0])*(-1. + 2.*n[0])*n[1]*(1. + n[1])*(-1. + 2.*n[1]));
  double complex out         = numerator/denominator * eft_mat_J(n[0], n[1]);

  return out;
}

static double complex eft_mat_I1301p3101(const double complex * const n)
{
  double complex numerator   = ((-7. + 9.*n[0])*ctan(n[0]*_PI_));
  double complex denominator = (336.*n[0]*(n[0] + 1.)*(n[0] - 1.)*(n[0] - 2.)*(n[0] - 3.)*_PI_);
  double complex out         = numerator/denominator;

  return out;
}

static double complex eft_mat_Idelta201(const double complex * const n)
{
  double complex numerator   = ((-3. + 2.*n[0] + 2.*n[1])*(-8. + 7.*n[0] + 7.*n[1]));
  double complex denominator = (28.*n[0]*n[1]);
  double complex out         = numerator/denominator * eft_mat_J(n[0], n[1]);

  return out;
}

static double complex eft_mat_IG201(const double complex * const n)
{
  double complex numerator   = -((-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(-2. + 7.*n[0] + 7.*n[1]));
  double complex denominator = (56.*n[0]*(1. + n[0])*n[1]*(1. + n[1]));
  double complex out         = numerator/denominator * eft_mat_J(n[0], n[1]);

  return out;
}

static double complex eft_mat_FG201(const double complex * const n)
{
  double complex numerator   =  -15. * ctan(n[0]*_PI_);
  double complex denominator = 112.*n[0]*(n[0] + 1.)*(n[0] - 1.)*(n[0] - 2.)*(n[0] - 3.)*_PI_;
  // double complex denominator = 112.*n[0]*(-6. + 5.*n[0] + 5.*n[0]*n[0] - 5.* cpow(n[0],3.) + cpow(n[0],4.))*_PI_;
  double complex out         = numerator/denominator;

  return out;
}

static double complex eft_mat_J12101(const double complex * const n)
{
  double complex numerator   = (9.*ctan(n[0]*_PI_));
  //double complex denominator = (224.*n[0]*(-6. + 11.*n[0] - 6.*n[0]*n[0] + cpow(n[0],3))*_PI_);
  double complex denominator = 224.*n[0]*(n[0] - 1.)*(n[0] - 2.)*(n[0] - 3.)*_PI_;
  double complex out         = numerator/denominator;

  return out;
}

static double complex eft_mat_J21101(const double complex * const n)
{
  double complex numerator   = ((-3. + n[0] + n[1])*(-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(-5. + n[0]*(-4. + 7.*n[0] + 7.*n[1])));
  double complex denominator = (14.*n[0]*(1. + n[0])*(-3. + 2.*n[0])*(-1. + 2.*n[0])*n[1]);
  double complex out         = numerator/denominator * eft_mat_M1(n[0], n[1]);

  return out;
}

static double complex eft_mat_Jdelta201(const double complex * const n)
{
  double complex numerator   = ((-3. + n[0] + n[1])*(-3. + 2.*n[0] + 2.*n[1]));
  double complex denominator = (n[0]*(-3. + 2.*n[0]));
  double complex out         = numerator/denominator * eft_mat_M1(n[0], n[1]);

  return out;
}

static double complex eft_mat_JG201(const double complex * const n)
{
  double complex numerator   = -0.5*((-3. + n[0] + n[1])*(-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1]));
  double complex denominator = (n[0]*(1. + n[0])*(-3. + 2.*n[0])*n[1]);
  double complex out         = numerator/denominator * eft_mat_M1(n[0], n[1]);

  return out;
}

/*
      _________________ 2-nd Moment _______________________
 */

static double complex eft_mat_J12102x(const double complex * const n)
{
  double complex numerator   = (-9.*ctan(n[0]*_PI_));
  double complex denominator = (224.*n[0]*(n[0] + 1.)*(n[0] - 1.)*(n[0] - 2.)*(n[0] - 3.)*_PI_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_J12102y(const double complex * const n)
{
  // double complex summand1   = (3.*(3. + 2.*(-2. + n[0])*n[0])*ctan(n[0]*_PI_))/(224.*(-3. + n[0])*(-2. + n[0])*(-1. + n[0])*n[0]*(1. + n[0])*_PI_);
  // //double complex summand2 = (3.*(1. - 2.*n[0])*ctan(n[0]*_PI_))/(224.*n[0]*(2. - n[0] - 2.*n[0]*n[0] + cpow(n[0],3.))*_PI_);
  // double complex summand2 = (3.*(1. - 2.*n[0])*ctan(n[0]*_PI_))/(224.*n[0]*(n[0] + 1.)*(n[0] - 1.)*(n[0] - 2.)*_PI_);
  // double complex out = summand1 + summand2;
  double complex numerator   = (9.*ctan(n[0]*_PI_));
  double complex denominator = (224.*(n[0] + 1.)*(n[0] - 1.)*(n[0] - 2.)*(n[0] - 3.)*_PI_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_J21102x(const double complex * const n)
{
  double complex numerator   = ((-1. + 2.*n[0])*(-1. + 2.*n[1])*(-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(6. + 7.*n[0] + 7.*n[1])*cGamma(-2.*n[0])*
                                cGamma(-2.*n[1])*cGamma(2.*(-2. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  double complex denominator = (28.*(1. + n[0])*(1. + n[1])*_PI_CUBED_);
  double complex out = - numerator/denominator;

  return out;
}

static double complex eft_mat_J21102y(const double complex * const n)
{
  double complex numerator   = ((-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(34. + n[0] + n[1] + 56.*n[0]*n[0]*n[0]*n[1] - 54.*n[1]*n[1] +
                                2.*n[0]*n[0]*(-27. - 2.*n[1] + 56.*n[1]*n[1]) + 4.*n[0]*n[1]*(-21. + n[1]*(-1. + 14.*n[1])))*cGamma(-2.*n[0])*cGamma(-2.*n[1])*
                                cGamma(2.*(-2. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  double complex denominator = (28.*(1. + n[0])*(1. + n[1])*_PI_CUBED_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_Jdelta202x(const double complex * const n)
{
  double complex numerator   = -0.25*((-3. + 2.*n[0] + 2.*n[1])*cGamma(2. - 2.*n[0])*cGamma(2. - 2.*n[1])*cGamma(2.*(-2. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  double complex denominator = (n[0]*n[1]*_PI_CUBED_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_Jdelta202y(const double complex * const n)
{
  double complex numerator   = ((-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*cGamma(2. - 2.*n[0])*cGamma(2. - 2.*n[1])*cGamma(2.*(-2. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  double complex denominator = (4.*n[0]*n[1]*_PI_CUBED_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_JG202x(const double complex * const n)
{
  double complex numerator   = ((-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*cGamma(2. - 2.*n[0])*cGamma(2. - 2.*n[1])*cGamma(2.*(-2. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  double complex denominator = (4.*n[0]*(1. + n[0])*n[1]*(1. + n[1])*_PI_CUBED_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_JG202y(const double complex * const n)
{
  double complex numerator   = -0.25*((1. + n[0] + n[1])*(-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*cGamma(2. - 2.*n[0])*cGamma(2. - 2.*n[1])*cGamma(2.*(-2. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  double complex denominator = (n[0]*(1. + n[0])*n[1]*(1. + n[1])*_PI_CUBED_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_I2211(const double complex * const n)
{
  double complex numerator   =((-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(50. + 98.*n[0]*n[0]*n[0]*n[1] - n[1]*(9. + 35.*n[1]) +
    7.*n[0]*n[0]*(-5. + 2.*n[1]*(-9. + 14.*n[1])) + n[0]*(-9. + 2.*n[1]*(-33. + 7.*n[1]*(-9. + 7.*n[1])))));
  double complex denominator = (392.*n[0]*(1. + n[0])*(-1. + 2.*n[0])*n[1]*(1. + n[1])*(-1. + 2.*n[1]));
  double complex out = numerator/denominator * eft_mat_J(n[0],n[1]);

  return out;
}

static double complex eft_mat_I1311(const double complex * const n)
{
  double complex numerator   = ((-5. + 3.*n[0])*ctan(n[0]*_PI_));
  double complex denominator = (224.*n[0]*(n[0] + 1.)*(n[0] - 1.)*(n[0] - 2.)*(n[0] - 3.)*_PI_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_J12111(const double complex * const n)
{
  double complex numerator   = (9.*ctan(n[0]*_PI_));
  double complex denominator = (224.*n[0]*(-6. + 11.*n[0] - 6.*n[0]*n[0] + n[0]*n[0]*n[0])*_PI_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_J21111(const double complex * const n)
{
  double complex numerator   = ((-3. + n[0] + n[1])*(-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(-3. + n[0]*(-8. + 7.*n[0] + 7.*n[1])));
  double complex denominator = (14.*n[0]*(1. + n[0])*(-3. + 2.*n[0])*(-1. + 2.*n[0])*n[1]);
  double complex out = numerator/denominator * eft_mat_M1(n[0],n[1]);

  return out;
}

static double complex eft_mat_N11x(const double complex * const n)
{
  double complex numerator   = (((1. + n[0])*(-1. + 2.*n[1])*cGamma(2. - 2.*n[0])*cGamma(1. - 2.*n[1]) + 2.*n[0]*n[1]*(-3. + 2.*n[1])*cGamma(-2.*n[0])*cGamma(2. - 2.*n[1]))*
                                cGamma(2.*(-1. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  double complex denominator = (8.*n[0]*(1. + n[0])*n[1]*(-2. + n[0] + n[1])*_PI_CUBED_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_N11y(const double complex * const n)
{
  double complex numerator   = ((-1. + 2.*n[0] + 2.*n[1])*(((-3. + n[0] + n[1])*(-2. + n[0] + n[1])*(-2. + n[0] + n[1])*(-3. + 2.*n[0] + 2.*n[1])*pow(_PI_,1.5)*
      ((cGamma(2.5 - n[0])*cGamma(1.5 - n[1])*cGamma(-2.5 + n[0] + n[1]))/(cGamma(-1. + n[0])*cGamma(4. - n[0] - n[1])*cGamma(n[1])) +
        (cGamma(1.5 - n[0])*(-((cGamma(2.5 - n[1])*cGamma(-2.5 + n[0] + n[1]))/(cGamma(4. - n[0] - n[1])*cGamma(-1. + n[1]))) +
              (cGamma(1.5 - n[1])*cGamma(-1.5 + n[0] + n[1]))/(cGamma(3. - n[0] - n[1])*cGamma(n[1]))))/cGamma(n[0])))/
    (n[0]*(-3. + 2.*n[0])*n[1]*(-1. + 2.*n[1])) + (2.*(-3. + 2.*n[0])*cGamma(2. - 2.*n[0])*cGamma(1. - 2.*n[1])*cGamma(2.*(-1. + n[0] + n[1]))*
      csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_))/(n[0]*n[1]) +
    (4.*(-1. + 2.*n[0])*cGamma(-2.*n[0])*cGamma(2. - 2.*n[1])*cGamma(2.*(-1. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_))/
    (1. + n[0])));
  double complex denominator = (16.*(-2. + n[0] + n[1])*_PI_CUBED_);
  double complex out = numerator/denominator;

  return out;
}


/*
      _________________ 3-rd Moment _______________________
 */

static double complex eft_mat_J21112x(const double complex * const n)
{
  double complex numerator   = -((-1. + 2.*n[0])*(-1. + 2.*n[1])*(-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(-2. + 7.*n[0] + 7.*n[1])*cGamma(-2.*n[0])*
    cGamma(-2.*n[1])*cGamma(2.*(-2. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  double complex denominator = 28.*(1. + n[0])*(1. + n[1])*_PI_CUBED_;
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_J21112y(const double complex * const n)
{
  double complex numerator   = ((-3. + 2.*n[0] + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(26. + 56.*n[0]*n[0]*n[0]*n[1] + (9. - 38.*n[1])*n[1] +
    2.*n[0]*n[0]*(-19. + 2.*n[1]*(-9. + 28.*n[1])) + n[0]*(9. + 4.*n[1]*(-21. + n[1]*(-9. + 14.*n[1]))))*cGamma(-2.*n[0])*cGamma(-2.*n[1])*
    cGamma(2.*(-2. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  double complex denominator = (28.*(1. + n[0])*(1. + n[1])*_PI_CUBED_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_J12112x(const double complex * const n)
{
  double complex numerator   = (-9.*ctan(n[0]*_PI_));
  double complex denominator = (224.*n[0]*(n[0] + 1.)*(n[0] - 1.)*(n[0] - 2.)*(n[0] - 3.)*_PI_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_J12112y(const double complex * const n)
{
  double complex numerator   = (9.*ctan(n[0]*_PI_));
  double complex denominator = (224.*(n[0] + 1.)*(n[0] - 1.)*(n[0] - 2.)*(n[0] - 3.)*_PI_);
  double complex out = numerator/denominator;

  return out;
}

static double complex eft_mat_N12x(const double complex * const n)
{
  double complex n12 = n[0] + n[1];

  // double complex numerator1   = (3.*(-1. + 2.*n[0] + 2.*n[1])*cGamma(2. - 2.*n[0])*cGamma(-2.*n[1])*cGamma(2.*(-1. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator1 = (8.*n[0]*(1. + n[0])*(-2. + n[0] + n[1])*_PI_CUBED_);

  // double complex numerator2   = -(3.*n[1]*(-1. + 2.*n[0] + 2.*n[1])*cGamma(2. - 2.*n[0])*cGamma(-2.*n[1])*cGamma(2.*(-1. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator2 = (4.*n[0]*(1. + n[0])*(-2. + n[0] + n[1])*_PI_CUBED_);

  // double complex numerator3   = (cGamma(-2.*n[0])*cGamma(2. - 2.*n[1])*cGamma(2.*(n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator3 = (8.*(1. + n[0])*n[1]*(-1. + n[0] + n[1])*_PI_CUBED_);

  double complex numerator1   = (3.*(-1. + 2.*n12)*cGamma(2. - 2.*n[0])*cGamma(-2.*n[1])*cGamma(2.*(-1. + n12))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator1 = (8.*n[0]*(1. + n[0])*(-2. + n12)*_PI_CUBED_);

  double complex numerator2   = -(3.*n[1]*(-1. + 2.*n12)*cGamma(2. - 2.*n[0])*cGamma(-2.*n[1])*cGamma(2.*(-1. + n12))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator2 = (4.*n[0]*(1. + n[0])*(-2. + n[0] + n[1])*_PI_CUBED_);

  double complex numerator3   = (cGamma(-2.*n[0])*cGamma(2. - 2.*n[1])*cGamma(2.*n12)*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator3 = (8.*(1. + n[0])*n[1]*(-1. + n12)*_PI_CUBED_);

  double complex out = numerator1/denominator1 + numerator2/denominator2 + numerator3/denominator3;

  return out;
}

static double complex eft_mat_N12y(const double complex * const n)
{
  double complex n12 = n[0] + n[1];

  // double complex numerator1   = -0.125*((-3. + 2.*n[0])*(-1. + 2.*n[0] + 2.*n[1])*(1. + 2.*n[0] + 2.*n[1])*cGamma(2. - 2.*n[0])*cGamma(-2.*n[1])*cGamma(2.*(-1. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator1 = (n[0]*(1. + n[0])*(-2. + n[0] + n[1])*_PI_CUBED_);

  // double complex numerator2   = ((-1. + 2.*n[0])*(1. + 2.*n[0] + 2.*n[1])*cGamma(-2.*n[0])*cGamma(-2.*n[1])*cGamma(2.*(n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator2 = (4.*(1. + n[0])*(-1. + n[0] + n[1])*_PI_CUBED_);

  double complex numerator1   = -0.125*((-3. + 2.*n[0])*(-1. + 2.*n12)*(1. + 2.*n12)*cGamma(2. - 2.*n[0])*cGamma(-2.*n[1])*cGamma(2.*(-1. + n12))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator1 = (n[0]*(1. + n[0])*(-2. + n12)*_PI_CUBED_);

  double complex numerator2   = ((-1. + 2.*n[0])*(1. + 2.*n12)*cGamma(-2.*n[0])*cGamma(-2.*n[1])*cGamma(2.*n12)*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator2 = (4.*(1. + n[0])*(-1. + n12)*_PI_CUBED_);

  double complex out = numerator1/denominator1 + numerator2/denominator2;

  return out;
}

/*
      _________________ 4-th Moment _______________________
 */

static double complex eft_mat_N22x(const double complex * const n)
{
  double complex n12 = n[0] + n[1];

  // double complex numerator1   = (-3.*(-1. + 2.*n[0])*(-1. + 2.*n[0] + 2.*n[1])*cGamma(-2.*n[0])*cGamma(-2.*(1 + n[1]))*cGamma(2.*(-1. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator1 = (4.*(1. + n[0])*(-2. + n[0] + n[1])*_PI_CUBED_);

  // double complex numerator2   = (3.*(-1. + 2.*n[0])*n[1]*n[1]*(-1. + 2.*n[0] + 2.*n[1])*cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2.*(-1. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator2 = ((1. + n[0])*(-2. + n[0] + n[1])*_PI_CUBED_);

  double complex numerator1   = (-3.*(-1. + 2.*n[0])*(-1. + 2.*n12)*cGamma(-2.*n[0])*cGamma(-2.*(1 + n[1]))*cGamma(2.*(-1. + n12))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator1 = (4.*(1. + n[0])*(-2. + n12)*_PI_CUBED_);

  double complex numerator2   = (3.*(-1. + 2.*n[0])*n[1]*n[1]*(-1. + 2.*n12)*cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2.*(-1. + n12))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator2 = ((1. + n[0])*(-2. + n12)*_PI_CUBED_);

  double complex out = numerator1/denominator1 + numerator2/denominator2;

  return out;
}

static double complex eft_mat_N22y(const double complex * const n)
{
  double complex n12 = n[0] + n[1];

  // double complex numerator1   = (3.*(-3. + 2.*n[0])*(-1. + 2.*n[0])*(1. + 2.*n[1])*(-1. + 2.*n[0] + 2.*n[1])*(1. + 2.*n[0] + 2.*n[1])*cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*
  //   cGamma(2.*(-1. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator1 = (2.*(1. + n[0])*(-2. + n[0] + n[1])*_PI_CUBED_);

  // double complex numerator2   = - (3.*(-1. + 2.*n[0])*(1. + 2.*n[0] + 2.*n[1])*cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2.*(n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*
  //   csin((n[0] + n[1])*_PI_));
  // double complex denominator2 = (2.*(1. + n[0])*(-1. + n[0] + n[1])*_PI_CUBED_);

  // double complex numerator3   = - (3.*(-1. + 2.*n[0])*n[1]*(1. + 2.*n[0] + 2.*n[1])*cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2.*(n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*
  //   csin((n[0] + n[1])*_PI_));
  // double complex denominator3 = ((1. + n[0])*(-1. + n[0] + n[1])*_PI_CUBED_);

  // double complex numerator4   = (cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2.*(1. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator4 = (4.*(1. + n[0])*(n[0] + n[1])*_PI_CUBED_);

  // double complex numerator5   = (n[1]*cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2*(1. + n[0] + n[1]))*csin(n[0]*_PI_)*
  //   csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator5 = (2.*(1. + n[0])*(n[0] + n[1])*_PI_CUBED_);

  double complex numerator1   = (3.*(-3. + 2.*n[0])*(-1. + 2.*n[0])*(1. + 2.*n[1])*(-1. + 2.*n12)*(1. + 2.*n12)*cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*
    cGamma(2.*(-1. + n12))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator1 = (2.*(1. + n[0])*(-2. + n12)*_PI_CUBED_);

  double complex numerator2   = - (3.*(-1. + 2.*n[0])*(1. + 2.*n12)*cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2.*n12)*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator2 = (2.*(1. + n[0])*(-1. + n12)*_PI_CUBED_);

  double complex numerator3   = - (3.*(-1. + 2.*n[0])*n[1]*(1. + 2.*n12)*cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2.*n12)*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator3 = ((1. + n[0])*(-1. + n12)*_PI_CUBED_);

  double complex numerator4   = (cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2.*(1. + n12))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator4 = (4.*(1. + n[0])*n12*_PI_CUBED_);

  double complex numerator5   = (n[1]*cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2*(1. + n12))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator5 = (2.*(1. + n[0])*n12*_PI_CUBED_);

  double complex out = numerator1/denominator1 + numerator2/denominator2 + numerator3/denominator3 + numerator4/denominator4 + numerator5/denominator5;

  return out;
}

static double complex eft_mat_N22z(const double complex * const n)
{
  double complex n12 = n[0] + n[1];

  // double complex numerator1   = ((-5. + 2.*n[0])*(-3. + 2.*n[0])*(-1. + 2.*n[0])*(-1. + 2.*n[0] + 2.*n[1])*(1. + 2.*n[0] + 2.*n[1])*(3. + 2.*n[0] + 2.*n[1])*cGamma(-2.*n[0])*
  // cGamma(-2.*(1. + n[1]))*cGamma(2.*(-1. + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator1 = (4.*(1. + n[0])*(-2. + n[0] + n[1])*_PI_CUBED_);

  // double complex numerator2   = - ((-3. + 2.*n[0])*(-1. + 2.*n[0])*(1. + 2.*n[0] + 2.*n[1])*(3. + 2.*n[0] + 2.*n[1])*
  // cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2.*(n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator2 = (2.*(1. + n[0])*(-1. + n[0] + n[1])*_PI_CUBED_);

  // double complex numerator3   = ((-1. + 2.*n[0])*(3. + 2.*n[0] + 2.*n[1])*cGamma(-2.*n[0])*cGamma(-2.*(1 + n[1]))*
  // cGamma(2.*(1 + n[0] + n[1]))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin((n[0] + n[1])*_PI_));
  // double complex denominator3 = (4.*(1. + n[0])*(n[0] + n[1])*_PI_CUBED_);

  double complex numerator1   = ((-5. + 2.*n[0])*(-3. + 2.*n[0])*(-1. + 2.*n[0])*(-1. + 2.*n12)*(1. + 2.*n12)*(3. + 2.*n12)*cGamma(-2.*n[0])*
    cGamma(-2.*(1. + n[1]))*cGamma(2.*(-1. + n12))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator1 = (4.*(1. + n[0])*(-2. + n12)*_PI_CUBED_);

  double complex numerator2   = - ((-3. + 2.*n[0])*(-1. + 2.*n[0])*(1. + 2.*n12)*(3. + 2.*n12)*
    cGamma(-2.*n[0])*cGamma(-2.*(1. + n[1]))*cGamma(2.*n12)*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator2 = (2.*(1. + n[0])*(-1. + n12)*_PI_CUBED_);

  double complex numerator3   = ((-1. + 2.*n[0])*(3. + 2.*n12)*cGamma(-2.*n[0])*cGamma(-2.*(1 + n[1]))*
    cGamma(2.*(1 + n12))*csin(n[0]*_PI_)*csin(n[1]*_PI_)*csin(n12*_PI_));
  double complex denominator3 = (4.*(1. + n[0])*n12*_PI_CUBED_);

  double complex out = numerator1/denominator1 + numerator2/denominator2 + numerator3/denominator3;

  return out;
}


/** static array holding function pointers to the individual moments */
const static double complex (*loop_mat[NUM_MOMENTS])(const double complex * const n)      \
  = {eft_mat_I2200, eft_mat_I1300, eft_mat_Idelta200, eft_mat_IG200, eft_mat_Idelta2delta200, eft_mat_IG2G200, eft_mat_Idelta2G200, eft_mat_FG200,  \
     eft_mat_I2201, eft_mat_Idelta201, eft_mat_IG201, eft_mat_J21101, eft_mat_Jdelta201, eft_mat_JG201, eft_mat_FG201, eft_mat_I1301p3101, eft_mat_J12101, NULL,  \
     eft_mat_J21102x, eft_mat_J21102y, eft_mat_Jdelta202x, eft_mat_Jdelta202y, eft_mat_JG202x, eft_mat_JG202y, eft_mat_I2211, eft_mat_J21111, eft_mat_N11x, eft_mat_N11y, eft_mat_J12102x, eft_mat_J12102y, eft_mat_I1311, eft_mat_J12111, NULL,  \
     eft_mat_J21112x, eft_mat_J21112y, eft_mat_N12x, eft_mat_N12y, eft_mat_J12112x, eft_mat_J12112y,                                \
     eft_mat_N22x, eft_mat_N22y, eft_mat_N22z};


int eft_compute_loop_matrices(struct eft * peft) {

  int index_moment, it1, it2, use_tracer;

  double complex n[2];
  double * const n1_real = (double *)&(n[0]);
  double * const n1_imag = (double *)&(n[0]) + 1;
  double * const n2_real = (double *)&(n[1]);
  double * const n2_imag = (double *)&(n[1]) + 1;

  for (index_moment = 0; index_moment < peft->index_num; index_moment++) {
    use_tracer = peft->use_tracer[index_moment];  // must be set for each individual moment
    /** - generate matrices for different symmetry types in LAPACK-compatible storage schemes */
    switch (peft->symmetry[index_moment])
    {
    case no_finite_part:
      break;

    case sym_vec:
      *n1_real = -0.5 * peft->hp->bias[use_tracer];
      for (it1 = 0; it1 < peft->hp->fourier_coeff_size; it1++) {
        *n1_imag = -0.5 * peft->fourier_frequencies[use_tracer][it1];
        peft->loop_matrices[index_moment][it1] = (*loop_mat[index_moment])((const double complex * const)&n);
      }
      break;

    case sym_mat_none:
      *n1_real = -0.5 * peft->hp->bias[use_tracer];
      *n2_real = -0.5 * peft->hp->bias[use_tracer];
      for (it2 = 0; it2 < peft->hp->fourier_coeff_size; it2++) {  /** - row-major array is best filled in this order */
        for (it1 = 0; it1 < peft->hp->fourier_coeff_size; it1++) {
          *n1_imag = -0.5 * peft->fourier_frequencies[use_tracer][it1];
          *n2_imag = -0.5 * peft->fourier_frequencies[use_tracer][it2];
          peft->loop_matrices[index_moment][it1 + it2*peft->hp->fourier_coeff_size] = (*loop_mat[index_moment])((const double complex * const)&n);
        }
      }
      break;

    case sym_mat_symmetric:
      *n1_real = -0.5 * peft->hp->bias[use_tracer];
      *n2_real = -0.5 * peft->hp->bias[use_tracer];
      for (it2 = 0; it2 < peft->hp->fourier_coeff_size; it2++) {
        for (it1 = 0; it1 <= it2; it1++) {
          *n1_imag = -0.5 * peft->fourier_frequencies[use_tracer][it1];
          *n2_imag = -0.5 * peft->fourier_frequencies[use_tracer][it2];
          /** - upper-triangular packed storage */
          peft->loop_matrices[index_moment][it1 + it2*(it2+1)/2] = (*loop_mat[index_moment])((const double complex * const)&n);
        }
      }
      break;

    default:
       break;
    }
  }

  return _SUCCESS_;
}
