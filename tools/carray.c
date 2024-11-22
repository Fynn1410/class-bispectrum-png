
#include "carray.h"


#ifdef __FAST_MATH__
#warning Kahan-Neumaier summation is incompatible with -ffast-math. Defaulting to uncompensated summation.

  inline double d2sum(double a, double b, double *__restrict__ t) {
    double s;
    s = a + b;
    *t = 0.;
    return s;
  }

/*
  inline double dcompsum(double * summands, int size, int stride) {
    int i;
    double s;

    s = 0.;
    for (i = 0; i < size; i++) {
      s += summands[i];
    }
    return s;
  }
*/

/*
  inline class_complex c2sum(class_complex a, class_complex b, class_complex * __restrict__ t) {
    class_complex s;
    s = a + b;
    *t = class_complex(0., 0.);
    return s;
  }
*/

  inline class_complex ccompsum(class_complex * summands, int size, int stride) {
    int i;
    class_complex s;

    s = 0.;
    for (i = 0; i < size; i++) {
      s += summands[i];
    }
    return s;
  }

#else

  /**
   * @brief Computes the floating-point error t of the sum a + b, such that a + b = s + t.
   * @param a    Input: first summand
   * @param b    Input: second summand
   * @param t    Output: floating-point error of the operation a + b
   * @return the sum s = float(a + b)
   */
  inline double d2sum(double a, double b, double * __restrict__ t) {
    double s, ap, bp, da, db;
    s = a + b;
    ap = s - b;
    bp = s - ap;
    da = a - ap;
    db = b - bp;
    *t = da + db;
    return s;
  }

  /**
   * @brief Kahan-Babushka-Neumaier compensated summation for an array of doubles.
   * @param summands  Input: array of summand values
   * @param size      Input: number of summands to use
   * @param stride    Input: stride in the summands array
   * @return compensated sum of the input array
   */
/*
  inline double dcompsum(double * summands, int size, int stride) {
    int i;
    double s, c, t, in;

    s = 0.; c = 0.;
    for (i = 0; i < size; i++) {
      in = summands[i];
      t = s + in;
      if (fabs(s) >= fabs(in)) {
        c += (s - t) + in;
      }
      else {
        c += (in - t) + s;
      }
      s = t;
    }

    return s + c;
  }
*/

  /**
   *  @brief Computes the floating-point error t of the complex sum a + b, such that a + b = s + t.
   *  @param a    Input: first summand
   *  @param b    Input: second summand
   *  @param t    Output: floating-point error of the operation a + b
   *  @return the complex sum s = float(a + b)
   */
/*
  inline class_complex c2sum(class_complex a, class_complex b, class_complex * __restrict__ t) {
    class_complex s, ap, bp, da, db;
    s = a + b;
    ap = s - b;
    bp = s - ap;
    da = a - ap;
    db = b - bp;
    *t = da + db;
    return s;
  }
*/

  /**
   * @brief Kahan-Babushka-Neumaier compensated summation for an array of complex doubles.
   * @param summands  Input: array of summand values
   * @param size      Input: number of summands to use
   * @param stride    Input: stride in the summands array
   * @return compensated sum of the input array
   */
  inline class_complex ccompsum(class_complex * summands, int size, int stride) {
    int i;
    double sr, si, cr, ci, tr, ti, inr, ini;

    sr = 0.; si = 0.; cr = 0.; ci = 0.;
    for (i = 0; i < size; i++) {
      inr = creal(summands[i]);
      ini = cimag(summands[i]);
      tr = sr + inr;
      ti = si + ini;
      if (fabs(sr) >= fabs(inr)) {
        cr += (sr - tr) + inr;
      }
      else {
        cr += (inr - tr) + sr;
      }
      if (fabs(si) >= fabs(ini)) {
        ci += (si - ti) + ini;
      }
      else {
        ci += (ini - ti) + si;
      }
      sr = tr;
      si = ti;
    }

    return class_complex(sr + cr, si + ci);
  }

#endif

int array_integrate_internal_exponential(
                                         double * x0,
                                         int x_size,
                                         int x_stride,
                                         double * y0,
                                         int y_stride,
                                         double * ddy0,
                                         int ddy_stride,
                                         class_complex phase,
                                         class_complex * result,
                                         ErrorMsg errmsg) {

  int index_x;
  double h1, h2;                     /**< distance between control points i & (i-1) and (i+1) & i respectively */
  class_complex phase_invpow[4];    /**< contains powers of phase as phase_pow[n] = phase^(-n-1) */
  class_complex sum;

  for (index_x = 0; index_x < 4; index_x++) phase_invpow[index_x] = cpow(phase, -(index_x+1));

  h1 = x0[1*x_stride] - x0[0*x_stride];
  h2 = x0[(x_size-1)*x_stride] - x0[(x_size-2)*x_stride];

  *result = ( phase_invpow[0] * y0[0*y_stride]
            + phase_invpow[1] * ((y0[1*y_stride] - y0[0*y_stride])/h1
                                - h1/6.*(2*ddy0[0*ddy_stride] + ddy0[0*ddy_stride]))
            + phase_invpow[2] * ddy0[0*ddy_stride]
            + phase_invpow[3] * (ddy0[1*ddy_stride] - ddy0[0*ddy_stride])/h1    \
            ) * cexp(-phase * x0[0*x_stride]);
  *result -=( phase_invpow[0] * y0[(x_size-1)*y_stride]
            + phase_invpow[1] * ((y0[(x_size-1)*y_stride] - y0[(x_size-2)*y_stride])/h2
                                + h2/6.*(2*ddy0[(x_size-1)*ddy_stride] + ddy0[(x_size-2)*ddy_stride]))
            + phase_invpow[2] * ddy0[(x_size-1)*ddy_stride]
            + phase_invpow[3] * (ddy0[(x_size-1)*ddy_stride] - ddy0[(x_size-2)*ddy_stride])/h2    \
            ) * cexp(-phase * x0[(x_size-1)*x_stride]);

  sum = 0.;
  for (index_x = 1; index_x < x_size-2; index_x++) {
    h1 = x0[index_x*x_stride] - x0[(index_x-1)*x_stride];
    h2 = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];

    sum += ( (ddy0[(index_x+1)*ddy_stride] - ddy0[index_x*ddy_stride]) / h2
            -(ddy0[index_x*ddy_stride] - ddy0[(index_x-1)*ddy_stride]) / h1       \
            ) * cexp(-phase * x0[index_x*x_stride]);
  }

  *result += phase_invpow[3] * sum;

  return _SUCCESS_;
}















/**
 * @brief Computes the squared spline integral.
 *        dI = dx S(x)^2
 * @param x           Input: contains x-values of the integration range
 * @param x_size      Input: size of x-array to be used
 * @param y_array     Input: contains y-values of the integrand with elements
 *                            y_array[index_x*y_size + index_y]
 * @param y_size      Input: number of columns in y-array
 * @param ddy_array   Input: contains y''-values of the integrand obtained from splining with elements
 *                            ddy_array[index_x*y_size + index_y]
 * @param result      Output: integration result I with elements[index_y]
 * @param errmsg
 *
 * @return the error status
 */
int array_square_integrate_all_spline_table_lines(
		      double * x,
			    int x_size,
			    double * y_array,
			    int y_size,
			    double * ddy_array,
			    double * result,
          ErrorMsg errmsg) {

  int index_y;

  for (index_y = 0; index_y < y_size; index_y++) {
    array_square_integrate_internal(x,
                                    x_size,
                                    1,
                                    y_array + index_y,
                                    y_size,
                                    ddy_array + index_y,
                                    y_size,
                                    result + index_y,
                                    errmsg);
  }

  return _SUCCESS_;
}

/**
 * @brief Computes the squared spline integral with an exponential bias a and arbitrary derivative order m.
 *        dI = dx (exp(a x) S^(m)(x))^2
 * @param x           Input: contains x-values of the integration range
 * @param x_size      Input: size of x-array to be used
 * @param y_array     Input: contains y-values of the integrand with elements
 *                            y_array[index_x*y_size + index_y]
 * @param y_size      Input: number of columns in y-array
 * @param ddy_array   Input: contains y''-values of the integrand obtained from splining with elements
 *                            ddy_array[index_x*y_size + index_y]
 * @param bias        Input: exponential bias a
 * @param derivative_order  Input: order of the spline derivative m (there are 3 non-zero derivatives)
 * @param result      Output: integration result I with elements[index_y]
 * @param errmsg
 *
 * @return the error status
 */
int array_square_integrate_exponential_all_spline_table_lines(
		      double * x,
			    int x_size,
			    double * y_array,
			    int y_size,
			    double * ddy_array,
          double bias,
          int derivative_order,
			    double * result,
          ErrorMsg errmsg) {

  int index_y, num_threads = 1;

  #ifdef _OPENMP
  num_threads = omp_get_max_threads();
  #endif

  for (index_y = 0; index_y < y_size; index_y++) {
    array_square_integrate_exponential_internal(x,
                                                x_size,
                                                1,
                                                y_array + index_y,
                                                y_size,
                                                ddy_array + index_y,
                                                y_size,
                                                bias,
                                                derivative_order,
                                                result + index_y,
                                                num_threads,
                                                errmsg);
  }

  return _SUCCESS_;
}






/**
 * @brief Computes the spline integral with a specified (complex) Exponential window exactly.
 *        Useful for computing Fourier coefficients.
 *        dI(p) = dx S(x) * exp(-p*x)
 * @param array       Input: contains x, y and y'' values retrieved from splining
 * @param n_columns   Input: number of columns in array, indexed by index_x/y/ddy
 * @param n_lines     Input: number of used control points
 * @param index_x     Input: index for x-values (->class_define_index)
 * @param index_y     Input: index for y-values (->class_define_index)
 * @param index_ddy   Input: index for y''-values (->class_define_index)
 * @param phase       Input: phase factor p which may be complex
 * @param result      Output: complex integration result I(p)
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_spline_exponential(
          double * array,
          int n_columns,
          int n_lines,
          int index_x,   /** from 0 to (n_columns-1) */
          int index_y,
          int index_ddy,
          class_complex phase,
          class_complex * result,
          ErrorMsg errmsg) {


  class_test(n_lines<2,
             errmsg,
             "integral is zero with less than one spline segment");

  array_integrate_internal_exponential(array + index_x,
                                       n_lines,
                                       n_columns,
                                       array + index_y,
                                       n_columns,
                                       array + index_ddy,
                                       n_columns,
                                       phase,
                                       result,
                                       errmsg);

  return _SUCCESS_;
}

/**
 * @brief Computes the spline integral with a specified (complex) Exponential window exactly.
 *        Useful for computing Fourier coefficients.
 *        dI(p) = dx S(x) * exp(-p*x)
 * @param x           Input: contains x-values of the integration range
 * @param x_size      Input: size of x-array to be used
 * @param y_array     Input: contains y-values of the integrand with elements
 *                            y_array[index_x*y_size + index_y]
 * @param y_size      Input: number of columns in y-array
 * @param ddy_array   Input: contains y''-values of the integrand obtained from splining with elements
 *                            ddy_array[index_x*y_size + index_y]
 * @param phase       Input: phase factor p which may be complex
 * @param result      Output: complex integration result I(p)
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_spline_table_lines_exponential(
          double * x,
		      int x_size,
		      double * y_array,
		      int y_size,
		      double * ddy_array,
          class_complex phase,
          class_complex * result,
          ErrorMsg errmsg) {

  int i, index_x, index_y;
  double h1, h2;                     /**< distance between control points i & (i-1) and (i+1) & i respectively */
  class_complex phase_invpow[4];    /**< contains powers of phase as phase_pow[n] = phase^(-n-1) */
  class_complex sum;

  class_test(x_size<2,
             errmsg,
             "integral is zero with less than one spline segment");

  for (index_y = 0; index_y < y_size; index_y++) {
    array_integrate_internal_exponential(x,
                                         x_size,
                                         1,
                                         y_array + index_y,
                                         y_size,
                                         ddy_array + index_y,
                                         y_size,
                                         phase,
                                         result + index_y,
                                         errmsg);
  }

  return _SUCCESS_;
}

int array_integrate_internal_exponential_pure_phase(
                        double * x0,
                        int x_size,
                        int x_stride,
                        double * y0,
                        int y_stride,
                        double * ddy0,
                        int ddy_stride,
                        double phase,
                        class_complex * result,
                        class_complex * condition_number,
                        ErrorMsg errmsg) {

  int index_x;
  double h1, h2;              /**< distance between control points i & (i-1) and (i+1) & i respectively */
  double phase_invpow[4];     /**< contains powers of phase as phase_pow[n] = phase^(-n-1) */
  class_complex sum;
  class_complex sum_abs;
  double spline_factor, real, imag;

  for (index_x = 0; index_x < 4; index_x++) phase_invpow[index_x] = pow(phase, -(index_x+1));

  h1 = x0[1*x_stride] - x0[0*x_stride];
  h2 = x0[(x_size-1)*x_stride] - x0[(x_size-2)*x_stride];

  *result = class_create_complex(- phase_invpow[1] * ((y0[1*y_stride] - y0[0*y_stride])/h1                 \
                                      - h1/6.*(2*ddy0[0*ddy_stride] + ddy0[1*ddy_stride]))  \
                  + phase_invpow[3] * (ddy0[1*ddy_stride] - ddy0[0*ddy_stride])/h1 ,        \
                  - phase_invpow[0] * y0[0*y_stride]          \
                  + phase_invpow[2] * ddy0[0*ddy_stride] )    \
            * class_create_complex( cos(phase * x0[0*x_stride]), -sin(phase * x0[0*x_stride]) );
  *condition_number = class_create_complex( fabs(creal(*result)), fabs(cimag(*result)) );

  sum    = -class_create_complex(- phase_invpow[1] * ((y0[(x_size-1)*y_stride] - y0[(x_size-2)*y_stride])/h2                 \
                                      + h2/6.*(2*ddy0[(x_size-1)*ddy_stride] + ddy0[(x_size-2)*ddy_stride]))  \
                  + phase_invpow[3] * (ddy0[(x_size-1)*ddy_stride] - ddy0[(x_size-2)*ddy_stride])/h2 ,        \
                  - phase_invpow[0] * y0[(x_size-1)*y_stride]         \
                  + phase_invpow[2] * ddy0[(x_size-1)*ddy_stride] )   \
            * class_create_complex( cos(phase * x0[(x_size-1)*x_stride]), -sin(phase * x0[(x_size-1)*x_stride]) );
  *result += sum;
  *condition_number += class_create_complex( fabs(creal(sum)), fabs(cimag(sum)) );

  sum = 0.;
  sum_abs = 0.;
  for (index_x = 1; index_x < x_size-1; index_x++) {
    h1 = x0[index_x*x_stride] - x0[(index_x-1)*x_stride];
    h2 = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
    spline_factor = (ddy0[(index_x+1)*ddy_stride] - ddy0[index_x*ddy_stride]) / h2       \
                   -(ddy0[index_x*ddy_stride] - ddy0[(index_x-1)*ddy_stride]) / h1;
    real = cos(phase * x0[index_x*x_stride]);
    imag = -sin(phase * x0[index_x*x_stride]);

    sum += spline_factor * class_create_complex(real, imag);
    sum_abs += fabs(spline_factor) * class_create_complex(fabs(real), fabs(imag));
  }

  *result += phase_invpow[3] * sum;
  *condition_number += phase_invpow[3] * sum_abs;
  *condition_number = class_create_complex( creal(*condition_number) / fabs(creal(*result)),       \
                             cimag(*condition_number) / fabs(cimag(*result)) );

  return _SUCCESS_;
}

/** - assumes periodic splines with S(x0) = S(xn), S'(x0) = S'(xn), S''(x0) = S''(xn) */
int array_integrate_internal_exponential_pure_phase_compensated(
                        double * x0,
                        int x_size,
                        int x_stride,
                        double * y0,
                        int y_stride,
                        double * ddy0,
                        int ddy_stride,
                        double phase,
                        class_complex * result,
                        double condition_number_threshold,
                        class_complex * condition_number,
                        ErrorMsg errmsg) {

  int index_x;
  double h1, h2, dM1, dM2, t1, t2;    /**< distances between control points, spline moment differences and error variable */
  double phase_invquart;              /**< contains phase^(-4) */
  class_complex summands[x_size-1];  /**< summands of the Spline Fourier sum */
  class_complex sum, sum_abs;        /**< sum and sum of absolute values of summands */

  phase_invquart = pow(phase, -4.);

  /** - compute the summands of the Fourier sum */
  h1 = x0[1*x_stride] - x0[0*x_stride];
  h2 = x0[(x_size-1)*x_stride] - x0[(x_size-2)*x_stride];
  dM1 = d2sum(ddy0[1*x_stride], -ddy0[0*x_stride], &t1);
  dM2 = d2sum(ddy0[(x_size-1)*x_stride], -ddy0[(x_size-2)*x_stride], &t2);
  summands[0] = ((dM1/h1 - dM2/h2) + (t1/h1 - t2/h2))     \
                * cexp(-_Complex_I * phase * x0[0*x_stride]);
  // printf("%.16e        %.16e \n", creal(summands[0]), cimag(summands[0]));

  for (index_x = 1; index_x < x_size-1; index_x++) {
    h1 = x0[index_x*x_stride] - x0[(index_x-1)*x_stride];
    h2 = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
    dM1 = d2sum(ddy0[index_x*x_stride], -ddy0[(index_x-1)*x_stride], &t1);
    dM2 = d2sum(ddy0[(index_x+1)*x_stride], -ddy0[index_x*x_stride], &t2);

    summands[index_x] = ((dM2/h2 - dM1/h1) + (t2/h2 - t1/h1))     \
                        * cexp(-_Complex_I * phase * x0[index_x*x_stride]);
    // printf("%.16e        %.16e \n", creal(summands[index_x]), cimag(summands[index_x]));
  }

  /** - sum them up naively and compute the condition number */
  sum = 0.; sum_abs = 0.;
  for (index_x = 0; index_x < x_size-1; index_x++) {
    sum += summands[index_x];
    #ifndef __FAST_MATH__
    sum_abs += class_create_complex( fabs(creal(summands[index_x])),
                      fabs(cimag(summands[index_x])) );
    #endif
  }
  *condition_number = class_create_complex( creal(sum_abs)/fabs(creal(sum)),
                          cimag(sum_abs)/fabs(cimag(sum)) );

  /** - if the condition number is higher than a threshold, recompute the sum with compensated summation */
  if ((creal(*condition_number) > condition_number_threshold)     \
      || (cimag(*condition_number) > condition_number_threshold)) {
    sum = ccompsum(summands, x_size-1, 1);
  }

  /** - multiply by an overall factor of phase^(-4) */
  *result = phase_invquart * sum;

  return _SUCCESS_;
}

/**
 * @brief Computes the spline integral with a specified exponential phase factor exactly.
 *        Useful for computing Fourier coefficients.
 *        dI(p) = dx S(x) * exp(-i p*x)
 * @param array       Input: contains x, y and y'' values retrieved from splining
 * @param n_columns   Input: number of columns in array, indexed by index_x/y/ddy
 * @param n_lines     Input: number of used control points
 * @param index_x     Input: index for x-values (->class_define_index)
 * @param index_y     Input: index for y-values (->class_define_index)
 * @param index_ddy   Input: index for y''-values (->class_define_index)
 * @param phase       Input: phase factor p
 * @param result      Output: complex integration result I(p)
 * @param condition_num_threshold   Input: threshold for condition number above which the sum is computed with compensation
 * @param condition_num   Output: condition number of the underlying sum
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_spline_fourier_compensated(
          double * array,
          int n_columns,
          int n_lines,
          int index_x,   /** from 0 to (n_columns-1) */
          int index_y,
          int index_ddy,
          double phase,
          class_complex * result,
          double condition_num_threshold,
          class_complex * condition_num,
          ErrorMsg errmsg) {

  class_test(n_lines<2,
             errmsg,
             "integral is zero with less than one spline segment");

  return array_integrate_internal_exponential_pure_phase_compensated(
                                            array + index_x,
                                            n_lines,
                                            n_columns,
                                            array + index_y,
                                            n_columns,
                                            array + index_ddy,
                                            n_columns,
                                            phase,
                                            result,
                                            condition_num_threshold,
                                            condition_num,
                                            errmsg);
}

/**
 * @brief Computes the spline integral with a specified (complex) Exponential window exactly.
 *        Useful for computing Fourier coefficients.
 *        dI(p) = dx S(x) * exp(-i p*x)
 * @param x           Input: contains x-values of the integration range
 * @param x_size      Input: size of x-array to be used
 * @param y_array     Input: contains y-values of the integrand with elements
 *                            y_array[index_x*y_size + index_y]
 * @param y_size      Input: number of columns in y-array
 * @param ddy_array   Input: contains y''-values of the integrand obtained from splining with elements
 *                            ddy_array[index_x*y_size + index_y]
 * @param phase       Input: phase factor p which may be complex
 * @param result      Output: complex integration result I(p)
 * @param condition_num_threshold   Input: threshold for condition number above which the sum is computed with compensation
 * @param condition_num   Output: condition number of the underlying sum with elements[index_y]
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_spline_table_lines_fourier_compensated(
          double * x,
		      int x_size,
		      double * y_array,
		      int y_size,
		      double * ddy_array,
          double phase,
          class_complex * result,
          double condition_num_threshold,
          class_complex * condition_num,
          ErrorMsg errmsg) {

  int index_y;

  class_test(x_size<2,
             errmsg,
             "integral is zero with less than one spline segment");

  for (index_y = 0; index_y < y_size; index_y++) {
    array_integrate_internal_exponential_pure_phase_compensated(
                                                    x,
                                                    x_size,
                                                    1,
                                                    y_array + index_y,
                                                    y_size,
                                                    ddy_array + index_y,
                                                    y_size,
                                                    phase,
                                                    result + index_y,
                                                    condition_num_threshold,
                                                    condition_num + index_y,
                                                    errmsg);
  }

  return _SUCCESS_;
}

// /**
//  * @brief Computes the value of a power series with coefficients pcoeff, starting from x^(min_pow)
//  *        with coeff_size terms
//  * @param pcoeff      Input: Coefficients of the power series
//  * @param coeff_size  Input: Number of terms in the sum
//  * @param min_pow     Input: Minimal exponent
//  * @param x           Input: function argument
//  *
//  * @return function value at x
//  */
// inline double array_polynomial(double * pcoeff,
//                         int coeff_size,
//                         int min_pow,
//                         double x) {
//   int j;
//   double x_min_pow = pow(x, min_pow);
//   double x_pow = 1., sum;

//   sum = pcoeff[0];
//   for (j = 1; j < coeff_size; j++) {
//     x_pow *= x;
//     sum += pcoeff[j] * x_pow;
//   }

//   return x_min_pow * sum;
// }

// /**
//  * @brief Returns the function value in the trigonometric ring corresponding to index n.
//  *        n = 0 is set to the sine function.
//  * @param x   Input: function argument
//  * @param n   Input: index in the ring
//  *
//  * @return function value at x
//  */
// inline double array_trigonometric_ring(double x,
//                                         int n) {
//   int m = (((n % 4) + 4) % 4);
//   switch (m)
//   {
//   case 0:
//     return sin(x); break;
//   case 1:
//     return cos(x); break;
//   case 2:
//     return -sin(x);break;
//   case 3:
//     return -cos(x);break;

//   default:  /* can never happen */
//     return 0.;break;
//   }
// }

// int array_bessel_series_expansion_integral(
//           double * pcoeff,
//           int n,
//           double x,
//           double * result,
//           ErrorMsg errmsg) {

//   int m, k, ksup;
//   double prefactor = 1., x_pow = 1., term, sum, polysum;
//   double x2_half = x*x/2;

//   /** - compute the prefactor to the series */
//   for (m = 0; m <= n; m++) {
//     prefactor *= x/(2*m+1.);
//   }

//   /** - compute the series expansion of the integral */
//   m = 0;
//   term = sum = 1.;
//   for (k = 1; fabs(term) > __DBL_EPSILON__; k++) {
//     term *= -x2_half / (k * (2*k + 2*n + 1.));
//     sum += (n+m+1.)/(n+2*k+m+1.) * term;
//   }
//   polysum = pcoeff[0] / (n+m+1.) * sum;
//   ksup = k;
//   /** m > 0 with the same number of terms */
//   for (m = 1; m <= 3; m++) {
//     term = sum = 1.;
//     x_pow *= x;
//     for (k = 1; k < ksup; k++) {
//       term *= -x2_half / (k * (2*k + 2*n + 1.));
//       sum += (n+m+1.)/(n+2*k+m+1.) * term;
//     }
//     polysum += pcoeff[m] / (n+m+1.) * sum * x_pow;
//   }

//   *result = prefactor * polysum;

//   return _SUCCESS_;
// }

// /**
//  * @brief Computes the spline integral with a specified spherical Bessel function J_n(x) exactly.
//  *        dI(n, c, s) = dx S(x) * j_n( c(x - s) )
//  * @param array       Input: contains x, y and y'' values retrieved from splining
//  * @param n_columns   Input: number of columns in array, indexed by index_x/y/ddy
//  * @param n_lines     Input: number of used control points
//  * @param index_x     Input: index for x-values (->class_define_index)
//  * @param index_y     Input: index for y-values (->class_define_index)
//  * @param index_ddy   Input: index for y''-values (->class_define_index)
//  * @param n           Input: Bessel function order
//  * @param coeff       Input: Coefficient c in the Bessel function
//  * @param shift       Input: shift variable s in the Bessel function
//  * @param result      Output: Integration result I(n, c, s)
//  * @param errmsg
//  *
//  * @return the error status
//  */
// int array_integrate_all_spline_spherical_bessel_J(
//           double * array,
//           int n_columns,
//           int n_lines,
//           int index_x,   /** from 0 to (n_columns-1) */
//           int index_y,
//           int index_ddy,
//           int n,
//           double coeff,
//           double shift,
//           double * result,
//           ErrorMsg errmsg) {

//   int i, j, k, m, min_pow, factor;
//   double h;                          /**< distance between control points */
//   double a, b, Ms, dM, ys, dy;
//   double prt_res, prt_sum, apow, bpow, pcoeff[4], bessel_j_a[n+1], bessel_j_b[n+1];
//   double intg_func_b, intg_func_a;
//   double coeff_border_ab[4] = {9.9123380e-1, 2.1567720e-2, -3.1745158e-5, -3.6693021e-8};
//   double border = array_polynomial(coeff_border_ab, 4, 0, (double)n);
//   ErrorMsg internal_error;

//   class_test( coeff * (array[0*n_columns+index_x] - shift) < 0.,
//               errmsg,
//               "%s(L:%d) First shifted abscissa is too close to zero for this integration method: %.5e",
//               __func__, __LINE__, coeff * (array[0*n_columns+index_x] - shift));


//   for (i = 0; i < n_lines-1; i++) {
//     h = array[(i+1)*n_columns+index_x] - array[i*n_columns+index_x];
//     a = coeff * (array[i*n_columns+index_x] - shift);
//     b = coeff * (array[(i+1)*n_columns+index_x] - shift);
//     Ms = array[i*n_columns+index_ddy] + array[(i+1)*n_columns+index_ddy];
//     dM = array[i*n_columns+index_ddy] - array[(i+1)*n_columns+index_ddy];
//     ys = array[i*n_columns+index_y] + array[(i+1)*n_columns+index_y];
//     dy = array[i*n_columns+index_y] - array[(i+1)*n_columns+index_y];

//     min_pow = 0;
//     /** - pcoeff initially contains f_{n}(x) */
//     pcoeff[0] = array[i*n_columns+index_y] + a/coeff * (dy/h + h*(3.*Ms + dM)/12.)
//                   +(2.*a*a*a*dM + 3.*a*a*coeff*h*(Ms + dM)) / (12.*coeff*coeff*coeff*h);
//     pcoeff[1] = -1./coeff * (dy/h + h*(3.*Ms + dM)/12.)
//                   -(a*a*dM + a*coeff*h*(Ms + dM)) / (2.*coeff*coeff*coeff*h);
//     pcoeff[2] = array[i*n_columns+index_ddy] / (2.*coeff*coeff)
//                   + a*dM / (2.*coeff*coeff*coeff*h);
//     pcoeff[3] = -dM / (6.*coeff*coeff*coeff*h);

//     if (log(b) < border) {
//       /** - use the series expansion */
//       class_call(array_bessel_series_expansion_integral(pcoeff, n, a, &intg_func_a, internal_error),
//                   internal_error, errmsg);
//       class_call(array_bessel_series_expansion_integral(pcoeff, n, b, &intg_func_b, internal_error),
//                   internal_error, errmsg);

//       prt_sum = intg_func_b - intg_func_a;
//       *result += prt_sum / coeff;
//     }
//     else {
//       /** - compute the bessel functions at a and b up to order n-1 */
//       class_call(spherical_bessel_j(n, a, bessel_j_a, internal_error),
//                   internal_error, errmsg);
//       class_call(spherical_bessel_j(n, b, bessel_j_b, internal_error),
//                   internal_error, errmsg);

//       prt_sum = 0.;
//       prt_res = 0.;
//       for (j = 0; j < n; j++) {
//         /** - compute boundary terms from partial integration */
//         prt_res += array_polynomial(pcoeff, 4, min_pow, a) * bessel_j_a[n-j-1]
//                   -array_polynomial(pcoeff, 4, min_pow, b) * bessel_j_b[n-j-1];
//         /** - compute f_{n-j-1}(x) = f'_{n-j}(x) + (n-j-1)*f_{n-j}(x)/x */
//         min_pow--;
//         for (k = -(j+1); k < 3-j; k++) {
//           pcoeff[k+j+1] *= (n-j+k);
//         }
//       }
//       prt_sum = prt_res;

//       /** - pcoeff now contains f_{0}(x)
//        *    for which we compute the integral dI0 =  dx f_{0}(x) * j_{0}(x) = dx f_{0}(x)/x * sin(x) */
//       min_pow--; /** < makes pcoeff hold f_{0}(x)/x */
//       for (j = min_pow; j < min_pow + 4 && j < 0; j++) {
//         /** - negative powers: use partial integration to reduce to Sine/CosineIntegral function */
//         factor = -1;
//         apow = pow(a, j);
//         bpow = pow(b, j);
//         prt_res = 0.;
//         for (k = 1; k <= -(j+1); k++) {
//           apow *= a; bpow *= b;
//           factor *= -(j+k);
//           prt_res += ( bpow * array_trigonometric_ring(b, k-1)
//                       -apow * array_trigonometric_ring(a, k-1)) / (double)factor;
//         }
//         /* k = -(j+1) => k+j = -1 */
//         factor *= -1;
//         m = (((-(j+1) % 4) + 4) % 4);
//         switch (m)
//         {
//         case 0:
//           class_call(sine_integral(a, &intg_func_a, internal_error),
//                       internal_error, errmsg);
//           class_call(sine_integral(b, &intg_func_b, internal_error),
//                       internal_error, errmsg);
//           break;
//         case 1:
//           class_call(cosine_integral(a, &intg_func_a, internal_error),
//                       internal_error, errmsg);
//           class_call(cosine_integral(b, &intg_func_b, internal_error),
//                       internal_error, errmsg);
//           break;
//         case 2:
//           class_call(sine_integral(a, &intg_func_a, internal_error),
//                       internal_error, errmsg);
//           class_call(sine_integral(b, &intg_func_b, internal_error),
//                       internal_error, errmsg);
//           factor *= -1;
//           break;
//         case 3:
//           class_call(cosine_integral(a, &intg_func_a, internal_error),
//                       internal_error, errmsg);
//           class_call(cosine_integral(b, &intg_func_b, internal_error),
//                       internal_error, errmsg);
//           factor *= -1;
//           break;
//         default: /** can never happen */
//           break;
//         }
//         prt_res += (intg_func_b - intg_func_a) / (double)factor;
//         prt_sum += pcoeff[j-min_pow] * prt_res;

//       }
//       if (j < min_pow + 4) {
//         /** - constant term */
//         prt_res = cos(a) - cos(b);
//         prt_sum += pcoeff[j-min_pow] * prt_res;
//       }
//       j++;
//       for (j; j < min_pow + 4; j++) {
//         /** - positive powers: use partial integration -> sequence terminates automatically
//          *                                                (result consists only of boundary terms) */
//         factor = 1;
//         apow = pow(a, j);
//         bpow = pow(b, j);
//         prt_res =  bpow * array_trigonometric_ring(b, -1)
//                   -apow * array_trigonometric_ring(a, -1);
//         for (k = 1; k <= j; k++) {
//           apow *= 1/a; bpow *= 1/b;
//           factor *= -(j-k+1);
//           prt_res += (double)factor * ( bpow * array_trigonometric_ring(b, -(k+1))
//                               -apow * array_trigonometric_ring(a, -(k+1)) );
//         }

//         prt_sum += pcoeff[j-min_pow] * prt_res;
//       }

//       *result += prt_sum / coeff;
//     }
//   }

//   return _SUCCESS_;
// }
