/**
 * module with tools for manipulating arrays
 * Julien Lesgourgues, 18.04.2010
 * Christian Radermacher, 2023/24
 */

#include "arrays.h"
# include "parallel.h"


/**
 * @brief Internal function for searching supremum and infimum indices
 *        for value in array using bisection
 *
 * @param array         Input: pointer to sorted array indexed as array[i*array_stride] with i < array_size
 * @param array_stride  Input: stride in array
 * @param array_size    Input: size of array
 * @param value         Input: value to search for
 * @param inf_index     In/Output: infimum index (previous must be given)
 * @param sup_index     Output: supremum index
 *
 * @return the error status
 */
int array_search_bisect_internal(
                const double * const __restrict__ array0,
                int array_stride,
                const int array_size,
                const double value,
                int * const __restrict__ inf_index,
                int * const __restrict__ sup_index,
                ErrorMsg errmsg) {

  int inf, sup, mid, reversed = _FALSE_;
  const double * array = array0;

  inf = 0;
  sup = array_size-1;

  if (array[0*array_stride] > array[(array_size-1)*array_stride]) {  /** - descending values */
    /** reverse the array indexing */
    array_stride *= -1;
    array += array_size-1;
    reversed = _TRUE_;
  }

  class_test((value < array[0*array_stride]) || (value > array[(array_size-1)*array_stride]), errmsg, \
            "value=%e is outside of the range of array values (y_min=%e, y_max=%e)",                  \
            value, array[0*array_stride], array[(array_size-1)*array_stride]);

  while (sup-inf > 1) {
    mid = (int)(0.5*(inf+sup));
    if (value < array[mid]) {sup=mid;}
    else {inf=mid;}
  }

  if (reversed) { /** - descending values */
    *inf_index = (array_size-1) - sup;
    *sup_index = (array_size-1) - inf;
  }
  else {          /** - ascending values */
    *inf_index = inf;
    *sup_index = sup;
  }

  return _SUCCESS_;
}

/**
 * @brief Internal function for searching supremum and infimum indices
 *        for value in array in the neighborhood of a previous infimum
 *
 * @param array         Input: pointer to sorted array indexed as array[i*array_stride] with i < array_size
 * @param array_stride  Input: stride in array
 * @param array_size    Input: size of array
 * @param value         Input: value to search for
 * @param inf_index     In/Output: infimum index (previous must be given)
 * @param sup_index     Output: supremum index
 *
 * @return the error status
 */
int array_search_closeby_internal(
                const double * const __restrict__ array0,
                int array_stride,
                const int array_size,
                const double value,
                int * const __restrict__ inf_index,
                int * const __restrict__ sup_index,
                ErrorMsg errmsg) {

  int inf, sup, reversed;
  const double * array = array0;

  if (array[0*array_stride] > array[(array_size-1)*array_stride]) {  /** - descending values */
    /** reverse the array indexing */
    array_stride *= -1;
    array += array_size-1;
    inf = (array_size-1) - (*inf_index + 1);
    reversed = _TRUE_;
  }
  else {                                                             /** - ascending values */
    inf = *inf_index;
    reversed = _FALSE_;
  }

  /** - search for spline interval containing x close to last inf */
  // if (inf<0 || inf>(x_size-1)) return _FAILURE_;  handled by caller
  while (value < array[inf*array_stride]) {
    inf--;
    class_test((inf < 0), errmsg, \
              "value=%e is outside of the range of array values (y_min=%e, y_max=%e)", \
              value, array[0*array_stride], array[(array_size-1)*array_stride]);
  }
  sup = inf + 1;
  while (value > array[sup*array_stride]) {
    sup++;
    class_test((sup > (array_size-1)), errmsg, \
              "value=%e is outside of the range of array values (y_min=%e, y_max=%e)", \
              value, array[0*array_stride], array[(array_size-1)*array_stride]);
  }
  inf = sup - 1;

  if (reversed) { /** - descending values */
    *inf_index = (array_size-1) - sup;
    *sup_index = (array_size-1) - inf;
  }
  else {          /** - ascending values */
    *inf_index = inf;
    *sup_index = sup;
  }

  return _SUCCESS_;
}

/**
 * Called by thermodynamics_init(); perturbations_sources().
 */
int array_derive(
		 double * array,
		 int n_columns,
		 int n_lines,
		 int index_x,   /** from 0 to (n_columns-1) */
		 int index_y,
		 int index_dydx,
		 ErrorMsg errmsg) {

  int i;

  double dx1,dx2,dy1,dy2,weight1,weight2;

  class_test((index_dydx == index_x) || (index_dydx == index_y),
	     errmsg,
	     "output column %d must differ from input columns %d and %d",index_dydx,index_x,index_y);

  dx2=array[1*n_columns+index_x]-array[0*n_columns+index_x];
  dy2=array[1*n_columns+index_y]-array[0*n_columns+index_y];

  for (i=1; i<n_lines-1; i++) {

    dx1 = dx2;
    dy1 = dy2;
    dx2 = array[(i+1)*n_columns+index_x]-array[i*n_columns+index_x];
    dy2 = array[(i+1)*n_columns+index_y]-array[i*n_columns+index_y];
    class_test((dx1 == 0) || (dx2 == 0),
	       errmsg,
	       "stop to avoid division by zero");
    weight1 = dx2*dx2;
    weight2 = dx1*dx1;
    array[i*n_columns+index_dydx] = (weight1*dy1+weight2*dy2) / (weight1*dx1+weight2*dx2);

    if (i == 1)
      array[(i-1)*n_columns+index_dydx] = 2.*dy1/dx1 - array[i*n_columns+index_dydx];

    if (i == n_lines-2)
      array[(i+1)*n_columns+index_dydx] = 2.*dy2/dx2 - array[i*n_columns+index_dydx];
  }

  return _SUCCESS_;
}

int array_derive_spline(
		 double * x_array,
		 int n_lines,
		 double * array,
		 double * array_splined,
		 int n_columns,
		 int index_y,
		 int index_dydx,
		 ErrorMsg errmsg) {

  int i;

  double h;

  class_test(index_dydx == index_y,
	     errmsg,
	     "Output column %d must differ from input columns %d",
	     index_dydx,
	     index_y);

  class_test(n_lines<2,
	     errmsg,
	     "no possible derivation with less than two lines");

  for (i=0; i<n_lines-1; i++) {

    h = x_array[i+1] - x_array[i];
    if (h == 0) {
      class_sprintf(errmsg,"%s(L:%d) h=0, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    array[i*n_columns+index_dydx] =
      (array[(i+1)*n_columns+index_y] - array[i*n_columns+index_y])/h
      - h / 6. * (array_splined[(i+1)*n_columns+index_y] + 2. * array_splined[i*n_columns+index_y]);

  }

  h = x_array[n_lines-1] - x_array[n_lines-2];

  array[(n_lines-1)*n_columns+index_dydx] =
    (array[(n_lines-1)*n_columns+index_y] - array[(n_lines-2)*n_columns+index_y])/h
    + h / 6. * (2. * array_splined[(n_lines-1)*n_columns+index_y] + array_splined[(n_lines-2)*n_columns+index_y]);

  return _SUCCESS_;
}

int array_derive_spline_table_line_to_line(
					   double * x_array,
					   int n_lines,
					   double * array,
					   int n_columns,
					   int index_y,
					   int index_ddy,
					   int index_dy,
					   ErrorMsg errmsg) {

  int i;

  double h;

  class_test(index_ddy == index_y,
	     errmsg,
	     "Output column %d must differ from input columns %d",
	     index_ddy,
	     index_y);

  class_test(index_ddy == index_dy,
	     errmsg,
	     "Output column %d must differ from input columns %d",
	     index_ddy,
	     index_dy);

  class_test(n_lines<2,
	     errmsg,
	     "no possible derivation with less than two lines");

  for (i=0; i<n_lines-1; i++) {

    h = x_array[i+1] - x_array[i];
    if (h == 0) {
      class_sprintf(errmsg,"%s(L:%d) h=0, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    array[i*n_columns+index_dy] =
      (array[(i+1)*n_columns+index_y] - array[i*n_columns+index_y])/h
      - h / 6. * (array[(i+1)*n_columns+index_ddy] + 2. * array[i*n_columns+index_ddy]);

  }

  h = x_array[n_lines-1] - x_array[n_lines-2];

  array[(n_lines-1)*n_columns+index_dy] =
    (array[(n_lines-1)*n_columns+index_y] - array[(n_lines-2)*n_columns+index_y])/h
    + h / 6. * (2. * array[(n_lines-1)*n_columns+index_ddy] + array[(n_lines-2)*n_columns+index_ddy]);

  return _SUCCESS_;
}

int array_derive1_order2_table_line_to_line(
				       double * x_array,
				       int n_lines,
				       double * array,
				       int n_columns,
				       int index_y,
				       int index_dy,
				       ErrorMsg errmsg) {

  int i=1;
  double dxp,dxm,dyp,dym;

  if (n_lines < 2) {
    class_sprintf(errmsg,"%s(L:%d) routine called with n_lines=%d, should be at least 2",__func__,__LINE__,n_lines);
    return _FAILURE_;
  }

  dxp = x_array[2] - x_array[1];
  dxm = x_array[0] - x_array[1];
  dyp = *(array+2*n_columns+index_y) - *(array+1*n_columns+index_y);
  dym = *(array+0*n_columns+index_y) - *(array+1*n_columns+index_y);

  if ((dxp*dxm*(dxm-dxp)) == 0.) {
    class_sprintf(errmsg,"%s(L:%d) stop to avoid division by zero",__func__,__LINE__);
    return _FAILURE_;
  }

  *(array+1*n_columns+index_dy) = (dyp*dxm*dxm-dym*dxp*dxp)/(dxp*dxm*(dxm-dxp));

  *(array+0*n_columns+index_dy) = *(array+1*n_columns+index_dy)
    - (x_array[1] - x_array[0]) * 2.*(dyp*dxm-dym*dxp)/(dxp*dxm*(dxp-dxm));

  for (i=2; i<n_lines-1; i++) {

    dxp = x_array[i+1] - x_array[i];
    dxm = x_array[i-1] - x_array[i];
    dyp = *(array+(i+1)*n_columns+index_y) - *(array+i*n_columns+index_y);
    dym = *(array+(i-1)*n_columns+index_y) - *(array+i*n_columns+index_y);

    if ((dxp*dxm*(dxm-dxp)) == 0.) {
      class_sprintf(errmsg,"%s(L:%d) stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    *(array+i*n_columns+index_dy) = (dyp*dxm*dxm-dym*dxp*dxp)/(dxp*dxm*(dxm-dxp));

  }

  *(array+(n_lines-1)*n_columns+index_dy) = *(array+(n_lines-2)*n_columns+index_dy)
    + (x_array[n_lines-1] - x_array[n_lines-2]) * 2.*(dyp*dxm-dym*dxp)/(dxp*dxm*(dxp-dxm));

  return _SUCCESS_;

}

int array_derive2_order2_table_line_to_line(
				       double * x_array,
				       int n_lines,
				       double * array,
				       int n_columns,
				       int index_y,
				       int index_dy,
				       int index_ddy,
				       ErrorMsg errmsg) {

  int i;
  double dxp,dxm,dyp,dym;

  for (i=1; i<n_lines-1; i++) {

    dxp = x_array[i+1] - x_array[i];
    dxm = x_array[i-1] - x_array[i];
    dyp = *(array+(i+1)*n_columns+index_y) - *(array+i*n_columns+index_y);
    dym = *(array+(i-1)*n_columns+index_y) - *(array+i*n_columns+index_y);

    if ((dxp*dxm*(dxm-dxp)) == 0.) {
      class_sprintf(errmsg,"%s(L:%d) stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    *(array+i*n_columns+index_dy) = (dyp*dxm*dxm-dym*dxp*dxp)/(dxp*dxm*(dxm-dxp));
    *(array+i*n_columns+index_ddy) = 2.*(dyp*dxm-dym*dxp)/(dxp*dxm*(dxp-dxm));

  }

  *(array+0*n_columns+index_dy) = *(array+1*n_columns+index_dy)
    - (x_array[1] - x_array[0]) * *(array+1*n_columns+index_ddy);
  *(array+0*n_columns+index_ddy) = *(array+1*n_columns+index_ddy);

  *(array+(n_lines-1)*n_columns+index_dy) = *(array+(n_lines-2)*n_columns+index_dy)
    + (x_array[n_lines-1] - x_array[n_lines-2]) * *(array+(n_lines-2)*n_columns+index_ddy);
  *(array+(n_lines-1)*n_columns+index_ddy) = *(array+(n_lines-2)*n_columns+index_ddy);

  return _SUCCESS_;

}

int array_integrate_spline_table_line_to_line(
					      double * x_array,
					      int n_lines,
					      double * array,
					      int n_columns,
					      int index_y,
					      int index_ddy,
					      int index_inty,
					      ErrorMsg errmsg) {

  int i;

  double h;

  *(array+0*n_columns+index_inty)  = 0.;

  for (i=0; i < n_lines-1; i++) {

    h = (x_array[i+1]-x_array[i]);

    *(array+(i+1)*n_columns+index_inty) = *(array+i*n_columns+index_inty) +
      (array[i*n_columns+index_y]+array[(i+1)*n_columns+index_y])*h/2. -
      (array[i*n_columns+index_ddy]+array[(i+1)*n_columns+index_ddy])*h*h*h/24.;

  }

  return _SUCCESS_;
}


 /**
 * Not called.
 */
int array_derive_two(
		     double * array,
		     int n_columns,
		     int n_lines,
		     int index_x,   /** from 0 to (n_columns-1) */
		     int index_y,
		     int index_dydx,
		     int index_ddydxdx,
		     ErrorMsg errmsg) {

  int i;

  double dx1,dx2,dy1,dy2,weight1,weight2;

  if ((index_dydx == index_x) || (index_dydx == index_y)) {
    class_sprintf(errmsg,"%s(L:%d) : Output column %d must differ from input columns %d and %d",__func__,__LINE__,index_dydx,index_x,index_y);
    return _FAILURE_;
  }

  dx2=*(array+1*n_columns+index_x)-*(array+0*n_columns+index_x);
  dy2=*(array+1*n_columns+index_y)-*(array+0*n_columns+index_y);

  for (i=1; i<n_lines-1; i++) {

    dx1 = dx2;
    dy1 = dy2;
    dx2 = *(array+(i+1)*n_columns+index_x)-*(array+i*n_columns+index_x);
    dy2 = *(array+(i+1)*n_columns+index_y)-*(array+i*n_columns+index_y);
    weight1 = dx2*dx2;
    weight2 = dx1*dx1;

    if ((dx1 == 0.) && (dx2 == 0.)) {
      class_sprintf(errmsg,"%s(L:%d) stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    *(array+i*n_columns+index_dydx) = (weight1*dy1+weight2*dy2) / (weight1*dx1+weight2*dx2);
    *(array+i*n_columns+index_ddydxdx) = (dx2*dy1-dx1*dy2) / (weight1*dx1+weight2*dx2);

    if (i == 1) {
      *(array+(i-1)*n_columns+index_dydx) = 2.*dy1/dx1 - *(array+i*n_columns+index_dydx);
      *(array+(i-1)*n_columns+index_ddydxdx) = *(array+i*n_columns+index_ddydxdx);
    }

    if (i == n_lines-2) {
      *(array+(i+1)*n_columns+index_dydx) = 2.*dy2/dx2 - *(array+i*n_columns+index_dydx);
      *(array+(i+1)*n_columns+index_dydx) = *(array+i*n_columns+index_ddydxdx);
    }
  }

  return _SUCCESS_;
}


/**
 * @brief Thomas algorithm for tridiagonal systems of linear equations.
 *        Stable if the coefficient matrix is diagonally dominant or symmetric positive-definite.
 *
 * @param size          Input: size of the coefficient matrix [size*size = (n+1)*(n+1)]
 * @param diag          Input: array of diagonal elements [size] with elements mu_i = diag[i]
 * @param superdiag     Input: array of super-diagonal elements [size-1] with elements lambda_i = superdiag[i]
 * @param subdiag       Input: array of sub-diagonal elements [size-1] with elements kappa_i = subdiag[i-1]
 * @param constant      Input: constant vector [size] with elements b_i = constant[i]
 * @param solution      Output: solution vector [size] with elements x_i = solution[i]; externally allocated
 * @param superdiag_f   In/Output: array of forward super-diagonal elements [size-1] with elements lambda'_i = superdiag_f[i]; externally allocated
 * @param compute_forward_superdiagonal   Input: overwrite the forward super-diagonal elements in superdiag_f
 * @param constant_f   Output: forward constant vector [size] with elements b'_i = constant_f[i]; externally allocated
 *
 * @return the error status
 */
int arrays_tridiagonal_solve(const int size,
                            const double * const diag,
                            const double * const superdiag,
                            const double * const subdiag,
                            const double * const constant,
                            double * solution,
                            double * superdiag_f,
                            const short compute_forward_superdiagonal,
                            double * constant_f) {

  int it;
  double denom;

  /** - forward sweep */
  if (compute_forward_superdiagonal) {  /** - compute the forward-substituted superdiagonal coefficients */
    superdiag_f[0] = superdiag[0] / diag[0];  // lambda'_0 = lambda_0 / mu_0
    constant_f[0]  = constant[0] / diag[0];   // b'_0 = b_0 / mu0
    for (it=1; it < size-1; it++) { // from 1 to n-1
      denom = diag[it] - subdiag[it-1] * superdiag_f[it-1]; // (mu_i - kappa_i * lambda'_{i-1})
      superdiag_f[it] = superdiag[it] / denom;  // lambda'_i = lambda_i / (mu_i - kappa_i * lambda'_{i-1})
      constant_f[it]  = (constant[it] - subdiag[it-1] * constant_f[it-1]) / denom;  // b'_i = (b_i - kappa_i * b'_{i-1})/(mu_i - kappa_i * lambda'_{i-1})
    }
    denom = diag[size-1] - subdiag[size-2] * superdiag_f[size-2]; // (mu_n - kappa_n * lambda'_{n-1})
    constant_f[size-1] = (constant[size-1] - subdiag[size-2] * constant_f[size-2]) / denom; // b'_n = (b_n - kappa_n * b'_{n-1})/(mu_n - kappa_n * lambda'_{n-1})
  }
  else {  /** - use predefined forward-substituted superdiagonal coefficients */
    constant_f[0]  = constant[0] / diag[0];   // b'_0 = b_0 / mu0
    for (it=1; it < size-1; it++) { // from 1 to n-1
      denom = diag[it] - subdiag[it-1] * superdiag_f[it-1]; // (mu_i - kappa_i * lambda'_{i-1})
      constant_f[it]  = (constant[it] - subdiag[it-1] * constant_f[it-1]) / denom;  // b'_i = (b_i - kappa_i * b'_{i-1})/(mu_i - kappa_i * lambda'_{i-1})
    }
    denom = diag[size-1] - subdiag[size-2] * superdiag_f[size-2]; // (mu_n - kappa_n * lambda'_{n-1})
    constant_f[size-1] = (constant[size-1] - subdiag[size-2] * constant_f[size-2]) / denom; // b'_n = (b_n - kappa_n * b'_{n-1})/(mu_n - kappa_n * lambda'_{n-1})
  }

  /** - backward substitution */
  solution[size-1] = constant_f[size-1];  // x_n = b'_n
  for (it=size-2; it >= 0; it--) {
    solution[it] = constant_f[it] - superdiag_f[it] * solution[it+1]; // x_i = b'_i - lambda'_i * x_{i+1}
  }

  return _SUCCESS_;
}


#ifdef __FAST_MATH__
#warning Kahan-Neumaier summation is incompatible with -ffast-math. Defaulting to uncompensated summation.

static inline double d2sum(const double a, const double b, double *__restrict__ t) {
  double s;
  s = a + b;
  *t = 0.;
  return s;
}

static inline double dcompsum(const double * const summands, const int size, const int stride) {
  int i;
  double s;

  s = 0.;
  for (i = 0; i < size; i++) {
    s += summands[i];
  }
  return s;
}

#else
static inline double d2sum(const double a, const double b, double * __restrict__ t) {
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
static inline double dcompsum(const double * const summands, const int size, const int stride) {
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
#endif

/**
 * @brief Perform internal splining for the given input arrays (can be made identical) which are indexed as array[i*array_stride].
 *        Needs pre-allocated working arrays for superdiagonal and constant entries in the tridiagonal system of equations.
 *        Implements the Thomas algorithm which is stable for all valid splines.
 *        Sets natural boundary conditions S''(x0) = S''(xn) = 0.
 *
 * @param x0          Input: pointer to first x value
 * @param x_stride    Input: stride in the x-array to get to the next value
 * @param y0          Input: pointer to first y-value
 * @param y_stride    Input: stride in the y-array
 * @param ddy0        Output: pointer to first moment of the spline
 * @param ddy_stride  Input: stride in the ddy-array
 * @param size        Input: number of equations to solve (#points)
 * @param super       Output: iterative super-diagonal coefficients (externally allocated: at least (size-1)*sizeof(double))
 * @param constants   Output: iterative constant terms / RHS (externally allocated: at least (size-1)*sizeof(double))
 *
 * @return the error status
 */
int array_spline_internal_natural(
              double * x0,
              int x_stride,
              double * y0,
              int y_stride,
              double * ddy0,
              int ddy_stride,
              int size,
              double * super,
              double * constants
              ) {

  int it;
  double h1, h2;
  double diag, denom; /**< 3 * current diagonal entry, 3 * current denominator */

  /** - set boundary condition at x0 */
  super[0] = 0.; // lambda'0 = lambda0/mu0
  constants[0] = 0.; // b0/mu0

  /** - forward sweep */
  for (it=1; it < size-1; it++) {
    h1 = x0[it*x_stride] - x0[(it-1)*x_stride]; // h_i
    h2 = x0[(it+1)*x_stride] - x0[it*x_stride]; // h_{i+1}
    diag = h1+h2; // 3*mu_i
    denom = diag - 0.5*h1*super[it-1];  // 3*(mu_i - kappa_i * lambda'_{i-1})
    super[it] = 0.5*h2/denom; // lambda'_i = lambda_i / (mu_i - kappa_i * lambda'_{i-1})
    constants[it] = (3.*((y0[(it+1)*y_stride]-y0[it*y_stride])/h2   \
                        -(y0[it*y_stride]-y0[(it-1)*y_stride])/h1)  \
                    -0.5*h1*constants[it-1]) / denom; // b'_i = (b_i - kappa_i * b'_{i-1})/(mu_i - kappa_i * lambda'_{i-1})
  }

  /** - set boundary condition at x_n */
  ddy0[(size-1)*ddy_stride] = 0.; // M_n = b'_{n} = 0 since b_n = kappa_n = 0 and mu_n = 1

  /** - backward substitution */
  for (it=size-2; it >= 0; it--) {
    ddy0[it*ddy_stride] = constants[it] - super[it] * ddy0[(it+1)*ddy_stride];
                        // M_i = b'_i - lambda'_i * M_{i+1}
  }

  return _SUCCESS_;
}

/**
 * @brief Perform internal splining for the given input arrays (can be made identical) which are indexed as array[i*array_stride].
 *        Needs pre-allocated working arrays for superdiagonal and constant entries in the tridiagonal system of equations.
 *        Implements the Thomas algorithm which is stable for all valid splines.
 *        Sets hermite boundary conditions S'(x0) = dy_first and S'(xn) = dy_last which can be approximated from the input.
 *
 * @param x0          Input: pointer to first x value
 * @param x_stride    Input: stride in the x-array to get to the next value
 * @param y0          Input: pointer to first y-value
 * @param y_stride    Input: stride in the y-array
 * @param ddy0        Output: pointer to first moment of the spline
 * @param ddy_stride  Input: stride in the ddy-array
 * @param size        Input: number of equations to solve (#points)
 * @param super       Output: iterative super-diagonal coefficients (externally allocated: at least (size-1)*sizeof(double))
 * @param constants   Output: iterative constant terms / RHS (externally allocated: at least (size-1)*sizeof(double))
 * @param use_approx  Input: use 3-point approximation of the derivative at the boundaries
 * @param dy_first    In/Output: derivative of y at x0
 * @param dy_last     In/Output: derivative of y at xn
 *
 * @return the error status
 */
int array_spline_internal_hermite(
              double * x0,
              int x_stride,
              double * y0,
              int y_stride,
              double * ddy0,
              int ddy_stride,
              int size,
              double * super,
              double * constants,
              short use_approx,
              double * dy_first,
              double * dy_last
              ) {

  int it;
  double h1, h2;
  double diag, denom; /**< 3 * current diagonal entry, 3 * current denominator */

  /** - compute the approximate derivative at x0 using a 3-point rule */
  if (use_approx)
  {
    h1 = x0[1*x_stride] - x0[0*x_stride]; // h1
    h2 = x0[2*x_stride] - x0[1*x_stride]; // h2
    *dy_first = ((h1+h2)*(h1+h2)*(y0[1*y_stride]-y0[0*y_stride])    \
                          -h1*h1*(y0[2*y_stride]-y0[0*y_stride]))   \
                /((h1+h2)*h1*h2);

    h1 = x0[(size-2)*x_stride] - x0[(size-3)*x_stride]; // h_{n-1}
    h2 = x0[(size-1)*x_stride] - x0[(size-2)*x_stride]; // h_n
    *dy_last  = ((h1+h2)*(h1+h2)*(y0[(size-1)*y_stride]-y0[(size-2)*y_stride])   \
                          -h2*h2*(y0[(size-1)*y_stride]-y0[(size-3)*y_stride]))  \
                /((h1+h2)*h1*h2);
  }

  /** - set boundary condition at x0 */
  h1 = x0[1*x_stride] - x0[0*x_stride]; // h1
  super[0] = 0.5; // lambda0/mu0
  constants[0] = (3./h1)*((y0[1*y_stride]-y0[0*y_stride])/h1  \
                            - *dy_first); // b0/mu0

  /** - forward sweep */
  for (it=1; it < size-1; it++) {
    h1 = x0[it*x_stride] - x0[(it-1)*x_stride]; // h_i
    h2 = x0[(it+1)*x_stride] - x0[it*x_stride]; // h_{i+1}
    diag = h1+h2; // 3*mu_i
    denom = diag - 0.5*h1*super[it-1];  // 3*(mu_i - kappa_i * lambda'_{i-1})
    super[it] = 0.5*h2/denom; // lambda'_i = lambda_i / (mu_i - kappa_i * lambda'_{i-1})
    constants[it] = (3.*((y0[(it+1)*y_stride]-y0[it*y_stride])/h2   \
                        -(y0[it*y_stride]-y0[(it-1)*y_stride])/h1)  \
                    -0.5*h1*constants[it-1]) / denom; // b'_i = (b_i - kappa_i * b'_{i-1})/(mu_i - kappa_i * lambda'_{i-1})
  }

  /** - set boundary condition at x_n */
  // h2 is still h_n
  ddy0[(size-1)*ddy_stride] = ((3./h2)*(*dy_last - (y0[(size-1)*y_stride]-y0[(size-2)*y_stride])/h2)   \
                                -0.5*constants[size-2]) / (1. - 0.5*super[size-2]); // M_n = b'_n = (b_n - kappa_n * b'_{n-1})/(mu_n - kappa_n * lambda'_{n-1})

  /** - backward substitution */
  for (it=size-2; it >= 0; it--) {
    ddy0[it*ddy_stride] = constants[it] - super[it] * ddy0[(it+1)*ddy_stride];
                        // M_i = b'_i - lambda'_i * M_{i+1}
  }

  return _SUCCESS_;
}

/**
 * @brief Perform internal splining for the given input arrays (can be made identical) which are indexed as array[i*array_stride].
 *        Needs pre-allocated working arrays for superdiagonal and constant entries in the tridiagonal system of equations.
 *        Implements the Thomas algorithm which is stable for all valid splines.
 *        Sets hermite boundary conditions S'(x0) = dy_first and S'(xn) = dy_last which can be approximated from the input.
 *
 * @param x0            Input: pointer to first x value
 * @param x_stride      Input: stride in the x-array to get to the next value
 * @param y0            Input: pointer to first y-value
 * @param y_stride      Input: stride in the y-array
 * @param ddy0          Output: pointer to first moment of the spline
 * @param ddy_stride    Input: stride in the ddy-array
 * @param size          Input: number of equations to solve (#points - 1)
 * @param super         Output: iterative super-diagonal coefficients (externally allocated: at least (size-1)*sizeof(double))
 * @param constants     Output: iterative constant terms / RHS (externally allocated: at least (size-1)*sizeof(double))
 * @param constants_aux Output: iterative constant terms / RHS for auxiliary problem (externally allocated: at least (size-1)*sizeof(double))
 * @param sol_aux       Output: solution of the auxiliary problem ~A.sol = u (externally allocated: at least size*sizeof(double))
 *
 * @return the error status
 */
int array_spline_internal_periodic(
              double * x0,
              int x_stride,
              double * y0,
              int y_stride,
              double * ddy0,
              int ddy_stride,
              int size,
              double * super,
              double * constants,
              double * constants_aux,
              double * sol_aux
              ) {

  int it;
  double h1, h2;
  double diag, denom; /**< 3 * current diagonal entry, 3 * current denominator */

  /** - set boundary condition at x0 */
  h1 = x0[1*x_stride] - x0[0*x_stride]; // h1
  h2 = x0[size*x_stride] - x0[(size-1)*x_stride]; // hn
  diag = h1 + 0.5*h2; // 3*mu0
  super[0] = 0.5*h1/diag; // lambda'0 = lambda0/mu0
  constants[0] = (3./diag) *((y0[1*y_stride]-y0[0*y_stride])/h1  \
                            -(y0[size*y_stride]-y0[(size-1)*y_stride])/h2); // b0/mu0
  constants_aux[0] = 0.5*h2/diag; // u'0 = u0/mu0

  /** - forward sweep */
  for (it=1; it < size-1; it++) {
    h1 = x0[it*x_stride] - x0[(it-1)*x_stride]; // h_i
    h2 = x0[(it+1)*x_stride] - x0[it*x_stride]; // h_{i+1}
    diag = h1+h2; // 3*mu_i
    denom = diag - 0.5*h1*super[it-1];  // 3*(mu_i - kappa_i * lambda'_{i-1})
    super[it] = 0.5*h2/denom; // lambda'_i = lambda_i / (mu_i - kappa_i * lambda'_{i-1})
    constants[it] = (3.*((y0[(it+1)*y_stride]-y0[it*y_stride])/h2   \
                        -(y0[it*y_stride]-y0[(it-1)*y_stride])/h1)  \
                    -0.5*h1*constants[it-1]) / denom; // b'_i = (b_i - kappa_i * b'_{i-1})/(mu_i - kappa_i * lambda'_{i-1})
    constants_aux[it] = -0.5*h1*constants_aux[it-1] / denom;  // u'_i = - kappa_i * u'_{i-1}/(mu_i - kappa_i * lambda'_{i-1})
  }

  /** - set boundary condition at x_{n-1} */
  h1 = x0[(size-1)*x_stride] - x0[(size-2)*x_stride]; // h_{n-1}
  h2 = x0[size*x_stride] - x0[(size-1)*x_stride]; // h_n
  diag = h1 + 0.5*h2; // 3*mu_{n-1}
  denom = diag - 0.5*h1*super[size-2];  // 3*(mu_{n-1} - kappa_{n-1} * lambda'_{n-2})
  ddy0[(size-1)*ddy_stride] = (3.*((y0[size*y_stride]-y0[(size-1)*y_stride])/h2      \
                                       -(y0[(size-1)*y_stride]-y0[(size-2)*y_stride])/h1) \
                                    -0.5*h1*constants[size-2]) / denom; // m_{n-1} = b'_{n-1} = (b_{n-1} - kappa_{n-1} * b'_{n-2})/(mu_{n-1} - kappa_{n-1} * lambda'_{n-2})
  sol_aux[size-1] = 0.5*(h2 - h1*constants_aux[size-2]) / denom;  // z_{n-1} = u'_{n-1} = (u_{n-1} - kappa_{n-1} * u'_{n-2})/(mu_{n-1} - kappa_{n-1} * lambda'_{n-2})

  /** - backward substitution */
  for (it=size-2; it >= 0; it--) {
    ddy0[it*ddy_stride] = constants[it] - super[it] * ddy0[(it+1)*ddy_stride];
                        // m_i = b'_i - lambda'_i * m_{i+1}
    sol_aux[it] = constants_aux[it] - super[it] * sol_aux[it+1];
                        // z_i = u'_i - lambda'_i * z_{i+1}
  }

  /** - compute scalar products */
  diag = ddy0[0*ddy_stride] + ddy0[(size-1)*ddy_stride];  // v.m = m_0 + m_{n-1}
  denom = 1. + sol_aux[0] + sol_aux[size-1];  // 1 + v.z = 1 + z_0 + z_{n-1}
  /** - correct output according to Sherman-Morrison formula */
  for (it=0; it < size; it++) {
    ddy0[it*ddy_stride] -= (diag/denom) * sol_aux[it];
  }

  /** - last moment is the same as the first one */
  ddy0[size*ddy_stride] = ddy0[0*ddy_stride];


  return _SUCCESS_;
}

int array_spline(
		  double * array,
		  int n_columns,
		  int n_lines,
		  int index_x,   /** from 0 to (n_columns-1) */
		  int index_y,
		  int index_ddydx2,
		  short spline_mode,
		  ErrorMsg errmsg) {

  double *super, *constants;
  double *constants_aux, *sol_aux;

  class_test(n_lines < 3, errmsg, "%s(L:%d) there is no spline with less than 3 points", __func__, __LINE__);

  switch (spline_mode)
  {
  case _SPLINE_NATURAL_:
    class_alloc(super, (n_lines-1)*sizeof(double), errmsg);
    class_alloc(constants, (n_lines-1)*sizeof(double), errmsg);

    array_spline_internal_natural(array+index_x, n_columns, array+index_y, n_columns, array+index_ddydx2, n_columns, \
                                  n_lines, super, constants);
    break;

  case _SPLINE_EST_DERIV_:
    class_alloc(super, (n_lines-1)*sizeof(double), errmsg);
    class_alloc(constants, (n_lines-1)*sizeof(double), errmsg);
    double dy_first, dy_last;

    array_spline_internal_hermite(array+index_x, n_columns, array+index_y, n_columns, array+index_ddydx2, n_columns, \
                                  n_lines, super, constants, _TRUE_, &dy_first, &dy_last);
    break;

  case _SPLINE_PERIODIC_:
    class_alloc(super, (n_lines-2)*sizeof(double), errmsg);
    class_alloc(constants, (n_lines-2)*sizeof(double), errmsg);
    class_alloc(constants_aux, (n_lines-2)*sizeof(double), errmsg);
    class_alloc(sol_aux, (n_lines-1)*sizeof(double), errmsg);

    array_spline_internal_periodic(array+index_x, n_columns, array+index_y, n_columns, array+index_ddydx2, n_columns, \
                                    n_lines-1, super, constants, constants_aux, sol_aux);

    free(constants_aux);
    free(sol_aux);
    break;

  default:
    class_stop(errmsg, "%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
    break;
  }

  free(super);
  free(constants);

  return _SUCCESS_;
}

int array_spline_table_line_to_line(
                double * x, /* vector of size x_size */
                int n_lines,
                double * array,
                int n_columns,
                int index_y,
                int index_ddydx2,
                short spline_mode,
                ErrorMsg errmsg
                ) {

  double *super, *constants;
  double *constants_aux, *sol_aux;

  class_test(n_lines < 3, errmsg, "%s(L:%d) there is no spline with less than 3 points", __func__, __LINE__);

  switch (spline_mode)
  {
  case _SPLINE_NATURAL_:
    class_alloc(super, (n_lines-1)*sizeof(double), errmsg);
    class_alloc(constants, (n_lines-1)*sizeof(double), errmsg);

    array_spline_internal_natural(x, 1, array+index_y, n_columns, array+index_ddydx2, n_columns, \
                                  n_lines, super, constants);
    break;

  case _SPLINE_EST_DERIV_:
    class_alloc(super, (n_lines-1)*sizeof(double), errmsg);
    class_alloc(constants, (n_lines-1)*sizeof(double), errmsg);
    double dy_first, dy_last;

    array_spline_internal_hermite(x, 1, array+index_y, n_columns, array+index_ddydx2, n_columns, \
                                  n_lines, super, constants, _TRUE_, &dy_first, &dy_last);
    break;

  case _SPLINE_PERIODIC_:
    class_alloc(super, (n_lines-2)*sizeof(double), errmsg);
    class_alloc(constants, (n_lines-2)*sizeof(double), errmsg);
    class_alloc(constants_aux, (n_lines-2)*sizeof(double), errmsg);
    class_alloc(sol_aux, (n_lines-1)*sizeof(double), errmsg);

    array_spline_internal_periodic(x, 1, array+index_y, n_columns, array+index_ddydx2, n_columns, \
                                    n_lines-1, super, constants, constants_aux, sol_aux);

    free(constants_aux);
    free(sol_aux);
    break;

  default:
    class_stop(errmsg, "%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
    break;
  }

  free(super);
  free(constants);

  return _SUCCESS_;
}

int array_spline_table_lines(
			     double * x, /* vector of size x_size */
			     int x_size,
			     double * y_array, /* array of size x_size*y_size with elements
						  y_array[index_x*y_size+index_y] */
			     int y_size,
			     double * ddy_array, /* array of size x_size*y_size */
			     short spline_mode,
			     ErrorMsg errmsg
			     ) {

  int index_y;
  double *super, *constants;
  double *constants_aux, *sol_aux;

  class_test(x_size < 3, errmsg, "%s(L:%d) there is no spline with less than 3 points", __func__, __LINE__);

  switch (spline_mode)
  {
  case _SPLINE_NATURAL_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_natural(x, 1, y_array+index_y, y_size, ddy_array+index_y, y_size, x_size, super, constants);
    }
    break;

  case _SPLINE_EST_DERIV_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);
    double dy_first, dy_last;

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_hermite(x, 1, y_array+index_y, y_size, ddy_array+index_y, y_size, x_size, super, constants, \
                                    _TRUE_, &dy_first, &dy_last);
    }
    break;

  case _SPLINE_PERIODIC_:
    class_alloc(super, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants_aux, (x_size-2)*sizeof(double), errmsg);
    class_alloc(sol_aux, (x_size-1)*sizeof(double), errmsg);

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_periodic(x, 1, y_array+index_y, y_size, ddy_array+index_y, y_size, x_size-1, super, constants, \
                                      constants_aux, sol_aux);
    }

    free(constants_aux);
    free(sol_aux);
    break;

  default:
    class_stop(errmsg, "%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
    break;
  }

  free(super);
  free(constants);

  return _SUCCESS_;
}

int array_spline_table_lines_parallel(
			     double * x, /* vector of size x_size */
			     int x_size,
			     double * y_array, /* array of size x_size*y_size with elements
						  y_array[index_x*y_size+index_y] */
			     int y_size,
			     double * ddy_array, /* array of size x_size*y_size */
			     short spline_mode,
			     ErrorMsg errmsg
			     ) {

  int index_y;
  double *super, *constants;
  double *constants_aux, *sol_aux;

  class_test(x_size < 3, errmsg, "%s(L:%d) there is no spline with less than 3 points", __func__, __LINE__);

  switch (spline_mode)
  {
  case _SPLINE_NATURAL_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);
    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_natural(x, 1, y_array+index_y, y_size, ddy_array+index_y, y_size, x_size, super, constants);
    }

    free(super);
    free(constants);
    break;

  case _SPLINE_EST_DERIV_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);
    double dy_first, dy_last;

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_hermite(x, 1, y_array+index_y, y_size, ddy_array+index_y, y_size, x_size, super, constants, \
                                    _TRUE_, &dy_first, &dy_last);
    }

    free(super);
    free(constants);
    break;

  case _SPLINE_PERIODIC_:
    class_alloc(super, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants_aux, (x_size-2)*sizeof(double), errmsg);
    class_alloc(sol_aux, (x_size-1)*sizeof(double), errmsg);

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_periodic(x, 1, y_array+index_y, y_size, ddy_array+index_y, y_size, x_size-1, super, constants, \
                                     constants_aux, sol_aux);
    }

    free(super);
    free(constants);
    free(constants_aux);
    free(sol_aux);
    break;

  default:
    ErrorMsg errmsg_mode;
    class_sprintf(errmsg_mode, "Spline mode not identified: %d", spline_mode);
    class_build_error_string(errmsg, "error; %s", errmsg_mode);
    break;
  }

  return _SUCCESS_;
}

int array_spline_table_columns(
		       double * x, /* vector of size x_size */
		       int x_size,
		       double * y_array, /* array of size x_size*y_size with elements
					  y_array[index_y*x_size+index_x] */
		       int y_size,
		       double * ddy_array, /* array of size x_size*y_size */
		       short spline_mode,
		       ErrorMsg errmsg
		       ) {

  int index_y;
  double *super, *constants;
  double *constants_aux, *sol_aux;

  class_test(x_size < 3, errmsg, "%s(L:%d) there is no spline with less than 3 points", __func__, __LINE__);

  switch (spline_mode)
  {
  case _SPLINE_NATURAL_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_natural(x, 1, y_array+index_y*x_size, 1, ddy_array+index_y*x_size, 1, x_size, super, constants);
    }
    break;

  case _SPLINE_EST_DERIV_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);
    double dy_first, dy_last;

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_hermite(x, 1, y_array+index_y*x_size, 1, ddy_array+index_y*x_size, 1, x_size, super, constants, \
                                    _TRUE_, &dy_first, &dy_last);
    }
    break;

  case _SPLINE_PERIODIC_:
    class_alloc(super, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants_aux, (x_size-2)*sizeof(double), errmsg);
    class_alloc(sol_aux, (x_size-1)*sizeof(double), errmsg);

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_periodic(x, 1, y_array+index_y*x_size, 1, ddy_array+index_y*x_size, 1, x_size-1, super, constants, \
                                      constants_aux, sol_aux);
    }

    free(constants_aux);
    free(sol_aux);
    break;

  default:
    class_stop(errmsg, "%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
    break;
  }

  free(super);
  free(constants);

  return _SUCCESS_;
}

int array_spline_table_columns_parallel(
			     double * x, /* vector of size x_size */
		       int x_size,
		       double * y_array, /* array of size x_size*y_size with elements
					  y_array[index_y*x_size+index_x] */
		       int y_size,
		       double * ddy_array, /* array of size x_size*y_size */
		       short spline_mode,
		       ErrorMsg errmsg
		       ) {

  int index_y;
  double *super, *constants;
  double *constants_aux, *sol_aux;

  class_test(x_size < 3, errmsg, "%s(L:%d) there is no spline with less than 3 points", __func__, __LINE__);

  switch (spline_mode)
  {
  case _SPLINE_NATURAL_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_natural(x, 1, y_array+index_y*x_size, 1, ddy_array+index_y*x_size, 1, x_size, super, constants);
    }

    free(super);
    free(constants);
    break;

  case _SPLINE_EST_DERIV_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);
    double dy_first, dy_last;

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_hermite(x, 1, y_array+index_y*x_size, 1, ddy_array+index_y*x_size, 1, x_size, super, constants, \
                                    _TRUE_, &dy_first, &dy_last);
    }

    free(super);
    free(constants);
    break;

  case _SPLINE_PERIODIC_:
    class_alloc(super, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants_aux, (x_size-2)*sizeof(double), errmsg);
    class_alloc(sol_aux, (x_size-1)*sizeof(double), errmsg);

    for (index_y=0; index_y < y_size; index_y++) {
      array_spline_internal_periodic(x, 1, y_array+index_y*x_size, 1, ddy_array+index_y*x_size, 1, x_size-1, super, constants, \
                                     constants_aux, sol_aux);
    }

    free(super);
    free(constants);
    free(constants_aux);
    free(sol_aux);
    break;

  default:
    ErrorMsg errmsg_mode;
    class_sprintf(errmsg_mode, "Spline mode not identified: %d", spline_mode);
    class_build_error_string(errmsg, "error; %s", errmsg_mode);
    break;
  }

  return _SUCCESS_;
}

int array_spline_table_one_column(
		       double * x, /* vector of size x_size */
		       int x_size,
		       double * y_array, /* array of size x_size*y_size with elements
					  y_array[index_y*x_size+index_x] */
		       int index_y,
		       double * ddy_array, /* array of size x_size*y_size */
		       short spline_mode,
		       ErrorMsg errmsg
		       ) {

  double *super, *constants;
  double *constants_aux, *sol_aux;

  class_test(x_size < 3, errmsg, "%s(L:%d) there is no spline with less than 3 points", __func__, __LINE__);

  switch (spline_mode)
  {
  case _SPLINE_NATURAL_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);

    array_spline_internal_natural(x, 1, y_array+index_y*x_size, 1, ddy_array+index_y*x_size, 1, x_size, super, constants);
    break;

  case _SPLINE_EST_DERIV_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);
    double dy_first, dy_last;

    array_spline_internal_hermite(x, 1, y_array+index_y*x_size, 1, ddy_array+index_y*x_size, 1, x_size, super, constants, \
                                  _TRUE_, &dy_first, &dy_last);
    break;

  case _SPLINE_PERIODIC_:
    class_alloc(super, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants_aux, (x_size-2)*sizeof(double), errmsg);
    class_alloc(sol_aux, (x_size-1)*sizeof(double), errmsg);

    array_spline_internal_periodic(x, 1, y_array+index_y*x_size, 1, ddy_array+index_y*x_size, 1, x_size-1, super, constants, \
                                    constants_aux, sol_aux);

    free(constants_aux);
    free(sol_aux);
    break;

  default:
    class_stop(errmsg, "%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
    break;
  }

  free(super);
  free(constants);

  return _SUCCESS_;
}

// DEPRECATED
int array_logspline_table_lines(
			     double * x, /* vector of size x_size */
			     int x_size,
			     double * y_array, /* array of size x_size*y_size with elements
						  y_array[index_x*y_size+index_y] */
			     int y_size,
			     double * ddlny_array, /* array of size x_size*y_size */
			     short spline_mode,
			     ErrorMsg errmsg
			     ) {

  int index_x, index_y;
  double *ln_x, *ln_y;
  double *super, *constants;
  double *constants_aux, *sol_aux;

  class_test(x_size < 3, errmsg, "%s(L:%d) there is no spline with less than 3 points", __func__, __LINE__);

  class_alloc(ln_x, x_size*sizeof(double), errmsg);
  class_alloc(ln_y, x_size*sizeof(double), errmsg);

  switch (spline_mode)
  {
  case _SPLINE_NATURAL_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);

    for (index_y=0; index_y < y_size; index_y++) {
      for (index_x=0; index_x < x_size; index_x++) {
        ln_x[index_x] = log( x[index_x] );
        ln_y[index_x] = log( y_array[index_x*y_size + index_y] );
      }

      array_spline_internal_natural(ln_x, 1, ln_y, 1, ddlny_array+index_y, y_size, x_size, super, constants);
    }
    break;

  case _SPLINE_EST_DERIV_:
    class_alloc(super, (x_size-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-1)*sizeof(double), errmsg);
    double dy_first, dy_last;

    for (index_y=0; index_y < y_size; index_y++) {
      for (index_x=0; index_x < x_size; index_x++) {
        ln_x[index_x] = log( x[index_x] );
        ln_y[index_x] = log( y_array[index_x*y_size + index_y] );
      }

      array_spline_internal_hermite(ln_x, 1, ln_y, 1, ddlny_array+index_y, y_size, x_size, super, constants, \
                                    _TRUE_, &dy_first, &dy_last);
    }
    break;

  case _SPLINE_PERIODIC_:
    class_alloc(super, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants, (x_size-2)*sizeof(double), errmsg);
    class_alloc(constants_aux, (x_size-2)*sizeof(double), errmsg);
    class_alloc(sol_aux, (x_size-1)*sizeof(double), errmsg);

    for (index_y=0; index_y < y_size; index_y++) {
      for (index_x=0; index_x < x_size; index_x++) {
        ln_x[index_x] = log( x[index_x] );
        ln_y[index_x] = log( y_array[index_x*y_size + index_y] );
      }

      array_spline_internal_periodic(ln_x, 1, ln_y, 1, ddlny_array+index_y, y_size, x_size-1, super, constants, \
                                      constants_aux, sol_aux);
    }

    free(constants_aux);
    free(sol_aux);
    break;

  default:
    class_stop(errmsg, "%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
    break;
  }

  free(super);
  free(constants);
  free(ln_x);
  free(ln_y);

  return _SUCCESS_;
}

// DEPRECATED
int array_logspline_table_one_column(
		       double * x, /* vector of size x_size */
		       int x_size,
		       int x_stop,
		       double * y_array, /* array of size x_size*y_size with elements
					  y_array[index_y*x_size+index_x] */
		       int y_size,
		       int index_y,
		       double * ddlny_array, /* array of size x_size*y_size */
		       short spline_mode,
		       ErrorMsg errmsg
		       ) {

  int index_x;
  double *ln_x, *ln_y;
  double *super, *constants;
  double *constants_aux, *sol_aux;

  class_test((x_size < 3) || (x_stop < 3), errmsg, "%s(L:%d) there is no spline with less than 3 points", __func__, __LINE__);

  class_alloc(ln_x, x_stop*sizeof(double), errmsg);
  class_alloc(ln_y, x_stop*sizeof(double), errmsg);

  switch (spline_mode)
  {
  case _SPLINE_NATURAL_:
    class_alloc(super, (x_stop-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_stop-1)*sizeof(double), errmsg);

    for (index_x=0; index_x < x_stop; index_x++) {
      ln_x[index_x] = log( x[index_x] );
      ln_y[index_x] = log( y_array[index_y*x_size + index_x] );
    }
    array_spline_internal_natural(ln_x, 1, ln_y, 1, ddlny_array+index_y*x_size, 1, x_stop, super, constants);
    break;

  case _SPLINE_EST_DERIV_:
    class_alloc(super, (x_stop-1)*sizeof(double), errmsg);
    class_alloc(constants, (x_stop-1)*sizeof(double), errmsg);
    double dy_first, dy_last;

    for (index_x=0; index_x < x_stop; index_x++) {
      ln_x[index_x] = log( x[index_x] );
      ln_y[index_x] = log( y_array[index_y*x_size + index_x] );
    }
    array_spline_internal_hermite(ln_x, 1, ln_y, 1, ddlny_array+index_y*x_size, 1, x_stop, super, constants, \
                                  _TRUE_, &dy_first, &dy_last);
    break;

  case _SPLINE_PERIODIC_:
    class_alloc(super, (x_stop-2)*sizeof(double), errmsg);
    class_alloc(constants, (x_stop-2)*sizeof(double), errmsg);
    class_alloc(constants_aux, (x_stop-2)*sizeof(double), errmsg);
    class_alloc(sol_aux, (x_stop-1)*sizeof(double), errmsg);

    for (index_x=0; index_x < x_stop; index_x++) {
      ln_x[index_x] = log( x[index_x] );
      ln_y[index_x] = log( y_array[index_y*x_size + index_x] );
    }
    array_spline_internal_periodic(ln_x, 1, ln_y, 1, ddlny_array+index_y*x_size, 1, x_stop-1, super, constants, \
                                    constants_aux, sol_aux);

    free(constants_aux);
    free(sol_aux);
    break;

  default:
    class_stop(errmsg, "%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
    break;
  }

  free(super);
  free(constants);

  return _SUCCESS_;
}

double dcompsumv(const double * const summands, const int size, const int num_accumulators) {
  /** - TODO */
  return 0.;
}

int array_integrate_internal(
                const double * const x0,
                const int x_size,
                const int x_stride,
                const double * const y0,
                const int y_stride,
                const double * const ddy0,
                const int ddy_stride,
                double * result,
                const double condition_number_threshold,
                double * condition_number) {

  int index_x;
  double h, sy, ty, sM, tM, sum, sum_abs;
  double summands[x_size-1];

  /** - compute the individual summands with compensation */
  for (index_x = 0; index_x < x_size-1; index_x++) {
    h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
    sy = d2sum(y0[(index_x+1)*y_stride], y0[index_x*y_stride], &ty);
    sM = d2sum(ddy0[(index_x+1)*ddy_stride], ddy0[index_x*ddy_stride], &tM);
    summands[index_x] = h/2. * ((sy - h*h/12. * sM) + (ty - h*h/12. * tM));
  }

  /** - sum them up and compute the condition number */
  sum = 0.; sum_abs = 0.;
  for (index_x = 0; index_x < x_size-1; index_x++) {
    sum += summands[index_x];
    #ifndef __FAST_MATH__
    sum_abs += fabs(summands[index_x]);
    #endif
  }
  *condition_number = sum_abs / fabs(sum);

  /** - if condition number is higher than the threshold, recompute using compensated summation */
  if (*condition_number > condition_number_threshold) {
    sum = dcompsum(summands, x_size-1, 1);
  }

  *result = sum;

  return _SUCCESS_;
}

int array_integrate_internal_partial_range(
                const double * const x0,
                const int x_size,
                const int x_stride,
                const double * const y0,
                const int y_stride,
                const double * const ddy0,
                const int ddy_stride,
                const double a,
                const double b,
                const int index_a,
                const int index_b,
                double * result) {

  int index_x;
  double h, b1, b2, sy, sM, dy, dM, sum;

  index_x = index_a;
  h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
  b1 = (a - x0[index_x*x_stride])/h;
  sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
  dy = y0[index_x*y_stride] - y0[(index_x+1)*y_stride];
  sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
  dM = ddy0[index_x*ddy_stride] - ddy0[(index_x+1)*ddy_stride];

  sum = -h/2. * ( ( b1 * sy - b1*(b1 - 1.) * dy )   \
                  -h*h/12. * ( b1*b1*(3. - 2.*b1) * sM + b1*b1*(b1*(b1 - 2.) + 1.) * dM ) );

  if (index_b > index_a) {
    sum += h/2. * ( sy - h*h/12. * sM );

    for (index_x = index_a+1; index_x < index_b; index_x++) {
      h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
      sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
      sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
      sum += h/2. * ( sy - h*h/12. * sM );
    }

    index_x = index_b;
    h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
    b2 = (b - x0[index_x*x_stride])/h;
    sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
    dy = y0[index_x*y_stride] - y0[(index_x+1)*y_stride];
    sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
    dM = ddy0[index_x*ddy_stride] - ddy0[(index_x+1)*ddy_stride];
  }

  sum += h/2. * ( ( b2 * sy - b2*(b2 - 1.) * dy )   \
                  -h*h/12. * ( b2*b2*(3. - 2.*b2) * sM + b2*b2*(b2*(b2 - 2.) + 1.) * dM ) );


  *result = sum;

  return _SUCCESS_;
}

int array_square_integrate_internal(
                const double * const x0,
                const int x_size,
                const int x_stride,
                const double * const y0,
                const int y_stride,
                const double * const ddy0,
                const int ddy_stride,
                double * result) {

  int index_x;
  double h, sy, sM, dM;
  register double sum;

  sum = 0.;
  for (index_x = 0; index_x < x_size-1; index_x++) {
    h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
    sy = y0[(index_x+1)*y_stride] + y0[index_x*y_stride];
    sM = ddy0[(index_x+1)*ddy_stride] + ddy0[index_x*ddy_stride];

    sum += h/3. * (sy*sy - y0[index_x*y_stride]*y0[(index_x+1)*y_stride]                                                                                \
                  -2*h*h/15. * (sM*sy - 0.125*(ddy0[index_x*ddy_stride]*y0[(index_x+1)*y_stride] + ddy0[(index_x+1)*ddy_stride]*y0[index_x*y_stride]))  \
                  +2*h*h*h*h/315. * (sM*sM - 0.0625*ddy0[index_x*ddy_stride]*ddy0[(index_x+1)*ddy_stride]) );
  }
  *result = sum;

  return _SUCCESS_;
}

#define _SPL_SQUARE_EXP_SERIES_THRESHOLD_ 0.001
int array_square_integrate_exponential_internal(
                const double * const x0,
                const int x_size,
                const int x_stride,
                const double * const y0,
                const int y_stride,
                const double * const ddy0,
                const int ddy_stride,
                const double bias,
                const int derivative_order,
                double * result,
                int num_threads) {

  int index_x;
  double h, sy, dy, sM, dM, sx, bias_h;
  register double sum;

  #pragma omp parallel shared(x0, x_size, x_stride, y0, y_stride, ddy0, ddy_stride, bias, derivative_order, result, sum),   \
                       private(index_x, h, sy, dy, sM, dM, sx, bias_h), default(none),                                    \
                       if((num_threads > 1)), num_threads(num_threads)
  {

  sum = 0.;

  switch (derivative_order)
  {
  case 0:
    /** - spline function itself */
    if (fabs(bias) < _SPL_SQUARE_EXP_SERIES_THRESHOLD_) {
      #pragma omp for schedule(static), reduction(+:sum)
      for (index_x = 0; index_x < x_size-1; index_x++) {
        sx = x0[(index_x+1)*x_stride] + x0[index_x*x_stride];
        h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
        sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
        dy = y0[index_x*y_stride] - y0[(index_x+1)*y_stride];
        sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
        dM = ddy0[index_x*ddy_stride] - ddy0[(index_x+1)*ddy_stride];

        sum += h*((2520*dy*dy - 84*dM*dy*h*h + h*h*h*h*(dM*dM + 63*sM*sM) - 1260*h*h*sM*sy + 7560*sy*sy)/30240.    \
                  + bias*(2520*dy*dy*sx + h*h*h*h*(-6*dM*sM*h + dM*dM*sx + 63*sM*sM*sx)    \
                          + 84*dy*h*(3*h*h*sM - dM*h*sx - 60*sy) + 84*h*h*(dM*h - 15*sM*sx)*sy + 7560*sx*sy*sy)/30240.      \
                  + bias*bias*(1512*dy*dy*(3*h*h + 5*sx*sx) + h*h*h*h*(-36*dM*sM*h*sx + dM*dM*(h*h + 3*sx*sx) + 27*sM*sM*(h*h + 7*sx*sx))     \
                               - 252*h*h*(3*h*h*sM - 2*dM*h*sx + 15*sM*sx*sx)*sy + 7560*(h*h + 3*sx*sx)*sy*sy - 36*dy*h*(3*dM*h*h*h - 42*h*h*sM*sx + 7*dM*h*sx*sx + 840*sx*sy))/181440.);
      }
    }
    else {
      #pragma omp for schedule(static), reduction(+:sum)
      for (index_x = 0; index_x < x_size-1; index_x++) {
        sx = x0[(index_x+1)*x_stride] + x0[index_x*x_stride];
        h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
        sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
        dy = y0[index_x*y_stride] - y0[(index_x+1)*y_stride];
        sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
        dM = ddy0[index_x*ddy_stride] - ddy0[(index_x+1)*ddy_stride];
        bias_h = bias * h;

        sum += exp(bias * sx) / (288.*pow(bias, 5)*bias_h*bias_h)   \
              * ( -3*bias_h*(3*dM*dM*(10. + bias_h*bias_h) + 2*dM*sM*bias_h*(15. + bias_h*bias_h) + 9*bias_h*bias_h*sM*sM   \
                              + 48*bias*bias*bias*bias*dy*(dy + bias_h*sy) + 12*bias*bias*bias_h*sM*(3*dy + bias_h*sy)      \
                              + 4*bias*bias*dM*(dy*(12. + bias_h*bias_h) + 3*bias_h*sy)) * cosh(bias_h)
                + (dM*dM*(90. + 39*bias_h*bias_h + bias_h*bias_h*bias_h*bias_h)       \
                    + 6*dM*(3*bias_h*(5. + 2*bias_h*bias_h)*sM + 2*bias*bias*(dy*(12. + 5*bias_h*bias_h) + bias_h*(3. + bias_h*bias_h)*sy))   \
                    + 9*(bias_h*bias_h*(3. + bias_h*bias_h)*sM*sM + 4*bias*bias*bias_h*sM*(dy*(3. + bias_h*bias_h) + bias_h*sy) + 8*bias*bias*bias*bias*(dy*dy*(2. + bias_h*bias_h) + 2*dy*sy*bias_h + bias_h*bias_h*sy*sy))) * sinh(bias_h));
      }
    }
    break;

  case 1:
    /** - first derivative */
    if (fabs(bias) < _SPL_SQUARE_EXP_SERIES_THRESHOLD_) {
      #pragma omp for schedule(static), reduction(+:sum)
      for (index_x = 0; index_x < x_size-1; index_x++) {
        sx = x0[(index_x+1)*x_stride] + x0[index_x*x_stride];
        h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
        // sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
        dy = y0[index_x*y_stride] - y0[(index_x+1)*y_stride];
        sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
        dM = ddy0[index_x*ddy_stride] - ddy0[(index_x+1)*ddy_stride];

        sum += dy*dy/h + (h*h*h*h*h*(dM*dM + 15*sM*sM))/720.      \
              + bias*(-dy*h*h*h*sM/6. + dy*dy*sx/h + h*h*h*h*h*(-2*dM*h*sM + sM*dM*sx + 15*sM*sM*sx))/720.    \
              + bias*bias*(dy*h*h*h*(dM*h - 15*sM*sx))/90. + (dy*dy*(h*h + 3*sx*sx)) / (6.*h)                 \
                          + (h*h*h*h*h*(-168*dM*h*sM*sx + 63*sM*sM*(3*h*h + 5*sx*sx) + dM*dM*(11*h*h + 21*sx*sx)))/30240.;
      }
    }
    else {
      #pragma omp for schedule(static), reduction(+:sum)
      for (index_x = 0; index_x < x_size-1; index_x++) {
        sx = x0[(index_x+1)*x_stride] + x0[index_x*x_stride];
        h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
        // sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
        dy = y0[index_x*y_stride] - y0[(index_x+1)*y_stride];
        sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
        dM = ddy0[index_x*ddy_stride] - ddy0[(index_x+1)*ddy_stride];
        bias_h = bias * h;

        sum += exp(bias * sx) / (144.*pow(bias, 5)*bias_h*bias_h)   \
              * (-6*bias_h*bias_h*(12*bias*bias*bias*dy*(dM + bias_h*sM) + bias_h*(dM*dM*(9. + bias_h*bias_h)     \
                                    + dM*sM*bias_h*(9. + bias_h*bias_h) + 3*bias_h*bias_h*sM*sM)) * cosh(bias_h)
                + (144*bias*bias*bias*bias*bias*bias*dy*dy + 24*bias*bias*bias*dy*bias_h*(dM*(3. + bias_h*bias_h)     \
                    + 3*bias_h*sM) + bias_h*bias_h*(dM*dM*(54. + 24*bias_h*bias_h + bias_h*bias_h*bias_h*bias_h)      \
                    + 6*dM*sM*bias_h*(9. + 4*bias_h*bias_h) + 9*bias_h*bias_h*(2. + bias_h*bias_h)*sM*sM)) * sinh(bias_h));
      }
    }
    break;

  case 2:
    /** - second derivative */
    if (fabs(bias) < _SPL_SQUARE_EXP_SERIES_THRESHOLD_) {
      #pragma omp for schedule(static), reduction(+:sum)
      for (index_x = 0; index_x < x_size-1; index_x++) {
        sx = x0[(index_x+1)*x_stride] + x0[index_x*x_stride];
        h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
        // sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
        // dy = y0[index_x*y_stride] - y0[(index_x+1)*y_stride];
        sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
        dM = ddy0[index_x*ddy_stride] - ddy0[(index_x+1)*ddy_stride];

        sum += h*((dM*dM + 3*sM*sM)/12.       \
                + bias*(-2*dM*h*sM + dM*dM*sx + 3*sM*sM*sx)/12.       \
                + bias*bias*(-20*dM*h*sM*sx + 5*sM*sM*(h*h + 3*sx*sx) + dM*dM*(3*h*h + 5*sx*sx))/120.);
      }
    }
    else {
      #pragma omp for schedule(static), reduction(+:sum)
      for (index_x = 0; index_x < x_size-1; index_x++) {
        sx = x0[(index_x+1)*x_stride] + x0[index_x*x_stride];
        h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
        // sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
        // dy = y0[index_x*y_stride] - y0[(index_x+1)*y_stride];
        sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
        dM = ddy0[index_x*ddy_stride] - ddy0[(index_x+1)*ddy_stride];
        bias_h = bias * h;

        sum += exp(bias * sx) / (4.*bias*bias_h*bias_h)     \
              * (-2*dM*bias_h*(dM + bias_h*sM) * cosh(bias_h)     \
                  + (dM*dM*(2. + bias_h*bias_h) + 2*dM*sM*bias_h + bias_h*bias_h*sM*sM) * sinh(bias_h));
      }
    }
    break;

  case 3:
    /** - third derivative */
    if (fabs(bias) < _SPL_SQUARE_EXP_SERIES_THRESHOLD_) {
      #pragma omp for schedule(static), reduction(+:sum)
      for (index_x = 0; index_x < x_size-1; index_x++) {
        sx = x0[(index_x+1)*x_stride] + x0[index_x*x_stride];
        h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
        // sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
        // dy = y0[index_x*y_stride] - y0[(index_x+1)*y_stride];
        // sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
        dM = ddy0[index_x*ddy_stride] - ddy0[(index_x+1)*ddy_stride];

        sum += 1./h*(dM*dM + bias*dM*dM*sx + bias*bias*dM*dM*(h*h + 3*sx*sx)/6.);
      }
    }
    else {
      #pragma omp for schedule(static), reduction(+:sum)
      for (index_x = 0; index_x < x_size-1; index_x++) {
        sx = x0[(index_x+1)*x_stride] + x0[index_x*x_stride];
        h = x0[(index_x+1)*x_stride] - x0[index_x*x_stride];
        // sy = y0[index_x*y_stride] + y0[(index_x+1)*y_stride];
        // dy = y0[index_x*y_stride] - y0[(index_x+1)*y_stride];
        // sM = ddy0[index_x*ddy_stride] + ddy0[(index_x+1)*ddy_stride];
        dM = ddy0[index_x*ddy_stride] - ddy0[(index_x+1)*ddy_stride];
        bias_h = bias * h;

        sum += bias*dM*dM * exp(bias * sx) * sinh(bias_h) / (bias_h*bias_h);
      }
    }
    break;

  default:
    /** - higher derivative vanish */
    break;
  }

  } /** - end of parallel region */

  *result = sum;

  return _SUCCESS_;
}

/**
 * @brief Computes the spline integral.
 *        dI = dx S(x)
 * @param array       Input: contains x, y and y'' values retrieved from splining
 * @param n_columns   Input: number of columns in array, indexed by index_x/y/ddy
 * @param n_lines     Input: number of used control points
 * @param index_x     Input: index for x-values (->class_define_index)
 * @param index_y     Input: index for y-values (->class_define_index)
 * @param index_ddy   Input: index for y''-values (->class_define_index)
 * @param result      Output: integration result I
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_spline(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y,
		   int index_ddy,
		   double * result,
		   ErrorMsg errmsg) {

  double cond_num;

  return array_integrate_internal(array + index_x,
                                  n_lines,
                                  n_columns,
                                  array + index_y,
                                  n_columns,
                                  array + index_ddy,
                                  n_columns,
                                  result,
                                  _SPL_INTEGRATION_DEF_COND_THRESHOLD_,
                                  &cond_num);
}

/**
 * @brief Computes the spline integral.
 *        dI = dx S(x)
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
int array_integrate_all_spline_table_lines(
		      double * x,
			    int x_size,
			    double * y_array,
			    int y_size,
			    double * ddy_array,
			    double * result,
          ErrorMsg errmsg) {

  int index_y;
  double cond_num[y_size];

  for (index_y = 0; index_y < y_size; index_y++) {
    array_integrate_internal(x,
                             x_size,
                             1,
                             y_array + index_y,
                             y_size,
                             ddy_array + index_y,
                             y_size,
                             result + index_y,
                             _SPL_INTEGRATION_DEF_COND_THRESHOLD_,
                             cond_num + index_y);
  }

  return _SUCCESS_;
}

/**
 * @brief Computes the spline integral in a partial range.
 *        dI = dx S(x) on [a, b]
 * @param x           Input: contains x-values of the integration range
 * @param x_size      Input: size of x-array to be used
 * @param y_array     Input: contains y-values of the integrand with elements
 *                            y_array[index_x*y_size + index_y]
 * @param y_size      Input: number of columns in y-array
 * @param ddy_array   Input: contains y''-values of the integrand obtained from splining with elements
 *                            ddy_array[index_x*y_size + index_y]
 * @param a           Input: lower limit of integration region
 * @param b           Input: upper limit of integration region
 * @param result      Output: integration result I with elements[index_y]
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_spline_partial_range_table_lines(
		      double * x,
			    int x_size,
			    double * y_array,
			    int y_size,
			    double * ddy_array,
          double a,
          double b,
			    double * result,
          ErrorMsg errmsg) {

  int index_y, index_a, index_b, sup;
  ErrorMsg err_search;

  class_call(array_search_bisect_internal(x,
                                          1,
                                          x_size,
                                          a,
                                          &index_a,
                                          &sup,
                                          err_search),
              err_search, errmsg);

  class_call(array_search_bisect_internal(x,
                                          1,
                                          x_size,
                                          b,
                                          &index_b,
                                          &sup,
                                          err_search),
              err_search, errmsg);

  for (index_y = 0; index_y < y_size; index_y++) {
    array_integrate_internal_partial_range(x,
                                           x_size,
                                           1,
                                           y_array + index_y,
                                           y_size,
                                           ddy_array + index_y,
                                           y_size,
                                           a,
                                           b,
                                           index_a,
                                           index_b,
                                           result + index_y);
  }

  return _SUCCESS_;
}

/**
 * @brief Computes the spline integral in a partial range.
 *        dI = dx S(x) on [a, b]
 * @param x           Input: contains x-values of the integration range
 * @param x_size      Input: size of x-array to be used
 * @param y_array     Input: contains y-values of the integrand with elements
 *                            y_array[index_x*y_size + index_y]
 * @param y_size      Input: number of columns in y-array
 * @param ddy_array   Input: contains y''-values of the integrand obtained from splining with elements
 *                            ddy_array[index_x*y_size + index_y]
 * @param a           Input: lower limit of integration region
 * @param b           Input: upper limit of integration region
 * @param last_index  In/Output: close to infimum index of a in x, overwritten with infimum index of b in x
 * @param result      Output: integration result I with elements[index_y]
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_spline_partial_range_table_lines_closeby(
		      double * x,
			    int x_size,
			    double * y_array,
			    int y_size,
			    double * ddy_array,
          double a,
          double b,
          int * last_index,
			    double * result,
          ErrorMsg errmsg) {

  int index_y, index_a, index_b, sup;
  ErrorMsg err_search;

  index_a = *last_index;
  class_call(array_search_closeby_internal(x,
                                          1,
                                          x_size,
                                          a,
                                          &index_a,
                                          &sup,
                                          err_search),
              err_search, errmsg);

  index_b = sup;
  class_call(array_search_closeby_internal(x,
                                          1,
                                          x_size,
                                          b,
                                          &index_b,
                                          &sup,
                                          err_search),
              err_search, errmsg);

  for (index_y = 0; index_y < y_size; index_y++) {
    array_integrate_internal_partial_range(x,
                                           x_size,
                                           1,
                                           y_array + index_y,
                                           y_size,
                                           ddy_array + index_y,
                                           y_size,
                                           a,
                                           b,
                                           index_a,
                                           index_b,
                                           result + index_y);
  }

  *last_index = index_b;

  return _SUCCESS_;
}

/**
 * @brief Computes the spline integral.
 *        dI = dx S(x)
 * @param x           Input: contains x-values of the integration range
 * @param x_size      Input: size of x-array to be used
 * @param y_array     Input: contains y-values of the integrand with elements
 *                            y_array[index_x*y_size + index_y]
 * @param y_size      Input: number of columns in y-array
 * @param ddy_array   Input: contains y''-values of the integrand obtained from splining with elements
 *                            ddy_array[index_x*y_size + index_y]
 * @param result      Output: integration result I with elements[index_y]
 * @param condition_num_threshold   Input: threshold for condition number above which the sum is computed with compensation
 * @param condition_num   Output: condition number of the underlying sum with elements[index_y]
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_spline_table_lines_compensated(
		      double * x,
			    int x_size,
			    double * y_array,
			    int y_size,
			    double * ddy_array,
			    double * result,
          const double condition_num_threshold,
          double * condition_num,
          ErrorMsg errmsg) {

  int index_y;

  for (index_y = 0; index_y < y_size; index_y++) {
    array_integrate_internal(x,
                             x_size,
                             1,
                             y_array + index_y,
                             y_size,
                             ddy_array + index_y,
                             y_size,
                             result + index_y,
                             condition_num_threshold,
                             condition_num + index_y);
  }

  return _SUCCESS_;
}


/**
 * @brief Computes the spline integral with a specified gaussian window exactly.
 *        dI(m, s) = dx S(x) * exp(-(x-m)^2/(2 s^2)) / (sqrt(2 pi) s)
 * @param array       Input: contains x, y and y'' values retrieved from splining
 * @param n_columns   Input: number of columns in array, indexed by index_x/y/ddy
 * @param n_lines     Input: number of used control points
 * @param index_x     Input: index for x-values (->class_define_index)
 * @param index_y     Input: index for y-values (->class_define_index)
 * @param index_ddy   Input: index for y''-values (->class_define_index)
 * @param mean        Input: mean of the gaussian distribution
 * @param stddev      Input: standard deviation
 * @param result      Output: integration result I(m, s)
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_spline_gaussian_window(
          double * array,
          int n_columns,
          int n_lines,
          int index_x,   /** from 0 to (n_columns-1) */
          int index_y,
          int index_ddy,
          double mean,
          double stddev,
          double * result,
          ErrorMsg errmsg) {

  int i;
  double h;                          /**< distance between control points */
  double a, b, g, d, Ms, dM, ys, dy; /**< alpha, beta, gamma, delta, sum/difference of y and y'' */

  *result = 0;

  for (i=0; i < n_lines-1; i++) {

    h = array[(i+1)*n_columns+index_x] - array[i*n_columns+index_x];
    a = (array[i*n_columns+index_x] - mean) / (_SQRT2_ * stddev);
    b = (array[(i+1)*n_columns+index_x] - mean) / (_SQRT2_ * stddev);
    g = (a + b)/2;
    d = h / (2.*_SQRT2_* stddev);
    Ms = array[i*n_columns+index_ddy] + array[(i+1)*n_columns+index_ddy];
    dM = array[i*n_columns+index_ddy] - array[(i+1)*n_columns+index_ddy];
    ys = array[i*n_columns+index_y] + array[(i+1)*n_columns+index_y];
    dy = array[i*n_columns+index_y] - array[(i+1)*n_columns+index_y];

    *result += 1./(2.*d) * (
      (erf(b) - erf(a))/2. * ( dy*g + ys*d + stddev*stddev/3. * (3.*Ms*d*(g*g - d*d + 0.5) + dM*g*(g*g - d*d + 1.5)) ) \
     -1./(2.*_SQRT_PI_) * ( exp(-(g-d)*(g-d))*( dy + stddev*stddev/3. * (3.*Ms*d*(g+d) + dM*(g*(g+d) + 1.)) ) \
                           -exp(-(g+d)*(g+d))*( dy + stddev*stddev/3. * (3.*Ms*d*(g-d) + dM*(g*(g-d) + 1.)) ) ) \
                            );

  }

  return _SUCCESS_;
}


/**
 * @brief Computes the spline integral with separate x-array.
 *        dI = dx S(x)
 * @param x_array     Input: contains x values
 * @param n_lines     Input: number of used control points
 * @param array       Input: contains y and y'' values retrieved from splining
 * @param n_columns   Input: number of columns in array, indexed by index_y/ddy
 * @param index_y     Input: index for y-values (->class_define_index)
 * @param index_ddy   Input: index for y''-values (->class_define_index)
 * @param result      Output: integration result I
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_spline_table_line_to_line(
                  double * x_array,
                  int n_lines,
                  double * array,
                  int n_columns,
                  int index_y,
                  int index_ddy,
                  double * result,
                  ErrorMsg errmsg){

  int i;
  double h;

  *result = 0;

  for (i=0; i < n_lines-1; i++) {

    h = (x_array[i+1]-x_array[i]);

    *result +=
      (array[i*n_columns+index_y]+array[(i+1)*n_columns+index_y])*h/2. -
      (array[i*n_columns+index_ddy]+array[(i+1)*n_columns+index_ddy])*h*h*h/24.;

  }

  return _SUCCESS_;
}

/**
 * @brief Computes the integral using trapezoidal and spline integration.
 *        dI = dx y or dI = dx S(x)
 * @param array       Input: contains x, y and y'' values retrieved from splining
 * @param n_columns   Input: number of columns in array, indexed by index_x/y/ddy
 * @param n_lines     Input: number of used control points
 * @param index_start_spline  Input: use spline integration from this index onwards
 * @param index_x     Input: index for x-values (->class_define_index)
 * @param index_y     Input: index for y-values (->class_define_index)
 * @param index_ddy   Input: index for y''-values (->class_define_index)
 * @param result      Output: integration result I
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_all_trapzd_or_spline(
        double * array,
        int n_columns,
        int n_lines,
        int index_start_spline,
        int index_x,   /** from 0 to (n_columns-1) */
        int index_y,
        int index_ddy,
        double * result,
        ErrorMsg errmsg) {

  int i;
  double h;

  if ((index_start_spline<0) || (index_start_spline>=n_lines)) {
    class_sprintf(errmsg,"%s(L:%d) index_start_spline outside of range",__func__,__LINE__);
    return _FAILURE_;
  }

  *result = 0;

  /* trapezoidal integration till given index */

  for (i=0; i < index_start_spline; i++) {

    h = (array[(i+1)*n_columns+index_x]-array[i*n_columns+index_x]);

    *result +=
      (array[i*n_columns+index_y]+array[(i+1)*n_columns+index_y])*h/2.;

  }

  /* then, spline integration */

  for (i=index_start_spline; i < n_lines-1; i++) {

    h = (array[(i+1)*n_columns+index_x]-array[i*n_columns+index_x]);

    *result +=
      (array[i*n_columns+index_y]+array[(i+1)*n_columns+index_y])*h/2. -
      (array[i*n_columns+index_ddy]+array[(i+1)*n_columns+index_ddy])*h*h*h/24.;

  }

  return _SUCCESS_;
}

 /**
 * @brief Computes the integral using trapezoidal rule.
 *        dI = dx y
 * @param array       Input: contains x, y and I
 * @param n_columns   Input: number of columns in array, indexed by index_x/y/int_y_dx
 * @param n_lines     Input: number of used control points
 * @param index_x     Input: index for x-values (->class_define_index)
 * @param index_y     Input: index for y-values (->class_define_index)
 * @param index_int_y_dx  Input: index for integration result
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y,
		   int index_int_y_dx,
		   ErrorMsg errmsg) {

  int i;
  double sum;

  if ((index_int_y_dx == index_x) || (index_int_y_dx == index_y)) {
    class_sprintf(errmsg,"%s(L:%d) : Output column %d must differ from input columns %d and %d",__func__,__LINE__,index_int_y_dx,index_x,index_y);
    return _FAILURE_;
  }

  sum=0.;
  *(array+0*n_columns+index_int_y_dx)=sum;

  for (i=1; i<n_lines; i++) {

    sum += 0.5 * (*(array+i*n_columns+index_y) + *(array+(i-1)*n_columns+index_y))
               * (*(array+i*n_columns+index_x) - *(array+(i-1)*n_columns+index_x));

    *(array+i*n_columns+index_int_y_dx)=sum;
  }


  return _SUCCESS_;
}

 /**
 * @brief Computes the integral of a ratio using trapezoidal rule.
 *        dI = dx y1/y2
 *        Called by thermodynamics_init()
 * @param array       Input: contains x, y1, y2 and I
 * @param n_columns   Input: number of columns in array, indexed by index_x/y1/y2/int_y1_over_y2_dx
 * @param n_lines     Input: number of used control points
 * @param index_x     Input: index for x-values (->class_define_index)
 * @param index_y1    Input: index for y1-values (->class_define_index)
 * @param index_y2    Input: index for y2-values (->class_define_index)
 * @param index_int_y1_over_y2_dx  Input: index for integration result
 * @param errmsg
 *
 * @return the error status
 */
int array_integrate_ratio(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y1,
		   int index_y2,
		   int index_int_y1_over_y2_dx,
		   ErrorMsg errmsg) {

  int i;
  double sum;

  if ((index_int_y1_over_y2_dx == index_x) || (index_int_y1_over_y2_dx == index_y1) || (index_int_y1_over_y2_dx == index_y2)) {
    class_sprintf(errmsg,"%s(L:%d) : Output column %d must differ from input columns %d, %d and %d",__func__,__LINE__,index_int_y1_over_y2_dx,index_x,index_y1,index_y2);
    return _FAILURE_;
  }

  sum=0.;

  *(array+0*n_columns+index_int_y1_over_y2_dx)=sum;

  for (i=1; i<n_lines; i++) {

    sum += 0.5 * (*(array+i*n_columns+index_y1) / *(array+i*n_columns+index_y2)
		  + *(array+(i-1)*n_columns+index_y1) / *(array+(i-1)*n_columns+index_y2))
      * (*(array+i*n_columns+index_x) - *(array+(i-1)*n_columns+index_x));

    *(array+i*n_columns+index_int_y1_over_y2_dx)=sum;
  }


  return _SUCCESS_;
}

int array_interpolate_internal(
                               const double * const x0,
                               const int x_size,
                               const int x_stride,
                               const double * const y0,
                               const int y_stride,
                               const double * const ddy0,
                               const int ddy_stride,
                               const short derivative_order,
                               const double x,
                               const int inf,
                               const int sup,
                               double * result) {

  int i;
  double h,a,b;

  /** - negative deriavtive order would require integration (non-local) */
  // if (derivative_order < 0) return _FAILURE_;  handled by caller

  h = x0[sup*x_stride] - x0[inf*x_stride];
  b = (x - x0[inf*x_stride]) / h;
  a = 1 - b;

  switch (derivative_order)
  {
  case 0:
    /* spline value */
    *result = y0[inf*y_stride]*a + y0[sup*y_stride]*b         \
              + h*h/6. * (ddy0[inf*ddy_stride]*a*(a*a - 1.) + ddy0[sup*ddy_stride]*b*(b*b - 1.));
    break;
  case 1:
    /* first derivative */
    *result = (y0[sup*y_stride] - y0[inf*y_stride]) / h       \
              + h/6. * (ddy0[sup*ddy_stride]*(3*b*b - 1.) - ddy0[inf*ddy_stride]*(3*a*a - 1.));
    break;
  case 2:
    /* second derivative */
    *result = ddy0[inf*ddy_stride]*a + ddy0[sup*ddy_stride]*b;
    break;
  case 3:
    /* third derivative (not continuous) */
    *result = (ddy0[sup*ddy_stride] - ddy0[inf*ddy_stride]) / h;
    break;

  default:
    *result = 0.;
    break;
  }

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   double x,
		   int * last_index,
		   double * result,
		   int result_size, /** from 1 to n_columns */
		   ErrorMsg errmsg) {

  int inf,sup,mid,i;
  double weight;

  inf=0;
  sup=n_lines-1;

  if (*(array+inf*n_columns+index_x) < *(array+sup*n_columns+index_x)){

    if (x < *(array+inf*n_columns+index_x)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,*(array+inf*n_columns+index_x));
      return _FAILURE_;
    }

    if (x > *(array+sup*n_columns+index_x)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,*(array+sup*n_columns+index_x));
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < *(array+mid*n_columns+index_x)) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < *(array+sup*n_columns+index_x)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,*(array+sup*n_columns+index_x));
      return _FAILURE_;
    }

    if (x > *(array+inf*n_columns+index_x)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,*(array+inf*n_columns+index_x));
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > *(array+mid*n_columns+index_x)) {sup=mid;}
      else {inf=mid;}

    }

  }

  *last_index = inf;

  weight=(x-*(array+inf*n_columns+index_x))/(*(array+sup*n_columns+index_x)-*(array+inf*n_columns+index_x));

  for (i=0; i<result_size; i++)
    *(result+i) = *(array+inf*n_columns+i) * (1.-weight)
      + weight * *(array+sup*n_columns+i);

  *(result+index_x) = x;

  return _SUCCESS_;
}


 /**
  * interpolate to get y(x), when x_i,y_i and ddy_i are all columns of the same array
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_spline_transposed(double * array,
                                        int x_size,
                                        int y_size,
                                        int index_x,
                                        int index_y,
                                        int index_ddy,
                                        double x,
                                        int * last_index,
                                        double * result,
                                        ErrorMsg errmsg) {

  int inf,sup,mid;
  double h,a,b;

  inf=0;
  sup=x_size-1;

  if (array[inf*y_size+index_x] < array[sup*y_size+index_x]){

    if (x < array[inf*y_size+index_x]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,array[inf*y_size+index_x]);
      return _FAILURE_;
    }

    if (x > array[sup*y_size+index_x]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,array[sup*y_size+index_x]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < array[mid*y_size+index_x]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < array[sup*y_size+index_x]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,array[sup*y_size+index_x]);
      return _FAILURE_;
    }

    if (x > array[inf*y_size+index_x]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,array[inf*y_size+index_x]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > array[mid*y_size+index_x]) {sup=mid;}
      else {inf=mid;}

    }

  }

  *last_index = inf;

  h = array[sup*y_size+index_x]-array[inf*y_size+index_x];
  b = (x - array[inf*y_size+index_x])/h;
  a = 1.0 - b;

  *result= (a*array[inf*y_size+index_y]+b*array[sup*y_size+index_y]
            + ((a*a*a-a)* array[inf*y_size+index_ddy] + (b*b*b-b)* array[sup*y_size+index_ddy])*h*h/6.);

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are in different arrays
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_spline(
                             double * __restrict__ x_array,
                             int n_lines,
                             double * __restrict__ array,
                             double * __restrict__ array_splined,
                             int n_columns,
                             double x,
                             int * __restrict__ last_index,
                             double * __restrict__ result,
                             int result_size, /** from 1 to n_columns */
                             ErrorMsg errmsg) {

  int inf,sup,mid,i;
  double h,a,b;

  inf=0;
  sup=n_lines-1;

  if (x_array[inf] < x_array[sup]){

    if (x < x_array[inf]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    if (x > x_array[sup]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < x_array[sup]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    if (x > x_array[inf]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  for (i=0; i<result_size; i++)
    *(result+i) =
      a * *(array+inf*n_columns+i) +
      b * *(array+sup*n_columns+i) +
      ((a*a*a-a)* *(array_splined+inf*n_columns+i) +
       (b*b*b-b)* *(array_splined+sup*n_columns+i))*h*h/6.;

  return _SUCCESS_;
}

int array_interpolate_spline_derivative(
                const double * const __restrict__ x_array,
                const int n_lines,
                const double * const __restrict__ array,
                const double * const __restrict__ array_splined,
                const int n_columns,
                const double x,
                const short derivative_order,
                int * const __restrict__ last_index,
                double * __restrict__ result,
                const int result_size, /** from 1 to n_columns */
                ErrorMsg errmsg) {

  int inf, sup, i;

  class_call(array_search_bisect_internal(x_array,
                                          1,
                                          n_lines,
                                          x,
                                          &inf,
                                          &sup,
                                          errmsg),
              errmsg, errmsg);

  *last_index = inf;

  for (i = 0; i < result_size; i++) {
    array_interpolate_internal(x_array,
                              n_lines,
                              1,
                              array,
                              n_columns,
                              array_splined,
                              n_columns,
                              derivative_order,
                              x,
                              inf,
                              sup,
                              result + i);
  }

  return _SUCCESS_;
}

int array_interpolate_spline_derivative_closeby(
                const double * const __restrict__ x_array,
                const int n_lines,
                const double * const __restrict__ array,
                const double * const __restrict__ array_splined,
                const int n_columns,
                const double x,
                const short derivative_order,
                int * const __restrict__ last_index,
                double * __restrict__ result,
                const int result_size, /** from 1 to n_columns */
                ErrorMsg errmsg) {

  int inf, sup, i;

  inf = *last_index;

  class_call(array_search_closeby_internal(x_array,
                                           1,
                                           n_lines,
                                           x,
                                           &inf,
                                           &sup,
                                           errmsg),
              errmsg, errmsg);

  *last_index = inf;

  for (i = 0; i < result_size; i++) {
    array_interpolate_internal(x_array,
                              n_lines,
                              1,
                              array,
                              n_columns,
                              array_splined,
                              n_columns,
                              derivative_order,
                              x,
                              inf,
                              sup,
                              result + i);
  }

  return _SUCCESS_;
}

 /**
  * Get the y[i] for which y[i]>c
  *
  * Called by fourier_HMcode()
  */
int array_search_bisect(
                        int n_lines,
                        double * __restrict__ array,
                        double c,
                        int * __restrict__ last_index,
                        ErrorMsg errmsg) {

  int inf,sup,mid;

  inf=0;
  sup=n_lines-1;

  if (array[inf] < array[sup]){

    if (c < array[inf]) {
      class_sprintf(errmsg,"%s(L:%d) : c=%e < y_min=%e",__func__,__LINE__,c,array[inf]);
      return _FAILURE_;
    }

    if (c > array[sup]) {
      class_sprintf(errmsg,"%s(L:%d) : c=%e > y_max=%e",__func__,__LINE__,c,array[sup]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (c < array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (c < array[sup]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,c,array[sup]);
      return _FAILURE_;
    }

    if (c > array[inf]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,c,array[inf]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (c > array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  *last_index = inf;

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are in different arrays
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_linear(
			     double * x_array,
			     int n_lines,
			     double * array,
			     int n_columns,
			     double x,
			     int * last_index,
			     double * result,
			     int result_size, /** from 1 to n_columns */
			     ErrorMsg errmsg) {

  int inf,sup,mid,i;
  double h,a,b;

  inf=0;
  sup=n_lines-1;

  if (x_array[inf] < x_array[sup]){

    if (x < x_array[inf]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    if (x > x_array[sup]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < x_array[sup]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    if (x > x_array[inf]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  for (i=0; i<result_size; i++)
    *(result+i) =
      a * *(array+inf*n_columns+i) +
      b * *(array+sup*n_columns+i);

  return _SUCCESS_;
}


/**
 * interpolate to get y_i(x), when x and y_i are in different arrays
 *
 * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
 */
int array_interpolate_logspline(
							 double * x_array,
							 int n_lines,
							 double * array,
							 double * array_logsplined,
							 int n_columns,
							 double x,
							 int * last_index,
							 double * result,
							 int result_size, /** from 1 to n_columns */
							 ErrorMsg errmsg) {

	int inf,sup,mid,i;
	double h,a,b;

	inf=0;
	sup=n_lines-1;

	if (x_array[inf] < x_array[sup]){

		if (x < x_array[inf]) {
			class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[inf]);
			return _FAILURE_;
		}

		if (x > x_array[sup]) {
			class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[sup]);
			return _FAILURE_;
		}

		while (sup-inf > 1) {

			mid=(int)(0.5*(inf+sup));
			if (x < x_array[mid]) {sup=mid;}
			else {inf=mid;}

		}

	}

	else {

		if (x < x_array[sup]) {
			class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[sup]);
			return _FAILURE_;
		}

		if (x > x_array[inf]) {
			class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[inf]);
			return _FAILURE_;
		}

		while (sup-inf > 1) {

			mid=(int)(0.5*(inf+sup));
			if (x > x_array[mid]) {sup=mid;}
			else {inf=mid;}

		}

	}

	*last_index = inf;

	h = log(x_array[sup]) - log(x_array[inf]);
	b = (log(x)-log(x_array[inf]))/h;
	a = 1-b;

	for (i=0; i<result_size; i++)
		*(result+i) = exp(
		a * log(array[inf*n_columns+i]) +
		b * log(array[sup*n_columns+i]) +
		((a*a*a-a)* array_logsplined[inf*n_columns+i] +
		 (b*b*b-b)* array_logsplined[sup*n_columns+i])*h*h/6.);

	return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are in different arrays
  *
  *
  */
int array_interpolate_spline_one_column(
					double * x_array,
					int x_size,
					double * y_array, /* array of size x_size*y_size with elements
							   y_array[index_y*x_size+index_x] */
					int y_size,
					int index_y,
					double * ddy_array, /* array of size x_size*y_size */
					double x,   /* input */
					double * y, /* output */
					ErrorMsg errmsg
					) {


  int inf,sup,mid;
  double h,a,b;

  inf=0;
  sup=x_size-1;

  if (x_array[inf] < x_array[sup]){

    if (x < x_array[inf]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    if (x > x_array[sup]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < x_array[sup]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    if (x > x_array[inf]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  *y =
    a * y_array[index_y * x_size + inf] +
    b * y_array[index_y * x_size + sup] +
    ((a*a*a-a)* ddy_array[index_y * x_size + inf] +
     (b*b*b-b)* ddy_array[index_y * x_size + sup])*h*h/6.;

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are in different arrays
  *
  *
  */
int array_interpolate_extrapolate_spline_one_column(
						    double * x_array,
						    int x_size,
						    double * y_array, /* array of size x_size*y_size with elements
									 y_array[index_y*x_size+index_x] */
						    int y_size,
						    int index_y,
						    double * ddy_array, /* array of size x_size*y_size */
						    double x,   /* input */
						    double * y, /* output */
						    ErrorMsg errmsg
						    ) {


  int inf,sup,mid;
  double h,a,b;

  if (x > x_array[x_size-2] || x < x_array[0]) {

    /*interpolate/extrapolate linearly y as a function of x*/

    h = x_array[x_size-1] - x_array[x_size-2];
    b = (x-x_array[x_size-2])/h;
    a = 1-b;

    *y = a * y_array[index_y * x_size + (x_size-2)] +
	     b * y_array[index_y * x_size + (x_size-1)];


  }

  else {

    /*interpolate y as a function of x with a spline*/

    inf=0;
    sup=x_size-1;

    if (x_array[inf] < x_array[sup]){

      if (x < x_array[inf]) {
	class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[inf]);
	return _FAILURE_;
      }

      if (x > x_array[sup]) {
	class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[sup]);
	return _FAILURE_;
      }

      while (sup-inf > 1) {

	mid=(int)(0.5*(inf+sup));
	if (x < x_array[mid]) {sup=mid;}
	else {inf=mid;}

      }

    }

    else {

      if (x < x_array[sup]) {
	class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[sup]);
	return _FAILURE_;
      }

      if (x > x_array[inf]) {
	class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[inf]);
	return _FAILURE_;
      }

      while (sup-inf > 1) {

	mid=(int)(0.5*(inf+sup));
	if (x > x_array[mid]) {sup=mid;}
	else {inf=mid;}

      }

    }

    h = x_array[sup] - x_array[inf];
    b = (x-x_array[inf])/h;
    a = 1-b;

    *y =
      a * y_array[index_y * x_size + inf] +
      b * y_array[index_y * x_size + sup] +
      ((a*a*a-a)* ddy_array[index_y * x_size + inf] +
       (b*b*b-b)* ddy_array[index_y * x_size + sup])*h*h/6.;

  }

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are in different arrays
  *
  *
  */
int array_interpolate_extrapolate_logspline_loglinear_one_column(
								 double * x_array,
								 int x_size,
								 int x_stop,
								 double * y_array, /* array of size x_size*y_size with elements
										      y_array[index_y*x_size+index_x] */
								 int y_size,
								 int index_y,
								 double * ddlogy_array, /* array of size x_size*y_size */
								 double x,   /* input */
								 double * y, /* output */
								 ErrorMsg errmsg
								 ) {


  int inf,sup,mid;
  double h,a,b;

  if (x > x_array[x_stop-1]) {

    /*interpolate/extrapolate linearly ln(y) as a function of ln(x)*/

    h = log(x_array[x_stop-1]) - log(x_array[x_stop-2]);
    b = (log(x)-log(x_array[x_stop-2]))/h;
    a = 1-b;

/*     *y = exp(a * log(y_array[index_y * x_size + (x_stop-2)]) + */
/* 	     b * log(y_array[index_y * x_size + (x_stop-1)])); */

    *y = exp(log(y_array[index_y * x_size + (x_stop-1)])
	     +(log(x)-log(x_array[x_stop-1]))
	     *((log(y_array[index_y * x_size + (x_stop-1)])-log(y_array[index_y * x_size + (x_stop-2)]))/h
	       +h/6.*(ddlogy_array[index_y * x_size + (x_stop-2)]+2.*ddlogy_array[index_y * x_size + (x_stop-1)])));


  }

  else {

    /*interpolate ln(y) as a function of ln(x) with a spline*/

    inf=0;
    sup=x_stop-1;

    if (x_array[inf] < x_array[sup]){

      if (x < x_array[inf]) {
	class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[inf]);
	return _FAILURE_;
      }

      if (x > x_array[sup]) {
	class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[sup]);
	return _FAILURE_;
      }

      while (sup-inf > 1) {

	mid=(int)(0.5*(inf+sup));
	if (x < x_array[mid]) {sup=mid;}
	else {inf=mid;}

      }

    }

    else {

      if (x < x_array[sup]) {
	class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[sup]);
	return _FAILURE_;
      }

      if (x > x_array[inf]) {
	class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[inf]);
	return _FAILURE_;
      }

      while (sup-inf > 1) {

	mid=(int)(0.5*(inf+sup));
	if (x > x_array[mid]) {sup=mid;}
	else {inf=mid;}

      }

    }

    h = log(x_array[sup]) - log(x_array[inf]);
    b = (log(x)-log(x_array[inf]))/h;
    a = 1-b;

    *y = exp(a * log(y_array[index_y * x_size + inf]) +
	     b * log(y_array[index_y * x_size + sup]) +
	     ((a*a*a-a)* ddlogy_array[index_y * x_size + inf] +
	      (b*b*b-b)* ddlogy_array[index_y * x_size + sup])*h*h/6.);

  }

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array, x is arranged in growing order, and the point x is presumably close to the previous point x from the last call of this function.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_growing_closeby(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   double x,
		   int * last_index,
		   double * result,
		   int result_size, /** from 1 to n_columns */
		   ErrorMsg errmsg) {

  int inf,sup,i;
  double weight;

  inf = *last_index;
  sup = *last_index+1;

  while (x < *(array+inf*n_columns+index_x)) {
    inf--;
    if (inf < 0) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,
	      x,array[index_x]);
      return _FAILURE_;
    }
  }
  sup = inf+1;
  while (x > *(array+sup*n_columns+index_x)) {
    sup++;
    if (sup > (n_lines-1)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
	      x,array[(n_lines-1)*n_columns+index_x]);
      return _FAILURE_;
    }
  }
  inf = sup-1;

  *last_index = inf;

  weight=(x-*(array+inf*n_columns+index_x))/(*(array+sup*n_columns+index_x)-*(array+inf*n_columns+index_x));

  for (i=0; i<result_size; i++)
    *(result+i) = *(array+inf*n_columns+i) * (1.-weight)
      + weight * *(array+sup*n_columns+i);

  *(result+index_x) = x;

  return _SUCCESS_;
}

/**
  * interpolate to get y(x), when x and y are two columns of the same array, x is arranged in growing order, and the point x is presumably close to the previous point x from the last call of this function.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_one_growing_closeby(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   double x,
		   int * last_index,
       int index_y,
		   double * result,
		   ErrorMsg errmsg) {

  int inf,sup;
  double weight;

  inf = *last_index;
  sup = *last_index+1;

  while (x < *(array+inf*n_columns+index_x)) {
    inf--;
    if (inf < 0) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,
	      x,array[index_x]);
      return _FAILURE_;
    }
  }
  sup = inf+1;
  while (x > *(array+sup*n_columns+index_x)) {
    sup++;
    if (sup > (n_lines-1)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
	      x,array[(n_lines-1)*n_columns+index_x]);
      return _FAILURE_;
    }
  }
  inf = sup-1;

  *last_index = inf;

  weight=(x-*(array+inf*n_columns+index_x))/(*(array+sup*n_columns+index_x)-*(array+inf*n_columns+index_x));

  *result = *(array+inf*n_columns+index_y) * (1.-weight) + *(array+sup*n_columns+index_y) * weight;

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array, x is arranged in growing order, and the point x is presumably very close to the previous point x from the last call of this function.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_spline_growing_closeby(
					     double * x_array,
					     int n_lines,
					     double * array,
					     double * array_splined,
					     int n_columns,
					     double x,
					     int * last_index,
					     double * result,
					     int result_size, /** from 1 to n_columns */
					     ErrorMsg errmsg) {

  int inf,sup,i;
  double h,a,b;

  /*
  if (*last_index < 0) {
    class_sprintf(errmsg,"%s(L:%d) problem with last_index =%d < 0",__func__,__LINE__,*last_index);
    return _FAILURE_;
  }
  if (*last_index > (n_lines-1)) {
    class_sprintf(errmsg,"%s(L:%d) problem with last_index =%d > %d",__func__,__LINE__,*last_index,n_lines-1);
    return _FAILURE_;
  }
  */

  inf = *last_index;
  class_test(inf<0 || inf>(n_lines-1),
	     errmsg,
	     "*lastindex=%d out of range [0:%d]\n",inf,n_lines-1);
  while (x < x_array[inf]) {
    inf--;
    if (inf < 0) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,
	      x,x_array[0]);
      return _FAILURE_;
    }
  }
  sup = inf+1;
  while (x > x_array[sup]) {
    sup++;
    if (sup > (n_lines-1)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
	      x,x_array[n_lines-1]);
      return _FAILURE_;
    }
  }
  inf = sup-1;

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  for (i=0; i<result_size; i++)
    *(result+i) =
      a * *(array+inf*n_columns+i) +
      b * *(array+sup*n_columns+i) +
      ((a*a*a-a)* *(array_splined+inf*n_columns+i) +
       (b*b*b-b)* *(array_splined+sup*n_columns+i))*h*h/6.;

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array, x is arranged in growing order, and the point x is presumably close (but maybe not so close) to the previous point x from the last call of this function.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_spline_growing_hunt(
					     double * x_array,
					     int n_lines,
					     double * array,
					     double * array_splined,
					     int n_columns,
					     double x,
					     int * last_index,
					     double * result,
					     int result_size, /** from 1 to n_columns */
					     ErrorMsg errmsg) {

  int inf,sup,mid,i,inc;
  double h,a,b;

  inc=1;

  if (x >= x_array[*last_index]) {
    if (x > x_array[n_lines-1]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
	      x,x_array[n_lines-1]);
      return _FAILURE_;
    }
    /* try closest neighboor upward */
    inf = *last_index;
    sup = inf + inc;
    if (x > x_array[sup]) {
      /* hunt upward */
      while (x > x_array[sup]) {
	inf = sup;
	inc += 1;
	sup += inc;
	if (sup > n_lines-1) {
	  sup = n_lines-1;
	}
      }
      /* bisect */
      while (sup-inf > 1) {
	mid=(int)(0.5*(inf+sup));
	if (x < x_array[mid]) {sup=mid;}
	else {inf=mid;}
      }
    }
   }
  else {
    if (x < x_array[0]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,
	      x,x_array[0]);
      return _FAILURE_;
    }
    /* try closest neighboor downward */
    sup = *last_index;
    inf = sup - inc;
    if (x < x_array[inf]) {
      /* hunt downward */
      while (x < x_array[inf]) {
	sup = inf;
	inc += 1;
	inf -= inc;
	if (inf < 0) {
	  inf = 0;
	}
      }
      /* bisect */
      while (sup-inf > 1) {
	mid=(int)(0.5*(inf+sup));
	if (x < x_array[mid]) {sup=mid;}
	else {inf=mid;}
      }
    }
  }

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  for (i=0; i<result_size; i++)
    *(result+i) =
      a * *(array+inf*n_columns+i) +
      b * *(array+sup*n_columns+i) +
      ((a*a*a-a)* *(array_splined+inf*n_columns+i) +
       (b*b*b-b)* *(array_splined+sup*n_columns+i))*h*h/6.;

  return _SUCCESS_;
}


#define COEFF_MIN 1e-20
/**
 * @brief Solves the equation S(x) - y = 0 for distinct non-multiple roots.
 *        There needs to be at least one control point higher and one lower than y.
 * @param x_array     Input: contains x-values
 * @param x_size      Input: number of used control points
 * @param y_array     Input: contains y-values indexed as y_array[index_x] with size x_size
 * @param ddy_array   Input: contains y''-values retrieved from splining with size x_size
 * @param y           Input: y-value for which to find a solution x
 * @param x           In/Output: solution x (> or <) x0, should be initialised to a starting value x0
 * @param ascending   Input: if _TRUE_, search in ascending x-direction, else descending
 * @param last        Output: index of the start of the spline segment with y_{i} <= y <= y_{i+1}
 * @param errmsg
 *
 * @return the error status
 */
int array_spline_solve_table_lines(
                      double * x_array,
                      const int x_size,
                      double * y_array,
                      double * ddy_array,
                      const double y,
                      double * x,
                      const int ascending,
                      int * last,
                      ErrorMsg errmsg) {

  int index_y_inf, index_x_inf = 0, k;
  double h, a, b, c, d, sol;
  ErrorMsg err_hunt;

  /** - get the index of x0 in x_array */
  class_call(array_hunt_ascending(x_array, x_size, *x, &index_x_inf, err_hunt),
              err_hunt, errmsg);

  /** - get the index of y in y_array */
  if (ascending) {
    index_y_inf = index_x_inf;
    /** - find the sign change ascending from x0 */
    while ((y_array[index_y_inf] - y)*(y_array[index_y_inf+1] - y) > 0) {
      index_y_inf++;
      if (index_y_inf > x_size-2) {
        class_sprintf(errmsg, "%s(L:%d) no sign change in spline segments with x > %.4e detected",__func__,__LINE__, x_array[index_x_inf]);
        return _FAILURE_;
      }
    }
  } else {
    index_y_inf = index_x_inf + 1;
    /** - find the sign change descending from x0 */
    while ((y_array[index_y_inf] - y)*(y_array[index_y_inf+1] - y) > 0) {
      index_y_inf--;
      if (index_y_inf < 0) {
        class_sprintf(errmsg, "%s(L:%d) no sign change in spline segments with x < %.4e detected",__func__,__LINE__, x_array[index_x_inf+1]);
        return _FAILURE_;
      }
    }
  }

  /** - compute coefficients of the polynomial */
  h = x_array[index_y_inf+1] - x_array[index_y_inf];
  a = h*h/6. * (ddy_array[index_y_inf+1] - ddy_array[index_y_inf]);
  b = h*h/2. * ddy_array[index_y_inf];
  c = y_array[index_y_inf+1] - y_array[index_y_inf] - h*h/6. * (ddy_array[index_y_inf+1] + 2.*ddy_array[index_y_inf]);
  d = y_array[index_y_inf] - y;

  if (fabs(a) > COEFF_MIN) {
    /** - cubic polynomial transform to depressed form t^3 + pt + q = 0 */
    double p = (3.*a*c - b*b) / (3.*a*a);
    double q = (2.*b*b*b - 9.*a*b*c + 27.*a*a*d) / (27.*a*a*a);
    double delta = p*p*p/27. + q*q/4.; /**< discriminant */
    double rootp = sqrt(abs(p)/3.);

    if (delta <= 0) {
      /** - trigonometric root solutions */
      double arccos = acos(3.*q/(2.*p*rootp))/3.;
      for (k = 0; k < 3; k++) {
        sol = 2.*rootp * cos(arccos - 2.*_PI_*k/3.) - b/(3.*a);
        if ((sol >= 0) && (sol <= 1)) {
          break;
        }
      }
    }
    else {
      /** - hyperbolic root solutions */
      if (p < 0) {
        sol = -2. * ((q>0.) - (q<0.)) * rootp * cosh( acosh(-3.*abs(q)/(2.*p*rootp))/3. ) - b/(3.*a);
      }
      else {
        sol = -2.*rootp * sinh( asinh(3.*q/(2.*p*rootp))/3. ) - b/(3.*a);
      }
    }
  }
  else if (fabs(b) > COEFF_MIN) {
    /** - quadratic polynomial */
    sol = -(c + ((c>=0.) - (c<0.)) * sqrt(c*c - 4.*b*d))/(2.*b);
    if ((sol < 0) || (sol > 1)) {
      if (fabs(d) > COEFF_MIN) {
        sol = d / (b*sol);  /**< Vietas formula */
      }
      else {
        sol = 0.; /**< second root is zero for bx^2 + c*x = 0 */
      }
    }
  }
  else {
    /** - linear solution */
    sol = -d/c;
  }
  *x = x_array[index_y_inf] + sol * h;
  *last = index_y_inf;


  return _SUCCESS_;
}


// [NS]
/**
 * Get the index in the array, and the relative offset,
 *  but do not yet actually interpolate
 */
int array_spline_hunt(double* x_array,
                       int x_size,
                       double x,
                       int* last,
                       double* h,
                       double* a,
                       double* b,
                       ErrorMsg errmsg){
  /* Define local quantities */
  int inf,sup,mid,inc;
  int last_index = *last;

  /* Error checking */
  if(last_index>=x_size-1){last_index=x_size-2;}
  if(last_index<0){last_index=0;}

  /* Hunt ! */
  inc=1;
  if (x >= x_array[last_index]) {
    if (x > x_array[x_size-1]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
        x,x_array[x_size-1]);
      return _FAILURE_;
    }
    /* try closest neighboor upward */
    inf = last_index;
    sup = inf + inc;
    if (x > x_array[sup]) {
      /* hunt upward */
      while (x > x_array[sup]) {
        inf = sup;
        inc += 1;
        sup += inc;
        if (sup > x_size-1) {
          sup = x_size-1;
        }
      }
      /* bisect */
      while (sup-inf > 1) {
        mid=(int)(0.5*(inf+sup));
        if (x < x_array[mid]) {sup=mid;}
        else {inf=mid;}
      }
    }
  }
  else {
    if (x < x_array[0]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%.20e < x_min=%.20e",__func__,__LINE__,
        x,x_array[0]);
      return _FAILURE_;
    }
    /* try closest neighboor downward */
    sup = last_index;
    inf = sup - inc;
    if (x < x_array[inf]) {
      /* hunt downward */
      while (x < x_array[inf]) {
        sup = inf;
        inc += 1;
        inf -= inc;
        if (inf < 0) {
          inf = 0;
        }
      }
      /* bisect */
      while (sup-inf > 1) {
        mid=(int)(0.5*(inf+sup));
        if (x < x_array[mid]) {sup=mid;}
        else {inf=mid;}
      }
    }
  }

  /* We have found the prey */
  last_index = inf;
  *last = last_index;
  *h = x_array[sup] - x_array[inf];
  *b = (x-x_array[inf])/(*h);
  *a = 1.0-(*b);

  return _SUCCESS_;
}

/**
 * interpolate linearily to get y_i(x), when x and y_i are in two different arrays
 *
 * Called by transfer_interpolate_sources(); transfer_functions_at_k(); perturbations_sources_at_eta().
 */
int array_interpolate_two(
		   double * array_x,
		   int n_columns_x,
		   int index_x,   /** from 0 to (n_columns_x-1) */
		   double * array_y,
		   int n_columns_y,
		   int n_lines,  /** must be the same for array_x and array_y */
		   double x,
		   double * result,
		   int result_size, /** from 1 to n_columns_y */
		   ErrorMsg errmsg) {

  int inf,sup,mid,i;
  double weight;

  inf=0;
  sup=n_lines-1;

  if (array_x[inf*n_columns_x+index_x] < array_x[sup*n_columns_x+index_x]){

    if (x < array_x[inf*n_columns_x+index_x]) {

      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,array_x[inf*n_columns_x+index_x]);
      return _FAILURE_;
    }

    if (x > array_x[sup*n_columns_x+index_x]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,array_x[sup*n_columns_x+index_x]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < array_x[mid*n_columns_x+index_x]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < *(array_x+sup*n_columns_x+index_x)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,*(array_x+sup*n_columns_x+index_x));
      return _FAILURE_;
    }

    if (x > *(array_x+inf*n_columns_x+index_x)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,*(array_x+inf*n_columns_x+index_x));
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > *(array_x+mid*n_columns_x+index_x)) {sup=mid;}
      else {inf=mid;}

    }

  }

  weight=(x-*(array_x+inf*n_columns_x+index_x))/(*(array_x+sup*n_columns_x+index_x)-*(array_x+inf*n_columns_x+index_x));

  for (i=0; i<result_size; i++)
    *(result+i) = *(array_y+i*n_lines+inf) * (1.-weight)
      + weight * *(array_y+i*n_lines+sup) ;

  return _SUCCESS_;
}

/**
 * Same as array_interpolate_two, but with order of indices exchanged in array_y
 */
int array_interpolate_two_bis(
		   double * array_x,
		   int n_columns_x,
		   int index_x,   /** from 0 to (n_columns_x-1) */
		   double * array_y,
		   int n_columns_y,
		   int n_lines,  /** must be the same for array_x and array_y */
		   double x,
		   double * result,
		   int result_size, /** from 1 to n_columns_y */
		   ErrorMsg errmsg) {

  int inf,sup,mid,i;
  double weight;

  inf=0;
  sup=n_lines-1;

  if (array_x[inf*n_columns_x+index_x] < array_x[sup*n_columns_x+index_x]){

    if (x < array_x[inf*n_columns_x+index_x]) {

      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,array_x[inf*n_columns_x+index_x]);
      return _FAILURE_;
    }

    if (x > array_x[sup*n_columns_x+index_x]) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,array_x[sup*n_columns_x+index_x]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < array_x[mid*n_columns_x+index_x]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < *(array_x+sup*n_columns_x+index_x)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,*(array_x+sup*n_columns_x+index_x));
      return _FAILURE_;
    }

    if (x > *(array_x+inf*n_columns_x+index_x)) {
      class_sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,*(array_x+inf*n_columns_x+index_x));
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > *(array_x+mid*n_columns_x+index_x)) {sup=mid;}
      else {inf=mid;}

    }

  }

  weight=(x-*(array_x+inf*n_columns_x+index_x))/(*(array_x+sup*n_columns_x+index_x)-*(array_x+inf*n_columns_x+index_x));

  for (i=0; i<result_size; i++)
    *(result+i) = *(array_y+inf*n_columns_y+i) * (1.-weight)
      + weight * *(array_y+sup*n_columns_y+i) ;

  return _SUCCESS_;
}


/**
 * interpolate linearily to get y_i(x), when x and y_i are in two different arrays
 *
 * Called by transfer_interpolate_sources(); transfer_functions_at_k(); perturbations_sources_at_eta().
 */
int array_interpolate_two_arrays_one_column(
					    double * array_x, /* assumed to be a vector (i.e. one column array) */
					    double * array_y,
					    int n_columns_y,
					    int index_y, /* between 0 and (n_columns_y-1) */
					    int n_lines,  /** must be the same for array_x and array_y */
					    double x,
					    double * result,
					    ErrorMsg errmsg) {

  int inf,sup,mid;
  double weight;
  double epsilon=1e-9;

  inf=0;
  sup=n_lines-1;

  if (array_x[inf] < array_x[sup]){

    class_test(x < array_x[inf]-epsilon,
	       errmsg,
	       "x=%e < x_min=%e",x,array_x[inf]);

    class_test(x > array_x[sup]+epsilon,
	       errmsg,
	       "x=%e > x_max=%e",x,array_x[sup]);

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < array_x[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    class_test(x < array_x[sup]-epsilon,
	       errmsg,
	       "x=%e < x_min=%e",x,array_x[sup]);

    class_test(x > array_x[inf]+epsilon,
	       errmsg,
	       "x=%e > x_max=%e",x,array_x[inf]);

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > array_x[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  weight=(x-array_x[inf])/(array_x[sup]-array_x[inf]);

  *result = array_y[index_y*n_lines+inf] * (1.-weight)
    + weight * array_y[index_y*n_lines+sup];

  return _SUCCESS_;
}

/**
 * Called by transfer_solve().
 */
int array_interpolate_equal(
			    double * array,
			    int n_columns,
			    int n_lines,
			    double x,
			    double x_min,
			    double x_max,
			    double * result,
			    ErrorMsg errmsg) {

  int index_minus,i;
  double x_step,x_minus,weight;

  if (x < x_min) {
    class_sprintf(errmsg,"%s(L:%d) : x out of bounds: x=%e,x_min=%e",__func__,__LINE__,x,x_min);
    return _FAILURE_;
  }

  if (x > x_max) {
    class_sprintf(errmsg,"%s(L:%d) : x out of bounds: x=%e,x_max=%e",__func__,__LINE__,x,x_max);
    return _FAILURE_;
  }

  x_step = (x_max-x_min)/(n_lines-1);
  index_minus = (int)((x-x_min)/x_step);
  x_minus = index_minus * x_step;
  weight = (x-x_minus) / x_step;

  for (i=0; i<n_columns; i++)
    result[i] = *(array+n_columns*index_minus+i)*(1.-weight)
      + *(array+n_columns*(index_minus+1)+i)*weight;

  return _SUCCESS_;

}

/**
 * cubic interpolation of array with equally space abscisses
 */

int array_interpolate_cubic_equal(
				  double x0,
				  double dx,
				  double *yarray,
				  int Nx,
				  double x,
				  double * result,
				  ErrorMsg errmsg) {

  int i;
  double frac;

  class_test((dx > 0 && (x<x0 || x>x0+dx*(Nx-1))),
	     errmsg,
	     "x=%e out of range [%e %e]",x,x0,x0+dx*(Nx-1));

  class_test((dx < 0 && (x>x0 || x<x0+dx*(Nx-1))),
	     errmsg,
	     "x=%e out of range [%e %e]",x,x0+dx*(Nx-1),x0);

  i = (int)floor((x-x0)/dx);
  if (i<1) i=1;
  if (i>Nx-3) i=Nx-3;
  frac = (x-x0)/dx-i;
  yarray += i-1;

  *result=-yarray[0]*frac*(1.-frac)*(2.-frac)/6.
    +yarray[1]*(1.+frac)*(1.-frac)*(2.-frac)/2.
    +yarray[2]*(1.+frac)*frac*(2.-frac)/2.
    +yarray[3]*(1.+frac)*frac*(frac-1.)/6.;

  return _SUCCESS_;
}

int array_interpolate_parabola(double x1,
			       double x2,
			       double x3,
			       double x,
			       double y1,
			       double y2,
			       double y3,
			       double * y,
			       double * dy,
			       double * ddy,
			       ErrorMsg errmsg) {

  double a,b,c;

  /*
    a x_i**2 + b x_i + c = y_i

    a (x1**2-x2**2) + b (x1-x2) = y1-y2
    a (x3**2-x2**2) + b (x3-x2) = y3-y2

    a (x1**2-x2**2)(x3**2-x2**2) + b (x1-x2)(x3**2-x2**2) = (y1-y2)(x3**2-x2**2)
    a (x3**2-x2**2)(x1**2-x2**2) + b (x3-x2)(x1**2-x2**2) = (y3-y2)(x1**2-x2**2)

    b = [(y1-y2)(x3**2-x2**2) - (y3-y2)(x1**2-x2**2)]/(x1-x2)(x3-x2)(x3-x1)

  */

  b = ((y1-y2)*(x3-x2)*(x3+x2) - (y3-y2)*(x1-x2)*(x1+x2))/(x1-x2)/(x3-x2)/(x3-x1);

  a = (y1-y2-b*(x1-x2))/(x1-x2)/(x1+x2);

  c = y2 - b*x2 - a*x2*x2;

  *y = a*x*x + b*x + c;
  *dy = 2.*a*x + b;
  *ddy = 2.*a;

  return _SUCCESS_;

}

int array_convert_spline_table_columns_to_local_power_basis(
              const double * const x, /* vector of size x_size */
		          const int x_size,
		          const double * const y_array, /* array of size x_size*y_size with elements
					                                      y_array[index_y*x_size+index_x] */
		          const int y_size,
		          const double * const ddy_array, /* array of size x_size*y_size */
              double * coefficients,  /* array of coefficients of size y_size*(x_size-1)*4, C-order */
              double * breakpoints) { /* array of breakpoints of size y_size*x_size */

  int index_x, index_y;
  double h;

  memcpy(breakpoints, x, x_size*sizeof(double));

  for (index_x = 0; index_x < x_size-1; index_x++) {
    h = x[index_x + 1] - x[index_x];
    for (index_y = 0; index_y < y_size; index_y++) {
      /** - coefficients in order of decreasing powers */
      coefficients[0*(x_size-1)*y_size + index_x*y_size + index_y] = (ddy_array[index_y*x_size + (index_x+1)] - ddy_array[index_y*x_size + index_x])/(6.*h);
      coefficients[1*(x_size-1)*y_size + index_x*y_size + index_y] = ddy_array[index_y*x_size + index_x]/2.;
      coefficients[2*(x_size-1)*y_size + index_x*y_size + index_y] = (y_array[index_y*x_size + (index_x+1)] - y_array[index_y*x_size + index_x])/h     \
                                                                    - h/6.*(2.*ddy_array[index_y*x_size + index_x] + ddy_array[index_y*x_size + (index_x+1)]);
      coefficients[3*(x_size-1)*y_size + index_x*y_size + index_y] = y_array[index_y*x_size + index_x];
    }
  }

  return _SUCCESS_;
}

/**
 * Called by transfer_solve().
 */
int array_integrate_all(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y,
		   double *result) {

  int i;
  double sum;

  sum=0.;

  for (i=1; i<n_lines; i++) {

    sum += 0.5 * (*(array+i*n_columns+index_y) + *(array+(i-1)*n_columns+index_y))
               * (*(array+i*n_columns+index_x) - *(array+(i-1)*n_columns+index_x));

  }

  *result = sum;

  return _SUCCESS_;

}

int array_smooth_trg(double * array,
		     int k_size,
		     int starting_k,
		     int eta_size,
		     int index_eta,
		     int radius, /*3, 5 or 7 */
		     ErrorMsg errmsg) {

  double * smooth;
  int i,j,jmin,jmax;
  double weigth;
  double *coeff;

  smooth=(double*)malloc(k_size*sizeof(double));
  if (smooth == NULL) {
    class_sprintf(errmsg,"%s(L:%d) Cannot allocate smooth",__func__,__LINE__);
    return _FAILURE_;
  }

  class_calloc(coeff,2*radius+1,sizeof(double),errmsg);

  switch(radius){
  case 3:
    weigth = 21;

    coeff[0] = -2;
    coeff[1] = 3;
    coeff[2] = 6;
    coeff[3] = 7;
    coeff[4] = 6;
    coeff[5] = 3;
    coeff[6] = -2;

    break;
  case 4:
    weigth = 231;

    coeff[0] = -21;
    coeff[1] = 14;
    coeff[2] = 39;
    coeff[3] = 54;
    coeff[4] = 59;
    coeff[5] = 54;
    coeff[6] = 39;
    coeff[7] = 14;
    coeff[8] = -21;

    break;
  case 5:
    weigth = 429;

    coeff[0] = -36;
    coeff[1] = 9;
    coeff[2] = 44;
    coeff[3] = 69;
    coeff[4] = 84;
    coeff[5] = 89;
    coeff[6] = 84;
    coeff[7] = 69;
    coeff[8] = 44;
    coeff[9] = 9;
    coeff[10] = -36;

    break;
  case 6:
    weigth = 143;

    coeff[0] = -11;
    coeff[1] = 0;
    coeff[2] = 9;
    coeff[3] = 16;
    coeff[4] = 21;
    coeff[5] = 24;
    coeff[6] = 25;
    coeff[7] = 24;
    coeff[8] = 21;
    coeff[9] = 16;
    coeff[10] = 9;
    coeff[11] = 0;
    coeff[12] = -11;

    break;
  case 7:
    weigth = 1105;

    coeff[0] = -78;
    coeff[1] = -13;
    coeff[2] = 42;
    coeff[3] = 87;
    coeff[4] = 122;
    coeff[5] = 147;
    coeff[6] = 162;
    coeff[7] = 167;
    coeff[8] = 162;
    coeff[9] = 147;
    coeff[10] = 122;
    coeff[11] = 87;
    coeff[12] = 42;
    coeff[13] = -13;
    coeff[14] = -78;

    break;

/*   case 8: */


  default:
    class_stop(errmsg,"Non valid radius %d: please chose between 3 4 5 or 6\n",radius);
    weigth=0;
    break;
  }

  for (i=starting_k; i<k_size-radius; i++) {
      smooth[i]=0.;
      jmin = MAX(i-radius,0);
      jmax = MIN(i+radius,k_size-1);
      for (j=jmin; j <= jmax; j++) {
	smooth[i] += coeff[j-jmin]*array[j+k_size*index_eta];
      }
      smooth[i] /= weigth;
  }

  for (i=starting_k; i<k_size-radius; i++)
    array[i+k_size*index_eta] = smooth[i];

  free(smooth);
  free(coeff);

  return _SUCCESS_;

}

int array_smooth(double * array,
		 int n_columns,
		 int n_lines,
		 int index, /** from 0 to (n_columns-1) */
		 int radius,
		 ErrorMsg errmsg) {

  double * smooth;
  int i,j,jmin,jmax;
  double weigth;

  smooth=(double*)malloc(n_lines*sizeof(double));
  if (smooth == NULL) {
    class_sprintf(errmsg,"%s(L:%d) Cannot allocate smooth",__func__,__LINE__);
    return _FAILURE_;
  }

  for (i=0; i<n_lines; i++) {
    smooth[i]=0.;
    weigth=0.;
    jmin = MAX(i-radius,0);
    jmax = MIN(i+radius,n_lines-1);
    for (j=jmin; j <= jmax; j++) {
      smooth[i] += array[j*n_columns+index];
      weigth += 1.;
    }
    smooth[i] /= weigth;
  }

  for (i=0; i<n_lines; i++)
    array[i*n_columns+index] = smooth[i];

  free(smooth);

  return _SUCCESS_;

}

/**
 * Compute quadrature weights for the trapezoidal integration method, xhen x is in gorwing order.
 *
 * @param x                     Input: Grid points on which f() is known.
 * @param n                     Input: number of grid points.
 * @param w_trapz               Output: Weights of the trapezoidal method.
 * @return the error status
 */

int array_trapezoidal_weights(
                              double * __restrict__ x,
                              int n,
                              double * __restrict__ w_trapz,
                              ErrorMsg errmsg
                              ) {
  int i;

  /* Case with just one point, w would normally be 0. */
  if (n==1){
    w_trapz[0] = 0.0;
  }
  else if (n>1){
    //Set edgeweights:
    w_trapz[0] = 0.5*(x[1]-x[0]);
    w_trapz[n-1] = 0.5*(x[n-1]-x[n-2]);
    //Set inner weights:
    for (i=1; i<(n-1); i++){
      w_trapz[i] = 0.5*(x[i+1]-x[i-1]);
    }
  }
  return _SUCCESS_;
}

/**
 * Compute quadrature weights for the trapezoidal integration method, when x is in decreasing order.
 *
 * @param x                     Input: Grid points on which f() is known.
 * @param n                     Input: number of grid points.
 * @param w_trapz               Output: Weights of the trapezoidal method.
 * @return the error status
 */

int array_trapezoidal_mweights(
                              double * __restrict__ x,
                              int n,
                              double * __restrict__ w_trapz,
                              ErrorMsg errmsg
                              ) {
  int i;

  /* Case with just one point. */
  if (n==1){
    w_trapz[0] = 1.0;
  }
  else if (n>1){
    //Set edgeweights:
    w_trapz[0] = 0.5*(x[0]-x[1]);
    w_trapz[n-1] = 0.5*(x[n-2]-x[n-1]);
    //Set inner weights:
    for (i=1; i<(n-1); i++){
      w_trapz[i] = 0.5*(x[i-1]-x[i+1]);
    }
  }
  return _SUCCESS_;
}

/**
 * Compute integral of function using trapezoidal method.
 *
 * @param integrand             Input: The function we are integrating.
 * @param n                     Input: Compute integral on grid [0;n-1].
 * @param w_trapz               Input: Weights of the trapezoidal method.
 * @param I                     Output: The integral.
 * @return the error status
 */

int array_trapezoidal_integral(
                                  double * __restrict__ integrand,
                                  int n,
                                  double * __restrict__ w_trapz,
                                  double * __restrict__ I,
                                  ErrorMsg errmsg
                                  ) {
  int i;
  double res=0.0;
  for (i=0; i<n; i++){
    res += integrand[i]*w_trapz[i];
  }
  *I = res;
  return _SUCCESS_;
}

/**
 * Compute convolution integral of product of two functions using trapezoidal method.
 *
 * @param integrand1            Input: Function 1.
 * @param integrand2            Input: Function 2.
 * @param n                     Input: Compute integral on grid [0;n-1].
 * @param w_trapz               Input: Weights of the trapezoidal method.
 * @param I                     Output: The integral.
 * @return the error status
 */

int array_trapezoidal_convolution(
                                     double * __restrict__ integrand1,
                                     double * __restrict__ integrand2,
                                     int n,
                                     double * __restrict__ w_trapz,
                                     double * __restrict__ I,
                                     ErrorMsg errmsg
                                     ) {
  int i;
  double res=0.0;
  for (i=0; i<n; i++){
    res += integrand1[i]*integrand2[i]*w_trapz[i];
  }
  *I = res;
  return _SUCCESS_;
}

/**
 * In general, to obtain a least-squared fit to N data points,
 * the matrix average(x^(i+j)) * a_i = average(x^j y) has to be solved,
 * where i,j are the individual matrix indices as powers, and
 * the averages are over all N data points.
 * Here we implement the case for 3 coefficients a_i explicitly,
 * using matrix calculations according to Cramer's rule.
 * (Instead of a full LU decomposition)
 *
 **/
int array_extrapolate_quadratic(double* x, double* y, double xnew, int x_size, double* ynew, double* dynew, ErrorMsg errmsg){

  int i;

  double * xarr = x;
  double * yarr = y;

  double av1=x_size;
  double avx=0.0, avxx=0.0, avxxx=0.0, avxxxx=0.0;
  double avy=0.0, avyx=0.0, avyxx=0.0;

  double div,a,b,c,xval,yval;

  /*
   * We pivot around the zero-th element
   * This removes natural offsets and scales
   * (i.e. transforms it to around unity for x and y)
   * This usually prevents numerical cancelation (offsets) and over/underflow (scales)
   */
  for(i=0;i<x_size;++i){
    xval = (xarr[i]-xarr[0])/xarr[0];
    yval = (yarr[i]-yarr[0])/yarr[0];
    avx += xval;
    avxx += xval*xval;
    avxxx += xval*xval*xval;
    avxxxx += xval*xval*xval*xval;
    avy += yval;
    avyx += yval*xval;
    avyxx += yval*xval*xval;
  }

  div = avxxxx*(avxx*av1-avx*avx) + avxxx*(avxx*avx-avxxx*av1) + avxx*(avxxx*avx-avxx*avxx);
  class_test(div == 0.0,
             errmsg,
             "Cannot extrapolate at x = %g for the given data set",xnew);
  a = avyxx*(avxx*av1-avx*avx) + avyx*(avxx*avx-avxxx*av1) + avy*(avxxx*avx-avxx*avxx);
  b = avyxx*(avxx*avx-avxxx*av1) + avyx*(avxxxx*av1-avxx*avxx) + avy*(avxxx*avxx-avx*avxxxx);
  c = avyxx*(avxxx*avx-avxx*avxx) + avyx*(avxxx*avxx-avxxxx*avx) + avy*(avxxxx*avxx-avxxx*avxxx);

  a/=div; b/=div; c/=div;

  xval = (xnew-xarr[0])/xarr[0];

  *ynew = yarr[0]+yarr[0]*(a*xval*xval + b*xval + c);
  *dynew = yarr[0]*(2.*a/xarr[0]*xval + b/xarr[0]);

  return _SUCCESS_;
}

/**
 * Assuming that array is a vector with elements arranged in descending
 * order, find i such that array[i] > value > array[i+1].
 */

int array_hunt_descending(
                          double * array,
                          int size,
                          double value,
                          int * index,
                          ErrorMsg errmsg
                          ) {
  int i_inf,i_sup,i_mid;

  i_inf=0;
  i_sup=size-1;

  /* checks */
  if (array[i_inf] < array[i_sup]) {
    class_sprintf(errmsg,"%s(L:%d) array is not in descending order (checked only the boundaries)",__func__,__LINE__);
    return _FAILURE_;
  }
  if ((value > array[i_inf]) || (value < array[i_sup])) {
    class_sprintf(errmsg,"%s(L:%d) %e is outside the range [%e, %e]",__func__,__LINE__,value,array[size-1],array[0]);
    return _FAILURE_;
  }

  /* bisection */
  while (i_sup-i_inf>1) {
    i_mid = (i_sup+i_inf)/2;
    if (value > array[i_mid])
      i_sup=i_mid;
    else
      i_inf=i_mid;
  }

  /* result */
  *index = i_inf;

  /* check for debug */
  /* routine never used and tested before, be careful */
  fprintf(stderr,
          "Check that array[%d]=%e>%e>array[%d]=%e\n",
          *index,
          array[*index],
          value,
          *index+1,
          array[*index+1]);
  return _FAILURE_;

  return _SUCCESS_;
}

/**
 * Assuming that array is a vector with elements arranged in ascending
 * order, find i such that array[i] < value < array[i+1].
 */

int array_hunt_ascending(
                          double * array,
                          int size,
                          double value,
                          int * index,
                          ErrorMsg errmsg
                          ) {
  int i_inf,i_sup,i_mid;

  i_inf=0;
  i_sup=size-1;

  /* checks */
  if (array[i_inf] > array[i_sup]) {
    class_sprintf(errmsg,"%s(L:%d) array is not in ascending order (checked only the boundaries)",__func__,__LINE__);
    return _FAILURE_;
  }
  if ((value < array[i_inf]) || (value > array[i_sup])) {
    class_sprintf(errmsg,"%s(L:%d) %e is outside the range [%e, %e]",__func__,__LINE__,value,array[0],array[size-1]);
    return _FAILURE_;
  }

  /* bisection */
  while (i_sup-i_inf>1) {
    i_mid = (i_sup+i_inf)/2;
    if (value > array[i_mid])
      i_inf=i_mid;
    else
      i_sup=i_mid;
  }

  /* result */
  *index = i_inf;

  /* check for debug */
  /*
  fprintf(stderr,
          "Check that array[%d]=%e<%e<array[%d]=%e\n",
          *index,
          array[*index],
          value,
          *index+1,
          array[*index+1]);
  return _FAILURE_;
  */

  return _SUCCESS_;
}
