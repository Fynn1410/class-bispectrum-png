
/** @file utilities.c Documented basic utility functions used by other modules of the code.
 *
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 *
 * In summary, the following functions can be called from other modules:
 * -# make_1Darray()  			dynamically allocates memory to a 1d array
 * -# make_2Darray()  			dynamically allocates memory to a 2d array
 * -# free_2Darray()  			free the memory allocated to a 2d array
 * -# init_1Darray()  			initialize a 1d array with linear spacing
 * -# loginit_1Darray()         initialize a 1d array with natural-log spacing
 * -# log10init_1Darray()       initialize a 1d array with log10 spacing
 * -# count_lines_in_file()     count the number of lines of a file
 * -# count_cols_in_file()      count number of columns of a file
 * -# return_arr()
 */


#include "header.h"
// #include <complex.h>
// #include <stdio.h>
// #include <math.h>
// #include <stdlib.h>


int _errno_util = 0;
double complex _err_last_arg = CMPLX(0., 0.);

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
double complex cGamma(const double complex z)
{
  const double g = 6.0246800407767295;
  const double co[LANCZOS_GAMMA_N] = {2.506628274631000,
      589.5106040667278, -888.0253935502019,
      395.8387847487511, -53.21395931507683, 
      1.277182848200161, -4.046172558017667e-4, 
      -7.347584327845915e-6, 8.208805790146655e-6, 
      -5.159543403041225e-6, 2.319631454949221e-6, 
      -6.671246136975432e-7, 9.060393467651553e-8};
  register double complex Lg;
  double complex t;
  double r;
  
  r = creal(z);
  if (r < 0.5)
  {
    /** - check the distance to the singularity */
    if (cabs(z - rint(r)) < _EPSILON_) {  /** - uses current rounding-mode, which defaults to round-to-nearest */
      _errno_util |= _ERR_RES_OUT_OF_RANGE_;
      _err_last_arg = z;
      #ifdef NAN
      return CMPLX(nan(""), nan(""));   /** Gamma(z) is undefined */
      #else
      return CMPLX(0., 0.);
      #endif
    }
    /** - use mirror identity */
    return _PI_ / (csin(_PI_ * z) * cGamma(1 - z));
  }

  Lg = co[0];
  for (int k = 1; k < LANCZOS_GAMMA_N; k++)
  { Lg += co[k] / (z + k - 1); }  /** Partial fraction sum */
  t = z + g - 0.5;
  return cpow(t, z - 0.5) * cexp(-t) * Lg;
} 


/**
 * @brief Real Gamma function approximated using the
 *        Lanczos method with g ~ 6 and N = 13.
 * 
 * @param x     Input: Real argument
 * @return      value of Gamma(x) / NaN at a singularity
*/
double rGamma(const double x)
{
  const double g = 6.0246800407767295;
  const double co[LANCZOS_GAMMA_N] = {2.506628274631000,
      589.5106040667278, -888.0253935502019,
      395.8387847487511, -53.21395931507683, 
      1.277182848200161, -4.046172558017667e-4, 
      -7.347584327845915e-6, 8.208805790146655e-6, 
      -5.159543403041225e-6, 2.319631454949221e-6, 
      -6.671246136975432e-7, 9.060393467651553e-8};
  register double Lg;
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




/**
 * Allocate memory to a 1d array of type double and length size
 *
 * @param size      Input: length of the array
 *
 * @return a pointer to a 1d array
 */
double *make_1Darray(long size)
{
	double *array;
	long i;

	array = (double*) malloc(size*sizeof(double));
	if (array==NULL) exit(2);
	for(i=0L;i<size;i++) array[i]= 0.0;
	return array;

}

/**
 * Allocate memory to a 1d array of type complex and length size
 *
 * @param size      Input: length of the array
 *
 * @return a pointer to a 1d array
 */
double complex *make_1D_c_array(long size)
{
	double complex *array;
	long i;

	array = (double complex*) malloc(size*sizeof(complex));
	if (array==NULL) exit(2);
	for(i=0L;i<size;i++) array[i]= 0.0;
	return array;

}

/**
 * Allocate memory to a 1d array of type integer and length size
 *
 * @param size      Input: length of the arrat
 *
 * @return a pointer to an integer type 1d array
 */
int *make_1D_int_array(long size)
{
	int *array;
	long i;

	array = (int*) malloc(size*sizeof(int));
	if (array==NULL) exit(2);
	for(i=0L;i<size;i++) array[i]= 0.0;

	return array;

}

/**
 * Allocate memory to a 2d array of type double
 *
 * @param nrows      Input: number of rows of the output array
 * @param ncols      Input: number of columns of the output array
 *
 * @return a double pointer to a double type 2d array
 */
 double complex **make_2D_c_array(long nrows, long ncolumns)
{
	long size;
	long i,j;
	double complex **array;

	size = nrows*ncolumns;

	array = (double complex **) malloc(nrows*sizeof(double complex *));
	if(array==NULL) exit(2);

	array[0] = (double complex *) malloc(size*sizeof(double complex));
	if(array[0]==NULL) exit(2);

	for(i=1L;i<nrows;i++)
	{
		array[i] = array[0] + i * ncolumns;
	}


	for(i=0L; i<nrows; i++)
	{
		for(j=0L;j<ncolumns;j++)
		{
			array[i][j] = 0.0;
		}
	}

	return array;
}

/**
 * Allocate memory to a 2d array of type double
 *
 * @param nrows      Input: number of rows of the output array
 * @param ncols      Input: number of columns of the output array
 *
 * @return a double pointer to a double type 2d array
 */
 double **make_2Darray(long nrows, long ncolumns)
{
	long size;
	long i,j;
	double **array;

	size = nrows*ncolumns;

	array = (double **) malloc(nrows*sizeof(double *));
	if(array==NULL) exit(2);

	array[0] = (double *) malloc(size*sizeof(double));
	if(array[0]==NULL) exit(2);

	for(i=1L;i<nrows;i++)
	{
		array[i] = array[0] + i * ncolumns;
	}


	for(i=0L; i<nrows; i++)
	{
		for(j=0L;j<ncolumns;j++)
		{
			array[i][j] = 0.0;
		}
	}

	return array;
}


/**
 * Free the memory allocated to a 2d array
 *
 * @param m      Input: double pointer to the elements of 2d array
 *
 * @return void
 */
void free_2Darray(double ** m)
{

	free(m[0]);
	free(m);

	return;
}

/**
 * Free the memory allocated to a 2d array
 *
 * @param m      Input: double pointer to the elements of 2d array
 *
 * @return void
 */
void free_2D_c_array(double complex ** m)
{

	free(m[0]);
	free(m);

	return;
}


/**
 * initialize a 1d array, with values in the range of [xmin,xmax] and evenely-space on linear scale
 *
 * @param n     	Input: number of elements
 * @param xmin      Input: start point
 * @param xmiax     Input: end point
 *
 * @return a pointer to a double type 1d array, with values initialized
 */
double *init_1Darray(long n,double xmin,double xmax)
{
	double inc;
	double *A;
	long i;

	A = make_1Darray(n+1L);

	inc = (xmax-xmin)/n;

	for(i=0L;i<n+1L;i++)
	{
		A[i] = xmin + inc * i;
		//printf("%12.6e\n",A[i]);
	}

	return A;
	}


/**
 * initialize a 1d array, with values in the range of [xmin,xmax] and evenely-space on natural-log scale
 *
 * @param n     	Input: number of elements
 * @param xmin      Input: start point
 * @param xmiax     Input: end point
 *
 * @return a pointer to a double type 1d array, with values initialized
 */
double *loginit_1Darray(long n,double xmin,double xmax)
{

	double inc;
	double *A;
	long i;

	A = make_1Darray(n);

	xmin = log(xmin);
	xmax = log(xmax);
	inc  = (xmax-xmin)/n;

	for(i=0L;i<n;i++)
	{
		A[i] = xmin + inc * i;
		A[i] = exp(A[i]);
	}

	return A;
}


/**
 * initialize a 1d array, with values in the range of [xmin,xmax] and evenely-space on log10 scale
 *
 * @param n     	Input: number of elements
 * @param xmin      Input: start point
 * @param xmiax     Input: end point
 *
 * @return a pointer to a double type 1d array, with values initialized
 */
double *log10init_1Darray(long n, double inc, double xmin)
{

	double *A, *logA;
	long i;

	double logxmin;
	logA = make_1Darray(n);
	A = make_1Darray(n);

	logxmin = log10(xmin);

	for(i=0L;i<n;i++)
	{
		logA[i] = logxmin + inc * i;
		A[i] = pow(10,logA[i]);
	}

	free(logA);

	return A;
}

/**
 * Count the number of lines of a file
 *
 * @param fname     	Input: filename
 * @return long integer value of nlines
 */
long count_lines_in_file(char *fname)
{
	FILE* ifp;
	long numlines = 0L;
	char c;
	int verbose = 1;
	char line[MAXL];

	ifp = fopen(fname, "r");
	if(ifp == NULL)
	{
    	fprintf( stderr, "Unable to open file \n" );
    	exit( 1 );
	}

	int numcoms = 0;
	while((c=fgetc(ifp)) != EOF)
	{
			if(c =='\n')
       			numlines++;

       		if(c =='#')
       			numcoms++;

	}
	if (verbose > 3) printf("Number of lines: %ld Number of comment lines: %d\n",numlines,numcoms);

	int data_len = 	numlines-numcoms;

	fclose(ifp);

	return data_len;
}


/**
 * Count the number of columns of a file
 *
 * @param fname     	Input: filename
 * @return long integer value of ncols
 */
long count_cols_in_file(char *fname)
{
	FILE* ifp;
	long nrows = 0L, ncols =1;
	char c;
	int verbose = 1;
	char line[MAXL];

	ifp = fopen(fname, "r");
	if(ifp == NULL)
	{
    	fprintf( stderr, "Unable to open file \n" );
    	exit( 1 );
	}

	int count = 0;
	while(fgets(line, sizeof line, ifp) )
	{
	 	if(*line == '#')  continue;
	 	while((c=fgetc(ifp)) != EOF)
		{
			if(c =='\n')
	       		nrows++;
			if (nrows == 2 && c== '\t')
				ncols ++;
		}
	}

	if (verbose > 3) printf("num_rows: %ld \t num_cols =%ld\n",nrows,ncols);
	fclose(ifp);

	return ncols;
}


void np_mat_fill(double complex (mat_func)(double complex, double complex), struct fft_struct *fft_input, long hm_switch, double complex **matrix)
{
	int Nmax = fft_input -> nfft;

	int i,j;
	double complex nu1, nu2;

    /* number of threads (always one if no openmp) */
    int number_of_threads=1;
    /* index of the thread (always 0 if no openmp) */
    int thread=0;

	#ifdef _OPENMP

	#pragma omp parallel
			{
			number_of_threads = omp_get_num_threads();
			}
	#endif

	#pragma omp parallel \
		shared(Nmax, hm_switch, fft_input, matrix)\
		private(i, j, nu1, nu2) \
		num_threads(number_of_threads)

		{
		#pragma omp for schedule (dynamic)
		for (i=0; i<Nmax+1; i++){
			for (j=0; j<Nmax+1; j++){
					if (hm_switch == MATTER){
						nu1 = -0.5 * fft_input->etam_m[i];
						nu2 = -0.5 * fft_input->etam_m[j];
					}
					else {
						nu1 = -0.5 * fft_input->etam_g[i];
						nu2 = -0.5 * fft_input->etam_g[j];
					}
					matrix[i][j] = mat_func(nu1,nu2);
				}
			}
		}
}

void p_mat_fill(double complex (mat_func)(double complex), struct fft_struct *fft_input, long hm_switch, double complex *matrix)
{
	int Nmax = fft_input -> nfft;

	int i;
	double complex nu1;
	    /* number of threads (always one if no openmp) */
    int number_of_threads=1;
    /* index of the thread (always 0 if no openmp) */
    int thread=0;

	#ifdef _OPENMP

	#pragma omp parallel
			{
			number_of_threads = omp_get_num_threads();
			}
	#endif

	#pragma omp parallel \
		shared(Nmax, hm_switch, fft_input, matrix)\
		private(i, nu1) \
		num_threads(number_of_threads)

		{
		#pragma omp for schedule (dynamic)
		for (i=0; i<Nmax+1; i++){
			if (hm_switch == MATTER){
				nu1 = -0.5 * fft_input->etam_m[i];
			}
			else {
				nu1 = -0.5 * fft_input->etam_g[i];
			}

			matrix[i] = mat_func(nu1);
		}
	}
}

void vec_fill(struct fft_struct *fft_input, double k, long hm_switch, double complex *vec)
{
	int Nmax = fft_input -> nfft;

	int i,j;
	for (i=0; i<Nmax+1; i++){
		if (hm_switch == MATTER){
			vec[i] = fft_input->cmsym_m[i] * cpow(k,fft_input->etam_m[i]);
		}
		else {
			vec[i] = fft_input->cmsym_g[i] * cpow(k,fft_input->etam_g[i]);
		}
	}
}

void c_matmul(double complex** matrix, double complex* arr, int rows, int cols, double complex* result)
{
	int i,j;
	for (i=0; i<rows; i++){
		for (j=0; j<cols; j++){
			result[i] += matrix[i][j] * arr[j];
		}
	}
}

void c_matadd(double complex** matrix1, double complex** matrix2, int rows, int cols, double complex** result){
	int i,j;
	for (i=0; i<rows; i++){
		for (j=0; j<cols; j++){
			result[i][j] = matrix1[i][j] + matrix2[i][j];
		}
	}
}

void c_matadd3(double complex** matrix1, double complex** matrix2, double complex** matrix3, int rows, int cols, double complex** result){
	int i,j;
	for (i=0; i<rows; i++){
		for (j=0; j<cols; j++){
			result[i][j] = matrix1[i][j] + matrix2[i][j] + matrix3[i][j];
		}
	}
}

void vecmult(double* arr1, double* arr2, int rows, double *result)
{
	int i;

	for (i=0; i<rows; i++){
		result[i] = arr1[i] * arr2[i];
	}
}

void c_dot(double complex* arr1, double complex* arr2, int rows, double *result)
{
	int i;
	double complex res = 0.;

	for (i=0; i<rows; i++){
		res += arr1[i] * arr2[i];
	}

	*result = creal(res);
}

void dot(double* arr1, double* arr2, int rows, double *result)
{
	int i;

	for (i=0; i<rows; i++){
		*result += arr1[i] * arr2[i];
	}
}

void c_nonprop(double complex* arr1, double complex** matrix, double complex* arr2, int rows, double *result)
{
	double complex* arr3;
	arr3 = make_1D_c_array(rows);
	c_matmul(matrix, arr2, rows, rows, arr3);
	c_dot(arr1, arr3, rows, result);
    free(arr3);
}
