//JL#include <stdlib.h>
//JL#include <stdio.h>
//JL#include <math.h>
//JL#include <string.h>

	int *make_1D_int_array(long size);
	double 	*make_1Darray(long size);
	double complex 	*make_1D_c_array(long size);
	double  **make_2Darray(long nrows, long ncolumns);
	double complex **make_2D_c_array(long nrows, long ncolumns);
	double 	*init_1Darray(long n,double xmin,double xmax);
	double 	*loginit_1Darray(long n,double xmin,double xmax);
	long 	count_lines_in_file(char *fname);
	void 	return_arr(double a, double b, long n,double ** arr1);
	void 	free_2Darray(double ** m);
	void 	free_2D_c_array(double complex ** m);
    long count_cols_in_file(char *fname);

	double *log10init_1Darray(long n, double inc, double xmin);

	void np_mat_fill(double complex (mat_func)(double complex, double complex), struct fft_struct *fft_input, double k, double complex **matrix);
	void np_mu_mat_fill(double complex (mat_func)(double complex, double complex, double), struct fft_struct *fft_input, double k, double mu, double complex **matrix);
	void p_mat_fill(double complex (mat_func)(double complex), struct fft_struct *fft_input, double k, double complex *matrix);
	void p_mu_mat_fill(double complex (mat_func)(double complex, double), struct fft_struct *fft_input, double k, double mu, double complex *matrix);
	void vec_fill(struct fft_struct *fft_input, double k, double complex *vec);

	void c_nonprop(double complex* arr1, double complex** matrix, double complex* arr2, int rows, double complex *result);
	void c_matmul(double complex** matrix, double complex* arr, int rows, int cols, double complex* result);
	void c_dot(double complex* arr1, double complex* arr2, int rows, double complex * result);
	void dot(double* arr1, double* arr2, int rows, double *result);